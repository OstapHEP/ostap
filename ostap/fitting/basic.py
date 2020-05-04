#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostap/fitting/basic.py
#  Set of useful basic utilities to build various fit models 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
# =============================================================================
"""Set of useful basic utilities to build various fit models"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'PDF'           , ## useful base class for 1D-models
    'MASS'          , ## useful base class to create "signal" PDFs for mass-fits
    'RESOLUTION'    , ## useful base class to create "resolution" PDFs
    ##
    'Fit1D'         , ## the basic compound 1D-fit model 
    ##
    'Flat1D'        , ## trivial 1D-pdf: constant 
    'Generic1D_pdf' , ## wrapper over imported RooFit (1D)-pdf
    'Sum1D'         , ## wrapper for RooAddPdf 
    'H1D_pdf'       , ## convertor of 1D-histo to RooHistPdf 
    ##
    )
# =============================================================================
import ROOT, math,  random
import ostap.fitting.roofit 
import ostap.fitting.variables
import ostap.fitting.roocollections 
from   builtins                import range
from   ostap.core.core         import cpp , Ostap , VE , hID , dsID , rootID, valid_pointer
from   ostap.math.base         import iszero , frexp10 
from   ostap.core.ostap_types  import ( is_integer     , string_types   , 
                                        integer_types  , num_types      ,
                                        list_types     , dictlike_types ) 
from   ostap.fitting.roofit    import SETVAR, FIXVAR, PDF_fun
from   ostap.logger.utils      import roo_silent   , rootWarning
from   ostap.fitting.utils     import ( RangeVar   , MakeVar  , numcpu , 
                                        fit_status , cov_qual , H1D_dset , get_i  )
from   ostap.fitting.funbasic  import FUNC 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.basic' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
## @class PDF
#  The helper base class for implementation of various PDF-wrappers 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-08-21
class PDF (FUNC) :
    """Useful helper base class for implementation of various PDF-wrappers 
    """
    def __init__ ( self , name ,  xvar , special = False ) :

        FUNC.__init__  ( self , name , xvar = xvar ) 

        logger.error('I AM PDF-init')

        self.__signals               = ROOT.RooArgList ()
        self.__backgrounds           = ROOT.RooArgList ()
        self.__components            = ROOT.RooArgList ()
        self.__crossterms1           = ROOT.RooArgSet  () 
        self.__crossterms2           = ROOT.RooArgSet  () 
        self.__combined_signals      = ROOT.RooArgList ()
        self.__combined_backgrounds  = ROOT.RooArgList ()
        self.__combined_components   = ROOT.RooArgList ()
        ## take care about sPlots 
        self.__splots          = []
        self.__histo_data      = None
        self.__fit_options     = () ## predefined fit options for this PDF
        self.__special         = True if special else False
        
        self.__alist1     = ROOT.RooArgList()
        self.__alist2     = ROOT.RooArgList()
        self.__alist3     = []

        self.__pdf        = None
        
        self.config = { 'name' : self.name , 'xvar' : self.xvar ,  'special' : self.special }
        
    ## conversion to string 
    def __str__ (  self ) :
        return '%s(%s,xvar=%s)' % ( self.__class__.__name__ , self.name , self.xvar.name )
    __repr__ = __str__ 
    
    @property
    def pdf  ( self ) :
        """The actual PDF (ROOT.RooAbsPdf)"""
        return self.fun 
    @pdf.setter
    def pdf  ( self , value ) :
        if value is None :
            self.fun = value
            return        
        if not self.special :
            assert isinstance ( value , ROOT.RooAbsPdf ) , "``pdf'' is not ROOT.RooAbsPdf"
        self.fun = value
        
    @property
    def pdf_name ( self ) :
        """``pdf_name'' : get the name of the underlying RooAbsPdf"""
        return  self.fun_name 
        
    @property
    def value ( self ) :
        """``value''  :  get the value of PDF"""
        v = float ( self )
        if self.fit_result :
            e = self.pdf.getPropagatedError ( self.fit_result )
            if 0 <= e : return  VE ( v ,  e * e )           ## RETURN
        return  v
    
    @property
    def special ( self ) :
        """``special'' : is this PDF ``special''   (does nor conform some requirements)?"""
        return self.__special
    
    @property
    def title ( self ) :
        """``title'' : get the title for RooAbsPdf"""
        return self.pdf.title if self.pdf else self.name
    
    @property
    def alist1 ( self ) :
        """list/RooArgList of PDF components for compound PDF"""
        return self.__alist1
    @alist1.setter
    def alist1 ( self , value ) :
        assert isinstance ( value , ROOT.RooArgList ) , "Value must be RooArgList, %s/%s is  given" % ( value , type(value) )
        self.__alist1 = value
        
    @property
    def alist2 ( self ) :
        """list/RooArgList of PDF  component's fractions (or yields for exteded fits) for compound PDF"""        
        return self.__alist2
    @alist2.setter
    def alist2 ( self , value ) :
        assert isinstance ( value , ROOT.RooArgList ) , "Value must be RooArgList, %s/%s is  given" % ( value , type(value) )
        self.__alist2 = value                    

    @property
    def alist3 ( self ) :
        """``alist3'' : list of components as PDFs"""
        if len (  self.alist1 ) != len ( self.__alist3 ) : self.error ( "Mismatch for ``alist1/alis3''" )
        return self.__alist3
    @alist3.setter
    def alist3 ( self , value ) :
        nv = [] 
        for v in value :
            assert isinstance ( v , PDF ) , 'Invalid component type %s/%s' % ( v , type ( v ) )
            nv.append ( v )
        if len (  self.alist1 ) != len ( self.__alist3 ) : self.error ( "Mismatch for ``alist1/alis3''" )
        self.__alist2 = nv 
        
    @property
    def signals     ( self ) :
        """The list/ROOT.RooArgList of all ``signal'' components,
        e.g. for visualization"""
        return self.__signals
    @property
    def backgrounds ( self ) :
        """The list/ROOT.RooArgList of all ``background'' components,
        e.g. for visualization"""
        return self.__backgrounds 
    @property
    def components  ( self ) :
        """The list/ROOT.RooArgList of all ``other'' components,
        e.g. for visualization"""
        return self.__components

    @property
    def combined_signals     ( self ) :
        """The list/ROOT.RooArgList of all combined ``signal'' components,
        e.g. for visualization"""
        return self.__combined_signals
    @property
    def combined_backgrounds ( self ) :
        """The list/ROOT.RooArgList of all combined ``background'' components,
        e.g. for visualization"""
        return self.__combined_backgrounds 
    @property
    def combined_components  ( self ) :
        """The list/ROOT.RooArgList of all combined ``other'' components,
        e.g. for visualization"""
        return self.__combined_components

    @property 
    def crossterms1 ( self ) :
        """``cross-terms'': cross-components for multidimensional PDFs e.g.
        - Signal(x)*Background(y)           for 2D-fits,
        - Signal(x)*Signal(y)*Background(z) for 3D-fits, etc...         
        """        
        return self.__crossterms1
    @property
    def crossterms2 ( self ) :
        """``cross-terms'': cross-components for multidimensional PDFs e.g.
        - Signal(y)*Background(x)               for 2D-fits,
        - Signal(x)*Background(y)*Background(z) for 3D-fits, etc...         
        """        
        return self.__crossterms2
    
    @property
    def histo_data  ( self ):
        """Histogram representation as DataSet (RooDataSet)"""
        return self.__histo_data
    
    @histo_data.setter
    def  histo_data ( self  , value ) :
        if   value is None :
            self.__histo_data = value 
        elif hasattr ( value , 'dset' ) and isinstance ( value.dset , ROOT.RooDataHist ) :
            self.__histo_data = value 
        else :
            raise AttributeError("``histo_data'' has invalid type %s/%s" % (   value , type(value) ) )
    @property
    def fit_options ( self ) :
        """``fit_options'' : the predefined ``fitTo''-options for this PDF
        - tuple of ROOT.RooArgCmd
        pdf = ...
        pdf.fit_options = ROOT.RooFit.Optimize ( 1 )
        pdf.fit_options = ROOT.RooFit.Optimize ( 1 ) , ROOT.RooFit.PrintEvalError ( 2 ) 
        """
        return self.__fit_options
    @fit_options.setter
    def fit_options ( self , value )  :
        if isinstance ( value , ROOT.RooCmdArg ) : value = value , 
        assert isinstance ( value , list_types ), 'Invalid fitTo-options %s' % value 
        _opts = []
        for v in value :
            assert isinstance ( v , ROOT.RooCmdArg ), 'Invalid fitTo-option %s' % v
            _opts.append ( v )
        self.__fit_options = tuple ( _opts ) 
            
    # =========================================================================
    ## make the actual fit (and optionally draw it!)
    #  @code
    #  r,f = model.fitTo ( dataset )
    #  r,f = model.fitTo ( dataset , weighted = True )    
    #  r,f = model.fitTo ( dataset , ncpu     = 10   )    
    #  r,f = model.fitTo ( dataset , draw = True , nbins = 300 )    
    #  @endcode 
    def fitTo ( self           ,
                dataset        ,
                draw   = False ,
                nbins  = 100   ,
                silent = False ,
                refit  = False ,
                timer  = False , 
                args   = ()    , **kwargs ) :
        """
        Perform the actual fit (and draw it)
        >>> r,f = model.fitTo ( dataset )
        >>> r,f = model.fitTo ( dataset , weighted = True )    
        >>> r,f = model.fitTo ( dataset , ncpu     = 10   )    
        >>> r,f = model.fitTo ( dataset , draw = True , nbins = 300 )    
        """
        if timer :
            from ostap.utils.timing import timing 
            with timing ( name = "`fitTo'" , logger = self.logger , format =  "Timing %-18s %7.1fs") :
                if not silent : args = args + ( ROOT.RooFit.Timer ( True ) , ) 
                return self.fitTo ( dataset = dataset ,
                                    draw    = draw    ,
                                    nbins   = nbins   ,
                                    silent  = silent  ,
                                    refit   = refit   ,
                                    timer   = False   , ## NB
                                    args    =  args   , **kwargs )
            
        
        if   isinstance ( dataset , H1D_dset ) : dataset = dataset.dset        
        elif isinstance ( dataset , ROOT.TH1 ) :
            density = kwargs.pop ( 'density' , False  )
            chi2    = kwargs.pop ( 'chi2'    , False  )
            return self.fitHisto ( dataset           ,
                                   draw    = draw    ,
                                   silent  = silent  ,
                                   density = density ,
                                   nbins   = nbins   , 
                                   chi2    = chi2    , args = args , **kwargs ) 
        #
        ## treat the arguments properly
        #
        opts = self.fit_options + ( ROOT.RooFit.Save () , ) + args 
        opts = self.parse_args ( dataset , *opts , **kwargs )
        if not silent and opts : self.info ('fitTo options: %s ' % list ( opts ) )

        ## play a bit with the binning cache for convolutions 
        if self.xvar.hasBinning ( 'cache' ) :
            nb1 = self.xvar.getBins( 'cache' ) 
            xv  = getattr ( dataset , self.xvar.name , None )
            if   xv and xv.hasBinning ( 'cache' ) :
                nb2 = xv.getBins('cache')
                if  nb1 != nb2 :
                    xv.setBins ( max (  nb1 , nb2 ) , 'cache' )
                    self.info ('Adjust binning cache %s->%s for variable %s in dataset' % ( nb2 , nb1 , xv.name ) )
            elif xv :
                xv.setBins (        nb1         , 'cache' )
                self    .info ('Set binning cache %s for variable %s in dataset' %  ( nb1 , xv.name )  )
                
        #
        ## define silent context
        with roo_silent ( silent ) :
            self.fit_result = None
            result          = self.pdf.fitTo ( dataset , *opts ) 
            self.fit_result = result 
            if hasattr ( self.pdf , 'setPars' ) : self.pdf.setPars() 

        if not valid_pointer (  result ) :
            self.fatal ( "fitTo: RooFitResult is invalid. Check model&data" )
            self.fit_result = None             
            return None , None
        
        st = result.status()
        ## if 0 != st and silent :
        ##     self.warning ( 'fitTo: status is %s. Refit in non-silent regime ' % fit_status ( st ) )
        ##     return self.fitTo ( dataset ,
        ##                         draw   = draw  ,
        ##                         nbins  = nbins ,
        ##                         silent = False ,
        ##                         refit  = refit ,
        ##                         args   = args  , **kwargs )
        
        for_refit = False
        if 0 != st   :
            for_refit = 'status' 
            self.warning ( 'fitTo: Fit status is %s ' % fit_status ( st ) )
        #
        qual = result.covQual()
        if   -1 == qual and dataset.isWeighted() : pass
        elif  3 != qual :
            for_refit = 'covariance'
            self.warning ( 'fitTo: covQual    is %s ' % cov_qual ( qual ) )

        #
        ## check the integrals (when possible)
        #
        if hasattr ( self , 'yields' ) and self.yields  :
            
            nsum = VE()            
            for i in self.yields:
                nsum += i.value
                if i.minmax() :
                    imn , imx = i.minmax()
                    idx = imx - imn
                    iv  = i.getVal()
                    ie  = i.error if hasattr ( iv  , 'error' ) else 0 
                    if    iv > imx - 0.05 * idx : 
                        self.warning ( "fitTo: variable ``%s'' == %s [very close (>95%%) to maximum %s]"
                                       % ( i.GetName() , i.value , imx ) )
                    elif  0  < ie and iv < imn + 0.1 * ie :
                        self.warning ( "fitTo: variable ``%s'' == %s [very close (<0.1sigma) to minimum %s]"
                                       % ( i.GetName() , i.value , imn ) )                        
                    elif  iv < imn + 0.01 * idx : 
                        self.debug   ( "fitTo: variable ``%s'' == %s [very close (< 1%%) to minimum %s]"
                                       % ( i.GetName() , i.value , imn ) )
                        
            if not dataset.isWeighted () :

                sums = [ nsum ]
                if 2 <= len ( self.yields ) : sums.append ( result.sum ( *self.yields ) )

                for ss in sums :
                    if 0 >= ss.cov2() : continue 
                    nl = ss.value() - 0.50 * ss.error() 
                    nr = ss.value() + 0.50 * ss.error()
                    if not nl <= len ( dataset ) <= nr :
                        self.warning ( 'fitTo: fit is problematic: ``sum'' %s != %s [%+.5g/%+.5g]' % ( ss , len( dataset ) , nl , nr ) )
                        for_refit = 'integral'
        #
        ## call for refit if needed
        #
        if refit and for_refit :
            self.info ( 'fitTo: call for refit:  %s/%s'  % ( for_refit , refit ) ) 
            if   is_integer ( refit ) : refit -= 1
            else                      : refit  = False
            return  self.fitTo ( dataset         ,
                                 draw   = draw   ,
                                 nbins  = nbins  ,
                                 silent = silent ,
                                 refit  = refit  ,
                                 args   = args   , **kwargs ) 

        frame = None
        
        ## draw it if requested
        if draw :  
            from ostap.plotting.fit_draw import draw_options
            draw_opts = draw_options ( **kwargs )
            if draw_opts and not draw     : draw = draw_opts
            if isinstance ( draw , dict ) : draw_opts.update( draw )            
            frame = self.draw ( dataset , nbins = nbins , silent = silent , **draw_opts ) 
                        
        if hasattr ( self.pdf , 'setPars' ) : self.pdf.setPars()
            
        for s in self.components  : 
            if hasattr ( s , 'setPars' ) : s.setPars()
        for s in self.backgrounds :  
            if hasattr ( s , 'setPars' ) : s.setPars() 
        for s in self.signals     : 
            if hasattr ( s , 'setPars' ) : s.setPars() 

        ## 
        return result, frame 

    ## helper method to draw set of components 
    def _draw ( self , what , frame , options , style = None , args = () ) :
        """ Helper method to draw set of components
        """

        from ostap.plotting.fit_draw import Styles, Style
        
        if isinstance ( options , ROOT.RooCmdArg ) : options = options, 
        elif not options                           : options = ()

        if   isinstance ( style , Styles     ) : pass
        elif isinstance ( style , Style      ) : style = Styles ( [ style ] )
        elif isinstance ( style , list_types ) : style = Styles (   style   )   

        if args :
            self.error ( "___DRAW: " + str ( args ) )
        
        for i , cmp in enumerate ( what ) :

            st         = style  ( i ) if style and callable  ( style ) else ()

            cmps       = ROOT.RooArgSet         ( cmp  )             
            components = ROOT.RooFit.Components ( cmps )
            
            command    = ROOT.RooLinkedList()
            command.add ( components )
            
            for s in st         : command.add ( s )
            for o in options    : command.add ( o ) 
            for a in args       : command.add ( a )
 
            self.pdf .plotOn ( frame , command )
            
            ncmps = [ c.GetName() for c in cmps ]
            if 1 == len ( ncmps )  :  ncmps = ncmps[0]
            self.debug ("draw ``%s'' with %s" % ( ncmps , st + options ) )
            
            del command

                                            
    # ================================================================================
    ## draw fit results
    #  @code
    #  r,f = model.fitTo ( dataset )
    #  model.draw ( dataset , nbins = 100 ) 
    #  @endcode
    #  @param dataset  dataset to be drawn 
    #  @param nbins    binning scheme for frame/RooPlot 
    #  @param silent   silent mode ?
    #  @param data_options          drawing options for dataset
    #  @param signal_options        drawing options for `signal'        components    
    #  @param background_options    drawing options for `background'    components 
    #  @param crossterm1_options    drawing options for `crossterm-1'   components 
    #  @param crossterm2_options    drawing options for `crossterm-2'   components 
    #  @param background2D_options  drawing options for `background-2D' components 
    #  @param component_options     drawing options for 'other'         components 
    #  @param fit_options           drawing options for fit curve    
    #  @param signal_style          style(s) for signal components 
    #  @param background_style      style(s) for background components
    #  @param component_style       style(s) for other components
    #  @param crossterm1_style      style(s) for ``crossterm-1''   components
    #  @param crossterm2_style      style(s) for ``crossterm-2''   components
    #  @param background2D_style    style(s) for ``background-2D'' components
    #  @see ostap.plotting.fit_draw
    #
    #  Drawing options can be specified as keyword arguments:
    #  @code
    #  fit_curve = ROOT.RooFit.LineColor ( ROOT.kRed ) , ROOT.RooFit.LineWidth ( 3 )
    #  f = pdf.draw ( ... , total_fit_options = fit_curve  , )
    #  @endcode
    #  When options are not provided explicitly, the options defined in the PDF are looked for:
    #  @code
    #  fit_curve = ROOT.RooFit.LineColor ( ROOT.kRed ) , ROOT.RooFit.LineWidth ( 3 )
    #  pdf.draw_opptions['total_fit_options'] = fit_curve 
    #  f = pdf.draw ( ...)
    #  @endcode
    #  Otherwise the default options,  defined in ostap.plotting.fit_draw module, are used 
    #  @see ostap.plotting.fit_draw
    def draw ( self                         ,
               dataset               = None ,
               nbins                 = 100  ,   ## Frame binning
               silent                = True ,   ## silent mode ?
               style                 = None ,   ## use another style ?
               drawvar               = None ,   ## drawvar 
               args                  = ()   , 
               **kwargs                     ) :
        """  Visualize the fits results
        >>> r,f = model.draw ( dataset )
        >>> model.draw ( dataset , nbins = 100 )
        >>> model.draw ( dataset , base_signal_color  = ROOT.kGreen+2 )
        >>> model.draw ( dataset , data_options = (ROOT.RooFit.DrawOptions('zp'),) )

        Produce also residual & pull:
        
        >>> f,r,h = model.draw ( dataset , nbins = 100 , residual = 'P' , pull = 'P')
        
        Drawing options:
        - data_options            ## drawing options for dataset  
        - background_options      ## drawing options for background    component(s)
        - crossterm1_options      ## drawing options for crossterm1    component(s)
        - crossterm2_options      ## drawing options for crossterm2    component(s)
        - signal_options          ## drawing options for signal        component(s)
        - component_options       ## drawing options for other         component(s)
        - background2D_options    ## drawing options for 2D-background component(s)
        - total_fit_options       ## drawing options for the total fit curve
        
        Style&Colors:                  
        - background_style        ## style(s) for background component(s)
        - crossterm1_style        ## style(s) for crossterm1 component(s)
        - crossterm2_style        ## style(s) for crossterm2 component(s)
        - signal_style            ## style(s) for signal     component(s)
        - component_style         ## style(s) for other      component(s)
        - background2D_style      ## style(s) for background component(s)

        Other options:
        -  residual               ## make also residual frame
        -  pull                   ## make also residual frame
        
        For default values see ostap.plotting.fit_draw
        
        - Drawing options can be specified as keyword arguments:
        
        >>> fit_curve = ROOT.RooFit.LineColor ( ROOT.kRed ) , ROOT.RooFit.LineWidth ( 3 )
        >>> f = pdf.draw ( ... , total_fit_options = fit_curve  , )
        
        - when options are not provided explicitly, the options defined in the PDF are looked for:
        
        >>> fit_curve = ROOT.RooFit.LineColor ( ROOT.kRed ) , ROOT.RooFit.LineWidth ( 3 )
        >>> pdf.draw_options['total_fit_options'] = fit_curve 
        >>> f = pdf.draw ( ...)
        
        - otherwise the default options, defined in ostap.plotting.fit_draw module, are used 

        """
        #
        
        from ostap.plotting.style import useStyle 

        #
        ## again the context
        # 
        with roo_silent ( silent ) , useStyle ( style ) :

            drawvar = drawvar if drawvar else ( self.draw_var if self.draw_var else self.xvar )  

            binned = dataset and isinstance ( dataset , ROOT.RooDataHist )

            if nbins :  frame = drawvar.frame ( nbins )            
            else     :  frame = drawvar.frame ()

            #
            ## draw invizible data (for normalzation of fitting curves)
            #
            data_options = self.draw_option ( 'data_options' , **kwargs )
            kwargs.pop ( 'data_options' , () )
            if dataset and dataset.isWeighted() and dataset.isNonPoissonWeighted() : 
                data_options = data_options + ( ROOT.RooFit.DataError( ROOT.RooAbsData.SumW2 ) , )
                
            if dataset and binned and nbins :
                data_options = data_options + ( ROOT.RooFit.Binning ( nbins ) , )

            if dataset :
                command  = ROOT.RooLinkedList()
                for o in data_options : command.add ( o )
                for a in args         : command.add ( o )
                invisible = ROOT.RooFit.Invisible()  
                command.add ( invisible ) 
                dataset .plotOn ( frame , command )
                
            ## draw various ``background'' terms
            boptions     = self.draw_option ( 'background_options' , **kwargs ) 
            bbstyle      = self.draw_option (   'background_style' , **kwargs )
            self._draw( self.backgrounds , frame , boptions , bbstyle )
            kwargs.pop ( 'background_options' , () )
            kwargs.pop ( 'background_style'   , () )

            ## draw combined ``background'' components 
            if self.combined_backgrounds :
                drawit   = self.draw_option ( 'draw_combined_background'    , **kwargs )
                doptions = self.draw_option ( 'combined_background_options' , **kwargs ) 
                dstyle   = self.draw_option ( 'combined_background_style'   , **kwargs )
                if drawit : self._draw ( self.combined_backgrounds , frame , doptions , dstyle , args )
                
            kwargs.pop ( 'combined_background_options' , ()   )
            kwargs.pop ( 'combined_background_style'   , ()   )
            kwargs.pop ( 'draw_combined_backgrounds'   , True )
            
            ## ugly :-(
            ct1options   = self.draw_option ( 'crossterm1_options' , **kwargs )
            ct1bstyle    = self.draw_option ( 'crossterm1_style'   , **kwargs ) 
            if hasattr ( self , 'crossterms1' ) and self.crossterms1 : 
                self._draw( self.crossterms1 , frame , ct1options , ct1bstyle , args )
            kwargs.pop ( 'crossterm1_options' , () )
            kwargs.pop ( 'crossterm1_style' , () )

            ## ugly :-(
            ct2options   = self.draw_option ( 'crossterm2_options' , **kwargs )
            ct2bstyle    = self.draw_option ( 'crossterm2_style'   , **kwargs ) 
            if hasattr ( self , 'crossterms2' ) and self.crossterms2 :
                self._draw( self.crossterms2 , frame , ct2options , ct2bstyle , args )
            kwargs.pop ( 'crossterm2_options' , () )
            kwargs.pop ( 'crossterm2_style'   , () )

            ## draw ``other'' components
            coptions     = self.draw_option ( 'component_options' , **kwargs )
            cbstyle      = self.draw_option ( 'component_style'   , **kwargs )
            self._draw( self.components , frame , coptions , cbstyle , args )
            kwargs.pop ( 'component_options' , () )
            kwargs.pop ( 'component_style'   , () )

            ## draw combined ``other'' components 
            if self.combined_components :
                drawit   = self.draw_option ( 'draw_combined_component'    , **kwargs )
                doptions = self.draw_option ( 'combined_component_options' , **kwargs ) 
                dstyle   = self.draw_option ( 'combined_component_style'   , **kwargs )
                if drawit : self._draw ( self.combined_components , frame , doptions , dstyle , args )
                
            kwargs.pop ( 'combined_component_options' , ()   )
            kwargs.pop ( 'combined_component_style'   , ()   )
            kwargs.pop ( 'draw_combined_component'    , True )
                    
            ## draw ``signal'' components
            soptions     = self.draw_option (    'signal_options'  , **kwargs )
            sbstyle      = self.draw_option (      'signal_style'  , **kwargs ) 
            self._draw( self.signals , frame , soptions , sbstyle , args )
            kwargs.pop ( 'signal_options' , () )
            kwargs.pop ( 'signal_style'   , () )

            ## draw combined ``signals'' components 
            if self.combined_signals :
                drawit    = self.draw_option ( 'draw_combined_signal'     , **kwargs )
                doptions  = self.draw_option ( 'combined_signal_options'  , **kwargs ) 
                dstyle    = self.draw_option (   'combined_signal_style'  , **kwargs )
                if drawit : self._draw ( self.combined_signals , frame , doptions , dstyle , args )
                
            kwargs.pop ( 'combined_signal_options' , ()   )
            kwargs.pop ( 'combined_signal_style'   , ()   )
            kwargs.pop ( 'draw_combined_signal'    , True )
            
            #
            ## the total fit curve
            #
            totoptions   = self.draw_option (  'total_fit_options' , **kwargs )
            self.pdf .plotOn ( frame , *totoptions )
            kwargs.pop ( 'total_fit_options' , () )            
            #
            ## draw data once more
            #
            if dataset :
                command    = ROOT.RooLinkedList()
                for o in data_options : command.add ( o )
                for a in args         : command.add ( a ) 
                dataset .plotOn ( frame , command )
                del command

            #
            ## suppress ugly axis labels
            #
            if not kwargs.get ( 'draw_axis_title' , False ) :  
                frame.SetXTitle ( '' )
                frame.SetYTitle ( '' )
                frame.SetZTitle ( '' )
                
            #
            ## Draw the frame!
            #
            if not ROOT.gROOT.IsBatch() :
                with rootWarning (): frame.draw( kwargs.pop ( 'draw_options','' ) )
            
            residual =  kwargs.pop ( 'residual' , False )
            if residual and not  dataset :
                self.warning("draw: can't produce residual without data")
                residual = False
                
            pull     =  kwargs.pop ( 'pull'     , False ) 
            if pull     and not  dataset :
                self.warning("draw: can't produce residual without data")
                residual = False

            if kwargs :
                self.warning("draw: ignored unknown options: %s" % list( kwargs.keys() ) ) 

            ## calculate chi2/ndf
            frame.chi2dnf = None 
            if dataset and not silent :             
                pars          = self.pdf.getParameters ( dataset )
                frame.chi2ndf = frame.chiSquare ( len ( pars ) )
                binw          = -1 
                if nbins and isinstance ( nbins , integer_types ) and 1 < nbins :
                    if hasattr ( drawvar , 'xminmax' ) and drawvar.xminmax () :
                        xmn , xmx =  drawvar.xminmax()
                        binw = ( xmx - xmn ) / float ( nbins )
                if 0 < binw : self.info ( 'chi2/ndf: %s, binwidth: %s' %  ( frame.chi2ndf , binw ) )
                else        : self.info ( 'chi2/ndf: %s' %                  frame.chi2ndf          )
                
            if not residual and not pull:
                return frame

            rframe =  None 
            if residual  :
                if   residual is True               : residual =      "P" ,
                elif isinstance  ( residual , str ) : residual = residual ,
                rframe  = frame.emptyClone ( rootID ( 'residual_' ) )
                rh      = frame.residHist()
                rframe.addPlotable ( rh , *residual ) 
                if not kwargs.get( 'draw_axis_title' , False ) :  
                    rframe.SetXTitle ( '' )
                    rframe.SetYTitle ( '' )
                    rframe.SetZTitle ( '' )
                    
            pframe = None 
            if pull      : 
                if   pull is True               : pull =   "P",
                elif isinstance  ( pull , str ) : pull = pull ,
                pframe  = frame.emptyClone ( rootID ( 'pull_' ) )
                ph      = frame.pullHist()
                pframe.addPlotable ( ph , *pull ) 
                if not kwargs.get( 'draw_axis_title' , False ) :  
                    pframe.SetXTitle ( '' )
                    pframe.SetYTitle ( '' )
                    pframe.SetZTitle ( '' )

            return frame, rframe, pframe  

    # =========================================================================
    ## fit the 1D-histogram (and draw it)
    #  @code
    #  histo = ...
    #  r,f = model.fitHisto ( histo , draw = True ) 
    #  @endcode 
    def fitHisto ( self            ,
                   histo           ,
                   draw    = False ,
                   silent  = False ,
                   density = False ,
                   chi2    = False ,
                   nbins   = None  , 
                   args    = () , **kwargs ) :
        """Fit the histogram (and draw it)

        >>> histo = ...
        >>> r,f = model.fitHisto ( histo , draw = True ) 
        
        """
        with RangeVar( self.xvar , *(histo.xminmax()) ) : 

            hdata = getattr ( self , 'histo_data' , None )
            if hdata and isinstance ( hdata , H1D_dset ) and \
                   hdata.histo      is histo             and \
                   hdata.density    == density           and \
                   hdata.histo_hash == hash ( histo ) :
                ## reuse the existing dataset
                self.debug ('Reuse the existing H1D_dset') 
                data = hdata.dset
            else :                
                ## convert it! 
                self.debug ('Create new H1D_dset'        ) 
                self.histo_data = H1D_dset ( histo , self.xvar , density , silent )
                data            = self.histo_data.dset

            if  self.xminmax() :
                
                mn  , mx  = self .xminmax ()
                hmn , hmx = histo.xminmax ()

                vmn = max ( mn , hmn )
                vmx = min ( mx , hmx )

                args = list  ( args )
                ## args.append  ( ROOT.RooFit.Range ( vmn , vmx ) )  
                args = tuple ( args )
                
            if chi2 : return self.chi2fitTo ( data               ,
                                              draw    = draw     ,
                                              silent  = silent   ,
                                              density = density  ,
                                              nbins   = nbins    ,
                                              args    = args     , **kwargs )
            else    : return self.fitTo     ( data               ,
                                              draw    = draw     ,
                                              nbins   = nbins    , 
                                              silent  = silent   ,
                                              args    = args     , **kwargs )

    # =========================================================================
    ## make chi2-fit for binned dataset or histogram
    #  @code
    #  histo = ...
    #  r,f = model.chi2FitTo ( histo , draw = True ) 
    #  @endcode
    #  @todo add proper parsing of arguments for RooChi2Var 
    def chi2fitTo ( self            ,
                    dataset         ,
                    draw    = False ,
                    silent  = False ,
                    density = False ,
                    nbins   = None  , 
                    args    = ()    , **kwargs ) :
        """ Chi2-fit for binned dataset or histogram
        >>> histo = ...
        >>> result , frame  = model.chi2FitTo ( histo , draw = True ) 
        """
        hdataset = dataset
        histo    = None 
        if isinstance  ( dataset , ROOT.TH1 ) :
            # if histogram, convert it to RooDataHist object:
            xminmax = dataset.xminmax() 
            with RangeVar( self.xvar , *xminmax ) :                
                self.histo_data = H1D_dset ( dataset , self.xvar , density , silent )
                hdataset        = self.histo_data.dset 
                histo           = dataset 
                
        with roo_silent ( silent ) : 

            
            lst1 = self.fit_options + ( ROOT.RooFit.Save () , ) + args 
            lst1 = list ( self.parse_args ( hdataset , *args , **kwargs ) )
            lst2 = []
            
            if       self.pdf.mustBeExtended () : lst2.append ( ROOT.RooFit.Extended ( True  ) )
            elif not self.pdf.canBeExtended  () : lst2.append ( ROOT.RooFit.Extended ( False ) )
            
            if not silent : lst2.append ( ROOT.RooFit.Verbose  () )
            if histo :
                if histo.natural() : lst2.append ( ROOT.RooFit.DataError ( ROOT.RooAbsData.Poisson ) )
                else               : lst2.append ( ROOT.RooFit.DataError ( ROOT.RooAbsData.SumW2   ) )  

            args_ = tuple ( lst2 + lst1  )
            #
            chi2 = ROOT.RooChi2Var ( rootID ( "chi2_" ) , "chi2(%s)" % self.name  , self.pdf , hdataset , *args_ )
            m    = ROOT.RooMinuit  ( chi2 ) 
            m.migrad   () 
            m.hesse    ()
            result = m.save ()
            ## save fit results 
            self.fit_result = result 

        if not draw :
            return result, None 
        
        from ostap.plotting.fit_draw import draw_options         
        draw_opts = draw_options ( **kwargs )
        if isinstance ( draw , dict ) : draw_opts.update( draw )

        return result, self.draw ( hdataset , nbins = nbins , silent = silent , **draw_opts )




    # =========================================================================
    ## draw/prepare NLL or LL-profiles for selected variable
    #  @code
    #  model.fitTo ( dataset , ... )
    #  nll  , f1 = model.draw_nll ( 'B' ,  dataset )
    #  prof , f2 = model.draw_nll ( 'B' ,  dataset , profile = True )
    #  @endcode    
    def draw_nll ( self            ,
                   var             ,
                   dataset         ,
                   profile = False ,
                   draw    = True  ,
                   silent  = True  , 
                   args    = ()    , **kwargs ) :
        """Draw/prepare NLL or LL-profile for seleted variable:        
        >>> model.fitTo ( dataset , ... )
        >>> nll  , f1 = model.draw_nll ( 'B' ,  dataset )
        >>> prof , f2 = model.draw_nll ( 'B' ,  dataset , profile = True )
        """

        # if histogram, convert it to RooDataHist object:
        if isinstance  ( dataset , ROOT.TH1 ) :
            # if histogram, convert it to RooDataHist object:
            xminmax = dataset.xminmax() 
            with RangeVar( self.xvar , *xminmax ) :
                density = kwargs.pop ( 'density' , False )
                silent  = kwargs.pop ( 'silent'  , True  )                
                self.histo_data   = H1D_dset ( dataset , self.xvar , density , silent )
                hdataset          = self.histo_data.dset
                kwargs [ 'ncpu' ] = 2   
                return self.draw_nll ( var     = var      ,
                                       dataset = hdataset ,
                                       profile = profile  ,
                                       draw    = draw     ,
                                       args    = args     , **kwargs )
            
        ## convert if needed 
        if not isinstance ( dataset , ROOT.RooAbsData ) and hasattr ( dataset , 'dset' ) :
            dataset = dataset.dset
            
        ## get all parametrs
        pars = self.pdf.getParameters ( dataset ) 
        assert var in pars , "Variable %s is not a parameter"   % var
        if not isinstance ( var , ROOT.RooAbsReal ) : var = pars[ var ]
        del pars 
        ##
        fargs = []
        ##
        bins  = kwargs.pop ( 'nbins' , 25 if profile else 200 )
        if bins   : fargs.append ( ROOT.RooFit.Bins      ( bins  ) )
        ## 
        rng   = kwargs.pop ( 'range' , None )
        if rng    : fargs.append ( ROOT.RooFit.Range     ( *rng  ) ) 
        ##
        fargs = tuple ( fargs )
        ##
        largs = [ i for i in  args ]
        ## 
        from ostap.utils.cidict import select_keys
        dopts = select_keys ( kwargs , ( 'line_color' , 'line_width' , 'line_style' , 
                                         'color'      , 'width'      , 'style'      ) , 
                              transform = lambda s : s.lower().replace('_','') )
        for color in  ( 'color' , 'line_color' ) : 
            if color in dopts : 
                largs.append ( ROOT.RooFit.LineColor ( dopts.pop ( color ) ) )                 
        for width in  ( 'width' , 'line_width' ) : 
            if width in dopts : 
                largs.append ( ROOT.RooFit.LineWidth ( dopts.pop ( width ) ) )
        for style in  ( 'style' , 'line_style' ) : 
            if style in dopts : 
                largs.append ( ROOT.RooFit.LineStyle ( dopts.pop ( style ) ) )
        ##
        largs.append  ( ROOT.RooFit.ShiftToZero() ) 
        largs  = tuple ( largs ) 
        
        ## create NLL 
        nll , sf = self.nll ( dataset , silent = silent , **kwargs )  

        result  = nll

        ## make profile? 
        if profile :
            avar    = ROOT.RooArgSet    (  var ) 
            profile = nll.createProfile ( avar )
            result  = profile

        from ostap.fitting.variables import KeepBinning
        
        with KeepBinning ( var ) :

            if bins :
                var.bins = bins
            
            self.debug ( 'draw_nll: frame  args: %s'% list ( fargs ) )
            ## prepare the frame & plot 
            frame = var.frame ( *fargs )
            
            self.debug ( 'draw_nll: plotOn args: %s'% list ( largs ) )
            result.plotOn ( frame , *largs  )
            
            import ostap.histos.graphs
            
            ## remove a bit strange drawing artefacts (the first and the last points )
            if var.minmax() :
                vmn , vmx = var.minmax()
                graph   = frame.getObject ( 0 )
                if graph.remove ( remove = lambda x,y : not vmn <= x <= vmx ) :
                    logger.debug ('draw_nll: remove drawing artefacts  at the first and the last points' ) 
                                        
            ## scale it if needed
            if 1 != sf :
                logger.info ('draw_nll: apply scale factor of %s due to dataset weights' % sf )
                graph  = frame.getObject ( 0 )
                graph *= sf 
                
            gr      = frame.getObject ( 0 )
            mn , mx = gr.minmax()
            m  , e  = frexp10 ( mx )
            m *= 10
            e -= 1
            mx  = math.floor ( m ) + 1
            mx *= 10**e 
            
            frame.SetMinimum ( 0  )
            frame.SetMaximum ( mx )
            
            if not kwargs.get('draw_axis_title' , False ) : 
                frame.SetXTitle  ( '' )
                frame.SetYTitle  ( '' )
                frame.SetZTitle  ( '' )
                
        

        ## draw it! 
        if not ROOT.gROOT.IsBatch() :
            with rootWarning ():
                if draw : frame.draw ( kwargs.get('draw_options', '' ) )

        return result , frame
            
    # =========================================================================
    ## create NLL
    #  @code
    #  model.fitTo ( dataset , ... )
    #  nll, sfactor  = model.nll ( dataset )
    #  @endcode
    #  @see RooAbsPdf::createNLL
    #  @attention Due to the bug/typo in<c> RooAbsPdf.createNLL</c>, line 817 
    #  <c>CloneData</c> depends on <c>Optimize</c>
    #  @todo report problem to RooFit and fix it! 
    def nll ( self            ,
              dataset         ,
              silent  = True  ,
              args    = ()    , **kwargs ) :
        """Create NLL object from the pdf
        >>> model.fitTo ( dataset , ... )
        >>> nll, sf = model.nll ( dataset )
        - see RooAbsPdf::createNLL 
        """
        
        ## convert if needed 
        if not isinstance ( dataset , ROOT.RooAbsData ) and hasattr ( dataset , 'dset' ) :
            dataset = dataset.dset 

        clone = kwargs.pop ( 'clone' , False )
        kwargs [ 'clone' ] = clone 
        if not clone : kwargs['optimize'] = False 

        opts = self.parse_args ( dataset , *args , **kwargs )

        ## skip some artifacts from MakeVars.parse_args 
        ok = []
        for o in opts :
            if o.name == 'SumW2Error'      : continue
            if o.name == 'AsymptoticError' : continue
            if o.name == 'PrintLevel'      : continue
            if o.name == 'PrintEvalErrors' : continue
            ok.append ( o )
            
        opts = tuple ( ok )        
        if not silent and opts : self.info ('NLL options: %s ' % list ( opts ) )
        
        ## get s-Factor 
        sf   = dataset.sFactor() 

        self.debug ( 'nll: createNLL args: %s'% list ( opts ) )            
        return self.pdf.createNLL ( dataset , *opts ) , sf 

    # =========================================================================
    ## get NLL/profile-graph for the variable, using the specified bscissas
    #  @code
    #  pdf   = ...
    #  graph = pdf.graph_nll ( 'S'               ,
    #                          range ( 0 , 100 ) ,
    #                          dataset           )
    #  @endcode
    def graph_nll ( self            ,
                    variable        , 
                    values          ,
                    dataset         ,
                    silent  = True  ,
                    args    = ()    , **kwargs ) :
        """Get NLL/profile-graph for the variable, using the specified abscissas
        >>> pdf   = ...
        >>> graph = pdf.graph_nll ( 'S'               ,
        ...                          range ( 0 , 100 ) ,
        ...                          dataset           )
        """

        ## 1) create NLL 
        nLL, sf = self.nll ( dataset , silent = silent ,  args = args , **kwargs )

        ## get the parametrs
        var  = variable 
        pars = self.pdf.getParameters ( dataset ) 
        assert var in pars , "Variable %s is not a parameter"   % var
        if not isinstance ( var , ROOT.RooAbsReal ) : var = pars[ var ]
        del pars 

        ## 2) collect NLL values 
        results = []
        vmin    = None 
        with SETVAR  ( var ) :
            from ostap.utils.progress_bar import progress_bar 
            for v in progress_bar  ( values , silent = silent ) :
                var.setVal ( v )
                n   = nLL.getVal() 
                res = v , n
                results.append ( res )
                vmin = n if vmin is None else min ( vmin , n ) 
        
        ## 3) create graph
        import ostap.histos.graphs
        graph = ROOT.TGraph ( len ( results ) )
        results.sort () 
        for i , point in enumerate ( results ) :
            x , y = point 
            if vmin is None : graph [ i ] = x , y
            else            : graph [ i ] = x , y - vmin 
            
        ## scale it if needed
        if 1 != sf :
            logger.info ('graph_nll: apply scale factor of %s due to dataset weights' % sf )
            graph *= sf 
            
        return graph 

    # =========================================================================
    ## get NLL-profile-graph for the variable, using the specified abscissas
    #  @code
    #  pdf   = ...
    #  graph = pdf.graph_profile ( 'S'                       ,
    #                              vrange ( 0 , 12.5 , 10  ) ,
    #                              dataset                   )
    #  @endcode
    def graph_profile ( self            ,
                        variable        , 
                        values          ,
                        dataset         ,
                        fix     = []    ,
                        silent  = True  ,
                        args    = ()    , **kwargs ) :
        """Get profile-graph for the variable, using the specified abscissas
        >>> pdf   = ...
        >>> graph = pdf.graph_profile ( 'S'                     ,
        ...                             range ( 0 , 12.5 , 20 ) ,
        ...                             dataset                 )
        """

        ## 1) create NLL 
        nLL , sf = self.nll ( dataset , silent = silent ,  args = args , **kwargs )

        ## get the parametrs
        var  = variable 
        pars = self.pdf.getParameters ( dataset ) 
        assert var in pars , "Variable %s is not a parameter"   % var
        if not isinstance ( var , ROOT.RooAbsReal ) : var = pars[ var ]

        vars = ROOT.RooArgSet ( var )
        for f in fix :
            fv = f
            if not isinstance ( fv , ROOT.RooAbsReal ) : fv = pars [ fv ]
            vars.add ( fv ) 
            
        ## 2)  create profile LL
        pLL = ROOT.RooProfileLL ( 'pLL_%s_%s'         % ( var.name , self.name ) ,
                                  'LL-profile(%s,%s)' % ( var.name , self.name ) ,
                                  nLL , vars )

        ## 2) collect pLL values 
        results = [] 
        vmin    = None 
        with SETVAR  ( var ) :
            from ostap.utils.progress_bar import progress_bar 
            for  v in progress_bar ( values , silent = silent )  :
                var.setVal ( v )
                p   = pLL.getVal() 
                res = v , p 
                results.append ( res )
                vmin = p if vmin is None else min ( vmin , p ) 
             
        ## 3) create graph 
        import ostap.histos.graphs
        graph = ROOT.TGraph ( len ( results ) )
        results.sort ()
        for i , point  in enumerate ( results ) :
            x , y = point 
            if vmin is None : graph [ i ] = x , y
            else            : graph [ i ] = x , y - vmin 
            
        ## scale it if needed
        if 1 != sf :
            logger.info ('graph_profile: apply scale factor of %s due to dataset weights' % sf )
            graph *= sf 
            
        return graph 
        
    # ========================================================================
    ## evaluate "significance" using Wilks' theorem via NLL
    #  @code
    #  data = ...
    #  pdf  = ...
    #  pdf.fitTo ( data , ... )
    #  sigmas = pdf.wilks ( 'S' , data )
    #  @endcode
    def wilks ( self                     ,
                var                      ,
                dataset                  ,
                range    = ( 0 , None )  ,
                silent   = True          ,
                args     = () , **kwargs ) :
        """Evaluate ``significance'' using Wilks' theorem via NLL
        >>> data = ...
        >>> pdf  = ...
        >>> pdf.fitTo ( data , ... )
        >>> sigmas = pdf.wilks ( 'S' , data )
        """
        # if histogram, convert it to RooDataHist object:
        if isinstance  ( dataset , ROOT.TH1 ) :
            # if histogram, convert it to RooDataHist object:
            xminmax = dataset.xminmax() 
            with RangeVar( self.xvar , *xminmax ) :
                density = kwargs.pop ( 'density' , False )
                silent  = kwargs.pop ( 'silent'  , True  )                
                self.histo_data = H1D_dset ( dataset , self.xvar , density , silent )
                hdataset        = self.histo_data.dset
                kwargs['ncpu']  = 1 
                return self.wilks ( var     = var      ,
                                    dataset = hdataset ,
                                    range   = range    ,
                                    silent  = silent   , 
                                    args    = args     , **kwargs )
        ## convert if needed 
        if not isinstance ( dataset , ROOT.RooAbsData ) and hasattr ( dataset , 'dset' ) :
            dataset = dataset.dset 
                          
        ## get all parameters
        pars = self.pdf.getParameters ( dataset ) 
        assert var in pars , "Variable %s is not a parameter"   % var
        if not isinstance ( var , ROOT.RooAbsReal ) : var = pars[ var ]
        del pars
        
        ## unpack the range 
        minv , maxv = range        
        if maxv is None : maxv = var.value

        error = 0 
        if isinstance ( maxv , VE ) :
            if 0 < maxv.cov2 () : error = maxv.error() 
            maxv = maxv.value ()
            
        with roo_silent ( silent ) :
            
            nll , sf = self.nll ( dataset         ,
                                  silent = silent ,
                                  clone  = False  ,
                                  args   = args   , **kwargs ) 
            

            vv = var.getVal()
            mn = min ( minv , maxv - error )
            with RangeVar ( var , mn , maxv + error ) :

                with SETVAR ( var ) : 
                    var.setVal ( minv )
                    val_minv = nll.getVal ()

                with SETVAR ( var ) : 
                    var.setVal ( maxv )
                    val_maxv = nll.getVal ()
                    
                dnll     = val_minv -  val_maxv
                
                if 0 < error :

                    with SETVAR ( var ) : 
                        var.setVal ( maxv + error ) 
                        val_maxvp = nll.getVal()

                    with SETVAR ( var ) :  
                        var.setVal ( maxv - error )
                        val_maxvm = nll.getVal()
                        
                    dnll = VE ( dnll , 0.25 * (val_maxvp - val_maxvm )**2 )
                    
            ## apply scale factor
            if 1 != sf :  logger.info ('Scale factor of %s is applied' % sf )
            dnll *= sf            
                
            ## convert the difference in likelihoods into sigmas 
            result = 2.0 * abs ( dnll )
            result = result**0.5

            del nll
            
        return result if 0<=dnll else -1*result 

    # ========================================================================
    ## evaluate "significance" using Wilks' theorem via NLL
    #  @code
    #  data = ...
    #  pdf  = ...
    #  pdf.fitTo ( data , ... )
    #  sigmas = pdf.wilks2 ( 'S' , data , fix = [ 'mean' , 'gamma' ] )
    #  @endcode
    def wilks2 ( self                           ,
                 var                            ,
                 dataset                        ,
                 fix                            , ## variables to fix 
                 range          = ( 0 , None )  ,
                 silent         = True          ,
                 args           = () , **kwargs ) :
        """Evaluate ``significance'' using Wilks' theorem via NLL
        >>> data = ...
        >>> pdf  = ...
        >>> pdf.fitTo ( data , ... )
        >>> sigmas = pdf.wilks2 ( 'S' , data , fix = [ 'mean' , 'gamma'] )
        """
        # if histogram, convert it to RooDataHist object:
        if isinstance  ( dataset , ROOT.TH1 ) :
            # if histogram, convert it to RooDataHist object:
            xminmax = dataset.xminmax() 
            with RangeVar( self.xvar , *xminmax ) :
                density = kwargs.pop ( 'density' , False )
                silent  = kwargs.pop ( 'silent'  , True  )                
                self.histo_data = H1D_dset ( dataset , self.xvar , density , silent )
                hdataset        = self.histo_data.dset
                kwargs['ncpu']  = 1 
                return self.wilks2 ( var            = var             ,
                                     dataset        = hdataset        ,
                                     fix            = fix             ,
                                     range          = range           , 
                                     silent         = silent          ,
                                     args           = args , **kwargs )
        ## convert if needed 
        if not isinstance ( dataset , ROOT.RooAbsData ) and hasattr ( dataset , 'dset' ) :
            dataset = dataset.dset 
                          
        ## get all parameters
        pars = self.pdf.getParameters ( dataset )
        
        assert var in pars , "Variable %s is not a parameter/1"   % var
        if not isinstance ( var , ROOT.RooAbsReal ) : var = pars[ var ]
        
        if   isinstance ( fix , string_types    ) : fix = [ fix ] 
        elif isinstance ( fix , ROOT.RooAbsReal ) : fix = [ fix ] 
        
        fixed = []
        for f in fix :
            assert f in pars , "Variable %s is not a parameter/2"   % f            
            fixed.append ( f if isinstance ( f , ROOT.RooAbsReal ) else pars [ f ] )

        del pars

        logger.info ( "Wilks: fixed variables: %s" % [ f.GetName()  for f in fixed] )

        ## unpack the range 
        minv , maxv = range        
        if maxv is None : maxv = var.value
        
        error = 0 
        if isinstance ( maxv , VE ) :
            if 0 < maxv.cov2 () : error = maxv.error() 
            maxv = maxv.value ()

        vname = var.GetName() 
        with roo_silent ( silent ) :

            ## fix "fixed" variables and redefine range for main variable
            mn = min ( minv , maxv - error )
            with SETVAR ( var ) , RangeVar ( var , mn , maxv + error ) :

                ## create NLL 
                nLL , sf = self.nll ( dataset         ,
                                      silent = silent ,
                                      clone  = False  ,
                                      offset = False  , 
                                      args   = args   , **kwargs )

                ## create profile
                vars = ROOT.RooArgSet()
                vars.add ( var )
                for f in  fixed : vars.add ( f )
                
                pLL = ROOT.RooProfileLL ( 'pLL_%s'          % self.name ,
                                          'LL-profile (%s)' % self.name , nLL , vars )
                
                ## 1st value:
                var.setVal ( maxv )
                nll_max = pLL.getVal()

                if 0 < error :
                    
                    var.setVal ( maxv + error )
                    ve1 = pLL.getVal()
                    
                    var.setVal ( maxv - error )
                    ve2 = pLL.getVal()

                    nll_max = VE ( nll_max , 0.25 * ( ve1 - ve2 ) ** 2 ) 
                    
                ## 2nd value:
                var.setVal ( minv )
                nll_min = pLL.getVal() 
                
            dnll = nll_min - nll_max

            ## apply scale factor
            if 1 != sf :  logger.info ('Scale factor of %s is applied' % sf )
            dnll *= sf            
        
            ## convert the difference in likelihoods into sigmas/significance
            result = 2.0 * abs ( dnll )
            result = result ** 0.5
            
            del pLL , nLL
            
        return result if 0 <= dnll else -1 * result 
                
    # ========================================================================
    ## get the actual minimizer for the explicit manipulations
    #  @code
    #  data = ...
    #  pdf  = ...
    #  m    = pdf.minuit  ( data )
    #  m.migrad()
    #  m.hesse ()
    #  m.minos ( param )
    #  @endcode
    #  @see RooMinimizer
    def minuit ( self                   ,
                 dataset         = None ,
                 max_calls       = -1   ,
                 max_iterations  = -1   , 
                 opt_const       = True , ## optimize const 
                 strategy        = None ,
                 nLL             = None , ## nLL  
                 args            =   () , **kwargs  ):
        """Get the actual minimizer for the explicit manipulations
        >>> data = ...
        >>> pdf  = ...
        >>> m    = pdf.minuit ( data )
        >>> m.migrad()
        >>> m.hesse ()
        >>> m.minos ( param )
        - see ROOT.RooMinimizer
        """

        ## parse the arguments 
        ## opts = self.parse_args    ( dataset ,
        ##                            ROOT.RooFit.Offset    ( True  ) ,
        ##                            ROOT.RooFit.CloneData ( False ) , *args , **kwargs )
        ## nll  = self.pdf.createNLL ( dataset , *opts )

        assert nLL or dataset , "minuit: nLL or dataset *must* be specified!" 
        if not nLL : 
            nLL , _ = self.nll ( dataset , args = args , **kwargs )
        
        m    = ROOT.RooMinimizer ( nLL )
        if isinstance  ( opt_const  , bool ) : m.optimizeConst ( opt_const ) 
        if isinstance  ( max_calls      , integer_types ) and 1 < max_calls :
            m.setMaxFunctionCalls ( max_calls )
        if isinstance  ( max_iterations , integer_types ) and 1 < max_iterations :
            m.setMaxIterations    ( max_iterations  )
        if isinstance  ( strategy , integer_types       ) and 0 <= strategy <= 2 :
            m.setStrategy ( strategy )
            
        return m  

    # =========================================================================
    ## perform sPlot-analysis 
    #  @code
    #  r,f = model.fitTo ( dataset )
    #  model.sPlot ( dataset ) 
    #  @endcode 
    def sPlot ( self , dataset , silent = True ) : 
        """ Make sPlot analysis
        >>> r,f = model.fitTo ( dataset )
        >>> model.sPlot ( dataset ) 
        """
        assert self.alist2,\
               "PDF(%s) has empty ``alist2''/(list of components)" + \
               "no sPlot is possible" % self.name 
        
        with roo_silent ( silent ) :
            
            splot = ROOT.RooStats.SPlot ( rootID( "sPlot_" ) ,
                                          "sPlot"            ,
                                          dataset            ,
                                          self.pdf           ,
                                          self.alist2        )
        
            self.__splots += [ splot ]            
            return splot 
    # =========================================================================
    ## generate toy-sample according to PDF
    #  @code
    #  model  = ....
    #  data   = model.generate ( 10000 ) ## generate dataset with 10000 events
    #  varset = ....
    #  data   = model.generate ( 100000 , varset )
    #  data   = model.generate ( 100000 , varset , sample = True )     
    #  @endcode
    def generate ( self             ,
                   nEvents          ,
                   varset   = None  ,
                   binning  = None  ,
                   sample   = False , 
                   args     = ()    ) :
        """Generate toy-sample according to PDF
        >>> model  = ....
        >>> data   = model.generate ( 10000 ) ## generate dataset with 10000 events
        
        >>> varset = ....
        >>> data   = model.generate ( 100000 , varset )
        >>> data   = model.generate ( 100000 , varset , sample = True )
        """
        nEvents = self.gen_sample ( nEvents ) if sample else nEvents 
        assert 0 <= nEvents , 'Invalid number of Events %s' % nEvents  

        args = args + ( ROOT.RooFit.Name ( dsID() ) , ROOT.RooFit.NumEvents ( nEvents ) )

        if   isinstance ( binning , integer_types ) and 0 < binning :
            args    = args + ( ROOT.RooFit.AllBinned () , ) 
        elif binning is True :
            args    = args + ( ROOT.RooFit.AllBinned () , ) 
            binning = {}
            
        if   not varset :
            varset = ROOT.RooArgSet ( self.xvar )
        elif isinstance ( varset , ROOT.RooAbsReal ) :
            varset = ROOT.RooArgSet ( varset    )
            
        if not self.xvar in varset :
            vs = ROOT.RooArgSet()
            vs . add ( self.xvar )
            for  v in varset : vs.add ( v )
            varset = vs  

        from ostap.fitting.variables import KeepBinning

            
        with KeepBinning ( self.xvar ) : 

            if isinstance ( binning , dict ) :
                binning = binning.get ( self.xvar.name , None )                 
            if binning : self.xvar.bins = binning
            
            return self.pdf.generate (  varset , *args )

    # =========================================================================
    ## simple 'function-like' interface 
    def __call__ ( self , x , error = False , normalized = True ) :
        """ PDF as a ``function''
        >>> pdf  = ...
        >>> x = 1
        >>> y = pdf ( x ) 
        """
        return FUNC.__call__ ( self , x , error = error , normalized = normalized ) 

    # ========================================================================
    ## check minmax of the PDF using the random shoots
    #  @code
    #  pdf     = ....
    #  mn , mx = pdf.minmax()            
    #  @endcode 
    def minmax ( self , nshoots =  50000 ) :
        """Check min/max for the PDF using  random shoots 
        >>> pdf     = ....
        >>> mn , mx = pdf.minmax()        
        """
        ## try to get minmax directly from pdf/function 
        if self.tricks and hasattr ( self.pdf , 'function' ) :
            if hasattr ( self.pdf , 'setPars' ) : self.pdf.setPars() 
            f = self.pdf.function()
            if hasattr ( f , 'minmax' ) :
                try :
                    mn , mx = f.minmax()
                    if  0<= mn and mn <= mx and 0 < mx :   
                        return mn , mx
                except :
                    pass
            if hasattr ( f , 'max' ) :
                try :
                    mx = f.max()
                    if 0 < mx : return 0 , mx
                except :
                    pass
            if hasattr ( f , 'mode' ) :
                try :
                    mode = f.mode()
                    if mode in self.xvar :
                        mx = f ( mode ) 
                        if 0 < mx : return 0 , mx
                except :
                    pass 
                    

        ## check RooAbsReal functionality
        code = self.pdf.getMaxVal( ROOT.RooArgSet ( self.xvar ) )
        if 0 < code :
            mx = self.pdf.maxVal ( code )
            if 0 < mx : return 0 , mx
            
                 
        mn , mx = -1 , -10
        if hasattr ( self.pdf , 'min' ) : mn = self.pdf.min()
        if hasattr ( self.pdf , 'max' ) : mx = self.pdf.max()
        if 0 <= mn and mn <= mx and 0 < mx : return mn , mx
        
        ## now try to use brute force and random shoots 
        if not self.xminmax() : return ()
        
        mn  , mx = -1 , -10
        xmn , xmx = self.xminmax()
        for i in range ( nshoots ) : 
            xx = random.uniform ( xmn , xmx )
            with SETVAR ( self.xvar ) :
                self.xvar.setVal ( xx )
                vv = self.pdf.getVal()
                if mn < 0 or vv < mn : mn = vv
                if mx < 0 or vv > mx : mx = vv 
        return mn , mx 

    # ========================================================================
    ## clean some stuff 
    def clean ( self ) :
        self.__splots     = []
        self.__histo_data = None 
        self.__fit_result = None
        
    # ========================================================================
    ## get the effective RMS 
    def rms ( self , **kwargs ) :
        """Get the effective RMS
        >>>  fun = ...
        >>>  print 'RMS: %s ' % fun.rms()
        """
        fun  = self.fun 
        if   hasattr ( fun , 'rms'            ) : return fun.rms()        
        elif hasattr ( fun , 'Rms'            ) : return fun.Rms()        
        elif hasattr ( fun , 'RMS'            ) : return fun.RMS()        
        elif self.tricks and hasattr ( fun , 'function' ) :
            
            if   hasattr ( fun , 'setPars'    ) : fun.setPars()
            
            ff = fun.function()            
            if   hasattr ( ff , 'rms'        ) : return ff.rms()
            elif hasattr ( ff , 'variance'   ) : return ff.variance   ()**0.5  
            elif hasattr ( ff , 'dispersion' ) : return ff.dispersion ()**0.5 
            
        from ostap.stats.moments import rms as _rms
        return  self._get_stat_ ( _rms , **kwargs )

    # =========================================================================
    ## get the effective Skewness
    def skewness ( self , **kwargs ) :
        """Get the effective Skewness
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'SKEWNESS: %s ' % pdf.skewness()
        """
        ## use generic machinery 
        from ostap.stats.moments import skewness as _skewness
        return self._get_stat_ ( _skewness , **kwargs )

    # =========================================================================
    ## get the effective Kurtosis
    def kurtosis ( self , **kwargs ) :
        """Get the effective Kurtosis
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'KURTOSIS: %s ' % pdf.kurtosis()
        """
        ## use generic machinery 
        from ostap.stats.moments import kurtosis as _kurtosis
        return self._get_stat_ ( _kurtosis , **kwargs )

    # =========================================================================
    ## get the effective median
    def median ( self , **kwargs ) :
        """Get the effective median
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'MEDIAN: %s ' % pdf.median()
        """
        from ostap.stats.moments import median as _median
        return self._gets_stat_ ( _median , **kwargs )

    # =========================================================================
    ## get the effective mean
    def get_mean ( self , **kwargs ) :
        """Get the effective Mean
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'MEAN: %s ' % pdf.get_mean()
        """
        from ostap.stats.moments import mean as _mean
        return self._get_stat_ ( _mean , **kwargs )
    
    # =========================================================================
    ## get the effective moment for the distribution
    def moment ( self , N , **kwargs ) :
        """Get the effective moment
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'MOMENT: %s ' % pdf.moment( 10 )
        """
        ## use generic machinery 
        from ostap.stats.moments import moment as _moment
        return self._get_stat_ ( _moment , N , **kwargs ) 
    
    # =========================================================================
    ## get the effective central moment for the distribution
    def central_moment ( self , N , **kwargs ) :
        """Get the effective central moment
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'MOMENT: %s ' % pdf.moment( 10 )
        """
        from ostap.stats.moments import central_moment as _moment
        return self._get_stat_ ( _moment , N , **kwargs ) 

    # =========================================================================
    ## get the effective quantile 
    def quantile ( self , prob  , **kwargs ) :
        """Get the effective quantile
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'QUANTILE: %s ' % pdf.quantile ( 0.10 )
        """
        from ostap.stats.moments import quantile as _quantile
        return self._get_stat_ ( quantile , prob , **kwargs ) 

    # =========================================================================
    ## get the symmetric confidence interval 
    def cl_symm ( self , prob , x0 =  None , **kwargs ) :
        """Get the symmetric confidence interval 
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'CL :  ',  pdf.cl_symm ( 0.10 )
        """
        from ostap.stats.moments import cl_symm as _cl
        return self._get_stat_ ( _cl , prob , x0 , **kwargs ) 

    # =========================================================================
    ## get the asymmetric confidence interval 
    def cl_asymm ( self , prob , **kwargs ) :
        """Get the asymmetric confidence interval 
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'CL :  ',  pdf.cl_asymm ( 0.10 )
        """
        from ostap.stats.moments import cl_asymm as _cl
        return self._get_stat_ ( _cl , prob , **kwargs )
    
    # =========================================================================
    ## get the integral between xmin and xmax 
    def integral ( self , xmin , xmax , nevents = True ) :
        """Get integral between xmin and xmax
        >>> pdf = ...
        >>> print pdf.integral ( 0 , 10 )
        """
        ## check limits
        if self.xminmax() :
            mn , mx = self.xminmax() 
            xmin = max ( xmin , mn )  
            xmax = max ( xmax , mn )

        ## initialize the value and the flag 
        value , todo = 0 , True
        
        ## 1) make a try to use ``analytical'' integral 
        if self.tricks :
            try:
                if hasattr ( self.pdf , 'setPars'  ) : self.pdf.setPars() 
                fun          = self.pdf.function()
                value , todo = fun.integral ( xmin , xmax ) , False 
            except:
                pass

        ## 2) use numerical integration
        from ostap.math.integral import integral as _integral

        extended =  self.pdf.canBeExtended() or isinstance ( self.pdf , ROOT.RooAddPdf )
        
        if   todo and extended : value = _integral ( self , xmin , xmax )

        elif todo :
                        
            ## use unormalized PDF here to speed up the integration 
            ifun   = lambda x :  self ( x , error = False , normalized = False )
            value  = _integral ( ifun , xmin , xmax )
            norm   = self.pdf.getNorm ( self.vars )
            value /= norm

        if nevents and self.pdf.mustBeExtended () :
            evts = self.pdf.expectedEvents( self.vars )
            if evts <= 0 or iszero ( evts ) :
                self.warning ( "integral: expectedEvents is %s" % evts )
            value *= evts 

        return value

    # ==========================================================================
    ## convert PDF into TF1 object, e.g. to profit from TF1::Draw options
    #  @code
    #  pdf = ...
    #  tf1 = pdf.tf()
    #  tf1.Draw('colz')
    #  @endcode
    def tf ( self , xmin = None , xmax = None ) :
        """Convert PDF  to TF1 object, e.g. to profit from TF1::Draw options
        >>> pdf = ...
        >>> tf2 = pdf.tf()
        >>> tf1.Draw('colz')
        """
        def _aux_fun_ ( x , pars = [] ) :
            return self ( x[0] , error = False )
        
        if xmin == None and self.xminmax() : xmin = self.xminmax()[0]
        if xmax == None and self.xminmax() : xmax = self.xminmax()[1]

        if xmin == None : xmin = 0.0
        if xmax == None : xmin = 1.0
        
        from ostap.core.core import fID
        return ROOT.TF1 ( fID() , _aux_fun_ , xmin , xmax ) 

    # =========================================================================
    ## helper function to create the histogram
    def make_histo ( self  ,
                     nbins    = 100   , xmin = None , xmax = None ,
                     hpars    = ()    , 
                     histo    = None  ) :
        """Create the histogram according to specifications
        """
        
        import ostap.histos.histos
        # histogram is provided 
        if histo :
            
            assert isinstance ( histo , ROOT.TH1 ) and not isinstance ( histo , ROOT.TH2 ) , \
                   "Illegal type of ``histo''-argument %s" % type( histo )
            
            histo = histo.clone()
            histo.Reset()

        # arguments for the histogram constructor 
        elif hpars :
            
            histo = ROOT.TH1F ( hID() , 'PDF%s' % self.name , *hpars  )
            if not histo.GetSumw2() : histo.Sumw2()

        # explicit construction from (#bins,min,max)-triplet  
        else :
            
            assert is_integer ( nbins ) and 0 < nbins, \
                   "Wrong ``nbins''-argument %s" % nbins 
            if xmin == None and self.xminmax() : xmin = self.xminmax()[0]
            if xmax == None and self.xminmax() : xmax = self.xminmax()[1]
            
            histo = ROOT.TH1F ( hID() , 'PDF%s' % self.name , nbins , xmin , xmax )
            if not histo.GetSumw2() : histo.Sumw2()

        return histo 
    
    # ==========================================================================
    ## Convert PDF to the 1D-histogram  in correct way.
    #  @code
    #  pdf = ...
    #  h1  = pdf.histo ( 100 , -1 , 10 ) ## specify histogram parameters
    #  histo_template = ...
    #  h2  = pdf.histo ( histo = histo_template ) ## use historgam template
    #  h3  = pdf.histo ( ... , integral = True  ) ## use PDF integral within the bin  
    #  h4  = pdf.histo ( ... , density  = True  ) ## convert to "density" histogram 
    #  @endcode
    #  @see PDF.roo_histo
    #  Unlike  <code>PDF.roo_histo</code> method, PDF is integrated within the bin
    def histo ( self             ,
                nbins    = 100   , xmin = None , xmax = None ,
                hpars    = ()    , 
                histo    = None  ,
                integral = True  ,
                errors   = False ) :
        """Convert PDF to the 1D-histogram in correct way
        - Unlike  <code>PDF.roo_histo</code> method, PDF is integrated within the bin
        >>> pdf = ...
        >>> h1  = pdf.histo ( 100 , 0. , 10. ) ## specify histogram parameters
        >>> histo_template = ...
        >>> h2  = pdf.histo ( histo = histo_template ) ## use historgam template
        >>> h3  = pdf.histo ( ... , integral = True  ) ## use PDF integral within the bin  
        >>> h4  = pdf.histo ( ... , density  = True  ) ## convert to 'density' histogram 
        """
        
        histo = self.make_histo ( nbins = nbins ,
                                  xmin  = xmin  ,
                                  xmax  = xmax  ,
                                  hpars = hpars ,
                                  histo = histo )

        # loop over the historgam bins 
        for i , x , y in histo.items() :

            xv , xe = x.value() , x.error()
            
            # value at the bin center 
            c = self ( xv , error = errors ) 

            if not integral : 
                histo[i] = c
                continue

            # integral over the bin 
            v  = self.integral( xv - xe , xv + xe )
            
            if errors :
                if    0 == c.cov2 () : pass
                elif  0 != c.value() and 0 != v : 
                    v = c * ( v / c.value() )
                    
            histo[i] = v 

        return histo

    # ==========================================================================
    ## Convert PDF to the 1D-histogram, taking PDF-values at bin-centres
    #  @code
    #  pdf = ...
    #  h1  = pdf.roo_histo ( 100 , -1 , 10 ) ## specify histogram parameters
    #  histo_template = ...
    #  h2  = pdf.roo_histo ( histo = histo_template ) ## use historgam template
    #  @endcode
    #  @see RooAbsPdf::createHistogram
    #  @see RooAbsPdf::fillHistogram
    #  @see PDF.histo
    def roo_histo ( self             ,
                    nbins    = 100   , xmin = None , xmax = None ,
                    hpars    = ()    , 
                    histo    = None  ,
                    events   = True  ) : 
        """Convert PDF to the 1D-histogram, taking PDF-values at bin-centres
        - see RooAbsPdf::createHistogram
        - see RooAbsPdf::fillHistogram
        - see PDF.histo
        >>> pdf = ...
        >>> h1  = pdf.roo_histo ( 100 , 0. , 10. ) ## specify histogram parameters
        >>> histo_template = ...
        >>> h2  = pdf.roo_histo ( histo = histo_template ) ## use histogram template
        """

        histo = self.make_histo ( nbins = nbins ,
                                  xmin  = xmin  ,
                                  xmax  = xmax  ,
                                  hpars = hpars ,
                                  histo = histo )

        hh = self.pdf.createHistogram (
            hID()     ,
            self.xvar ,
            self.binning ( histo.GetXaxis() , 'histo1x' ) ,
            ROOT.RooFit.Extended ( False ) ,
            ROOT.RooFit.Scaling  ( False ) ,            
            )
        
        for i in hh : hh.SetBinError ( i , 0 ) 
        
        if events and self.pdf.mustBeExtended() :
            
            for i , x , y in hh.items() :
                volume  = 2*x.error() 
                hh[i]  *= volume
                
            hh *= self.pdf.expectedEvents ( self.vars ) / hh.sum() 
                
        histo += hh
        
        return histo 

    # ==========================================================================
    ## create the residual histogram  :  (data - pdf)
    #  @param   data_histo  the data histogram
    #  @return  the residual histogram 
    #  @code
    #  data = ...
    #  pdf  = ...
    #  pdf.fitTo ( data )
    #
    #  histo = ..
    #  data.project ( histo , 'var1' )
    #
    #  residual = pdf.residual_histo ( histo )
    #  @endcode 
    def residual_histo  ( self , data_histo ) :
        """Create the residual histogram   (data - fit)
        >>> data = ... 
        >>> pdf  = ...
        >>> pdf.fitTo ( data )
        
        >>> histo = ..
        >>> data.project ( histo , 'var1' )
        
        >>> residual = pdf.residual_histo ( histo )
        """
        
        hpdf = self.histo ( histo = data_histo )

        for i in hpdf :
            
            d       = data_histo [i]
            v       = hpdf       [i]
            hpdf[i] = d - v.value()    ## data - pdf 
            
        return hpdf 

    # ==========================================================================
    ## create the pull histogram  : (data - pdf)/data_error
    #  @param   data_histo  the data histogram
    #  @return  the pull histogram 
    #  @code
    #  data = ...
    #  pdf  = ...
    #  pdf.fitTo ( data )
    #
    #  histo = ..
    #  data.project ( histo , 'var1' )
    #
    #  pull = pdf.pull_histo ( histo )
    #  @endcode 
    def pull_histo  ( self , data_histo ) :
        """Create the residual histogram   (data - fit)
        >>> data = ... 
        >>> pdf  = ...
        >>> pdf.fitTo ( data )
        
        >>> histo = ..
        >>> data.project ( histo , 'var1' )
        
        >>> pull = pdf.pull_histo ( histo )
        """
        h = self.residual_histo ( data_histo = data_histo )
        
        for i in h :
            v    = h         [i]
            e    = data_histo[i].error()
            if 0 < e : h[i] = v / e   ## (data-pdf)/data_error

        return h 

    # ==========================================================================
    ## get the residual histogram : (data-fit) 
    #  @see PDF.histo
    #  @see PDF.residual_histo
    #  @see PDF.make_histo
    #  @code
    #  data = ...
    #  pdf  = ...
    #  pdf.fitTo ( data )
    #  residual = pdf.residual ( data , nbins = 100 ) 
    #  @endcode 
    def residual ( self  , dataset , **kwargs ) :
        """Get the residual histogram
        - see PDF.histo
        - see PDF.residual_histo
        - see PDF.make_histo

        >>> data = ...
        >>> pdf  = ...
        >>> pdf.fitTo ( data )
        >>> residual = pdf.residual ( data , nbins = 100 ) 
        """
        hdata = self.make_histo ( **kwargs )
        dataset.project ( hdata , self.xvar.name )
        return self.residual_histo ( hdata ) 
        

    # ==========================================================================
    ## get the pull histogram : (data-fit)/data_error 
    #  @see PDF.histo
    #  @see PDF.residual_histo
    #  @see PDF.make_histo
    #  @code
    #  data = ...
    #  pdf  = ...
    #  pdf.fitTo ( data )
    #  residual = pdf.pull ( data , nbins = 100 ) 
    #  @endcode 
    def pull ( self  , dataset , **kwargs ) :
        """Get the pull  histogram: (data-fit)/data_error
        - see PDF.histo
        - see PDF.residual_histo
        - see PDF.make_histo

        >>> data = ...
        >>> pdf  = ...
        >>> pdf.fitTo ( data )
        >>> residual = pdf.residual ( data , nbins = 100 ) 
        """
        hdata = self.make_histo ( **kwargs )
        dataset.project ( hdata , self.xvar.name )
        return self.pull_histo ( hdata ) 
        

        
    # ==========================================================================
    # Several purely technical methods 
    # ==========================================================================
        
    # ==========================================================================
    ## Pure technical methods to make a sum of PDFs with *constant* equal fractions
    #  @code 
    #  sum = self.make_sum ( 'A' , 'B' , pdf1 , pdf2 , pdf3 )
    #  @endcode
    #  It is very useful for e.g. creation of ""symmetrized" PDFs
    #  f(x,y) = 0.5 * [f1(x)*f2(y)] + 0.5* [ f2(x)*f1(y) ]    
    def make_sum ( self , name , title , pdf1 , pdf2 , *pdfs ) :
        """ Pure technical methods to make a sum of PDFs with *constant* equal fractions
        >>> sum = self.make_sum ( 'A' , 'B' , pdf1 , pdf2 , pdf3 )
        It is very useful for e.g. creation of 'symmetrized'PDF:
        f(x,y) = 0.5 * [f1(x)*f2(y)] + 0.5* [ f2(x)*f1(y) ]    
        """
        if self.name : name = name + '_' + self.name 
        _pdfs  = ROOT.RooArgList()
        _pdfs.add ( pdf1 )
        _pdfs.add ( pdf2 )
        for p in pdfs : _pdfs.add ( p )
        n = len(_pdfs)
        _fracs = [] 
        for i in  range(n) :
            f = ROOT.RooConstVar ( "Fraction%d_%s"     % ( i+1 , name         ) ,
                                   "fraction%d(%s,%s)" % ( i+1 , name , title ) , 1.0 / n )
            _fracs.append ( f )
            
        _rlst = ROOT.RooArgList()
        for f in _fracs : _rlst.add ( f ) 
        ## create PDF 
        result = ROOT.RooAddPdf ( name , title , _pdfs , _rlst , False )
        ##
        self.aux_keep.append ( _pdfs  )
        self.aux_keep.append ( _fracs )
        self.aux_keep.append ( _rlst  )
        #
        return result
    
    # =============================================================================
    ## make list of variables/fractions for compound PDFs 
    def make_fracs ( self             ,
                     N                ,
                     pname            ,
                     ptitle           ,
                     fractions = True ,
                     recursive = True ,
                     fracs     = []   )  :
        """Make list of variables/fractions for compound PDF
        """
        assert is_integer ( N ) and 2 <= N , \
               "PDF.make_fracs: Invalid N=%s/%s" % ( N, type ( N ) )
        ##
        ufracs   = [] 
        n        =  ( N - 1 ) if fractions else N
        NN       = n + 1
        vminmax  =  ( 0 , 1 ) if fractions else ( 0 , 1.e+7 )
        value    = 1 
        prod     = 1.0
        for i in range ( 0 , n ) :            
            if not fractions :                
                fv = 1.0 / NN
                if recursive :    
                    fv   /=  prod
                    prod *= ( 1.0 - fv )
                value = fv 
            ## finally create the fraction
            fi = get_i ( fracs , i , None )
            
            if isinstance ( fi , num_types ) and not vminmax[0] <= fi <= vminmax[1] :
                self.warning ("make_fracs: fraction %s is outside the interval %s, ignore" % ( fi , vminmax ) )
                fi = None
                
            f  = self.make_var ( fi , pname % i , ptitle % i , None , value , *vminmax ) 
            ufracs.append ( f )
            
        return ufracs
    
    # =============================================================================
    ## helper function to build composite (non-extended) PDF from components 
    def add_pdf ( self             ,
                  pdflist          ,
                  name             ,
                  title            ,
                  fname            ,
                  ftitle           ,
                  recursive = True ,
                  fractions = []   ) :
        """Helper function to build composite (non-extended) PDF from components 
        """
        ##
        pdfs   = ROOT.RooArgList() 
        for pdf in pdflist : pdfs.add  ( pdf )
        fs     = self.make_fracs ( len ( pdfs )          ,
                                   fname                 ,
                                   ftitle                ,
                                   fractions = True      ,
                                   recursive = recursive ,
                                   fracs     = fractions )
        fracs  = ROOT.RooArgList()
        for f in fs : fracs.add ( f ) 
        pdf    = ROOT.RooAddPdf ( name , title , pdfs , fracs , recursive )
        ##
        self.aux_keep.append ( pdf   )
        self.aux_keep.append ( pdfs  )
        self.aux_keep.append ( fracs )
        ##
        return pdf , fracs , pdfs

    # =========================================================================
    ## create popular 1D ``background''  function
    #  @param bkg  the type of background function/PDF
    #  @param name the name of background function/PDF
    #  @param xvar the observable
    #  Possible values for <code>bkg</code>:
    #  - None or 0          : <code>Flat1D</code>
    #  - positive integer N : <code>Bkg_pdf(power=N)</code> 
    #  - negative integer K : <code>PolyPos_pdf(power=abs(K))</code> 
    #  - any Ostap/PDF      : PDF will be copied or cloned  
    #  - any RooAbsPdf    P : <code>Generic1D_pdf(pdf=P)</code> 
    #  - any RooAbsReal   V : <code>Bkg_pdf(power=0,tau=V)</code> 
    #  - math.exp           : <code>Bkg_pdf(power=0)</code>
    #  - ''  , 'const', 'constant' , 'flat' , 'uniform' : <code>Flat1D</code>
    #  - 'p0', 'pol0' , 'poly0' : <code>Flat1D</code>
    #  - 'e' , 'exp'  , 'expo'  : <code>Bkg_pdf(power=0)</code>
    #  - 'e+', 'exp+' , 'expo+' : <code>Bkg_pdf(power=0)</code> with tau>0
    #  - 'e-', 'exp-' , 'expo-' : <code>Bkg_pdf(power=0)</code> with tau<0
    #  - 'e0', 'exp0' , 'expo0' : <code>Bkg_pdf(power=0)</code>
    #  - 'eN', 'expN' , 'expoN' : <code>Bkg_pdf(power=N)</code>
    #  - 'pN', 'polN' , 'polyN' : <code>PolyPos_pdf(power=N)</code>
    #  - 'iN', 'incN' , 'incrN','increasingN' : <code>Monotonic_pdf(power=N,increasing=True)</code>
    #  - 'dN', 'decN' , 'decrN','decreasingN' : <code>Monotonic_pdf(power=N,increasing=False)</code>     
    #  @see ostap.fitting.backrgound.make_bkg 
    def make_bkg ( self , bkg , name , xvar , **kwargs ) :
        """Create popular 1D ``background''  function.
        
        Possible values for ``bkg'':
        
        - None or 0                               : Flat1D
        - positive integer ``N''                  : Bkg_pdf(power=N)
        - negative integer ``K''                  : PolyPos_pdf(power=abs(K))
        - any Ostap-PDF                           : PDF will be copied or cloned  
        - RooAbsPdf      ``pdf''                  : Generic1D_pdf(pdf=pdf)
        - RooAbsReal     ``var''                  : Bkg_pdf(power=0,tau=var)
        - math.exp                                : Bkg_pdf(power=0)
        - 'const' or 'constant'                   : Flat1D
        - '' , 'flat' or 'uniform'                : Flat1D
        - 'e' , 'exp'  or 'expo'                  : Bkg_pdf(power=0)
        - 'e+', 'exp+' or 'expo+'                 : Bkg_pdf(power=0) with tau>0 (increasing)
        - 'e-', 'exp-' or 'expo-'                 : Bkg_pdf(power=0) with tau<0 (decreasing)
        - 'e0', 'exp0' or 'expo0'                 : Bkg_pdf(power=0) 
        - 'eN', 'expN' or 'expoN'                 : Bkg_pdf(power=N)
        - 'p0', 'pol0' or 'poly0'                 : Flat1D
        - 'pN', 'polN' or 'polyN'                 : PolyPos_pdf(power=N)
        - 'iN', 'incN' , 'incrN' or 'increasingN' : Monotonic_pdf(power=N,increasing=True)
        - 'dN', 'decN' , 'decrN' or 'decreasingN' : Monotonic_pdf(power=N,increasing=False)
        For more information see 
        see Ostap.FitBkgModels.make_bkg 
        """
        from ostap.fitting.background import make_bkg as bkg_make
        return bkg_make ( bkg    = bkg         ,
                          name   = name        ,
                          xvar   = xvar        ,
                          logger = self.logger , **kwargs ) 

    # =========================================================================
    ## Load parameters from:
    #    - external dictionary <code>{ name : value }</code>
    #    - sequence of <code>RooAbsReal</code> object
    #    - <code>ROOT.RooFitResult</code> object 
    #  @code
    #  pdf     = ...
    #  dataset = ...
    #  params  = { 'A' : 10 , 'B' : ... }
    #  pdf.load_params ( dataset , params ) 
    #  params  = ( A , B , C , ... )
    #  pdf.load_params ( dataset , params )  
    #  @endcode 
    def load_params ( self , dataset = None , params = {} , silent = False  ) :
        """Load parameters from
        - external dictionary `{ name : value }`
        - sequence of `RooAbsReal` objects
        - `RooFitResult` object
        
        >>> pdf      = ...
        >>> dataset = ... 
        >>> params = { 'A' : 10 , 'B' : ... }
        >>> pdf.load_params ( dataset , params ) 
        >>> params = ( A , B , C , ... )
        >>> pdf.load_params ( dataset , params )  
        """
        ## nothing to load 
        if not params : return 

        if isinstance ( params , ROOT.RooFitResult ) :
            params = params.dct_params () 
        
        ## get the list of the actual parameters 
        pars = self.pdf.getParameters ( dataset )

        table = [] 
        if isinstance ( params , dictlike_types ) :
            keys   = set () 
            for key in params :
                for p in pars :
                    if not hasattr ( p  , 'setVal' ) : continue
                    if p.name != key                 : continue
                    
                    v  = params[key]
                    vv = float ( v  )
                    pv = p.getVal ()   
                    if vv != pv : 
                        p.setVal   ( vv )
                        item = p.name , "%-14.6g" % pv , "%-+14.6g" % vv 
                        table.append ( item ) 
                    keys.add ( key )

            not_used = set ( params.keys() ) - keys 

        ## list of objects 
        else :
            
            keys = set()        
            for i , pp in enumerate ( params ) :  
                if not isinstance ( pp , ROOT.RooAbsReal ) : continue
                for p in pars :
                    if not hasattr ( p  , 'setVal' )       : continue
                    if p.name != pp.name                   : continue
                    
                    vv = float ( pp )
                    pv = p.getVal () 
                    if vv != pv :
                        p.setVal   ( vv )
                        item = p.name , "%-14.6g" % pv , "%-+14.6g" % vv 
                        table.append ( item ) 
                    keys.add  ( i )

            not_used = []
            for i , pp in enumerate ( params ) :  
                if i in keys : continue
                not_used.append ( pp )

        if not silent :
            
            table.sort()
            npars = len ( table )
            
            if npars :            
                title = 'Parameters loaded: %s' % npars 
                table = [ ('Parameter' ,'old value' , 'new value' ) ] + table
                import ostap.logger.table
                table = ostap.logger.table.table ( table , title , prefix = "# " )
                
            self.info ( "%s parameters loaded:\n%s" % ( npars , table ) ) 
            
            not_used = list ( not_used )
            not_used.sort() 
            if not_used :
                self.warning ("Following keys are unused %s" % not_used ) 
        
        return 

    # =========================================================================
    ## get all parameters/variables in form of dictionary
    #  @code
    #  pdf    = ...
    #  params = pdf.parameters ( dataset ) 
    #  @endcode
    def parameters ( self , dataset = None ) :
        """ Get all parameters/variables in form of dictionary
        >>> pdf    = ...
        >>> params = pdf.parameters ( dataset ) 
        """
        
        ## get the list of the actual parameters 
        pars = self.pdf.getParameters ( dataset )

        tmp    = {}
        for p in pars :
            if not isinstance ( p, ROOT.RooAbsCategory ) :
                tmp [ p.name ] = p.value
                
        keys   = tmp.keys()
        result = {} 
        for key in sorted ( keys ) : result [ key ] = tmp [ key ] 
            
        return result 

    # ========================================================================
    ## get the parameter value by name
    #  @code
    #  pdf = ...
    #  p   = pdf.parameter  ( 'A' )
    #  @endcode
    def parameter ( self , param , dataset = None ) :
        """Get the parameter value by name
        >>> pdf = ...
        >>> p   = pdf.parameter  ( 'A' )
        """
        ## get the list of the actual parameters 
        pars = self.pdf.getParameters ( dataset )

        for p in pars :
            if p.name == param : return p
            
        self.error ( "No parameter %s defined" % param )
        raise KeyError ( "No parameter %s defined" % param )

    # ==========================================================================
    ## get parameter by name 
    #  @code
    #  pdf = ...
    #  a   = pdf['A']
    #  @endcode
    def __getitem__ ( self , param ) :
        """Get parameter by name 
        >>> pdf = ...
        >>> a   = pdf['A']
        """
        ## get the list of the actual parameters 
        pars = self.pdf.getParameters ( None )
        for p in pars :
            if p.name == param : return p
        raise KeyError ( "No parameter %s defined" % param )
        
# =============================================================================
##  helper utilities to imlement resolution models.
# =============================================================================
class _CHECKMEAN(object) :
    check = True
def checkMean() :
    return True if  _CHECKMEAN.check else False
class Resolution(object) :    
    def __init__  ( self , resolution = True ) :
        self.check = False if resolution else True 
    def __enter__ ( self ) :
        self.old         = _CHECKMEAN.check 
        _CHECKMEAN.check =  self.check
    def __exit__  ( self , *_ ) :
        _CHECKMEAN.check =  self.old 
# =============================================================================
## helper base class for implementation  of various helper pdfs 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class MASS(PDF) :
    """Helper base class for implementation of various pdfs
    It is useful for ``peak-like'' distributions, where one can talk about
    - ``mean/location''
    - ``sigma/width/scale'' 
    """
    def __init__ ( self            ,
                   name            ,
                   xvar            ,
                   mean     = None ,
                   sigma    = None ) : 
        
        m_name  = "m_%s"     % name
        m_title = "mass(%s)" % name

        if   isinstance ( xvar , ROOT.TH1   ) :
            m_title = xvar.GetTitle ()            
            xvar    = xvar.xminmax  ()
        elif isinstance ( xvar , ROOT.TAxis ) :
            xvar    = xvar.GetXmin() , mass.GetXmax()

        ## create the variable 
        if isinstance ( xvar , tuple ) and 2 == len(xvar) :  
            xvar = self.make_var ( xvar       , ## var 
                                   m_name     , ## name 
                                   m_title    , ## title/comment
                                   None       , ## fix ? 
                                   *xvar      ) ## min/max 
        elif isinstance ( xvar , ROOT.RooAbsReal ) :
            xvar = self.make_var ( xvar       , ## var 
                                   m_name     , ## name 
                                   m_title    , ## title/comment
                                   fix = None ) ## fix ? 
        else :
            raise AttributeError("MASS: Unknown type of ``xvar'' parameter %s/%s" % ( type ( xvar ) , xvar ) )

        ## intialize the base 
        PDF.__init__ ( self , name , xvar = xvar )
        

        limits_mean  = ()
        limits_sigma = ()
        
        if   self.xminmax() :            
            mn , mx = self.xminmax()
            dm      =  mx - mn
            limits_mean  = mn - 0.2 * dm , mx + 0.2 * dm
            sigma_max    =  2 * dm / math.sqrt(12)  
            limits_sigma = 1.e-3 * sigma_max , sigma_max 
        #
        ## mean-value
        #
        self.__mean = self.make_var ( mean              ,
                                      "mean_%s"  % name ,
                                      "mean(%s)" % name , mean , *limits_mean )
        ## 
        if checkMean () and self.xminmax() : 
            mn , mx = self.xminmax() 
            dm      =  mx - mn
            if   self.mean.isConstant() :
                if not mn <= self.mean.getVal() <= mx : 
                    raise AttributeError ( 'MASS(%s): Fixed mass %s is not in mass-range (%s,%s)' % ( name , self.mean.getVal() , mn , mx  ) )
            elif self.mean.minmax() :
                mmn , mmx = self.mean.minmax()
                self.mean.setMin ( max ( mmn , mn - 0.1 * dm ) )
                self.mean.setMax ( min ( mmx , mx + 0.1 * dm ) )
                self.debug ( 'mean range is redefined to be %s' % list( self.mean.minmax() ) )
        #
        ## sigma
        #
        self.__sigma = self.make_var ( sigma               ,
                                       "sigma_%s"   % name ,
                                       "#sigma(%s)" % name , sigma , *limits_sigma )
        
        ## save the configuration
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'mean'  : self.mean  ,
            'sigma' : self.sigma
            }
            

    @property 
    def mass ( self ) :
        """``mass''-variable (the same as ``x'' or ``xvar'')"""
        return self.xvar
    
    @property
    def mean ( self ):
        """``mean/location''-variable (the same as ``location'')"""
        return self.__mean
    @mean.setter
    def mean ( self , value ) :
        value =  float ( value )
        if self.xminmax() : 
            mn , mx = self.xminmax()
            dm = mx - mn
            m1 = mn - 1.0 * dm
            m2 = mx + 1.0 * dm
            if not m1 <= value <= m2 :
                self.warning ("``mean'' %s is outside the interval  %s,%s"  % ( value , m1 , m2 ) )                
        self.mean.setVal ( value )
        
    @property
    def location ( self ):
        """``location/mean''-variable (the same as ``mean'')"""
        return self.mean
    @location.setter
    def location ( self , value ) :
        self.mean =  value 
    
    @property
    def sigma ( self ):
        """``sigma/width/scale''-variable"""
        return self.__sigma
    @sigma.setter
    def sigma ( self , value ) :
        value =   float ( value )
        if self.xminmax() : 
            mn , mx = self.xminmax()
            dm = mx - mn
            smax = 2 * dm / math.sqrt ( 12 ) 
            smin = 2.e-5 * smax  
            if not smin <= value <= smax :
                self.warning ("``sigma'' %s is outside the interval (%s,%s)" % ( value , smin , smax ) )
        self.sigma.setVal ( value )


# =============================================================================
## @class RESOLUTION
#  helper base class  to parameterize the resolution
#  - It allows setting of the <code>mean</code> to zero,
#  - It containg "fudge-factor" for the resolution parameter <code>sigma</code>
#  - It simplify creation of the soft/gaussian constraint for the "fudge-factor"
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2017-07-13
class RESOLUTION(MASS) :
    """Helper base class  to parameterize the resolution
    - It allows setting of the ``mean'' to zero,
    - It containg ``fudge-factor'' for the resolution parameter ``sigma''
    - It simplify creation of the soft/gaussian constraint for the ``fudge-factor''    
    """
    ## constructor
    #  @param name   the name of PDF
    #  @param xvar   the variable/observable
    #  @param sigma  sigma/resoltuion parameter 
    #  @param mean   "mean"-variable
    #  @param fudge  "fudge-factor" to be aplied to sigma
    def __init__ ( self            ,
                   name            ,
                   xvar     = None ,
                   sigma    = None , 
                   mean     = None ,
                   fudge    = 1.0  ) :
        
        with Resolution() :
            super(RESOLUTION,self).__init__ ( name  = name  ,
                                              xvar  = xvar  ,
                                              sigma = sigma ,
                                              mean  = mean  )
            
        self.__fudge            = fudge
        
        self.__fudge_constraint = None
        
        if isinstance ( fudge , VE ) :
            
            assert 0 < fudge.value() and 0 < fudge.cov2(),\
                   "Invalid value for ``fudge-factor'': %s" % s 
            
            value  = fudge.value()
            error  = fudge.error()
            vmin   = max ( 0 , value - 5 * error )
            vmax   =           value + 5 * error
            
            ## make fudge-factor to be a variable 
            self.__fudge = self.make_var ( value ,
                                           'fudge_factor_%s'  % self.name ,
                                           'fudge_factor(%s)' % self.name ,
                                           None                           ,
                                           value , vmin , vmax            )
            
            ## create soft/gaussian constraint for fudge-factor
            self.__fudge_constraint = self.soft_constraint (
                self.fudge ,
                fudge      ,
                name  = 'Fudge_constraint_%s'  % self.name ,
                title = 'Fudge_constraint(%s)' % self.name )

        if isinstance ( self.fudge , num_types ) and 1 == fudge : 
            ## corrected sigma is trivial 
            self.__sigma_corr   =    self.sigma
        else : 
            ## create the corrected sigma 
            self.__sigma_corr = self.vars_multiply ( self.sigma ,
                                                     self.fudge ,
                                                     'Corrected_%s' % self.sigma.name )
            
        ## save the configuration
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'mean'  : self.mean  ,
            'sigma' : self.sigma ,
            'fudge' : self.fudge ,
            }

    @property 
    def fudge ( self ) :
        """``fudge'' : fudge factor for resolution"""
        return self.__fudge
    @property 
    def fudge_constraint ( self ) :
        """``fudge_constraint'' : constraint for fudge factor for the resolution"""
        return self.__fudge_constraint
    @property 
    def sigma_corr ( self ) :
        """``sigma_corr'' : the corrected sigma parameter: sigma*fudge """
        return self.__sigma_corr

    
# =============================================================================
## @class Flat1D
#  The most trivial 1D-model - constant
#  @code 
#  pdf = Flat1D( 'flat' , xvar = ... )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
class Flat1D(PDF) :
    """The most trival 1D-model - constant
    >>> pdf = Flat1D ( 'flat' , xvar = ... )
    """
    def __init__ ( self , xvar , name = 'Flat1D' , title = '' ) :
        
        PDF.__init__ ( self  , name , xvar ) 
        
        if not title : title = 'flat1(%s)' % name 
        self.pdf = Ostap.Models.Uniform ( name , title , self.xvar )
        assert 1 == self.pdf.dim() , 'Flat1D: wrong dimensionality!'
        
        ## save configuration
        self.config = {
            'xvar'     : self.xvar ,
            'name'     : self.name ,            
            'title'    : title     
            }
        
# =============================================================================
## @class Generic1D_pdf
#  "Wrapper" over generic RooFit (1D)-pdf
#  @code
#  raw_pdf = RooGaussian  ( ...     )
#  pdf     = Generic1D_pdf ( raw_pdf , xvar = x )  
#  @endcode 
#  If more functionality is required , more actions are possible:
#  @code
#  ## for sPlot 
#  pdf.alist2 = ROOT.RooArgList ( n1 , n2 , n3 ) ## for sPlotting 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-03-29
class Generic1D_pdf(PDF) :
    """Wrapper for generic RooFit pdf
    >>> raw_pdf = RooGaussian   ( ...     )
    >>> pdf     = Generic1D_pdf ( raw_pdf , xvar = x )
    """
    ## constructor 
    def __init__ ( self , pdf , xvar = None  ,
                   name           = ''    ,
                   special        = False ,
                   add_to_signals = True  ,
                   prefix         = ''    ,
                   suffix         = ''    ) :
        """Wrapper for generic RooFit pdf        
        >>> raw_pdf = RooGaussian   ( ...     )
        >>> pdf     = Generic1D_pdf ( raw_pdf , xvar = x )
        """
        assert xvar and isinstance ( xvar , ROOT.RooAbsReal ) , "``xvar'' must be ROOT.RooAbsReal"
        assert pdf  and isinstance ( pdf  , ROOT.RooAbsReal ) , "``pdf'' must be ROOT.RooAbsReal"
        
        name = name if name else prefix + pdf.GetName () + suffix 
        ## initialize the base 
        PDF . __init__ ( self , name , xvar , special = special )
        ##
        if not self.special :
            assert isinstance ( pdf  , ROOT.RooAbsPdf ) , "``pdf'' must be ROOT.RooAbsPdf"

        ## PDF itself 
        self.pdf  = pdf

        if not self.xvar in self.pdf.getParameters ( 0 ) : 
            self.warning ("Function/PDF does not depend on xvar=%s" % self.xvar.name )
            
        ## add it to the list of signal components ?
        self.__add_to_signals = True if add_to_signals else False
        
        if self.add_to_signals :
            self.signals.add ( self.pdf )
        
        ## save the configuration
        self.config = {
            'pdf'            : self.pdf            ,
            'xvar'           : self.xvar           ,
            'name'           : self.name           , 
            'special'        : self.special        , 
            'add_to_signals' : self.add_to_signals ,
            'prefix'         : prefix              ,
            'suffix'         : suffix              ,            
            }

        self.checked_keys.add  ( 'pdf'     )
        self.checked_keys.add  ( 'xvar'    )
        self.checked_keys.add  ( 'special' )
        
    @property
    def add_to_signals ( self ) :
        """``add_to_signals'' : should PDF be added into list of signal components?"""
        return self.__add_to_signals 
        
    
# =============================================================================
## @class Sum1D
#  Non-extended sum of two PDFs
#  @code
#  pdf1 = ...
#  pdf2 = ...
#  sum  = Sum1D ( pdf1 , pdf2 ) 
#  @endcode
#  It is just a small wrapper for <code>ROOT.RooAddPdf</code>
#  @see RooAddPdf 
class Sum1D(PDF) :
    """Non-extended sum of two PDFs:
    
    It is just a small wrapper for <code>ROOT.RooAddPdf</code>
    - see RooAddPdf 
    
    pdf1 = ...
    pdf2 = ...
    sum  = Sum1D ( pdf1 , pdf2 ) 

    """
    def __init__ ( self            ,
                   pdf1            ,
                   pdf2            ,  
                   xvar     = None , 
                   name     = ''   , 
                   fraction = None ) :

        if   isinstance ( pdf1 , PDF ) :
            assert ( not xvar ) or xvar is pdf1.xvar, "Invalid xvar/pdf1.xvar: %s/%s" %( xvar , pdf1.xvar )             
            xvar = pdf1.xvar        
        elif isinstance ( pdf1 , ROOT.RooAbsPdf ) and xvar and isinstance ( xvar , ROOT.RooAbsReal ) :            
            pdf1 = Generic1D_pdf ( pdf1 , xvar )
        else :
            raise TypeError ( "Invalid type: pdf1, xvar %s/%s , %s,%s" % ( pdf1, type(pdf1) , xvar , type(xvar) ) )

        if   isinstance ( pdf2 , PDF ) and xvar in pdf2.vars : pass 
        elif isinstance ( pdf2 , ROOT.RooAbsPdf ) and xvar and isinstance ( xvar , ROOT.RooAbsReal ) :            
            pdf2 = Generic1D_pdf ( pdf2 , xvar )
        else :
            raise TypeError ( "Invalid type: pdf1, xvar %s/%s , %s,%s" % ( pdf2, type(pdf2) , xvar , type(xvar) ) )

        name = name if name else 'Sum_%s_%s' % (  pdf1.name , pdf2.name ) 

        ## initialize the base class
        PDF.__init__ ( self , name , xvar )

        self.__pdf1     = pdf1
        self.__pdf2     = pdf2
        
        self.__fraction = self.make_var ( fraction ,
                                          'f_%s_%s'            % ( pdf1.name , pdf2.name ) ,
                                          'Fraction:(%s)+(%s)' % ( pdf1.name , pdf2.name ) ,
                                          fraction , 0 , 1 )
        
        self.alist1     = ROOT.RooArgList ( self.pdf1.pdf ,
                                            self.pdf2.pdf )
        self.alist2     = ROOT.RooArgList ( self.fraction )
        
        self.pdf = ROOT.RooAddPdf ( name , '(%s)+(%s)' % (  pdf1.name , pdf2.name ) ,
                                    self.pdf1.pdf ,
                                    self.pdf2.pdf ,
                                    self.fraction )
        
        if self.pdf1.pdf.canBeExtended() : self.error ("``pdf1'' can be extended!") 
        if self.pdf2.pdf.canBeExtended() : self.error ("``pdf2'' can be extended!") 

        self.config = {
            'pdf1'     : self.pdf1 ,
            'pdf2'     : self.pdf2 ,
            'xvar'     : self.xvar ,
            'name'     : self.name , 
            'fraction' : self.fraction 
            }

    @property
    def pdf1 ( self ) :
        """``pdf1'' : the first PDF"""
        return self.__pdf1
    
    @property
    def pdf2 ( self ) :
        """``pdf2'' : the second PDF"""
        return self.__pdf2

    @property
    def fraction ( self ) :
        """``fraction'' : the fraction of the first PDF in the sum"""
        return self.__fraction
    @fraction.setter
    def fraction ( self , value ) :
        val = float ( value )
        self.__fraction.setVal ( val )

    @property
    def F  ( self ) :
        """``F'' : the fratcion of the first PDF in the sum (the same  as ``fraction'')"""
        return self.__fraction
    @F.setter
    def F ( self , value ) :
        self.fraction = value 
        

# =============================================================================
## simple convertor of 1D-histogram into PDF
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H1D_pdf(H1D_dset,PDF) :
    """Simple convertor of 1D-histogram into PDF
    """
    def __init__ ( self            ,
                   name            ,
                   histo           ,
                   xvar    = None  ,
                   density = False ,
                   silent  = False ) :
        
        H1D_dset.__init__ ( self , histo , xvar , density , silent )
        PDF     .__init__ ( self , name  , self.xaxis ) 
        
        with roo_silent ( silent ) : 
            #
            ## finally create PDF :
            self.__vset = ROOT.RooArgSet  ( self.xvar )        
            self.pdf    = ROOT.RooHistPdf (
                'hpdf_%s'             % name ,
                'Histo1PDF(%s/%s/%s)' % ( name , histo.GetName() , histo.GetTitle() ) , 
                self.__vset , 
                self.dset   )
            
        ## and declare it be be a "signal"
        self.signals.add ( self.pdf ) 

        ## save the configuration
        self.config = {
            'name'    : self.name    , 
            'histo'   : self.histo   , 
            'xvar'    : self.xvar    , 
            'density' : self.density , 
            'silent'  : self.silent  ,             
            }
        
# =============================================================================
## @class Fit1D
#  The actual model for 1D-mass fits
#  @param signal              PDF for 'signal'     component                 (Ostap/PDF or RooAbsPdf)
#  @param background          PDF for 'background' component                 (Ostap/PDF or RooAbsPdf)
#  @param othersignals        list of PDFs for other 'signal' components     (Ostap/PDF or RooAbsPdf)
#  @param otherbackgrouds     list of PDFs for other 'background' components (Ostap/PDF or RooAbsPdf)
#  @param others              list of 'other' components                     (Ostap/PDF or RooAbsPdf)
#  @param name                the name of compound PDF 
#  @param suffix              ... add this  suffix for the PDF name
#  @param extended            build 'extended' PDF
#  @param combine_signals     combine all signal components into single SIGNAL?
#  @param combine_backgrounds combine all background components into single BACKGROUND?
#  @param combine_others      combine all other components into single COMPONENT?
#  @param recirsive           use recursive fractions for compound PDF
#  @param xvar                the fitting variable, must be specified if components are given as RooAbsPdf
#  @code 
#  gauss = Gauss_pdf( ... ) 
#  pdf   = Fit1D ( signal = gauss , background = 0 ) ## Gauss as signal and exponent as background
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-08-02
class Fit1D (PDF) :
    """The actual fit-model for generic 1D-fits
    Parameters
    - signal              : PDF for 'signal'     component                 (Ostap/PDF or RooAbsPdf)
    - background          : PDF for 'background' component                 (Ostap/PDF or RooAbsPdf)
    - othersignals        : list of PDFs for other 'signal' components     (Ostap/PDF or RooAbsPdf)
    - otherbackgrouds     : list of PDFs for other 'background' components (Ostap/PDF or RooAbsPdf)
    - others              : list of 'other' components                     (Ostap/PDF or RooAbsPdf)
    - name                : The name of compound PDF 
    - suffix              : ... add this  suffix for the PDF name
    - extended            : build 'extended' PDF
    - combine_signals     : combine all signal components into single SIGNAL?
    - combine_backgrounds : combine all background components into single BACKGROUND?
    - combine_others      : combine all other components into single COMPONENT?
    - recursive           : use recursive fractions for compound PDF
    - xvar                : the fitting variable, must be specified if components are given as RooAbsPdf

    >>> gauss = Gauss_pdf( ... ) 
    >>> pdf   = Fit1D ( signal = gauss , background = 0 ) ## Gauss as signal and exponent as background 
    """
    def __init__ ( self                          , 
                   signal              = None    ,    ## the main signal 
                   background          = None    ,    ## the main background 
                   othersignals        = []      ,    ## additional signal         components
                   otherbackgrounds    = []      ,    ## additional background     components
                   others              = []      ,    ## additional non-classified components
                   signals             = []      ,    ## alternative : all signals 
                   backgrounds         = []      ,    ## alternative : all backgrounds  
                   suffix              = ''      ,    ## the suffix 
                   name                = ''      ,    ## the name 
                   extended            = True    ,    ## extended fits ?
                   combine_signals     = False   ,    ## combine signal PDFs into single "SIGNAL"     ? 
                   combine_backgrounds = False   ,    ## combine signal PDFs into single "BACKGROUND" ?            
                   combine_others      = False   ,    ## combine signal PDFs into single "COMPONENT"  ?             
                   recursive           = True    ,    ## recursive fractions for NON-extended models?
                   xvar                = None    ,
                   S                   = []      ,    ## yields for ``signals''
                   B                   = []      ,    ## yields for ``background''
                   C                   = []      ,    ## yields for ``components''
                   F                   = []      ,    ## fractions for noin-extended fit 
                   fS                  = []      ,    ## fraction for combined signal
                   fB                  = []      ,    ## fraction for combined background
                   fC                  = []      ) :  ## fraction for combined components
        
        ##  save all arguments 
        self.__args = {
            'signal'              : signal           ,
            'othersignals'        : othersignals     ,
            'signals'             : signals          ,
            'background'          : background       ,
            'otherbackgrounds'    : otherbackgrounds ,
            'backgrounds'         : backgrounds      ,
            'others'              : others           ,
            'extended'            : extended         ,
            ##
            'suffix'              : suffix           ,
            'name'                : name             , 
            ##
            'combine_signals'     : combine_signals     ,
            'combine_backgrounds' : combine_backgrounds ,
            'combine_others'      : combine_others      ,
            'recursive'           : recursive           ,
            ##
            'xvar'                : xvar  ,
            ## 
            'S'                   : S  , ## signal     yields 
            'B'                   : B  , ## background yields       
            'C'                   : C  , ## component  yields 
            'F'                   : F  , ## fractions ( for non-extended fits)
            'fS'                  : fS , ## fractions for combined signal
            'fB'                  : fB , ## fractions for combined background
            'fC'                  : fC , ## fractions for combined component
            }
        
        self.__suffix              = suffix
        self.__recursive           = recursive
        self.__extended            = True if extended            else False
        self.__combine_signals     = True if combine_signals     else False
        self.__combine_backgrounds = True if combine_backgrounds else False
        self.__combine_others      = True if combine_others      else False

        self.__args_S = S
        self.__args_B = B
        self.__args_C = C
        self.__args_F = F
        
        assert signal or signals , "Fit1D:``signal'' or ``signals'' must be specified!"
        
        if signals :
            if signal       : raise AttributeError ( "Fit1D: ``signal''       specified for valid ``signals''!" )
            if othersignals : raise AttributeError ( "Fit1D: ``othersignals'' specified for valid ``signals''!" )
            signal       = signals [ 0   ]
            othersignals = signals [ 1 : ]
            self.debug ( 'Split signals %s into %s and %s' % ( signals , signal , othersignals  ) )                                                                   
            signals     = [] 
            
        if backgrounds :
            if background       : raise AttributeError ( "Fit1D: ``background''       specified for valid ``backgrounds''!")
            if otherbackgrounds : raise AttributeError ( "Fit1D: ``otherbackgrounds'' specified for valid ``backgrounds''!")
            background       = backgrounds [ 0   ]
            otherbackgrounds = backgrounds [ 1 : ]
            self.debug ( 'Split backgrounds %s into %s and %s' % ( backgrounds , background , otherbackgrounds ) )                                                                   
            backgrounds      = [] 
            
        ## wrap signal if needed 
        if   isinstance ( signal , PDF )                     : self.__signal = signal ## .clone() 
        ## if bare RooFit pdf,  fit variable must be specified
        elif isinstance ( signal , ROOT.RooAbsPdf ) and xvar :
            self.__signal = Generic1D_pdf ( signal , xvar , prefix = 'S_' , suffix = suffix )
        else :
            raise AttributeError ( "Fit1D:Invalid type for ``signal'': %s/%s"  % (  signal , type( signal ) ) )
        
        if not name :
            name = '%s' % self.__signal.name 
            if suffix : name += '_' + suffix 

        ## Init base class
        PDF.__init__ ( self , name , self.__signal.xvar )             
        
        ## create the background component 
        self.__background = self.make_bkg ( background , 'Bkg_' + self.name , self.xvar )

        ##  keep the lists of signals and backgrounds 
        self.signals     .add ( self.__signal     .pdf )
        self.backgrounds .add ( self.__background .pdf )

        #
        ## treat additional signals
        #        
        self.__more_signals       = [] 
        for i , c in enumerate ( othersignals ) :
            if   isinstance ( c , PDF            ) : cc = c 
            elif isinstance ( c , ROOT.RooAbsPdf ) : cc = Generic1D_pdf ( c ,  self.xvar , prefix = 'S%d_' % i , suffix = suffix ) 
            else :
                self.error ('unknown signal component %s/%s, skip it!' % ( c , type ( c ) ) )
                continue  
            self.__more_signals.append ( cc     )
            self.signals.add           ( cc.pdf ) 
        #
        ## treat additional backgounds 
        #
        self.__more_backgrounds   = [] 
        for i, c in enumerate ( otherbackgrounds ) :
            if   isinstance ( c , PDF            ) : cc = c  
            elif isinstance ( c , ROOT.RooAbsPdf ) : cc = Generic1D_pdf ( cs ,  self.xvar , prefix = 'B%d_' % i , suffix = suffix ) 
            else :
                self.error ('unknown background component %s/%s, skip it!' % ( c , type ( c ) ) )
                continue  
            self.__more_backgrounds.append ( cc     )
            self.backgrounds.add           ( cc.pdf ) 
        #
        ## treat additional components
        #
        self.__more_components    = []
        for i , c in enumerate ( others ) : 
            if   isinstance ( c , PDF            ) : cc = c  
            elif isinstance ( c , ROOT.RooAbsPdf ) : cc = Generic1D_pdf ( cs ,  self.xvar , prefix = 'C%d_' % i , suffix = suffix ) 
            else :
                self.error ("unknown ``other''component %s/%s, skip it!" % ( c , type ( c ) ) )
                continue  
            self.__more_components.append ( cc     )
            self.components.add           ( cc.pdf ) 

        # =====================================================================
        ## build PDF
        # =====================================================================
        
        self.__all_signals     = self.signals     
        self.__all_backgrounds = self.backgrounds 
        self.__all_components  = self.components  
        
        self.__save_signal     = self.__signal
        self.__save_background = self.__background
        
        ## combine signal components into single signal  (if needed)
        self.__signal_fractions = ()  
        if combine_signals and 1 < len( self.signals ) :            
            sig , fracs , sigs = self.add_pdf ( self.signals          ,
                                                'signal_'    + suffix ,
                                                'signal(%s)' % suffix ,
                                                'fS%s_%%d'   % suffix ,
                                                'fS%s_%%d'   % suffix ,
                                                recursive = True      ,
                                                fractions = fS        )
            ## new signal
            self.__signal      = Generic1D_pdf   ( sig , self.xvar , prefix = 'SIGNAL_' , suffix = suffix )
            self.__all_signals = ROOT.RooArgList ( sig )
            self.__sigs        = sigs 
            self.__signal_fractions = fracs
            self.verbose('%2d signals     are combined into single SIGNAL'     % len ( sigs ) )
            self.combined_signals.add ( sig ) 

        ## combine background components into single backhround (if needed ) 
        self.__background_fractions = () 
        if combine_backgrounds and 1 < len( self.backgrounds ) :            
            bkg , fracs , bkgs = self.add_pdf ( self.backgrounds          ,
                                                'background_'    + suffix ,
                                                'background(%s)' % suffix ,
                                                'fB%s_%%d'       % suffix ,
                                                'fB%s_%%d'       % suffix ,
                                                recursive = True          ,
                                                fractions = fB            )
            ## new background
            self.__background      = Generic1D_pdf   ( bkg , self.xvar , prefix = 'BACKGROUND_' , suffix =  suffix )
            self.__all_backgrounds = ROOT.RooArgList ( bkg )
            self.__bkgs            = bkgs 
            self.__background_fractions = fracs 
            self.verbose ('%2d backgrounds are combined into single BACKGROUND' % len ( bkgs ) ) 
            self.combined_backgrounds.add ( bkg ) 

        ## combine other components into single component (if needed ) 
        self.__components_fractions = () 
        if combine_others and 1 < len( self.components ) :
            
            cmp , fracs , cmps = self.add_pdf ( self.components      ,
                                                'other_'    + suffix ,
                                                'other(%s)' % suffix ,
                                                'fC%s_%%d'  % suffix ,
                                                'fC%s_%%d'  % suffix ,
                                                recursive = True     ,
                                                fractions = fC       ) 
            ## save old background
            self.__other          = Generic1D_pdf   ( cmp , self.xvar , prefix = 'COMPONENT_' , suffix = suffix )
            self.__all_components = ROOT.RooArgList ( cmp )
            self.__components_fractions = fracs 
            self.verbose('%2d components  are combined into single COMPONENT'    % len ( cmps ) )
            self.combined_components.add ( cmp ) 


        self.__nums_signals     = [] 
        self.__nums_backgrounds = [] 
        self.__nums_components  = []
        self.__nums_fractions   = []
        
        ## build models 
        if self.extended :

            if F : self.warning("Non empty list of ``fractions'' is specified: %s, ignore" % F ) 

            ns = len ( self.__all_signals )
            if 1 == ns :
                sf = self.make_var  ( get_i ( S , 0 ) , "S"+suffix , "Signal"     + suffix , None , 1 , 0 , 1.e+7 )
                self.alist1    .add ( self.__all_signals[0]  )
                self.__nums_signals.append ( sf ) 
            elif 2 <= ns : 
                fis = self.make_fracs ( ns , 'S%s_%%d' % suffix ,  'S%s_%%d'  % suffix , fractions  = False , fracs = S )
                for s in self.__all_signals : self.alist1.add ( s )
                for f in fis                : self.__nums_signals.append ( f ) 

            nb = len ( self.__all_backgrounds )
            if 1 == nb :
                bf = self.make_var ( get_i ( B , 0 ) , "B"+suffix , "Background" + suffix , None , 1 , 0 , 1.e+7 )
                self.alist1.add ( self.__all_backgrounds[0]  )
                self.__nums_backgrounds.append ( bf ) 
            elif 2 <= nb :
                fib = self.make_fracs ( nb , 'B%s_%%d' % suffix ,  'B%s_%%d'  % suffix , fractions  = False , fracs = B )
                for b in self.__all_backgrounds : self.alist1.add ( b )
                for f in fib                    : self.__nums_backgrounds.append ( f ) 

            nc = len ( self.__all_components )
            if 1 == nc :
                cf = self.make_var ( get_i ( C , 0 )  , "C"+suffix , "Component" + suffix , None , 1 , 0 , 1.e+7 )
                self.alist1.add  ( self.__all_components[0]  )
                self.__nums_components.append ( cf ) 
            elif 2 <= nc : 
                fic = self.make_fracs ( nc , 'C%s_%%d' % suffix ,  'C%s_%%d'  % suffix , fractions  = False , fracs = C )
                for c in self.__all_components : self.alist1.add ( c )
                for f in fic                   : self.__nums_components.append ( f )

            for s in self.__nums_signals     : self.alist2.add ( s ) 
            for b in self.__nums_backgrounds : self.alist2.add ( b ) 
            for c in self.__nums_components  : self.alist2.add ( c ) 
                    
        else :

            if S : self.warning("Non empty list of ``signals''     is specified: %s, ignore" % S ) 
            if C : self.warning("Non empty list of ``components''  is specified: %s, ignore" % C ) 
            if B : self.warning("Non empty list of ``backgrounds'' is specified: %s, ignore" % B ) 

            ns = len ( self.__all_signals     )
            nb = len ( self.__all_backgrounds )
            nc = len ( self.__all_components  )
            
            for s in self.__all_signals     : self.alist1.add ( s )
            for b in self.__all_backgrounds : self.alist1.add ( b )
            for c in self.__all_components  : self.alist1.add ( c )

            fic = self.make_fracs ( ns + nb + nc ,
                                    'f%s_%%d' % suffix          ,
                                    'f%s_%%d' % suffix          ,
                                    fractions  = True           ,
                                    recursive  = self.recursive ,
                                    fracs      = F              )
                
            for f in fic                    : self.__nums_fractions.append ( f )   
            for f in self.__nums_fractions  : self.alist2.add ( f ) 


        self.__nums_signals     = tuple ( self.__nums_signals     )
        self.__nums_backgrounds = tuple ( self.__nums_backgrounds ) 
        self.__nums_components  = tuple ( self.__nums_components  ) 
        self.__nums_fractions   = tuple ( self.__nums_fractions   ) 

        #
        ## The final PDF
        #       

        pdfname  = "Fit1D_"    + self.name
        pdftitle = "Fit1D(%s)" % self.name
        pdfargs  = pdfname , pdftitle , self.alist1 , self.alist2
        
        if not self.extended :
            pdfargs = pdfargs + ( True if recursive else False , ) ## RECURSIVE ? 
        self.pdf = ROOT.RooAddPdf ( *pdfargs )
        
        if self.extended : 
            self.debug ( "extended     model ``%s'' with %s/%s components"  % ( self.pdf.GetName() , len( self.alist1) , len(self.alist2) ) )
        else : 
            self.debug ( "non-extended model ``%s'' with %s/%s components"  % ( self.pdf.GetName() , len( self.alist1) , len(self.alist2) ) )

        ## save the configuration
        self.config = {
            'signal'              : self.save_signal         ,
            'background'          : self.save_background     ,
            'othersignals'        : self.more_signals        ,
            'otherbackgrounds'    : self.more_backgrounds    ,
            'others'              : self.more_components     ,
            'suffix'              : self.suffix              ,
            'name'                : self.name                ,
            'extended'            : self.extended            ,
            'combine_signals'     : self.combine_signals     ,
            'combine_backgrounds' : self.combine_backgrounds ,
            'combine_others'      : self.combine_others      ,
            'recursive'           : self.recursive           ,
            'xvar'                : self.xvar                ,
            'S'                   : self.S                   ,
            'B'                   : self.B                   ,
            'C'                   : self.C                   ,
            'F'                   : self.F                   ,            
            'fS'                  : self.fS                  ,
            'fB'                  : self.fB                  ,
            'fC'                  : self.fC                  ,
            }

        self.checked_keys.add  ( 'xvar' )
        
    @property
    def extended ( self ) :
        """``extended'': build extended PDF?"""
        return  self.__extended 
    @property
    def suffix   ( self ) :
        """``suffix'' : append the names  with the specified suffix"""
        return self.__suffix
    @property
    def recursive ( self ) :
        """``recursive'':  use recursive fitfractions?"""
        return  self.__recursive
    @property
    def combine_signals ( self ) :
        """Combine all ``signal''-components into single ``signal'' componet?"""
        return self.__combine_signals
    @property
    def combine_backgrounds ( self ) :
        """Combine all ``background''-components into single ``background'' componet?"""
        return self.__combine_backgrounds
    @property
    def combine_others ( self ) :
        """Combine all ``others''-components into single ``other'' componet?"""
        return self.__combine_others 
                
    @property
    def signal (  self ) :
        """The main ``signal'' component (possible compound)"""
        return self.__signal
    
    @property
    def background (  self ) :
        """The main ``background'' component (possible compound)"""
        return self.__background

    @property
    def save_signal (  self ) :
        """The original ``signal'' component (possible compound)"""
        return self.__save_signal
    
    @property
    def save_background (  self ) :
        """The original ``background'' component (possible compound)"""
        return self.__save_background

    @property
    def more_signals     ( self ) :
        """Additional ``signal'' components"""
        return tuple( self.__more_signals )
    
    @property
    def more_backgrounds ( self ) :
        """additional ``background'' components"""
        return tuple( self.__more_backgrounds  )
    
    @property
    def more_components ( self ) :
        """additional ``other'' components"""
        return tuple( self.__more_components  )
    
    @property
    def fS ( self  ) :
        """(Recursive) fractions for the compound signal components (empty for simple signal) """
        lst = [ i for i in self.__signal_fractions ]
        return tuple ( lst )
    @fS.setter
    def fS ( self , value ) :

        if   isinstance ( value , num_types          ) : value = [ value           ]
        elif isinstance ( value , VE                 ) : value = [ value.value()   ]
        elif isinstance ( value , ROOT.RooAbsReal    ) : value = [ float ( value ) ] 
        elif isinstance ( value , list_types         ) : pass
        elif isinstance ( value , ROOT.RooArgList    ) : pass

        for f , v in zip ( self.__signal_fractions , value ) :
            vv = float ( v )
            if f.minmax() and not vv in f :
                logger.error ("Value %s is outside the allowed region %s"  % ( vv , f.minmax() ) ) 
            f.setVal   ( vv ) 
            
    @property
    def fB ( self  ) :
        """(Recursive) fractions for the compound background components (empty for simple background)"""
        lst = [ i for i in self.__background_fractions ]
        return tuple ( lst )
    @fB.setter
    def fB ( self , value ) :

        if   isinstance ( value , num_types          ) : value = [ value           ]
        elif isinstance ( value , VE                 ) : value = [ value.value()   ]
        elif isinstance ( value , ROOT.RooAbsReal    ) : value = [ float ( value ) ] 
        elif isinstance ( value , list_types         ) : pass
        elif isinstance ( value , ROOT.RooArgList    ) : pass

        for f , v in zip ( self.__background_fractions , value ) :
            vv = float ( v )
            if f.minmax() and not vv in f :
                logger.error ("Value %s is outside the allowed region %s"  % ( vv , f.minmax() ) ) 
            f.setVal   ( vv ) 
                
    @property
    def fC ( self  ) :
        """(Recursive) fractions for the compound ``other'' components (empty for no additional commponents case)"""
        lst = [ i for i in self.__components_fractions ]
        return tuple ( lst )
    @fC.setter
    def fC ( self , value ) :

        if   isinstance ( value , num_types          ) : value = [ value           ]
        elif isinstance ( value , VE                 ) : value = [ value.value()   ]
        elif isinstance ( value , ROOT.RooAbsReal    ) : value = [ float ( value ) ] 
        elif isinstance ( value , list_types         ) : pass
        elif isinstance ( value , ROOT.RooArgList    ) : pass

        for f , v in zip ( self.__components_fractions , value ) :
            vv = float ( v )
            if f.minmax() and not vv in f :
                logger.error ("Value %s is outside the allowed region %s"  % ( vv , f.minmax() ) ) 
            f.setVal   ( vv ) 
            
    @property
    def S ( self ) :
        """Get the  yields of signal component(s) (empty for non-extended fits)
        For single signal component:
        >>> print pdf.S          ## read the single single component 
        >>> pdf.S = 100          ## assign to it
        For multiple signal components:
        >>> print pdf.S[4]       ## read the 4th signal component 
        >>> pdf.S = (1,2,3,4,5,6)## assign to it 
        ... or, alternatively:
        >>> print pdf.S[4]       ## read the 4th signal component 
        >>> pdf.S[4].value = 100 ## assign to it         
        """
        lst = [ i for i in self.__nums_signals ]
        if not lst          : return ()     ## extended fit? 
        elif  1 == len(lst) : return lst[0] ## simple signal?
        return tuple ( lst )
    @S.setter
    def S (  self , value ) :
        
        ns = len ( self.__nums_signals )
        assert 1 <= ns , "No signals are defined, assignement is impossible"
        
        ##
        if   isinstance ( value , num_types          ) : value = [ value           ]
        elif isinstance ( value , VE                 ) : value = [ value.value()   ]
        elif isinstance ( value , ROOT.RooAbsReal    ) : value = [ float ( value ) ] 
        elif isinstance ( value , list_types         ) : pass
        elif isinstance ( value , ROOT.RooArgList    ) : pass

        ss = [ self.S ] if 1 == ns else self.S

        for s , v in zip ( ss , value ) :

            vv = float ( v  )
            if s.minmax() and not vv in s :
                logger.error ("Value %s is outside the allowed region %s"  % ( vv , s.minmax() ) ) 
            s.setVal   ( vv ) 
    
    @property
    def B ( self ) :
        """Get the  yields of background  component(s) (empty for non-extended fits)
        For single background component:
        >>> print pdf.B          ## read the single background component 
        >>> pdf.B = 100          ## assign to it 
        For multiple background components:
        >>> print pdf.B[4]            ## read the 4th background component 
        >>> pdf.B = ( 1, 2, 3, 4, 5 ) ## assign to it 
        ... or, alternatively:
        >>> print pdf.B[4]       ## read the 4th background component 
        >>> pdf.B[4].value = 100 ## assign to it 
        """
        lst = [ i for i in self.__nums_backgrounds ]
        if not lst          : return ()     ## extended fit? 
        elif  1 == len(lst) : return lst[0] ## simple background?
        return tuple ( lst )
    @B.setter
    def B (  self , value ) :
        
        nb = len ( self.__nums_backgrounds )
        assert 1 <= nb , "No backgrounds are defined, assignement is impossible"

        if   isinstance ( value , num_types          ) : value = [ value           ]
        elif isinstance ( value , VE                 ) : value = [ value.value()   ]
        elif isinstance ( value , ROOT.RooAbsReal    ) : value = [ float ( value ) ] 
        elif isinstance ( value , list_types         ) : pass
        elif isinstance ( value , ROOT.RooArgList    ) : pass

        ss = [ self.B ] if 1 == nb else self.B

        for s , v in zip ( ss , value ) :

            vv = float ( v  )
            if s.minmax() and not vv in s :
                logger.error ("Value %s is outside the allowed region %s"  % ( vv  , s.minmax() ) ) 
            s.setVal   ( vv ) 

    @property
    def C ( self ) :
        """Get the  yields of ``other'' component(s) (empty for non-extended fits)
        For single ``other'' component:
        >>> print pdf.C           ## read the single ``other'' component 
        >>> pdf.C = 100           ## assign to it 
        For multiple ``other'' components:
        >>> print pdf.C[4]            ## read the 4th ``other'' component 
        >>> pdf.C = ( 1, 2, 3, 4, 5 ) ## assign to it 
        ... or, alternatively:
        >>> print pdf.C[4]        ## read the 4th ``other'' component 
        >>> pdf.C[4].value 100    ## assign to it         
        """
        lst = [ i for i in self.__nums_components ]
        if not lst          : return ()     ## extended fit? no other components?
        elif  1 == len(lst) : return lst[0] ## single component?
        return tuple ( lst )
    @C.setter
    def C (  self , value ) :
        
        nc = len ( self.__nums_components )
        assert 1 <= nc , "No ``other'' components are defined, assignement is impossible"

        if   isinstance ( value , num_types          ) : value = [ value           ]
        elif isinstance ( value , VE                 ) : value = [ value.value()   ]
        elif isinstance ( value , ROOT.RooAbsReal    ) : value = [ float ( value ) ] 
        elif isinstance ( value , list_types         ) : pass
        elif isinstance ( value , ROOT.RooArgList    ) : pass

        ss = [ self.C ] if 1 == nc else self.C

        for s , v in zip ( ss , value ) :

            vv = float ( v  )
            if s.minmax() and not vv in s :
                logger.error("Value %s is outside the allowed region %s"  % ( vv , s.minmax() ) ) 
            s.setVal   ( vv ) 

    @property 
    def F ( self ) :
        """Get fit fractions for non-expended fits (empty for extended fits)
        For single fraction (2 fit components):
        >>> print pdf.F           ## read the single fraction 
        >>> pdf.F = 0.1           ## assign to it 
        For multiple fractions (>2 fit components):
        >>> print pdf.F[4]        ## read the 4th fraction
        >>> pdf.F = (0.1,0.2,0.3,0.4,0.6) ## assign to it 
        ... or, alternatively:
        >>> print pdf.F[4]        ## read the 4th fraction
        >>> pdf.F[4].value = 0.1  ## assign to it         
        """
        lst = [ i for i in self.__nums_fractions ]
        if not lst          : return ()     ## extended fit? 
        elif  1 == len(lst) : return lst[0] ## simple two component fit ?
        return tuple ( lst )
    @F.setter
    def F (  self , value ) :
        nf = len ( self.__nums_fractions )
        assert 1 <= nf , "No fractions are defined, assignement is impossible"

        ss = [ self.F ] if 1 == nf else self.F
        if   isinstance ( value , num_types          ) : value = [ value           ]
        elif isinstance ( value , VE                 ) : value = [ value.value()   ]
        elif isinstance ( value , ROOT.RooAbsReal    ) : value = [ float ( value ) ] 
        elif isinstance ( value , list_types         ) : pass
        elif isinstance ( value , ROOT.RooArgList    ) : pass

        for s , v in zip ( ss , value ) :

            vv = float ( v  )
            if s.minmax() and not vv in s :
                logger.error ("Value %s is outside the allowed region %s"  % ( vv , s.minmax() ) ) 
            s.setVal   ( vv ) 

    @property
    def  yields    ( self ) :
        """The list/tuple of the yields of all numeric components (empty for non-extended fit)"""
        return tuple ( [ i for i in  self.alist2 ] ) if     self.extended else ()
    
    def total_yield ( self ) :
        """``total_yield''' : get the total yield"""
        if not self.extended    : return None 
        if not self.fit_result                                 : return None
        if not valid_pointer ( self.fit_result )               : return None
        yields = self.yields
        if not yields                                          : return None
        if 1 ==  len ( yields )                                : return yields[0].value  
        return self.fit_result.sum ( *yields ) 
 
    @property
    def  fractions ( self ) :
        """The list/tuple of fit fractions of all numeric components (empty for extended fit)"""
        return tuple ( [ i for i in  self.alist2 ] ) if not self.extended else () 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
# The END 
# =============================================================================
