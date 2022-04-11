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
    'MASSMEAN'      , ## useful base class to create "signal" PDFs for mass-fits
    'MASS'          , ## useful base class to create "signal" PDFs for mass-fits
    'RESOLUTION'    , ## useful base class to create "resolution" PDFs
    ##
    'Fit1D'         , ## the basic compound 1D-fit model 
    ##
    'Flat1D'        , ## trivial 1D-pdf: constant 
    'Generic1D_pdf' , ## wrapper over imported RooFit (1D)-pdf
    'Sum1D'         , ## wrapper for RooAddPdf 
    'H1D_pdf'       , ## convertor of 1D-histo to RooHistPdf
    'Shape1D_pdf'   , ## simple PDF from C++ shape 
    'make_pdf'      , ## helper function to make PDF
    'all_args'      , ## check that all arguments has correct type 
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
                                        list_types     , all_numerics   ) 
from   ostap.fitting.roofit    import SETVAR, FIXVAR, PDF_fun
from   ostap.logger.utils      import roo_silent   , rootWarning
from   ostap.fitting.utils     import ( RangeVar   , MakeVar   , numcpu   ,
                                        Phases     , make_name ,  
                                        fit_status , cov_qual  , H1D_dset , get_i )
from   ostap.utils.utils       import make_iterable 
from   ostap.fitting.funbasic  import FUNC,  SETPARS 
from   ostap.utils.cidict      import select_keys
from   ostap.fitting.roocmdarg import check_arg , nontrivial_arg , flat_args , command  
import ostap.histos.histos 
from   ostap.core.meta_info    import root_info
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.basic' )
else                       : logger = getLogger ( __name__              )
# =============================================================================

# =============================================================================
## @var arg_types
#  list of "good" argument  types 
arg_types = num_types + ( VE , ROOT.RooAbsReal )
# =============================================================================
## are all args of "good" type?
#  - ROOT.RooAbsReal
#  - numeric type
#  - VE
#  - tuple of 1-3 arguments of numeric values 
def all_args ( *args ) :
    """Are all arguments of ``good'' type?
    - ROOT.RooAbsReal
    - numeric type
    - VE
    - tuple of 1-3 arguments of numeric values 
    """
    ## try to find 
    for a in args :
        if   isinstance ( a , arg_types ) : pass
        elif isinstance ( a , tuple     ) and \
             1 <= len ( a ) <= 3          and \
             all_numerics ( *a )          : pass
        else                              : return False 
        
    return True 

# =============================================================================
## @class PDF
#  The helper base class for implementation of various PDF-wrappers 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-08-21
class PDF (FUNC) :
    """Useful helper base class for implementation of various PDF-wrappers 
    """
    def __init__ ( self , name ,  xvar , special = False , **kwargs ) :

        FUNC.__init__  ( self , name , xvar = xvar , **kwargs ) 

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
        elif hasattr ( value , 'dset' ) and isinstance ( value.dset , ROOT.RooAbsData ) :
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
        opts     = self.fit_options + ( ROOT.RooFit.Save () , ) + args 
        opts     = self.parse_args ( dataset , *opts , **kwargs )
        
        if silent :
            pl = check_arg ('PrintLevel'      , *opts ) 
            if not pl : opts = opts + ( ROOT.RooFit.PrintLevel      ( -1    ) , )
            vl = check_arg ('Verbose'         , *opts )
            if not vl : opts = opts + ( ROOT.RooFit.Verbose         ( False ) , )
            pe = check_arg ('PrintEvalErrors' , *opts )
            if not pe : opts = opts + ( ROOT.RooFit.PrintEvalErrors ( 0     ) , )
                
        weighted = dataset.isWeighted() if dataset else False
        if weighted :
            sw2 = check_arg ( 'SumW2Error'      , *opts )
            aer = check_arg ( 'AsymptoticError' , *opts )
            if not sw2 and not aer :
                self.warning ( "fitTo: Neither ``SumW2Error'' and ``AsymptoticError'' are specified for weighted dataset!" )

        if 1 < len ( self.vars ) :
            rng = check_arg ( 'Range' , *opts ) 
            if rng : self.warning ( 'fitTo: %s is specified for >1D function - it is ambuguous!' % rng )

        ## check fit ranges 
        rng = check_arg ( 'RangeByName' , *opts )
        ok  = self.check_ranges ( dataset , rng.getString(0) if rng else '' )
        if not ok : self.warning ( 'fitTo: ranges are not OK' ) 

        if not silent and opts and nontrivial_arg ( ( 'Save' , 'NumCPU' ) , *opts ) :
            self.info ('fitTo options: %s ' % list ( flat_args ( *opts ) ) ) 

        ## play a bit with the binning cache for 1D-convolutions 
        if self.xvar.hasBinning ( 'cache' ) :
            nb1 = self.xvar.getBins( 'cache' ) 
            xv  = getattr ( dataset , self.xvar.name , None )
            if   xv and xv.hasBinning ( 'cache' ) :
                nb2 = xv.getBins('cache')
                if  nb1 != nb2 :
                    xv.setBins ( max (  nb1 , nb2 ) , 'cache' )
                    if not silent :
                        self.info ('Adjust binning cache %s->%s for variable %s in dataset' % ( nb2 , nb1 , xv.name ) )
            elif xv :
                xv.setBins (        nb1         , 'cache' )
                if not silent : 
                    self    .info ('Set binning cache %s for variable %s in dataset' %  ( nb1 , xv.name )  )

        ## define silent context
        with roo_silent ( silent ) :
            self.fit_result = None            
            result          = self.fit_to ( self.pdf , dataset , *opts ) 
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
        ## check the integrals (if/when possible)
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
                        
            if hasattr ( self , 'natural' ) and self.natural and not dataset.isWeighted () :

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

        cov2_good = ( qual == 3 ) or ( dataset.isWeighted() and qual == -1 )
        
        if   result and 0 == result.status() and not silent :
            logger.info     ( "Fit result is\n%s" % result.table ( prefix = "# " ) ) 
        elif result and ( not cov2_good ) and not silent : 
            logger.warning  ( "Fit result is\n%s" % result.table ( prefix = "# " ) ) 
        elif result and not silent :
            logger.warning  ( "Fit result is\n%s" % result.table ( prefix = "# " ) )

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

        ## ##
                         
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

        for i , cmp in enumerate ( what ) :

            st         = style  ( i ) if style and callable  ( style ) else ()
            
            component  = ROOT.RooFit.Components ( cmp.name )

            atup = args + tuple ( options ) + tuple ( st ) 
            self.debug   ( 'drawing component %s with options %s' % ( cmp.name , ( component, ) + atup ) )             
            self.plot_on ( self.pdf , frame , component , *atup  )
            
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
        used_options = set() 
        
        with roo_silent ( silent ) , useStyle ( style ) :

            drawvar = drawvar if drawvar else ( self.draw_var if self.draw_var else self.xvar )  

            binned = dataset and isinstance ( dataset , ROOT.RooDataHist )

            if nbins :  frame = drawvar.frame ( nbins )            
            else     :  frame = drawvar.frame ()

            #
            ## draw invizible data (for normalzation of fitting curves)
            #
            data_options = self.draw_option ( 'data_options' , **kwargs )
            used_options.add ( 'data_options' ) 
            if dataset and dataset.isWeighted() and dataset.isNonPoissonWeighted() : 
                data_options = data_options + ( ROOT.RooFit.DataError( ROOT.RooAbsData.SumW2 ) , )
                
            if dataset and binned and nbins :
                data_options = data_options + ( ROOT.RooFit.Binning ( nbins ) , )

            if dataset :

                commands = data_options 
                commands = data_options + args +  ( ROOT.RooFit.Invisible() , ) 
                self.plot_on ( dataset , frame , *commands ) 
                
            ## draw various ``background'' terms
            boptions     = self.draw_option ( 'background_options' , **kwargs ) 
            bbstyle      = self.draw_option (   'background_style' , **kwargs )
            self._draw( self.backgrounds , frame , boptions , bbstyle )
            used_options.add ( 'background_options' ) 
            used_options.add ( 'background_style'   ) 

            ## draw combined ``background'' components 
            if self.combined_backgrounds :
                
                drawit   = self.draw_option ( 'draw_combined_background'    , **kwargs )
                doptions = self.draw_option ( 'combined_background_options' , **kwargs ) 
                dstyle   = self.draw_option ( 'combined_background_style'   , **kwargs )
                
                if drawit : self._draw ( self.combined_backgrounds , frame , doptions , dstyle , args )
                
            used_options.add ( 'draw_combined_background'    ) 
            used_options.add ( 'combined_background_options' ) 
            used_options.add ( 'combined_background_style  ' )
                
            ## ugly :-(
            ct1options   = self.draw_option ( 'crossterm1_options' , **kwargs )
            ct1bstyle    = self.draw_option ( 'crossterm1_style'   , **kwargs ) 
            if hasattr ( self , 'crossterms1' ) and self.crossterms1 : 
                self._draw( self.crossterms1 , frame , ct1options , ct1bstyle , args )
                
            used_options.add ( 'crossterm1_options' ) 
            used_options.add ( 'crossterm1_style'   ) 

            ## ugly :-(
            ct2options   = self.draw_option ( 'crossterm2_options' , **kwargs )
            ct2bstyle    = self.draw_option ( 'crossterm2_style'   , **kwargs )
            
            if hasattr ( self , 'crossterms2' ) and self.crossterms2 :
                self._draw( self.crossterms2 , frame , ct2options , ct2bstyle , args )
                
            used_options.add ( 'crossterm2_options' ) 
            used_options.add ( 'crossterm2_style'   ) 

            ## draw ``other'' components
            coptions     = self.draw_option ( 'component_options' , **kwargs )
            cbstyle      = self.draw_option ( 'component_style'   , **kwargs )
            self._draw( self.components , frame , coptions , cbstyle , args )

            used_options.add ( 'component_options' ) 
            used_options.add ( 'component_style'   ) 

            ## draw combined ``other'' components 
            if self.combined_components :
                
                drawit   = self.draw_option ( 'draw_combined_component'    , **kwargs )
                doptions = self.draw_option ( 'combined_component_options' , **kwargs ) 
                dstyle   = self.draw_option ( 'combined_component_style'   , **kwargs )
                
                if drawit : self._draw ( self.combined_components , frame , doptions , dstyle , args )
                
            used_options.add ( 'draw_combined_component'    ) 
            used_options.add ( 'combined_component_options' ) 
            used_options.add ( 'combined_component_style'   )
            
            ## draw ``signal'' components
            soptions     = self.draw_option (    'signal_options'  , **kwargs )
            sbstyle      = self.draw_option (      'signal_style'  , **kwargs ) 
            self._draw( self.signals , frame , soptions , sbstyle , args )

            used_options.add ( 'signal_options' ) 
            used_options.add ( 'signal_style'   )

            ## draw combined ``signals'' components 
            if self.combined_signals :
                drawit    = self.draw_option ( 'draw_combined_signal'     , **kwargs )
                doptions  = self.draw_option ( 'combined_signal_options'  , **kwargs ) 
                dstyle    = self.draw_option (   'combined_signal_style'  , **kwargs )
                if drawit : self._draw ( self.combined_signals , frame , doptions , dstyle , args )
                
            used_options.add ( 'draw_combined_signal'    ) 
            used_options.add ( 'combined_signal_options' ) 
            used_options.add ( 'combined_combined_style' )
            
            #
            ## the total fit curve
            #
            totoptions   = self.draw_option (  'total_fit_options' , **kwargs )
            self.plot_on ( self.pdf , frame , *totoptions ) 
            used_options.add ( 'total_fit_options'    ) 

            #
            ## draw data once more
            #
            if dataset :
                commands = data_options + args
                self.plot_on ( dataset , frame , *commands )

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
            groot = ROOT.ROOT.GetROOT()
            if not groot.IsBatch() :
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
                used = set() 
                from ostap.plotting.fit_draw import key_compare as draw_key_compare 
                for k in kwargs :
                    for key in used_options :
                        if draw_key_compare ( k , key ) :
                            used.add ( k ) 
                extra = list ( sorted ( set ( kwargs.keys() ) - used ) )
                if extra : self.warning ( "draw: ignored unknown options: %s" % extra ) 

            ## calculate chi2/ndf
            frame.chi2dnf = None 
            if dataset and not silent :             
                pars          = self.params ( dataset )
                frame.chi2ndf = frame.chiSquare ( len ( pars ) )
                binw          = -1 
                if nbins and isinstance ( nbins , integer_types ) and 1 < nbins :
                    if hasattr ( drawvar , 'xminmax' ) and drawvar.xminmax () :
                        xmn , xmx =  drawvar.xminmax()
                        binw = ( xmx - xmn ) / float ( nbins )
                if 0 < binw : self.info ( 'chi2/ndf: %.3f, binwidth: %s' %  ( frame.chi2ndf , binw ) )
                else        : self.info ( 'chi2/ndf: %.3f' %                  frame.chi2ndf          )
                
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
                self.histo_data = H1D_dset ( histo = histo , xaxis = self.xvar , density = density , silent = silent )
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
                self.histo_data = H1D_dset ( histo = dataset , xaxis = self.xvar , density = density , silent = silent )
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
                self.histo_data   = H1D_dset ( histo = dataset , xaxis = self.xvar , density = density , silent = silent )
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
        pars = self.params ( dataset )
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
        
        with SETPARS ( self , dataset ) , KeepBinning ( var ) :

            if bins :
                var.bins = bins
            
            self.debug ( 'draw_nll: frame  args: %s'% list ( fargs ) )
            ## prepare the frame & plot 
            frame = var.frame ( *fargs )
            
            self.debug ( 'draw_nll: plotOn args: %s'% list ( largs ) )
            self.plot_on ( result , frame , *largs )
            
            import ostap.histos.graphs
            
            ## remove a bit strange drawing artefacts (the first and the last points )
            if var.minmax() :
                vmn , vmx = var.minmax()
                graph   = frame.getObject ( 0 )
                if graph.remove ( remove = lambda x,y : not vmn <= x <= vmx ) :
                    logger.debug ('draw_nll: remove drawing artefacts  at the first and the last points' ) 
                                        
            ## scale it if needed
            if 1 != sf :
                logger.info ('draw_nll: apply scale factor of %.4g due to dataset weights' % sf )
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
        groot = ROOT.ROOT.GetROOT()        
        if not groot.IsBatch() :
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
        if len ( opts ) < 8 : 
            return self.pdf.createNLL ( dataset , *opts             ) , sf
        else :
            return self.pdf.createNLL ( dataset , command ( *opts ) ) , sf
            

    # =========================================================================
    ## get NLL/profile-graph for the variable, using the specified abscissas
    #  @code
    #  pdf   = ...
    #  graph = pdf.graph_nll ( 'S'                      ,
    #                          vrange ( 0 , 100 , 100 ) ,
    #                          dataset                  )
    #  @endcode
    def graph_nll ( self             ,
                    variable         , 
                    values           ,
                    dataset          ,
                    silent   = True  ,
                    draw     = False ,
                    subtract = True  , 
                    args     = ()    , **kwargs ) :
        """Get NLL/profile-graph for the variable, using the specified abscissas
        >>> pdf   = ...
        >>> graph = pdf.graph_nll ( 'S'                     ,
        ...                          vrange ( 0 , 100 , 100 ) ,
        ...                          dataset                )
        """

        ## 1) create NLL 
        nLL, sf = self.nll ( dataset , silent = silent ,  args = args , **kwargs )

        ## get the parametrs
        var  = variable 
        pars = self.params ( dataset ) 
        assert var in pars , "Variable %s is not a parameter"   % var
        if not isinstance ( var , ROOT.RooAbsReal ) : var = pars[ var ]
        del pars 

        import ostap.histos.graphs
        ## 2) create graph if drawing reqested 
        graph = ROOT.TGraph () if draw else None 
            
        ## 3) collect NLL values 
        results   = []
        vmin      = None
        with SETPARS ( self , dataset ) , SETVAR  ( var ) :
            from ostap.utils.progress_bar import progress_bar 
            for v in progress_bar  ( values , silent = silent ) :
                var.setVal ( v )
                n   = nLL.getVal() 
                res = v , n
                results.append ( res )
                vmin = n if vmin is None else min ( vmin , n )
                if draw :
                    graph.SetPoint ( len ( graph ) , v , n ) ## add the point 
                    if 1 == len ( graph ) : graph.draw ( "ap" )  

        ## 3) create graph
        graph = ROOT.TGraph ( len ( results ) )
        results.sort ()
        ymin   = None 
        for i , point in enumerate ( results ) :
            x , y = point 
            graph [ i ]  = x , y 
            ymin = y if ymin is None else min ( ymin , y )
            
       ## subtract the minimum
        if subtract :
            logger.info ( "graph_nll: minimal value of %.5g is subtracted" % ymin ) 
            graph -= ymin 

        ## scale it if needed
        if 1 != sf :
            logger.info ('graph_nll: apply scale factor of %.5g due to dataset weights' % sf )
            graph *= sf 
            
        if draw : graph.draw ('ap')

        return graph 

    # =========================================================================
    ## get NLL-profile-graph for the variable, using the specified abscissas
    #  @code
    #  pdf   = ...
    #  graph = pdf.graph_profile ( 'S'                       ,
    #                              vrange ( 0 , 12.5 , 10  ) ,
    #                              dataset                   )
    #  @endcode
    def graph_profile ( self             ,
                        variable         , 
                        values           ,
                        dataset          ,
                        fix      = []    ,
                        silent   = True  ,
                        draw     = False ,
                        subtract = True  , 
                        args     = ()    , **kwargs ) :
        """Get profile-graph for the variable, using the specified abscissas
        >>> pdf   = ...
        >>> graph = pdf.graph_profile ( 'S'                      ,
        ...                             vrange ( 0 , 12.5 , 20 ) ,
        ...                             dataset                  )
        """
        
        vals = [ v for v in values ]
        assert vals, 'graph_profile: no points are specified!'
                
        vmin = min ( vals )
        vmax = max ( vals )

        ## 1) create NLL 
        nLL , sf = self.nll ( dataset , silent = silent ,  args = args , **kwargs )

        ## get the parametrs
        var  = variable 
        pars = self.params ( dataset ) 
        assert var in pars , "Variable %s is not a parameter"   % var
        if not isinstance ( var , ROOT.RooAbsReal ) : var = pars[ var ]

        if var.minmax () :            
            minv = min ( var.getMin () , vmin )
            maxv = max ( var.getMax () , vmax )

        vars = ROOT.RooArgSet ( var )
        for f in fix :
            fv = f
            if not isinstance ( fv , ROOT.RooAbsReal ) : fv = pars [ fv ]
            vars.add ( fv ) 
            
        ## 2)  create profile LL
        pLL = ROOT.RooProfileLL ( 'pLL_%s_%s'         % ( var.name , self.name ) ,
                                  'LL-profile(%s,%s)' % ( var.name , self.name ) ,
                                  nLL , vars )

        ## 2) create graph if requested 
        import ostap.histos.graphs
        graph = ROOT.TGraph () if draw else None 
                        

        ## 3) collect pLL values 
        results = [] 
        with SETPARS ( self , dataset ) , RangeVar ( var , minv , maxv ) , SETVAR  ( var ) :
            from ostap.utils.progress_bar import progress_bar 
            for i , v in enumerate ( progress_bar ( vals , silent = silent )  ) :
                var.setVal ( v )
                p   = pLL.getVal() 
                res = v , p 
                results.append ( res )
                if draw :
                    graph.SetPoint ( len ( graph ) , v , p ) ## add the point 
                    if 1 == len ( graph ) : graph.draw ("ap")  
                    
        ## 4) re-create the graph
        graph   = ROOT.TGraph ( len ( results ) )
        results.sort ()
        ymin    = None 
        for i , point  in enumerate ( results ) :
            x , y = point
            graph [ i ]  = x , y 
            ymin = y if ymin is None else min ( ymin , y )
                    
        ## subtract the minimum
        if subtract :
            logger.info ( "graph_profile: minimal value of %.5g is subtracted" % ymin ) 
            graph -= ymin 

        ## scale it if needed
        if 1 != sf :
            logger.info ('graph_profile: apply scale factor of %.5g due to dataset weights' % sf )
            graph *= sf 

        if draw : graph.draw ('ap')
        
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
                self.histo_data = H1D_dset ( histo = dataset , xaxis = self.xvar , density = density , silent = silent )
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
        pars = self.params ( dataset ) 
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
            
        with SETPARS ( self , dataset ) , roo_silent ( silent ) :
            
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
            if 1 != sf :  logger.info ('Scale factor of %.4g is applied' % sf )
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
                self.histo_data = H1D_dset ( histo = dataset , xaxis = self.xvar , density = density , silnet = silent )
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
        pars = self.params ( dataset ) 
        
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
        with SETPARS ( self , dataset ) , roo_silent ( silent ) :

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
            if 1 != sf :  logger.info ('Scale factor of %.4g is applied' % sf )
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
    def minuit ( self                    ,
                 dataset         = None  ,
                 max_calls       = -1    ,
                 max_iterations  = -1    , 
                 opt_const       = True  , ## optimize const 
                 strategy        = None  ,
                 silent          = False ,
                 print_level     = 0     , 
                 nLL             = None  , ## nLL
                 scale           = True  , ## scale weighted dataset ?
                 offset          = True  , ## offset the FCN values? 
                 args            =   ()  , **kwargs  ):
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

        scale_errdef = 1 
        if not nLL :
            nLL , sf = self.nll ( dataset , args = args , **kwargs )
            if dataset.isWeighted() and 1 != sf :
                if scale : scale_errdef = sf
                else     : self.warning("minuit: no FCN scaling for the weighted dataset is defined!")
            self.aux_keep.append ( nLL )
            
        m    = ROOT.RooMinimizer ( nLL )

        if  silent : m.setPrintLevel ( -1 ) 
        else       : m.setPrintLevel ( print_level )   
        
        if isinstance  ( opt_const  , bool ) : m.optimizeConst ( opt_const ) 
        if isinstance  ( max_calls      , integer_types ) and 1 < max_calls :
            m.setMaxFunctionCalls ( max_calls )
        if isinstance  ( max_iterations , integer_types ) and 1 < max_iterations :
            m.setMaxIterations    ( max_iterations  )
        if isinstance  ( strategy , integer_types       ) and 0 <= strategy <= 2 :
            m.setStrategy ( strategy )

        ## offset ? 
        m.setOffsetting ( offset )
        
        if 1 != scale_errdef :
            old_errdef = nLL.defaultErrorLevel ()
            new_errdef = old_errdef / scale_errdef
            logger.info ("minuit: Error Level is redefined from %.1f to %.4g" %  ( old_errdef ,
                                                                                   new_errdef ) )
            m.setErrorLevel ( new_errdef ) 

        return m  

    # =========================================================================
    ## perform sPlot-analysis 
    #  @code
    #  r,f = model.fitTo ( dataset )
    #  model.sPlot ( dataset ) 
    #  @endcode 
    def sPlot ( self , dataset , silent = False ) : 
        """ Make sPlot analysis
        >>> r,f = model.fitTo ( dataset )
        >>> model.sPlot ( dataset ) 
        """
        assert self.alist2,\
               "PDF(%s) has empty ``alist2''/(list of components)" + \
               "no sPlot is possible" % self.name 


        vars = set ( ( v.name for v in dataset.varset() ) ) 
        with roo_silent ( True ) :
            
            splot = ROOT.RooStats.SPlot ( rootID( "sPlot_" ) ,
                                          "sPlot"            ,
                                          dataset            ,
                                          self.pdf           ,
                                          self.alist2        )
        
            self.__splots += [ splot ]

        if not silent :
            vars = set ( ( v.name for v in dataset.varset() ) ) - vars
            if vars :  self.info ( 'sPlot: %d variables are added to dataset:\n%s' % (
                len ( vars ) ,
                dataset.table ( title     = 'Variables added by the sPlot' ,
                                prefix    = '# ' ,
                                variables = list ( vars ) ) ) )
            
        return splot 
    # =========================================================================
    ## generate toy-sample according to PDF
    #  @code
    #  model  = ....
    #  data   = model.generate ( 10000 ) ## generate dataset
    #  varset = ....
    #  data   = model.generate ( 100000 , varset , sample = False )
    #  data   = model.generate ( 100000 , varset , sample = True  )     
    #  @endcode
    def generate ( self             ,
                   nEvents          ,
                   varset   = None  ,
                   binning  = None  ,
                   sample   = True  , ##  sample number of events ?  
                   args     = ()    ) :
        """Generate toy-sample according to PDF
        >>> model  = ....
        >>> data   = model.generate ( 10000 ) ## generate dataset 
        
        >>> varset = ....
        >>> data   = model.generate ( 100000 , varset , sample = False )
        >>> data   = model.generate ( 100000 , varset , sample = True  )
        """
        
        ## sample number of events in dataset ?
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
        
        from ostap.stats.moments import rms      as sp_rms
        from ostap.stats.moments import variance as sp_variance 
        
        fun   = self.fun
        ftype = type ( fun ) 
        if   hasattr ( ftype , 'rms' ) and not ftype.rms is sp_rms :
            return fun.rms()        
        elif hasattr ( ftype , 'Rms' ) :
            return fun.Rms()        
        elif hasattr ( ftype , 'RMS' ) :
            return fun.RMS()        
        elif self.tricks and hasattr ( fun , 'function' ) :

            if   hasattr ( fun , 'setPars'    ) : fun.setPars()
            
            ff = fun.function()
            ftype = type ( ff ) 
            if   hasattr ( ftype , 'rms'        ) and not ftype.rms        is sp_rms      :
                return ff.rms()
            elif hasattr ( ftype , 'variance'   ) and not ftype.variance   is sp_variance :
                return ff.variance   ()**0.5  
            elif hasattr ( ftype , 'dispersion' ) and not ftype.dispersion is sp_variance :
                return ff.dispersion ()**0.5 

        return  self._get_stat_ ( sp_rms , **kwargs )

    # =========================================================================
    ## get the effective Skewness
    def skewness ( self , **kwargs ) :
        """Get the effective Skewness
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'SKEWNESS: %s ' % pdf.skewness()
        """
        ## use generic machinery 
        from ostap.stats.moments import skewness as sp_skewness
        return self._get_stat_ ( sp_skewness , **kwargs )

    # =========================================================================
    ## get the effective Kurtosis
    def kurtosis ( self , **kwargs ) :
        """Get the effective Kurtosis
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'KURTOSIS: %s ' % pdf.kurtosis()
        """
        ## use generic machinery 
        from ostap.stats.moments import kurtosis as sp_kurtosis
        return self._get_stat_ ( sp_kurtosis , **kwargs )

    # =========================================================================
    ## get the effective median
    def median ( self , **kwargs ) :
        """Get the effective median
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'MEDIAN: %s ' % pdf.median()
        """
        from ostap.stats.moments import median as _median
        return self._get_stat_ ( _median , **kwargs )

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
    ## get moment using <code>RooAbsPdf::moment</code> method
    #  @see RooAbsPdf::moment
    #  @code
    #  pdf = ...
    #  v4  = pdf.sroo_moment ( 5 , central = True )
    #  @endcode
    def roo_moment ( self , order , central ) :
        """Get moment using <code>RooAbsPdf::moment</code> method
        >>> pdf = ...
        >>> v5  = pdf.sroo_moment ( 5 , central = True )
        - see `ROOT.RooAbsPdf.moment`
        """
        assert isinstance ( order , integer_types ) and 0<= order , \
               'roo_moment: invalid moment order %s !' % order

        from ostap.logger.utils import rootWarning, roo_silent 
        with rootWarning() , roo_silent ( True ) : 
            mom    = self.pdf.moment ( self.xvar , order , central , False  )
            result = mom.getVal ()
            del mom
            return result 

    # =========================================================================
    ## get mean using RooAbdPdf method
    #  @see RooAbdPdf::moment 
    def roo_mean ( self ) :
        """get mean using RooAbdPdf method
         -see `ROOT.RooAbdPdf.moment`
         """
        return self.roo_moment ( 1 , central = False )

    # =========================================================================
    ## get variance using RooAbdPdf method
    #  @see RooAbdPdf::moment 
    def roo_variance ( self ) :
        """get variance  using RooAbdPdf method
         -see `ROOT.RooAbdPdf.moment`
         """
        return self.roo_moment ( 2 , central = True )

    # =========================================================================
    ## get dispersion  using RooAbdPdf method
    #  @see RooAbdPdf::moment 
    def roo_dispersion  ( self ) :
        """get dispersion using RooAbdPdf method
         -see `ROOT.RooAbdPdf.moment`
         """
        return self.roo_variance () 

    # =========================================================================
    ## get RMS using RooAbdPdf method
    #  @see RooAbdPdf::moment 
    def roo_rms ( self ) :
        """get RMS using RooAbdPdf method
         -see `ROOT.RooAbdPdf.moment`
         """
        v2 = self.roo_variance() 
        return math.sqrt ( v2 )

    # =========================================================================
    ## get skewness using RooAbdPdf method
    #  @see RooAbdPdf::moment 
    def roo_skewness ( self ) :
        """get skewness using RooAbdPdf method
        -see `ROOT.RooAbdPdf.moment`
        """
        m2 = self.roo_moment ( 2 , central = True )
        m3 = self.roo_moment ( 3 , central = True )        
        return m3/ ( m2 ** ( 3.0 / 2 ) )

    # =========================================================================
    ## get kurtosis  using RooAbdPdf method
    #  @see RooAbdPdf::moment 
    def roo_kurtosis ( self ) :
        """get kurtosis using RooAbdPdf method
        -see `ROOT.RooAbdPdf.moment`
        """
        m2 = self.roo_moment ( 2 , central = True )
        m4 = self.roo_moment ( 4 , central = True )    
        return m4 / ( m2 * m2 ) - 3.0 
    

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
        - Unlike  `PDF.roo_histo` method, PDF is integrated within the bin
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
    ## make 2D-cpontours
    # ==========================================================================
    def contours ( self              ,
                   var1              ,
                   var2              ,
                   dataset           ,
                   levels  = ( 1 , ) ,
                   npoints = 100     ,
                   **kwargs          ) :
        
        ## create the minuit 
        mn = self.minuit ( dataset , **kwargs )
        
        ## get the parametrs
        pars = self.params ( dataset ) 
        assert var1 in pars , "Variable %s is not a parameter"   % var1
        if not isinstance ( var1 , ROOT.RooAbsReal ) : var1 = pars [ var1 ]
        assert var2 in pars , "Variable %s is not a parameter"   % var2
        if not isinstance ( var2 , ROOT.RooAbsReal ) : var2 = pars [ var2 ]
        del pars 

        import ostap.fitting.roofitresult
        status = mn.migrad( tag = 'contours' )

        return mn.contour ( var1 , var2 , npoints , *levels ) 

        
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
    ## Create list of variables that can be used as ``fractions'' for N-component fit
    #  @code
    #  MV = ...
    #  fractions = MV.make_fractions ( 5 )
    #  fractions = MV.make_fractions ( 5 , name = 'F%d_A' , suffix = 'D' )
    #  fractions = MV.make_fractions ( 5 , name = 'F%d_B' , recursive = False )    
    #  fractions = MV.make_fractions ( 5 , name = 'F%d_C' , fractions = (0.4, 0.1 , 0.3 , 0.1 ) )   
    #  @endcode 
    def make_fractions ( self              ,
                         N                 ,
                         name      = 'f%d' , ## pattern to contruct the fraction name from index 
                         title     = ''    , 
                         recursive = True  ,
                         fractions = ()    ) :
        
        assert is_integer ( N ) and 2 <= N ,\
               "make_fractions: there must be at least two components!"

        my_fractions = make_iterable ( fractions , None )
            
        fracs = []
        value = 1.0
        prod  = 1.0


        for i , ff in zip ( range ( N - 1  ) , my_fractions ) :

            value = 1.0 / N
            
            if recursive :
                
                value /= prod
                prod  *= ( 1.0 - value ) 
                
            if   isinstance  ( ff , num_types       ) and not 0 <= ff <= 1 :
                self.error   ("make_fractions: fraction %s is outside [0,1] interval, ignore it!" %  ff )
                ff = value
                
            elif isinstance  ( ff , ROOT.RooAbsReal ) and not 0 <= float ( ff ) <= 1 :
                self.warning ("make_fractions: fraction %s is outside [0,1] interval" %  ff )
                
            if ff is None : ff = value

            fname = name if ( 2 == N and  '%s' not in name and '%d' not in name ) else name % i  
            
            tit   =  title if title else 'Fraction #%d: %s %s' % ( i , fname , self.name )                 
            fvar  = self.make_var ( ff   , fname  , tit , None , value , 0 , 1 )
            
            fracs.append ( fvar  )
            
        return tuple ( fracs )

    # =============================================================================
    ## create a list of variables that can be used as ``yields'' for N-component fit
    #  @code
    #  M = ...
    #  nums = M.make_yields ( 5 , name = 'S_%d_A' , minmax = (0, 1000) )
    #  nums = M.make_yields ( 5 , name = 'S_%d_B' , minmax = (0, 1000) ,  yields =  ( 1 , 100 , 50 , 10 , 10 ) )
    #  @endcode     
    def make_yields ( self                   ,
                      N                      ,
                      name  = 'N%d'          , ## pattern to construct yield names from index
                      title = ''             ,
                      minmax = ( 0 , 1.e+6 ) , 
                      yields = ()            ) :
        """create a list of variables that can be used as ``yields'' for N-component fit
        >>> M = ...
        >>> nums = M.make_yields ( 5 , name = 'S_%d_A' , minmax = (0, 1000) )
        >>> nums = M.make_yields ( 5 , name = 'S_%d_B' , minmax = (0, 1000) ,  yields =  ( 1 , 100 , 50 , 10 , 10 ) )
        """
        my_yields = make_iterable ( yields , None )

        nums = []        
        for i , n in zip ( range ( N ) , my_yields ) : 

            fname = ( name  % i )  if 1 != N else name           ## ATTENTION!  
            tit   = title if title else 'Yield #%d: %s %s' % ( i , fname , self.name )                 
            if not n is None : nvar  = self.make_var ( n    , fname  , tit , None , n , *minmax )
            else             : nvar  = self.make_var ( n    , fname  , tit , None ,     *minmax )
            
            nums.append ( nvar )

        return tuple ( nums ) 
            
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
        pdf    = ROOT.RooAddPdf ( self.roo_name ( name ) , title , pdfs , fracs , recursive )
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

        
# =============================================================================
##  helper utilities to imlement resolution models.
# =============================================================================
class _CHECKMEAN(object) :
    check = True
def checkMean() :
    return True if  _CHECKMEAN.check else False
# =============================================================================
## @class CheckMean 
#  Helper contex manager to enable/disable check for the mean/location-values
class CheckMean(object) :
    """Helper contex manager to enable/disable check for the mean/location-values
    """
    def __init__  ( self , check ) :
        self.__check = True if check else False 
    def __enter__ ( self ) :
        self.__old       = _CHECKMEAN.check 
        _CHECKMEAN.check =  self.__check
    def __exit__  ( self , *_ ) :
        _CHECKMEAN.check =  self.__old
    @property
    def check ( self ) :
        """``check''  : check the mean/location?"""
        return self.__check
    
# =============================================================================
## helper base class for implementation  of various helper pdfs
#  - it defines alias <code>mass</code> for <code>xvar</code>
#  - it defiens a variable <code>mean</code> alias <code>location</code>
#  - optionally it checks that this variable is withing the specified range  
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class MASSMEAN(PDF) :
    """Helper base class for implementation of various pdfs
    It is useful for ``peak-like'' distributions, where one can talk about
    - ``mean/location''
    - it defines alias `mass` for `xvar`
    - it defiens a variable `mean` (alias `location`)
    - optionally it checks that this variable is withing the specified range  
    """
    def __init__ ( self              ,
                   name              ,
                   xvar              ,
                   mean      = None  ,
                   mean_name  = ''   , 
                   mean_title = ''   ) : 
        
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
            raise AttributeError("MASSMEAN: Unknown type of ``xvar'' parameter %s/%s" % ( type ( xvar ) , xvar ) )

        ## intialize the base 
        PDF.__init__ ( self , name , xvar = xvar )

        ## check mean/location values ? 
        self.__check_mean = self.xminmax () and checkMean () 
        
        self.__limits_mean  = ()
        if  self.check_mean and self.xminmax () and not isinstance ( mean , ROOT.RooAbsReal ) :      
            mn , mx = self.xminmax()
            dm      =  mx - mn
            self.__limits_mean  = mn - 0.2 * dm , mx + 0.2 * dm

        ## mean-value
        m_name  = mean_name  if mean_name  else "mean_%s"  % name
        m_title = mean_title if mean_title else "mean(%s)" % name
        self.__mean = self.make_var ( mean , m_name , m_title , mean , *self.limits_mean )
        
        ##
        if self.limits_mean :  
            mn , mx = self.limits_mean  
            dm      =  mx - mn
            if   self.mean.isConstant() :
                if not mn <= self.mean.getVal() <= mx : 
                    self.error ( 'MASSMEAN(%s): Fixed mass %s is not in mass-range (%s,%s)' % ( name , self.mean.getVal() , mn , mx  ) )
            elif self.mean.minmax() :
                mmn , mmx = self.mean.minmax()
                self.mean.setMin ( max ( mmn , mn ) )
                self.mean.setMax ( min ( mmx , mx ) )
                self.debug ( 'mean range is adjusted  to be %s' % list ( self.mean.minmax() ) )

        ## save the configuration
        self.config = {
            'name'        : self.name  ,
            'xvar'        : self.xvar  ,
            'mean'        : self.mean  ,
            'mean_name'   : mean_name  ,
            'mean_title'  : mean_title ,
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
        mn , mx = self.mean.minmax()
        if not mn <= value <= mx :
            self.warning( "``%s'': %s is outside the interval (%s,%s)/1" % ( self.mean.name , value , mn , mx ) )
        if self.check_mean and self.limits_mean  :  
            mn , mx = self.limits_mean 
            if not mn <= value <= mx :
                self.error ("``%s'': %s is outside the interval (%s,%s)/2"  % ( self.mean.name , value , mn , mx ) )                
        self.mean.setVal ( value )
        
    @property
    def location ( self ):
        """``location/mean''-variable (the same as ``mean'')"""
        return self.mean
    @location.setter
    def location ( self , value ) :
        self.mean =  value

    @property
    def check_mean ( self ) :
        """``check_mean'' : Is mean-value to be checked?"""
        return self.__check_mean
    @check_mean.setter
    def check_mean ( self, value ) :
        self.__check_mean = True if  value else False
        
    @property
    def limits_mean ( self ) :
        """``limits_mean'' : reasonable limits for mean/location"""
        return self.__limits_mean
    
# =============================================================================
## @class MASS
#  helper base class for implementation  of various helper pdfs 
#  - mean/location
#  - sigma/width/scale
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class MASS(MASSMEAN) :
    """Helper base class for implementation of various pdfs
    It is useful for ``peak-like'' distributions, where one can talk about
    - ``mean/location''
    - ``sigma/width/scale'' 
    """
    def __init__ ( self               ,
                   name               ,
                   xvar               ,
                   mean        = None ,
                   sigma       = None , 
                   mean_name   = ''   , 
                   mean_title  = ''   ,
                   sigma_name  = ''   , 
                   sigma_title = ''   ) : 
            
        ## base class 
        MASSMEAN.__init__ ( self                    ,
                            name       = name       ,
                            xvar       = xvar       , 
                            mean       = mean       ,
                            mean_name  = mean_name  ,
                            mean_title = mean_title )
        
        self.__limits_sigma = ()        
        if  self.xminmax() and not isinstance ( sigma , ROOT.RooAbsReal ) :            
            mn , mx   = self.xminmax()
            dm        =  mx - mn
            sigma_max =  2 * dm / math.sqrt(12)  
            self.__limits_sigma = 1.e-4 * sigma_max , sigma_max 

        ## sigma
        s_name  = sigma_name  if sigma_name  else "sigma_%s"   % name
        s_title = sigma_title if sigma_title else "#sigma(%s)" % name
        #
        self.__check_sigma = True 
        self.__sigma = self.make_var ( sigma  , s_name , s_title , sigma , *self.limits_sigma )
        
        ## save the configuration
        self.config = {
            'name'        : self.name   ,
            'xvar'        : self.xvar   ,
            'mean'        : self.mean   ,
            'sigma'       : self.sigma  ,
            'mean_name'   : mean_name   ,
            'mean_title'  : mean_title  ,
            'sigma_name'  : sigma_name  ,
            'sigma_title' : sigma_title ,
            }
            
    @property
    def sigma ( self ):
        """``sigma/width/scale/spread''-variable"""
        return self.__sigma
    @sigma.setter
    def sigma ( self , value ) :
        value =   float ( value )
        mn , mx = self.sigma.minmax()
        if not mn <= value <= mx :
            self.warning ("``%s'': %s is outside the interval (%s,%s)/1" % ( self.sigma.name , value , mn , mx ) )
        if self.limits_sigma and self.check_sigma  : 
            mn , mx = self.limits_sigma 
            if not mn <= value <= mx :
                self.error ("``%s'': %s is outside the interval (%s,%s)/2" % ( self.sigma.name , value , mn , mx ) )
        self.sigma.setVal ( value )

    @property
    def check_sigma ( self ) :
        """``check_mean'' : Is mean-value to be checked?"""
        return self.__check_sigma 
    @check_sigma.setter
    def check_sigma ( self, value ) :
        self.__check_sigma = True if  value else False 
    
    @property
    def limits_sigma ( self ) :
        """``limits_sigma'' : reasonable limits for sigma/width"""
        return self.__limits_sigma

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
    - It contains ``fudge-factor'' for the resolution parameter ``sigma''
    - It simplify creation of the soft/gaussian constraint for the ``fudge-factor''    
    """
    ## constructor
    #  @param name   the name of PDF
    #  @param xvar   the variable/observable
    #  @param sigma  sigma/resoltuion parameter 
    #  @param mean   "mean"-variable
    #  @param fudge  "fudge-factor" to be aplied to sigma
    def __init__ ( self               ,
                   name               ,
                   xvar        = None ,
                   sigma       = None , 
                   mean        = None ,
                   fudge       = 1.0  ,
                   mean_name   = ''   ,
                   mean_title  = ''   ,
                   sigma_name  = ''   ,
                   sigma_title = ''   ) :
        
        ## mean-value
        if mean is None :
            mean = ROOT.RooRealConstant.value ( 0 ) 
            
        with CheckMean ( False ) :
            super(RESOLUTION,self).__init__ ( name        = name        ,
                                              xvar        = xvar        ,
                                              sigma       = sigma       ,
                                              mean        = mean        ,
                                              mean_name   = mean_name   ,
                                              mean_title  = mean_title  ,
                                              sigma_name  = sigma_name  ,
                                              sigma_title = sigma_title )
            
        self.__fudge            = fudge        
        self.__fudge_constraint = None
        self.__sigma_corr       = None
        
        if isinstance ( fudge , VE ) :
            
            assert 0 < fudge.value() and 0 < fudge.cov2(),\
                   "Invalid value for ``fudge-factor'': %s" % s 
            
            value  = fudge.value()
            error  = fudge.error()
            vmin   = max ( 1.e-3 , value - 10 * error )
            vmax   =               value + 10 * error
            
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
            
        elif isinstance ( fudge , ROOT.RooAbsReal ) :
            
            ## make fudge-factor to be a variable 
            self.__fudge = self.make_var ( fudge  ,
                                           'fudge_factor_%s'  % self.name ,
                                           'fudge_factor(%s)' % self.name , None )
            
        elif isinstance ( fudge , num_types ) and 1 == fudge :

            ## fudge is trivial 
            self.__fudge = ROOT.RooFit.RooConst ( fudge )

            ## corrected sigma is trivial 
            self.__sigma_corr   =    self.sigma
            
        elif isinstance ( fudge , num_types ) :
            
            ## fudge is trivial 
            self.__fudge = ROOT.RooFit.RooConst ( fudge )

        else :
                
            ## make fudge-factor to be a variable 
            self.__fudge = self.make_var ( fudge  ,
                                           'fudge_factor_%s'  % self.name ,
                                           'fudge_factor(%s)' % self.name ,
                                           None                           , fudge , 0.01 , 10 ) 
            
        ## create corrected sigma 
        if self.__sigma_corr is None :            
            ## corrected sigma 
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
    def __init__ ( self , xvar , name = '' , title = '' ) :
        
        name = name if name else self.generate_name ( prefix = 'flat1D_')
        PDF.__init__ ( self  , name , xvar ) 
        
        if not title : title = 'flat1(%s)' % name
        
        self.pdf = Ostap.Models.Uniform ( self.roo_name ( 'flat_' ) , title , self.xvar )
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

        name = name if name else self.generate_name ( prefix = prefix + '%s_' % pdf.GetName() , suffix = suffix )
        
        ## initialize the base 
        PDF . __init__ ( self , name , xvar , special = special )
        ##
        if not self.special :
            assert isinstance ( pdf  , ROOT.RooAbsPdf ) , "``pdf'' must be ROOT.RooAbsPdf"

        ## Does PDF depends on XVAR ?
        if not pdf.depends_on ( xvar ) :
            self.warning ( "PDF/%s does not depend on %s!" % ( pdf.name , xvar.name ) ) 

        ## PDF itself 
        self.pdf  = pdf

        if not self.xvar in self.params () : 
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
## @class Combine1D
#  Non-extended sum of several PDFs
#  It is just a small wrapper for <code>ROOT.RooAddPdf</code>
#  @see RooAddPdf 
class Combine1D (PDF) :
    """Non-extended sum of several PDFs:
    
    It is just a small wrapper for <code>ROOT.RooAddPdf</code>
    - see RooAddPdf 
    
    >>> sum  = Combine1D ( [ pdf1 , pdf2 , pdf3 ]  ) 
    
    """
    def __init__ ( self             ,
                   pdfs             , ## input list of PDFs  
                   xvar      = None , 
                   name      = ''   ,
                   recursive = True ,
                   prefix    = 'f'  , ## prefix for fraction names 
                   suffix    = ''   , ## suffix for fraction names 
                   fractions = None ) :

        assert 2 <= len ( pdfs ) , 'Combine1D: at least two PDFs are needed!'

        pdf_list = []        
        for i , p in enumerate ( pdfs ) :

            if isinstance ( p , PDF ) :
                
                assert ( not xvar ) or xvar is p.xvar, "Invalid xvar/pdf%d.xvar: %s/%s" % ( i , xvar , p.xvar ) 
                xvar = p.xvar
                pdf_list.append ( p )
                
            elif isinstance ( p  , ROOT.RooAbsPdf ) and xvar and isinstance ( xvar , ROOT.RooAbsReal ) :
                
                pdf_list.append ( Generic1D_pdf ( p , xvar ) )
                
            else :
                raise TypeError ( "Invalid type: pdf%d xvar %s/%s xvar=%s" % ( i , p , type(p) , xvar ) )


        ## check the name 
        name = name if name else self.generate_name ( prefix = 'sum1' )
        
        ## ininialize the base class
        PDF.__init__ ( self , name , xvar ) 

        for i , p in enumerate ( pdf_list )  :
            if p.pdf.canBeExtended() : self.warning ("``pdf%f'' can be extended!" % i ) 
                
        while prefix.endswith  ('_') : prefix = prefix[:-1]
        while suffix.startswith('_') : suffix = suffix[1:]
        
        self.__prefix    = prefix if prefix else 'f'
        self.__suffix    = suffix
        self.__recursive = True if recursive else False 

        N = len ( pdf_list )

        if 2 < len ( pdf_list ) : 
            if self.prefix and self.suffix : fr_name = '%s_%%d_%s' % ( self.prefix , self.suffix )
            else                           : fr_name = '%s_%%d'    %   self.prefix
        else :
            if self.prefix and self.suffix : fr_name = '%s_%s'     % ( self.prefix , self.suffix )
            else                           : fr_name =                 self.prefix
    
        ## make list of fractions
        fraction_list = self.make_fractions  ( N                          ,
                                               name      = fr_name        , 
                                               recursive = self.recursive ,
                                               fractions = fractions      )

        ## keep them
        self.__pdfs      = tuple ( pdf_list )
        self.__fractions = tuple ( fraction_list ) 
        
        for p in self.__pdfs      : self.alist1.add ( p.pdf )
        for f in self.__fractions : self.alist2.add ( f     )
        
        ## finally build PDF 
        self.pdf = ROOT.RooAddPdf ( self.roo_name ( 'combine1' ) ,
                                    ' + '.join ( '(%s)' % p.name for p in self.pdfs  ) ,
                                    self.alist1    ,
                                    self.alist2    ,
                                    self.recursive )
        
        self.config = {
            'pdfs'      : self.pdfs      ,
            'xvar'      : self.xvar      ,
            'name'      : self.name      , 
            'prefix'    : self.prefix    ,
            'suffix'    : self.suffix    ,
            'fractions' : self.fractions ,
            'recursive' : self.recursive        
            }
        
    @property
    def prefix ( self ) :
        """``prefix'' : prefix for fraction names"""
        return self.__prefix

    @property
    def suffix ( self ) :
        """``suffix'' : suffix for fraction names"""
        return self.__suffix

    @property
    def recursive ( self ) :
        """``recursive'' : recursive fractions?"""
        return self.__recursive
    
    @property
    def pdfs ( self ) :
        """``pdfs'' : get list/tuple of involved PDFs (same as ``components'')"""
        return self.__pdfs
    @property
    def components ( self ) :
        """``components'' : get list/tuple of involved PDFs (same as ``pdfs'')"""
        return self.pdfs
        
    @property
    def fractions ( self ) :
        """``fractions'' : get involved fractions (same as ``F'')"""
        return self.component_getter ( self.__fractions )
    @fractions.setter
    def fractions ( self , values ) :
        self.component_setter ( self.__fractions , values )

    @property
    def F         ( self ) :
        """``F'' : get involved fractions (same as ``fractions'')"""
        return self.component_getter ( self.__fractions )
    @F.setter
    def F         ( self , values ) :
        self.fractions = values 

        
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
class Sum1D(Combine1D) :
    """Non-extended sum of two PDFs:    
    It is just a small wrapper for `ROOT.RooAddPdf`    
    >>> pdf1 = ...
    >>> pdf2 = ...
    >>> sum  = Sum1D ( pdf1 , pdf2 ) 
    - see `ROOT.RooAddPdf` 
    """
    def __init__ ( self             ,
                   pdf1             ,
                   pdf2             ,  
                   xvar      = None ,                   
                   name      = ''   , 
                   prefix    = 'f'  ,
                   suffix    = ''   ,
                   fraction  = None ,
                   others    = []   ,
                   recursive = True ) :                    
        
        ## check the name 
        name = name if name else self.generate_name ( prefix = 'sum1' )
        
        ## initialize the base class 
        Combine1D.__init__ ( self                  ,
                             name      = name      , 
                             pdfs      = [ pdf1 , pdf2 ] + others ,
                             xvar      = xvar      ,
                             recursive = recursive ,         
                             prefix    = prefix    ,                             
                             suffix    = suffix    ,
                             fractions = fraction  )

        self.config = {
            'pdf1'      : self.pdf1      ,
            'pdf2'      : self.pdf2      ,
            'xvar'      : self.xvar      ,
            'name'      : self.name      ,
            'prefix'    : self.prefix    ,
            'suffix'    : self.suffix    , 
            'fraction'  : self.fraction  ,
            'others'    : self.others    ,            
            'recursive' : self.recursive ,
            }
        
    @property
    def pdf1 ( self ) :
        """``pdf1'' : the first PDF"""
        return self.pdfs[0]
    
    @property
    def pdf2 ( self ) :
        """``pdf2'' : the second PDF"""
        return self.pdfs[1]
    
    @property 
    def others ( self ) :
        """``others'' : other PDFs (if any)"""
        return self.pdfs[2:]
    
    @property
    def fraction ( self ) :
        """``fraction'' : the fraction of the first PDF in the sum (same as ``fractions'')"""
        return self.fractions 
    @fraction.setter
    def fraction ( self , value ) :
        self.fractions = value 

# =============================================================================
    
# =============================================================================
## Helper function to create the PDF/PDF2/PDF3
#  @param pdf   input pdf of funcntion   <code>RooAbsReal</code> or <code>RooAbsPdf</code>
#  pa
def make_pdf ( pdf , args , name = '' ) :
    """Helper function to create the PDF/PDF2/PDF3
    """
    
    assert pdf and isinstance ( pdf , ROOT.RooAbsReal ), \
           'make_pdf: Invalid type %s' % type ( pdf )
    
    name = name if name else "PDF_from_%s" % pdf.name
    
    if not isinstance ( pdf , ROOT.RooAbsPdf ) :
        if (6,20) <= root_info : 
            pdf = ROOT.RooWrapperPdf  ( name , 'PDF from %s' % pdf.name , pdf )
        else :
            raise TypeError("make_pdf: RooWrapperPdf is not available for ROOT %s" % root_version_int )
        
    num = len ( args )
    if   1 == num :
        return Generic1D_pdf ( pdf , name = name , *args )
    elif 2 == num :
        from ostap.fitting.fit2d import Generic2D_pdf 
        return Generic2D_pdf ( fun , name = name , *args )
    elif 3 == num :
        from ostap.fitting.fit3d import Generic3D_pdf 
        return Generic3D_pdf ( fun , name = name , *args )
    
    raise TypeError ( "Invalid length of arguments %s " % num ) 

# =============================================================================
## Generic 1D-shape from C++ callable
#  @see Ostap::Models:Shape1D
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2020-07-20
class Shape1D_pdf(PDF) :
    """ Generic 1D-shape from C++ callable
    - see Ostap::Models:Shape1D
    """
    
    def __init__ ( self , name , shape , xvar ) :


        if isinstance ( shape , ROOT.TH1 ) and not isinstance ( shape , ROOT.TH2 ) and not xvar :
            xvar = shape.xminmax() 

        if isinstance ( shape , ROOT.TH1 ) and not isinstance ( shape , ROOT.TH2 ) :
            self.histo = shape
            shape      = Ostap.Math.Histo1D ( shape )

        ##  iniialize the base 
        PDF.__init__ ( self , name , xvar ) 
        
        self.__shape = shape

        if isinstance ( self.shape , Ostap.Math.Histo1D ) :
        
            ## create the actual pdf
            self.pdf = Ostap.Models.Histo1D ( self.roo_name ( 'histo1_' ) , 
                                              "Histo-1D %s" % self.name   ,
                                              self.xvar                   ,
                                              self.shape                  )
            
        else :
            
            ## create the actual pdf
            self.pdf = Ostap.Models.Shape1D.create  (
                self.roo_name ( 'shape1_' ) , 
                "Shape-1D %s" % self.name ,
                self.xvar                 ,
                self.shape                ) 

        ## save the configuration
        self.config = {
            'name'    : self.name    , 
            'shape'   : self.shape   , 
            'xvar'    : self.xvar    , 
            }
        
    @property
    def shape  ( self ) :
        """``shape'': the actual C++ callable shape"""
        return self.__shape 
            
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
                   order   = 0     , ## interpolation order 
                   silent  = False ) :
        
        H1D_dset.__init__ ( self , histo = histo , xaxis = xvar , density = density , silent = silent )
        PDF     .__init__ ( self , name  , self.xaxis ) 

        assert isinstance ( order, integer_types ) and 0 <= order ,\
               'Invalid interpolation order: %s/%s' % ( order , type ( order ) )

        with roo_silent ( silent ) : 
            #
            ## finally create PDF :
            self.__vset = ROOT.RooArgSet  ( self.xvar )        
            self.pdf    = ROOT.RooHistPdf (
                self.roo_name ( 'histo1_' ) ,
                'Histo-1D PDF: %s/%s' % ( histo.GetName() , histo.GetTitle() ) , 
                self.__vset , 
                self.dset   ,
                order       )
            
        ## and declare it be be a "signal"
        self.signals.add ( self.pdf ) 
        
        ## save the configuration
        self.config = {
            'name'    : self.name    , 
            'histo'   : self.histo   , 
            'xvar'    : self.xvar    , 
            'density' : self.density , 
            'silent'  : self.silent  ,
            'order'   : self.order   ,
            }
        
    @property
    def order  ( self ) :
        """``order'': interpolation order"""
        return self.pdf.getInterpolationOrder () 
    @order.setter
    def order  ( self , value ) :
        assert isinstance ( value , integer_types ) and 0 <= value,\
               'Invalid interpolation order %s/%s' % ( value , type ( value ) )
        self.pdf.setInterpolationOrder ( value )
        
# =============================================================================
## @class Fit1D
#  The actual model for 1D-mass fits
#  @param signal                PDF for 'signal'     component                 (Ostap/PDF or RooAbsPdf)
#  @param background            PDF for 'background' component                 (Ostap/PDF or RooAbsPdf)
#  @param othersignals          list of PDFs for other 'signal' components     (Ostap/PDF or RooAbsPdf)
#  @param otherbackgrouds       list of PDFs for other 'background' components (Ostap/PDF or RooAbsPdf)
#  @param others                list of 'other' components                     (Ostap/PDF or RooAbsPdf)
#  @param signals               list 'signal'     components
#  @param backgrounds           list of 'background' component               
#  @param suffix                ... add this  suffix for the PDF name
#  @param name                  the name of compound PDF 
#  @param xvar                  the fitting variable, must be specified if components are given as RooAbsPdf
#  @param extended              build 'extended' PDF
#  @param combine_signals       combine all signal components into single SIGNAL?
#  @param combine_backgrounds   combine all background components into single BACKGROUND?
#  @param combine_others        combine all other components into single COMPONENT?
#  @param recursive             use recursive fractions for compound PDF
#  @param recursive_signals     use recursive fractions for compound signal
#  @param recursive_backgriunds use recursive fractions for compound background
#  @param recursive_others      use recursive fractions for compound others
#  @param S                     yields of signal components 
#  @param B                     yields of background components 
#  @param C                     yields of others components
#  @param F                     component fractions for non-extended fit
#  @param fS                    fractions for compound signal
#  @param fB                    fractions for compound background
#  @param fC                    fractions for compound others 
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
                   signal                = None    ,    ## the main signal 
                   background            = None    ,    ## the main background 
                   othersignals          = ()      ,    ## additional signal         components
                   otherbackgrounds      = ()      ,    ## additional background     components
                   others                = ()      ,    ## additional non-classified components
                   signals               = ()      ,    ## alternative : all signals 
                   backgrounds           = ()      ,    ## alternative : all backgrounds  
                   suffix                = ''      ,    ## the suffix 
                   name                  = ''      ,    ## the name
                   xvar                  = None    ,    ## x-variable 
                   extended              = True    ,    ## extended fits ?
                   combine_signals       = False   ,    ## combine signal     components into single "SIGNAL"     ? 
                   combine_backgrounds   = False   ,    ## combine background components into single "BACKGROUND" ?            
                   combine_others        = False   ,    ## combine other      components into single "COMPONENT"  ?             
                   recursive             = True    ,    ## recursive fractions for NON-extended models?
                   recursive_signals     = True    ,    ## recursive fractions for combined signal?
                   recursive_backgrounds = True    ,    ## recursive fractions for combined background?
                   recursive_others      = True    ,    ## recursive fractions for combined other components?
                   S                     = ()      ,    ## yields for ``signals''
                   B                     = ()      ,    ## yields for ``background''
                   C                     = ()      ,    ## yields for ``components''
                   F                     = ()      ,    ## fractions for noin-extended fit 
                   fS                    = ()      ,    ## fraction for combined signal
                   fB                    = ()      ,    ## fraction for combined background
                   fC                    = ()      ,    ## fraction for combined components
                   **kwargs                        ) :  ## other arguments (e.g. drawing)
        

        self.__suffix                = suffix
        self.__extended              = True if extended              else False
        
        self.__combine_signals       = True if combine_signals       else False
        self.__combine_backgrounds   = True if combine_backgrounds   else False
        self.__combine_others        = True if combine_others        else False

        self.__recursive             = True if recursive             else False
        self.__recursive_signals     = True if recursive_signals     else False
        self.__recursive_backgrounds = True if recursive_backgrounds else False
        self.__recursive_others      = True if recursive_others      else False

        self.__signal_components     = ()
        self.__background_components = () 
        self.__other_components      = ()
        
        # =====================================================================
        ## Signals
        # =====================================================================
        
        assert       signals   or       signal         , "Fit1D:``signals'' or ``signal'' must be specified!"        
        assert ( not signals ) or ( not signal       ) , "Fit1D:``signals'' and ``signal'' are mutually exclusive!"        
        assert ( not signals ) or ( not othersignals ) , "Fit1D:``signals'' and ``othersignals'' are mutually exclusive!"  
        
        sig_lst = list ( othersignals ) + list ( signals ) 
        if signal : sig_lst.insert ( 0 , signal )

        pdfs = [] 
        for i , signal  in enumerate ( sig_lst ) :

            if isinstance ( signal , PDF ) :
                
                assert ( not xvar ) or xvar is signal.xvar, \
                       "Invalid ``signal'' component #%d/xvar/xvar %s/%s" % ( i , xvar , signal.xvar )
                
                pdfs.append ( signal )
                
                xvar = signal.xvar
                
            elif isinstance ( signal , ROOT.RooAbsPdf ) and xvar and isinstance ( xvar , ROOT.RooAbsReal ) :

                pdfs.append ( Generic1D_pdf ( signal , xvar = xvar , prefix = "Sig%d_" % i , suffix = self.suffix ) )
                
            else :
                
                raise AttributeError ( "Fit1D:Invalid ``signal'' #%d/xvar: %s/%s %s"  % ( i , sig , type( signal ) , xvar ) )

        ## sinal components 
        self.__signal_components  = tuple ( pdfs )
        
        # =====================================================================
        ## initialize the base class
        # =====================================================================
        name = name if name else self.generate_name ( prefix = 'Fit%s' % self.signal_components[0].name , suffix = self.suffix ) 
        PDF.__init__ ( self , name , xvar , **kwargs ) 
                            
        # =====================================================================
        ## Backgrounds
        # =====================================================================

        assert ( not backgrounds ) or ( not otherbackgrounds ) ,\
               "Fit1D:``backgrounds'' and ``otherbackgrounds'' are mutually exclusive!"        
        assert ( not backgrounds ) or ( not background ) ,\
               "Fit1D:``backgrounds'' and ``backgrounds'' are mutually exclusive!"        
        
        if not backgrounds :
            ## create background
            background = self.make_bkg ( background , name = 'Background' , xvar = self.xvar ) 
            
        bkg_lst  = list ( otherbackgrounds ) + list ( backgrounds )
        if background : bkg_lst.insert ( 0 , background ) 

        pdfs = []
        for i , bkg in enumerate ( bkg_lst ) :

            if isinstance ( bkg , PDF ) :
                
                assert self.xvar is bkg.xvar, \
                       "Invalid ``background'' component #%d/xvar/xvar %s/%s" % ( i , self.xvar , bkg.xvar )
                
                pdfs.append ( bkg )
                
            elif isinstance ( bkg  , ROOT.RooAbsPdf ) :
                
                pdfs.append ( Generic1D_pdf ( bkg , xvar = xvar , prefix = "Bkg%d_" % i , suffix = self.suffix ) )
                
            else :
                
                raise AttributeError ( "Fit1D:Invalid ``background''#%d: %s/%s"  % ( i , bkg , type( bkg )  ) )

        ## background components 
        self.__background_components = tuple ( pdfs ) 
        
        # =====================================================================
        ## Other fit components
        # =====================================================================
        
        pdfs = [] 
        for i, cmp in enumerate ( others ) :
            
            if isinstance ( cmp , PDF ) :
                
                assert self.xvar is cmp.xvar, \
                       "Invalid ``other'' component #%d/xvar/xvar %s/%s" % ( i , self.xvar , cmp.xvar )
                
                pdfs.append ( cmp )
                
            elif isinstance ( cmp  , ROOT.RooAbsPdf ) :
                
                cmp_pdfs.append ( Generic1D_pdf ( cmp , xvar = xvar , prefix = "Cmp%d_" % i , suffix = self.suffix ) )
                
            else :
                
                raise AttributeError ( "Fit1D:Invalid ``other''#%d: %s/%s"  % ( i , cmp , type( cmp ) ) )

        ## all other componnets 
        self.__other_components = tuple ( pdfs ) 
        
        # =====================================================================
        ## Merge them if requested 
        # =====================================================================
        self.__combined_signal     = None
        self.__combined_background = None
        self.__combined_others     = None

        sigs = list ( self.signal_components     )
        bkgs = list ( self.background_components ) 
        cmps = list ( self.other_components      ) 
        
        if combine_signals     and 1 < len ( self.signal_components     ) :
            
            combined = Combine1D ( self.signal_components    ,
                                   prefix    = 'fS'          ,
                                   suffix    = suffix        ,
                                   recursive = recursive     ,                                              
                                   fractions = fS            ) ## read  fS from arguments
            
            self.__combined_signal = combined 
            sigs = [ combined ] 
            
        if combine_backgrounds and 1 < len ( self.background_components ) :
                
            combined = Combine1D ( self.background_components ,
                                   prefix    = 'fB'           ,
                                   suffix    = suffix         ,
                                   recursive = recursive      ,                                             
                                   fractions = fB             ) ## read fB from arguments 
            
            self.__combined_background = combined
            bkgs = [ combined ]  
            
        if combine_others      and 1 < len ( self.other_components           ) :            

            combined = Combine1D ( self.other_components     ,
                                   prefix    = 'fC'          ,
                                   suffix    = suffix        ,
                                   recursive = recursive     ,                                              
                                   fractions = fC            ) ## read fC from arguments 
            
            self.__combined_others = combined
            cmps = [ combined ] 

        self.__fit_signals     = tuple ( sigs )
        self.__fit_backgrounds = tuple ( bkgs )
        self.__fit_others      = tuple ( cmps )
        
        ## final list of fit components 
        self.__fit_components  = self.fit_signals + self.fit_backgrounds + self.fit_others 
        for p in self.fit_components  : self.alist1.add ( p.pdf )  

        ## Yields/fractions 
        
        self.__S = () 
        self.__B = () 
        self.__C = ()
        self.__F = ()
        
        if self.extended :
            
            ns        = len  ( self.fit_signals ) 
            fname     = make_name ( 'S' , '%d' if 1 != ns else '' , suffix )
            title     = "Yield(s) for ``signal'' component(s)/%s" % self.name
            title     = title if 1 != ns else title.replace ( '(s)', '' )

            self.__S  = self.make_yields ( ns , fname , title , yields = S ) ## read S from arguments
            
            nb        = len ( self.fit_backgrounds ) 
            fname     = make_name ( 'B' , '%d' if 1 != nb else '' , suffix )
            title     = "Yield(s) for ``background'' component(s)/%s" % self.name
            title     = title if 1 != nb else title.replace ( '(s)', '' )             
            self.__B  = self.make_yields ( nb , fname , title , yields = B ) ## read S from arguments

            nc        = len ( self.fit_others  ) 
            fname     = make_name ( 'C' , '%d' if 1 != nc else '' , suffix )
            title     = "Yield(s) for ``other'' component(s)/%s" % self.name
            title     = title if 1 != nc else title.replace ( '(s)', '' )             
            self.__C  =  self.make_yields ( nc , fname , title , yields = C ) ## read S from arguments

            for y in self.yields : self.alist2.add ( y )   

            assert len ( self.alist2 ) == len ( self.alist1 ) ,\
                   'Fit1D: inconsistent parameters for ROOT.RooAddPdf' 

        else :
            
            assert 1 < len ( self.fit_components ) ,\
                   'Fot1D: At least two components are required to build proper non-extended PDF!'

            nt       = len ( self.fit_components )
            nf       = nt - 1 
            fname    = make_name ( 'F' , '%d' if 1 != nf else '' , suffix )
            title    = "Fraction(s) for various component(s)/%s" % self.name
            title    = title if 1 != nf else title.replace ( '(s)', '' )                         
            self.__F = self.make_fractions ( nt , fname ,  title , fractions = F ) ## read F from arguments
                        
            for f in self.__F : self.alist2.add ( f )

            assert len ( self.alist2 ) + 1 == len ( self.alist1 ) ,\
                   'Fit1D: inconsistent parameters for ROOT.RooAddPdf' 

        
            
        ## now we finally can create PDF
            
        pdf_name  = self.roo_name ( 'fit1d_' ) 
        pdf_title = "Fit1D %s" % self.name
        pdf_args  = pdf_name , pdf_title , self.alist1 , self.alist2

        if not self.extended :
            pdf_args = pdf_args + ( True if recursive else False , ) ## RECURSIVE ?
            
        self.pdf = ROOT.RooAddPdf ( *pdf_args )

        ## sanity checks

        ## drawing stuff
        if self.combined_background         : self.combined_backgrounds.add ( self.combined_background.pdf ) ## for drawing 
        if self.combined_signal             : self.combined_signals    .add ( self.combined_signal    .pdf ) ## for drawing 
        if self.combined_others             : self.combined_components .add ( self.combined_others    .pdf ) ## for drawing 

        for p in self.signal_components     : self.signals    .add ( p.pdf ) 
        for p in self.background_components : self.backgrounds.add ( p.pdf ) 
        for p in self.other_components      : self.components .add ( p.pdf ) 

        ## save the configuration
        self.config = {
            ## 
            'signals'               : self.signal_components     ,
            'backgrounds'           : self.background_components ,
            'others'                : self.other_components      ,
            ##
            'suffix'                : self.suffix                ,
            'name'                  : self.name                  ,            
            'extended'              : self.extended              ,
            'xvar'                  : self.xvar                  ,
            ## 
            'combine_signals'       : self.combine_signals       ,
            'combine_backgrounds'   : self.combine_backgrounds   ,
            'combine_others'        : self.combine_others        ,
            ## 
            'recursive'             : self.recursive             ,
            'recursive_signals'     : self.recursive_signals     ,
            'recursive_backgrounds' : self.recursive_backgrounds ,
            'recursive_others'      : self.recursive_others      ,
            ##
            'fS'                    : self.fS                    ,
            'fB'                    : self.fB                    ,
            'fC'                    : self.fC                    ,
            ##                      
            'S'                     : self.S                     ,
            'B'                     : self.B                     ,
            'C'                     : self.C                     ,
            'F'                     : self.F                     ,
            ##
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
    def recursive ( self ) :
        """``recursive'':  use recursive fit fractions fro non-extended fit?"""
        return  self.__recursive
    @property
    def recursive_signals ( self ) :
        """``recursive_signals'':  use recursive fractions for combined signal?"""
        return  self.__recursive_signals
    @property
    def recursive_backgrounds ( self ) :
        """``recursive_backgrounds'':  use recursive fractions for combined background?"""
        return  self.__recursive_backgrounds
    @property
    def recursive_others      ( self ) :
        """``recursive_backgrounds'':  use recursive fractions for combined other components?"""
        return  self.__recursive_others

    @property
    def signal_components ( self )  :
        """``signal_components'' : all ``signal'' components"""
        return self.__signal_components         
    @property
    def background_components ( self )  :
        """``background_components'' : all ``background'' components"""
        return self.__background_components 
    @property
    def other_components ( self )  :
        """``other_components'' : all ``other'' components"""
        return self.__other_components 

    @property
    def combined_signal     ( self ) :
        """``combined_signal'' : PDF for combined ``signal'' component"""
        return self.__combined_signal
    @property
    def combined_background ( self ) :
        """``combined_background'' : PDF for combined ``background'' component"""
        return self.__combined_background
    @property
    def combined_others    ( self ) :
        """``combined_background'' : PDF for combined ``others'' component"""
        return self.__combined_others
    
    @property
    def fS ( self ) :
        """``fS'' : fractions (possible recursive) of components in combined signal"""
        return () if not self.combined_signal else self.combined_signal.F
    @fS.setter
    def fS ( self , value ) :
        assert  self.combined_signal, "``fS'': no combined signal is defined!"
        self.combined_signal.F = value

    @property
    def fB ( self ) :
        """``fB'' : fractions (possible recursive) of components in combined background"""
        return () if not self.combined_background else self.combined_background.F
    @fB.setter
    def fB ( self , value ) :
        assert  self.combined_background, "``fB'': no combined background is defined!"
        self.combined_background.F = value

    @property
    def fC ( self ) :
        """``fC'' : fractions (possible recursive) of components in combined ``others''"""
        return () if not self.combined_others else self.combined_others.F
    @fC.setter
    def fC ( self , value ) :
        assert  self.combined_others, "``fC'': no combined ``others'' is defined!"
        self.combined_others.F = value

    @property
    def fit_components  ( self ) :
        """``fit_components'' : list of fit components"""
        return self.__fit_components
    @property
    def fit_signals      ( self ) :
        """``fit_signals'': list of (the 1st order) ``signal'' components in the model"""
        return self.__fit_signals 
    @property
    def fit_backgrounds  ( self ) :
        """``fit_backgrounds'': list of (the 1st order) ``background'' components in the model"""
        return self.__fit_backgrounds
    @property
    def fit_others       ( self ) :
        """``fit_others'': list of (the 1st order) ``others'' components in the model"""
        return self.__fit_others 

    @property
    def signal ( self ) :
        """``signal'' : get ``signal'' PDF (``combined_signal'' or the 1st from ``signal_components'')"""
        if self.__combined_signal : return self.__combined_signal
        assert self.signal_components, "signal: empty lst of ``signal'' components!"    
        if 1 != len ( self.signal_components ) :
            logger.warning ("signal: get the 1st ``signal'' component")
        return self.signal_components[0]
    @property
    def background ( self ) :
        """``background'' : get ``background'' PDF (``combined_background'' or the 1st from ``background_components'')"""
        if self.combined_background : return selfcombined_background
        assert self.background_components, "background: empty lst of ``background'' components!"    
        if 1 != len ( self.background_components ) :
            logger.warning ("background: get the 1st ``background'' component")
        return self.background_components[0]
        
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
        return () if not self.extended else self.component_getter ( self.__S )    
    @S.setter
    def S (  self , value ) :
        assert self.extended, "``S'' cannot be set for non-exteded model!"
        self.component_setter ( self.__S , value )

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
        return () if not self.extended else self.component_getter ( self.__B ) 
    @B.setter
    def B (  self , value ) :
        assert self.extended, "``B'' cannot be set for non-extended model!"
        self.component_setter ( self.__B , value )

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
        return () if not self.extended else self.component_getter ( self.__C )     
    @C.setter
    def C (  self , value ) : 
        assert self.extended, "``C'' cannot be set for non-extended model!"
        self.component_setter ( self.__C , value )

    @property
    def yields ( self ) :
        """``yields'' : yields for ``all'' components, same as ``S+B+C'', emtpy for non-extended fits"""
        return () if not self.extended else self.__S + self.__B + self.__C 

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
        return () if self.extended else self.component_getter ( self.__F )     
    @F.setter
    def F (  self , value ) :
        assert not self.extended, "``F'' cannot be set for extended model!"        
        self.component_setter ( self.__F , value )

    @property
    def natural ( self ) :
        """Are all yields natural? """
        if not self.yeilds : return False
        for y in self.yields :
            if not isinstance  ( y , ROOT.RooRealVar ) : return False 
        return True 
        
    def total_yield ( self ) :
        """``total_yield''' : get the total yield if/when possible"""
        if not self.extended                                   : return None 
        if not self.fit_result                                 : return None
        if not valid_pointer ( self.fit_result )               : return None
        yields = self.yields
        if not yields                                          : return None
        if not self.natural                                    : return None 
        if 1 ==  len ( yields )                                : return yields[0].value  
        return self.fit_result.sum ( *yields ) 
 
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END 
# =============================================================================
