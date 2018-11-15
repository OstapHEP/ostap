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
    'H1D_pdf'       , ## convertor of 1D-histo to RooHistPdf 
    ##
    )
# =============================================================================
import ROOT, math,  random
import ostap.fitting.roofit 
from   ostap.core.core      import cpp , Ostap , VE , hID , dsID , rootID, valid_pointer
from   ostap.fitting.roofit import SETVAR, PDF_fun
from   ostap.logger.utils   import roo_silent   , rootWarning , numcpu 
from   ostap.fitting.utils  import ( RangeVar   , fitArgs  , MakeVar  , 
                                     fit_status , cov_qual , Adjust1D , H1D_dset ) 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.basic' )
else                       : logger = getLogger ( __name__              )
# =============================================================================        
## @class PDF
#  The helper base class for implementation of various PDF-wrappers 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-08-21
class PDF (MakeVar) :
    """Useful helper base class for implementation of various PDF-wrappers 
    """
    __name = 'Unnamed'
    def __init__ ( self , name ,  xvar = None , special = False ) :
        
        self.__name         = name        
        self.__signals      = ROOT.RooArgList ()
        self.__backgrounds  = ROOT.RooArgList ()
        self.__components   = ROOT.RooArgList ()
        self.__crossterms1  = ROOT.RooArgSet  () 
        self.__crossterms2  = ROOT.RooArgSet  () 
        ## take care about sPlots 
        self.__splots       = []
        self.__histo_data   = None
        self.__draw_var     = None
        self.__special      = True if special else False 
        self.__fit_result   = None
        self.__vars         = ROOT.RooArgSet  ()
        self.__tricks       = True
        
        if   isinstance ( xvar , ROOT.TH1   ) : xvar = xvar.xminmax()
        elif isinstance ( xvar , ROOT.TAxis ) : xvar = xvar.GetXmin() , xvar.GetXmax()
            
        self.__xvar = None 
        ## create the variable 
        if isinstance ( xvar , tuple ) and 2 == len(xvar) :  
            self.__xvar = self.make_var ( xvar               , ## var 
                                          'x'                , ## name 
                                          'x-variable(mass)' , ## title/comment
                                          None               , ## fix ? 
                                          *xvar              ) ## min/max 
        elif isinstance ( xvar , ROOT.RooAbsReal ) :
            self.__xvar = self.make_var ( xvar               , ## var 
                                          'x'                , ## name 
                                          'x-variable/mass'  , ## title/comment
                                          fix = None         ) ## fix ? 
        else :
            self.warning('PDF : ``x-variable''is not specified properly %s/%s' % ( xvar , type ( xvar ) ) )
            self.__xvar = self.make_var( xvar , 'x' , 'x-variable' )
            
        self.__alist1     = ROOT.RooArgList()
        self.__alist2     = ROOT.RooArgList()
        self.__adjustment = None
        self.__config     = {}
        self.__pdf        = None

        self.vars.add ( self.__xvar ) 

        self.config = { 'name' : self.name , 'xvar' : self.xvar }
        
    def xminmax ( self ) :
        """Min/max values for x-variable"""
        return self.__xvar.minmax()

    @property 
    def vars ( self ) :
        """``vars'' : fitting variables/observables (as ROOT.RooArgSet)"""
        return self.__vars
    
    @property
    def value ( self ) :
        """``value''  :  get the value of PDF"""
        v = float ( self )
        if self.fit_result :
            e = self.pdf.getPropagatedError ( self.fit_result )
            if 0<= e : return  VE ( v ,  e * e )           ## RETURN
        return  v
 
    @property 
    def xvar ( self ) :
        """``x''-variable for the fit (same as ``x'')"""
        return self.__xvar

    @property 
    def x    ( self ) :
        """``x''-variable for the fit (same as ``xvar'')"""
        return self.__xvar

    @property
    def config ( self ) :
        """The full configuration info for the PDF"""
        conf = {}
        conf.update ( self.__config )
        return conf
    @config.setter
    def config ( self , value ) :
        conf = {}
        conf.update ( value )
        self.__config = conf

    @property
    def special ( self ) :
        """``special'' : is this PDF ``special''   (does nor conform some requirements)?"""
        return self.__special

    @property
    def tricks ( self ) :
        """``tricks'' : flag to allow some  tricks&shortcuts """
        return self.__tricks
    @tricks.setter
    def tricks ( self , value ) :
        val = True if value else False 
        if val and not self.__tricks :
            raise ValueError("Can't allow tricks&shortcuts!")
        self.__tricks = val

    @property
    def pdf  ( self ) :
        """The actual PDF (ROOT.RooAbsPdf)"""
        return self.__pdf
    @pdf.setter
    def pdf  ( self , value ) :
        assert isinstance ( value , ROOT.RooAbsReal ) , "``pdf'' is not ROOT.RooAbsReal"
        if not self.special :
            assert isinstance ( value , ROOT.RooAbsPdf ) , "``pdf'' is not ROOT.RooAbsPdf"
        self.__pdf = value 

    @property
    def fit_result ( self ) :
        """``fit_result'' : (the latest) fit resut (TFitResult)"""
        return self.__fit_result
    @fit_result.setter
    def fit_result ( self , value ) :
        assert value is None or isinstance ( value , ROOT.RooFitResult ) , \
               "Invalid value: %s/%s" % ( value , type ( value ) )
        self.__fit_result = None
        if isinstance ( value , ROOT.RooFitResult ) and valid_pointer ( value ) : 
            self.__fit_result = value
    
    @property
    def name ( self ) :
        """The name of the PDF"""
        return self.__name
    @name.setter
    def name ( self , value ) :
        assert isinstance ( value , str ) , "``name'' must  be a string, %s/%s is given" % ( value , type(value) ) 
        self.__name = value
    
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
    def signals     ( self ) :
        """The list/ROOT.RooArgList of all ``signal'' components, e.g. for visualization"""
        return self.__signals
    @property
    def backgrounds ( self ) :
        """The list/ROOT.RooArgList of all ``background'' components, e.g. for visualization"""
        return self.__backgrounds 
    @property
    def components  ( self ) :
        """The list/ROOT.RooArgList of all ``other'' components, e.g. for visualization"""
        return self.__components      
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
    def  adjustment ( self ) :
        """``adjustement'' object for the pdf (``None'' if no adjustment was performed)"""
        return self.__adjustment 

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
    def draw_var ( self ) :
        """``draw_var''  :  varibale to be drawn if not specified explicitely"""
        return self.__draw_var
    @draw_var.setter
    def draw_var ( self , value ) :
        assert value is None or isinstance ( value , ROOT.RooAbsReal ) , \
               "``draw_var'' has invalid type %s/%s" % ( value , type(value) )
        self.__draw_var = value 
        
    # =========================================================================
    ## make a clone for the given PDF with optional  replacement of certain parameters
    #  @code
    #  >>> xpdf = ...
    #  >>> ypdf = xpdf.clone ( xvar = yvar ,  name = 'PDFy' ) 
    #  @endcode 
    def clone ( self , **kwargs ) :
        """Make a clone for the given PDF with the optional replacement of the certain parameters
        >>> xpdf = ...
        >>> ypdf = xpdf.clone ( xvar = yvar ,  name = 'PDFy' ) 
        """

        ## get config 
        conf = {}
        conf.update ( self.config ) 
        
        ## modify the name if the name is in config  
        if conf.has_key ('name' ) : conf['name'] += '_copy'
            
        ## update (if needed)
        conf.update ( kwargs )

        KLASS = self.__class__
        cloned = KLASS ( **conf )
        
        ## if PDF is adjusted, adjust it!
        if self.adjustment :
            cloned.adjust ( self.adjustment.fraction )
            
        return cloned 

    # =========================================================================
    ## make a copy/clone for the given PDF 
    #  @code
    #  >>> import copy 
    #  >>> xpdf = ...
    #  >>> ypdf = copy.copy ( xpdf ) 
    #  @endcode 
    def __copy__ ( self ) :
        """Make a copy/clone for the given PDF 
        >>> import copy 
        >>> xpdf = ...
        >>> ypdf = copy.copy ( xpdf ) 
        """
        return self.clone()

    # =========================================================================
    ## adjust PDF a little bit to avoid zeroes
    #  A tiny  ``flat'' component is added and the orginal PDF is replaced by a new compound PDF.
    #  The fraction of added  component is fixed and defined by ``value''
    #  @code
    #  >>> pdf = ...
    #  >>> pdf.adjust ( 1.e-6 )
    #  @endcode
    #  The  fraction can be changed and/or relesed:
    #  @code
    #  >>> pdf.adjustment.fraction = 1.e-4    ## release it
    #  >>> pdf.adjustment.fraction.release()  ## allow to  vary in the fit 
    #  @endcode 
    #  The original PDF is stored as:
    #  @code
    #  >>> orig_pdf = pdf.adjustment.old_pdf 
    #  @endcode
    def adjust ( self , value =  1.e-5 ) :
        """``adjust'' PDF a little bit to avoid zeroes
        A tiny  ``flat'' component is added and the orginal PDF is replaced by a new compound PDF.
        The fraction of added  component is fixed and defined by ``value''

        >>> pdf = ...
        >>> pdf.adjust ( 1.e-6 )

        The  fraction can be changed/relesed

        >>> pdf.adjustment.fraction = 1.e-4    ## change the value 
        >>> pdf.adjustment.fraction.release()  ## release it, allow to vary in the fit 
        
        The original PDF is stored as:
        
        >>> orig_pdf = pdf.adjustment.old_pdf 
        
        """
        if self.adjustment :
            self.warning ( "PDF is already adjusted, skip it!")
            return

        ## create adjustment object and  use it to adjust PDF:
        self.__adjustment = Adjust1D ( self.name , self.xvar , self.pdf , value )
        ## replace original PDF  with  adjusted one:
        self.pdf          = self.__adjustment.pdf

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
                args   = ()    , **kwargs ) :
        """
        Perform the actual fit (and draw it)
        >>> r,f = model.fitTo ( dataset )
        >>> r,f = model.fitTo ( dataset , weighted = True )    
        >>> r,f = model.fitTo ( dataset , ncpu     = 10   )    
        >>> r,f = model.fitTo ( dataset , draw = True , nbins = 300 )    
        """
        if isinstance ( dataset , ROOT.TH1 ) :
            density = kwargs.pop ( 'density' , True   )
            chi2    = kwargs.pop ( 'chi2'    , False  )
            return self.fitHisto ( dataset           ,
                                   draw    = draw    ,
                                   silent  = silent  ,
                                   density = density ,
                                   chi2    = chi2    , 
                                   args    = args    , **kwargs ) 
        #
        ## treat the arguments properly
        #
        opts = fitArgs ( "PDF(%s).fitTo:" % self.name , dataset , *args , **kwargs )

        #
        ## define silent context
        with roo_silent ( silent ) :
            self.fit_result = None             
            result          =  self.pdf.fitTo ( dataset , ROOT.RooFit.Save () , *opts ) 
            self.fit_result = result 
            if hasattr ( self.pdf , 'setPars' ) : self.pdf.setPars() 

        if not valid_pointer (  result ) :
            self.fatal ( "fitTo: RooFitResult is invalid. Check model&data" )
            self.fit_result = None             
            return None , None
        
        st = result.status()
        if 0 != st and silent :
            self.warning ( 'fitTo: status is %s. Refit in non-silent regime ' % fit_status ( st ) )
            return self.fitTo ( dataset ,
                                draw   = draw  ,
                                nbins  = nbins ,
                                silent = False ,
                                refit  = refit ,
                                args   = args  , **kwargs )
        
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
            if isinstance ( refit , ( int , long ) )  : refit -= 1
            else                                      : refit  = False
            return  self.fitTo ( dataset         ,
                                 draw   = draw   ,
                                 nbins  = nbins  ,
                                 silent = silent ,
                                 refit  = refit  ,
                                 args   = args   , **kwargs ) 

        ## draw it if requested
        from ostap.plotting.fit_draw import draw_options
        draw_opts = draw_options ( **kwargs )
        if draw_opts and not draw     : draw = draw_opts
        if isinstance ( draw , dict ) : draw_opts.update( draw )
        
        frame = self.draw ( dataset , nbins = nbins , silent = silent , **draw_opts ) if draw else None 

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
    def _draw ( self , what , frame , options , base_color , step_color = 1 ) :
        """ Helper method to draw set of components
        """
        i = 0
        for cmp in what : 
            cmps = ROOT.RooArgSet( cmp )
            if 0 <= base_color : 
                self.pdf .plotOn (
                    frame ,
                    ROOT.RooFit.Components ( cmps                        ) ,
                    ROOT.RooFit.LineColor  ( base_color + i * step_color ) , *options )
            else :
                self.pdf .plotOn (
                    frame ,
                    ROOT.RooFit.Components ( cmps ) , *options )
                
            i += 1
                                            
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
    #  @param signal_options        drawing options for `signal'     components    
    #  @param background_options    drawing options for `background' components 
    #  @param component_options     drawing options for 'other'      components 
    #  @param fit_options           drawing options for fit curve    
    #  @param base_signal_color     base color for signal components 
    #  @param base_background_color base color for background components
    #  @param base_component_color  base color for other components 
    #  @param data_overlay          draw points atop of fitting curves?  
    #  @see ostap.plotting.fit_draw 
    def draw ( self ,
               dataset               = None ,
               nbins                 = 100  ,   ## Frame binning
               silent                = True ,   ## silent mode ?               
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
        - background_options      ## drawing options for background component(s)
        - crossterm1_options      ## drawing options for crossterm1 component(s)
        - crossterm2_options      ## drawing options for crossterm2 component(s)
        - signal_options          ## drawing options for signal     component(s)
        - component_options       ## drawing options for other      component(s)
        - total_fit_options       ## drawing options for the total fit curve
        
        Colors:                  
        - base_background_color   ## base color for background component(s)
        - base_crossterm1_color   ## base color for crossterm1 component(s)
        - base_crossterm2_color   ## base color for crossterm2 component(s)
        - base_signal_color       ## base color for signal     component(s)
        - base_component_color    ## base color for other      component(s)

        Other options:
        -  residual               ## make also residual frame
        -  pull                   ## make also residual frame
        
        For default values see ostap.plotting.fit_draw 
        """
        #
        import ostap.plotting.fit_draw as FD
        #
        ## again the context
        # 
        with roo_silent ( silent ) :

            drawvar = self.draw_var if self.draw_var else self.xvar  

            if nbins :  frame = drawvar.frame ( nbins )
            else     :  frame = drawvar.frame ()
            
            #
            ## draw invizible data (for normalzation of fitting curves)
            #
            data_options = kwargs.pop ( 'data_options' , FD.data_options )
            if dataset and dataset.isWeighted() and dataset.isNonPoissonWeighted() : 
                data_options = data_options + ( ROOT.RooFit.DataError( ROOT.RooAbsData.SumW2 ) , )
                
            if dataset : dataset .plotOn ( frame , ROOT.RooFit.Invisible() , *data_options )

            ## draw various ``background'' terms
            boptions = kwargs.pop (     'background_options' , FD.   background_options )
            bbcolor  = kwargs.pop (  'base_background_color' , FD.base_background_color )
            bbcstep  = kwargs.pop (  'background_step_color' , 1 )
            self._draw( self.backgrounds , frame , boptions , bbcolor , bbcstep )

            ## ugly :-(
            ct1options = kwargs.pop (     'crossterm1_options' , FD.   crossterm1_options )
            ct1bcolor  = kwargs.pop (  'base_crossterm1_color' , FD.base_crossterm1_color ) 
            ct1cstep   = kwargs.pop (  'crossterm1_step_color' , 1 )
            if hasattr ( self , 'crossterms1' ) and self.crossterms1 : 
                self._draw( self.crossterms1 , frame , ct1options , ct1bcolor , ct1cstep )

            ## ugly :-(
            ct2options = kwargs.pop (     'crossterm2_options' , FD.   crossterm1_options )
            ct2bcolor  = kwargs.pop (  'base_crossterm2_color' , FD.base_crossterm1_color )
            ct2cstep   = kwargs.pop (  'crossterm2_step_color' , 1 )         
            if hasattr ( self , 'crossterms2' ) and self.crossterms2 :
                self._draw( self.crossterms2 , frame , ct2options , ct2bcolor , ct2cstep )

            ## draw ``other'' components
            coptions   = kwargs.pop (      'component_options' , FD.    component_options  )
            cbcolor    = kwargs.pop (   'base_component_color' , FD. base_component_color  ) 
            ccstep     = kwargs.pop (   'component_step_color' , 1 )         
            self._draw( self.components , frame , coptions , cbcolor , ccstep )

            ## draw ``signal'' components
            soptions   = kwargs.pop (         'signal_options' , FD.      signal_options  )
            sbcolor    = kwargs.pop (      'base_signal_color' , FD.   base_signal_color  ) 
            scstep     = kwargs.pop (      'signal_step_color' , 1 )         
            self._draw( self.signals , frame , soptions , sbcolor , scstep )

            #
            ## the total fit curve
            #
            self.pdf .plotOn ( frame , *kwargs.pop ( 'total_fit_options' , FD. total_fit_options  ) )
            
            #
            ## draw data once more
            #
            if dataset : dataset  .plotOn ( frame , *data_options )            

            #
            ## suppress ugly axis labels
            #
            if not kwargs.get( 'draw_axis_title' , False ) :  
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
                self.warning("draw: ignored unknown options: %s" % kwargs.keys() )

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
    ## fit the histogram (and draw it)
    #  @code
    #  histo = ...
    #  r,f = model.fitHisto ( histo , draw = True ) 
    #  @endcode 
    def fitHisto ( self , histo , draw = False , silent = False , density = True , chi2 = False , args = () , **kwargs ) :
        """Fit the histogram (and draw it)

        >>> histo = ...
        >>> r,f = model.fitHisto ( histo , draw = True ) 
        
        """
        with RangeVar( self.xvar , *(histo.xminmax()) ) : 
            
            ## convert it! 
            self.histo_data = H1D_dset ( histo , self.xvar , density , silent )
            data            = self.histo_data.dset
            
            if chi2 : return self.chi2fitTo ( data               ,
                                              draw    = draw     ,
                                              silent  = silent   ,
                                              density =  density ,
                                              args    = args     , **kwargs )
            else    : return self.fitTo     ( data ,
                                              draw    = draw     ,
                                              nbins   = None     , 
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
                    density = True  ,
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
                hdataset        = self.__histo_data.dset 
                histo           = dataset 
                
        with roo_silent ( silent ) : 

            lst1 = list ( fitArgs ( self.name , hdataset , *args , **kwargs ) )
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

        return result, self.draw ( hdataset , nbins = None , silent = silent , **draw_opts )

    # =========================================================================
    ## draw NLL or LL-profiles for selected variable
    #  @code
    #  model.fitTo ( dataset , ... )
    #  nll  , f1 = model.draw_nll ( 'B' ,  dataset )
    #  prof , f2 = model.draw_nll ( 'B' ,  dataset , profile = True )
    #  @endcode    
    def draw_nll ( self            ,
                   var             ,
                   dataset         ,
                   profile = False ,
                   args    = ()    , **kwargs ) :

        ## get all parametrs
        pars = self.pdf.getParameters ( dataset ) 
        assert var in pars , "Variable %s is not a parameter"   % var
        if not isinstance ( var , ROOT.RooAbsReal ) : var = pars[ var ]
        del pars 
        ##
        fargs = []
        ##
        bins  = kwargs.pop ( 'nbins' , 100 )
        if bins   : fargs.append ( ROOT.RooFit.Bins      ( bins  ) ) 
        ## 
        rng   = kwargs.pop ( 'range' , None )
        if rng    : fargs.append ( ROOT.RooFit.Range     ( *rng  ) ) 
        ##
        fargs = tuple ( fargs )
        ##
        largs = [ i for i in  args ] 
        color = kwargs.pop ( 'color' , None )
        if color  : largs.append ( ROOT.RooFit.LineColor ( color ) ) 
        ##
        style = kwargs.pop ( 'style' , None )
        if style  : largs.append ( ROOT.RooFit.LineStyle ( style ) )
        ##
        width = kwargs.pop ( 'width' , None )        
        if width  : largs.append ( ROOT.RooFit.LineWidth ( width ) ) 
        ##
        if kwargs : self.warning("draw_nll: unknown parameters, ignore: %s"    % kwargs)
        ##
        largs.append  ( ROOT.RooFit.ShiftToZero() ) 
        largs = tuple ( largs ) 
    
        nll    = self.pdf.createNLL ( dataset , ROOT.RooFit.NumCPU (  numcpu() ) ) 
        result = nll

        ## make profile? 
        if profile :
            avar    = ROOT.RooArgSet    (  var ) 
            profile = nll.createProfile ( avar )
            result  = profile
            
        ## prepare the  frame & plot 
        frame  = var.frame ( *fargs )
        result.plotOn ( frame , *largs  )

        frame.SetMinimum ( 0  )

        if not kwargs.get('draw_axis_title' , False ) : 
            frame.SetXTitle  ( '' )
            frame.SetYTitle  ( '' )
            frame.SetZTitle  ( '' )
        
        ## draw it! 
        if not ROOT.gROOT.IsBatch() :
            with rootWarning ():
                frame.draw ( kwargs.get('draw_options', '' ) )
                        
        return result , frame
    
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
        assert self.alist2, "PDF(%s) has empty ``alist2''/(list of components), no sPlot is possible" % self.name 
        
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
    #  data   = model.generate ( 100000 , varset , extended = True )     
    #  @endcode
    def generate ( self ,  nEvents , varset = None , extended = False ,  *args ) :
        """Generate toy-sample according to PDF
        >>> model  = ....
        >>> data   = model.generate ( 10000 ) ## generate dataset with 10000 events
        
        >>> varset = ....
        >>> data   = model.generate ( 100000 , varset )
        >>> data   = model.generate ( 100000 , varset , extended =  =   True )
        """
        args = args + ( ROOT.RooFit.Name ( dsID() ) , ROOT.RooFit.NumEvents ( nEvents ) )
        if  extended :
            args = args + ( ROOT.RooFit.Extended () , )
        if   not varset :
            varset = ROOT.RooArgSet( self.xvar )
        elif isinstance ( varset , ROOT.RooAbsReal ) :
            varset = ROOT.RooArgSet( varser )

        if not self.xvar in varset :
            vs = ROOT.RooArgSet()
            vs . add ( self.xvar )
            for  v in varset : vs.add ( v )
            varset = vs  

        return self.pdf.generate (  varset , *args )


    # =========================================================================
    ## simple 'function-like' interface 
    def __call__ ( self , x , error = False ) :
        """ PDF as a ``function''
        >>> pdf  = ...
        >>> x = 1
        >>> y = pdf ( x ) 
        """
        if isinstance ( self.xvar , ROOT.RooRealVar ) :
            mn,mx = self.xminmax()
            if mn <= x <= mx :
                with SETVAR( self.xvar ) :
                    self.xvar.setVal ( x )
                    v = self.pdf.getVal ( self.vars )  
                    if error and self.fit_result :
                        e = self.pdf.getPropagatedError ( self.fit_result )
                        if 0<= e : return  VE ( v ,  e * e )
                    return v 
            else :
                return 0.0                    ## RETURN 
            
        raise AttributeError, 'something wrong goes here'

    # ========================================================================
    ## convert to float 
    def __float__ ( self ) :
        """Convert to float
        >>> pdf = ...
        >>> v = float ( pdf )
        """
        return self.pdf.getVal ( self.vars ) 
       
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
        for i in xrange ( nshoots ) : 
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
    # some generic stuff 
    # ========================================================================
    ## helper  function to implement some math stuff 
    def _get_stat_ ( self , funcall , *args , **kwargs ) :
        """Helper  function to implement some math stuff 
        """
        pdf         = self.pdf
        xmin , xmax = self.xminmax()
        
        if self.tricks and hasattr ( pdf , 'function' ) :    
            fun = pdf.function()
            if   hasattr ( pdf  , 'setPars'   ) : pdf.setPars()             
        else :
            fun = PDF_fun ( pdf , self.xvar , xmin , xmax )
            
        return funcall (  fun , xmin , xmax , *args , **kwargs ) 
        
    # ========================================================================
    ## get the effective RMS 
    def rms ( self ) :
        """Get the effective RMS
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'RMS: %s ' % pdf.rms()
        """
        pdf  = self.pdf 
        if   hasattr ( pdf , 'rms'            ) : return pdf.rms()        
        elif hasattr ( pdf , 'Rms'            ) : return pdf.Rms()        
        elif hasattr ( pdf , 'RMS'            ) : return pdf.RMS()        
        elif self.tricks and hasattr ( pdf , 'function' ) :
            
            fun = pdf.function()
            if   hasattr ( pdf , 'setPars'    ) : pdf.setPars() 
            if   hasattr ( fun , 'rms'        ) : return fun.rms()
            elif hasattr ( fun , 'variance'   ) : return fun.variance   ()**0.5  
            elif hasattr ( fun , 'dispersion' ) : return fun.dispersion ()**0.5 
            
        from ostap.stats.moments import rms as _rms
        return  self._get_stat_ ( _rms )

    # ========================================================================
    ## get the effective Full Width at Half Maximum
    def fwhm ( self ) :
        """Get the effective Full Width at  Half Maximum
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'FWHM: %s ' % pdf.fwhm()
        """
        ## use generic machinery 
        from ostap.stats.moments import width as _width
        w = self._get_stat_ ( _width )
        return w[1]-w[0]

    # =========================================================================
    ## get the effective Skewness
    def skewness ( self ) :
        """Get the effective Skewness
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'SKEWNESS: %s ' % pdf.skewness()
        """
        ## use generic machinery 
        from ostap.stats.moments import skewness as _skewness
        return self._get_stat_ ( _skewness )

    # =========================================================================
    ## get the effective Kurtosis
    def kurtosis ( self ) :
        """Get the effective Kurtosis
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'KURTOSIS: %s ' % pdf.kurtosis()
        """
        ## use generic machinery 
        from ostap.stats.moments import kurtosis as _kurtosis
        return self._get_stat_ ( _kurtosis )

    # =========================================================================
    ## get the effective mode 
    def mode ( self ) :
        """Get the effective mode
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'MODE: %s ' % pdf.mode()
        """
        from ostap.stats.moments import mode as _mode
        return self._get_stat_ ( _mode )

    # =========================================================================
    ## get the effective median
    def median ( self ) :
        """Get the effective median
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'MEDIAN: %s ' % pdf.median()
        """
        from ostap.stats.moments import median as _median
        return self._gets_stat_ ( _median )

    # =========================================================================
    ## get the effective mean
    def get_mean ( self ) :
        """Get the effective Mean
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'MEAN: %s ' % pdf.get_mean()
        """
        from ostap.stats.moments import mean as _mean
        return self._get_stat_ ( _mean )
    
    # =========================================================================
    ## get the effective moment for the distribution
    def moment ( self , N ) :
        """Get the effective moment
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'MOMENT: %s ' % pdf.moment( 10 )
        """
        ## use generic machinery 
        from ostap.stats.moments import moment as _moment
        return self._get_stat_ ( _moment , N ) 
    
    # =========================================================================
    ## get the effective central moment for the distribution
    def central_moment ( self , N ) :
        """Get the effective central moment
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'MOMENT: %s ' % pdf.moment( 10 )
        """
        from ostap.stats.moments import central_moment as _moment
        return self._get_stat_ ( _moment , N ) 

    # =========================================================================
    ## get the effective quantile 
    def quantile ( self , prob  ) :
        """Get the effective quantile
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'QUANTILE: %s ' % pdf.quantile ( 0.10 )
        """
        from ostap.stats.moments import quantile as _quantile
        return self._get_stat_ ( quantile , prob ) 

    # =========================================================================
    ## get the symmetric confidence interval 
    def cl_symm ( self , prob , x0 =  None ) :
        """Get the symmetric confidence interval 
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'CL :  ',  pdf.cl_symm ( 0.10 )
        """
        from ostap.stats.moments import cl_symm as _cl
        return self._get_sstat_ ( _cl , prob , x0 ) 

    # =========================================================================
    ## get the asymmetric confidence interval 
    def cl_asymm ( self , prob ) :
        """Get the asymmetric confidence interval 
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'CL :  ',  pdf.cl_asymm ( 0.10 )
        """
        from ostap.stats.moments import cl_asymm as _cl
        return self._get_sstat_ ( _cl , prob )
    
    # =========================================================================
    ## get the integral between xmin and xmax 
    def integral ( self , xmin , xmax ) :
        """Get integral between xmin and xmax
        >>> pdf = ...
        >>> print pdf.integral ( 0 , 10 )
        """
        ## check limits
        if self.xminmax() :
            mn , mx = self.xminmax() 
            xmin = max ( xmin , mn )  
            xmax = max ( xmax , mn )
 
        ## make a try to use analytical integral 
        if self.tricks and hasattr ( self , 'pdf' ) :
            _pdf = self.pdf 
            if hasattr ( _pdf , 'setPars'  ) : _pdf.setPars() 
            try: 
                if hasattr ( _pdf , 'function' ) :
                    _func = _pdf.function() 
                    if hasattr ( _func , 'integral' ) :
                        return _func.integral ( xmin , xmax )
            except:
                pass
            
        ## use numerical ingeration
        from ostap.math.integral import integral as _integral
        return _integral ( self , xmin , xmax )

    # =========================================================================
    ## get the derivative at  point x 
    def derivative ( self , x ) :
        """Get derivative at point x 
        >>> pdf = ...
        >>> print pdf.derivative ( 0 ) 
        """
        ## check limits 
        if hasattr ( self.xvar , 'getMin' ) and x < self.xvar.getMin() : return 0.
        if hasattr ( self.xvar , 'getMax' ) and x > self.xvar.getMax() : return 0.

        ## make a try to use analytical derivatives 
        if self.tricks  and hasattr ( self , 'pdf' ) :
            _pdf = self.pdf 
            if hasattr ( _pdf , 'setPars'  ) : _pdf.setPars() 
            try: 
                if hasattr ( _pdf , 'function' ) :
                    _func = _pdf.function() 
                    if hasattr ( _func , 'integral' ) :
                        return _func.derivative ( x )
            except:
                pass
            
        ## use numerical derivatives 
        from ostap.math.derivative import derivative as _derivatve
        return _derivative ( self , x )

    # ==========================================================================
    ## get a minimum of PDF for certain interval
    #  @code
    #  pdf = ...
    #  x   = pdf.minimum() 
    #  @endcode 
    def minimum ( self , xmin = None , xmax = None , x0 = None ) :
        """Get a minimum of PDF for certain interval
        >>> pdf = ...
        >>> x = pdf.minimum()
        """
        if xmin is None : xmin = self.xminmax()[0]
        if xmax is None : xmax = self.xminmax()[1]
        if self.xminmax() :
            xmin =  max ( xmin , self.xminmax()[0] )
            xmax =  min ( xmax , self.xminmax()[1] )
            
        if x0 is None           : x0 = 0.5 * ( xmin + xmax )
        
        if not xmin <= x0 <= xmax :
            logger.error("Wrong xmin/x0/xmax: %s/%s/%s"   % ( xmin , x0 , xmax ) )
        
        from ostap.math.minimize import sp_minimum_1D
        return sp_minimum_1D (  self , xmin , xmax , x0 )

    # ==========================================================================
    ## get a maximum of PDF for certain interval
    #  @code
    #  pdf = ...
    #  x   = pdf.maximum() 
    #  @endcode 
    def maximum ( self , xmin = None , xmax = None , x0 = None ) :
        """Get a maximum of PDF for certain interval
        >>> pdf = ...
        >>> x = pdf.maximum()
        """
        if xmin is None : xmin = self.xminmax()[0]
        if xmax is None : xmax = self.xminmax()[1]
        if self.xminmax() :
            xmin =  max ( xmin , self.xminmax()[0] )
            xmax  = min ( xmax , self.xminmax()[1] )
            
        if x0 is None           : x0 = 0.5 * ( xmin + xmax )

        if not xmin <= x0 <= xmax :
            logger.error("Wrong xmin/x0/xmax: %s/%s/%s"   % ( xmin , x0 , xmax ) )
        
        from ostap.math.minimize import sp_maximum_1D
        return sp_maximum_1D (  self , xmin , xmax , x0 )

        
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

    # ==========================================================================
    ## Convert PDF to the 1D-histogram
    #  @code
    #  pdf = ...
    #  h1  = pdf.histo ( 100 , -1 , 10 ) ## specify histogram parameters
    #  histo_template = ...
    #  h2  = pdf.histo ( histo = histo_template ) ## use historgam template
    #  h3  = pdf.histo ( ... , integral = True  ) ## use PDF integral within the bin  
    #  h4  = pdf.histo ( ... , density  = True  ) ## convert to "density" histogram 
    #  @endcode
    def histo ( self             ,
                nbins    = 100   , xmin = None , xmax = None ,
                hpars    = ()    , 
                histo    = None  ,
                intergal = False ,
                errors   = False , 
                density  = False ) :
        """Convert PDF to the 1D-histogram
        >>> pdf = ...
        >>> h1  = pdf.histo ( 100 , 0. , 10. ) ## specify histogram parameters
        >>> histo_template = ...
        >>> h2  = pdf.histo ( histo = histo_template ) ## use historgam template
        >>> h3  = pdf.histo ( ... , integral = True  ) ## use PDF integral within the bin  
        >>> h4  = pdf.histo ( ... , density  = True  ) ## convert to 'density' histogram 
        """
        
        import ostap.histos.histos

        # histogram is provided 
        if histo :
            
            assert isinstance ( histo  , ROOT.TH1 ) and not isinstance ( ROOT.TH2 ) , \
                   "Illegal type of ``histo''-argument %s" % type( histo )
            
            histo = histo.clone()
            histo.Reset()

        # arguments for the histogram constructor 
        elif hpars :
            
            from ostap.core.core import hID
            histo = ROOT.TH1F ( hID() , 'PDF%s' % self.name , *hpars  )
            if not histo.GetSumw2() : histo.Sumw2()

        # explicit construction from (#bins,min,max)-triplet  
        else :
            
            assert isinstance ( nbins , ( int , long) ) and 0 < nbins, \
                   "Wrong ``nbins''-argument %s" % nbins 
            if xmin == None and self.xminmax() : xmin = self.xminmax()[0]
            if xmax == None and self.xminmax() : xmax = self.xminmax()[1]
            
            from ostap.core.core import hID
            histo = ROOT.TH1F ( hID() , 'PDF%s' % self.name , nbins , xmin , xmax )
            if not histo.GetSumw2() : histo.Sumw2()


        # loop over the historgam bins 
        for i,x,y in histo.iteritems() :

            xv , xe = x.value() , x.error()
            
            # value at the bin center 
            c = self ( xv , error = errors ) 

            if not integral : 
                histo[i] = c
                continue

            # integral over the bin 
            v  = self.integral( xv -  xe , xv + xe )
            
            if errors :
                if    0 == c.cov2 () : pass
                elif  0 != c.value() and 0 != v : 
                    v = c * ( v / c.value() )
                    
            histo[i] = v 

        ## coovert to density historgam, if requested 
        if density : histo =  histo.density()
        
        return histo
    
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
    def make_fracs ( self , N , pname , ptitle , fractions = True , recursive = True )  :
        """Make list of variables/fractions for compound PDF
        """
        assert isinstance ( N , ( int , long ) ) and 2 <= N , \
               "PDF.make_fracs: Invalid N=%s/%s" % ( N, type ( N ) )
        ##
        fracs    = [] 
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
            f = self.make_var ( None , pname % i , ptitle % i , None , value , *vminmax ) 
            fracs.append ( f )
            
        return fracs
    
    # =============================================================================
    ## helper function to build composite (non-extended) PDF from components 
    def add_pdf ( self , pdflist , name , title , fname , ftitle , recursive = True ) :
        """Helper function to build composite (non-extended) PDF from components 
        """
        ##
        pdfs   = ROOT.RooArgList() 
        for pdf in pdflist : pdfs.add  ( pdf )
        fs     = self.make_fracs ( len ( pdfs ) , fname , ftitle ,
                                  fractions = True , recursive = recursive )
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
            mn, mx = self.xminmax()
            dm     =  mx - mn
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
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2017-07-13
class RESOLUTION(MASS) :
    """Helper base class  to parameterize the resolution
    """
    def __init__ ( self            ,
                   name            ,
                   xvar     = None ,
                   sigma    = None , 
                   mean     = None ) : 
        with Resolution() :
            super(RESOLUTION,self).__init__ ( name  = name  ,
                                              xvar  = xvar  ,
                                              sigma = sigma ,
                                              mean  = mean  )
        ## save the configuration
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'mean'  : self.mean  ,
            'sigma' : self.sigma
            }
            
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
        
    ## simple conversion to string    
    def __str__ ( self ) :
        """Simple conversion to string"""        
        return "Flat1D(name='%s',xvar='%s')" % ( self.name , self.xvar.name )
  
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
    def __init__ ( self , pdf , xvar ,
                   name           = None  ,
                   special        = False ,
                   add_to_signals = True  ) :
        """Wrapper for generic RooFit pdf        
        >>> raw_pdf = RooGaussian   ( ...     )
        >>> pdf     = Generic1D_pdf ( raw_pdf , xvar = x )
        """
        assert isinstance ( xvar , ROOT.RooAbsReal ) , "``xvar'' must be ROOT.RooAbsReal"
        assert isinstance ( pdf  , ROOT.RooAbsReal ) , "``pdf'' must be ROOT.RooAbsReal"
        if not name : name = pdf.GetName()
        
        ## initialize the base 
        PDF . __init__ ( self , name , xvar , special = special )
        ##
        if not self.special :
            assert isinstance ( pdf  , ROOT.RooAbsPdf ) , "``pdf'' must be ROOT.RooAbsPdf"

        ## PDF itself 
        self.pdf  = pdf

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
            }
        
    @property
    def add_to_signals ( self ) :
        """``add_to_signals'' : shodul PDF be added into list of signal components?"""
        return self.__add_to_signals 
        
    ## redefine the clone method, allowing only the name to be changed
    #  @attention redefinition of parameters and variables is disabled,
    #             since it can't be done in a safe way                  
    def clone ( self , name = '' , xvar = None ) :
        """Redefine the clone method, allowing only the name to be changed
         - redefinition of parameters and variables is disabled,
         since it can't be done in a safe way          
        """
        if xvar and not xvar is self.xvar :
            raise AttributeError("Generic1D_pdf can not be cloned with different ``xvar''")
        return PDF.clone ( self , name = name ) if name else PDF.clone ( self )
    
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
                   density = True  ,
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
#  @param signal             PDF for 'signal'     component                 (Ostap/PDF or RooAbsPdf)
#  @param background         PDF for 'background' component                 (Ostap/PDF or RooAbsPdf)
#  @param othersignals       list of PDFs for other 'signal' components     (Ostap/PDF or RooAbsPdf)
#  @param otherbackgrouds    list of PDFs for other 'background' components (Ostap/PDF or RooAbsPdf)
#  @param others             list of 'other' components                     (Ostap/PDF or RooAbsPdf)
#  @param name               the name of compound PDF 
#  @param suffix             ... add this  suffix for the PDF name
#  @param extended           build 'extended' PDF
#  @param combine_signals    combine all signal components into single SIGNAL?
#  @param combine_background combine all background components into single BACKGROUND?
#  @param combine_others     combine all other components into single COMPONENT?
#  @param recirsive          use recursive fractions for compound PDF
#  @param xvar               the fitting variable, must be specified if components are given as RooAbsPdf
#  @code 
#  gauss = Gauss_pdf( ... ) 
#  pdf   = Fit1D ( signal = gauss , background = 0 ) ## Gauss as signal ans exponent as background
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-08-02
class Fit1D (PDF) :
    """The actual fit-model for generic 1D-fits
    Parameters
    - signal             : PDF for 'signal'     component                 (Ostap/PDF or RooAbsPdf)
    - background         : PDF for 'background' component                 (Ostap/PDF or RooAbsPdf)
    - othersignals       : list of PDFs for other 'signal' components     (Ostap/PDF or RooAbsPdf)
    - otherbackgrouds    : list of PDFs for other 'background' components (Ostap/PDF or RooAbsPdf)
    - others             : list of 'other' components                     (Ostap/PDF or RooAbsPdf)
    - name               : The name of compound PDF 
    - suffix             : ... add this  suffix for the PDF name
    - extended           : build 'extended' PDF
    - combine_signals    : combine all signal components into single SIGNAL?
    - combine_background : combine all background components into single BACKGROUND?
    - combine_others     : combine all other components into single COMPONENT?
    - recirsive          : use recursive fractions for compound PDF
    - xvar               : the fitting variable, must be specified if components are given as RooAbsPdf

    >>> gauss = Gauss_pdf( ... ) 
    >>> pdf   = Fit1D ( signal = gauss , background = 0 ) ## Gauss as signal ans exponent as background 
    """
    def __init__ ( self                          , 
                   signal                        ,  ## the main signal 
                   background          = None    ,  ## the main background 
                   othersignals        = []      ,  ## additional signal         components
                   otherbackgrounds    = []      ,  ## additional background     components
                   others              = []      ,  ## additional non-classified components 
                   suffix              = ''      ,  ## the suffix 
                   name                = ''      ,  ## the name 
                   extended            = True    ,  ## extended fits ?
                   combine_signals     = False   ,  ## combine signal PDFs into single "SIGNAL"     ? 
                   combine_backgrounds = False   ,  ## combine signal PDFs into single "BACKGROUND" ?            
                   combine_others      = False   ,  ## combine signal PDFs into single "COMPONENT"  ?             
                   recursive           = True    ,  ## recursive fractions for NON-extended models?
                   xvar                = None    ) : 

        ##  save all arguments 
        self.__args = {
            'signal'     : signal     , 'othersignals'     : othersignals     ,
            'background' : background , 'otherbackgrounds' : otherbackgrounds ,
            'others'     : others     ,
            'extended'   : extended   ,
            'suffix'     : suffix     , 'name'             : name             , 
            ##
            'combine_signals'     : combine_signals     ,
            'combine_backgrounds' : combine_backgrounds ,
            'combine_others'      : combine_others      ,
            'recursive'           : recursive           ,        
            }
        
        self.__suffix              = suffix
        self.__recursive           = recursive
        self.__extended            = True if extended            else False
        self.__combine_signals     = True if combine_signals     else False
        self.__combine_backgrounds = True if combine_backgrounds else False
        self.__combine_others      = True if combine_others      else False
        
        ## wrap signal if needed 
        if   isinstance ( signal , PDF )                     : self.__signal = signal
        ## if bare RooFit pdf,  fit variable must be specified
        elif isinstance ( signal , ROOT.RooAbsPdf ) and xvar :
            self.__signal = Generic1D_pdf (  signal , xvar )
        else :
            raise AttributeError("Invalid type for ``signal'': %s/%s"  % (  signal , type( signal ) ) )
        
        if not name :
            name = '%s' % self.__signal.name 
            if suffix : name += '_' + suffix

        ## Init base class
        PDF.__init__ ( self , name + suffix , self.__signal.xvar ) 
            
        
        ## create the background component 
        self.__background = self.make_bkg ( background , 'Background' + suffix , self.xvar )

        ##  keep the lists of signals and backgrounds 
        self.signals     .add ( self.__signal     .pdf )
        self.backgrounds .add ( self.__background .pdf )

        #
        ## treat additional signals
        #        
        self.__more_signals       = [] 
        for c in othersignals :
            if   isinstance ( c , PDF            ) : cc = c 
            elif isinstance ( c , ROOT.RooAbsPdf ) : cc = Generic1D_pdf ( c ,  self.xvar ) 
            else :
                self.error ('unknown signal component %s/%s, skip it!' % ( c , type(c) ) )
                continue  
            self.__more_signals.append ( cc     )
            self.signals.add           ( cc.pdf ) 
        #
        ## treat additional backgounds 
        #
        self.__more_backgrounds   = [] 
        for c in otherbackgrounds :
            if   isinstance ( c , PDF            ) : cc = c  
            elif isinstance ( c , ROOT.RooAbsPdf ) : cc = Generic1D_pdf ( cs ,  self.xvar ) 
            else :
                self.error ('unknown background component %s/%s, skip it!' % ( cc , type(cc) ) )
                continue  
            self.__more_backgrounds.append ( cc     )
            self.backgrounds.add           ( cc.pdf ) 
        #
        ## treat additional components
        #
        self.__more_components    = []
        for c in others : 
            if   isinstance ( c , PDF            ) : cc = c  
            elif isinstance ( c , ROOT.RooAbsPdf ) : cc = Generic1D_pdf ( cs ,  self.xvar ) 
            else :
                self.error ("unknown ``other''component %s/%s, skip it!" % ( cc , type(cc) ) )
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
                                                'fS_%%d%s'   % suffix ,
                                                'fS(%%d)%s'  % suffix , recursive = True )
            ## new signal
            self.__signal      = Generic1D_pdf   ( sig , self.xvar , 'SIGNAL_' + suffix )
            self.__all_signals = ROOT.RooArgList ( sig )
            self.__sigs        = sigs 
            self.__signal_fractions = fracs 
            self.verbose('%2d signals     are combined into single SIGNAL'     % len ( sigs ) )

        ## combine background components into single backhround (if needed ) 
        self.__background_fractions = () 
        if combine_backgrounds and 1 < len( self.backgrounds ) :            
            bkg , fracs , bkgs = self.add_pdf ( self.backgrounds          ,
                                                'background_'    + suffix ,
                                                'background(%s)' % suffix ,
                                                'fB_%%d%s'       % suffix ,
                                                'fB(%%d)%s'      % suffix , recursive = True )
            ## new background
            self.__background      = Generic1D_pdf   ( bkg , self.xvar , 'BACKGROUND_' + suffix )
            self.__all_backgrounds = ROOT.RooArgList ( bkg )
            self.__bkgs            = bkgs 
            self.__background_fractions = fracs 
            self.verbose ('%2d backgrounds are combined into single BACKGROUND' % len ( bkgs ) ) 

        ## combine other components into single component (if needed ) 
        self.__components_fractions = () 
        if combine_others and 1 < len( self.components ) :
            
            cmp , fracs , cmps = self.add_pdf ( self.components      ,
                                                'other_'    + suffix ,
                                                'other(%s)' % suffix ,
                                                'fC_%%d%s'  % suffix ,
                                                'fC(%%d)%s' % suffix )
            ## save old background
            self.__other          = Generic1D_pdf   ( cmp , self.xvar , 'COMPONENT_' + suffix )
            self.__all_components = ROOT.RooArgList ( cmp )
            self.__components_fractions = fracs 
            self.verbose('%2d components  are combined into single COMPONENT'    % len ( cmps ) )


        self.__nums_signals     = [] 
        self.__nums_backgrounds = [] 
        self.__nums_components  = []
        self.__nums_fractions   = []
        
        ## build models 
        if self.extended :
            
            ns = len ( self.__all_signals )
            if 1 == ns :
                sf = self.make_var ( None , "S"+suffix , "Signal"     + suffix , None , 1 , 0 , 1.e+7 )
                self.alist1.add ( self.__all_signals[0]  )
                self.__nums_signals.append ( sf ) 
            elif 2 <= ns : 
                fis = self.make_fracs ( ns , 'S{%%d}%s' % suffix ,  'S(%%d)%s'  % suffix , fractions  = False )
                for s in self.__all_signals : self.alist1.add ( s )
                for f in fis                : self.__nums_signals.append ( f ) 

            nb = len ( self.__all_backgrounds )
            if 1 == nb :
                bf = self.make_var ( None , "B"+suffix , "Background" + suffix , None , 1 , 0 , 1.e+7 )
                self.alist1.add ( self.__all_backgrounds[0]  )
                self.__nums_backgrounds.append ( bf ) 
            elif 2 <= nb :
                fib = self.make_fracs ( nb , 'B{%%d}%s' % suffix ,  'B(%%d)%s'  % suffix , fractions  = False )
                for b in self.__all_backgrounds : self.alist1.add ( b )
                for f in fib                    : self.__nums_backgrounds.append ( f ) 

            nc = len ( self.__all_components )
            if 1 == nc :
                cf = self.make_var ( None , "C"+suffix , "Component" + suffix , None , 1 , 0 , 1.e+7 )
                self.alist1.add  ( self.__all_components[0]  )
                self.__nums_components.append ( cf ) 
            elif 2 <= nc : 
                fic = self.make_fracs ( nc , 'C{%%d}%s' % suffix ,  'C(%%d)%s'  % suffix , fractions  = False )
                for c in self.__all_components : self.alist1.add ( c )
                for f in fic                   : self.__nums_components.append ( f )

            for s in self.__nums_signals     : self.alist2.add ( s ) 
            for b in self.__nums_backgrounds : self.alist2.add ( b ) 
            for c in self.__nums_components  : self.alist2.add ( c ) 
                    
        else :

            ns = len ( self.__all_signals     )
            nb = len ( self.__all_backgrounds )
            nc = len ( self.__all_components  )
            
            for s in self.__all_signals     : self.alist1.add ( s )
            for b in self.__all_backgrounds : self.alist1.add ( b )
            for c in self.__all_components  : self.alist1.add ( c )
            
            fic = self.make_fracs ( ns + nb + nc , 'f{%%d}%s' % suffix , 'f(%%d)%s'  % suffix ,
                                    fractions  = True , recursive = self.recursive )
            
            for f in fic                    : self.__nums_fractions.append ( f )   
            for f in self.__nums_fractions  : self.alist2.add ( f ) 


        self.__nums_signals     = tuple ( self.__nums_signals     )
        self.__nums_backgrounds = tuple ( self.__nums_backgrounds ) 
        self.__nums_components  = tuple ( self.__nums_components  ) 
        self.__nums_fractions   = tuple ( self.__nums_fractions   ) 

        #
        ## The final PDF
        #       
        if   self.name       : title = "model(%s)"     % self.name
        else                 : title = "model(Fit1D)"
        
        name     = name if name else ( 'model_' + self.name )
        
        pdfargs  = name , title , self.alist1 , self.alist2
        if not self.extended :
            pdfargs = pdfargs + ( True if recursive else False , ) ## RECURSIVE ? 
        self.pdf = ROOT.RooAddPdf ( *pdfargs )
        
        if self.extended : 
            self.debug ( "extended     model ``%s'' with %s/%s components"  % ( self.pdf.GetName() , len( self.alist1) , len(self.alist2) ) )
        else : 
            self.debug ( "non-extended model ``%s'' with %s/%s components"  % ( self.pdf.GetName() , len( self.alist1) , len(self.alist2) ) )


        ## save the configurtaion
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
            'combine_backrgounds' : self.combine_backgrounds ,
            'combine_others'      : self.combine_others      ,
            'recursive'           : self.recursive           ,
            'xvar'                : self.xvar                ,
            }
        
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

    @property
    def fB ( self  ) :
        """(Recursive) fractions for the compound background components (empty for simple background)"""
        lst = [ i for i in self.__background_fractions ]
        return tuple ( lst )
    
    @property
    def fC ( self  ) :
        """(Recursive) fractions for the compound ``other'' components (empty for no additional commponents case)"""
        lst = [ i for i in self.__components_fractions ]
        return tuple ( lst )

    @property
    def S ( self ) :
        """Get the  yields of signal component(s) (empty for non-extended fits)
        For single signal component:
        >>> print pdf.S         ## read the single single component 
        >>> pdf.S = 100         ## assign to it
        For multiple signal components:
        >>> print pdf.S[4]      ## read the 4th signal component 
        >>> pdf.S = 4,100       ## assign to it 
        ... or, alternatively:
        >>> print pdf.S[4]      ## read the 4th signal component 
        >>> pdf.S[4].value 100  ## assign to it         
        """
        lst = [ i for i in self.__nums_signals ]
        if not lst          : return ()     ## extended fit? 
        elif  1 == len(lst) : return lst[0] ## simple signal?
        return tuple ( lst )
    @S.setter
    def S (  self , value ) :
        _n = len ( self.__nums_signals )
        assert 1 <= _n , "No signals are defined, assignement is impossible"
        if 1 ==  _n :
            _s    = self.S 
            value = float ( value )
        else :
            index = value [0]
            assert isinstance ( index , int ) and 0 <= index < _n, "Invalid signal index %s/%d" % ( index , _n ) 
            value = float ( value[1] )
            _s    = self.S[index]
        ## assign 
        assert value in _s , "Value %s is outside the allowed region %s"  % ( value , _s.minmax() )
        _s.setVal ( value )
    
    @property
    def B ( self ) :
        """Get the  yields of background  component(s) (empty for non-extended fits)
        For single background component:
        >>> print pdf.B           ## read the single background component 
        >>> pdf.B = 100           ## assign to it 
        For multiple background components:
        >>> print pdf.B[4]        ## read the 4th background component 
        >>> pdf.B = 4,100         ## assign to it 
        ... or, alternatively:
        >>> print pdf.B[4]        ## read the 4th background component 
        >>> pdf.B[4].value 100    ## assign to it 
        """
        lst = [ i for i in self.__nums_backgrounds ]
        if not lst          : return ()     ## extended fit? 
        elif  1 == len(lst) : return lst[0] ## simple background?
        return tuple ( lst )
    @B.setter
    def B (  self , value ) :
        _n = len ( self.__nums_backgrounds )
        assert 1 <= _n , "No backgrounds are defined, assignement is impossible"
        if 1 ==  _n :
            _b    = self.B 
            value = float ( value )
        else : 
            index = value [0]
            assert isinstance ( index , int ) and 0 <= index < _n, "Invalid background index %s/%d" % ( index , _n ) 
            value = float ( value[1] )
            _b    = self.B[index]
        ## assign 
        assert value in _b , "Value %s is outside the allowed region %s"  % ( value , _b.minmax() )
        _b.setVal ( value )

    @property
    def C ( self ) :
        """Get the  yields of ``other'' component(s) (empty for non-extended fits)
        For single ``other'' component:
        >>> print pdf.C           ## read the single ``other'' component 
        >>> pdf.C = 100           ## assign to it 
        For multiple ``other'' components:
        >>> print pdf.C[4]        ## read the 4th ``other'' component 
        >>> pdf.C = 4,100         ## assign to it 
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
        _n = len ( self.__nums_components )
        assert 1 <= _n , "No ``other'' components are defined, assignement is impossible"
        if 1 ==  _n :
            _c    = self.C 
            value = float ( value )
        else : 
            index = value [0]
            assert isinstance ( index , int ) and 0 <= index < _n, "Invalid ``other'' index %s/%d" % ( index , _n ) 
            value = float ( value[1] )
            _c    = self.C[index]
        ## assign 
        assert value in _c , "Value %s is outside the allowed region %s"  % ( value , _c.minmax() )
        _c.setVal ( value )

    @property 
    def F ( self ) :
        """Get fit fractions for non-expended fits (empty for extended fits)
        For single fraction (2 fit components):
        >>> print pdf.F           ## read the single fraction 
        >>> pdf.F = 0.1           ## assign to it 
        For multiple fractions (>2 fit components):
        >>> print pdf.F[4]        ## read the 4th fraction
        >>> pdf.F = 4,0.1         ## assign to it 
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
        _n = len ( self.__nums_fractions )
        assert 1 <= _n , "No fractions are defined, assignement is impossible"
        if 1 ==  _n :
            _f    = self.F 
            value = float ( value )
        else : 
            index = value [0]
            assert isinstance ( index , int ) and 0 <= index < _n, "Invalid fraction index %s/%d" % ( index , _n ) 
            value = float ( value[1] )
            _f    = self.F[index]
        ## assign 
        assert value in _f , "Value %s is outside the allowed region %s"  % ( value , _f.minmax() )
        _f.setVal ( value )

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
