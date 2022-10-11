#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostap/fitting/pdfbasic.py
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
    'APDF1'         , ## useful base class for 1D-models
    'APDF2'         , ## useful base class for 2D-models
    'APDF3'         , ## useful base class for 3D-models
    ##
    'PDF1'          , ## useful base class for 1D-models
    'PDF2'          , ## useful base class for 2D-models
    'PDF3'          , ## useful base class for 3D-models
    ##
    'Generic1D_pdf' , ## wrapper over imported RooFit (1D)-pdf
    'Generic2D_pdf' , ## wrapper over imported RooFit (2D)-pdf
    'Generic3D_pdf' , ## wrapper over imported RooFit (3D)-pdf
    ##
    'make_pdf'      , ## helper function to make PDF
    'all_args'      , ## check that all arguments has correct type 
    ##
    )
# =============================================================================
import ostap.fitting.roofit 
import ostap.fitting.variables
import ostap.fitting.roocollections 
from   builtins                 import range
from   ostap.core.core          import ( Ostap , VE , hID , dsID , rootID,
                                         valid_pointer ,
                                         roo_silent    , rootWarning  )
from   ostap.math.base          import iszero , frexp10 
from   ostap.core.ostap_types   import ( is_integer     , string_types   , 
                                         integer_types  , num_types      ,
                                         list_types     , all_numerics   ) 
from   ostap.fitting.roofit     import SETVAR
from   ostap.fitting.utils      import ( RangeVar   , numcpu     ,
                                         make_name  , fit_status ,
                                         cov_qual   , get_i      )
from   ostap.fitting.fithelpers import H1D_dset, H2D_dset, H3D_dset , SETPARS  
from   ostap.fitting.funbasic   import FUN1, FUN2, FUN3, Fun1D , Fun2D , Fun3D 
from   ostap.utils.cidict       import select_keys
from   ostap.fitting.roocmdarg  import check_arg , nontrivial_arg , flat_args , command  
from   ostap.core.meta_info     import root_info
import ostap.histos.histos 
import ROOT, math,  random
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.basic' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
## @var bkg_types
#  shortcuts for predefined backgroud types 
bkg_types = string_types + integer_types + ( None , ) 
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
    """Are all arguments of 'good' type?
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
## @class APDF1
#  The helper MIXIN class for implementation of various PDF-wrappers
#  - it relies on <code>xvar</code> method
#  - it relies on <code>parse_args</code> method
#  - it relies on <code>warning</code> method
#  - it relies on <code>info</code> method
#  - it relies on <code>error</code> method
#  - it relies on <code>fun</code> attribute 

#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-08-21
class APDF1 ( object ) :
    """Useful helper base class for implementation of various PDF-wrappers 
    - it relies on `xvar`       method
    - it relies on `parse_args` method
    - it relies on `warning`    method
    - it relies on `info`       method
    - it relies on `error`      method
    - it relies on `fun`        attribute 
    """
    def __init__ ( self ) :

        self.__signals               = ROOT.RooArgList ()
        self.__backgrounds           = ROOT.RooArgList ()
        self.__components            = ROOT.RooArgList ()
        self.__crossterms1           = ROOT.RooArgSet  () 
        self.__crossterms2           = ROOT.RooArgSet  () 
        self.__combined_signals      = ROOT.RooArgList ()
        self.__combined_backgrounds  = ROOT.RooArgList ()
        self.__combined_components   = ROOT.RooArgList ()

        ## take care about sPlots 
        self.__splots                = []
        self.__histo_data            = None
        self.__fit_options           = () ## predefined fit options for this PDF
        
        self.__alist1                = ROOT.RooArgList()
        self.__alist2                = ROOT.RooArgList()
        
        self.__fit_result   = None

    @property
    def pdf  ( self ) :
        """The actual PDF (ROOT.RooAbsPdf)"""
        return self.fun 
    @pdf.setter
    def pdf  ( self , value ) :
        if value is None : self.fun = value
        else : 
            assert value and isinstance ( value , ROOT.RooAbsPdf ) , "'pdf' is not ROOT.RooAbsPdf"
            self.fun = value
        
    @property
    def pdf_name ( self ) :
        """'pdf_name' : get the name of the underlying RooAbsPdf (same as 'fun_name') """
        return  self.fun_name 

    @property
    def value ( self ) :
        """'value'  :  get the value of PDF"""
        v = float ( self )
        if self.fit_result :
            e = self.pdf.getPropagatedError ( self.fit_result )
            if 0 <= e : return  VE ( v ,  e * e )           ## RETURN
        return  v
    
    @property
    def special ( self ) :
        """'special' : is this PDF 'special' (does nor conform some requirements)?"""
        return ( not self.__pdf ) or not isinstance ( self.__pdf , ROOT.RooAbsPdf )
    
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
    def fit_result ( self ) :
        """'fit_result' : result of the latest call to `fitTo` method (or `None`"""
        return self.__fit_result
    @fit_result.setter
    def fit_result ( self , value ) :
        assert ( value is None ) or ( value and isinstance ( value , ROOT.RooFitResult ) and valid_pointer ( value ) ) , \
               "Invalid value for 'fit-result' object! %s/%s" %( value , type ( value ) ) 
        self.__fit_result = value
            
    @property
    def signals     ( self ) :
        """The list/ROOT.RooArgList of all 'signal' components,
        e.g. for visualization"""
        return self.__signals
    @property
    def backgrounds ( self ) :
        """The list/ROOT.RooArgList of all 'background' components,
        e.g. for visualization"""
        return self.__backgrounds 
    @property
    def components  ( self ) :
        """The list/ROOT.RooArgList of all 'other' components,
        e.g. for visualization"""
        return self.__components

    @property
    def combined_signals     ( self ) :
        """The list/ROOT.RooArgList of all combined 'signal' components,
        e.g. for visualization"""
        return self.__combined_signals
    @property
    def combined_backgrounds ( self ) :
        """The list/ROOT.RooArgList of all combined 'background' components,
        e.g. for visualization"""
        return self.__combined_backgrounds 
    @property
    def combined_components  ( self ) :
        """The list/ROOT.RooArgList of all combined 'other' components,
        e.g. for visualization"""
        return self.__combined_components

    @property 
    def crossterms1 ( self ) :
        """'cross-terms2': cross-components for multidimensional PDFs e.g.
        - Signal(x)*Background(y)           for 2D-fits,
        - Signal(x)*Signal(y)*Background(z) for 3D-fits, etc...         
        """        
        return self.__crossterms1
    
    @property
    def crossterms2 ( self ) : 
        """'cross-terms2': cross-components for multidimensional PDFs e.g.
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
            raise AttributeError("'histo_data' has invalid type %s/%s" % (   value , type(value) ) )
    @property
    def fit_options ( self ) :
        """'fit_options' : the predefined 'fitTo'-options for this PDF
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
                self.warning ( "fitTo: Neither 'SumW2Error' and 'AsymptoticError' are specified for weighted dataset!" )

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
            if result and valid_pointer ( result ) : self.fit_result = result 
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
                        self.warning ( "fitTo: variable '%s' == %s [very close (>95%%) to maximum %s]"
                                       % ( i.GetName() , i.value , imx ) )
                    elif  0  < ie and iv < imn + 0.1 * ie :
                        self.warning ( "fitTo: variable '%s' == %s [very close (<0.1sigma) to minimum %s]"
                                       % ( i.GetName() , i.value , imn ) )                        
                    elif  iv < imn + 0.01 * idx : 
                        self.debug   ( "fitTo: variable '%s' == %s [very close (< 1%%) to minimum %s]"
                                       % ( i.GetName() , i.value , imn ) )
                        
            if hasattr ( self , 'natural' ) and self.natural and not dataset.isWeighted () :

                sums = [ nsum ]
                    
                if 2 <= len ( self.yields ) : sums.append ( result.sum ( *self.yields ) )

                for ss in sums :
                    if 0 >= ss.cov2() : continue
                    nl = ss.value() - 0.50 * ss.error() 
                    nr = ss.value() + 0.50 * ss.error()
                    if not nl <= len ( dataset ) <= nr :
                        self.warning ( "fitTo: fit is problematic: 'sum' %s != %s [%+.5g/%+.5g]" % ( ss , len( dataset ) , nl , nr ) )
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
            self.info      ( "Fit result is\n%s" % result.table ( prefix = "# " ) ) 
        elif result and ( not cov2_good ) and not silent : 
            self.warning   ( "Fit result is\n%s" % result.table ( prefix = "# " ) ) 
        elif result and not silent :
            self.warning  ( "Fit result is\n%s" % result.table ( prefix = "# " ) )

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

    # =========================================================================
    ## invoke <code>model.fitTo  ( data  , *options)</code> command
    #  - merge arguments using <code>RooFit::MultiArg</code> to shorted list
    def fit_to ( self , model , data , *options ) :
        """Invoke `model.fitTo ( data , *options)` command
        - merge arguments using `ROOT.RooFit::MultiArg` to shorted list
        """
        
        NARGS = 8
        
        assert all ( isinstance ( o , ROOT.RooCmdArg ) for o in options  ), \
               "fit_to: invalid argument types: %s" % list ( options  ) 

        ## for "small" number of arguments use the standard function 
        if len ( options ) <= NARGS :
            return model.fitTo ( data , *options )
        
        from ostap.fitting.roocmdarg import command 
        cmd = command ( *options )
        
        return model.fitTo ( frame , cmd  )


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
    #  @param crossterm1_style      style(s) for "crossterm-1"   components
    #  @param crossterm2_style      style(s) for "crossterm-2"   components
    #  @param background2D_style    style(s) for "background-2D" components
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
                
            ## draw various "background" terms
            boptions     = self.draw_option ( 'background_options' , **kwargs ) 
            bbstyle      = self.draw_option (   'background_style' , **kwargs )
            self._draw( self.backgrounds , frame , boptions , bbstyle )
            used_options.add ( 'background_options' ) 
            used_options.add ( 'background_style'   ) 

            ## draw combined "background" components 
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

            ## draw "other" components
            coptions     = self.draw_option ( 'component_options' , **kwargs )
            cbstyle      = self.draw_option ( 'component_style'   , **kwargs )
            self._draw( self.components , frame , coptions , cbstyle , args )

            used_options.add ( 'component_options' ) 
            used_options.add ( 'component_style'   ) 

            ## draw combined "other" components 
            if self.combined_components :
                
                drawit   = self.draw_option ( 'draw_combined_component'    , **kwargs )
                doptions = self.draw_option ( 'combined_component_options' , **kwargs ) 
                dstyle   = self.draw_option ( 'combined_component_style'   , **kwargs )
                
                if drawit : self._draw ( self.combined_components , frame , doptions , dstyle , args )
                
            used_options.add ( 'draw_combined_component'    ) 
            used_options.add ( 'combined_component_options' ) 
            used_options.add ( 'combined_component_style'   )
            
            ## draw "signal" components
            soptions     = self.draw_option (    'signal_options'  , **kwargs )
            sbstyle      = self.draw_option (      'signal_style'  , **kwargs ) 
            self._draw( self.signals , frame , soptions , sbstyle , args )

            used_options.add ( 'signal_options' ) 
            used_options.add ( 'signal_style'   )

            ## draw combined "signals" components 
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

            self.fit_result = None 

            args_ = tuple ( lst2 + lst1  )
            #
            chi2 = ROOT.RooChi2Var ( rootID ( "chi2_" ) , "chi2(%s)" % self.name  , self.pdf , hdataset , *args_ )
            m    = ROOT.RooMinuit  ( chi2 ) 
            m.migrad   () 
            m.hesse    ()
            result = m.save ()
            ## save fit results
            if result and valid_pointer  ( result ) : 
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
                    self.debug ('draw_nll: remove drawing artefacts  at the first and the last points' ) 
                                        
            ## scale it if needed
            if 1 != sf :
                self.info ('draw_nll: apply scale factor of %.4g due to dataset weights' % sf )
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

        ## skip some artifacts from VarMaker.parse_args 
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
            self.info ( "graph_nll: minimal value of %.5g is subtracted" % ymin ) 
            graph -= ymin 

        ## scale it if needed
        if 1 != sf :
            self.info ('graph_nll: apply scale factor of %.5g due to dataset weights' % sf )
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
            self.info ( "graph_profile: minimal value of %.5g is subtracted" % ymin ) 
            graph -= ymin 

        ## scale it if needed
        if 1 != sf :
            self.info ('graph_profile: apply scale factor of %.5g due to dataset weights' % sf )
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
        """Evaluate 'significance' using Wilks' theorem via NLL
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
            if 1 != sf :  self.info ('Scale factor of %.4g is applied' % sf )
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
        """Evaluate 'significance' using Wilks' theorem via NLL
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

        self.info ( "Wilks: fixed variables: %s" % [ f.GetName()  for f in fixed] )

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
            if 1 != sf :  self.info ('Scale factor of %.4g is applied' % sf )
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
            self.info ("minuit: Error Level is redefined from %.1f to %.4g" %  ( old_errdef ,
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
               "PDF(%s) has empty 'alist2'/(list of components)" + \
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

    # ========================================================================
    ## clean some stuff 
    def clean ( self ) :
        self.__splots     = []
        self.__histo_data = None 
        self.__fit_result = None
        
    # ==========================================================================
    ## Convert PDF to the 1D-histogram  in correct way.
    #  @code
    #  pdf = ...
    #  h1  = pdf.histo ( 100 , -1 , 10 ) ## specify histogram parameters
    #  histo_template = ...
    #  h2  = pdf.histo ( histo = histo_template ) ## use histogram template
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
        >>> h2  = pdf.histo ( histo = histo_template ) ## use histogram template
        >>> h3  = pdf.histo ( ... , integral = True  ) ## use PDF integral within the bin  
        >>> h4  = pdf.histo ( ... , density  = True  ) ## convert to 'density' histogram 
        """
        
        histo = self.make_histo ( nbins = nbins ,
                                  xmin  = xmin  ,
                                  xmax  = xmax  ,
                                  hpars = hpars ,
                                  histo = histo )

        # loop over the histogram bins 
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
    #  h2  = pdf.roo_histo ( histo = histo_template ) ## use histogram template
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
            hID ()    ,
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

            

    # =========================================================================
    ## create popular 1D "background"  function
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
        """Create popular 1D 'background'  function.
        
        Possible values for 'bkg':
        
        - None or 0                               : Flat1D
        - positive integer  'N'                   : Bkg_pdf(power=N)
        - negative integer  'K'                   : PolyPos_pdf(power=abs(K))
        - any Ostap-PDF                           : PDF will be copied or cloned  
        - RooAbsPdf        'pdf'                  : Generic1D_pdf(pdf=pdf)
        - RooAbsReal       'var'                  : Bkg_pdf(power=0,tau=var)
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


    # ================================================================================
    ## Check the ranges for variables  in dataset 
    def check_ranges ( self , dataset , range = '' ) :
        """Check the ranges for varibales in dataset 
        """

        import ostap.trees.cuts
        
        cuts = '' 
        for v in self.vars :
            
            ## has range? 
            if ( hasattr ( v , 'hasMin' ) and not v.hasMin() ) and \
               ( hasattr ( v , 'hasMax' ) and not v.hasMax() ) : continue
            
            if not v in dataset                                : continue
            
            vv_minmax = v.minmax ()
            if not vv_minmax : continue

            ##variable in dataset 
            vd = getattr ( dataset , v.name , None )
            if vd is None    : continue

            vd_minmax = vd.minmax()
            if not vd_minmax : continue

            vcut1 = ROOT.TCut ( '%s<%.17g'  % ( vd.name      , vd_minmax[0] ) )
            vcut2 = ROOT.TCut ( '%.17g<=%s' % ( vd_minmax[1] , vd.name      ) )
            if   cuts : cuts  = cuts | ( vcut1 | vcut2 )
            else      : cuts  =          vcut1 | vcut2 


        if dataset and not cuts : return  True

        has_entry  = dataset.hasEntry ( cuts , range ) if range else dataset.hasEntry ( cuts )
        
        return not has_entry 
    
    # =========================================================================
    ## Make PDF1 object
    #  @code
    #  pdf = ...
    #  pdf2 , xvar = pdf.make_PDF( ... , xvar = .. , )
    #  @endcode
    def make_PDF1 ( self , pdf , xvar = None , *args , **kwargs ) :
        """Make PDF1 object
        >>> pdf = ...
        >>> pdf2 , xvar = pdf.make_PDF( ... , xvar = .. , )        
        """
        if   isinstance ( pdf , PDF1 ) :
            
            assert ( not xvar ) or ( xvar is pdf.xvar ) , \
                       "make_PDF1: Invalid setting of xvar %s vs %s" % ( xvar , pdf.xvar )
                
            return pdf , pdf.xvar 
        
        elif xvar and isinstance ( pdf , ROOT.RooAbsPdf ) :
            
            return Generic1D_pdf ( pdf , xvar = xvar , *args , **kwargs ) , xvar 
        
        elif xvar and isinstance ( pdf , bkg_types ) :
            
            return self.make_bkg ( pdf , xvar = xvar , **kwargs ) , xvar 

        raise TypeError( "make_PDF1: invalid pdf/xvar %s/%s" % ( pdf , xvar ) )

    # =========================================================================
    ## helper functon to make a raw product of PDFs or RooAbsPDF objects
    def raw_product ( self , *pdfs ) :
        """Make a raw product of PDFs or RooAbsPDF objects
        """
        lpdfs = [] 
        for i , p in enumerate ( pdfs ) :
            if   p and isinstance ( p , APDF1          ) : lpdfs.append ( p.pdf )
            elif p and isinstance ( p , ROOT.RooAbsPdf ) : lpdfs.append ( p     )
            else : raise TypeError ( "Invalid type for %s component %s/%s" % ( i , p , type ( p ) ) ) 

        assert 2 <= len ( lpdfs ) , 'raw_product: there should be at leats two elements in the PDF list!'

            
        name  = self.new_name ( '*'.join ( '(%s)' % p.name for p in pdfs ) ) 
        title = 'product: ' + ( '*'.join ( '(%s)' % p.name for p in pdfs ) )

        self.aux_keep.append ( lpdfs ) 
        if 2 == len ( lpdf ) : return ROOT.RooProdPdf ( name , title , *lpdfs )
        
        plst = ROOT.RooArgList()
        for p in lpdfs : plst.add ( p )
        
        self.aux_keep.append ( plst ) 
        return ROOT.RooProdPdf ( name , title , plst  )
        
# =============================================================================
## @class PDF1
#  The main helper base class for implementation of various 1D PDF-wrappers 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-08-21
class PDF1(APDF1,FUN1) :
    """The main helper base class for implementation of various 1D PDF-wrappers
    """
    def __init__ ( self , name ,  xvar , tricks = True , **kwargs ) :

        FUN1  .__init__ ( self , name = name , xvar = xvar  , tricks = tricks ,
                          fun = kwargs.pop ( 'pdf' , None ) , **kwargs )
        APDF1 .__init__ ( self )
        
        self.config   = { 'name'    : self.name    ,
                          'xvar'    : self.xvar    ,
                          'tricks'  : self.tricks  ,
                          'pdf'     : self.pdf     }
        self.config.update ( kwargs )

        self.__call_OK = isinstance ( self.xvar , ROOT.RooAbsRealLValue ) 
        

    # =========================================================================
    ## simple 'function-like' interface 
    def __call__ ( self , x , error = False , normalized = True  ) :
        """ Function as a 'function'
        >>> fun  = ...
        >>> x = 1
        >>> y = fun ( x ) 
        """
        
        assert self.__call_OK , "Invalid type for xvar!"
        
        if error and not normalized :
            self.error("Can't get error for non-normalized call" )
            error = False
            
        xmnmx = self.xminmax()
        if xmnmx :
            xmn , xmx = xmnmx 
            if not xmn <= x <= xmx : return 0
        
        with SETVAR( self.xvar ) :
            
            self.xvar.setVal ( x )
            
            v = self.fun.getVal ( self.vars ) if normalized else self.fun.getVal ()  
            
            if error and self.fit_result :
                e = self.pdf.getPropagatedError ( self.fit_result )
                if 0<= e : v= VE ( v ,  e * e )

            return v

    # ========================================================================
    ## convert to float 
    def __float__ ( self ) :
        """Convert to float
        >>> pdf = ...
        >>> v = float ( pdf )
        """
        return self.fun.getVal ( self.vars ) 

    # =========================================================================
    ## make a product of two PDFs
    #  @code
    #  pdf1 = ...
    #  pdf2 = ...
    #  pdf  = pdf1 * pdf2 
    #  @endcode
    #  Rules: 
    #  - PDF3 ( x , y , z ) * PDF3 ( x , y , z ) -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF2 ( x , y )     -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF2 ( x , z )     -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF2 ( y , z )     -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF1 ( x )         -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF1 ( y )         -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF1 ( z )         -> PDF3 ( x , y , z ) 
    #  - PDF2 ( x , y )     * PDF3 ( ... )       -> process as PDF3 (...) * PDF2 ( x , y )
    #  - PDF2 ( x , y )     * PDF2 ( x , y )     -> PDF2 ( x , y )
    #  - PDF2 ( x , y )     * PDF1 ( x )         -> PDF2 ( x , y )
    #  - PDF2 ( x , y )     * PDF1 ( y )         -> PDF2 ( x , y )
    #  - PDF1 ( x )         * PDF3 ( ... )       -> process as PDF3 (...) * PDF ( x )
    #  - PDF1 ( x )         * PDF2 ( ... )       -> process as PDF2 (...) * PDF ( x )
    #  - PDF1 ( x )         * PDF1 ( x )         -> PDF1 ( x )
    #  - PDF1 ( x )         * PDF1 ( y )         -> PDF2 ( x , y )
    def __mul__ ( self , other ) :
        """Make a product of two PDFs
        
        >>> pdf1 = ...
        >>> pdf2 = ...
        >>> pdf  = pdf1 * pdf2

        Rules:
        - PDF3 ( x , y , z ) * PDF3 ( x , y , z ) -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF2 ( x , y )     -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF2 ( x , z )     -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF2 ( y , z )     -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF1 ( x )         -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF1 ( y )         -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF1 ( z )         -> PDF3 ( x , y , z ) 
        - PDF2 ( x , y )     * PDF3 ( ... )       -> process as PDF3 (...) * PDF2 ( x , y )
        - PDF2 ( x , y )     * PDF2 ( x , y )     -> PDF2 ( x , y )
        - PDF2 ( x , y )     * PDF1 ( x )         -> PDF2 ( x , y )
        - PDF2 ( x , y )     * PDF1 ( y )         -> PDF2 ( x , y )
        - PDF1 ( x )         * PDF3 ( ... )       -> process as PDF3 (...) * PDF ( x )
        - PDF1 ( x )         * PDF2 ( ... )       -> process as PDF2 (...) * PDF ( x )
        - PDF1 ( x )         * PDF1 ( x )         -> PDF1 ( x )
        - PDF1 ( x )         * PDF1 ( y )         -> PDF2 ( x , y )
        """
        from ostap.fitting.pdf_ops import pdf1_product        
        return pdf1_product ( self , other )

    # =========================================================================
    ## make a product of two PDFs
    #  @code
    #  pdf1 = ...
    #  pdf2 = ...
    #  pdf  = pdf1 * pdf2 
    #  @endcode
    #  Rules 
    #  - PDF3 ( x , y , z ) * PDF3 ( x , y , z ) -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF2 ( x , y )     -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF2 ( x , z )     -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF2 ( y , z )     -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF1 ( x )         -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF1 ( y )         -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF1 ( z )         -> PDF3 ( x , y , z ) 
    #  - PDF2 ( x , y )     * PDF3 ( ... )       -> process as PDF3 (...) * PDF2 ( x , y )
    #  - PDF2 ( x , y )     * PDF2 ( x , y )     -> PDF2 ( x , y )
    #  - PDF2 ( x , y )     * PDF1 ( x )         -> PDF2 ( x , y )
    #  - PDF2 ( x , y )     * PDF1 ( y )         -> PDF2 ( x , y )
    #  - PDF1 ( x )         * PDF3 ( ... )       -> process as PDF3 (...) * PDF ( x )
    #  - PDF1 ( x )         * PDF2 ( ... )       -> process as PDF2 (...) * PDF ( x )
    #  - PDF1 ( x )         * PDF1 ( x )         -> PDF1 ( x )    
    #  - PDF1 ( x )         * PDF1 ( y )         -> PDF2 ( x , y )
    def __rmul__ ( self , other ) :
        """Make a product of two PDFs
        
        >>> pdf1 = ...
        >>> pdf2 = ...
        >>> pdf  = pdf1 * pdf2
        
        Rules:
        - PDF3 ( x , y , z ) * PDF3 ( x , y , z ) -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF2 ( x , y )     -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF2 ( x , z )     -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF2 ( y , z )     -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF1 ( x )         -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF1 ( y )         -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF1 ( z )         -> PDF3 ( x , y , z ) 
        - PDF2 ( x , y )     * PDF3 ( ... )       -> process as PDF3 (...) * PDF2 ( x , y )
        - PDF2 ( x , y )     * PDF2 ( x , y )     -> PDF2 ( x , y )
        - PDF2 ( x , y )     * PDF1 ( x )         -> PDF2 ( x , y )
        - PDF2 ( x , y )     * PDF1 ( y )         -> PDF2 ( x , y )
        - PDF1 ( x )         * PDF3 ( ... )       -> process as PDF3 (...) * PDF ( x )
        - PDF1 ( x )         * PDF2 ( ... )       -> process as PDF2 (...) * PDF ( x )
        - PDF1 ( x )         * PDF1 ( x )         -> PDF1 ( x )
        - PDF1 ( x )         * PDF1 ( y )         -> PDF2 ( x , y )
        """
        from ostap.fitting.pdf_ops import pdf2_product, pdf3_product
        if   isinstance ( other , PDF3 ) : return pdf3_product ( other , self )
        elif isinstance ( other , PDF2 ) : return pdf2_product ( other , self )
        return NotImplemented 

    # =========================================================================
    ## make a non-extender sum of 1D PDFs
    #  @code
    #  pdf1 = ...
    #  pdf2 = ...
    #  pdf  = pdf1 + pdf2     
    #  @endcode
    def __add__ ( self , other ) :
        """Make a no-extended sum of 1D PDFs
        >>> pdf1 = ...
        >>> pdf2 = ...
        >>> pdf  = pdf1 + pdf2     
        """
        from ostap.fitting.pdf_ops import pdf1_sum        
        return pdf1_sum ( self , other )

    # =========================================================================
    ## make a non-extended sum of 1D PDFs
    #  @code
    #  pdf1 = ...
    #  pdf2 = ...
    #  pdf  = pdf1 + pdf2     
    #  @endcode
    def __radd__ ( self , other ) :
        """Make a no-extended sum of 1D PDFs
        >>> pdf1 = ...
        >>> pdf2 = ...
        >>> pdf  = pdf1 + pdf2     
        """
        from ostap.fitting.pdf_ops import pdf1_sum        
        return pdf1_sum ( self , other )

    #  ========================================================================
    ## Make a convolution for two PDFs
    #  @code 
    #  pdf    = ...
    #  other  =  ...
    #  result = pdf % other     ## Python 2 & 3 
    #  # result = pdf @ other   ## Python 3 only
    #  @endcode
    #  <code>pther</code> can be
    #  - fully configured Convolution object
    #
    #  It also can be 
    #  - resolution  PDF, it will be treated as resoltuoon function 
    #  - <code>RooAbsPdf</code>, it wil lbe treated as resoltuion function
    #  - <code>RooAbsReal</code>, it wil lbe treated as sigma for Gaussian resolution function
    #  - positive constant, it wil lbe treated as sigma for Gaussian resolution function
    #  - 2 or 3-tuple: it wil lbe treated as sigma for Gaussian resolution function
    #  
    #  The configuration can bve specified via <code>ConvolutionConfig</code>
    #  context manager:
    #  @code 
    #  pdf    = ...
    #  other  =  ...
    #  @code
    #  with ConvolutionConfig ( buffer = 0.25 , nbins = 1000 ):  
    #  ... result = pdf % other     ## python 2 & 3 
    #  ... ## result = pdf @ other  ## python 3 only 
    #  @endcode
    def __mod__ ( self , other ) :
        """ Make a convolution for two PDFs
        >>> pdf    = ...
        >>> other  =  ...
        >>> result = pdf % other     ## python 2&3
        >>> ## result = pdf % other  ## python 3 only 
        
        `Other` can be
        - fully configured Convolution object

        It also can be 
        - resolution  PDF   , it will be treated as resoltuoon function 
        - `RooAbsPdf`       , it wil lbe treated as resoltuion function
        - `RooAbsReal`      , it wil lbe treated as sigma for Gaussian resolution function
        - positive constant , it wil lbe treated as sigma for Gaussian resolution function
        - 2 or 3-tuple: it wil lbe treated as sigma for Gaussian resolution function
        
        The configuration can bve specified via `ConvolutionConfig`
        context manager:
        >>> pdf    = ...
        >>> other  =  ...
        >>> with ConvolutionConfig ( buffer = 0.25 , nbins = 1000 ):  
        >>> ... result = pdf % other    ## python 2 an d3 
        >>> ... ## result = pdf @ other ## python 3 only 
        """
        from ostap.fitting.pdf_ops import pdf_convolution
        return pdf_convolution ( self , other )

    
    def __rmod__ ( self , other ) :
        from ostap.fitting.pdf_ops import pdf_convolution
        return pdf_convolution ( self , other )
        
    __matmul__  = __mod__
    __rmatmul__ = __rmod__

    # =========================================================================
    ## Convert PDF into simple function
    #  @code
    #  pdf = ...
    #  fun = pdf.as_FUN () 
    #  @endcode
    def as_FUN ( self , name = '' ) : 
        """Convert PDF into simple function
        >>> pdf = ...
        >>> fun = pdf.as_FUN ()
        """
        return Fun1D ( self.pdf  ,
                       xvar = self.xvar  ,
                       name = name if name else self.new_name ( 'fun1' ) ) 
    
    
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
class Generic1D_pdf(PDF1) :
    """Wrapper for generic RooFit pdf
    >>> raw_pdf = RooGaussian   ( ...     )
    >>> pdf     = Generic1D_pdf ( raw_pdf , xvar = x )
    """
    ## constructor 
    def __init__ ( self , pdf , xvar = None  ,
                   name           = ''    ,
                   add_to_signals = True  ,
                   prefix         = ''    ,
                   suffix         = ''    ) :
        """Wrapper for generic RooFit pdf        
        >>> raw_pdf = RooGaussian   ( ...     )
        >>> pdf     = Generic1D_pdf ( raw_pdf , xvar = x )
        """
        assert xvar and   isinstance ( xvar , ROOT.RooAbsReal ) , "'xvar' must be ROOT.RooAbsReal"
        assert pdf  and ( isinstance ( pdf  , ROOT.RooAbsPdf  ) or \
                          ( isinstance ( pdf  , ROOT.RooAbsReal ) ) ) , \
                          "Invalid `pdf'' type`"
        
        name = ( prefix + name + suffix ) if name \
               else self.generate_name ( prefix = prefix , suffix = suffix , name = pdf.GetName() )
        
        ## initialize the base 
        PDF1 . __init__ ( self , name , xvar )
        ##

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
            'add_to_signals' : self.add_to_signals ,
            'prefix'         : prefix              ,
            'suffix'         : suffix              ,            
            }

        self.checked_keys.add  ( 'pdf'     )
        self.checked_keys.add  ( 'xvar'    )

    @property
    def add_to_signals ( self ) :
        """'add_to_signals' : should PDF be added into list of signal components?"""
        return self.__add_to_signals 

        

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
            pdf = ROOT.RooWrapperPdf        ( name , 'PDF from %s' % pdf.name , pdf )
        else :
            pdf = Ostap.MorERooFit.WrapPdf  ( name , 'PDF from %s' % pdf.name , pdf )
        
    num = len ( args )
    if   1 == num : return Generic1D_pdf ( pdf , name = name , *args )
    elif 2 == num : return Generic2D_pdf ( fun , name = name , *args )
    elif 3 == num : return Generic3D_pdf ( fun , name = name , *args )
    
    raise TypeError ( "Invalid length of arguments %s " % num ) 


# =============================================================================
## 2D Stuff goes here 
# =============================================================================
# @class APDF2
# The helper MIXIN class for implementation of 2D-pdfs 
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date 2014-08-21
class APDF2 (APDF1) :
    """ Useful helper MIXIN class for implementation of PDFs for 2D-fit
    """
    def __init__ ( self ) : 
        
        APDF1.__init__ ( self )

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
                nbins  =    50 ,
                ybins  =  None , 
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
        if   isinstance ( dataset , H2D_dset ) : dataset = dataset.dset        
        elif isinstance ( dataset , ROOT.TH2 ) :
            density = kwargs.pop ( 'density' , False ) 
            chi2    = kwargs.pop ( 'chi2'    , False ) 
            return self.fitHisto ( dataset   ,
                                   draw    = draw    ,
                                   silent  = silent  ,
                                   density = density ,
                                   chi2    = chi2    , args = args , **kwargs )
        
        ## play a bit with binning cache for convolutions 
        if self.yvar.hasBinning ( 'cache' ) :
            nb1 = self.yvar.getBins( 'cache' ) 
            yv  = getattr ( dataset , self.yvar.name , None )
            if   yv and yv.hasBinning ( 'cache' ) :
                nb2 = yv.getBins('cache')
                if  nb1 != nb2 :
                    yv.setBins ( max (  nb1 , nb2 ) , 'cache' )
                    self.info ('Adjust binning cache %s->%s for variable %s in dataset' % ( nb2 , nb1 , yv.name ) )
            elif yv :
                yv.setBins (        nb1         , 'cache' )
                self    .info ('Set binning cache %s for variable %s in dataset' %  ( nb1 , yv.name )  )
                                
        result , f = APDF1.fitTo ( self            ,
                                   dataset         ,
                                   draw   = False  , ## false here!
                                   nbins  = nbins  ,
                                   silent = silent ,
                                   refit  = refit  ,
                                   timer  = timer  , 
                                   args   = args   , **kwargs ) 
        if not draw :
            return result , None

        
        ## 2D 
        if 1 < nbins and isinstance ( ybins , integer_types ) and 1 < ybins :
            return result, self.draw ( None , dataset , nbins , ybins , silent = silent )
        
        if isinstance ( draw , str ) :
            if   draw.upper() in ( '1' , 'X' ) :
                return result, self.draw1 ( dataset , nbins = nbins , silent = silent )
            elif draw.upper() in ( '2' , 'Y' ) :
                return result, self.draw2 ( dataset , nbins = nbins , silent = silent )

        ## return 2D 
        return result, self.draw ( None , dataset , silent = silent )
    
    # =========================================================================
    ## draw the projection over 1st variable
    #
    #  @code
    #  r,f = model.fitTo ( dataset ) ## fit dataset
    #  fx  = model.draw1 ( dataset , nbins = 100 ) ## draw results
    #
    #  f1  = model.draw1 ( dataset , nbins = 100 , in_range = (2,3) ) ## draw results
    #
    #  model.yvar.setRange ( 'QUQU2' , 2 , 3 ) 
    #  f1  = model.draw1 ( dataset , nbins = 100 , in_range = 'QUQU2') ## draw results
    #
    #  @endcode 
    def draw1 ( self            ,
                dataset  = None ,
                nbins    = 100  ,
                silent   = True ,
                in_range = None ,
                args     = ()   , **kwargs ) :
        """ Draw the projection over 1st variable
        
        >>> r,f = model.fitTo ( dataset ) ## fit dataset
        >>> fx  = model.draw1 ( dataset , nbins = 100 ) ## draw results
        
        >>> f1  = model.draw1 ( dataset , nbins = 100 , in_range = (2,3) ) ## draw results

        >>> model.yvar.setRange ( 'QUQU2' , 2 , 3 ) 
        >>> f1  = model.draw1 ( dataset , nbins = 100 , in_range = 'QUQU2') ## draw results
        
        """
        if in_range and isinstance ( in_range , tuple ) and 2 == len ( in_range ) :
            range_name = 'aux2_rng2_%s' % self.name 
            with roo_silent ( silent , 3 ) : 
              self.yvar.setRange ( range_name , in_range[0] , in_range[1] )
              if dataset:
                dataset.get_var(self.yvar.GetName()).setRange ( range_name , in_range[0] , in_range[1] )

            in_range = range_name 

        return self.draw ( drawvar  = self.xvar , 
                           dataset  = dataset   ,
                           nbins    = nbins     ,
                           ybins    = 20        , ## fake 
                           silent   = silent    ,
                           in_range = in_range  ,
                           args     = args      , **kwargs )
    
    # =========================================================================
    ## draw the projection over 2nd variable
    #
    #  @code
    #  r,f = model.fitTo ( dataset ) ## fit dataset
    #  fy  = model.draw2 ( dataset , nbins = 100 ) ## draw results
    #
    #  f2  = model.draw2 ( dataset , nbins = 100 , in_range = (2,3) ) ## draw results
    #
    #  model.xvar.setRange ( 'QUQU1' , 2 , 3 ) 
    #  f2  = model.draw2 ( dataset , nbins = 100 , in_range = 'QUQU1') ## draw results
    #
    #  @endcode 
    def draw2 ( self            ,
                dataset  = None ,
                nbins    = 100  ,
                silent   = True ,
                in_range = None ,
                args     = ()   , **kwargs ) :
        """
        Draw the projection over 2nd variable
        
        >>> r,f = model.fitTo ( dataset ) ## fit dataset
        >>> fy  = model.draw2 ( dataset , nbins = 100 ) ## draw results
        
        >>> f2  = model.draw2 ( dataset , nbins = 100 , in_range = (2,3) ) ## draw results

        >>> model.xvar.setRange ( 'QUQU1' , 2 , 3 ) 
        >>> f2  = model.draw2 ( dataset , nbins = 100 , in_range = 'QUQU1') ## draw results

        """
        if in_range and isinstance ( in_range , tuple ) and 2 == len ( in_range ) :
            range_name = 'aux2_rng1_%s' % self.name 
            with roo_silent ( silent , 3 ) : 
                self.xvar.setRange ( range_name , in_range[0] , in_range[1] )
                if dataset:
                    dataset.get_var(self.xvar.GetName()).setRange ( range_name , in_range[0] , in_range[1] )
                    
            in_range = range_name

        return self.draw ( drawvar  = self.yvar ,
                           dataset  = dataset   ,
                           nbins    = nbins     ,
                           ybins    = 20        , ## fake
                           silent   = silent    ,
                           in_range = in_range  ,
                           args     = args      , **kwargs )

    # =========================================================================
    ## draw as 2D-histograms 
    def draw_H2D ( self           ,
                   dataset = None ,  
                   xbins   = 20   ,
                   ybins   = 20   ) :
        """
        Make/draw 2D-histograms 
        """
        
        _xbins = ROOT.RooFit.Binning ( xbins ) 
        _ybins = ROOT.RooFit.Binning ( ybins ) 
        _yvar  = ROOT.RooFit.YVar    ( self.yvar , _ybins )
        _clst  = ROOT.RooLinkedList  ()
        hdata  = self.pdf.createHistogram ( hID() , self.xvar , _xbins , _yvar )
        hpdf   = self.pdf.createHistogram ( hID() , self.xvar , _xbins , _yvar )
        hdata.SetTitle(';;;')
        hpdf .SetTitle(';;;')
        _lst   = ROOT.RooArgList ( self.xvar , self.yvar )  
        if dataset : dataset.fillHistogram( hdata , _lst ) 
        self.pdf.fillHistogram  ( hpdf , _lst )

        groot = ROOT.ROOT.GetROOT()
        if not groot.IsBatch() :
            with rootWarning ():
                hdata.lego ()
                hpdf .Draw ( 'same surf')
        
        return hpdf , hdata 
    
    # =========================================================================
    ## make 1D-plot
    def draw ( self                         ,
               drawvar               = None ,
               dataset               = None ,
               nbins                 =  100 ,
               ybins                 =   20 ,
               silent                = True ,
               in_range              = None ,
               args                  = ()   , 
               **kwargs                     ) : 
        """
        Make 1D-plot:
        """
        if   drawvar in ( 'x'  , 'X' , '1' , 1 , self.xvar.name ) : drawvar = self.xvar
        elif drawvar in ( 'y'  , 'Y' , '2' , 2 , self.yvar.name ) : drawvar = self.yvar

        # 
        ## special case:  do we need it? 
        #

        if drawvar is None : return self.draw_H2D( dataset , nbins , ybins )
                
        newargs = kwargs.copy ()
        
        if in_range and isinstance ( in_range , list_types ) and 2 == len ( in_range ) :
            low  = in_range [ 0 ]
            high = in_range [ 1 ]
            if isinstance ( low , num_types ) and isinstance ( high , num_types ) and low < high :
                range_name = 'aux2_range_%s' % self.name 
                with roo_silent ( true , 3 ) : 
                    drawvar.setRange ( range_name , low , high )
                    if dataset : 
                        dataset.get_var(drawvar.GetName()).setRange ( range_name , low , high )
                    in_range = range_name
                    
  #      if in_range and not isinstance ( in_range , list_types ) :
   #         in_range = in_range ,
        
        if in_range :
            options_cut =tuple ( [ ROOT.RooFit.CutRange ( in_range ) , ] )
            newargs [ 'data_options' ] = self.draw_option ( 'data_options' , **newargs ) + options_cut
            
        if in_range : 
            options_project =  tuple (  [ROOT.RooFit.ProjectionRange ( in_range ) ,] )
            for key in  ( 'total_fit_options'           ,
                          #
                          'signal_options'              ,
                          'background_options'          ,
                          'component_options'           ,
                          'crossterm1_options'          ,
                          'crossterm2_options'          ,
                          #
                          'combined_signal_options'     ,
                          'combined_background_options' ,
                          'combined_component_options'  ) :
                newargs [ key ] =  self.draw_option ( key , **newargs ) + options_project
                
        #
        ## redefine the drawing variable:
        #

        self.draw_var = drawvar
        ##
        if in_range and root_info < ( 6 , 24 ) :
            from itertools import chain 
            for p in chain ( self.signals               ,
                             self.backgrounds           ,
                             self.components            ,
                             self.crossterms1           ,
                             self.crossterms2           ,
                             self.combined_signals      ,
                             self.combined_backgrounds  ,
                             self.combined_components   ) :
                
                if   isinstance ( p , ( ROOT.RooHistPdf , ROOT.RooParamHistFunc ) ) :
                    self.warning ("'in_range' is specified, it does not work properly for ROOT<6.24")
                elif isinstance ( p , ( ROOT.RooAddPdf , ROOT.RooProdPdf ) ) : 
                    for pp in p.pdfList() :
                        if isinstance  ( pp , ( ROOT.RooHistPdf , ROOT.RooParamHistFunc ) ) :
                            self.warning ("'in_range' is specified, it does not work properly for ROOT<6.24")
                            
        #
        ## delegate the actual drawing to the base class
        #
        
        result = APDF1.draw ( self            ,
                              dataset         ,
                              nbins  = nbins  ,
                              silent = silent ,
                              args   = args   , **newargs )
        
        self.draw_var = None
        return result 
    
    # =========================================================================
    ## fit the 2D-histogram (and draw it)
    #
    #  @code
    #
    #  histo = ...
    #  r,f = model.fitHisto ( histo )
    #
    #  @endcode
    def fitHisto ( self            ,
                   histo           ,
                   draw    = False ,
                   silent  = False ,
                   density = False ,
                   chi2    = False ,
                   args    = ()    , **kwargs ) :
        """Fit the histogram (and draw it)
        
        >>> histo = ...
        >>> r,f = model.fitHisto ( histo , draw = True )
        
        """

        xminmax = histo.xminmax()
        yminmax = histo.yminmax()        
        with RangeVar( self.xvar , *xminmax ) , RangeVar ( self.yvar , *yminmax ):
            
            hdata = getattr ( self , 'histo_data' , None )
            if hdata and isinstance ( hdata  , H2D_dset ) and \
                   hdata.histo      is histo              and \
                   hdata.density    == density            and \
                   hdata.histo_hash == hash ( histo ) :
                ## reuse the existing dataset
                self.debug ('Reuse the existing H2D_dset') 
                data = hdata.dset
            else :                
                ## convert it!
                self.debug ('Create new H2D_dset'        ) 
                self.histo_data = H2D_dset ( histo , self.xvar , self.yvar  , density , silent )
                data = self.histo_data.dset 
            
            ## fit it!!
            if chi2 : return self.chi2fitTo ( data                     ,
                                              draw    = draw           ,
                                              silent  = False          ,
                                              density = density        ,
                                              args    = args           , **kwargs )
            else     : return self.fitTo    ( data                     ,
                                              draw    = draw           ,
                                              nbins   = histo.nbinsx() ,
                                              ybins   = histo.nbinsy() ,
                                              silent  = silent         ,
                                              args    = args           , **kwargs )
            
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
                   binning  = {}    ,
                   sample   = True  , 
                   args     = ()    ) :
        """Generate toy-sample according to PDF
        >>> model  = ....
        >>> data   = model.generate ( 10000 ) ## generate dataset
        
        >>> varset = ....
        >>> data   = model.generate ( 100000 , varset , sample = False )
        >>> data   = model.generate ( 100000 , varset , sample = True  )
        """
        nEvents = self.gen_sample ( nEvents ) if sample else nEvents 
        assert 0 <= nEvents , 'Invalid number of Events %s' % nEvents  
        
        args = args + ( ROOT.RooFit.Name ( dsID() ) , ROOT.RooFit.NumEvents ( nEvents ) )
        
        if binning is True :
            args    = args + ( ROOT.AllBinned() , ) 
            binning = {}

        if   not varset :
            varset = ROOT.RooArgSet( self.xvar , self.yvar )
        elif isinstance ( varset , ROOT.RooAbsReal ) :
            varset = ROOT.RooArgSet( varset )

        if not self.xvar in varset :
            vs = ROOT.RooArgSet()
            vs . add ( self.xvar )
            for  v in varset : vs.add ( v )
            varset = vs

        if not self.yvar in varset :
            vs = ROOT.RooArgSet()
            vs . add ( self.yvar )
            for  v in varset : vs.add ( v )
            varset = vs

        from ostap.fitting.variables import KeepBinning        
        with KeepBinning ( self.xvar ) , KeepBinning ( self.yvar ) : 

            if binning :
                
                xbins = binning.get ( self.xvar.name , None )
                ybins = binning.get ( self.yvar.name , None )

                if xbins : self.xvar.bins = xbins
                if ybins : self.yvar.bins = ybins
                                        
            return self.pdf.generate ( varset , *args )


    # ========================================================================
    ## check minmax of the PDF using the random shoots
    #  @code
    #  pdf     = ....
    #  mn , mx = pdf.minmax()            
    #  @endcode 
    def minmax ( self , nshoots =  100000 ) :
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

        ## check RooAbsReal functionality
        code = self.pdf.getMaxVal( ROOT.RooArgSet ( self.xvar , self.yvar ) )
        if 0 < code :
            mx = self.pdf.maxVal ( code )
            if 0 < mx : return 0 , mx
            
        ## not try  to use random
                
        mn , mx = -1 , -10
        if hasattr ( self.pdf , 'min' ) : mn = self.pdf.min()
        if hasattr ( self.pdf , 'max' ) : mx = self.pdf.max()
        if 0 <= mn and mn <= mx and 0 < mx : return mn , mx
        
        if not self.xminmax() : return ()
        if not self.yminmax() : return ()
        
        mn  , mx = -1 , -10
        xmn , xmx = self.xminmax()
        ymn , ymx = self.yminmax()
        for i in range ( nshoots ) : 
            xx = random.uniform ( xmn , xmx )
            yy = random.uniform ( ymn , ymx )
            with SETVAR ( self.xvar ) , SETVAR ( self.yvar ) :
                self.xvar.setVal ( xx )
                self.yvar.setVal ( yy )
                vv = self.pdf.getVal()
                if mn < 0 or vv < mn : mn = vv
                if mx < 0 or vv > mx : mx = vv
                    
        return mn , mx 
        

    # =========================================================================
    ## get integral over (xmin,xmax,ymin,ymax) region
    #  @code
    #  pdf = ...
    #  print ( pdf.integral( 0,1,0,2) ) 
    #  @endcode
    def integral ( self, xmin , xmax , ymin , ymax , nevents = True ) :
        """Get integral over (xmin,xmax,ymin,ymax) region
        >>> pdf = ...
        >>> print ( pdf.integral( 0,1,0,2) ) 
        """
        if self.xminmax() :
            xmn , xmx = self.xminmax()
            xmin = max ( xmin , xmn )
            xmax = min ( xmax , xmx )

        if self.yminmax() : 
            ymn , ymx = self.yminmax()            
            ymin = max ( ymin , ymn )
            ymax = min ( ymax , ymx )

        value , todo  = 0 , True 
        
        ## 1) make a try to use analytical integral (could be fast)
        if self.tricks :
            try:
                if hasattr ( self.pdf , 'setPars'  ) : self.pdf.setPars() 
                fun          = self.pdf.function()
                value , todo = fun.integral ( xmin , xmax , ymin , ymax ) , False 
            except:
                pass

        ## use numerical integration 
        from ostap.math.integral import integral2 as _integral2

        extended =  self.pdf.canBeExtended() or isinstance ( self.pdf , ROOT.RooAddPdf )

        if   todo and extended : value   = _integral2 ( self , xmin , xmax , ymin , ymax )
        elif todo  :
            
            ## use unormalized PDF here to speed up the integration 
            ifun   = lambda x, y  :  self ( x , y , error = False , normalized = False )
            value  = _integral2 ( ifun , xmin , xmax , ymin , ymax )
            norm   = self.pdf.getNorm ( self.vars )
            value /= norm

        if nevents and self.pdf.mustBeExtended () :
            evts = self.pdf.expectedEvents( self.vars )
            if evts  <= 0 or iszero ( evts ) :
                self.warning ( "integral: expectedEvents is %s" % evts )
            value *= evts 

        return value


    # ==========================================================================
    ## get a minimum of PDF for certain interval
    #  @code
    #  pdf2 = ...
    #  x ,y = pdf2.minimum() 
    #  @endcode 
    def minimum ( self ,
                  xmin = None , xmax = None ,
                  ymin = None , ymax = None , x0 = () ) :
        """Get a minimum of PDF for certain interval
        >>> pdf2 = ...
        >>> x, y = pdf2.minimum()
        """
        
        if xmin is None : xmin = self.xminmax()[0]
        if xmax is None : xmax = self.xminmax()[1]
        if self.xminmax() :
            xmin =  max ( xmin , self.xminmax()[0] )
            xmax =  min ( xmax , self.xminmax()[1] )

        if ymin is None : ymin = self.yminmax()[0]
        if ymax is None : ymax = self.yminmax()[1]
        if self.yminmax() :
            ymin =  max ( ymin , self.yminmax()[0] )
            ymax =  min ( ymax , self.yminmax()[1] )
            
        if not x0 : x0 = 0.5 * ( xmin + xmax ) , 0.5 * ( ymin + ymax )
        
        if not xmin <= x0[0] <= xmax :
            self.error("Wrong xmin/x0[0]/xmax: %s/%s/%s"   % ( xmin , x0[0] , xmax ) )

        if not ymin <= x0[1] <= ymax : 
            self.error("Wrong ymin/x0[1]/ymax: %s/%s/%s"   % ( ymin , x0[1] , ymax ) )
        
        from ostap.math.minimize import sp_minimum_2D
        return sp_minimum_2D (  self ,
                                xmin , xmax ,
                                ymin , ymax , x0 )

    # ==========================================================================
    ## get a maximum of PDF for certain interval
    #  @code
    #  pdf2 = ...
    #  x,y  = pdf2.maximum() 
    #  @endcode 
    def maximum ( self , xmin = None , xmax = None , x0 = None ) :
        """Get a maximum of PDF for certain interval
        >>> pdf2  = ...
        >>> x , y  = pdf2.maximum()
        """
        if xmin is None : xmin = self.xminmax()[0]
        if xmax is None : xmax = self.xminmax()[1]
        if self.xminmax() :
            xmin =  max ( xmin , self.xminmax()[0] )
            xmax =  min ( xmax , self.xminmax()[1] )

        if ymin is None : ymin = self.yminmax()[0]
        if ymax is None : ymax = self.yminmax()[1]
        if self.yminmax() :
            ymin =  max ( ymin , self.yminmax()[0] )
            ymax =  min ( ymax , self.yminmax()[1] )
            
        if not x0 : x0 = 0.5 * ( xmin + xmax ) , 0.5 * ( ymin + ymax )

        if not xmin <= x0[0] <= xmax :
            self.error("Wrong xmin/x0[0]/xmax: %s/%s/%s"   % ( xmin , x0[0] , xmax ) )

        if not ymin <= x0[1] <= ymax : 
            self.error("Wrong ymin/x0[1]/ymax: %s/%s/%s"   % ( ymin , x0[1] , ymax ) )

        from ostap.math.minimize import sp_maximum_2D
        return sp_maximum_2D (  self ,
                                xmin , xmax ,
                                ymin , ymax , x0 )
    
    # ==========================================================================
    ## convert PDF into TF2 object, e.g. to profit from TF2::Draw options
    #  @code
    #  pdf = ...
    #  tf2 = pdf.tf()
    #  tf2.Draw('colz')
    #  @endcode
    def tf ( self , xmin = None , xmax = None , ymin = None , ymax = None ) :
        """Convert PDF to TF2 object, e.g. to profit from TF2::Draw options
        >>> pdf = ...
        >>> tf2 = pdf.tf()
        >>> tf1.Draw('colz')
        """
        def _aux_fun_ ( x , pars = [] ) :
            return self ( x[0] , x[1] , error = False )

        if xmin == None and self.xminmax() : xmin = self.xminmax()[0]
        if xmax == None and self.xminmax() : xmax = self.xminmax()[1]
        if ymin == None and self.yminmax() : ymin = self.yminmax()[0]
        if ymax == None and self.yminmax() : ymax = self.yminmax()[1]
        
        if xmin == None : xmin = 0.0
        if xmax == None : xmin = 1.0
        if ymin == None : ymin = 0.0
        if ymax == None : ymin = 1.0
        
        from ostap.core.core import fID
        return ROOT.TF2 ( fID() , _aux_fun_ , xmin , xmax , ymin , ymax ) 

    
    # ==========================================================================
    ## Create the histo according to specifications 
    def make_histo ( self , 
                     xbins    = 20    , xmin = None , xmax = None ,
                     ybins    = 20    , ymin = None , ymax = None ,
                     hpars    = ()    , 
                     histo    = None  ) :
        """Create the histogram accordig to specifications
        """
        
        import ostap.histos.histos

        # histogram is provided 
        if histo :
            
            assert isinstance ( histo , ROOT.TH2 ) and not isinstance ( histo , ROOT.TH3 ) , \
                   "Illegal type of 'histo'-argument %s" % type( histo )
            
            histo = histo.clone()
            histo.Reset()

        # arguments for the histogram constructor 
        elif hpars :
            
            histo = ROOT.TH2F ( hID () , 'PDF%s' % self.name , *hpars  )
            if not histo.GetSumw2() : histo.Sumw2()

        # explicit construction from (#bins,min,max)-triplet  
        else :
            
            assert isinstance ( xbins , integer_types ) and 0 < xbins, \
                   "Wrong 'xbins'-argument %s" % xbins 
            assert isinstance ( ybins , integer_types ) and 0 < ybins, \
                   "Wrong 'ybins'-argument %s" % ybins 
            if xmin == None and self.xminmax() : xmin = self.xminmax()[0]
            if xmax == None and self.xminmax() : xmax = self.xminmax()[1]
            if ymin == None and self.yminmax() : ymin = self.yminmax()[0]
            if ymax == None and self.yminmax() : ymax = self.yminmax()[1]
            
            histo = ROOT.TH2F ( hID() , 'PDF%s' % self.name ,
                                xbins , xmin , xmax ,
                                ybins , ymin , ymax )
            if not histo.GetSumw2() : histo.Sumw2()

        return histo 
                     
    # ==========================================================================
    ## Convert PDF to the 2D-histogram
    #  @code
    #  pdf = ...
    #  h1  = pdf.histo ( 100 , 0. , 10. , 20 , 0. , 10 ) ## specify histogram parameters
    #  histo_template = ...
    #  h2  = pdf.histo ( histo = histo_template ) ## use histogram template
    #  h3  = pdf.histo ( ... , integral = True  ) ## use PDF integral within the bin  
    #  h4  = pdf.histo ( ... , density  = True  ) ## convert to "density" histogram 
    #  @endcode
    def histo ( self             ,
                xbins    = 20    , xmin = None , xmax = None ,
                ybins    = 20    , ymin = None , ymax = None ,
                hpars    = ()    , 
                histo    = None  ,
                integral = False ,
                errors   = False , 
                density  = False ) :
        """Convert PDF to the 2D-histogram
        >>> pdf = ...
        >>> h1  = pdf.histo ( 100 , 0. , 10. , 20 , 0. , 10 ) ## specify histogram parameters
        >>> histo_template = ...
        >>> h2  = pdf.histo ( histo = histo_template ) ## use histogram template
        >>> h3  = pdf.histo ( ... , integral = True  ) ## use PDF integral within the bin  
        >>> h4  = pdf.histo ( ... , density  = True  ) ## convert to 'density' histogram 
        """
        
        
        histos = self.make_histo ( xbins = xbins , xmin = xmin , xmax = xmax ,
                                   ybins = ybins , ymin = ymin , ymax = ymax ,
                                   hpars = hpars ,
                                   histo = histo )

        # loop over the histogram bins 
        for ix,iy,x,y,z in histo.items() :

            xv , xe = x.value() , x.error()
            yv , ye = y.value() , y.error()
            
            # value at the bin center 
            c = self ( xv , yv , error = errors ) 

            if not integral : 
                histo[ix,iy] = c
                continue

            # integral over the bin 
            v  = self.integral( xv - xe , xv + xe , yv - ye , yv + ye )
            
            if errors :
                if    0 == c.cov2 () : pass
                elif  0 != c.value() and 0 != v : 
                    v = c * ( v / c.value() )
                    
            histo[ix,iy] = v 

        ## coovert to density histogram, if requested 
        if density : histo =  histo.density()
        
        return histo


    # ==========================================================================
    ## Convert PDF to the 2D-histogram, taking taking PDF-values at bin-centres
    #  @code
    #  pdf = ...
    #  h1  = pdf.roo_histo ( 100 , 0. , 10. , 20 , 0. , 10 ) 
    #  histo_template = ...
    #  h2  = pdf.roo_histo ( histo = histo_template ) ## use histogram template
    #  h3  = pdf.roo_histo ( ... , density  = True  ) ## convert to "density" histogram 
    #  @endcode
    def roo_histo ( self             ,
                   xbins   = 20    , xmin = None , xmax = None ,
                   ybins   = 20    , ymin = None , ymax = None ,
                   hpars   = ()    , 
                   histo   = None  , 
                   events  = True  ) : 
        """Convert PDF to the 2D-histogram, taking PDF-values at bin-centres
        >>> pdf = ...
        >>> h1  = pdf.as_histo ( 100 , 0. , 10. , 20 , 0. , 10 ) 
        >>> histo_template = ...
        >>> h2  = pdf.as_histo ( histo = histo_template ) ## use histogram template
        >>> h3  = pdf.as_histo ( ... , density  = True  ) ## convert to 'density' histogram 
        """
        
        histo = self.make_histo ( xbins = xbins , xmin = xmin , xmax = xmax ,
                                  ybins = ybins , ymin = ymin , ymax = ymax ,
                                  hpars = hpars ,
                                  histo = histo )
        
        hh = self.pdf.createHistogram (
            hID()     ,
            self.xvar ,                    self.binning ( histo.GetXaxis() , 'histo2x' )   ,
            ROOT.RooFit.YVar ( self.yvar , self.binning ( histo.GetYaxis() , 'histo2y' ) ) , 
            ROOT.RooFit.Scaling  ( False ) , 
            ROOT.RooFit.Extended ( False ) ) 

        for i in hh : hh.SetBinError ( i , 0 ) 
        
        if events and self.pdf.mustBeExtended() :
            
            for ix , iy , x , y , z in hh.items() :
                volume          = 4 * x.error() * y.error() 
                hh [ ix , iy ] *= volume
                
            hh *= self.pdf.expectedEvents ( self.vars ) / hh.sum() 
                
        histo += hh
        
        return histo 

    
    # ==========================================================================
    ## get the residual histogram : (data-fit) 
    #  @see PDF.as_histo
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
        - see PDF.as_histo
        - see PDF.residual_histo
        - see PDF.make_histo

        >>> data = ...
        >>> pdf  = ...
        >>> pdf.fitTo ( data )
        >>> residual = pdf.residual ( data , nbins = 100 ) 
        """
        hdata = self.make_histo ( **kwargs )
        dataset.project ( hdata , ( self.yvar.name , self.xvar.name )  )
        return self.residual_histo ( hdata ) 
        
    # ==========================================================================
    ## get the pull histogram : (data-fit)/data_error 
    #  @see PDF.as_histo
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
        - see PDF.as_histo
        - see PDF.residual_histo
        - see PDF.make_histo

        >>> data = ...
        >>> pdf  = ...
        >>> pdf.fitTo ( data )
        >>> residual = pdf.residual ( data , nbins = 100 ) 
        """
        hdata = self.make_histo ( **kwargs )
        dataset.project ( hdata , ( self.yvar.name , self.xvar.name ) ) 
        return self.pull_histo ( hdata ) 
        
    ## conversion to string 
    def __str__ (  self ) :
        return '%s(%s,xvar=%s,yvar=%s)' % (
            self.__class__.__name__ , self.name , self.xvar.name , self.yvar.name )
    __repr__ = __str__ 

    ## Make PDF2 object 
    def make_PDF2 ( self , pdf , xvar = None , yvar = None , *args , **kwargs ) :
        """Make PDF1 object
        """
        if   isinstance  ( pdf , PDF2 ) :
            
            assert ( not xvar ) or ( xvar in pdf.vars ) , \
                   "make_PDF2: Invalid setting of xvar %s vs %s" % ( xvar , pdf.xvar )
            assert ( not yvar ) or ( yvar in pdf.vars ) , \
                   "make_PDF2: Invalid setting of yvar %s vs %s" % ( yvar , pdf.yvar )
   
            return pdf, pdf.xvar, pdf.yvar
        
        elif isinstance ( pdf , ROOT.RooAbsPdf   ) and xvar and yvar :
            
            return Generic2D_pdf ( pdf , xvar = xvar , yvar = yvar , *args , **kwargs ) , xvar , yvar

        raise TypeError( "make_PDF2: invalid pdf/xvar %s/%s" % ( pdf , xvar ) )

# =============================================================================
## @class PDF2
#  The main helper base class for implementation of various 1D PDF-wrappers 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-08-21
class PDF2(APDF2,FUN2) :
    """The main helper base class for implementation of various 1D PDF-wrappers
    """
    def __init__ ( self , name ,  xvar , yvar , tricks = True , **kwargs ) :
        
        FUN2  .__init__ ( self , name = name , xvar = xvar , yvar = yvar , tricks = tricks ,
                          fun = kwargs.pop ( 'pdf' , None ) , **kwargs )
        APDF2 .__init__ ( self )
        
        self.config   = { 'name'    : self.name    ,
                          'xvar'    : self.xvar    ,
                          'yvar'    : self.yvar    ,
                          'tricks'  : self.tricks  ,
                          'pdf'     : self.pdf     }
        self.config.update ( kwargs )
        
        self.__call_OK = isinstance ( self.xvar , ROOT.RooAbsRealLValue ) and \
                         isinstance ( self.yvar , ROOT.RooAbsRealLValue ) 
        

    # =========================================================================
    ## simple 'function-like' interface 
    def __call__ ( self , x , y , error = False , normalized = True ) :
        """ Simple  function-like interface
        >>>  pdf = ...
        >>>  print ( pdf(0.1,0.5) ) 
        """
        assert self.__call_OK , "Invalid types for xvar/yvar!"

        xmnmx = self.xminmax()
        if xmnmx :
            xmn , xmx = xmnmx 
            if not xmn <= x <= xmx : return 0

        ymnmx = self.yminmax()
        if ymnmx :
            ymn , ymx = ymnmx 
            if not ymn <= y <= ymx : return 0

        with SETVAR ( self.xvar ) , SETVAR( self.yvar ) :
            self.xvar.setVal ( x )
            self.yvar.setVal ( y )
            
            v = self.pdf.getVal ( self.vars ) if normalized else self.pdf.getValV ()
            
            if error and self.fit_result :
                e = self.pdf.getPropagatedError ( self.fit_result )
                if 0<= e : v = VE ( v ,  e * e )

            return v 
            
    # ========================================================================
    ## convert to float 
    def __float__ ( self ) :
        """Convert to float
        >>> fun = ...
        >>> v   = float ( fun )
        """
        return self.fun.getVal ( self.vars ) 


    # =========================================================================
    ## make a product of two PDFs
    #  @code
    #  pdf1 = ...
    #  pdf2 = ...
    #  pdf  = pdf1 * pdf2 
    #  @endcode
    #  Rules: 
    #  - PDF3 ( x , y , z ) * PDF3 ( x , y , z ) -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF2 ( x , y )     -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF2 ( x , z )     -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF2 ( y , z )     -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF1 ( x )         -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF1 ( y )         -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF1 ( z )         -> PDF3 ( x , y , z ) 
    #  - PDF2 ( x , y )     * PDF3 ( ... )       -> process as PDF3 (...) * PDF2 ( x , y )
    #  - PDF2 ( x , y )     * PDF2 ( x , y )     -> PDF2 ( x , y )
    #  - PDF2 ( x , y )     * PDF1 ( x )         -> PDF2 ( x , y )
    #  - PDF2 ( x , y )     * PDF1 ( y )         -> PDF2 ( x , y )
    #  - PDF1 ( x )         * PDF3 ( ... )       -> process as PDF3 (...) * PDF ( x )
    #  - PDF1 ( x )         * PDF2 ( ... )       -> process as PDF2 (...) * PDF ( x )
    #  - PDF1 ( x )         * PDF1 ( x )         -> PDF1 ( x )
    #  - PDF1 ( x )         * PDF1 ( y )         -> PDF2 ( x , y )
    def __mul__ ( self , other ) :
        """Make a product of two PDFs
        
        >>> pdf1 = ...
        >>> pdf2 = ...
        >>> pdf  = pdf1 * pdf2

        Rules:
        - PDF3 ( x , y , z ) * PDF3 ( x , y , z ) -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF2 ( x , y )     -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF2 ( x , z )     -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF2 ( y , z )     -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF1 ( x )         -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF1 ( y )         -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF1 ( z )         -> PDF3 ( x , y , z ) 
        - PDF2 ( x , y )     * PDF3 ( ... )       -> process as PDF3 (...) * PDF2 ( x , y )
        - PDF2 ( x , y )     * PDF2 ( x , y )     -> PDF2 ( x , y )
        - PDF2 ( x , y )     * PDF1 ( x )         -> PDF2 ( x , y )
        - PDF2 ( x , y )     * PDF1 ( y )         -> PDF2 ( x , y )
        - PDF1 ( x )         * PDF3 ( ... )       -> process as PDF3 (...) * PDF ( x )
        - PDF1 ( x )         * PDF2 ( ... )       -> process as PDF2 (...) * PDF ( x )
        - PDF1 ( x )         * PDF1 ( x )         -> PDF1 ( x )
        - PDF1 ( x )         * PDF1 ( y )         -> PDF2 ( x , y )
        """
        from ostap.fitting.pdf_ops import pdf2_product        
        return pdf2_product ( self , other )

    # =========================================================================
    ## make a product of two PDFs
    #  @code
    #  pdf1 = ...
    #  pdf2 = ...
    #  pdf  = pdf1 * pdf2 
    #  @endcode
    #  Rules 
    #  - PDF3 ( x , y , z ) * PDF3 ( x , y , z ) -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF2 ( x , y )     -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF2 ( x , z )     -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF2 ( y , z )     -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF1 ( x )         -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF1 ( y )         -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF1 ( z )         -> PDF3 ( x , y , z ) 
    #  - PDF2 ( x , y )     * PDF3 ( ... )       -> process as PDF3 (...) * PDF2 ( x , y )
    #  - PDF2 ( x , y )     * PDF2 ( x , y )     -> PDF2 ( x , y )
    #  - PDF2 ( x , y )     * PDF1 ( x )         -> PDF2 ( x , y )
    #  - PDF2 ( x , y )     * PDF1 ( y )         -> PDF2 ( x , y )
    #  - PDF1 ( x )         * PDF3 ( ... )       -> process as PDF3 (...) * PDF ( x )
    #  - PDF1 ( x )         * PDF2 ( ... )       -> process as PDF2 (...) * PDF ( x )
    #  - PDF1 ( x )         * PDF1 ( x )         -> PDF1 ( x )    
    #  - PDF1 ( x )         * PDF1 ( y )         -> PDF2 ( x , y )
    def __rmul__ ( self , other ) :
        """Make a product of two PDFs
        
        >>> pdf1 = ...
        >>> pdf2 = ...
        >>> pdf  = pdf1 * pdf2
        
        Rules:
        - PDF3 ( x , y , z ) * PDF3 ( x , y , z ) -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF2 ( x , y )     -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF2 ( x , z )     -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF2 ( y , z )     -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF1 ( x )         -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF1 ( y )         -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF1 ( z )         -> PDF3 ( x , y , z ) 
        - PDF2 ( x , y )     * PDF3 ( ... )       -> process as PDF3 (...) * PDF2 ( x , y )
        - PDF2 ( x , y )     * PDF2 ( x , y )     -> PDF2 ( x , y )
        - PDF2 ( x , y )     * PDF1 ( x )         -> PDF2 ( x , y )
        - PDF2 ( x , y )     * PDF1 ( y )         -> PDF2 ( x , y )
        - PDF1 ( x )         * PDF3 ( ... )       -> process as PDF3 (...) * PDF ( x )
        - PDF1 ( x )         * PDF2 ( ... )       -> process as PDF2 (...) * PDF ( x )
        - PDF1 ( x )         * PDF1 ( x )         -> PDF1 ( x )
        - PDF1 ( x )         * PDF1 ( y )         -> PDF2 ( x , y )
        """
        from ostap.fitting.pdf_ops import pdf2_product, pdf3_product
        if   isinstance ( other , PDF3 ) : return pdf3_product ( other , self  )
        elif isinstance ( other , PDF1 ) : return pdf2_product ( self  , other )
        return NotImplemented 

    # =========================================================================
    ## make a non-extender sum of 1D  PDFs
    #  @code
    #  pdf1 = ...
    #  pdf2 = ...
    #  pdf  = pdf1 + pdf2     
    #  @endcode
    def __add__ ( self , other ) :
        """Make a no-extended sum of 1D PDFs
        >>> pdf1 = ...
        >>> pdf2 = ...
        >>> pdf  = pdf1 + pdf2     
        """
        from ostap.fitting.pdf_ops import pdf2_sum        
        return pdf2_sum ( self , other )

    # =========================================================================
    ## make a non-extended sum of 1D  PDFs
    #  @code
    #  pdf1 = ...
    #  pdf2 = ...
    #  pdf  = pdf1 + pdf2     
    #  @endcode
    def __radd__ ( self , other ) :
        """Make a no-extended sum of 1D PDFs
        >>> pdf1 = ...
        >>> pdf2 = ...
        >>> pdf  = pdf1 + pdf2     
        """
        from ostap.fitting.pdf_ops import pdf2_sum        
        return pdf2_sum ( self , other )

    # =========================================================================
    ## Convert PDF into simple function
    #  @code
    #  pdf = ...
    #  fun = pdf.as_FUN () 
    #  @endcode
    def as_FUN ( self , name = '' ) : 
        """Convert PDF into simple function
        >>> pdf = ...
        >>> fun = pdf.as_FUN () 
        """
        return Fun2D ( self.pdf  ,
                       xvar = self.xvar  ,
                       yvar = self.yvar  ,
                       name = name if name else self.new_name ( 'fun2' ) ) 
    

# =============================================================================
## @class Generic2D_pdf
#  "Wrapper" over generic RooFit (2D)-pdf
#  @code
#  raw_pdf = 
#  pdf     = Generic2D_pdf ( raw_pdf , xvar = ... , yvar = ... )  
#  @endcode 
#  If more functionality is required , more actions are possible:
#  @code
#  ## for sPlot 
#  pdf.alist2 = ROOT.RooArgList ( n1 , n2 , n3 ) ## for sPlotting 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-03-29
class Generic2D_pdf(PDF2) :
    """ Wrapper for generic (2D) RooFit pdf    
    >>> raw_pdf = 
    >>> pdf     = Generic2D_pdf ( raw_pdf , xvar = ... , yvar = ... )
    """
    ## constructor 
    def __init__ ( self , pdf , xvar , yvar ,
                   name           = None    ,
                   add_to_signals = True    ,
                   prefix         = ''      ,
                   suffix         = ''      ) :

        assert isinstance ( xvar , ROOT.RooAbsReal ) , "'xvar' must be ROOT.RooAbsReal"
        assert isinstance ( yvar , ROOT.RooAbsReal ) , "'yvar' must be ROOT.RooAbsReal"        
        assert isinstance ( pdf  , ROOT.RooAbsReal ) , "'pdf' must be ROOT.RooAbsReal"
        
        name = name if name else self.generate_name ( prefix = prefix + '%s_' % pdf.GetName() , suffix = suffix ) 
        PDF2  . __init__ ( self , name , xvar , yvar )

        ## PDF! 
        self.pdf = pdf

        if not self.xvar in self.params () : 
            self.warning ( "Function/PDF does not depend on xvar=%s" % self.xvar.name )
        if not self.yvar in self.params () : 
            self.warning ( "Function/PDF does not depend on yvar=%s" % self.yvar.name )

        ## add it to the list of signal components ?
        self.__add_to_signals = True if add_to_signals else False
        
        if self.add_to_signals :
            self.signals.add ( self.pdf )
       
        ## save the configuration
        self.config = {
            'pdf'            : self.pdf            ,
            'xvar'           : self.xvar           ,
            'yvar'           : self.yvar           ,
            'name'           : self.name           ,
            'add_to_signals' : self.add_to_signals , 
            'prefix'         : prefix              ,
            'suffix'         : suffix              ,            
            }

        self.checked_keys.add ( 'pdf'     )
        self.checked_keys.add ( 'xvar'    )
        self.checked_keys.add ( 'yvar'    )
        
    @property
    def add_to_signals ( self ) :
        """'add_to_signals' : should PDF be added into list of signal components?"""
        return self.__add_to_signals 

# =============================================================================
## 3D Stuff goes here 
# =============================================================================

# =============================================================================
# @class APDF3
# The helper MIXIN class for implementation of 3D-pdfs 
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date 2017-11-11
class APDF3 (APDF2) :
    """ Useful helper MIXIN class for implementation of PDFs for 3D-fit
    """
    def __init__ ( self ) : 
        
        APDF2.__init__ ( self )
        
    # =========================================================================
    ## make the actual fit 
    #  @code
    #  r,f = model.fitTo ( dataset )
    #  r,f = model.fitTo ( dataset , weighted = True )    
    #  r,f = model.fitTo ( dataset , ncpu     = 10   )    
    #  r,f = model.fitTo ( dataset )    
    #  @endcode 
    def fitTo ( self           , 
                dataset        ,
                silent = False ,
                refit  = False ,
                timer  = False ,
                draw   = False , 
                args   = ()    , **kwargs ) :
        """
        Perform the actual fit (and draw it)
        >>> r,f = model.fitTo ( dataset )
        >>> r,f = model.fitTo ( dataset , weighted = True )    
        >>> r,f = model.fitTo ( dataset , ncpu     = 10   )    
        >>> r,f = model.fitTo ( dataset )    
        """
        if   isinstance ( dataset , H3D_dset ) : dataset = dataset.dset        
        elif isinstance ( dataset , ROOT.TH3 ) :
            density = kwargs.pop ( 'density' , False ) 
            chi2    = kwargs.pop ( 'chi2'    , False ) 
            return self.fitHisto ( dataset   ,
                                   draw    = draw    ,
                                   silent  = silent  ,
                                   density = dentity ,
                                   chi2    = chi2    , args = args , **kwargs )

        ## play a bit with binning cache for convolutions 
        if self.zvar.hasBinning ( 'cache' ) :
            nb1 = self.zvar.getBins( 'cache' ) 
            zv  = getattr ( dataset , self.zvar.name , None )
            if   zv and zv.hasBinning ( 'cache' ) :
                nb2 = zv.getBins('cache')
                if  nb1 != nb2 :
                    zv.setBins ( max (  nb1 , nb2 ) , 'cache' )
                    self.info ('Adjust binning cache %s->%s for variable %s in dataset' % ( nb2 , nb1 , zv.name ) )
            elif zv :
                zv.setBins (        nb1         , 'cache' )
                self    .info ('Set binning cache %s for variable %s in dataset' %  ( nb1 , zv.name )  )
                                
        
        result , f2 = APDF2.fitTo ( self    ,
                                    dataset = dataset ,
                                    draw    = False   , ## False here!
                                    nbins   = 50      , ## fake  here!
                                    ybins   = 20      , ## fake  here!
                                    silent  = silent  ,
                                    refit   = refit   ,
                                    timer   = timer   ,
                                    args    = args    , **kwargs )

        if   draw and draw in ( 1 , '1' , 'x' , 'X' , self.xvar.name ) :
            f = self.draw1 ( daatset , silent = silent , args = args , **kwargs )
            return result, f
        elif draw and draw in ( 2 , '2' , 'y' , 'Y' , self.yvar.name ) :
            f = self.draw2 ( daatset , silent = silent , args = args , **kwargs )
            return result, f
        elif draw and draw in ( 3 , '3' , 'z' , 'Z' , self.zvar.name ) :
            f = self.draw3 ( daatset , silent = silent , args = args , **kwargs )
            return result, f

        return result
    
    # =========================================================================
    ## draw the projection over 1st variable
    #
    #  @code
    #  r,f = model.fitTo ( dataset ) ## fit dataset
    #  fx  = model.draw1 ( dataset , nbins = 100 ) ## draw results
    #
    #  fx  = model.draw1 ( dataset , nbins = 100 , in_range2 = (2,3) ) ## draw results
    #
    #  model.yvar.setRange ( 'QUQU2' , 2 , 3 ) 
    #  fx  = model.draw1 ( dataset , nbins = 100 , in_range2 = 'QUQU2') ## draw results
    #
    #  @endcode 
    def draw1 ( self             ,
                dataset   = None ,
                nbins     = 100  ,
                silent    = True ,
                in_range2 = None ,
                in_range3 = None ,
                args      = ()   , **kwargs ) :
        """ Draw the projection over 3rd variable
        
        >>> r,f = model.fitTo ( dataset ) ## fit dataset
        >>> fx  = model.draw1 ( dataset , nbins = 100 ) ## draw results
        
        >>> fx  = model.draw1 ( dataset , nbins = 100 , in_range2 = (2,3) ) ## draw results

        >>> model.yvar.setRange ( 'QUQU2' , 2 , 3 ) 
        >>> fx  = model.draw1 ( dataset , nbins = 100 , in_range2 = 'QUQU2') ## draw results
        
        """
        in_range = None
        if in_range2        and  not in_range3 :
            in_range = 'aux3_rng12_%s'  % self.name
        elif not in_range2  and      in_range3 :
            in_range = 'aux3_rng13_%s'  % self.name
        elif in_range2  and in_range3:
            if  not isinstance ( in_range2 , tuple ):
                in_range = in_range2
            elif not isinstance ( in_range3 , tuple ):
                in_range = in_range3
            else:
                in_range = 'aux3_rng123_%s' % self.name


        if in_range2 and isinstance ( in_range2 , tuple ) and 2 == len ( in_range2 ) :
            with roo_silent ( silent , 3 ) :
                self.yvar.setRange ( in_range, in_range2[0] , in_range2[1] )
                if dataset:
                    dataset.get_var(self.yvar.GetName()).setRange ( in_range, in_range2[0] , in_range2[1] )
#            in_range2  = range_name 

        if in_range3 and isinstance ( in_range3 , tuple ) and 2 == len ( in_range3 ) :
            with roo_silent ( silent , 3 ) : 
                self.zvar.setRange ( in_range , in_range3[0] , in_range3[1] )
                if dataset:
                    dataset.get_var(self.zvar.GetName()).setRange ( in_range, in_range3[0] , in_range3[1] )
                #    in_range3  = range_name 
                
                
        #  if in_range2 and in_range3: 
        #     in_range= in_range3 
        # elif in_range2:
        #     in_range=  in_range2 
        # elif in_range3:
        #     in_range= in_range3 
        
        return self.draw ( drawvar  = self.xvar , 
                           dataset  = dataset   ,
                           nbins    = nbins     ,
                           ybins    = 20        , ## fake 
                           silent   = silent    ,
                           in_range = in_range  ,
                           args     = args      , **kwargs )


    # =========================================================================
    ## draw the projection over 2nd variable
    #
    #  @code
    #  r,f = model.fitTo ( dataset ) ## fit dataset
    #  fy  = model.draw1 ( dataset , nbins = 100 ) ## draw results
    #
    #  fy  = model.draw1 ( dataset , nbins = 100 , in_range1 = (2,3) ) ## draw results
    #
    #  model.xvar.setRange ( 'QUQU1' , 2 , 3 ) 
    #  fy  = model.draw1 ( dataset , nbins = 100 , in_range1 = 'QUQU1') ## draw results
    #
    #  @endcode 
    def draw2 ( self            ,
                dataset   = None ,
                nbins     = 100  ,
                silent    = True ,
                in_range1 = None ,
                in_range3 = None ,
                args      = ()   , **kwargs ) :
        """ Draw the projection over 2nd variable
        
        >>> r,f = model.fitTo ( dataset ) ## fit dataset
        >>> fy  = model.draw2 ( dataset , nbins = 100 ) ## draw results
        
        >>> fx  = model.draw2 ( dataset , nbins = 100 , in_range1 = (2,3) ) ## draw results

        >>> model.xvar.setRange ( 'QUQU1' , 2 , 3 ) 
        >>> fx  = model.draw2 ( dataset , nbins = 100 , in_range1 = 'QUQU1') ## draw results
        
        """
        in_range=None

        if in_range1  and  not in_range3:
            in_range = 'aux3_rng21_%s' % self.name
        elif not in_range1  and in_range3:
            in_range = 'aux3_rng23_%s' % self.name
        elif in_range1  and in_range3:
            if  not isinstance ( in_range1 , tuple ):
                in_range = in_range1
            elif  not isinstance ( in_range3 , tuple ):
                in_range = in_range3
            else:
                in_range = 'aux3_rng213_%s' % self.name

        if in_range1 and isinstance ( in_range1 , tuple ) and 2 == len ( in_range1 ) :
            with roo_silent ( silent , 3 ) : 
                self.xvar.setRange ( in_range , in_range1[0] , in_range1[1] )  
                if dataset:
                    dataset.get_var(self.xvar.GetName()).setRange ( in_range , in_range1[0] , in_range1[1] )


        if in_range3 and isinstance ( in_range3 , tuple ) and 2 == len ( in_range3 ) :
            with roo_silent ( silent , 3 ) : 
                self.zvar.setRange ( in_range, in_range3[0] , in_range3[1] )
                if dataset:
                    dataset.get_var(self.zvar.GetName()).setRange (in_range , in_range3[0] , in_range3[1] )

        
        return self.draw ( drawvar  = self.yvar , 
                           dataset  = dataset   ,
                           nbins    = nbins     ,
                           ybins    = 20        , ## fake 
                           silent   = silent    ,
                           in_range = in_range  ,
                           args     = args      , **kwargs )


    # =========================================================================
    ## draw the projection over 3rd variable
    #
    #  @code
    #  r,f = model.fitTo ( dataset ) ## fit dataset
    #  fz  = model.draw3 ( dataset , nbins = 100 ) ## draw results
    #
    #  fz  = model.draw3 ( dataset , nbins = 100 , in_range2 = (2,3) ) ## draw results
    #
    #  model.yvar.setRange ( 'QUQU2' , 2 , 3 ) 
    #  f  = model.draw3 ( dataset , nbins = 100 , in_range2 = 'QUQU2') ## draw results
    #  @endcode 
    def draw3 ( self             ,
                dataset   = None ,
                nbins     = 100  ,
                silent    = True ,
                in_range1 = None ,
                in_range2 = None ,
                args      = ()   ,  **kwargs ) :
        """ Draw the projection over 3rd variable
        
        >>> r,f = model.fitTo ( dataset ) ## fit dataset
        >>> fx  = model.draw3 ( dataset , nbins = 100 ) ## draw results
        
        >>> fx  = model.draw3 ( dataset , nbins = 100 , in_range2 = (2,3) ) ## draw results

        >>> model.yvar.setRange ( 'QUQU2' , 2 , 3 ) 
        >>> fx  = model.draw3 ( dataset , nbins = 100 , in_range2 = 'QUQU2') ## draw results
        
        """
        in_range=None

        if in_range1  and  not in_range2:
            in_range = 'aux3_rng31_%s' % self.name
        elif not in_range1  and in_range2:
            in_range = 'aux3_rng32_%s' % self.name
        elif in_range1  and in_range2:
            if  not isinstance ( in_range1 , tuple ):
                in_range = in_range1
            elif  not isinstance ( in_range2 , tuple ):
                in_range = in_range2
            else:
                in_range = 'aux3_rng312_%s' % self.name

        if in_range1 and isinstance ( in_range1 , tuple ) and 2 == len ( in_range1 ) :
            with roo_silent ( silent , 3 ) : 
                self.xvar.setRange ( in_range ,  in_range1[0] , in_range1[1] )       
                if dataset:
                    dataset.get_var(self.xvar.GetName()).setRange ( in_range , in_range1[0] , in_range1[1] )

        if in_range2 and isinstance ( in_range2 , tuple ) and 2 == len ( in_range2 ) :
            with roo_silent ( silent , 3 ) : 
                self.yvar.setRange ( in_range , in_range2[0] , in_range2[1] )    
                if dataset:
                    dataset.get_var(self.yvar.GetName()).setRange ( in_range , in_range2[0] , in_range2[1] )

        
        return self.draw ( drawvar  = self.zvar , 
                           dataset  = dataset   ,
                           nbins    = nbins     ,
                           ybins    = 20        , ## fake 
                           silent   = silent    ,
                           in_range = in_range  ,
                           args     = args      , **kwargs )
    
    
    # =========================================================================
    ## make 1D-plot
    def draw ( self            ,
               drawvar  = None ,
               dataset  = None ,
               nbins    =  100 ,
               silent   = True ,
               in_range = None ,
               args     = ()   , 
               **kwargs        ) : 
        """
        Make 1D-plot:
        """
        if drawvar in ( 'z'  , 'Z' , '3' , 3 , self.zvar.name ) :
            drawvar = self.zvar
            
        return APDF2.draw  ( self ,
                             drawvar  = drawvar  ,
                             dataset  = dataset  ,
                             nbins    = nbins    ,
                             silent   =  silent  ,
                             in_range = in_range ,
                             args     = args     , **kwargs )
    
    # =========================================================================
    ## fit the 3D-histogram (and draw it)
    #
    #  @code
    #
    #  histo = ...
    #  r,f = model.fitHisto ( histo )
    #
    #  @endcode
    def fitHisto ( self            ,
                   histo           ,
                   draw    = False ,
                   silent  = False ,
                   density = False ,
                   chi2    = False , 
                   args    = ()    , **kwargs ) :
        """Fit the histogram (and draw it)
        
        >>> histo = ...
        >>> r,f = model.fitHisto ( histo , draw = True )
        
        """
        
        xminmax = histo.xminmax()
        yminmax = histo.yminmax()
        zminmax = histo.zminmax()
        
        with     RangeVar ( self.xvar , *xminmax ) , \
                 RangeVar ( self.yvar , *yminmax ) , \
                 RangeVar ( self.xvar , *zminmax ): 

            hdata = getattr ( self , 'histo_data' , None )
            if hdata and isinstance ( hdata  , H3D_dset ) and \
                   hdata.histo      is histo              and \
                   hdata.density    == density            and \
                   hdata.histo_hash == hash ( histo ) :
                ## reuse the existing dataset
                self.debug ('Reuse the existing H3D_dset') 
                data = hdata.dset
            else :
                ## convert it! 
                self.debug ('Create new H3D_dset'        ) 
                self.histo_data = H3D_dset ( histo , self.xvar , self.yvar  , self.zvar ,
                                             density , silent )
                data = self.histo_data
                
            if chi2 : return self.chi2fitTo ( data              ,
                                              draw    = draw    ,
                                              silent  = False   ,
                                              density = density ,
                                              args    = args    , **kwargs )
            else    : return self.fitTo     ( data              ,
                                              silent  = silent  ,
                                              args    =  args   , **kwargs ) 

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
                   binning  = {}    ,
                   sample   = True  , 
                   args     = ()    ) :
        """Generate toy-sample according to PDF
        >>> model  = ....
        >>> data   = model.generate ( 10000 ) ## generate dataset
        
        >>> varset = ....
        >>> data   = model.generate ( 100000 , varset , sample = False )
        >>> data   = model.generate ( 100000 , varset , sample = True  )
        """
        nEvents = self.gen_sample ( nEvents ) if sample else nEvents 
        assert 0 <= nEvents , 'Invalid number of Events %s' % nEvents  

        args = args + ( ROOT.RooFit.Name ( dsID() ) , ROOT.RooFit.NumEvents ( nEvents ) )

        if binning is True :
            args    = args + ( ROOT.AllBinned() , ) 
            binning = {}

        if   not varset :
            varset = ROOT.RooArgSet( self.xvar , self.yvar , self.zvar )
        elif isinstance ( varset , ROOT.RooAbsReal ) :
            varset = ROOT.RooArgSet( varset )

        if not self.xvar in varset :
            vs = ROOT.RooArgSet()
            vs . add ( self.xvar )
            for  v in varset : vs.add ( v )
            varset = vs

        if not self.yvar in varset :
            vs = ROOT.RooArgSet()
            vs . add ( self.yvar )
            for  v in varset : vs.add ( v )
            varset = vs

        if not self.zvar in varset :
            vs = ROOT.RooArgSet()
            vs . add ( self.zvar )
            for  v in varset : vs.add ( v )
            varset = vs
            
        from ostap.fitting.variables import KeepBinning        
        with KeepBinning ( self.xvar ) , KeepBinning ( self.yvar ), KeepBinning ( self.zvar ) : 
            
            if binning :
                
                xbins = binning.get ( self.xvar.name , None )
                ybins = binning.get ( self.yvar.name , None )
                zbins = binning.get ( self.zvar.name , None )
                
                if xbins : self.xvar.bins = xbins
                if ybins : self.yvar.bins = ybins
                if zbins : self.zvar.bins = zbins
   
        return self.pdf.generate ( varset , *args )


    # ========================================================================
    ## check minmax of the PDF using the random shoots
    #  @code
    #  pdf     = ....
    #  mn , mx = pdf.minmax()            
    #  @endcode 
    def minmax ( self , nshoots = 200000 ) :
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

        ## check RooAbsReal functionality
        code = self.pdf.getMaxVal( ROOT.RooArgSet ( self.xvar , self.yvar , self.zvar ) )
        if 0 < code :
            mx = self.pdf.maxVal ( code )
            if 0 < mx : return 0 , mx
            
        ## not try  to use random
                
        mn , mx = -1 , -10
        if hasattr ( self.pdf , 'min' ) : mn = self.pdf.min()
        if hasattr ( self.pdf , 'max' ) : mx = self.pdf.max()
        if 0 <= mn and mn <= mx and 0 < mx : return mn , mx
        
        if not self.xminmax() : return ()
        if not self.yminmax() : return ()
        if not self.zminmax() : return ()
        
        mn  , mx = -1 , -10
        xmn , xmx = self.xminmax()
        ymn , ymx = self.yminmax()
        zmn , zmx = self.zminmax()
        for i in range ( nshoots ) : 
            xx = random.uniform ( xmn , xmx )
            yy = random.uniform ( ymn , ymx )
            zz = random.uniform ( zmn , zmx )
            with SETVAR ( self.xvar ) , SETVAR ( self.yvar ) , SETVAR ( self.zvar ) :
                self.xvar.setVal ( xx )
                self.yvar.setVal ( yy )
                self.zvar.setVal ( zz )
                vv = self.pdf.getVal()
                if mn < 0 or vv < mn : mn = vv
                if mx < 0 or vv > mx : mx = vv
                        
        return mn , mx 

    # =========================================================================
    ## get integral over (xmin,xmax,ymin,ymax,zmin,zmax) region
    #  @code
    #  pdf = ...
    #  print ( pdf.integral( 0,1,0,2,0,5) ) 
    #  @endcode
    def integral ( self, xmin , xmax , ymin , ymax , zmin , zmax , nevents = True ) :
        """Get integral over (xmin,xmax,ymin,ymax,zmin,zmax) region
        >>> pdf = ...
        >>> print ( pdf.integral( 0,1,0,2,0,5) ) 
        """
        if self.xminmax() :            
            xmn , xmx = self.xminmax()
            xmin = max ( xmin , xmn )
            xmax = min ( xmax , xmx )
            
        if self.yminmax() : 
            ymn , ymx = self.yminmax()
            ymin = max ( ymin , ymn )
            ymax = min ( ymax , ymx )
            
        if self.zminmax() :
            zmn , zmx = self.zminmax()
            zmin = max ( zmin , zmn )
            zmax = min ( zmax , zmx )

        
        value , todo = 0 , True 
        
         ## 1) make a try to use analytical integral (could be fast)
        if self.tricks :
            try:
                if hasattr ( self.pdf , 'setPars'  ) : self.pdf.setPars() 
                fun          = self.pdf.function()
                value , todo = fun.integral ( xmin , xmax ,
                                              ymin , ymax ,
                                              zmin , zmax ) , False 
            except:
                pass


        ## for numerical integration 
        from ostap.math.integral import integral3 as _integral3
        
        extended = self.pdf.canBeExtended() or isinstance ( self.pdf , ROOT.RooAddPdf )
        if todo  and extended  :
            value = _integral3 ( self , xmin , xmax , ymin , ymax , zmin , zmax )
        elif todo : 
                        
            ## use unormalized PDF here to speed up the integration 
            ifun   = lambda x , y , z : self ( x , y , z , error = False , normalized = False )
            value  = _integral3 ( ifun , xmin , xmax , ymin , ymax , zmin , zmax )
            norm   = self.pdf.getNorm ( self.vars )
            value /= norm
            
        if nevents and self.pdf.mustBeExtended () :
            evts = self.pdf.expectedEvents( self.vars )
            if evts  <= 0 or iszero ( evts ) :
                self.warning ( "integral: expectedEvents is %s" % evts )
            value *= evts 
                
        return value
    
    # ==========================================================================
    ## get a minimum of PDF for certain interval
    #  @code
    #  pdf2 = ...
    #  x , y , z = pdf3.minimum() 
    #  @endcode 
    def minimum ( self ,
                  xmin = None , xmax = None ,
                  ymin = None , ymax = None ,
                  zmin = None , zmax = None , x0 = () ) :
        """Get a minimum of PDF for certain interval
        >>> pdf3 = ...
        >>> x, y , z = pdf3.minimum()
        """
        
        if xmin is None : xmin = self.xminmax()[0]
        if xmax is None : xmax = self.xminmax()[1]
        if self.xminmax() :
            xmin =  max ( xmin , self.xminmax()[0] )
            xmax =  min ( xmax , self.xminmax()[1] )

        if ymin is None : ymin = self.yminmax()[0]
        if ymax is None : ymax = self.yminmax()[1]
        if self.yminmax() :
            ymin =  max ( ymin , self.yminmax()[0] )
            ymax =  min ( ymax , self.yminmax()[1] )

        if zmin is None : zmin = self.zminmax()[0]
        if zmax is None : zmax = self.zminmax()[1]
        if self.zminmax() :
            zmin =  max ( zmin , self.zminmax()[0] )
            zmax =  min ( zmax , self.zminmax()[1] )
            
        if not x0 : x0 = 0.5 * ( xmin + xmax ) , 0.5 * ( ymin + ymax ) , 0.5 * ( zmin + zmax )

        if not xmin <= x0[0] <= xmax :
            self.error("Wrong xmin/x0[0]/xmax: %s/%s/%s"   % ( xmin , x0[0] , xmax ) )

        if not ymin <= x0[1] <= ymax : 
            self.error("Wrong ymin/x0[1]/ymax: %s/%s/%s"   % ( ymin , x0[1] , ymax ) )
        
        if not zmin <= x0[2] <= zmax : 
            self.error("Wrong zmin/x0[2]/zmax: %s/%s/%s"   % ( zmin , x0[2] , zmax ) )

        from ostap.math.minimize import sp_minimum_3D
        return sp_minimum_3D (  self ,
                                xmin , xmax ,
                                ymin , ymax ,
                                zmin , zmax , x0 )
    

    # ==========================================================================
    ## get a maximum of PDF for certain interval
    #  @code
    #  pdf2 = ...
    #  x , y , z = pdf3.maximum() 
    #  @endcode 
    def minimum ( self ,
                  xmin = None , xmax = None ,
                  ymin = None , ymax = None ,
                  zmin = None , zmax = None , x0 = () ) :
        """Get a maximum of PDF for certain interval
        >>> pdf3 = ...
        >>> x, y , z = pdf3.maximum()
        """
        
        if xmin is None : xmin = self.xminmax()[0]
        if xmax is None : xmax = self.xminmax()[1]
        if self.xminmax() :
            xmin =  max ( xmin , self.xminmax()[0] )
            xmax =  min ( xmax , self.xminmax()[1] )

        if ymin is None : ymin = self.yminmax()[0]
        if ymax is None : ymax = self.yminmax()[1]
        if self.yminmax() :
            ymin =  max ( ymin , self.yminmax()[0] )
            ymax =  min ( ymax , self.yminmax()[1] )

        if zmin is None : zmin = self.zminmax()[0]
        if zmax is None : zmax = self.zminmax()[1]
        if self.zminmax() :
            zmin =  max ( zmin , self.zminmax()[0] )
            zmax =  min ( zmax , self.zminmax()[1] )
            
        if not x0 : x0 = 0.5 * ( xmin + xmax ) , 0.5 * ( ymin + ymax ) , 0.5 * ( zmin + zmax )

        if not xmin <= x0[0] <= xmax :
            self.error("Wrong xmin/x0[0]/xmax: %s/%s/%s"   % ( xmin , x0[0] , xmax ) )

        if not ymin <= x0[1] <= ymax : 
            self.error("Wrong ymin/x0[1]/ymax: %s/%s/%s"   % ( ymin , x0[1] , ymax ) )
        
        if not zmin <= x0[2] <= zmax : 
            self.error("Wrong zmin/x0[2]/zmax: %s/%s/%s"   % ( zmin , x0[2] , zmax ) )

        from ostap.math.minimize import sp_maximum_3D
        return sp_maximum_3D ( self ,
                               xmin , xmax ,
                               ymin , ymax ,
                               zmin , zmax , x0 )
    

    # ==========================================================================
    ## convert PDF into TF2 object, e.g. to profit from TF3::Draw options
    #  @code
    #  pdf = ...
    #  tf3 = pdf.tf()
    #  tf3.Draw( options )
    #  @endcode
    def tf ( self                      ,
             xmin = None , xmax = None ,
             ymin = None , ymax = None ,
             zmin = None , zmax = None ) :
        """Convert PDF  to TF3 object, e.g. to profit from TF3::Draw options
        >>> pdf = ...
        >>> tf3 = pdf.tf()
        >>> tf3.Draw('colz')
        """
        def _aux_fun_ ( x , pars = [] ) :
            return self ( x[0] , x[1] , x[2] , error = False )

        if xmin == None and self.xminmax() : xmin = self.xminmax()[0]
        if xmax == None and self.xminmax() : xmax = self.xminmax()[1]
        if ymin == None and self.yminmax() : ymin = self.yminmax()[0]
        if ymax == None and self.yminmax() : ymax = self.yminmax()[1]
        if zmin == None and self.zminmax() : zmin = self.zminmax()[0]
        if zmax == None and self.zminmax() : zmax = self.zminmax()[1]
        
        if xmin == None : xmin = 0.0
        if xmax == None : xmin = 1.0
        if ymin == None : ymin = 0.0
        if ymax == None : ymin = 1.0
        if zmin == None : zmin = 0.0
        if zmax == None : zmin = 1.0
        
        from ostap.core.core import fID
        return ROOT.TF3 ( fID() , _aux_fun_ , xmin , xmax , ymin , ymax , zmin , zmax ) 


    # ==========================================================================
    ## create the histogram accoring to specifications 
    def make_histo ( self ,
                     xbins    = 10    , xmin = None , xmax = None ,
                     ybins    = 10    , ymin = None , ymax = None ,
                     zbins    = 10    , zmin = None , zmax = None ,
                     hpars    = ()    , 
                     histo    = None  ) :
        """Create the histogram accoring to specifications"""
        
        import ostap.histos.histos

        # histogram is provided 
        if histo :
            
            assert isinstance ( histo  , ROOT.TH3 ), \
                   "Illegal type of 'histo'-argument %s" % type( histo )
            
            histo = histo.clone()
            histo.Reset()

        # arguments for the histogram constructor 
        elif hpars :
            
            from ostap.core.core import hID
            histo = ROOT.TH3F ( hID () , 'PDF%s' % self.name , *hpars  )
            if not histo.GetSumw2() : histo.Sumw2()

        # explicit contruction from (#bins,min,max)-triplet  
        else :
            
            assert isinstance ( xbins , integer_types ) and 0 < xbins, \
                   "Wrong 'xbins'-argument %s" % xbins 
            assert isinstance ( ybins , integer_types ) and 0 < ybins, \
                   "Wrong 'ybins'-argument %s" % ybins 
            assert isinstance ( zbins , integer_types ) and 0 < zbins, \
                   "Wrong 'zbins'-argument %s" % zbins 
            if xmin == None and self.xminmax() : xmin = self.xminmax()[0]
            if xmax == None and self.xminmax() : xmax = self.xminmax()[1]
            if ymin == None and self.yminmax() : ymin = self.yminmax()[0]
            if ymax == None and self.yminmax() : ymax = self.yminmax()[1]
            if zmin == None and self.zminmax() : zmin = self.zminmax()[0]
            if zmax == None and self.zminmax() : zmax = self.zminmax()[1]
            
            from ostap.core.core import hID
            histo = ROOT.TH3F ( hID() , 'PDF%s' % self.name ,
                                xbins , xmin , xmax ,
                                ybins , ymin , ymax ,
                                zbins , zmin , zmax )
            if not histo.GetSumw2() : histo.Sumw2()

        return histo 


    # ==========================================================================
    ## Convert PDF to the 3D-histogram in correct  way 
    #  @code
    #  pdf = ...
    #  h1  = pdf.histo ( 10 , 0. , 10. , 10 , 0. , 4. , 10 , 0. , 3 ) ## specify histogram parameters
    #  histo_template = ...
    #  h2  = pdf.histo ( histo = histo_template ) ## use histogram template
    #  h3  = pdf.histo ( ... , integral = True  ) ## use PDF integral within the bin  
    #  @endcode
    def histo ( self             ,
                xbins    = 10    , xmin = None , xmax = None ,
                ybins    = 10    , ymin = None , ymax = None ,
                zbins    = 10    , zmin = None , zmax = None ,
                hpars    = ()    , 
                histo    = None  ,
                intergal = True  ,
                errors   = False ) :
        """Convert PDF to the 3D-histogram in correct way
        >>> pdf = ...
        >>> h1  = pdf.histo ( 10 , 0. , 10. , 10 , 0. , 4. , 10 , 0. , 3 ) ## specify histogram parameters
        >>> histo_template = ...
        >>> h2  = pdf.histo ( histo = histo_template ) ## use histogram template
        >>> h3  = pdf.histo ( ... , integral = True  ) ## use PDF integral within the bin  
        """


        histo = self.make_histo ( xbins = xbins , xmin = xmin , xmax = xmax ,
                                  ybins = ybins , ymin = ymin , ymax = ymax ,
                                  zbins = zbins , zmin = zmin , zmax = zmax ,
                                  hpars = hpars ,
                                  histo = histo )
        
        
        # loop over the histogram bins 
        for ix , iy , iz , x , y , z , w in histo.items() :

            xv , xe = x.value() , x.error()
            yv , ye = y.value() , y.error()
            zv , ze = z.value() , z.error()
            
            # value at the bin center 
            c = self ( xv , yv , zv , error = errors ) 

            if not integral : 
                histo[ix,iy,iz] = c
                continue

            # integral over the bin 
            v  = self.integral( xv - xe , xv + xe ,
                                yv - ye , yv + ye ,
                                zv - ze , zv + ze )
            
            if errors :
                if    0 == c.cov2 () : pass
                elif  0 != c.value() and 0 != v : 
                    v = c * ( v / c.value() )
                    
            histo[ix,iy,iz] = v 

        return histo
    
    # ==========================================================================
    ## Convert PDF to the 3D-histogram, taking PDF-values at bin-centres
    #  @code
    #  pdf = ...
    #  h1  = pdf.roo_histo ( 10 , 0. , 10. , 10 , 0. , 4. , 10 , 0. , 3 ) 
    #  histo_template = ...
    #  h2  = pdf.roo_histo ( histo = histo_template ) ## use histogram template
    #  @endcode
    def roo_histo ( self            ,
                   xbins    = 10    , xmin = None , xmax = None ,
                   ybins    = 10    , ymin = None , ymax = None ,
                   zbins    = 10    , zmin = None , zmax = None ,
                   hpars    = ()    , 
                   histo    = None  ,
                   events   = True) :
        """Convert PDF to the 3D-histogram, taking PDF-values at bin-centres
        >>> pdf = ...
        >>> h1  = pdf.roo_histo ( 10 , 0. , 10. , 10 , 0. , 4. , 10 , 0. , 3 )
        >>> histo_template = ...
        >>> h2  = pdf.roo_histo ( histo = histo_template ) ## use histogram template
        """
        
        histo = self.make_histo ( xbins = xbins , xmin = xmin , xmax = xmax ,
                                  ybins = ybins , ymin = ymin , ymax = ymax ,
                                  zbins = zbins , zmin = zmin , zmax = zmax ,
                                  hpars = hpars ,
                                  histo = histo )

        hh = self.pdf.createHistogram (
            hID()     ,
            self.xvar ,                    self.binning ( histo.GetXaxis() , 'histo3x' )   ,
            ROOT.RooFit.YVar ( self.yvar , self.binning ( histo.GetYaxis() , 'histo3y' ) ) , 
            ROOT.RooFit.ZVar ( self.zvar , self.binning ( histo.GetZaxis() , 'histo3z' ) ) , 
            ROOT.RooFit.Scaling  ( False ) , 
            ROOT.RooFit.Extended ( False ) )
        
        for i in hh : hh.SetBinError ( i , 0 ) 
        
        if events and self.pdf.mustBeExtended() :            
            for ix , iy , iz , x , y , z , v  in hh.items() :
                volume               = 8 * x.error()  * y.error() * z.error() 
                hh [ iz , iy , iz ] *= volume
                
            hh *= self.pdf.expectedEvents ( self.vars ) / hh.sum() 
                
        histo += hh
            
        return histo 

    # ==========================================================================
    ## get the residual histogram : (data-fit) 
    #  @see PDF.as_histo
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
        - see PDF.as_histo
        - see PDF.residual_histo
        - see PDF.make_histo

        >>> data = ...
        >>> pdf  = ...
        >>> pdf.fitTo ( data )
        >>> residual = pdf.residual ( data , nbins = 100 ) 
        """
        hdata = self.make_histo ( **kwargs )
        dataset.project ( hdata , ( self.xvar.name , self.yvar.name , self.xvar.name )  )
        return self.residual_histo ( hdata ) 
        
    # ==========================================================================
    ## get the pull histogram : (data-fit)/data_error 
    #  @see PDF.as_histo
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
        - see PDF.as_histo
        - see PDF.residual_histo
        - see PDF.make_histo

        >>> data = ...
        >>> pdf  = ...
        >>> pdf.fitTo ( data )
        >>> residual = pdf.residual ( data , nbins = 100 ) 
        """
        hdata = self.make_histo ( **kwargs )
        dataset.project ( hdata , ( self.zvar.name , self.yvar.name , self.xvar.name ) ) 
        return self.pull_histo ( hdata ) 

    ## conversion to string 
    def __str__ (  self ) :
        return '%s(%s,xvar=%s,yvar=%s,zvar=%s)' % (
            self.__class__.__name__ , self.name ,
            self.xvar.name , self.yvar.name , self.zvar.name )
    __repr__ = __str__ 

    ## Make PDF3 object 
    def make_PDF3 ( self , pdf , xvar = None , yvar = None , zvar = None , *args , **kwargs ) :
        """Make PDF1 object
        """
        if isinstance ( pdf , PDF3 ) :
            
            assert ( not xvar ) or ( xvar in pdf.xvar ) , \
                   "make_PDF3: Invalid setting of xvar %s vs %s" % ( xvar , pdf.xvar )
            assert ( not yvar ) or ( yvar in pdf.xvar ) , \
                   "make_PDF3: Invalid setting of xvar %s vs %s" % ( yvar , pdf.yvar )
            assert ( not zvar ) or ( zvar in pdf.xvar ) , \
                   "make_PDF3: Invalid setting of xvar %s vs %s" % ( xvar , pdf.zvar )
          
            return pdf, pdf.xvar, pdf.yvar, pdf.zvar 
        
        elif isinstance ( pdf , ROOT.RooAbsPdf ) and xvar and yvar and zvar :
            
            return Generic3D_pdf ( pdf , xvar = xvar , yvar = yvar , zvar = zvar , *args , **kwargs ) , xvar , yvar , zvar 

        raise TypeError( "make_PDF2: invalid pdf/xvar %s/%s" % ( pdf , xvar ) )


# =============================================================================
## @class PDF3
#  The main helper base class for implementation of various 3D PDF-wrappers 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-08-21
class PDF3(APDF3,FUN3) :
    """The main helper base class for implementation of various 1D PDF-wrappers
    """
    def __init__ ( self , name , xvar , yvar , zvar , tricks = True , **kwargs ) :

        FUN3  .__init__ ( self            ,
                          name   = name   ,
                          xvar   = xvar   ,
                          yvar   = yvar   ,
                          zvar   = zvar   ,
                          tricks = tricks ,
                          fun    = kwargs.pop ( 'pdf' , None ) ,
                          **kwargs )
        APDF3 .__init__ ( self )
        
        self.config   = { 'name'    : self.name    ,
                          'xvar'    : self.xvar    ,
                          'yvar'    : self.yvar    ,
                          'zvar'    : self.yvar    ,
                          'tricks'  : self.tricks  ,
                          'pdf'     : self.pdf     }
        self.config.update ( kwargs )
        
        self.__call_OK = isinstance ( self.xvar , ROOT.RooAbsRealLValue ) and \
                         isinstance ( self.yvar , ROOT.RooAbsRealLValue ) and \
                         isinstance ( self.zvar , ROOT.RooAbsRealLValue )
        
    # ====================================================================================
    ## simple 'function-like' interface 
    def __call__ ( self , x , y , z , error = False , normalized = True ) :
        
        assert self.__call_OK , "Invalid types for xvar/yvar/zvar!"
        
        xmnmx = self.xminmax()
        if xmnmx :
            xmn , xmx = xmnmx 
            if not xmn <= x <= xmx : return 0

        ymnmx = self.yminmax()
        if ymnmx :
            ymn , yxmx = ymnmx 
            if not ymn <= y <= ymx : return 0

        zmnmx = self.zminmax()
        if zmnmx :
            zmn , zxmx = zmnmx 
            if not zmn <= z <= zmx : return 0

        with SETVAR ( self.xvar ) , SETVAR ( self.yvar ) , SETVAR ( self.zvar ) :
            self.xvar.setVal ( x )
            self.yvar.setVal ( y )
            self.zvar.setVal ( z )
            
            v = self.pdf.getVal ( self.vars ) if normalized else self.pdf.getValV ()
            
            if error and self.fit_result :
                e = self.pdf.getPropagatedError ( self.fit_result )
                if 0<= e : v = VE ( v ,  e * e )

            return v

    # ========================================================================
    ## convert to float 
    def __float__ ( self ) :
        """Convert to float
        >>> fun = ...
        >>> v   = float ( fun )
        """
        return self.fun.getVal ( self.vars ) 
    
    # =========================================================================
    ## make a product of two PDFs
    def __mul__ ( self , other ) :
        """Make a product of two PDFs
        """
        from ostap.fitting.pdf_ops import pdf2_product        
        return pdf2_product ( self , other )

    # =========================================================================
    ## make a product of two PDFs
    #  @code
    #  pdf1 = ...
    #  pdf2 = ...
    #  pdf  = pdf1 * pdf2 
    #  @endcode
    #  Rules: 
    #  - PDF3 ( x , y , z ) * PDF3 ( x , y , z ) -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF2 ( x , y )     -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF2 ( x , z )     -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF2 ( y , z )     -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF1 ( x )         -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF1 ( y )         -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF1 ( z )         -> PDF3 ( x , y , z ) 
    #  - PDF2 ( x , y )     * PDF3 ( ... )       -> process as PDF3 (...) * PDF2 ( x , y )
    #  - PDF2 ( x , y )     * PDF2 ( x , y )     -> PDF2 ( x , y )
    #  - PDF2 ( x , y )     * PDF1 ( x )         -> PDF2 ( x , y )
    #  - PDF2 ( x , y )     * PDF1 ( y )         -> PDF2 ( x , y )
    #  - PDF1 ( x )         * PDF3 ( ... )       -> process as PDF3 (...) * PDF ( x )
    #  - PDF1 ( x )         * PDF2 ( ... )       -> process as PDF2 (...) * PDF ( x )
    #  - PDF1 ( x )         * PDF1 ( x )         -> PDF1 ( x )
    #  - PDF1 ( x )         * PDF1 ( y )         -> PDF2 ( x , y )
    def __mul__ ( self , other ) :
        """Make a product of two PDFs
        
        >>> pdf1 = ...
        >>> pdf2 = ...
        >>> pdf  = pdf1 * pdf2

        Rules:
        - PDF3 ( x , y , z ) * PDF3 ( x , y , z ) -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF2 ( x , y )     -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF2 ( x , z )     -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF2 ( y , z )     -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF1 ( x )         -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF1 ( y )         -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF1 ( z )         -> PDF3 ( x , y , z ) 
        - PDF2 ( x , y )     * PDF3 ( ... )       -> process as PDF3 (...) * PDF2 ( x , y )
        - PDF2 ( x , y )     * PDF2 ( x , y )     -> PDF2 ( x , y )
        - PDF2 ( x , y )     * PDF1 ( x )         -> PDF2 ( x , y )
        - PDF2 ( x , y )     * PDF1 ( y )         -> PDF2 ( x , y )
        - PDF1 ( x )         * PDF3 ( ... )       -> process as PDF3 (...) * PDF ( x )
        - PDF1 ( x )         * PDF2 ( ... )       -> process as PDF2 (...) * PDF ( x )
        - PDF1 ( x )         * PDF1 ( x )         -> PDF1 ( x )
        - PDF1 ( x )         * PDF1 ( y )         -> PDF2 ( x , y )
        """
        from ostap.fitting.pdf_ops import pdf3_product        
        return pdf3_product ( self , other )

    # =========================================================================
    ## make a product of two PDFs
    #  @code
    #  pdf1 = ...
    #  pdf2 = ...
    #  pdf  = pdf1 * pdf2 
    #  @endcode
    #  Rules 
    #  - PDF3 ( x , y , z ) * PDF3 ( x , y , z ) -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF2 ( x , y )     -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF2 ( x , z )     -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF2 ( y , z )     -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF1 ( x )         -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF1 ( y )         -> PDF3 ( x , y , z )
    #  - PDF3 ( x , y , z ) * PDF1 ( z )         -> PDF3 ( x , y , z ) 
    #  - PDF2 ( x , y )     * PDF3 ( ... )       -> process as PDF3 (...) * PDF2 ( x , y )
    #  - PDF2 ( x , y )     * PDF2 ( x , y )     -> PDF2 ( x , y )
    #  - PDF2 ( x , y )     * PDF1 ( x )         -> PDF2 ( x , y )
    #  - PDF2 ( x , y )     * PDF1 ( y )         -> PDF2 ( x , y )
    #  - PDF1 ( x )         * PDF3 ( ... )       -> process as PDF3 (...) * PDF ( x )
    #  - PDF1 ( x )         * PDF2 ( ... )       -> process as PDF2 (...) * PDF ( x )
    #  - PDF1 ( x )         * PDF1 ( x )         -> PDF1 ( x )    
    #  - PDF1 ( x )         * PDF1 ( y )         -> PDF2 ( x , y )
    def __rmul__ ( self , other ) :
        """Make a product of two PDFs
        
        >>> pdf1 = ...
        >>> pdf2 = ...
        >>> pdf  = pdf1 * pdf2
        
        Rules:
        - PDF3 ( x , y , z ) * PDF3 ( x , y , z ) -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF2 ( x , y )     -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF2 ( x , z )     -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF2 ( y , z )     -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF1 ( x )         -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF1 ( y )         -> PDF3 ( x , y , z )
        - PDF3 ( x , y , z ) * PDF1 ( z )         -> PDF3 ( x , y , z ) 
        - PDF2 ( x , y )     * PDF3 ( ... )       -> process as PDF3 (...) * PDF2 ( x , y )
        - PDF2 ( x , y )     * PDF2 ( x , y )     -> PDF2 ( x , y )
        - PDF2 ( x , y )     * PDF1 ( x )         -> PDF2 ( x , y )
        - PDF2 ( x , y )     * PDF1 ( y )         -> PDF2 ( x , y )
        - PDF1 ( x )         * PDF3 ( ... )       -> process as PDF3 (...) * PDF ( x )
        - PDF1 ( x )         * PDF2 ( ... )       -> process as PDF2 (...) * PDF ( x )
        - PDF1 ( x )         * PDF1 ( x )         -> PDF1 ( x )
        - PDF1 ( x )         * PDF1 ( y )         -> PDF2 ( x , y )
        """
        from ostap.fitting.pdf_ops import pdf3_product
        if   isinstance ( other , PDF2 ) : return pdf3_product ( self , other )
        elif isinstance ( other , PDF1 ) : return pdf3_product ( self , other )
        return NotImplemented 

    # =========================================================================
    ## make a non-extender sum of 3D PDFs
    #  @code
    #  pdf1 = ...
    #  pdf2 = ...
    #  pdf  = pdf1 + pdf2     
    #  @endcode
    def __add__ ( self , other ) :
        """Make a no-extended sum of 3D PDFs
        >>> pdf1 = ...
        >>> pdf2 = ...
        >>> pdf  = pdf1 + pdf2     
        """
        from ostap.fitting.pdf_ops import pdf3_sum        
        return pdf3_sum ( self , other )

    # =========================================================================
    ## make a non-extended sum of 3D PDFs
    #  @code
    #  pdf1 = ...
    #  pdf2 = ...
    #  pdf  = pdf1 + pdf2     
    #  @endcode
    def __radd__ ( self , other ) :
        """Make a no-extended sum of 3D PDFs
        >>> pdf1 = ...
        >>> pdf2 = ...
        >>> pdf  = pdf1 + pdf2     
        """
        from ostap.fitting.pdf_ops import pdf3_sum        
        return pdf3_sum ( self , other )

    # =========================================================================
    ## Convert PDF into simple function
    #  @code
    #  pdf = ...
    #  fun = pdf.as_FUN () 
    #  @endcode
    def as_FUN ( self , name = '' ) : 
        """Convert PDF into simple function
        >>> pdf = ...
        >>> fun = pdf.as_FUN () 
        """
        return Fun3D ( self.pdf  ,
                       xvar = self.xvar  ,
                       yvar = self.yvar  ,
                       zvar = self.zvar  ,
                       name = name if name else self.new_name ( 'fun3' ) ) 
    

# =============================================================================
## @class Generic3D_pdf
#  "Wrapper" over generic RooFit (3D)-pdf
#  @code
#     
#  raw_pdf = 
#  pdf     = Generic3D_pdf ( raw_pdf )  
# 
#  @endcode 
#  If more functionality is required , more actions are possible:
#  @code
#  ## for sPlot 
#  pdf.alist2 = ROOT.RooArgList ( n1 , n2 , n3 ) ## for sPlotting 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-03-29
class Generic3D_pdf(PDF3) :
    """ Wrapper for generic (3D) RooFit pdf:    
    >>> raw_pdf =
    >>> x,y,z   = 
    >>> pdf     = Generic3D_pdf ( raw_pdf ,  xvar = x ,   yvar = y , zvar =  z) 
    """
    ## constructor 
    def __init__ ( self , pdf , xvar , yvar , zvar ,
                   name           = None    ,
                   add_to_signals = True    ,
                   prefix         = ''      ,
                   suffix         = ''      ) :
        
        assert isinstance ( xvar , ROOT.RooAbsReal ) , "'xvar' must be ROOT.RooAbsReal"
        assert isinstance ( yvar , ROOT.RooAbsReal ) , "'yvar' must be ROOT.RooAbsReal"
        assert isinstance ( zvar , ROOT.RooAbsReal ) , "'zvar' must be ROOT.RooAbsReal"
        assert isinstance ( pdf  , ROOT.RooAbsReal ) , "'pdf'  must be ROOT.RooAbsReal"

        name = name if name else self.generate_name ( prefix , suffix , pdf.GetName() )
        
        PDF3  . __init__ ( self , name , xvar , yvar , zvar )

        ## PDF! 
        self.pdf = pdf

        if not self.xvar in self.params () : 
            self.warning ( "Function/PDF does not depend on xvar=%s" % self.xvar.name )
        if not self.yvar in self.params () : 
            self.warning ( "Function/PDF does not depend on yvar=%s" % self.yvar.name )
        if not self.zvar in self.params () :
            self.warning ( "Function/PDF does not depend on zvar=%s" % self.zvar.name )

        ## add it to the list of signal components ?
        self.__add_to_signals = True if add_to_signals else False
        
        if self.add_to_signals :
            self.signals.add ( self.pdf )

        ## save the configuration
        self.config = {
            'pdf'            : self.pdf            ,
            'xvar'           : self.xvar           ,
            'yvar'           : self.yvar           ,
            'zvar'           : self.zvar           ,
            'name'           : self.name           ,
            'add_to_signals' : self.add_to_signals ,
            'prefix'         : prefix              ,
            'suffix'         : suffix              ,                        
            }

        self.checked_keys.add ( 'pdf'    )
        self.checked_keys.add ( 'xvar'   )
        self.checked_keys.add ( 'yvar'   )
        self.checked_keys.add ( 'zvar'   )
    
    @property
    def add_to_signals ( self ) :
        """'add_to_signals' : should PDF be added into list of signal components?"""
        return self.__add_to_signals 


# ===========================================================================
## @class Flat3D
#  The most trivial 3D-model - constant
#  @code 
#  pdf = Flat3D( 'flat' , xvar = ...  , yvar = ... , zvar = ... )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
class Flat3D(PDF3) :
    """The most trival 3D-model - constant
    >>> pdf = Flat3D( 'flat' , xvar = ...  , yvar = ... , zvar = ... )
    """
    def __init__ ( self , xvar , yvar , zvar , name = ''  , title = '' ) :

        name = name if name else self.generate_name ( prefix = 'Flat3D_')                            
        PDF3.__init__ ( self  , name , xvar , yvar , zvar ) 
        
        if not title : title = 'flat3(%s)' % name 
        self.pdf = Ostap.Models.Uniform ( name , title , self.xvar , self.yvar , self.zvar )
        assert 3 == self.pdf.dim() , 'Flat3D: wrong dimensionality!'
        
        ## save configuration
        self.config = {
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'zvar'     : self.zvar ,
            'name'     : self.name ,            
            'title'    : title     ,             
            }
        

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END 
# =============================================================================
