#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostap/fitting/pdfbasic.py
#  Set of useful basic utilities to build various fit models 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
# =============================================================================
""" Set of useful basic utilities to build various fit models"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'APDF1'            , ## useful base class for 1D-models
    'APDF2'            , ## useful base class for 2D-models
    'APDF3'            , ## useful base class for 3D-models
    ##
    'PDF1'             , ## useful base class for 1D-models
    'PDF2'             , ## useful base class for 2D-models
    'PDF3'             , ## useful base class for 3D-models
    ##
    'Generic1D_pdf'    , ## wrapper over imported RooFit (1D)-pdf
    'Generic2D_pdf'    , ## wrapper over imported RooFit (2D)-pdf
    'Generic3D_pdf'    , ## wrapper over imported RooFit (3D)-pdf
    ##
    'Shape1D_pdf'      , ## 1D-Generic/fixed shape from (C++) callable 
    'Shape2D_pdf'      , ## 2D-Generic/fixed shape from (C++) callable 
    'Shape3D_pdf'      , ## 3D-Generic/fixed shape from (C++) callable 
    ##
    'Histo1D_pdf'      , ## simple PDF from 1D-histogram
    'Histo2D_pdf'      , ## simple PDF from 2D-histogram
    'Histo3D_pdf'      , ## simple PDF from 3D-histogram
    ##
    'Constrained'      , ## mixin for creation of constrained PDFs
    'Constrained1D'    , ## 1D-PDF with constraints
    'Constrained2D'    , ## 1D-PDF with constraints
    'Constrained3D'    , ## 1D-PDF with constraints
    ##
    'make_pdf'         , ## helper function to make PDF
    'all_args'         , ## check that all arguments has correct type 
    ##
    'H1D_pdf'          , ## convertor of 1D-histo to RooHistPdf
    'H2D_pdf'          , ## convertor of 2D-histo to RooHistPdf
    'H3D_pdf'          , ## convertor of 3D-histo to RooHistPdf
    ##
    'Flat1D'           , ## uniform in 1D 
    'Flat2D'           , ## uniform in 2D 
    'Flat3D'           , ## uniform in 3D
    ## 
    'Sum1D'            , ## helper pdf: a non-extensive sum of 1D PDFs 
    'Sum2D'            , ## helper pdf: a non-extensive sum of 2D PDFs 
    'Sum3D'            , ## helper pdf: a non-extensive sum of 3D PDFs 
    ## 
)
# =============================================================================
from   ostap.core.meta_info     import root_info, python_info 
from   ostap.core.ostap_types   import ( is_integer     , string_types   , 
                                         integer_types  , num_types      ,
                                         list_types     , all_numerics   ) 
from   ostap.math.base          import ( iszero , isfinite , isequal , frexp10 ,  
                                         vct1_call_method ,
                                         vct2_call_method ,
                                         vct3_call_method ) 

from   ostap.core.core          import ( Ostap , VE , hID , dsID , rootID   ,
                                         valid_pointer , 
                                         rootException , 
                                         roo_silent    , rootWarning  )
from   ostap.fitting.utils      import ( RangeVar   , numcpu     ,
                                         make_name  , fit_status ,
                                         cov_qual   , get_i      )
from   ostap.fitting.fithelpers import ( H1D_dset, H2D_dset, H3D_dset ,
                                         SETPARS , Fractions          )
from   ostap.fitting.funbasic   import FUN1, FUN2, FUN3, Fun1D , Fun2D , Fun3D 
from   ostap.fitting.variables  import SETVAR
from   ostap.utils.cidict       import select_keys
from   ostap.utils.basic        import typename 
from   ostap.fitting.roocmdarg  import check_arg , nontrivial_arg , flat_args , command
from   ostap.logger.pretty      import pretty_float 
# 
import ostap.fitting.roofit 
import ostap.fitting.dataset 
import ostap.histos.histos
import ostap.fitting.roocollections 
# 
import ROOT, numpy, math,  random, sys, abc, warnings   
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
    """ Are all arguments of 'good' type?
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
## @class Components
#  Helper base class that keeps the major fit-pdf components
#
#  - signals
#  - backgrounds 
#  - other components 
#  - cross-terms1
#  - cross-terms2
#  - combined signals
#  - combined backgrounds 
#  - combined components 
# 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2023-04-18
class Components ( object ) :
    """ Helper base class that keeps the major fit-pdf components
    - signals
    - backgrounds 
    - other components 
    - cross-terms1
    - cross-terms2
    - combined signals
    - combined backgrounds 
    - combined components 
    """
    def __init__ ( self ) :

        ## major components 
        self.__signals               = ROOT.RooArgList ()
        self.__backgrounds           = ROOT.RooArgList ()
        self.__components            = ROOT.RooArgList ()
        self.__crossterms1           = ROOT.RooArgSet  () 
        self.__crossterms2           = ROOT.RooArgSet  () 
        self.__combined_signals      = ROOT.RooArgList ()
        self.__combined_backgrounds  = ROOT.RooArgList ()
        self.__combined_components   = ROOT.RooArgList ()

        ## needed for sPlot
        self.__alist1                = ROOT.RooArgList () 
        self.__alist2                = ROOT.RooArgList ()

    # =========================================================================
    ## copy (some) structural elements 
    def copy_structures ( self , source ) :
        """ Copy (some) structural elements"""
        
        for c in source.signals :
            if not c in self.signals :
                self.signals              .add ( c ) 
        for c in source.backgrounds :
            if not c in self.backgrounds :
                self.backgrounds          .add ( c ) 
        for c in source.components :
            if not c in self.components  :
                self.components           .add ( c ) 
        for c in source.crossterms1 :
            if not c in self.crossterms1 :
                self.crossterms1          .add ( c ) 
        for c in source.crossterms2 :
            if not c in self.crossterms2 :
                self.crossterms2          .add ( c ) 
        for c in source.combined_signals :
            if not c in self.combined_signals :
                self.combined_signals     .add ( c ) 
        for c in source.combined_backgrounds :
            if not c in self.combined_backgrounds :
                self.combined_backgrounds .add ( c ) 
        for c in source.combined_components  :
            if not c in self.combined_components  :
                self.combined_components  .add ( c ) 
        ##
        for c in source.alist1 :
            if not c in self.alist1  : self.alist1.add ( c ) 
        for c in source.alist2 :
            if not c in self.alist2  : self.alist2.add ( c ) 
        
    @property
    def signals     ( self ) :
        """ The list/ROOT.RooArgList of all 'signal' components,
        e.g. for visualization"""
        return self.__signals
    @property
    def backgrounds ( self ) :
        """ The list/ROOT.RooArgList of all 'background' components,
        e.g. for visualization"""
        return self.__backgrounds 
    @property
    def components  ( self ) :
        """ The list/ROOT.RooArgList of all 'other' components,
        e.g. for visualization"""
        return self.__components

    @property
    def combined_signals     ( self ) :
        """ The list/ROOT.RooArgList of all combined 'signal' components,
        e.g. for visualization"""
        return self.__combined_signals
    @property
    def combined_backgrounds ( self ) :
        """ The list/ROOT.RooArgList of all combined 'background' components,
        e.g. for visualization"""
        return self.__combined_backgrounds 
    @property
    def combined_components  ( self ) :
        """ The list/ROOT.RooArgList of all combined 'other' components,
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
        - Signal(y)*Background(x)     
          for 2D-fits,
        - Signal(x)*Background(y)*Background(z) for 3D-fits, etc...         
        """        
        return self.__crossterms2
        
    @property
    def alist1 ( self ) :
        """ The list/RooArgList of PDF components for compound PDF"""
        return self.__alist1
    @alist1.setter
    def alist1 ( self , value ) :
        assert isinstance ( value , ROOT.RooArgList ) , "Value must be RooArgList, %s/%s is  given" % ( value , type(value) )
        self.__alist1 = value
        
    @property
    def alist2 ( self ) :
        """ The list/RooArgList of PDF component's fractions (or yields/fractions for exteded fits) for compound PDF"""        
        return self.__alist2
    @alist2.setter
    def alist2 ( self , value ) :
        assert isinstance ( value , ROOT.RooArgList ) , "Value must be RooArgList, %s/%s is  given" % ( value , type(value) )
        self.__alist2 = value                    

    # =========================================================================
    ## List/tuple of structural components from `self.alist1` 
    @abc.abstractmethod
    def cmp_alist ( self ) :
        """ Generator of structural components from `self.alist` 
        """
        pass 

    # =========================================================================
    ## A list of structural components from `self.alist1`
    @property
    def alist_cmp ( self ) :
        """`alist_cmps` : a list of structural components from `self.alist1`
        """
        return self.cmp_alist()

    # ======================================================================
    ## add few rows into the table 
    def tab_rows ( self , dataset = None ) :
        """ Add few rows into the table 
        """
        ## get the tablle rows: 
        rows = super().tab_rows ( dataset )
        rows = list ( rows )

        if self.alist1 :
            row = '#components' , '%d' % len ( self.alist1 )
            rows.append ( row )
            for c in self.alist1 :
                n    = c.GetName() 
                what = '' 
                if   c in self.signals     : what = 'signal'
                elif c in self.backgrounds : what = 'bkg.'    
                elif c in self.crossterms1 : what = 'X-term-1'
                elif c in self.crossterms2 : what = 'X-term-2'
                elif c in self.components  : what = 'cmp.'
                ## 
                row = '   ' + n , typename ( c ) , what 
                rows.append ( row )
                
        return rows
    
    # =========================================================================
    ## Get all parameters from all signal components
    #  @code
    #  data  = ...
    #  model = ...    
    #  signal_pars = model.signal_parameters ( data  )
    #  @endcode 
    def signal_parameters ( self , dataset = None ) :
        """ Get all parameters from all signal components

        >>> data  = ...
        >>> model = ...    
        >>> signal_pars = model.signal_parameters ( data  )
        """        
        sigpars = ROOT.RooArgSet()
        
        if dataset is None or not dataset :
            dataset = self.vars if hasattr ( self , 'vars' ) else ROOT.nullptr             

        assert dataset is ROOT.nullptr or \
            ( isinstance ( dataset , ( ROOT.RooAbsData , ROOT.RooArgSet ) ) and dataset ) ,\
            "Invalid dataset/varset: %s" % typename ( dataset )
        
        for s in self.signals :
            pars = s.getParameters ( dataset )
            for p in pars :
                if not p in sigpars  : sigpars.add ( p ) 

        return sigpars 

    # =========================================================================
    ## Get all parameters from all background components
    #  @code
    #  data  = ...
    #  model = ...    
    #  bkg_pars = model.background_parameters ( data  )
    #  @endcode 
    def background_parameters ( self , dataset = None ) :
        """ Get all parameters from all signal components

        >>> data  = ...
        >>> model = ...    
        >>> bkg_pars = model.background_parameters ( data  )
        """        
        bkgpars = ROOT.RooArgSet()
        
        if dataset is None or not dataset :
            dataset = self.vars if hasattr ( self , 'vars' ) else ROOT.nullptr             

        assert dataset is ROOT.nullptr or \
            ( isinstance ( dataset , ( ROOT.RooAbsData , ROOT.RooArgSet ) ) and dataset ) ,\
            "Invalid dataset :%s" % typename ( dataset )
        
        for s in self.backgrounds :
            pars = s.getParameters ( dataset )
            for p in pars :
                if not p in bkgpars  : bkgpars.add ( p ) 

        return bkgpars 
    
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
class APDF1 ( Components ) :
    """ Useful helper base class for implementation of various PDF-wrappers 
    - it relies on `xvar`       method
    - it relies on `parse_args` method
    - it relies on `warning`    method
    - it relies on `info`       method
    - it relies on `error`      method
    - it relies on `fun`        attribute 
    """
    def __init__ ( self ) :

        ## initialize the base class 
        Components.__init__ ( self )
        
        ## take care about sPlots 
        self.__splots       = []
        self.__histo_data   = None        
        self.__fit_result   = None

    @property
    def pdf  ( self ) :
        """ The actual PDF (ROOT.RooAbsPdf)"""
        return self.fun 
    @pdf.setter
    def pdf  ( self , value ) :
        if value is None : self.fun = value
        else : 
            assert value and isinstance ( value , ROOT.RooAbsPdf ) , "'pdf' is not ROOT.RooAbsPdf"
            self.fun = value

    @property
    def roo_pdf ( self ) :
        """'roo_pdf' : get the underlying RooFit `RooAbsPdf` object (same as `pdf` here)"""
        return self.pdf

    @property
    def pdf_name ( self ) :
        """'pdf_name' : get the name of the underlying `RooAbsPdf` (same as 'fun_name') """
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
    def fit_result ( self ) :
        """'fit_result' : result of the latest call to `fitTo` method (or `None`"""
        return self.__fit_result
    @fit_result.setter
    def fit_result ( self , value ) :
        assert ( value is None ) or ( value and isinstance ( value , ROOT.RooFitResult ) and valid_pointer ( value ) ) , \
               "Invalid value for 'fit-result' object! %s/%s" %( value , type ( value ) ) 
        self.__fit_result = value
            
    @property
    def histo_data  ( self ):
        """ Histogram representation as DataSet (RooDataSet)"""
        return self.__histo_data
    @histo_data.setter
    def  histo_data ( self  , value ) :
        if   value is None :
            self.__histo_data = value 
        elif hasattr ( value , 'dset' ) and isinstance ( value.dset , ROOT.RooAbsData ) :
            self.__histo_data = value 
        else :
            raise AttributeError("'histo_data' has invalid type %s/%s" % (   value , type(value) ) )
        
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
        """ Perform the actual fit (and draw it)
        >>> r , f = model.fitTo ( dataset )
        >>> r , f = model.fitTo ( dataset , weighted = True )    
        >>> r , f = model.fitTo ( dataset , ncpu     = 10   )    
        >>> r , f = model.fitTo ( dataset , draw = True , nbins = 300 )    
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
            
        
        if   isinstance ( dataset , ( H1D_dset , H2D_dset , H3D_dset ) ) : dataset = dataset.dset        
        elif isinstance ( dataset , ROOT.TH1 )  :
            assert 1 == dataset.GetDimension () , 'Invalid histogram dimension: %s' % dataset.GetDimension() 
            density = kwargs.pop ( 'density' , False  )
            chi2    = kwargs.pop ( 'chi2'    , False  )
            return self.fitHisto ( dataset           ,
                                   draw    = draw    ,
                                   silent  = silent  ,
                                   density = density ,
                                   nbins   = nbins   ,
                                   refit   = refit   , 
                                   chi2    = chi2    , args = args , **kwargs )
        #
        ## treat the arguments properly
        #
        opts     = self.fit_options + ( ROOT.RooFit.Save () , ) + args 
        opts     = self.parse_args ( dataset , *opts , **kwargs )
        
        if silent :
            pl = check_arg ( 'PrintLevel'      , *opts ) 
            if not pl : opts = opts + ( ROOT.RooFit.PrintLevel      ( -1    ) , )
            vl = check_arg ( 'Verbose'         , *opts )
            if not vl : opts = opts + ( ROOT.RooFit.Verbose         ( False ) , )
            pe = check_arg ( 'PrintEvalErrors' , *opts )
            if not pe : opts = opts + ( ROOT.RooFit.PrintEvalErrors ( 0     ) , )

        if self.vars and 1 < len ( self.vars ) :
            rng = check_arg ( 'Range' , *opts ) 
            if rng : self.warning ( 'fitTo: %s is specified for >1D function - it is ambuguous!' % rng )

        ## check sumw2/asymptoticerorr flags 
        weighted    = dataset.isWeighted() if dataset else False
        weighted_ok = weighted 
        if weighted :
            sw2 = check_arg ( 'SumW2Error'      , *opts )
            aer = check_arg ( 'AsymptoticError' , *opts )
            if not sw2 and not aer :
                weighted_ok = False 
                self.warning ( "fitTo: Neither 'SumW2Error' and 'AsymptoticError' are specified for weighted dataset!" )

        ## check fit ranges 
        rng = check_arg ( 'RangeByName' , *opts )
        ok  = self.check_ranges ( dataset , rng.getString ( 0 ) if rng else '' )
        if not ok : self.warning ( 'fitTo: ranges are not OK (Dataset has entries outside fitting range)' ) 
    
        ## #
        ## ## check the limits/ranges 
        ## #
        ## if dataset and self.vars and \
        ##        ( not check_arg ( 'Range'       , *opts ) ) and \
        ##        ( not check_arg ( 'RangeByName' , *opts ) ) :
        ##     vnames = [ v.name for v in self.vars if v.xminmax() ]
        ##     if vnames :
        ##         stat = dataset.statVars ( vnames )
        ##         for v in self.vars :
        ##             if not v.name in vnames : continue 
        ##             vmn , vmx = v.xminmax()
        ##             mnv , mxv = stat [ v.name ].minmax()
        ##             if mnv < vmn or vmx < mxv :
        ##                 self.warning ( 'Limits for %s variable [%+.3g,%+.3g], while the range is [%+.3g,%+.3g]' % ( v.name    ,
        ##                                                                                                             vmn , vmx ,
        ##                                                                                                             mnv , mxv ) )

        if opts  and ( not silent and not 'quiet' in kwargs ) and nontrivial_arg ( ( 'Save' , 'NumCPU' ) , *opts ) :
            rows = [ ( 'Option' , ) ]
            for o in opts :
                row = str ( o ) ,
                rows.append ( row )
            import ostap.logger.table as T
            title = 'fitTo: Fit-Options'
            table = T.table ( rows , title = 'FitOptions' , prefix = '# ' )
            self.info ( '%s:\n%s' % ( title , table ) )
                        
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

        if not silent and not 'quiet' in kwargs  :
            params = self.params ( dataset )
            rows = [  ( 'Parameter' , 'Value' , 'Factor' ) ]
            for p in params :
                name = p.name 
                v    = float ( p * 1 )
                v , expo = pretty_float ( v , precision = 4 , width = 6 )
                if expo : row = name , v , '10^%+d' % expo  
                else    : row = name , v , '' 
                rows.append ( row )
            import ostap.logger.table as T
            rows  = T.remove_empty_columns ( rows ) 
            title = 'Before fit'
            table = T.table ( rows , title = title , prefix = '# ' , alignment = 'rl' )            
            self.info ( '%s:\n%s' % ( title , table ) ) 
            
        ## define silent context
        with roo_silent ( silent ) :            
            if self.fit_result : self.fit_result = None
            ## actual fit! 
            result = self.fit_to ( self.pdf , dataset , *opts )
            if result and valid_pointer ( result ) : self.fit_result = result
            else :
                self.fatal ( "fitTo: RooFitResult object is invalid! Check model&data!" )
                self.fit_result = None             
                return None , None
            
            if hasattr ( self.pdf , 'setPars' ) : self.pdf.setPars() 

        ## result = Ostap.MoreRooFit.delete_result ( result ) 
        ## return result, None
    
        st = result.status()
        
        for_refit = '' 
        if 0 != st   :
            for_refit = 'status is %s' % fit_status ( st ) 
            ## if not silent : self.warning ( 'fitTo: Fit status is %s ' % fit_status ( st ) )
        #
        qual = result.covQual()
        cov2_good = ( qual == 3 ) or ( weighted_ok and qual == -1 )
        if not cov2_good :
            fr2 = 'covariance is %s' % cov_qual ( qual )  
            for_refit = fr2 if not for_refit else '%s & %s' % ( for_refit , fr2 )
            ## if not silent : self.warning ( 'fitTo: covQual    is %s ' % cov_qual ( qual ) )

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
                        for_refit = 'integral' if not for_refit else  '%s & %s' % ( for_refit , 'integral' )
        #
        ## call for refit if needed
        #
        if for_refit and 0 < refit :
            if not silent : self.warning ( 'fitTo: call for refit #%-2d reason: %s' % ( refit , for_refit ) )
            refit -= 1 
            return  self.fitTo ( dataset         ,
                                 draw   = draw   ,
                                 nbins  = nbins  ,
                                 silent = silent ,
                                 refit  = refit  ,
                                 args   = args   , **kwargs )
            
        if   result and  0 == result.status () and not silent :
            self.info    ( "Fit result is\n%s" % result.table ( prefix = "# " ) ) 
        elif result and  0 != result.status () and not silent : 
            self.warning ( "Fit result is\n%s" % result.table ( prefix = "# " ) ) 
        elif result and  ( not cov2_good )     and not silent : 
            self.warning ( "Fit result is\n%s" % result.table ( prefix = "# " ) ) 
        elif result and not silent :
            self.warning ( "Fit result is\n%s" % result.table ( prefix = "# " ) )

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

    # =========================================================================
    ## invoke <code>model.fitTo  ( data  , *options)</code> command
    #  - merge arguments using <code>RooFit::MultiArg</code> to shorted list
    def fit_to ( self , model , data , *options ) :
        """ Invoke `model.fitTo ( data , *options)` command
        - merge arguments using `ROOT.RooFit::MultiArg` to shorted list
        """
        assert all ( isinstance ( o , ROOT.RooCmdArg ) for o in options  ), \
            "fit_to: invalid argument types: %s" % list ( options  ) 
        
        if 'fitrange' in model.all_string_attributes() :
            model.removeStringAttribute ( 'fitrange' )

        if ( 3 , 11 ) <= python_info : kw = { 'category' : RuntimeWarning }
        else                         : kw = {}
        
        # ================================================================
        ## No need to restructure the options
        if  ( 6 , 32 ) <= root_info :            
            # ============================================================
            try : # ======================================================
                # ========================================================
                ## convert ROOT errors to exceptions and ROOT warnings to warnings  
                with rootException() , warnings.catch_warnings ( **kw ) :
                    warnings.simplefilter ( 'always' , RuntimeWarning ) 
                    return Ostap.MoreRooFit.fitTo ( model , data , *options )                    
                # ========================================================                
            except Exception : # =========================================
                # ========================================================
                self.error ( "Exception from %s/%s fit_tp\n%s" % ( typename ( self  ) , 
                                                                   typename ( model ) , 
                                                                   self               ) , exc_info = True )
                raise
            
        """
        # ================================================================
        ## special treatment for old ROOT: it does not like too many arguments 
        elif root_info < ( 6 , 29 ) : # ==================================
            # ============================================================
            NARGMAX = 8
            # ============================================================
            ##  for "small" number of arguments use the standard function 
            if len ( options ) <= NARGMAX  :
                # ========================================================
                try : # ==================================================
                    # ====================================================
                    ## convert ROOT errors to exceptions and ROOT warnings to warnings  
                    with rootException() , warnings.catch_warnings ( **kw ) :
                        warnings.simplefilter ( 'always' , RuntimeWarning ) 
                        return model.fitTo ( data , *options )
                    # ====================================================
                except Exception : # =====================================
                    # ====================================================
                    self.error ( "Exception from %s/%s fit_to\n%s" % ( typename ( self  ) ,
                                                                       typename ( model ) , 
                                                                       self             ) , exc_info = True )
                    raise
        """
                
        # ================================================================
        ## merge arguments into linked list 
        from ostap.fitting.roocmdarg import command 
        cmd = command ( *options )
        # ================================================================
        try : # ==========================================================
            # ============================================================
            ## convert ROOT errors to exceptions and ROOT warnings to warnings  
            with rootException() , warnings.catch_warnings ( **kw ) :
                warnings.simplefilter ( 'default' , RuntimeWarning ) 
                return Ostap.MoreRooFit.fitTo ( model , data , cmd  )
            # ============================================================                        
        except Exception : # =============================================
            # ============================================================
            self.error ( "Exception from %s/%s fit_to\n%s" % ( typename ( self  ) , 
                                                               typename ( model ) , 
                                                           self             ) , exc_info = True )
            raise
    
    # ===============================================================================
    ## refit up-to N-times till good fit result/cov matrix quality 
    def reFit ( self , *args , **kwargs ) :
        """ efit up-to N-times till good fit result/cov nmarix quality 
        """
        N           = kwargs.pop ( "repeat"      , 1 )
        good_status = kwargs.pop ( "good_status" , 0 )
        good_cov    = kwargs.pop ( "good_cov"    , ( 3 , -1 ) ) 
        ## 
        if isinstance ( good_status , int ) : good_status = good_status ,
        if isinstance ( good_cov    , int ) : good_cov    = good_cov    ,
        
        assert isinstance ( N , it ) and 1 <= N , "reFit: Invalid N=%s" % N
        assert all ( isinstance ( a , int ) and  0 <= a <=5 for a in good_status ) , \
            "Invalid `good_status' list: %s" % str ( good_status )
        assert all ( isinstance ( a , int ) and -1 <= a <=3 for a in good_covs   ) , \
            "Invalid `good_cov' list: %s" % str ( good_cov )
        
        for n in range ( N ) :
            
            r , f = self.fitTo ( *args , **kwargs )
            s = r.status()
            c = r.covQual() 
            if s in good_status and c in good_covs : return r , f
            
        self.warning ( "reFit: status:%s covQual:%s" % ( fit_status ( s ) , cov_qual   ( c ) ) )
        return r , f 
        
    # ================================================================================
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
    ## Helper method to draw total fit curve 
    def _draw_total  ( self  , frame  , *args , **kwargs ) :
        """ Helper method to draw total fit curve 
        """
        totoptions   = self.draw_option (  'total_fit_options' , **kwargs )
        self.plot_on ( self.pdf , frame , *totoptions ) 
        
    # ================================================================================
    ## helper method to draw "signal-like" components
    def _draw_signals ( self , frame  , *args , **kwargs ) :
        """ Helper method to draw `signal-like' components
        """

        if self.combined_signals or self.signals :
            drawit1   = self.draw_option ( 'draw_combined_signal'    , **kwargs )
            drawit2   = self.draw_option ( 'draw_signals'            , **kwargs )
            if ( not drawit1 ) and ( not  drawit2 ) : 
                self.warning ( "Plotting of `signal-like' components is disabled!")
                
        ## draw "combined-signal"
        if self.combined_signals and drawit1 :
            doptions  = self.draw_option ( 'combined_signal_options'  , **kwargs ) 
            dstyle    = self.draw_option ( 'combined_signal_style'    , **kwargs )
            self._draw ( self.combined_signals , frame , doptions , dstyle , args )            
            
        ## draw "signal" components
        if self.signals and drawit2 : 
            soptions     = self.draw_option (    'signal_options'  , **kwargs )
            sbstyle      = self.draw_option (      'signal_style'  , **kwargs ) 
            self._draw ( self.signals , frame , soptions , sbstyle , args )
            
    # ================================================================================
    ## helper method to draw "crossterm-like" components
    def _draw_crossterms ( self , frame , *args , **kwargs ) :
        """ Helper method to draw `crossterm-like' components
        """

        if self.crossterms1 or self.crossterms2 : 
            drawit1   = self.draw_option ( 'draw_crossterm1'    , **kwargs )
            drawit2   = self.draw_option ( 'draw_crossterm2'    , **kwargs )
            if ( not drawit1 ) and ( not  drawit2 ) : 
                self.warning ( "Plotting of `crossterm-like' components is disabled!")
                
        if self.crossterms1 and drawit1 : 
            ct1options   = self.draw_option ( 'crossterm1_options' , **kwargs )
            ct1bstyle    = self.draw_option ( 'crossterm1_style'   , **kwargs ) 
            self._draw ( self.crossterms1 , frame , ct1options , ct1bstyle , args )
            
        if self.crossterms2 and drawit2 :
            ct2options   = self.draw_option ( 'crossterm2_options' , **kwargs )
            ct2bstyle    = self.draw_option ( 'crossterm2_style'   , **kwargs )        
            self._draw ( self.crossterms2 , frame , ct2options , ct2bstyle , args )
                
    # ================================================================================
    ## helper method to draw "other" components
    def _draw_components ( self , frame , *args , **kwargs ) :
        """ Helper method to draw `other' components
        """

        if self.combined_components or self.components :
            drawit1   = self.draw_option ( 'draw_combined_component'    , **kwargs )
            drawit2   = self.draw_option ( 'draw_components'            , **kwargs )
            if ( not drawit1 ) and ( not  drawit2 ) : 
                self.warning ( "Plotting of `other' components is disabled!")
        
        ## draw combined "other" components 
        if self.combined_components  and drawit1 :
            doptions = self.draw_option ( 'combined_component_options' , **kwargs ) 
            dstyle   = self.draw_option ( 'combined_component_style'   , **kwargs )                
            self._draw ( self.combined_components , frame , doptions , dstyle , args )
            
        ## draw "other" components
        if self.components and drawit2 : 
            coptions     = self.draw_option ( 'component_options' , **kwargs )
            cbstyle      = self.draw_option ( 'component_style'   , **kwargs )
            self._draw ( self.components , frame , coptions , cbstyle , args )

    # ================================================================================
    ## helper method to draw "background-like" components
    def _draw_backgrounds ( self , frame , *args , **kwargs ) :
        """ Helper method to draw `background-like' components
        """
        
        if self.combined_backgrounds or self.backgrounds :
            drawit1   = self.draw_option ( 'draw_combined_background'   , **kwargs )
            drawit2   = self.draw_option ( 'draw_backgrounds'           , **kwargs )
            if ( not drawit1 ) and ( not  drawit2 ) : 
                self.warning ( "Plotting of `background-like' components is disabled!")
        
        ## draw combined "background" components 
        if self.combined_backgrounds and drawit1 :
            doptions = self.draw_option ( 'combined_background_options' , **kwargs ) 
            dstyle   = self.draw_option ( 'combined_background_style'   , **kwargs )                
            self._draw ( self.combined_backgrounds , frame , doptions , dstyle , args )

        ## draw various "background" terms
        if self.backgrounds and drawit2 :
            boptions     = self.draw_option ( 'background_options' , **kwargs ) 
            bbstyle      = self.draw_option (   'background_style' , **kwargs )
            self._draw ( self.backgrounds , frame , boptions , bbstyle , args )
            
    # ================================================================================
    ## draw fit results
    #  @code
    #  r,f = model.fitTo ( dataset )
    #  model.draw ( dataset , nbins = 100 ) 
    #  @endcode
    #  @param dataset  dataset to be drawn 
    #  @param nbins    binning scheme for frame/RooPlot 
    #  @param silent   silent mode ?
    #  @param data_options             drawing options for dataset
    #  @param signal_options           drawing options for `signal'        components    
    #  @param background_options       drawing options for `background'    components 
    #  @param crossterm1_options       drawing options for `crossterm-1'   components 
    #  @param crossterm2_options       drawing options for `crossterm-2'   components 
    #  @param background2D_options     drawing options for `background-2D' components 
    #  @param component_options        drawing options for 'other'         components 
    #  @param fit_options              drawing options for fit curve    
    #  @param signal_style             style(s) for signal components 
    #  @param background_style         style(s) for background components
    #  @param component_style          style(s) for other components
    #  @param crossterm1_style         style(s) for "crossterm-1"   components
    #  @param crossterm2_style         style(s) for "crossterm-2"   components
    #  @param background2D_style       style(s) for "background-2D" components
    #  @param draw_order               order of drawing components
    #  @param draw_combined_signal     flag to draw combined signal component
    #  @param draw_signals             flag to draw signal components 
    #  @param draw_crossterm1          flag to draw crossterm-like components
    #  @param draw_crossterm2          flag to draw crossterm-like components
    #  @param draw_combined_component  flag to draw combined 'other' component
    #  @param draw_components          flag to draw 'other' components
    #  @param draw_combined_background flag to draw combined 'background-like' component
    #  @param draw_background          flag to draw 'background-like' components
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
        """ Visualize the fits results
        >>> r,f = model.draw ( dataset )
        >>> model.draw ( dataset , nbins = 100 )
        >>> model.draw ( dataset , base_signal_color  = ROOT.kGreen+2 )
        >>> model.draw ( dataset , data_options = (ROOT.RooFit.DrawOptions('zp'),) )

        Produce also residual & pull:
        
        >>> f,r,p = model.draw ( dataset , nbins = 100 , residual = 'P' , pull = 'P')
        
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
        -  draw_order             ## order of drawing components

        Flags:
        
        - draw_combined_signal     ## flag to draw combined signal component
        - draw_signals             ## flag to draw signal components 
        - draw_crossterm1          ## flag to draw crossterm-like components
        - draw_crossterm2          ## flag to draw crossterm-like components
        - draw_combined_component  ## flag to draw combined 'other' component
        - draw_components          ## flag to draw 'other' components
        - draw_combined_background ## flag to draw combined 'background-like' component
        - draw_background          ## flag to draw 'background-like' components


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

            binned  = dataset and isinstance ( dataset , ROOT.RooDataHist )

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
                commands = data_options + args +  ( ROOT.RooFit.Invisible() , ) 
                self.plot_on ( dataset , frame , *commands ) 

            ## drawing order
            import ostap.plotting.fit_draw as FD 
            draw_order = self.draw_option ( 'draw_order' , **kwargs )
            draw_order = FD.get_draw_order ( draw_order )

            if ( self.signals or self.combined_signals         ) and not 'S' in draw_order :
                self.info ( "Plotting of `signal-like'    components is omitted" )  
            if ( self.crossterms1 or self.crossterms2          ) and not 'X' in draw_order :
                self.info ( "Plotting of `crossterm-like' components is omitted" )  
            if ( self.components or self.combined_components   ) and not 'C' in draw_order :
                self.info ( "Plotting of `other'          components is omitted" )  
            if ( self.backgrounds or self.combined_backgrounds ) and not 'B' in draw_order :
                self.info ( "Plotting of `other'          components is omitted" )  
            if not 'T' in draw_order :
                self.info ( "Plotting of `total' curve               is omitted" )  
            if dataset  and not 'D' in draw_order :
                self.info ( "Plotting of `data'                      is omitted" )  

            ## now draw the classified components
            drargs = args

            for c in draw_order :
                if   'S' == c :
                    self._draw_signals     ( frame , *drargs , **kwargs )
                    used_options.add ( 'draw_combined_signal'        )
                    used_options.add ( 'draw_signals'                )
                    used_options.add ( 'combined_signal_options'     )
                    used_options.add ( 'combined_signal_style'       )
                    used_options.add ( 'signal_options'              )
                    used_options.add ( 'signal_style'                )
                elif 'X' == c :
                    self._draw_crossterms  ( frame , *drargs , **kwargs )
                    used_options.add ( 'draw_crossterm1'             )
                    used_options.add ( 'draw_crossterm2'             )
                    used_options.add ( 'crossterm1_options'          )
                    used_options.add ( 'crossterm1_style'            )
                    used_options.add ( 'crossterm2_options'          )
                    used_options.add ( 'crossterm2_style'            )                    
                elif 'C' == c :
                    self._draw_components  ( frame , *drargs , **kwargs )
                    used_options.add ( 'draw_combined_component'     )
                    used_options.add ( 'draw_components'             )
                    used_options.add ( 'combined_component_options'  )
                    used_options.add ( 'combined_component_style'    )
                    used_options.add ( 'component_options'           )
                    used_options.add ( 'component_style'             )                                        
                elif 'B' == c :
                    self._draw_backgrounds ( frame , *drargs , **kwargs )
                    used_options.add ( 'draw_combined_background'    )
                    used_options.add ( 'draw_backgrounds'            )
                    used_options.add ( 'combined_background_options' )
                    used_options.add ( 'combined_background_style'   )
                    used_options.add ( 'background_options'          )
                    used_options.add ( 'background_style'            )                                                            
                elif 'T' == c :
                    self._draw_total       ( frame , *drargs , **kwargs )
                    used_options.add ( 'total_fit_options'    ) 
                elif 'D' == c :
                    if dataset : 
                        ## draw data once more
                        commands = data_options + args
                        self.plot_on ( dataset , frame , *commands )
                else :
                    self.warning ( "Invalid draw-order component `%s', skip it!" % c )
                                
            for k in FD.keys : used_options.add ( k ) 

            #
            ## suppress ugly axis labels
            #                
            if frame and not kwargs.get ( 'draw_axis_title' , False ) :  
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
            frame.chi2ndf = None 
            if dataset and not silent :
                if self.fit_result :
                    npf = len ( self.fit_result.floatParsFinal() )
                    frame.chi2ndf = frame.chiSquare ( npf )
                else : 
                    npf = len ( self.params ( dataset ) ) 
                    frame.chi2ndf = frame.chiSquare ( npf )
                ## 
                binw          = -1 
                if nbins and isinstance ( nbins , integer_types ) and 1 < nbins :
                    if hasattr ( drawvar , 'xminmax' ) and drawvar.xminmax () :
                        xmn , xmx =  drawvar.xminmax()
                        binw = ( xmx - xmn ) / float ( nbins )
                if 0 < binw : self.info ( 'chi2/ndf: %.3f, binwidth: %.4g' % ( frame.chi2ndf , binw ) )
                else        : self.info ( 'chi2/ndf: %.3f' %                   frame.chi2ndf          )


            if not residual and not pull:
                frame = frame.copy () ## a bit strange action but it helps to avoid decolorization/reset for the last created frame
                if frame and not kwargs.get ( 'draw_axis_title' , False ) : 
                    frame.SetXTitle ( '' )
                    frame.SetYTitle ( '' )
                    frame.SetZTitle ( '' )                    
                return frame

            rframe =  None 
            if residual and frame :
                if   residual is True               : residual =      "P" ,
                elif isinstance  ( residual , str ) : residual = residual ,
                rframe  = frame.emptyClone ( rootID ( 'residual_' ) )
                rh      = frame.residHist()
                rframe.addPlotable ( rh , *residual )
                
                rframe = rframe.copy()                
                if not kwargs.get( 'draw_axis_title' , False ) :  
                    rframe.SetXTitle ( '' )
                    rframe.SetYTitle ( '' )
                    rframe.SetZTitle ( '' )
                    
            pframe = None 
            if pull  and frame : 
                if   pull is True               : pull =   "P",
                elif isinstance  ( pull , str ) : pull = pull ,
                pframe  = frame.emptyClone ( rootID ( 'pull_' ) )
                ph      = frame.pullHist()
                pframe.addPlotable ( ph , *pull )
                
                pframe = pframe.copy()                
                if not kwargs.get( 'draw_axis_title' , False ) :  
                    pframe.SetXTitle ( '' )
                    pframe.SetYTitle ( '' )
                    pframe.SetZTitle ( '' )
                                    
            ## a bit strange action but it helps to avoid decolorization/reset for the last created frame
            frame  =  frame.copy () if  frame else  frame 
            rframe = rframe.copy () if rframe else rframe 
            pframe = pframe.copy () if pframe else pframe 
            
            if frame and not kwargs.get ( 'draw_axis_title' , False ) : 
                frame.SetXTitle ( '' )
                frame.SetYTitle ( '' )
                frame.SetZTitle ( '' )                    
                
            if rframe and not kwargs.get ( 'draw_axis_title' , False ) : 
                rframe.SetXTitle ( '' )
                rframe.SetYTitle ( '' )
                rframe.SetZTitle ( '' )
                
            if pframe and not kwargs.get ( 'draw_axis_title' , False ) : 
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
                   refit   = False , 
                   args    = () , **kwargs ) :
        """ Fit the 1D-histogram (and draw it)

        >>> histo = ...
        >>> r,f = model.fitHisto ( histo , draw = True ) 
        
        """
        with RangeVar ( self.xvar , *(histo.xminmax()) ) :
            
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
                                              refit   = False    ,  
                                              args    = args     , **kwargs )
            else    : return self.fitTo     ( data               ,
                                              draw    = draw     ,
                                              nbins   = nbins    , 
                                              silent  = silent   ,
                                              refit   = refit    , 
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
                    refit   = False , 
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
            
            if silent :
                pl = check_arg ('PrintLevel'      , *lst1 ) 
                if not pl : lst1.append ( ROOT.RooFit.PrintLevel      ( -1    ) ) 
                vl = check_arg ('Verbose'         , *lst1 )
                if not vl : lst1.append ( ROOT.RooFit.Verbose         ( False ) )
                pe = check_arg ('PrintEvalErrors' , *lst1 )
                if not pe : lst1.append ( ROOT.RooFit.PrintEvalErrors ( 0     ) )

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
            ## 
            if root_info < ( 6 , 28 ) :  m = ROOT.RooMinuit    ( chi2 )
            else                      :  m = ROOT.RooMinimizer ( chi2 )
            ##
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
        """ Draw/prepare NLL or LL-profile for seleted variable:        
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
                self.histo_data   = H1D_dset ( histo = dataset , xaxis = self.xvar , density = density , silent = silent )
                hdataset          = self.histo_data.dset
                kwargs [ 'ncpu' ] = 2   
                return self.draw_nll ( var     = var      ,
                                       dataset = hdataset ,
                                       profile = profile  ,                                       
                                       draw    = draw     ,
                                       silent  = silent   , 
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
        nbins  = kwargs.pop ( 'nbins' , 25 if profile else 200 )
        if nbins  : fargs.append ( ROOT.RooFit.Bins      ( nbins  ) )
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
        """ Create NLL object from the pdf
        >>> model.fitTo ( dataset , ... )
        >>> nll, sf = model.nll ( dataset )
        - see RooAbsPdf::createNLL 
        """
        
        ## convert if needed 
        if not isinstance ( dataset , ROOT.RooAbsData ) and hasattr ( dataset , 'dset' ) :
            dataset = dataset.dset 

        if root_info < ( 6 , 28 ) : 
            clone = kwargs.pop ( 'clone' , False )
            kwargs [ 'clone' ] = clone 
            if not clone : kwargs [ 'optimize' ] = False 

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
        
        if not silent and opts and nontrivial_arg ( ( 'Save' , 'NumCPU' ) , *opts ) :
            rows = [ ( 'Option' , ) ]
            for o in opts :
                row = str ( o ) ,
                rows.append ( row )
            import ostap.logger.table as T
            title = 'nll: createNLL options'
            table = T.table ( rows , title = 'CreateNLL options' , prefix = '# ' )
            self.info ( '%s:\n%s' % ( title , table ) )

        ## get s-Factor 
        sf   = dataset.sFactor() 

        if ( 6 , 32 ) <= root_info :
            return Ostap.MoreRooFit.createNLL ( self.pdf , dataset , *opts ), sf 
            
        if len ( opts ) < 8 and root_info < ( 6 , 29 ) : 
            return self.pdf.createNLL ( dataset , *opts   ) , sf
        else :
            cmdlist =  command ( *opts )
            return self.pdf.createNLL ( dataset , cmdlist  ) , sf
            
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
        """ Get NLL/profile-graph for the variable, using the specified abscissas
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
            for v in progress_bar  ( values , silent = silent , description = 'Points:' ) :
                var.setVal ( v )
                n   = nLL.getVal() 
                res = v , n
                results.append ( res )
                vmin = n if vmin is None else min ( vmin , n )
                if draw and graph :
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
            
        if draw : graph.draw ('apl')

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
                        silent   = False ,
                        draw     = False ,
                        subtract = True  , 
                        args     = ()    , **kwargs ) :
        """ Get profile-graph for the variable, using the specified abscissas
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
        pname  = self.roo_name ( 'pLL_%s_%s'         % ( var.name , self.name ) )
        ptitle =                 'LL-profile(%s,%s)' % ( var.name , self.name ) 

        if subtract : 
            pLL = ROOT.RooProfileLL          ( pname , ptitle , nLL , vars )
        else :
            pLL = Ostap.MoreRooFit.ProfileLL ( pname , ptitle , nLL , vars )
            
        ## 2) create graph if requested 
        import ostap.histos.graphs
        graph = ROOT.TGraph () if draw else None 
                        
        ## 3) collect pLL values 
        results = [] 
        with SETPARS ( self , dataset ) , RangeVar ( var , minv , maxv ) , SETVAR  ( var ) :
            from ostap.utils.progress_bar import progress_bar 
            for i , v in enumerate ( progress_bar ( vals , silent = silent , description = 'Points:' )  ) :
                var.setVal ( v )
                p   = pLL.getVal() 
                res = v , p 
                results.append ( res )
                if draw :
                    graph.SetPoint ( len ( graph ) , v , p ) ## add the point 
                    if   1 == len ( graph ) : graph.draw ("ap")
                    elif 1 <  len ( graph ) : 
                        cnv = Ostap.Utils.get_canvas()
                        if cnv : cnv.Update()
                        
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

        if draw : graph.draw ('apl')
        
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
        """ Evaluate 'significance' using Wilks' theorem via NLL
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
                        var.setVal ( max ( minv , maxv + error ) )
                        val_maxvp = nll.getVal()

                    with SETVAR ( var ) :  
                        var.setVal ( max ( minv , maxv - error ) ) 
                        val_maxvm = nll.getVal()
                        
                    dnll = VE ( dnll , 0.25 * (val_maxvp - val_maxvm ) ** 2 )
                    
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
        """ Evaluate 'significance' using Wilks' theorem via NLL
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
                self.histo_data = H1D_dset ( histo = dataset , xaxis = self.xvar , density = density , silent = silent )
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

        if fixed and not silent :
            self.info ( "Wilks2: fixed variables: %s" % [ f.GetName()  for f in fixed] )

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
                    
                    var.setVal ( max ( minv , maxv + error ) )
                    ve1 = pLL.getVal()
                    
                    var.setVal ( max ( minv , maxv - error ) )
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
        """ Get the actual minimizer for the explicit manipulations
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

        ## create the the minimizer 
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
            self.info ("minuit: Error Level is redefined from %.3f to %.4g" %  ( old_errdef ,
                                                                                 new_errdef ) )
            m.setErrorLevel ( new_errdef ) 

        return m  

    # =========================================================================
    ## perform sPlot-analysis 
    #  @code
    #  r,f = model.fitTo ( dataset )
    #  model.sPlot ( dataset ) 
    #  @endcode
    #  @attention: since internally it performs the fit, all additional fit options,
    #              like ranges, constraints, etc. must be specified!
    def sPlot ( self                            ,
                dataset                         ,
                silent       = False            ,
                proj_deps    = ()               , 
                use_weights  = True             ,
                *args       , **kwargs          ) :
        """ Make sPlot analysis
        >>> r,f = model.fitTo ( dataset )
        >>> model.sPlot ( dataset )
        
        Attention: since internally it performs the fit, all addtional fit options,
        like ranges, constraints, etc. must be specified!
        """
        assert self.alist2,\
               "PDF(%s) has empty 'alist2'/(list of components)" + \
               "no sPlot is possible" % self.name 
        
        projdeps = proj_deps if proj_deps else ROOT.RooArgSet()
        if dataset.isWeighted() and not silent : 
            if use_weights : logger.info( 'sPlot: internal weights will be used!') 
            else           : logger.info( 'sPlot: internal weights will *not* be used!')

        ## parse additional options 
        opts = self.parse_args ( dataset , *args , **kwargs )

        if opts :
            ## those can be duplicated 
            vok  = ( 'Extended' , 'SumW2Error' , 'PrintLevel' , 'PrintEvalLevel' ) 
            opts = tuple ( o for o in opts if not o.name in vok  )
                                   
        nopts = len ( opts )
        if   nopts <= 4  : copts = opts
        elif nopps <= 8  : copts = ROOT.RooFit.MultiArg ( *opts ) ,  ## note comma...
        elif nopts <= 32 :
            from ostap.utils.utils import divide
            split = divide ( 4 , opts )
            copts = tuple ( ROOT.RooFit.MultiArg  ( *[ a for a in s ] ) for s in split )
        else :
            raise TypeError( "Too many RooCmdArgs (>32)! Merge them with ROOT.RooFit.MultiArg" ) 

        if opts and not silent :
            logger.info ( 'sPlot options: %s' % str ( copts ) )  

        
        vars = set ( ( v.name for v in dataset.varset() ) )    
        with roo_silent ( True ) :

            splot = ROOT.RooStats.SPlot ( rootID ( "sPlot_" ) ,
                                          "sPlot"             ,
                                          dataset             ,
                                          self.pdf            ,
                                          self.alist2         ,
                                          projdeps            ,
                                          use_weights         ,
                                          False               , ## clone data 
                                          ''                  , ## neew name 
                                          *copts              )
        
            ## self.__splots += [ splot ]

        if not silent :
            vars = set ( ( v.name for v in dataset.varset() ) ) - vars
            if vars :  self.info ( 'sPlot: %d variables are added to dataset:\n%s' % (
                len ( vars ) ,
                dataset.table ( title     = 'Variables added by sPlot' ,
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
                   sample   = True  , ## sample number of events ?
                   silent   = True  , ## silent processing ? 
                   storage  = None  , ## storage type for dataset
                   args     = ()    ) :
        """ Generate toy-sample according to PDF
        >>> model  = ....
        >>> data   = model.generate ( 10000 ) ## generate dataset 
        
        >>> varset = ....
        >>> data   = model.generate ( 100000 , varset , sample = False )
        >>> data   = model.generate ( 100000 , varset , sample = True  )
        """

        ## sample number of events in dataset ?
        nEvents = self.gen_sample ( nEvents , sample ) 
        assert 0 <= nEvents , 'Invalid number of Events %s' % nEvents  

        args = args + ( ROOT.RooFit.Name ( dsID() ) , ROOT.RooFit.NumEvents ( nEvents ) )

        if  silent : args = args + ( ROOT.RooFit.Verbose ( False ) , )
        else       : args = args + ( ROOT.RooFit.Verbose ( True  ) , )
                                                  
        if   isinstance ( binning , integer_types ) and 0 < binning :
            args    = args + ( ROOT.RooFit.AllBinned () , ) 
        elif binning is True :
            args    = args + ( ROOT.RooFit.AllBinned () , ) 
            binning = {}
            
        if   not varset :
            varset = ROOT.RooArgSet ( self.vars )
        elif isinstance  ( varset , ROOT.RooAbsData ) :
            vs  = varset.get()
            vs2 = ROOT.RooArgSet ()
            for v in vs :
                if v in self.vars : vs2.add ( v )
            varset = vs2 
        elif isinstance ( varset , ROOT.RooAbsReal ) :
            varset = ROOT.RooArgSet ( varset    )
            
        for v in self.vars :
            if not v in varset :
                vs = ROOT.RooArgSet()
                vs . add ( v )
                for vv in varset : vs.add ( vv )
                varset = vs

        from ostap.fitting.variables import KeepBinning

        if args and not silent :
            rows = [ ( 'Option' , ) ]
            for a in args : 
                row = str ( a ) ,
                rows.append ( row )
            title = 'generate: Gen-Options'
            table = T.table ( rows , title = 'GenOptions', prefix = '# ' )
            self.info ( '%s:\n%s' % ( title , table ) )
         
        with KeepBinning ( self.xvar ) : 

            if isinstance ( binning , dict ) :
                binning = binning.get ( self.xvar.name , None )                 
            if binning : self.xvar.bins = binning
            
            if storage in ( ROOT.RooAbsData.Tree , ROOT.RooAbsData.Vector ) :
                from ostap.fitting.dataset import useStorage
                with useStorage ( storage ) : 
                    return self.pdf.generate (  varset , *args )

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
                nbins    = 100   , 
                hpars    = ()    , 
                histo    = None  ,
                integral = True  ,
                events   = True  , 
                errors   = False , **kwargs ) :
        """ Convert PDF to the 1D-histogram in correct way
        - Unlike  `PDF.roo_histo` method, PDF is integrated within the bin
        >>> pdf = ...
        >>> h1  = pdf.histo ( 100 , 0. , 10. ) ## specify histogram parameters
        >>> histo_template = ...
        >>> h2  = pdf.histo ( histo = histo_template ) ## use histogram template
        >>> h3  = pdf.histo ( ... , integral = True  ) ## use PDF integral within the bin  
        >>> h4  = pdf.histo ( ... , density  = True  ) ## convert to 'density' histogram 
        """
        
        histo = self.make_histo ( nbins = nbins ,
                                  hpars = hpars ,
                                  histo = histo , **kwargs )

        # loop over the histogram bins 
        for i , x , y in histo.items() :

            xv , xe = x.value() , x.error()
            
            # value at the bin center 
            c = self ( xv , error = errors ) 

            if not integral : 
                histo [ i ] = c
                continue

            # integral over the bin 
            v       = self.integral( xv - xe , xv + xe , nevents = events )

            # scale it by the bin volume 
            volume  = 2 * x.error() 
            v      /= volume
            
            if errors :
                if    0 == c.cov2 () : pass
                elif  0 != c.value() and 0 != v : 
                    v = c * ( v / c.value() )
                    
                    histo [ i ] = v
                    
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
                    nbins    = 100   , 
                    hpars    = ()    , 
                    histo    = None  ,
                    events   = True  , **kwargs ) : 
        """ Convert PDF to the 1D-histogram, taking PDF-values at bin-centres
        - see RooAbsPdf::createHistogram
        - see RooAbsPdf::fillHistogram
        - see PDF.histo
        >>> pdf = ...
        >>> h1  = pdf.roo_histo ( 100 , 0. , 10. ) ## specify histogram parameters
        >>> histo_template = ...
        >>> h2  = pdf.roo_histo ( histo = histo_template ) ## use histogram template
        """
        
        histo = self.make_histo ( nbins = nbins ,
                                  hpars = hpars ,
                                  histo = histo , **kwargs )
        
        with rootException() , warnings.catch_warnings() :
            warnings.simplefilter ( 'ignore' , RuntimeWarning ) 
            hh = self.pdf.createHistogram (
                hID ()    ,
                self.xvar ,
                self.binning ( histo.GetXaxis() , 'histo1x' ) ,
                ROOT.RooFit.Extended ( False ) ,
                ROOT.RooFit.Scaling  ( False ) ,            
            )
            
        ## nullify errors 
        for i in hh : hh.SetBinError ( i , 0 ) 
    
        if events and self.pdf.mustBeExtended() :            
            for i , x , y in hh.items() :
                volume  = 2*x.error() 
                hh [ i ]  *= volume
                
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
        """ Create the residual histogram   (data - fit)
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
        """ Create the pull histogram   (data - fit)/data_error
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
        """ Get the residual histogram
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
        """ Get the pull  histogram: (data-fit)/data_error
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
    ## make 2D-contours
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
        """ Create popular 1D 'background'  function.
        
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
    ## Check the ranges for variables in dataset:
    #
    #  @param dataset dataset to check
    #  @return True if dataset has no  entries outside the fitting range
    #
    #  @code
    #  data = ...
    #  pdf  = ...
    #  ok   = pdf.check_ranges ( data ) 
    #  @endcode
    # 
    def check_ranges ( self , dataset , cut_range = '' ) :
        """ Check the ranges for variables in dataset:  
        - return True if dataset has no entries outside the fitting range  

        >>> data = ...
        >>> pdf  = ...
        >>> ok   = pdf.check_ranges ( data ) 
        """

        import ostap.trees.cuts

        cuts = ROOT.TCut() 
        for v in self.vars :
            
            ## check that variable in the dataset 
            if dataset and not v in dataset : continue
            
            if hasattr ( v , 'hasMin' ) and hasattr ( v , 'getMin' ) and v.hasMin ( cut_range ) :
                vmin  = v.getMin ( cut_range )
                if isfinite ( vmin ) : cuts |= ( ROOT.TCut ( v.name ) < vmin ) 
                
            if hasattr ( v , 'hasMax' ) and hasattr ( v , 'getMax' ) and v.hasMax ( cut_range ) :
                vmax  = v.getMax ( cut_range )
                if isfinite ( vmax ) : cuts |= ( ROOT.TCut ( v.name ) > vmax )

        has_entry = dataset.hasEntry ( cuts = cuts , cut_range = cut_range ) 
        return not has_entry 

    # =========================================================================
    ## Make PDF1 object
    #  @code
    #  pdf = ...
    #  pdf2 , xvar = pdf.make_PDF( ... , xvar = .. , )
    #  @endcode
    def make_PDF1 ( self , pdf , xvar = None , *args , **kwargs ) :
        """ Make PDF1 object
        >>> pdf = ...
        >>> pdf2 , xvar = pdf.make_PDF( ... , xvar = .. , )        
        """
        if   isinstance ( pdf , PDF1 ) :
            
            assert ( not xvar ) or ( xvar in pdf.vars ) , \
                       "make_PDF1: Invalid setting of xvar %s vs %s" % ( xvar , pdf.vars )
                
            return pdf , pdf.xvar 
        
        elif xvar and isinstance ( pdf , ROOT.RooAbsPdf ) :
            
            return Generic1D_pdf ( pdf , xvar = xvar , *args , **kwargs ) , xvar 
        
        elif xvar and isinstance ( pdf , bkg_types ) :

            name = kwargs.pop ( 'name' , self.generate_name ( 'BB' ) )
            return self.make_bkg ( pdf , name = name , xvar = xvar , **kwargs ) , xvar 

        raise TypeError( "make_PDF1: invalid pdf/xvar %s/%s" % ( pdf , xvar ) )

    # =========================================================================
    ## helper functon to make a raw product of PDFs or RooAbsPDF objects
    def raw_product ( self , *pdfs ) :
        """ Make a raw product of PDFs or RooAbsPDF objects
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

    # ========================================================================
    ## create constrained PDF
    #  @see Constrained
    #  @see Constrained1D    
    def make_constrained ( self , *constraints ) :
        """ Create constrained PDF
        - see Constrained
        - see Constrained1D
        """
        return Constrained1D ( self , *constraints )

    # ========================================================================== 
    ## List/tuple of structural components from `self.alist1` 
    def cmp_alist ( self )  :
        """ List/tuple of structural components from `self.alist1`"""
        return tuple ( Generic1D_pdf ( p , xvar = self.xvar ) for p in self.alist1 )
        
# =============================================================================
## @class PDF1
#  The main helper base class for implementation of various 1D PDF-wrappers 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-08-21
class PDF1(APDF1,FUN1) :
    """ The main helper base class for implementation of various 1D PDF-wrappers
    """
    def __init__ ( self , name ,  xvar , tricks = True , **kwargs ) :

        ## initialize the base 
        FUN1  .__init__ ( self , name = name , xvar = xvar  , tricks = tricks ,
                          fun = kwargs.pop ( 'pdf' , None ) , **kwargs )

        ## initialize the base 
        APDF1 .__init__ ( self )
        
        self.config   = { 'name'    : self.name    ,
                          'xvar'    : self.xvar    ,
                          'tricks'  : self.tricks  ,
                          'pdf'     : self.pdf     }
        self.config.update ( kwargs )

        self.__call_OK = isinstance ( self.xvar , ROOT.RooAbsRealLValue ) 

    # =========================================================================
    ## simple 'function-like' interface
    @vct1_call_method
    def __call__ ( self , x , error = False , normalized = True  ) :
        """ Function as a 'function'
        >>> fun  = ...
        >>> x = 1
        >>> y = fun ( x ) 
        """
        assert self.__call_OK , "Invalid type for xvar!"
        
        if error and not normalized :
            self.error ( "Can't get error for non-normalized call" )
            error = False
            
        if error and not self.fit_result :
            error = False 
            
        ## min-max, if defined 
        xmnmx = self.xminmax()

        if xmnmx :    
            xmn , xmx = xmnmx
            if not xmn <= x <= xmx : return 0

        ## ensure the value does not change after the call 
        with SETVAR ( self.xvar ) :
            
            ## change x-value
            self.xvar.setVal ( x )
            
            ## evaluate the function
            v = self.fun.getVal ( self.vars ) if normalized else self.fun.getVal ( ROOT.nullptr )
            
            ## get ucertainties if/when available 
            if error and self.fit_result :
                e = self.pdf.getPropagatedError ( self.fit_result )
                if 0 <= e : v = VE ( v ,  e * e )
                
        return v
        
    # ========================================================================
    ## convert to float 
    def __float__ ( self ) :
        """ Convert to float
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
        """ Make a product of two PDFs
        
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
        """ Make a product of two PDFs
        
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
    ## make a non-extended sum of 1D PDFs
    #  @code
    #  pdf1 = ...
    #  pdf2 = ...
    #  pdf  = pdf1 + pdf2     
    #  @endcode
    def __add__ ( self , other ) :
        """ Make a nom-extended sum of 1D PDFs
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
        """ Make a non-extended sum of 1D PDFs
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
        >>> ## result = pdf @ other  ## python 3 only 
        
        `Other` can be
        - fully configured Convolution object

        It also can be 
        - resolution  PDF   , it will be treated as resolution function 
        - `RooAbsPdf`       , it will be treated as resolution function
        - `RooAbsReal`      , it will be treated as sigma for Gaussian resolution function
        - positive constant , it will be treated as sigma for Gaussian resolution function
        - 2 or 3-tuple: it will be treated as sigma for Gaussian resolution function
        
        The configuration can be specified via `ConvolutionConfig`
        context manager:
        >>> pdf    = ...
        >>> other  =  ...
        >>> with ConvolutionConfig ( buffer = 0.25 , nbins = 1000 ):  
        >>> ... result = pdf % other    ## python 2 an d3 
        >>> ... ## result = pdf @ other ## python 3 only 
        """
        from ostap.fitting.pdf_ops import pdf_convolution
        return pdf_convolution ( self , other )

    __matmul__  = __mod__

    # =========================================================================
    ## Convert PDF into simple function
    #  @code
    #  pdf = ...
    #  fun = pdf.as_FUN () 
    #  @endcode
    def as_FUN ( self , name = '' ) : 
        """ Convert PDF into simple function
        >>> pdf = ...
        >>> fun = pdf.as_FUN ()
        """
        return Fun1D ( self.pdf  ,
                       xvar = self.xvar  ,
                       name = name if name else self.new_name ( 'fun1' ) ) 
    
    # =========================================================================
    ## Get the 1D-CDF from 1D-PDF
    #  @see RooAbdPdf::createCdf
    #  @code
    #  pdf = ...
    #  cdf = pdf.cdf () 
    #  @endcode
    def cdf ( self ) :
        """ Get the 1D-CDF from 1D-PDF
        - see `ROOT.RooAbdPdf.createCdf`
        >>> pdf = ...
        >>> cdf = pdf.cdf () 
        """
        assert 1 == len ( self.vars ) , 'One can make CDF only for 1D-PDF!'
        fun  = self.pdf.createCdf ( self.vars  )
        name = self.generate_name ( prefix = 'cdf_' , name = self.name )  
        return Fun1D ( fun , xvar = self.xvar , name = name )
    
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
    """ Wrapper for generic RooFit pdf
    >>> raw_pdf = RooGaussian   ( ...     )
    >>> pdf     = Generic1D_pdf ( raw_pdf , xvar = x )
    """
    ## constructor 
    def __init__ ( self           , pdf   ,
                   xvar           = None  ,
                   name           = ''    ,
                   add_to_signals = True  ,
                   prefix         = ''    ,
                   suffix         = ''    ) :
        """ Wrapper for generic RooFit pdf        
        >>> raw_pdf = RooGaussian   ( ...     )
        >>> pdf     = Generic1D_pdf ( raw_pdf , xvar = x )
        """
        assert xvar and   isinstance ( xvar  , ROOT.RooAbsReal ) , "'xvar' must be ROOT.RooAbsReal"
        assert pdf  and ( isinstance ( pdf   , ROOT.RooAbsPdf  ) or \
                          ( isinstance ( pdf , ROOT.RooAbsReal ) ) ) , \
                          "Invalid `pdf' type"
        
        name = ( prefix + name + suffix ) if name \
            else self.generate_name ( prefix = prefix , suffix = suffix , name = pdf.GetName() )
        
        ## initialize the base 
        PDF1 . __init__ ( self , name , xvar )
        ##

        ## Does PDF depends on XVAR ?
        if not pdf.depends_on ( self.xvar ) :
            self.warning ( "PDF/%s does not depend on %s!" % ( pdf.name , self.xvar.name ) ) 

        ## PDF itself 
        self.pdf  = pdf

        ## get some structure 
        if isinstance ( self.pdf , ROOT.RooAddPdf ) :
            for p in self.pdf.pdfList    ()    : self.alist1.add ( p )
            for f in self.pdf.orig_fracs ()[0] : self.alist2.add ( f )
            
        if isinstance ( self.xvar , ROOT.RooAbsRealLValue ) and not self.pdf.dependsOn ( self.xvar ) : 
            self.warning ("Function/PDF does not depend on xvar=%s" % self.xvar.name )

        ## add it to the list of signal components ?
        self.__add_to_signals = True if add_to_signals else False
        
        if self.add_to_signals : self.signals.add ( self.pdf )
        
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
    """ Helper function to create the PDF/PDF2/PDF3
    """
    
    assert pdf and isinstance ( pdf , ROOT.RooAbsReal ), \
           'make_pdf: Invalid type %s' % type ( pdf )
    
    name = name if name else "PDF_from_%s" % pdf.name
    
    if not isinstance ( pdf , ROOT.RooAbsPdf ) :
        pdf = ROOT.RooWrapperPdf        ( name , 'PDF from %s' % pdf.name , pdf )
        
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
    #  r , f = model.fitTo ( dataset )
    #  r , f = model.fitTo ( dataset , weighted = True )    
    #  r , f = model.fitTo ( dataset , ncpu     = 10   )    
    #  r , f = model.fitTo ( dataset , draw = True , nbins = 300 )    
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
        """ Perform the actual fit (and draw it)
        >>> r , f = model.fitTo ( dataset )
        >>> r , f = model.fitTo ( dataset , weighted = True )    
        >>> r , f = model.fitTo ( dataset , ncpu     = 10   )    
        >>> r , f = model.fitTo ( dataset , draw = True , nbins = 300 )    
        """
        if   isinstance ( dataset , H2D_dset ) : dataset = dataset.dset        
        elif isinstance ( dataset , ROOT.TH2 ) :
            assert 2 == dataset.GetDimension() , 'Invalid histogram dimension: %s' % dataset.GetDimension() 
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
                if dataset :
                    dsyvar = dataset.get_var ( self.yvar.name ) 
                    if dsyvar : dsyvar.setRange ( range_name , in_range[0] , in_range[1] )
                
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
        """ Draw the projection over 2nd variable
        
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
                if dataset :
                    dsxvar = dataset.get_var ( self.xvar.name )
                    if dsxvar : dsxvar.setRange ( range_name , in_range[0] , in_range[1] )
                    
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
        """ Make/draw 2D-histograms 
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
        if dataset : dataset.fillHistogram ( hdata , _lst ) 
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
        """ Make 1D-plot:
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
            options_cut = ROOT.RooFit.CutRange ( in_range ) , 
            newargs [ 'data_options' ] = self.draw_option ( 'data_options' , **newargs ) + options_cut
            
        if in_range : 
            options_project =  ROOT.RooFit.ProjectionRange ( in_range ) , 
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
    ## fit the 2D-histogram
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
        """ Fit the 2D-histogram        
        >>> histo = ...
        >>> r,f = model.fitHisto ( histo )        
        """

        xminmax = histo.xminmax()
        yminmax = histo.yminmax()        
        with RangeVar ( self.xvar , *xminmax ) , RangeVar ( self.yvar , *yminmax ):
            
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
                self.histo_data = H2D_dset ( histo , self.xvar , self.yvar , density = density , silent = silent )
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
                   silent   = True  , ## silent processins?
                   storage  = None  ,  
                   args     = ()    ) :
        """ Generate toy-sample according to PDF
        >>> model  = ....
        >>> data   = model.generate ( 10000 ) ## generate dataset
        
        >>> varset = ....
        >>> data   = model.generate ( 100000 , varset , sample = False )
        >>> data   = model.generate ( 100000 , varset , sample = True  )
        """
        nEvents = self.gen_sample ( nEvents , sample ) 
        assert 0 <= nEvents , 'Invalid number of Events %s' % nEvents  
        
        args = args + ( ROOT.RooFit.Name ( dsID() ) , ROOT.RooFit.NumEvents ( nEvents ) )
        
        if  silent : args = args + ( ROOT.RooFit.Verbose ( False ) , )
        else       : args = args + ( ROOT.RooFit.Verbose ( True  ) , )
                                                  
        if binning is True :
            args    = args + ( ROOT.AllBinned() , ) 
            binning = {}

        if   not varset :
            varset = ROOT.RooArgSet( self.xvar , self.yvar )
        elif isinstance  ( varset , ROOT.RooAbsData ) :
            vs  = varset.get()
            vs2 = ROOT.RooArgSet ()
            for v in vs :
                if v in self.vars : vs2.add ( v )
            varset = vs2             
        elif isinstance ( varset , ROOT.RooAbsReal ) :
            varset = ROOT.RooArgSet( varset )

        for v in self.vars :
            if not v in varset :
                vs = ROOT.RooArgSet()
                vs . add ( v )
                for vv in varset : vs.add ( vv )
                varset = vs

        if args and not silent :
            rows = [ ( 'Option' , ) ]
            for a in args : 
                row = str ( a ) ,
                rows.append ( row )
            title = 'generate: Gen-Options'
            table = T.table ( rows , title = 'GenOptions', prefix = '# ' )
            self.info ( '%s:\n%s' % ( title , table ) )
    
        from ostap.fitting.variables import KeepBinning        
        with KeepBinning ( self.xvar ) , KeepBinning ( self.yvar ) : 

            if binning :
                
                xbins = binning.get ( self.xvar.name , None )
                ybins = binning.get ( self.yvar.name , None )

                if xbins : self.xvar.bins = xbins
                if ybins : self.yvar.bins = ybins
                
            if storage in ( ROOT.RooAbsData.Tree , ROOT.RooAbsData.Vector ) :
                from ostap.fitting.dataset import useStorage
                with useStorage ( storage ) : 
                    return self.pdf.generate (  varset , *args )

            return self.pdf.generate ( varset , *args )

    # ========================================================================
    ## check minmax of the PDF using the random shoots
    #  @code
    #  pdf     = ....
    #  mn , mx = pdf.minmax()            
    #  @endcode 
    def minmax ( self , nshoots =  100000 ) :
        """ Check min/max for the PDF using  random shoots 
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
        """ Get integral over (xmin,xmax,ymin,ymax) region
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
        """ Get a minimum of PDF for certain interval
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
        """ Get a maximum of PDF for certain interval
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
        """ Convert PDF to TF2 object, e.g. to profit from TF2::Draw options
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
                     xbins    = 20    , 
                     ybins    = 20    ,  
                     hpars    = ()    , 
                     histo    = None  , **kwargs ) :
        """ Create the 2D histogram according to specifications
        """
        
        from ostap.histos.histos import histo_book 

        # histogram is provided 
        if histo :
            
            assert isinstance ( histo , ROOT.TH2 ) and 2 == histo.GetDimension() , \
                "Illegal type of 'histo'-argument %s" % type( histo )
            
            histo = histo.clone()
            histo.Reset()

        # arguments for the histogram constructor 
        elif hpars :
            
            histo = ROOT.TH2d ( hID () , 'PDF%s' % self.name , *hpars  )
            if not histo.GetSumw2() : histo.Sumw2()

        # explicit construction from (#bins,min,max)-triplet  
        else :
            
            ranges = [ ( self.xvar.name , self.xminmax() ) ,
                       ( self.yvar.name , self.yminmax() ) ] 
            histo  = histo_book ( ranges , xbins = xbins , ybins = ybins , **kwargs )

        histo.SetDirectory ( ROOT.nullptr ) 
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
                events   = True  , 
                density  = False ) :
        """ Convert PDF to the 2D-histogram
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
        for ix , iy , x , y , z in histo.items() :

            xv , xe = x.value() , x.error()
            yv , ye = y.value() , y.error()
            
            # value at the bin center 
            c = self ( xv , yv , error = errors ) 

            if not integral : 
                histo[ix,iy] = c
                continue

            # integral over the bin 
            v  = self.integral( xv - xe , xv + xe , yv - ye , yv + ye , nevents = events )

            # scale it by the bin volume 
            volume  = 4 * x.error() * y.error() 
            v      /= volume

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
    def roo_histo ( self           ,
                    xbins   = 20    , 
                    ybins   = 20    , 
                    hpars   = ()    , 
                    histo   = None  , 
                    events  = True  , **kwargs ) : 
        """ Convert PDF to the 2D-histogram, taking PDF-values at bin-centres
        >>> pdf = ...
        >>> h1  = pdf.roo_histo ( 100 , 0. , 10. , 20 , 0. , 10 ) 
        >>> histo_template = ...
        >>> h2  = pdf.roo_histo ( histo = histo_template ) ## use histogram template
        >>> h3  = pdf.roo_histo ( ... , density  = True  ) ## convert to 'density' histogram 
        """
        
        histo = self.make_histo ( xbins = xbins , 
                                  ybins = ybins , 
                                  hpars = hpars ,
                                  histo = histo , **kwargs )
        
        with rootException() , warnings.catch_warnings() : 
            warnings.simplefilter ( 'ignore' , RuntimeWarning )         
            hh = self.pdf.createHistogram (
                hID()     ,
                self.xvar ,                    self.binning ( histo.GetXaxis() , 'histo2x' )   ,
                ROOT.RooFit.YVar ( self.yvar , self.binning ( histo.GetYaxis() , 'histo2y' ) ) , 
                ROOT.RooFit.Scaling  ( False ) , 
                ROOT.RooFit.Extended ( False ) ) 

        ## nullify errors 
        for i , j in hh : hh.SetBinError ( i , j , 0 ) 
        
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
        """ Get the residual histogram
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
        """ Get the pull  histogram: (data-fit)/data_error
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
        return '%s(%s,xvar=%s,yvar=%s)' % ( typename ( self ) ,
                                            self.name         ,
                                            self.xvar.name    ,
                                            self.yvar.name    )
    __repr__ = __str__ 

    # ===========================================================================
    ## Make PDF2 object 
    def make_PDF2 ( self , pdf , xvar = None , yvar = None , *args , **kwargs ) :
        """ Make PDF2 object
        """
        if   isinstance  ( pdf , PDF2 ) :
            
            assert ( not xvar ) or ( xvar in pdf.vars ) , \
                   "make_PDF2: Invalid setting of xvar %s vs %s" % ( xvar , pdf.vars )
            assert ( not yvar ) or ( yvar in pdf.vars ) , \
                   "make_PDF2: Invalid setting of yvar %s vs %s" % ( yvar , pdf.vars )
   
            return pdf, pdf.xvar, pdf.yvar
        
        elif isinstance ( pdf , ROOT.RooAbsPdf   ) and xvar and yvar :
            
            return Generic2D_pdf ( pdf , xvar = xvar , yvar = yvar , *args , **kwargs ) , xvar , yvar

        raise TypeError( "make_PDF2: invalid pdf/xvar %s/%s" % ( pdf , xvar ) )

    # ========================================================================
    ## create constrained PDF
    #  @see Constrained
    #  @see Constrained1D    
    def make_constrained ( self , *constraints ) :
        """ Create constrained PDF
        - see Constrained
        - see Constrained1D
        """
        return Constrained2D ( self , *constraints )
    
    # ========================================================================== 
    ## List/tuple of structural components from `self.alist1` 
    def cmp_alist ( self )  :
        """ List/tuple of structural components from `self.alist1`"""
        return tuple ( Generic2D_pdf ( p ,
                                       xvar = self.xvar ,
                                       yvar = self.yvar ) for p in self.alist1 )
                
# =============================================================================
## @class PDF2
#  The main helper base class for implementation of various 1D PDF-wrappers 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-08-21
class PDF2(APDF2,FUN2) :
    """ The main helper base class for implementation of various 1D PDF-wrappers
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
    @vct2_call_method
    def __call__ ( self , x , y , error = False , normalized = True ) :
        """  Simple  function-like interface
        >>>  pdf = ...
        >>>  print ( pdf(0.1,0.5) ) 
        """
        assert self.__call_OK , "Invalid types for xvar/yvar!"

        if error and not normalized :
            self.error ( "Can't get error for non-normalized call" )
            error = False
            
        if error and not self.fit_result :
            error = False 
       
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
            
            v = self.pdf.getVal ( self.vars ) if normalized else self.pdf.getVal ( ROOT.nullptr )
            
            if error and self.fit_result :
                e = self.pdf.getPropagatedError ( self.fit_result )
                if 0 <= e : v = VE ( v , e * e )

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
        """ Make a product of two PDFs
        
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
        """ Make a product of two PDFs
        
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
        """ Make a no-extended sum of 1D PDFs
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

        ## Does PDF depends on XVAR ?
        if not pdf.depends_on ( self.xvar ) :
            self.warning ( "PDF/%s does not depend on %s!" % ( pdf.name , self.xvar.name ) ) 
        ## Does PDF depends on YVAR ?
        if not pdf.depends_on ( self.yvar ) :
            self.warning ( "PDF/%s does not depend on %s!" % ( pdf.name , self.xvar.name ) ) 
            
        ## PDF! 
        self.pdf = pdf

        if isinstance ( self.xvar , ROOT.RooAbsRealLValue ) and not self.pdf.dependsOn ( self.xvar ) : 
            self.warning ("Function/PDF does not depend on xvar=%s" % self.xvar.name )
        if isinstance ( self.yvar , ROOT.RooAbsRealLValue ) and not self.pdf.dependsOn ( self.yvar ) : 
            self.warning ("Function/PDF does not depend on yvar=%s" % self.yvar.name )

        ## get some structure 
        if isinstance ( self.pdf , ROOT.RooAddPdf ) :
            for p in self.pdf.pdfList    ()       : self.alist1.add ( p )
            for f in self.pdf.orig_fracs () [ 0 ] : self.alist2.add ( f )

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
    #  r , f = model.fitTo ( dataset )
    #  r , f = model.fitTo ( dataset , weighted = True )    
    #  r , f = model.fitTo ( dataset , ncpu     = 10   )    
    #  r , f = model.fitTo ( dataset )    
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
            assert 3 == dataset.GetDimension () , 'Invalid histogram dimension: %s' % dataset.GetDimension() 
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

        return result , None 
    
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
        """ Draw the projection over 1st variable
        
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
        in_range = None

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
        """ Fit the 3D histogram
        
        >>> histo = ...
        >>> r,f = model.fitHisto ( histo )
        
        """
        
        xminmax = histo.xminmax()
        yminmax = histo.yminmax()
        zminmax = histo.zminmax()
        
        with     RangeVar ( self.xvar , *xminmax ) , \
                 RangeVar ( self.yvar , *yminmax ) , \
                 RangeVar ( self.xvar , *zminmax ) : 

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
                self.histo_data = H3D_dset ( histo     ,
                                             self.xvar ,
                                             self.yvar ,
                                             self.zvar ,
                                             density = density ,
                                             silent  = silent   )
                data = self.histo_data
                
            if chi2 : return self.chi2fitTo ( data              ,
                                              draw    = draw    ,
                                              silent  = False   ,
                                              density = density ,
                                              args    = args    , **kwargs )
            else    : return self.fitTo     ( data              ,
                                              silent  = silent  ,
                                              args    = args    , **kwargs ) 

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
                   silent   = True  , ## silent processing 
                   storage  = None  ,
                   args     = ()    ) :
        """ Generate toy-sample according to PDF
        >>> model  = ....
        >>> data   = model.generate ( 10000 ) ## generate dataset
        
        >>> varset = ....
        >>> data   = model.generate ( 100000 , varset , sample = False )
        >>> data   = model.generate ( 100000 , varset , sample = True  )
        """
        nEvents = self.gen_sample ( nEvents , sample ) 
        assert 0 <= nEvents , 'Invalid number of Events %s' % nEvents  

        args = args + ( ROOT.RooFit.Name ( dsID() ) , ROOT.RooFit.NumEvents ( nEvents ) )
        
        if  silent : args = args + ( ROOT.RooFit.Verbose ( False ) , )
        else       : args = args + ( ROOT.RooFit.Verbose ( True  ) , )
                                                      
        if binning is True :
            args    = args + ( ROOT.AllBinned() , ) 
            binning = {}

        if   not varset :
            varset = ROOT.RooArgSet( self.xvar , self.yvar , self.zvar )
        elif isinstance  ( varset , ROOT.RooAbsData ) :
            vs  = varset.get()
            vs2 = ROOT.RooArgSet ()
            for v in vs :
                if v in self.vars : vs2.add ( v )
            varset = vs2                         
        elif isinstance ( varset , ROOT.RooAbsReal ) :
            varset = ROOT.RooArgSet( varset )

        for v in self.vars :
            if not v in varset :
                vs = ROOT.RooArgSet()
                vs . add ( v )
                for vv in varset : vs.add ( vv )
                varset = vs
                
        if args and not silent :
            rows = [ ( 'Option' , ) ]
            for a in args : 
                row = str ( a ) ,
                rows.append ( row )
            title = 'generate: Gen-Options'
            table = T.table ( rows , title = 'GenOptions', prefix = '# ' )
            self.info ( '%s:\n%s' % ( title , table ) )
        
        from ostap.fitting.variables import KeepBinning        
        with KeepBinning ( self.xvar ) , KeepBinning ( self.yvar ), KeepBinning ( self.zvar ) : 
            
            if binning :
                
                xbins = binning.get ( self.xvar.name , None )
                ybins = binning.get ( self.yvar.name , None )
                zbins = binning.get ( self.zvar.name , None )
                
                if xbins : self.xvar.bins = xbins
                if ybins : self.yvar.bins = ybins
                if zbins : self.zvar.bins = zbins

            if storage in ( ROOT.RooAbsData.Tree , ROOT.RooAbsData.Vector ) :
                from ostap.fitting.dataset import useStorage
                with useStorage ( storage ) : 
                    return self.pdf.generate (  varset , *args )

            return self.pdf.generate ( varset , *args )

    # ========================================================================
    ## check minmax of the PDF using the random shoots
    #  @code
    #  pdf     = ....
    #  mn , mx = pdf.minmax()            
    #  @endcode 
    def minmax ( self , nshoots = 200000 ) :
        """ Check min/max for the PDF using  random shoots 
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
        """ Get integral over (xmin,xmax,ymin,ymax,zmin,zmax) region
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
        """ Get a minimum of PDF for certain interval
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
        """ Convert PDF  to TF3 object, e.g. to profit from TF3::Draw options
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
                     xbins    = 10   , 
                     ybins    = 10   , 
                     zbins    = 10   , 
                     hpars    = ()   , 
                     histo    = None , **kwargs ) :
        """ Create the histogram accoring to specifications"""
        
        from ostap.histos.histos import histo_book 

        # histogram is provided 
        if histo :
            
            assert isinstance ( histo , ROOT.TH3 ) and 3 == histo.GetDimension() , \
                "Illegal type of 'histo'-argument %s" % typename ( histo )
            
            histo = histo.clone()
            histo.Reset()

        # arguments for the histogram constructor 
        elif hpars :
            
            from ostap.core.core import hID
            histo = ROOT.TH3F ( hID () , 'PDF%s' % self.name , *hpars  )
            if not histo.GetSumw2() : histo.Sumw2()

        # explicit contruction from (#bins,min,max)-triplet  
        else :

            ranges = [ ( self.xvar.name , self.xminmax() ) ,
                       ( self.yvar.name , self.yminmax() ) , 
                       ( self.zvar.name , self.zminmax() ) ] 
            histo  = histo_book ( ranges        ,
                                  xbins = xbins ,
                                  ybins = ybins ,
                                  zbins = xbins , **kwargs )

        histo.SetDirecotry ( ROOT.nullptr )
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
                events   = True  , 
                errors   = False ) :
        """ Convert PDF to the 3D-histogram in correct way
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
                                zv - ze , zv + ze , nevents = events )

            # scale it by the bin volume 
            volume  = 8 * x.error() * y.error() * z.error() 
            v      /= volume

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
    def roo_histo ( self             ,
                    xbins    = 10    , 
                    ybins    = 10    , 
                    zbins    = 10    , 
                    hpars    = ()    , 
                    histo    = None  ,
                    events   = True  , **kwargs ) :
        """ Convert PDF to the 3D-histogram, taking PDF-values at bin-centres
        >>> pdf = ...
        >>> h1  = pdf.roo_histo ( 10 , 0. , 10. , 10 , 0. , 4. , 10 , 0. , 3 )
        >>> histo_template = ...
        >>> h2  = pdf.roo_histo ( histo = histo_template ) ## use histogram template
        """
        
        histo = self.make_histo ( xbins = xbins , 
                                  ybins = ybins , 
                                  zbins = zbins , 
                                  hpars = hpars ,
                                  histo = histo , **kwargs )
        
        with rootException() , warnings.catch_warnings() : 
            warnings.simplefilter ( 'ignore' , RuntimeWarning )         
            hh = self.pdf.createHistogram (
                hID()     ,
                self.xvar ,                    self.binning ( histo.GetXaxis() , 'histo3x' )   ,
                ROOT.RooFit.YVar ( self.yvar , self.binning ( histo.GetYaxis() , 'histo3y' ) ) , 
                ROOT.RooFit.ZVar ( self.zvar , self.binning ( histo.GetZaxis() , 'histo3z' ) ) , 
                ROOT.RooFit.Scaling  ( False ) , 
                ROOT.RooFit.Extended ( False ) )
        
        for i,j,k in hh : hh.SetBinError ( i , j , k , 0 ) 
        
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
        """ Get the residual histogram
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
        """ Get the pull  histogram: (data-fit)/data_error
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
        return '%s(%s,xvar=%s,yvar=%s,zvar=%s)' % ( typename ( self ) ,
                                                    self.name         ,
                                                    self.xvar.name    ,
                                                    self.yvar.name    ,
                                                    self.zvar.name    )
    __repr__ = __str__ 

    # =========================================================================
    ## Make PDF3 object 
    def make_PDF3 ( self , pdf , xvar = None , yvar = None , zvar = None , *args , **kwargs ) :
        """ Make PDF1 object
        """
        if isinstance ( pdf , PDF3 ) :

            assert ( not xvar ) or ( xvar in pdf.vars ) , \
                   "make_PDF3: Invalid setting of xvar %s vs %s" % ( xvar , pdf.vars )
            assert ( not yvar ) or ( yvar in pdf.vars ) , \
                   "make_PDF3: Invalid setting of xvar %s vs %s" % ( yvar , pdf.vars )
            assert ( not zvar ) or ( zvar in pdf.vars ) , \
                   "make_PDF3: Invalid setting of xvar %s vs %s" % ( xvar , pdf.vars )
          
            return pdf, pdf.xvar, pdf.yvar, pdf.zvar 
        
        elif isinstance ( pdf , ROOT.RooAbsPdf ) and xvar and yvar and zvar :
            
            return Generic3D_pdf ( pdf , xvar = xvar , yvar = yvar , zvar = zvar , *args , **kwargs ) , xvar , yvar , zvar 

        raise TypeError( "make_PDF3: invalid pdf/xvar %s/%s" % ( pdf , xvar ) )

    # ========================================================================
    ## create constrained PDF
    #  @see Constrained
    #  @see Constrained1D    
    def make_constrained ( self , *constraints ) :
        """ Create constrained PDF
        - see Constrained
        - see Constrained1D
        """
        return Constrained3D ( self , *constraints )
    
    # ========================================================================== 
    ## List/tuple of structural components from `self.alist1` 
    def cmp_alist ( self )  :
        """ List/tuple of structural components from `self.alist1`"""
        return tuple ( Generic3D_pdf ( p ,
                                       xvar = self.xvar ,
                                       yvar = self.yvar ,
                                       zvar = self.zvar ) for p in self.alist1 )
        
# =============================================================================
## @class PDF3
#  The main helper base class for implementation of various 3D PDF-wrappers 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-08-21
class PDF3(APDF3,FUN3) :
    """ The main helper base class for implementation of various 1D PDF-wrappers
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
    @vct3_call_method
    def __call__ ( self , x , y , z , error = False , normalized = True ) :
        """ Simple function-like interface, converting PDF3 to callable 
        """
        assert self.__call_OK , "Invalid types for xvar/yvar/zvar!"
        
        if error and not normalized :
            self.error ( "Can't get error for non-normalized call" )
            error = False
            
        if error and not self.fit_result :
            error = False 
            
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
            
            v = self.pdf.getVal ( self.vars ) if normalized else self.pdf.getVal ( ROOT.nullpr )
            
            if error and self.fit_result :
                e = self.pdf.getPropagatedError ( self.fit_result )
                if 0 <= e : v = VE ( v , e * e )
                
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
        """ Make a product of two PDFs
        
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
        """ Make a product of two PDFs
        
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
    ## make a non-extended sum of 3D PDFs
    #  @code
    #  pdf1 = ...
    #  pdf2 = ...
    #  pdf  = pdf1 + pdf2     
    #  @endcode
    def __add__ ( self , other ) :
        """ Make a non-extended sum of 3D PDFs
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
        """ Make a no-extended sum of 3D PDFs
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
        """ Convert PDF into simple function
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

        ## Does PDF depend on XVAR ?
        if not pdf.depends_on ( self.xvar ) :
            self.warning ( "PDF/%s does not depend on %s!" % ( pdf.name , self.xvar.name ) ) 
        ## Does PDF depend on YVAR ?
        if not pdf.depends_on ( self.yvar ) :
            self.warning ( "PDF/%s does not depend on %s!" % ( pdf.name , self.xvar.name ) ) 
        ## Does PDF depend on ZVAR ?
        if not pdf.depends_on ( self.zvar ) :
            self.warning ( "PDF/%s does not depend on %s!" % ( pdf.name , self.zvar.name ) ) 

        ## PDF! 
        self.pdf = pdf
        
        if isinstance ( self.xvar , ROOT.RooAbsRealLValue ) and not self.pdf.dependsOn ( self.xvar ) : 
            self.warning ("Function/PDF does not depend on xvar=%s" % self.xvar.name )
        if isinstance ( self.yvar , ROOT.RooAbsRealLValue ) and not self.pdf.dependsOn ( self.yvar ) : 
            self.warning ("Function/PDF does not depend on yvar=%s" % self.yvar.name )
        if isinstance ( self.zvar , ROOT.RooAbsRealLValue ) and not self.pdf.dependsOn ( self.zvar ) : 
            self.warning ("Function/PDF does not depend on zvar=%s" % self.zvar.name )

        ## get some structure 
        if isinstance ( self.pdf , ROOT.RooAddPdf ) :
            for p in self.pdf.pdfList    ()    : self.alist1.add ( p )
            for f in self.pdf.orig_fracs ()[0] : self.alist2.add ( f )
            

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
    """ The most trival 3D-model - constant
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
## @class Constrained
#  Helper mixin/base class for creation of constrained PDFs
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date 2023-04-18
class Constrained(object) :
    """ Helper mixin/base class for creation of constrained PDFs
    """
    def __init__  ( self , original , *constraints ) :

        assert isinstance  ( original , APDF1 ) , '"original" must be APDF1!'
        assert constraints , '"constraints" are not specified!'
        
        self.__original_pdf     = original        
        self.__constraints      = () 
        
        ## add constraints 
        self.add_constraints ( *constraints )
        
    # =========================================================================
    ## add more constraints  
    def add_constraints ( self , *constraints ) :
        """ Add more constraints"""
        
        assert constraints , '"Constraints" are not specified!'

        cnts1 = [ c         for c in cnts if isinstance ( c , ROOT.RooAbsPdf ) ]
        cnts2 = [ c.roo_pdf for c in cnts if isinstance ( c , APDF1          ) ]
        assert len ( cnts1 ) + len ( cnts2 ) == len ( cnts ) , \
               'Invalid constraints are specified %s' % str ( constraints ) 
        
        ## safe/preserve the previous PDF 
        if self.pdf       : self.aux_keep.append ( self.pdf       ) 

        
        self.__constraints     = set   ( self.constraints ) + set ( cnts1 + cnts2 )
        self.__constraints     = tuple ( self.constraints ) 
        
        ## create the actual PDF
        from ostp.fitting.consraint import make_constrained as mk_constrained 
        self.pdf = mk_constrained ( self.original_pdf.roo_pdf , *self.constraints )
        
    @property
    def original_pdf ( self ) :
        """'original_pdf' : get the original PDF (same as `unconstrained`)"""
        return self.__original_pdf

    @property
    def unconstrained ( self ) :
        """'unconstrained' : get the unconstrained PDF (same as `original_pdf`)"""
        return self.__original_pdf
    
    @property
    def constraints  ( self ) :
        """'constraints' : get the constrainsts (as list of RooAbsPdf)"""
        return self.__constraints


# =============================================================================
## @class Constrained1D
#  PDF1 with constraints
#  @code
#  opdf1 = ...
#  constrains = cpdf1, cpdf2, cpdf3
#  cpdf = Constrained1D ( opdf1 , consraints )
#  @endcode
class Constrained1D(PDF1,Constrained) :
    """ PDF1 with constraints
    >>> opdf1 = ...
    >>> constrains = cpdf1, cpdf2, cpdf3
    >>> cpdf = Constrained1D ( opdf1 , consraints )
    """    
    def __init__ ( self          ,
                   original      ,
                   *constraints  ) :

        assert isinstance  ( original , PDF1 ) , '"original" must be PDF1!'
        assert constraints , '"constraints" are not specified!'
        
        name = 'Constrained_%s' % original.name

        ## initialize the 1st base 
        PDF1.__init__        ( self , name , original.xvar )
        
        ## initialize the 2nd base
        Constrained.__init__ ( self , original , *constraints ) 
        
        ## copy structural elements 
        self.copy_structures     ( self.original_pdf )
        ## copy drawing options 
        self.draw_options.update ( self.original_pdf.draw_options )
        ## copy fit options 
        self.fit_options = self.original_pdf.fit_options 
        
        ## save the configuration
        self.config = {
            'name'        : self.name         ,
            'original'    : self.original_pdf ,
            'constraints' : self.constraints 
            }
        
    ## access to attributes of original PDF
    def __getattr__ ( self , attr ) :
        """ Delegate the access to missing attributes to the original (unconstrained) PDF"""
        return getattr ( self.original_pdf , attr ) 

# =============================================================================
## @class Constrained2D
#  PDF2 with constraints
#  @code
#  opdf2 = ...
#  constrains = cpdf1, cpdf2, cpdf3
#  cpdf = Constrained2D ( opdf2 , consraints )
#  @endcode
class Constrained2D(PDF2,Constrained) :
    """ PDF2 with constraints
    >>> opdf2 = ...
    >>> constrains = cpdf1, cpdf2, cpdf3
    >>> cpdf = Constrained2D ( opdf2 , consraints )
    """
    
    def __init__ ( self          ,
                   original      ,
                   *constraints  ) :

        assert isinstance  ( original , PDF2 ) , '"original" must be PDF2!'
        assert constraints , '"constraints" are not specified!'
        
        name = 'Constrained_%s' % original.name

        ## initialize the 1st base 
        PDF2.__init__ ( self , name , original.xvar , original.yvar )
        
        ## initialize the 2nd base
        Constrained.__init__ ( self , original , *constraints ) 
        
        ## copy structural elements 
        self.copy_structures     ( self.original_pdf )
        ## copy drawing options 
        self.draw_options.update ( self.original_pdf.draw_options )
        ## copy fit options 
        self.fit_options = self.original_pdf.fit_options 
        
        ## save the configuration
        self.config = {
            'name'        : self.name         ,
            'original'    : self.original_pdf ,
            'constraints' : self.constraints 
            }

    ## access to  attributes of original PDF 
    def __getattr__ ( self , attr ) :
        """Delegate the access to missing attribite to the original (unconstrained) PDF"""
        return getattr ( self.__original_pdf , attr ) 

# =============================================================================
## @class Constrained3D
#  PDF3 with constraints
#  @code
#  opdf3 = ...
#  constrains = cpdf1, cpdf2, cpdf3
#  cpdf = Constrained3D ( opdf3 , consraints )
#  @endcode
class Constrained3D(PDF3,Constrained) :
    """ PDF3 with constraints
    >>> opdf3 = ...
    >>> constrains = cpdf1, cpdf2, cpdf3
    >>> cpdf = Constrained3D ( opdf3 , consraints )
    """    
    def __init__ ( self         ,
                   original     ,
                   *constraints ) : 

        assert isinstance  ( original , PDF3 ) , '"original" must be PDF3!'
        assert constraints , '"constraints" are not specified!'
        
        name = 'Constrained_%s' % original.name

        ## initialize the 1st base base 
        PDF3.__init__ ( self , name , original.xvar , original.yvar , original.zvar ) 

        ## initialize the 2nd base
        Constrained.__init__ ( self , original , *constraints ) 
        
        ## copy structural elements 
        self.copy_structures     ( self.original_pdf )
        ## copy drawing options 
        self.draw_options.update ( self.original_pdf.draw_options )
        ## copy fit options 
        self.fit_options = self.original_pdf.fit_options 
        
        ## save the configuration
        self.config = {
            'name'        : self.name         ,
            'original'    : self.original_pdf ,
            'constraints' : self.constraints  
            }

    ## access to  attributes of original PDF 
    def __getattr__ ( self , attr ) :
        """ Delegate the access to missing attribite to the original (unconstrained) PDF"""
        return getattr ( self.__original_pdf , attr ) 

# =============================================================================
## generic shapes
# =============================================================================
# =============================================================================
## Generic 1D-shape from C++ callable
#  @see Ostap::Models:Shape1D
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2020-07-20
class Shape1D_pdf(PDF1) :
    """ Generic 1D-shape from C++ callable
    - see Ostap::Models:Shape1D
    """
    
    def __init__ ( self , name , shape , xvar , tag = 0 ) :

        
        th1 = None 
        if   isinstance ( shape , ROOT.TH1           ) and 1 == shape.dim() : th1 = shape
        elif isinstance ( shape , Ostap.Math.Histo3D )                      : th1 = shape.h ()
        
        if th1 and not xvar : xvar = th1.xminmax()
        
        if isinstance ( shape , ROOT.TH1 ) and 1 == shape.dim() :
            
            shape      = Ostap.Math.Histo1D     ( shape )
            tag        = Ostap.Utils.hash_histo ( shape )
            
        elif hasattr  ( shape , 'tag' ) and not tag : 

            tag = shape.tag() 
            
        ##  initialize the base 
        PDF1.__init__ ( self , name , xvar ) 
        
        self.__shape = shape
        self.__tag   = tag
        
        if isinstance ( self.shape , Ostap.Math.Histo1D ) :
            
            ## create the actual pdf
            self.pdf = Ostap.Models.Histo1D (
                self.roo_name ( 'histo1_' ) , 
                "Histo-1D %s" % self.name   ,
                self.xvar                   ,
                self.shape                  )
            
        else :

            ## create the actual pdf
            self.pdf = Ostap.Models.Shape1D ( ## create ? 
                self.roo_name ( 'shape1_' ) , 
                "Shape-1D %s" % self.name   ,
                self.xvar                   ,
                self.shape                  ,
                self.tag                    ) 
            
        ## save the configuration
        self.config = {
            'name'    : self.name    , 
            'shape'   : self.shape   , 
            'xvar'    : self.xvar    ,
            'tag'     : self.tag     ,
            }
        
    @property
    def shape  ( self ) :
        """'shape': the actual C++ callable shape"""
        return self.__shape 
    @property
    def tag   ( self ) :
        """'tag' : unique tag used for cache-integration"""
        return self.__tag

# ============================================================================= 
## Generic 2D-shape from C++ callable
#  @see Ostap::Models:Shape2D
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2020-07-20
class Shape2D_pdf(PDF2) :
    """ Generic 2D-shape from C++ callable
    - see Ostap::Models:Shape2D
    """
    
    def __init__ ( self , name , shape , xvar , yvar , tag = 0 ) :
        
        
        th2 = None 
        if   isinstance ( shape , ROOT.TH2           ) and 2 == shape.dim() : th2 = shape
        elif isinstance ( shape , Ostap.Math.Histo3D )                      : th2 = shape.h ()
        
        if th2 and not xvar : xvar = th2.xminmax()
        if th2 and not yvar : yvar = th2.yminmax()

        if isinstance ( shape , ROOT.TH2 ) and 2 == shape.dim () :
            
            self.histo_obj = shape
            shape          = Ostap.Math.Histo2D     ( shape )
            tag            = Ostap.Utils.hash_histo ( shape ) 
            
        elif hasattr ( shape , 'tag' ) and not tag :
            
            tag = shape.tag() 
            
        ##  iniialize the base 
        PDF2.__init__ ( self , name , xvar , yvar ) 
        
        self.__shape = shape
        self.__tag   = tag
        
        if isinstance ( self.shape , Ostap.Math.Histo2D ) :
            
            ## create the actual pdf
            self.pdf = Ostap.Models.Histo2D (
                self.roo_name ( 'histo2_' ) , 
                "Histo-2D %s" % self.name   ,
                self.xvar                   ,
                self.yvar                   ,
                self.shape                  )            
        else :

            ## create the actual pdf
            self.pdf = Ostap.Models.Shape2D (
                self.roo_name  ( 'shape2_' ) , 
                "Shape-2D %s" % self.name    ,
                self.xvar                    ,
                self.yvar                    ,
                self.shape                   ,
                self.tag                     )
            
        ## save the configuration
        self.config = {
            'name'    : self.name    , 
            'shape'   : self.shape   , 
            'xvar'    : self.xvar    , 
            'yvar'    : self.yvar    ,
            'tag'     : self.tag     , 
            }
        
    @property
    def shape  ( self ) :
        """'shape' : the actual C++ callable shape"""
        return self.__shape 
    @property
    def tag   ( self ) :
        """'tag' : unique tag used for cache-integration"""
        return self.__tag 

# ============================================================================= 
## Generic 3D-shape from C++ callable
#  @see Ostap::Models:Shape3D
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2020-07-20
class Shape3D_pdf(PDF3) :
    """ Generic 3D-shape from C++ callable
    - see Ostap::Models:Shape3D
    """
    
    def __init__ ( self , name , shape , xvar , yvar , zvar , tag = 0 ) :

        th3 = None 
        if   isinstance ( shape , ROOT.TH3           ) and 3 == shape.dim() : th3 = shape
        elif isinstance ( shape , Ostap.Math.Histo3D )                      : th3 = shape.h ()
        
        if th3 and not xvar : xvar = th3.xminmax()
        if th3 and not yvar : yvar = th3.yminmax()
        if th3 and not zvar : zvar = th3.zminmax()
                
        if isinstance ( shape , ROOT.TH3 ) :
            
            self.histo = shape
            shape      = Ostap.Math.Histo3D     ( shape )
            tag        = Ostap.Utils.hash_histo ( shape ) 
            
        elif hasattr ( shape , 'tag' ) and not tag : 
            tag = shape.tag() 
            
        ##  iniialize the base 
        PDF3.__init__ ( self , name , xvar , yvar , zvar ) 
        
        self.__shape = shape
        self.__tag   = tag
        
        if isinstance ( self.shape , Ostap.Math.Histo2D ) :
            
            ## create the actual pdf
            self.pdf = Ostap.Models.Histo3D ( self.roo_name ( 'histo3_' ) , 
                                              "Histo-3D %s" % self.name   ,
                                              self.xvar                   ,
                                              self.yvar                   ,
                                              self.zvar                   ,
                                              self.shape                  )
        else :
            
            ## create the actual pdf
            self.pdf = Ostap.Models.Shape3D (
                self.roo_name ( 'shape3_' ) , 
                "Shape-3D %s" % self.name   ,
                self.xvar                   ,
                self.yvar                   ,
                self.zvar                   ,
                self.shape                  ,
                self.tag                    )
            
        ## save the configuration
        self.config = {
            'name'    : self.name    , 
            'shape'   : self.shape   , 
            'xvar'    : self.xvar    , 
            'yvar'    : self.yvar    , 
            'zvar'    : self.zvar    , 
            'tag'     : self.tag     , 
            }
            
    @property
    def shape  ( self ) :
        """'shape' : the actual C++ callable shape"""
        return self.__shape  
    @property
    def tag   ( self ) :
        """'tag' : unique tag used for cache-integration"""
        return self.__tag
    

# =============================================================================
## Histo*D PDFs
# =============================================================================

# =============================================================================
## Use histogram as PDF
#  @see Ostap::Models::Histo1D
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2020-07-20
class Histo1D_pdf(PDF1) :
    """ Use histgram as PDF
    - see `Ostap.Models.Histo1D
    """
        
    def __init__ ( self , name , histo , xvar ) :

        assert isinstance ( histo , Ostap.Math.Histo1D ) or  \
               ( isinstance ( histo , ROOT.TH1 ) and 1 == histo.dim()  ) , 'Invalid histogram object'
        
        name = name if name else VarMaker.generate_name ( prefix = 'h1d' , suffix = 'pdf' , name = histo.name ) 

        th1 = histo if isinstance ( histo , ROOT.TH1 ) else histo.h ()         
        if not xvar : xvar = th1.xminmax ()

        
        ##  initialize the base 
        PDF1.__init__ ( self , name , xvar ) 

        shape = histo 
        if not isinstance ( shape , Ostap.Math.Histo1D ) :
            shape = Ostap.Math.Histo1D ( shape )

        self.__shape = shape 
        
        ## create the actual pdf
        self.pdf = Ostap.Models.Histo1D (
            self.roo_name ( 'histo1_' ) , 
            "Histo-1D %s" % self.name   ,
            self.xvar                   ,
            self.shape                  )            
        
        ## save the configuration
        self.config = {
            'name'    : self.name    , 
            'histo'   : self.shape   , 
            'xvar'    : self.xvar    , 
            }
            
    @property
    def shape  ( self ) :
        """'shape' : the actual C++ callable shape for TH1"""
        return self.__shape

# =============================================================================
## Use histogram as PDF
#  @see Ostap::Models::Histo2D
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2020-07-20
class Histo2D_pdf(PDF2) :
    """Use histgram as PDF
    - see `Ostap.Models.Histo2D
    """
        
    def __init__ ( self , name , histo , xvar , yvar ) :

        assert isinstance ( histo , Ostap.Math.Histo2D ) or  \
               ( isinstance ( histo , ROOT.TH2 ) and 2 == histo.dim()  ) , 'Invalid histogram object'
        
        name = name if name else VarMaker.generate_name ( prefix = 'h2d' , suffix = 'pdf' , name = histo.name ) 

        th2 = histo if isinstance ( histo , ROOT.TH2 ) else histo.h ()         
        if not xvar : xvar = th2.xminmax ()
        if not yvar : yvar = th2.yminmax ()
        
        ##  initialize the base 
        PDF2.__init__ ( self , name , xvar , yvar ) 

        shape = histo 
        if not isinstance ( shape , Ostap.Math.Histo2D ) :
            shape = Ostap.Math.Histo2D ( shape )

        self.__shape = shape 
        
        ## create the actual pdf
        self.pdf = Ostap.Models.Histo2D (
            self.roo_name ( 'histo2_' ) , 
            "Histo-2D %s" % self.name   ,
            self.xvar                   ,
            self.yvar                   ,
            self.shape                  )            
        
        ## save the configuration
        self.config = {
            'name'    : self.name    , 
            'histo'   : self.shape   , 
            'xvar'    : self.xvar    , 
            'yvar'    : self.yvar    ,
            }
        
    @property
    def shape  ( self ) :
        """'shape' : the actual C++ callable shape for TH2"""
        return self.__shape

# =============================================================================
## Use histogram as PDF
#  @see Ostap::Models::Histo3D
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2020-07-20
class Histo3D_pdf(PDF3) :
    """Use histogram as PDF
    - see `Ostap.Models.Histo3D
    """
        
    def __init__ ( self , name , histo , xvar , yvar , zvar ) :

        assert isinstance ( histo , Ostap.Math.Histo3D ) or  \
               ( isinstance ( histo , ROOT.TH3 ) and 3 == histo.dim()  ) , 'Invalid histogram object'
        
        name = name if name else VarMaker.generate_name ( prefix = 'h3d' , suffix = 'pdf' , name = histo.name ) 

        th3 = histo if isinstance ( histo , ROOT.TH3 ) else histo.h ()         
        if not xvar : xvar = th3.xminmax ()
        if not yvar : yvar = th3.yminmax ()
        if not zvar : zvar = th3.zminmax ()
        
        ##  initialize the base 
        PDF3.__init__ ( self , name , xvar , yvar , zvar ) 

        shape = histo 
        if not isinstance ( shape , Ostap.Math.Histo3D ) :
            shape = Ostap.Math.Histo3D ( shape )

        self.__shape = shape 
        
        ## create the actual pdf
        self.pdf = Ostap.Models.Histo3D (
            self.roo_name ( 'histo3_' ) , 
            "Histo-3D %s" % self.name   ,
            self.xvar                   ,
            self.yvar                   ,
            self.zvar                   ,
            self.shape                  )            
        
        ## save the configuration
        self.config = {
            'name'    : self.name    , 
            'histo'   : self.shape   , 
            'xvar'    : self.xvar    , 
            'yvar'    : self.yvar    ,
            'zvar'    : self.zvar    ,
            }
            
    @property
    def shape  ( self ) :
        """'shape' : the actual C++ callable shape for TH3"""
        return self.__shape
            

# =============================================================================
## H*D_pdf
# =============================================================================

# =============================================================================
## simple convertor of 1D-histogram into PDF
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H1D_pdf(PDF1) :
    """ Simple convertor of 1D-histogram into PDF
    """
    def __init__ ( self            ,
                   name            ,
                   histo           ,
                   xvar    = None  ,
                   density = False ,
                   order   = 0     , ## interpolation order 
                   silent  = False ) :
        
        assert isinstance ( order, integer_types ) and 0 <= order ,\
               'Invalid interpolation order: %s/%s' % ( order , type ( order ) )

        name = name if name else VarMaker.generate_name ( prefix = 'h1d' , suffix = 'pdf' , name = histo.name ) 

        self.__ds = H1D_dset ( histo             ,
                               xaxis   = xvar    ,
                               density = density ,
                               silent  = silent  )
        
        PDF1    .__init__ ( self , name  , self.ds.xaxis ) 

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
        ## self.signals.add ( self.pdf ) 
        
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
    def ds ( self ) :
        """'ds' : the H1D_dset object"""
        return self.__ds 
    @property     
    def xaxis  ( self ) :
        """The histogram x-axis variable (same as xvar)"""
        return self.xvar
    @property
    def histo ( self ) :
        """The  histogram itself"""
        return self.ds.histo    
    @property
    def density( self ) :
        """Treat the histo as 'density' histogram?"""
        return self.ds.density    
    @property
    def skip_zero ( self ) :
        """'skip_zero' : skip zero bins for weighted dataset in histo?"""
        return self.ds.skip_zero    
    @property
    def silent( self ) :
        """Use the silent mode?"""
        return self.ds.silent
    @property
    def dset ( self ) :
        """'dset' : ROOT.RooDataHist object"""
        return self.ds.dset
    @property
    def histo_hash ( self ) :
        """Hash value for the histogram"""
        return self.ds.histo_hash
    @property
    def weight ( self ) :
        """'weight' : get weight variable if defined, None otherwise"""
        return self.ds.wvar

    @property
    def order  ( self ) :
        """'order': interpolation order"""
        return self.pdf.getInterpolationOrder () 
    @order.setter
    def order  ( self , value ) :
        assert isinstance ( value , integer_types ) and 0 <= value,\
               'Invalid interpolation order %s/%s' % ( value , type ( value ) )
        self.pdf.setInterpolationOrder ( value )

# ===================================================it==========================
## simple convertor of 2D-histogram into PDF
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H2D_pdf(PDF2) :
    """Simple convertor of 2D-histogram into PDF 
    """
    def __init__ ( self            ,
                   name            ,
                   histo           ,
                   xvar    = None  , 
                   yvar    = None  ,
                   density = False ,
                   order   = 0     , ## interpolation order 
                   silent  = False ) :

        assert isinstance ( order, integer_types ) and 0 <= order ,\
               'Invalid interpolation order: %s/%s' % ( order , type ( order ) )

        name = name if name else VarMaker.generate_name ( prefix = 'h1d' , suffix = 'pdf' , name = histo.name ) 

        self.__ds = H2D_dset ( histo             ,
                               xaxis   = xvar    ,
                               yaxis   = yvar    , 
                               density = density ,  
                               silent  = silent  )
        
        PDF2    .__init__ ( self , name  , self.ds.xaxis , self.ds.yaxis ) 
        
        self.__vset  = ROOT.RooArgSet  ( self.xvar , self.yvar )

        #
        ## finally create PDF :
        #
        with roo_silent ( silent ) : 
            self.pdf    = ROOT.RooHistPdf (
                self.roo_name  ( 'histo2_' ) , 
                'Histo-2D PDF: %s/%s' % ( self.histo.GetName() , self.histo.GetTitle() ) , 
                self.__vset , 
                self.dset   ,
                order       )
            
        ## and declare it be be a "signal"
        ## self.signals.add ( self.pdf ) 

        ## save the configuration
        self.config = {
            'name'    : self.name    , 
            'histo'   : self.histo   , 
            'xvar'    : self.xvar    , 
            'yvar'    : self.yvar    , 
            'density' : self.density , 
            'silent'  : self.silent  ,             
            'order'   : self.order   ,             
            }
        
    @property
    def ds ( self ) :
        """'ds' : the H2D_dset object"""
        return self.__ds 
    @property     
    def xaxis  ( self ) :
        """The histogram x-axis variable (same as xvar)"""
        return self.xvar
    @property     
    def yaxis  ( self ) :
        """The histogram y-axis variable (same as yvar)"""
        return self.yvar
    @property
    def histo ( self ) :
        """The  histogram itself"""
        return self.ds.histo    
    @property
    def density( self ) :
        """Treat the histo as 'density' histogram?"""
        return self.ds.density    
    @property
    def skip_zero ( self ) :
        """'skip_zero' : skip zero bins for weighted dataset in histo?"""
        return self.ds.skip_zero    
    @property
    def silent( self ) :
        """Use the silent mode?"""
        return self.ds.silent
    @property
    def dset ( self ) :
        """'dset' : ROOT.RooDataHist object"""
        return self.ds.dset
    @property
    def histo_hash ( self ) :
        """Hash value for the histogram"""
        return self.ds.histo_hash
    @property
    def weight ( self ) :
        """'weight' : get weight variable if defined, None otherwise"""
        return self.ds.wvar

    @property
    def order  ( self ) :
        """'order' : interpolation order"""
        return self.pdf.getInterpolationOrder () 
    @order.setter
    def order  ( self , value ) :
        assert isinstance ( value , integer_types ) and 0 <= value,\
               'Invalid interpolation order %s/%s' % ( value , type ( value ) )
        self.pdf.setInterpolationOrder ( value )


# =============================================================================
## simple convertor of 3D-histogram into PDF
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H3D_pdf(PDF3) :
    """Simple convertor of 3D-histogram into PDF 
    """
    def __init__ ( self            ,
                   name            ,
                   histo           ,
                   xvar    = None  , 
                   yvar    = None  ,
                   zvar    = None  ,
                   density = False ,
                   order   = 0     , 
                   silent  = False ) :
        
        assert isinstance ( order, integer_types ) and 0 <= order ,\
               'Invalid interpolation order: %s/%s' % ( order , type ( order ) )

        name = name if name else VarMaker.generate_name ( prefix = 'h1d' , suffix = 'pdf' , name = histo.name ) 

        self.__ds = H3D_dset ( histo ,
                               xaxis   = xvar    ,
                               yaxis   = yvar    ,
                               zaxis   = zvar    ,
                               density = density ,
                               silent  = silent  )
        
        PDF3    .__init__ ( self , name = name  ,
                            xvar = self.ds.xvar ,
                            yvar = self.ds.yvar ,
                            zvar = self.ds.zvar ) 
        
        self.__vset  = ROOT.RooArgSet  ( self.xvar , self.yvar , self.zvar )
        
        #
        ## finally create PDF :
        #
        with roo_silent ( silent ) : 
            self.pdf    = ROOT.RooHistPdf (
                self.roo_name ( 'histo3_' ) , 
                'Histo-3D PDF: %s/%s' % ( histo3.GetName() , histo2.GetTitle() ) , 
                self.__vset  , 
                self.dset    ,
                order        )

        ## and declare it be be a "signal"
        ## self.signals.add ( self.pdf ) 
            
        ## save the configuration
        self.config = {
            'name'    : self.name    , 
            'histo'   : self.histo   , 
            'xvar'    : self.xvar    , 
            'yvar'    : self.yvar    , 
            'zvar'    : self.zvar    , 
            'density' : self.density , 
            'silent'  : self.silent  ,             
            'order'   : self.order   ,             
            }

    @property
    def ds ( self ) :
        """'ds' : the H3D_dset object"""
        return self.__ds 
    @property     
    def xaxis  ( self ) :
        """The histogram x-axis variable (same as xvar)"""
        return self.xvar
    @property     
    def yaxis  ( self ) :
        """The histogram y-axis variable (same as yvar)"""
        return self.yvar
    @property     
    def zaxis  ( self ) :
        """The histogram z-axis variable (same as zvar)"""
        return self.zvar
    
    @property
    def histo ( self ) :
        """The  histogram itself"""
        return self.ds.histo    
    @property
    def density( self ) :
        """Treat the histo as 'density' histogram?"""
        return self.ds.density    
    @property
    def skip_zero ( self ) :
        """'skip_zero' : skip zero bins for weighted dataset in histo?"""
        return self.ds.skip_zero    
    @property
    def silent( self ) :
        """Use the silent mode?"""
        return self.ds.silent
    @property
    def dset ( self ) :
        """'dset' : ROOT.RooDataHist object"""
        return self.ds.dset
    @property
    def histo_hash ( self ) :
        """Hash value for the histogram"""
        return self.ds.histo_hash
    @property
    def weight ( self ) :
        """'weight' : get weight variable if defined, None otherwise"""
        return self.ds.wvar
        
    @property
    def order  ( self ) :
        """'order' : interpolation order"""
        return self.pdf.getInterpolationOrder () 
    @order.setter
    def order  ( self , value ) :
        assert isinstance ( value , integer_types ) and 0 <= value,\
               'Invalid interpolation order %s/%s' % ( value , type ( value ) )
        self.pdf.setInterpolationOrder ( value )


# ============================================================================
## Flat*D
# ============================================================================

# ============================================================================
## @class Flat1D
#  The most trivial 1D-model - constant
#  @code 
#  pdf = Flat1D( 'flat' , xvar = ... )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
class Flat1D(PDF1) :
    """ The most trival 1D-model - constant
    >>> pdf = Flat1D ( 'flat' , xvar = ... )
    """
    def __init__ ( self , xvar , name = '' , title = '' ) :
        
        name = name if name else self.generate_name ( prefix = 'flat1D_')
        PDF1.__init__ ( self  , name , xvar ) 
        
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
## @class Flat2D
#  The most trivial 2D-model - constant
#  @code 
#  pdf = Flat2D( 'flat' , xvar = ...  , yvar = ... )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
class Flat2D(PDF2) :
    """ The most trival 2D-model - constant
    >>> pdf = Flat2D( 'flat' , xvar = ...  , yvar = ... )
    """
    def __init__ ( self , xvar , yvar , name = '' ,  title = '' ) :

        name = name if name else self.generate_name ( prefix = 'flat2D_')                            
        PDF2.__init__ ( self  , name , xvar , yvar ) 
                        
        if not title : title = 'flat2(%s)' % name 
        self.pdf = Ostap.Models.Uniform ( name , title , self.xvar , self.yvar )
        assert 2 == self.pdf.dim() , 'Flat2D: wrong dimensionality!'
        
        ## save configuration
        self.config = {
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'name'     : self.name ,            
            'title'    : title     ,             
            }                   

# ===========================================================================
## @class Flat3D
#  The most trivial 3D-model - constant
#  @code 
#  pdf = Flat3D( 'flat' , xvar = ...  , yvar = ... , zvar = ... )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
class Flat3D(PDF3) :
    """ The most trival 3D-model - constant
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
## Sum*D
# =============================================================================

# =============================================================================
## @class Sum1D
#  Non-extended sum of several PDFs
#  It is just a small wrapper for <code>ROOT.RooAddPdf</code>
#  @see RooAddPdf 
class Sum1D (PDF1,Fractions) :
    """ Non-extended sum of several PDFs:    
    It is just a small wrapper for <code>ROOT.RooAddPdf</code>
    - see RooAddPdf     
    >>> sum  = Sum1D ( [ pdf1 , pdf2 , pdf3 ]  )     
    """
    def __init__ ( self              ,
                   pdfs              , ## input list of PDFs  
                   xvar      = None  , 
                   name      = ''    ,
                   recursive = True  ,
                   prefix    = 'f'   , ## prefix for fraction names 
                   suffix    = ''    , ## suffix for fraction names 
                   fractions = None  ,
                   fix_norm  = False ) :

        assert 2 <= len ( pdfs ) , 'Sum1D: at least two PDFs are needed!'
        pdf_list = []           
        for i , p in enumerate ( pdfs ) :
            cmp , xvar = self.make_PDF1 ( p , xvar = xvar ) 
            pdf_list.append ( cmp )
            
        ## generic name 
        patname =  '+'.join ( '(%s)' % p.name for p in pdf_list )
        ## check the instance name 
        name    = name if name else self.new_name ( patname ) 
        
        ## ininialize the base class
        PDF1.__init__ ( self , name , xvar ) 
        Fractions.__init__ ( self , pdf_list       ,
                             prefix    = prefix    ,
                             suffix    = suffix    ,
                             recursive = recursive ,
                             fractions = fractions ) 

        for p in self.pdfs      : self.alist1.add ( p.pdf )
        for f in self.frac_list : self.alist2.add ( f     )
        
        ## finally build PDF
        self.pdf = ROOT.RooAddPdf ( self.new_roo_name ( patname , suffix ) , 
                                    patname        ,
                                    self.alist1    ,
                                    self.alist2    ,
                                    self.recursive )
        
        ## attention!
        if fix_norm : self.pdf.SetCoefNormalization ( self.vars )
            
        self.config = {
            'pdfs'      : self.pdfs      ,
            'xvar'      : self.xvar      ,
            'name'      : self.name      , 
            'prefix'    : self.prefix    ,
            'suffix'    : self.suffix    ,
            'fractions' : self.fractions ,
            'recursive' : self.recursive ,
            'fix_norm'  : self.fix_norm  
            }

    # ==========================================================================
    ## Get the parameter by name or component by oindex 
    #  @code
    #  pdf = ...
    #  a   = pdf['A']
    #  c   = pdf[2] 
    #  @endcode
    def __getitem__ ( self , index ) :
        """ Get the component by index (or parameter by name) 
        >>> pdf = ...
        >>> a   = pdf['A']
        >>> c   = pdf[2] 
        """
        if isinstance ( index , integer_types  ) or isinstance ( index , slice ) :
            return self.pdfs [ index ] 
        return super().__getitem__ ( index ) 

    @property
    def fix_norm ( self ) :
        """`fix-norm`: 
        - see `ROOT.RooAbsPdf.SetCoefNormalization`
        - see `ROOT.RooAbsPdf.getCoefNormalization`
        """
        pars = self.pdf.getCoefNormalization()
        return True if pars else False 
    
# =============================================================================
## @class Sum2D
#  Non-extended sum of several PDFs
#  It is just a small wrapper for <code>ROOT.RooAddPdf</code>
#  @see RooAddPdf 
class Sum2D (PDF2,Fractions) :
    """ Non-extended sum of several PDFs:
    
    It is just a small wrapper for <code>ROOT.RooAddPdf</code>
    - see RooAddPdf 
    
    >>> sum  = Sum2D ( [ pdf1 , pdf2 , pdf3 ]  ) 

    """
    def __init__ ( self              ,
                   pdfs              , ## input list of PDFs  
                   xvar      = None  , 
                   yvar      = None  , 
                   name      = ''    ,
                   recursive = True  ,
                   prefix    = 'f'   , ## prefix for fraction names 
                   suffix    = ''    , ## suffix for fraction names 
                   fractions = None  ,
                   fix_norm  = False ) :
        
        assert 2 <= len ( pdfs ) , 'Sum2D: at least two PDFs are needed!'
        
        pdf_list = []           
        for i , p in enumerate ( pdfs ) :
            cmp , xvar , yvar = self.make_PDF2 ( p , xvar = xvar , yvar = yvar ) 
            pdf_list.append ( cmp )

        ## generic name 
        patname =  '+'.join ( '(%s)' % p.name for p in pdf_list )

        ## check the instance name 
        name    = name if name else self.new_name ( patname ) 
        
        ## ininialize the base classes 
        PDF2.     __init__ ( self , name , xvar , yvar )        
        Fractions.__init__ ( self , pdf_list       , 
                             prefix    = prefix    ,
                             suffix    = suffix    ,
                             recursive = recursive ,
                             fractions = fractions )

        for p in self.pdfs      : self.alist1.add ( p.pdf )
        for f in self.frac_list : self.alist2.add ( f     )
        
        ## finally build PDF
        self.pdf = ROOT.RooAddPdf ( self.new_roo_name ( patname , suffix ) , 
                                    patname        ,
                                    self.alist1    ,
                                    self.alist2    ,
                                    self.recursive )

        ## attention!
        if fix_norm : self.pdf.SetCoefNormalization ( self.vars )
        
        self.config = {
            'pdfs'      : self.pdfs      ,
            'xvar'      : self.xvar      ,
            'yvar'      : self.yvar      ,
            'name'      : self.name      , 
            'prefix'    : self.prefix    ,
            'suffix'    : self.suffix    ,
            'fractions' : self.fractions ,
            'recursive' : self.recursive ,        
            'fix_norm'  : self.fix_norm  
            }
  
    # ==========================================================================
    ## Get the parameter by name or component by oindex 
    #  @code
    #  pdf = ...
    #  a   = pdf['A']
    #  c   = pdf[2] 
    #  @endcode
    def __getitem__ ( self , index ) :
        """ Get the component by index (or parameter by name) 
        >>> pdf = ...
        >>> a   = pdf['A']
        >>> c   = pdf[2] 
        """
        if isinstance ( index , integer_types  ) or isinstance ( index , slice ) :
            return self.pdfs [ index ] 
        return super().__getitem__ ( index )
    
    @property
    def fix_norm ( self ) :
        """`fix-norm`: 
        - see `ROOT.RooAbsPdf.SetCoefNormalization`
        - see `ROOT.RooAbsPdf.getCoefNormalization`
        """
        pars = self.pdf.getCoefNormalization()
        return True if pars else False 
    
        
# =============================================================================
## @class Sum3D
#  Non-extended sum of several PDFs
#  It is just a small wrapper for <code>ROOT.RooAddPdf</code>
#  @see RooAddPdf 
class Sum3D (PDF3,Fractions) :
    """ Non-extended sum of several PDFs:
    
    It is just a small wrapper for <code>ROOT.RooAddPdf</code>
    - see RooAddPdf 
    
    >>> sum  = Sum3D ( [ pdf1 , pdf2 , pdf3 ]  ) 
    
    """
    def __init__ ( self              ,
                   pdfs              , ## input list of PDFs  
                   xvar      = None  , 
                   yvar      = None  , 
                   zvar      = None  , 
                   name      = ''    ,
                   recursive = True  ,
                   prefix    = 'f'   , ## prefix for fraction names 
                   suffix    = ''    , ## suffix for fraction names 
                   fractions = None  ,
                   fix_norm  = False ) :

        assert 2 <= len ( pdfs ) , 'Sum3D: at least two PDFs are needed!'

        pdf_list = []           
        for i , p in enumerate ( pdfs ) :
            cmp , xvar , yvar , zvar = self.make_PDF3 ( p , xvar = xvar , yvar = yvar , zvar = zvar ) 
            pdf_list.append ( cmp )
            
        ## generic name 
        patname =  '+'.join ( '(%s)' % p.name for p in pdf_list )
        ## check the instance name 
        name    = name if name else self.new_name ( patname ) 

        ## initialize the base class
        PDF3.     __init__ ( self , name , xvar , yvar , zvar ) 
        Fractions.__init__ ( self , pdf_list ,
                             prefix    = prefix    ,
                             suffix    = suffix    ,
                             recursive = recursive ,
                             fractions = fractions ) 

        for p in self.pdfs      : self.alist1.add ( p.pdf )
        for f in self.frac_list : self.alist2.add ( f     )
        
        ## finally build PDF
        self.pdf = ROOT.RooAddPdf ( self.new_roo_name ( patname , suffix ) , 
                                    patname        ,
                                    self.alist1    ,
                                    self.alist2    ,
                                    self.recursive )

        ## attention!
        if fix_norm : self.pdf.SetCoefNormalization ( self.vars )

        self.config = {
            'pdfs'      : self.pdfs      ,
            'xvar'      : self.xvar      ,
            'yvar'      : self.yvar      ,
            'zvar'      : self.zvar      ,
            'name'      : self.name      , 
            'prefix'    : self.prefix    ,
            'suffix'    : self.suffix    ,
            'fractions' : self.fractions ,
            'recursive' : self.recursive , 
            'fix_norm'  : self.fix_norm  
            }

    # ==========================================================================
    ## Get the parameter by name or component by oindex 
    #  @code
    #  pdf = ...
    #  a   = pdf['A']
    #  c   = pdf[2] 
    #  @endcode
    def __getitem__ ( self , index ) :
        """ Get the component by index (or parameter by name) 
        >>> pdf = ...
        >>> a   = pdf['A']
        >>> c   = pdf[2] 
        """
        if isinstance ( index , integer_types  ) or isinstance ( index , slice ) :
            return self.pdfs [ index ] 
        return super().__getitem__ ( index ) 

    @property
    def fix_norm ( self ) :
        """`fix-norm`: 
        - see `ROOT.RooAbsPdf.SetCoefNormalization`
        - see `ROOT.RooAbsPdf.getCoefNormalization`
        """
        pars = self.pdf.getCoefNormalization()
        return True if pars else False 
        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END 
# =============================================================================
