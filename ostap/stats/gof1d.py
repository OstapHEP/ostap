#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/gof_1d.py
#  Set of utulities for goodness-of-fit studies for 1D-fits
#  @see https://doi.org/10.1111/1467-9868.00337 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2024-09-16
# =============================================================================
""" Simple utulities for goodness-of-fit studies for 1D-fits  
- see https://doi.org/10.1111/1467-9868.00337
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2024-09-16"
__all__     = (
    'kolmogorov_smirnov' , ## Kolmogorov-Sminov GoF estimator 
    'anderson_darling'   , ## Anderson-Darling  GoF estimator 
    'cramer_von_mises'   , ## Cramer-von Mises  GoF estimator 
    'kuiper'             , ## Kuiper            GoF estimator 
    'ZK'                 , ## ZK                GoF estimator
    'ZA'                 , ## ZA                GoF estimator
    'ZC'                 , ## ZC                GoF estimator
    'GoF1D'              , ## helper utility for GoF estimate 
    'GoF1DToys'          , ## helper utility for GoF estimate with toys 
    )
# =============================================================================
from   collections              import defaultdict, namedtuple
from   ostap.core.meta_info     import root_info 
from   ostap.fitting.funbasic   import AFUN1
from   ostap.fitting.pdfbasic   import PDF1
from   ostap.core.core          import VE, Ostap
from   ostap.math.base          import doubles, axis_range, np2raw    
from   ostap.math.models        import f1_draw
from   ostap.utils.cidict       import cidict_fun
from   ostap.utils.basic        import numcpu, loop_items, typename   
from   ostap.stats.counters     import SE, EffCounter 
from   ostap.logger.pretty      import pretty_float
from   ostap.math.ve            import fmt_pretty_ve
from   ostap.math.math_ve       import significance
from   ostap.logger.symbols     import plus_minus, times
from   ostap.logger.colorized   import infostr 
import ostap.logger.table       as     T
import ostap.fitting.ds2numpy 
import ostap.fitting.roofit
import ROOT, math, numpy 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.gof1d' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Simple utilities for goodness-of-1D-fit studies' )
# =============================================================================
if  ( 6 , 32 ) <= root_info : data2vct = lambda s : s
else                        : data2vct = lambda s : doubles ( s ) 
# =============================================================================
## @var NL
#  use C++ if length od data exceeds NL, otherwise Python is OK 
NL  = 100 
# =============================================================================
## Get Kolmogorov-Smirnov statistics KS
#  @code
#  cdf_data =...
#  ks2  = kolmogorov_smirnov ( cdf_data )
#  @endcode
#  @param cdf_data sorted array of F0(X_i) - values of CDF at X data points
#  @return Kolmogorov-Smirnov statistics KS
def kolmogorov_smirnov ( cdf_data ) :
    """ Get Kolmogorov-Smirnov statistics  KS
    - `cdf_data` : sorted array of F0(X_i) - values of CDF at (sorted) X data points
    >>> cdf_data =...
    >>> ks2  = kolmogorov_smirnov ( cdf_data )
    """

    n = len ( cdf_data )
    
    ## Long numpy arrays 
    if NL < n and numpy and np2raw and isinstance ( cdf_data , numpy.ndarray ) :
        if cdf_data.dtype in ( numpy.float32 , numpy.float64 ) :
            raw_buffer , n = np2raw ( cdf_data )
            the_buffer = Ostap.Utils.make_buffer ( raw_buffer , n )
            return Ostap.GoF.kolmogorov_smirnov ( the_buffer )

    ## Long arrays to be cnverted to numpy 
    if NL < n and numpy : 
        data = numpy.asarray ( cdf_data , dtype = numpy.float64 ) 
        return kolmogorov_smirnov ( data ) 

    ## for short arraya python is OK 
    result = max ( max ( ( i + 1.0 ) / n - Fi , Fi - float ( i ) / n ) for ( i , Fi )  in enumerate ( cdf_data )  ) ## ** 2 
    return result

# =============================================================================
## Get Anderson-Darling  statistiscs AD^2
#  @code
#  cdf_data =...
#  ad2      = anderson_darling ( cdf_data )
#  @endcode
#  @param cdf_data sorted array of F0(X_i) - values of CDF at (sorted) X data points
#  @return Anderson-Darling statistics AD^2
def anderson_darling ( cdf_data ) :
    """ Get Anderson-Darling statistiscs AD^2
    - `cdf_data` : sorted array of F0(X_i) - values of CDF at X data points
    >>> cdf_data =...
    >>> ad2      = anderson_darling ( cdf_data )    
    """
    n       = len ( cdf_data )
    
    ## Long numpy arrays 
    if NL < n and numpy and np2raw and isinstance ( cdf_data , numpy.ndarray ) :
        if cdf_data.dtype in ( numpy.float32 , numpy.float64 ) :
            raw_buffer , n = np2raw ( cdf_data )
            the_buffer = Ostap.Utils.make_buffer ( raw_buffer , n )
            return Ostap.GoF.anderson_darling ( the_buffer )
        
    ## Long arrays to be cnverted to numpy 
    if NL < n and numpy : 
        data = numpy.asarray ( cdf_data , dtype = numpy.float64 ) 
        return anderson_darling ( data ) 

    ## for short array pure python is OK     
    flog    = math.log
    result  = sum ( ( i + 0.5 ) * flog ( Fi ) + ( n - i -  0.5 ) * flog ( 1 - Fi ) for ( i , Fi )  in enumerate ( cdf_data ) ) 
    ## 
    return -2.0 * result / n - n 
# =============================================================================
## Get Cramer-von Mises statistics CM^2
#  @code
#  cdf_data = ...
#  cm2      = cramer_von_mises ( cdf_data )
#  @endcode
#  @param cdf_data sorted array of F0(X_i) - values of CDF at X data points
#  @return Cramer-von Mises statistics CM^2
def cramer_von_mises ( cdf_data  ) :
    """ Get Cramer-von Mises statistics CM^2
    - `cdf_data` : sorted array of F0(X_i) - values of CDF at (sorted) X data points
    >>> cdf_data =...
    >>> cm2     = cramer_von_mises ( cdf_data )    
    """
    n       = len ( cdf_data )

    ## Long numpy arrays 
    if NL < n and numpy and np2raw and isinstance ( cdf_data , numpy.ndarray ) :
        if cdf_data.dtype in ( numpy.float32 , numpy.float64 ) :
            raw_buffer , n = np2raw ( cdf_data )
            the_buffer = Ostap.Utils.make_buffer ( raw_buffer , n )
            return Ostap.GoF.cramer_von_mises  ( the_buffer )
        
    ## Long arrays to be cnverted to numpy 
    if NL < n and numpy : 
        data = numpy.asarray ( cdf_data , dtype = numpy.float64 ) 
        return cramer_von_mises ( data ) 
    
    ## for short array pure python is OK     
    result  = sum ( ( Fi - ( i + 0.5 ) / n ) ** 2 for ( i, Fi ) in enumerate ( cdf_data ) ) 
    return result + 1 / ( 12.0 * n ) 
# =============================================================================
## Get Kuiper's statistis K
#  @code
#  cdf_data =...
#  k        = kuiper ( cdf_data )
#  @endcode
#  @param cdf_data sorted array of F0(X_i) - values of CDF at X data points
#  @return Kolmogorov-Smirnov statistics K
def kuiper ( cdf_data ) :
    """ Get Kuiper statistis  KS
    - `cdf_data` : sorted array of F0(X_i) - values of CDF at (sorted) X data points
    >>> cdf_data =...
    >>> ks2  = kolmogorov_smirnov ( cdf_data )
    """
    n       = len ( cdf_data )
    
    ## Long numpy arrays 
    if NL < n and numpy and np2raw and isinstance ( cdf_data , numpy.ndarray ) :
        if cdf_data.dtype in ( numpy.float32 , numpy.float64 ) :
            raw_buffer , n = np2raw ( cdf_data )
            the_buffer = Ostap.Utils.make_buffer ( raw_buffer , n )
            return Ostap.GoF.kuiper  ( the_buffer )
        
    ## Long arrays to be cnverted to numpy 
    if NL < n and numpy : 
        data = numpy.asarray ( cdf_data , dtype = numpy.float64 ) 
        return kuiper ( data ) 
    
    ## for short array pure python is OK     
    d_plus  = max ( ( i + 1.0 ) / n - Fi for ( i, Fi ) in enumerate ( cdf_data ) )
    d_minus = max ( Fi - ( i + 1.0 ) / n for ( i, Fi ) in enumerate ( cdf_data ) )
    return d_plus + d_minus  
# =============================================================================
## Get ZA statististics
#  @code
#  cdf_data = ...
#  za       = ZA ( cdf_data )
#  @endcode
#  @see https://doi.org/10.1111/1467-9868.00337
#  @param cdf_data sorted array of F0(X_i) - values of CDF at X data points
#  @return ZA statistics ZA 
def ZA  ( cdf_data ) :
    """ Get ZA statististics 
    - `cdf_data` : sorted array of F0(X_i) - values of CDF at (sorted) X data points
    >>> cdf_data =...
    >>> za       = ZA ( cdf_data )    
    """
    n       = len ( cdf_data )
    
    ## Long numpy arrays 
    if NL < n and numpy and np2raw and isinstance ( cdf_data , numpy.ndarray ) :
        if cdf_data.dtype in ( numpy.float32 , numpy.float64 ) :
            raw_buffer , n = np2raw ( cdf_data )
            the_buffer = Ostap.Utils.make_buffer ( raw_buffer , n )
            return Ostap.GoF.ZA  ( the_buffer )
        
    ## Long arrays to be cnverted to numpy 
    if NL < n and numpy : 
        data = numpy.asarray ( cdf_data , dtype = numpy.float64 ) 
        return ZA ( data ) 
    
    flog    = math.log
    result  = sum ( flog ( Fi ) / ( n - i - 0.5 ) + flog ( 1 - Fi ) / ( i + 0.5 ) for ( i , Fi )  in enumerate ( cdf_data ) )
    return -1 * result    
# =============================================================================
## Get Zhang's ZK statististics
#  @code
#  cdf_data = ...
#  zk       = ZK ( cdf_data )
#  @endcode
#  @see https://doi.org/10.1111/1467-9868.00337
#  @param cdf_data sorted array of F0(X_i) - values of CDF at X data points
#  @return ZK statistics ZK 
def ZK  ( cdf_data ) :
    """ Get Zhang's ZK statististics 
    - `cdf_data` : sorted array of F0(X_i) - values of CDF at X data points
    >>> cdf_data =...
    >>> zk       = ZK ( cdf_data )    
    """
    n      = len ( cdf_data )

    ## Long numpy arrays 
    if NL < n and numpy and np2raw and isinstance ( cdf_data , numpy.ndarray ) :
        if cdf_data.dtype in ( numpy.float32 , numpy.float64 ) :
            raw_buffer , n = np2raw ( cdf_data )
            the_buffer = Ostap.Utils.make_buffer ( raw_buffer , n )
            return Ostap.GoF.ZK  ( the_buffer )
        
    ## Long arrays to be cnverted to numpy 
    if NL < n and numpy : 
        data = numpy.asarray ( cdf_data , dtype = numpy.float64 ) 
        return ZK ( data ) 
    
    flog   = math.log 
    result = max ( ( i     + 0.5 ) * flog ( ( i + 0.5     ) / ( n *       Fi   ) ) +
                   ( n - i - 0.5 ) * flog ( ( n - i - 0.5 ) / ( n * ( 1 - Fi ) ) ) for ( i , Fi ) in enumerate ( cdf_data ) )
    return result 
# =============================================================================
## Get Zhang's ZC statististics
#  @code
#  cdf_data = ...
#  zc       = ZC ( cdf_data )
#  @endcode
#  @see https://doi.org/10.1111/1467-9868.00337
#  @param cdf_data sorted array of F0(X_i) - values of CDF at X data points
#  @return ZC statistics ZC
def ZC  ( cdf_data ) :
    """ Get Zhang's ZC statististics 
    - `cdf_data` : sorted array of F0(X_i) - values of CDF at (sorted) X data points
    >>> cdf_data =...
    >>> zc       = ZC ( cdf_data )    
    """
    n      = len ( cdf_data )
    
    ## Long numpy arrays 
    if NL < n and numpy and np2raw and isinstance ( cdf_data , numpy.ndarray ) :
        if cdf_data.dtype in ( numpy.float32 , numpy.float64 ) :
            raw_buffer , n = np2raw ( cdf_data )
            the_buffer = Ostap.Utils.make_buffer ( raw_buffer , n )
            return Ostap.GoF.ZC  ( the_buffer )
        
    ## Long arrays to be cnverted to numpy 
    if NL < n and numpy : 
        data = numpy.asarray ( cdf_data , dtype = numpy.float64 ) 
        return ZC ( data ) 
    
    flog   = math.log
    result = sum ( ( flog ( ( 1.0 / Fi - 1 ) / ( ( n - 0.5 ) / ( i + 0.25 ) - 1 ) ) ) ** 2 for ( i , Fi ) in enumerate ( cdf_data ) )
    return result

# =============================================================================
## Short labels for various statitical estimators 
Labels = {
    'KS' : 'Kolmogorov-Smirnov' ,
    'K'  : 'Kuiper'             ,
    'AD' : 'Anderson-Darling'   ,
    'CM' : 'Cramer-von Mises'   ,
    'ZK' : 'Zhang/ZK'           ,
    'ZA' : 'Zhang/ZA'           ,
    'ZC' : 'Zhang/ZC'           ,        
}
## lower-case shortcuts:
KS_keys = 'ks' , 'kolmogorov' , 'kolmogorovsmirnov' 
K_keys  = 'k'  , 'kuiper'  
AD_keys = 'ad' , 'anderson'   , 'andersondarling' 
CM_keys = 'cm' , 'cramer'     , 'cramervonmises' 
ZK_keys = 'zk' , 'zhangk'     , 'zhangzk'
ZA_keys = 'za' , 'zhanga'     , 'zhangza'
ZC_keys = 'zc' , 'zhangc'     , 'zhangzc'
# =========================================================================
## Clip input CDF arrays 
def vct_clip ( input , silent = True ) :
    """ Clip input CDF arrays"""
    vmin , vmax = 1.e-12 , 1 - 1.e-10
    if silent and ( numpy.min ( input ) < vmin or vmax < numpy.max ( input ) ) :
        logger.warning ( 'Adjust CDF to be %s<cdf<%s' % ( vmin , fmax ) ) 
    return numpy.clip ( input , a_min = vmin , a_max = vmax )

# ==============================================================================
## @class GoF1D
#  Goodness of 1D-fits 
#  @see https://doi.org/10.1111/1467-9868.00337
class GoF1D(object) :
    """ Goodness-of-fit 1D-case 
    - see https://doi.org/10.1111/1467-9868.00337
    """
    def __init__ ( self              ,
                   pdf               ,
                   dataset           ,
                   cdf        = None ,
                   parameters = {}   ) :
        
        assert isinstance ( pdf     , PDF1            ) , 'Invalid type of `pdf`:%s'     % typename ( pdf     )
        assert isinstance ( dataset , ROOT.RooDataSet ) , 'Invalid type of `dataset`:%s' % typename ( dataset )
        
        vars = pdf.vars
        assert 1 == len ( vars )   , 'GoF1D: Only 1D-PDFs are allowed!'
        assert pdf.xvar in dataset , 'GoF1D: `xvar`:%s is not in dataset!' % ( self.xvar.name ) 

        original_cdf = cdf
        
        cdf_ok = cdf and callable ( cdf )

        self.__store = () 
        if not cdf_ok :
            cdf = pdf.cdf ()
            ## make a try to get the underlying c++ cdf-function from PDF 
            if hasattr ( pdf.pdf , 'setPars'  ) : pdf.pdf.setPars() 
            if hasattr ( pdf.pdf , 'function' ) :
                if hasattr ( pdf.pdf , 'setPars'  ) : pdf.pdf.setPars() 
                fun = pdf.pdf.function()
                if hasattr ( fun , 'cdf' ) :
                    # =========================================================
                    try : # ===================================================
                        # =====================================================
                        a = fun.cdf ( 0.0 )
                        def cdf ( x ) : return fun.cdf ( x ) 
                        self.__store = pdf , cdf , fun
                        # =====================================================
                    except TypeError : # ======================================
                        # =====================================================
                        pass 

        ## store PDF 
        self.__pdf        = pdf        
        self.__parameters = parameters if parameters else pdf.parameters ( dataset ) 

        ## store CDF 
        self.__cdf   = cdf

        ## vectorized version of CDF 
        self.__vct_cdf = numpy.vectorize ( self.__cdf )
        
        ## data in a form of numpy sructured array
        varname = pdf.xvar.name 
        data    = dataset.tonumpy ( varname ) [ varname ] 

        ## sorted data 
        data    = numpy.sort  ( data )

        ## empirical CDF function 
        self.__ecdf = Ostap.Math.ECDF ( data2vct ( data ) )

        ## evaluate CDF for sorted data 
        cdf_data = vct_clip ( self.__vct_cdf ( data ) ) 
        
        self.__estimators = {
            'KS'  : kolmogorov_smirnov ( cdf_data ) , 
            'K'   : kuiper             ( cdf_data ) , 
            'AD'  : anderson_darling   ( cdf_data ) , 
            'CM'  : cramer_von_mises   ( cdf_data ) , 
            'ZK'  : ZK                 ( cdf_data ) , 
            'ZA'  : ZA                 ( cdf_data ) , 
            'ZC'  : ZC                 ( cdf_data ) ,
        }
        
    ## serialize the object 
    def __getstate__ ( self ) :
        """ Serialize the object"""
        self.pdf.load_params ( self.__parameters , silent = True )
        return { 'pdf'        : self.pdf          ,
                 'parameters' : self.__parameters , 
                 'cdf'        : self.cdf          ,
                 'ecdf'       : self.ecdf         ,
                 'estimators' : self.estimators   , 
                 'store'      : self.__store      }
    
    ## De-serialize the object 
    def __setstate__ ( self , state ) :
        """ De-serialize the object """
        print ( 'I AM SET STATE' , state.keys() ) 
        self.__pdf        = state.pop ( 'pdf'        ) 
        self.__parameters = state.pop ( 'parameters' , {}  ) 
        self.__cdf        = state.pop ( 'cdf'        )
        self.__ecdf       = state.pop ( 'ecdf'       )
        self.__estimators = state.pop ( 'estimators' ) 
        self.__store      = state.pop ( 'store' , () )
        ## vectorized form of CDF 
        self.__vct_cdf    = numpy.vectorize ( self.cdf )
        self.pdf.load_params ( self.__parameters , silent = True )
                               
    # =========================================================================
    @property
    def estimators ( self ) :
        """`estimators` : get a dictionary of all estimators"""
        return self.__estimators
    
    # =========================================================================
    ## size of dataset
    @property 
    def N  ( self ) :
        """`N` : size of dataset"""
        return len ( self.__ecdf ) 
    
    # =========================================================================
    ## fit-PDF 
    @property
    def pdf ( self ) :
        """`pdf` : PDF for fit-function"""
        return self.__pdf
    
    # =========================================================================
    ## fit-CDF
    @property
    def cdf ( self ) :
        """`cdf` : CDF for fit-function"""
        return self.__cdf

    # =========================================================================
    ## vectorized CDF
    @property 
    def vcdf ( self ) :
        """`vcdf` : vectorized fit-CDF"""
        return self.__vct_cdf 
    
    # =========================================================================
    ## empirical CDF for data
    @property 
    def ecdf ( self ) :
        """`ecdf` : empirical CDF for data"""
        return self.__ecdf

    # =========================================================================
    @property
    def parameters ( self ) :
        """ PDF parameters """
        return self.__parameters
    
    # =========================================================================
    ## Get Kolmogorov-Smirnov statistics 
    @property 
    def kolmogorov_smirnov_estimator ( self ) :
        """ Get Kolmogorov-Smirnov' statistics KS
        """
        return self.__estimators.get('KS',None) 

    # =========================================================================
    ## Get Kuiper  statististics 
    @property 
    def kuiper_estimator ( self ) :
        """ Get Kuiper' statistics K 
        """        
        return self.__estimators.get('K', None ) 
                
    # =========================================================================
    ## Get Anderson-Darling  statistiscs 
    @property 
    def anderson_darling_estimator ( self ) :
        """ Get Anderson-Darling statistiscs AD
        """
        return self.__estimators.get ( 'AD' , None ) 

    # =========================================================================
    ## Get Cramer-von Mises statistics 
    @property 
    def cramer_von_mises_estimator ( self ) :
        """ Get Cramer-von Mises statistics 
        """
        return self.estimators.get ( 'CM' , None ) 
    
    # =========================================================================
    ## Get ZK statististics 
    @property 
    def ZK_estimator  ( self ) :
        """ Get ZK statistics
        """
        return self.estimators.get ( 'ZK' , None ) 
        
    # =========================================================================
    ## Get ZA statististics 
    @property 
    def ZA_estimator  ( self ) :
        """ Get ZA statistics
        """        
        return self.estimators.get( 'ZA' , None ) 
    
    # =========================================================================
    ## Get ZC statististics 
    @property 
    def ZC_estimator ( self ) :
        """ Get ZC statistics
        """        
        return self.__estimators.get ( 'ZC' , None ) 

    # =========================================================================
    ## Print the summary as Table
    def table ( self             , * , 
                title     = ''   ,
                prefix    = ''   ,
                precision = 4    , 
                width     = 6    ,
                style     = None ) :
        """ Print the summary as Table
        """
        ##
        header = ( 'Statistics' , 'Value' , '' ) 
        rows   = [] 
        
        for label , value  in loop_items ( self.estimators ) :
            
            the_label = Labels.get ( label , label )

            result , expo = pretty_float ( value , width = width , precision = precision )

            if expo : row = the_label , result , '10^%+d' % expo
            else    : row = the_label , result 
            rows.append ( row )

        rows = [ header ] + sorted ( rows ) 
        title = title if title else 'Goodness of 1D Fit' 
        rows  = T.remove_empty_columns ( rows )
        return T.table ( rows , title = title , prefix = prefix , alignment = 'lcl' , style = style  )

    ## print estimator as table 
    __repr__ = table
    __str__  = table

    # =========================================================================
    ## Draw fit CDF & empirical ECDF 
    def draw  ( self , opts = '' , *args , **kwargs ) :
        """ Draw fit CDF & empirical CDF
        """
        cdf  = self.__cdf
        ecdf = self.__ecdf
        
        xmin, xmax =  ecdf.xmin () , ecdf.xmax ()
        
        if hasattr ( cdf , 'xminmax' ) :
            xmn , xmx =  cdf.xminmax()
            xmin = min ( xmin , xmn )
            xmax = max ( xmax , xmx )
        else :
            if hasattr ( cdf , 'xmin' ) : xmin = min ( xmin , cdf.xmin () )
            if hasattr ( cdf , 'xmax' ) : xmax = max ( xmax , cdf.xmax () )

        xmin , xmax = axis_range ( xmin , xmax , delta = 0.20 )
        
        xmin = kwargs.pop ( 'xmin' , xmin )
        xmax = kwargs.pop ( 'xmax' , xmax )

        opts    = opts.strip() 
        optsame = 'same' if not opts else '%s %s' % ( 'same' , opts ) 

        if isinstance ( cdf , AFUN1 ) :
            self.__frame = cdf.draw ()
            self.__frame.draw ( opts )
            r1 = self.__frame 
        else : 
            self.__draw_fun = lambda x : cdf ( x ) 
            r1 = f1_draw   ( self.__draw_fun , opts , color = ROOT.kOrange + 1 , linewidth = 3 , xmin = xmin , xmax = xmax , **kwargs ) 
        
        r2 = ecdf.draw ( optsame  , color = 2 , linewidth = 3 , xmin = xmin , xmax = xmax , **kwargs )
        return r1 , r2
    
# =============================================================================
## @class GoF1DToys
#  Check Goodness of 1D-Fits using toys 
class GoF1DToys(GoF1D) :
    """ Check Goodness-of-fit with toys (1D-case) 
    """
    ## result of GoF-toys 
    Result = namedtuple ( 'Result' , 'statistics counter pvalue nsigma' )
    # =========================================================================
    ## Initialize GoF1D toys object :
    #  @code
    #  gof  = GoF1D     ( ... ) 
    #  toys = GoF1DToys ( gof )
    #  @endcode 
    def __init__ ( self , gof ) :
        """ Initialize GoF1D toys object :
        >>> gof  = GoF1D     ( ... ) 
        >>> toys = GoF1DToys ( gof ) 
        """
        assert isinstance ( gof , GoF1D ) , "Invalid `gof`-parameter"

        ## mimic the copy-constructor for the base class 
        state = GoF1D.__getstate__ ( gof ) 
        GoF1D.__setstate__ ( self , state )

        self.__counters = defaultdict(SE) 
        self.__ecdfs    = {}        
        self.__nToys    = 0
        
    ## serialize the object 
    def __getstate__ ( self ) :
        """ Serialize the object 
        """
        #
        ## (1) serialize the base 
        state = GoF1D.__getstate__ ( self )
        # 
        state [ 'counters' ] = self.__counters
        state [ 'ecdfs'    ] = self.__ecdfs 
        state [ 'nToys'    ] = self.__nToys
        # 
        return state 
    
    ## De-serialize the object 
    def __setstate__ ( self , state ) :
        """ De-serialize the object """
        
        ## (1) de-serialize the base 
        GoF1D.__setstate__ ( self , state )
        # 
        self.__counters   = state.pop ( 'counters'   )
        self.__ecdfs      = state.pop ( 'ecdfs'      )
        self.__nToys      = state.pop ( 'nToys'    , 0  )

    # ===============================================================================
    ## run toys 
    def run ( self , nToys = 1000 , parallel = False , silent = False , nSplit = 0 ) :
        """ Run toys 
        """
        assert isinstance ( nToys , int ) and 0 < nToys , "Invalid `nToys` argument!"

        if parallel :
            from ostap.parallel.parallel_gof1d import parallel_gof1dtoys as parallel_toys 
            self += parallel_toys ( gof      = self       ,
                                    nToys    = nToys      ,
                                    nSplit   = nSplit     ,
                                    silent   = True       ,
                                    progress = not silent )
            return self 
        
        varname  = self.pdf.xvar.name        
        results  = defaultdict(list)
        counters = self.counters 
        vct_cdf  = self.vcdf

        from ostap.utils.progress_bar import progress_bar

        cnt = SE() 
        for i in progress_bar ( nToys , silent = silent , description = 'Toys:') :

            self.pdf.load_params ( self.parameters , silent = True )            
            dset     = self.pdf.generate ( self.N  , sample = True )                
            data     = dset.tonumpy ( varname ) [ varname ] 
            data     = numpy.sort ( data )

            ## CLIP... does one need it? 
            cdf_data = vct_clip ( vct_cdf ( data ) )
            
            ks       = kolmogorov_smirnov ( cdf_data )
            k        = kuiper             ( cdf_data )
            ad       = anderson_darling   ( cdf_data )
            cm       = cramer_von_mises   ( cdf_data )
            zk       = ZK                 ( cdf_data )
            za       = ZA                 ( cdf_data )
            zc       = ZC                 ( cdf_data )

            counters [ 'KS' ] += ks
            counters [ 'K'  ] += k
            counters [ 'AD' ] += ad
            counters [ 'CM' ] += cm
            counters [ 'ZK' ] += zk
            counters [ 'ZA' ] += za
            counters [ 'ZC' ] += zc
            
            results  [ 'KS' ].append ( ks )    
            results  [ 'K'  ].append ( k  )
            results  [ 'AD' ].append ( ad ) 
            results  [ 'CM' ].append ( cm ) 
            results  [ 'ZK' ].append ( zk ) 
            results  [ 'ZA' ].append ( za ) 
            results  [ 'ZC' ].append ( zc ) 

            cnt += ks 
            ## delete data
            if isinstance ( dset , ROOT.RooDataSet ) :
                dset.clear () 
                dset = Ostap.MoreRooFit.delete_data ( dset )
                
            del dset
            del data
            del cdf_data            

        ## accumulate number of toys 
        self.__nToys += nToys 

        ECDF = Ostap.Math.ECDF 
        for key in results :
            data = results [ key ]
            if not data : continue
            if not key in self.__ecdfs : self.__ecdfs [ key ] = ECDF ( data , True ) ## complementary ECDF!
            else                       : self.__ecdfs [ key ]  .add  ( data2vct ( data ) ) 

        del results 
        return self

    # =========================================================================
    ## number of toys 
    @property
    def nToys ( self ) :
        """`nToys` : number of toys"""
        return self.__nToys
    
    # =========================================================================
    ## ECDFs
    @property 
    def ecdfs ( self ) :
        """`ecdfs` : toy results as empirical cumulative distribution functions"""
        return self.__ecdfs
    # =========================================================================
    ## Counters  
    @property 
    def counters ( self ) :
        """`counters` : toy results as counters"""
        return self.__counters

    # =========================================================================
    ## Get Kolmogorov-Smirnov statistics 
    @property 
    def kolmogorov_smirnov ( self ) :
        """ Get Kolmogorov-Smirnov statistiscs KS
        """
        return self.result ( 'KS' ) 
    
    # ===============================================
    ## Get Anderson-Darling  statistics
    @property 
    def anderson_darling ( self ) :
        """ Get Anderson-Darling statistiscs 
        """
        return self.result ( 'AD' ) 
        
    # =========================================================================
    ## Get Cramer-von Mises statistics 
    @property 
    def cramer_von_mises ( self ) :
        """ Get Cramer-von Mises statistics 
        """
        return self.result ( 'CM' ) 
        
    # =========================================================================
    ## Get ZK statistics
    @property 
    def ZK  ( self ) :
        """ Get ZK statistics 
        """
        return self.result ( 'ZK' ) 
        
    # =========================================================================
    ## Get ZA statistics
    @property 
    def ZA  ( self ) :
        """ Get ZA statistics
        """        
        return self.result ( 'ZA' ) 
    
    # =========================================================================
    ## Get ZC statistics
    @property 
    def ZC  ( self ) :
        """ Get ZC statistics 
        """        
        return self.result ( 'ZC' ) 

    # =========================================================================
    ## Get Kuiper statistics
    @property 
    def kuiper ( self ) :
        """ Get Kuiper statistics 
        """        
        return self.result ( 'K' ) 

    # =========================================================================
    ## merge two objects:
    def merge ( self , other ) :
        self += other
        return self
    
    ## merge two objects:
    def __iadd__ ( self , other ) :
        """ Merge two GoF-toys objects """        
        if not isinstance ( other , GoF1DToys ) : return NotImplemented 

        ecdfs = other.ecdfs        
        for key, ecdf    in loop_items ( ecdfs   ) : 
            if key in self.__ecdfs    : self.__ecdfs    [ key ] += ecdf
            else                      : self.__ecdfs    [ key ]  = ecdf
            
        counters = other.counters
        for key, counter in loop_items ( counters ) : 
            if key in self.__counters : self.__counters [ key ] += counter 
            else                      : self.__counters [ key ]  = counter 

        self.__nToys += other.nToys
        ##

        return self 

    # ========================================================================
    ## Helper method to get the result
    def result ( self , label ) :
        """ Helper method to get the result 
        """
        if not label in self.estimators : return None
        if not label in self.ecdfs      : return None
        if not label in self.counters   : return None
        ##
        value   = self.estimators   [ label ]
        ecdf    = self.ecdfs        [ label ] 
        counter = self.counters     [ label ] 
        ##
        pvalue = ecdf. estimate ( value  ) ## estimate the p-value
        #
        clip = 0.5 * pvalue.error()
        pv   = pvalue 
        if   1 <= pvalue.value() : pv = VE ( 1 - clip , pv.cov2() )
        elif 0 >= pvalue.value() : pv = VE (     clip , pv.cov2() )
        nsigma = significance ( pv ) ## convert  it to significace
        
        return self.Result ( value   ,
                             counter ,
                             pvalue  ,
                             nsigma  )
    
    # =========================================================================
    ## format a row in the summary table
    def row  ( self , what , result , width = 6 , precision = 4 ) :
        """ Format a row in the sumamry table
        """
        value      = result.statistics
        counter    = result.counter
        pvalue     = result.pvalue
        nsigma     = result.nsigma
        
        mean       = counter.mean   ()
        rms        = counter.rms    () 
        vmin, vmax = counter.minmax () 
        
        mxv = max ( abs ( value        ) ,
                    abs ( mean.value() ) ,
                    mean.error()         , rms ,
                    abs ( vmin )  , abs ( vmax ) ) 
        
        fmt, fmtv , fmte , expo = fmt_pretty_ve ( VE ( mxv ,  mean.cov2() ) ,
                                                  width       = width       ,
                                                  precision   = precision   , 
                                                  parentheses = False       )
        
        if expo : scale = 10**expo
        else    : scale = 1
        
        fmt2 = '%s/%s' % ( fmtv , fmtv ) 

        vs  = value / scale
        vm  = mean  / scale
        vr  = rms   / scale
        vmn = vmin  / scale
        vmx = vmax  / scale
        
        pvalue = str ( ( 100 * pvalue ) .toString ( '%% 5.2f %s %%-.2f' % plus_minus ) )
        sigma  = str ( nsigma.toString ( '%%.2f %s %%-.2f' % plus_minus ) if float ( nsigma ) < 1000 else '+inf' ) 

        return ( what  ,
                 fmtv  % vs ,
                 fmt   % ( vm.value() , vm.error() ) ,
                 fmtv  % vr                          ,
                 fmt2  %  ( vmn , vmx )              ,
                 ( '%s10^%+d' %  ( times , expo )  if expo else '' ) , pvalue , sigma )

    # =========================================================================
    ## Make a summary table
    def table ( self , title = '' , prefix = '' , width = 6 , precision = 4 , style = None ) :
        """ Make a summary table
        """
        import ostap.logger.table  as     T                 
        header = ( 'Statistics' , 'value' , 'mean' , 'rms' , 'min/max' , 'factor' , 'p-value [%]' , '#sigma' )
        rows   = [] 
        
        for label in self.ecdfs :
            
            result  = self.result ( label )
            if not result : continue

            the_label = Labels.get ( label , label )
            row = self.row ( the_label , result , width = width , precision = precision )
            rows.append ( row ) 
                    
        if rows : rows = T.remove_empty_columns ( rows )
        if   not title and self.nToys : title = 'Goodness of 1D-fit with #%d toys' % self.nToys  
        elif not title                : title = 'Goodness of 1D-fit'

        rows = [ header ] + sorted ( rows ) 
        rows = T.remove_empty_columns ( rows ) 
        return T.table ( rows , title = title , prefix = prefix , alignment = 'lccccccc' , style = style )

    __repr__ = table
    __str__  = table

    # =========================================================================
    ## Draw ECDF for toys & statistical estimator 
    def draw  ( self , what , opts = '' , *args , **kwargs ) :
        """ Draw ECDF for toys & statistical estgimator 
        """
        key = cidict_fun ( what ) 
        if   key in KS_keys and 'KS' in self.ecdfs :             
            result = self.result ( 'KS' )
            ecdf   = self.ecdfs  [ 'KS' ]
            logger.info ( 'Toy resuls for Kolmogorov-Smirnov estimate' ) 
        elif key in K_keys  and 'K'  in self.ecdfs : 
            result = self.result ( 'K' )
            ecdf   = self.ecdfs  [ 'K' ]
            logger.info ( 'Toy resuls for Kuiper estimate' ) 
        elif key in AD_keys and 'AD' in self.ecdfs :             
            result = self.result ( 'AD' )
            ecdf   = self.ecdfs  [ 'AD' ]
            logger.info ( 'Toy resuls for Anderson-Darling estimate' ) 
        elif key in CM_keys  and 'CM' in self.ecdfs : 
            result = self.result  ( 'CM' )
            ecdf   = self.ecdfs   [ 'CM' ]
            logger.info ( 'Toy resuls for Cramer-von Mises  estimate' ) 
        elif key in ZK_keys  and 'ZK' in self.ecdfs : 
            result = self.result  ( 'ZK' )
            ecdf   = self.ecdfs   [ 'ZK' ]
            logger.info ( 'Toy resuls for Zhang/ZK estimate' ) 
        elif key in ZA_keys  and 'ZA' in self.ecdfs :  
            result = self.result  ( 'ZA' )
            ecdf   = self.ecdfs   [ 'ZA' ]
            logger.info ( 'Toy resuls for Zhang/ZA estimate' ) 
        elif key in ZC_keys and 'ZC' in self.ecdfs : 
            result = self.result  ( 'ZC' )
            ecdf   = self.ecdfs   [ 'ZC' ]
            logger.info ( 'Toy resuls for Zhang/ZC estimate' ) 
        else :
            raise KeyError (  "draw: Invalid `what`:%s" % what )
            
        xmin , xmax = ecdf.xmin () , ecdf.xmax ()
        value     = result.statistics
        xmin      = min ( xmin , value )
        xmax      = max ( xmax , value )
        xmin , xmax  = axis_range ( xmin , xmax , delta = 0.20 )

        kwargs [ 'xmin' ] = kwargs.get ( 'xmin' , xmin ) 
        kwargs [ 'xmax' ] = kwargs.get ( 'xmax' , xmax )

        result    = ecdf.draw  ( opts , *args , **kwargs ) 
        line      = ROOT.TLine ( value , 1e-3 , value , 1 - 1e-3 )
        ## 
        line.SetLineWidth ( 4 ) 
        line.SetLineColor ( 8 ) 
        line.draw ( 'same' )
        ##
        self._line = line
        return result, line  

# =============================================================================
## Goodness-of-Fit for SimFit
# =============================================================================

# =============================================================================
## @class GoFSimFit
#  Goodness-of-fit for simultaneous 1D-fits
#  - All components of the simultaneous fit must be 1D-components 
#  - GoF is estimated for each 1D-component
class GoFSimFit(object) :
    """ Goodness-of-fit for Simultaneous 1D-fits
    - All components of the Simultaneous fit must be 1D-components!
    - GoF is estimated for each 1D-component
    """
    def __init__ ( self               ,
                   pdf                ,
                   dataset            ,
                   parameters  = None ) :

        from ostap.fitting.simfit import SimFit 
        assert isinstance ( pdf     , SimFit          ) , 'Invalid type of `pdf`:%s'     % typename ( pdf     )
        assert isinstance ( dataset , ROOT.RooDataSet ) , 'Invalid type of `dataset`:%s' % typename ( dataset )

        assert pdf.sample in dataset , "Sample category `%s` not in dataset!" % pdf.sample.name 

        self.__pdf        = pdf
        self.__parameters = parameters if parameters else pdf.parameters ( dataset )

        ## load parameters here
        self.__pdf.load_params ( self.__parameters , silent = True ) 
        
        self.__gofs = {}
        self.__N    = {}

        name = self.sample.name
        for key , cmp  in pdf.categories.items ()  :
            
            assert isinstance ( cmp , PDF1 ) , "Component `%s` is not PDF1` %s" % ( key , typename ( cmp ) )
            obs      = cmp.pdf.getObservables ( dataset )
            category = '%s==%s::%s' % ( name , name , key ) 
            ds       = dataset.subset ( variables = obs ,cuts = category )
            gof      = GoF1D ( cmp , ds )
            
            self.__gofs [ key ] = gof 
            self.__N    [ key ] = len ( ds ) 
            
    @property
    def pdf ( self ) :
        """`pdf`: SimFit/PDF for simultaneous fit
        """
        return self.__pdf
    
    @property
    def sample ( self ) :
        """sample`: sample/category  variable for simultaneous fit
        """
        return self.pdf.sample 

    @property
    def parameters ( self ) :
        """`parameters' : fit parameters, e.g. fit-resutl or dictinoary or ...
        """
        return self.__parameters 
        
    @property
    def gofs ( self ) :
        """`gofs` : individual GoF estimators for Simfit components
        """
        return self.__gofs 

    @property
    def N ( self ) :
        """`N` : dictionary { sample : #events}
        """
        return { key : v.N  for key , v in self.__gofs.items () }
    
    # =========================================================================
    ## all estimatrs togather
    @property 
    def estimators ( self ) :
        """`estimators` : get all statistical estimators
        """
        return { key : v.estimators for key , v in self.__gofs.items () }
        
    @property
    def kolmogorov_smirnov_estimator ( self ) :
        """ Get Kolmogorov-Smirnov statistics' KS
        """
        return { key : v.kolmogorov_smirnov_estimator for key , v in self.__gofs.items () }
    
    @property
    def kuiper_estimator ( self ) :
        """ Get Kuiper' statistics K
        """
        return { key : v.kuiper_estimator for key , v in self.__gofs.items () }

    @property
    def anderson_darling_estimator ( self ) :
        """ Get Anderson-Darling' statistics AD
        """
        return { key : v.anderson_darling_estimator for key , v in self.__gofs.items () }
    
    @property
    def cramer_von_mises_estimator ( self ) :
        """ Get Cramer-on Mises' statistics CM
        """
        return { key : v.cramer_von_mises_estimator for key , v in self.__gofs.items () }
    
    @property
    def ZK_estimator ( self ) :
        """ Get Zhang's statistics ZK
        """
        return { key : v.ZK_estimator for key , v in self.__gofs.items () }
    
    @property
    def ZA_estimator ( self ) :
        """ Get Zhang's statistics ZA
        """
        return { key : v.ZA_estimator for key , v in self.__gofs.items () }
    
    @property
    def ZC_estimator ( self ) :
        """ Get Zhang's statistics ZC
        """
        return { key : v.ZC_estimator for key , v in self.__gofs.items () }
   
    # ====================================================================================
    ## Print the summary as Table  (for simfit)
    def table ( self             , * ,
                title     = ''   ,
                prefix    = ''   ,
                precision = 4    , 
                width     = 6    ,
                style     = None ) :
        """ Print the summary as Table  (for SimFit)
        """
        ##
        
        keys = tuple ( self.gofs.keys() )
        
        header = "Statistics",  
        for k in keys : header += ( str ( k ) , '' )
        
        rows = []

        ## get the estimators from the fist GOF 
        estimators = self.gofs [ keys [ 0 ] ].estimators.keys()
        for label in estimators : 
            
            the_label = Labels.get ( label , label )

            if   'KS' == label : values = self.kolmogorov_smirnov_estimator
            elif 'K'  == label : values = self.kuiper_estimator
            elif 'AD' == label : values = self.anderson_darling_estimator
            elif 'CM' == label : values = self.cramer_von_mises_estimator
            elif 'ZK' == label : values = self.ZK_estimator
            elif 'ZA' == label : values = self.ZA_estimator
            elif 'ZC' == label : values = self.ZC_estimator
            else               : continue
            
            row = the_label ,
            
            for k, v in values.items () :
                result , expo = pretty_float ( v  , width = width , precision = precision )                
                if expo : row += ( result , '10^%+d' % expo )
                else    : row += ( result , ''              )
                
            rows.append ( row )

        rows  = [ header ] + sorted ( rows ) 
        title = title if title else 'Goodness of 1D (Sim)Fit'
        rows  = T.remove_empty_columns ( rows )
        return T.table ( rows , title = title , prefix = prefix , alignment = 'lccccccccccccccccccc' , style = style  )

    ## print estimators as table 
    __repr__ = table
    __str__  = table

    ## Draw fit CDF & empirical ECDF 
    def draw  ( self , sample , opts = '' , *args , **kwargs ) :
        """ Draw fit CDF & empirical CDF
        """
        gof = self.__gofs.get ( sample , None )
        if gof is None : raise KeyError ( "Invaild sample `%s`" % sample )
        return gof.draw ( opts , *args , **kwargs )
    
    ## serialize the object 
    def __getstate__ ( self ) :
        """ Serialize the object"""
        self.__pdf.load_params ( self.__parameters , silent = True )
        return { 'pdf'        : self.pdf    ,
                 'parameters' : self.__parameters , 
                 'gofs'       : self.__gofs }
    
    ## De-serialize the object 
    def __setstate__ ( self , state ) :
        """ De-serialize the object """         
        self.__pdf        = state.pop ( 'pdf'  )
        self.__parameters = state.pop ( 'parameters' , {} ) 
        self.__gofs       = state.pop ( 'gofs' )
        self.__pdf.load_params ( self.__parameters , silent = True )
        
# =============================================================================
## @class GoFSimFitToys
#  Check Goodness of 1D (Sim)Fits using toys 
class GoFSimFitToys(GoFSimFit) :
    """ Check Goodness-of-Fit with toys (Simfit caase)
    """
    ## result of GoF-toys 
    Result = namedtuple ( 'Result' , 'statistics counter pvalue nsigma' )
    # =========================================================================
    ## Initialize GoF1D toys object :
    #  @code
    #  gof  = GoFSimfit     ( ... ) 
    #  toys = GoFSimFitToys ( gof )
    #  @endcode 
    def __init__ ( self , gof ) :
        """ Initialize GoF1D toys object :
        >>> gof  = GoFSimFit     ( ... ) 
        >>> toys = GoFSimFitToys ( gof ) 
        """
        assert isinstance ( gof , GoFSimFit ) , "Invalid `gof`-parameter"

        ## mimic the copy-constructor for the base class 
        state = GoFSimFit.__getstate__ ( gof ) 
        GoFSimFit.__setstate__ ( self , state )

        self.__counters = { k : defaultdict(SE) for k in self.gofs }
        self.__ecdfs    = { k : {}              for k in self.gofs }
        self.__total    = defaultdict(EffCounter) 
        self.__nToys    = 0

    ## serialize the object 
    def __getstate__ ( self ) :
        """ Serialize the object 
        """
        #
        ## (1) serialize the base 
        state = GoFSimFit.__getstate__ ( self )
        # 
        state [ 'counters' ] = self.__counters
        state [ 'ecdfs'    ] = self.__ecdfs 
        state [ 'total'    ] = self.__total 
        state [ 'nToys'    ] = self.__nToys
        # 
        return state 
    
    ## De-serialize the object 
    def __setstate__ ( self , state ) :
        """ De-serialize the object """
        
        ## (1) de-serialize the base 
        GoFSimFit.__setstate__ ( self , state )
        # 
        self.__counters   = state.pop ( 'counters'  )
        self.__ecdfs      = state.pop ( 'ecdfs'     )
        self.__total      = state.pop ( 'total'     )
        self.__nToys      = state.pop ( 'nToys' , 0 )

    # ===============================================================================
    ## run toys 
    def run ( self , nToys = 1000 , parallel = False , silent = False , nSplit = 0 ) :
        """ Run toys 
        """
        assert isinstance ( nToys , int ) and 0 < nToys , "Invalid `nToys` argument!"

        ## 
        ## if parallel :
        ##    from ostap.parallel.parallel_gof1d import parallel_gof1dtoys as parallel_toys 
        ##    self += parallel_toys ( gof      = self       ,
        ##                            nToys    = nToys      ,
        ##                            nSplit   = nSplit     ,
        ##                            silent   = True       ,
        ##                            progress = not silent )
        ##    return self 

        
        results  = { k : defaultdict(list) for k in self.gofs } 
        counters = self.counters 

        from ostap.utils.progress_bar import progress_bar

        for i in progress_bar ( nToys , silent = silent , description = 'Toys:') :

            self.pdf.load_params ( self.parameters , silent = True )            
            dset     = self.pdf.generate ( self.N  , sample = True )
            sample   = self.pdf.sample

            total_KS = True  
            total_K  = True 
            total_AD = True 
            total_CM = True 
            total_ZK = True 
            total_ZA = True 
            total_ZC = True 
                      
            for key , gof  in self.gofs.items ()  :

                obs      = gof.pdf.pdf.getObservables ( dset )
                category = '%s==%s::%s' % ( sample.name , sample.name , key )
                ds       = dset.subset ( variables = obs , cuts = category )
                varname  = gof.pdf.xvar.name

                ## convert to numpy & sort it! 
                data     = ds.tonumpy ( varname ) [ varname ]
                data     = numpy.sort ( data )
                
                vct_cdf  = gof.vcdf

                ## CLIP... does one need it? 
                cdf_data = vct_clip ( vct_cdf ( data ) )
            
                ks       = kolmogorov_smirnov ( cdf_data )
                k        = kuiper             ( cdf_data )
                ad       = anderson_darling   ( cdf_data )
                cm       = cramer_von_mises   ( cdf_data )
                zk       = ZK                 ( cdf_data )
                za       = ZA                 ( cdf_data )
                zc       = ZC                 ( cdf_data )
                
                total_KS = total_KS and gof .kolmogorov_smirnov_estimator <= ks
                total_K  = total_K  and gof .            kuiper_estimator <= k
                total_AD = total_AD and gof .  anderson_darling_estimator <= ad
                total_CM = total_CM and gof .  cramer_von_mises_estimator <= cm
                total_ZK = total_ZK and gof .                ZK_estimator <= zk 
                total_ZA = total_ZK and gof .                ZA_estimator <= za 
                total_ZC = total_ZK and gof .                ZC_estimator <= zc 
                
                cnts = counters [ key ]                
                cnts [ 'KS' ] += ks
                cnts [ 'K'  ] += k
                cnts [ 'AD' ] += ad
                cnts [ 'CM' ] += cm
                cnts [ 'ZK' ] += zk
                cnts [ 'ZA' ] += za
                cnts [ 'ZC' ] += zc

                res  = results [ key ]
                res [ 'KS'  ].append ( ks )    
                res [ 'K'   ].append ( k  )
                res [ 'AD'  ].append ( ad ) 
                res [ 'CM'  ].append ( cm ) 
                res [ 'ZK'  ].append ( zk ) 
                res [ 'ZA'  ].append ( za ) 
                res [ 'ZC'  ].append ( zc ) 

                ## delete data
                if isinstance ( ds , ROOT.RooDataSet ) : ds.clear () 
                del ds
                del data
                del cdf_data
                
            self.__total [ 'KS' ] += total_KS
            self.__total [ 'K'  ] += total_K
            self.__total [ 'AD' ] += total_AD
            self.__total [ 'CM' ] += total_CM
            self.__total [ 'ZK' ] += total_ZK
            self.__total [ 'ZA' ] += total_ZA
            self.__total [ 'ZC' ] += total_ZC
                
            ## delete data
            if isinstance ( dset , ROOT.RooDataSet ) : dset.clear () 
            del dset
            
        ## accumulate number of toys 
        self.__nToys += nToys 

        ECDF = Ostap.Math.ECDF 
        for key, vv  in results.items()  :
            for e , data in vv.items() : 
                if not data : continue
                ecdfs = self.__ecdfs [ key ]
                if not e in ecdfs : ecdfs [ e ] = ECDF ( data , True ) ## complementary ECDF!
                else              : ecdfs [ e ] .add   ( data2vct ( data ) ) 
                
        del results 
        return self
    
    # =========================================================================
    ## number of toys 
    @property
    def nToys ( self ) :
        """`nToys` : number of toys"""
        return self.__nToys
    
    # =========================================================================
    ## ECDFs
    @property 
    def ecdfs ( self ) :
        """`ecdfs` : toy results as empirical cumulative distribution functions"""
        return self.__ecdfs
    # =========================================================================
    ## Counters  
    @property 
    def counters ( self ) :
        """`counters` : toy results as counters"""
        return self.__counters

    # =========================================================================
    ## Total/global counters
    @property 
    def total ( self ) :
        """`total` : total/global counters
        """
        return self.__total
    
    # =========================================================================
    ## Get Kolmogorov-Smirnov statistics 
    @property 
    def kolmogorov_smirnov ( self ) :
        """ Get Kolmogorov-Smirnov statistiscs KS
        """
        return self.result ( 'KS' ) 
    
    # ===============================================
    ## Get Anderson-Darling  statistics
    @property 
    def anderson_darling ( self ) :
        """ Get Anderson-Darling statistiscs 
        """
        return self.result ( 'AD' ) 
        
    # =========================================================================
    ## Get Cramer-von Mises statistics 
    @property 
    def cramer_von_mises ( self ) :
        """ Get Cramer-von Mises statistics 
        """
        return self.result ( 'CM' ) 
        
    # =========================================================================
    ## Get ZK statistics
    @property 
    def ZK  ( self ) :
        """ Get ZK statistics 
        """
        return self.result ( 'ZK' ) 
        
    # =========================================================================
    ## Get ZA statistics
    @property 
    def ZA  ( self ) :
        """ Get ZA statistics
        """        
        return self.result ( 'ZA' ) 
    
    # =========================================================================
    ## Get ZC statistics
    @property 
    def ZC  ( self ) :
        """ Get ZC statistics 
        """        
        return self.result ( 'ZC' ) 

    # =========================================================================
    ## Get Kuiper statistics
    @property 
    def kuiper ( self ) :
        """ Get Kuiper statistics 
        """        
        return self.result ( 'K' ) 

    # ========================================================================
    ## Helper method to get the result
    def result ( self , sample , label ) :
        """ Helper method to get the result 
        """
        ##
        value   = self.estimators   [ sample ] [ label ]
        ecdf    = self.ecdfs        [ sample ] [ label ] 
        counter = self.counters     [ sample ] [ label ] 
        ##
        pvalue = ecdf. estimate ( value  ) ## estimate the p-value
        #
        clip = 0.5 * pvalue.error()
        pv   = pvalue 
        if   1 <= pvalue.value() : pv = VE ( 1 - clip , pv.cov2() )
        elif 0 >= pvalue.value() : pv = VE (     clip , pv.cov2() )
        nsigma = significance ( pv ) ## convert  it to significace
        
        return self.Result ( value   ,
                             counter ,
                             pvalue  ,
                             nsigma  )

    # =========================================================================
    ## format a row in the summary table
    def row  ( self , sample , what , result , width = 6 , precision = 4 ) :
        """ Format a row in the sumamry table
        """
        value      = result.statistics
        counter    = result.counter
        pvalue     = result.pvalue
        nsigma     = result.nsigma
        
        mean       = counter.mean   ()
        rms        = counter.rms    () 
        vmin, vmax = counter.minmax () 
        
        mxv = max ( abs ( value        ) ,
                    abs ( mean.value() ) ,
                    mean.error()         , rms ,
                    abs ( vmin )  , abs ( vmax ) ) 
        
        fmt, fmtv , fmte , expo = fmt_pretty_ve ( VE ( mxv ,  mean.cov2() ) ,
                                                  width       = width       ,
                                                  precision   = precision   , 
                                                  parentheses = False       )
        
        if expo : scale = 10**expo
        else    : scale = 1
        
        fmt2 = '%s/%s' % ( fmtv , fmtv ) 

        vs  = value / scale
        vm  = mean  / scale
        vr  = rms   / scale
        vmn = vmin  / scale
        vmx = vmax  / scale
        
        pvalue = str ( ( 100 * pvalue ) .toString ( '%% 5.2f %s %%-.2f' % plus_minus ) )
        sigma  = str ( nsigma.toString ( '%%.2f %s %%-.2f' % plus_minus ) if float ( nsigma ) < 1000 else '+inf' ) 

        return ( what ,
                 sample ,
                 fmtv   % vs ,
                 fmt    % ( vm.value() , vm.error() ) ,
                 fmtv   % vr                          ,
                 fmt2   %  ( vmn , vmx )              ,
                 ( '%s10^%+d' %  ( times , expo )  if expo else '' )   ,                  
                 pvalue , sigma )
    
    # =========================================================================
    ## Make a summary table
    def table ( self , title = '' , prefix = '' , width = 6 , precision = 4 , style = None ) :
        """ Make a summary table
        """
        import ostap.logger.table  as     T                 
        header = ( 'Statistics'  ,
                   'Sample'      , 
                   'value'       , 
                   'mean'        ,
                   'rms'         ,
                   'min/max'     ,
                   'factor'      ,
                   'p-value [%]' ,
                   '#sigma'      ) 
        
        rows = []
        for sample , ecdfs in self.ecdfs.items()  :
            for label in ecdfs :
                result  = self.result ( sample , label )
                if not result : continue
                the_label = Labels.get ( label , label )
                row = self.row ( sample , the_label , result , width = width , precision = precision )
                rows.append ( row ) 
                    
        if   not title and self.nToys : title = 'Goodness of 1D-fit with #%d toys' % self.nToys  
        elif not title                : title = 'Goodness of 1D-fit'

        for e , cnt in self.total.items() :

            ## get the binomial efficiency 
            pvalue = cnt.efficiency

            clip   = 0.5 * pvalue.error()
            pv     = pvalue 
            if   1 <= pvalue.value() : pv = VE ( 1 - clip , pv.cov2() )
            elif 0 >= pvalue.value() : pv = VE (     clip , pv.cov2() )
            nsigma = significance ( pv ) ## convert  it to significace
            
            pvalue = str ( ( 100 * pvalue ) .toString ( '%% 5.2f %s %%-.2f' % plus_minus ) )
            sigma  = str ( nsigma.toString ( '%%.2f %s %%-.2f' % plus_minus ) if float ( nsigma ) < 1000 else '+inf' ) 

            label = Labels.get( e , e ) 
            row   = label , infostr ( '*' ) , '' , '' , '' , '' , '' , pvalue , sigma 
            rows.append ( row ) 


        rows  = [ header ] + sorted ( rows )
        rows  = T.remove_empty_columns ( rows ) 
        return T.table ( rows , title = title , prefix = prefix , alignment = 'lcccccccc' , style = style )
        
    __repr__ = table
    __str__  = table

    # =========================================================================
    ## merge two objects:
    def merge ( self , other ) :
        self += other
        return self
    
    ## merge two objects:
    def __iadd__ ( self , other ) :
        """ Merge two GoF-toys objects """        
        if not isinstance ( other , GoFSimFitToys ) : return NotImplemented 

        ## (1) merge ECDFs 
        for key, content in loop_items ( other.ecdfs   ) :
            if not key in self.__ecdfs : self.__ecdfs [ key ] = content
            else :
                ecdfs = self.__ecdfs [ key ]
                for e , ecdf in content.items () :
                    if   e in ecdfs : ecdfs [ e ] += ecdf
                    else            : ecdfs [ e ]  = ecdf
                    
        ## (2) merge counters 
        for key, content in loop_items ( other.counters ) :
            if not key in self.__counters : self.__counters [ key ] = content
            else :
                counters = self.__counters [ key ]
                for e , cnt in content.items () : counters [ e ] += cnt

        ## (3) merge total
        for key, content in loop_items ( other.total ) :
            if not key in self.__total : self.__total [ key ]  = content
            else                       : self.__total [ key ] += content

        return self 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================


