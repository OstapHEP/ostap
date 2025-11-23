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
from   ostap.fitting.funbasic   import AFUN1
from   ostap.fitting.pdfbasic   import PDF1
from   ostap.core.core          import VE, Ostap
from   ostap.math.base          import axis_range, np2raw    
from   ostap.math.models        import f1_draw
from   ostap.utils.cidict       import cidict_fun
from   ostap.utils.basic        import loop_items, typename   
from   ostap.stats.counters     import SE, EffCounter 
from   ostap.logger.pretty      import pretty_float
from   ostap.math.ve            import fmt_pretty_ve
from   ostap.math.math_ve       import significance
from   ostap.logger.symbols     import plus_minus, times, greek_lower_sigma
from   ostap.logger.colorized   import infostr
from   ostap.stats.gof_utils    import Labels, Keys, clip_pvalue, data2vct  
from   collections              import defaultdict, namedtuple
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
## @var NL
#  use C++ if length of data exceeds NL, otherwise Python is OK 
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
    if NL < n and np2raw and isinstance ( cdf_data , numpy.ndarray ) :
        if cdf_data.dtype in ( numpy.float32 , numpy.float64 ) :
            raw_buffer , n = np2raw ( cdf_data )
            the_buffer = Ostap.Utils.make_buffer ( raw_buffer , n )
            return Ostap.GoF.kolmogorov_smirnov ( the_buffer )

    ## Long arrays to be converted to numpy 
    if NL < n : return kolmogorov_smirnov ( numpy.asarray ( cdf_data , dtype = numpy.float64 ) )

    ## for short arrays plain vanilla python is OK 
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
    if NL < n and np2raw and isinstance ( cdf_data , numpy.ndarray ) :
        if cdf_data.dtype in ( numpy.float32 , numpy.float64 ) :
            raw_buffer , n = np2raw ( cdf_data )
            the_buffer = Ostap.Utils.make_buffer ( raw_buffer , n )
            return Ostap.GoF.anderson_darling ( the_buffer )
        
    ## Long arrays to be converted to numpy 
    if NL < n : return anderson_darling ( numpy.asarray ( cdf_data , dtype = numpy.float64 ) )

    ## for short arrays plain vanilla python is OK 
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
    if NL < n and np2raw and isinstance ( cdf_data , numpy.ndarray ) :
        if cdf_data.dtype in ( numpy.float32 , numpy.float64 ) :
            raw_buffer , n = np2raw ( cdf_data )
            the_buffer = Ostap.Utils.make_buffer ( raw_buffer , n )
            return Ostap.GoF.cramer_von_mises  ( the_buffer )
        
    ## Long arrays to be converted to numpy 
    if NL < n : return cramer_von_mises ( numpy.asarray ( cdf_data , dtype = numpy.float64 ) )
    
    ## for short arrays plain vanilla python is OK 
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
    if NL < n and np2raw and isinstance ( cdf_data , numpy.ndarray ) :
        if cdf_data.dtype in ( numpy.float32 , numpy.float64 ) :
            raw_buffer , n = np2raw ( cdf_data )
            the_buffer = Ostap.Utils.make_buffer ( raw_buffer , n )
            return Ostap.GoF.kuiper  ( the_buffer )
        
    ## Long arrays to be converted to numpy 
    if NL < n : return kuiper ( numpy.asarray ( cdf_data , dtype = numpy.float64 ) )
    
    ## for short arrays plain vanilla python is OK 
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
    if NL < n and np2raw and isinstance ( cdf_data , numpy.ndarray ) :
        if cdf_data.dtype in ( numpy.float32 , numpy.float64 ) :
            raw_buffer , n = np2raw ( cdf_data )
            the_buffer = Ostap.Utils.make_buffer ( raw_buffer , n )
            return Ostap.GoF.ZA  ( the_buffer )
        
    ## Long arrays to be converted to numpy 
    if NL < n : return ZA ( numpy.asarray ( cdf_data , dtype = numpy.float64 ) )
    
    ## for short arrays plain vanilla python is OK 
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
    if NL < n and np2raw and isinstance ( cdf_data , numpy.ndarray ) :
        if cdf_data.dtype in ( numpy.float32 , numpy.float64 ) :
            raw_buffer , n = np2raw ( cdf_data )
            the_buffer = Ostap.Utils.make_buffer ( raw_buffer , n )
            return Ostap.GoF.ZK  ( the_buffer )
        
    ## Long arrays to be converted to numpy 
    if NL < n : return ZK ( numpy.asarray ( cdf_data , dtype = numpy.float64 ) ) 
    
    ## for short arrays plain vanilla python is OK 
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
    if NL < n and np2raw and isinstance ( cdf_data , numpy.ndarray ) :
        if cdf_data.dtype in ( numpy.float32 , numpy.float64 ) :
            raw_buffer , n = np2raw ( cdf_data )
            the_buffer = Ostap.Utils.make_buffer ( raw_buffer , n )
            return Ostap.GoF.ZC  ( the_buffer )
        
    ## Long arrays to be converted to numpy 
    if NL < n : return ZC ( numpy.asarray ( cdf_data , dtype = numpy.float64 ) ) 
    
    ## for short arrays plain vanilla python is OK 
    flog   = math.log
    result = sum ( ( flog ( ( 1.0 / Fi - 1 ) / ( ( n - 0.5 ) / ( i + 0.25 ) - 1 ) ) ) ** 2 for ( i , Fi ) in enumerate ( cdf_data ) )
    return result

# =========================================================================
## Clip input CDF arrays 
def vct_clip ( input , silent = True ) :
    """ Clip input CDF arrays"""
    vmin , vmax = 1.e-12 , 1 - 1.e-10
    if not silent and ( numpy.min ( input ) < vmin or vmax < numpy.max ( input ) ) :
        logger.warning ( 'Adjust CDF to be %s<cdf<%s' % ( vmin , vmax ) ) 
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

        self.__original_cdf = cdf if cdf and callable ( cdf ) else None 
        
        ## store PDF 
        self.__pdf        = pdf        
        self.__parameters = parameters if parameters else pdf.parameters ( dataset )
        
        ## re-load parameters
        self.__pdf.load_params ( self.__parameters , silent = True )

        ## get/construct/store  CDF 
        self.__cdf = self.get_cdf ( cdf ) 
        
        ## get the vectorized version of CDF 
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
        """ Serialize the object
        """
        self.pdf.load_params ( self.__parameters , silent = True )
        return { 'pdf'        : self.pdf            ,
                 'parameters' : self.__parameters   , 
                 'cdf'        : self.__original_cdf ,
                 'ecdf'       : self.ecdf           ,
                 'estimators' : self.estimators     } 
    
    ## De-serialize the object 
    def __setstate__ ( self , state ) :
        """ De-serialize the object
        """
        self.__pdf          = state.pop ( 'pdf'        ) 
        self.__parameters   = state.pop ( 'parameters' , {}  ) 
        self.__original_cdf = state.pop ( 'cdf'        )
        self.__ecdf         = state.pop ( 'ecdf'       )
        self.__estimators   = state.pop ( 'estimators' ) 

        ## (1) re-load parameters 
        self.pdf.load_params ( self.__parameters , silent = True )
        ## (2) re-reconstruct CDF 
        self.__cdf     = self.get_cdf ( self.__original_cdf )         
        ## vectorized form of CDF 
        self.__vct_cdf = numpy.vectorize ( self.cdf )
                
    # =========================================================================
    ## Get/Construct CDF 
    def get_cdf ( self , cdf ) :
        """ Get/Construct CDF
        """
        ## CDF is OK 
        if cdf and callable ( cdf ) : return cdf

        ## get CDF from PDF:
        pdf = self.__pdf

        ## consruct the CDF based on RooFit machinery 
        cdf = pdf.cdf()
        
        ## make a try to get the underlying c++ cdf-function from PDF 
        if hasattr ( pdf.pdf , 'setPars'  ) : pdf.pdf.setPars() 
        if hasattr ( pdf.pdf , 'function' ) :
            if hasattr ( pdf.pdf , 'setPars' ) : pdf.pdf.setPars()
            ## get the underlying C++ function 
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
        ## 
        return cdf
    
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
    def run ( self ,
              nToys    = 1000   , * ,
              parallel = False  ,
              silent   = False  ,
              nSplit   = 0      ) :
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
                ## dset = Ostap.MoreRooFit.delete_data ( dset )
                
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
        pv     = clip_pvalue ( pvalue , 0.5 ) 
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
        nsigma = str ( nsigma.toString ( '%%.2f %s %%-.2f' % plus_minus ) if float ( nsigma ) < 1000 else '+inf' ) 

        return ( what  ,
                 fmtv  % vs ,
                 fmt   % ( vm.value() , vm.error() ) ,
                 fmtv  % vr                          ,
                 fmt2  %  ( vmn , vmx )              ,
                 ( '%s10^%+d' %  ( times , expo )  if expo else '' ) , pvalue , nsigma )

    # =========================================================================
    ## Make a summary table
    def table ( self , title = '' , prefix = '' , width = 6 , precision = 4 , style = None ) :
        """ Make a summary table
        """
        import ostap.logger.table  as     T                 
        header = ( 'Statistics'        ,
                   'value'             ,
                   'mean'              , 
                   'rms'               ,
                   'min/max'           ,
                   'factor'            ,
                   'p-value [%]'       ,
                   '#%s' % greek_lower_sigma )
        rows   = [] 
        
        for label in self.ecdfs :
            
            result  = self.result ( label )
            if not result : continue

            the_label = Labels.get ( label , label )
            row = self.row ( the_label , result , width = width , precision = precision )
            rows.append ( row ) 
                    
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
        """ Draw ECDF for toys & statistical estimator 
        """
        key = cidict_fun ( what ) 
        if   key in Keys [ 'KS' ] and 'KS' in self.ecdfs :             
            result = self.result ( 'KS' )
            ecdf   = self.ecdfs  [ 'KS' ]
            ## logger.info ( 'Toy results for Kolmogorov-Smirnov estimate' ) 
        elif key in Keys [ 'K'  ] and 'K'  in self.ecdfs : 
            result = self.result ( 'K' )
            ecdf   = self.ecdfs  [ 'K' ]
            ## logger.info ( 'Toy results for Kuiper estimate' ) 
        elif key in Keys [ 'AD' ] and 'AD' in self.ecdfs :             
            result = self.result ( 'AD' )
            ecdf   = self.ecdfs  [ 'AD' ]
            ## logger.info ( 'Toy results for Anderson-Darling estimate' ) 
        elif key in Keys [ 'CM' ] and 'CM' in self.ecdfs : 
            result = self.result  ( 'CM' )
            ecdf   = self.ecdfs   [ 'CM' ]
            ## logger.info ( 'Toy results for Cramer-von Mises  estimate' ) 
        elif key in Keys [ 'ZK' ]  and 'ZK' in self.ecdfs : 
            result = self.result  ( 'ZK' )
            ecdf   = self.ecdfs   [ 'ZK' ]
            ## logger.info ( 'Toy results for Zhang/ZK estimate' ) 
        elif key in Keys [ 'ZA' ]   and 'ZA' in self.ecdfs :  
            result = self.result  ( 'ZA' )
            ecdf   = self.ecdfs   [ 'ZA' ]
            ## logger.info ( 'Toy results for Zhang/ZA estimate' ) 
        elif key in Keys [ 'ZC' ] and 'ZC' in self.ecdfs : 
            result = self.result  ( 'ZC' )
            ecdf   = self.ecdfs   [ 'ZC' ]
            ## logger.info ( 'Toy results for Zhang/ZC estimate' ) 
        else :
            raise KeyError (  "draw: Invalid `what`:%s" % what )
            
        xmin , xmax = ecdf.xmin () , ecdf.xmax ()
        value       = result.statistics
        xmin        = min ( xmin , value )
        xmax        = max ( xmax , value )
        xmin , xmax = axis_range ( xmin , xmax , delta = 0.20 )

        kwargs [ 'xmin' ] = kwargs.get ( 'xmin' , xmin ) 
        kwargs [ 'xmax' ] = kwargs.get ( 'xmax' , xmax )

        ## draw ECDF 
        result    = ecdf.draw  ( opts , *args , **kwargs )
        
        ## vertical line 
        line1    = ROOT.TLine ( value , 1e-3 , value , 1 - 1e-3 )
        ##

        ## horisontal line 
        xmin      = kwargs['xmin']
        xmax      = kwargs['xmax']
        dx        = ( xmax - xmin ) / 100 
        e         = ecdf ( value )
        line2     = ROOT.TLine ( xmin + dx , e , xmax - dx , e )
        ## 
        line2.SetLineWidth ( 2 ) 
        line2.SetLineColor ( 4 ) 
        line2.SetLineStyle ( 9 ) 
        ##
        line1.SetLineWidth  ( 4 ) 
        line1.SetLineColor  ( 8 )
        ##
        line2.draw ( 'same' )
        line1.draw ( 'same' )
        ##
        self._line1 = line1
        self._line2 = line2
        ##
        return result, line1, line2   

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================


