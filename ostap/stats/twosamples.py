#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/twosamples.py
#  Two-Sample Tests
#  @see https://www.jstor.org/stable/25471118
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2024-10-16
# =============================================================================
""" Two-Sample Tests
 - see https://www.jstor.org/stable/25471118
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2024-09-16"
__all__     = (
)
# =============================================================================
from   ostap.core.meta_info   import root_info 
from   ostap.core.ostap_types import listlike_types, sequence_types  
from   ostap.core.core        import SE, VE, Ostap
from   ostap.utils.basic      import loop_items, typename 
from   ostap.math.base        import doubles, axis_range, numpy   
from   ostap.math.models      import f1_draw
from   ostap.utils.cidict     import cidict_fun
from   ostap.logger.pretty    import pretty_float
from   ostap.stats.gof_utils  import Labels, Keys, clip_pvalue 
from   ostap.math.math_ve     import significance
from   ostap.math.ve          import fmt_pretty_ve
from   ostap.logger.symbols   import plus_minus, times, greek_lower_sigma
from   collections            import defaultdict, namedtuple 
import ostap.logger.table     as     T
import ostap.fitting.ds2numpy 
import ostap.fitting.roofit
import ROOT, math, array   
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.twosamples' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Two-sample test' )
# =============================================================================
if  ( 6 , 32 ) <= root_info : data2vct = lambda s : s
else                        : data2vct = lambda s : doubles ( s )
# =============================================================================
## transform data into empirical CDF
#  @code
#  data = ...
#  cdf  = ecdf_from_data ( data ) 
#  @endcode 
#  @param data input dataset 
#  @return empirical cumulative distribution function
#  @see Ostap::Math::ECDF 
def ecdf_from_data  ( data ) :
    """ Transform data into empirical CDF 
    - data : input dataset 
    >>> data = ...
    >>> cdf  = ecdf_from_data ( data )
    - see `Ostap.Math.ECDF`
    """
    if   isinstance  ( data , Ostap.Math.ECDF        ) : return data 
    elif numpy and isinstance ( data , numpy.ndarray ) : return Ostap.Math.ECDF ( data2vct ( data ) ) 
    elif isinstance  ( data  , array.array           ) : return Ostap.Math.ECDF ( data2vct ( data ) ) 
    elif isinstance  ( data1 , sequence_types        ) : return Ostap.Math.ECDF ( doubles  ( data ) ) 
    ## 
    raise TypeError ( "ecdf_from_data: Unsupported `data' type: %s" % typename ( data1 ) )
# ===============================================================================
## Prepare data
#  @code
#  data1 = ...
#  data2 = ...
#  ecdf1 , ecdf2 , ecdf = prepare_data1 ( data1 , data2 ) 
#  @endcode
#  @param data1 the 1st    dataset
#  @param data2 the second dataset
#  @param pooled (optional) the CDF for pooled data 
#  @return triplet of empiricla CDFs for the 1st, 2nd and pooled datasets 
def prepare_data1 ( data1 , data2 , pooled = None ) :
    """ Prepare data
    >>> data1 = ...
    >>> data2 = ...
    >>> ecdf1 , ecdf2 , ecdf = prepare_data1 ( data1 , data2 ) 
    - data1  : the 1st    dataset
    - data2  : the second dataset
    - pooled : (optiomal) the CDF for pooled data 
    - return a triplet of empirical CDFs for the 1st, 2nd and pooled datasets 
    """
    ecdf1 = ecdf_from_data ( data1 )
    ecdf2 = ecdf_from_data ( data2 )
    ##
    n1    = len ( ecdf1 )
    n2    = len ( ecdf2 )
    ##
    ecdf  = pooled if ( isinstance ( pooled , Ostap.Math.ECDF ) and n1 + n2 == len ( pooled ) ) else ecdf1 + ecdf2
    ## 
    return ecdf1 , ecdf2 , ecdf
# =============================================================================
## Prepare data
#  @code
#  data1 = ...
#  data2 = ...
#  ecdf1 , ecdf2 , ecdf , ranks1 , ranks2 = prepare_data2 ( data1 , data2 ) 
#  @endcode
#  @param data1 the 1st    dataset
#  @param data2 the second dataset
#  @param pooled (optional) the CDF for pooled data
#  @param ranks1 (optional) array of ranks of elements from the pooled sample in the 1st dataset
#  @param ranks2 (optional) array of ranks of elements from the pooled sample in the 2nd dataset
#  @return quintet of empirical CDFs for the 1st, 2nd and pooled datasets and two ransk arrays 
def prepare_data2 ( data1 , data2 , pooled = None , ranks1 = None , ranks2 = None ) :
    """ Prepare data
    >>> data1 = ...
    >>> data2 = ...
    >>> ecdf1 , ecdf2 , ecdf , ranks1 , ranks2 = prepare_data2 ( data1 , data2 ) 
    - data1  : the 1st    dataset
    - data2  : the second dataset
    - pooled : (optional) the CDF for pooled data
    - ranks1 : (optional) array of ranks of elements from the pooled sample in the 1st dataset
    - ranks2 : (optional) array of ranks of elements from the pooled sample in the 2nd dataset
    - return quintet of empirical CDFs for the 1st, 2nd and pooled datasets and two ransk arrays 
    """
    ##
    ecdf1 , ecdf2 , ecdf = prepare_data1 ( data1 , data2 , pooled )
    ##
    n1 = len ( ecdf1 )
    n2 = len ( ecdf2 )
    ## 
    ok1 = isinstance ( ranks1 , listlike_types ) and n1 + n2 == len ( ranks1 )
    if not ok1 :
        r1     = ecdf1.ranks ( ecdf )            
        ranks1 = numpy.asarray  ( r1 , dtype = numpy.dtype ( 'uint16' ) ) if numpy else array.array ( 'I' , r1 )
        ## 
    ok2 = isinstance ( ranks2 , listlike_types ) and n1 + n2 == len ( ranks2 )
    if not ok2 :
        r2     = ecdf2.ranks ( ecdf )            
        ranks2 = numpy.asarray  ( r2 , dtype = numpy.dtype ( 'uint16' ) ) if numpy else array.array ( 'I' , r2 )
    ## 
    return ecdf1, ecdf2, ecdf, ranks1 , ranks2

# =============================================================================
## Get Kolmogorov-Smirnov statistics KS
#  @code
#  ecdf1 =...
#  ecdf2 =...
#  ks = kolmogorov_smirnov ( ecdf1, ecdf2 )
#  @endcode
#  @param data1  the 1st dataset or empirical CDF
#  @param data2  the 2nd dataset or empirical CDF 
#  @param pooled (optional) pooled dataset or empirical CDF
#  @return Kolmogorov-Smirnov statistics KS
def kolmogorov_smirnov ( data1 , data2 , pooled = None ) :
    """ Get Kolmogorov-Smirnov statistics  KS
    - data1  : the 1st dataset or empirical CDF 
    - data2  : the 2nd dataset or empirical CDF 
    - pooled : (optional) pooled dataset or CDF
    """
    ecdf1, ecdf2, ecdf = prepare_data1 ( data1 , data2 , pooled ) 
    return max ( abs ( ecdf1 ( x ) - ecdf2 ( x ) ) for x in ecdf )
# ================================================================================
## Get Kuiper statistics K
#  @code
#  ecdf1 =...
#  ecdf2 =...
#  k  = kuiper ( ecdf1, ecdf2 )
#  @endcode
#  @param data1  the 1st dataset or empirical CDF
#  @param data2  the 2nd dataset or empirical CDF 
#  @param pooled (optional) pooled dataset or empirical CDF
#  @return Kuiper statistics K
def kuiper ( data1 , data2 , pooled = None  ) :
    """ Get kuiper statistics K
    - data1  : the 1st dataset or empirical CDF 
    - data2  : the 2nd dataset or empirical CDF 
    - pooled : (optional) pooled dataset or CDF
    """
    ecdf1, ecdf2, ecdf = prepare_data1 ( data1 , data2 , pooled ) 
    return max ( ecdf1 ( x ) - ecdf2 ( x ) for x in ecdf ) - min ( ecdf1 ( x ) - ecdf2 ( x ) for x in ecdf )
# =============================================================================
## Get Anderson-Darling statistics AD
#  @code
#  ecdf1 =...
#  ecdf2 =...
#  zk  = anderson_darling ( ecdf1, ecdf2 )
#  @endcode
#  @param data1  the 1st dataset or empirical CDF
#  @param data2  the 2nd dataset or empirical CDF 
#  @param pooled (optional) pooled dataset or empirical CDF
#  @param ranks1 (optional) ranks of the elements from the pooled data in the 1st dataset
#  @param ranks2 (optional) ranks of the elements from the pooled data in the 2nd dataset
#  @return Anderspon Darling statistics AD 
def anderson_darling  ( data1 , data2 , pooled = None , ranks1 = None  , ranks2 = None  ) :
    """ Get Anderson-Darling  statistics AD 
    - data1  : the 1st dataset or empirical CDF 
    - data2  : the 2nd dataset or empirical CDF 
    - pooled : (optional) pooled dataset or CDF
    - ranks1 : (optional) ranks of the elements from the pooled data in the 1st dataset
    - ranks2 : (optional) ranks of the elements from the pooled data in the 2nd dataset
    """
    ## transform input (if needed) 
    ecdf1, ecdf2, ecdf, ranks1, ranks2 = prepare_data2 ( data1 , data2 , pooled , ranks1 , ranks2 ) 
    
    n1   = len ( ecdf1 )
    n2   = len ( ecdf2 )
    n    = n1 + n2

    def term1 ( j ) :
        r1  = int ( ranks1 [ j ] ) 
        return ( n * r1 - j * n1 ) ** 2  / ( ( j + 1 ) *  ( n - j ) )
    def term2 ( j ) :
        r2  = int ( ranks2 [ j ] ) 
        return ( n * r2 - j * n2 ) ** 2  / ( ( j + 1 ) *  ( n - j ) )
        
    s1 = sum ( term1 ( j ) for j in range ( n ) ) / n1
    s2 = sum ( term2 ( j ) for j in range ( n ) ) / n2
    
    return ( s1 + s2 ) / n 

# =============================================================================
## Get Cramers-von Mises statistics CM
#  @code
#  ecdf1 =...
#  ecdf2 =...
#  zk  = cramer_von_mises ( ecdf1, ecdf2 )
#  @endcode
#  @param data1  the 1st dataset or empirical CDF
#  @param data2  the 2nd dataset or empirical CDF 
#  @param pooled (optional) pooled dataset or empirical CDF
#  @return Anderspon Darling statistics AD 
def cramer_von_mises ( data1 , data2 , pooled = None  ) :
    """ Get Crsmer-von Misesstatistics CM
    - data1  : the 1st dataset or empirical CDF 
    - data2  : the 2nd dataset or empirical CDF 
    - pooled : (optional) pooled dataset or CDF
    """
    ## transform input (if needed) 
    ecdf1, ecdf2, ecdf = prepare_data1 ( data1 , data2 , pooled ) 
    
    n1   = len ( ecdf1 )
    n2   = len ( ecdf2 )
    n    = n1 + n2
    
    ranks1 = ecdf.ranks ( ecdf1 )
    ranks2 = ecdf.ranks ( ecdf2 )
    
    u  = n1 * sum ( ( int ( ranks1 [ i ] ) - i ) ** 2 for i in range ( n1 ) ) 
    u += n2 * sum ( ( int ( ranks2 [ i ] ) - i ) ** 2 for i in range ( n2 ) ) 
    
    return u / ( n1 * n2 * n ) - ( 4 * n1 * n2 - 1 ) / ( 6.0 * n ) 

# =============================================================================
## Get ZK statistics ZK
#  @see https://www.jstor.org/stable/25471118
#  @code
#  ecdf1 =...
#  ecdf2 =...
#  zk  = ZK ( ecdf1, ecdf2 )
#  @endcode
#  @param data1  the 1st dataset or empirical CDF
#  @param data2  the 2nd dataset or empirical CDF 
#  @param pooled (optional) pooled dataset or empirical CDF
#  @param ranks1 (optional) ranks of the elements from the pooled data in the 1st dataset
#  @param ranks2 (optional) ranks of the elements from the pooled data in the 2nd dataset
#  @return ZK statistics ZK 
def ZK  ( data1 , data2 , pooled = None , ranks1 = None  , ranks2 = None  ) :
    """ Get ZK statistics 
    - data1  : the 1st dataset or empirical CDF 
    - data2  : the 2nd dataset or empirical CDF 
    - pooled : (optional) pooled dataset or CDF
    - ranks1 : (optional) ranks of the elements from the pooled data in the 1st dataset
    - ranks2 : (optional) ranks of the elements from the pooled data in the 2nd dataset
    - see https://www.jstor.org/stable/25471118
    """
    ## transform input (if needed) 
    ecdf1, ecdf2, ecdf, ranks1, ranks2 = prepare_data2 ( data1 , data2 , pooled , ranks1 , ranks2 ) 
    
    n1   = len ( ecdf1 )
    n2   = len ( ecdf2 )
    n    = n1 + n2

    flog = math.log 
    def term ( k ) :

        r1  = int ( ranks1 [ k ] ) 
        r2  = int ( ranks2 [ k ] ) 
        
        r1  = max ( 0.5 , min ( n1 - 0.5 , r1 ) ) ## r1 + 0.5 ??
        r2  = max ( 0.5 , min ( n2 - 0.5 , r2 ) ) ## r2 + 0.5 ?? 
        
        fk  = 0.5 / n if not k else float ( k ) / n
        
        f1k = float ( r1 ) / n1
        f2k = float ( r2 ) / n2
        
        rr  = n1 * ( f1k * flog ( f1k / fk ) + ( 1 - f1k ) * flog ( ( 1. - f1k ) / ( 1. - fk ) ) )
        rr += n2 * ( f2k * flog ( f2k / fk ) + ( 1 - f2k ) * flog ( ( 1. - f2k ) / ( 1. - fk ) ) ) 
        
        return rr
    
    return max ( term ( k ) for k in range ( n ) )

# =============================================================================
## Get ZA statistics ZA
#  @see https://www.jstor.org/stable/25471118
#  @code
#  ecdf1 =...
#  ecdf2 =...
#  za  = ZA ( ecdf1, ecdf2 )
#  @endcode
#  @param data1  the 1st dataset or empirical CDF
#  @param data2  the 2nd dataset or empirical CDF 
#  @param pooled (optional) pooled dataset or empirical CDF
#  @param ranks1 (optional) ranks of the elements from the pooled data in the 1st dataset
#  @param ranks2 (optional) ranks of the elements from the pooled data in the 2nd dataset
#  @return ZA statistics ZA 
def ZA ( data1 , data2 , pooled = None , ranks1 = None , ranks2 = None ) :
    """ Get ZA statistics 
    - data1  : the 1st dataset or empirical CDF 
    - data2  : the 2nd dataset or empirical CDF 
    - pooled : (optional) pooled dataset or CDF
    - ranks1 : (optional) ranks of the elements from the pooled data in the 1st dataset
    - ranks2 : (optional) ranks of the elements from the pooled data in the 2nd dataset
    - see https://www.jstor.org/stable/25471118
    """
    ## transform input (if needed) 
    ecdf1, ecdf2, ecdf, ranks1, ranks2 = prepare_data2 ( data1 , data2 , pooled , ranks1 , ranks2 ) 

    n1   = len ( ecdf1 )
    n2   = len ( ecdf2 )
    n    = n1 + n2
    
    flog = math.log 
    def term ( k ) :

        r1  = int ( ranks1 [ k ] ) 
        r2  = int ( ranks2 [ k ] ) 
        
        r1  = max ( 0.5 , min ( n1 - 0.5 , r1 ) ) ## r1 + 0.5 ??
        r2  = max ( 0.5 , min ( n2 - 0.5 , r2 ) ) ## r2 + 0.5 ?? 
        
        fk  = 0.5 / n if not k else float ( k ) / n  
        f1k = float ( r1 ) / n1
        f2k = float ( r2 ) / n2

        rr  = n1 * ( f1k * flog ( f1k ) + ( 1 - f1k ) * flog ( 1. - f1k ) ) / ( ( k + 0.5 ) * ( n - k - 0.5 ) )
        rr += n2 * ( f2k * flog ( f2k ) + ( 1 - f2k ) * flog ( 1. - f2k ) ) / ( ( k + 0.5 ) * ( n - k - 0.5 ) )
        
        return rr
    
    return sum ( term ( k ) for k in range ( n ) )

# =============================================================================
## Get ZC statistics ZC
#  @see https://www.jstor.org/stable/25471118
#  @code
#  ecdf1 =...
#  ecdf2 =...
#  zc  = ZC ( ecdf1, ecdf2 )
#  @endcode
#  @param data1  the 1st dataset or empirical CDF
#  @param data2  the 2nd dataset or empirical CDF 
#  @param pooled (optional) pooled dataset or empirical CDF
#  @param ranks1 (optional) ranks of the elements from the pooled data in the 1st dataset
#  @param ranks2 (optional) ranks of the elements from the pooled data in the 2nd dataset
#  @return ZA statistics ZA 
def ZC ( data1 , data2 , pooled = None , ranks1 = None , ranks2 = None ) :
    """ Get ZC statistics 
    - data1  : the 1st dataset or empirical CDF 
    - data2  : the 2nd dataset or empirical CDF 
    - pooled : (optional) pooled dataset or CDF
    - ranks1 : (optional) ranks of the elements from the pooled data in the 1st dataset
    - ranks2 : (optional) ranks of the elements from the pooled data in the 2nd dataset
    - see https://www.jstor.org/stable/25471118
    """
    ## transform input (if needed) 
    ecdf1, ecdf2, ecdf, ranks1, ranks2 = prepare_data2 ( data1 , data2 , pooled , ranks1 , ranks2 ) 

    n1  = len ( ecdf1 )
    n2  = len ( ecdf2 )
    n   = n1 + n2
    
    flog = math.log 
    def term1 ( j ) :
        r1  = int ( ranks1 [ j ] ) 
        return flog ( n1 / ( j + 0.5 ) - 1 ) * flog ( n / ( r1 + 0.5 ) - 1 )
    
    def term2 ( j ) :
        r2  = int ( ranks2 [ j ] ) 
        return flog ( n2 / ( j + 0.5 ) - 1 ) * flog ( n / ( r2 + 0.5 ) - 1 )
    
    t1 = sum ( term1 ( j ) for j in range ( n1 ) )
    t2 = sum ( term2 ( j ) for j in range ( n2 ) )

    return ( t1 + t2 ) / n 

# =============================================================================
## @class TSTest
#  Two-Sample Test
class TSTest(object):
    """ Two-sample test 
    """
    def __init__ ( self  ,
                   data1 ,
                   data2 ) :
        
        ## transform input data 
        ecdf1, ecdf2, ecdf, ranks1, ranks2 = prepare_data2 ( data1 , data2 ) 
        
        self.__ecdf1  = ecdf1 
        self.__ecdf2  = ecdf2
        self.__ecdf   = ecdf 
        self.__ranks1 = ranks1 
        self.__ranks2 = ranks2 
        
        self.__estimators = {            
            'KS' : kolmogorov_smirnov ( self.ecdf1 , self.ecdf2 , self.ecdf ) , 
            'K'  : kuiper             ( self.ecdf1 , self.ecdf2 , self.ecdf ) , 
            'AD' : anderson_darling   ( self.ecdf1 , self.ecdf2 , self.ecdf , self.ranks1 , self.ranks2 ) , 
            'CM' : cramer_von_mises   ( self.ecdf1 , self.ecdf2 , self.ecdf ) , 
            'ZK' : ZK                 ( self.ecdf1 , self.ecdf2 , self.ecdf , self.ranks1 , self.ranks2 ) , 
            'ZA' : ZA                 ( self.ecdf1 , self.ecdf2 , self.ecdf , self.ranks1 , self.ranks2 ) , 
            'ZC' : ZC                 ( self.ecdf1 , self.ecdf2 , self.ecdf , self.ranks1 , self.ranks2 ) , 
        }

    @property
    def ecdf1 ( self ) :
        """`ecdf1` : empirical CDF for the 1st dataset
        """
        return self.__ecdf1 
    @property
    def ecdf2 ( self ) :
        """`ecdf2` : empirical CDF for the 2nd dataset
        """
        return self.__ecdf2 
    @property
    def ecdf  ( self ) :
        """`ecdf` : empirical CDF for combined dataset
        """
        return self.__ecdf
    @property
    def ranks1 ( self ) :
        """ ranks of elements from pooled sample in the 1st dataset"""
        return self.__ranks1
    @property
    def ranks2 ( self ) :
        """ ranks of elements from pooled sample in the 2nd dataset"""
        return self.__ranks2

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
        title = title if title else 'Two samples Test' 
        rows  = T.remove_empty_columns ( rows )
        return T.table ( rows , title = title , prefix = prefix , alignment = 'lcl' , style = style  )

    ## print estimator as table 
    __repr__ = table
    __str__  = table
    
    # =========================================================================
    @property
    def estimators ( self ) :
        """`estimators` : get a dictionary of all estimators"""
        return self.__estimators

    # =========================================================================
    ## Get Kolmogorov-Smirnov statistiscs 
    @property 
    def kolmogorov_smirnov_estimator ( self ) :
        """ Get Kolmogorov-Smirnov statistiscs KS
        """
        return self.estimators.get ( 'KS' , None ) 

    # =========================================================================
    ## Get Kuiper  statististics 
    @property 
    def kuiper_estimator ( self ) :
        """ Get Kuiperstatistics
        """        
        return self.estimators.get ( 'K' , None ) 

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
        return self.estimators.get ( 'ZA' , None ) 
    
    # =========================================================================
    ## Get ZC statististics 
    @property 
    def ZC_estimator ( self ) :
        """ Get ZC statistics
        """        
        return self.estimators.get ( 'ZC' , None ) 
    
    # =========================================================================
    def ks_acccept ( self , alpha ) :

        assert 0 < alpha < 1 , 'Invalid "alpha"%s' % alpha
        c_alpha = math.sqrt ( -1 * math.log ( alpha / 2. ) / 2. ) 
        n   = len ( self.ecdf1 )
        m   = len ( self.ecdf2 )
        d   = self.__KS_val
        return d > c_alpha * math.sqrt ( ( n + m ) / ( n * m ) ) 

    # =========================================================================
    def ks_alpha ( self ) :
        
        n = len ( self.ecdf1 )
        m = len ( self.ecdf2 )
        d = self.__KS_val
        c = d / math.sqrt ( ( n + m ) / ( n * m ) ) 
        return 2 * math.exp ( -2 * c * c )

    # =========================================================================
    ## draw all involved empirical CDFS
    def draw  ( self , opts = '' , *args , **kwargs ) :
        """ Draw fit CDF & empirical CDF
        """

        ecdf1 = self.ecdf1
        ecdf2 = self.ecdf2
        ecdf  = self.ecdf

        xmin  , xmax  = ecdf .xmin (), ecdf .xmax()
        xmin1 , xmax1 = ecdf1.xmin (), ecdf1.xmax()
        xmin2 , xmax2 = ecdf2.xmin (), ecdf2.xmax()

        xmin = min ( xmin , xmin1 , xmin2 )
        xmax = max ( xmax , xmax1 , xmax2 )

        xmin , xmax = axis_range ( xmin , xmax , delta = 0.15 )

        xmin = kwargs.pop ( 'xmin' , xmin )
        xmax = kwargs.pop ( 'xmax' , xmax )
        
        opts    = opts.strip() 
        optsame = 'same' if not opts else '%s %s' % ( 'same' , opts ) 
        r  = ecdf .draw ( opts    , color = 8 , linewidth = 3 , xmin = xmin , xmax = xmax , **kwargs )
        r1 = ecdf1.draw ( optsame , color = 2 , linewidth = 3 , xmin = xmin , xmax = xmax , **kwargs )
        r2 = ecdf2.draw ( optsame , color = 4 , linewidth = 3 , xmin = xmin , xmax = xmax , **kwargs )
        
        return r, r1, r2

# =============================================================================
## @class TSToys
#  Two-Sample Test
class TSToys(TSTest):
    """ Two-sample test 
    """
    ## result of 2-samples test
    Result = namedtuple ( 'Result' , 'statistics counter pvalue nsigma' )
    # =========================================================================
    ## Initialize TSToys object 
    def __init__ ( self           ,
                   data1          ,
                   data2          ,
                   nToys  = 100   , 
                   silent = False ) :
        
        ## initialize the base
        TSTest.__init__ ( self , data1 , data2 )

        self.__counters = defaultdict(SE) 
        self.__ecdfs    = {} 
        self.__silent   = True if silent else False
        
        self.__nToys  = 0 
        if 0 < nToys : self.run ( nToys , silent = silent ) 
    
    # =========================================================================
    ## number of toys/permutations 
    @property
    def nToys ( self ) :
        """`nToys` : number of toys/permutations"""
        return self.__nToys
    
    # =========================================================================
    ## ECDFs
    @property 
    def ecdfs ( self ) :
        """`ecdfs` : toys/permutations results as empirical cumulative distribution functions"""
        return self.__ecdfs
    
    # =========================================================================
    ## Counters  
    @property 
    def counters ( self ) :
        """`counters` : toy results as counters"""
        return self.__counters
    
    # ===============================================================================
    ## run toys 
    def run ( self , nToys = 1000 , silent = False ) :
        """ Run toys 
        """ 
        assert isinstance ( nToys , int ) and 0 < nToys , "Invalid `nToys` argument!"

        n1   = len ( self.ecdf1 )
        n2   = len ( self.ecdf1 )
        
        data = numpy.asarray ( self.ecdf.data () , dtype = float )
        
        results  = defaultdict(list)
        counters = self.counters
        
        from ostap.utils.progress_bar import progress_bar 
        for i in progress_bar ( nToys , silent = silent , description = 'Permutations:') :

            numpy.random.shuffle ( data )
            data1 = data [    : n1 ]
            data2 = data [ n1 :    ] 
            
            ## transform input data 
            ecdf1, ecdf2, ecdf, ranks1, ranks2 = prepare_data2 ( data1 , data2 ) 
            
            ks = kolmogorov_smirnov ( ecdf1 , ecdf2 , ecdf )
            k  = kuiper             ( ecdf1 , ecdf2 , ecdf )
            ad = anderson_darling   ( ecdf1 , ecdf2 , ecdf , ranks1 , ranks2 )            
            cm = cramer_von_mises   ( ecdf1 , ecdf2 , ecdf )             
            zk = ZK                 ( ecdf1 , ecdf2 , ecdf , ranks1 , ranks2 )
            za = ZA                 ( ecdf1 , ecdf2 , ecdf , ranks1 , ranks2 )
            zc = ZC                 ( ecdf1 , ecdf2 , ecdf , ranks1 , ranks2 )
            
            counters [ 'KS'  ] += ks 
            counters [ 'K'   ] += k 
            counters [ 'AD'  ] += ad
            counters [ 'CM'  ] += cm
            counters [ 'ZK'  ] += zk
            counters [ 'ZA'  ] += za
            counters [ 'ZC'  ] += zc
            
            results  [ 'KS'  ] .append ( ks )    
            results  [ 'K'   ] .append ( k  ) 
            results  [ 'AD'  ] .append ( ad ) 
            results  [ 'CM'  ] .append ( cm ) 
            results  [ 'ZK'  ] .append ( zk ) 
            results  [ 'ZA'  ] .append ( za ) 
            results  [ 'ZC'  ] .append ( zc ) 
            
        ## accumulate number of toys 
        self.__nToys += nToys 

        ECDF = Ostap.Math.ECDF 
        for key in results :
            data = results [ key ]
            if not data : continue
            if not key in self.__ecdfs : self.__ecdfs [ key ] = ECDF ( data , True ) ## complementary ECDF!
            else                       : self.__ecdfs [ key ]  .add  ( data2vct ( data ) ) 
            
        del results

    # =========================================================================
    ## Get Kolmogorov-Smirnov statistiscs 
    @property 
    def kolmogorov_smirnov ( self ) :
        """ Get Kolmogorov-Smirnov statistiscs KS
        """
        return self.result ( 'KS' ) 
    
    # =========================================================================
    ## Get Kuiper statististics 
    @property 
    def kuiper ( self ) :
        """ Get Kuiper statistics 
        """        
        return self.result ( 'K' ) 
           
    # =========================================================================
    ## Get ZK statististics 
    @property 
    def ZK  ( self ) :
        """ Get ZK statistics 
        """
        return self.result ( 'ZK' ) 
        
    # =========================================================================
    ## Get ZA statististics 
    @property 
    def ZA  ( self ) :
        """ Get ZA statistics
        """        
        return self.result ( 'ZA' ) 
    
    # =========================================================================
    ## Get ZC statististics 
    @property 
    def ZC  ( self ) :
        """ Get ZC statistics 
        """        
        return self.result ( 'ZC' ) 

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
        pv = clip_pvalue ( pvalue , 0.5 ) 
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
        header = ( 'Statistics'   ,
                   'value'       ,
                   'mean'        ,
                   'rms'         ,
                   'min/max'     ,
                   'factor'      ,
                   'p-value [%]' ,
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
        """ Draw ECDF for toys & statistical estgimator 
        """
        key = cidict_fun ( what ) 
        if   key in Keys [ 'KS' ] and 'KS' in self.ecdfs :             
            result = self.result ( 'KS' )
            ecdf   = self.ecdfs  [ 'KS' ]
            logger.info ( 'Toy results for Kolmogorov-Smirnov estimate' ) 
        elif key in Keys [ 'K'  ] and 'K'  in self.ecdfs : 
            result = self.result ( 'K' )
            ecdf   = self.ecdfs  [ 'K' ]
            logger.info ( 'Toy results for Kuiper estimate' ) 
        elif key in Keys [ 'AD' ] and 'AD' in self.ecdfs :             
            result = self.result ( 'AD' )
            ecdf   = self.ecdfs  [ 'AD' ]
            logger.info ( 'Toy results for Anderson-Darling estimate' ) 
        elif key in Keys [ 'CM' ] and 'CM' in self.ecdfs : 
            result = self.result  ( 'CM' )
            ecdf   = self.ecdfs   [ 'CM' ]
            logger.info ( 'Toy results for Cramer-von Mises  estimate' ) 
        elif key in Keys [ 'ZK' ]  and 'ZK' in self.ecdfs : 
            result = self.result  ( 'ZK' )
            ecdf   = self.ecdfs   [ 'ZK' ]
            logger.info ( 'Toy results for Zhang/ZK estimate' ) 
        elif key in Keys [ 'ZA' ]   and 'ZA' in self.ecdfs :  
            result = self.result  ( 'ZA' )
            ecdf   = self.ecdfs   [ 'ZA' ]
            logger.info ( 'Toy results for Zhang/ZA estimate' ) 
        elif key in Keys [ 'ZC' ] and 'ZC' in self.ecdfs : 
            result = self.result  ( 'ZC' )
            ecdf   = self.ecdfs   [ 'ZC' ]
            logger.info ( 'Toy results for Zhang/ZC estimate' ) 
        else :
            raise KeyError (  "draw: Invalid `what`:%s" % what )

        xmin , xmax = ecdf.xmin () , ecdf.xmax ()
        value       = result.statistics
        xmin        = min ( xmin , value )
        xmax        = max ( xmax , value )
        xmin , xmax = axis_range ( xmin , xmax , delta = 0.20 )

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
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================




