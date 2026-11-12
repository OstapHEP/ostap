#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/twosamples.py
#  Two-Sample Tests for 1D non-weighted datasets 
#  @see https://www.jstor.org/stable/25471118
#  For dataset with weights, see 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2024-10-16
# =============================================================================
""" Two-Sample Tests for 1D non-weighted datasets 
 - see https://www.jstor.org/stable/25471118
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2024-09-16"
__all__     = (
    'kolmogorov_smirnov' , ## Kolmogorov-smirnov statistics
    'kuiper'             , ## Kuiper             statistics
    'anderson_darling'   , ## Abderson-Darling   statistics
    'cramer_von_mises'   , ## Cramer-von Mises   statistics
    'berk_jones'         , ## Berk-Jones         statistics
    'ZA'                 , ## Zhang ZA           statistics
    'ZC'                 , ## Zhang ZC           statistics
    'ZK'                 , ## Zhang ZK           statistics
    'wasserstein1'       , ## Wasserstein/1      statistics
    'wasserstein2'       , ## Wasserstein/2      statistics
    'rosenblatt'         , ## Rosenblatt         statistics
    'watson_u2'          , ## Watson U2          statistics
    'mann_whitney'       , ## Mann-Whitney       statistics
    'savage'             , ## symmetric Savage   statistics
    ##
    "TwoSamplesTest"     , ## Two-Samples Test
    "TwoSamplesToys"     , ## toys for Two-samples Test 
    ## 
)
# =============================================================================
from   ostap.core.ostap_types import listlike_types, sequence_types  
from   ostap.core.core        import VE, Ostap
from   ostap.utils.core       import typename 
from   ostap.utils.basic      import safe_loop 
from   ostap.math.math_base   import doubles, axis_range, data2vct 
from   ostap.math.models      import f1_draw
from   ostap.utils.cidict     import cidict_fun
from   ostap.logger.pretty    import pretty_float, format_pow10 
from   ostap.stats.gof_utils  import Labels, Keys, clip_pvalue 
from   ostap.stats.utils      import full_header, t_header  
from   ostap.math.math_ve     import significance
from   ostap.math.ve          import fmt_pretty_ve
from   ostap.stats.counters   import SE, ECDF
from   ostap.stats.pvalue     import permutation_symbol
from   ostap.logger.symbols   import plus_minus, times, greek_lower_sigma, infinity_pos 
from   collections            import defaultdict, namedtuple 
import ostap.logger.table     as     T
import ostap.fitting.ds2numpy 
import ostap.fitting.roofit
import ROOT, math, array, numpy    
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
logger = getLogger ( 'ostap.stats.twosamples' )
# =============================================================================
logger.debug ( 'Two-Samples Test for 1D non-weighted datasets' )
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
    if   isinstance ( data , ECDF             ) : return data 
    elif isinstance ( data , numpy.ndarray    ) : return ECDF ( data2vct ( data ) ) 
    elif isinstance ( data , array.array      ) : return ECDF ( data2vct ( data ) ) 
    elif isinstance ( data , sequence_types   ) : return ECDF ( doubles  ( data ) ) 
    ## 
    raise TypeError ( "ecdf_from_data: Unsupported `data' type: %s" % typename ( data1 ) )

# ===============================================================================
## Prepare data
#  @code
#  data1 = ...
#  data2 = ...
#  ecdf1 , ecdf2 = prepare_data1 ( data1 , data2 ) 
#  @endcode
#  @param data1 the 1st    dataset
#  @param data2 the second dataset
#  @return pair of empirical CDFs for the 1st and 2nd second daatsets 
def prepare_data ( data1 , data2 ) : 
    """ Prepare data
    >>> data1 = ...
    >>> data2 = ...
    >>> ecdf1 , ecdf2 , ecdf = prepare_data1 ( data1 , data2 ) 
    - data1  : the 1st    dataset
    - data2  : the second dataset
    - return a pair of of empirical CDFs for the 1st and 2nd datasets 
    """
    ecdf1 = ecdf_from_data ( data1 )
    ecdf2 = ecdf_from_data ( data2 )
    ##
    return ecdf1 , ecdf2

# =============================================================================
## Get Kolmogorov-Smirnov statistics KS
#  @code
#  ecdf1 =...
#  ecdf2 =...
#  ks = kolmogorov_smirnov ( ecdf1, ecdf2 )
#  @endcode
#  @param data1  the 1st dataset or empirical CDF
#  @param data2  the 2nd dataset or empirical CDF 
#  @return Kolmogorov-Smirnov statistics KS
def kolmogorov_smirnov ( data1 , data2 ) :
    """ Get Kolmogorov-Smirnov statistics  KS
    - data1  : the 1st dataset or empirical CDF 
    - data2  : the 2nd dataset or empirical CDF 
    """
    ecdf1, ecdf2 = prepare_data ( data1 , data2  ) 
    return Ostap.Math.kolmogorov_smirnov ( ecdf1 , ecdf2 )

# ================================================================================
## Get Kuiper statistics K
#  @code
#  ecdf1 =...
#  ecdf2 =...
#  k  = kuiper ( ecdf1, ecdf2 )
#  @endcode
#  @param data1  the 1st dataset or empirical CDF
#  @param data2  the 2nd dataset or empirical CDF 
#  @return Kuiper statistics K
def kuiper ( data1 , data2 ) :
    """ Get kuiper statistics K
    - data1  : the 1st dataset or empirical CDF 
    - data2  : the 2nd dataset or empirical CDF 
    """
    ecdf1, ecdf2 = prepare_data ( data1 , data2 )
    return Ostap.Math.kuiper( ecdf1 , ecdf2 )

# =============================================================================
## Get Anderson-Darling statistics AD
#  @code
#  ecdf1 =...
#  ecdf2 =...
#  ad  = anderson_darling ( ecdf1, ecdf2 )
#  @endcode
#  @param data1  the 1st dataset or empirical CDF
#  @param data2  the 2nd dataset or empirical CDF 
#  @return Anderson-Darling statistics AD 
def anderson_darling  ( data1 , data2  ) : 
    """ Get Anderson-Darling  statistics AD 
    - data1  : the 1st dataset or empirical CDF 
    - data2  : the 2nd dataset or empirical CDF 
    """
    ecdf1, ecdf2 = prepare_data ( data1 , data2 ) 
    return Ostap.Math.anderson_darling ( ecdf1 , ecdf2 )

# =============================================================================
## Get Cramers-von Mises statistics CM
#  @code
#  ecdf1 =...
#  ecdf2 =...
#  cm  = cramer_von_mises ( ecdf1, ecdf2 )
#  @endcode
#  @param data1  the 1st dataset or empirical CDF
#  @param data2  the 2nd dataset or empirical CDF 
#  @return Cramer-von Mises statistics CM 
def cramer_von_mises ( data1 , data2 ) :
    """ Get Cramer-von Mises statistics CM
    - data1  : the 1st dataset or empirical CDF 
    - data2  : the 2nd dataset or empirical CDF 
    """
    ## transform input (if needed) 
    ecdf1, ecdf2 = prepare_data ( data1 , data2 ) 
    return Ostap.Math.cramer_von_mises ( ecdf1 , ecdf2 )

# =============================================================================
## Get Berk-Jones statistics BJ
#  @code
#  ecdf1 =...
#  ecdf2 =...
#  bj  = cramer_von_mises ( ecdf1, ecdf2 )
#  @endcode
#  @param data1  the 1st dataset or empirical CDF
#  @param data2  the 2nd dataset or empirical CDF 
#  @return Berk-Jones statistics BJ 
def berk_jones ( data1 , data2 ) :
    """ Get Berk-Jones statistics BJ
    - data1  : the 1st dataset or empirical CDF 
    - data2  : the 2nd dataset or empirical CDF 
    """
    ## transform input (if needed) 
    ecdf1, ecdf2 = prepare_data ( data1 , data2 ) 
    return Ostap.Math.berk_jones ( ecdf1 , ecdf2 )

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
#  @return ZK statistics ZK 
def ZK  ( data1 , data2 ) :
    """ Get ZK statistics 
    - data1  : the 1st dataset or empirical CDF 
    - data2  : the 2nd dataset or empirical CDF 
    - see https://www.jstor.org/stable/25471118
    """
    ecdf1, ecdf2 = prepare_data ( data1 , data2 ) 
    return Ostap.Math.ZK ( ecdf1 , ecdf2 )

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
#  @return ZA statistics ZA 
def ZA ( data1 , data2 ) : 
    """ Get ZA statistics 
    - data1  : the 1st dataset or empirical CDF 
    - data2  : the 2nd dataset or empirical CDF 
    - see https://www.jstor.org/stable/25471118
    """
    ecdf1, ecdf2 = prepare_data ( data1 , data2 ) 
    return Ostap.Math.ZA ( ecdf1 , ecdf2 )

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
#  @return ZA statistics ZA 
def ZC ( data1 , data2 ) :
    """ Get ZC statistics 
    - data1  : the 1st dataset or empirical CDF 
    - data2  : the 2nd dataset or empirical CDF 
    - see https://www.jstor.org/stable/25471118
    """
    ecdf1, ecdf2 = prepare_data ( data1 , data2  ) 
    return Ostap.Math.ZC ( ecdf1 , ecdf2 )

# =============================================================================
## Wasserstein statistics W1
#  @code
#  ecdf1 =...
#  ecdf2 =...
#  w1    = wasserstein1( ecdf1, ecdf2 )
#  @endcode
#  @param data1  the 1st dataset 
#  @param data2  the 2nd dataset 
#  @return Wasserstein statistics W1
def wasserstein1 ( data1 , data2 ) : 
    """ Get Wasserstein statistics W1
    - data1  : the 1st dataset 
    - data2  : the 2nd dataset 
    """
    ## transform input (if needed) 
    ecdf1, ecdf2 = prepare_data ( data1 , data2 ) 
    return Ostap.Math.wasserstein1 ( ecdf1 , ecdf2 )

# =============================================================================
## Wasserstein statistics W2
#  @code
#  ecdf1 =...
#  ecdf2 =...
#  w2    = wasserstein1( ecdf1, ecdf2 )
#  @endcode
#  @param data1  the 1st dataset 
#  @param data2  the 2nd dataset 
#  @return Wasserstein statistics W2
def wasserstein2 ( data1 , data2 ) : 
    """ Get Wasserstein statistics W2
    - data1  : the 1st dataset 
    - data2  : the 2nd dataset 
    """
    ## transform input (if needed) 
    ecdf1, ecdf2 = prepare_data ( data1 , data2 ) 
    return Ostap.Math.wasserstein2 ( ecdf1 , ecdf2 )

# =============================================================================
## Rosenblatt statistics R 
#  @code
#  ecdf1 =...
#  ecdf2 =...
#  r     = rosenblatt ( ecdf1, ecdf2 )
#  @endcode
#  @param data1  the 1st dataset 
#  @param data2  the 2nd dataset 
#  @return Rosenblatt statistics W1
def rosenblatt ( data1 , data2 ) : 
    """ Get Rosenblatt statistics R
    - data1  : the 1st dataset 
    - data2  : the 2nd dataset 
    """
    ## transform input (if needed) 
    ecdf1, ecdf2 = prepare_data ( data1 , data2 ) 
    return Ostap.Math.rosenblatt  ( ecdf1 , ecdf2 )

# =============================================================================
## Watson statistics U2
#  @code
#  ecdf1 =...
#  ecdf2 =...
#  u2    = watson_u2 ( ecdf1, ecdf2 )
#  @endcode
#  @param data1  the 1st dataset 
#  @param data2  the 2nd dataset 
#  @return Watson statistics U2 
def watson_u2 ( data1 , data2 ) :
    """ Get Watson statistics U2
    - data1  : the 1st dataset 
    - data2  : the 2nd dataset 
    """
    ## transform input (if needed) 
    ecdf1, ecdf2 = prepare_data ( data1 , data2 ) 
    return Ostap.Math.watson_u2 ( ecdf1 , ecdf2 )

# =============================================================================
## Matt-Whitney statistics 
#  @code
#  ecdf1 =...
#  ecdf2 =...
#  mw    = mann_whitney ( ecdf1, ecdf2 )
#  @endcode
#  @param data1  the 1st dataset 
#  @param data2  the 2nd dataset 
#  @return Matt-Whitney statistics 
def mann_whitney ( data1 , data2 ) :
    """ Get Mann-Whitney statistics 
    - data1  : the 1st dataset 
    - data2  : the 2nd dataset 
    """
    ## transform input (if needed) 
    ecdf1, ecdf2 = prepare_data ( data1 , data2 ) 
    return Ostap.Math.mann_whitney ( ecdf1 , ecdf2 )

# =============================================================================
## Symmetrized Savage statistics 
#  @code
#  ecdf1 =...
#  ecdf2 =...
#  ss    = savage ( ecdf1, ecdf2 )
#  @endcode
#  @param data1  the 1st dataset 
#  @param data2  the 2nd dataset 
#  @return Savage statistics 
def savage ( data1 , data2   )
    """ Get Savage statistics 
    - data1  : the 1st dataset 
    - data2  : the 2nd dataset 
    """
    ## transform input (if needed) 
    ecdf1, ecdf2 = prepare_data ( data1 , data2 ) 
    return Ostap.Math.savage ( ecdf1 , ecdf2 )

# =============================================================================
## @class ThoSamplesTest 
#  Two-Samples Test
class TwoSamplesTest ( object ) :
    """ Two-sample test 
    """
    def __init__ ( self  ,
                   data1 ,
                   data2 ) :
        
        ## transform input data 
        ecdf1, ecdf2 = prepare_data ( data1 , data2 )
        
        self.__ecdf1 = ecdf1
        self.__ecdf2 = ecdf2
        self.__ecdf  = ecdf1.merge ( ecdf2 ) 
        
        self.__estimators = {            
            'KS' : kolmogorov_smirnov ( self.ecdf1 , self.ecdf2 ) , 
            'K'  : kuiper             ( self.ecdf1 , self.ecdf2 ) , 
            'AD' : anderson_darling   ( self.ecdf1 , self.ecdf2 ) , 
            'CM' : cramer_von_mises   ( self.ecdf1 , self.ecdf2 ) , 
            'BJ' : berk_jones         ( self.ecdf1 , self.ecdf2 ) , 
            'ZK' : ZK                 ( self.ecdf1 , self.ecdf2 ) , 
            'ZA' : ZA                 ( self.ecdf1 , self.ecdf2 ) , 
            'ZC' : ZC                 ( self.ecdf1 , self.ecdf2 ) , 
            'W1' : wasserstein1       ( self.ecdf1 , self.ecdf2 ) , 
            'W2' : wasserstein2       ( self.ecdf1 , self.ecdf2 ) , 
            'R'  : rosenblatt         ( self.ecdf1 , self.ecdf2 ) , 
            'U2' : watson_u2          ( self.ecdf1 , self.ecdf2 ) , 
            'MW' : mann_whitney       ( self.ecdf1 , self.ecdf2 ) , 
            'S'  : savage             ( self.ecdf1 , self.ecdf2 )
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
    def ecdf ( self ) :
        """`ecdf` : total/combined CDF for both datasets
        """
        return self.__ecdf
    
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
        header = ( 'Statistics' , )  + t_header 
        rows   = [] 
        
        for label , value  in self.estimators.items () :
            
            the_label = Labels.get ( label , label )
            
            result , expo = pretty_float ( value , width = width , precision = precision )            
            row = the_label , result , ( '[%s]' % format_pow10 ( expo ) )  if expo else ''
            
            rows.append ( row )
            
        rows  = [ header ] + sorted ( rows ) 
        rows  = T.remove_empty_columns ( rows )
        title = title if title else 'Two-Samples Test'    
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
    def draw  ( self , option = '' , *options , **kwargs ) :
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

        xmin , xmax = axis_range ( xmin , xmax , delta = 0.10 )

        xmin = kwargs.pop ( 'xmin' , xmin )
        xmax = kwargs.pop ( 'xmax' , xmax )
        
        option  = option.strip() 
        optsame = 'same' if not option else '%s %s' % ( 'same' , option ) 
        r  = ecdf .draw ( option  , *options , color = 8 , linewidth = 3 , xmin = xmin , xmax = xmax , **kwargs )
        r1 = ecdf1.draw ( optsame , *options , color = 2 , linewidth = 3 , xmin = xmin , xmax = xmax , **kwargs )
        r2 = ecdf2.draw ( optsame , *options , color = 4 , linewidth = 3 , xmin = xmin , xmax = xmax , **kwargs )
        
        return r, r1, r2

# =============================================================================
## @class TwoSamplesToys
#  Toys for Two-Samples Test
class TwoSamplesToys(TwoSamplesTest):
    """ Toys for Two-sample test 
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
        super().__init__ ( data1 , data2 )

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
        
        data1    = numpy.asarray     ( self.ecdf1.data () , dtype = float )
        data2    = numpy.asarray     ( self.ecdf1.data () , dtype = float )
        data     = numpy.concatenate ( ( data1 , data2 ) )
                
        results  = defaultdict(list)
        counters = self.counters
        
        from ostap.utils.progress_bar import progress_bar 
        for i in safe_loop ( progress_bar ( nToys , silent = silent , description = permutation_symbol ) ) : 

            numpy.random.shuffle ( data )
            data1 = data [    : n1 ]
            data2 = data [ n1 :    ] 
            
            ## transform input data 
            ecdf1, ecdf2 = prepare_data ( data1 , data2 ) 
            
            ks = kolmogorov_smirnov ( ecdf1 , ecdf2 )
            k  = kuiper             ( ecdf1 , ecdf2 )
            ad = anderson_darling   ( ecdf1 , ecdf2 )            
            cm = cramer_von_mises   ( ecdf1 , ecdf2 )             
            bj = berk_jones         ( ecdf1 , ecdf2 )             
            zk = ZK                 ( ecdf1 , ecdf2 )
            za = ZA                 ( ecdf1 , ecdf2 )
            zc = ZC                 ( ecdf1 , ecdf2 )
            w1 = wasserstein1       ( ecdf1 , ecdf2 )             
            w2 = wasserstein2       ( ecdf1 , ecdf2 )             
            r  = rosenblatt         ( ecdf1 , ecdf2 )             
            u2 = watson_u2          ( ecdf1 , ecdf2 )             
            mw = mann_whitney       ( ecdf1 , ecdf2 )             
            s  = savage             ( ecdf1 , ecdf2 )             
            
            counters [ 'KS'  ] += ks 
            counters [ 'K'   ] += k 
            counters [ 'AD'  ] += ad
            counters [ 'CM'  ] += cm
            counters [ 'BJ'  ] += bj
            counters [ 'ZK'  ] += zk
            counters [ 'ZA'  ] += za
            counters [ 'ZC'  ] += zc
            counters [ 'W1'  ] += w1 
            counters [ 'W2'  ] += w2 
            counters [ 'R'   ] += r
            counters [ 'U2'  ] += u2
            counters [ 'MW'  ] += mw
            counters [ 'S'   ] += s 
                        
            results  [ 'KS'  ] .append ( ks )    
            results  [ 'K'   ] .append ( k  ) 
            results  [ 'AD'  ] .append ( ad ) 
            results  [ 'CM'  ] .append ( cm ) 
            results  [ 'BJ'  ] .append ( bj ) 
            results  [ 'ZK'  ] .append ( zk ) 
            results  [ 'ZA'  ] .append ( za ) 
            results  [ 'ZC'  ] .append ( zc ) 
            results  [ 'W1'  ] .append ( w1 ) 
            results  [ 'W2'  ] .append ( w2 ) 
            results  [ 'R'   ] .append ( r  ) 
            results  [ 'U2'  ] .append ( u2 ) 
            results  [ 'MW'  ] .append ( mw ) 
            results  [ 'S'   ] .append ( s  ) 
            
        ## accumulate number of toys 
        self.__nToys += nToys 

        ECDF = Ostap.Math.ECDF 
        for key in results :
            data = results [ key ]
            if not data : continue
            if not key in self.__ecdfs : self.__ecdfs [ key ] = ECDF ( data , True ) ## complementary ECDF!
            else                       : self.__ecdfs [ key ]  .add  ( data2vct ( data ) ) 
            
        del results

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
        sigma  = str ( nsigma.toString ( '%%.2f %s %%-.2f' % plus_minus ) if float ( nsigma ) < 1000 else infinity_pos  ) 

        return ( what  ,
                 fmtv  % vs ,
                 fmt   % ( vm.value() , vm.error() ) ,
                 fmtv  % vr                          ,
                 fmtv  % vmn                         ,
                 fmtv  % vmx                         ,
                 ( '[%s]' % format_pow10 ( expo ) )  if expo else '' , pvalue , sigma )
    
    # =========================================================================
    ## Make a summary table
    def table ( self , title = '' , prefix = '' , width = 6 , precision = 4 , style = None ) :
        """ Make a summary table
        """
        import ostap.logger.table  as     T                 
        header = ( 'Statistics'  , ) + full_header
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
    def draw  ( self , what , option = '' , *options , **kwargs ) :
        """ Draw ECDF for toys & statistical estgimator 
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
        elif key in Keys [ 'BJ' ] and 'BJ' in self.ecdfs : 
            result = self.result  ( 'BJ' )
            ecdf   = self.ecdfs   [ 'BJ' ]
            ## logger.info ( 'Toy results for Berk-Jones  estimate' )             
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
        elif key in Keys [ 'W1' ] and 'W1' in self.ecdfs : 
            result = self.result  ( 'W1' )
            ecdf   = self.ecdfs   [ 'W1' ]
            ## logger.info ( 'Toy results for Wasserstein/1 estimate' ) 
        elif key in Keys [ 'W2' ] and 'W2' in self.ecdfs : 
            result = self.result  ( 'W2' )
            ecdf   = self.ecdfs   [ 'W2' ]
            ## logger.info ( 'Toy results for Wasserstein/2 estimate' ) 
        elif key in Keys [ 'R' ] and 'R' in self.ecdfs : 
            result = self.result  ( 'R' )
            ecdf   = self.ecdfs   [ 'R' ]
            ## logger.info ( 'Toy results for Rosenblatt estimate' ) 
        elif key in Keys [ 'U2' ] and 'U2' in self.ecdfs : 
            result = self.result  ( 'U2' )
            ecdf   = self.ecdfs   [ 'U2' ]
            ## logger.info ( 'Toy results for Watson U2 estimate' ) 
        elif key in Keys [ 'MW' ] and 'MW' in self.ecdfs : 
            result = self.result  ( 'MW' )
            ecdf   = self.ecdfs   [ 'MW' ]
            ## logger.info ( 'Toy results for Mann-Whitney estimate' ) 
        elif key in Keys [ 'S' ] and 'S' in self.ecdfs : 
            result = self.result  ( 'S' )
            ecdf   = self.ecdfs   [ 'S' ]
            ## logger.info ( 'Toy results for Savage estimate' ) 
        else :
            raise KeyError (  "draw: Invalid `what`:%s" % what )

        xmin , xmax = ecdf.xmin () , ecdf.xmax ()
        value       = result.statistics
        xmin        = min ( xmin , value )
        xmax        = max ( xmax , value )
        xmin , xmax = axis_range ( xmin , xmax , delta = 0.20 )

        kwargs [ 'xmin' ] = kwargs.get ( 'xmin' , xmin ) 
        kwargs [ 'xmax' ] = kwargs.get ( 'xmax' , xmax )

        result    = ecdf.draw  ( option , *options , **kwargs ) 
        line      = ROOT.TLine ( value  , 1e-3 , value , 1 - 1e-3 )
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




