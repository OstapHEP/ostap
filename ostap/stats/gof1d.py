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
from   collections            import defaultdict, namedtuple
from   ostap.core.meta_info   import root_info 
from   ostap.fitting.funbasic import AFUN1
from   ostap.fitting.pdfbasic import PDF1
from   ostap.core.core        import SE, VE, Ostap, cidict_fun 
from   ostap.math.base        import doubles, axis_range  
from   ostap.math.models      import f1_draw
from   ostap.utils.basic      import numcpu, loop_items  
from   ostap.stats.gof_utils  import Estimators,Summary
import ostap.fitting.ds2numpy 
import ostap.fitting.roofit
import ROOT, math  
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import numpy as np
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    np = None
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
    n      = len ( cdf_data ) 
    result = max ( max ( ( i + 1.0 ) / n - Fi , Fi - float ( i ) / n ) for ( i, Fi )  in enumerate ( cdf_data )  ) ** 2  
    return math.sqrt ( result )
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
    flog    = math.log
    result  = sum ( ( i + 0.5 ) * flog ( Fi ) + ( n - i -  0.5 ) * flog ( 1 - Fi ) for ( i , Fi )  in enumerate ( cdf_data ) ) 
    result *= -2.0 / n
    result -= n
    return result 
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
    result  = sum ( ( Fi - ( i + 0.5 ) / n ) ** 2 for ( i, Fi ) in enumerate ( cdf_data ) ) 
    result += 12.0 / n
    return result#
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
    d_plus  = max ( ( i + 1.0 ) / n - Fi for ( i, Fi ) in enumerate ( cdf_data ) )
    d_minus = max ( Fi - ( i + 1.0 ) / n for ( i, Fi ) in enumerate ( cdf_data ) )
    result  = d_plus + d_minus  
    return result 
# =============================================================================
## Get ZK statististics
#  @code
#  cdf_data = ...
#  zk       = ZK ( cdf_data )
#  @endcode
#  @see https://doi.org/10.1111/1467-9868.00337
#  @param cdf_data sorted array of F0(X_i) - values of CDF at X data points
#  @return ZK statistics ZK 
def ZK  ( cdf_data ) :
    """ Get ZK statististics 
    - `cdf_data` : sorted array of F0(X_i) - values of CDF at X data points
    >>> cdf_data =...
    >>> zk       = ZK ( cdf_data )    
    """
    n      = len ( cdf_data ) 
    flog   = math.log 
    result = max ( ( i     + 0.5 ) * flog ( ( i + 0.5     ) / ( n *       Fi   ) ) +
                   ( n - i - 0.5 ) * flog ( ( n - i - 0.5 ) / ( n * ( 1 - Fi ) ) ) for ( i , Fi ) in enumerate ( cdf_data ) )
    return result 
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
    flog    = math.log
    result  = sum ( flog ( Fi ) / ( n - i - 0.5 ) + flog ( 1 - Fi ) / ( i + 0.5 ) for ( i , Fi )  in enumerate ( cdf_data ) )
    result *= -1 
    return result    
# =============================================================================
## Get ZC statististics
#  @code
#  cdf_data = ...
#  zc       = ZC ( cdf_data )
#  @endcode
#  @see https://doi.org/10.1111/1467-9868.00337
#  @param cdf_data sorted array of F0(X_i) - values of CDF at X data points
#  @return ZC statistics ZC
def ZC  ( cdf_data ) :
    """ Get ZC statististics 
    - `cdf_data` : sorted array of F0(X_i) - values of CDF at (sorted) X data points
    >>> cdf_data =...
    >>> zc       = ZC ( cdf_data )    
    """
    n      = len ( cdf_data ) 
    flog   = math.log
    result = sum ( ( flog ( ( 1.0 / Fi - 1 ) / ( ( n - 0.5 ) / ( i + 0.25 ) - 1 ) ) ) ** 2 for ( i , Fi ) in enumerate ( cdf_data ) )
    return result 
# ==============================================================================
## @class GoF1D
#  Goodness of 1D-fits 
#  @see https://doi.org/10.1111/1467-9868.00337
class GoF1D(Estimators) :
    """ Goodness-of-fit 1D-case 
    - see https://doi.org/10.1111/1467-9868.00337
    """
    def __init__ ( self    ,
                   pdf     ,
                   dataset ) :
        
        assert isinstance ( pdf     , PDF1            ) , 'Invalid type of `pdf`:%s'     % type ( pdf     )
        assert isinstance ( dataset , ROOT.RooDataSet ) , 'Invalid type of `daatset`:%s' % type ( dataset )
        
        vars = pdf.vars
        assert 1 == len ( vars )   , 'GoF1D: Only 1D-pdfs are allowed!'
        assert pdf.xvar in dataset , 'GoF1D: `xvar`:%s is not in dataset!' % ( self.xvar.name ) 
        
        cdf = pdf.cdf ()
        if hasattr ( pdf.pdf , 'setPars'  ) : pdf.pdf.setPars() 
        if hasattr ( pdf.pdf , 'function' ) :
            if hasattr ( pdf.pdf , 'setPars'  ) : pdf.pdf.setPars() 
            fun = pdf.pdf.function()
            if hasattr ( fun , 'cdf' ) :
                try :                    
                    a = fun.cdf ( 0.0 )
                    self.__store = pdf, fun 
                    def cdf ( x ) : return fun.cdf ( x ) 
                except TypeError : pass 

        self.__cdf   = cdf
        self.__xmnmx = pdf.xminmax()
        
        ## vectorized version of CDF 
        self.__vct_cdf = np.vectorize ( cdf )
        
        ## data in form of numpy sructured array
        varname = pdf.xvar.name 
        data    = dataset.tonumpy ( varname ) [ varname ] 

        ## sorted data 
        self.__data     = np.sort  ( data )

        ## empirical CDF function 
        self.__ecdf     = Ostap.Math.ECDF ( data2vct ( self.__data ) )

        ## evalute CDF for sorted data 
        self.__cdf_data = self.__vct_cdf ( self.__data )

        self.__estimators = {
            'KS'  : kolmogorov_smirnov ( self.__cdf_data ) , 
            'K'   : kuiper             ( self.__cdf_data ) , 
            'AD'  : anderson_darling   ( self.__cdf_data ) , 
            'CM'  : cramer_von_mises   ( self.__cdf_data ) , 
            'ZK'  : ZK                 ( self.__cdf_data ) , 
            'ZA'  : ZA                 ( self.__cdf_data ) , 
            'ZC'  : ZC                 ( self.__cdf_data ) , 
        }
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
        return len ( self.__data ) 
    
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
    def ecdf ( self ) :
        """`ecdf` : empirical CDF for data"""
        return self.__ecdf
    
    @property
    def data ( self ) :
        """`data` : sorted input data (as numpy array)"""
        return self.__data

    @property
    def cdf_data ( self ) :
        """`cdata` : vector of CDF(x) data (as numpy array)"""
        return self.__cdf_data 

    # =========================================================================
    ## Get Kolmogorov-Smirnov statistiscs 
    @property 
    def kolmogorov_smirnov_estimator ( self ) :
        """ Get Kolmogorov-Smirnov statistiscs KS
        """
        return self.__estimators.get('KS',None) 

    # =========================================================================
    ## Get Kuiper  statististics 
    @property 
    def kuiper_estimator ( self ) :
        """ Get Kuiperstatistics
        """        
        return self.__estimators.get('K', None ) 
                
    # =========================================================================
    ## Get Anderson-Darling  statistiscs 
    @property 
    def anderson_darling_estimator ( self ) :
        """ Get Anderson-Darling statistiscs 
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
    
    __repr__ = Estimators.table
    __str__  = Estimators.table
    
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

        xmin , xmax = axis_range ( xmin , xmax , delta = 0.10 )
        
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
#  Check Goodness of 1D-fits using toys 
class GoF1DToys(GoF1D,Summary) :
    """ Check Goodness-of-fit with toys (1D-case) 
    """
    ## result of GoF-toys 
    Result = namedtuple ( 'Result' , 'statistics counter pvalue nsigma' )
    # =========================================================================
    def __init__ ( self              ,
                   pdf               ,
                   dataset           ,
                   nToys     = 1000  ,
                   parallel  = False ,
                   nSplit    = 0     , 
                   silent    = False ) :
        
        assert isinstance ( nToys , int ) and 0 < nToys , "Invalid `nToys` argument!"

        ## initialize the base
        GoF1D.__init__ ( self , pdf , dataset ) 

        self.__pdf      = pdf
        self.__dataset  = dataset
        
        self.__counters = defaultdict(SE) 
        self.__ecdfs    = {}
        
        self.__nToys    = 0
        
        if 0 < nToys :
            self.run ( nToys , parallel = parallel , silent = silent , nSplit = nSplit ) 

    # ===============================================================================
    ## run toys 
    def run ( self , nToys = 1000 , parallel = False , silent = False , nSplit = 0 ) :
        """ Run toys 
        """
        assert isinstance ( nToys , int ) and 0 < nToys , "Invalid `nToys` argument!"

        if parallel :
            from ostap.parallel.parallel_gof1d import parallel_gof1dtoys as parallel_toys 
            self += parallel_toys ( self.__pdf            ,
                                    self.__dataset        ,
                                    nToys    = nToys      ,
                                    nSplit   = nSplit     ,
                                    silent   = True       ,
                                    progress = not silent )
            return 

        varname = self.__pdf.xvar.name        
        results  = defaultdict(list)
        counters = self.counters 
        vct_cdf  = self.vcdf
        
        from ostap.utils.progress_bar import progress_bar 
        for i in progress_bar ( nToys , silent = silent , description = 'Toys:') :
        
            dset     = self.__pdf.generate ( self.N  , sample = True )
            data     = dset.tonumpy ( varname ) [ varname ] 
            data     = np.sort ( data )
            cdf_data = vct_cdf ( data )
            
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

            ## delete data
            if isinstance ( dset , ROOT.RooDataSet ) : 
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
    ## Get Kolmogorov-Smirnov statistiscs 
    @property 
    def kolmogorov_smirnov ( self ) :
        """ Get Kolmogorov-Smirnov statistiscs KS
        """
        return self.result ( 'KS' ) 
    
    # ===============================================
    ## Get Anderson-Darling  statistiscs 
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

    # =========================================================================
    ## Get Kuiper statististics 
    @property 
    def kuiper ( self ) :
        """ Get Kuiper statistics 
        """        
        return self.result ( 'K' ) 

    table    = Summary.table 
    __repr__ = Summary.table
    __str__  = Summary.table

    # =========================================================================
    ## merge two objects:
    def merge ( self , other ) :
        self += other
        return self
    
    ## merge two objects:
    def __iadd__ ( self , other ) :
        """ Merge two objects """        
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
        return self 

    plot = Summary.draw 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================


