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
from   ostap.core.core          import SE, VE, Ostap
from   ostap.math.base          import doubles, axis_range, numpy, np2raw    
from   ostap.math.models        import f1_draw
from   ostap.utils.basic        import numcpu, loop_items  
from   ostap.stats.gof_utils    import Estimators,Summary
import ostap.fitting.ds2numpy 
import ostap.fitting.roofit
import ROOT, math


from ostap.utils.timing import timing 
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
# ==============================================================================
## @class GoF1D
#  Goodness of 1D-fits 
#  @see https://doi.org/10.1111/1467-9868.00337
class GoF1D(Estimators) :
    """ Goodness-of-fit 1D-case 
    - see https://doi.org/10.1111/1467-9868.00337
    """
    def __init__ ( self           ,
                   pdf            ,
                   dataset        ,
                   cdf     = None ) :
        
        assert isinstance ( pdf     , PDF1            ) , 'Invalid type of `pdf`:%s'     % type ( pdf     )
        assert isinstance ( dataset , ROOT.RooDataSet ) , 'Invalid type of `dataset`:%s' % type ( dataset )

        
        vars = pdf.vars
        assert 1 == len ( vars )   , 'GoF1D: Only 1D-pdfs are allowed!'
        assert pdf.xvar in dataset , 'GoF1D: `xvar`:%s is not in dataset!' % ( self.xvar.name ) 

        original_cdf = cdf
        
        cdf_ok = cdf and callable ( cdf )

        self.__store = () 
        if not cdf_ok :
            cdf = pdf.cdf ()
            ## make atry to get underlying c++ cdf-function from PDF 
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
        self.__parameters = pdf.parameters ( dataset ) 

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

        ## evalute CDF for sorted data 
        cdf_data = self.clip ( self.__vct_cdf ( data ) ) 
        
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
    
    ## De-serialize th eobject 
    def __setstate__ ( self , state ) :
        """ De-serialize the object """         
        self.__pdf        = state.pop ( 'pdf'        ) 
        self.__parameters = state.pop ( 'parameters' ) 
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
    ## Clip input CDF arrays 
    def clip ( self , input ) :
        """ Clip input CDF arrays"""
        vmin , vmax = 1.e-12 , 1 - 1.e-10
        ## if numpy.min ( input ) < vmin or vmax < numpy.max ( input ) :
        ## logger.warning ( 'Adjust CDF' ) 
        return numpy.clip ( input , a_min = vmin , a_max = vmax )
        ## return input 
    
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
    ## Initialize Gof!D toys object :
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
            cdf_data = self.clip ( vct_cdf ( data ) )
            
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

    plot = Summary.draw 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================


