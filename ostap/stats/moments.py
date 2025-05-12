#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  moments.py 
#  Utilities to get moments for various functions/distributions/pdfs
#  - Moment
#  - CentralMoment
#  - StdlMoment
#  - Mean
#  - Variance 
#  - RMS 
#  - Skewness
#  - Kurtosis
#  - Median
#  - Quantile 
#  - Mode
#  - Width
#  - Mode
#  - symmetric and asymmetric "confidence intervals"
#  For these quantities numerical integration and root-findinng are  used.
#
#  All objects exists as classes/functors and as standalone simple functions
#  - moment
#  - central_moment
#  - std_moment
#  - mean
#  - variance 
#  - rms 
#  - skewness
#  - kurtosis
#  - median
#  - quantile 
#  - mode
#  - width
#  - FWHM
#  - cl_symm and sl_asymm 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06  
# =============================================================================
"""Utilities to get moments for various functions/distributions/pdfs
- Moment
- CentralMoment
- StdMoment
- Mean
- Variance 
- RMS 
- Skewness
- Kurtosis
- Median
- Quantile 
- Mode
- Width
- symmetric and asymmetric ``confidence intervals''
For these quantities numerical integration and root-finding are used.

All objects exists as classes/functors and as standalone simlpe functions
- moment
- central_moment
- std_moment
- mean
- variance 
- rms 
- skewness
- kurtosis
- median
- quantile 
- mode
- width
- fwhm
- cl_symm and sl_asymm 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-06-06"
__all__     = (
    ##
    ## stat-quantities 
    "Moment"        , ## calculate N-th moment of functions/distributions, etc 
    "CentralMoment" , ## calculate N-th central moment of functions/distributions
    "StdMoment"     , ## calculate N-th standartized moment of functions/distributions
    "Mean"          , ## "mean"     for functions/distributions, etc 
    "Variance"      , ## calculate "variance" for functions/distributions, etc 
    "RMS"           , ## calculate "RMS"      for functions/distributions, etc 
    "Skewness"      , ## calculate "skewness" for functions/distributions, etc 
    "Kurtosis"      , ## calculate "kurtosis" for functions/distributions, etc 
    "Median"        , ## calculate "median"   for functions/distributions, etc 
    "Quantile"      , ## calculate "quantile" for functions/distributions, etc 
    "Mode"          , ## calculate "mode"     for functions/distributions, etc 
    "Width"         , ## calculate "width"    for functions/distributions, etc 
    "CL_symm"       , ## calcualte symmetrical confidence intervals            
    "CL_asymm"      , ## calcualte asymmetrical confidence intervals           
    ##
    ## stat-quantities   
    "moment"        , ## calculate N-th moment of functions/distributions, etc 
    "central_moment", ## calculate N-th moment of functions/distributions, etc 
    "mean"          , ## calculate "mean"     for functions/distributions, etc 
    "variance"      , ## calculate "variance" for functions/distributions, etc 
    "rms"           , ## calculate "RMS"      for functions/distributions, etc 
    "skewness"      , ## calculate "skeness"  for functions/distributions, etc 
    "kurtosis"      , ## calculate "kurtosis" for functions/distributions, etc 
    "median"        , ## calculate "median"   for functions/distributions, etc 
    "quantile"      , ## calculate "quantile" for functions/distributions, etc 
    "mode"          , ## calculate "mode"     for functions/distributions, etc 
    "width"         , ## calculate "width"    for functions/distributions, etc
    "fwhm"          , ## calculate "fwhm"     for functions/distributions, etc    
    "cl_symm"       , ## calculate symmetrical  confidence intervals            
    "cl_asymm"      , ## calculate asymmetrical confidence intervals           
    ##
    ) 
# =============================================================================
from ostap.core.ostap_types import integer_types, num_types
from ostap.math.base        import pos_infinity, neg_infinity
from ostap.core.core        import Ostap
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.moments' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
## @class BaseMoment
#  Base class for calculation of variosu moments 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class BaseMoment(object) :
    """Calculate the N-th moment for the distribution
    
    >>> xmin,xmax = 0,math.pi 
    >>> mean  = Moment(1,xmin,xmax)  ## specify min/max
    >>> value = mean  ( math.sin )
    """
    ## constructor
    def __init__ ( self , N , xmin , xmax , err = False , *args ) :
        """Contructor 
        """
        assert isinstance ( N , integer_types ) and  0 <= N , \
               'BaseMoment: illegal order: %s' % N 
        
        self.__N     = N 
        self.__xmin  = float ( xmin ) if isinstance ( xmin , integer_types ) else xmin 
        self.__xmax  = float ( xmax ) if isinstance ( xmax , integer_types ) else xmax 
        self.__err   = err
        self.__args  = args

    # =========================================================================
    ## get the normalized moment
    def moment ( self , K , func , center , *args ) :
        """Get the normalized momentum
        """
        assert isinstance ( K , integer_types ) and 0 <= K , \
               'Invalid type/valeu for K=%s' %k
        
        m0 = self.normalization ( func , *args )
        mK = self.umoment ( K , func , center = center , *args )
        return mK / float ( m0 ) 

    # =========================================================================
    ## get the unormalized moment
    def umoment ( self , K , func , center  , *args ) :
        """Get the unormalized moment
        """
        if 0 == K : return self.normalization ( func , *args ) 

        x0  = float ( center ) 
        fnK = lambda x , *a  : func ( x , *a ) * ( ( x - x0 ) ** K )
        ##
        return self.integral ( fnK , self.xmin , self.xmax , *args )

    # =========================================================================
    ## get the normalization integral within the specified range 
    def normalization ( self, func , *args ) :
        """get the normalization integral within the specified range 
        """        
        return self.integral ( func , xmin = self.xmin , xmax = self.xmax , *args )

    # =========================================================================
    ## get the mean value
    def mean ( self , func , *args ) :
        """Get the mean value 
        """
        return self.moment ( 1 , func , center = 0.0 , *args )
    
    # =========================================================================
    ## get the  central moment
    def central_moment ( self , K , func , *args ) :
        """Get the central moment
        """
        ##
        if   0 == K : return 1.0 
        elif 1 == K : return 0.0
        ##
        m0 = self.normalization (     func ,                *args ) 
        m1 = self.umoment       ( 1 , func , center = 0.0 , *args )
        mu = float ( m1 ) / float ( m0 ) 
        mK = self.umoment       ( K , func , center = mu  , *args )        
        ##
        return mK / float ( m0 )

    # =========================================================================
    ## get the standartizeds central moment
    def std_moment ( self , K , func , *args ) :
        """Get the standartized central moment
        """
        ##
        if   0 == K : return 1.0 
        elif 1 == K : return 0.0
        elif 2 == K : return 1.0
        ##
        m0 = self.normalization (     func ,                *args ) 
        m1 = self.umoment       ( 1 , func , center = 0.0 , *args )
        mu = float ( m1 ) / float ( m0 ) 
        mK = self.umoment       ( K , func , center = mu  , *args )        
        m2 = self.umoment       ( 2 , func , center = mu  , *args )        
        ##
        return mK / pow ( m2 , 0.5 * K ) 
    
    # =========================================================================
    ## get the variance 
    def variance      ( self , func , *args ) :
        """Get the  variance"""
        return self.central_moment ( 2 , func , *args )

    # =========================================================================
    ## get the RMS
    def rms          ( self , func , *args ) :
        """Get the  RMS"""
        return self.variance ( func , *args ) **0.5 

    # =========================================================================
    ## get the skewness 
    def skewness      ( self , func , *args ) :
        """Get the skewness"""
        return self.std_moment ( 3 , func , *args )

    # =========================================================================
    ## get the (excess) kurtosis 
    def kurtosis      ( self , func , *args ) :
        """Get the (excess) kurtosis """
        return self.std_moment ( 4 , func , *args ) - 3.0 
    
    # =========================================================================
    ## integrate the function between xmin and xmax 
    def integral ( self , func , xmin , xmax , *args ) :
        """Integrate the function between xmin and xmax"""
        from ostap.math.integral import IntegralCache 
        integrator = IntegralCache ( func , xmin , err = self.err , args = args )
        return integrator ( xmax , *args )

    # ==========================================================================
    ## calculate the median
    def median ( self , func , *args ) :

        ## need to know the integral
        from ostap.math.integral import Integral

        iint   = Integral      ( func ,  self.xmin , err = False ,  args = args )
        half   = 2.0 / iint    ( self.xmax ) 

        ifun   = lambda x : iint( x ) * half - 1.0

        from ostap.math.rootfinder import findroot
        
        ## @see https://en.wikipedia.org/wiki/Median#Inequality_relating_means_and_medians
        try: 
            meanv = self.mean ( func , *args )
            sigma = self.rms  ( func , *args )
            import math
            xmn   = meanv - 1.5 * sigma ## use 1.5 instead of 1 
            xmx   = meanv + 1.5 * sigma ## use 1.5 instead of 1
            #
            if isinstance ( self.xmin , float ) : xmn = max ( xmn , self.xmin ) 
            if isinstance ( self.xmax , float ) : xmx = min ( xmx , self.xmax )
            #
            result = findroot ( ifun , xmn       , xmx       , maxiter = 500 )
        except :
            result = findroot ( ifun , self.xmin , self.xmax , maxiter = 500 )
            
        return result

    @property
    def args ( self ) :
        "``args'' - other arguments for the function call"
        return self.__args
    
    @property
    def  err( self ) :
        "``err''- evaluate the error/uncertanty?"
        return self.__err

    @property
    def  xmin ( self ) :
        "``xmin''- low edge of the interval"
        return self.__xmin
    @property
    def  xmax ( self ) :
        "``xmax''- high edge of the interval"
        return self.__xmax
    
    @property
    def N ( self ) :
        "``N''- the moment order"
        return self.__N

# =============================================================================    
## @class Moment
#  Calculate the N-th moment for the distribution 
#  @code
#   xmin,xmax = 0,math.pi 
#   mean  = Moment(1,xmin,xmax)  ## specify min/max
#   value = mean  ( math.sin )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class Moment(BaseMoment) :
    """Calculate the N-th moment for the distribution
    
    >>> xmin,xmax = 0,math.pi 
    >>> mean  = Moment(1,xmin,xmax)  ## specify min/max
    >>> value = mean  ( math.sin )
    """
    ## constructor
    def __init__ ( self , N , xmin , xmax , err = False , center = 0.0 , *args ) :
        """Contructor 
        """
        BaseMoment.__init__ ( self , N = N , xmin = xmin , xmax = xmax , err = err , *args )
        self.__center = float ( center )   

    # =========================================================================
    ## The main method: get the moment of the function/distribution 
    def __call__ ( self , func , *args ) :
        ##
        args   = args if args else self.args
        ##
        return self.moment ( self.N , func , self.__center , *args )

    ## print it!
    def __str__ ( self ) :
        return "Moment(%d,%s,%s,%s,%s)" % ( self.N      ,
                                            self.xmin   ,
                                            self.xmax   ,
                                            self.err    ,
                                            self.center )
    @property
    def x0 ( self ) :
        "``x0''- the center"
        return self.__center 
    @property
    def center ( self ) :
        "``center''- the center (same as ``x0'')"
        return self.__center
    
# =============================================================================
## @class CentralMoment
#  Calculate the N-th central moment for the distribution 
#  @code
#   xmin,xmax = 0,math.pi 
#   mc        = CentralMoment(1,xmin,xmax)  ## specify min/max
#   value     = mome  ( math.sin )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class CentralMoment(BaseMoment) :
    """Calculate the N-th central moment for the distribution
    
    >>> xmin,xmax = 0,math.pi 
    >>> mc        = CentralMoment(1,xmin,xmax)  ## specify min/max
    >>> value     = mc  ( math.sin )    
    """
    ## constructor
    def __init__ ( self , N , xmin , xmax , err = False, *args ) :
        """Contructor 
        """
        BaseMoment.__init__ ( self        ,
                              N    = N    ,
                              xmin = xmin ,
                              xmax = xmax ,
                              err  = err  , *args ) 
        
    ## calculate the central moment
    def __call__ ( self , func , *args ) :
        ##
        args   = args if args else self.args
        ##
        if   0 == self.N : return 1.0
        elif 1 == self.N : return 0.0
        ##
        return self.central_moment ( self.N , func , *args )

    ## print it!
    def __str__ ( self ) :
        return "CentralMoment(%d,%s,%s,%s)" % ( self.N    ,
                                                self.xmin ,
                                                self.xmax ,
                                                self.err  )

# =============================================================================
## @class StdMoment
#  Calculate the N-th standartizez moment for the distribution 
#  @code
#   xmin,xmax = 0,math.pi 
#   mc        = StdMoment(1,xmin,xmax)  ## specify min/max
#   value     = mome  ( math.sin )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class StdMoment(CentralMoment) :
    """Calculate the N-th central moment for the distribution
    
    >>> xmin,xmax = 0,math.pi 
    >>> mc        = CentralMoment(1,xmin,xmax)  ## specify min/max
    >>> value     = mc  ( math.sin )    
    """
    ## constructor
    def __init__ ( self , N , xmin , xmax , err = False, *args ) :
        """Contructor 
        """
        CentralMoment.__init__ ( self , N = N , xmin = xmin , xmax = xmax , err = err , *args ) 

    ## calculate the central moment
    def __call__ ( self , func , *args ) :
        ##
        args   = args if args else self.args
        ##
        if   0 == self.N : return 1.0
        elif 1 == self.N : return 0.0
        elif 2 == self.N : return 1.0 
        ##
        return self.std_moment ( self.N , func , *args )

    ## print it!
    def __str__ ( self ) :
        return "StdMoment(%d,%s,%s,%s)" % ( self.N    ,
                                            self.xmin ,
                                            self.xmax ,
                                            self.err  )
    
# =============================================================================
## @class Mean
#  Calculate the mean-value for the distribution 
#  @code
#   xmin,xmax = 0,math.pi 
#   mean  = Mean(xmin,xmax)  ## specify min/max
#   value = mean  ( math.sin )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class Mean(Moment) :
    """Calculate the N-th moment for the distribution or function 
    
    >>> xmin,xmax = 0,math.pi 
    >>> mean  = Mean ( xmin , xmax )  ## specify min/max
    >>> value = mean ( math.sin    )
    
    """
    
    def __init__ ( self , xmin , xmax , err = False , *args ) :
        Moment.__init__ ( self , N = 1 , xmin = xmin , xmax = xmax , err = err , center = 0.0 , *args )

    def __str__ ( self ) :
        return "Mean(%s,%s,%s)" % ( self.xmin ,
                                    self.xmax ,
                                    self.err  )                                            

# =============================================================================
## @class Variance
#  Calculate the variance for the distribution or function  
#  @code
#   xmin,xmax = 0,math.pi 
#   variance  = Variance ( xmin,xmax )  ## specify min/max
#   value     = variance ( math.sin  )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class Variance(CentralMoment) :
    """Calculate the variance for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> variance  = Variance ( xmin,xmax )  ## specify min/max
    >>> value     = variance ( math.sin  )
    """
    def __init__ ( self , xmin , xmax , err = False , *args ) :
        CentralMoment.__init__ ( self , N = 2 , xmin = xmin , xmax = xmax , err = err , *args )

    def __str__ ( self ) :
        return "Variance(%s,%s,%s)" % ( self.xmin ,
                                        self.xmax ,
                                        self.err  )                                            

# =============================================================================
## @class RMS
#  Calculate the variance for the distribution or function  
#  @code
#   xmin,xmax = 0,math.pi 
#   rms       = RMS ( xmin,xmax )  ## specify min/max
#   value     = rms ( math.sin  )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06
class RMS(Variance) :
    """Calculate the RMS for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> rms       = RMS ( xmin,xmax )  ## specify min/max
    >>> value     = rms ( math.sin  )
    """
    def __init__ ( self , xmin , xmax , err = False , *args ) :
        Variance.__init__ ( self , xmin = xmin , xmax = xmax , err = err , *args )
        
    ## calculate the variance 
    def __call__ ( self , func , *args ) :
        """Calculate the RMS for the distribution or function          
        """
        ##
        args   = args if args else self.args
        ##
        var2 = Variance.__call__ ( self , func , *args )
        import ostap.math.math_ve as ME 
        return ME.sqrt ( var2 ) 

    def __str__ ( self ) :
        return "RMS(%s,%s,%s)" % ( self.xmin ,
                                   self.xmax ,
                                   self.err  )                                            

# =============================================================================
## @class Skewness
#  Calculate the skewness for the distribution or function  
# \f$ k = \frac{\mu_3}{\sigma^3} \f$
#  @code
#   xmin,xmax = 0,math.pi 
#   skew     = Skewness ( xmin,xmax )  ## specify min/max
#   value     = skew ( math.sin  )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-08-09
class Skewness(StdMoment) :
    """
    Calculate the variance for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> skew      = Skewness ( xmin,xmax )  ## specify min/max
    >>> value     = skew     ( math.sin  )
    """
    def __init__ ( self , xmin , xmax , err = False , *args ) :
        StdMoment.__init__ ( self , N = 3 ,  xmin = xmin , xmax = xmax , err = err , *args )
        
    def __str__ ( self ) :
        return "Skewness(%s,%s,%s)" % ( self.xmin ,
                                        self.xmax ,
                                        self.err  )                                            

# =============================================================================
## @class  Kurtosis
#  Calculate the (excess) kurtosis for the distribution or function
#  \f$ k = \frac{\mu_4}{\sigma_2}-3\f$
#  @code
#   xmin,xmax = 0,math.pi 
#   kurt      = Kurtosis ( xmin,xmax )  ## specify min/max
#   value     = kurt ( math.sin  )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-08-09
class Kurtosis(StdMoment) :
    """Calculate the variance for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> kurt      = Kurtosis ( xmin,xmax )  ## specify min/max
    >>> value     = kurt     ( math.sin  )
    """
    def __init__ ( self , xmin , xmax , err = False , *args ) :
        StdMoment.__init__ ( self , N = 4 , xmin = xmin , xmax = xmax , err = err , *args )
        
    ## calculate the kurtosis
    def __call__ ( self , func , *args ) :
        ##
        return StdMoment.__call__ ( self , func , *args ) - 3.0 

    def __str__ ( self ) :
        return "Kurtosis(%s,%s,%s)" % ( self.xmin ,
                                        self.xmax ,
                                        self.err  )                                            


# =============================================================================
## @class Median 
#  Calculate the median for the distribution or function  
#  @code
#  xmin,xmax = 0,math.pi 
#  median    = Median ( xmin,xmax )  ## specify min/max
#  value     = median ( math.sin  )
#  @endcode 
#  @see https://en.wikipedia.org/wiki/Median#Inequality_relating_means_and_medians
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-12
class Median(RMS) :
    """Calculate median for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> median    = Median ( xmin,xmax )  ## specify min/max
    >>> value     = median ( math.sin  )
    """
    def __init__ ( self , xmin , xmax ) :
        RMS.__init__ ( self , xmin , xmax , err = False )
        
    ## calculate the median 
    def __call__ ( self , func , *args ) :
        return self.median ( func ,  *args )

    def __str__ ( self ) :
        return "Median(%s,%s)" % ( self.xmin , self.xmax )
    
# =============================================================================
## get the quantile
#  Calculate the quantile for the distribution or function  
#  @code
#  xmin,xmax = 0,math.pi 
#  quantile  = Quantile ( 0.1 , xmin,xmax )  ## specify min/max
#  value     = quantile ( math.sin  )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-12
class Quantile(Median) :
    """Calculate quantiles for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> quantile  = Quantile ( 0.1 , xmin,xmax )  ## specify min/max
    >>> value     = quantile ( math.sin  )
    """
    def __init__ ( self , Q , xmin , xmax ) :
        Median.__init__ ( self , xmin , xmax )
        #
        assert 0 <= Q <=  1, 'Quantile is invalid %s' % Q 

        self.__Q = float( Q ) 
        
    def __str__ ( self ) :
        return "Quantile(%s,%s,%s)" % ( self.Q , self.xmin , self.xmax )

    ## calculate the median 
    def __call__ ( self , func , *args ) :
        ##

        if    0.5 == self.Q : return self.median ( func , *args ) 
        elif  0.0 == self.Q : return self.xmin
        elif  1.0 == self.Q : return self.xmax

        ## need to know the integral
        from ostap.math.integral import IntegralCache
        iint = IntegralCache ( func, self.xmin, err = False ,  args = args )
        quan = 1.0 / iint    (  self.xmax ) / self.Q 
        
        ifun   = lambda x : iint( x ) * quan - 1.0

        xmn = self.xmin
        xmx = self.xmax

        p   = 0.5
        l   = 0.5

        ## make some bracketing before next step 
        while ( not isinstance ( xmn , float ) ) or ( not isinstance ( xmx , float ) ) or l>0.1 :   
        
            l /= 2
            mm = Median ( xmn , xmx )
            m  = mm .median ( func , *args )
            
            if   self.Q < p :
                xmn   = xmn 
                xmx   = float( m ) 
                p    -= l 
            elif self.Q > p :
                xmn   = float ( m ) 
                xmx   = xmx
                p    +=  l  
            else : return m               ## RETURN 

        ## finally, calculate quantile
        from ostap.math.rootfinder import findroot 
        result = findroot ( ifun , xmn , xmx )
            
        return result

    @property
    def  Q ( self ) :
        "``Q''-quantile level"
        return self.__Q
    

# =============================================================================
## @class Mode
#  Calculate the mode for the distribution or function  
#  @code
#  xmin,xmax = 0,math.pi 
#  mode      = Mode ( xmin,xmax )  ## specify min/max
#  value     = mode ( math.sin  )
#  @endcode 
#  @see https://en.wikipedia.org/wiki/Mode_(statistics)#Comparison_of_mean,_median_and_mode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-12
#  @attention it is not very efficient, could be (and should be) improved 
class Mode(Median) :
    """ Calculate the mode for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> mode      = Mode ( xmin,xmax )  ## specify min/max
    >>> value     = mode ( math.sin  )
    """
    def __init__ ( self , xmin , xmax ) :
        Median.__init__ ( self , xmin , xmax )
        
    ## calculate the mode 
    def __call__ ( self , func , *args ) :

        ## median 
        _median   = self.median ( func , *args )

        ## mean
        _mean     = self.mean   ( func , *args )

        ## rms
        _sigma    = self.rms    ( func , *args )

        _imin = _median - 1.75 * _sigma
        _imax = _median + 1.75 * _sigma
        
        mn = max ( _imin , self.xmin )
        mx = min ( _imax , self.xmax )

        mm = Mode ( mn , mx )
        _median   = mm.median ( func , *args )
        _sigma    = mm.rms    ( func , *args )
        
        ## readjust the interval 
        mn        = max ( mn , _median - 1.75 * _sigma )
        mx        = min ( mx , _median + 1.75 * _sigma )

        mm        = Mode ( mn , mx ) 
        _mean     = mm.mean   ( func , *args )
        _median   = mm.median ( func , *args )

        ## use the empirical approximation 
        _mode = 3.0 * _median - 2.0 * _mean
        _fnm  = func ( _mode ,  *args )

        ## make the interval small enough
        for i in range ( 10 ) : 

            if abs ( mx - mn ) < 6 * _sigma : break
            
            
            if abs ( mx - mn ) * 20 < abs ( self.xmax -  self.xmin ) : break
            
            if mn < _mode < mx and func ( mn , *args ) < _fnm and func ( mx , *args ) < _fnm :
                
                x1 = 0.5 * ( mn + _mode )
                x2 = 0.5 * ( mx + _mode )
                
                f1 = func ( x1 , *args )
                f2 = func ( x2 , *args )
                
                if f1 < _fnm and f2 < _fnm :
                    mn = x1
                    mx = x2
                elif f1 < _fnm : mn = x1
                elif f2 < _fnm : mx = x2
                else :
                    break

                mm = Mode ( mn , mx )
                _mean     = mm.mean   ( func , *args )
                _median   = mm.median ( func , *args )
                _new_mode = 3.0 * _median - 2.0 * _mean

                _new_fnm  = func ( _new_mode , *args ) 
                if  mn  < _new_mode < mx and func ( mn , *args ) < _new_fnm and func ( mx , *args ) < _new_fnm :
                    _mode = _new_mode
                    _fnm  = _new_fnm
                    
            else :
                
                break

            
        ifun = lambda x,*a : -1.0 * float( func ( x , *a ) )
        
        from ostap.math.minimize import minimize_scalar as _ms 
        result = _ms (
            ifun                     , 
            ## x0     = float ( m0 ) ,
            method = 'bounded'       , 
            bounds = [ mn , mx ]     ,
            args   = args )
        
        return result.x
    
    def __str__ ( self ) :
        return "Mode(%s,%s)" % ( self.xmin , self.xmax )

# =============================================================================
## @class Width
#  Calculate the full width at half height for the distribution or function  
#  @code
#  xmin,xmax = 0,math.pi 
#  width     = Width ( xmin,xmax )  ## specify min/max
#  x1 , x2   = width ( math.sin  )
#  fwhm      = x2 - x1
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-12
class Width(Mode) :
    """Calculate the mode for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> width     = Width ( xmin,xmax )  ## specify min/max
    >>> x1 , x2   = width ( math.sin )
    >>> fwhm      = x2 - x1
    """
    def __init__ ( self , xmin , xmax , height_factor = 0.5 ) :
        Mode.__init__ ( self , xmin , xmax )
        self._hfactor = height_factor
        
    ## calculate the width
    def __call__ ( self , func , mode = None , *args ) :
        ##

        ## mode is specified 
        if isinstance ( mode , float )        and \
           self.xmin < mode < self.xmax       and \
           func ( self.xmin ) < func ( mode ) and \
           func ( self.xmax ) < func ( mode ) :  m0 = mode
        else :
            ## mode needs to be calculated
            m0 = Mode.__call__ ( self , func , *args )
            
        ## function  value at the maximum
        v0      = float ( func ( m0 , *args ) ) 

        ## half height 
        vheight = 1.0 * v0 * self._hfactor
        
        ifun = lambda x,*a : float ( func ( x , *a ) ) - vheight

        from ostap.math.rootfinder import findroot
        x1 = findroot ( ifun , self.xmin , m0        , args = args )
        x2 = findroot ( ifun , m0        , self.xmax , args = args ) 
        
        return x1 , x2

    def __str__ ( self ) :
        return "Width(%s,%s,%s)" % ( self.xmin , self.xmax , self._hfactor)
    
# =============================================================================
## @class CL_symm
#  Calcualate symmetic confidence interval around x0
#  for the given function on (xmin,xmax) interval
#  function is assumed to be zero outside the interval
#  @code
#  fun = lambda x : exp( -0.5*x*x)
#  reg = CL_symm ( 0.68 , -10 , 10 )
#  print reg ( fun )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-08-03
class CL_symm(object) :
    """Calculate symmetic confidence interval around x0
    for the given function on (xmin,xmax) interval
    function is assumed to be zero outside the interval
    >>> fun = lambda x : exp( -0.5*x*x)
    >>> reg = CL_symm ( 0.68 , -10 , 10 )
    >>> print reg ( fun )
    """
    def __init__ ( self , prob ,  xmin , xmax , x0 = None , *args ) :
        
        if not 0.0 < prob < 1.0 :
            raise AttributeError ("Invalid value of prob/CL=%g" % prob)

        self.__prob = prob 
        self.__xmin = float ( xmin ) if isinstance ( xmin , integer_types ) else xmin 
        self.__xmax = float ( xmax ) if isinstance ( xmax , integer_types ) else xmax 
        self.__x0   = float ( x0   ) if isinstance ( x0   , integer_types ) else x0  
        self.__args = args
        
    def __str__ ( self ) :
        return "CL_sym(%s,%s,%s,%s)" % ( self.prob ,
                                         self.xmin ,
                                         self.xmax ,
                                         self.x0   )

    def __call__ ( self , func , *args ) :

        ## additional arguments
        args   = args if args else self.args

        #
        ## define integration rules
        #
        if hasattr ( func , 'integral' ) :
            _integral_ = lambda f , low , high : f.integral (      low , high , *args )
        else                             :
            from ostap.math.integral import integral 
            _integral_ = lambda f , low , high :   integral ( f  , low , high , *args )

        #
        ## xmin/max
        #
        xmin,xmax = self.xmin, self.xmax
        
        #
        ## calculate x0 as "mean"-value
        #
        x0 = self.x0 
        if x0 is None :
            if hasattr ( func , 'mean' ) : x0 = func.mean()
            else                         :
                m  = Mean ( xmin , xmax , False )
                x0 = m    ( func , *args )  

        #
        ## check validity of x0
        #
        if not xmin <= x0 <= xmax :
            raise AttributeError ("Invalid x0 value %s<=%s<=%s" % ( xmin , x0 , xmax ) )

        #
        ## get the normalization
        #
        norm  = _integral_ ( func , xmin , xmax )
        if 0 >= norm :
            raise AttributeError ("Normalization integral is not positive %s" % norm )

        #
        ## Equation:  ifun(x) \equiv \int_{x0-x}^{x0+x}f(t)dt - N*prob = 0   
        #
        yval  = self.prob * norm
        def ifun ( x ) :
            if 0 >= x : return -yval 
            return _integral_ ( func , max ( xmin , x0 - x ) , min ( xmax , x0 + x )  ) - yval  
        
        from ostap.math.rootfinder import findroot
        s = findroot (  ifun , 0 , max ( xmax - x0 , x0 - xmin ) )

        from ostap.math.ve import VE 
        return VE ( x0 , s * s )

    @property
    def args ( self ) :
        "``args'' - other arguments for the function call"
        return self.__args
    @property
    def  xmin ( self ) :
        "``xmin''- low edge of the interval"
        return self.__xmin
    @property
    def  xmax ( self ) :
        "``xmax''- high edge of the interval"
        return self.__xmax
    @property
    def prob ( self ) :
        "``prop'' - confidence level"
        return self.__prob
    @property
    def x0 ( self ) :
        "``x0''- the center"
        return self.__x0

# =============================================================================
## @class CL_asymm
#  Calcualate asymmetic confidence interval (x1,x2) 
#  for the given function on (xmin,xmax) interval,
#  such as 
#  \f$ \begin{array}{l} f(x_1)=f(x_2)            \\
#  \int_{x_1}^{x_2}f(t)dr = p \int_{x_{min}}^{x_{max}}f(t)ft \end{array}\f$
#  Assumtions on function
#  - nonnegative
#  - unimodal
#  - positive integral between (xmin,xmax)
#  - zero outside interval     (xmin,xmax)
#  @code
#  fun   = lambda x : exp( -0.5*x*x)
#  reg   = CL_asymm ( 0.68 , -10 , 10 )
#  print reg ( fun )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-08-03
class CL_asymm(object) :
    """Calculate asymmetic confidence interval around x0
    for the given function on (xmin,xmax) interval
    function is assumed to be zero outside the interval
    >>> fun = lambda x : exp( -0.5*x*x)
    >>> reg = CL_asymm ( 0.68 , -10 , 10 )
    >>> print reg ( fun )
    """
    def __init__ ( self , prob ,  xmin , xmax , *args ) :
        
        if not 0.0 < prob < 1.0 :
            raise AttributeError ("Invalid value of prob/CL=%g" % prob)

        self.__prob = prob 
        self.__xmin = float ( xmin ) if isinstance ( xmin , integer_types ) else xmin 
        self.__xmax = float ( xmax ) if isinstance ( xmax , integer_types ) else xmax 
        self.__args = args

    def __str__ ( self ) :
        return "CL_asymm(%s,%s,%s)" % ( self.prob , self.xmin , self.xmax )

    ## solve equation f(x)=a 
    def _solve_  ( self , func , fval , xmn , xmx , *args ) :

        ifun = lambda x,*a : func(x,*a) - fval

        from ostap.math.rootfinder import findroot
        return findroot (  ifun , xmn , xmx , args = args )
                   
    def __call__ ( self , func , *args ) :

        ## additional arguments
        args   = args if args else self.args
        
        #
        ## define integration rules
        #
        if hasattr ( func , 'integral' ) :
            _integral_ = lambda f , low , high : f.integral (      low , high , *args )
        else                             :
            from ostap.math.integral import integral
            _integral_ = lambda f , low , high :   integral ( f  , low , high , *args )

        #
        ## xmin/max
        #
        xmin,xmax = self.xmin, self.xmax
        
        #
        # calculate mode
        if hasattr ( func , 'mode' ) : xmode = func.mode()
        else :
            md    = Mode ( xmin , xmax )
            xmode = md   ( func , *args )

        if not xmin <= xmode <= xmax :
            raise AttributeError ("Invalid mode value %s<=%s<=%s" % ( xmin , xmode , xmax ) )

        #
        ## get the normalization
        #
        norm  = _integral_ ( func , xmin , xmax )
        if 0 >= norm :
            raise AttributeError ("Normalization integral is not positive %s" % norm )

        normL = _integral_ ( func , xmin  , xmode )
        normR = _integral_ ( func , xmode , xmax  )

        from ostap.math.base import isequal, iszero
        
        ## solve equation f(x)=a 
        def _solve_  ( func , fval , xmn , xmx , *args ) :
            ##
            if isequal (  xmn , xmx   )  : return xmn
            ## 
            ifun = lambda x,*a : func(x,*a) - fval
            ## 
            fmn = ifun  ( xmn )
            if iszero  ( ifun ( xmn ) )  : return xmn
            fmx = ifun  ( xmx )
            if iszero  ( ifun ( xmx ) )  : return xmx 
            ##
            if 0 < fmx * fmn : ## more or less arbitrary choice 
                return xmx if abs ( fmx ) <= abs ( fmn ) else xmn 
            #
            from ostap.math.rootfinder import findroot 
            return findroot (  ifun , xmn , xmx , args = args )

        yval = self.prob * norm
        fm   = func ( xmode ) 
        def iifun ( f ) :

            if   iszero  ( f      ) : x1 , x2 = xmin,xmax
            elif isequal ( f , fm ) : return -yval 
            else : 
                x1 = _solve_ ( func ,  f , xmin  , xmode )
                x2 = _solve_ ( func ,  f , xmode , xmax  )

            return _integral_( func , x1 , x2 ) - yval 
            
        from ostap.math.rootfinder import findroot 
        l  = findroot (  iifun , 0 , func ( xmode ) )
        x1 = _solve_ ( func ,  l , xmin  , xmode )
        x2 = _solve_ ( func ,  l , xmode , xmax  )

        return x1 , x2 

    @property
    def args ( self ) :
        "``args'' - other arguments for the function call"
        return self.__args
    @property
    def  xmin ( self ) :
        "``xmin''- low edge of the interval"
        return self.__xmin
    @property
    def  xmax ( self ) :
        "``xmax''- high edge of the interval"
        return self.__xmax
    @property
    def prob ( self ) :
        "``prop'' - confidence level"
        return self.__prob

# =============================================================================
## calculate some statistical quantities of variable,
#  considering function to be PDF 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def sp_action ( func , actor , xmin = None , xmax = None , *args , **kwargs ) :
    """ Calculate some statistical quantities of variable, considering function to be PDF 
    """
    ##
    if   isinstance  ( xmin , num_types ) : xmn =  float ( xmin            ) 
    elif hasattr     ( func ,'GetXmin'  ) : xmn =  float ( func.GetXmin () )
    elif hasattr     ( func ,'xmin'     ) : xmn =  float ( func.xmin    () ) 
    else                                  : xmn =  neg_infinity 
    ##
    if   isinstance  ( xmax , num_types ) : xmx =  float ( xmax            )
    elif hasattr     ( func ,'GetXmax'  ) : xmx =  float ( func.GetXmax () ) 
    elif hasattr     ( func ,'xmax'     ) : xmx =  float ( func.xmax    () )
    else                                  : xmx =  pos_infinity 
    ##
    xmn = float ( xmn ) if isinstance ( xmn , integer_types ) else xmn 
    xmx = float ( xmx ) if isinstance ( xmx , integer_types ) else xmx
    #
    ## instantiate calculator and use it 
    calc = actor ( xmn , xmx )
    ##
    return calc  ( func , *args , **kwargs )

# =============================================================================
## get the N-moment of variable, considering function to be PDF 
# @code 
# >>> fun  = ...
# >>> m5   = moment( fun , 5 , xmin = 10 , xmax = 50 )
# @endcode
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date 2015-07-11
def moment ( func , N , xmin = None , xmax = None , err = False , x0 = 0 ) :
    """ Get the moment for the distribution
    >>> fun  = ...
    >>> mom5 = moment ( fun , 5 , xmin = 10 , xmax = 50 )
    """
    ## get the functions from ostap.stats.moments 
    actor = lambda x1,x2 : Moment ( N , x1 , x2 , err , x0 ) 
    return sp_action ( func , actor , xmin , xmax )

# =============================================================================
## get the N-th centralmoment of variable, considering function to be PDF 
# @code 
# >>> fun  = ...
# >>> m5   = central_moment( fun , 5 , xmin = 10 , xmax = 50 )
# @endcode
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date 2015-07-11
def central_moment ( func , N , xmin = None , xmax = None , err = False ) :
    """ Get the central moment for the distribution using 
    >>> fun  = ...
    >>> mom5 = central_moment ( fun , 5 , xmin = 10 , xmax = 50 )
    """
    ## get the functions from ostap.stats.moments
    if   0 == N : return 1
    elif 1 == N : return 0
    actor = lambda x1,x2 : CentralMoment ( N , x1 , x2 , err ) 
    return sp_action ( func , actor , xmin , xmax )

# =============================================================================
## get the N-th standartized moment of variable, considering function to be PDF 
# @code 
# >>> fun  = ...
# >>> m5   = std_moment( fun , 5 , xmin = 10 , xmax = 50 )
# @endcode
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date 2022-07-13
def std_moment ( func , N , xmin = None , xmax = None , err = False ) :
    """ Get the standartized moment for the distribution using 
    >>> fun  = ...
    >>> mom5 = std_moment ( fun , 5 , xmin = 10 , xmax = 50 )
    """
    ## get the functions from ostap.stats.moments
    if   0 == N : return 1
    elif 1 == N : return 0
    elif 2 == N : return 1
    actor = lambda x1,x2 : StdMoment ( N , x1 , x2 , err ) 
    return sp_action ( func , actor , xmin , xmax )

# =============================================================================
## get the mean value of variable, considering function to be PDF 
#  @code 
#  >>> fun = ...
#  >>> m   = mean( fun , xmin = 10 , xmax = 50 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def mean ( func , xmin = None , xmax = None , err = False ) :
    """ Get the mean-value for the distribution
    >>> fun = ...
    >>> m   = mean( fun , xmin = 10 , xmax = 50 )
    """
    ## get the functions from ostap.stats.moments 
    actor = lambda x1,x2 : Mean ( x1 , x2 , err ) 
    ## use it! 
    return sp_action ( func , actor , xmin , xmax )

# =============================================================================
## get the variance of the variable, considering function to be PDF 
#  @code 
#  >>> fun = ...
#  >>> v   = variance( fun , xmin = 10 , xmax = 50 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def variance ( func , xmin = None , xmax = None , err = False ) :
    """ Get the variance for the distribution using
    >>> fun = ...
    >>> v   = variance( fun , xmin = 10 , xmax = 50 )
    """
    ##
    ## get the functions from ostap.stats.moments 
    actor = lambda x1,x2 : Variance ( x1 , x2 , err ) 
    ## use it! 
    return sp_action ( func , actor  , xmin , xmax )

# =============================================================================
## get the rms of the variable, considering function to be PDF 
#  @code 
#  >>> fun = ...
#  >>> v   = rms( fun , xmin = 10 , xmax = 50 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def rms ( func , xmin = None , xmax = None , err = False ) :
    """ Get RMS for the distribution using
    >>> fun = ...
    >>> v   = rms( fun , xmin = 10 , xmax = 50 )
    """
    ##
    ## get the functions from ostap.stats.moments 
    actor = lambda x1,x2 : RMS ( x1 , x2 , err ) 
    ## use it! 
    return sp_action ( func , actor  , xmin , xmax )

# =============================================================================
## get the skewness of the variable, considering function to be PDF 
#  @code 
#  >>> fun  = ...
#  >>> skew = skewness ( fun , xmin = -10 , xmax = 10 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def skewness ( func , xmin = None , xmax = None , err = False ) :
    """ Get the skewness for the distribution using
    >>> fun = ...
    >>> v   = skewness ( fun , xmin = -10 , xmax = 10 )
    """
    ##
    ## get the functions from ostap.stats.moments 
    actor = lambda x1,x2 : Skewness ( x1 , x2 , err ) 
    ## use it! 
    return sp_action ( func , actor  , xmin , xmax )


# =============================================================================
## get the excessive kurtosis of the variable, considering function to be PDF
#  \f$ k = \frac{\mu_4}{\sigma^4}-3\f$ 
#  @code 
#  >>> fun  = ...
#  >>> kurt = kurtosis ( fun , xmin = -10 , xmax = 10 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def kurtosis ( func , xmin = None , xmax = None , err = False ) :
    """ Get the (exessive) kurtosis for the distribution
    >>> fun  = ...
    >>> kurt = kurtosis ( fun , xmin = 10 , xmax = 50 )
    """
    ##
    ## get the functions from ostap.stats.moments 
    actor = lambda x1,x2 : Kurtosis ( x1 , x2 , err ) 
    ## use it! 
    return sp_action ( func , actor  , xmin , xmax )

# =============================================================================
## get the median the variable, considering function to be PDF 
#  @code 
#  >>> fun = ...
#  >>> med = median ( fun , xmin = 10 , xmax = 50 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def median ( func , xmin = None , xmax = None ) :
    """ Get the median for the distribution using
    >>> fun = ...
    >>> v   = median( fun , xmin = 10 , xmax = 50 )
    """
    ## get the functions from ostap.stats.moments 
    actor = lambda x1,x2 : Median ( x1 , x2 ) 
    ##
    return sp_action ( func , actor , xmin , xmax )

# =============================================================================
## get the quantile of variable, considering function to be PDF 
#  @code 
#  >>> fun  = ...
#  >>> quan = quantile( fun , 0.1 , xmin = 10 , xmax = 50 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def quantile ( func , Q , xmin = None , xmax = None , err = False , x0 = 0 ) :
    """ Get quantile for the distribution
    >>> fun  = ...
    >>> quan = quantile ( fun , 0.1 , xmin = 10 , xmax = 50 )
    """
    ## get the functions from ostap.stats.moments 
    actor = lambda x1,x2 : Quantile ( Q , x1 , x2 ) 
    return sp_action ( func , actor , xmin , xmax )

# =============================================================================
## get the mode the variable, considering function to be PDF 
#  @code 
#  >>> fun = ...
#  >>> m   = mode( fun ,  xmin = 10 , xmax = 50 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def mode ( func , xmin = None , xmax = None ) :
    """ Get the mode for the distribution
    >>> fun = ...
    >>> v   = mode( fun ,  xmin = 10 , xmax = 50 )
    """
    ## get the functions from ostap.stats.moments 
    ## use it! 
    ## get the functions from ostap.stats.moments 
    actor = lambda x1,x2 : Mode  ( x1 , x2 ) 
    return sp_action ( func , actor , xmin , xmax )

# =============================================================================
## get the width, considering function to be PDF 
#  @code 
#  >>> fun   = ...
#  >>> x1,x2 = width( fun ,  xmin = 10 , xmax = 50 )
#  >>> fwhm  = x2-x1
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def width ( func , xmin = None , xmax = None , height_factor = 0.5 , mode = None ) :
    """ Get the width for the distribution
    >>> fun   = ...
    >>> x1,x2 = width ( fun ,  xmin = 10 , xmax = 50 )
    >>> fwhm  = x2-x1   
    """
    ## get the functions from ostap.stats.moments 
    ## use it! 
    ## get the functions from ostap.stats.moments
    actor = lambda x1,x2 : Width  ( x1 , x2 , height_factor ) 
    return sp_action ( func , actor , xmin , xmax , mode = mode )

# =============================================================================
## get the FWHM, considering function to be PDF 
#  @code 
#  >>> fun   = ...
#  >>> width = fwhm ( fun , xmin = 10 , xmax = 50 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def fwhm ( func , xmin = None , xmax = None , mode = None ) :
    """ Get the width for the distribution
    >>> fun   = ...
    >>> x1,x2 = width ( fun ,  xmin = 10 , xmax = 50 )
    >>> fwhm  = x2-x1   
    """
    x1 , x2 = width ( func , xmin = xmin , xmax = xmax , height_factor = 0.5 , mode = mode )
    return  x2 - x1

# =============================================================================
## get the symmetric confidence interval around x0 for (xmin,xmax) interval 
#  @code 
#  fun  = lambda x : exp( - 0.5 * x * x ) 
#  x_1  = cl_symm ( fun , 0.68 , -10 , 10 )
#  print x_1 
#  x_2  = cl_symm ( fun , 0.68 ,   0 , 10 , x0 = 0 )
#  print x_2 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-08-03
def cl_symm ( func , prob , xmin = None , xmax = None , x0 = None ) :
    """ Symmetric confidence interval around x0    
    >>> fun  = lambda x : exp( - 0.5 * x * x )
    >>> x_1  = cl_symm ( fun , 0.68 , -10 , 10 )
    >>> print x_1 
    >>> x_2  = cl_symm ( fun , 0.68 ,   0 , 10 , x0 = 0 )
    >>> print x_2 
    """
    ## get the functions
    actor = lambda x1,x2 : CL_symm ( prob , x1 , x2 , x0 = x0 ) 
    ## and use it! 
    return sp_action ( func , actor , xmin , xmax )

# =============================================================================
## get the asymmetric confidence interval around x0 for (xmin,xmax) interval 
#  @code 
#  fun  = lambda x : exp( - 0.5 * x * x ) 
#  x_1,x_2  = cl_asymm ( fun , 0.68 , -10 , 10 )
#  print x_1 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-08-03
def cl_asymm ( func , prob , xmin = None , xmax = None ) :
    """ Asymmetric confidence interval around x0    
    >>> fun  = lambda x : exp( - 0.5 * x * x )
    >>> x_1,x_2  = cl_asymm ( fun , 0.68 , -10 , 10 )
    >>> print x_1,x_2 
    """
    ## get the functions
    actor = lambda x1,x2 : CL_asymm ( prob , x1 , x2 ) 
    ## and use it! 
    return sp_action ( func , actor , xmin , xmax )


                                
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
