#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  moments.py 
#  Utilities to get moments for various functions/distributions/pdfs
#  - Moment
#  - CentralMoment
#  - Mean
#  - Variance 
#  - RMS 
#  - Skewness
#  - Kurtosis
#  For these quantities scipy.integrate is used integration engine, in case scipy
#  is not available, a hand-made replacement is used
#
#  With help of scipy.optimize.brent additional quantities can be calculated
#  - Median
#  - Quantile 
#  - Mode
#  - Width
#  - symmetric and asymmetric "confidence intervals"
#
#  All objects exists as classes/functors and as standalone simple functions
#  - moment
#  - centralMoment
#  - mean
#  - variance 
#  - rms 
#  - skewness
#  - kurtosis
#  - median
#  - quantile 
#  - mode
#  - width
#  - cl_symm and sl_asymm 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-06  
# =============================================================================
"""Utilities to get moments for various functions/distributions/pdfs
- Moment
- CentralMoment
- Mean
- Variance 
- RMS 
- Skewness
- Kurtosis
For these quantities scipy.integrate is used integration engine, in case scipy
is not available, a hand-made replacement is used

Also it calculates with help of scipy.optimize.brent following quantities
- Median
- Quantile 
- Mode
- Width
- symmetric and asymmetric ``confidence intervals''

All objects exists as classes/functors and as standalone simlpe functions
- moment
- centralMoment
- mean
- variance 
- rms 
- skewness
- kurtosis
- median
- quantile 
- mode
- width
- cl_symm and sl_asymm 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-06-06"
__all__     = (
    ##
    ## stat-quantities 
    "Moment"        , ## calculate N-th moment of functions/distribitions, etc 
    "CentralMoment" , ## calculate N-th central moment of functions/distribitions
    "Mean"          , ## calculate "mean"     for functions/distribitions, etc 
    "Variance"      , ## calculate "variance" for functions/distribitions, etc 
    "RMS"           , ## calculate "RMS"      for functions/distribitions, etc 
    "Skewness"      , ## calculate "skewness" for functions/distribitions, etc 
    "Kurtosis"      , ## calculate "kurtosis" for functions/distribitions, etc 
    "Median"        , ## calculate "median"   for functions/distribitions, etc (brentq)
    "Quantile"      , ## calculate "quantile" for functions/distribitions, etc (brentq)
    "Mode"          , ## calculate "mode"     for functions/distribitions, etc (brentq)
    "Width"         , ## calculate "width"    for functions/distribitions, etc (brentq)
    "CL_symm"       , ## calcualte symmetrical confidence intervals            (brentq)
    "CL_asymm"      , ## calcualte asymmetrical confidence intervals           (brentq)
    ##
    ## stat-quantities   
    "moment"        , ## calculate N-th moment of functions/distribitions, etc 
    "central_moment", ## calculate N-th moment of functions/distribitions, etc 
    "mean"          , ## calculate "mean"     for functions/distribitions, etc 
    "variance"      , ## calculate "variance" for functions/distribitions, etc 
    "rms"           , ## calculate "RMS"      for functions/distribitions, etc 
    "skewness"      , ## calculate "skeness"  for functions/distribitions, etc 
    "kurtosis"      , ## calculate "kurtosis" for functions/distribitions, etc 
    "median"        , ## calculate "median"   for functions/distribitions, etc (brentq)
    "quantile"      , ## calculate "quantile" for functions/distribitions, etc (brentq)
    "mode"          , ## calculate "mode"     for functions/distribitions, etc (brentq)
    "width"         , ## calculate "width"    for functions/distribitions, etc (brentq
    "cl_symm"       , ## calculate symmetrical confidence intervals            (brentq)
    "cl_asymm"      , ## calculate asymmetrical confidence intervals           (brentq)
    ##
    ) 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.moments' )
else                       : logger = getLogger ( __name__              )
# =============================================================================

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
class Moment(object) :
    """Calculate the N-th moment for the distribution
    
    >>> xmin,xmax = 0,math.pi 
    >>> mean  = Moment(1,xmin,xmax)  ## specify min/max
    >>> value = mean  ( math.sin )
    """
    ## constructor
    def __init__ ( self , N , xmin , xmax , err = False , x0 = 0 , *args ) :
        """Contructor 
        """
        if not isinstance ( N , ( int , long ) ) or 0 > N  :
            raise TypeError('Moment: illegal order')
        
        self._N    = N 
        self._xmin = float ( xmin ) if isinstance ( xmin , ( int , long ) ) else xmin 
        self._xmax = float ( xmax ) if isinstance ( xmax , ( int , long ) ) else xmax 
        self._x0   = x0  
        self._err  = err
        self._args = args
        self._moms = {} 

    ## make an integral 
    def _integral_ ( self , func , xmn , xmx , *args ) :
        from ostap.math.integral import Integral
        integrator = Integral ( func , xmn , err = self._err )
        return integrator ( xmx , *args )
    
    ## calculate un-normalized 0-moment  
    def _moment0_ ( self , func , *args ) :
        """Calculate un-normalized 0-moment
        """
        return self._integral_ ( func , self._xmin , self._xmax , *args )
    
    ## calculate un-normalized k-moment  
    def _momentK_ ( self , k , func , mu = None , *args ) :
        """Calculate unnormalized k-moment
        """
        x0     = self._x0 if mu is None else mu 
        func_N = lambda x,*a : func( x , *a ) * ( ( x - x0 ) ** k  )
        return self._integral_ ( func_N , self._xmin , self._xmax , *args )
    
    ## calculate the moment 
    def __call__ ( self , func , *args ) :
        ## 
        args   = args if args else self._args
        ##
        n0  = self._moment0_ (            func ,            *args ) 
        nN  = self._momentK_ ( self._N  , func , self._x0 , *args ) 
        ##
        return nN/n0

    ## print it!
    def __str__ ( self ) :
        return "Moment(%d,%s,%s,%s,%s)" % ( self._N    ,
                                            self._xmin ,
                                            self._xmax ,
                                            self._err  ,
                                            self._x0   )                                            
                                            

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
class CentralMoment(Moment) :
    """Calculate the N-th central moment for the distribution
    
    >>> xmin,xmax = 0,math.pi 
    >>> mc        = CentralMoment(1,xmin,xmax)  ## specify min/max
    >>> value     = mc  ( math.sin )    
    """
    ## constructor
    def __init__ ( self , N , xmin , xmax , err = False, *args ) :
        """Contructor 
        """
        Moment.__init__ ( self , N , xmin , xmax , err , 0.0 , *args ) 

    ## calculate the central moment
    def __call__ ( self , func , *args ) :
        ## 
        args   = args if args else self._args
        ##
        n0  = self._moment0_ (     func ,             *args ) 
        n1  = self._momentK_ ( 1 , func , mu = 0.0  , *args )
        ## get mean
        mu  = float(n1/n0)
        ## use it 
        nN  = self._momentK_ ( self._N  , func , mu , *args ) 
        ##
        return nN/n0

    ## print it!
    def __str__ ( self ) :
        return "CentralMoment(%d,%s,%s,%s)" % ( self._N    ,
                                                self._xmin ,
                                                self._xmax ,
                                                self._err  )
                                                
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
    def __init__ ( self , xmin , xmax , err = False ) :
        Moment.__init__ ( self , 1 , xmin , xmax , err )

    def __str__ ( self ) :
        return "Mean(%s,%s,%s)" % ( self._xmin ,
                                    self._xmax ,
                                    self._err  )                                            

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
class Variance(Mean) :
    """Calculate the variance for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> variance  = Variance ( xmin,xmax )  ## specify min/max
    >>> value     = variance ( math.sin  )
    """
    def __init__ ( self , xmin , xmax , err = False ) :
        Mean.__init__ ( self , xmin , xmax , err )
    ## calculate the variance 
    def __call__ ( self , func , *args ) :
        ## 
        args   = args if args else self._args
        ##
        n0 = self._moment0_ (     func ,            *args ) ## moment-0
        n1 = self._momentK_ ( 1 , func , mu = 0.0 , *args ) ## moment-1 
        ##
        mu = float(n1/n0)                        ## mean-value 
        ## central moment 
        m2 = self._momentK_ ( 2 , func , mu , *args ) 
        ##
        return m2/n0
    
    def __str__ ( self ) :
        return "Variance(%s,%s,%s)" % ( self._xmin ,
                                        self._xmax ,
                                        self._err  )                                            

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
    def __init__ ( self , xmin , xmax , err = False ) :
        Variance.__init__ ( self , xmin , xmax , err )
        
    ## calculate the variance 
    def __call__ ( self , func , *args ) :
        """Calculate the RMS for the distribution or function          
        """
        ##
        args   = args if args else self._args
        ##
        var2 = Variance.__call__ ( self , func , *args )
        import ostap.math.math_ve as ME 
        return ME.sqrt ( var2 ) 

    def __str__ ( self ) :
        return "RMS(%s,%s,%s)" % ( self._xmin ,
                                   self._xmax ,
                                   self._err  )                                            

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
class Skewness(Variance) :
    """
    Calculate the variance for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> skew      = Skewness ( xmin,xmax )  ## specify min/max
    >>> value     = skew     ( math.sin  )
    """
    def __init__ ( self , xmin , xmax , err = False ) :
        Variance.__init__ ( self , xmin , xmax , err )
    ## calculate the variance 
    def __call__ ( self , func , *args ) :
        ## 
        args   = args if args else self._args
        ##
        n0 = self._moment0_ (     func ,          *args ) ## norm
        n1 = self._momentK_ ( 1 , func , mu = 0 , *args ) ## m1
        ## get mean-value 
        mu = float(n1/n0) ## mean-value
        ## 
        m2 = self._momentK_ ( 2 , func , mu     , *args ) ## mu2 
        m3 = self._momentK_ ( 3 , func , mu     , *args ) ## mu3 
        ##
        m2 /= n0 ## normalize 
        m3 /= n0 ## normalize
        ## 
        return m3/(m2**(3.0/2))
    
    def __str__ ( self ) :
        return "Skewness(%s,%s,%s)" % ( self._xmin ,
                                        self._xmax ,
                                        self._err  )                                            

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
class Kurtosis(Skewness) :
    """Calculate the variance for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> kurt      = Kurtosis ( xmin,xmax )  ## specify min/max
    >>> value     = kurt     ( math.sin  )
    """
    def __init__ ( self , xmin , xmax , err = False ) :
        Skewness.__init__ ( self , xmin , xmax , err )
    ## calculate the variance 
    def __call__ ( self , func , *args ) :
        ## 
        args   = args if args else self._args
        ##
        n0 = self._moment0_ (     func ,          *args ) ## norm
        n1 = self._momentK_ ( 1 , func , mu = 0 , *args ) ## m1
        ## get mean-value 
        mu = float(n1/n0) ## mean-value
        ## 
        m2 = self._momentK_ ( 2 , func , mu     , *args ) ## mu2 
        m4 = self._momentK_ ( 4 , func , mu     , *args ) ## mu3 
        ##
        m2 /= n0 ## normalize 
        m4 /= n0 ## normalize
        ## 
        return m4/(m2*m2)-3.0 
    
    def __str__ ( self ) :
        return "Kurtosis(%s,%s,%s)" % ( self._xmin ,
                                        self._xmax ,
                                        self._err  )                                            


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
#  @attention <code>scipy.optimize.brentq</code> is needed!
#  @see scipy.optimize.brentq 
#  @date   2015-07-12
class Median(RMS) :
    """Calculate median for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> median    = Median ( xmin,xmax )  ## specify min/max
    >>> value     = median ( math.sin  )
    - scipy.optimize.brentq is used 
    """
    def __init__ ( self , xmin , xmax ) :
        RMS.__init__ ( self , xmin , xmax , err = False )

    ## calculate he median
    def _median_ ( self , func , xmin , xmax , *args ) :
        ## need to know the integral
        from ostap.math.integral import IntegralCache
        iint   = IntegralCache ( func ,  xmin , False , *args )
        half   = 2.0 / iint    ( xmax ) 
        
        from scipy import optimize
        ifun   = lambda x : iint( x ) * half - 1.0

        ## @see https://en.wikipedia.org/wiki/Median#Inequality_relating_means_and_medians
        try: 
            meanv = Mean . __call__ ( self , func , *args )
            sigma = RMS  . __call__ ( self , func , *args )
            import math
            xmn   = meanv - 2 * sigma ## use 2 instead of 1 
            xmx   = meanv + 2 * sigma ## use 2 instead of 1
            #
            if isinstance ( xmin , float ) : xmn = max ( xmn , xmin ) 
            if isinstance ( xmax , float ) : xmx = min ( xmx , xmax )
            #
            result = optimize.brentq ( ifun , xmn , xmx )
        except :
            result = optimize.brentq ( ifun , xmin , xmax )
            
        return result

        
    ## calculate the median 
    def __call__ ( self , func , *args ) :
        return self._median_ ( func , self._xmin , self._xmax )

    def __str__ ( self ) :
        return "Median(%s,%s)" % ( self._xmin , self._xmax )
    
# =============================================================================
## get the quantile
#  Calculate the quantile for the distribution or function  
#  @code
#  xmin,xmax = 0,math.pi 
#  quantile  = Quantile ( 0.1 , xmin,xmax )  ## specify min/max
#  value     = quantile ( math.sin  )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @attention scipy.optmize.brentq is used 
#  @date   2015-07-12
class Quantile(Median) :
    """Calculate quantiles for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> quantile  = Quantile ( 0.1 , xmin,xmax )  ## specify min/max
    >>> value     = quantile ( math.sin  )
    - scipy.optmize.brentq is used 
    """
    def __init__ ( self , Q , xmin , xmax ) :
        Median.__init__ ( self , xmin , xmax )
        #
        if Q < 0 : raise ArrtibuteError ( 'Quantile is invalid %s' % Q )
        if Q > 1 : raise ArrtibuteError ( 'Quantile is invalid %s' % Q )
        self._Q = float( Q ) 
        
    def __str__ ( self ) :
        return "Quantile(%s,%s,%s)" % ( self._Q , self._xmin , self._xmax )

    ## calculate the median 
    def __call__ ( self , func , *args ) :
        ##

        if    0.5 == self._Q : return Median.__call__ ( self , func , *args ) 
        elif  0.0 == self._Q : return self._xmin
        elif  1.0 == self._Q : return self._xmax

        ## need to know the integral
        from ostap.math.integral import IntegralCache
        iint = IntegralCache ( func, self._xmin, False , *args )
        quan = 1.0 / iint    (  self._xmax ) / self._Q 
        
        from scipy import optimize
        ifun   = lambda x : iint( x ) * quan - 1.0

        xmn = self._xmin
        xmx = self._xmax

        p   = 0.5
        l   = 0.5

        ## make some bracketing before next step 
        while ( not isinstance ( xmn , float ) ) or ( not isinstance ( xmx , float ) ) or l>0.1 :   
        
            l /= 2            
            m = self._median_ ( func , xmn , xmx , *args ) 
            
            if   self._Q < p :
                xmn   = xmn 
                xmx   = float( m ) 
                p    -= l 
            elif self._Q > p :
                xmn   = float ( m ) 
                xmx   = xmx
                p    +=  l  
            else : return m               ## RETURN 

        ## finally, calculate quantile 
        result = optimize.brentq ( ifun , xmn , xmx )
            
        return result


# =============================================================================
## @class Mode
#  Calculate the mode for the distribution or function  
#  @code
#  xmin,xmax = 0,math.pi 
#  mode      = Mode ( xmin,xmax )  ## specify min/max
#  value     = mode ( math.sin  )
#  @endcode 
#  @attention scipy.optimize.brentq is used 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-12
class Mode(Median) :
    """Calculate the mode for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> mode      = Mode ( xmin,xmax )  ## specify min/max
    >>> value     = mode ( math.sin  )
    - scipy.optimize.brentq is used 
    """
    def __init__ ( self , xmin , xmax ) :
        Median.__init__ ( self , xmin , xmax )
        
    ## calculate the mode 
    def __call__ ( self , func , *args ) :
        ##
        
        ## use mean    as intial approximation for mode 
        m1     = Mean   .__call__ ( self , func , *args )
        
        ## use median as intial approximation for mode 
        ## m2     = Median.__call__ ( self , func , *args )
        
        ## use the point intermediate between mean and median as approximation 
        ## m0     = 0.5 * ( m1 + m2 )

        m0 = m1 
        ifun = lambda x,*a : -1.0 * float( func ( x , *a ) )
        
        from scipy import optimize
        result = optimize.minimize (
            ifun                                   , 
            x0     = float ( m0 )                  ,
            bounds = [ (self._xmin , self._xmax) ] ,
            args   = args )
        
        return result.x[0]
    
    def __str__ ( self ) :
        return "Mode(%s,%s)" % ( self._xmin , self._xmax )

# =============================================================================
## @class Width
#  Calculate the full width at half heigh for the distribution or function  
#  @code
#  xmin,xmax = 0,math.pi 
#  width     = Width ( xmin,xmax )  ## specify min/max
#  x1,x2     = width ( math.sin  )
#  fwhm      = x2-x1
#  @endcode 
#  @attention scipy.optimize.brentq is used 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-12
class Width(Mode) :
    """Calculate the mode for the distribution or function  
    >>> xmin,xmax = 0,math.pi 
    >>> width     = Width ( xmin,xmax )  ## specify min/max
    >>> x1,x2     = width ( math.sin )
    >>> fwhm      = x2-x1
    - scipy.optimize.brentq is used 
    """
    def __init__ ( self , xmin , xmax , height_factor = 0.5 ) :
        Mode.__init__ ( self , xmin , xmax )
        self._hfactor = height_factor
        
    ## calculate the width
    def __call__ ( self , func , *args ) :
        ##

        ## get the position of the mode
        m0  = Mode.__call__ ( self , func , *args )

        ## function  value at the maximum
        v0      = func ( m0 , *args )

        ## half height 
        vheight = 1.0 * v0 * self._hfactor

        
        ## use scipy to find solution 
        from scipy import optimize        
        ifun = lambda x,*a : float(func (x,*a))-vheight
        x1 = optimize.brentq ( ifun , self._xmin , m0         , args = args )
        x2 = optimize.brentq ( ifun , m0         , self._xmax , args = args ) 
        
        return x1,x2

    def __str__ ( self ) :
        return "Width(%s,%s,%s)" % ( self._xmin , self._xmax , self._hfactor)
    

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
#  @attention scipy.optimize.brentq is used 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-08-03
class CL_symm(object) :
    """Calculate symmetic confidence interval around x0
    for the given function on (xmin,xmax) interval
    function is assumed to be zero outside the interval
    >>> fun = lambda x : exp( -0.5*x*x)
    >>> reg = CL_symm ( 0.68 , -10 , 10 )
    >>> print reg ( fun )
    - scipy.optimize.brentq is used 
    """
    def __init__ ( self , prob ,  xmin , xmax , x0 = None , *args ) :
        
        if not 0.0 < prob < 1.0 :
            raise AttributeError ("Invalid value of prob/CL=%g" % prob)

        self._prob = prob 
        self._xmin = float ( xmin ) if isinstance ( xmin , ( int , long ) ) else xmin 
        self._xmax = float ( xmax ) if isinstance ( xmax , ( int , long ) ) else xmax 
        self._x0   = float ( x0   ) if isinstance ( x0   , ( int , long ) ) else x0  
        self._args = args
        
    def __str__ ( self ) :
        return "CL_sym(%s,%s,%s,%s)" % ( self._prob ,
                                         self._xmin , self._xmax , self._x0   )

    def __call__ ( self , func , *args ) :

        ## additional arguments
        args   = args if args else self._args
        
        #
        ## define integration rules
        #
        if hasattr ( func , 'integral' ) :
            _integral_ = lambda f , low , high : f.integral (      low , high , *args )
        else                             :
            from ostap.math.integral import integral 
            _integral_ = lambda f , low , high :   integral ( f  , low , high , False , *args )

        #
        ## xmin/max
        #
        xmin,xmax = self._xmin, self._xmax
        
        #
        ## calculate x0 as "mean"-value
        #
        x0 = self._x0 
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
        yval  = self._prob * norm
        def ifun ( x ) :
            if 0 >= x : return -yval 
            return _integral_ ( func , max ( xmin , x0 - x ) , min ( xmax , x0 + x )  ) - yval  
        
        ## use scipy to find solution 
        from scipy import optimize        
        s = optimize.brentq (  ifun , 0 , max ( xmax - x0 , x0 - xmin ) )

        from ostap.math.ve import VE 
        return VE ( x0 , s * s )


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
#  @attention scipy.optimize.brentq is used 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-08-03
class CL_asymm(object) :
    """Calculate asymmetic confidence interval around x0
    for the given function on (xmin,xmax) interval
    function is assumed to be zero outside the interval
    >>> fun = lambda x : exp( -0.5*x*x)
    >>> reg = CL_asymm ( 0.68 , -10 , 10 )
    >>> print reg ( fun )
    - scipy.optimize.brentq is used 
    """
    def __init__ ( self , prob ,  xmin , xmax , *args ) :
        
        if not 0.0 < prob < 1.0 :
            raise AttributeError ("Invalid value of prob/CL=%g" % prob)

        self._prob = prob 
        self._xmin = float ( xmin ) if isinstance ( xmin , ( int , long ) ) else xmin 
        self._xmax = float ( xmax ) if isinstance ( xmax , ( int , long ) ) else xmax 
        self._args = args

    def __str__ ( self ) :
        return "CL_asymm(%s,%s,%s)" % ( self._prob , self._xmin , self._xmax )

    ## solve equation f(x)=a 
    def _solve_  ( self , func , fval , xmn , xmx , *args ) :

        ifun = lambda x,*a : func(x,*a) - fval

        ## use scipy to find solution 
        from scipy import optimize        
        return optimize.brentq (  ifun , xmn , xmx , args = args )

                   
    def __call__ ( self , func , *args ) :

        ## additional arguments
        args   = args if args else self._args
        
        #
        ## define integration rules
        #
        if hasattr ( func , 'integral' ) :
            _integral_ = lambda f , low , high : f.integral (      low , high , *args )
        else                             :
            from ostap.math.integral import integral
            _integral_ = lambda f , low , high :   integral ( f  , low , high , False , *args )

        #
        ## xmin/max
        #
        xmin,xmax = self._xmin, self._xmax
        
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
            ## use scipy to find solution 
            from scipy import optimize
            return optimize.brentq (  ifun , xmn , xmx , args = args )

        yval = self._prob * norm
        fm   = func ( xmode ) 
        def iifun ( f ) :

            if   iszero  ( f      ) : x1 , x2 = xmin,xmax
            elif isequal ( f , fm ) : return -yval 
            else : 
                x1 = _solve_ ( func ,  f , xmin  , xmode )
                x2 = _solve_ ( func ,  f , xmode , xmax  )

            return _integral_( func , x1 , x2 ) - yval 
            
        from scipy import optimize
        l  = optimize.brentq (  iifun , 0 , func ( xmode ) )
        x1 = _solve_ ( func ,  l , xmin  , xmode )
        x2 = _solve_ ( func ,  l , xmode , xmax  )

        return x1 , x2 

 
# =============================================================================
## calculate some statistical quantities of variable,
#  considering function to be PDF 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def sp_action ( func , actor , xmin = None , xmax = None ) :
    """Calculate some statistical quantities of variable, considering function to be PDF 
    """
    ##
    import numpy
    ##
    if   isinstance  ( xmin , (float,int,long) ) : xmn =  float ( xmin            ) 
    elif hasattr     ( func ,'GetXmin' )         : xmn =  float ( func.GetXmin () )
    elif hasattr     ( func ,'xmin'    )         : xmn =  float ( func.xmin    () ) 
    else                                         : xmn = -numpy.inf
    ##
    if   isinstance  ( xmax , (float,int,long) ) : xmx =  float ( xmax            )
    elif hasattr     ( func ,'GetXmax' )         : xmx =  float ( func.GetXmax () ) 
    elif hasattr     ( func ,'xmax'    )         : xmx =  float ( func.xmax    () )
    else                                         : xmx = +numpy.inf
    ##
    xmn = float ( xmn ) if isinstance ( xmn , ( int , long ) ) else xmn 
    xmx = float ( xmx ) if isinstance ( xmx , ( int , long ) ) else xmx
    #
    ## instantiate calculator and use it 
    calc = actor ( xmn , xmx )
    ##
    return calc  ( func )

# =============================================================================
## get the N-moment of variable, considering function to be PDF 
# @code 
# >>> fun  = ...
# >>> m5   = moment( fun , 5 , xmin = 10 , xmax = 50 )
# @endcode
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date 2015-07-11
def moment ( func , N , xmin = None , xmax = None , err = False , x0 = 0 ) :
    """ Get the moment for the distribution using scipy/numpy
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
    """Get the central moment for the distribution using scipy/numpy
    >>> fun  = ...
    >>> mom5 = central_moment ( fun , 5 , xmin = 10 , xmax = 50 )
    """
    ## get the functions from ostap.stats.moments 
    actor = lambda x1,x2 : CentralMoment ( N , x1 , x2 , err ) 
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
    """Get the mean-value for the distribution using scipy/numpy
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
    """Get the variance for the distribution using scipy/numpy
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
    """Get RMS for the distribution using scipy/numpy
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
    """Get the skewness for the distribution using scipy/numpy
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
    """Get the (exessive) kurtosis for the distribution using scipy/numpy
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
    """Get the median for the distribution using scipy/numpy
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
    """Get quantile for the distribution using scipy/numpy
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
    """Get the mode for the distribution using scipy/numpy
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
def width ( func , xmin = None , xmax = None , height_factor = 0.5 ) :
    """ Get the width for the distribution using scipy/numpy
    >>> fun   = ...
    >>> x1,x2 = width ( fun ,  xmin = 10 , xmax = 50 )
    >>> fwhm  = x2-x1   
    """
    ## get the functions from ostap.stats.moments 
    ## use it! 
    ## get the functions from ostap.stats.moments
    actor = lambda x1,x2 : Width  ( x1 , x2 , height_factor ) 
    return sp_action ( func , actor , xmin , xmax )

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
    """Asymmetric confidence interval around x0    
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
# The END 
# =============================================================================
