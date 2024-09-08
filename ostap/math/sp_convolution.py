#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/sp_convolution.py
#  Module with some useful NumPy/SciPy-based utilities for convolution 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2020-02-28
# =============================================================================
"""NumPy/SciPy-based utilities for convolution 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2020-02-28"
__all__     = ()
# =============================================================================
import warnings 
from   ostap.math.operations import Function
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.sp_convolution' )
else                       : logger = getLogger ( __name__                    )
# =============================================================================
try : # =======================================================================
    # =========================================================================
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")        
        import scipy.signal as SS
    import numpy        as np
    # =========================================================================
    # =========================================================================
    ## simple class for scipy-based interpolation  
    class ArrayConvolution(object) :
        """ Simple class for SciPy-based convolution 
        >>> x = ... ##
        >>> y = ... ##
        >>> sp = ArrayConvolution ( x , y , ... )
        """
        # =====================================================================
        ## Create the convolution 
        #  @code
        #  >>> x  = ... ##
        #  >>> y  = ... ##
        #  >>> sp = ArrayConvolution (x , y )
        #  @endcode
        def __init__ ( self , x , y , mode = 'full' ) :
            """ Create the convolution 
            >>> x  = ... ##
            >>> y  = ... ##
            >>> sp = ArrayConvolution (x , y )
            """            
            self.__result = None 
            self.__result = SS.fftconvolve ( x , y , mode ) 
            
        @property
        def result( self ) :
            """`result' : result with array convolution"""
            return self.__result 
        
    # =========================================================================
    ## simple class for SciPy/FFT-based convolution 
    class GaussConvolution(ArrayConvolution,Function) :
        """ Simple class for SciPy/FFT-based convolution 
        >>> func  = ...
        >>> sigma = ...
        >>> sp    = GaussConvolution ( func = func , xmin = 0 , xmax = 1 , sigma = 1 , N = 2048 ) 
        """
        # =================================================================
        ## Create the convolution 
        #  @code
        #  >>> x  = ... ##
        #  >>> y  = ... ##
        #  >>> sp = ArrayConvolution ( x , y )
        #  @endcode
        def __init__ ( self           ,
                       func           ,
                       xmin           ,
                       xmax           ,
                       N     = 1024   ,
                       sigma = 1      ,
                       mode  = 'same' ) :
            """ Create the convolution 
            >>> x  = ... ##
            >>> y  = ... ##
            >>> sp = GaussConvolution ( x , y )
            """            
            from ostap.core.ostap_types import integer_types
            assert isinstance ( N , integer_types ) and 1 < N , 'Illegal N %s!' % N 
            
            self.__func   = func 
            self.__params = xmin , xmax , N , abs ( sigma ) 
            
            from ostap.math.operations import digitize
            dfunc  = digitize ( func , xmin , xmax , N )

            from ostap.math.math_ve import gauss_pdf         
            gpdf   = lambda x : gauss_pdf ( x , sigma = self.sigma )
            
            vgauss = np.vectorize ( gpdf )            
            dx = float ( xmax - xmin ) / N 
            a  = float ( xmax - xmin ) / 2 
            
            x  = np.arange ( -a , a , dx ) 
            dgauss = vgauss ( x )  
            
            ## make the real convolution
            super(GaussConvolution,self).__init__ ( dfunc , dgauss , mode )
            
            r  = self.result
            r *= dx   
            
            self.__spline = None

        # =====================================================================
        def xmin   ( self ) : return self.__params [ 0 ]
        def xmax   ( self ) : return self.__params [ 1 ]

        @property
        def N      ( self ) : return self.__params [ 2 ]
        @property
        def sigma  ( self ) : return self.__params [ 3 ]
        
        @property
        def params ( self ) : return self.__params
        
        @property
        def spline ( self ) :
            """ `spline' : get the spline function for the result of convolution"""
            if self.__spline  : return self.__spline
            from ostap.math.sp_interpolation import SplineInterpolator 
            x = np.linspace ( self.xmin() , self.xmax() , self.N )
            self.__spline = SplineInterpolator ( ( x , self.result ) , 3 )
            return self.__spline

        # =====================================================================
        def __call__ ( self , x  ) :
            if not self.__spline : self.spline
            return self.__spline ( x ) 
        
        def __str__ ( self ) :
            return "(%s{*}Gauss(%s))" % ( self.__func , self.sigma ) 
        
        __repr__ = __str__

        # =====================================================================
        ## simple class for SciPy/FFT-based convolution 
        class Convolution ( ArrayConvolution,Function) :
            """ Simple class for SciPy/FFT-based convolution 
            >>> func  = ...
            >>> reso  = ...
            >>> sp    = Convolution ( func1 = func , xmin1 = 0 , xmax1 = 1  , 
            >>>                       func2 = reso , xmin2 = 0 , xmax2 = 10 , N = 2048 ) 
            """
            # =================================================================
            ## Create the convolution 
            #  @code
            #  >>> x  = ... ##
            #  >>> y  = ... ##
            #  >>> sp = Convolution (x , y )
            #  @endcode
            def __init__ ( self           ,
                           func1          ,
                           xmin1          ,
                           xmax1          ,
                           func2          ,
                           xmin2          ,
                           xmax2          ,
                           N     = 1024   ,
                           mode  = 'same' ) :
                """ Create the convolution 
                >>> x  = ... ##
                >>> y  = ... ##
                >>> sp = Convolution (x , y )
                """                
                from ostap.core.ostap_types import integer_types
                assert isinstance ( N , integer_types ) and 1 < N , 'Illegal N %s!' % N 
                
                self.__func1   = func1
                self.__params1 = xmin1 , xmax1 , N
                
                self.__func2   = func2 
                self.__params2 = xmin2 , xmax2 
                
                from ostap.math.operations import digitize
                dfunc1 = digitize ( func1 , xmin1 , xmax1 , N )
                
                dx     = float ( xmax1 - xmin1 ) / N 
                
                dfunc2 = digitize ( func2 , xmin12 , xmax2 , dx )
                
                ## make the real convolution
                super(Convolution,self).__init__ ( dfunc1 , dfunc2 , mode )
                
                r  = self.result
                r *= dx 

                self.__spline = None
                
            def xmin    ( self ) : return self.__params1 [ 0 ]
            def xmax    ( self ) : return self.__params1 [ 1 ]
            def xmin1   ( self ) : return self.__params1 [ 0 ]
            def xmax1   ( self ) : return self.__params1 [ 1 ]
            def xmin2   ( self ) : return self.__params2 [ 0 ]
            def xmax2   ( self ) : return self.__params2 [ 1 ]
            
            @property
            def N       ( self ) : return self.__params1 [ 2 ]
            
            @property
            def params1 ( self ) : return self.__params1
            @property
            def params2 ( self ) : return self.__params2
            
            @property
            def spline ( self ) :
                """`spline' : get the spline function for the result of convolution"""
                if self.__spline  : return self.__spline
                from ostap.math.sp_interpolation import SplineInterpolator 
                x = np.linspace ( self.xmin() , self.xmax() , self.N )
                self.__spline = SplineInterpolator ( ( x , self.result ) , 3 )
                return self.__spline
            
            def __call__ ( self , x  ) :
                if not self.__spline : self.spline
                return self.__spline ( x )
            
            def __str__ ( self ) :
                return "(%s{*}%s)" % ( self.__func1 , self.__func2 ) 
            
            __repr__ = __str__

    __all__     = (
        'ArrayConvolution' , ## SciPy' FFT-based convolution for two arrays 
        'GaussConvolution' , ## SciPy' FFT-based convolution for function with a gaussian 
        'Convolution'      , ## SciPy' FFT-based convolution for functions 
        )

    # =========================================================================
except ImportError :
    # =========================================================================
    logger.warning ("No numpy/scipy is available: convolution is disabled")

    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    if not __all__ :
        logger.error ( "No convolution is available" )
        
# =============================================================================
##                                                                      The END 
# =============================================================================
