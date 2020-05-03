#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  combine.py
#  Few helper utilities for combining the correlated measurements,
#  aka ``BLUE'': Best Linear Unbiased Estimator
# 
#  @see P.Avery "Combining measurements with correlated errors", CBX 95 55
#  @see http://www.phys.ufl.edu/~avery/fitting/error_correl.ps.gz
#  @see http://www.researchgate.net.publication/2345482_Combining_Measurements_with_Correlated_Errors
#
#  @see Louis Lyons, Duncan Gibaut, Peter Clifford, 
#       "How to combine correlated estimates of a single physical quantity",
#        Nuclear Instruments and Methods in Physics Research Section A:
#        Accelerators, Spectrometers, Detectors and Associated Equipment
#        Volume 270, Issue 1, 1 July 1988, Pages 110-117
#  @see https://doi.org/10.1016/0168-9002(88)90018-6
# 
#  @see Ostap::Math::Combine
#  @code
#  x = VE ( 0.95 , 0.08**2 )        ## the first  measurement (with stat uncertainty)
#  y = VE ( 1.08 , 0.08**2 )        ## the second measurement (with stat uncertainty)
#  syst = Ostap.Math.SymMatrix2x2() ## systematic uncertainty (correlated and uncorrelated)
#  syst [0,0] = 0.10**2
#  syst [0,1] = 0.08**2
#  syst [1,1] = 0.10**2
#  combiner = Combine ( [x,y] , syst )
#  print combiner.result(), combiner.errComponents()
#  @endcode 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-09-28
#  
# =============================================================================
""" Few helper utilities for combining the correlated measurements
- see P.Avery ``Combining measurements with correlated errors'', CBX 95 55
- see http://www.phys.ufl.edu/~avery/fitting/error_correl.ps.gz
- see http://www.researchgate.net.publication/2345482_Combining_Measurements_with_Correlated_Errors

- see Ostap::Math::Combine

>>> x = VE ( 0.95 , 0.08**2 ) ## the first measurement   (with stat uncertainty)
>>> y = VE ( 1.08 , 0.08**2 ) ## the second measurement  (with stat uncertainty)
>>> syst = Ostap.Math.SymMatrix2x2() ## syst uncertainty (correlated and uncorrelated)
>>> syst [0,0] = 0.10**2
>>> syst [0,1] = 0.08**2
>>> syst [1,1] = 0.10**2
>>> combiner = Combine ( [x,y] , syst )
>>> print combiner.result(), combiner.errComponents()

"""
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2015-09-29"
__version__ = "$Revision$"
# =============================================================================
__all__     = (
    'Combine',
    ) 
# =============================================================================
import ROOT
import ostap.math.linalg
from   ostap.math.ve import VE,Ostap  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.combine' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
## @class Combine
#  Helper class to combine correlated measurements
#  @see Ostap::Math::Combine
#  @see P.Avery "Combining measurements with correlated errors", CBX 95 55
#  @see http://www.phys.ufl.edu/~avery/fitting/error_correl.ps.gz
#  @see http://www.researchgate.net.publication/2345482_Combining_Measurements_with_Correlated_Errors
#  @see Ostap::Math::Combine
#  @code
#  x = VE ( 0.95 , 0.08**2 ) ## the first  measurement (with stat uncertainty)
#  y = VE ( 1.08 , 0.08**2 ) ## the second measurement (with stat uncertainty)
#  ## systematic uncertainty (correlated and uncorrelated)
#  syst = Ostap.Math.SymMatrix2x2()
#  syst [0,0] = 0.10**2
#  syst [0,1] = 0.08**2
#  syst [1,1] = 0.10**2
#  combiner = Combine ( [x,y] , syst )
#  print combiner.result(), combiner.errComponents()
#  @endcode 
#  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
#  @date 2015-09-29
class Combine(object) :
    """Helper class to combine the correlated measurements
    - see P.Avery ``Combining measurements with correlated errors'', CBX 95 55
    - see http://www.phys.ufl.edu/~avery/fitting/error_correl.ps.gz
    - see http://www.researchgate.net.publication/2345482_Combining_Measurements_with_Correlated_Errors
    - see Ostap::Math::Combine
    
    >>> x = VE ( 0.95 , 0.08**2 ) ## the first measurement   (with stat uncertainty)
    >>> y = VE ( 1.08 , 0.08**2 ) ## the second measurement  (with stat uncertainty)
    >>> syst = Ostap.Math.SymMatrix2x2() ## syst uncertainty (correlated and uncorrelated)
    >>> syst [0,0] = 0.10**2
    >>> syst [0,1] = 0.08**2
    >>> syst [1,1] = 0.10**2
    >>> combiner = Combine ( [x,y] , syst )
    >>> print combiner.result(), combiner.errComponents()
    """
    def __init__ ( self , data , cov2 , *args ) :
        """Create the combiner
        - data : vector of data (floats or VE-objects)
        - cov2 : covariance matrix
        - args : additional covariance matrices for different uncertainty components
        """

        N        = len ( data )
        COMBINER = Ostap.Math.Combine ( N , 'double' )
        DATA     = Ostap.Vector       ( N ) 
        COV2     = Ostap.SymMatrix    ( N ) 

        self.__cov2 = None  
        if isinstance  ( data , DATA ) : self.__data = data
        else :
            self.__data = DATA() 
            for i in range ( N )  :
                e  = data[i]
                if hasattr ( e , 'cov2' ) and 0 < e.cov2() :
                    if not self.__cov2 : self.__cov2 = COV2() 
                    self.__cov2[i,i] = e.cov2()  
                self.__data[i] = data[i]

        self.__covs = [] 
        ## add covariance matrix
        if not self.__cov2 : self.__cov2 =        COV2 (             )
        else               : self.__covs.append ( COV2 ( self.__cov2 )  )

        _covs = ( cov2, ) + args
        for c in _covs :

            c1 = COV2()
            c1 = c1 + c 

            self.__cov2 = self.__cov2 + c1               
            self.__covs.append ( COV2 ( c1 ) )

        self.__covs   = tuple ( self.__covs ) 
        self.combiner = COMBINER( self.data , self.cov2 )

    @property
    def data ( self )  :
        """``data'' : input data/measurements"""
        return self.__data

    @property
    def covs ( self )  :
        """``covs'' : covariance matrices"""
        return self.__covs
    
    @property
    def cov2 ( self )  :
        """``cov2'' : overall covariance matrix"""
        return self.__cov2 

    @property
    def result ( self ) :
        """``result'' : the final  result"""
        return self.combiner.result  ()
    
    @property 
    def weights ( self ) :
        """``weights'' : the weights"""
        return self.combiner.weights ()

    # =========================================================================
    ## get the decomposition of the final variances 
    #  @code
    #  combiner = Combiner ( ... ) 
    #  print c.result(), c.covComponents()
    #  @endcode
    def covComponents ( self ) :
        """Get the decomposition of the final variances
        >>> combiner = Combiner ( ... ) 
        >>> print c.result(), c.covComponents() 
        """
        r = []
        w = self.weights
        for c in self.covs :
            a = c.sim ( w )
            r.append  ( a )
            
        return tuple ( r ) 

    # =========================================================================
    ## get the decomposition of the final uncertainty
    #  @code 
    #  combiner = ...
    #  print combiner.result(), combiner.errComponents()
    #  @endcode 
    def errComponents ( self ) :
        """Get the decomposition of the final uncertainty
        >>> combiner = ...
        >>> print combiner.result(), combiner.errComponents()
        """
        
        import math
        def _sqrt_ ( c ) :
            if 0 <= c : return math.sqrt ( c )
            return -1.0 * math.sqrt( abs  ( c ) )
        
        cc = self.covComponents()
        ce = [ _sqrt_ ( c ) for c in cc ]
        
        return tuple ( ce ) 

    @property
    def variances  ( self ) :
        """``variances'' : get variance components"""
        return self.covComponents()
    
    @property
    def errors  ( self ) :
        """``errors'' : get error components"""
        return self.errComponents()
    

        
# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    logger.info ( ' - Example from Avery:' ) 
    
    x = VE ( 0.95 , 0.08**2 )
    y = VE ( 1.08 , 0.08**2 )

    syst1 = Ostap.Math.SymMatrix2x2()
    syst1[0,0] = 0.08**2
    syst1[1,1] = 0.08**2

    syst2 = Ostap.Math.SymMatrix2x2()
    syst2[0,0] = 0.08**2
    syst2[0,1] = 0.08**2
    syst2[1,1] = 0.08**2

    c1 = Combine( [x,y] , syst1 )
    c2 = Combine( [x,y] , syst2 )

    r1 = c1.result
    e1 = c1.errors 

    r2 = c2.result
    e2 = c2.errors

    logger.info ( 'CORRELATED  : %.3f+-%.3f = (%.3f+-%.3f+-%.3f) ' % ( r2.value() , r2.error()  , r2.value() , e2[0] , e2[1] ) ) 
    logger.info ( 'UNCORRELATED: %.3f+-%.3f = (%.3f+-%.3f+-%.3f) ' % ( r1.value() , r1.error()  , r1.value() , e1[0] , e1[1] ) ) 
    
    logger.info ( ' - Lambda_b mass average:' ) 

    x = VE(5619.44 , 0.70**2 )
    y = VE(5619.44 , 0.13**2 )

    syst1 = Ostap.Math.SymMatrix2x2()
    syst1[0,0] = 0.30**2
    syst1[1,1] = 0.45**2
    
    
    syst2 = Ostap.Math.SymMatrix2x2()
    syst2[0,0] = 0.30**2
    syst2[1,1] = 0.45**2
    syst2[0,1] = 0.30*0.45

    c1 = Combine( [x,y] , syst1 )
    c2 = Combine( [x,y] , syst2 )

    r1 = c1.result
    e1 = c1.errors 

    r2 = c2.result
    e2 = c2.errors 

    logger.info ( 'CORRELATED  : %.3f+-%.3f = (%.3f+-%.3f+-%.3f) ' % ( r2.value() , r2.error()  , r2.value() , e2[0] , e2[1] ) ) 
    logger.info ( 'UNCORRELATED: %.3f+-%.3f = (%.3f+-%.3f+-%.3f) ' % ( r1.value() , r1.error()  , r1.value() , e1[0] , e1[1] ) ) 

    logger.info ( 80*'*' ) 
    
    
# =============================================================================
##                                                                      The END 
# =============================================================================
