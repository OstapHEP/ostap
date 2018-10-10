#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  combine.py
#  Few helper utilities for combining the correlated measurements
#  @see P.Avery "Combining measurements with correlated errors", CBX 95 55
#  @see http://www.phys.ufl.edu/~avery/fitting/error_correl.ps.gz
#  @see http://www.researchgate.net.publication/2345482_Combining_Measurements_with_Correlated_Errors
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

        self.cov2 = None  
        if isinstance  ( data , DATA ) : self.data = data
        else :
            self.data = DATA() 
            for i in range(N)  :
                e  = data[i]
                if hasattr ( e , 'cov2' ) and 0 < e.cov2() :
                    if not self.cov2 : self.cov2 = COV2() 
                    self.cov2[i,i] = e.cov2()  
                self.data[i] = data[i]

        self.covs = [] 
        ## add covariance matrix
        if not self.cov2 : self.cov2 = COV2()
        else             : self.covs.append ( COV2( self.cov2 )  )

        _covs = ( cov2, ) + args
        for c in _covs :

            c1 = COV2()
            c1.increment ( c )

            self.cov2.increment ( c1 )
            self.covs.append    ( COV2 ( c1 ) )

        self.combiner = COMBINER( self.data , self.cov2 )

    ## get the final result: combined measurement with total uncertainty
    #  @code 
    #  combiner = Combiner ( ... )
    #  print combiner.result()
    #  @endcode
    def result  ( self ) :
        """Get the final result: combined measurement with total combined uncertainty
        >>> combiner = Combiner ( ... )
        >>> print combiner.result() 
        """
        return self.combiner.result  ()

    ## get the weights used to get the result  
    #  @code
    #  combiner = Combiner ( ... )
    #  print combiner.weights ()
    #  @endcode 
    def weight  ( self ) :
        """Get the weights used to get the result
        >>> combiner = Combiner ( ... )
        >>> print combiner.weights ()         
        """
        return self.combiner.weights ()

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
        w = self.combiner.weights()
        for c in self.covs :
            a = c.sim ( w )
            r.append  ( a )
            
        return tuple ( r ) 
                    
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

    r1 = c1.result()
    e1 = c1.errComponents()

    r2 = c2.result()
    e2 = c2.errComponents()

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

    r1 = c1.result()
    e1 = c1.errComponents()

    r2 = c2.result()
    e2 = c2.errComponents()

    logger.info ( 'CORRELATED  : %.3f+-%.3f = (%.3f+-%.3f+-%.3f) ' % ( r2.value() , r2.error()  , r2.value() , e2[0] , e2[1] ) ) 
    logger.info ( 'UNCORRELATED: %.3f+-%.3f = (%.3f+-%.3f+-%.3f) ' % ( r1.value() , r1.error()  , r1.value() , e1[0] , e1[1] ) ) 

    logger.info ( 80*'*' ) 
    
    
# =============================================================================
# The END 
# =============================================================================
