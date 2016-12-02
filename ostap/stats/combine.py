#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  combine.py
#  Few helper utilities for combining correlated measurements
#  @see P.Avery "Combining measurements with correlated errors", CBX 95 55
#  @see http://www.phys.ufl.edu/~avery/fitting/error_correl.ps.gz
#  @see http://www.researchgate.net.publication/2345482_Combining_Measurements_with_Correlated_Errors
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-09-28
#  
# =============================================================================
""" Few helper utilities for combining correlated measurements
- see P.Avery ``Combining measurements with correlated errors'', CBX 95 55
- see http://www.phys.ufl.edu/~avery/fitting/error_correl.ps.gz
- see http://www.researchgate.net.publication/2345482_Combining_Measurements_with_Correlated_Errors
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
#  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
#  @date 2015-09-29
class Combine(object) :
    """Helper class to combine correlated measurements
    - see P.Avery ``Combining measurements with correlated errors'', CBX 95 55
    - see http://www.phys.ufl.edu/~avery/fitting/error_correl.ps.gz
    - see http://www.researchgate.net.publication/2345482_Combining_Measurements_with_Correlated_Errors
    """
    def __init__ ( self , data , cov2 , *args ) :
        ## types involved

        N        = len ( data )
        COMBINER = Ostap.Math.Combine( N , 'double' )
        DATA     = Ostap.Vector      ( N ) 
        COV2     = Ostap.SymMatrix   ( N ) 

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

    ##  get final result: combined measurement with total uncertainty 
    def result  ( self ) :
        """get final result: combined measurement with total uncertainty"""
        return self.combiner.result  ()

    ## get weigts used to get the result  
    def weight  ( self ) :
        """get weigts used to get the result"""
        return self.combiner.weights ()

    ## get the decomposition of the final variance 
    def covComponents ( self ) :
        """Get the decomposition of the final variance
        >>> c = ...
        >>> print c.result(), c.covComponents() 
        """
        r = []
        w = self.combiner.weights()
        for c in self.covs :
            a = c.sim ( w )
            r.append  ( a )
            
        return tuple ( r ) 
                    
    ## get the decomposition of the final variance 
    def errComponents ( self ) :
        """Get the decomposition of the final uncertainty
        >>> c = ...
        >>> print c.result(), c.errComponents() 
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

    logger.info ( 'Example from Avery:' ) 
    
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

    logger.info ( 'CORRELATED  : %s %s ' % ( c2.combiner.result(), c2.errComponents() ) ) 
    logger.info ( 'UNCORRELATED: %s %s ' % ( c1.combiner.result(), c1.errComponents() ) ) 
    
    logger.info ( 'Lb mass average:' ) 

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
    
    logger.info ( 'CORRELATED  : %s %s ' % ( c2.combiner.result(), c2.errComponents() ) ) 
    logger.info ( 'UNCORRELATED: %s %s ' % ( c1.combiner.result(), c1.errComponents() ) ) 

    
    logger.info ( 80*'*' ) 
    
    
# =============================================================================
# The END 
# =============================================================================
