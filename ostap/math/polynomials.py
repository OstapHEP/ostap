#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/polynomials.py
#  Module with some useful utilities for dealing with polymonials
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-12-01"
__all__     = ()
# =============================================================================
import  ROOT
from    ostap.math.base    import Ostap
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.polynomials' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================


# =============================================================================
## factory for deserialization of simple polynomians 
#  @see Ostap.Math.Chebyshev  
#  @see Ostap.Math.ChebyshevU  
#  @see Ostap.Math.Hermite     
#  @see Ostap.Math.Legendre
#  @see Ostap.Math.PLegendre
def pN_factory ( klass , *args ) :
    """Factory for deserialization of simple polynomians
    - see Ostap.Math.Chebyshev  
    - see Ostap.Math.ChebyshevU 
    - see Ostap.Math.Hermite    
    - see Ostap.Math.Legendre  
    - see Ostap.Math.PLegendre  
    """
    return klass ( *args ) 

# =============================================================================
## reduce simple polynomials
#  @see Ostap.Math.Chebyshev  
#  @see Ostap.Math.ChebyshevU  
#  @see Ostap.Math.Hermite     
#  @see Ostap.Math.Legendre   
def pN_reduce ( p ) :
    """Reduce simple polynomials
    - see Ostap.Math.Chebyshev  
    - see Ostap.Math.ChebyshevU 
    - see Ostap.Math.Hermite    
    - see Ostap.Math.Legendre  

    """
    return pN_factory ,  ( type ( p ) , p.degree() ) 

# =============================================================================
## reduce simple polynomials
#  @see Ostap.Math.PLegendre   
def pLM_reduce ( p ) :
    """Reduce simple polynomials
    - see Ostap.Math.PLegendre  

    """
    return pN_factory ,  ( type ( p ) , p.L() , p.M() ) 

for t in ( Ostap.Math.Chebyshev  ,
           Ostap.Math.ChebyshevU , 
           Ostap.Math.Hermite    , 
           Ostap.Math.Legendre   ) :
    
    t.__reduce__ = pN_reduce
    
Ostap.Math.PLegendre.__reduce__ = pLM_reduce 


# =============================================================================
## factory for deserisalization of polynomials with parameters
#  @see Ostap::Math::Polynomial
#  @see Ostap::Math::ChebyshevSum
#  @see Ostap::Math::LegendreSum
#  @see Ostap::Math::HermiteSum
#  @see Ostap::Math::Bernstein
#  @see Ostap::Math::BernsteinEven
def poly_factory ( klass , params , *args ) :
    """Factory for deserisalization of polynomials with parameters
    - see Ostap.Math.Polynomial
    - see Ostap.Math.ChebyshevSum
    - see Ostap.Math.LegendreSum
    - see Ostap.Math.HermiteSum
    - see Ostap.Math.Bernstein
    - see Ostap.Math.BernsteinEven
    """
    from ostap.math.base import doubles 
    return klass ( doubles ( params ) , *args ) 

# =============================================================================
## Reduce polynomials with parameters
#  @see Ostap::Math::Polynomial
#  @see Ostap::Math::ChebyshevSum
#  @see Ostap::Math::LegendreSum
#  @see Ostap::Math::HermiteSum
#  @see Ostap::Math::Bernstein
#  @see Ostap::Math::BernsteinEven
#  @see Ostap::Math::Positive 
def poly_reduce ( p ) : 
    """Reduce polynomials with parameters
    - see Ostap.Math.Polynomial
    - see Ostap.Math.ChebyshevSum
    - see Ostap.Math.LegendreSum
    - see Ostap.Math.HermiteSum
    - see Ostap.Math.Bernstein
    - see Ostap.Math.BernsteinEven
    - see Ostap.Math.Positive 
    - see Ostap.Math.PositiveEven 
    """
    return poly_factory , ( type ( p ) ,
                            tuple ( p.pars() ) ,
                            p.xmin () ,
                            p.xmax () )

    
for t in (  Ostap.Math.Polynomial     , 
            Ostap.Math.ChebyshevSum   , 
            Ostap.Math.LegendreSum    , 
            Ostap.Math.HermiteSum     , 
            Ostap.Math.Bernstein      , 
            Ostap.Math.BernsteinEven  ,
            Ostap.Math.Positive       , 
            Ostap.Math.PositiveEven   ) :
    
    t.__reduce__ = poly_reduce
    

# =============================================================================
## reduce monotonic polymomial
#  @see Ostap::Math::Monotonic 
def pm_reduce ( p ) :
    """reduce monotonic polymomial
    - see Ostap.Math.Monotonic
    """
    return poly_factory , ( type ( p ) ,
                            tuple ( p.pars() ) ,
                            p.xmin () ,
                            p.xmax () ,
                            True if p.increasing() else False )

# =============================================================================
## reduce convex polynomial
#  @see Ostap::Math::Convex  
def pc_reduce ( p ) :
    """reduce convex polymomial
    - see Ostap.Math.Convex
    """
    return poly_factory , ( type ( p ) ,
                            tuple ( p.pars() ) ,
                            p.xmin () ,
                            p.xmax () ,
                            True if p.increasing () else False , 
                            True if p.convex     () else False ) 

# =============================================================================
## reduce convex-only polynomial
#  @see Ostap::Math::ConvexOnly
def pco_reduce ( p ) :
    """reduce convex-only polynomial
    - see Ostap.Math.ConvexOnly
    """
    return poly_factory , ( type ( p ) ,
                            tuple ( p.pars() ) ,
                            p.xmin () ,
                            p.xmax () ,
                            True if p.convex     () else False ) 


Ostap.Math.Monotonic  .__reduce__ =  pm_reduce
Ostap.Math.Convex     .__reduce__ =  pc_reduce
Ostap.Math.ConvexOnly .__reduce__ = pco_reduce



# =============================================================================
## factory for deserisalization of splines 
#  @see Ostap::Math::BSPline 
#  @see Ostap::Math::PositiveSpline 
#  @see Ostap::Math::MonotonicSpline 
#  @see Ostap::Math::ConvexSpline 
#  @see Ostap::Math::ConvexOnlySpline 
def sp_factory ( klass , knots , pars , *args ) :
    """Factory for deserisalization of splines 
    - see Ostap.Math.BSPline 
    - see Ostap.Math.PositiveSpline 
    - see Ostap.Math.MonotonicSpline 
    - see Ostap.Math.ConvexSpline 
    - see Ostap.Math.ConvexOnlySpline 
    
    """
    from ostap.math.base import doubles 
    return klass ( doubles ( knots) , doubles ( pars ) , *args ) 


# =============================================================================
## factory for deserisalization of splines 
#  @see Ostap::Math::BSPline 
#  @see Ostap::Math::PositiveSpline 
def sp_reduce (  sp ) :
    """Factory for deserisalization of splines 
    - see Ostap.Math.BSPline 
    - see Ostap.Math.PositiveSpline 
    """
    return sp_factory , ( type  ( sp ) ,
                          tuple ( sp.knots() ) ,
                          tuple ( sp.pars () ) ) 

Ostap.Math.BSpline        . __reduce__ = sp_reduce 
Ostap.Math.PositiveSpline . __reduce__ = sp_reduce 


# =============================================================================
## factory for deserisalization of splines 
#  @see Ostap::Math::MonotonicSpline 
def spm_reduce (  sp ) :
    """Factory for deserisalization of splines 
    - see Ostap.Math.MonotonicSpline 
    """
    return sp_factory , ( type  ( sp ) ,
                          tuple ( sp.knots() ) ,
                          tuple ( sp.pars () ) , 
                          True if sp.increasing () else False )
                          
Ostap.Math.MonotonicSpline . __reduce__ = spm_reduce 

# =============================================================================
## factory for deserisalization of splines 
#  @see Ostap::Math::ConvexSpline 
def spc_reduce (  sp ) :
    """Factory for deserisalization of splines 
    - see Ostap.Math.ConvexSpline 
    """
    return sp_factory , ( type  ( sp ) ,
                          tuple ( sp.knots() ) ,
                          tuple ( sp.pars () ) , 
                          True if sp.increasing () else False , 
                          True if sp.convex     () else False )
                          
Ostap.Math.ConvexSpline . __reduce__ = spc_reduce 



# =============================================================================
## factory for deserisalization of splines 
#  @see Ostap::Math::ConvexOnlySpline 
def spco_reduce (  sp ) :
    """Factory for deserisalization of splines 
    - see Ostap.Math.ConvexOnlySpline 
    """
    return sp_factory , ( type  ( sp ) ,
                          tuple ( sp.knots() ) ,
                          tuple ( sp.pars () ) , 
                          True if p.convex     () else False )
                          
Ostap.Math.ConvexOnlySpline . __reduce__ = spco_reduce 


# =============================================================================
## decorated classes 
_decorated_classes_  = (
    ##
    Ostap.Math.Chebyshev        ,
    Ostap.Math.ChebyshevU       , 
    Ostap.Math.Hermite          , 
    Ostap.Math.Legendre         , 
    Ostap.Math.PLegendre        ,    
    ## 
    Ostap.Math.Polynomial       , 
    Ostap.Math.ChebyshevSum     , 
    Ostap.Math.LegendreSum      , 
    Ostap.Math.HermiteSum       , 
    Ostap.Math.Bernstein        , 
    Ostap.Math.BernsteinEven    ,
    ##
    Ostap.Math.Positive         ,
    Ostap.Math.PositiveEven     ,
    ## 
    Ostap.Math.Monotonic        ,
    Ostap.Math.Convex           ,
    Ostap.Math.ConvexOnly       ,
    ##
    Ostap.Math.BSpline          ,
    Ostap.Math.PositiveSpline   , 
    Ostap.Math.MonotonicSpline  , 
    Ostap.Math.ConvexSpline     , 
    Ostap.Math.ConvexOnlySpline ,
    ##
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END 
# =============================================================================

