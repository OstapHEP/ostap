#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/reduce.py
#  Module with some useful utilities for reducing some math objects 
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-12-01"
__all__     = ()
# =============================================================================
from    ostap.math.base        import Ostap, doubles 
from    ostap.core.ostap_types import sequence_types 
import  ROOT, array 
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.reduce' )
else                       : logger = getLogger ( __name__            )
# =============================================================================

# =============================================================================
## Simple (basic) polynomials 
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

for t in ( Ostap.Math.Chebyshev  ,
           Ostap.Math.ChebyshevU , 
           Ostap.Math.Hermite    , 
           Ostap.Math.Legendre   ) :
    
    t.__reduce__ = pN_reduce
    
Ostap.Math.PLegendre.__reduce__ = pLM_reduce 

# =============================================================================
## Regular polynomials  
# =============================================================================   

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
                            array.array ( 'd' ,  p.pars() ) ,
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
## Specific forms of Bernstein  polynomials 
# =============================================================================


# =============================================================================
## reduce monotonic polynomial
#  @see Ostap::Math::Monotonic 
def pm_reduce ( p ) :
    """reduce monotonic polynomial
    - see Ostap.Math.Monotonic
    """
    return poly_factory , ( type ( p ) ,
                            array.array ( 'd' , p.pars() ) ,
                            p.xmin () ,
                            p.xmax () ,
                            True if p.increasing() else False )

# =============================================================================
## reduce convex polynomial
#  @see Ostap::Math::Convex  
def pc_reduce ( p ) :
    """reduce convex polynomial
    - see Ostap.Math.Convex
    """
    return poly_factory , ( type ( p ) ,
                            array.array ( 'd' ,  p.pars() ) ,
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
                            array.array ( 'd' , p.pars() ) ,
                            p.xmin () ,
                            p.xmax () ,
                            True if p.convex     () else False ) 


Ostap.Math.Monotonic  .__reduce__ =  pm_reduce
Ostap.Math.Convex     .__reduce__ =  pc_reduce
Ostap.Math.ConvexOnly .__reduce__ = pco_reduce


# =============================================================================
## B-splines 
# =============================================================================

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
                          array.array ( 'd' , sp.knots() ) ,
                          array.array ( 'd' , sp.pars () ) ) 


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
                          array.array ( 'd' , sp.knots() ) ,
                          array.array ( 'd' , sp.pars () ) , 
                          True if sp.increasing () else False )
                          

# =============================================================================
## factory for deserisalization of splines 
#  @see Ostap::Math::ConvexSpline 
def spc_reduce (  sp ) :
    """Factory for deserisalization of splines 
    - see Ostap.Math.ConvexSpline 
    """
    return sp_factory , ( type  ( sp ) ,
                          array.array ( 'd' , sp.knots() ) ,
                          array.array ( 'd' , sp.pars () ) , 
                          True if sp.increasing () else False , 
                          True if sp.convex     () else False )
                          


# =============================================================================
## factory for deserisalization of splines 
#  @see Ostap::Math::ConvexOnlySpline 
def spco_reduce (  sp ) :
    """Factory for deserisalization of splines 
    - see Ostap.Math.ConvexOnlySpline 
    """
    return sp_factory , ( type  ( sp ) ,
                          array.array ( 'd' , sp.knots() ) ,
                          array.array ( 'd' , sp.pars () ) , 
                          True if p.convex     () else False )
                          
Ostap.Math.MonotonicSpline . __reduce__ = spm_reduce 
Ostap.Math.ConvexSpline . __reduce__ = spc_reduce 
Ostap.Math.ConvexOnlySpline . __reduce__ = spco_reduce 



# =============================================================================
## Interpolation stuff 
# =============================================================================



# ============================================================================
## factory for deserialisation of interpolation abscissas 
def abs_factory ( arg , *args ) :
    """Factory for deserialisation of interpolation abscissas
    """
    if isinstance ( arg , sequence_types ) :
        vals = doubles ( arg )
        return Ostap.Math.Interpolation.Abscissas ( vals , *args  )
    
    return Ostap.Math.Interpolation.Abscissas ( arg , *args ) 


# =============================================================================
## Reduce interpolation abscissas 
def abs_reduce ( a ) :
    """Reduce interpolation abscissas 
    """
    at = a.atype()
    if at in ( Ostap.Math.Interpolation.Abscissas.Uniform    ,
               Ostap.Math.Interpolation.Abscissas.Chebyshev  ,
               Ostap.Math.Interpolation.Abscissas.Chebyshev2 ) :
        return abs_factory, ( a.n () , a.xmin() , a.xmax () , int ( at ) )

    return abs_factory, ( array.array ('d' , a.x() ) , ) 

# ============================================================================
## the factory for serialisation of the interpolation table 
def tab_factory ( abscissas , values ) :
    """The factory for serialisation of the interpolation table
    """
    return Ostap.Math.Interpolation.Table ( abscissas , doubles ( values ) ) 
## ===========================================================================
## Reduce the interpolation table 
def tab_reduce ( table ) :
    """Reduce the interpolation table"""
    return tab_factory , ( table.abscissas ()              ,
                           array.array ( 'd' , table.values () ) )


Ostap.Math.Interpolation.Abscissas . __reduce__ = abs_reduce 
Ostap.Math.Interpolation.Table     . __reduce__ = tab_reduce


# ============================================================================
## the factory for serialisation of the interpolation objects 
def int_factory ( klass , abscissas , values , *args ) :
    """The factory for serialisation of the interpolation table
    """
    the_table = Ostap.Math.Interpolation.Table ( abscissas , doubles ( values ) )
    return klass ( the_table , *args )

## ===========================================================================
## Reduce the interpolation object 
def int_reduce ( table ) :
    """Reduce the interpolation object"""
    return int_factory , ( type ( table )                  ,
                           table.abscissas ()              ,
                           array.array ( 'd' , table.values () ) ) 

## ===========================================================================
## Reduce the interpolation Floater-Hormann interpolant 
def intfh_reduce ( table ) :
    """Reduce the Floater-Hormann interpolant"""
    return int_factory , ( type ( table )                  ,
                           table.abscissas ()              ,
                           array.array ( 'd' , table.values () ) , 
                           table.d ()                      ) 

for t in ( Ostap.Math.Neville     ,
           Ostap.Math.Lagrange    , 
           Ostap.Math.Newton      , 
           Ostap.Math.Barycentric , 
           Ostap.Math.Berrut1st   , 
           Ostap.Math.Berrut2nd   ) :
    
    t.__reduce__ = int_reduce 

Ostap.Math.FloaterHormann. __reduce__ = intfh_reduce 



# =============================================================================
## Dalitz' objects 
# =============================================================================   


# =============================================================================   
## Factory for deserialization of Dalitz' objects
def dalitz_factory ( klass , *params ) :
    """Factory for deserialization of `Dalitz` objects
    """
    return klass ( *params )

# ============================================================================
## Serialise class <code>Ostap::Kinematics::Dalitz0</code>
#  @see Ostap::Kinematcis.Dalitz0
def _dalitz0_reduce_ ( dalitz ) :
    """Serialise class `Ostap.Kinematics.Dalitz0`
    - see Ostap.Kinematics.Dalitz0
    """
    return dalitz_factory , ( type ( dalitz ) ,
                              dalitz.m1 () ,
                              dalitz.m2 () ,
                              dalitz.m3 () ) 

# ============================================================================
## Serialise class <code>Ostap::Kinematics::Dalitz</code>
#  @see Ostap::Kinematcis.Dalitz
def _dalitzm_reduce_ ( dalitz ) :
    """Serialise class `Ostap.Kinematics.Dalitz`
    - see Ostap.Kinematics.Dalitz
    """
    return dalitz_factory , ( type ( dalitz ) ,
                              dalitz.M  () , 
                              dalitz.m1 () ,
                              dalitz.m2 () ,
                              dalitz.m3 () ) 


Ostap.Kinematics.Dalitz0. __reduce__ = _dalitz0_reduce_ 
Ostap.Kinematics.Dalitz . __reduce__ = _dalitzm_reduce_   



# =============================================================================
## (Redefine standard constructor to allow usage of python lists&tuples)
#  Lists and tuples are converted on flight to :
# - std::vector<double> 
def _new_init_ ( t ,  *args )  :
    """(Redefine standard constructor to allow usage of python lists&tuples)
    Lists and tuples are  converted on flight to :
    - std::vector<double> 
    """
    from ostap.math.base        import doubles  , VCT_TYPES 
    from ostap.core.ostap_types import Generator, Sequence, list_types  
    
    largs = list (  args )

    for i , arg in enumerate ( largs ) :
        
        if   isinstance ( arg , VCT_TYPES  ) : continue 
        
        if   isinstance ( arg , Generator  ) : pass
        elif isinstance ( arg , Sequence   ) : pass
        elif isinstance ( arg , list_types ) : pass
        else :
            continue
        
        try: 
            _arg = doubles  ( arg  )
            largs [ i ] = _arg
            continue 
        except TypeError : pass
                
    targs = tuple ( largs )
        
    ## use old constructor 
    return t._old_init_ ( *targs ) 

# =============================================================================
for p in ( Ostap.Math.Polynomial    , 
           Ostap.Math.ChebyshevSum  , 
           Ostap.Math.LegendreSum   , 
           Ostap.Math.HermiteSum    ) :
    
    if not hasattr ( p , '_old_init_' ) :
        
        p._old_init_ = p.__init__
        ## Modifed constructor to allow python lists/tuples
        def _p_new_init_ ( s ,  *args ) :
            """Modifed constructor to allow python lists/tuples
            """
            _new_init_ ( s , *args )
            
        _p_new_init_.__doc__ += '\n' +   _new_init_.__doc__ 
        _p_new_init_.__doc__ += '\n' + p._old_init_.__doc__ 
        p.__init__ = _p_new_init_ 



# =============================================================================
## decorated classes 
_decorated_classes_  = (
    ## 
    Ostap.Math.Chebyshev               ,
    Ostap.Math.ChebyshevU              , 
    Ostap.Math.Hermite                 , 
    Ostap.Math.Legendre                , 
    Ostap.Math.PLegendre               ,
    ## 
    Ostap.Math.Polynomial              , 
    Ostap.Math.ChebyshevSum            , 
    Ostap.Math.LegendreSum             , 
    Ostap.Math.HermiteSum              , 
    Ostap.Math.Bernstein               , 
    Ostap.Math.BernsteinEven           ,
    Ostap.Math.Positive                , 
    Ostap.Math.PositiveEven            , 
    ##
    Ostap.Math.Monotonic               , 
    Ostap.Math.Convex                  , 
    Ostap.Math.ConvexOnly              , 
    ##
    Ostap.Math.BSpline                 , 
    Ostap.Math.PositiveSpline          ,
    ## 
    Ostap.Math.MonotonicSpline         ,
    Ostap.Math.ConvexSpline            ,
    Ostap.Math.ConvexOnlySpline        , 
    ##
    Ostap.Math.Interpolation.Abscissas , 
    Ostap.Math.Interpolation.Table     , 
    ##
    Ostap.Math.Neville                 ,
    Ostap.Math.Lagrange                , 
    Ostap.Math.Newton                  , 
    Ostap.Math.Barycentric             , 
    Ostap.Math.Berrut1st               , 
    Ostap.Math.Berrut2nd               , 
    Ostap.Math.FloaterHormann          , 
    ##
    Ostap.Kinematics.Dalitz0           , 
    Ostap.Kinematics.Dalitz            ,
    )

# =============================================================================
## new methdods 
_new_methods_       = (
    ##
    Ostap.Math.Chebyshev               . __reduce__  , 
    Ostap.Math.ChebyshevU              . __reduce__  , 
    Ostap.Math.Hermite                 . __reduce__  , 
    Ostap.Math.Legendre                . __reduce__  , 
    Ostap.Math.PLegendre               . __reduce__  , 
    ##
    Ostap.Math.Polynomial              . __reduce__  , 
    Ostap.Math.ChebyshevSum            . __reduce__  ,  
    Ostap.Math.LegendreSum             . __reduce__  ,  
    Ostap.Math.HermiteSum              . __reduce__  , 
    Ostap.Math.Bernstein               . __reduce__  , 
    Ostap.Math.BernsteinEven           . __reduce__  , 
    Ostap.Math.Positive                . __reduce__  , 
    Ostap.Math.PositiveEven            . __reduce__  , 
    ## 
    Ostap.Math.Monotonic               . __reduce__  , 
    Ostap.Math.Convex                  . __reduce__  ,
    Ostap.Math.ConvexOnly              . __reduce__  ,
    ##
    Ostap.Math.BSpline                 . __reduce__  , 
    Ostap.Math.PositiveSpline          . __reduce__  ,
    ## 
    Ostap.Math.MonotonicSpline         . __reduce__  , 
    Ostap.Math.ConvexSpline            . __reduce__  ,
    Ostap.Math.ConvexOnlySpline        . __reduce__  ,
    ##
    Ostap.Math.Interpolation.Abscissas . __reduce__  , 
    Ostap.Math.Interpolation.Table     . __reduce__  , 
    ##
    Ostap.Math.Neville                 . __reduce__  ,
    Ostap.Math.Lagrange                . __reduce__  ,  
    Ostap.Math.Newton                  . __reduce__  ,  
    Ostap.Math.Barycentric             . __reduce__  ,  
    Ostap.Math.Berrut1st               . __reduce__  ,  
    Ostap.Math.Berrut2nd               . __reduce__  , 
    Ostap.Math.FloaterHormann          . __reduce__  ,
    ##
    Ostap.Kinematics.Dalitz0           . __reduce__  ,
    Ostap.Kinematics.Dalitz            . __reduce__  ,
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END 
# =============================================================================

