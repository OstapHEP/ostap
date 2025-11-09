#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostath/integrator.py
#  Decoration modlel for C++ class Ostap::Math::Integrator
#  @see Ostap::Math::Integrator
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2022-12-18
# =============================================================================
""" Decoration modlel for C++ class Ostap::Math::Integrator
- see `Ostap.Math.Integrator`
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2022-12-18"
__all__     = (
    "Integrator"      , ## C++ numerical integrator
    'integral_ostap'  , ## 1D (C++) integrator 
    'integral2_ostap' , ## 2D (C++) integrator 
    'integral3_ostap' , ## 3D (C++) integrator 
)
# =============================================================================
from ostap.core.ostap_types import num_types 
from ostap.math.base        import Ostap, isfinite
from ostap.math.ve          import VE 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.integrator' )
else                       : logger = getLogger ( __name__                )
# =============================================================================
## actual C++ interator 
Integrator = Ostap.Math.Integrator
# =============================================================================
eps_abs = 1.49e-8
eps_rel = 1.49e-8
# =============================================================================
## Estimate 1D integrals using Ostap::Math::Integrator
#  @see Ostap::Math::Integrator
#  @see Ostap::Math::Integrator::integrate
#  @see Ostap::Math::Integrator::integrate_err
#  @see Ostap::Math::Integrator::integrate_infinity
#  @see Ostap::Math::Integrator::integrate_infinity_err
#  @see Ostap::Math::Integrator::integrate_to_infinity
#  @see Ostap::Math::Integrator::integrate_to_infinity_err
#  @see Ostap::Math::Integrator::integrate_from_infinity
#  @see Ostap::Math::Integrator::integrate_from_infinity_err
def integral_ostap ( fun               ,
                     xmin              ,
                     xmax              , * , 
                     args    = ()      ,
                     kwargs  = {}      ,
                     err     = False   ,
                     epsabs  = eps_abs ,
                     epsrel  = eps_rel ,
                     rule    = 0       ,
                     tag     = 0       ,
                     rescale = 0       ,
                     silent  = False   , **other ) :
    """ Estimate 1D integrals using Ostap::Math::Integrator
    - see Ostap.Math.Integrator
    - see Ostap.Math.Integrator.integrate
    - see Ostap.Math.Integrator.integrate_err
    - see Ostap.Math.Integrator.integrate_infinity
    - see Ostap.Math.Integrator.integrate_infinity_err
    - see Ostap.Math.Integrator.integrate_to_infinity
    - see Ostap.Math.Integrator.integrate_to_infinity_err
    - see Ostap.Math.Integrator.integrate_from_infinity
    - see Ostap.Math.Integrator.integrate_from_infinity_err
    """
    
    assert callable    ( fun ) , "Function  must be callable!"
    assert isinstance  ( xmin     , num_types )                 , "Invalid `xmin' %s" % xmin 
    assert isinstance  ( xmax     , num_types )                 , "Invalid `xmax' %s" % xmax  
    assert isinstance  ( epsabs   , float     ) and 0 <  epsabs , "Invalid 'epsabs'!"
    assert isinstance  ( epsrel   , float     ) and 0 <  epsrel , "Invalid 'epsrel'!"

    if other : logger.warning ( "Unknown extra  parameters are ignored: %s"  % [ k for k in other ] )

    ## wrap additional positional  and/or keyword arguments 
    if args or kwargs : func = lambda x : fun ( x , *args , **kwargs )
    else              : func = fun
    
    ## get the integrator 
    I     = Ostap.Math.Integrator ()
    
    finmin = isfinite ( xmin )
    finmax = isfinite ( xmax )
    
    if finmin and finmax :
        r , e = I.integrate_err ( func    , xmin    , xmax    , 
                                  tag     , rescale ,
                                  epsabs  , epsrel  ,
                                  rule    )
        
    elif finmin :

        if xmin < xmax : 
            r , e = I.integrate_to_infinity_err   ( func , xmin , tag , epsabs , epsrel )
        else : 
            r , e = I.integrate_from_infinity_err ( func , xmin , tag , epsabs , epsrel )
            r = -r
            
    elif finmax :
        
        if xmin < xmax :
            r , e = I.integrate_from_infinity_err ( func , xmax , tag , epsabs , epsrel )
        else :
            r , e = I.integrate_to_infinity_err   ( func , xmax , tag , epsabs , epsrel )
            r = -r 

    else :
    
        r , e = I.integrate_infinity_err ( func , tag , epsabs , epsrel )     
        r = r if xmin < xmax else -r
    
    ## requested precision
    rp = max ( epsabs , abs ( r ) * epsrel ) 
    if rp < e and not silent :
        logger.warning ( "Requested precision %.3g is not reached: %.3g" % ( rp , e ) ) 
    
    return VE ( r , e * e ) if err else r

# =================================================================================
## Estimate 2D integrals using Ostap::Math::Integrator
#  @see Ostap::Math::Integrator
#  @see Ostap::Math::Integrator::integrate2
#  @see Ostap::Math::Integrator::integrate2_err
def integral2_ostap ( fun               ,
                      xmin              ,
                      xmax              , 
                      ymin              ,
                      ymax              , * , 
                      args    = ()      ,
                      kwargs  = {}      ,
                      err     = False   ,
                      epsabs  = eps_abs ,
                      epsrel  = eps_rel ,
                      tag     = 0       ,
                      silent  = False   , **other ) :
    """ Estimate 2D integrals using Ostap::Math::Integrator
    - see Ostap.Math.Integrator
    - see Ostap.Math.Integrator.integrate2
    - see Ostap.Math.Integrator.integrate2_err
    """
    
    assert callable    ( fun ) , "Function  must be callable!"
    assert isinstance  ( xmin     , num_types ) and isfinite ( xmin ) , "Invalid `xmin' %s" % xmin 
    assert isinstance  ( xmax     , num_types ) and isfinite ( xmax ) , "Invalid `xmax' %s" % xmax  
    assert isinstance  ( ymin     , num_types ) and isfinite ( ymin ) , "Invalid `ymin' %s" % ymin 
    assert isinstance  ( ymax     , num_types ) and isfinite ( ymax ) , "Invalid `ymax' %s" % ymax  
    assert isinstance  ( epsabs   , float     ) and 0 <  epsabs       , "Invalid 'epsabs'!"
    assert isinstance  ( epsrel   , float     ) and 0 <  epsrel       , "Invalid 'epsrel'!"

    if other : logger.warning ( "Unknown extra  parameters are ignored: %s"  % [ k for k in other ] )

    ## wrap additional positional  and/or keyword arguments 
    if args or kwargs : func = lambda x : fun ( x , *args , **kwargs )
    else              : func = fun
    
    ## get the integrator 
    I     = Ostap.Math.Integrator ()
    
    r , e = I.integrate2_err ( func   ,
                               xmin   ,
                               xmax   , 
                               ymin   ,
                               ymax   , 
                               tag    , 
                               epsabs ,
                               epsrel )
    
    ## requested precision
    rp = max ( epsabs , abs ( r ) * epsrel ) 
    if rp < e and not silent :
        logger.warning ( "Requested precision %.3g is not reached: %.3g" % ( rp , e ) ) 
    
    return VE ( r , e * e ) if err else r

# =================================================================================
## Estimate 3D integrals using Ostap::Math::Integrator
#  @see Ostap::Math::Integrator
#  @see Ostap::Math::Integrator::integrate3
#  @see Ostap::Math::Integrator::integrate3_err
def integral3_ostap ( fun               ,
                      xmin              ,
                      xmax              , 
                      ymin              ,
                      ymax              , 
                      zmin              ,
                      zmax              , * , 
                      args    = ()      ,
                      kwargs  = {}      ,
                      err     = False   ,
                      epsabs  = eps_abs ,
                      epsrel  = eps_rel ,
                      tag     = 0       ,
                      silent  = False   , **other ) :
    """ Estimate 3D integrals using Ostap::Math::Integrator
    - see Ostap.Math.Integrator
    - see Ostap.Math.Integrator.integrate3
    - see Ostap.Math.Integrator.integrate3_err
    """
    
    assert callable    ( fun ) , "Function  must be callable!"
    assert isinstance  ( xmin     , num_types ) and isfinite ( xmin ) , "Invalid `xmin' %s" % xmin 
    assert isinstance  ( xmax     , num_types ) and isfinite ( xmax ) , "Invalid `xmax' %s" % xmax  
    assert isinstance  ( ymin     , num_types ) and isfinite ( ymin ) , "Invalid `ymin' %s" % ymin 
    assert isinstance  ( ymax     , num_types ) and isfinite ( ymax ) , "Invalid `ymax' %s" % ymax  
    assert isinstance  ( zmin     , num_types ) and isfinite ( zmin ) , "Invalid `zmin' %s" % zmin 
    assert isinstance  ( zmax     , num_types ) and isfinite ( zmax ) , "Invalid `zmax' %s" % zmax  
    assert isinstance  ( epsabs   , float     ) and 0 <  epsabs       , "Invalid 'epsabs'!"
    assert isinstance  ( epsrel   , float     ) and 0 <  epsrel       , "Invalid 'epsrel'!"

    if other : logger.warning ( "Unknown extra  parameters are ignored: %s"  % [ k for k in other ] )

    ## wrap additional positional  and/or keyword arguments 
    if args or kwargs : func = lambda x : fun ( x , *args , **kwargs )
    else              : func = fun
    
    ## get the integrator 
    I     = Ostap.Math.Integrator ()
    
    r , e = I.integrate3_err ( func   ,
                               xmin   ,
                               xmax   , 
                               ymin   ,
                               ymax   , 
                               zmin   ,
                               zmax   , 
                               tag    , 
                               epsabs ,
                               epsrel )
    
    ## requested precision
    rp = max ( epsabs , abs ( r ) * epsrel ) 
    if rp < e and not silent :
        logger.warning ( "Requested precision %.3g is not reached: %.3g" % ( rp , e ) ) 
    
    return VE ( r , e * e ) if err else r

                               
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                     The EMD   
# =============================================================================

