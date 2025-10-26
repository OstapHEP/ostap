#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
# @file ostap/histos/tests/test_histos_axes.py
# - Test TAxis decorators 
# ============================================================================= 
""" Test module for ostap/histos/axes.py
- It tests TAxis decorators 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import 
# ============================================================================= 
from   ostap.utils.root_utils import batch_env 
import ostap.histos.axes 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'ostap.test_histos_axes' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
logger.info ( 'Test decorators for ROOT.Taxis objects ')
# =============================================================================
## set batch form environment 
batch_env ( logger )
# =============================================================================


## Test decorators for TAxis objects 
def test_axes () : 
    """ Test decorators for TAxis objects 
    """
    logger.info ( 'Test decorators for ROOT.TAxis objetc' )

    nbins = 10
    xmin  = 0
    xmax  = 10
    
    ## (1) uniform axis 
    axis1 = ROOT.TAxis ( nbins , xmin  , xmax  ) 
    
    edges = [ random.uniform ( 0 , 10 ) for i in range ( nbins - 1 ) ]
    edges = sorted ( edges + [ xmin ] + [ xmax ] )
    axis2 = ROOT.TAxis.from_edges  ( edges , check = True ) 
    
    logger.info ( '1st axis    : %s' % axis1 )
    logger.info ( '2nd axis    : %s' % axis2 )

    logger.info ( '1st len     : %s' % len ( axis1 ) ) 
    logger.info ( '2nd len     : %s' % len ( axis2 ) ) 
    
    logger.info ( '1nd == 2nd? : %s ' %  ( axis1 == axis2 ) )
    logger.info ( '1nd != 2nd? : %s ' %  ( axis1 != axis2 ) )
    
    logger.info ( '1st uniform : %s' % axis1.uniform() )
    logger.info ( '2nd uniform : %s' % axis2.uniform() )



    logger.info ( '1st axis iterator: [%s] ' %  ( ', '.join ( str ( e ) for e in axis1 ) ) )
    logger.info ( '2nd axis iterator: [%s] ' %  ( ', '.join ( str ( e ) for e in axis2 ) ) )
    logger.info ( '1st axis reversed: [%s] ' %  ( ', '.join ( str ( e ) for e in reversed ( axis1 ) ) ) )
    logger.info ( '2nd axis reversed: [%s] ' %  ( ', '.join ( str ( e ) for e in reversed ( axis2 ) ) ) ) 
    
    logger.info ( '1st axis contains: %s ' %  [ ( i , i in axis1 ) for i in range ( -3 , nbins + 3 ) ] ) 
    logger.info ( '2nd axis contains: %s ' %  [ ( i , i in axis2 ) for i in range ( -3 , nbins + 3 ) ] ) 


    logger.info ( '1st axis[3],axis[7]: %s %s ' %  ( axis1[3] , axis1[7] ) )
    logger.info ( '2nd axis[3],axis[7]: %s %s ' %  ( axis2[3] , axis2[7] ) )

    logger.info ( '1st axis[3:7]     : %s ' %  axis1[3:7] )
    logger.info ( '2nd axis[3:7]     : %s ' %  axis2[3:7] )

    logger.info ( '1st axis[3:7:2]   : %s ' %  axis1[3:7:2] )
    logger.info ( '2nd axis[3:7:2]   : %s ' %  axis2[3:7:2] )

    logger.info ( '1st axis[4.1:]    : %s ' %  axis1[4.1:] )
    logger.info ( '2nd axis[4.1:]    : %s ' %  axis2[4.1:] )
     
    logger.info ( '1st axis[:4.3]    : %s ' %  axis1[:4.3] )
    logger.info ( '2nd axis[:4.3]    : %s ' %  axis2[:4.3] )
    
    logger.info ( '1st axis[1.2:5.2] : %s ' %  axis1[1.2:5.2] )
    logger.info ( '2nd axis[1.2:5.2] : %s ' %  axis2[1.2:5.2] )

    logger.info ( '1st axis bin_iterator: %s ' %  [ b for b in axis1.bin_iterator() ] )
    logger.info ( '2ns axis bin_iterator: %s ' %  [ b for b in axis2.bin_iterator() ] )
    
    logger.info ( '1st axis bin_edges   : %s ' %  [ b for b in axis1.bin_edges () ] )
    logger.info ( '2ns axis bin_edges   : %s ' %  [ b for b in axis2.bin_edges () ] )

    logger.info ( '1st axis items : %s ' %  [ b for b in axis1.items () ] )
    logger.info ( '2ns axis items : %s ' %  [ b for b in axis2.items () ] )


    logger.info ( '1st axis scale(10) : %s ' %  axis1.scale(10) )
    logger.info ( '2nd axis scale(10) : %s ' %  axis2.scale(10) )
    
    logger.info ( '1st axis *10 : %s ' %  ( axis1 * 10 ) )
    logger.info ( '2nd axis *10 : %s ' %  ( axis2 * 10 ) ) 
    
    logger.info ( '1st axis 10* : %s ' %  ( 10 * axis1 ) )
    logger.info ( '2nd axis 10* : %s ' %  ( 10 * axis2 ) ) 

    logger.info ( '1st split(2) : %s ' %  ( axis1.split ( 2 ) ) )
    logger.info ( '2nd split(2) : %s ' %  ( axis2.split ( 2 ) ) )
    
    logger.info ( '1st /2  : %s ' %  ( axis1 / 2 ) )
    logger.info ( '2nd /2  : %s ' %  ( axis2 / 2 ) )

    logger.info ( '1st //2 : %s ' %  ( axis1 // 2 ) )
    logger.info ( '2nd //2 : %s ' %  ( axis2 // 2 ) )

    logger.info ( '1st merge(3): %s ' %  axis1.merge ( 3 ) )
    logger.info ( '2nd merge(3): %s ' %  axis2.merge ( 3 ) )

    logger.info ( '1st merge(2): %s ' %  axis1.merge ( 2 ) )
    logger.info ( '2nd merge(2): %s ' %  axis2.merge ( 2 ) )

    logger.info ( '1st %%2 : %s ' %  ( axis1 % 2 ) )
    logger.info ( '2nd %%2 : %s ' %  ( axis2 % 2 ) )

    logger.info ( '1st %%3 : %s ' %  ( axis1 % 3 ) )
    logger.info ( '2nd %%3 : %s ' %  ( axis2 % 3 ) )

    ha = axis1 @ axis2
    hb = axis2 @ axis1

    logger.info ( '1st join(1,4) : %s ' % axis1.join ( 1,4) ) 
    logger.info ( '2nd join(1.4) : %s ' % axis2.join ( 1,4) ) 

    logger.info ( '1st range(1.0,4.0) : %s ' % axis1.range ( 1.,4.) ) 
    logger.info ( '2nd range(1.0,4.0) : %s ' % axis2.range ( 1.,4.) ) 

    logger.info ( '1st axis+0.1  : %s ' % ( axis1 + 0.1 ) ) 
    logger.info ( '2nd axis+0.1  : %s ' % ( axis2 + 0.1 ) ) 
    
    logger.info ( '1st -0.1 + axis : %s ' % ( -0.1 + axis1 ) ) 
    logger.info ( '2nd -0.1 + axis : %s ' % ( -0.1 + axis2 ) ) 

    logger.info ( '1st axis+(-0.1,-0.11) : %s ' % ( axis1 + (-0.1 ,-0.11 ) ) ) 
    logger.info ( '2nd axis+(-0.1,-0.11) : %s ' % ( axis2 + (-0.1 ,-0.11 ) ) ) 
     
    
# =============================================================================
if '__main__' == __name__ :

    test_axes () ## test decorators for ROOT.TAxis objexts 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
