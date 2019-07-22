#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_math_dalitz.py
#  Test module for the file ostap/math/dalitz.py
#  @see Ostap::Kinematics::Dalitz
# ============================================================================= 
""" Test module for ostap/math/dalitz.py
- see Ostap.Kinematics.Dalitz
"""
# ============================================================================= 
from __future__ import print_function
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_math_dalitz'   )
else                       : logger = getLogger ( __name__             )
# ============================================================================= 
import ROOT, random 
from   ostap.core.core import Ostap
import ostap.math.kinematic
import ostap.math.dalitz

LV = Ostap.LorentzVector
CT = Ostap.Kinematics.cos_theta

p1 = LV ( 0 ,  1 , 10 , 15 )
p2 = LV ( 1 ,  1 , 13 , 18 )
p3 = LV ( 0 , -2 , 15 , 25 )

P  = p1 
Q  = p1 + p2
D  = p1 + p2 + p3

s1 = ( p1 + p2 ).M2()
s2 = ( p2 + p3 ).M2()
s  = ( p1 + p2 + p3 ).M2()

d  = Ostap.Kinematics.Dalitz ( D.M () , p1.M () , p2.M () , p3.M () )


def test_dalitz1 () : 
    

    logger.info ( 'E1     : %s' , d.E1 ( s1 ,  s2 ) )
    logger.info ( 'E2     : %s' , d.E2 ( s1 ,  s2 ) )
    logger.info ( 'E3     : %s' , d.E3 ( s1 ,  s2 ) )
    
    logger.info ( 'P1     : %s' , d.P1 ( s1 ,  s2 ) )
    logger.info ( 'P2     : %s' , d.P2 ( s1 ,  s2 ) )
    logger.info ( 'P3     : %s' , d.P3 ( s1 ,  s2 ) )

    ct_12 = d.cos_12 (  s1 , s2 ) 
    logger.info ( 'cos_12 : %s [R=%s] ' % ( ct_12 , ct_12 / CT ( p1 , p2 , p1 + p2 + p3 ) ) ) 
    
    ct_23 = d.cos_23 (  s1 , s2 ) 
    logger.info ( 'cos_23 : %s [R=%s] ' % ( ct_23 , ct_23 / CT ( p2 , p3 , p1 + p2 + p3 ) ) )
    
    ct_31 = d.cos_31 (  s1 , s2 ) 
    logger.info ( 'cos_31 : %s [R=%s] ' % ( ct_31 , ct_31 / CT ( p3 , p1 , p1 + p2 + p3 ) ) ) 
    
    ## R12 

    logger.info ( 'E_R12  : %s' , d.E_R12  ( s1 ,  s2 ) )
    logger.info ( 'E1_R12 : %s' , d.E1_R12 ( s1 ,  s2 ) )
    logger.info ( 'E2_R12 : %s' , d.E2_R12 ( s1 ,  s2 ) )
    logger.info ( 'E3_R12 : %s' , d.E3_R12 ( s1 ,  s2 ) )
    logger.info ( 'P_R12  : %s' , d.P_R12  ( s1 ,  s2 ) )
    logger.info ( 'P1_R12 : %s' , d.P1_R12 ( s1 ,  s2 ) )
    logger.info ( 'P2_R12 : %s' , d.P2_R12 ( s1 ,  s2 ) )
    logger.info ( 'P3_R12 : %s' , d.P3_R12 ( s1 ,  s2 ) )
    
    ct_31  = d.cos_31_R12  ( s1 , s2 )
    st2_31 = d.sin2_31_R12 ( s1 , s2 )
    logger.info ( 'cos_31 : %s [R=%s] [r=%s]' % ( ct_31 ,
                                                  ct_31 / CT ( p3 , p1 , p1 + p2 ) ,
                                                  ct_31 / ( 1 - st2_31 ) ** 0.5  ) )
    
    ## R23 

    logger.info ( 'E_R23  : %s' , d. E_R23 ( s1 ,  s2 ) )
    logger.info ( 'E1_R23 : %s' , d.E1_R23 ( s1 ,  s2 ) )
    logger.info ( 'E2_R23 : %s' , d.E2_R23 ( s1 ,  s2 ) )
    logger.info ( 'E3_R23 : %s' , d.E3_R23 ( s1 ,  s2 ) )
    logger.info ( 'P_R23  : %s' , d. P_R23 ( s1 ,  s2 ) )
    logger.info ( 'P1_R23 : %s' , d.P1_R23 ( s1 ,  s2 ) )
    logger.info ( 'P2_R23 : %s' , d.P2_R23 ( s1 ,  s2 ) )
    logger.info ( 'P3_R23 : %s' , d.P3_R23 ( s1 ,  s2 ) )
    
    ct_12  = d. cos_12_R23 ( s1 , s2 )
    st2_12 = d.sin2_12_R23 ( s1 , s2 )
    logger.info ( 'cos_12 : %s [R=%s] [r=%s]' % ( ct_12 ,
                                                  ct_12 / CT ( p1 , p2 , p2 + p3 ) ,
                                                  ct_12 / ( 1 - st2_12 ) ** 0.5  ) )

    ## R31 

    logger.info ( 'E_R31  : %s' , d. E_R31 ( s1 ,  s2 ) )
    logger.info ( 'E1_R31 : %s' , d.E1_R31 ( s1 ,  s2 ) )
    logger.info ( 'E2_R31 : %s' , d.E2_R31 ( s1 ,  s2 ) )
    logger.info ( 'E3_R31 : %s' , d.E3_R31 ( s1 ,  s2 ) )
    logger.info ( 'P_R31  : %s' , d. P_R31 ( s1 ,  s2 ) )
    logger.info ( 'P1_R31 : %s' , d.P1_R31 ( s1 ,  s2 ) )
    logger.info ( 'P2_R31 : %s' , d.P2_R31 ( s1 ,  s2 ) )
    logger.info ( 'P3_R31 : %s' , d.P3_R31 ( s1 ,  s2 ) )
    
    ct_23  = d. cos_23_R31 ( s1 , s2 )
    st2_23 = d.sin2_23_R31 ( s1 , s2 )
    logger.info ( 'cos_23 : %s [R=%s] [r=%s]' % ( ct_23 ,
                                                  ct_23 / CT ( p2 , p3 , p1 + p3 ) ,
                                                  ct_23 / ( 1 - st2_23 ) ** 0.5  ) )
    

    logger.info ( 'density: %s '  % d.density      ( s1      , s2      ) )
    logger.info ( 'density: %s '  % d.density_mass ( s1**0.5 , s2**0.5 ) )


def test_dalitz2 () :

    gr21  = d.graph21 () 
    gr21m = d.graph21 ( masses = True ) 

    gr31  = d.graph31 () 
    gr31m = d.graph31 ( masses = True ) 

    gr32  = d.graph32 () 
    gr32m = d.graph32 ( masses = True ) 

    gr21.draw  ( 'al' )
    gr31.draw  ( 'al' )
    gr32.draw  ( 'al' )

    gr21m.draw ( 'al' )
    gr31m.draw ( 'al' )
    gr32m.draw ( 'al' )
    
# =============================================================================
if '__main__' == __name__ :

    test_dalitz1 ()
    test_dalitz2 ()

# =============================================================================
##                                                                      The END 
# ============================================================================= 
