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
import ROOT, random, time  
from   ostap.core.core        import Ostap
import ostap.math.kinematic
import ostap.math.dalitz
from   ostap.utils.utils      import wait
from   ostap.plotting.canvas  import use_canvas
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.test_math_dalitz'   )
else                       : logger = getLogger ( __name__                   )
# ============================================================================= 

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


d0 = Ostap.Kinematics.Dalitz ( 1 , 0   , 0   , 0   )
d1 = Ostap.Kinematics.Dalitz ( 1 , 0.2 , 0   , 0   )
d2 = Ostap.Kinematics.Dalitz ( 1 , 0   , 0.2 , 0   )
d3 = Ostap.Kinematics.Dalitz ( 1 , 0   , 0   , 0.2 )
d4 = Ostap.Kinematics.Dalitz ( 1 , 0.2 , 0.2 , 0   )
d5 = Ostap.Kinematics.Dalitz ( 1 , 0.2 , 0   , 0.2 )
d6 = Ostap.Kinematics.Dalitz ( 1 , 0   , 0.2 , 0.2 )
d7 = Ostap.Kinematics.Dalitz ( 1 , 0.2 , 0.2 , 0.2 )

plots = d0 , d1 , d2 , d3 , d4 , d5 , d6 , d7 

# =============================================================================
def test_dalitz1 () : 

    logger = getLogger  ('test_dalitz1' ) 

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


# =============================================================================
def test_dalitz2 () :

    logger = getLogger  ('test_dalitz2' )

    gr21  = d.graph21 () 
    gr21m = d.graph21 ( masses = True ) 

    gr31  = d.graph31 () 
    gr31m = d.graph31 ( masses = True ) 

    gr32  = d.graph32 () 
    gr32m = d.graph32 ( masses = True ) 
    
    with wait ( 2 ) , use_canvas( 'test_dalitz2, squared masses (21)' ) :
        gr21.draw  ( 'alcf' , linecolor = 2, linewidth = 2 , fillcolor = 2 )
    with wait ( 2 ) , use_canvas( 'test_dalitz2, squared masses (31)' ) :
        gr31.draw  ( 'alcf' , linecolor = 2, linewidth = 2 , fillcolor = 2 )
    with wait ( 2 ) , use_canvas( 'test_dalitz2, squared masses (32)' ) :
        gr32.draw  ( 'alcf' , linecolor = 2, linewidth = 2 , fillcolor = 2 )

    with wait ( 2 ) , use_canvas( 'test_dalitz2, masses (21)' ) :
        gr21m.draw ( 'alcf' , linecolor = 2, linewidth = 2 , fillcolor = 2 )
    with wait ( 2 ) , use_canvas( 'test_dalitz2, masses (31)' ) :
        gr31m.draw ( 'alcf' , linecolor = 2, linewidth = 2 , fillcolor = 2 )
    with wait ( 2 ) , use_canvas( 'test_dalitz2, masses (32)' ) :
        gr32m.draw ( 'alcf' , linecolor = 2, linewidth = 2 , fillcolor = 2 )


# =============================================================================
def test_dalitz3 () :

    for i , p  in enumerate ( plots ) : 
        
        with wait ( 2 ) , use_canvas( 'test_dalitz3 #%d'  % i ) :
            
            if   i == 0 : logger.info ( "All masses are     zero" )
            elif i <  4 : logger.info ( "Two masses are     zero" )
            elif i <  7 : logger.info ( "One mass   is      zero" )
            else        : logger.info ( "All masses are non-zero" )
            
            gr  = p.graph21 ( masses = False ) 
            gr.draw  ( 'alf' , linecolor = 2, linewidth = 2 , fillcolor=2)
        
# =============================================================================
if '__main__' == __name__ :

    test_dalitz1 ()
    test_dalitz2 ()
    test_dalitz3 ()

# =============================================================================
##                                                                      The END 
# ============================================================================= 
