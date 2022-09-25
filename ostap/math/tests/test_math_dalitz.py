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
from   ostap.core.core        import Ostap, hID 
import ostap.math.kinematic
import ostap.math.dalitz
import ostap.histos.histos 
import ostap.logger.table     as     T
from   ostap.utils.utils      import wait
from   ostap.plotting.canvas  import use_canvas
import ROOT, random, time  
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

    rows = [ ( 'Quantity' , 'Value' ) ]

    rows.append ( ( 'E1'            , '%+.7g' % d.E1     ( s1 ,  s2 ) ) )
    rows.append ( ( 'E2'            , '%+.7g' % d.E2     ( s1 ,  s2 ) ) )
    rows.append ( ( 'E3'            , '%+.7g' % d.E3     ( s1 ,  s2 ) ) )
    
    rows.append ( ( 'P1'            , '%+.7g' % d.P1     ( s1 ,  s2 ) ) )
    rows.append ( ( 'P2'            , '%+.7g' % d.P2     ( s1 ,  s2 ) ) )
    rows.append ( ( 'P3'            , '%+.7g' % d.P3     ( s1 ,  s2 ) ) )


    rows.append ( ( 'cos_12'        , '%+.7g' % d.cos_12 ( s1 ,  s2 ) ) )
    rows.append ( ( 'cos_23'        , '%+.7g' % d.cos_23 ( s1 ,  s2 ) ) )
    rows.append ( ( 'cos_31'        , '%+.7g' % d.cos_31 ( s1 ,  s2 ) ) )

    rows.append ( ( 'R:cos_12'      , '%+.7g' % ( d.cos_12 ( s1 ,  s2 ) / CT ( p1 , p2 , p1 + p2 + p3 ) ) ) ) 
    rows.append ( ( 'R:cos_23'      , '%+.7g' % ( d.cos_23 ( s1 ,  s2 ) / CT ( p2 , p3 , p1 + p2 + p3 ) ) ) )
    rows.append ( ( 'R:cos_31'      , '%+.7g' % ( d.cos_31 ( s1 ,  s2 ) / CT ( p3 , p1 , p1 + p2 + p3 ) ) ) )
    
    ## R12
    
    rows.append ( ( 'E_R12'         ,  '%+.7g'  % d.E_R12         ( s1 ,  s2 ) ) )
    rows.append ( ( 'E1_R12'        ,  '%+.7g'  % d.E1_R12        ( s1 ,  s2 ) ) )
    rows.append ( ( 'E2_R12'        ,  '%+.7g'  % d.E2_R12        ( s1 ,  s2 ) ) )
    rows.append ( ( 'E3_R12'        ,  '%+.7g'  % d.E3_R12        ( s1 ,  s2 ) ) )

    rows.append ( ( 'P_R12'         ,  '%+.7g'  % d.P_R12         ( s1 ,  s2 ) ) )
    rows.append ( ( 'P1_R12'        ,  '%+.7g'  % d.P1_R12        ( s1 ,  s2 ) ) )
    rows.append ( ( 'P2_R12'        ,  '%+.7g'  % d.P2_R12        ( s1 ,  s2 ) ) )
    rows.append ( ( 'P3_R12'        ,  '%+.7g'  % d.P3_R12        ( s1 ,  s2 ) ) )
 
    rows.append ( ( 'cos_31_R12'    ,  '%+.7g'  % d.cos_31_R12    ( s1 ,  s2 ) ) )
    rows.append ( ( 'sin2_31_R12'   ,  '%+.7g'  % d.sin2_31_R12   ( s1 ,  s2 ) ) )
    rows.append ( ( 'R1:cos_31_R12' ,  '%+.7g'  % ( d.cos_31_R12  ( s1 ,  s2 ) / CT ( p3 , p1 , p1 + p2 ) ) ) )
    rows.append ( ( 'R2:cos_31_R12' ,  '%+.7g'  % ( d.cos_31_R12  ( s1 ,  s2 ) / ( 1 - d.sin2_31_R12 ( s1 ,  s2 ) ) ** 0.5 ) ) ) 

    ## R23
    
    rows.append ( ( 'E_R23'         ,  '%+.7g'  % d.E_R23         ( s1 ,  s2 ) ) )
    rows.append ( ( 'E1_R23'        ,  '%+.7g'  % d.E1_R23        ( s1 ,  s2 ) ) )
    rows.append ( ( 'E2_R23'        ,  '%+.7g'  % d.E2_R23        ( s1 ,  s2 ) ) )
    rows.append ( ( 'E3_R23'        ,  '%+.7g'  % d.E3_R23        ( s1 ,  s2 ) ) )

    rows.append ( ( 'P_R23'         ,  '%+.7g'  % d.P_R23         ( s1 ,  s2 ) ) )
    rows.append ( ( 'P1_R23'        ,  '%+.7g'  % d.P1_R23        ( s1 ,  s2 ) ) )
    rows.append ( ( 'P2_R23'        ,  '%+.7g'  % d.P2_R23        ( s1 ,  s2 ) ) )
    rows.append ( ( 'P3_R23'        ,  '%+.7g'  % d.P3_R23        ( s1 ,  s2 ) ) )
 
    rows.append ( ( 'cos_12_R23'    ,  '%+.7g'  % d.cos_12_R23    ( s1 ,  s2 ) ) )
    rows.append ( ( 'sin2_12_R23'   ,  '%+.7g'  % d.sin2_12_R23   ( s1 ,  s2 ) ) )
    rows.append ( ( 'R1:cos_12_R23' ,  '%+.7g'  % ( d.cos_12_R23  ( s1 ,  s2 ) / CT ( p1 , p2 , p2 + p3 ) ) ) )
    rows.append ( ( 'R2:cos_12_R23' ,  '%+.7g'  % ( d.cos_12_R23  ( s1 ,  s2 ) / ( 1 - d.sin2_12_R23 ( s1 ,  s2 ) ) ** 0.5 ) ) ) 

    ## R31
    
    rows.append ( ( 'E_R31'         ,  '%+.7g'  % d.E_R31         ( s1 ,  s2 ) ) )
    rows.append ( ( 'E1_R31'        ,  '%+.7g'  % d.E1_R31        ( s1 ,  s2 ) ) )
    rows.append ( ( 'E2_R31'        ,  '%+.7g'  % d.E2_R31        ( s1 ,  s2 ) ) )
    rows.append ( ( 'E3_R31'        ,  '%+.7g'  % d.E3_R31        ( s1 ,  s2 ) ) )

    rows.append ( ( 'P_R31'         ,  '%+.7g'  % d.P_R31         ( s1 ,  s2 ) ) )
    rows.append ( ( 'P1_R31'        ,  '%+.7g'  % d.P1_R31        ( s1 ,  s2 ) ) )
    rows.append ( ( 'P2_R31'        ,  '%+.7g'  % d.P2_R31        ( s1 ,  s2 ) ) )
    rows.append ( ( 'P3_R32'        ,  '%+.7g'  % d.P3_R31        ( s1 ,  s2 ) ) )
 
    rows.append ( ( 'cos_23_R31'    ,  '%+.7g'  % d.cos_23_R31    ( s1 ,  s2 ) ) )
    rows.append ( ( 'sin2_23_R31'   ,  '%+.7g'  % d.sin2_23_R31   ( s1 ,  s2 ) ) )
    rows.append ( ( 'R1:cos_23_R31' ,  '%+.7g'  % ( d.cos_23_R31  ( s1 ,  s2 ) / CT ( p2 , p3 , p3 + p1 ) ) ) )
    rows.append ( ( 'R2:cos_23_R31' ,  '%+.7g'  % ( d.cos_23_R31  ( s1 ,  s2 ) / ( 1 - d.sin2_23_R31 ( s1 ,  s2 ) ) ** 0.5 ) ) ) 


    ## Misha's formalism
    
    rows.append ( ( 'cos_theta12'   ,  '%+.7g'  % d.cos_theta12   ( s1 ,  s2 ) ) )
    rows.append ( ( 'cos_theta23'   ,  '%+.7g'  % d.cos_theta23   ( s1 ,  s2 ) ) )
    rows.append ( ( 'cos_theta31'   ,  '%+.7g'  % d.cos_theta31   ( s1 ,  s2 ) ) )


    rows.append ( ( 'cos_zeta120'   ,  '%+.7g'  % d.cos_zeta120   ( s1 ,  s2 ) ) )
    rows.append ( ( 'cos_zeta230'   ,  '%+.7g'  % d.cos_zeta230   ( s1 ,  s2 ) ) )
    rows.append ( ( 'cos_zeta310'   ,  '%+.7g'  % d.cos_zeta310   ( s1 ,  s2 ) ) )


    rows.append ( ( 'cos_zeta131'   ,  '%+.7g'  % d.cos_zeta131   ( s1 ,  s2 ) ) )
    rows.append ( ( 'cos_zeta211'   ,  '%+.7g'  % d.cos_zeta211   ( s1 ,  s2 ) ) )
    rows.append ( ( 'cos_zeta212'   ,  '%+.7g'  % d.cos_zeta212   ( s1 ,  s2 ) ) )
    rows.append ( ( 'cos_zeta322'   ,  '%+.7g'  % d.cos_zeta322   ( s1 ,  s2 ) ) )
    rows.append ( ( 'cos_zeta323'   ,  '%+.7g'  % d.cos_zeta323   ( s1 ,  s2 ) ) )
    rows.append ( ( 'cos_zeta133'   ,  '%+.7g'  % d.cos_zeta133   ( s1 ,  s2 ) ) )
    rows.append ( ( 'cos_zeta231'   ,  '%+.7g'  % d.cos_zeta231   ( s1 ,  s2 ) ) )
    rows.append ( ( 'cos_zeta312'   ,  '%+.7g'  % d.cos_zeta312   ( s1 ,  s2 ) ) )
    rows.append ( ( 'cos_zeta123'   ,  '%+.7g'  % d.cos_zeta123   ( s1 ,  s2 ) ) )
    

    rows.append ( ( 'density'       , '%+.7g'  %  d.density      ( s1        , s2        ) ) ) 
    rows.append ( ( 'density mass'  , '%+.7g'  %  d.density      ( s1 ** 0.5 , s2 ** 0.5 ) ) ) 
    


    title = '%s' % d 
    table = T.table ( rows , title = title , alignment = 'll' , prefix = '# ')
    logger.info ( 'Test %s\n%s' % ( title , table ) )

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

    logger = getLogger  ('test_dalitz3' )

    for i , p  in enumerate ( plots ) : 
        
        with wait ( 2 ) , use_canvas( 'test_dalitz3 #%d'  % i ) :
            
            if   i == 0 : logger.info ( "All masses are     zero" )
            elif i <  4 : logger.info ( "Two masses are     zero" )
            elif i <  7 : logger.info ( "One mass   is      zero" )
            else        : logger.info ( "All masses are non-zero" )
            
            gr  = p.graph21 ( masses = False ) 
            gr.draw  ( 'alf' , linecolor = 2, linewidth = 2 , fillcolor=2)


# =============================================================================
def test_dalitz4 () :

    logger = getLogger  ('test_dalitz4' )


    histos = [] 
    for p in plots :
        
        h = ROOT.TH2D ( hID() , 'Dalitz plot %s' % p ,
                        25 , 0 , 1 ,
                        25 , 0 , 1 )
        
        for s1,s2,s3  in p.random ( 10000 ) :
            h.Fill ( s1 , s2 )
            
        with wait ( 2 ) , use_canvas( 'test_dalitz4 #%s'  % d ) :            
            gr  = p.graph21 ( masses = False ) 
            h.colz ()
            gr.draw  ( 'l' , linecolor = 2, linewidth = 5 )            
            histos.append ( ( h , gr ) )
            
            
# =============================================================================
if '__main__' == __name__ :
    
    test_dalitz1 ()
    test_dalitz2 ()
    test_dalitz3 ()
    test_dalitz4 ()

# =============================================================================
##                                                                      The END 
# ============================================================================= 
