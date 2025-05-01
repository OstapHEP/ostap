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
from   ostap.core.core          import Ostap, hID 
from   ostap.utils.root_utils   import batch_env 
from   ostap.plotting.canvas    import use_canvas
from   ostap.plotting.style     import useStyle as use_style
from   ostap.utils.progress_bar import progress_bar 
import ostap.logger.table       as     T
import ostap.math.kinematic
import ostap.math.dalitz
import ostap.histos.histos 
import ROOT, random, time, math  
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.test_math_dalitz'   )
else                       : logger = getLogger ( __name__                   )
# ============================================================================= 
## set batch from environment 
batch_env ( logger )
# =============================================================================
#

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

graphs = set() 
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
    
    with use_canvas( 'test_dalitz2, squared masses (21)' , wait = 1 ) :
        gr21.draw  ( 'alcf' , linecolor = 2, linewidth = 2 , fillcolor = 2 )
    with use_canvas( 'test_dalitz2, squared masses (31)' , wait = 1 ) :
        gr31.draw  ( 'alcf' , linecolor = 2, linewidth = 2 , fillcolor = 2 )
    with use_canvas( 'test_dalitz2, squared masses (32)' , wait = 1 ) :
        gr32.draw  ( 'alcf' , linecolor = 2, linewidth = 2 , fillcolor = 2 )

    with use_canvas( 'test_dalitz2, masses (21)' , wait = 1 ) :
        gr21m.draw ( 'alcf' , linecolor = 2, linewidth = 2 , fillcolor = 2 )
    with use_canvas( 'test_dalitz2, masses (31)' , wait = 1 ) :
        gr31m.draw ( 'alcf' , linecolor = 2, linewidth = 2 , fillcolor = 2 )
    with use_canvas( 'test_dalitz2, masses (32)' , wait = 1 ) :
        gr32m.draw ( 'alcf' , linecolor = 2, linewidth = 2 , fillcolor = 2 )


# =============================================================================
def test_dalitz3 () :

    logger = getLogger  ('test_dalitz3' )

    for i , p  in enumerate ( plots ) : 
        
        with use_canvas( 'test_dalitz3 #%d'  % i , wait = 1 ) :
            
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
            
        with use_canvas( 'test_dalitz4 #%s'  % d , wait = 1 ) :            
            gr  = p.graph21 ( masses = False ) 
            h.colz ()
            gr.draw  ( 'l' , linecolor = 2, linewidth = 5 )            
            histos.append ( ( h , gr ) )
            

# =============================================================================
def test_dalitz5 () :
    
    logger = getLogger  ('test_dalitz5' )

    dd = Ostap.Kinematics.Dalitz ( 1 , 0.20 , 0.15 , 0.10  )

    NBX, NBY = 15, 15
    
    h1 = ROOT.TH2F ( hID() , 'Dalitz plot in s2,s1 varibales' ,
                     NBX , dd.s2_min () ,  dd.s2_max () , 
                     NBY , dd.s1_min () ,  dd.s1_max () )
    h2 = ROOT.TH2F ( hID() , 'Dalitz plot in x2,x1 variables' ,
                     NBX , dd.x2_min () ,  dd.x2_max () , 
                     NBY , dd.x1_min () ,  dd.x1_max () )
    h3 = ROOT.TH2F ( hID() , 'Dalitz plot in x2,x1 variables (weighted)' ,
                     NBX , dd.x2_min () ,  dd.x2_max () , 
                     NBY , dd.x1_min () ,  dd.x1_max () )
    
    ## s1,s2,s3 -> x1,x2

    N = 1000000
    for s1,s2,s3 in progress_bar ( dd.random ( N ) , max_value = N ) :

        h1.Fill ( s2 , s1 )
        
        x1,x2 = dd.s2x ( s1 , s2 )
        J     = dd.J   ( s1 , s2 )

        h2.Fill ( x2 , x1 )        
        h3.Fill ( x2 , x1 , 1./J )

    with use_style ( 'Z' ) : 
        with use_canvas( 'test_dalitz5 Dalitz(s2,s1)'          , wait = 2 ) : h1.colz ()
        with use_canvas( 'test_dalitz5 Dalitz(x2,x1)'          , wait = 2 ) : h2.colz ()
        with use_canvas( 'test_dalitz5 Dalitz(x2,x1) weighted' , wait = 2 ) : h3.colz ()
        
    graphs.add ( h1 )
    graphs.add ( h2 )
    graphs.add ( h3 )


# =============================================================================
def test_dalitz6 () :
    
    logger = getLogger  ('test_dalitz6' )

    dd = Ostap.Kinematics.Dalitz ( 1 , 0.20 , 0.15 , 0.10  )

    NBX, NBY = 20, 20
    
    h1 = ROOT.TH2F ( hID() , 'Dalittz plot in x2,x1 variables (flat)'    ,
                     NBX , dd.x2_min () ,  dd.x2_max () , 
                     NBY , dd.x1_min () ,  dd.x1_max () )
    h2 = ROOT.TH2F ( hID() , 'Dalitz plot in x2,x1 variables (weighted)' ,
                     NBX , dd.x2_min () ,  dd.x2_max () , 
                     NBY , dd.x1_min () ,  dd.x1_max () )
    h3 = ROOT.TH2F ( hID() , 'Dalitz plot in s2,s1 variable' ,
                     NBX , dd.s2_min () ,  dd.s2_max () , 
                     NBY , dd.s1_min () ,  dd.s1_max () )
    h4 = ROOT.TH2F ( hID() , 'Dalitz plot in s2,s1 variable (weighted)' ,
                     NBX , dd.s2_min () ,  dd.s2_max () , 
                     NBY , dd.s1_min () ,  dd.s1_max () )
    
    ## x1,x2 -> s1 , s2 

    N = 2000000
    for x1,x2,J in progress_bar ( dd.random_x ( N ) , max_value = N ) :

        h1.Fill ( x2 , x1     )
        h2.Fill ( x2 , x1 , J )

        s1 , s2 = dd.x2s ( x1 , x2 ) 

        h3.Fill ( s2 , s1 )
        h4.Fill ( s2 , s1 , J )

    with use_style ( 'Z' ) : 
        with use_canvas( 'test_dalitz6 Dalitz(x1,x1)'          , wait = 1 ) : h1.colz ()
        with use_canvas( 'test_dalitz6 Dalitz(x2,x1)'          , wait = 1 ) : h2.colz ()
        with use_canvas( 'test_dalitz6 Dalitz(s2,s1)'          , wait = 1 ) : h3.colz ()
        with use_canvas( 'test_dalitz6 Dalitz(s2,s1) weighted' , wait = 1 ) : h4.colz ()
        
    graphs.add ( h1 )
    graphs.add ( h2 )
    graphs.add ( h3 )
    graphs.add ( h4 )



# =============================================================================
def test_dalitz7 () :
    
    logger = getLogger  ('test_dalitz7' )

    dd = Ostap.Kinematics.Dalitz ( 1 , 0.20 , 0.15 , 0.10  )

    NBX, NBY = 15, 15
    
    h1 = ROOT.TH2F ( hID() , 'Dalitz plot in s2,s1 varibales' ,
                     NBX , dd.s2_min () ,  dd.s2_max () , 
                     NBY , dd.s1_min () ,  dd.s1_max () )
    h2 = ROOT.TH2F ( hID() , 'Dalitz plot in z2,z1 variables' ,
                     NBX , dd.z2_min () ,  dd.z2_max () , 
                     NBY , dd.z1_min () ,  dd.z1_max () )
    h3 = ROOT.TH2F ( hID() , 'Dalitz plot in x2,x1 variables (weighted)' ,
                     NBX , dd.z2_min () ,  dd.z2_max () , 
                     NBY , dd.z1_min () ,  dd.z1_max () )
    
    ## s1,s2,s3 -> z1,z2

    N = 1000000
    for s1,s2,s3 in progress_bar ( dd.random ( N ) , max_value = N ) :

        h1.Fill ( s2 , s1 )

        z1,z2 = dd.s2z ( s1 , s2 ) 
        Jz    = dd.Jz  ( s1 , s2 )

        h2.Fill ( z2 , z1 )             
        h3.Fill ( z2 , z1 , 1./Jz )
        

    with use_style ( 'Z' ) : 
        with use_canvas( 'test_dalitz7 Dalitz(s2,s1)'          , wait = 1 ) : h1.colz ()
        with use_canvas( 'test_dalitz7 Dalitz(z2,z1)'          , wait = 1 ) : h2.colz ()
        with use_canvas( 'test_dalitz7 Dalitz(z2,z1) weighted' , wait = 1 ) : h3.colz ()
        
    graphs.add ( h1 )
    graphs.add ( h2 )
    graphs.add ( h3 )

# =============================================================================
def test_dalitz8 () :
    
    logger = getLogger  ('test_dalitz8' )

    dd = Ostap.Kinematics.Dalitz ( 1 , 0.20 , 0.15 , 0.10  )

    NBX, NBY = 20, 20
    
    h1 = ROOT.TH2F ( hID() , 'Dalittz plot in z2,z1 variables (flat)'    ,
                     NBX , dd.z2_min () ,  dd.z2_max () , 
                     NBY , dd.z1_min () ,  dd.z1_max () )
    h2 = ROOT.TH2F ( hID() , 'Dalitz plot in z2,z1 variables (weighted)' ,
                     NBX , dd.z2_min () ,  dd.z2_max () , 
                     NBY , dd.z1_min () ,  dd.z1_max () )
    h3 = ROOT.TH2F ( hID() , 'Dalitz plot in s2,s1 variable' ,
                     NBX , dd.s2_min () ,  dd.s2_max () , 
                     NBY , dd.s1_min () ,  dd.s1_max () )
    h4 = ROOT.TH2F ( hID() , 'Dalitz plot in s2,s1 variable (weighted)' ,
                     NBX , dd.s2_min () ,  dd.s2_max () , 
                     NBY , dd.s1_min () ,  dd.s1_max () )
    
    ## x1,x2 -> s1 , s2 

    N = 2000000
    for z1,z2,Jz in progress_bar ( dd.random_z ( N ) , max_value = N ) :

        h1.Fill ( z2 , z1     )
        h2.Fill ( z2 , z1 , Jz )

        s1 , s2 = dd.z2s ( z1 , z2 ) 

        h3.Fill ( s2 , s1 )
        h4.Fill ( s2 , s1 , Jz )
        
    with use_style ( 'Z' ) : 
        with use_canvas( 'test_dalitz8 Dalitz(z1,z1)'          , wait = 1 ) : h1.colz ()
        with use_canvas( 'test_dalitz8 Dalitz(z2,z1)'          , wait = 1 ) : h2.colz ()
        with use_canvas( 'test_dalitz8 Dalitz(s2,s1)'          , wait = 1 ) : h3.colz ()
        with use_canvas( 'test_dalitz8 Dalitz(s2,s1) weighted' , wait = 1 ) : h4.colz ()
        
    graphs.add ( h1 )
    graphs.add ( h2 )
    graphs.add ( h3 )
    graphs.add ( h4 )


# =============================================================================
def test_dalitz9 () :
    
    logger = getLogger  ('test_dalitz9' )

    dd = Ostap.Kinematics.Dalitz ( 1 , 0.20 , 0.15 , 0.10  )

    NBX, NBY = 50, 50

    s1mn, s1mx = dd.s1_min() , dd.s1_max()
    s2mn, s2mx = dd.s2_min() , dd.s2_max()
    s3mn, s3mx = dd.s3_min() , dd.s3_max()

    d1 = s1mx - s1mn
    d2 = s2mx - s2mn
    d3 = s3mx - s3mn

    s1mn += 0.1 * d1
    s1mx -= 0.1 * d1
    s2mn += 0.1 * d2
    s2mx -= 0.1 * d2
    s3mn += 0.1 * d3
    s3mx -= 0.1 * d3
    
    regions = [
        lambda s1 , s2 :  s1 < s1mn ,
        lambda s1 , s2 :  s1 > s1mx ,
        lambda s1 , s2 :  s2 < s2mn ,
        lambda s1 , s2 :  s2 > s2mx ,
        lambda s1 , s2 :  s3 < s3mn ,
        lambda s1 , s2 :  s3 > s3mx ,
        ]
    
    histos = tuple ( 
        ( ROOT.TH2F ( hID() , 'Dalitz plot in s2,s1 variables' ,
                      NBX , dd.s2_min () ,  dd.s2_max () , 
                      NBY , dd.s1_min () ,  dd.s1_max () ) ,
          ROOT.TH2F ( hID() , 'Dalitz plot in z2,z1 variables (weighted)' ,
                      NBX , dd.z2_min () ,  dd.z2_max () , 
                      NBY , dd.z1_min () ,  dd.z1_max () )
          ) for i in regions        
        )

    colors = ( 2 , 4 , 6, 7, 8, 5 )
    for hh , c in zip ( histos , colors ) :
        for h in hh :
            h.SetLineColor (c)
            h.SetFillColor (c)

    ## s1,s2,s3 -> z1,z2
    N = 1000000
    for s1,s2,s3 in progress_bar ( dd.random ( N ) , max_value = N ) :
        
        z1,z2 = dd.s2z ( s1 , s2 ) 
        Jz    = dd.Jz  ( s1 , s2 )

        wz    = 1/Jz

        for r, hh in zip ( regions , histos ) :
            h1, h2 = hh 
            if r ( s1 , s2 ) :
                h1.Fill ( s2 , s1      )
                h2.Fill ( z2 , z1 , wz )
                

    with use_style ( 'Z' ) :

        with use_canvas( 'test_dalitz9 Dalitz(s2,s1)'          , wait = 1 ) :

            histos[0][0].draw('box') 
            for h1,h2 in histos[1:] : 
                h1.draw ( 'box same' )

            g1  = dd.graph21 ()
            g2  = g1.T()
            g2.SetLineWidth(2)
            g2.draw('c')
                
        with use_canvas( 'test_dalitz9 Dalitz(z2,z1)'          , wait = 1 ) :

            histos[0][1].draw('box') 
            for h1,h2 in histos[1:] : 
                h2.draw ( 'box same' )
                
    graphs.add ( histos )
    graphs.add ( g1     )
    graphs.add ( g2     )

    
# =============================================================================
def test_dalitz10 () :
    
    logger = getLogger  ('test_dalitz10' )

    dd = Ostap.Kinematics.Dalitz ( 1 , 0.20 , 0.15 , 0.10  )

    NBX, NBY = 50, 50

    s1mn, s1mx = dd.s1_min() , dd.s1_max()
    s2mn, s2mx = dd.s2_min() , dd.s2_max()
    s3mn, s3mx = dd.s3_min() , dd.s3_max()

    d1 = s1mx - s1mn
    d2 = s2mx - s2mn
    d3 = s3mx - s3mn

    m1 , g1 = math.sqrt ( s1mn + 0.3 * d1 ) , math.sqrt ( d1 ) / 20 
    m2 , g2 = math.sqrt ( s1mx - 0.3 * d1 ) , math.sqrt ( d1 ) / 20 
    m3 , g3 = math.sqrt ( s2mn + 0.3 * d2 ) , math.sqrt ( d2 ) / 20 
    m4 , g4 = math.sqrt ( s2mx - 0.3 * d2 ) , math.sqrt ( d2 ) / 20 
    m5 , g5 = math.sqrt ( s3mn + 0.3 * d3 ) , math.sqrt ( d3 ) / 20 
    m6 , g6 = math.sqrt ( s3mx - 0.3 * d3 ) , math.sqrt ( d3 ) / 20 

    regions = [
        lambda s1 , s2 , s3 :  abs (  math.sqrt ( s1 ) - m1 ) < g1 ,
        lambda s1 , s2 , s3 :  abs (  math.sqrt ( s1 ) - m2 ) < g2 ,
        lambda s1 , s2 , s3 :  abs (  math.sqrt ( s2 ) - m3 ) < g3 ,
        lambda s1 , s2 , s3 :  abs (  math.sqrt ( s2 ) - m4 ) < g4 ,
        lambda s1 , s2 , s3 :  abs (  math.sqrt ( s3 ) - m5 ) < g5 ,
        lambda s1 , s2 , s3 :  abs (  math.sqrt ( s3 ) - m6 ) < g6 ,
        ]
    
    histos = tuple ( 
        ( ROOT.TH2F ( hID() , 'Dalitz plot in s2,s1 variables' ,
                      NBX , dd.s2_min () ,  dd.s2_max () , 
                      NBY , dd.s1_min () ,  dd.s1_max () ) ,
          ROOT.TH2F ( hID() , 'Dalitz plot in z2,z1 variables (weighted)' ,
                      NBX , dd.z2_min () ,  dd.z2_max () , 
                      NBY , dd.z1_min () ,  dd.z1_max () )
          ) for i in regions        
        )

    colors = ( 2 , 4 , 6, 7, 8, 5 )
    for hh , c in zip ( histos , colors ) :
        for h in hh :
            h.SetLineColor (c)
            h.SetFillColor (c)

    ## s1,s2,s3 -> z1,z2
    N = 1000000
    for s1,s2,s3 in progress_bar ( dd.random ( N ) , max_value = N ) :
        
        z1,z2 = dd.s2z ( s1 , s2 ) 
        Jz    = dd.Jz  ( s1 , s2 )

        wz    = 1/Jz

        for r, hh in zip ( regions , histos ) :
            h1, h2 = hh 
            if r ( s1 , s2 , s3 ) :
                h1.Fill ( s2 , s1      )
                h2.Fill ( z2 , z1 , wz )
                

    with use_style ( 'Z' ) :

        with use_canvas( 'test_dalitz10  Dalitz(s2,s1)'          , wait = 1 ) :

            histos[0][0].draw('box') 
            for h1,h2 in histos[1:] : 
                h1.draw ( 'box same' )

            g1  = dd.graph21 ()
            g2  = g1.T()
            g2.SetLineWidth(2)
            g2.draw('c')
                
        with use_canvas( 'test_dalitz10 Dalitz(z2,z1)'          , wait = 1 ) :

            histos[0][1].draw('box') 
            for h1,h2 in histos[1:] : 
                h2.draw ( 'box same' )
                
    graphs.add ( histos )
    graphs.add ( g1     )
    graphs.add ( g2     )

# =============================================================================
if '__main__' == __name__ :
    
    test_dalitz1  ()
    test_dalitz2  ()
    test_dalitz3  ()
    test_dalitz4  ()

    test_dalitz5  ()
    test_dalitz6  ()    
    test_dalitz7  ()
    test_dalitz8  ()

    test_dalitz9  ()
    test_dalitz10 ()


# =============================================================================
##                                                                      The END 
# ============================================================================= 
