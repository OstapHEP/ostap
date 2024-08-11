#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file ostap/trees/tests/test_trees_param.py
# - It tests unbinned parameterization
# @see Ostap::Math::LegendreSum
# @see Ostap::Math::LegendreSum2
# @see Ostap::Math::LegendreSum3
# @see Ostap::Math::LegendreSum4
# @see Ostap::DataParam
# ============================================================================= 
""" Test unbinned parameterizations
- see Ostap::Math::LegendreSum
- see Ostap::Math::LegendreSum2
- see Ostap::Math::LegendreSum3
- see Ostap::Math::LegendreSum4
- see Ostap::DataParam
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ostap.core.pyrouts
import ostap.trees.param
import ostap.math.models
from   ostap.core.core          import hID, SE, Ostap
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait
from   ostap.utils.timing       import timing
from   ostap.utils.progress_bar import progress_bar, ProgressBar 
import ostap.logger.table       as     T 
import ROOT, os,  random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_trees_param' )
else : 
    logger = getLogger ( __name__           )
# =============================================================================
from ostap.utils.cleanup import CleanUp
data_file = CleanUp.tempfile ( suffix = '.root' , prefix = 'ostap-test-trees-param-' ) 

xmin, xmax = -2, 2 
ymin, ymax = -2, 2 
zmin, zmax = -4, 3 
umin, umax = -4, 6 

if not os.path.exists( data_file ) :
    
    import random
    
    N = 100000
    
    with timing ('Prepare input ROOT file with data %s' % data_file , logger = logger ) :
        with ROOT.TFile.Open( data_file ,'recreate') as test_file:
            ## test_file.cd()
            tree = ROOT.TTree('S','signal     tree')
            tree.SetDirectory ( test_file ) 
            
            from array import array 
            var1 = array ( 'd', [0] )
            var2 = array ( 'd', [0] )
            var3 = array ( 'd', [0] )
            var4 = array ( 'd', [0] )
            var5 = array ( 'd', [0] )
            
            tree .Branch ( 'x' , var1 , 'x/D' )
            tree .Branch ( 'y' , var2 , 'y/D' )
            tree .Branch ( 'z' , var3 , 'z/D' )
            tree .Branch ( 'u' , var4 , 'u/D' )
            tree .Branch ( 'v' , var5 , 'v/D' )
            
            I   = 0
            bar = ProgressBar ( N ) 
            while bar : 
                
                x =      random.uniform     ( xmin , xmax )
                y = -5 + random.expovariate ( 1/5.0 )
                z =      random.gauss       (   .0 , 4.0 )
                u =      random.gauss       (  2.0 , 5.0 )
                v =      random.gauss       (    0 , 1.0 )
            
                var1[0] =  x + 0.5 * y
                
                var2[0] =  x - 0.5 * y
                
                var3[0] = -x +       z + 0.5 * y + random.uniform ( -5 , 5 )
                
                var4[0] = -z + u - y   + 1.5 * x + random.uniform ( -5 , 5 )
                
                var5[0] =  v 

                if not xmin <= var1 [0] <= xmax : continue
                if not ymin <= var2 [0] <= ymax : continue
                if not zmin <= var3 [0] <= zmax : continue
                if not umin <= var4 [0] <= umax : continue

                bar += 1 
                tree.Fill()
                
            test_file.Write()
            del bar 
            test_file.ls   ()

def diff ( a , b ) :
    aa = float ( a  )
    bb = float ( b  )
    if aa == bb : return 0 
    return ( aa - bb ) / ( abs ( aa ) + abs ( bb ) ) 

# =============================================================================
## 1D parameterizations
# =============================================================================
def test_parameterize_1D () :

    logger  = getLogger("test_parameterize_1D")
    logger.info ( 'Test 1D parameterisations' ) 

    with ROOT.TFile.Open(data_file,'READ') as f :
        
        tree = f.S
        N    = 6
        
        with timing ( "Unbinned 1D-Legendre  parameterization" , logger = logger ) :
            
            lx = Ostap.Math.LegendreSum  ( N , xmin , xmax )
            ly = Ostap.Math.LegendreSum  ( N , ymin , ymax )
            lz = Ostap.Math.LegendreSum  ( N , zmin , zmax )
            lu = Ostap.Math.LegendreSum  ( N , umin , umax )
            
            lx.parameterize ( tree , 'x' )
            ly.parameterize ( tree , 'y' )
            lz.parameterize ( tree , 'z' )
            lu.parameterize ( tree , 'u' )
            
        with timing ( "Unbinned 1D-Chebyshev parameterization" , logger = logger ) :
            
            cx = Ostap.Math.ChebyshevSum ( N , xmin , xmax )
            cy = Ostap.Math.ChebyshevSum ( N , ymin , ymax )
            cz = Ostap.Math.ChebyshevSum ( N , zmin , zmax )
            cu = Ostap.Math.ChebyshevSum ( N , umin , umax )

            cx.parameterize ( tree , 'x' )
            cy.parameterize ( tree , 'y' )
            cz.parameterize ( tree , 'z' )
            cu.parameterize ( tree , 'u' )
            
        with timing ( "Unbinned 1D-Bernstein parameterization" , logger = logger ) :
            
            bx = Ostap.Math.Bernstein    ( N , xmin , xmax )
            by = Ostap.Math.Bernstein    ( N , ymin , ymax )
            bz = Ostap.Math.Bernstein    ( N , zmin , zmax )
            bu = Ostap.Math.Bernstein    ( N , umin , umax )
            
            bx.parameterize ( tree , 'x' )
            by.parameterize ( tree , 'y' )
            bz.parameterize ( tree , 'z' )
            bu.parameterize ( tree , 'u' )
            
        with timing ( "Prepare 1D-histos" , logger = logger ) :

            NX , NY , NZ , NU = 200 , 200 , 200 , 200 
            hx = ROOT.TH1D ( hID() , '' , NX , xmin , xmax )
            hy = ROOT.TH1D ( hID() , '' , NY , ymin , ymax )
            hz = ROOT.TH1D ( hID() , '' , NZ , zmin , zmax )
            hu = ROOT.TH1D ( hID() , '' , NU , umin , umax )
            
            tree.project ( hx , 'x' )
            tree.project ( hy , 'y' )
            tree.project ( hz , 'z' )
            tree.project ( hu , 'u' )

            for fx in ( lx, cx, bx ) : fx *= (xmax-xmin)/NX 
            for fy in ( ly, cy, by ) : fy *= (ymax-ymin)/NY
            for fz in ( lz, cz, bz ) : fz *= (zmax-zmin)/NZ
            for fu in ( lu, cu, bu ) : fu *= (umax-umin)/NU
                
            hx.SetMinimum(0)
            hy.SetMinimum(0)
            hz.SetMinimum(0)
            hu.SetMinimum(0)

        with use_canvas ( 'test_trees_param:X-variable' ) , wait ( 1 ) , timing ( 'draw X-variable' , logger = logger ) :
            
            with timing ( 'draw histogram' , logger= logger ) :
                hx.draw()        
            with timing ( 'draw Legendre'  , logger= logger ) :
                lx.draw('same', linecolor=2)
            with timing ( 'draw Chebyshev' , logger= logger ) :
                cx.draw('same', linecolor=4)
            with timing ( 'draw Bernstein' , logger= logger ) :
                bx.draw('same', linecolor=8)
            
        with use_canvas ( 'test_trees_param:Y-variable' ) , wait ( 1 ) , timing ( 'draw Y-variable' , logger = logger ) : 

            with timing ( 'draw histogram' , logger= logger ) :
                hy.draw()            
            with timing ( 'draw Legendre'  , logger= logger ) :
                ly.draw('same', linecolor=2)
            with timing ( 'draw Chebyshev' , logger= logger ) :
                cy.draw('same', linecolor=4)
            with timing ( 'draw Bernstein' , logger= logger ) :
                by.draw('same', linecolor=8)
            
        with use_canvas ( 'test_trees_param:Z-variable' ) , wait ( 1 ) , timing ( 'draw Z-variable' , logger = logger ) :

            with timing ( 'draw histogram' , logger= logger ) :
                hz.draw()
            with timing ( 'draw Legendre'  , logger= logger ) :
                lz.draw('same', linecolor=2)
            with timing ( 'draw Chebyshev' , logger= logger ) :
                cz.draw('same', linecolor=4)
            with timing ( 'draw Bernstein' , logger= logger ) :
                bz.draw('same', linecolor=8)
            
        with use_canvas ( 'test_trees_param:U-variable' ) , wait ( 1 ) , timing ( 'draw U-variable' , logger = logger ) :
            
            with timing ( 'draw histogram' , logger= logger ) :
                hu.draw()
            with timing ( 'draw Legendre'  , logger= logger ) :
                lu.draw('same', linecolor=2)
            with timing ( 'draw Chebyshev' , logger= logger ) :
                cu.draw('same', linecolor=4)
            with timing ( 'draw Bernstein' , logger= logger ) :
                bu.draw('same', linecolor=8)

        dhx, dhy, dhz, dhu = SE() , SE() , SE() , SE()
        dcx, dcy, dcz, dcu = SE() , SE() , SE() , SE()
        dbx, dby, dbz, dbu = SE() , SE() , SE() , SE()

        for i in progress_bar ( range ( 1000 ) ):
            
            x   = random.uniform ( xmin , xmax )
            y   = random.uniform ( ymin , ymax )
            z   = random.uniform ( zmin , zmax )
            u   = random.uniform ( umin , umax )

            dhx  += diff ( lx ( x ) , hx ( x ) )
            dhy  += diff ( ly ( y ) , hy ( y ) )
            dhz  += diff ( lz ( z ) , hz ( x ) )
            dhu  += diff ( lu ( u ) , hu ( u ) )

            dcx  += diff ( lx ( x ) , cx ( x ) )
            dcy  += diff ( ly ( y ) , cy ( y ) )
            dcz  += diff ( lz ( z ) , cz ( x ) )
            dcu  += diff ( lu ( u ) , cu ( u ) )

            dbx  += diff ( lx ( x ) , bx ( x ) )
            dby  += diff ( ly ( y ) , by ( y ) )
            dbz  += diff ( lz ( z ) , bz ( x ) )
            dbu  += diff ( lu ( u ) , bu ( u ) )

        rows = [ ( '' , 'Variable'  , 'mean[%]' , 'rms[%]' , 'min[%]' , 'max[%]') ]
        
        for  d, v in zip ( ( dhx , dhy , dhz , dhu ) , ( 'X' , 'Y' , 'Z' , 'U' ) ) :  
            row = 'Legendre vs histogram'  , v , \
                  '%+6.3f' % float ( d.mean() * 100 ) , \
                  '%6.3f'  % float ( d.rms () * 100 ) , \
                  '%+6.3f' % float ( d.min () * 100 ) , \
                  '%+6.3f' % float ( d.max () * 100 ) , 
            rows.append ( row )
            
        for  d, v in zip ( ( dcx , dcy , dcz , dcu ) , ( 'X' , 'Y' , 'Z' , 'U' ) ) :  
            row = 'Legendre vs Chebyshev'  , v , \
                  '%+6.3f' % float ( d.mean() * 100 ) , \
                  '%6.3f'  % float ( d.rms () * 100 ) , \
                  '%+6.3f' % float ( d.min () * 100 ) , \
                  '%+6.3f' % float ( d.max () * 100 ) , 
            rows.append ( row )

        for  d, v in zip ( ( dbx , dby , dbz , dbu ) , ( 'X' , 'Y' , 'Z' , 'U' ) ) :  
            row = 'Legendre vs Bernstein'  , v , \
                  '%+6.3f' % float ( d.mean() * 100 ) , \
                  '%6.3f'  % float ( d.rms () * 100 ) , \
                  '%+6.3f' % float ( d.min () * 100 ) , \
                  '%+6.3f' % float ( d.max () * 100 ) , 
            rows.append ( row )

        title = '1D parameterisation'
        table = T.table ( rows , title = title , prefix = '# ' , alignment =  'llcc' )
        logger.info ( '%s\n%s' % ( title , table ) )
        

# =============================================================================
## 2D parameterizations
# =============================================================================
def test_parameterize_2D () :
    
    logger  = getLogger("test_parameterize_2D")
    logger.info ( 'Test 2D parameterisations' ) 

    I1, I2 = 10 , 10
    with ROOT.TFile.Open(data_file,'READ') as f :
        
        tree = f.S

        with timing ( "2D-Legendre  prameterization" , logger = logger ) :
            
            lxy = Ostap.Math.LegendreSum2 ( I1 , I2 , xmin , xmax , ymin , ymax )
            lzu = Ostap.Math.LegendreSum2 ( I1 , I2 , zmin , zmax , umin , umax )
            lxu = Ostap.Math.LegendreSum2 ( I1 , I2 , xmin , xmax , umin , umax )
            
            lxy.parameterize ( tree , 'x' , 'y' )
            lzu.parameterize ( tree , 'z' , 'u' )
            lxu.parameterize ( tree , 'x' , 'u' )

        with timing ( "2D-Bernstein prameterization" , logger = logger ) :
            
            bxy = Ostap.Math.Bernstein2D ( I1 , I2 , xmin , xmax , ymin , ymax )
            bzu = Ostap.Math.Bernstein2D ( I1 , I2 , zmin , zmax , umin , umax )
            bxu = Ostap.Math.Bernstein2D ( I1 , I2 , xmin , xmax , umin , umax )
            
            bxy.parameterize ( tree , 'x' , 'y' )
            bzu.parameterize ( tree , 'z' , 'u' )
            bxu.parameterize ( tree , 'x' , 'u' )
                 
        with timing ( "2D-histograms" , logger = logger ) :

            N1 , N2 = 25 , 25 
            hxy = ROOT.TH2F(hID() , '', N1 , xmin , xmax , N2 , ymin , ymax )
            hzu = ROOT.TH2F(hID() , '', N1 , zmin , zmax , N2 , umin , umax )
            hxu = ROOT.TH2F(hID() , '', N1 , xmin , xmax , N2 , umin , umax )
            
            tree.project ( hxy , 'y:x' )
            tree.project ( hzu , 'u:z' )
            tree.project ( hxu , 'u:x' )

        lxy *= ( xmax - xmin ) * ( ymax - ymin ) / ( N1 * N2 ) 
        lzu *= ( zmax - zmin ) * ( umax - umin ) / ( N1 * N2 ) 
        lxu *= ( xmax - xmin ) * ( umax - umin ) / ( N1 * N2 ) 

        bxy *= ( xmax - xmin ) * ( ymax - ymin ) / ( N1 * N2 ) 
        bzu *= ( zmax - zmin ) * ( umax - umin ) / ( N1 * N2 ) 
        bxu *= ( xmax - xmin ) * ( umax - umin ) / ( N1 * N2 ) 

        d1h , d2h, d3h = SE() , SE() , SE()
        d1b , d2b, d3b = SE() , SE() , SE()
        
        for i in progress_bar ( range ( 10000 ) ) :
            
            x   = random.uniform ( xmin , xmax )
            y   = random.uniform ( ymin , ymax )
            z   = random.uniform ( zmin , zmax )
            u   = random.uniform ( umin , umax )

            d1h  += diff ( lxy ( x , y ) , hxy ( x , y ) )
            d2h  += diff ( lzu ( z , u ) , hzu ( z , u ) )
            d3h  += diff ( lxu ( x , u ) , hxu ( x , u ) )
            
            d1b  += diff ( lxy ( x , y ) , bxy ( x , y ) )
            d2b  += diff ( lzu ( z , u ) , bzu ( z , u ) )
            d3b  += diff ( lxu ( x , u ) , bxu ( x , u ) )


        rows = [ ( '' , 'Variable'  , 'mean[%]' , 'rms[%]' , 'min[%]' , 'max[%]') ]        
        for  d, v in zip ( ( d1h , d2h , d3h ) , ( 'X:Y' , 'Z:U' , 'X:U' ) ) :  
            row = 'Legendre vs histogram'  , v , \
                  '%+6.3f' % float ( d.mean() * 100 ) , \
                  '%6.3f'  % float ( d.rms () * 100 ) , \
                  '%+6.3f' % float ( d.min () * 100 ) , \
                  '%+6.3f' % float ( d.max () * 100 ) , 
            rows.append ( row )
            
        for  d, v in zip ( ( d1b , d2b , d3b ) , ( 'X:Y' , 'Z:U' , 'X:U' ) ) :  
            row = 'Legendre vs Bernstein'  , v , \
                  '%+6.3f' % float ( d.mean() * 100 ) , \
                  '%6.3f'  % float ( d.rms () * 100 ) , \
                  '%+6.3f' % float ( d.min () * 100 ) , \
                  '%+6.3f' % float ( d.max () * 100 ) , 
            rows.append ( row )
            
        title = '2D parameterisation'
        table = T.table ( rows , title = title , prefix = '# ' , alignment =  'llcc' )
        logger.info ( '%s\n%s' % ( title , table ) )
        
# =============================================================================
## 3D parameterizations
# =============================================================================
def test_parameterize_3D () :
    
    logger  = getLogger("test_parameterize_3D")
    logger.info ( 'Test 3D parameterisations' ) 

    with ROOT.TFile.Open(data_file,'READ') as f :
        
        tree = f.S
        
        I1 , I2 , I3  = 9 , 9 , 9 
        with timing ( "3D-Legendre  parameterization" , logger = logger ) :
            
            lx = Ostap.Math.LegendreSum3 ( I1 , I2 , I3 , ymin , ymax , zmin , zmax , umin , umax )
            ly = Ostap.Math.LegendreSum3 ( I1 , I2 , I3 , xmin , xmax , zmin , zmax , umin , umax )
            lz = Ostap.Math.LegendreSum3 ( I1 , I2 , I3 , xmin , xmax , ymin , ymax , umin , umax )
            lu = Ostap.Math.LegendreSum3 ( I1 , I2 , I3 , xmin , xmax , ymin , ymax , zmin , zmax )
            
            lx.parameterize ( tree , 'y' , 'z' , 'u' )
            ly.parameterize ( tree , 'x' , 'z' , 'u' )
            lz.parameterize ( tree , 'x' , 'y' , 'u' )
            lu.parameterize ( tree , 'x' , 'y' , 'z' )

        with timing ( "3D-Bernstein parameterization" , logger = logger ) :

            bx = Ostap.Math.Bernstein3D  ( I1 , I2 , I3 , ymin , ymax , zmin , zmax , umin , umax )
            by = Ostap.Math.Bernstein3D  ( I1 , I2 , I3 , xmin , xmax , zmin , zmax , umin , umax )
            bz = Ostap.Math.Bernstein3D  ( I1 , I2 , I3 , xmin , xmax , ymin , ymax , umin , umax )
            bu = Ostap.Math.Bernstein3D  ( I1 , I2 , I3 , xmin , xmax , ymin , ymax , zmin , zmax )
            
            bx.parameterize ( tree , 'y' , 'z' , 'u' )
            by.parameterize ( tree , 'x' , 'z' , 'u' )
            bz.parameterize ( tree , 'x' , 'y' , 'u' )
            bu.parameterize ( tree , 'x' , 'y' , 'z' )

        with timing ( "3D-histograms" , logger = logger ) :

            N1 , N2 , N3 = 30 , 30 , 30 
            hx = ROOT.TH3F ( hID() , '' , N1 , ymin , ymax , N2 , zmin, zmax , N3 , umin, umax )
            hy = ROOT.TH3F ( hID() , '' , N1 , xmin , xmax , N2 , zmin, zmax , N3 , umin, umax )
            hz = ROOT.TH3F ( hID() , '' , N1 , xmin , xmax , N2 , ymin, ymax , N3 , umin, umax )
            hu = ROOT.TH3F ( hID() , '' , N1 , xmin , xmax , N2 , ymin, ymax , N3 , zmin, zmax )
            
            tree.project ( hx , 'u:z:y' )
            tree.project ( hy , 'u:z:x' )
            tree.project ( hz , 'u:y:x' )
            tree.project ( hu , 'z:y:x' )
            
            lx *= ( ymax - ymin ) * ( zmax - zmin ) * ( umax - umin ) / ( N1 * N2 * N3 )
            ly *= ( xmax - xmin ) * ( zmax - zmin ) * ( umax - umin ) / ( N1 * N2 * N3 )
            lz *= ( xmax - xmin ) * ( ymax - ymin ) * ( umax - umin ) / ( N1 * N2 * N3 )
            lu *= ( xmax - xmin ) * ( ymax - ymin ) * ( zmax - zmin ) / ( N1 * N2 * N3 )

            bx *= ( ymax - ymin ) * ( zmax - zmin ) * ( umax - umin ) / ( N1 * N2 * N3 )
            by *= ( xmax - xmin ) * ( zmax - zmin ) * ( umax - umin ) / ( N1 * N2 * N3 )
            bz *= ( xmax - xmin ) * ( ymax - ymin ) * ( umax - umin ) / ( N1 * N2 * N3 )
            bu *= ( xmax - xmin ) * ( ymax - ymin ) * ( zmax - zmin ) / ( N1 * N2 * N3 )

        dhx, dhy, dhz, dhu = SE() , SE() , SE() , SE()
        dbx, dby, dbz, dbu = SE() , SE() , SE() , SE()

        for i in progress_bar ( range ( 100000 ) ) :
            
            x   = random.uniform ( xmin , xmax )
            y   = random.uniform ( ymin , ymax )
            z   = random.uniform ( zmin , zmax )
            u   = random.uniform ( umin , umax )

            dhx += diff ( lx ( y , z , u  ) , hx ( y , z , u ) )
            dhy += diff ( ly ( x , z , u  ) , hy ( x , z , u ) )
            dhz += diff ( lz ( x , y , u  ) , hz ( x , y , u ) )
            dhu += diff ( lu ( x , y , z  ) , hu ( x , y , z ) )

            dbx += diff ( lx ( y , z , u  ) , bx ( y , z , u ) )
            dby += diff ( ly ( x , z , u  ) , by ( x , z , u ) )
            dbz += diff ( lz ( x , y , u  ) , bz ( x , y , u ) )
            dbu += diff ( lu ( x , y , z  ) , bu ( x , y , z ) )

        rows = [ ( '' , 'Variable'  , 'mean[%]' , 'rms[%]' , 'min[%]' , 'max[%]') ]
        for  d, v in zip ( ( dhx , dhy , dhz, dhu ) , ( 'Y:Z:U' , 'X:Z:U' , 'X:Y:U'  , 'X:Y:Z') ) :  
            row = 'Legendre vs histogram'  , v , \
                  '%+6.3f' % float ( d.mean() * 100 ) , \
                  '%6.3f'  % float ( d.rms () * 100 ) , \
                  '%+6.3f' % float ( d.min () * 100 ) , \
                  '%+6.3f' % float ( d.max () * 100 ) , 
            rows.append ( row )
            
        for  d, v in zip ( ( dbx , dby , dbz, dbu ) , ( 'Y:Z:U' , 'X:Z:U' , 'X:Y:U'  , 'X:Y:Z') ) :  
            row = 'Legendre vs Bernstein'  , v , \
                  '%+6.3f' % float ( d.mean() * 100 ) , \
                  '%6.3f'  % float ( d.rms () * 100 ) , \
                  '%+6.3f' % float ( d.min () * 100 ) , \
                  '%+6.3f' % float ( d.max () * 100 ) , 
            rows.append ( row )
            
        title = '3D parameterisation'
        table = T.table ( rows , title = title , prefix = '# ' , alignment =  'llcc' )
        logger.info ( '%s\n%s' % ( title , table ) )
        


# =============================================================================
## 4D parameterizations
# =============================================================================
def test_parameterize_4D () :
    
    logger  = getLogger("test_parameterize_4D")
    logger.info ( 'Test 4D parameterisations' ) 

    with ROOT.TFile.Open(data_file,'READ') as f :
        
        tree = f.S
        
        l  = Ostap.Math.LegendreSum4 ( 10 , 10 , 10 , 10 ,
                                       xmin , xmax ,
                                       ymin , ymax ,
                                       zmin , zmax ,
                                       umin , umax ) 
        
        with timing ( "4D-Legendre  parameterization" , logger = logger ) :
            l.parameterize ( tree , 'x' , 'y' , 'z' , 'u' )
            
        lxy = l.integralU().integralZ() 
        lzu = l.integralX().integralX() 
        lxu = l.integralY().integralY() 
        
        hxy = ROOT.TH2F(hID() , '', 20 , xmin , xmax , 20 , ymin , ymax )
        hzu = ROOT.TH2F(hID() , '', 20 , zmin , zmax , 20 , umin , umax )
        hxu = ROOT.TH2F(hID() , '', 20 , xmin , xmax , 20 , umin , umax )
        
        with timing ( "Histogram projections" , logger = logger ) :
            tree.project ( hxy , 'y:x' )
            tree.project ( hzu , 'u:z' )
            tree.project ( hxu , 'u:x' )
        
        lxy *= ( ( xmax - xmin ) / 20 ) * ( ( ymax - ymin ) / 20 ) 
        lzu *= ( ( zmax - zmin ) / 20 ) * ( ( umax - umin ) / 20 ) 
        lxu *= ( ( xmax - xmin ) / 20 ) * ( ( umax - umin ) / 20 ) 
        
        d1 = SE()
        d2 = SE()
        d3 = SE() 
        
        for i in progress_bar ( 100000 ) :
            
            x   =  random.uniform ( xmin , xmax )
            y   =  random.uniform ( ymin , ymax )
            z   =  random.uniform ( zmin , zmax )
            u   =  random.uniform ( umin , umax )
        
            d1 += diff ( lxy ( x , y ) , hxy ( x , y ) )
            d2 += diff ( lzu ( z , u ) , hzu ( z , u ) ) 
            d3 += diff ( lxu ( x , u ) , hxu ( x , u ) ) 

            
        rows = [ ( '' , 'Variable'  , 'mean[%]' , 'rms[%]' , 'min[%]' , 'max[%]') ]
        for  d, v in zip ( ( d1 , d2 , d3 ) , ( 'X:Y' , 'Z:U' , 'X:U' ) ) :  
            row = 'Legendre vs histogram'  , v , \
                  '%+6.3f' % float ( d.mean() * 100 ) , \
                  '%6.3f'  % float ( d.rms () * 100 ) , \
                  '%+6.3f' % float ( d.min () * 100 ) , \
                  '%+6.3f' % float ( d.max () * 100 ) , 
            rows.append ( row )
            
        title = '4D parameterisation+projections'
        table = T.table ( rows , title = title , prefix = '# ' , alignment =  'llcc' )
        logger.info ( '%s\n%s' % ( title , table ) )
        
# =============================================================================

# =============================================================================
## 1D statistics 
# =============================================================================
def test_statistics_1D () :
    
    logger  = getLogger("test_statistics_1D")
    logger.info ( 'Test 1D statistics' ) 
    
    with ROOT.TFile.Open(data_file,'READ') as f :
        
        tree = f.S

        rows = [  ('Parameter' , 'Value' ) ]

        for i in range ( 5 ) :            
            v = tree.get_moment ( i , 0.0 , 'v' )
            row = "moment[%d,0.0,'v']" % i , '%+.4f' % v 
            rows.append ( row )

        for i in range ( 5 ) :            
            v = tree.moment ( i , 'v' )
            row = "moment[%d,'v']" % i , v.toString ( '%+.4f +/- %-.4f' ) 
            rows.append ( row )

        for i in range ( 5 ) :            
            v = tree.central_moment ( i , 'v' )
            row = "central moment[%d,'v']" % i , v.toString ( '%+.4f +/- %-.4f' ) 
            rows.append ( row )

            
        vv  = 'abs(1+0.01*x)'
        
        v   = tree.harmonic_mean ( vv )
        row = 'harmonic   mean' , '%+.6f' % v.value()
        rows.append ( row )
        
        v   = tree.geometric_mean ( vv )
        row = 'geometric  mean' , '%+.6f' % v.value() 
        rows.append ( row )
        
        v   = tree.arithmetic_mean ( vv )             
        row = 'arithmetic mean' , '%+.6f' % v.value() 
        rows.append ( row )

        for p in range ( -2 , 6 ) :
            v   = tree.power_mean ( p , vv  )
            row = 'power [%+d] mean' % p , '%+.6f' % v.value() 
            rows.append ( row )

        for p in range ( -2 , 6 ) :
            v   = tree.lehmer_mean ( p , vv  )
            row = 'Lehmer[%+d] mean' % p , '%+.6f' % v.value() 
            rows.append ( row )
            
        
            
        title = '1D statistics'
        table = T.table ( rows , title = title , prefix = '# ' , alignment =  'lc' )
        logger.info ( '%s\n%s' % ( title , table ) )

        print ( tree.statVar ( 'v' ) ) 
        

# =============================================================================
if '__main__' == __name__ :

    test_parameterize_1D () 
    test_parameterize_2D () 
    test_parameterize_3D () 
    test_parameterize_4D () 
    
    test_statistics_1D () 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
