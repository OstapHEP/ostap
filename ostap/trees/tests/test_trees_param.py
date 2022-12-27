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
from   ostap.utils.progress_bar import progress_bar 
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

if not os.path.exists( data_file ) :
    
    import random
    
    N = 100000
        
    logger.info('Prepare input ROOT file with data %s' % data_file )
    with ROOT.TFile.Open( data_file ,'recreate') as test_file:
        ## test_file.cd()
        tree = ROOT.TTree('S','signal     tree')
        tree.SetDirectory ( test_file ) 
        
        from array import array 
        var1 = array ( 'd', [0] )
        var2 = array ( 'd', [0] )
        var3 = array ( 'd', [0] )
        var4 = array ( 'd', [0] )
        
        tree .Branch ( 'x' , var1 , 'x/D' )
        tree .Branch ( 'y' , var2 , 'y/D' )
        tree .Branch ( 'z' , var3 , 'z/D' )
        tree .Branch ( 'u' , var4 , 'u/D' )

        for i in progress_bar ( range  ( N ) ) : 
            
            x =      random.uniform     ( -4.0 , 4.0 )
            y = -5 + random.expovariate ( 1/5.0 )
            z =      random.gauss       (   .0 , 4.0 )
            u =      random.gauss       (  2.0 , 5.0 )
            
            var1[0] =  x + 0.5 * y  
            var2[0] =  x - 0.5 * y  
            var3[0] = -x +       z + 0.5 * y + random.uniform ( -5 , 5 ) 
            var4[0] = -z + u - y   + 1.5 * x + random.uniform ( -5 , 5 ) 

            if not -2 <= var1[0] <= 2 : continue
            if not -2 <= var2[0] <= 2 : continue
            if not -4 <= var3[0] <= 4 : continue
            if not -4 <= var4[0] <= 6 : continue
                        
            tree.Fill()
            
        test_file.Write()
        test_file.ls()

cut1 = ROOT.TCut ( '-2<=x&&x<=2' )
cut2 = ROOT.TCut ( '-2<=y&&y<=2' )
cut3 = ROOT.TCut ( '-4<=z&&z<=4' )
cut4 = ROOT.TCut ( '-4<=u&&u<=6' )
cuts = cut1&cut2&cut4&cut4

# =============================================================================
## 1D parameterizations
# =============================================================================
def test_parameterize_1D () :

    logger  = getLogger("test_parameterize_1D")
    logger.info ( 'Test 1D parameterisations' ) 

    with ROOT.TFile.Open(data_file,'READ') as f :
        
        tree = f.S
        
        lx = Ostap.Math.LegendreSum ( 4 , -2 , 2 )
        ly = Ostap.Math.LegendreSum ( 4 , -2 , 2 )
        lz = Ostap.Math.LegendreSum ( 4 , -4 , 4 )
        lu = Ostap.Math.LegendreSum ( 4 , -4 , 6 )

        with timing ( "LegendreSums" , logger = logger ) : 
            lx.parameterize ( tree , 'x' , cuts )
            ly.parameterize ( tree , 'y' , cuts )
            lz.parameterize ( tree , 'z' , cuts )
            lu.parameterize ( tree , 'u' , cuts )
            
        cx = Ostap.Math.ChebyshevSum ( 4 , -2 , 2 )
        cy = Ostap.Math.ChebyshevSum ( 4 , -2 , 2 )
        cz = Ostap.Math.ChebyshevSum ( 4 , -4 , 4 )
        cu = Ostap.Math.ChebyshevSum ( 4 , -4 , 6 )

        with timing ( "ChebyshevSums" , logger = logger ) : 
            cx.parameterize ( tree , 'x' , cuts )
            cy.parameterize ( tree , 'y' , cuts )
            cz.parameterize ( tree , 'z' , cuts )
            cu.parameterize ( tree , 'u' , cuts )
            
        bx = Ostap.Math.Bernstein  ( 4 , -2 , 2 )
        by = Ostap.Math.Bernstein  ( 4 , -2 , 2 )
        bz = Ostap.Math.Bernstein  ( 4 , -4 , 4 )
        bu = Ostap.Math.Bernstein  ( 4 , -4 , 6 )

        with timing ( "Bernstein" , logger = logger ) : 
            bx.parameterize ( tree , 'x' , cuts )
            by.parameterize ( tree , 'y' , cuts )
            bz.parameterize ( tree , 'z' , cuts )
            bu.parameterize ( tree , 'u' , cuts )
            
        hx = ROOT.TH1D(hID(),'',100,-2,2)
        hy = ROOT.TH1D(hID(),'',100,-2,2)
        hz = ROOT.TH1D(hID(),'',100,-4,4)
        hu = ROOT.TH1D(hID(),'',100,-4,6)
        
        with timing ( "Histos" , logger = logger ) : 
            tree.project ( hx , 'x' , cuts )
            tree.project ( hy , 'y' , cuts )
            tree.project ( hz , 'z' , cuts )
            tree.project ( hu , 'u' , cuts )
            
        hx.SetMinimum(0)
        hy.SetMinimum(0)
        hz.SetMinimum(0)
        hu.SetMinimum(0)

        with use_canvas ( 'X-variable' ) , wait ( 1 ) , timing ( 'X-varibale' , logger = logger ) : 
            hx.draw()        
            lx *= 0.04   ## bin-width
            lx.draw('same', linecolor=2)
            cx *= 0.04   ## bin-width
            cx.draw('same', linecolor=4)
            bx *= 0.04   ## bin-width
            bx.draw('same', linecolor=8)
            
        with use_canvas ( 'Y-variable' ) , wait ( 1 ) , timing ( 'Y-variable' , logger = logger ) : 
            hy.draw()            
            ly *= 0.04   ## bin-width
            ly.draw('same', linecolor=2)
            cy *= 0.04   ## bin-width
            cy.draw('same', linecolor=4)
            by *= 0.04   ## bin-width
            by.draw('same', linecolor=8)
            
        with use_canvas ( 'Z-variable' ) , wait ( 1 ) , timing ( 'Z-variable' , logger = logger ) : 
            hz.draw()
            lz *= 0.08   ## bin-width
            lz.draw('same', linecolor=2)
            cz *= 0.08   ## bin-width
            cz.draw('same', linecolor=4)
            bz *= 0.08   ## bin-width
            bz.draw('same', linecolor=8)
            
        with use_canvas ( 'U-variable' ) , wait ( 1 ) , timing ( 'U-variable' , logger = logger ) :
            hu.draw()
            lu *= 0.10   ## bin-width
            lu.draw('same', linecolor=2)
            cu *= 0.10   ## bin-width
            cu.draw('same', linecolor=4)
            bu *= 0.10   ## bin-width
            bu.draw('same', linecolor=8)

        d1  = SE()
        d2  = SE()
        d3  = SE()
        d4  = SE()
        dp1 = SE()
        dp2 = SE()
        dp3 = SE()
        dp4 = SE()
        db1 = SE()
        db2 = SE()
        db3 = SE()
        db4 = SE()
        
        for i in progress_bar ( range ( 1000 ) ):
            
            x   = random.uniform ( -2 , 2 )
            y   = random.uniform ( -2 , 2 )
            z   = random.uniform ( -4 , 4 )
            u   = random.uniform ( -4 , 6 )
            
            d1  += ( hx ( x ) - lx ( x ) ) / max ( hx ( x ) , lx ( x ) )
            d2  += ( hy ( y ) - ly ( y ) ) / max ( hy ( y ) , ly ( y ) )
            d3  += ( hz ( z ) - lz ( z ) ) / max ( hz ( z ) , lz ( z ) )
            d4  += ( hu ( u ) - lu ( u ) ) / max ( hu ( u ) , lu ( u ) )

            dp1 += ( cx ( x ) - lx ( x ) ) / max ( cx ( x ) , lx ( x ) )
            dp2 += ( cy ( y ) - ly ( y ) ) / max ( cy ( y ) , ly ( y ) )
            dp3 += ( cz ( z ) - lz ( z ) ) / max ( cz ( z ) , lz ( z ) )
            dp4 += ( cu ( u ) - lu ( u ) ) / max ( cu ( u ) , lu ( u ) )

            db1 += ( bx ( x ) - lx ( x ) ) / max ( cx ( x ) , lx ( x ) )
            db2 += ( by ( y ) - ly ( y ) ) / max ( cy ( y ) , ly ( y ) )
            db3 += ( bz ( z ) - lz ( z ) ) / max ( cz ( z ) , lz ( z ) )
            db4 += ( bu ( u ) - lu ( u ) ) / max ( cu ( u ) , lu ( u ) )

        logger.info ( '1D-(x)      DIFFERENCES are %s ' % d1  ) 
        logger.info ( '1D-(y)      DIFFERENCES are %s ' % d2  ) 
        logger.info ( '1D-(z)      DIFFERENCES are %s ' % d3  ) 
        logger.info ( '1D-(u)      DIFFERENCES are %s ' % d4  )
        
        logger.info ( '1D-(x) (L/C)DIFFERENCES are %s ' % dp1 ) 
        logger.info ( '1D-(y) (L/C)DIFFERENCES are %s ' % dp2 ) 
        logger.info ( '1D-(z) (L/C)DIFFERENCES are %s ' % dp3 ) 
        logger.info ( '1D-(u) (L/C)DIFFERENCES are %s ' % dp4 ) 

        logger.info ( '1D-(x) (B/L)DIFFERENCES are %s ' % db1 ) 
        logger.info ( '1D-(y) (L/L)DIFFERENCES are %s ' % db2 ) 
        logger.info ( '1D-(z) (B/L)DIFFERENCES are %s ' % db3 ) 
        logger.info ( '1D-(u) (B/L)DIFFERENCES are %s ' % db4 ) 


# =============================================================================
## 2D parameterizations
# =============================================================================
def test_parameterize_2D () :
    
    logger  = getLogger("test_parameterize_2D")
    logger.info ( 'Test 2D parameterisations' ) 

    with ROOT.TFile.Open(data_file,'READ') as f :
        
        tree = f.S

        with timing ( "2D-Legendre prameterization" , logger = logger ) :
            
            lxy = Ostap.Math.LegendreSum2 ( 12 , 12 , -2 , 2 , -2 , 2 )
            lzu = Ostap.Math.LegendreSum2 ( 12 , 12 , -4 , 4 , -4 , 6 )
            lxu = Ostap.Math.LegendreSum2 ( 12 , 12 , -2 , 2 , -4 , 6 )
            
            lxy.parameterize ( tree , 'x' , 'y' ,  cuts )
            lzu.parameterize ( tree , 'z' , 'u' ,  cuts )
            lxu.parameterize ( tree , 'x' , 'u' ,  cuts )

        with timing ( "2D-Bernstein prameterization" , logger = logger ) :
            
            bxy = Ostap.Math.Bernstein2D ( 12 , 12 , -2 , 2 , -2 , 2 )
            bzu = Ostap.Math.Bernstein2D ( 12 , 12 , -4 , 4 , -4 , 6 )
            bxu = Ostap.Math.Bernstein2D ( 12 , 12 , -2 , 2 , -4 , 6 )
            
            bxy.parameterize ( tree , 'x' , 'y' ,  cuts )
            bzu.parameterize ( tree , 'z' , 'u' ,  cuts )
            bxu.parameterize ( tree , 'x' , 'u' ,  cuts )
                 
        with timing ( "2D-histograms" , logger = logger ) :
            
            hxy = ROOT.TH2F(hID() , '', 20 , -2, 2, 20 , -2 , 2 )
            hzu = ROOT.TH2F(hID() , '', 20 , -4, 4, 20 , -4 , 6 )
            hxu = ROOT.TH2F(hID() , '', 20 , -2, 2, 20 , -4 , 6 )
            
            tree.project ( hxy , 'y:x' , cuts )
            tree.project ( hzu , 'u:z' , cuts )
            tree.project ( hxu , 'u:x' , cuts )

        lxy *= ( 4.0 / 20 ) * (  4.0 / 20 ) 
        lzu *= ( 8.0 / 20 ) * ( 10.0 / 20 ) 
        lxu *= ( 4.0 / 20 ) * ( 10.0 / 20 ) 

        bxy *= ( 4.0 / 20 ) * (  4.0 / 20 ) 
        bzu *= ( 8.0 / 20 ) * ( 10.0 / 20 ) 
        bxu *= ( 4.0 / 20 ) * ( 10.0 / 20 ) 

        d1 = SE()
        d2 = SE()
        d3 = SE() 
        
        for i in range ( 1000 ) :
            
            x   =  random.uniform ( -2 , 2 )
            y   =  random.uniform ( -2 , 2 )
            z   =  random.uniform ( -4 , 4 )
            u   =  random.uniform ( -4 , 6 )
            
            d1 += ( lxy ( x , y ) - hxy ( x , y ) ) / max ( lxy ( x , y ) , hxy ( x , y ) )
            d2 += ( lzu ( z , u ) - hzu ( z , u ) ) / max ( lzu ( z , u ) , hzu ( z , u ) )
            d3 += ( lxu ( x , u ) - hxu ( x , u ) ) / max ( lxu ( x , u ) , hxu ( x , u ) )
            
        logger.info ( '2D-(xy)-DIFFERENCES are %s ' % d1 ) 
        logger.info ( '2D-(zu)-DIFFERENCES are %s ' % d2 ) 
        logger.info ( '2D-(xu)-DIFFERENCES are %s ' % d3 ) 


# =============================================================================
## 3D parameterizations
# =============================================================================
def test_parameterize_3D () :
    
    logger  = getLogger("test_parameterize_3D")
    logger.info ( 'Test 3D parameterisations' ) 

    with ROOT.TFile.Open(data_file,'READ') as f :
        
        tree = f.S

        with timing ( "3D-Legendre prameterization" , logger = logger ) :
            
            lx = Ostap.Math.LegendreSum3 ( 8 , 8 , 8 , -2 , 2 , -4 , 4 , -4 , 6 )
            ly = Ostap.Math.LegendreSum3 ( 8 , 8 , 8 , -2 , 2 , -4 , 4 , -4 , 6 )
            lz = Ostap.Math.LegendreSum3 ( 8 , 8 , 8 , -2 , 2 , -2 , 2 , -4 , 6 )
            lu = Ostap.Math.LegendreSum3 ( 8 , 8 , 8 , -2 , 2 , -2 , 2 , -4 , 4 )
            
            lx.parameterize ( tree , 'y' , 'z' , 'u' ,  cuts )
            ly.parameterize ( tree , 'x' , 'z' , 'u' ,  cuts )
            lz.parameterize ( tree , 'x' , 'y' , 'u' ,  cuts )
            lu.parameterize ( tree , 'x' , 'y' , 'z' ,  cuts )

        with timing ( "3D-Bernsteinprameterization" , logger = logger ) :
            
            bx = Ostap.Math.Bernstein3D ( 8 , 8 , 8 , -2 , 2 , -4 , 4 , -4 , 6 )
            by = Ostap.Math.Bernstein3D ( 8 , 8 , 8 , -2 , 2 , -4 , 4 , -4 , 6 )
            bz = Ostap.Math.Bernstein3D ( 8 , 8 , 8 , -2 , 2 , -2 , 2 , -4 , 6 )
            bu = Ostap.Math.Bernstein3D ( 8 , 8 , 8 , -2 , 2 , -2 , 2 , -4 , 4 )
            
            bx.parameterize ( tree , 'y' , 'z' , 'u' ,  cuts )
            by.parameterize ( tree , 'x' , 'z' , 'u' ,  cuts )
            bz.parameterize ( tree , 'x' , 'y' , 'u' ,  cuts )
            bu.parameterize ( tree , 'x' , 'y' , 'z' ,  cuts )

        lxy = lz.integralZ()
        lzu = lx.integralX()
        lxu = ly.integralY()
        
        with timing ( "2D-histograms" , logger = logger ) :
            
            hxy = ROOT.TH2F(hID() , '', 20 , -2, 2, 20 , -2 , 2 )
            hzu = ROOT.TH2F(hID() , '', 20 , -4, 4, 20 , -4 , 6 )
            hxu = ROOT.TH2F(hID() , '', 20 , -2, 2, 20 , -4 , 6 )
            
            tree.project ( hxy , 'y:x' , cuts )
            tree.project ( hzu , 'u:z' , cuts )
            tree.project ( hxu , 'u:x' , cuts )
            
        lxy *= ( 4.0 / 20 ) * (  4.0 / 20 ) 
        lzu *= ( 8.0 / 20 ) * ( 10.0 / 20 ) 
        lxu *= ( 4.0 / 20 ) * ( 10.0 / 20 ) 
        
        d1 = SE()
        d2 = SE()
        d3 = SE() 
        
        for i in range ( 100 ) :
            
            x   =  random.uniform ( -2 , 2 )
            y   =  random.uniform ( -2 , 2 )
            z   =  random.uniform ( -4 , 4 )
            u   =  random.uniform ( -4 , 6 )
            
            d1 += ( lxy ( x , y ) - hxy ( x , y ) ) / max ( lxy ( x , y ) , hxy ( x , y ) )
            d2 += ( lzu ( z , u ) - hzu ( z , u ) ) / max ( lzu ( z , u ) , hzu ( z , u ) )
            d3 += ( lxu ( x , u ) - hxu ( x , u ) ) / max ( lxu ( x , u ) , hxu ( x , u ) )
            
        logger.info ( '3D-(xy)-DIFFERENCES are %s ' % d1 ) 
        logger.info ( '3D-(zu)-DIFFERENCES are %s ' % d2 ) 
        logger.info ( '3D-(xu)-DIFFERENCES are %s ' % d3 ) 



# =============================================================================
## 4D parameterizations
# =============================================================================
def test_parameterize_4D () :
    
    logger  = getLogger("test_parameterize_4D")
    logger.info ( 'Test 4D parameterisations' ) 

    with ROOT.TFile.Open(data_file,'READ') as f :
        
        tree = f.S
        
        l  = Ostap.Math.LegendreSum4 ( 8 , 8 , 8 , 8 ,
                                       -2 , 2 ,
                                       -2 , 2 ,
                                       -4 , 4 , -4 , 6 )
        
        l.parameterize ( tree , 'x' , 'y' , 'z' , 'u' ,  cuts )
        
        lxy = l.integralU().integralZ() 
        lzu = l.integralX().integralX() 
        lxu = l.integralY().integralY() 
        
        hxy = ROOT.TH2F(hID() , '', 20 , -2, 2, 20 , -2 , 2 )
        hzu = ROOT.TH2F(hID() , '', 20 , -4, 4, 20 , -4 , 6 )
        hxu = ROOT.TH2F(hID() , '', 20 , -2, 2, 20 , -4 , 6 )
        
        tree.project ( hxy , 'y:x' , cuts )
        tree.project ( hzu , 'u:z' , cuts )
        tree.project ( hxu , 'u:x' , cuts )
        
        lxy *= ( 4.0 / 20 ) * (  4.0 / 20 ) 
        lzu *= ( 8.0 / 20 ) * ( 10.0 / 20 ) 
        lxu *= ( 4.0 / 20 ) * ( 10.0 / 20 ) 
        
        d1 = SE()
        d2 = SE()
        d3 = SE() 
        
        for i in range ( 1000 ) :
            
            x   =  random.uniform ( -2 , 2 )
            y   =  random.uniform ( -2 , 2 )
            z   =  random.uniform ( -4 , 4 )
            u   =  random.uniform ( -4 , 6 )
        
            d1 += ( lxy ( x , y ) - hxy ( x , y ) ) / max ( lxy ( x , y ) , hxy ( x , y ) )
            d2 += ( lzu ( z , u ) - hzu ( z , u ) ) / max ( lzu ( z , u ) , hzu ( z , u ) )
            d3 += ( lxu ( x , u ) - hxu ( x , u ) ) / max ( lxu ( x , u ) , hxu ( x , u ) )
            
        logger.info ( '4D-(xy)-DIFFERENCES are %s ' % d1 ) 
        logger.info ( '4D-(zu)-DIFFERENCES are %s ' % d2 ) 
        logger.info ( '4D-(xu)-DIFFERENCES are %s ' % d3 ) 
        

    
# =============================================================================
if '__main__' == __name__ :

    test_parameterize_1D() 
    test_parameterize_2D() 
    test_parameterize_3D() 
    test_parameterize_4D() 
    
# =============================================================================
# The END 
# =============================================================================
