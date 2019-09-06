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
import ROOT, os,  random
import ostap.core.pyrouts
import ostap.trees.param
import ostap.math.models
from   ostap.core.core    import hID, SE, Ostap 
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
data_file = CleanUp.tempfile ( suffix = '.root' , prefix = 'test_trees_param_' ) 

if not os.path.exists( data_file ) :
    
    import random
    
    N =  200000
        
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

        I = 0
        while I < N : 
            
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
            I+= 1
            
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
    
    with ROOT.TFile.Open(data_file,'READ') as f :
        
        tree = f.S
        
        lx = Ostap.Math.LegendreSum ( 4 , -2 , 2 )
        ly = Ostap.Math.LegendreSum ( 4 , -2 , 2 )
        lz = Ostap.Math.LegendreSum ( 4 , -4 , 4 )
        lu = Ostap.Math.LegendreSum ( 4 , -4 , 6 )
        
        lx.parameterize ( tree , 'x' , cuts )
        ly.parameterize ( tree , 'y' , cuts )
        lz.parameterize ( tree , 'z' , cuts )
        lu.parameterize ( tree , 'u' , cuts )
        
        hx = ROOT.TH1D(hID(),'',100,-2,2)
        hy = ROOT.TH1D(hID(),'',100,-2,2)
        hz = ROOT.TH1D(hID(),'',100,-4,4)
        hu = ROOT.TH1D(hID(),'',100,-4,6)
        
        tree.project ( hx , 'x' , cuts )
        tree.project ( hy , 'y' , cuts )
        tree.project ( hz , 'z' , cuts )
        tree.project ( hu , 'u' , cuts )
        
        hx.SetMinimum(0)
        hx.draw()
        lx *= 0.04   ## bin-width
        lx.draw('same', linecolor=2)
        
        hy.SetMinimum(0)
        hy.draw()
        ly *= 0.04   ## bin-width
        ly.draw('same', linecolor=2)

        hz.SetMinimum(0)
        hz.draw()
        lz *= 0.08   ## bin-width
        lz.draw('same', linecolor=2)
        
        hu.SetMinimum(0)
        hu.draw()
        lu *= 0.10   ## bin-width
        lu.draw('same', linecolor=2)

        d1 = SE()
        d2 = SE()
        d3 = SE()
        d4 = SE()

        for i in range ( 1000 ) :
            
            x = random.uniform ( -2 , 2 )
            y = random.uniform ( -2 , 2 )
            z = random.uniform ( -4 , 4 )
            u = random.uniform ( -4 , 6 )

            d1 += ( hx ( x ) - lx ( x ) ) / max ( hx ( x ) , lx ( x ) )
            d2 += ( hy ( y ) - ly ( y ) ) / max ( hy ( y ) , ly ( y ) )
            d3 += ( hz ( z ) - lz ( z ) ) / max ( hz ( z ) , lz ( z ) )
            d4 += ( hu ( u ) - lu ( u ) ) / max ( hu ( u ) , lu ( u ) )
            
        logger.info ( '1D-(x)-DIFFERENCES are %s ' % d1 ) 
        logger.info ( '1D-(y)-DIFFERENCES are %s ' % d2 ) 
        logger.info ( '1D-(z)-DIFFERENCES are %s ' % d3 ) 
        logger.info ( '1D-(u)-DIFFERENCES are %s ' % d4 ) 


# =============================================================================
## 2D parameterizations
# =============================================================================
def test_parameterize_2D () :
    
    with ROOT.TFile.Open(data_file,'READ') as f :
        
        tree = f.S
        
        lxy = Ostap.Math.LegendreSum2 ( 12 , 12 , -2 , 2 , -2 , 2 )
        lzu = Ostap.Math.LegendreSum2 ( 12 , 12 , -4 , 4 , -4 , 6 )
        lxu = Ostap.Math.LegendreSum2 ( 12 , 12 , -2 , 2 , -4 , 6 )
        
        lxy.parameterize ( tree , 'x' , 'y' ,  cuts )
        lzu.parameterize ( tree , 'z' , 'u' ,  cuts )
        lxu.parameterize ( tree , 'x' , 'u' ,  cuts )
        
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
            
        logger.info ( '2D-(xy)-DIFFERENCES are %s ' % d1 ) 
        logger.info ( '2D-(zu)-DIFFERENCES are %s ' % d2 ) 
        logger.info ( '2D-(xu)-DIFFERENCES are %s ' % d3 ) 


# =============================================================================
## 3D parameterizations
# =============================================================================
def test_parameterize_3D () :
    
    with ROOT.TFile.Open(data_file,'READ') as f :
        
        tree = f.S
        
        lx = Ostap.Math.LegendreSum3 ( 8 , 8 , 8 , -2 , 2 , -4 , 4 , -4 , 6 )
        ly = Ostap.Math.LegendreSum3 ( 8 , 8 , 8 , -2 , 2 , -4 , 4 , -4 , 6 )
        lz = Ostap.Math.LegendreSum3 ( 8 , 8 , 8 , -2 , 2 , -2 , 2 , -4 , 6 )
        lu = Ostap.Math.LegendreSum3 ( 8 , 8 , 8 , -2 , 2 , -2 , 2 , -4 , 4 )
        
        lx.parameterize ( tree , 'y' , 'z' , 'u' ,  cuts )
        ly.parameterize ( tree , 'x' , 'z' , 'u' ,  cuts )
        lz.parameterize ( tree , 'x' , 'y' , 'u' ,  cuts )
        lu.parameterize ( tree , 'x' , 'y' , 'z' ,  cuts )
        
        lxy = lz.integralZ()
        lzu = lx.integralX()
        lxu = ly.integralY()
        
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
