#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_shapes.py
# ============================================================================= 
""" Test module for ostap/fitting/models..py
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.meta_info     import root_info 
from   ostap.core.core          import Ostap, hID
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait, batch_env  
import ostap.fitting.models     as     Models
import ostap.histos.histos
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_shapes' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================
#

def gauss ( x , mu = 0 , sigma = 1 ) :
    return Ostap.Math.gauss_pdf ( x , mu , sigma ) 

# ============================================================================
def test_shapes_1d() :
    
    logger = getLogger("test_shapes_1d")


    if root_info < ( 6 , 18 ) :
        logger.info ( "Test is distabled for ROOT %s" % str ( root_info ) )
        return 

    ## C++ callable as shape

    xvar = ROOT.RooRealVar('x', '', -5, 5)
    
    s1 = Models.Shape1D_pdf ( 'S1'                                 ,
                              shape = Ostap.Math.BifurcatedGauss() ,
                              xvar  = xvar                        )

    with wait ( 1 ) , use_canvas ( "shape1d : C++ functor" ) :        
        s1.draw()

    ## histogram as shape
        
    h2  = ROOT.TH1D ( hID() , '' , 50 , -5 , 5 )
    h2 += lambda x : 100 * gauss ( x )
        
    s2 = Models.Shape1D_pdf ( 'S2'         ,
                              shape = h2   ,
                              xvar  = xvar ) 

    with wait ( 1 ) , use_canvas ( "shape1d : histogram " ) :        
        s2.draw()

# ============================================================================
def test_shapes_2d() :
    
    logger = getLogger("test_shapes_2d")

    if root_info < ( 6 , 18 ) :
        logger.info ( "Test is distabled for ROOT %s" % str ( root_info ) )
        return 

    ## histogram as shape
        
    h2  = ROOT.TH2D ( hID() , '' , 40 , -5 , 5 , 40 , -10 , 10 )
    
    h2 += lambda x,y : 100 * gauss ( x ) * gauss ( y , sigma = 2 ) 
    
    s2 = Models.Shape2D_pdf ( 'S3'               ,
                              shape = h2         ,
                              xvar  = ( -5  ,  5 ) ,  
                              yvar  = ( -10 , 10 ) )

    with wait ( 1 ) , use_canvas ( "shape2d : histogram/x" ) :        
        s2.draw1 ()
    with wait ( 1 ) , use_canvas ( "shape2d : histogram/x in y-range" ) :        
        s2.draw1 ( in_range = (-6,6) )

    with wait ( 1 ) , use_canvas ( "shape2d : histogram/y" ) :        
        s2.draw2 ()
    with wait ( 1 ) , use_canvas ( "shape2d : histogram/y in x-range" ) :        
        s2.draw2 ( in_range = (-3,3) )


# ============================================================================
def test_shapes_3d() :
    
    logger = getLogger("test_shapes_3d")

    if root_info < ( 6 , 18 ) :
        logger.info ( "Test is distabled for ROOT %s" % str ( root_info ) )
        return 

    ## histogram as shape
        
    h2  = ROOT.TH3D ( hID() , '' , 30 , -5 , 5 , 30 , -10 , 10 , 30 , -15 , 15 )
    
    h2 += lambda x,y,z : 100 * gauss ( x ) * gauss ( y , sigma = 2 ) * gauss ( z , sigma = 3 ) 
    
    s2 = Models.Shape3D_pdf ( 'S4'               ,
                              shape = h2         ,
                              xvar  = ( -5  , 5  ) ,  
                              yvar  = ( -10 , 10 ) ,
                              zvar  = ( -15 , 15 ) )

    with wait ( 1 ) , use_canvas ( "shape3d : histogram/x" ) :        
        s2.draw1 ()
    with wait ( 1 ) , use_canvas ( "shape3d : histogram/x in y-range" ) :        
        s2.draw1 ( in_range2 = (-3,3) )
    with wait ( 1 ) , use_canvas ( "shape3d : histogram/x in y&z-range" ) :        
        s2.draw1 ( in_range2 = (-3,3) , in_range3 = (-3,3) )


    with wait ( 1 ) , use_canvas ( "shape3d : histogram/y" ) :        
        s2.draw2 ()
    with wait ( 1 ) , use_canvas ( "shape3d : histogram/y in x-range" ) :        
        s2.draw2 ( in_range1 = (-3,3) )
    with wait ( 1 ) , use_canvas ( "shape3d : histogram/y in x&z-range" ) :        
        s2.draw2 ( in_range1 = (-3,3) , in_range3 = (-3,3) )


    with wait ( 1 ) , use_canvas ( "shape3d : histogram/z" ) :        
        s2.draw3 ()
    with wait ( 1 ) , use_canvas ( "shape3d : histogram/z in x-range" ) :        
        s2.draw3 ( in_range1 = (-3,3) )
    with wait ( 1 ) , use_canvas ( "shape3d : histogram/y in x&y-range" ) :        
        s2.draw3 ( in_range1 = (-3,3) , in_range2 = (-3,3) )


# =============================================================================
if '__main__' == __name__ :

    test_shapes_1d () 
    test_shapes_2d () 
    test_shapes_3d () 

# =============================================================================
##                                                                      The END 
# ============================================================================= 
