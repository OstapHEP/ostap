#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_plotting_sum.py
# Test sum of RooPlto objects  
# ============================================================================= 
""" Test sum of RoPlot objects 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.root_utils   import batch_env 
import ostap.fitting.models     as     Models
import ostap.histos.graphs
import ROOT 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_plotting_sum' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
batch_env ( logger ) 
# =============================================================================


x   = ROOT.RooRealVar  ( 'x' , 'x-obsevable' , 0 , 10 )
g   = Models.Gauss_pdf ( 'G' , xvar = x , mean = 5 , sigma = 1 )
ds1 = g.generate ( 1000 )
ds2 = g.generate ( 1000 )

title = 'The first plot'
with use_canvas ( title , wait = 2 ) :
    r1 , f1 = g.fitTo ( ds1 , draw = True , nbins = 100 , silent = True )
    f1.draw()
    logger.info ( '%s:\n%s' % ( title , f1.table ( prefix = '# ' ) ) )  

                 
title = 'The second plot'
with use_canvas ( title , wait = 2 ) :
    r2 , f2 = g.fitTo ( ds2 , draw = True , nbins = 100 , silent = True )
    f2.draw()
    logger.info ( '%s:\n%s' % ( title , f2.table ( prefix = '# ' ) ) )  

title = 'The sum'
with use_canvas ( title , wait = 5 ) :
    fsum = f1 + f2
    fsum.draw()
    logger.info ( '%s:\n%s' % ( title , fsum.table ( prefix = '# ' ) ) )  


assert len ( f1 ) == len ( f2 ) and len ( f1 ) == len ( fsum ) , \
       'Invalid dimentions/1!'


for a, b, c in zip ( f1 , f2  , fsum ) :
    
    if not  isinstance ( a , ROOT.RooHist ) : continue
    
    assert type ( a ) is type ( b )  and type ( b ) is type  ( c  ) , \
           "Invalid typeS!"
    
    assert len ( a ) == len ( b ) and len ( b ) == len ( c ) , \
           'Invalid dimentions/1!'
    
    for i, j, k in zip ( a , b , c ) :
        
        v1x, v1y = a [ i ]
        v2x, v2y = b [ i ]
        vsx, vsy = c [ i ]
        
        assert v1x.value == v2x.value and v2x.value == vsx.value, \
              'Invalid x: %s %s %s ' %  ( v1x , v2x , vsx )
        
        
        assert v1y.value + v2y.value == vsy.value, \
              'Invalid y: %s %s %s ' %  ( v1y , v2y , vsy ) 
        

## h1 = f1.getObject(0)
## h2 = f2.getObject(0)

## print ( 'nominal bin widths h1' , h1.getNominalBinWidth() )
## print ( 'nominal bin widths h2' , h1.getNominalBinWidth() )


## h1 = f1.getObject(2)
## h2 = f2.getObject(2)

## print ( 'nominal bin widths h1' , h1.getNominalBinWidth() )
## print ( 'nominal bin widths h2' , h1.getNominalBinWidth() )



    
# =============================================================================
if '__main__' == __name__  :
    pass

# =============================================================================
##                                                                      The END 
# =============================================================================
