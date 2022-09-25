#!/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright Ostap developers
# =============================================================================
## @file examples/histos/histos_ex001_basic.py
#  Simple example of basic operations with ROOT histograms
# =============================================================================
"""Simple example of basic operations with ROOT histograms
"""
# =============================================================================
from   __future__ import print_function
from   ostap.histos.histos import VE
import ROOT

h1 = ROOT.TH1F('h1','title',20,0,1)

## 1) set values for some bins
for i in  h1 :
    h1[i] = VE(i,i)
    
## 2) dump it! 
print( 'Histogram: %s'  % h1.dump(50,20)  )

## 3) basic properties:

print ( 'MEAN              %s' % h1.mean     () ) 
print ( 'RMS               %s' % h1.rms      () ) 
print ( 'SKEWNESS          %s' % h1.skewness () ) 
print ( 'KURTOSIS          %s' % h1.kurtosis () ) 
for i in range(1,10) :
    print ( 'MOMENT         %2d %s' % ( i , h1.moment ( i ) ) )
for i in range(1,10) :
    print ( 'CENTRAL MOMENT %2d %s' % ( i , h1.moment ( i ) ) ) 
print ( 'MIN/MAX           %s' % list ( h1. minmax() ) ) 
print ( 'XMIN/XMAX         %s' % list ( h1.xminmax() ) ) 
print ( 'YMIN/YMAX         %s' % list ( h1.yminmax() ) ) 
print ( '#BINS             %d' % len ( h1 )            )

# =============================================================================
##                                                                      The END 
# =============================================================================
