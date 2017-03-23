#!/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright Ostap developers
# =============================================================================
"""Simple example of operations with ROOT histograms
"""
# =============================================================================
import ROOT,random 
from   ostap.histos import VE 

h1 = ROOT.TH1F('h1','title',10,0,10)

## 1) set values for some bins
for i in  h1 :
    h1[i] = VE(i*i,i*i)

    
## 2) use historgam as function

print 'Histogram as function: '
for i in range(10) :

    x = random.uniform ( *h1.xminmax() )

    v  = h1 ( x )                        ## use default interpolation 
    v0 = h1 ( x , interpolate = False )  ## no        interpolation
    v1 = h1 ( x , interpolate = 1     )  ## linear    interpolation
    v2 = h1 ( x , interpolate = 2     )  ## parabolic interpolation 
    v3 = h1 ( x , interpolate = 3     )  ## cubic     interpolation 

    print  'x=%s \tv=%s \tv0/v1/v2/v3=%s/%s/%s/%s ' % ( x , v , v0 , v1 , v2 , v3 ) 
 
# =============================================================================
# The END 
# =============================================================================
