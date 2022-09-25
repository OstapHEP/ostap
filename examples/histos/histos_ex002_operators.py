#!/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright Ostap developers
# =============================================================================
## @file examples/histos/histos_ex002_operators.py
#  Simple example of operations and operators with ROOT histograms
# =============================================================================
"""Simple example of operations and operators with ROOT histograms
"""
# =============================================================================
from   __future__ import print_function
from   ostap.histos.histos import VE 
import ROOT,random 

h1 = ROOT.TH1F('h1','title',50,0,20)
h2 = ROOT.TH1F('h2','title',50,0,20)
h3 = ROOT.TH1F('h3','title',80,0,40) ## different bining! 

## 1) set values for some bins
for i in  h1 :
    h1[i] = VE(i,i)

## 2) set values for some bins
for i in  h2 :
    h2[i] = VE(i*i,i*i)

## 3) dump them 
print ( 'Histogram H1: %s ' % h1.dump ( 50 , 20 ) )
print ( 'Histogram H2: %s ' % h1.dump ( 50 , 20 ) )

## 4) binary operations with histograms 
print ( 'H1+H2           %s ' % ( h1 + h2 ).dump(50,20) ) 
print ( 'H1-H2           %s ' % ( h1 - h2 ).dump(50,20) ) 
print ( 'H1*H2           %s ' % ( h1 * h2 ).dump(50,20) ) 
print ( 'H1/H2           %s ' % ( h1 / h2 ).dump(50,20) ) 
print ( 'H2/H1           %s ' % ( h2 / h1 ).dump(50,20) )
print ( 'H1/(H1+H2)      %s ' % ( h1.frac ( h2 ) ).dump(50,20) )
print ( '(H1-H2)/(H1+H2) %s ' % ( h1.asym ( h2 ) ).dump(50,20) )
 

## 5) binary operations with historgams (different binning)

print ( 'H1+H3           %s' % (h1+h3).dump(50,20)   ) 
print ( 'H1-H3           %s' % (h1-h3).dump(50,20)   ) 
print ( 'H1*H3           %s' % (h1*h3).dump(50,20)   )
print ( 'H1/H3           %s' % (h1/h3).dump(50,20)   )
print ( 'H1/(H1+H3)      %s' % (h1.frac(h3)).dump(50,20) )
print ( '(H1-H3)/(H1+H3) %s' % (h1.asym(h3)).dump(50,20) ) 

## 6) binary operations with constants

c = 0.5 
print ( 'H1+1/2   %s' % (h1+c).dump (50,20) ) 
print ( 'H1-1/2   %s' % (h1-c).dump (50,20) ) 
print ( 'H1*1/2   %s' % (h1*c).dump (50,20) ) 
print ( 'H1/0.5   %s' % (h1/c).dump (50,20) ) 
print ( '1/2+H1   %s' % (c+h1).dump (50,20) ) 
print ( '1/2-H1   %s' % (c-h1).dump (50,20) ) 
print ( '1/2*H1   %s' % (c*h1).dump (50,20) ) 
print ( '0.5/H1   %s' % (c/h1).dump (50,20) ) 

## 7) binary operations with functions
fun = lambda x : x
print ( 'H1+x     %s' % (h1+fun).dump(50,20) )
print ( 'H1-x     %s' % (h1-fun).dump(50,20) )
print ( 'H1*x     %s' % (h1*fun).dump(50,20) )
print ( 'H1/x     %s' % (h1/fun).dump(50,20) )
print ( 'x+H1     %s' % (fun+h1).dump(50,20) )
print ( 'x-H1     %s' % (fun-h1).dump(50,20) )
print ( 'x*H1     %s' % (fun*h1).dump(50,20) )
print ( 'x/H1     %s' % (fun/h1).dump(50,20) )

## 8) math
from ostap.math.math_ve import sqrt, log  
print ( 'H1**2    %s ' % (h1**2).dump(50,20)   )
print ( 'sqrt(H2) %s ' % sqrt(h2).dump(50,20)  )
print ( 'log (H2) %s ' % log(h2).dump(50,20)   )


# =============================================================================
##                                                                      The END 
# =============================================================================
