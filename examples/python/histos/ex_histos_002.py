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

h1 = ROOT.TH1F('h1','title',50,0,20)
h2 = ROOT.TH1F('h2','title',50,0,20)

## 1) set values for some bins
for i in  h1 :
    h1[i] = VE(i,i)

## 2) set values for some bins
for i in  h2 :
    h2[i] = VE(i*i,i*i)

## 3) dump them 
print 'Histogram H1:' , h1.dump(50,20) 
print 'Histogram H2:' , h1.dump(50,20) 

## 4) operations
print 'H1+H2' , (h1+h2).dump(50,20)
print 'H1-H2' , (h1-h2).dump(50,20)
print 'H1*H2' , (h1*h2).dump(50,20)
print 'H1/H2' , (h1/h2).dump(50,20)
print 'H2/H1' , (h2/h1).dump(50,20)

## 5) operations with constants 
c = 0.5 
print 'H1+1/2' , (h1+c).dump(50,20)
print 'H1-1/2' , (h1-c).dump(50,20)
print 'H1*1/2' , (h1*c).dump(50,20)
print 'H1/0.5' , (h1/c).dump(50,20)
print '1/2+H1' , (c+h1).dump(50,20)
print '1/2-H1' , (c-h1).dump(50,20)
print '1/2*H1' , (c*h1).dump(50,20)
print '0.5/H1' , (c/h1).dump(50,20)

## 6) operations with functions
fun = lambda x : x
print 'H1+x' , (h1+fun).dump(50,20)
print 'H1-x' , (h1-fun).dump(50,20)
print 'H1*x' , (h1*fun).dump(50,20)
print 'H1/x' , (h1/fun).dump(50,20)
print 'x+H1' , (fun+h1).dump(50,20)
print 'x-H1' , (fun-h1).dump(50,20)
print 'x*H1' , (fun*h1).dump(50,20)
print 'x/H1' , (fun/h1).dump(50,20)

## 7) math
from ostap.math.math_ve import sqrt, log  
print 'H1**2'    ,  (h1**2).dump(50,20)
print 'sqrt(H2)' , sqrt(h2).dump(50,20)
print 'log (H2)' , log(h2).dump(50,20)


# =============================================================================
# The END 
# =============================================================================
