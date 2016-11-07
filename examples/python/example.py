#!/usr/bin/env python

## import cppyy
## cppyy.load_reflection_info("libostap")
import ROOT

## cpp    = cppyy.gbl
Ostap  = ROOT.Ostap

ve = Ostap.Math.ValueWithError(2, 2)
print("value=%f, cov2=%f, error=%f"  % (ve.value(), ve.cov2(), ve.error()))

#print Ostap.Math.lomont_compare_double(1.011, 1.01,2)
#print Ostap.math.next_double(0.1,1000)
#print Ostap.Math.Equal_To(float)(1)(1,2)
#print Ostap.Math.knuth_equal_to_double(1,2)
