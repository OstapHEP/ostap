#!/usr/bin/env python

import ROOT
ROOT.gSystem.Load("libostap")

ve = ROOT.ostap.math.ValueWithError(2, 1)
print("value=%f, cov2=%f, error=%f"  % (ve.value(), ve.cov2(), ve.error()))


print ROOT.ostap.math.lomont_compare_double(1.011, 1.01,2)
print ROOT.ostap.math.next_double(0.1,1000)
print ROOT.ostap.math.absMin(1.0,2.0)
print ROOT.ostap.math.Equal_To(float)(1)