#!/usr/bin/env python

import cppyy
cppyy.load_reflection_info("libostap")
import ROOT

ve = cppyy.gbl.ostap.math.ValueWithError(2, 1)
print("value=%f, cov2=%f, error=%f"  % (ve.value(), ve.cov2(), ve.error()))


print cppyy.gbl.ostap.math.lomont_compare_double(1.011, 1.01,2)
print cppyy.gbl.ostap.math.next_double(0.1,1000)
# print ROOT.ostap.math.absMin(1.0, 2.0)
print cppyy.gbl.ostap.math.Equal_To(float)(1)(1,2)
print cppyy.gbl.ostap.math.knuth_equal_to_double(1,2)
# print ROOT.ostap.math.pow("int")(1,2)