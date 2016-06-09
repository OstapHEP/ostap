#!/usr/bin/env python

import ROOT
ROOT.gSystem.Load("libostap")

ve = ROOT.ostap.ValueWithError(2, 1)
print("value=%f, cov2=%f, error=%f"  % (ve.value(), ve.cov2(), ve.error()))