#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
# =============================================================================
## @file fixes.py 
#  Couple of minor fixes for Ostap
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-02-23
# =============================================================================
"""Couple o fminot fixes for Ostap
"""
# =============================================================================
__version__ = '$Revision$'
__author__  = 'Vanya BELYAEV Ivan.Belyaev@itep.ru'
__date__    = '2016-02-23'
__all__     = () ## noting to import 
# =============================================================================
import ROOT 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.fixes.fixes')
else                      : logger = getLogger ( __name__           ) 
# =============================================================================
from ostap.logger.utils import mute
with mute() :
    v = ROOT.RooRealVar()
    del v

# =============================================================================
# The EDN
# =============================================================================
