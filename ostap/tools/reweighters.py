#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/reweghters.py
#  Density ratio reweighters & iterative reweighters based on LightGBM, XGBoost, CatBoost, SkLearn and HepML
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @Date  2026-07-18
# =============================================================================
""" Density ratio & iterative reweighters based on
 - LightGBM
 - XGBoost
 - CatBoost
 - SkLearn 
 - HepML.
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-06-07"
__all__     = (
    'LightGBMDensityReweighter'      , # LightGBM-based density reweighter
    'XGBoostDensityReweighter'       , # XGBoost-based density reweighter
    'CatBoostDensityReweighter'      , # CatBoost-based density reweighter
    'PyTorchDensityReweighter'       , # PyTorch-based density reweighter
    'LogRegressionDensityReweighter' , # Density reweighter based on Logistic Regression
    'GBReweighter'                   , # Reweighter based on GBReweighter from hep_ml    
) 
# =============================================================================
from  ostap.tools.dreweighters import *
from  ostap.tools.ireweighters import *
# =============================================================================
# Logging setup
# =============================================================================
from ostap.logger.logger import getLogger, logAttention 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.tools.reweighters' )
else                       : logger = getLogger( __name__ )
# =============================================================================

# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
#                                                                       The END 
# =============================================================================
