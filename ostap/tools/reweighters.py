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
    ## 
    'LightGBMDensityReweighter'       , # LightGBM-based density reweighter
    'XGBoostDensityReweighter'        , # XGBoost-based density reweighter
    'CatBoostDensityReweighter'       , # CatBoost-based density reweighter
    'PyTorchAda,DensityReweighter'    , # PyTorch-based density reweighter
    'PyTorchLBFGSDensityReweighter'   , # PyTorch-based density reweighter
    'LogRegressionDensityReweighter'  , # Density reweighter based on Logistic Regression
    ##
    'LightGBMIterativeReweighter'     , # LightGBM-based iterative density reweighter
    'XGBoostIterativeReweighter'      , # XGBoost-based  iterative density reweighter
    'CatBoostIterativeReweighter'     , # CatBoost-based iterative density reweighter
    'PyTorchAdamIterativeReweighter'  , # PyTorch-based  iterative density reweighter
    'PyTorchLBFGSIterativeReweighter' , # PyTorch-based  iterative density reweighter
    ## 
    'GBReweighter'                   , # Reweighter based on GBReweighter from hep_ml by Alex Rogozhnikov
    ## 
) # ===========================================================================
# =============================================================================
from  ostap.tools.dreweighters import  ( LightGBMDensityReweighter       , 
                                         XGBoostDensityReweighter        , 
                                         CatBoostDensityReweighter       , 
                                         PyTorchAdamDensityReweighter    , 
                                         PyTorchLbfgsDensityReweighter   , 
                                         LogRegressionDensityReweighter  )
from  ostap.tools.ireweighters import  ( LightGBMIterativeReweighter     , 
                                         XGBoostIterativeReweighter      , 
                                         CatBoostIterativeReweighter     ,
                                         PyTorchAdamIterativeReweighter  ,
                                         PyTorchLBFGSIterativeReweighter )
from  ostap.tools.greweighters import GBReweighter
# =============================================================================
# Logging setup
# =============================================================================
from ostap.logger.logger import getLogger, logAttention 
logger = getLogger ( 'ostap.tools.reweighters' )
# =============================================================================

# ============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    from ostap.stats.tools import ( hasLightGBM ,
                                    hasXGBoost  ,
                                    hasCatBoost ,
                                    hasPyTorch  ,
                                    hasSkLearn  , 
                                    hasHepML    )

    if not hasLightGBM ( False ) : logger.warning  ( "No LightGBM available!" ) 
    if not hasXGBoost  ( False ) : logger.warning  ( "No XGBoost  available!" ) 
    if not hasCatBoost ( False ) : logger.warning  ( "No CatBoost available!" ) 
    if not hasPyTorch  ( False ) : logger.warning  ( "No PyTorch  available!" ) 
    if not hasSkLearn  ( False ) : logger.warning  ( "No SkLearn  available!" ) 
    if not hasHepML    ( False ) : logger.warning  ( "No HepML    available!" ) 

# =============================================================================
#                                                                       The END 
# =============================================================================
