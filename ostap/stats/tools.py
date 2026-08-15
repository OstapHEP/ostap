#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/tools.py
#  Advanced Tools used for statistics 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2023-12-06
# =============================================================================
""" Simple utulities for goodness-of-fit studies 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2023-12-06"
__all__     = (
    # =========================================================================
    'useLightGBM'        , ## Are LigthGBM  library and classificator available?
    'useXGBoost'         , ## Are XGBoost   library and claffificators available?
    'useCatBoost'        , ## Are CatBoost  library and claffificators available?
    'usePyTorch'         , ## Are (Py)Torch library and claffificators available?
    'useKeras'           , ## Are Keras     library and claffificators available?
    'useSkLearn'         , ## Are sckearn   library and claffificators available?
    'useHepML'           , ## Are HepML tools available?
    # =========================================================================
)
# =============================================================================
from   packaging.version        import Version
import os 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.tools' )
else                       : logger = getLogger( __name__ )
# =============================================================================
## use LightGBM ?
# - Are LightGBM library & classificators available? 
# - There is some mess with lightgbm&narwhals installation 
def useLightGBM ( silent = True ) :
    """ Use LightGBM ?
    - Are LightGBM library & classificators available? 
    - There is some mess with LightBM&Narwhals installation 
    """
    # ==========================================================================
    try : # ====================================================================
        # ======================================================================
        import lightgbm
        if not silent : logger.info ( 'LightGBM version: %s' % lightgbm.__version__ ) 
        if Version ( lightgbm.__version__ ) <  Version ( "4.7.0"  ) : return True
        import narwhals
        if not silent : logger.info ( 'Narwhals version: %s' % narwhals.__version__ ) 
        return Version ( "2.0" ) <= Version ( narwhals.__version__ )
        # ======================================================================
    except ImportError : # =====================================================
        # ======================================================================
        return False 
    
# ==============================================================================
## use XGBoost ?
# - Are XGBoost library & classificators available? 
def useXGBoost ( silent = True ) : 
    """ Use XGBoost
     - Are XGBoost library & classificators available? 
    """
    # ==========================================================================
    try : # ====================================================================
        # ======================================================================
        import xgboost        
        if not silent : logger.info ( 'XGBoost  version: %s' % xgboost .__version__ )
        return Version ( "1.0" ) <= Version ( xgboost.__version__ )
        # ======================================================================
    except ImportError : # =====================================================
        # ======================================================================
        return False 

# ===============================================================================
## use CatBoost ?
# - Are CatBoost library & classificators available? 
def useCatBoost ( silent = True ) : 
    """ Use CatBoost
    - Are CatBoost library & classificators available? 
    """
    from ostap.core.cpu_info import HAS_AVX2
    if not HAS_AVX2 : return  False 
    # ==========================================================================
    try : # ====================================================================
        # ======================================================================
        import catboost
        if not silent : logger.info ( 'CatBoost version: %s' % catboost .__version__ )
        return True 
        # ======================================================================
    except ImportError : # =====================================================
        # ======================================================================
        return False

# ==============================================================================
## use PyTorch ?
#  Are (Py)Torch library and claffificators available?
def usePyTorch ( silent = True ) :
    """ Use PyTorch ?
    - Are (Py)Torch library and claffificators available?
    """
    # ==========================================================================
    try : # ====================================================================
        # ======================================================================
        import torch
        if not silent : logger.info ( 'PyTorch  version: %s' %    torch .__version__ )
        return Version ( "1.10" ) <= Version ( torch.__version__ )
        # ======================================================================
    except ImportError : # =====================================================
        # ======================================================================
        return False 
    
# ==============================================================================
## use Keras  ?
#  - Are Keras    library and claffificators available?
def useKeras ( silent = True ) : 
    """ Use Keras
    - Are Keras    library and claffificators available?
    """
    from ostap.core.cpu_info import HAS_AVX2
    if not HAS_AVX2 : return  False 
    # ==========================================================================
    try : # ====================================================================
        # ======================================================================
        ## silence TensorFlow & oneDNN
        os.environ [ 'TF_CPP_MIN_LOG_LEVEL'  ] = '2'        
        os.environ [ 'TF_ENABLE_ONEDNN_OPTS' ] = '0'
        # ======================================================================
        import keras 
        import torch 
        if not silent :
            logger.info ( 'PyTorch  version: %s' %    torch .__version__ )
            logger.info ( 'Keras    version: %s' %    keras .__version__ )
        return  ( Version ( "3.0" ) <= Version ( keras.__version__ ) and
                  Version ( "2.0" ) <= Version ( torch.__version__ ) )
        # ======================================================================
    except ImportError : # =====================================================
        # ======================================================================
        return False

# =============================================================================
## Use sklearn?
# - Are sklearn library & classificators available? 
def useSkLearn ( silent = True ) :
    """ Use sklearn?
    - Are sklearn library & classificators available? 
    """
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        import sklearn 
        from   sklearn.ensemble import HistGradientBoostingClassifier as _HGBC 
        from   sklearn.ensemble import     GradientBoostingClassifier as _GBC        
        from   sklearn.ensemble import         RandomForestClassifier as _RFC
        # ====================================================================
        if not silent : logger.info ( 'sklearn  version: %s' %  sklearn.__version__ )
        return True 
        # ====================================================================
    except ImportError : # ===================================================
        # ====================================================================
        return False 
                
# ============================================================================
## Use hep_ml
# - Are hep_ml tools available? 
def useHepML ( silent = True ) :
    """ Use sklearn?
    - Are hep_ml tools available? 
    """
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        import hep_ml 
        from hep_ml.reweight import      GBReweighter as _GBRW
        from hep_ml.reweight import FoldingReweighter as _FRW
        # ====================================================================
        if not silent : logger.info ( 'hep_ml   version: %s' % hep_ml.__version__ )
        return True 
        # ====================================================================
    except ImportError : # ===================================================
        # ====================================================================
        return False 
    
# ============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    if not useLightGBM ( False ) : logger.warning  ( "No LightGBM available!" ) 
    if not useXGBoost  ( False ) : logger.warning  ( "No XGBoost  available!" ) 
    if not useCatBoost ( False ) : logger.warning  ( "No CatBoost available!" ) 
    if not useSkLearn  ( False ) : logger.warning  ( "No scikit-learn available!" ) 
    if not usePyTorch  ( False ) : logger.warning  ( "No PyTorch available!" ) 
    if not useKeras    ( False ) : logger.warning  ( "No Keras   available!" ) 

# =============================================================================
##                                                                      The END 
# =============================================================================

    
