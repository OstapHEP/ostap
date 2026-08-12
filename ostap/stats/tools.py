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
def useLightGBM () :
    """ Use LightGBM ?
    - Are LightGBM library & classificators available? 
    - There is some mess with LightBM&Narwhals installation 
    """
    # ============================================================================
    try : # ======================================================================
        # ========================================================================
        import lightgbm
        logger.info ( 'LightGBM version : %s' % lightgbm.__version__ ) 
        if Version ( lightgbm.__version__ ) <  Version ( "4.7.0"  ) : return True
        import narwhals
        logger.info ( 'Narwhals version : %s' % narwhals.__version__ ) 
        return Version ( "2.0" ) <= Version ( narwhals.__version__ )
        # ========================================================================
    except ImportError : # =======================================================
        # ========================================================================
        return False 
    
# ===============================================================================
## use XGBoost ?
# - Are XGBoost library & classificators available? 
def useXGBoost () : 
    """ Use XGBoost
     - Are XGBoost library & classificators available? 
    """
    # ==========================================================================
    try : # ====================================================================
        # ======================================================================
        import xgboost        
        return Version ( "1.0" ) <= Version ( xgboost.__version__ )
        # ======================================================================
    except ImportError : # =====================================================
        # ======================================================================
        return False 

# ===============================================================================
## use CatBoost ?
# - Are CatBoost library & classificators available? 
def useCatBoost () : 
    """ Use CatBoost
    - Are CatBoost library & classificators available? 
    """
    from ostap.core.cpu_info import HAS_AVX2
    if not HAS_AVX2 : return  False 
    # ==========================================================================
    try : # ====================================================================
        # ======================================================================
        import catboost
        return True 
        # ======================================================================
    except ImportError : # =====================================================
        # ======================================================================
        return False

# ==============================================================================
## use PyTorch ?
#  Are (Py)Torch library and claffificators available?
def usePyTorch() :
    """ Use PyTorch ?
    - Are (Py)Torch library and claffificators available?
    """
    # ==========================================================================
    try : # ====================================================================
        # ======================================================================
        import torch
        return Version ( "1.10" ) <= Version ( torch.__version__ )
        # ======================================================================
    except ImportError : # =====================================================
        # ======================================================================
        return False 
    
# ==============================================================================
## use Keras  ?
#  - Are Keras    library and claffificators available?
def useKeras() : 
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
        return  ( Version ( "3.0" ) <= Version ( keras.__version__ ) and
                  Version ( "2.0" ) <= Version ( torch.__version__ ) )
        # ======================================================================
    except ImportError : # =====================================================
        # ======================================================================
        return False

# =============================================================================
## Use sklearn?
# - Are sklearn library & classificators available? 
def useSkLearn () :
    """ Use sklearn?
    - Are sklearn library & classificators available? 
    """
    # ===========================================================================
    try : # =====================================================================
        # =======================================================================
        import sklearn.ensemble 
        from   sklearn.ensemble import HistGradientBoostingClassifier as _HGBC 
        from   sklearn.ensemble import     GradientBoostingClassifier as _GBC        
        from   sklearn.ensemble import         RandomForestClassifier as _RFC
        # ======================================================================
        return True 
        # ======================================================================
    except ImportError : # =====================================================
        # ======================================================================
        return False 
                
# =============================================================================
## Use hep_ml
# - Are hep_ml tools available? 
def useHepML () :
    """ Use sklearn?
    - Are hep_ml tools available? 
    """
    # ===========================================================================
    try : # =====================================================================
        # =======================================================================
        from hep_ml.reweight import      GBReweighter as _GBRW
        from hep_ml.reweight import FoldingReweighter as _FRW
        # ======================================================================
        return True 
        # ======================================================================
    except ImportError : # =====================================================
        # ======================================================================
        return False 
    
# ==============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================

    
