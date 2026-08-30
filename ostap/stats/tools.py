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
    'hasLightGBM'        , ## Are LigthGBM  library and classificator available?
    'hasXGBoost'         , ## Are XGBoost   library and claffificators available?
    'hasCatBoost'        , ## Are CatBoost  library and claffificators available?
    'hasPyTorch'         , ## Are (Py)Torch library and claffificators available?
    'hasKeras'           , ## Are Keras     library and claffificators available?
    'hasSkLearn'         , ## Are sckearn   library and claffificators available?
    'hasTensorFlow'      , ## Is  TensorFlow library available? 
    'hasHepML'           , ## Are HepML tools available?
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
##  Has  LightGBM ?
# - Are LightGBM library & classificators available? 
# - There is some mess with lightgbm&narwhals installation 
def hasLightGBM ( silent = True ) :
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
##  Has XGBoost ?
# - Are XGBoost library & classificators available? 
def hasXGBoost ( silent = True ) : 
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
##  Has CatBoost ?
# - Are CatBoost library & classificators available? 
def hasCatBoost ( silent = True ) : 
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

# =============================================================================
##  Has sklearn?
# - Has sklearn library & classificators available? 
def hasSkLearn ( silent = True ) :
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
        ## 
        import warnings
        warnings.filterwarnings ( "ignore", category = UserWarning , module = "sklearn.utils.parallel" )
        ## 
        return True 
        # ====================================================================
    except ImportError : # ===================================================
        # ====================================================================
        return False 
                    
# ==============================================================================
## Has PyTorch ?
#  Are (Py)Torch library and claffificators available?
def hasPyTorch ( silent = True ) :
    """ Use PyTorch ?
    - Are (Py)Torch library and claffificators available?
    """
    # ==========================================================================
    from ostap.core.cpu_info import HAS_AVX2
    if not HAS_AVX2 : return  False
    # ==========================================================================
    try : # ====================================================================
        # ======================================================================
        import torch
        if not silent : logger.info ( 'PyTorch  version: %s' %    torch .__version__ )
        return Version ( "2.1.0" ) <= Version ( torch.__version__ )
        # ======================================================================
    except ImportError : # =====================================================
        # ======================================================================
        return False 

# ==============================================================================
## Has TensorFlow ?
#  Is TensorFlow library available?
def hasTensorFlow ( silent = True ) :
    """ Has TensorFlow  ?
    - Is TensorFlow library available?
    """
    # ========================================================================== 
    from ostap.core.cpu_info import HAS_AVX2
    if not HAS_AVX2 : return  False
    # ==========================================================================
    try : # ====================================================================
        # ======================================================================
        ## silence TensorFlow & oneDNN
        os.environ [ 'TF_CPP_MIN_LOG_LEVEL'  ] = '2'        
        os.environ [ 'TF_ENABLE_ONEDNN_OPTS' ] = '0'
        # =====================================================================
        import tensorflow
        if not silent : logger.info ( 'TensorFlow version: %s' % tensorflow .__version__ )
        tensorflow.get_logger().setLevel ( 'ERROR' )
        return Version ( "2.16.1" ) <= Version ( tensorflow.__version__ )
        # ======================================================================
    except ImportError : # =====================================================
        # ======================================================================
        return False 
    
# ==============================================================================
##   Has Keras  ?
#  - Are Keras    library and claffificators available?
def hasKeras ( silent = True ) : 
    """ Use Keras
    - Are Keras    library and claffificators available?
    """
    # ========================================================================
    if not hasTensorFlow ( silent ) and not hasPyTorch ( silent ) : return False 
    # ========================================================================
    from ostap.core.cpu_info import HAS_AVX2
    if not HAS_AVX2 : return  False 
    # ========================================================================
    logger.warning ( 'Keras is temporarily(?) disabled!' )
    return False
    # ========================================================================
    try : # ==================================================================
        # ====================================================================
        import keras 
        if not silent : logger.info ( 'Keras    version: %s' %    keras .__version__ )
        return Version ( "3.0" ) <= Version ( keras.__version__ ) 
        # ====================================================================
    except ImportError : # ===================================================
        # ====================================================================
        return False

# ============================================================================
##  Has hep_ml
# - Are hep_ml tools available? 
def hasHepML ( silent = True ) :
    """ Use HepML?
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
    
    if not hasLightGBM ( False ) : logger.warning  ( "No LightGBM available!" ) 
    if not hasXGBoost  ( False ) : logger.warning  ( "No XGBoost  available!" ) 
    if not hasCatBoost ( False ) : logger.warning  ( "No CatBoost available!" ) 
    if not hasSkLearn  ( False ) : logger.warning  ( "No scikit-learn available!" ) 
    if not hasPyTorch  ( False ) : logger.warning  ( "No PyTorch  available!" ) 
    if not hasKeras    ( False ) : logger.warning  ( "No Keras    available!" ) 
    if not hasHepML    ( False ) : logger.warning  ( "No HepML    available!" ) 

# =============================================================================
##                                                                      The END 
# =============================================================================

    
