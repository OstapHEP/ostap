#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/tests/test_tools_reweight5.py
#  Test for nD-reweighting machinery
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2023-01-20
# =============================================================================
"""Test for nD-reweighting machinery
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2023-01-20"
__all__     = ()  ## nothing to be imported 
# =============================================================================
from   ostap.utils.core         import typename 
from   ostap.utils.timing       import timing
from   ostap.logger.colorized   import allright
from   ostap.utils.root_utils   import batch_env 
from   ostap.logger.symbols     import script_p
from   ostap.utils.memory       import memory_usage, delta_ram
from   ostap.utils.basic        import numcpu 
from   ostap.utils.progress_bar import progress_bar 
from   ostap.stats.tools        import ( hasLightGBM , hasXGBoost ,
                                         hasCatBoost , hasSkLearn ,
                                         hasHepML    ) 
import ostap.logger.table       as     T 
import ostap.io.root_file 
import ostap.parallel.kisa
import ROOT, random, math, os, numpy 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__'  == __name__ : 
    logger = getLogger ( 'ostap.test_tools_reweight5' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================    
logger.info ( 'Test for nD-Reweighting machinery')
# ============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================
## Generate two N-dimensional numpy datasets for density reweighting tests.
def make_datasets ( n_samples = 10000 ,
                    n_dim     =     3 ,
                    seed      = None  ) :
    """ Generate two N-dimensional numpy datasets for density reweighting tests.
    
    Parameters:
    -----------
    n_samples : int
        Number of events in each dataset.
    n_dim : int
        Feature space dimension (N).
    seed : int
        Random seed for reproducibility.
        
    Returns:
    --------
    X_nontrivial : numpy.ndarray, shape (n_samples, n_dim)
        Dataset with non-linear 2D/3D+ correlations.
    X_uniform : numpy.ndarray, shape (n_samples, n_dim)
        Uniformly distributed dataset covering the bounding box.
    """
    numpy.random.seed ( seed )
    
    ## 1. Nontrivial dataset with strong non-linear correlations
    X_nontrivial = numpy.zeros ( ( n_samples , n_dim ) )
    
    ## Primary base Gaussian variable x_0
    z0 = numpy.random.normal ( loc = 0.0 , scale = 1.0 , size = n_samples )
    X_nontrivial [ : , 0 ] = z0
    
    if n_dim > 1 :
        ## x_1: Quadratic ("banana") correlation with x_0
        X_nontrivial [ : , 1 ] = z0**2 - 1.0 + numpy.random.normal ( loc = 0.0 , scale = 0.3 , size = n_samples )
        
    if n_dim > 2 :
        ## x_2: Oscillating joint dependence on x_0 and x_1
        X_nontrivial [ : , 2 ] = (
            numpy.sin ( 1.8 * z0 ) 
            + 0.4 * X_nontrivial [ : , 1 ] 
            + numpy.random.normal ( loc = 0.0 , scale = 0.2 , size = n_samples )
        )
        
    ## Recursive construction for higher dimensions (d >= 3)
    for d in range ( 3 , n_dim ) :
        X_nontrivial [ : , d ] = (
            0.5 * X_nontrivial [ : , d - 1 ] 
            + numpy.cos ( z0 ) 
            + numpy.random.normal ( loc = 0.0 , scale = 0.3 , size = n_samples )
        )
        
    ## 2. Uniform dataset spanning the bounding box of nontrivial distribution
    mins = numpy.min ( X_nontrivial , axis = 0 )
    maxs = numpy.max ( X_nontrivial , axis = 0 )
    
    X_uniform = numpy.random.uniform (
        low  = mins , 
        high = maxs , 
        size = ( n_samples , n_dim )
    )
    
    return X_nontrivial , X_uniform

# =========================================================================
has_lightgbm  = hasLightGBM  ()
if has_lightgbm :  logger.attention ( 'USE LightGBM!'              )
else            :  logger.warning   ( 'LightGBM is not available!' )
            
has_xgboost   = hasXGBoost  ()
if has_xgboost  :  logger.attention ( 'USE XGBoost!'               )
else            :  logger.warning   ( 'XGBoost is not available!'  )

has_catboost  = hasCatBoost  ()
if has_catboost :  logger.attention ( 'USE CatBoost!'              )
else            :  logger.warning   ( 'CatBoost is not available!' )

has_sklearn   = hasSkLearn  ()
if has_sklearn  :  logger.attention ( 'USE SkLearn!'              )
else            :  logger.warning   ( 'SkLearn  is not available!' )

has_hepml     = hasHepML  ()
if has_hepml    :  logger.attention ( 'USE HepML!'                 )
else            :  logger.warning   ( 'HepML    is not available!' )

# ==============================================================================
## Compare datasets using several methods 
# ==============================================================================
from ostap.stats.gof_np       import  ( MIXnp           as COMPARATOR4 , 
                                        KullbackLeibler as COMPARATOR3 , 
                                        Hotelling       as COMPARATOR2 , 
                                        Mahalanobis     as COMPARATOR1 ) 

from ostap.stats.data_compare import data_compare     
comparators = ( COMPARATOR1 ( parallel = True , nToys = 100 ) ,
                COMPARATOR2 ( parallel = True , nToys = 100 ) ,
                COMPARATOR3 ( parallel = True , nToys = 100 ) ,
                COMPARATOR4 ( parallel = True , nToys = 100 ) )

if has_lightgbm :  
    from ostap.stats.adval        import ADVAL_LGBM  as COMPARATOR5
    comparators += ( COMPARATOR5 ( parallel = True , nToys = 100 ) , ) 

if has_xgboost :  
    from ostap.stats.adval        import ADVAL_XGB  as COMPARATOR6
    comparators += ( COMPARATOR6 ( parallel = True , nToys = 100 ) , ) 

if has_catboost:  
    from ostap.stats.adval        import ADVAL_CATB  as COMPARATOR7
    comparators += ( COMPARATOR7 ( parallel = True , nToys = 100 ) , ) 

if False and has_sklearn:    
    from ostap.stats.adval        import ADVAL_HGBC  as COMPARATOR8
    comparators += ( COMPARATOR8 ( parallel = True , nToys = 100 ) , )
    
    from ostap.stats.adval        import ADVAL_GBC   as COMPARATOR9
    comparators += ( COMPARATOR9 ( parallel = True , nToys = 100 ) , ) 


# ============================
def run_reweight (  n_dim     = 3     ,
                    n_samples = 10000 ) : 
    
    ## generate samples 
    target , original = make_datasets ( n_dim = n_dim , n_samples = n_samples) 
    
    reweighters = [ None ]
    
    if has_hepml :         
        from ostap.tools.reweighters  import GBReweighter   as GBRW
        rw1 = GBRW ( target = target , original = original )
        reweighters.append ( rw1 )
        
    if has_lightgbm : 
        from ostap.tools.reweighters  import LightGBMDensityReweighter as  LGBM
        rw2 = LGBM ( target = target , original = original )
        reweighters.append ( rw2 )

    if has_xgboost : 
        from ostap.tools.reweighters  import XGBoostDensityReweighter as  XGB 
        rw3 = XGB  ( target = target , original = original )
        reweighters.append ( rw3 )

    if has_catboost : 
        from ostap.tools.reweighters  import CatBoostDensityReweighter as CATB
        rw4 = CATB ( target = target , original = original )
        reweighters.append ( rw4 )

    header = []
    rows   = []     
    for rw in reweighters :

        if rw is None : original_weight = None 
        else          : original_weight = rw.weights ( original )

        for c in comparators :
            
            tv , pv = c.pvalue ( data1 = target            ,
                                 data2 = original          , 
                                 weight2 = original_weight ) 
            
            header , row = c.the_row ()

            rw_type = '' if rw is None else typename ( rw ) 
            row     = ( rw_type , c.method ) + tuple ( row )
            rows.append ( row ) 

    header =  ( 'Reweighter' , 'Comparator' ) + header
    rows   = [ header ] + rows 
    title = 'Test results : %s-values [%%]' % script_p
    table = T.table ( rows , title = title , prefix = '# ' )
    logger.info ( '%s\n%s' % ( title , table ) )
  
# =============================================================================
if '__main__' == __name__ :

    run_reweight ( 3 , 1000 ) 
    
# =============================================================================

        
# =============================================================================
##                                                                      The END 
# =============================================================================
