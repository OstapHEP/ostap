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
## ATTENTION! 
import os 
os.environ [ "OMP_NUM_THREADS"      ]  = "1"
os.environ [ "MKL_NUM_THREADS"      ]  = "1"
os.environ [ "OPENBLAS_NUM_THREADS" ]  = "1"
# =============================================================================
from   ostap.utils.root_utils   import batch_env 
from   ostap.logger.symbols     import script_p
from   ostap.utils.basic        import numcpu
from   ostap.logger.pretty      import nice_print
from   ostap.stats.utils        import nEff 
from   ostap.stats.tools        import ( hasLightGBM , hasXGBoost ,
                                         hasCatBoost , hasSkLearn ,
                                         hasPyTorch  , hasHepML   ) 
import ostap.logger.table       as     T 
import ostap.io.root_file 
import ostap.parallel.kisa
import os, numpy 
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
def make_datasets ( n_samples = 5000 ,
                    n_dim     =    3 ,
                    seed      = None ) :
    """ Temporarily simplified datasets for baseline validation:
        smooth distributions with moderate shifts and no sharp boundaries.
    """
    numpy.random.seed ( seed )
    
    # 1. Original dataset: standard normal distribution
    original = numpy.random.normal ( loc = 0.0 , scale = 1.0 , size = ( n_samples , n_dim ) )
    
    # 2. Target dataset: shifted and mildly correlated normal distribution
    target = numpy.random.normal ( loc = 0.3 , scale = 1.1 , size = ( n_samples , n_dim ) )
    
    if n_dim > 1 :
        # Mild smooth correlation instead of complex trigonometric/recursive functions
        target [ : , 1 ] = target [ : , 1 ] + 0.3 * ( target [ : , 0 ] ** 2 )
        
    if n_dim > 2 :
        target [ : , 2 ] = target [ : , 2 ] - 0.2 * target [ : , 0 ]

    return target , original


def make_datasets2 ( n_samples = 5000 ,
                    n_dim     =    3 ,
                    seed      = None ) :
    """ Generate datasets by oversampling, filtering strict boundaries, 
        adding a narrow Gaussian peak in the first dimension, 
        and mixing a small uniform background to avoid zero-density gaps.
    """
    numpy.random.seed ( seed )
    
    low_bound, high_bound = -3.0, 3.0
    
    # --- ORIGINAL: Uniform distribution within strict bounds ---
    original = numpy.zeros ( ( n_samples , n_dim ) )
    for d in range ( n_dim ) :
        original [ : , d ] = numpy.random.uniform ( low = low_bound , high = high_bound , size = n_samples )
        
    # --- TARGET: Oversample to filter out-of-bounds points naturally ---
    oversample_factor = 4
    n_large = n_samples * oversample_factor
    
    # Decide which events belong to the main structure vs the narrow peak (10%)
    n_peak = int ( 0.10 * n_large )
    n_main = n_large - n_peak
    
    # 1. Main structure generation
    x0_main = numpy.random.uniform ( low = low_bound , high = high_bound , size = n_main )
    
    # 2. Narrow Gaussian peak generation in the first dimension (width = 0.5)
    x0_peak = numpy.random.normal ( loc = 1.0 , scale = 0.5 , size = n_peak )
    
    # Combine x_0 components
    x0_large = numpy.concatenate ( [ x0_main , x0_peak ] )
    numpy.random.shuffle ( x0_large )
    
    target_raw = numpy.zeros ( ( n_large , n_dim ) )
    target_raw [ : , 0 ] = x0_large
    
    if n_dim > 1 :
        target_raw [ : , 1 ] = 0.5 * ( x0_large ** 2 ) - 2.0 + numpy.random.normal ( loc = 0.0 , scale = 0.2 , size = n_large )
        
    if n_dim > 2 :
        target_raw [ : , 2 ] = numpy.sin ( numpy.pi * x0_large / 3.0 ) + 0.3 * target_raw [ : , 1 ] + numpy.random.normal ( loc = 0.0 , scale = 0.15 , size = n_large )
        
    for d in range ( 3 , n_dim ) :
        target_raw [ : , d ] = 0.4 * target_raw [ : , d - 1 ] + 0.5 * numpy.cos ( x0_large ) + numpy.random.normal ( loc = 0.0 , scale = 0.2 , size = n_large )

    # Filter points strictly within [low_bound, high_bound] for all dimensions
    mask = numpy.all ( ( target_raw >= low_bound ) & ( target_raw <= high_bound ) , axis = 1 )
    target_filtered = target_raw [ mask ]
    
    # Trim or fallback to n_samples
    if len ( target_filtered ) >= n_samples :
        target = target_filtered [ : n_samples ]
    else :
        target = target_filtered
        
    # Add a small uniform background component (50%) to eliminate exact zero-density pockets
    n_bg = int ( 0.50 * len ( target ) )
    if n_bg > 0 and len ( target ) > n_bg :
        bg_data = numpy.random.uniform ( low = low_bound , high = high_bound , size = ( n_bg , n_dim ) )
        target [ : n_bg ] = bg_data
        
    return target , original


# ==============================================================================
## Compare datasets using several methods 
# ==============================================================================
from ostap.stats.gof_np       import  ( KullbackLeibler as COMPARATOR3 , 
                                        Hotelling       as COMPARATOR2 , 
                                        Mahalanobis     as COMPARATOR1 ) 

from ostap.stats.data_compare import data_compare     
comparators = ( COMPARATOR1 ( parallel = True , nToys = 100 ) ,
                COMPARATOR2 ( parallel = True , nToys = 100 ) ,
                COMPARATOR3 ( parallel = True , nToys = 100 ) ) 

if hasLightGBM () :
    from ostap.stats.adval        import ADVAL_LGBM  as COMPARATOR5
    comparators += ( COMPARATOR5 ( parallel = True , nToys = 25 ) , ) 

if hasXGBoost  () :  
    from ostap.stats.adval        import ADVAL_XGB  as COMPARATOR6
    comparators += ( COMPARATOR6 ( parallel = True , nToys = 25 ) , ) 

if hasCatBoost () :  
    from ostap.stats.adval        import ADVAL_CATB  as COMPARATOR7
    comparators += ( COMPARATOR7 ( parallel = True , nToys = 25 ) , ) 

if hasSkLearn (): 
    from ostap.stats.adval        import ADVAL_HGBC  as COMPARATOR8
    comparators += ( COMPARATOR8 ( parallel = True , nToys = 25 ) , )
    
    from ostap.stats.adval        import ADVAL_GBC   as COMPARATOR9
    comparators += ( COMPARATOR9 ( parallel = True , nToys = 25 ) , ) 

# ============================
def run_reweight ( n_dim     = 3     ,
                   n_samples = 10000 ) : 
    
    ## generate samples 
    target , original = make_datasets ( n_dim = n_dim , n_samples = n_samples) 
    
    reweighters = [ None ]
    
    if hasHepML    () :         
        from ostap.tools.reweighters  import GBReweighter   as GBRW
        rw1 = GBRW ( target = target , original = original )
        reweighters.append ( rw1 )
        
    if hasLightGBM () : 
        from ostap.tools.reweighters  import LightGBMDensityReweighter as  LGBM
        rw2 = LGBM ( target = target , original = original )
        reweighters.append ( rw2 )

    if hasXGBoost  () : 
        from ostap.tools.reweighters  import XGBoostDensityReweighter as  XGB 
        rw3 = XGB  ( target = target , original = original )
        reweighters.append ( rw3 )

    if hasCatBoost () : 
        from ostap.tools.reweighters  import CatBoostDensityReweighter as CATB
        rw4 = CATB ( target = target , original = original )
        reweighters.append ( rw4 )
        
    if hasPyTorch () : 
        from ostap.tools.reweighters  import PyTorchDensityReweighter as TORCH
        rw5 = TORCH ( target = target , original = original , n_splits = 1 )
        reweighters.append ( rw5 )

    header = []
    rows   = []     
    for rw in reweighters :

        if rw is None : original_weight = None 
        else          : original_weight = rw.weights ( original )

        wmin, wmax = '' , '' 
        if not original_weight is None :
            wmin = float ( numpy.min ( original_weight ) ) 
            wmax = float ( numpy.max ( original_weight ) )
            wmin , wmax = '%.4g' % wmin , '%.4g' % wmax 
                               
        for c in comparators :
            
            tv , pv = c.pvalue ( data1   = target          ,
                                 data2   = original        , 
                                 weight2 = original_weight ) 
            
            header , row = c.the_row ()

            neff    = nEff ( original , original_weight ) 
            rw_type = '' if rw is None else rw.method            
            neff    = nice_print ( neff , precision = 3 ) 
            row     = ( rw_type , c.method , neff ) + tuple ( row ) + ( wmin , wmax )             
            rows.append ( row )
                        
    header =  ( 'Reweighter' , 'Comparator' , '#eff' ) + header + ( 'w-min' , 'w-max' ) 
    rows   = [ header ] + rows 
    title = 'Test results : %s-values [%%]' % script_p
    table = T.table ( rows , title = title , prefix = '# ' )
    logger.info ( '%s\n%s' % ( title , table ) )
  
# =============================================================================
if '__main__' == __name__ :

    run_reweight ( 3 , 3000 ) 
    
# =============================================================================

        
# =============================================================================
##                                                                      The END 
# =============================================================================
