#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/pvalues.py
#  Utility to combine several p-values into single p-value 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2026-07-18
# =============================================================================
""" Utility to combine several p-values into single p-value 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-06-07"
__all__     = (
    'combine_pvalues' , ## Combine several p-values into a single p-value with uncertainty handling
    'combine_methods' , ## all known valid combine methdos    
)
# =============================================================================
from   packaging.version  import Version
from   scipy.stats        import ( combine_pvalues as _combine_pvs_ ,
                                   chi2            as _chi2_        ,
                                   cauchy          as _cauchy_      ,
                                   norm            as _norm_        ) 
from   ostap.math.ve      import VE
from   ostap.stats.utils  import clip_pvalue 
import math, numpy, scipy
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
logger = getLogger ( 'ostap.stats.pvalues' )
# ============================================================================
## Combine several p-values into a single p-value 
# ============================================================================
SINGLE_P_METHODS      = ( "bonferroni" , "tippett"  , 'sidak' )
SCIPY_COMBINE_METHODS = ( "fisher"     , "stouffer" )
if Version ( "1.6" ) <= Version ( scipy.__version__ ) :
    SCIPY_COMBINE_METHODS = SCIPY_COMBINE_METHODS + ( "tippett" , "mudholkar_george" )
# ============================================================================
EXTRA_COMBINE_METHODS = ( "bonferroni"    ,
                          "simes"         ,
                          ## "pearson"    ,
                          "harmonic_mean" ,
                          "acat"          )
# ============================================================================
## all known valid combine methods 
COMBINE_METHODS = frozenset ( SCIPY_COMBINE_METHODS + EXTRA_COMBINE_METHODS )
## all known valid combine methods 
combine_methods = tuple ( sorted ( COMBINE_METHODS ) ) 
# =============================================================================
## Combine single-minimum p-value directly using VE error propagation
#  @param p_min_ve minimum p-value wrapped as VE
#  @param k        number of combined p-values
#  @param method   combination method name ("bonferroni" or "tippett")
def _combine_single_p_ve ( p_min_ve , k , method ) :
    """ Combine single-minimum p-value directly using VE error propagation
    """
    if "bonferroni" == method :
        # Bonferroni upper bound clamped at 1.0; VE propagates uncertainty
        pv    = VE  ( p_min_ve )
        kpv   = k * pv
        c2kpv = min ( max ( kpv.cov2() , 0 ) , 1 ) 
        return min ( VE ( 1.0 , c2kpv ) , kpv )

    if method in ( "tippett" , 'sidak' ) :
        # Tippett / Sidak combination: 1 - (1 - p_min)^k; VE propagates uncertainty
        one = VE ( 1.0 , 0.0 )
        return one - ( one - p_min_ve ) ** k

    raise ValueError ( "Not a single-p method: %s" % repr ( method ) )

# =============================================================================
## Compute combined p-value at a single point (no toy sampling)
#  @param p_arr  1D numpy array of p-values
#  @param method combination method name
def _combine_pvalues_point ( p_arr , method ) :
    """ Compute combined p-value at a single point (no toy sampling)
    """
    k = p_arr.size
    if 0 == k : raise ValueError ( "Empty sequence of p-values" )

    if "bonferroni" == method :
        # Multiply minimum p-value by number of tests, capped at 1.0
        return float ( min ( 1.0 , k * float ( numpy.min ( p_arr ) ) ) )

    if "harmonic_mean" == method :
        # Wilson harmonic mean p-value approximation
        hmp = float ( k ) / float ( numpy.sum ( 1.0 / p_arr ) )
        return float ( min ( 1.0 , hmp ) )

    if "acat" == method :
        # Aggregated Cauchy Association Test (Liu et al., 2019)
        stat = float ( numpy.mean ( numpy.tan ( ( 0.5 - p_arr ) * math.pi ) ) )
        return float ( _cauchy_.sf ( stat ) )

    if "simes" == method :
        # Simes global test across ordered p-values
        p_sorted = numpy.sort   ( p_arr )
        idx      = numpy.arange ( 1 , k + 1 , dtype = numpy.float64 )
        return float ( min ( 1.0 , float ( numpy.min ( k * p_sorted / idx ) ) ) )

    if "pearson" == method :
        # Pearson chi-square test on (1 - p) (more sensitive to large p-values)
        eps         = numpy.finfo ( float ).tiny
        one_minus_p = numpy.clip  ( 1.0 - p_arr , eps , 1.0 )
        stat        = -2.0 * float ( numpy.sum ( numpy.log ( one_minus_p ) ) )
        return float ( _chi2_.sf ( stat , df = 2.0 * k ) )

    # Delegate built-in methods to scipy.stats.combine_pvalues
    if Version ( "1.10" ) <= Version ( scipy.__version__ ) :
        return float ( _combine_pvs_ ( p_arr , method = method ).pvalue )
    
    return float ( _combine_pvs_ ( p_arr , method = method ) [ 1 ] )

# =============================================================================
## Propagate input p-value uncertainties via vectorized toys
#  @param p_vals array of nominal p-values
#  @param p_stds array of input standard deviations (sqrt of variance)
#  @param method combination method name
#  @param tol    boundary clipping tolerance away from 0 and 1
#  @param N      number of toy draws
#  @param rng    numpy random generator instance
def _combine_pvalues_toys_vec ( p_vals , p_stds , method , tol , N , rng ) :
    """ Propagate input p-value uncertainties via vectorized toys
    """
    k       = p_vals.size
    lo , hi = tol , 1.0 - tol
    has_err = p_stds > 0.0

    # No non-zero input uncertainties -> zero output variance
    if not numpy.any ( has_err ) : return 0.0

    # Sample perturbations only for entries with non-zero error
    samples = numpy.broadcast_to ( p_vals , ( N , k ) ).copy ()
    idx_err = numpy.flatnonzero ( has_err )
    draws   = rng.normal ( loc   = p_vals [ idx_err ] ,
                          scale = p_stds [ idx_err ] ,
                          size  = ( N , idx_err.size ) )
    samples [ : , idx_err ] = numpy.clip ( draws , lo , hi )

    if "fisher" == method :
        
        # Fisher chi-square combination across toys
        stat  = -2.0 * numpy.sum ( numpy.log ( samples ) , axis = 1 )
        toy_p = _chi2_.sf ( stat , df = 2.0 * k )
        
    elif "pearson" == method :
        # Pearson inverse-Fisher combination across toys
        one_m = numpy.clip ( 1.0 - samples , numpy.finfo ( float ).tiny , 1.0 )
        stat  = -2.0 * numpy.sum ( numpy.log ( one_m ) , axis = 1 )
        toy_p = _chi2_.sf ( stat , df = 2.0 * k )
        
    elif "harmonic_mean" == method :
        
        # Harmonic mean p-value across toys
        toy_p = numpy.minimum ( 1.0 , k / numpy.sum ( 1.0 / samples , axis = 1 ) )
        
    elif "acat" == method :
        
        # Cauchy combination statistic across toys
        stat  = numpy.mean ( numpy.tan ( ( 0.5 - samples ) * math.pi ) , axis = 1 )
        toy_p = _cauchy_.sf ( stat )
        
    elif "simes" == method :
        
        # Simes ordered p-value test across toys
        sorted_s = numpy.sort ( samples , axis = 1 )
        denom    = numpy.arange ( 1 , k + 1 , dtype = numpy.float64 ) [ None , : ]
        toy_p    = numpy.minimum ( 1.0 , numpy.min ( k * sorted_s / denom , axis = 1 ) )
        
    elif "stouffer" == method :
        
        # Stouffer Z-score combination across toys
        z      = _norm_.ppf ( 1.0 - samples )
        z_comb = numpy.sum ( z , axis = 1 ) / math.sqrt ( k )
        toy_p  = _norm_.sf ( z_comb )
        
    elif "mudholkar_george" == method :
        
        # Mudholkar-George logit sum approximated by Student-t distribution
        from scipy.stats import t as student_t
        logits = numpy.log ( samples / ( 1.0 - samples ) )
        scale  = math.sqrt ( 3.0 * ( 5.0 * k + 4.0 ) / ( k * ( 5.0 * k + 2.0 ) * ( math.pi ** 2 ) ) )
        stat   = -numpy.sum ( logits , axis = 1 ) * scale
        toy_p  = student_t.sf ( stat , df = 5 * k + 4 )
        
    else :
        
        # Fallback for any remaining method
        toy_p = numpy.array ( [ _combine_pvalues_point ( row , method = method ) for row in samples ] ,
                              dtype = numpy.float64 )

    # Return variance (rms^2) of combined p-values across toys
    return float ( numpy.var ( toy_p , ddof = 0 ) )

# =============================================================================
## Combine several p-values into a single p-value with uncertainty handling
#  @param pvalues      sequence of p-values (floats or VE instances)
#  @param method       combination method name
#  @param tol          boundary tolerance away from 0 and 1
#  @param N            number of toy draws for multi-p VE propagation
#  @param random_state seed or generator for toy sampling
def combine_pvalues ( pvalues      ,
                      method       ,
                      tol          = 1.e-9 ,
                      N            = 400   ,
                      random_state = None  ) :
    """ Combine several p-values into a single p-value with uncertainty handling
    """
    ## 
    method   = method.strip().lower()
    if not method in COMBINE_METHODS : raise ValueError ( "Unknown combination method: `%s'" % method )  
    if 'sidak' == method : method = 'tippett'


    if any ( not 0.0 <= float ( p ) <= 1.0 for p in pvalues ) :
        raise ValueError ( "Invalid input p-value outside [0,1] !" )  

    raw_list = list ( pvalues )
    k        = len ( raw_list )
    if 0 == k : raise ValueError ( "Empty sequence of p-values" )

    # Check if any input is a VE with non-zero variance inside [0, 1]
    has_ve = any ( isinstance ( p , VE ) and 0.0 < p.cov2 () and 0.0 <= float ( p ) <= 1.0
                   for p in raw_list )

        
    if has_ve :
        # Clip input VEs into safe probability range before combining
        clipped_ves = [
            clip_pvalue ( VE ( p ) if isinstance ( p , VE ) else VE ( float ( p ) , 0.0 ) , 0.5 )
            for p in raw_list
        ]

        # Direct VE path for methods depending solely on minimum p-value
        if method in SINGLE_P_METHODS :
            p_min_ve = min ( clipped_ves , key = float )
            return _combine_single_p_ve ( p_min_ve , k = k , method = method )

        # Multi-p path: nominal point value plus variance from vectorized toys
        vals     = numpy.array ( [ min ( max ( tol , float ( v ) ) , 1.0 - tol ) for v in clipped_ves ] ,
                                 dtype = numpy.float64 )
        vars_in  = numpy.array ( [ max ( 0.0 , float ( v.cov2 () ) ) for v in clipped_ves ] ,
                                 dtype = numpy.float64 )
        mean_val = _combine_pvalues_point ( vals , method = method )
        rng      = numpy.random.default_rng ( random_state )
        var_out  = _combine_pvalues_toys_vec ( p_vals = vals ,
                                              p_stds = numpy.sqrt ( vars_in ) ,
                                              method = method ,
                                              tol    = tol ,
                                              N      = N ,
                                              rng    = rng )
        return VE ( mean_val , var_out )

    # Plain float path without error propagation
    p_arr = numpy.array ( [ min ( max ( tol , float ( p ) ) , 1.0 - tol ) for p in raw_list ] ,
                          dtype = numpy.float64 )
    
    return _combine_pvalues_point ( p_arr , method = method )


# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
        

# =============================================================================
##                                                                      The END 
# =============================================================================
