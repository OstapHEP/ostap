#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
# @file
# Test module for binomial intervals 
# ============================================================================= 
""" Test module binomial intervals 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import 
# ============================================================================= 
from   ostap.core.core        import Ostap 
from   ostap.utils.root_utils import batch_env
from   ostap.logger.symbols   import efficiency as eff_symbol 
import ostap.logger.table     as     T 
import ROOT, random, math 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'ostap.test_stats_binomial' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
logger.info ( 'Test for binomial intervals ')
# =============================================================================
## set batch form environment 
batch_env ( logger )
# =============================================================================
#  the width of the +-1 sigma confidence interval 
one_sigma   = Ostap.Math.gauss_cdf ( 1.0 ) - Ostap.Math.gauss_cdf ( -1.0 )

# =============================================================================
## test binomial intervals
def test_binomial () :
    """ Test obinomial intervals
    """

    NN = 12
    
    for i in range ( 50 ) :

        accepted = random.randint ( 0 , NN )
        rejected = random.randint ( 0 , NN )
        
        total    = accepted + rejected
        
        eff  = 50. if not total else accepted * 100. / total

        A , B = accepted , rejected
        
        rows = [ ( 'Method' , 'A[%]' , 'B[%' ) ]
        
        a , b = Ostap.Math.wald_interval ( A , B , one_sigma )
        a , b = 100 * a , 100 * b 
        row   = 'Wald' , '%.2f' %a , '%.2f' % b
        rows.append ( row )
        
        a , b = Ostap.Math.wilson_score_interval ( A , B , one_sigma )
        a , b = 100 * a , 100 * b         
        row  = 'Wilson score' , '%.2f' %a , '%.2f' % b
        rows.append ( row )

        a , b = Ostap.Math.wilson_score_continuity_interval ( A , B , one_sigma ) 
        a , b = 100 * a , 100 * b 
        row  = 'Wilson score continuity' , '%.2f' %a , '%.2f' % b
        rows.append ( row )

        a , b = Ostap.Math.arcsin_interval ( A , B , one_sigma ) 
        a , b = 100 * a , 100 * b 
        row  = 'Arcsin' , '%.2f' %a , '%.2f' % b
        rows.append ( row )
        
        a , b = Ostap.Math.agresti_coull_interval ( A , B , one_sigma ) 
        a , b = 100 * a , 100 * b 
        row  = 'Agresti-Coull' , '%.2f' %a , '%.2f' % b
        rows.append ( row )

        a , b = Ostap.Math.jeffreys_interval ( A , B , one_sigma ) 
        a , b = 100 * a , 100 * b 
        row  = 'Jeffreys' , '%.2f' %a , '%.2f' % b
        rows.append ( row )

        a , b = Ostap.Math.clopper_pearson_interval ( A , B , one_sigma ) 
        a , b = 100 * a , 100 * b 
        row  = 'Clopper-Pearson' , '%.2f' %a , '%.2f' % b
        rows.append ( row )

        a , b = Ostap.Math.bayes_interval ( A , B , one_sigma ) 
        a , b = 100 * a , 100 * b 
        row  = 'Bayes' , '%.2f' %a , '%.2f' % b
        rows.append ( row )
        
        title = "(%d) Accepted/rejected=%d/%d, %s=%.2f%% @%.1f%%CL" % ( i , accepted , rejected , eff_symbol , eff , one_sigma * 100 ) 
        table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lcc' )
        logger.info ( '%s:\n%s' % ( title , table ) ) 
    
        
        

    
# =============================================================================
if '__main__' == __name__ :

    test_binomial ()  
    
# =============================================================================
##                                                                      The END 
# =============================================================================
