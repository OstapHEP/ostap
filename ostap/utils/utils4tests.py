#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
# @file ostap/utils/utils4tests.py
# Collection of simple utilities for Ostap tests
# ============================================================================= 
""" Collection of simple utilities for Ostap tests
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = (
    'batch_env'   , ## test in batch mode
    'light_tests' , ## Should we run "light" version of tests ?
    'accept'      , ## accept tiwh some prbabilty? 
    'active_test' , ## should we run this test ? 
)
# =============================================================================
from   ostap.utils.core       import typename 
from   ostap.utils.basic      import numcpu , is_ci
from   ostap.utils.root_utils import batch_env
import ostap.utils.env        as     E 
import random 
# =============================================================================
## cache soem local information
_light_tests_ = is_ci  () or ( ( v := E.get_env ( E.OSTAP_LIGHT_TESTS , '' ).strip() ) and v.lower() in E.boolean_true_values ) 
_ncpus_       = numcpu ()

# ==============================================================================
## Should we run "light" version of the tests ?
#  - Continious integration (GitHub, GitLab, Azurre, ...)
#  - check OSTAP_LIGHT_TESTS  environment varibale
#  - number of CPUS i excess of mincpu
#  @code
#  light , weak = light_tests () 
#  @endcode
def light_tests ( mincpu = 8 ) :
    """ Should we run "light" version of tests ?
    - Continious integration (GitHub, GitLab, Azurre, ...)
    - check OSTAP_LGHT_TESTS  environment varibale
    - number of CPUS i excess of mincpu

    >>> light , weak = light_tests () 
    """
    light = _light_tests_
    weak  = _ncpus_ <= mincpu
    
    return light or weak , weak 

# =============================================================================
## Get randm boolean according to "accept_fraction 
def accept ( fraction = 0.10 ) :
    """ Get randm boolean according to "accept_fraction
    """
    if  isinstance ( fraction , float ) :
        
        if   1 == fraction    : True
        elif 0 == fraction    : False
        elif 0 < fraction < 1 : return random.uniform ( 0 , 1 ) < fraction 

    elif isinstance ( fraction , int ) :
        
        if   1 == fraction    : True
        elif 0 == fraction    : False
        
    else : raise typeError ( "Invalid type for `fraction' : %s" % typename ( fraction ) )
    
    raise ValueError ( "Invalid Value for `fraction' : %s" % fraction ) 


# =============================================================================    
## Should we run this prescaled test ?
def active_test ( fraction = 0.10 , mincpu = 8 ) :
    """ Should we run this prescaled test ? 
    """
    light , _ = light_tests  ( mincpu )
    if not light  : return True
    return accept ( fraction )

  
# =============================================================================
if '__main__' == __name__ :
    
    # =========================================================================
    from ostap.logger.logger import getLogger
    logger = getLogger ( 'ostap.utils.utils4test' )
    # =========================================================================
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    

    light, weak = light_tests ( 8 )
    logger.info  ( 'LIGHT tests   : %s' % light )
    logger.info  ( 'WEAK  machine : %s' % weak  )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
