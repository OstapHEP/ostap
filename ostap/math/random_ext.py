#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
# =============================================================================
# @file random_ext.py
# The simple extention for the standard python module random
# @author Vanya BELYAEV
# @date   2012-04-28
# =============================================================================
""" The simple extension for the standard python module random
- bifurcated gaussian
- gaussian using Ostap.Math.ValueWithError as  argument
- poisson (missing in python random module)
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__version__ = "$Revision$"
__date__    = "2012-04-28"
# =============================================================================
__all__ = (
    'bifur'   ,  ## bifurcated gaussian
    've_gau'  ,  ## gaussian using ValueWithError construction
    'poisson' ,  ## poisson (missing in python random module) 
    )
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.math.random_ext')
else                      : logger = getLogger ( __name__               ) 

# =============================================================================
## generate bifurcated gaussian
#  @code
#  value = bifur ( 0 , -1 , +2 ) 
#  @endcode 
def _bifur_ ( self , mu , sigma1 , sigma2 ) :
    """Generate the bifurcated gaussian
    >>> value = bifur ( 0 , -1 , +2 )    
    """
    if   sigma1 * sigma2 > 0.0 :
        raise ValueError, 'Lower and upper errors must have opposite signs'
        
    _as1  =  abs ( float ( sigma1 ) ) 
    _as2  =  abs ( float ( sigma2 ) ) 
    _frac = _as1 / (  _as1  +  _as2 )
    
    _aux  =       self.random ()
    _gau  = abs ( self.gauss  ( 0 , 1 ) ) 
    
    if _aux <= _frac : return mu + sigma1 * _gau
    else             : return mu + sigma2 * _gau 

# =============================================================================
## generate bifurcated gaussian using Value
#  @see Ostap::Math::ValueWithError
def _ve_gau_ ( self , val ) :
    """Generate the gaussian accoridng to Ostap.Math.ValueWithError
    >>> ve    = VE ( 1 , 2 ) 
    >>> value = ve_gau ( ve  )
    """
    mean  = val.value ()
    sigma = val.error ()
    return self.gauss ( mean , sigma )

# =============================================================================
try :
    from scipy.random import poisson as _poisson
    def _poisson_ ( self , mu ) : return _poisson ( mu )
    logger.debug ('use scipy.random.poisson')
except ImportError :
    logger.debug ('scipy.random.poisson is not accessible, use hand-made replacement')
    _STEP = 500.0
    _MAX  =  30.0
    import math
    _sqrt  = math.sqrt
    _exp   = math.exp
    import ROOT,cppyy 
    _round = cppyy.gbl.Ostap.Math.round
    ## hand-made replacement for poisson random number generator  
    def _poisson_ ( self , mu ) :
        mu = float ( mu ) 
        if _MAX <= mu :
            r = _gauss ( mu , _sqrt( mu ) ) 
            return max ( _round ( r ) , 0 )
        x  = 0
        p  = _exp ( -mu )
        s  = p
        u  = self.uniform ( 0 , 1 )
        while s < u :
            x += 1
            p *= mu / x
            s += p
        return x

import random 
if not hasattr ( random.Random , 'bifur'   ) : random.Random.bifur   = _bifur_
if not hasattr ( random        , 'bifur'   ) : random.bifur          = random._inst.bifur
if not hasattr ( random.Random , 've_gau'  ) : random.Random.ve_gau  = _ve_gau_
if not hasattr ( random        , 've_gau'  ) : random.ve_gau         = random._inst.ve_gau
if not hasattr ( random.Random , 'poisson' ) : random.Random.poisson = _poisson_
if not hasattr ( random        , 'poisson' ) : random.poisson        = random._inst.poisson

bifur   = random.bifur
ve_gau  = random.ve_gau
poisson = random.poisson

# =============================================================================
if '__main__' == __name__  :
    
    from ostap.logger.line import line 
    logger.info ( __file__  + '\n' + line  ) 
    logger.info ( 80*'*'   )
    logger.info ( __doc__  )
    logger.info ( 80*'*' )
    logger.info ( ' Author  : %s' %         __author__    ) 
    logger.info ( ' Version : %s' %         __version__   ) 
    logger.info ( ' Date    : %s' %         __date__      )
    logger.info ( ' Symbols : %s' %  list ( __all__     ) )
    logger.info ( 80*'*' ) 

    from ostap.stats.counters import SE
    cnt = SE()
    mu  = 0.4
    for i in xrange(10000) :
        cnt += poisson(0.4)

    logger.info ( 'Poisson(mu=%s) : %s' % ( mu , cnt ) )
        
    logger.info ( 80*'*' ) 

    
# =============================================================================
# The END
# =============================================================================
