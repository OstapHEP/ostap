#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/karlin.py
#  Implemenetation of positive polynomials
#  - A non-negative polynomial on the [a,b] interval 
#  @see S.Karlin, L.S.Shapley, "Geometry of Moment Space", 
#       Mem. Am. Math. Soc., 12, 1953 
#  - A non-negative polynomial on the [ x_0,+infinityL interval
#  @see S.Karlin, W.J.Studden, "Tchebysheff systems: with application
#       in analysis and statistics", 
#
# @see https://bookstore.ams.org/view?ProductCode=MEMO/1/12 
# @see https://www.scirp.org/(S(351jmbntvnsjt1aadkposzje))/reference/ReferencesPapers.aspx?ReferenceID=1762140
# @see Ostap::Math::KarlinShapley
# @see Ostap::Math::KarlinStudden
#
# @author Vanya BELYAEV Ivan.Belayev@itep.ru
# @date   2023-01-16
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-12-01"
__all__     = (
    )
# =============================================================================
from ostap.core.core      import Ostap
from ostap.math.bernstein import _p_set_par_, _p_get_par_, _p_new_init_ 
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.karlin' )
else                       : logger = getLogger ( __name__            )
# =============================================================================
## Get parameters for <code>Ostap::Math::KarlinShapley</code> and
#                     <code>Ostap::Math::KarlinStudden</code> objects
# @code
# poly = ...
# pars = poly.pars() 
# @endcode 
# @see Ostap::Math::KarlinShapley
# @see Ostap::Math::KarlinStudden
def _ks_pars_ ( ks ) :
    """Get parameters for `Ostap.Math.KarlinShapley`
    and `Ostap.Math.KarlinStudden` objects  
    >>> poly = ...
    >>> pars = poly.pars() 
    - see `Ostap.Math.KarlinShapley`
    - see `Ostap.Math.KarlinStudden`
    """
    np = ks.npars()
    return tuple ( ks.par ( i ) for i in range ( np ) )

# =============================================================================
## Get phase parameters for <code>Ostap::Math::KarlinShapley</code> and
#                     <code>Ostap::Math::KarlinStudden</code> objects
# @code
# poly   = ...
# phases = poly.phases() 
# @endcode 
# @see Ostap::Math::KarlinShapley
# @see Ostap::Math::KarlinStudden
# @see Ostap::Math::KarlinShapley::phases1
# @see Ostap::Math::KarlinStudden::phases1
# @see Ostap::Math::KarlinShapley::phases2
# @see Ostap::Math::KarlinStudden::phases2
def _ks_phases_ ( ks ) :
    """Get phase parameters for `Ostap.Math.KarlinShapley`
    and `Ostap.Math.KarlinStudden` objects  
    >>> poly   = ...
    >>> phases = poly.phases() 
    - see `Ostap.Math.KarlinShapley`
    - see `Ostap.Math.KarlinStudden`
    - see `Ostap.Math.KarlinShapley.phases1`
    - see `Ostap.Math.KarlinStudden.phases1`
    - see `Ostap.Math.KarlinShapley.phases2`
    - see `Ostap.Math.KarlinStudden.phases2`
    """
    np = ks.npars()
    return tuple ( ks.par ( i ) for i in range ( 1 , np ) )

## decorate them 
for p in  ( Ostap.Math.KarlinShapley ,
            Ostap.Math.KarlinStudden ) :
    
    p.pars        = _ks_pars_
    p.phases      = _ks_phases_
    p.__getitem__ = _p_get_par_
    p.__setitem__ = _p_set_par_
    
    if not hasattr ( p , '_old_init_' ) :
        p._old_init_ = p.__init__
        ## Modifed constructor to allow python lists/tuples
        def _pp_new_init_ ( s ,  *args ) :
            """Modifed constructor to allow python lists/tuples
            """
            _p_new_init_ ( s , *args )
            
        _pp_new_init_.__doc__ += '\n' + _p_new_init_.__doc__ 
        _pp_new_init_.__doc__ += '\n' + p._old_init_.__doc__ 
        p.__init__ = _pp_new_init_ 
        
# =============================================================================
_decorated_classes_ = ( Ostap.Math.KarlinShapley ,
                        Ostap.Math.KarlinStudden )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
