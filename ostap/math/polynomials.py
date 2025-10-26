#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/polynomials.py
#  Redefine C++ operator for polynoial functions
#  - essentially we need to return NotImplemnted when TypeError occurs 
# =============================================================================
""" Redefine C++ operator for polynoial functions   
- redefine C++ opoerators, 
- enhance them 
- return `NotImplemented` if/when TypeError occurs 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-12-01"
__all__     = ()
# =============================================================================
from    ostap.core.ostap_types import num_types 
from    ostap.math.base        import Ostap
from    ostap.stats.funstats   import FunMNMX
import  ostap.math.reduce
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.polynomials' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================
## Modified iadd
def _poly_iadd_  ( poly , other  ) : 
    """ Modified __iadd__ 
    """
    try : return poly.iadd ( other )   
    except TypeError : pass     
    return NotImplemented 
# =============================================================================
## Modified isub
def _poly_isub_  ( poly , other  ) : 
    """ Modified __isub__ 
    """ 
    try : return poly.isub ( other )   
    except TypeError : pass     
    return NotImplemented 
# =============================================================================
## Modified imul
def _poly_imul_  ( poly , other  ) : 
    """ Modified __imul__ 
    """
    try : return poly.imul ( other ) 
    except TypeError : pass     
    return NotImplemented 
# =============================================================================
## Modified idiv
def _poly_idiv_  ( poly , other  ) : 
    """ Modified __idiv__ 
    """ 
    if isinstance ( other , num_types ) and not other  : 
        raise ZeroDivisionError ( "One cannot divide by zero!" )
    try: return poly.idiv ( other )   
    except TypeError : pass     
    return NotImplemented 
# =============================================================================
## Modified ipow
def _poly_ipow_  ( poly , other  ) : 
    """ Modified __ipow__ 
    """    
    try :
        if isinstance ( other , int ) and 0 <= other : return poly.ipow ( other ) 
    except TypeError : pass     
    return NotImplemented 
# ==============================================================================
## Modified add
def _poly_add_  ( poly , other  ) : 
    """ Modified __add__ 
    """
    try : return poly.add ( other )   
    except TypeError : pass     
    return NotImplemented 
# =============================================================================
## Modified isub
def _poly_sub_  ( poly , other  ) : 
    """ Modified __sub__ 
    """ 
    try : return poly.sub ( other )   
    except TypeError : pass     
    return NotImplemented 
# =============================================================================
## Modified mul
def _poly_mul_  ( poly , other  ) : 
    """ Modified __mul__ 
    """
    try : return poly.mul ( poly ) 
    except TypeError : pass     
    return NotImplemented 
# =============================================================================
def _poly_div_  ( poly , other  ) : 
    """ Modified __div__ 
    """ 
    if isinstance ( other ,  num_types ) and not other  : 
        raise ZeroDivisionError ( "One cannot divide by zero!" )
    try: return poly.div ( other )   
    except TypeError : pass     
    return NotImplemented 
# =============================================================================
## Modified pow 
def _poly_pow_  ( poly , other  ) : 
    """ Modified __pow__ 
    """
    try :
        if isinstance ( other , int ) and 0 <= other : return poly.pow ( other ) 
    except TypeError : pass     
    return NotImplemented 
# =============================================================================
## rigth subtraction
def _poly_rsub_  ( poly , other  ) : 
    """ Modified __rsub__ 
    """ 
    try: return poly.rsub ( other )   
    except TypeError :pass     
    return NotImplemented 

# =============================================================================
## Convert to "density" if/when possible
def _poly_density_( poly , *limits ) :
    """ Convert to `density' if/when possible
    """
    ## get the integral
    i = poly.integral ( *limits ) 
    
    if i <= 0 or not i : 
        logger.error ( "density: integral is non-positive, skip")
        return poly  
    
    return poly / i 

# =============================================================================
## copy/clone polynomials
def _poly_copy_  ( poly ) :
    """ Copy/clone polynomial (using copy contructor)  
    """
    ptype = type ( poly ) 
    return ptype ( poly )

# ============================================================================
## estimate min/max values of polynomials
#  @attentino - these are not real minmax, 
#  but just such pair of values that gurantes 
#  \f$ v_{min} \le p(x) \le v_{max} \f$
#     min <= Ppoly <= max 
# The etimate can be rather crude
def _poly_minmax_  ( poly , raw = True ) :
    """ Estimate min/max values of polynomials
    - these are not real minmax, but just 
    - such pair of values that gurantees min <= poly <= max
    """
    if isinstance   ( poly , Ostap.Math.Monotonic ) : 
        
        parss = poly.bernstein ().pars()  
        vmin , vmax = min ( pars ) , max ( pars )
        ## the result is exact 
        return vmin , vmax
    
    elif isinstance ( poly , Ostap.Math.Bernstein ) :
        
        pars = poly.pars() 
        vmin , vmax = min ( pars ) , max ( pars )
        if raw or 1 >= poly.degree () : return vmin , vmax
    
    elif hasattr  ( poly , 'bernstein' )  : 
        return _poly_minmax_  ( poly.bernstein () , raw = raw )
    
    elif raw : 
        b = Ostap.Math.Bernstein ( poly ) 
        return _poly_minmax_ ( b , raw = raw ) 

    mnmx = FunMNMX ( poly.xmin () , poly.xmax () , N = poly.degree() * 2 )   
    return mnmx ( poly )

# ============================================================================
## (re)decorate polynomials
# ============================================================================
for poly  in ( Ostap.Math.Polynomial   , 
               Ostap.Math.ChebyshevSum , 
               Ostap.Math.LegendreSum  , 
               Ostap.Math.Bernstein    ) : 
        
    if hasattr ( poly , 'iadd' ) : poly.__iadd__ = _poly_iadd_ 
    if hasattr ( poly , 'isub' ) : poly.__isub__ = _poly_isub_ 
    if hasattr ( poly , 'imul' ) : poly.__imul__ = _poly_imul_ 
    if hasattr ( poly , 'idiv' ) : 
        poly.__idiv__     = _poly_idiv_ 
        poly.__itruediv__ = _poly_idiv_ 
    
    if hasattr ( poly , 'add' ) :  
        poly. __add__ = _poly_add_ 
        poly.__radd__ = _poly_add_ 
        
    if hasattr ( poly , 'sub' ) :  poly.__sub__  = _poly_sub_ 
    
    if hasattr ( poly , 'mul' ) : 
        poly. __mul__ = _poly_mul_ 
        poly.__rmul__ = _poly_mul_
         
    if hasattr ( poly , 'div' ) : 
        poly. __div__    = _poly_div_ 
        poly.__truediv__ = _poly_div_ 

    if hasattr ( poly , 'rsub' ) :  poly.__rsub__  = _poly_rsub_
    
    if hasattr ( poly , 'ipow' ) :  poly.__ipow__  = _poly_ipow_ 
    if hasattr ( poly , 'pow'  ) :  poly.__pow__   = _poly_pow_ 

    if not hasattr ( poly , 'density'  ) : poly.density  = _poly_density_ 
    if not hasattr ( poly , 'clone'    ) : poly.clone    = _poly_copy_ 
    if not hasattr ( poly , 'copy'     ) : poly.copy     = _poly_copy_ 
    if not hasattr ( poly , '__copy__' ) : poly.__copy__ = _poly_copy_ 
    if not hasattr ( poly , 'minmax'   ) : poly.minimax  = _poly_minmax_ 

Ostap.Math.Monotonic.minmax  = _poly_minmax_         
# =============================================================================
## decorated classes 
_decorated_classes_  = ()

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================

