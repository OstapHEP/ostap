#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/polynomials.py
#  Module with some useful utilities for dealing with polynonials
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-12-01"
__all__     = ()
# =============================================================================
import  ostap.math.base  
import  ostap.math.reduce
import  ROOT
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.polynomials' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================


# =============================================================================
## (Redefine standard constructor to allow usage of python lists&tuples)
#  Lists and tuples are converted on flight to :
# - std::vector<double> 
def _new_init_ ( t ,  *args )  :
    """(Redefine standard constructor to allow usage of python lists&tuples)
    Lists and tuples are  converted on flight to :
    - std::vector<double> 
    """
    from ostap.math.base        import doubles  , VCT_TYPES 
    from ostap.core.ostap_types import Generator, Sequence, list_types  
    
    largs = list (  args )

    for i , arg in enumerate ( largs ) :
        
        if   isinstance ( arg , VCT_TYPES  ) : continue 
        
        if   isinstance ( arg , Generator  ) : pass
        elif isinstance ( arg , Sequence   ) : pass
        elif isinstance ( arg , list_types ) : pass
        else :
            continue
        
        try: 
            _arg = doubles  ( arg  )
            largs [ i ] = _arg
            continue 
        except TypeError : pass
                
    targs = tuple ( largs )
        
    ## use old constructor 
    return t._old_init_ ( *targs ) 

# =============================================================================
for p in ( ROOT.Ostap.Math.Polynomial    , 
           ROOT.Ostap.Math.ChebyshevSum  , 
           ROOT.Ostap.Math.LegendreSum   , 
           ROOT.Ostap.Math.HermiteSum    ) :
    
    if not hasattr ( p , '_old_init_' ) :
        
        p._old_init_ = p.__init__
        ## Modifed constructor to allow python lists/tuples
        def _p_new_init_ ( s ,  *args ) :
            """Modifed constructor to allow python lists/tuples
            """
            _new_init_ ( s , *args )
            
        _p_new_init_.__doc__ += '\n' +   _new_init_.__doc__ 
        _p_new_init_.__doc__ += '\n' + p._old_init_.__doc__ 
        p.__init__ = _p_new_init_ 

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

