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
from    ostap.core.ostap_types import num_types 
from    ostap.math.base        import Ostap 
import  ostap.math.reduce
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.polynomials' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================


for poly  in ( Ostap.Math.Polynomial, ) : 
    
    def _poly_add_  ( self , other  ) : 
        """ Modified __add__ 
        """
        if   isinstance ( other , num_types ) : return self.add ( other ) 
        elif isinstance ( other , poly      ) : return self.add ( other )
        
        return NotImplemented 
    
    def _poly_mul_  ( self , other  ) : 
        """ Modified __mul__ 
        """
        if   isinstance ( other , num_types ) : return self.add ( other ) 
        ## elif isinstance ( other , poly      ) : return self.add ( other )
                
        return NotImplemented 
    
    def _poly_sub_  ( self , other  ) : 
        """ Modified __sub__ 
        """
        if   isinstance ( other , num_types ) : return self.sub ( other ) 
        elif isinstance ( other , poly      ) : return self.sub ( other )

        return NotImplemented 
        
    def _poly_div_  ( self , other  ) : 
        """ Modified __div__ 
        """
        if   isinstance ( other , num_types ) : 
            if not other  : raise ZeroDivisionError()
            return self.div ( other ) 
        
        return NotImplemented 
        
    def _poly_rsub_  ( self , other  ) : 
        """ Modified __rsub__ 
        """
        if   isinstance ( other , num_types ) : return self.rsub ( other ) 
        
        return NotImplemented 
        
    ## redefine C++ __add__ operation 
    poly.__add__      = _poly_add_     
    poly.__mul__      = _poly_mul_     
    poly.__sub__      = _poly_sub_ 
    poly.__div__      = _poly_div_ 
    poly.__truediv__  = _poly_div_     
    poly.__rsub__     = _poly_rsub_
    
    poly.__rdiv__     = lambda x,*s : NotImplemented
    poly.__rtruediv__ = lambda x,*s : NotImplemented

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

