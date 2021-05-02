#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/bootstrap.py
#  Primitive bootstrap generator 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2021-05-01
# =============================================================================
"""Trivial bootstrap generator
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'bootstrap'  , ## primitive bootstrap generator 
    )
# =============================================================================
from   builtins          import range
# =============================================================================
try :
    
    import numpy as np
    
    # ============================================================================
    ## generate bootstrapping data
    #  @code
    #  data = ...
    #  for ds = bootstrap ( data , size = 100 ) :
    #  ...
    #  @endcode
    #  Generate indices for bootstrapped samples:
    #  @code
    #  data = ...
    #  N    = len ( data )
    #  for indices  in bootstrap ( range ( N ) , size = 100  ) :
    #      sample = data [ indices ]
    #  @endcode
    def bootstrap ( data , size = 100 ) :
        """Generate bootstrapping data
        >>> data = ...
        >>> for ds in bootstrap ( data , size = 100 ) :
        >>> ,,,
        - Generate indices for bootstrapped samples:
        >>> data = ...
        >>> N    = len ( data )
        >>> for indices  in bootstrap ( range ( N ) , size = 100  ) :
        >>> ... sample = data [ indices ]
        """
        N = len ( data )
        for i in range ( size ) :
            yield np.random.choice ( data , size = N )

except ImportError :

    import random 
    # ============================================================================
    ## generate bootstrapping data
    #  @code
    #  data = ...
    #  for ds = bootstrap ( data , size = 100 ) :
    #  ...
    #  @endcode
    #  Generate indices for bootstrapped samples:
    #  @code
    #  data = ...
    #  N    = len ( data )
    #  for indices  in bootstrap ( range ( N ) , size = 100  ) :
    #      sample = data [ indices ]
    #  @endcode
    def bootstrap ( data , size = 100 ) :
        """Generate bootstrapping data
        >>> data = ...
        >>> for ds in bootstrap ( data , size = 100 ) :
        >>> ,,,
        - Generate indices for bootstrapped samples:
        >>> data = ...
        >>> N    = len ( data )
        >>> for indices  in bootstrap ( range ( N ) , size = 100  ) :
        >>> ... sample = data [ indices ]
        """
        N = len ( data )
        for i in range ( size ) :
            yield tuple ( random.choices ( data , k = N ) ) 


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.logger.logger import getLogger 
    logger = getLogger( 'ostap.stats.bootstrap' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
