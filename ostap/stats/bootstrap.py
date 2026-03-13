#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/bootstrap.py
#  Primitive bootstrap generator 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2021-05-01
# =============================================================================
""" Trivial bootstrap generator
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'bootstrap'                   , ## primitive bootstrap generator 
    'extended_bootstrap'          , ## extedend  bootstrap generator 
    'bootstrap_indices'           , ## primitive bootstrap indices generator 
    'extended_bootstrap_indices'  , ## extedend  bootstrap indices generator 
    )
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import numpy as np
    # =========================================================================
    if hasattr ( np.random , 'default_rng' ) : np_rng = np.random.default_rng () 
    else                                     : np_rng = np.random
    # =========================================================================
    ## generate bootstrap samples 
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
    def bootstrap ( data , size = 100 , sort = False ) :
        """ Generate bootstrap samples
        >>> data = ...
        >>> for ds in bootstrap ( data , size = 100 ) :
        >>> ,,,
        - Generate indices for bootstrapped samples:
        >>> data = ...
        >>> N    = len ( data )
        >>> for indices  in bootstrap ( range ( N ) , size = 100  ) :
        >>> ... sample = data [ indices ]
        """
        N   = len ( data )
        for i in range ( size ) :
            indices = np_rng.choice ( data , size = N )
            if sort : indices.sort () 
            yield indices 

    # =========================================================================
    ## Generate "extended" bootstrap
    #  @code
    #  data = ...
    #  for ds = extended_bootstrap ( data , size = 25 ) :
    #  @endcode
    #  Generate indices for bootstrapped samples:
    #  @code
    #  data = ...
    #  N    = len ( data )
    #  for indices  in extended_bootstrap ( range ( N ) , size = 100  ) :
    #      sample = data [ indices ]
    #  @endcode        
    def extended_bootstrap ( data , size , sort = False  )  :
        """ Generate "extended" bootstrap
        >>> data = ...
        >>> for ds = extended_bootstrap ( data , size = 25 ) :
        - Generate indices for bootstrap samples:
        >>> data = ...
        >>> N    = len ( data )
        >>> for indices  in extended_bootstrap ( range ( N ) , size = 100  ) :
        >>> ... sample = data [ indices ]
        """
        N    = len ( data )
        for i in range ( size ) :
            indices = np_rng.choice ( data , size = np_rng.poisson ( N ) )
            if sort : indices.sort () 
            yield indices 
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    import random
    from   ostap.utils.utils import choices    
    # =========================================================================
    ## Generate bootstrap samples 
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
    def bootstrap ( data , size = 100 , sort = False ) :
        """ Generate bootstrap samples 
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
            indices = choices ( data , k = N )
            if sort : indices = sorted ( indices ) 
            yield tuple ( indices ) 

    from ostap.math.random_ext import poisson

    # =========================================================================
    ## Generate "extended" bootstrap
    #  @code
    #  data = ...
    #  for ds = extended_bootstrap ( data , size = 25 ) :
    #  @endcode
    #  Generate indices for bootstrapped samples:
    #  @code
    #  data = ...
    #  N    = len ( data )
    #  for indices  in extended_bootstrap ( range ( N ) , size = 100  ) :
    #      sample = data [ indices ]
    #  @endcode        
    def extended_bootstrap ( data , size , sort = False ) :
        """ Generate "extended" bootstrap
        >>> data = ...
        >>> for ds = extended_bootstrap ( data , size = 25 ) :
        - Generate indices for bootstrap samples:
        >>> data = ...
        >>> N    = len ( data )
        >>> for indices  in extended_bootstrap ( range ( N ) , size = 100  ) :
        >>> ... sample = data [ indices ]
        """        
        N = len ( data )
        for i in range ( size ) :
            indices = choices ( data , k = poisson ( N ) )
            if sort : indices = sorted ( indices )             
            yield tuple ( indices )
            
# =============================================================================
## Generate indices for bootstrap samples:
#  @code
#  data = ...
#  for indices  in bootstrap_indices ( N , size = 100  ) :
#      sample = data [ indices ]
#  @endcode
def bootstrap_indices ( N , size = 100 , sort = False ) :
    """ Generate indices for bootstrap  samples:
    >>> data = ...
    >>> for indices  in bootstrap_indices ( N , size = 100  ) :
    >>> ... sample = data [ indices ]
    """
    for indices in bootstrap ( range ( N ) , size , sort = sort ) :
        yield indices

# =============================================================================
## Generate indices for extedend bootstrap samples:
#  @code
#  data = ...
#  for indices  in extended_bootstrap_indices ( N , size = 100  ) :
#      sample = data [ indices ]
#  @endcode
def extended_bootstrap_indices ( N , size = 100 , sort = False ) :
    """ Generate indices for bootstrap  samples:
    >>> data = ...
    >>> for indices  in extended_bootstrap_indices ( N , size = 100  ) :
    >>> ... sample = data [ indices ]
    """
    for indices in extended_bootstrap ( range ( N ) , size , sort = sort ) :
        yield indices
  
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.logger.logger import getLogger 
    logger = getLogger( 'ostap.stats.bootstrap' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
