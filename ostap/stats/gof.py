#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/gof.py
#  Set of utulities for goodness-of-fit studies 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2023-12-06
# =============================================================================
""" Simple utulities for goodness-of-fit studies 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2023-12-06"
__all__     = (
    'AGoF'      , ## an abstract base cladd for Goodness-of-Fit estimators
    'AGoFnp'    , ## an abstract base cladd for Goodness-of-Fit estimators
    )
# =============================================================================
from   ostap.core.core        import VE, Ostap
from   ostap.core.ostap_types import string_types
import abc  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.gof' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Simple utilities for goodness-of-fit studies' )
# =============================================================================
## @class AGoF
#  An absract base class for family of methods to probe Goodness-of-Fit.
# - `RooFit` interface 
#  There are two abstract methods
#  - <code>__call__</code> to evaluate t-value, the value of GoF estimator 
#  - <code>pvalue</code> to evaluate (t,p)-values for two datasets
#  @code
#  gof  = ...
#  pdf  = ...
#  data = ...
#  pdf.fitTo ( data ) 
#  t_value           = gof        ( pdf , data )
#  t_value, p_value  = gof.pvalue ( pdf , data )
#  @endcode
class AGoF(object) :
    """ An abstract base class for family of methods to probe Goodness-of-Git
    There are two abstract methods
    - `__call__` to evaluate t-value, the value of GoF estimator 
    - `pvalue` to evaluate (t,p)-vaues
    
    >>> gof  = ...
    >>> pdf  = 
    >>> data = ...
    >>> pdf.fitTo ( data , ... ) 
    >>> t_value            = gof        ( pdf , data )
    >>> t_value , p_value  = gof.pvalue ( pdf , data )
    """
    # =========================================================================
    ## Calculate T-value for Goodness-of-Git 
    #  @code
    #  gof   = ...
    #  pdf   = ...  
    #  data  = ... 
    #  t_value = gof ( pdf , data ) 
    #  @endcode
    @abc.abstractmethod
    def __call__ ( self , pdf , data ) :
        """ Calculate T-value for Goodness-of-Fit
        >>> gof   = ...
        >>> pdf   = ... 
        >>> data  = ... 
        >>> t_value = gof ( pdf , data ) 
        """
        return NotImplemented 
    # =========================================================================
    ## Calculate the t & p-values
    #  @code
    #  gof  = ...
    #  pdf  = ...
    #  data = ... 
    #  t_value , p_value = gof.pvalue ( pdf , data )
    #  @endcode 
    @abc.abstractmethod
    def pvalue ( self , pdf , data ) :
        """ Calculate the t & p-values
        >>> gof  = ...
        >>> pdf  = ... 
        >>> data = ... 
        >>> t_value , p_value = gof.pvalue ( pdf , data ) 
        """
        return NotImplemented

# =============================================================================
## @class AGoFnp
#  An absract base class for numpy-related family of methods to probe goodness-of fit
#
#  There are two abstract methods
#  - <code>__call__</code> to evaluate t-value for two datasets 
#  - <code>pvalue</code> to evaluate (t,p)-vaues for two datasets
#  @code
#  gof = ...
#  ds1, ds2 = ...
#  t    = god        ( ds1 , ds2 , normalize = True )
#  t,p  = god.pvalue ( ds1 , ds2 , normalize = True )
#  @endcode
class AGoFnp(object) :
    """ An absract base class for numpy-related family of methods to probe goodness-of fit
    
    There are two abstract methods
    - `__call__` to evaluate t-value for two datasets 
    - `pvalue` to evaluate (t,p)-vaues for two datasets
    
    >>> gof = ...
    >>> ds1, ds2 = ...
    >>> t    = gof ( ds1 , ds2 , normalize = True )
    >>> t,p  = gof ( ds1 , ds2 , normalize = True )
    """
    # =========================================================================
    ## Calculate T-value for two datasets 
    #  @code
    #  data1 = ... ## the first  data set 
    #  data2 = ... ## the second data set
    #  gof   = ...
    #  t = gof ( data1 , data1 , normalize = False ) 
    #  t = gof ( data1 , data1 , normalize = True  ) 
    #  @endcode
    @abc.abstractmethod
    def __call__ ( self , data1 , data2 , normalize = True ) :
        """ Calculate T-value for two data sets 
        >>> gof    = ...
        >>> data1 = ... ## the first  data set 
        >>> data2 = ... ## the second data set
        >>> t = gof ( data1 , data1 , normalize = False ) 
        >>> t = gof ( data1 , data1 , normalize = True  ) 
        """
        return NotImplemented 
    # =========================================================================
    ## Calculate the t & p-values
    #  @code
    #  gof = ...
    #  ds1 , ds2 = ...
    #  t , p = gof.pvalue ( ds1 , ds2 , normalize = True ) 
    #  @endcode 
    @abc.abstractmethod
    def pvalue ( self , data1 , data2 , normalize = True ) :
        """ Calculate the t & p-values
        >>> gof = ...
        >>> ds1 , ds2 = ...
        >>> t , p = gof.pvalue ( ds1 , ds2 , normalize = True ) 
        """
        return NotImplemented


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================


    
