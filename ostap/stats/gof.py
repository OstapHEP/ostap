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
import ROOT, abc  
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
## class AGoF(object) :
class AGoF(abc.ABC) :
    """ An abstract base class for family of methods to probe Goodness-of-Fit
    
    There are three abstract methods
    - `__call__` to evaluate t-value, the value of GoF estimator 
    - `tvalue` to evaluate t-value 
    - `pvalue` to evaluate (t,p)-values
    
    >>> gof  = ...
    >>> pdf  = 
    >>> data = ...
    >>> pdf.fitTo ( data , ... ) 
    >>> t_value            = gof        ( pdf , data )
    >>> t_value            = gof.tvalue ( pdf , data )
    >>> t_value , p_value  = gof.pvalue ( pdf , data )
    
    """
    # =========================================================================
    ## Calculate T-value for Goodness-of-Fit test
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
    ## Calculate T-value for Goodness-of-Git
    #  @code
    #  gof    = ...
    #  pdf    = ...
    #  data   = ...
    #  tvalue = gof.tvalue ( pdf , data )
    #  @endcode
    @abc.abstractmethod
    def tvalue ( self , pdf , data ) :
        """ Calculate T-value for Goodness-of-Fit
        >>> gof    = ...
        >>> pdf    = ...
        >>> data   = ...
        >>> tvalue = gof.tvalue ( pdf , data )
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
    
    # ==========================================================================
    ## Get GoF results in form of the table
    #  @code
    #  gof = ....
    #  table = gof.table  ()
    #  table = gof.report () ## ditto
    #  @endcode 
    def report ( self             ,
                 tvalue    = None ,
                 pvalue    = None ,
                 ecdf      = None ,
                 counter   = None ,
                 precision = 4    ,
                 width     = 6    , 
                 title     = ''   ,
                 prefix    = ''   ,
                 style     = ''   ) :
        """ Get results in form of the table
        >>> gof = ....
        >>> table = gof.report () ## ditto
        
        """
        from ostap.stats.gof_utils import format_table
        from ostap.utils.basic     import typename 
        return format_table  ( tvalue    = tvalue    ,
                               pvalue    = pvalue    ,
                               ecdf      = ecdf      ,
                               counter   = counter   ,
                               precision = precision ,
                               width     = width     ,
                               title     = title if title else '%s report' % typename ( self ) , 
                               prefix    = prefix    ,
                               style     = style     )
    
    # ==========================================================================
    ## Get results in form of the row in the summary table
    #  @code
    #  gof = ...
    #  header , row = gof.the_row ( ... ) 
    #  @endcode `
    def the_row ( self             ,
                  tvalue    = None ,
                  pvalue    = None ,
                  ecdf      = None ,
                  counter   = None ,                  
                  precision = 4    ,
                  width     = 6    ) : 
        
        """ Get results in form of the table 
        >>> gof = ...
        >>> header , row = gof.the_row ( ... ) 
        """
        from ostap.stats.gof_utils import format_row 
        return format_row ( tvalue    = tvalue    ,
                            pvalue    = pvalue    ,
                            ecdf      = ecdf      ,  
                            counter   = counter   ,
                            precision = precision ,
                            width     = width     )

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
#  t    = gof        ( ds1 , ds2 , normalize = True )
#  t    = gof.tvalue ( ds1 , ds2 , normalize = True )
#  t,p  = gof.pvalue ( ds1 , ds2 , normalize = True )
#  @endcode
## class AGoF(abc.ABC) :
class AGoFnp(abc.ABC) :
    """ An abstract base class for numpy-related family of methods to probe goodness-of fit
    
    There are two abstract methods
    - `__call__` to evaluate t-value for two datasets 
    - `pvalue` to evaluate (t,p)-vaues for two datasets
    
    >>> gof = ...
    >>> ds1, ds2 = ...
    >>> t    = gof ( ds1 , ds2 , normalize = True )
    >>> t,p  = gof ( ds1 , ds2 , normalize = True )
    """
    # =========================================================================
    ## Calculate t-value for two datasets 
    #  @code
    #  data1 = ... ## the first  data set 
    #  data2 = ... ## the second data set
    #  gof   = ...
    #  t = gof ( data1 , data2 , normalize = False ) 
    #  t = gof ( data1 , data2 , normalize = True  ) 
    #  @endcode
    @abc.abstractmethod
    def __call__ ( self      ,
                   data1     ,
                   data2     , * ,
                   weight1   = None ,
                   weight2   = None ,                   
                   normalize = True ) :
        """ Calculate T-value for two data sets 
        >>> gof    = ...
        >>> data1 = ... ## the first  data set 
        >>> data2 = ... ## the second data set
        >>> t = gof ( data1 , data2 , normalize = False ) 
        >>> t = gof ( data1 , data2 , normalize = True  ) 
        """
        return NotImplemented 
            
    # =========================================================================
    ## Calculate t-value 
    #  @code
    #  gof = ...
    #  ds1 , ds2 = ...
    #  t   = gof.tvalue ( ds1 , ds2 , normalize = True ) 
    #  @endcode 
    @abc.abstractmethod
    def tvalue ( self      ,
                 data1     ,
                 data2     , *    ,
                 weight1   = None ,
                 weight2   = None ,
                 normalize = True ) : 
        """ Calculate the t-value
        >>> gof = ...
        >>> ds1 , ds2 = ...
        >>> t   = gof.tvalue ( ds1 , ds2 , normalize = True ) 
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
    def pvalue ( self      ,
                 data1     ,
                 data2     , *    ,
                 weight1   = None ,
                 weight2   = None ) : 
        """ Calculate the t & p-values
        >>> gof = ...
        >>> ds1 , ds2 = ...
        >>> t , p = gof.pvalue ( ds1 , ds2 , normalize = True ) 
        """
        return NotImplemented

    # =========================================================================
    ## Are weights supported by this GoF estimator?
    @property 
    @abc.abstractmethod
    def weights_supported ( self ) :
        """`weghts_supported`: Are weights supported by this estimator?
        """
        return NotImplemented

    # =========================================================================
    ## Good for two-samples comparison?
    #  Can this estimator be used for comparison of two samples?
    @property 
    @abc.abstractmethod
    def two_samples ( self ) :
        """`two_samples`: Can this estimator be used for comparison of two samples?
        """
        return NotImplemented
    
    # ==========================================================================
    ## Get results in form of the table 
    def report ( self             ,
                 tvalue    = None ,
                 pvalue    = None ,
                 ecdf      = None ,
                 counter   = None ,
                 precision = 4    ,
                 width     = 6    , 
                 title     = ''   ,
                 prefix    = ''   ,
                 style     = ''   ) :
        """ Get results in form of the table 
        """
        from ostap.stats.gof_utils import format_table
        from ostap.utils.basic     import typename 
        return format_table  ( tvalue    = tvalue    ,
                               pvalue    = pvalue    ,
                               ecdf      = ecdf      ,
                               counter   = counter   ,
                               precision = precision ,
                               width     = width     ,
                               title     = title if title else '%s report' % typename ( self ) , 
                               prefix    = prefix    ,
                               style     = style     )

    # ==========================================================================
    ## Get results in form of the row in the table
    #  @code
    #  gof = ...
    #  header , row = gof.the_row ( ... ) 
    #  @endcode `
    def the_row ( self             ,
                  tvalue    = None ,
                  pvalue    = None ,
                  ecdf      = None ,
                  counter   = None ,                  
                  precision = 4    ,
                  width     = 6    ) : 
        
        """ Get results in form of the table 
        >>> gof = ...
        >>> header , row = gof.the_row ( ... ) 
        """
        from ostap.stats.gof_utils import format_row 
        return format_row ( tvalue    = tvalue    ,
                            pvalue    = pvalue    ,
                            ecdf      = ecdf      ,  
                            counter   = counter   ,
                            precision = precision ,
                            width     = width     )
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================


    
