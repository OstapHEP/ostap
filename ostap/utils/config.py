#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Base class for configration holders 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2021-09-22
# =============================================================================
""" Base class for configration holders 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-06-07"
__all__     = (
    'ConfigBase' , ## the abstract base clss for advanced reweighters 
) 
# =============================================================================
from   ostap.utils.basic     import typename 
import numpy, os, abc 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils.config' )
else                       : logger = getLogger( __name__ )
# =============================================================================
## @class Config
#  The base class for configuration holders
#  It provides:
#   - property <code>params</code>
#   - property <code>silent</code>
#   - method <code>table</code>
#   - prints ...
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2021-09-22 
class Config ( object ) :
    """ The base class for configuration holders
    It provides:
    - property `params`
    - property `silent`    
    - method `table`
    - prints ...
    """
    def __init__ ( self , silent = False , **params  ) :
        """ The base class for configuration holders
        """
        self.__params = params 
        self.__silent = True if silent else False 
        
        if not silent : 
            title = '%s configuration' % typename ( self ) 
            table = self.table ( title = title  , prefix = '# ' )
            logger.info ( '%s:\n%s' %  ( title , table ) )
            
    @property
    def silent ( self ) :
        """`silent` : silent processing?"""
        return self.__silent
    @silent.setter
    def silent ( self , value ) :
        self.__silent = True if value else False 
    
    @property
    def params ( self ) :
        """`params`: actual configuration used for creation"""
        return self.__params

    @property 
    def config ( self ) :
        """`config` : overridable configuraton"""
        conf = {}
        conf.update ( self.params )
        conf [ 'type'   ] = typename ( self ) 
        conf [ 'silent' ] = self.silent 
        return conf 

    # =========================================================================
    ## self-print get the configuration 
    def table (  self , title = '' , prefix = '# ') : 
        """ print configuration """
        from ostap.logger.utils import map2table_ex
        title = title if title else "%s configuration " % typename ( self ) 
        return map2table_ex ( self.config , 
                              header      = ( 'Parameter' , 'type' , 'value' ) ,
                              alignment   = 'rcw'  , 
                              prefix      = prefix ,
                              title       = title  )
        
    def __str__  ( self ) : return self.table ( prefix = '' )
    def __repr__ ( self ) : return self.__str__ () 
    
# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
        
# =============================================================================
##                                                                      The END 
# =============================================================================
