#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
# =============================================================================
# @file shelve_ext
# 
# This is small extension for standars shelve module to allow more
# flexible treatment for the data base 
#
# @code
#
# >>> db1 = shelve.open ('$HOME/newdb'  , 'n')    ## create new DB
# >>> db2 = shelve.open ('../dir/newdb' , 'n')    ## create new DB
# ...
# >>> db1.ls()  ## list the content 
#
# @endcode 
#
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date   2012-05-13
# 
#                    $Revision$
#  Last Modification $Date$
#                 by $Author$
# =============================================================================
"""Small extension for standars shelve module to allow more flexible treatment for the data base 

>>> db1 = shelve.open ('$HOME/newdb'  , 'n')    ## create new DB
>>> db2 = shelve.open ('../dir/newdb' , 'n')    ## create new DB
...
>>> db1.ls()  ## list the content 

"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2012-05-13"
__version__ = "$Revision$" 
__all__     = ()   ## nothing to import 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.io.shelve_ext' )
else                       : logger = getLogger( __name__ )
# =============================================================================
# 
import shelve
#
_old_shelve_open_ = shelve.open
#
# ========================================================
## A bit extended version of shelve.open
def _new_shelve_open_ ( filename , *kargs , **kwargs ) :
    """A bit extended version of shelve.open:
    
    >>> db1 = open ('$HOME/a.db')
    >>> db2 = open ('../anotherdir/b.db')
    
    """
    import os
    filename = os.path.expandvars ( filename )
    filename = os.path.expanduser ( filename )
    filename = os.path.expandvars ( filename )
    filename = os.path.expandvars ( filename )
    filename = os.path.abspath    ( filename )
    #
    return _old_shelve_open_ ( filename , *kargs , **kwargs ) 

_new_shelve_open_ .__doc__ += '\n' + _old_shelve_open_ .__doc__ 

# =============================================================================
## List DB-keys 
def _ls_ ( self )  :
    """List DB-keys :
    
    >>> db = ...
    >>> db.ls() 
    
    """
    keys = self.keys()
    keys.sort()
    for k in keys : print k

# =============================================================================
## replace the method in module, if not done yet
if not hasattr  ( shelve , '_new_open_' ) :
    shelve._new_open_ = _new_shelve_open_
    shelve.     open  = _new_shelve_open_
    logger.debug ( 'Decorate shelve.open method') 

# =============================================================================
## add method to Shelve, if not done yet
if not hasattr ( shelve.Shelf, 'ls' ) : 
    shelve.Shelf.ls = _ls_
    logger.debug ( "Add 'ls' method for shelve.Shelf class") 

# =============================================================================
## add context-manager functionality to shelve 
def _shelf_enter_ ( self ) :
    """ Context-manager
    
    ## with shelve.open('A.db') as db :
    ##    ....
    
    """
    return self
def _shelf_exit_ ( self , *_ ) :
    """ Context-manager
    
    ## with shelve.open('A.db') as db :
    ##    ....
    
    """
    try :
        self.close()
    except: pass
    
    return self

# =============================================================================
## add method to Shelve, if not done yet
if hasattr ( shelve.Shelf, '__enter__' ) and hasattr ( shelve.Shelf, '__exit__' ) : pass 
else :    
    shelve.Shelf.__enter__ = _shelf_enter_
    shelve.Shelf.__exit__  = _shelf_exit_
    logger.debug ( "Add 'enter/exit' methods for shelve.Shelf class") 

# =============================================================================
if '__main__' == __name__ :

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
    
# =============================================================================
# The END 
# =============================================================================
