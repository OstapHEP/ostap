#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with some simple but useful utilities
#   - suppression of stdout/stderr 
#   - dumpting of stdout/stderr into file 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
#  
# =============================================================================
""" Module with some simple but useful utilities
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    #
    'tee_py'             , ## tee for Python's printouts
    'tee_cpp'            , ## tee for C++'s    printouts
    'output'             , ## redirect stdout/stderr into the file 
    'mute_py'            , ## suppress stdout/strerr Python printout 
    'silence_py'         , ## ditto 
    'mute'               , ## context manager to suppress stdout/strerr printout 
    'silence'            , ## ditto 
    'NoContext'          , ## empty context manager
    ## logging   
    'logColor'           , ## switch on  locally the colored logging
    'logNoColor'         , ## switch off locally the colored logging
    'logVerbose'         , ## redefine (locally) the logging level
    'logDebug'           , ## redefine (locally) the logging level
    'logInfo'            , ## redefine (locally) the logging level
    'logWarning'         , ## redefine (locally) the logging level
    'logError'           , ## redefine (locally) the logging level
    ##
    'multicolumn'        , ## format the list of strings into multicolumn block
    'DisplayTree'        , ## display tree-like structures
    ##
    'print_args'         , ## print all aguments as 3-column table
    'map2table'          , ## Format map/dict-like object  as a 2-column table
    'map2table_ex'       , ## Format map/dict-like object  as a 3-column table
    # =========================================================================
) # ===========================================================================
# =============================================================================
from   ostap.core.ostap_types import integer_types, num_types, string_types  
from   ostap.logger.logger    import logVerbose,  logDebug, logInfo, logWarning, logError
from   ostap.math.base        import isfinite, iszero, frexp10
from   ostap.logger.mute      import ( mute   , mute_py ,
                                       tee_py , tee_cpp ,
                                       output , silence , silence_py ,
                                       MuteC  , MutePy  ,
                                       TeeCpp , TeePy   , OutputC    )
from   ostap.utils.basic      import NoContext, typename, prntrf, items_loop  
from   ostap.logger.symbols   import plus_minus, times  
import time, os, sys, math  ## attention here!!
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger    import getLogger, logColor, logNoColor 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.logger.utils' )
else                       : logger = getLogger( __name__ )
del getLogger 

# =============================================================================
## format list of strings into multicolumn string
#  @code 
#  >>> strings =  ....
#  >>> table   = multicolumn (  strings , indent = 2 )
#  >>> print table
#  @endcode 
def multicolumn ( lines , term_width = None , indent = 0 , pad = 2 ):
    """ Format list of strings into multicolumn string
    >>> strings =  ....
    >>> table   = multicolumn (  strings , indent = 2 )
    >>> print table 
    """
    n_lines = len(lines)
    if n_lines == 0:
        return
    
    if not term_width :
        from ostap.utils.basic import terminal_size 
        term_width , h = terminal_size()
        
    col_width = max(len(line) for line in lines)
    n_cols = int((term_width + pad - indent)/(col_width + pad))
    n_cols = min(n_lines, max(1, n_cols))
    
    col_len = int(n_lines/n_cols) + (0 if n_lines % n_cols == 0 else 1)
    if (n_cols - 1) * col_len >= n_lines:
        n_cols -= 1
        
    cols = [lines[i*col_len : i*col_len + col_len] for i in range(n_cols)]
    
    rows        = list(zip(*cols))
    rows_missed = list(zip(*[col[len(rows):] for col in cols[:-1]]))
    rows.extend(rows_missed)

    result = []
    for row in rows:
        line = " "*indent + (" "*pad).join(line.ljust(col_width) for line in row)
        result.append ( line )
    return '\n'.join ( result )

# =======================================================================================
## @class DisplayTree
#  Very simple class to allow a rendering of Tree-like objects,
#  in particular the structure and the content of directories 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2020-01-029
class DisplayTree ( object ) :
    """ Very simple class to allow a rendering of Tree-like objects,
    in particular the structure and the content of directories 
    """    
    prefix_item_middle   = '\033(0\x74\033(B' + '\033(0\x71\033(B' + '\033(0\x71\033(B'
    prefix_item_last     = '\033(0\x6d\033(B' + '\033(0\x71\033(B' + '\033(0\x71\033(B'
    prefix_parent_middle = '    '  
    prefix_parent_last   = '\033(0\x78\033(B' +  '   '


    def __init__  ( self            ,
                    payload         ,   ## the object 
                    parent  = None  ,   ## parent if any  
                    display = None  ,   ## display function: how to display it? 
                    isdir   = False ,   ## directory object ?
                    last    = False ) : ## last in the list? 
        
        self.__payload = payload
        self.__parent  = parent

        if   display :         ## use display function if specified 
            self.__display = display
        elif parent  :         ## otherwise pick if from the parent 
            self.__display = parent.__display
            
        self.__last    = last
        self.__isdir   = isdir        
        self.__depth   = self.parent.depth + 1 if parent else 0 

    @property
    def payload ( self ) :
        """`payload'  : the actual object, hold by the tree leaf"""
        return self.__payload
    @property
    def parent ( self ) :
        """`parent' : the parent element in the tree"""
        return self.__parent
    @property
    def last   (  self ) :
        """`last'  : is it a last object oin the list? """
        return self.__last
    @property
    def isdir  (  self ) :
        """``isdir''  : is it `directory-like' item? """
        return self.__isdir    
    @property
    def depth  (  self ) :
        """`depth'  : how deep is our tree? """
        return self.__depth 
    
    @property 
    def itemname ( self ) :
        """`itemname' : how the name of the leaf is displayed?"""
        return self.__display ( self.payload , self.isdir )

    # ==========================================================================
    ##  show the (sub)tree  
    def showme ( self , prefix = '' ) :
        """ Show the constructed (sub)tree  
        """        
        if self.parent is None:
            return self.itemname

        symbols = self.prefix_item_last if self.last else self.prefix_item_middle

        graph   = [ prefix + '{!s} {!s}'.format ( symbols , self.itemname ) ]

        parent = self.parent        
        while parent and parent.parent is not None:
            pp     = self.prefix_parent_middle if parent.last else self.prefix_parent_last
            graph.append  ( pp )
            parent = parent.parent

        return ''.join ( reversed ( graph ) )

    def __str__ (  self ) :
        return 'Leaf(%s,isdir=%s,last=%s)' % ( self.itemname ,
                                               self.isdir    ,
                                               self.last     ) 
    __repr__ = __str__
    
# ============================================================================
## Print all arguments as table
#  Two keyword argument are not considered as "arguments" and just ignored
#  (used of only for formatting the output table)
#  - prefix 
#  - title
def print_args ( *args , title = '' , prefix = '' , **kwargs ) :
    """ Print all arguments as table 
    Two keyword argument are not considered as "arguments" and  justr ignored 
    (used of only for formatting the output table)
    - `prefix` 
    - `title` 
    """
    rows  = [  ( 'Argument' , 'type' , 'value' ) ]

    for i , a in enumerate ( args ) : 
        row  = '#%d' % i , typename ( a ) , prntrf ( a )
        rows.append ( row )

    ## if 'prefix' in kwargs : kwargs [ '*prefix' ] = kwargs.pop ( 'prefix' )
    ## if 'title'  in kwargs : kwargs [ '*title'  ] = kwargs.pop ( 'title'  ) 
    
    for key in sorted ( kwargs ) :
        v    = kwargs [ key ] 
        row  = '%s' % key , typename ( v ) , prntrf  ( v )
        rows.append ( row  )

    ## prefix = kwargs.pop ( '*prefix' , '' )
    ## if not prefix or not isinstance ( prefix , string_types ) : prefix = ''
    ## title  = kwargs.pop ( '*title' , 'Arguments' )
    ## if not title  or not isinstance ( title  , string_types ) : title = 'Arguments'
    
    import ostap.logger.table as T
    return T.table ( rows      ,
                     title     = title  if title else 'Arguments',
                     prefix    = prefix ,
                     alignment = 'rww'  ) 

# =============================================================================
## Format map/dict-like object  as a table
#  @code
#  map = ...
#  title = 'My dictionary'
#  print ( map2table ( map , title = title ) ) 
#  @endcode 
def map2table ( dct                             ,
                header    = ( 'Key' , 'Value' ) ,
                prefix    = ''                  ,
                title     = ''                  ,
                alignment = 'rl'                ,
                style     = ''                  ) :
    """ Format map/dict-like object  as a table
    >>> map = ...
    >>> title = 'My dictionary'
    >>> print ( map2table ( map , title = title ) ) 
    """
    rows = [ header ]
    map  = {}
    map.update ( dct )
    for key in sorted ( map ) :
        value = map [ key ]
        row   = str ( key ) , str ( value )
        rows.append ( row )
        
    import ostap.logger.table as T
    return T.table ( rows                  ,
                     title     = title     ,
                     prefix    = prefix    ,
                     alignment = alignment ,
                     style     = style     ) 

# =============================================================================
## Format map/dict-like object  as a "extended" table (with type for values)
#  @code
#  map = ...
#  title = 'My dictionary'
#  print ( map2table_ex ( map , title = title ) ) 
#  @endcode 
def map2table_ex ( dct                                    ,
                   header    = ( 'Key' , 'Type' , 'Value' ) ,
                   prefix    = ''                  ,
                   title     = ''                  ,
                   alignment = 'rcc'               ,
                   style     = ''                  ) :
    """ Format map/dict-like object  as a table
    >>> map = ...
    >>> title = 'My dictionary'
    >>> print ( map2table ( map , title = title ) ) 
    """
    rows = [ header ]
    map  = {}
    map.update ( dct )
    for key in sorted ( map ) :
        value = map [ key ]
        row   = str ( key ) , typename ( value ) , str ( value )
        rows.append ( row )
        
    import ostap.logger.table as T
    return T.table ( rows                  ,
                     title     = title     ,
                     prefix    = prefix    ,
                     alignment = alignment ,
                     style     = style     ) 

# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
        
# =============================================================================
##                                                                     The END 
# =============================================================================
