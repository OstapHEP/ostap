#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/roocollections.py
#  Module with decoration for  RooArgSet, RooArgList and similar objects 
#  @see RooArgSet 
#  @see RooArgList 
#  @see RooLinkedList 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Module with decoration for  RooArgSet, RooArgList and similar objects
- see RooArgSet 
- see RooArgList 
- see RooLinkedList 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'KeepArgs'  , ## context manager to preserve the content of the list
    ) 
# =============================================================================
import ROOT, sys, random
from   ostap.core.core         import Ostap
from   ostap.core.ostap_types  import string_types 
import ostap.fitting.variables
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.roocollections' )
else                       : logger = getLogger( __name__ )
# =============================================================================
if (3,3) < sys.version_info : from collections.abc import Container
else                        : from collections     import Container
# =============================================================================
_new_methods_ = []
# =============================================================================


# =============================================================================
## iterator for RooArgList
#  @code
#  lst = ...
#  for a in lst : print a 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2011-06-07
def _ral_iter_ ( self ) :
    """Iterator for RooArgList:
    >>> arg_list = ...
    >>> for p in arg_list : print p    
    """
    l = len ( self )
    for i in range ( 0 , l )  :
        yield self[i]

# =============================================================================
## check presence of element or index in the list
def _ral_contains_ ( self , i ) :
    """Check the presence of element or index in the list
    """
    if   isinstance ( i , int            ) : return 0<= i < len(self)
    elif isinstance ( i , string_types   ) : return self.find     ( i ) 
    elif isinstance ( i , ROOT.RooAbsArg ) : return self.contains ( i ) 
    obj = self.find ( i )
    return True if obj else False  
    ## return  0 <= self.index ( i )

# =============================================================================
##  get item form the list 
def _ral_getitem_ ( self , index ) :
    """Get item from the list
    >>>  lst  = ...
    >>>  item = lst[3]
    Slice notation is also supported (note that  retubned type is python tuple)
    >>>  res = lst[3:5:2]
    """
    if isinstance ( index , slice ) :
        l = len ( self )
        indices  = index.indices ( l )
        return tuple ( [ self[i] for i in range( *indices ) ] )
    
    if not isinstance ( index , int ) or not index in self :
        raise IndexError('List Index %s is out of the range [%d,%d)' % ( index , 0 , len(self) ) )
    return self.at ( index  )

# =============================================================================
## some decoration over RooArgList 
ROOT.RooArgList . __len__       = lambda s   : s.getSize()
ROOT.RooArgList . __contains__  = _ral_contains_ 
ROOT.RooArgList . __iter__      = _ral_iter_
ROOT.RooArgList . __nonzero__   = lambda s   : 0 != len ( s ) 
ROOT.RooArgList . __getitem__   = _ral_getitem_
ROOT.RooArgList . __setitem__   = lambda s,*_ : NotImplemented 

_new_methods_ += [
    ROOT.RooArgList. __len__       ,
    ROOT.RooArgList. __contains__  ,
    ROOT.RooArgList. __iter__      ,
    ROOT.RooArgList. __nonzero__   ,
    ROOT.RooArgList. __getitem__   ,
    ROOT.RooArgList. __setitem__   ,
]    

# =============================================================================
## helper function to print collection
def _rs_list_ ( self ) :
    _l = []
    for i in self :
        
        if   hasattr  ( i , 'GetName' ) and hasattr ( i , 'getVal' ) :
            _l.append ( i.GetName() + ":%s" % i.getVal() )
        elif hasattr  ( i , 'GetName' ) :
            _l.append ( i.GetName()   )
        elif hasattr  ( i , 'getVal'  ) :
            _l.append ( "%s" % i.getVal ()  )
        else :
            _l.append (  str ( i )    )
            
    return _l ;

# =============================================================================
## Get tuple of names for the objects in <code>RooAbsCollection</code>
#  @code
#  collection = ...
#  names = collection.names ()  
#  @endcode
#  @see RooArgList
#  @see RooArgSet
#  @see RooAbsCollection
def _rac_names_ ( self ) :
    """ Get tuple of names for the objects in `ROOT.RooAbsCollection`
    >>> collection = ...
    >>> names = collection.names ()  
    - see `ROOT.RooArgList`
    - see `ROOT.RooArgSet`
    - see `ROOT.RooAbsCollection`
    """
    return tuple ( [ i.name for i in self ] ) 

ROOT.RooAbsCollection . names = _rac_names_
ROOT.RooArgList.names         = _rac_names_
ROOT.RooArgSet .names         = _rac_names_

# =============================================================================
## printout for RooArgList 
ROOT.RooArgList . __str__       = lambda s : str ( _rs_list_ ( s ) )  
ROOT.RooArgList . __repr__      = lambda s : str ( _rs_list_ ( s ) )  


# =============================================================================
## iterator for RooArgSet
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _ras_iter_ ( self ) :
    """Simple iterator for RootArgSet:
    >>> arg_set = ...
    >>> for i in arg_set : print i    
    """    
    it  = Ostap.Utils.Iterator ( self )
    val = it.Next()
    while val :
        yield val 
        val = it.Next()
        
    del it

# =============================================================================
## get the attibute for RooArgSet 
def _ras_getattr_ ( self , aname ) :
    """Get the attibute from RooArgSet
    >>> aset = ...
    >>> print aset.pt    
    """
    _v = self.find ( aname )
    if not _v : raise  AttributeError
    return _v 

# =============================================================================
## get the attibute for RooArgSet 
def _ras_getitem_ ( self , aname ) :
    """Get the attibute from RooArgSet
    >>> aset = ...
    >>> print aset['pt']    
    """
    _v = self.find ( aname )
    if not _v : raise  IndexError
    return _v 

# =============================================================================
## check the presence of variable in set 
def _ras_contains_ ( self , aname ) :
    """Check the presence of variable in set 
    """
    _v = self.find ( aname )
    if not _v : return False 
    return             True 

# =============================================================================
## some decoration over RooArgSet 
ROOT.RooArgSet . __len__           = lambda s   : s.getSize()
ROOT.RooArgSet . __iter__          = _ras_iter_ 
ROOT.RooArgSet . __getattr__       = _ras_getattr_ 
ROOT.RooArgSet . __getitem__       = _ras_getitem_ 
ROOT.RooArgSet . __contains__      = _ras_contains_ 
ROOT.RooArgSet . __nonzero__       = lambda s   : 0 != len ( s ) 
        
ROOT.RooArgSet     . __str__       = lambda s : str ( tuple ( _rs_list_ ( s ) ) )  
ROOT.RooArgSet     . __repr__      = lambda s : str ( tuple ( _rs_list_ ( s ) ) )  
ROOT.RooLinkedList . __repr__      = lambda s : str (  _rs_list_ ( s ) )

ROOT.RooAbsCollection.__iter__     = _ras_iter_
ROOT.RooAbsCollection.__len__      = lambda s   : s.getSize()
ROOT.RooAbsCollection. __nonzero__ = lambda s   : 0 != len ( s ) 


# =============================================================================
## iterator for class RooLinkedList
#  @code
#  lst = ...
#  for l in lst : print l 
#  @endcode 
def _rll_iter_  ( self ) :
    """Iterator over RooLinekdList
    >>> lst = ...
    >>> for l in lst : print l     
    """
    l = len ( self )  
    for i in range ( l ) : yield self.At ( i )
        
ROOT.RooLinkedList . __iter__  = _rll_iter_ 

_new_methods_ += [
    ROOT.RooArgSet . __len__      ,
    ROOT.RooArgSet . __iter__     ,
    ROOT.RooArgSet . __getattr__  , 
    ROOT.RooArgSet . __setattr__  ,
    ROOT.RooArgSet . __contains__ ,
    ROOT.RooArgSet . __nonzero__  ,
    ROOT.RooArgSet . __str__      ,
    ROOT.RooArgSet . __repr__     ,
    ]

ROOT.RooLinkedList.add     = ROOT.RooLinkedList.Add
ROOT.RooLinkedList.__len__ = ROOT.RooLinkedList.GetSize

_new_methods_ += [
    ROOT.RooLinkedList . __len__  ,
    ROOT.RooLinkedList . __repr__ ,
    ROOT.RooLinkedList . __iter__ ,
    ROOT.RooLinkedList . add 
    ]
    
# =============================================================================
## add more data into list/set
def _ral_iadd_ ( self , other ) :
    """Update/increment collections
    >>> lst = ....
    >>> lst += another_lst
    """
    _CNT = Container
    _RAC = ROOT.RooAbsCollection
    _RAA = ROOT.RooAbsArg      
    if not isinstance ( other , ( _CNT, _RAC , _RAA ) ) : return NotImplemented
    if     isinstance ( other , str )                   : return NotImplemented

    ##
    if isinstance ( other , _RAA ) and not isinstance ( other , _RAC ) : other = [ other ]

    for o in other : self.add ( o )
    return self

# =============================================================================
## add more data into list/set
def _ral_add_ ( self , other ) :
    """Make a sum of two lists/sets/collections
    >>> lst1 = ...
    >>> set2 = ...
    >>> lst2 = lst1 + set2 
    """
    _CNT = Container
    _RAC = ROOT.RooAbsCollection
    _RAA = ROOT.RooAbsArg  
    if not isinstance ( other , ( _CNT, _RAC , _RAA ) ) : return NotImplemented
    if     isinstance ( other , str )                   : return NotImplemented
    _clone  = self.clone('')
    _clone += other
    return _clone

# =============================================================================
## add two list/sets 
def _ral_radd_ ( self , other ) : 
    """Make a sum of two lists/sets/collections
    >>> lst1 = ...
    >>> set2 = ...
    >>> lst2 = lst1 + set2 
    """
    return self + other

# ============================================================================
def _ral_clone_  ( self , name = '' ) :
    return self.Clone( name )

# ============================================================================
def _ral_clear_ ( self ) :
    """Clear the container
    """
    return self.removeAll() 
# =============================================================================
for t in ( ROOT.RooArgList , ROOT.RooArgSet , ROOT.RooLinkedList ) :
    t. clone    =  _ral_clone_ 
    t. clear    =  _ral_clear_ 
    t. __add__  =  _ral_add_
    t.__iadd__  = _ral_iadd_
    t.__radd__  = _ral_radd_
    t.append    = _ral_iadd_

# ==============================================================================
## get an intersection of  two sets e
#  @code
#  set1 = ...
#  set2 = ...
#  set3 = set1.intersection ( set2 ) 
#  @endcode 
def _ras_intersection_ ( self , another ) :
    """Get an intersection of two sets
    >>> set1 = ...
    >>> set2 = ...
    >>> set3 = set1.intersection ( set2 ) 
    """
    Type   = type ( self ) 
    result = Type ()
    for arg in self :
        if arg in another : result.add ( arg )
    return result

# ==============================================================================
## get an union of two sets 
#  @code
#  set1 = ...
#  set2 = ...
#  set3 = set1.union  ( set2 ) 
#  @endcode 
def _ras_union_ ( self , another ) :
    """Get an union of two sets
    >>> set1 = ...
    >>> set2 = ...
    >>> set3 = set1.union ( set2 ) 
    """
    Type   = type ( self ) 
    result = Type ()
    for arg in self    : result.add ( arg )
    for arg in another :
        if not arg in self : result.add ( arg ) 
    return result

# ================================================================================
## get a difference fot two sets
#  @code
#  set1 = ...
#  set2 = ...
#  set3 = set1.difference( set2 )
#  set3 = set1 - set2 
#  @endcode 
def _ras_difference_ (  self , another ) :
    """Get a difference fot two sets
    >>> set1 = ...
    >>> set2 = ...
    >>> set3 = set1.difference( set2 )
    >>> set3 = se1 - set2
    """
    Type   = type ( self ) 
    result = Type ()
    for arg in self    :
        if not arg in another : result.add ( arg )
    return result

# ================================================================================
## get a symmetric difference fot two sets
#  @code
#  set1 = ...
#  set2 = ...
#  set3 = set1.symmetric_difference( set2 ) 
#  @endcode 
def _ras_symmetric_difference_ (  self , another ) :
    """Get a difference fot two sets
    >>> set1 = ...
    >>> set2 = ...
    >>> set3 = set1.symmetric_difference( set2 )
    """
    Type   = type ( self ) 
    result = Type ()
    for arg in self    :
        if not arg in another : result.add ( arg )
    for arg in another :
        if not arg in self    : result.add ( arg )        
    return result

ROOT.RooArgSet . intersection         = _ras_intersection_
ROOT.RooArgSet . union                = _ras_union_
ROOT.RooArgSet . difference           = _ras_difference_
ROOT.RooArgSet . symmetric_difference = _ras_symmetric_difference_
ROOT.RooArgSet . __sub__              = _ras_difference_ 
ROOT.RooArgSet . __or__               = ROOT.RooArgSet.__add__
ROOT.RooArgSet . __ior__              = ROOT.RooArgSet.__iadd__

_new_methods_ += [
    ROOT.RooArgList    . clone        ,
    ROOT.RooArgSet     . clone        ,
    ROOT.RooLinkedList . clone        ,
    ROOT.RooArgList    . __add__      ,
    ROOT.RooArgSet     . __add__      ,
    ROOT.RooLinkedList . __add__      ,
    ROOT.RooArgList    . __iadd__     ,
    ROOT.RooArgSet     . __iadd__     ,
    ROOT.RooLinkedList . __iadd__     ,
    ROOT.RooArgList    . __radd__     ,
    ROOT.RooArgSet     . __radd__     ,
    ROOT.RooLinkedList . __radd__     ,
    ROOT.RooArgList    . append       ,
    ROOT.RooArgSet     . append       ,
    ROOT.RooLinkedList . append       ,
    ##
    ROOT.RooArgSet     . intersection         ,
    ROOT.RooArgSet     . union                ,
    ROOT.RooArgSet     . difference           ,
    ROOT.RooArgSet     . symmetric_difference ,
    ROOT.RooArgSet     . __sub__  , 
    ROOT.RooArgSet     . __or__   , 
    ROOT.RooArgSet     . __ior__  ,
    ##
    ROOT.RooAbsCollection . __iter__  ,
    ROOT.RooAbsCollection . __len__   ,
    ROOT.RooAbsCollection . names     ,
    ##
    ROOT.RooArgList       . names     ,
    ROOT.RooArgSet        . names     ,    
    ##
    ]

# =============================================================================
## @class KeepArg
#  Simple contect manager for temporary redefiniiton of some mutable collection
#  @code
#  signals = ...
#  new_signals = ...
#  print signals 
#  with KeepArgs ( signals , new_signals ) :
#      print signals
#  print signals 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
class  KeepArgs(object) :
    """Simple contect manager for temporary redefitnnnonn of some mutbale collection
    
    >>> signals = ...
    >>> new_signals = ...
    >>> print signals 
    >>> with KeepArgs ( signals , new_signals ) :
    ...     print signals
    >>> print signals 
    """
    def __init__ ( self , old_list , new_list ) :

        self.__old_list = old_list
        self.__new_list = new_list
        self.__content = []
        
    def __enter__ ( self ) :

        ## preserve the content of the old list 
        self.__content = [ i for i in self.__old_list ]
        ## clear the old list 
        self.__old_list.clear ()
        ## fill it with the content of new list 
        for i in self.__new_list : self.__old_list.add ( i )
        
        return self

    def __exit__ ( self , *_ ) :

        ## clear old list 
        self.__old_list.clear()
        ##  restore it content 
        for i in self.__content : self.__old_list.add ( i )
        ##
        self.__old_list = None 
    
# =============================================================================
_decorated_classes_ = (
    ROOT.RooArgSet        , 
    ROOT.RooArgList       , 
    ROOT.RooLinkedList    , 
    ROOT.RooAbsCollection , 
    )

_new_methods_ = tuple ( _new_methods_ ) 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
