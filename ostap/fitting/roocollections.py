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
""" Module with decoration for  RooArgSet, RooArgList and similar objects
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
from   collections.abc         import Container
from   ostap.core.meta_info    import root_info 
from   ostap.core.core         import std, Ostap, valid_pointer
from   ostap.utils.basic       import typename 
from   ostap.core.ostap_types  import string_types , integer_types
import ostap.fitting.variables
import ROOT, sys, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.roocollections' )
else                       : logger = getLogger( __name__ )
# =============================================================================
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
    """ Iterator for RooArgList:
    >>> arg_list = ...
    >>> for p in arg_list : print p    
    """
    l = len ( self )
    for i in range ( 0 , l )  :
        yield self[i]

# =============================================================================
## check presence of element or index in the list
def _ral_contains_ ( self , i ) :
    """ Check the presence of element or index in the list
    """
    if   isinstance ( i , int            ) : return 0<= i < len(self)
    elif isinstance ( i , string_types   ) : return self.find     ( i ) 
    elif isinstance ( i , ROOT.RooAbsArg ) : return self.contains ( i ) 
    obj = self.find ( i )
    return True if obj else False  
    ## return  0 <= self.index ( i )

# =============================================================================
##  get item from the list 
def _ral_getitem_ ( self , index ) :
    """ Get item from the list
    >>>  lst  = ...
    >>>  item = lst[3]
    >>>  iast = lst[-1]
    >>>  res  = lst[3:5:2]
    """
    if isinstance ( index , slice ) :
        l = len ( self )
        indices  = index.indices ( l )
        
        result   = ROOT.RooArgList ()
        for i in range ( *indices  ) : result.add ( self [i] )
        return result 
        ## return tuple ( [ self[i] for i in range( *indices ) ] )

    l = len ( self )
    
    ## allow slightly negative indices 
    if index < 0 : index += l 

    if not isinstance ( index , integer_types ) or l <= index :
        raise IndexError('List Index %s is out of the range [%d,%d)' % ( index , 0 , l ) )
        
    if not 0 <= index < l :
        raise IndexError('List Index %d is out of the range [%d,%d)' % ( index , 0 , l ) )
    
    return self.at ( index  )

# =============================================================================
## some decoration over RooArgList 
ROOT.RooArgList . __len__       = lambda s   : s.getSize() if valid_pointer ( s ) else 0 
ROOT.RooArgList . __contains__  = _ral_contains_ 
ROOT.RooArgList . __nonzero__   = lambda s   : valid_pointer ( s ) and 0 < s.getSize()  
ROOT.RooArgList . __bool__      = lambda s   : valid_pointer ( s ) and 0 < s.getSize()
ROOT.RooArgList . __getitem__   = _ral_getitem_
ROOT.RooArgList . __setitem__   = lambda s,*_ : NotImplemented 

_new_methods_ += [
    ROOT.RooArgList. __len__       ,
    ROOT.RooArgList. __contains__  ,
    ROOT.RooArgList. __nonzero__   ,
    ROOT.RooArgList. __bool__      ,
    ROOT.RooArgList. __getitem__   ,
    ROOT.RooArgList. __setitem__   ,
]

if not hasattr ( ROOT.RooAbsCollection , 'assign' ) :
    ROOT.RooAbsCollection.assign = ROOT.RooAbsCollection.assignValueOnly 
    _new_methods_ += [
        ROOT.RooAbsCollection.assign ,
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
ROOT.RooArgList . __str__       = lambda s : str ( _rs_list_ ( s ) ) if s else '[]' 
ROOT.RooArgList . __repr__      = lambda s : str ( _rs_list_ ( s ) ) if s else '[]'  


# =============================================================================
## iterator for RooArgSet
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
if root_info < ( 6 , 31 ) : 
    def _ras_iter_ ( self ) :
        """ Simple iterator for RooArgSet:
        >>> arg_set = ...
        >>> for i in arg_set : print i    
        """    
        it  = Ostap.Utils.Iterator ( self ) ## only for ROOT < 6.31 
        val = it.Next()
        while val :
            yield val 
            val = it.Next()    
        del it        
else :
    def _ras_iter_ ( self ) :
        """ Simple iterator for RooArgSet:
        >>> arg_set = ...
        >>> for i in arg_set : print i    
        """    
        cnt = self.get()
        for v in cnt : yield v 
            
# =============================================================================
## get the attibute for RooArgSet 
def _ras_getattr_ ( self , aname ) :
    """Get the attibute from RooArgSet
    >>> aset = ...
    >>> print aset.pt    
    """
    _v = self.find ( aname )
    if not _v :
        raise  AttributeError("%s: invalid attribute `%s'" % ( type ( self ) , aname ) )
    return _v 

# =============================================================================
## get the item for RooArgSet
def _ras_getitem_ ( self , aname ) :
    """ Get the attibute from RooArgSet
    >>> aset = ...
    >>> print ( aset[' pt']  )    
    >>> print ( aset[  0  ]  ) ## the first element
    >>> print ( aset[ -1 ]   ) ## the last element
    >>> print ( aset[ 3,5,2] ) ## get the range 
    """
    
    if   isinstance ( aname , slice ) :
        l = len ( self )
        indices = aname.indices ( l )
        result  = ROOT.RooArgSet()
        for i in range (*indices) : result.add ( self[i] ) 
        return result
    elif isinstance ( aname , integer_types ) :
        l = len ( self ) 
        ## allow "slightly negative values 
        if aname < 0 : aname += l 
        if 0 <= aname < l : 
            return ROOT.RooAbsCollection.__getitem__ ( self , aname )
        raise IndexError('Invalid index!')
    
    ## normal lookup 
    _v = self.find ( aname )
    if not _v : raise  IndexError('Invalid nidex/name!')
    return _v

# =============================================================================
## get the first element in collection 
def _rac_front_ ( self ) :
    """ Get the first element in collection
    >>> lst   = ...
    >>> first = lst.front ()
    """
    assert 0 < len ( self ) , 'front: collection is empty!'
    return self [ 0 ]
# =============================================================================
## get the last element in collection 
def _rac_back_ ( self ) :
    """ Get the last element in collection
    >>> lst = ...
    >>> last = lst.back() 
    """
    assert 0 < len ( self ) , 'back: collection is empty!'
    return self [ -1 ]
# ==============================================================================
## Remove and return item at index (default last).
#  @code
#  lst  = ...
#  item = lst.pop( )   ## get and remove from collection the last   element 
#  item = lst.pop(2)   ## get and remove from collection the second element 
#  item = lst.pop('a') ## get and remove from collection the element with name 'a'
#  @endcode 
def _rac_pop_ ( self , item = -1 ) :
    """ Remove and return item at index( or name) (default last).
    >>> lst  = ...
    >>> item = lst.pop( )   ## get and remove from collection the last   element 
    >>> item = lst.pop(2)   ## get and remove from collection the second element 
    >>> item = lst.pop('a') ## get and remove from collection the element with name 'a'
    """
    entry = self [ item ]
    self.remove ( entry ) 
    return entry 
    
for r in ( ROOT.RooArgSet  ,
           ROOT.RooArgList ) :
    r.front = _rac_front_
    r.back  = _rac_back_
    r.pop   = _rac_pop_
    
_new_methods_ += [
    ROOT.RooArgSet  . front ,
    ROOT.RooArgSet  . back  ,
    ROOT.RooArgSet  . pop   ,
    ROOT.RooArgList . front ,
    ROOT.RooArgList . back  ,
    ROOT.RooArgList . pop   ,
    ]

# =============================================================================
## check the presence of variable/name/index in set 
def _ras_contains_ ( self , aname ) :
    """ Check the presence of variable/name/index in set 
    """
    if isinstance ( aname , integer_types ) :
        return  0 <= aname < len ( self )

    if ( 6 , 31 ) <= root_info and isinstance ( aname , std.string ) : 
        _v = self.find ( str ( aname ) )
    else  : 
        _v = self.find (       aname   )
        
    if not _v : return False 
    return             True

# =============================================================================
## Check the presence of object/name/index in RooLinkedList
def _rll_contains_ ( self , what ) :
    """ Check the presence of object/name/index in RooLinkedList
    """ 
    if isinstance ( what , ROOT.RooAbsArg ) and valid_pointer ( what ) :
        return True if self.findArg    ( what   ) else False
    
    elif isinstance ( what , ROOT.TObject    ) and valid_pointer ( what ) :
        return True if self.FindObject ( what   ) else False
    
    elif isinstance ( what , string_types  ) :
        return True if self.find ( str ( what ) ) else False
    
    elif isinstance ( what , integer_types ) :
        return 0 <= what < len ( self )
    
    return False

# =============================================================================
## Get the object by index or name 
def _rll_getitem_ ( self , what ) :
    """ Get the object by index or name 
    """
    if   isinstance ( what , integer_types ) :
        l = len ( self )
        ## allow "slightly negative" indices:
        if what < 0 : what += l
        if 0 <= what < l : return self.At ( what )         
        raise IndexError ( 'List index is out of the range')
    elif isinstance ( what , string_types ) :
        obj = self.find ( str ( what ) )
        if valid_pointer (  obj  )  : return obj
        raise KeyError   ( "Key '%s' is not found")
    
    raise TypeError ( 'Unknown index type %s' % typename ( what ) )
    
# =============================================================================
## some decoration over RooArgSet 
ROOT.RooArgSet . __len__           = lambda s   : s.getSize() if valid_pointer ( s ) else 0 
ROOT.RooArgSet . __getattr__       = _ras_getattr_ 
ROOT.RooArgSet . __getitem__       = _ras_getitem_ 
ROOT.RooArgSet . __contains__      = _ras_contains_ 
ROOT.RooArgSet . __nonzero__       = lambda s : valid_pointer ( s ) and 0 < s.getSize()  
ROOT.RooArgSet . __bool__          = lambda s : valid_pointer ( s ) and 0 < s.getSize() 
        
ROOT.RooArgSet     . __str__       = lambda s : str ( set   ( _rs_list_ ( s ) ) ) if s else '{}'
ROOT.RooArgSet     . __repr__      = lambda s : str ( set   ( _rs_list_ ( s ) ) ) if s else '{}'

ROOT.RooLinkedList . __str__       = lambda s : str ( tuple ( _rs_list_ ( s ) ) ) if s else '[]'
ROOT.RooLinkedList . __repr__      = lambda s : str ( tuple ( _rs_list_ ( s ) ) ) if s else '[]'
ROOT.RooLinkedList . __contains__  = _rll_contains_ 
ROOT.RooLinkedList . __getitem__   = _rll_getitem_ 

ROOT.RooAbsCollection.__len__      = lambda s   : s.getSize() if valid_pointer ( s ) else 0 
ROOT.RooAbsCollection. __nonzero__ = lambda s   : valid_pointer ( s ) and 0 < s.getSize()
ROOT.RooAbsCollection. __bool__    = lambda s   : valid_pointer ( s ) and 0 < s.getsize()  

# ========================================================================================
## check it the item is in the list 
def _stl_contains_ ( lst , item ) :
    """ Check if the item is in the list 
    """
    if   isinstance ( item , integer_types  ) : return 0 <= item < len ( lst ) 
    elif isinstance ( item , ROOT.RooAbsArg ) :
        return any ( ( i is item ) or ( i == item ) or ( i.name == item.name ) for i in lst )
    elif isinstance ( item , string_types ) :
        return any (  i.name == item  for i in lst )    
    return False 


_STLList = ROOT.RooSTLRefCountList(ROOT.RooAbsArg)
_STLList .__str__      = lambda s : str ( tuple ( _rs_list_ ( s ) ) ) if s else '[]'
_STLList .__repr__     = lambda s : str ( tuple ( _rs_list_ ( s ) ) ) if s else '[]'
_STLList .__contains__ = _stl_contains_ 

if not hasattr ( ROOT.RooArgList , '__iter__' ) or root_info < ( 6 , 31 ) : 
    ROOT.RooArgList . __iter__      = _ral_iter_
    _new_methods_ += [ ROOT.RooArgList. __iter__ ]
    
if not hasattr ( ROOT.RooArgSet , '__iter__' ) or root_info < ( 6 , 31 ) : 
    ROOT.RooArgSet . __iter__          = _ras_iter_ 
    _new_methods_ += [ ROOT.RooArgList. __iter__ ]

if not hasattr ( ROOT.RooAbsCollection , '__iter__' ) or root_info < ( 6 , 31 ) : 
    ROOT.RooAbsCollection.__iter__     = _ras_iter_
    _new_methods_ += [ ROOT.RooAbsCollection. __iter__ ]

# =============================================================================
## iterator for class RooLinkedList
#  @code
#  lst = ...
#  for l in lst : print l 
#  @endcode 
def _rll_iter_  ( self ) :
    """ Iterator over RooLinekdList
    >>> lst = ...
    >>> for l in lst : print l     
    """
    l = len ( self )  
    for i in range ( l ) : yield self.At ( i )

if not hasattr ( ROOT.RooLinkedList , '__iter__' ) or root_info < ( 6, 31 ) :
    ROOT.RooLinkedList . __iter__  = _rll_iter_ 
    _new_methods_ += [ ROOT.RooLinkedList . __iter__ ]
        
_new_methods_ += [
    ROOT.RooArgSet . __len__      ,
    ROOT.RooArgSet . __getattr__  , 
    ROOT.RooArgSet . __setattr__  ,
    ROOT.RooArgSet . __contains__ ,
    ROOT.RooArgSet . __nonzero__  ,
    ROOT.RooArgSet . __bool__     ,
    ROOT.RooArgSet . __str__      ,
    ROOT.RooArgSet . __repr__     ,
    ]

ROOT.RooLinkedList.add     = ROOT.RooLinkedList.Add
ROOT.RooLinkedList.__len__ = ROOT.RooLinkedList.GetSize

_new_methods_ += [
    ROOT.RooLinkedList . __len__  ,
    ROOT.RooLinkedList . __repr__ ,
    ROOT.RooLinkedList . add 
    ]
    
# =============================================================================
## add more data into list/set
def _ral_iadd_ ( self , other ) :
    """ Update/increment collections
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
    """ Make a sum of two lists/sets/collections
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
    t. clone    = _ral_clone_ 
    t. clear    = _ral_clear_ 
    t. __add__  = _ral_add_
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
    """ Get an intersection of two sets
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
    """ Get an union of two sets
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
    """ Get a difference fot two sets
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
    """ Get a difference fot two sets
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
## Same content of two containers (by names!)
def _rac_same_ ( cnt , other ) :
    """ same content of two containers (by names!)
    """

    if other is cnt : return True
    assert isinstance ( other . ROOT.RooAbsCollection ), \
        'Invaild sequence type'
    
    a1 = set ( v.name for v in cnt   )
    a2 = set ( v.name for v in other )
    diff = a1 ^ a2
    
    return True if not diff else False

# =============================================================================
## Non-Equality of content of two containers (by names!)
def _rac_diff_ ( cnt , other ) :
    """ Non-Equality of content of two containers (by names!)
    """

    if other is cnt : return False 
    assert isinstance ( other . ROOT.RooAbsCollection ), \
        'Invaild sequence type'

    a1 = set ( v.name for v in cnt   )
    a2 = set ( v.name for v in other )
    diff = a1 ^ a2
    
    return True if diff else False

ROOT.RooAbsCollection. same      = _rac_same_ 
ROOT.RooAbsCollection. different = _rac_diff_

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
    """ Simple context manager for temporary redefitnnnonn of some mutbale collection
    
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
## attention 
if not hasattr ( ROOT.RooAbsCollection , 'assign' ) :
    ## attention !!!
    ROOT.RooAbsCollection.assign = ROOT.RooAbsCollection.assignValueOnly 

# =============================================================================
__item_store = set() 
# =============================================================================
## Unpickle RooArgList/RooArgSet
#  @attention it stores the constructore arguments  in the local store 
#  @see RooArgList
def _rac_factory ( klass , *args ) :
    """ Unpickle RooArgList/RooArgSet  instance
    0 attention: it stores the constructore arguments  in the local store 
    """
    c = klass () 
    for a in args : c.add ( a )
    __item_store.add ( args ) 
    return c
# =============================================================================
## reduce `RooArgList`
def _ral_reduce_ ( rac ) :
    """ Reduce `RooArgList` instances"""    
    return _rac_factory, ( ROOT.RooArgList, ) + tuple ( a for a in rac ) 

# =============================================================================
## reduce `RooArgSet`
def _ras_reduce_ ( rac ) :
    """ Reduce `RooArgSet` instances"""    
    return _rac_factory, ( ROOT.RooArgSet, ) + tuple ( a for a in rac )  

ROOT.RooArgSet  .__reduce__ = _ras_reduce_
ROOT.RooArgList .__reduce__ = _ral_reduce_

_new_methods_ += [
    ROOT.RooArgSet  .__reduce__ , 
    ROOT.RooArgList .__reduce__ ,
    ]

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
