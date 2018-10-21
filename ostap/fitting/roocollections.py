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
    ) 
#__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    ) 
# =============================================================================
import ROOT, random
from   ostap.core.core         import Ostap
import ostap.fitting.variables
# =============================================================================
# logging \
# =============================================================================
from ostap.logger.logger import getLogger , allright,  attention
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.roocollections' )
else                       : logger = getLogger( __name__ )
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
    if isinstance ( i , int ) : return 0<= i < len(self)
    return  0 <= self.index ( i )

# =============================================================================
##  get item form the list 
def _ral_getitem_ ( self , index ) :
    """Get item form the list
    >>>  lst  = ...
    >>>  item = lst[3]
    """
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

## ============================================================================
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
## printout for   RooArgList 
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
ROOT.RooArgSet . __len__       = lambda s   : s.getSize()
ROOT.RooArgSet . __iter__      = _ras_iter_ 
ROOT.RooArgSet . __getattr__   = _ras_getattr_ 
ROOT.RooArgSet . __getitem__   = _ras_getitem_ 
ROOT.RooArgSet . __contains__  = _ras_contains_ 
ROOT.RooArgSet . __nonzero__   = lambda s   : 0 != len ( s ) 
        
ROOT.RooArgSet     . __str__   = lambda s : str ( tuple ( _rs_list_ ( s ) ) )  
ROOT.RooArgSet     . __repr__  = lambda s : str ( tuple ( _rs_list_ ( s ) ) )  
ROOT.RooLinkedList . __repr__  = lambda s : str (  _rs_list_ ( s ) )
ROOT.RooLinkedList . __iter__  = _ras_iter_ 

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

ROOT.RooLinkedList.add = ROOT.RooLinkedList.Add

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
    from collections import Container as _CNT
    _RAC = ROOT.RooAbsCollection
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
    from collections import Container as _CNT
    _RAC = ROOT.RooAbsCollection
    _RAA = ROOT.RooAbsArg  
    if not isinstance ( other , ( _CNT, _RAC , _RAA ) ) : return NotImplemented
    if     isinstance ( other , str )                   : return NotImplemented
    _clone = self.clone('')
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
    return self.Clone(name)


# =============================================================================
for t in ( ROOT.RooArgList , ROOT.RooArgSet , ROOT.RooLinkedList ) :
    t. clone    =  _ral_clone_ 
    t. __add__  =  _ral_add_
    t.__iadd__  = _ral_iadd_
    t.__radd__  = _ral_radd_
    t.append    = _ral_iadd_

_new_methods_ += [
    ROOT.RooArgList    . clone    ,
    ROOT.RooArgSet     . clone    ,
    ROOT.RooLinkedList . clone    ,
    ROOT.RooArgList    . __add__  ,
    ROOT.RooArgSet     . __add__  ,
    ROOT.RooLinkedList . __add__  ,
    ROOT.RooArgList    . __iadd__ ,
    ROOT.RooArgSet     . __iadd__ ,
    ROOT.RooLinkedList . __iadd__ ,
    ROOT.RooArgList    . __radd__ ,
    ROOT.RooArgSet     . __radd__ ,
    ROOT.RooLinkedList . __radd__ ,
    ROOT.RooArgList    . append   ,
    ROOT.RooArgSet     . append   ,
    ROOT.RooLinkedList . append   ,
    ]


# =============================================================================
_decorated_classes_ = (
    ROOT.RooArgSet     , 
    ROOT.RooArgList    , 
    ROOT.RooLinkedList , 
    )

_new_methods_ = tuple ( _new_methods_ ) 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
