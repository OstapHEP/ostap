#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/printable.py
#  Module with decoration for RooPrintable 
#  @see RooPrintable 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
""" Module with decoration for RooPrintable 
- see RooPrintable 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    ) 
# =============================================================================
from   ostap.core.core        import Ostap, valid_pointer  
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger    import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.printable' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug( 'Some useful decorations for RooPrintable objects')
# =============================================================================
_new_methods_ = []
# =============================================================================
##  Use RooPrintable::printMultiline function
#   @see RooPrintable 
#   @see Ostap::Utils::print_printable   
def print_multiline ( o , content = 1 , verbose = False , indent = '' ) :
    """ Use RooPrintable::printMultiline function
    - see RooPrintable 
    - see Ostap::Utils::print_printable1 
    """
    if not valid_pointer ( o ) : return 'Invalid object'
    return Ostap.Utils.print_printable1 ( o , content , verbose , indent )
# =============================================================================
##  Use RooPrintable::printStream function
#   @see RooPrintable 
#   @see Ostap::Utils::print_printable2   
def print_stream  ( o , content = 1 , style = 3 , indent = '' ) :
    """ Use RooPrintable::printStream function
    - see RooPrintable 
    - see Ostap::Utils::print_printable2
    """
    if not valid_pointer ( o ) : return 'Invalid object'
    return Ostap.Utils.print_printable2 ( o , content , style , indent )
# =============================================================================
##  Use RooPrintable::printTree function
#   @see RooPrintable 
#   @see Ostap::Utils::print_printable_tree   
def print_tree ( o , indent = '' ) :
    """ Use RooPrintable::printTree function
    - see RooPrintable 
    - see Ostap::Utils::print_printable_tree
    """
    if not valid_pointer ( o ) : return 'Invalid object'
    return Ostap.Utils.print_printable_tree ( o , indent )

# =============================================================================
## make easy print for RooPrintable 
def _rp_print_ ( obj , opts = 'vv' , *style ) :
    """Make easy print for RooPrintable
    >>> o = ...
    >>> print o 
    """
    return Ostap.Utils.print_printable ( obj , opts , *style )

ROOT.RooPrintable.print_multiline = print_multiline
ROOT.RooPrintable.print_stream    = print_stream 
ROOT.RooPrintable.print_tree      = print_tree 
ROOT.RooPrintable.print_printable = _rp_print_ 
ROOT.RooPrintable.__str__         = _rp_print_
ROOT.RooPrintable.__repr__        = _rp_print_

_new_methods_ += [
    ROOT.RooPrintable.print_printable ,
    ROOT.RooPrintable.print_multiline ,
    ROOT.RooPrintable.print_stream    ,
    ROOT.RooPrintable.print_tree      ,
    ROOT.RooPrintable.__str__         ,
    ROOT.RooPrintable.__repr__        ,
    ]
# =============================================================================
_decorated_classes_ = (
    ROOT.RooPrintable  ,
    )

_new_methods_ = tuple ( _new_methods_ ) 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
