#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with decoration of ROOT.TCut object
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
#
# =============================================================================
"""Decoration of some ROOT.TCut objects for efficient use in python"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = () ## nothing to import 
# =============================================================================
import ROOT
from   ostap.core.core import cpp, VE, hID, dsID   
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.trees.cuts' )
else                       : logger = getLogger( __name__           )
# =============================================================================
logger.debug( 'Some useful decorations for ROOT.TCut objects')
# =============================================================================

# =============================================================================
## minor extension for TCut
#  @code
#  cut       =
#  other_cut = 
#  new_cut1  = cut + other_cut ## create new cut  
#  new_cut2  = cut & other_cut ## ditto 
#  @endcode 
#  @see ROOT::TCut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_and_( c , other ) :
    """Add cuts using logical AND
    >>> cut       =
    >>> other_cut = 
    >>> new_cut1  = cut + other_cut ## create new cut  
    >>> new_cut2  = cut & other_cut ## ditto 
    """
    #
    nc  = ROOT.TCut( c.strip() )
    nc += other.strip()
    #
    return nc 

# =============================================================================
## minor extension for TCut
#  @code 
#  cut       =
#  other_cut = 
#  new_cut1  = cut * other_cut ## create new cut
#  @endcode 
#  @see ROOT::TCut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_mul_( c , other ) :
    """Multiply cuts values
    >>> cut       =
    >>> other_cut = 
    >>> new_cut1  = cut * other_cut ## create new cut  
    """
    #
    nc  = ROOT.TCut( c.strip() )
    nc *= other.strip()
    #
    return nc 

# =============================================================================
## minor extension for TCut
#  @see ROOT::TCut
#  cut       =
#  other_cut = 
#  new_cut1  = cut | other_cut ## create new cut  
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_or_ ( c , other ) :
    """Logical OR for TCut 
    >>> cut       =
    >>> other_cut = 
    >>> new_cut1  = cut | other_cut ## create new cut  
    """
    #
    nc   = ROOT.TCut ( c     .strip() )
    oc   = ROOT.TCut ( other .strip() )
    #
    ntit = nc.GetTitle()
    otit = oc.GetTitle()
    #
    if   ntit and otit : return ROOT.TCut ( "(%s)||(%s)" % ( ntit , otit ) )
    elif ntit          : return nc 
    elif otit          : return oc  
    ##
    return ROOT.TCut()

# =============================================================================
## minor extension for TCut
#  @code 
#  cut     =
#  new_cut = ~cut
#  @endcode 
#  @see ROOT::TCut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_invert_ ( c ) :
    """Invert the cut 
    >>> cut     =
    >>> new_cut = ~cut 
    """
    nc  = ROOT.TCut ( c .strip() )
    tit = nc.GetTitle() 
    if tit : return ROOT.TCut( '!(%s)' % tit )
    ##
    return ROOT.TCut()

# =============================================================================
## Remove leading/traling and excessive blanks from TCut
#  @code 
#  cut = ...
#  cut.strip()
#  @endcode
#  @see ROOT::TCut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_strip_ ( c ) :
    """Remove leading/trailing and excessive blanks from TCut
    >>> cut = ...
    >>> cut.strip() 
    """
    t = c.GetTitle()
    t = t.strip()
    while 0 <= t.find ( '  ' ) : t = t.replace ( '  ' , ' ' )
    c.SetTitle ( t ) 
    return c

# =============================================================================
## create new cut by replacing the expressions
#  @code 
#  oldcut = ...
#  newcut = oldcut.replace ( 'pt_Bu' , 'pt_Bc' )
#  @endcode 
#  @see ROOT::TCut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_replace_ ( c , oldexp , newexp ) :
    """Create new cut by replacing the expressions  
    >>> oldcut = ...
    >>> newcut = oldcut.replace ( 'pt_Bu' , 'pt_Bc' ) 
    """
    t = c.strip ()
    t = ROOT.TCut ( t.replace ( oldexp , newexp ) ) 
    return ROOT.TCut ( t.strip () ) 
                       
ROOT.TCut . __add__    = _tc_and_
ROOT.TCut . __and__    = _tc_and_
ROOT.TCut . __or__     = _tc_or_
ROOT.TCut . __mul__    = _tc_mul_
ROOT.TCut . __invert__ = _tc_invert_

ROOT.TCut . strip      = _tc_strip_
ROOT.TCut . replace    = _tc_replace_

ROOT.TCut . __rand__   = lambda s,o : ROOT.TCut(o)&s
ROOT.TCut . __radd__   = lambda s,o : ROOT.TCut(o)+s   
ROOT.TCut . __ror__    = lambda s,o : ROOT.TCut(o)|s  
ROOT.TCut . __rmul__   = lambda s,o : ROOT.TCut(o)*s   

ROOT.TCut . __repr__   = ROOT.TCut.__str__ 

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
