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
# ROOT.TCut
# =============================================================================
ROOT.TCut.__str__      = lambda s :       s.GetTitle().strip() 
ROOT.TCut.__repr__     = lambda s :       s.GetTitle().strip() 
ROOT.TCut.__nonzero__  = lambda s : bool( s.GetTitle().strip() )

# =============================================================================
## Remove leading/traling and excessive blanks from TCut
#  @code 
#  cut = ...
#  cut.strip()
#  @endcode
#  @see ROOT::TCut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_strip_ ( self ) :
    """Remove leading/trailing and excessive blanks from TCut
    >>> cut = ...
    >>> cut.strip() 
    """
    ## 1. remove  starting and trailing spaces 
    t = self.GetTitle().strip()
    ## 2. remove excessive internal spaces 
    while 0 <= t.find ( '  ' ) : t = t.replace ( '  ' , ' ' )
    ## 3. treat empty expressions 
    if t in  ( '()' , '( )' )  : t = ''
    ## ...
    self.SetTitle ( t ) 
    return self

# =============================================================================
## create new cut by replacing the expressions
#  @code 
#  oldcut = ...
#  newcut = oldcut.replace ( 'pt_Bu' , 'pt_Bc' )
#  @endcode 
#  @see ROOT::TCut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_replace_ ( self , oldexp , newexp ) :
    """Create new cut by replacing the expressions  
    >>> oldcut = ...
    >>> newcut = oldcut.replace ( 'pt_Bu' , 'pt_Bc' ) 
    """
    t = self.strip ()
    t = ROOT.TCut ( t.replace ( oldexp , newexp ) ) 
    return ROOT.TCut ( t.strip () ) 

ROOT.TCut . strip      = _tc_strip_
ROOT.TCut . replace    = _tc_replace_

# =============================================================================
## Logical *AND* for TCut objects 
#  @see ROOT::TCut
#  cut       = ...
#  other_cut = ...
#  cut      &= other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_iand_ ( self , other ) :
    """Logical *AND* for TCut objects 
    >>> cut        = ...
    >>> other_cut  = ... 
    >>> cut       &= other_cut 
    """
    if not isinstance ( other , ( str , ROOT.TCut ) ) : return NotImplemented
    ##
    self.strip() 
    ##
    other = ROOT.TCut ( other.strip() )
    if   other and self : 
        self.SetTitle("(%s)&&(%s)" % ( self , other ) )
        return self
    elif other :
        self.SetTitle ( "%s" % other )
        return self
    ## 
    return self 

# =============================================================================
## Logical *OR* for TCut objects 
#  @see ROOT::TCut
#  cut       = ...
#  other_cut = ...
#  cut      |= other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_ior_ ( self , other ) :
    """Logical *OR* for TCut objects 
    >>> cut        = ...
    >>> other_cut  = ... 
    >>> cut       |= other_cut ## create new cut 
    """
    if not isinstance ( other , ( str , ROOT.TCut ) ) : return NotImplemented
    ##
    self.strip() 
    ##
    other = ROOT.TCut ( other.strip() )
    ## 
    if   other and self : 
        self.SetTitle("(%s)||(%s)" % ( self , other ) )
        return self
    elif other :
        self.SetTitle ( "%s" % other )
        return self
    ## 
    return self 

# =============================================================================
## Multiplication for TCut objects 
#  @see ROOT::TCut
#  cut       = ...
#  other_cut = ...
#  cut      &= other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_imul_ ( self , other ) :
    """Multiplication for TCut objects 
    >>> cut    = ...
    >>> other  = ... 
    >>> cut   *= other 
    """
    ## 
    self.strip()
    if not self : return NotImplemented 
    ##
    if   isinstance ( other , ( float , int , long ) ) :
        self.SetTitle( "(%s)*%s"  %  ( self , other ) ) 
        return self
    elif isinstance ( other , ( str , ROOT.TCut ) ) :
        other = ROOT.TCut ( other.strip() )
        if not other : return NotImplemented 
    else :
        return NotImplemented 
                
    self.SetTitle("(%s)*(%s)" % ( self , other ) )
    return self

# =============================================================================
## Multiplication for TCut objects 
#  @see ROOT::TCut
#  cut       = ...
#  other_cut = ...
#  cut      /= other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_idiv_ ( self , other ) :
    """Multiplication for TCut objects 
    >>> cut    = ...
    >>> other  = ... 
    >>> cut   *= other 
    """
    ## 
    self.strip()
    if not self : return NotImplemented 
    ##
    if   isinstance ( other , ( float , int , long ) ) :
        self.SetTitle( "(%s)/%s"  %  ( self , other ) ) 
        return self
    elif isinstance ( other , ( str , ROOT.TCut ) ) :
        other = ROOT.TCut ( other.strip() )
        if not other : return NotImplemented 
    else :
        return NotImplemented 
                
    self.SetTitle("(%s)/(%s)" % ( self , other ) )
    return self

# =============================================================================
## Logical *AND* for TCut objects  
#  @see ROOT::TCut
#  cut       =
#  other_cut = 
#  new_cut   = cut & other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_and_ ( self , other ) :
    """Logical *AND* for TCut objects 
    >>> cut        = ...
    >>> other_cut  = ... 
    >>> new_cut    = cut & other_cut 
    """
    ##
    if not isinstance ( other , ( str , ROOT.TCut ) ) : return NotImplemented
    ## 
    new_cut  = ROOT.TCut ( self )
    new_cut &= other 
    return new_cut

# =============================================================================
## Logical *OR* for TCut objects 
#  @see ROOT::TCut
#  cut       =
#  other_cut = 
#  new_cut   = cut | other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_or_ ( self , other ) :
    """Logical *OR* for TCut objects 
    >>> cut        = ...
    >>> other_cut  = ... 
    >>> new_cut    = cut | other_cut 
    """
    ##
    if not isinstance ( other , ( str , ROOT.TCut ) ) : return NotImplemented
    ## 
    new_cut  = ROOT.TCut ( self )
    new_cut |= other 
    return new_cut
# =============================================================================
## Multiplication for TCut objects 
#  @see ROOT::TCut
#  cut       =
#  other_cut = 
#  new_cut   = cut & other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_mul_ ( self , other ) :
    """Multiplication for TCut objects 
    >>> cut        = ...
    >>> other_cut  = ... 
    >>> new_cut    = cut * other_cut 
    """
    ##
    if not isinstance ( other , ( str , ROOT.TCut , float , int , long  ) ) :
        return NotImplemented
    ## 
    new_cut  = ROOT.TCut ( self )
    new_cut *= other 
    return new_cut

# =============================================================================
## Division for TCut objects 
#  @see ROOT::TCut
#  cut       =
#  other_cut = 
#  new_cut   = cut / other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_div_ ( self , other ) :
    """Division for TCut objects 
    >>> cut        = ...
    >>> other_cut  = ... 
    >>> new_cut    = cut / other_cut 
    """
    ##
    if not isinstance ( other , ( str , ROOT.TCut , float , int , long  ) ) :
        return NotImplemented
    ## 
    new_cut  = ROOT.TCut ( self )
    new_cut /= other 
    return new_cut

# =============================================================================
## minor extension for TCut
#  @see ROOT::TCut
#  cut       = ...
#  new_cut   = 'pt>1' & other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_rand_ ( self , other ) :
    """Logical *AND* for TCut objects 
    >>> cut        = ...
    >>> new_cut    = 'pt>1' & other_cut 
    """
    if not isinstance ( other , ( str , ROOT.TCut ) ) : return NotImplemented
    return ROOT.TCut(other) & self  

# =============================================================================
## minor extension for TCut
#  @see ROOT::TCut
#  cut       = ...
#  new_cut   = 'pt>1' | other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_ror_ ( self , other ) :
    """Logical *AND* for TCut objects 
    >>> cut        = ...
    >>> new_cut    = 'pt>1' | other_cut 
    """
    if not isinstance ( other , ( str , ROOT.TCut ) ) : return NotImplemented
    return ROOT.TCut(other) | self  
# =============================================================================
## Multiplication for TCut objects 
#  @see ROOT::TCut
#  cut       =
#  other_cut = 
#  new_cut   = cut * other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_rmul_ ( self , other ) :
    """Multiplication for TCut objects 
    >>> cut        = ...
    >>> other_cut  = ... 
    >>> new_cut    = other_cut * cut
    """
    ##
    if isinstance ( other , ( str , ROOT.TCut ) ) :
        other = ROOT.TCut( other.strip() )
    elif isinstance ( other ,  ( float , int, long ) ) :
        other = ROOT.TCut ("%s"   % other )
    else :
        return NotImplemented

    if not other : return NotImplemented 
    if not self  : return NotImplemented
    
    return other * self 

# =============================================================================
## Division for TCut objects 
#  @see ROOT::TCut
#  cut       =
#  other_cut = 
#  new_cut   = other_cut /  div 
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_rdiv_ ( self , other ) :
    """Multiplication for TCut objects 
    >>> cut        = ...
    >>> other_cut  = ... 
    >>> new_cut    = other_cut / cut
    """
    ##
    if isinstance ( other , ( str , ROOT.TCut ) ) :
        other = ROOT.TCut( other.strip() )
    elif isinstance ( other ,  ( float , int, long ) ) :
        other = ROOT.TCut ("%s"   % other )
    else :
        return NotImplemented
    
    if not other : return NotImplemented 
    if not self  : return NotImplemented
    
    return other / self 

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
    return ROOT.TCut

ROOT.TCut. __and__    = _tc_and_
ROOT.TCut.__iand__    = _tc_iand_
ROOT.TCut.__rand__    = _tc_rand_
ROOT.TCut.  __or__    = _tc_or_
ROOT.TCut. __ior__    = _tc_ior_
ROOT.TCut. __ror__    = _tc_ror_
ROOT.TCut. __mul__    = _tc_mul_
ROOT.TCut.__imul__    = _tc_imul_
ROOT.TCut.__rmul__    = _tc_rmul_
ROOT.TCut. __div__    = _tc_div_
ROOT.TCut.__idiv__    = _tc_idiv_
ROOT.TCut.__rdiv__    = _tc_rdiv_
ROOT.TCut.__invert__  = _tc_invert_


# =============================================================================
_decorated_classes_ = (
    ROOT.TCut ,
    )
_new_methods_       = (
    #
    ROOT.TCut . __add__    ,
    ROOT.TCut . __and__    ,
    ROOT.TCut . __or__     ,
    ROOT.TCut . __mul__    ,
    ROOT.TCut . __invert__ ,
    #
    ROOT.TCut . strip      ,
    ROOT.TCut . replace    ,    
    )
# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
