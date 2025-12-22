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
""" Decoration of some ROOT.TCut objects for efficient use in python"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'expression_types' , ## valid types for expressions/selections/weights
    'vars_and_cuts'    , ## helper routibe to treat expressions
)
# =============================================================================
from   ostap.core.meta_info   import ostap_info
from   ostap.core.ostap_types import num_types, string_types, integer_types 
from   ostap.core.core        import cpp, VE, hID, dsID, Ostap
from   ostap.utils.strings    import split_string 
from   ostap.utils.utils      import balanced 
from   math                   import isfinite, isnan
import ROOT, ast, re
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.trees.cuts' )
else                       : logger = getLogger( __name__           )
# =============================================================================
logger.debug( 'Some useful decorations for ROOT.TCut objects')
# =============================================================================
## warning about the order of variables for `project`
order_warning = ostap_info < ( 2 , 0 )
if order_warning :
    import datetime
    order_warning = datetime.datetime.now () < datetime.datetime ( 2025, 7 , 1 )
# =============================================================================
## float aformat for conversion of floats into cut-expressions
fmt_float = '%.16g'
fmt_int   = '%d'
fmt_other = '%s'
# =============================================================================
## types for expressions and cuts 
expression_types  = string_types +  ( ROOT.TCut , )
# ==============================================================================
##  Is the expression/cut  really trivial ? 
#  @see Ostap::tivial 
def trivial ( expression ) :
    """ Is the expression/cut  really trivial ? 
    - see Ostap.tivial 
    """
    return expression and Ostap.trivial ( expression )

# =============================================================================
## Prepare the arguments: variable names and cuts 
#  @code
#  vars_lst, cuts, input_string = vars_and_cuts ( 'x', 'y<0' )
#  vars_lst, cuts, input_string = vars_and_cuts ( ('x','y','z') , 'y<0' )
#  @endcode
def vars_and_cuts ( expressions , cuts , allow_empty = False ) :
    """ Prepare the arguments: variable names and cuts 
    >>> vars_lst, cuts, input_string = vars_and_cuts ( 'x', 'y<0' )
    >>> vars_lst, cuts, input_string = vars_and_cuts ( ('x','y','z') , 'y<0' )
    """
    
    ## single string as input 
    input_string = False 
    if isinstance ( expressions , expression_types ) :
        expressions  = split_string ( str ( expressions ) , strip = True , respect_groups = True )
        input_string = 1 == len ( expressions )


    # ===========================================================================
    ## CUTS
    # ===========================================================================
    
    assert isinstance ( cuts , expression_types ) or not cuts , \
        'Invalid type of cuts: %s' % str ( cuts ) 
    
    cuts = str ( cuts ).strip() if cuts else ''

    ## if cut is trivial : make it really trivial! 
    if cuts and trivial ( cuts ) : cuts = ''

    # ============================================================================
    ## VARIABLES 
    # ============================================================================

    ## EMPTY ?
    if allow_empty and not expressions : return () , cuts , input_string ## ATTENTION!
        
    assert expressions and all ( s and isinstance ( s , expression_types ) for s in expressions ) , \
        "Invalid expression(s) : %s" % str ( expressions )

    full = []
    for subexpr in expressions :
        subexpr  = split_string ( str( subexpr ).strip()  , strip = True , respect_groups = True )
        if subexpr : full    += subexpr
    expressions  = tuple ( full ) 

    ## EMPTY ?
    if allow_empty and not expressions : return () , cuts , input_string ## ATTENTION!! 
    
    assert expressions and all ( s and isinstance ( s , expression_types ) for s in expressions ) , \
        "Invalid expression(s) : %s" % str ( expressions )
    
    exprs = tuple ( str(s).strip() for s in expressions )

    exprs = tuple ( s for s in exprs if s )

    ## EMPTY ? 
    if allow_empty and not exprs : return () , cuts , input_string ## ATTENTION!! 
        
    assert exprs and all ( exprs ) , "Invalid expressions: %s" % str ( exprs )
    
    return exprs , cuts, input_string 

# =============================================================================
# ROOT.TCut
# =============================================================================
symbols = '+-=*/)(&|^%!~:[]{}$<>/ '
## printout of TCut (add extra parentheses if/when needed
def _tc_str_ ( cut ) :
    """ {rintout of TCut (add extra parentheses if/when needed
    """
    cut.strip()
    t = cut.GetTitle ()
    if any ( s in t for s in symbols ) : return '(%s)' % t 
    return t

ROOT.TCut.__str__      = _tc_str_ 
ROOT.TCut.__repr__     = _tc_str_ 
ROOT.TCut.__nonzero__  = lambda s : bool ( s.GetTitle().strip() )
ROOT.TCut.__bool__     = lambda s : bool ( s.GetTitle().strip() )

# =============================================================================
## Parse the cut using ast and then unparse it back
#  @code
#  cut = ...
#  print ( cut.ast() )
#  @endcode
def _tc_ast_ ( cut ) :  
    """ Parse the cut expression using ast and the unparse it back
    >>> cut = ...
    >>> print ( cut.ast () )
    """
    cut.strip()
    t = cut.GetTitle ()
    t = t.replace ( '&&'       , ' and '    )
    t = t.replace ( '||'       , ' or '     )

    if '!' in t : 
        t = t.replace ( '!='            , ' __TCUT_NE__ ' ) ## ATTENTION 
        t = t.replace ( '!'             , '~'             ) ## ATTERNTION 
        t = t.replace ( ' __TCUT_NE__ ' , '!='            ) ## ATTENTION 

    t = t.replace ( '::' , '__DOUBLE_COLUMN__' )
        
    p = ast.parse   ( t )
    t = ast.unparse ( p )
    
    t = t.replace ( '~'        , '!'        )
    t = t.replace ( ' or '     , ' || '     )
    t = t.replace ( ' and '    , ' && '     )
    t = t.replace ( '__DOUBLE_COLUMN__' , '::' )  

    return t 
    
ROOT.TCut.ast       = _tc_ast_     

# =============================================================================
## modified constructor:
#  - strip
#  - check parentheses 
def _tc_new_init_ ( self , *args ) :
    """ Modified constructor:
    - strip
    - check parentheses 
    """
    self._old_init_ ( *args ) 
    ## strip it 
    cut = self.GetTitle().strip()
    ##
    while '  ' in cut : cut = cut.replace ( '  ' , ' ' )
    if cut in  ( '()' , '( )' )  : cut = ''
    ##
    if any ( s in cut for s in symbols ) :
        if cut.startswith('(') and cut.endswith(')') : pass
        else : cut = '(%s)' % cut
    ## 
    self.SetTitle ( cut  )
    ##
    assert balanced ( cut , left = '[(' , right = '])' ) , \
           'TCut: unbalanced parentheses/square brackets:"%s"' % cut
    
if not hasattr ( ROOT.TCut , '_old_init_' ) :
    ROOT.TCut._old_init_ =  ROOT.TCut.__init__
    ROOT.TCut._new_init_ =  _tc_new_init_ 
    ROOT.TCut.__init__   =  _tc_new_init_ 
    
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
    """ Remove leading/trailing and excessive blanks from TCut
    >>> cut = ...
    >>> cut.strip() 
    """
    t = self.GetTitle().strip()
    ## 
    while '  ' in t : t = t.replace ( '  ' , ' ' )
    if t in  ( '()' , '( )' )  : t = ''
    ## 
    self.SetTitle ( t )
    ##
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
    """ Logical *AND* for TCut objects 
    >>> cut        = ...
    >>> other_cut  = ... 
    >>> cut       &= other_cut 
    """
    self.strip() 
    if       isinstance ( other , string_types ) : other = ROOT.TCut ( other ) 
    elif     other is True                       : return self 
    elif     other is False                      : return ROOT.TCut ( '0' )
    ## 
    if not isinstance ( other , ROOT.TCut    ) : return NotImplemented
    ## 
    if   other and self  :
        self.SetTitle("%s&&%s" % ( self , other ) )
        return self
    elif other : self.SetTitle ( other.GetTitle () )
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
    """ Logical *OR* for TCut objects 
    >>> cut        = ...
    >>> other_cut  = ... 
    >>> cut       |= other_cut ## create new cut 
    """
    self.strip() 
    if       isinstance ( other , string_types ) : other = ROOT.TCut ( other )
    elif     other is True                       : return ROOT.TCut ( '1<2' )  
    elif     other is False                      : return self 
    ##  
    if not isinstance ( other , ROOT.TCut    ) : return NotImplemented
    ## 
    if   other and self  :
        self.SetTitle("%s||%s" % ( self , other ) )
        return self
    elif other : self.SetTitle ( other.GetTitle () )
    ## 
    return self 

# =============================================================================
## Multiplication for TCut objects 
#  @see ROOT::TCut
#  cut       = ...
#  other_cut = ...
#  cut      *= other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_imul_ ( self , other ) :
    """ Multiplication for TCut objects 
    >>> cut    = ...
    >>> other  = ... 
    >>> cut   *= other 
    """
    self.strip() 
    if       isinstance  ( other , string_types  ) : other = ROOT.TCut ( other )
    elif     isinstance  ( other , num_types     ) :
        ## 
        if not other :
            self.SetTitle ( '0' ) 
            return self 
        if 1 == other  : return self
        ##
        if   isinstance  ( other , integer_types  ) : other = ROOT.TCut ( fmt_int   % other )
        elif isinstance  ( other , float          ) :
            if not isfinite ( other ) : return NotImplemented 
            other = ROOT.TCut ( fmt_float % other )
        else                                        : other = ROOT.TCut ( fmt_other % other )
        ## 
    if not isinstance ( other , ROOT.TCut     ) : return NotImplemented
    ## 
    if   other and self  :
        self.SetTitle("%s*%s" % ( self , other ) )
        return self
    elif other : self.SetTitle ( other.GetTitle () )
    ## 
    return self 

# =============================================================================
## Addition for TCut objects 
#  @see ROOT::TCut
#  cut       = ...
#  other_cut = ...
#  cut      += other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_iadd_ ( self , other ) :
    """ Addition for TCut objects 
    >>> cut    = ...
    >>> other  = ... 
    >>> cut   += other 
    """
    self.strip() 
    if       isinstance  ( other , string_types  ) : other = ROOT.TCut ( other ) 
    elif     isinstance  ( other , num_types     ) :
        ##
        if not other : return self
        ## 
        if   isinstance  ( other , integer_types ) : other = ROOT.TCut ( fmt_int   % other )
        elif isinstance  ( other , float         ) :
            if not isfinite ( other ) : return NotImplemented             
            other = ROOT.TCut ( fmt_float % other ) 
        else                                       : other = ROOT.TCut ( fmt_other % other )  
    ## 
    if not isinstance ( other , ROOT.TCut     ) : return NotImplemented
    ## 
    if   other and self  :
        self.SetTitle("%s+%s" % ( self , other ) )
        return self
    elif other : self.SetTitle ( other.GetTitle () )
    ## 
    return self 

# =============================================================================
## Division for TCut objects 
#  @see ROOT::TCut
#  cut       = ...
#  other_cut = ...
#  cut      /= other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_idiv_ ( self , other ) :
    """ Multiplication for TCut objects 
    >>> cut    = ...
    >>> other  = ... 
    >>> cut   *= other 
    """
    self.strip() 
    if       isinstance  ( other , string_types ) : other = ROOT.TCut ( other )
    elif     isinstance  ( other , num_types    ) :
        ## 
        if not other  : raise ZeroDivisionError ( "ROOT.TCut(%s)/=%s  is illegal!" % ( self , other ) )
        ## 
        if 1 == other : return self
        ## 
        if   isinstance  ( other , integer_types ) : other = ROOT.TCut ( fmt_int   % other )
        elif isinstance  ( other , float         ) :
            if not isfinite ( other ) : return NotImplemented                         
            other = ROOT.TCut ( fmt_float % other ) 
        else                                       : other = ROOT.TCut ( fmt_other % other )  
    ##  
    if not isinstance ( other , ROOT.TCut    ) : return NotImplemented
    ## 
    if   other and self  :
        self.SetTitle("%s/%s" % ( self , other ) )
        return self
    elif other : return NotImplemented 
    ## 
    return self 

# =============================================================================
## Subtraction for TCut objects 
#  @see ROOT::TCut
#  cut       = ...
#  other_cut = ...
#  cut      -= other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_isub_ ( self , other ) :
    """ Subtraction for TCut objects 
    >>> cut    = ...
    >>> other  = ... 
    >>> cut   -= other 
    """
    self.strip() 
    if       isinstance  ( other , string_types  ) : other = ROOT.TCut ( other ) 
    elif     isinstance  ( other , num_types     ) :
        ##
        if not other : return self
        ## 
        if   isinstance  ( other , integer_types ) : other = ROOT.TCut ( fmt_int   % other )
        elif isinstance  ( other , float         ) :
            if not isfinite ( other ) : return NotImplemented                                                 
            other = ROOT.TCut ( fmt_float % other ) 
        else                                       : other = ROOT.TCut ( fmt_other % other )  
    ##        
    if not isinstance ( other , ROOT.TCut    ) : return NotImplemented
    ## 
    if   other and self  :
        self.SetTitle("%s-%s" % ( self , other ) )
        return self
    elif other : return NotImplemented 
    ## 
    return self 

# =============================================================================
## Module for TCut objects 
#  @see ROOT::TCut
#  cut       = ...
#  other_cut = ...
#  cut      %= other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_imod_ ( self , other ) :
    """ Multiplication for TCut objects 
    >>> cut    = ...
    >>> other  = ... 
    >>> cut   *= other 
    """
    self.strip()
    ## 
    if       isinstance ( other , string_types  ) : other = ROOT.TCut ( other )
    elif     isinstance ( other , num_types     ) :
        ## 
        if not other  : raise ZeroDivisionError ( "ROOT.TCut(%s)%=%s  is illegal!" % ( self , other ) )
        ## 
        if   isinstance  ( other , integer_types ) : other = ROOT.TCut ( fmt_int   % other )
        elif isinstance  ( other , float         ) :
            if not isfinite ( other ) : return NotImplemented                                                             
            other = ROOT.TCut ( fmt_float % other ) 
        else                                       : other = ROOT.TCut ( fmt_other % other )  
    ##  
    if not isinstance ( other , ROOT.TCut ) : return NotImplemented
    ## 
    if   other and self  :
        self.SetTitle("%s%%%s" % ( self , other ) )
        return self
    elif other : return NotImplemented 
    ## 
    return self 

# =============================================================================
## abs function for TCut object
#  @see ROOT:TCut
#  @code
#  cut = ...
#  new_cut = abs ( cut )
#  @endcode
def _tc_abs_ ( self ) :
    """ abs function for TCut object
    >>> cut = ...
    >>> new_cut = abs ( cut )
    """
    self.strip()
    if not self : return NotImplemented
    ##
    new_cut = ROOT.TCut ( self )
    new_cut.SetTitle ( 'abs(%s)' % self )
    return new_cut 

atypes = string_types +             ( ROOT.TCut , bool )
btypes = string_types + num_types + ( ROOT.TCut ,      )
# =============================================================================
## Logical *AND* for TCut objects  
#  @see ROOT::TCut
#  cut       =
#  other_cut = 
#  new_cut   = cut & other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_and_ ( self , other ) :
    """ Logical *AND* for TCut objects 
    >>> cut        = ...
    >>> other_cut  = ... 
    >>> new_cut    = cut & other_cut 
    """
    if not isinstance ( other , atypes ) : return NotImplemented 
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
    """ Logical *OR* for TCut objects 
    >>> cut        = ...
    >>> other_cut  = ... 
    >>> new_cut    = cut | other_cut 
    """
    if not isinstance ( other , atypes ) : return NotImplemented 
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
    if not isinstance ( other , btypes ) : return NotImplemented 
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
    """ Division for TCut objects 
    >>> cut        = ...
    >>> other_cut  = ... 
    >>> new_cut    = cut / other_cut 
    """
    ##
    if not isinstance ( other , btypes ) : return NotImplemented
    new_cut  = ROOT.TCut ( self )
    new_cut /= other 
    return new_cut

# =============================================================================
## Addition for TCut objects 
#  @see ROOT::TCut
#  cut       =
#  other_cut = 
#  new_cut   = cut + other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_add_ ( self , other ) :
    """ Addtion for TCut objects 
    >>> cut        = ...
    >>> other_cut  = ... 
    >>> new_cut    = cut + other_cut 
    """
    ##
    if not isinstance ( other , btypes ) : return NotImplemented
    new_cut  = ROOT.TCut ( self )
    new_cut += other 
    return new_cut

# =============================================================================
## Subtraction for TCut objects 
#  @see ROOT::TCut
#  cut       =
#  other_cut = 
#  new_cut   = cut - other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_sub_ ( self , other ) :
    """ Subtraction for TCut objects 
    >>> cut        = ...
    >>> other_cut  = ... 
    >>> new_cut    = cut - other_cut 
    """
    ##
    if not isinstance ( other , btypes ) : return NotImplemented
    new_cut  = ROOT.TCut ( self )
    new_cut -= other 
    return new_cut

# =============================================================================
## Mod=-operation for TCut objects 
#  @see ROOT::TCut
#  cut       =
#  other_cut = 
#  new_cut   = cut % other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_mod_ ( self , other ) :
    """ Mod- operation for TCut objects 
    >>> cut        = ...
    >>> other_cut  = ... 
    >>> new_cut    = cut % other_cut 
    """
    ##
    if not isinstance ( other , btypes ) : return NotImplemented
    new_cut  = ROOT.TCut ( self )
    new_cut %= other 
    return new_cut

# =============================================================================
## Exponentiation of certain expression## Exponentiation of certain expression
#  @code
#  cut       =     
#  other_cut = 
#  new_cut   = cut ** other_cut
#  @endcode
#  @see ROOT::TCut 
#
def _tc_pow_ ( self, other ) : 
    """ Exponentiation of certain expression
    """   
    if not isinstance ( other , btypes ) : return NotImplemented
    
    if   isinstance   ( other , string_types ) : other = ROOT.TCut ( other )
    elif isinstance   ( other , num_types    ) :
        ## 
        if   not other  : return ROOT.TCut ( '1'  )   
        elif 1 == other : return ROOT.TCut ( self )
        ## 
        if   isinstance  ( other , integer_types ) : other = ROOT.TCut ( fmt_int   % other )
        elif isinstance  ( other , float         ) :
            if not isfinite ( other ) : return NotImplemented                                                                         
            other = ROOT.TCut ( fmt_float % other ) 
        else                                       : other = ROOT.TCut ( fmt_other % other )  
    ## 
    if not isinstance ( other , ROOT.TCut    ) : return NotImplemented
    ## 
    if   not self  : return NotImplemented    
    elif not other : return NotImplemented  
    ##
    new_cut  = ROOT.TCut ( self )
    new_cut.SetTitle("%s**%s" % ( new_cut, other ) )
    ##
    return new_cut
    
# =============================================================================
## minor extension for TCut
#  @see ROOT::TCut
#  cut       = ...
#  new_cut   = 'pt>1' & other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_rand_ ( self , other ) :
    """ Logical *AND* for TCut objects 
    >>> cut        = ...
    >>> new_cut    = 'pt>1' & other_cut 
    """
    if not isinstance ( other , atypes ) : return NotImplemented
    return self & other 

# =============================================================================
## minor extension for TCut
#  @see ROOT::TCut
#  cut       = ...
#  new_cut   = 'pt>1' | other_cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_ror_ ( self , other ) :
    """ Logical *AND* for TCut objects 
    >>> cut        = ...
    >>> new_cut    = 'pt>1' | other_cut 
    """
    if not isinstance ( other , atypes ) : return NotImplemented
    return self | other 

# =============================================================================
## Multiplication for TCut objects 
#  @see ROOT::TCut
#  cut       =
#  other_cut = 
#  new_cut   = other_cut   *  cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_rmul_ ( self , other ) :
    """ Multiplication for TCut objects 
    >>> cut        = ...
    >>> other_cut  = ... 
    >>> new_cut    = other_cut * cut
    """
    if not isinstance ( other , btypes ) : return NotImplemented
    return self * other 

# =============================================================================
## Addition for TCut objects 
#  @see ROOT::TCut
#  cut       =
#  other_cut = 
#  new_cut   = other_cut + cut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_radd_ ( self , other ) :
    """ Addtition for TCut objects 
    >>> cut        = ...
    >>> other_cut  = ... 
    >>> new_cut    = other_cut + cut
    """
    if not isinstance ( other , btypes ) : return NotImplemented
    return self + other 

# =============================================================================
## Subtraction for TCut objects 
#  @see ROOT::TCut
#  cut       =
#  other_cut = 
#  new_cut   = other_cut -  cut 
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_rsub_ ( self , other ) :
    """ Subtraction for TCut objects 
    >>> cut        = ...
    >>> other_cut  = ... 
    >>> new_cut    = other_cut - cut
    """
    ##
    self.strip()
    ## 
    if   not self : return NotImplemented
    ## 
    elif isinstance ( other , ROOT.TCut    ) : return             other   - self  
    elif isinstance ( other , string_types ) : return ROOT.TCut ( other ) - self 
    elif isinstance ( other , num_types    ) : return ROOT.TCut ( '%s-%s' % ( other , self ) )
    ## 
    return NotImplemented
    
# =============================================================================
## Division for TCut objects 
#  @see ROOT::TCut
#  cut       =
#  other_cut = 
#  new_cut   = other_cut /  div 
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_rdiv_ ( self , other ) :
    """ Right division for TCut objects 
    >>> cut        = ...
    >>> other_cut  = ... 
    >>> new_cut    = other_cut / cut
    """
    ##
    self.strip ()
    ## 
    if   not self : return NotImplemented    
    elif isinstance ( other , ROOT.TCut    ) : return             other   / self  
    elif isinstance ( other , string_types ) : 
        if not other : return NotImplemented 
        return ROOT.TCut ( other ) / self 
    elif isinstance ( other , num_types    ) : 
        if not other : return ROOT.TCut ( '0' )
        if   isinstance  ( other , integer_types ) : return ROOT.TCut ( ( fmt_int % other ) + (  '/%s' % self ) )
        elif isinstance  ( other , float         ) :
            if not isfinite ( other ) : return NotImplemented 
            return ROOT.TCut ( ( fmt_float % other ) + (  '/%s' % self ) )
    ## 
    return NotImplemented

# =============================================================================
## Exponentiation of certain expression## Exponentiation of certain expression
#  @code
#  cut       =     
#  other_cut = 
#  new_cut   = cut ** other_cut
#  @endcode
#  @see ROOT::TCut 
#
def _tc_rpow_ ( self, other ) : 
    """ Exponentiation of certain expression
    """   
    if not isinstance ( other , btypes ) : return NotImplemented
    
    if   isinstance      ( other , string_types  ) : other = ROOT.TCut ( other ) 
    elif isinstance      ( other , num_types     ) :
        ## 
        if   1 == other  : return ROOT.TCut ( '1' )
        ## 
        elif isinstance  ( other , integer_types ) : other = ROOT.TCut ( fmt_int   % other )
        elif isinstance  ( other , float         ) :
            if not isfinite ( other ) :            return NotImplemented 
            other = ROOT.TCut ( fmt_float % other ) 
        else                                       : other = ROOT.TCut ( fmt_other % other )  
    ## 
    if   not isinstance ( other , ROOT.TCut    ) : return NotImplemented
    ## 
    if   not self  : return NotImplemented    
    elif not other : return NotImplemented  
    ##
    return other ** self

# =============================================================================
## (right) mod operator 
#  @code
#  cut       =     
#  other_cut = 
#  new_cut   = cut % other_cut
#  @endcode
#  @see ROOT::TCut 
#
def _tc_rmod_ ( self, other ) : 
    """ (right) mod operator 
    """   
    if not isinstance ( other , btypes ) :         return NotImplemented
    ## 
    if isinstance      ( other , num_types     ) :
        ##
        if   isinstance  ( other , integer_types ) : other = ROOT.TCut ( fmt_int   % other )
        elif isinstance  ( other , float         ) :
            if not isfinite ( other ) :            return NotImplemented 
            other = ROOT.TCut ( fmt_float % other ) 
        else                                       : other = ROOT.TCut ( fmt_other % other )  
    ## 
    if   not isinstance ( other , ROOT.TCut    ) : return NotImplemented
    ## 
    if   not self  : return NotImplemented    
    elif not other : return NotImplemented  
    ##
    return other % a 

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
    """ Invert the cut 
    >>> cut     =
    >>> new_cut = ~cut 
    """
    nc  = ROOT.TCut ( c .strip() )
    tit = nc.GetTitle() 
    if tit : return ROOT.TCut( '!(%s)' % tit )
    ##
    return ROOT.TCut

# ===============================================================================
## Less for TCut object(s)
#  @code 
#  cut     =
#  new_cut = cut < 9 
#  @endcode 
#  @see ROOT::TCut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_lt_ ( cut , other ) :
    """ Less for TCut object(s)
    >>> cut     =
    >>> new_cut = cut < 9 
    """
    if not cut : return NotImplemented
    ##
    if   isinstance  ( other , integer_types ) : other = ROOT.TCut ( fmt_int   % other )
    elif isinstance  ( other , float         ) : other = ROOT.TCut ( fmt_float % other )
    elif isinstance  ( other , string_types  ) : other = ROOT.TCut ( other )
    ##
    if isinstance    ( other , ROOT.TCut    ) :
        return ROOT.TCut ( '%s<%s' % ( cut , other ) ) if other else NotImplemented 
    ## 
    return NotImplemented 

# ===============================================================================
## Less-or-equal for TCut object(s)
#  @code 
#  cut     =
#  new_cut = cut <= 9 
#  @endcode 
#  @see ROOT::TCut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_le_ ( cut , other ) :
    """ Less-or-equal for TCut object(s)
    >>> cut     =
    >>> new_cut = cut < 9 
    """
    if not cut : return NotImplemented
    ## 
    if   isinstance  ( other , integer_types ) : other = ROOT.TCut ( fmt_int   % other )
    elif isinstance  ( other , float         ) : other = ROOT.TCut ( fmt_float % other )
    elif isinstance  ( other , string_types  ) : other = ROOT.TCut ( other )
    ##
    if isinstance    ( other , ROOT.TCut    ) :
        return ROOT.TCut ( '%s<=%s' % ( cut , other ) ) if other else NotImplemented 
    ## 
    return NotImplemented 

# ===============================================================================
## Greater for TCut object(s)
#  @code 
#  cut     =
#  new_cut = cut > 9 
#  @endcode 
#  @see ROOT::TCut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_gt_ ( cut , other ) :
    """ Greater for TCut object(s)
    >>> cut     =
    >>> new_cut = cut < 9 
    """
    if not cut : return NotImplemented
    ## 
    if   isinstance  ( other , integer_types ) : other = ROOT.TCut ( fmt_int   % other )
    elif isinstance  ( other , float         ) : other = ROOT.TCut ( fmt_float % other )
    elif isinstance  ( other , string_types  ) : other = ROOT.TCut ( other )
    ##
    if isinstance    ( other , ROOT.TCut    ) :
        return ROOT.TCut ( '%s>%s' % ( cut , other ) ) if other else NotImplemented 
    ## 
    return NotImplemented 

# ===============================================================================
## Greater-or-equal for TCut object(s)
#  @code 
#  cut     =
#  new_cut = cut >= 9 
#  @endcode 
#  @see ROOT::TCut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_ge_ ( cut , other ) :
    """ Greater-or-equal for TCut object(s)
    >>> cut     =
    >>> new_cut = cut >= 9 
    """
    if not cut : return NotImplemented
    ## 
    if   isinstance  ( other , integer_types ) : other = ROOT.TCut ( fmt_int   % other )
    elif isinstance  ( other , float         ) : other = ROOT.TCut ( fmt_float % other )
    elif isinstance  ( other , string_types  ) : other = ROOT.TCut ( other )
    ##
    if isinstance    ( other , ROOT.TCut    ) :
        return ROOT.TCut ( '%s>=%s' % ( cut , other ) ) if other else NotImplemented 
    ## 
    return NotImplemented 

# ===============================================================================
## Equality for TCut object(s)
#  @code 
#  cut     =
#  new_cut = cut == 9 
#  @endcode 
#  @see ROOT::TCut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_eq_ ( cut , other ) :
    """ Greater-or-equal for TCut object(s)
    >>> cut     =
    >>> new_cut = cut == 9 
    """
    if not cut : return NotImplemented
    ## 
    if   isinstance  ( other , integer_types ) : other = ROOT.TCut ( fmt_int  % other )
    elif isinstance  ( other , float         ) :
        if       isnan    ( other ) : return ROOT.TCut ( 'TMath::IsNaN(%s)'   % cut )
        elif not isfinite ( other ) : return ROOT.TCut ( '!TMath::Finite(%s)' % cut )
        other = ROOT.TCut ( fmt_float % other )        
    elif isinstance  ( other , string_types  ) : other = ROOT.TCut ( other )
    ##
    if isinstance    ( other , ROOT.TCut    ) and other :
        return ROOT.TCut ( '%s==%s' % ( cut , other ) ) 
    ## 
    return NotImplemented 

# ===============================================================================
## Non-equality for TCut object(s)
#  @code 
#  cut     =
#  new_cut = cut != 9 
#  @endcode 
#  @see ROOT::TCut
#  @author Vanya BELYAEV Ivan.Belyaev
#  @date   2014-08-31
def _tc_ne_ ( cut , other ) :
    """ Non-equalily for TCut object(s)
    >>> cut     =
    >>> new_cut = cut != 9 
    """
    if not cut : return NotImplemented
    ## 
    if   isinstance  ( other , integer_types ) : other = ROOT.TCut ( fmt_int % other )
    elif isinstance  ( other , float         ) :
        if       isnan  ( other ) : return ROOT.TCut ( '!TMath::IsNaN(%s)'   % cut )
        if not isfinite ( other ) : return ROOT.TCut ( 'TMath::Finite(%s)'   % cut )        
        other = ROOT.TCut ( fmt_float % other )
    elif isinstance  ( other , string_types  ) : other = ROOT.TCut ( other )
    ##
    if isinstance    ( other , ROOT.TCut ) and other :
        return ROOT.TCut ( '%s!=%s' % ( cut , other ) ) 
    ## 
    return NotImplemented 

ROOT.TCut. __and__      = _tc_and_
ROOT.TCut.  __or__      = _tc_or_

ROOT.TCut.__rand__      = _tc_rand_
ROOT.TCut. __ror__      = _tc_ror_

ROOT.TCut. __mul__      = _tc_mul_
ROOT.TCut. __add__      = _tc_add_
ROOT.TCut. __sub__      = _tc_sub_
ROOT.TCut. __div__      = _tc_div_
ROOT.TCut. __truediv__  = _tc_div_
ROOT.TCut. __pow__      = _tc_pow_
ROOT.TCut. __mod__      = _tc_mod_

ROOT.TCut. __ior__      = _tc_ior_
ROOT.TCut. __iand__     = _tc_iand_

ROOT.TCut.__imul__      = _tc_imul_
ROOT.TCut.__iadd__      = _tc_iadd_
ROOT.TCut.__isub__      = _tc_isub_
ROOT.TCut.__idiv__      = _tc_idiv_
ROOT.TCut.__itruediv__  = _tc_idiv_
ROOT.TCut.__imod__      = _tc_imod_

ROOT.TCut.__radd__      = _tc_radd_
ROOT.TCut.__rmul__      = _tc_rmul_
ROOT.TCut.__rdiv__      = _tc_rdiv_
ROOT.TCut.__rsub__      = _tc_rsub_
ROOT.TCut.__rtruediv__  = _tc_rdiv_
ROOT.TCut.__rpow__      = _tc_rpow_
ROOT.TCut.__rmod__      = _tc_rmod_

ROOT.TCut.__invert__    = _tc_invert_
ROOT.TCut.__abs__       = _tc_abs_

ROOT.TCut. __lt__       = _tc_lt_
ROOT.TCut. __le__       = _tc_le_
ROOT.TCut. __gt__       = _tc_gt_
ROOT.TCut. __ge__       = _tc_ge_
ROOT.TCut. __eq__       = _tc_eq_
ROOT.TCut. __ne__       = _tc_ne_

# =============================================================================
_decorated_classes_ = (
    ROOT.TCut ,
    )
_new_methods_       = (
    #
    ROOT.TCut .  __and__  ,
    ROOT.TCut . __iand__  ,
    ROOT.TCut . __rand__  ,
    #
    ROOT.TCut .  __or__   ,
    ROOT.TCut . __ior__   ,
    ROOT.TCut . __ror__   ,
    #
    ROOT.TCut .  __add__  ,
    ROOT.TCut . __iadd__  ,
    ROOT.TCut . __radd__  ,
    #
    ROOT.TCut .  __sub__  ,
    ROOT.TCut . __isub__  ,
    ROOT.TCut . __rsub__  ,
    #
    ROOT.TCut .  __mul__  ,
    ROOT.TCut . __imul__  ,
    ROOT.TCut . __rmul__  ,
    #
    ROOT.TCut .  __div__  ,
    ROOT.TCut . __idiv__  ,
    ROOT.TCut . __rdiv__  ,
    #
    ROOT.TCut .  __truediv__  ,
    ROOT.TCut . __itruediv__  ,
    ROOT.TCut . __rtruediv__  ,
    #
    ROOT.TCut. __pow__     , 
    ROOT.TCut. __rpow__    ,
    #
    ROOT.TCut. __mod__     , 
    ROOT.TCut. __rmod__    , 
    ROOT.TCut. __imod__    , 
    #
    ROOT.TCut . __abs__    ,
    ROOT.TCut . __invert__ ,
    #
    ROOT.TCut. __lt__      , 
    ROOT.TCut. __le__      , 
    ROOT.TCut. __gt__      , 
    ROOT.TCut. __ge__      , 
    ROOT.TCut. __eq__      , 
    ROOT.TCut. __ne__      , 
    #
    ROOT.TCut . strip      ,
    ROOT.TCut . replace    ,  
    #  
    ROOT.TCut . ast        ,    
    )
# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
