#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with some simple functions for strings 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
# =============================================================================
""" Module with some simple functions for strings 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    ##
    'var_separators'       , ## separators form the split  
    'split_string_respect' , ## split the string  according to separators 
    'split_string'         , ## split the string  according to separators
    ##
    'math_symbols'         , ## symbols for math-operations
    'has_symbol'           , ## Has any of the symbols?
    'is_formula'           , ## Is this string expression represent math formula?    
    # =========================================================================
) # ===========================================================================
# =============================================================================
## defalt separators for the string expressions
var_separators = ',:;'
# =============================================================================
## defalt separators for the string expressions
var_separators = ',:;'
## rx_separators  = re.compile ( r'[,:;]\s*(?![^()]*\))' )
## rx_separators  = re.compile ( '[ ,:;](?!(?:[^(]*\([^)]*\))*[^()]*\))')
# =============================================================================
## mark for double columns: double column is a special for C++ namespaces 
dc_mark = '_SSPPLLIITT_'
# =============================================================================
## split string using separators and respecting the (),[] and {} groups.
#  - group can be nested
def split_string_respect  ( text , separators = var_separators , strip = True ) :
    """ Split string using separators and respecting the (),[] and {} groups.
    - groups can be nested
    """
    assert isinstance ( text , str ) , "Invalid `text` %s" % type ( text ) 
    protected = False 
    if ':' in separators and '::' in text :
        text      = text.replace ( '::' , dc_mark ) 
        protected = True 
    
    flag1  = 0
    flag2  = 0
    flag3  = 0
    item   = ''
    items  = []
    for c in text:
        if   c == '(' : flag1 += 1
        elif c == ')' : flag1 -= 1
        elif c == '[' : flag2 += 1
        elif c == ']' : flag2 -= 1
        elif c == '{' : flag3 += 1
        elif c == '}' : flag3 -= 1
        elif 0 == flag1 and 0 == flag2 and 0 == flag3 and c in separators :
            items .append ( item )
            item = ''
            continue
        item += c
        
    if item : items.append ( item  )

    if protected :
        nlst = []
        for item in items :
            if dc_mark in item : nlst.append ( item.replace ( dc_mark , '::' ) )
            else               : nlst.appenf ( item ) 
        items = nlst 
                              
    ## strip items if required 
    if strip : items = [ item.strip() for item in items ] 
    ## remove empty items
    return tuple ( item for item in items if item  )

# =============================================================================
## split string using separators:
#  @code
#  split_string ( ' a b cde,fg;jq', ',;:' )
#  @endcode
def split_string ( line                            ,
                   separators     = var_separators ,
                   strip          = False          ,
                   respect_groups = False          ) :
    """ Split the string using separators
    >>> split_string ( ' a b cde,fg;jq', ',;:' )
    """
    assert isinstance ( line , str ) , "Invalid `line` %s" % type ( line )
    
    if respect_groups :
        return split_string_respect ( line                    ,
                                      separators = separators ,
                                      strip      = strip      )
    ##
    protected = False 
    if ':' in separators and '::' in line :
        line      = line.replace ( '::' , dc_mark ) 
        protected = True 
    
    items = [ line ]
    for s in separators :
        result = []
        for item in items :
            if s in item : result += item.split ( s )
            else         : result.append ( item ) 
        items = result

    if protected :
        nlst = []
        for item in items :
            if dc_mark in item : nlst.append ( item.replace ( dc_mark , '::' ) )
            else               : nlst.appenf ( item ) 
        items = nlst 
        
    ## strip items if required 
    if strip : items = [ i.strip() for i in items ] 
    ## remove empty items 
    return tuple ( item for item in items if item )

# =========================================================================
## Has any of the symbols?
def has_symbol ( expression , symbols ) :
    """ Has any of the symbols?"""
    assert isinstance ( expression , str ) , "Invalid `expression` %s" % type ( expression  )
    return any ( s in expression for s in symbols )

# =========================================================================
## markers of the math formulas 
math_symbols = ' +-*/=><()[]^%&|'
## Does this string expression represent math formula?
def is_formula ( expression , symbols = math_symbols ) :
    """ Does  this string expression represent math formula?
    """
    assert isinstance ( expression , str ) , "Invalid `expression` %s" % type ( expression  )
    return has_symbol ( expression.strip() , symbols ) 

# =============================================================================
if '__main__' == __name__ :

    from ostap.logger.logger import getLogger
    logger = getLogger ( 'ostap.utils.strings' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
