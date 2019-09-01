#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Format the list of rows as a  table.
#  When <code>terminaltables</code> module is accessible,
#  use it to get nice pretty tables, otherwise use some home-made primitives
#  @see https://pypi.org/project/terminaltables
#  @see https://github.com/Robpol86/terminaltables
# 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2019-08-30
# =============================================================================
""" Format the list of rows as a table.
When terminaltables module is accessible,
use it to get nice pretty tables, otherwise use some home-made primitives
- see https://pypi.org/project/terminaltables
- see https://github.com/Robpol86/terminaltables
""" 
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2019-08-31"
__all__     = (
    'table'       , ## format the list of rows as a  table.
    'the_table'   , ## format the list of rows as a  table (local version)
    'table_width' , ## true  width of the table  
    )
# =============================================================================
from ostap.logger.colorized import infostr, allright, decolorize        
# =============================================================================
try :
    import terminaltables
except ImportError :
    terminaltables = None
# =============================================================================
## Format the list of rows as a  table (home-made primitive) 
#  - Each row is a sequence of column cells.
#  - The first row defines the column headers.
#  @code
#  table_data = [
#     ( 'Name'  , 'Occupation' , 'Note' ) ,
#     ( 'Alice' , '?'          , '---'  ) ,
#     ( 'Bob'   , 'unemployed' , ''     ) ,
#     ]
#  t = the_table ( table_data , 'Title' )
#  print (t)
#  @endcode
def the_table ( rows , title = '' ) :
    """Format the list of rows as a  table (home-made primitive) 
    - Each row is a sequence of column cells.
    - The first row defines the column headers.
    >>> table_data = [
    ...   ( 'Name'  , 'Occupation' , 'Note' ) ,
    ...   ( 'Alice' , '?'          , '---'  ) ,
    ...   ( 'Bob'   , 'unemployed' , ''     ) ]
    >>> t = the_table ( table_data , 'Title' )
    >>> print (t)
    """
    ## calculate the number of columns 
    nc  = 0
    for row in rows : nc = max ( nc , len ( row ) )
        
    ## calculate the maximum width for columns 
    widths = {}
    for row in rows :
        cols = [ c for  c in row ]
        while len  ( cols ) < nc : cols.append ('')             
        for i , c in enumerate ( cols ) :
            if not i in widths : widths[i] = 1
            widths[i] = max ( widths[i] , len ( decolorize ( c ) ) )

    totwidth = 0
    for c in widths : totwidth += widths[c]
    totwidth += ( nc - 1 ) * 3
    if totwidth < len ( title )  :
        delta = 1 + ( len ( title ) - totwidth ) // nc
        for c in widths : widths[c] += delta 

    hformats = [  "{:^%d}"  % widths[c] for c in range ( nc ) ]
    rformats = [ " {:^%d} " % widths[c] for c in range ( nc ) ]

    # the first column is left  adjusted
    rformats[ 0] = rformats[ 0].replace ( '^' , '<' ) 
    # the last column  is right adjusted
    rformats[-1] = rformats[-1].replace ( '^' , '>' ) 

    seps     = [   '-' * ( widths[c] + 2 ) for c in range ( nc ) ]
    sepline  = '+' + "+".join (  seps ) +  '+'


    table    = []
    if title :
        sline = '+' + '-' * ( len ( sepline ) - 2 ) + '+'
        tfmt  = "{:^%d}"  % ( len ( sepline ) - 4 )
        t     = '| ' + allright ( tfmt.format ( decolorize ( title ) ) ) + ' |'
        table.append ( sline  )
        table.append ( t      )

    table.append ( sepline )    
    for i , row in enumerate ( rows ) :
        
        cols = [ c for  c in row ]
        while len  ( cols ) < nc : cols.append ('')             

        if 0 == i :
            table.append ( '| ' + ' | '.join ( [ infostr ( f.format ( decolorize ( i ) ) )  for ( f , i ) in zip ( hformats , cols ) ] ) + ' |' ) 
            table.append ( sepline )             
        else :
            table.append ( '|'  +  '|' .join ( [           f.format ( decolorize ( i ) )    for ( f , i ) in zip ( rformats , cols ) ] ) +  '|' ) 


    table.append ( sepline )

    return '\n'.join ( table )

# =============================================================================
## Format the list of rows as a  table.
#  - Each row is a sequence of column cells.
#  - The first row defines the column headers.
# 
#  When <code>terminaltables</code> module is accessible,
#  use it to get nice pretty tables, otherwise use some
#  home-made primitive replacement 
#  @see https://pypi.org/project/terminaltables
#  @see https://github.com/Robpol86/terminaltables
#  @code
#  table_data = [
#     ( 'Name'  , 'Occupation' , 'Note' ) ,
#     ( 'Alice' , '?'          , '---'  ) ,
#     ( 'Bob'   , 'unemployed' , ''     ) ,
#     ]
#  t = table ( table_data , 'Title' )
#  print (t)
#  @endcode
def table ( rows , title = '' ) :
    """Format the list of rows as a  table.
    - Each row is a sequence of column cells.
    - The first row defines the column headers.

    When terminaltables package is accessible,
    use it to get nice pretty tables, otherwise use some
    home-made primitive replacement 
    - see https://pypi.org/project/terminaltables
    - see https://github.com/Robpol86/terminaltables

    >>> table_data = [
    ...   ( 'Name'  , 'Occupation' , 'Note' ) ,
    ...   ( 'Alice' , '?'          , '---'  ) ,
    ...   ( 'Bob'   , 'unemployed' , ''     ) ]
    >>> t = table ( table_data , 'Title' )
    >>> print (t)
    """
    
    from ostap.utils.basic import isatty

    title = allright ( decolorize ( title ) )
    if rows :
        rows       = list  ( rows ) 
        header_row = rows[0]
        header_row = [ infostr ( decolorize ( c ) ) for c in header_row ]
        rows [0]   = header_row 
        rows = tuple ( rows )
        
    if terminaltables and isatty () :
        
        table_instance = terminaltables.SingleTable ( rows , title)
        table_instance.justify_columns[ 0] = 'left'
        table_instance.justify_columns[ 2] = 'right'
        return table_instance.table
    
    elif terminaltables :
        
         title = allright ( title ) 
         table_instance = terminaltables.AsciiTable ( rows , title)
         table_instance.justify_columns[ 0] = 'left'
         table_instance.justify_columns[ 2] = 'right'
         return table_instance.table

     ## use the local replacement 
    return the_table ( rows , title ) 

# =============================================================================
## get the true  table width 
def table_width  ( table ) :
    """Get the true width of the table
    """
    width = 0
    for row in table.split('\n') :
        width = max ( width , len ( decolorize ( row ) ) )
    return width 

# =============================================================================
if __name__ == '__main__' :

    from ostap.logger.logger import getLogger
    logger = getLogger ( 'ostap.logger.table')

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    if terminaltables :
        logger.info     ( "``terminaltables'' will be used for formatting" )
    else :
        logger.warning  ( "``terminaltables'' is not available, use local replacement" )
        
    table_data = [
        ( 'Name'  , 'Occupation' , 'Note' ) ,
        ( 'Alice' , '?'          , '---'  ) ,
        ( 'Bob'   , 'unemployed' , ''     ) ]
    
    logger.info ( 'The table is \n' + table ( table_data , 'Title' ) ) 

    
    
# =============================================================================
##                                                                      The END 
# =============================================================================
