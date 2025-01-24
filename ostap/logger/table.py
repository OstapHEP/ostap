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
    'table'                 , ## format the list of rows as a  table.
    'the_table'             , ## format the list of rows as a  table (local version)
    'table_width'           , ## true  width of the table
    'align_column'          , ## align the certain column of the table  
    'add_prefix'            , ## add the prefix to each row of the table
    'empty_columns'         , ## find empty columns in the table 
    'remove_empty_columns'  , ## remove empty columns from the table 
    )
# =============================================================================
from   ostap.core.ostap_types import string_types 
from   ostap.logger.colorized import infostr, allright, decolorize        
from   ostap.utils.basic      import terminal_size, zip_longest 
import textwrap, os, sys
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.logger.table' )
else                       : logger = getLogger( __name__             )
# =============================================================================
terminaltables = None
# =============================================================================
if ( 3 , 9 ) <= sys.version_info : # ==========================================
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        import terminaltables3 as terminaltables
        # =====================================================================
    except ImportError : # ====================================================
        # =====================================================================
        terminaltables = None
# =============================================================================
if not terminaltables : # =====================================================
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        import terminaltables
        # =====================================================================
    except ImportError : # ====================================================
        # =====================================================================
        terminaltables = None 

# =============================================================================
visible_width = None
# =============================================================================
if terminaltables and not visible_width : # ===================================
    # =========================================================================
    visible_width = terminaltables.width_and_alignment.visible_width
    # =========================================================================
# =============================================================================
if not visible_width : # ======================================================
    # =========================================================================
    import unicodedata 
    ## visible lenght of the string/expression
    def visible_width ( what ) :
        """ Visible length of the string/expression (local, use unicodedata) 
        """
        if not what : return 0
        item = decolorize ( what )
        if not item : return 0
        ##
        if  sys.version_info < ( 3 , 0 ) : item = item.decode("u8")
        ## 
        for char in item :
            if unicodedata.east_asian_width ( char ) in ( 'F' , 'W' ) : width += 2
            else                                                      : width += 1
        return width 
    
# =============================================================================
if terminaltables : # =========================================================
    ## Vaild `style` arguments (case insensitive) 
    table_styles = ( 'local'     , ## use local, handmade replacement
                     'ascii'     , ## use `AsciiTable` 
                     'single'    , ## use `SingleTable` 
                     'porcelain' , ## use `PorcelainTable` 
                     'github'    , ## use `GithubFlavoredMarkdownTable`
                     'markdown'  , ## use `GithubFlavoredMarkdownTable`
                     'double'    , ## use `DoubleTable` 
                     'default'   ) ## use `DoubleTable`
    # =========================================================================
else : # ======================================================================
    # =========================================================================
    table_styles = 'local' , 
    # =========================================================================

    
print ( 'TERMINALTABLES' , terminaltables )
print ( 'VISIBLE_WIDTH'  , visible_width  )
raise TypeError('HERE!!!') 


    
# =============================================================================
## default style
default_style = 'default' if terminaltables else 'local'

# =============================================================================
## environment variable
env_var = 'OSTAP_TABLE_STYLE'
# =============================================================================
if    not terminaltables : # ==================================================
    # =========================================================================
    ## only 'local' is available for this case 
    default_style = 'local'
    # =========================================================================
else : # ======================================================================
    # =========================================================================
    from ostap.utils.basic import has_env as ostap_hasenv 
    if  ostap_hasenv ( env_var ) :
        # =====================================================================
        from ostap.utils.basic import get_env as ostap_getenv 
        default_style = ostap_getenv ( env_var, default_style )
        # =====================================================================
    else : # ==================================================================
        # =====================================================================
        ## get the preferred table style from the configuration file(s)
        import ostap.core.config as OCC
        if 'STYLE' in OCC.tables :
            default_style = OCC.tables.get ( 'STYLE' , fallback = default_style )
# =============================================================================
## finally adjust the style
default_style = default_style.strip().lower()
# =============================================================================
## very light wrapper for tabulate tables
#  @see https://github.com/astanin/python-tabulate
# =============================================================================
try : # =======================================================================
    # =========================================================================
    #  @see https://github.com/astanin/python-tabulate
    import tabulate
    # =========================================================================
    ## available formats 
    tabulate_styles = tuple ( tabulate.tabulate_formats ) 
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    ## no formats are available
    tabulate         = None
    tabulate_styles  = () 
    pass


# =============================================================================
left        = '<' , 'l' , 'left'  
right       = '>' , 'r' , 'right' 
center      = '^' , 'c' , 'center' , '='
wrapped     = 'w' , 'p' , 'wrap'   , 'wrapped'
max_width   = 50
wrap_indent = ' '
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
def the_table ( rows                          ,
                title           = ''          ,
                prefix          = ''          ,
                alignment       = ()          ,
                wrap_width      = -1          ,
                colorize_header = True        , 
                indent          = wrap_indent ,
                style           = ''          ,
                maxwidth        = -1          ) :
    """ Format the list of rows as a  table (home-made primitive) 
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

    title = allright ( decolorize ( title ).strip() ) 
    
    rows = list ( rows )
    
    nc   = 0
    for row in rows : nc = max ( nc , len ( row ) )

    wraps = [] 
    for i , a in  zip_longest ( range ( nc ) , alignment ) :
        if a and isinstance ( a , str ) :
            al = a.lower() 
            if   al in left   : pass 
            elif al in right  : pass 
            elif al in center : pass 
            elif al in wrapped : wraps.append ( i ) 

    ## calculate the maximum width for columns 
    widths = {}
    for k , row in enumerate ( rows ) :
        cols = [ c for  c in row ]
        while len  ( cols ) < nc : cols.append ('')             
        for i , c in enumerate ( cols ) :
            if not i in widths : widths[i] = 1
            widths [ i ] = max ( widths[i] , len ( decolorize ( c ) ) )
        cols = tuple ( cols )
        rows [ k ] = cols 

    ## play with wrapped columns
    while  wraps :

        twidth = 1 + len ( prefix )
        
        for k in widths :
            if not k in wraps :
                twidth += widths [ k ] + 2
                twidth += nc + 1

        ww = 12
        
        w , _ = terminal_size()
        if w <= twidth : break
        
        nw = len ( wraps ) 
        ww = ( w - twidth ) - 2 * nw
        ww , _ = divmod ( ww  , nw )


        if 12 < wrap_width and wrap_width < ww :
            ww = wrap_width
            
        if ww < 15 : break

        lw = len ( wraps )                 
        wraps = [ i for i in wraps if ww <= widths [i] ]
        
        if len ( wraps ) == lw  : break 

    for i in wraps :
        widths [ i ] = max ( min ( ww  , widths [ i ] ) , 10 ) 
            
    hformats = [  "{:^%d}"  % widths [ c ] for c in range ( nc ) ]
    rformats = [ " {:^%d} " % widths [ c ] for c in range ( nc ) ]

    for i , a in  zip ( range ( nc ) , alignment ) :
        if a and isinstance ( a , str ) :
            al = a.lower() 
            if   al in left  or al in wrapped :
                hformats [ i ] = hformats [ i ].replace ( '^' , '<' )
                rformats [ i ] = rformats [ i ].replace ( '^' , '<' )
            elif al in right :
                hformats [ i ] = hformats [ i ].replace ( '^' , '>' )
                rformats [ i ] = rformats [ i ].replace ( '^' , '>' )

    if wraps :
        rows_ = rows
        rows  = [] 
        for row in rows_ :
            cells = [] 
            for i , c in enumerate ( row ) :
                if i in wraps and wrap_width < len ( c ) : cells.append ( textwrap.wrap ( indent + c , widths [ i ]  ) )
                else          : cells.append ( [ c ] )
            nr = 0 
            for c in cells : nr = max ( nr , len ( c ) )
            
            for l in cells :
                while len  (l )  < nr :
                    l.insert ( 0 , '' )
                    l.append (     '' )
                                        
            for r in range ( nr ) :
                new_row = [] 
                for i , c in enumerate ( cells ) :
                    lc = len ( c )
                    if r < lc : new_row.append ( c[r] )
                    else      : new_row.append ( ''   )
                rows.append ( new_row )
                
    totwidth = 0
    for c in widths : totwidth += widths[c]
    totwidth += ( nc - 1 ) * 3
    titwidth  = len ( decolorize ( title ) )
    if totwidth < titwidth :
        delta = 1 + ( titwidth - totwidth ) // nc
        for c in widths : widths[c] += delta 

    seps     = [   '-' * ( widths [ c ]  + 2 ) for c in range ( nc ) ]
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
            if colorize_header : 
                table.append ( '| ' + ' | '.join ( [ infostr ( f.format ( decolorize ( i ) ) )  for ( f , i ) in zip ( hformats , cols ) ] ) + ' |' )
            else :
                table.append ( '|'  +  '|' .join ( [           f.format ( decolorize ( i ) )    for ( f , i ) in zip ( rformats , cols ) ] ) +  '|' )                
            table.append ( sepline )             
        else :
            table.append ( '|'  +  '|' .join ( [           f.format ( decolorize ( i ) )    for ( f , i ) in zip ( rformats , cols ) ] ) +  '|' ) 

    table.append ( sepline )
    
    return prefix + ( '\n' + prefix ) .join ( table ) if prefix else '\n'.join ( table ) 

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
#  Vaild `style` arguments (case insensitive) 
#   - local :  use local, handmade repalcement
#   - ascii :  use AsciiTable 
#   - single :  use SingleTable 
#   - porcelain :  use PorcelainTable 
#   - github :  use GithubFlavoredMarkdownTable
#   - markdown :  use GithubFlavoredMarkdownTable
#   - double (default) : DoubleTable 
def table ( rows                          ,
            title           = ''          ,
            prefix          = ''          ,
            alignment       = ()          ,
            wrap_width      = -1          ,
            colorize_header = True        ,
            indent          = wrap_indent ,            
            style           = None        ,
            maxwidth        = -1          , **kwargs ) :
    """ Format the list of rows as a  table.
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
    
    Valid `style` arguments (case insensitive) 
    - `local`     : use local, handmade repalcement
    - `ascii`     : use `AsciiTable` 
    - `single`    : use `SingleTable` 
    - `porcelain` : use `PorcelainTable` 
    - `github`    : use `GithubFlavoredMarkdownTable`
    - `markdown`  : use `GithubFlavoredMarkdownTable`
    - `double`    : use `DoubleTable` 
    - (default)   : use `DoubleTable` 
    """

    from ostap.utils.basic import isatty

    title_bw = decolorize  ( title    ).strip ()
    title    = allright    ( title_bw )    
    twidth   = visible_width ( title  ) ## visible lenght of the title 
    pwidth   = visible_width ( prefix ) ## visible lenght of the prefix 

    if maxwidth <= pwidth : maxwidth = terminal_size() [ 0 ]

    rows = [ list ( row ) for row in preprocess_table ( rows ) ]
    
    if rows and colorize_header : 
        header_row = rows [ 0 ]
        header_row = [ infostr ( decolorize ( c ) ) for c in header_row ]
        rows [ 0 ] = header_row 
        
    ## if style is None : style = '%s' % default_style
    if not style : style = '%s' % default_style
    
    fmt = style.lower()
    
    if not fmt in table_styles and not fmt in tabulate_styles :
        fmt = 'local' 

    if 'local' == fmt or not terminaltables :
        if kwargs :
            keys = ', '.join ( key for key in kwargs ) 
            logger.warning ( 'Ignore keyword arguments: %s' % keys )  
        ## use the local replacement 
        return the_table ( rows                              ,
                           title           = title           , 
                           prefix          = prefix          ,
                           alignment       = alignment       ,
                           wrap_width      = wrap_width      , 
                           colorize_header = colorize_header ,
                           indent          = indent          ,
                           style           = style           ,
                           maxwidth        = maxwidth        ) 
    
    elif terminaltables and fmt in table_styles    :  pass         
    elif tabulate       and fmt in tabulate_styles :
        colalign = [] 
        if alignment and not 'colalign' in kwargs :
            for a in alignment :
                aL = a.lower()
                if   aL in left   : colalign.append ( 'left'   )
                elif aL in right  : colalign.append ( 'right'  )
                elif aL in center : colalign.append ( 'center' )
                else              : colalign.append ( None     )
        colalign = tuple ( colalign )
        if colalign  : result = tabulate.tabulate ( rows , tablefmt = fmt , colalign = colalign , **kwargs )
        else         : result = tabulate.tabulate ( rows , tablefmt = fmt                       , **kwargs )
        ##
        if title :
            tab_width = table_width ( result )
            if fmt.startswith( 'latex' ) : title = '%% ' + title             
            title , twidth = shorten_title ( title , tab_width  )
            result = title + '\n' + result 
            
        if prefix : result = add_prefix ( result , prefix )
        ##
        return result

    if kwargs :
        keys = ', '.join ( key for key in kwargs ) 
        logger.warning ( 'Ignore keyword arguments: %s' % keys )  
        
    if 'ascii' == fmt or not isatty() : 
        table_instance = terminaltables.AsciiTable                  ( rows , title )        
    elif 'single' == fmt : 
        table_instance = terminaltables.SingleTable                 ( rows , title )
    elif 'porcelain' == fmt :
        table_instance = terminaltables.PorcelainTable              ( rows )
    elif fmt in ( 'github' , 'markdown' ) :
        table_instance = terminaltables.GithubFlavoredMarkdownTable ( rows )
    elif 'double' == fmt : 
        table_instance = terminaltables.DoubleTable                 ( rows , title )
    else :
        table_instance = terminaltables.DoubleTable                 ( rows , title )
        
    cw    = table_instance.column_widths
    nc    = len ( cw )
    wraps = [ i for ( i , a ) in enumerate ( alignment ) if a in wrapped ]

    wrap_width = 12 
    if wraps :
        
        max_dimensions = terminaltables.width_and_alignment.max_dimensions
        widths         = max_dimensions ( table_instance.table_data    ,
                                          table_instance.padding_left  ,
                                          table_instance.padding_right ) [2]
        widths  = sum ( l  for ( i , l ) in enumerate ( widths ) if not i in wraps ) 
        widths += nc + 1 + pwidth + 4 + 2 * len ( wraps )        
        ww , _ = divmod ( maxwidth  - widths  , len ( wraps ) )        
        wrap_width = max ( ww , 10  )

    nw = len ( wraps )        
    for i, a in zip ( range ( nc ) , alignment ) :
        if a and isinstance ( a , str ) :
            al = a.lower() 
            if   al in left   : table_instance.justify_columns [ i ] = 'left'
            elif al in right  : table_instance.justify_columns [ i ] = 'right'
            elif al in center : table_instance.justify_columns [ i ] = 'center'
            elif al in wrapped :
                maxw  = table_instance.column_max_width ( i )
                if 15 < wrap_width * nw < maxw : maxw = ( wrap_width - 3 ) * nw if 1 < nw else wrap_width 
                if maxw < 15 :                   maxw = ( wrap_width - 3 ) * nw if 1 < nw else wrap_width 
                if maxw < 15 :                   maxw = ( max_width  - 3 ) * nw if 1 < nw else max_width
                width = maxw / nw if 1 < nw else maxw
                width = wrap_width 
                for l , line in enumerate ( table_instance.table_data ) :
                    if i < len ( line ) and width < len ( line [ i ] ) : 
                        table_instance.table_data[l][i] = textwrap. fill ( line [ i ] , width = wrap_width , initial_indent = indent )
                        
    if title and ( 'github' , 'markdown' ) :
        title   = '# ' + title
        twidth += 2

    tabwidth = table_instance.table_width
    if title and twidth and tabwidth < twidth + 2 :
        ## is there a room for stretch ?
        if tabwidth + pwidth <= maxwidth :            
            newwidth = min ( twidth + 2 + pwidth , maxwidth  )
            if tabwidth + pwidth <= newwidth : 
                a  , _  = divmod ( newwidth - tabwidth - pwidth , nc )
                a1 , r  = divmod ( a        , 2 )
                table_instance.padding_left  += a1 + r 
                table_instance.padding_right += a1  
                tabwidth = table_instance.table_width
                
    if tabwidth < twidth + 2 :
        title , twidth = shorten_title ( title , tabwidth - 2 )
        table_instance.title = title 
        
    ## get the final table 
    table  = table_instance.table

    # special final adjustment 
    if title  and fmt in ( 'porcelain' , ) :
        line  = twidth * '-'
        table = title + '\n' + line + '\n' + table 
    elif title and fmt in ( 'github' , 'markdown' ) :
        line  = '---'
        table = line  + '\n' + title + '\n' + line + '\n' + table 

    result = add_prefix ( table , prefix )

    if sys.version_info < ( 3 , 0 ) :
        if isinstance ( result , unicode ) :
            result = result.encode ('utf-8')
            
    return result 
    
# =============================================================================
## Shorten the title
def shorten_title ( title , size ) :
    """ Shorten the title
    """
    twidth = visible_width ( title )
    size   = max ( size , 5 ) 
    if twidth <= size : return title , twidth
    title  = decolorize ( title )

    if   5 == size : title = '%s<.>%s'   % ( title[0] , title[-1] )
    elif 6 == size : title = '%s<.,>%s'  % ( title[0] , title[-1] )
    elif 7 == size : title = '%s<...>%s' % ( title[0] , title[-1] )
    else :
        a , r = divmod ( size - 5 , 2 )
        head = title [      : a ] 
        tail = title [ -a-r :   ] if r else title [ -a : ]
        title = '%s<...>%s' % ( head , tail )
        
    title = allright ( title )
    return title , visible_width ( title ) 

# =============================================================================
plain_string = lambda s : isinstance ( s , string_types ) and not '\n' in s
def get_item ( item ) :
    if   plain_string ( item )                : yield item
    elif isinstance   ( item , string_types ) :
        for line in item.split ( '\n' )       : yield line 
    ## sequence ? 
    else :
        for line in item : yield line 
# =============================================================================
## Preprcess table
def preprocess_table ( rows ) :
    """ Preprocess table
    """
    make_item  = lambda s : s if s else '' 
    for row in rows :
        if   plain_string ( row )                          : yield ( row , ) 
        elif all ( plain_string ( item ) for item in row ) : yield tuple ( row )
        else : 
            generators = tuple ( get_item ( item ) for item in row )
            for nrow in zip_longest ( *generators ) :
                yield tuple ( make_item ( i ) for i in nrow )   
# =============================================================================
## get the true  table width 
def table_width  ( table ) :
    """ Get the true width of the table
    """
    if terminaltables and isinstance ( table , terminaltables.AsciiTable ) :
        return table.table_width
    
    width = 0
    for row in table.split('\n') :
        width = max ( width , visible_width ( row ) ) 
    return width 

# =============================================================================
## Add certain prefix to  each line of the table
#  @code
#  table = ...
#  table =  add_prefix ( table , '# ') 
#  @endcode
def add_prefix ( table , prefix = '' ) :
    """ Add certain prefix to  each line of the table
    >>> table = ...
    >>> table =  add_prefix ( table , '# ') 
    """
    return prefix + table.replace ( '\n' , '\n' + prefix ) if prefix else table 

# ==============================================================================
## Align the certain column of the table
#  @code
#  aligned = align_column ( table , 1 , 'left' ) 
#  @endcode 
def align_column ( table , index , align = 'left') :
    """ Aling the certain column of the table
    >>> aligned = align_column ( table , 1 , 'left' ) 
    """
    nrows = [ list ( row ) for row in table ]
    lmax  =  0 

    for row in nrows :
        if index <= len ( row ) :
            item = decolorize ( row [ index ] )
            lmax = max ( lmax , len ( item )  )

    if not lmax : return table 

    aleft   =               align.lower() in left 
    aright  = not aleft and align.lower() in right 

    new_table = []
    for row in nrows :
        if index <= len ( row ) :
            item   = decolorize ( row [ index ] )
            nspace = lmax - len ( item ) 
            if   aleft :
                item = row [ index ] + nspace * ' '
            elif aright:
                item = nspace * ' ' + row [ index ]
            else :
                sl = nspace / 2
                sr = nspace - sl 
                item = sl * ' ' + row [ index ] + sr * ' '
            row[ index ] = item                        
        new_table.append ( row )
            
    return [ tuple ( row ) for row in new_table ] 

# =============================================================================
## Get the list of empty columns in the table:
#  @code
#  table = ...
#  empty = empty_columns ( table ) 
#  @endcode
def empty_columns ( table ) :
    """ Get the list of empty columns in the table
    >>> table = ...
    >>> empty = empty_columns ( table ) 
    """
    if not table : return ()
    
    maxcols = -1 
    for row in table : maxcols = max ( maxcols , len ( row ) )
    empty_cols = maxcols * [ True ]
    
    ## skip the first row (assumed to be the header)
    for row in table [ 1 : ] :
        for col , item in enumerate ( row ) :
            if item :  empty_cols [ col ] = False

    return tuple ( i for i, col in enumerate ( empty_cols ) if col )  
        
# ============================================================================
## Remove empty columns from the table
#  @code
#  table = ...
#  new_table = remove_empty_colums ( table ) 
#  @endcode
def remove_empty_columns ( table ) :
    """ Remove empty columns from the table
    >>> table = ...
    >>> new_table = remove_empty_colums ( table ) 
    """
    empty = empty_columns ( table )
    if not empty : return table
    empty = tuple ( reversed ( empty ) )
    
    newtable = []
    for row in table :
        newrow = list ( row    )
        lenrow = len  ( newrow ) 
        for i in empty:
            if i < lenrow : del newrow [ i ] 
        newtable.append ( tuple ( newrow ) )
        
    return tuple ( newtable )

    
# =============================================================================
if __name__ == '__main__' :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    if terminaltables :
        logger.info     ( "`terminaltables/%s' will be used for formatting"  % terminaltables.__name__ )
    else :
        logger.warning  ( "`terminaltables' is not available, use local replacement" )

    logger.info         ( " Default style is: `%s`" % default_style )
    
    table_data = [
        ( 'Name'  , 'Occupation' , 'Note' ) ,
        ( 'Alice' , '?'          , '---'  ) ,
        ( 'Bob'   , 'unemployed' , ''     ) ]

    logger.info ( 'The table is \n%s' % table     ( table_data , 'Title' , alignment = 'rrr' ) )

    
    logger.info ( 'The table is \n%s' % table     ( table_data , 'Title' , alignment = 'lll' ) ) 
    logger.info ( 'The table is \n%s' % table     ( table_data , 'Title' , alignment = 'ccc' ) ) 
    logger.info ( 'The table is \n%s' % the_table ( table_data , 'Title' , alignment = 'rrr' ) ) 
    logger.info ( 'The table is \n%s' % the_table ( table_data , 'Title' , alignment = 'lll' ) ) 
    logger.info ( 'The table is \n%s' % the_table ( table_data , 'Title' , alignment = 'ccc' ) ) 
    logger.info ( 'The table with prefix is \n%s' %
                  table     ( table_data , 'Title' , prefix = '# ' ) ) 
    logger.info ( 'The table with prefix is \n%s' %
                  the_table ( table_data , 'Title' , prefix = '# ' ) ) 

    for fmt in tabulate_styles + table_styles : 
        result = table ( table_data , title = 'Title' , style  = fmt , prefix = '# ' )            
        logger.info ( 'Use the format="%s":\n%s' % ( fmt , result ) )

    logger.info ( 'Available styles: \n%s' % ( '\n'.join ( table_styles + tabulate_styles ) ) ) 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
