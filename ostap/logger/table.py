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
    try : # ===================================================================
        # =====================================================================
        visible_width = terminaltables.width_and_alignment.visible_width
        # =====================================================================
    except ( AttributeError,NameError ) : # ===================================
        # =====================================================================
        visible_width = None 

# =============================================================================
if not visible_width : # ======================================================
    # =========================================================================
    import unicodedata 
    ## visible length of the string/expression
    def visible_width ( what ) :
        """ Visible length of the string/expression (local, use unicodedata) 
        """
        if not what : return 0
        item = decolorize ( what )
        if not item : return 0
        ##
        if  sys.version_info < ( 3 , 0 ) : item = item.decode("u8")
        ##
        width = 0 
        for char in item :
            if unicodedata.east_asian_width ( char ) in ( 'F' , 'W' ) : width += 2
            else                                                      : width += 1
        return width 

# ==============================================================================
## available table styles 
table_styles    = 'local' , 'default'         ## use local, handmade replacement
ascii_styles    = 'local' , 
terminal_styles = () 
tabulate_styles = () 
# ==============================================================================
if terminaltables : # ==========================================================
    ## Vaild `style` arguments (case insensitive) 
    terminal_styles += ( 'ascii'     , ## use `AsciiTable` 
                         'single'    , ## use `SingleTable` 
                         'porcelain' , ## use `PorcelainTable` 
                         'github'    , ## use `GithubFlavoredMarkdownTable`
                         'markdown'  , ## use `GithubFlavoredMarkdownTable`
                         'double'    , ## use `DoubleTable` 
                         'default'   ) ## use `DoubleTable` 
    table_styles += terminal_styles
    ascii_styles += ( 'ascii' , ) 
# ==============================================================================
## very light wrapper for tabulate tables
#  @see https://github.com/astanin/python-tabulate
# =============================================================================
try : # =======================================================================
    # =========================================================================
    #  @see https://github.com/astanin/python-tabulate
    import tabulate
    # =========================================================================
    ## available formats 
    tabulate_styles  = tuple ( tabulate.tabulate_formats ) 
    table_styles    += tabulate_styles
    ascii_styles    += tuple ( s for s in sorted ( (
        'asciidoc'  , 'github'     , 'markdown'       , 'grid'            ,
        'jira'      , 'mediawiki'  , 'moinmoin'       , 'org_tbl'         ,
        'outline'   , 'pipe'       , 'plain'          , 'html'            ,
        'presto'    , 'pretty'     , 'psql'           , 'rst'             ,
        'simple'    , 'texile'     , 'tsv'            , 'unsafehtml'      ,
        'yourtrack' ,  'latex'     , 'latex_booktabs' , 'latex_longtable' ,
        'latex_raw' ) ) if s in tabulate_styles )
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    ## no formats are available
    tabulate         = None

# =============================================================================
## default style
default_style = 'default'
# =============================================================================
from ostap.utils.env import has_env, get_env, OSTAP_TABLE 
# =============================================================================
if  has_env ( OSTAP_TABLE ) :
    # =========================================================================
    default_style = get_env ( OSTAP_TABLE , default_style )
    # =========================================================================
else : # ======================================================================
    # =========================================================================
    ## get the preferred table style from the configuration file(s)
    import ostap.core.config as OCC
    if 'STYLE' in OCC.tables : default_style = OCC.tables.get ( 'STYLE' , fallback = default_style )
    
# =============================================================================
default_style = default_style.lower().strip()
if not default_style in table_styles : 
    if   isatty () and terminaltables : default_style = 'double'
    elif isatty () and tabulate       : default_style = 'fancy_grid'
    elif               terminaltables : default_style = 'ascii' 
    else                              : default_style = 'local' 
# =============================================================================

# =============================================================================
## finally adjust the style
# =============================================================================

# =============================================================================
columns_left    = '<' , 'l' , 'left'  
columns_right   = '>' , 'r' , 'right' 
columns_center  = '^' , 'c' , 'center' , '='
columns_wrapped = 'w' , 'p' , 'wrap'   , 'wrapped'
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
                maxwidth        = -1          , **kwargs ) :
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

    title_bw = decolorize    ( title    ).strip ()
    title    = allright      ( title_bw )    
    twidth   = visible_width ( title    ) ## visible lenght of the title 
    pwidth   = visible_width ( prefix   ) ## visible lenght of the prefix 

    if not maxwidth or maxwidth <= pwidth : maxwidth = terminal_size() [ 0 ]
    
    ## take care on "wrapped&wrapped-like" columns 
    rows = [ list ( row ) for row in preprocess_table ( rows ) ]
    
    if not rows : return ''

    # =================================================================
    ## Basis structure 
    # =================================================================

    ## (max)number of columns
    num_cols    = max ( len ( row ) for row in rows )        
    ## take care on "wrapped&wrapped-like" columns 
    rows        = [ list ( row ) for row in preprocess_table ( rows ) ]
    ## equalize columns: ensure that all rows have the same (max) length  
    rows        = [ row for row in equalize_columns ( rows , num_cols ) ]
    ## widths of preprocessed columns 
    widths      = column_widths ( rows )

    # =================================================================
    ## Column alignment:
    # =================================================================
    
    cols_wrap   = column_types ( alignment , columns_wrapped , num_cols ) 
    cols_left   = column_types ( alignment , columns_left    , num_cols )
    cols_right  = column_types ( alignment , columns_right   , num_cols )
    cols_center = column_types ( alignment , columns_center  , num_cols ) 

    if not maxwidth or maxwidth <= pwidth : maxwidth = terminal_size() [ 0 ]

    ## play with wrapped columns
    rows , widths = wrap_cells ( rows , widths , cols_wrap , maxwidth , indent , twidth, pwidth ) 

    ## actual "formatting"
    for r , row in enumerate ( rows ) :
        for c, cell in enumerate ( row ) :
            cell_width = visible_width ( cell )
            delta = widths [ c ] - cell_width
            if   delta < 0  : cell = widths [ c ] * '*'
            if   0 == delta : pass 
            if   c in cols_right           : cell = ( delta * ' ' ) + cell ## extend 
            elif c in cols_left            : cell = cell + ( delta * ' ' ) ## extend 
            elif c in cols_wrap and 0 != r : cell = cell + ( delta * ' ' ) ## extend 
            else :  
                d1 , d2 = divmod ( delta , 2 ) 
                cell = ( ( d1 + d2 )  * ' ' ) + cell + ( d1 * ' ' ) ## extend    
            rows [ r ] [ c ] = cell

    pad_left  = 1
    pad_right = 1
    
    tabwidth = sum ( widths ) + num_cols * ( pad_left + pad_right ) + ( num_cols + 1 )

    if pwidth + tabwidth < maxwidth and twidth and title :
        ## try to stretch the table, if possible 
        delta = ( pwidth + 2 + twidth  ) - ( pwidth + tabwidth )        
        if 0 < delta :
            d1 , d2    = divmod ( delta , 2 )
            pad_left  += d1 + d2
            pad_right += d1
        tabwidth = sum ( widths ) + num_cols * ( pad_left + pad_right ) + ( num_cols + 1 )
        
    glue = pad_left * ' '  + '|' + pad_right * ' '

    seps     = [   '-' * ( widths [ c ]  + pad_left + pad_right  ) for c in range ( num_cols ) ]
    sepline  = '+' + "+".join (  seps ) +  '+'
        
    table = [ sepline ]
    
    if title and twidth :

        if tabwidth < twidth + 4 :
            title , twidth = shorten_title ( title , tabwidth - 4 )

        d1 , d2 = divmod ( tabwidth - 4 - twidth , 2 )

        line = '| ' + allright ( ( d1 + d2 ) * ' ' + decolorize ( title ) + d1 * ' ' ) + ' |'

        table.append ( line    ) 
        table.append ( sepline )
                
    for r , row in enumerate ( rows ) :
        if colorize_header and 0 == r : row  = [ infostr ( decolorize ( cell ) ) for cell in row ]        
        line = '|' + pad_left * ' ' + glue.join ( row ) + pad_right * ' ' + '|'        
        table.append ( line )
        if 0 == r : table.append ( sepline )
        
    if table [ -1 ] != sepline : table.append ( sepline )

    
    table = '\n'.join ( table )    
    if sys.version_info < ( 3 , 0 ) :
        if isinstance ( table , unicode ) : table = table.encode ('utf-8')
            
    if kwargs : logger.warning ( 'Ignore keyword arguments: %s' %  ( ', '.join ( key for key in kwargs ) ) )             
    return add_prefix ( table , prefix ) if prefix else table 

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
    """

    from ostap.utils.basic import isatty

    title_bw = decolorize  ( title    ).strip ()
    title    = allright    ( title_bw )    
    twidth   = visible_width ( title  ) ## visible lenght of the title 
    pwidth   = visible_width ( prefix ) ## visible lenght of the prefix 

    if not maxwidth or maxwidth <= pwidth : maxwidth = terminal_size() [ 0 ]

    rows = [ list ( row ) for row in preprocess_table ( rows ) ]

    if rows and colorize_header : 
        header_row = rows [ 0 ]
        header_row = [ infostr ( decolorize ( c ) ) for c in header_row ]
        rows [ 0 ] = header_row 

    if not style : style = '%s' % default_style    
    fmt = style.lower()

    if   not  fmt in table_styles : fmt = 'local' ## switch to local

    if isatty()                   : pass
    elif not fmt in ascii_styles  : fmt = 'local' ## swicth to local 

    # =================================================================
    ## Basis structure 
    # =================================================================

    ## (max)number of columns
    num_cols    = max ( len ( row ) for row in rows )        
    ## take care on "wrapped&wrapped-like" columns 
    rows        = [ list ( row ) for row in preprocess_table ( rows ) ]
    ## equalize columns: ensure that all rows have the same (max) length  
    rows        = [ row for row in equalize_columns ( rows , num_cols ) ]
    ## widths of preprocessed columns    
    widths      = column_widths ( rows )

    # =================================================================
    ## Column alignment:
    # =================================================================
    
    cols_wrap   = column_types ( alignment , columns_wrapped , num_cols ) 
    cols_left   = column_types ( alignment , columns_left    , num_cols )
    cols_right  = column_types ( alignment , columns_right   , num_cols )
    cols_center = column_types ( alignment , columns_center  , num_cols ) 
    
    ## 
    tty = isatty ()
    
    # =========================================================================
    ## Terminal tables
    # =========================================================================
    if terminaltables and fmt in terminal_styles and ( tty or fmt in ascii_styles ) : 
                                 
        if   'ascii'  == fmt or not tty       : table_ = terminaltables.AsciiTable                  ( rows , title )        
        elif 'single' == fmt                  : table_ = terminaltables.SingleTable                 ( rows , title )
        elif 'porcelain' == fmt               : table_ = terminaltables.PorcelainTable              ( rows         )
        elif fmt in ( 'github' , 'markdown' ) : table_ = terminaltables.GithubFlavoredMarkdownTable ( rows         )
        elif 'double' == fmt                  : table_ = terminaltables.DoubleTable                 ( rows , title )
        else                                  : table_ = terminaltables.DoubleTable                 ( rows , title )
            
        cw       = table_.column_widths
        num_cols = len ( cw )
        
        wrap_width = 15
        if cols_wrap :
            max_dimensions = terminaltables.width_and_alignment.max_dimensions
            col_widths     = max_dimensions ( table_.table_data    ,
                                              table_.padding_left  ,
                                              table_.padding_right ) [2]
            tw  = sum ( l for ( i , l ) in enumerate ( col_widths ) if not i in cols_wrap )
            tw += num_cols + 1 ## vertical stulls             
            tw += ( table_.padding_left + table_.padding_right )  * len ( cols_wrap )
            tw += pwidth
            
            max_width = maxwidth 

            wrap_total     = max_width - tw 
            wrap_width , _ = divmod ( wrap_total , len ( cols_wrap  ) )            
            wrap_width     = max ( wrap_width , 15  )

            ## adjust wrapping
            
            narrow = 0 
            used   = 0 
            for c in cols_wrap :
                colw  = col_widths [ c ] 
                delta = wrap_width - colw
                if delta <= 0 : continue
                narrow += 1 
                used   += colw 

            if narrow and narrow < len ( cols_wrap )  and used :
                wrap_width , _ = divmod ( wrap_total - used  , len ( cols_wrap ) - narrow  )            
                wrap_width     = max ( wrap_width , 15  )
            
            for c in cols_wrap   :
                for l , line in enumerate ( table_.table_data ) :
                    cell  = decolorize ( line [ c ] )                    
                    width = min ( wrap_width , col_widths [ c ] )
                    table_.table_data [ l ] [ c ] = textwrap. fill ( cell , width = width , initial_indent = indent )

            if colorize_header and table_.table_data :
                header_row              = table_.table_data [ 0 ] 
                header_row              = [ infostr ( decolorize ( cell ) ) for cell in header_row ]
                table_.table_data [ 0 ] = header_row 
            
            wrap_width , _ = divmod ( max_width  - tw  , len ( cols_wrap  ) )        
            wrap_width     = max ( wrap_width , 15  )
        
        for i in cols_right  : table_.justify_columns [ i ] = 'right'
        for i in cols_left   : table_.justify_columns [ i ] = 'left'
        for i in cols_center : table_.justify_columns [ i ] = 'center'
                        
        if title and style in ( 'github' , 'markdown' ) :
            title   = '# ' + title
            twidth += 2

        tabwidth = table_.table_width
        if title and twidth and tabwidth < twidth + 2 :
            ## is there a room for stretch ?
            if tabwidth + pwidth <= maxwidth :            
                newwidth = min ( twidth + 2 + pwidth , maxwidth  )
                if tabwidth + pwidth <= newwidth : 
                    a  , _  = divmod ( newwidth - tabwidth - pwidth , num_cols )
                    a1 , r  = divmod ( a        , 2 )
                    table_.padding_left  += a1 + r 
                    table_.padding_right += a1  
                    tabwidth = table_.table_width
                    
        if tabwidth < twidth + 2 :
            title , twidth = shorten_title ( title , tabwidth - 2 )
            table_.title = title
                        
        ## get the final table 
        table = table_.table

        ## special final adjustment 
        if title  and fmt in ( 'porcelain' , ) :
            line  = twidth * '-'
            table = title + '\n' + line + '\n' + table 
        elif title and fmt in ( 'github' , 'markdown' ) :
            line  = '---'
            table = line  + '\n' + title + '\n' + line + '\n' + table 

        if sys.version_info < ( 3 , 0 ) :
            if isinstance ( table , unicode ) : table = table.encode ('utf-8')
        
        if kwargs : logger.warning ( 'Ignore keyword arguments: %s' %  ( ', '.join ( key for key in kwargs ) ) )             
        return add_prefix ( table  , prefix ) if prefix else table 

    # =========================================================================
    ## Tanulate
    # =========================================================================
    elif tabulate  and fmt in tabulate_styles and ( tty or  fmt in ascii_styles ) : 
            
        colalign = kwargs.get ( 'colalign' , [] ) 
        if alignment and not colalign :
            for i in range ( num_cols ) :
                if   i in cols_left   : colalign.append ( 'leff'   ) 
                elif i in cols_right  : colalign.append ( 'right'  ) 
                elif i in cols_center : colalign.append ( 'center' ) 
                else                  : colalign.append ( 'left'   ) # #ATTENTION 
        colalign = tuple ( colalign )
        
        ## play with wrapped columns
        rows , widths = wrap_cells ( rows , widths , cols_wrap , maxwidth , indent , twidth, pwidth ) 

        if colalign  : table = tabulate.tabulate ( rows , tablefmt = fmt , colalign = colalign , **kwargs )
        else         : table = tabulate.tabulate ( rows , tablefmt = fmt                       , **kwargs )
        
        if title :
            tabwidth = table_width ( table  )
            if tabwidth < twidth + 2 :                
                title , twidth = shorten_title ( title , tabwidth - 2 )
            if fmt.startswith( 'latex' ) : title = '%% ' + title
            d1 , d2  = divmod ( tabwidth  - 2 - twidth , 2 )              
            title    = allright ( ( d1 + d2 ) * ' ' + decolorize ( title ) + d1 * ' ' )                                 
            table    = title + '\n' + table  

        if   fmt.startswith( 'latex' ) : table = decolorize ( table )
        elif 'html' in fmt             : table = decolorize ( table )
        
        if sys.version_info < ( 3 , 0 ) :
            if isinstance ( table , unicode ) : table = table.encode ('utf-8')
        
        if kwargs : logger.warning ( 'Ignore keyword arguments: %s' %  ( ', '.join ( key for key in kwargs ) ) )             
        return add_prefix ( table  , prefix ) if prefix else table 

    # =========================================================================
    ## Local
    # =========================================================================    
    return the_table ( rows                              ,
                       title           = title           , 
                       prefix          = prefix          ,
                       alignment       = alignment       ,
                       wrap_width      = wrap_width      , 
                       colorize_header = colorize_header ,
                       indent          = indent          ,
                       style           = style           ,
                       maxwidth        = maxwidth        , **kwargs )

# =============================================================================
## Shorten the title
def shorten_title ( title , size ) :
    """ Shorten the title
    """
    twidth = visible_width ( title )
    from ostap.logger.symbols import ellipsis 
    inset  = '{%s}' % ellipsis 
    insize  = len ( inset ) 
    size   = max ( 5 , size , insize +  2  ) 
    if twidth <= size : return title , twidth
    title  = decolorize ( title )

    a , r = divmod ( size - 5 , 2 )
    head = title [      : a ] 
    tail = title [ -a-r :   ] if r else title [ -a : ]
    head = head.strip()
    tail = tail.strip()
    title = '%s%s%s' % ( head , inset , tail )
        
    title = allright ( title )
    return title , visible_width ( title ) 

# =============================================================================
plain_string = lambda s : isinstance ( s , string_types ) and not '\n' in s
def get_item ( item ) :
    if   plain_string ( item )                : yield item
    elif isinstance   ( item , string_types ) :
        ## spliting can hav ebad interference with colorization
        dcitem = decolorize ( item ) 
        for line in dcitem.split ( '\n' )       : yield line 
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
        if   plain_string ( row )                          : yield          row ,  
        elif all ( plain_string ( item ) for item in row ) : yield  tuple ( row )
        else :
            generators = tuple ( get_item ( item ) for item in row )
            for nrow in zip_longest ( *generators ) :
                yield tuple ( make_item ( i ) for i in nrow )

# =============================================================================
## "wrap" certan cells in the table uinsg `textwrap.fill`
def wrap_cells ( rows , widths , wrapped , maxwidth , indent = '' , twidth = 0 , pwidth = 0 ) :
    """ `Wrap' certan cells in the table uinsg `textwrap.fill`
    """

    if not wrapped : return rows, widths  ## no action 
    
    num_cols = len ( widths )
    
    ## total width of non-wrapped columns + 2 blanks per column plus (nc+1) vertical separators +prefix 
    tw   = sum ( widths [ i ] for i in range ( num_cols ) if not i in wrapped )
    
    tw  += 2 * num_cols ##  two padding blanks per each column 
    tw  += num_cols + 1 ##  vertical separators     
    tw  += pwidth       ##  prefix
    
    max_width = maxwidth 
    ## if twidth : max_width = min ( pwidth + 2 + twidth , max_width )

    wrap_total     = max_width - tw 
    wrap_width , _ = divmod ( wrap_total , len ( wrapped ) ) 
    wrap_width     = max    ( wrap_width , 15 )

    # adjust wrapping:
    
    used   = 0
    narrow = 0 
    for c in wrapped :
        colw  = widths [ c ] 
        delta = wrap_width - colw
        if delta <= 0 : continue
        narrow += 1 
        used   += colw 
        
    if narrow and narrow < len ( wrapped ) and used :
        wrap_width , _ = divmod ( wrap_total - used  , len ( wrapped ) - narrow )            
        wrap_width     = max ( wrap_width , 15  )
            
    ## make wrapping
    reprocess = False  
    for row in rows :
        for c in wrapped :
            if wrap_width < widths [ c ] :
                cell  = decolorize ( row [ c ] )
                width = min ( wrap_width , widths [ c ] )
                row [ c ] = textwrap.fill ( cell , width = width , initial_indent = indent )
                reprocess = True
                
    if reprocess : 
        ## re-process wrapped table 
        rows = [ list ( row ) for row in preprocess_table ( rows ) ]        
        ## adjust new width of preprocessed columns 
        widths      = column_widths ( rows )

    return rows , widths
    

# =============================================================================
## Generator that makes all rows of equal (max) length 
def equalize_columns ( rows , numcols = None  ) :
    """ Generator that makes all rows of equal (max) length 
    """
    if not numcols or numcols <= 0 : numcols  = max ( len ( row ) for row in rows )
    for row in rows :
        row = list ( row ) 
        while len  ( row ) < numcols : row.append ( '' ) ## add empty cell
        yield row

# =============================================================================
## Get indices of certain columns 
def column_types ( alignment , col_types , num_cols ) :
    """ Get indices of wrapped columns
    """
    return tuple ( i for ( i, a ) in enumerate ( alignment ) if ( i < num_cols and a and isinstance ( a , str ) and a.lower() in col_types ) ) if alignment else ()

# =============================================================================
## Get indices of wrapped columns 
def wrapped_columns ( alignment , num_cols ) :
    """ Get indices of wrapped columns
    """
    return column_types ( alignment , columns_wrapped , num_cols )

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
## widths of columns
#  - rows are assumed to be preprocessed!
def column_widths ( rows ) :
    """ Widths of columns in the table
    - rows are assumed to be preprocessed!
    """
    widths = max ( len ( row ) for row in rows ) * [ 0 ]
    for row in rows :
        for col , cell in enumerate ( row ) :
            widths [ col ] = max ( widths [ col ] , visible_width ( cell ) )            
    return tuple ( widths )
    
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
    """ Align the certain column of the table
    >>> aligned = align_column ( table , 1 , 'left' ) 
    """
    nrows = [ list ( row ) for row in table ]
    lmax  =  0 

    for row in nrows :
        if index <= len ( row ) :
            item = decolorize ( row [ index ] )
            lmax = max ( lmax , len ( item )  )

    if not lmax : return table 

    aleft   =               align.lower() in columns_left 
    aright  = not aleft and align.lower() in columns_right 

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
            cell = decolorize ( item ).strip() 
            if cell :  empty_cols [ col ] = False

    return tuple ( i for i, col in enumerate ( empty_cols ) if col )  
        
# ============================================================================
## Remove empty columns from the table
#  @code
#  table = ...
#  new_table = remove_empty_columns ( table ) 
#  @endcode
def remove_empty_columns ( table ) :
    """ Remove empty columns from the table
    >>> table = ...
    >>> new_table = remove_empty_columns ( table ) 
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

    all_styles  = terminal_styles + tabulate_styles
    all_styles  = tabulate_styles
    
    for fmt in sorted ( all_styles ) : 
        result = table ( table_data , title = 'Title' , style  = fmt , prefix = '# ' )            
        logger.info ( 'Use the format="%s":\n%s' % ( fmt , result ) )

    rows = [ ( 'Style' , 'From' ) ]
    
    for s in sorted ( all_styles ) :
        if    'local' == s                           : row = s , 'Local'
        elif s in terminal_styles and terminaltables : row = s , 'terminaltables'
        elif s in tabulate_styles and tabulate       : row = s , 'tabulate'
        else                                         : row = s , ''
        rows.append ( row ) 
        
    title = 'All styles'
    logger.info ( '%s: \n%s' % ( title , table ( rows              ,
                                                 title     = title ,
                                                 prefix    = '# '  ,
                                                 alignment = 'll'  ) ) ) 

    logger.info  ('Default style: %s' % default_style )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
