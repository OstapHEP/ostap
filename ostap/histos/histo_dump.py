#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/histos/histo_dump.py
#  ASCII/Unicode histogram renderer for Ostap framework.
#  Features out-of-frame Underflow (UF) / Overflow (OV) bins, automatic 
#  computation of tick marks and scaling factors (10^N), flexible bin edges,
#  and safe handling of invalid/missing numeric values and errors.
# 
#  The code was co-developed and optimized with Gemini (Google AI). 
# 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2026-08-19
# =============================================================================
""" ASCII/Unicode histogram renderer for Ostap framework

 - Out-of-frame Underflow (UF) and Overflow (OV) bins.
 - Automatic computation of clean tick marks and 10^N scale factors for X and Y axes.
 - Flexible edge input: full edge lists or simple (xmin, xmax) tuples.
 - Safe handling of None, NaN, Inf values and missing error parameters.

 - The code was co-developed and optimized with Gemini (Google AI). 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2026-08-19"
__all__     = (
    'data2text'  , ## dump the content as ASCII/UNICODE/pseudographics chart 
    'data2text_' , ## dump the content as ASCII/UNICODE/pseudographics chart 
)
# =============================================================================
from   ostap.logger.pretty    import format_pow10 
import ostap.logger.colorized as     C
import math
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.histos.histo_dump' )
else                       : logger = getLogger( __name__                  )
# =============================================================================
logger.debug ( "ASCII/Unicode 1D histogram renderer" )
# =============================================================================
CIRCLE        = "●"
BAR           = "│"
RECTANGLE_POS = "░"
RECTANGLE_NEG = "▒"
CIRCLE_COLORED        = C.colored_string ( CIRCLE        , foreground = C.RED   )
BAR_COLORED           = C.colored_string ( BAR           , foreground = C.RED   )
RECTANGLE_POS_COLORED = C.colored_string ( RECTANGLE_POS , foreground = C.GREEN )
RECTANGLE_NEG_COLORED = C.colored_string ( RECTANGLE_NEG , foreground = C.BLUE  )
# =============================================================================
## Check if a value is a valid finite float
def _valid_value ( val ) :
    """ Check if a value is a valid finite float.
    """
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        return False if val is None else math.isfinite ( float ( val ) )
        # =====================================================================
    except ( TypeError , ValueError ) : # =====================================
        # =====================================================================
        return False

# =============================================================================
## Safely convert error input to a list of valid float errors, defaulting to 0.0
def _normalize_errors ( err_list , expected_len ) :
    """ Safely convert error input to a list of valid float errors, defaulting to 0.0.
    """
    if not err_list : return tuple ( [ 0.0 ] * expected_len ) 

    result = []
    for i in range ( expected_len ) :
        if i < len ( err_list ) and _valid_value ( err_list [ i ] ) :
            result.append ( max ( 0.0 , float ( err_list [ i ] ) ) ) 
        else :
            result.append ( 0.0 )
    return result

# =============================================================================
## Resolve X-axis bin edges
def _resolve_edges ( edges , n_bins ) :
    """ Resolve X-axis bin edges.
        Accepts explicit edge array (len = n_bins + 1) or a (xmin, xmax) tuple/list.
    """
    if edges is None : return [ float ( i ) for i in range ( n_bins + 1 ) ]

    # Handle (xmin, xmax) pair
    if len ( edges ) == 2 and n_bins >= 1 :
        x_min, x_max = float ( edges [ 0 ] ) , float ( edges [ 1 ] )
        step = ( x_max - x_min ) / n_bins
        return [ x_min + i * step for i in range ( n_bins + 1 ) ]

    return tuple ( float ( e ) if _valid_value ( e ) else float ( i ) for i , e in enumerate ( edges ) )

# =============================================================================
## Calculate human-friendly tick step and nice values in range
def _calculate_nice_ticks ( val_min , val_max , target_ticks = 6 ) :
    """ Calculate human-friendly tick step and nice values in range.
    """
    val_range = abs ( val_max - val_min )
    if val_range == 0 : val_range = 1.0

    raw_delta = val_range / target_ticks
    exponent  = math.floor ( math.log10 ( raw_delta ) ) if raw_delta > 0 else 0
    fraction  = raw_delta / ( 10 ** exponent ) if raw_delta > 0 else 1.0

    nice_steps = [ 1.0 , 2.0 , 2.5 , 5.0 , 10.0 ]
    delta      = 10.0 * ( 10 ** exponent )
    for step in nice_steps :
        if step >= fraction :
            delta = step * ( 10 ** exponent )
            break

    first_idx = math.ceil ( val_min / delta )
    last_idx  = math.floor ( val_max / delta )
    nice_vals = [ i * delta for i in range ( first_idx , last_idx + 1 ) ]

    return delta , nice_vals

# =============================================================================
## Generate Y-axis row specifications with steps and tick positions
def _build_grid_rows ( y_min , y_max , max_height ) :
    """ Generate Y-axis row specifications with steps and tick positions.
    """
    delta_y , _ = _calculate_nice_ticks ( y_min , y_max , target_ticks = max_height )

    steps_up       = int ( math.ceil ( y_max / delta_y ) ) if y_max > 0 else 0
    steps_down     = int ( math.ceil ( abs ( y_min ) / delta_y ) ) if y_min < 0 else 0
    lines_per_tick = 4

    rows = []
    for step in range ( steps_up , 0 , -1 ) :
        is_tick = ( step % lines_per_tick == 0 )
        rows.append ( { "y_high"   : step * delta_y ,
                        "y_low"    : ( step - 1 ) * delta_y ,
                        "tick_val" : round ( step * delta_y , 6 ) if is_tick else None ,
                        "is_tick"  : is_tick ,
                       "is_zero"  : False   } )

    rows.append ( { "y_high"   : 0.0  ,
                    "y_low"    : 0.0  ,
                    "tick_val" : 0.0  ,
                    "is_tick"  : True ,
                    "is_zero"  : True } )

    for step in range ( 1 , steps_down + 1 ) :
        is_tick = ( step % lines_per_tick == 0 )
        rows.append (
            {
                "y_high"   : -( step - 1 ) * delta_y ,
                "y_low"    : -step * delta_y ,
                "tick_val" : round ( -step * delta_y , 6 ) if is_tick else None ,
                "is_tick"  : is_tick ,
                "is_zero"  : False ,
            }
        )

    return rows , delta_y , steps_up

# =============================================================================
## Determine the ASCII character representation for a single cell
def the_glyph ( val          , 
                e_low        , 
                e_high       , 
                r_low        , 
                r_high       , 
                default_char , 
                has_errors   , 
                steps_up     , 
                delta_y      ,
                use_color = True ) :
    """ Determine the ASCII character representation for a single cell.
    """
    if not _valid_value ( val ) : return default_char

    val = float ( val )

    if has_errors :
        y_err_min = val - e_low
        y_err_max = val + e_high

        in_center = r_low <= val < r_high or ( r_high == steps_up * delta_y and val == r_high )
        in_error  = y_err_min < r_high and y_err_max > r_low

        if in_center  : return CIRCLE_COLORED if use_color else CIRCLE 
        elif in_error : return BAR_COLORED    if use_color else BAR 
    else :
        if   val > 0 and r_low >= 0  and val > r_low  : return RECTANGLE_POS_COLORED if use_color else RECTANGLE_POS 
        elif val < 0 and r_high <= 0 and val < r_high : return RECTANGLE_NEG_COLORED if use_color else RECTANGLE_NEG 

    return default_char

# =============================================================================
## Format vertical multi-line X-axis labels strictly aligned to the bottom row
def _render_vertical_x_labels ( marked_cols_dict , n_bins , prefix_len ) :
    """ Format vertical multi-line X-axis labels strictly aligned to the bottom row.
    """
    underflow_label = "UFLW"
    overflow_label  = "OFLW"

    all_labels   = list ( marked_cols_dict.values () ) + [ underflow_label , overflow_label ]
    max_len      = max ( len ( s ) for s in all_labels )
    prefix_space = " " * prefix_len

    lines = []
    for char_row_idx in range ( max_len ) :
        x_chars = [ " " ] * n_bins

        # Align X-axis labels strictly to bottom row
        for b_idx , label in marked_cols_dict.items () :
            offset       = max_len - len ( label )
            idx_in_label = char_row_idx - offset
            if 0 <= idx_in_label < len ( label ) :
                x_chars [ b_idx ] = label [ idx_in_label ]

        # Align U/FLOW and O/FLOW strictly to bottom row
        uf_offset = max_len - len ( underflow_label )
        uf_idx    = char_row_idx - uf_offset
        uf_char   = underflow_label [ uf_idx ] if 0 <= uf_idx < len ( underflow_label ) else " "

        ov_offset = max_len - len ( overflow_label )
        ov_idx    = char_row_idx - ov_offset
        ov_char   = overflow_label [ ov_idx ] if 0 <= ov_idx < len ( overflow_label ) else " "

        lines.append ( f"{prefix_space} {uf_char} {''.join ( x_chars )} {ov_char}" )

    return lines

# =============================================================================
## Main function to render ASCII histogram with UF/OV placed outside plot frame
def data2text_ ( bins               ,
                 underflow   = 0    ,
                 overflow    = 0    ,
                 edges       = None ,
                 errors_low  = None ,
                 errors_high = None ,
                 max_height  = 30   ,
                 use_color   = True ) :
    """ Main function to render ASCII histogram with UF/OV placed outside the plot frame.
    """
    full_data = [ underflow ] + list ( bins ) + [ overflow ]
    num_cols  = len ( full_data )
    n_bins    = len ( bins )

    edges     = _resolve_edges ( edges , n_bins )

    e_low_list  = _normalize_errors ( errors_low  , num_cols )
    e_high_list = _normalize_errors ( errors_high , num_cols )
    has_errors  = any ( e > 0 for e in e_low_list ) or any ( e > 0 for e in e_high_list )

    low_bounds , high_bounds = [] , []
    for i in range ( num_cols ) :
        v = full_data [ i ]
        if _valid_value ( v ) :
            val = float ( v )
            low_bounds.append  ( val - e_low_list  [ i ] )
            high_bounds.append ( val + e_high_list [ i ] )

    if not low_bounds : return "Histogram data is empty or invalid."

    raw_y_min , raw_y_max = min ( 0.0 , min ( low_bounds ) ) , max ( 0.0 , max ( high_bounds ) )
    rows , delta_y , steps_up = _build_grid_rows ( raw_y_min , raw_y_max , max_height )

    max_abs_y   = max ( abs ( raw_y_max ) , abs ( raw_y_min ) )
    y_exp       = int ( math.floor ( math.log10 ( max_abs_y ) ) ) if max_abs_y > 0 else 0
    y_scale     = 10 ** y_exp if ( y_exp >= 4 or y_exp <= -2 ) else 1.0
    ## y_scale_str = f" [10^{y_exp}]" if y_scale != 1.0 else ""
    y_scale_str = '[%s]' % format_pow10 ( y_exp ) if y_exp else ""

    max_abs_x   = max ( abs ( edges [ 0 ] ) , abs ( edges [ -1 ] ) )
    x_exp       = int ( math.floor ( math.log10 ( max_abs_x ) ) ) if max_abs_x > 0 else 0
    x_scale     = 10 ** x_exp if ( x_exp >= 4 or x_exp <= -2 ) else 1.0
    ## x_scale_str = f" [10^{x_exp}]" if x_scale != 1.0 else ""
    x_scale_str = '[%s]' % format_pow10 ( x_exp ) if x_exp else ""

    def fmt_x ( v ) :
        s = round ( v / x_scale , 6 )
        return f"{s:+g}" if s != 0 else "+0"

    _ , nice_x_vals  = _calculate_nice_ticks ( edges [ 0 ] , edges [ -1 ] , target_ticks = 6 )
    marked_cols_set  = set ()
    marked_cols_dict = { 0 : fmt_x ( edges [ 0 ] ) , n_bins - 1 : fmt_x ( edges [ -1 ] ) }

    for xv in nice_x_vals :
        b_idx = min ( range ( n_bins ) , key = lambda i : abs ( xv - edges [ i ] ) )
        if 2 <= b_idx <= n_bins - 3 :
            marked_cols_set.add ( b_idx )
            marked_cols_dict [ b_idx ] = fmt_x ( xv )

    prefix_len    = 6
    output        = []
    header_prefix = " " * prefix_len + "   "
    output.append ( f"{header_prefix}▲ Y{y_scale_str}".rstrip () )

    for idx , row in enumerate ( rows ) :
        is_top , is_bottom = ( idx == 0 ) , ( idx == len ( rows ) - 1 )

        if is_top :
            top_cols = [ "┬" if b in marked_cols_set else "─" for b in range ( n_bins ) ]
            output.append ( f"{' ' * prefix_len}  ┌{''.join ( top_cols )}┐" )

        if row [ "is_zero" ] :
            y_label = f"{'+0':>{prefix_len}}"
        elif row [ "is_tick" ] and row [ "tick_val" ] is not None :
            s_val   = row [ "tick_val" ] / y_scale
            val_str = f"{s_val:+g}" if s_val != 0 else "+0"
            y_label = f"{val_str:>{prefix_len}}"
        else :
            y_label = " " * prefix_len

        left_b  = "┼" if row [ "is_zero" ] else ( "┤" if row [ "is_tick" ] else "│" )
        right_b = "┼" if row [ "is_zero" ] or row [ "is_tick" ] else "│"

        uf_char = the_glyph ( full_data   [ 0 ] ,
                               e_low_list  [ 0 ] ,
                               e_high_list [ 0 ] ,
                               row         [ "y_low"  ] ,
                               row         [ "y_high" ] ,
                               " "         ,
                               has_errors  ,
                               steps_up    ,
                               delta_y     ,
                               use_color = use_color )
        ov_char = the_glyph ( full_data   [ -1 ] ,
                              e_low_list  [ -1 ] ,
                              e_high_list [ -1 ] ,
                              row         [ "y_low"  ] ,
                              row         [ "y_high" ] ,
                              " "         ,
                              has_errors  ,
                              steps_up    ,
                              delta_y     , 
                              use_color = use_color )
        
        bin_chars = []
        for i in range ( 1 , n_bins + 1 ) :
            b_idx       = i - 1
            is_tick_col = b_idx in marked_cols_set
            def_c       = (
                "┼"
                if ( row [ "is_tick" ] and is_tick_col )
                else (
                    "⋯"
                    if row [ "is_tick" ]
                    else ( "┊" if is_tick_col else " " )
                )
            )
            if row [ "is_zero" ] : def_c = "┼" if is_tick_col else "─"

            c = the_glyph ( bins         [ b_idx ] ,
                            e_low_list   [ i ] ,
                            e_high_list  [ i ] ,
                            row          [ "y_low"  ] ,
                            row          [ "y_high" ] ,
                            def_c        ,
                            has_errors   ,
                            steps_up     ,
                            delta_y      , 
                            use_color = use_color )
            
            bin_chars.append ( c )

        x_arrow = f"► X{x_scale_str}" if row [ "is_zero" ] else ""
        output.append ( f"{y_label} {uf_char}{left_b}{''.join ( bin_chars )}{right_b}{ov_char}{x_arrow}" )

        if is_bottom :
            bot_cols = [ "┴" if b in marked_cols_set else "─" for b in range ( n_bins ) ]
            output.append ( f"{' ' * prefix_len}  └{''.join ( bot_cols )}┘" )

    output.extend ( _render_vertical_x_labels ( marked_cols_dict , n_bins , prefix_len ) )

    return "\n".join ( output )

# =============================================================================
## Wrapper function to handle symmetric and asymmetric error formats
def data2text ( bins               , * , 
                underflow   = 0    ,
                overflow    = 0    ,
                edges       = None ,
                errors      = None ,
                max_height  = 26   ,
                use_color   = True ) :
    """ Wrapper function to handle symmetric and asymmetric error formats.
    """
    return data2text_ (
        bins        = bins        ,
        underflow   = underflow   ,
        overflow    = overflow    ,
        edges       = edges       ,
        errors_low  = errors      ,
        errors_high = errors      ,
        max_height  = max_height  ,
        use_color   = use_color   
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ )
    
    import ROOT, random
    logger.info ( "=== Test 1: Standard Histogram with (xmin, xmax) tuple ===" )

    th1 = ROOT.TH1F ( "h1" , "Test Histogram" , 70 , 0  , 100  )
    for i in range ( 500 ) :
        v = random.gauss ( 0 , 100  )
        th1.Fill ( v , -1 if v >= 50  else 2 )

    nbins      = th1.GetNbinsX ()
    values     = tuple ( th1.GetBinContent ( i ) for i in range ( 1 , nbins + 1 ) )
    errors     = tuple ( th1.GetBinError   ( i ) for i in range ( 0 , nbins + 2 ) )
    underflow  = th1.GetBinContent ( 0         )
    overflow   = th1.GetBinContent ( nbins + 1 )
    axis       = th1.GetXaxis() 
    edges      = axis.GetXmin() , axis.GetXmax()
    logger.info ( 'Histogram with error bars:\n%s'    % data2text ( values    ,
                                                                    errors    = errors    ,
                                                                    underflow = underflow , 
                                                                    overflow  = overflow  , 
                                                                    edges     = edges     ,
                                                                    use_color = True      ) )
    logger.info ( 'Histogram without error bars:\n%s' % data2text ( values    ,
                                                                    underflow = underflow , 
                                                                    overflow  = overflow  , 
                                                                    edges     = edges     ,
                                                                    use_color = True      ) )
    

# =============================================================================
##                                                                      The END 
# =============================================================================
