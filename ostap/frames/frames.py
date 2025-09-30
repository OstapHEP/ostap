#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/frames/frames.py
#  Module with decoration of TDataFrame objects for efficient use in python
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2018-06-16
# =============================================================================
""" Decoration of Tree/Chain objects for efficient use in python"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'DataFrame'             , ## RDataFrame object
    ##
    'report_print_table'   , ## print the report 
    'report_as_table'      , ## print the report
    'report_print'         , ## print the report 
    ##
    'frame_columns'         , ## all defined columns 
    'frame_branches'        , ## all defined columns  
    ## 
    'frame_progress1'       , ## progress bar for frame (OK for ROOT<6.32)
    'frame_progress2'       , ## progress bar for frame (OK for ROOT>6.29) 
    'frame_progress'        , ## progress bar for frame
    ## 
    'frame_prescale'        , ## prescale frame (trivial filter)
    'frame_print'           , ## over-simplified print frame
    'frame_table'           , ## Print frame as detailed table 
    ## 
    'frame_length'          , ## length/size/#entries of the frame 
    'frame_size'            , ## length/size/#entries of the frame 
    ##
    'frame_statistic'       , ## statistics for the variable(s) 
    'frame_minmax'          , ## min/max for the variable(s) 
    'frame_range'           , ## ranges  for the variable(s)
    ##
    'frame_the_moment'      , ## get moments for the variable(s)
    ##
    'frame_nEff'            , ## number of effective entries 
    'frame_mean'            , ## mean value for variable(s)
    'frame_variance'        , ## variance   for variable(s)
    'frame_dispersion'      , ## dispersion for variable(s)
    'frame_rms'             , ## rms for variable(s)
    'frame_skewness'        , ## skewness for variable(s)
    'frame_kurtosis'        , ## skewness for variable(s)
    ##
    'frame_arithmetic_mean' , ## get the arithmetc mean for variables 
    'frame_harmonic_mean'   , ## get the harmonic  mean for variables 
    'frame_geometric_mean'  , ## get the geometric mean for variables 
    'frame_power_mean'      , ## get the power     mean for variables 
    'frame_lehmer_mean'     , ## get the Lehmer    mean for variables 
    ## 
    'frame_ECDF'            , ## get the empirical cumulative distribution functions for variable(s)
    'frame_slice'           , ## get the slice as numpy array 
    ## 
    'frame_project'         , ## project data frame to the (1D/2D/3D) histogram    
    'frame_param'           , ## parameterize data n flight     
    'frame_draw'            , ## draw variable from the frame
    ## 
    'frame_types'          , ## the basic DataFrame/FeameNode types 
    ) 
# =============================================================================
from   ostap.core.meta_info      import root_info 
from   ostap.core.ostap_types    import ( integer_types , dictlike_types , 
                                          num_types     , ordered_dict   )    
from   ostap.core.core           import cpp, Ostap, SE , WSE 
from   ostap.math.base           import isequal, iszero, axis_range                             
from   ostap.logger.utils        import multicolumn
from   ostap.utils.cidict        import cidict, cidict_fun      
from   ostap.utils.progress_conf import progress_conf 
from   ostap.utils.basic         import isatty, loop_items, typename 
from   ostap.frames.frame2histo  import ( DF_P2Model , DF_P1Model , 
                                          DF_H3Model , DF_H2Model , DF_H1Model ) 
import ostap.core.config         as     OCC 
import ostap.stats.statvars      as     SV
import ostap.logger.table        as     T
import ROOT, math
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.frames.frames' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Some useful decorations for ROOT::RDataFrame objects')
# =============================================================================
## @see Ostap/DataFrame.h 
DataFrame = Ostap.DataFrame
## @see Ostap/DataFrame.h 
FrameNode = Ostap.FrameNode
# =============================================================================
if    FrameNode is DataFrame : frame_types = DataFrame ,
else                         : frame_types = DataFrame , FrameNode
# =============================================================================
if not hasattr ( Ostap ,'DataFrame' ) :  Ostap.DataFrame = DataFrame
if not hasattr ( Ostap ,'FrameNode' ) :  Ostap.FrameNode = FrameNode 
# =============================================================================
if not hasattr ( ROOT.TTree  , '__len__') : ROOT.TTree .__len__ = lambda s : s.GetEntries() 
if not hasattr ( ROOT.TChain , '__len__') : ROOT.TChain.__len__ = lambda s : s.GetEntries()
std_move = ROOT.std.move
# =============================================================================
## The shortcuts for actions  
SA1  = ROOT.Detail.RDF.StatAction1
SA1w = ROOT.Detail.RDF.StatAction1w 
SA2  = ROOT.Detail.RDF.StatAction2
SA2w = ROOT.Detail.RDF.StatAction2w 
SA3  = ROOT.Detail.RDF.StatAction3
SA3w = ROOT.Detail.RDF.StatAction3w 
SA4  = ROOT.Detail.RDF.StatAction4
SA4w = ROOT.Detail.RDF.StatAction4w
# ================================================================================    
## type for the column names 
CNT  = DataFrame.ColumnNames_t 
# =============================================================================
## get frame-like stuff  as rnode 
def as_rnode ( frame ) :
    """ Get frame-like stuff  as rnode"""
    return ROOT.RDF.AsRNode ( frame ) 
# ==============================================================================
## Local helper function  to get the lazy results
def get_values ( results                               , * , 
                 frame     = None                      ,
                 transform = lambda s : s.GetValue ()  , 
                 report    = False                     ,
                 title     = 'DataFrame'               ) :
    """ Local helper function to get the Lazy results 
    """
    assert callable ( transform ) , "`transform` must be callable!"
    
    if not transform : 
        transform = lambda s : s
        
    if not isinstance ( results , dictlike_types ) :
        if frame and report :
            report = frame.Report()
            logger.info ( '%s\n%s' % ( title , report_print ( report , title = title , prefix = '# ') ) )
        return transform ( results ) 
    ## 
    values = {}
    for key, value in loop_items ( results ) : 
        values [ key ] = transform ( value ) 
    ## 
    if frame and report :
        report = frame.Report()
        logger.info ( '%s\n%s' % ( title , report_print ( report , title = title , prefix = '# ') ) )
    return values

# ==============================================================================
## Local helper method to generate new, unique name for the variable 
def var_name ( prefix , used_names , *garbage ) :
    """ Local helper method to generate new, unique name for the variable 
    """
    names = tuple ( n for n in used_names ) 
    name  =    prefix + '%x' % ( hash ( (        prefix , names , garbage ) ) % ( 2**32 ) ) 
    while name in used_names :
        name = prefix + '%x' % ( hash ( ( name , prefix , names , garbage ) ) % ( 2**32 ) )
    return name

# ===============================================================================
## Print the frame report data
def report_print_table ( report , title  = '' , prefix = '' , more_rows = [] ) :
    """ Print a frame report data 
    """
    from ostap.core.core import binomEff
    
    n0     = -1 
    lmax   =  5
    table  = []
    
    for name, passed, all in report :
        n0   = max ( n0 , all , passed )        
        eff1 = binomEff ( passed , all ) * 100        
        eff2 = binomEff ( passed ,  n0 ) * 100
        lmax = max ( len ( name ) , lmax , len ( 'Filter ' ) )     
        item = name ,  passed , all , eff1 , eff2 
        table.append ( item )
        
    lmax          =  max ( lmax + 2 , len ( 'Selection' ) + 2 )
    fmt_name      =  '%%-%ds ' % lmax 
    fmt_input     =  '%10d'
    fmt_passed    =  '%-10d'
    fmt_eff       =  '%8.3g +/- %-8.3g'
    fmt_cumulated =  '%8.3g +/- %-8.3g'
        
    header = ( ( '{:^%d}' % lmax ).format ( 'Filter'   ) ,               
               ( '{:>10}'        ).format ( '#input '  ) ,
               ( '{:<10}'        ).format ( '#passed'  ) ,
               ( '{:^20}'        ).format ( 'efficiency [%]' ) ,
               ( '{:^20}'        ).format ( 'cumulated efficiency [%]' ) )

    table_data = [ header ]
    for entry in table :
        n, p, a , e1 , e2 = entry
        table_data.append ( ( fmt_name      % n ,
                              fmt_input     % a ,
                              fmt_passed    % p  ,
                              fmt_eff       % ( e1.value () , e1.error () ) ,
                              fmt_cumulated % ( e2.value () , e2.error () ) ) )
    for row in more_rows :
        table_data.append ( row ) 
    
    import ostap.logger.table as T
    return T.table ( table_data , title = title , prefix = prefix , alignment = 'lcccc' )

# ===============================================================================
## Print the frame report
def report_as_table ( report ) :
    """ Print a frame report
    """
    table = []
    for c in report:
        name    = c.GetName ()
        passed  = c.GetPass ()
        all     = c.GetAll  ()
        table.append (  ( name , passed , all )  )

    return table 

# ===============================================================================
## Print the data frame report
def report_print ( report , title  = '' , prefix = '' , more_rows = [] ) :
    """ Print the data frame report
    """
    table = report_as_table ( report )    
    return report_print_table ( table , title, prefix , more_rows )

# ==============================================================================
## Is implicit MC globally enabled? 
mt_global = OCC.general.getboolean ( 'ImplicitMT' , fallback = True )    

# ==============================================================================
## Helper function that define expressions and cuts 
#  @code
#  frame = ...
#  currnt , vexpr , cexpr, input_string  = _fr_helper_ ( frame , 'x*x' , 'z<0' )
#  #endcode 
def _fr_helper_ ( frame , expressions , cuts = '' , progress = False ) :
    """ Helper function that define expressions and cuts 
    >>> frame = ...
    >>> current, vexpr , cexpr = _fr_helper_ ( frame , 'x*x' , 'z<0' )
    """
    
    if progress and isinstance ( frame , ROOT.TTree ) : progress = len ( frame )

    exprs, cuts, input_string = SV.vars_and_cuts ( expressions , cuts ) 
    
    ## Frame/Tree ?
    lenght = -1
    
    if   isinstance ( frame , ROOT.TTree  ) : node , length = DataFrame ( frame ) , len ( frame ) 
    elif isinstance ( frame , frame_types ) : node = frame 
    else                                    : node = as_rnode  ( frame ) 

    ## get the list of currently known names
    vars     = frame_columns ( node ) 
    all_vars = set ( vars ) 
    
    current  = node

    ## progress ? 
    if   progress is False               : pass
    elif progress is True and 0 < lenght : 
        current , cnt = frame_progress ( current , length ) 
    elif isinstance ( progress , integer_types ) and 1 < progress :
        current , cnt = frame_progress ( current , progress ) 
    elif progress :
        current , cnt = frame_progress ( current , lenght   ) 

    vnames   = ordered_dict ()
    added    = ordered_dict ()
    for i , expr in enumerate ( exprs , start = 1 ) :
        vname = expr
        if not expr in all_vars :
            used    = tuple ( all_vars | set ( frame_columns ( current ) ) ) 
            vn      = var_name ( 'var%d_' % i , used , expr , *vars )
            all_vars.add ( vn )
            current = current.Define ( vn , expr )
            vname   = vn
            added [ vn ] = expr 
        vnames [ expr ] = vname 

    cname = cuts
    if cuts and not cuts in all_vars :
        used    = tuple ( all_vars | set ( frame_columns ( current ) ) ) 
        cn      = var_name ( 'cut_' , used , cuts , *vars )
        all_vars.add ( cn )
        ncuts   = '1.0*(%s)' % cuts 
        current = current.Define ( cn , ncuts )
        cname   = cn
        added [ cn ] = ncuts 

    ## attach cuts as BOOLEAN filter 
    if cname :
        fname      = '(PRE)FILTER %s' % cuts
        the_filter = '(bool) %s' % cname 
        current    = current.Filter ( the_filter , fname ) 
        logger.debug ( 'Added pre-filter: %s' % the_filter ) 
                       
    if added :
        title = 'Added variables'
        rows  = [ ( 'Variable' , 'Expression' ) ] 
        for row in loop_items ( added ) : rows.append ( row )
        table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lr' )        
        logger.debug ( '%s:\n%s' % ( title , table ) )
        
    return current, vnames, cname, input_string

# ==================================================================================
## The second helper method to implement various "statistics"-related actions  
def _fr_helper2_ ( frame                                ,
                   creator                              ,
                   expressions                          , * ,  
                   cuts      = ''                       ,  
                   progress  = False                    ,
                   report    = False                    ,
                   lazy      = True                     , 
                   transform = lambda s : s.GetValue () ,  
                   title     = 'DataFrame'              ) : 
    """ The second helper method to implement various statistic-related actions  
    """

    current, var_names, cut_name, input_string = \
        _fr_helper_ ( frame , expressions , cuts , progress = progress )

    results = {}
    for expr, var_name in loop_items ( var_names ) : 
        results [ expr ] = creator ( current , var_name , cut_name ) 

    ## single variable as argument 
    if input_string and 1 == len ( results ) :
        _ , results = results.popitem ()

    ## Lazy processing ? 
    if lazy :
        return results , current

    ## get result/perform the loop 
    return get_values ( results               ,                        
                        frame     = current   ,
                        transform = transform , 
                        report    = report    ,
                        title     = title     ) 

# ==============================================================================
## modify constructor for RDataFrame to enable/disable Implicit multithreading
#  @code
#  f1 = DataFrame ( ... , enable = True  , progress = False ) ## default
#  f2 = DataFrame ( ... , enable = False , progress = False ) ## default
#  @endcode
#  @see ROOT::EnableImplicitMT 
#  @see ROOT::DisableImplicitMT
#  @see ROOT::IsImplicitMTEnabled
def _fr_new_init_ ( self , *args , **kwargs ) :
    """ Modify the DataFrame constuctor to allow (semi)automatic
    manipulations wth ROOT.ROOT.EnableImplicitMT/DisableImplicitMT
    - see ROOT.ROOT.EnableImplicitMT 
    - see ROOT.ROOT.DisableImplicitMT
    - see ROOT.ROOT.IsImplicitMTEnabled
    >>> f = DataFrame ( ... , enable = True  , progress = False ) ## default 
    >>> f = DataFrame ( ... , enable = False , progress = False )
    """
        
    mt        = kwargs.pop ( 'enable'   , mt_global ) and mt_global 
    progress  = kwargs.pop ( 'progress' , False     )
    
    if       mt and not ROOT.ROOT.IsImplicitMTEnabled() :
        ROOT.ROOT.EnableImplicitMT  ()
        logger.debug ( 'DataFrame:ImplicitMT is %s' % ( 'Enabled' if ROOT.ROOT.IsImplicitMTEnabled() else 'Disabled' ) )
    elif not mt and     ROOT.ROOT.IsImplicitMTEnabled() :
        ROOT.ROOT.DisableImplicitMT ()
        logger.info  ( 'DataFrame:ImplicitMT is %s' % ( 'Enabled' if ROOT.ROOT.IsImplicitMTEnabled() else 'Disabled' ) ) 
        
    self._fr_old_init_ ( *args , **kwargs ) 

    if progress and args :
        lenght = -1 
        if   isinstance ( args [ 0 ] , integer_types ) and 1 < args [ 0 ] : lenght = args [ 0 ]
        elif isinstance ( args [ 0 ] , ROOT.TTree    )                    : lenght = len ( args [ 0 ] )
        elif isinstance ( progress   , inetger_types ) and 1 < progress   : lenght = progress         
        _ , _ = frame_progress ( self , lenght ) 
        
    
if not hasattr ( DataFrame , '_fr_old_init_' ) :
    DataFrame._fr_old_init_ = DataFrame.__init__
    DataFrame.__init__      = _fr_new_init_

# =============================================================================
## get all column names
#  @code
#  cols = frame_colums   ( frame ) 
#  cols = frame_branches ( frame ) ## ditto
#  @endcode 
def frame_columns ( frame ) :
    """ Get all column names
    >>> cols = frame_colums   ( frame ) 
    >>> cols = frame_branches ( frame ) ## ditto 
    """
    names  = [ str ( c ) for c in frame.GetColumnNames()        ]
    names += [ str ( c ) for c in frame.GetDefinedColumnNames() ]            
    return tuple ( sorted ( set ( names ) ) )
    
# ==============================================================================
## Frame branches, same as frame columns
#  @see frame_columns 
frame_branches = frame_columns
# =============================================================================

# =============================================================================
## Display the progress bar for the DataFrame:
#  @code
#  f     = DataFrame ( ... )
#  f , c = frame_progress1 ( f , 1000 ) ## number of elements!
#  ....
#  @endcode 
#  @attention  active only for ROOT < 6.32
def frame_progress1 ( frame  , length ) :
    """ Draw (lazy) progress bar for the    DataFrame:
    >>> f      = DataFrame ( ... )
    >>> f , c  = frame_progress1 ( f , 1000 ) ## number of elements!
    >>> ....
    - active only for ROOT < 6.32
    """
    
    if   isinstance ( frame , frame_types ) : pass
    elif isinstance ( frame , ROOT.TTree  ) :
        length = len       ( frame )
        frame  = DataFrame ( frame ) 
    else                                    :
        frame = as_rnode   ( frame )
    
    cnt = frame.Count ()
    
    ## commented out
    if    not isatty () : return frame , cnt ## ATTENTION 
    elif  length <= 0   : return frame_progress2 ( frame , length ) 
        
    if   2048 < length  : nchunks = 2000 
    elif 1024 < length  : nchunks = 1000 
    elif  512 < length  : nchunks =  500 
    elif  101 < length  : nchunks =  100 
    else                : return frame, cnt     ## no progress bar for short frames 
    
    csize , rr = divmod ( length , nchunks )
    csize      = max    ( csize  , 1       )
    
    if rr : nchunks += 1 

    if ( 6 , 32 ) <= root_info :
        cnt = Ostap.Utils.add_progress_bar ( cnt , nchunks , csize , progress_conf () ) 
    else : 
        fun = Ostap.Utils.frame_progress ( nchunks , progress_conf () )
        cnt = cnt.OnPartialResultSlot  ( csize , fun )
        
    return frame , cnt

# =========================================================================
## make use of new `ROOT::RDF::Experimental::AddProgressbar` utility
#  @see ROOT::RDF::Experimental::AddProgressbar 
#  @attention no action for ROOT < 6.30 
def frame_progress2 ( frame , length = -1 ) :
    """ Make use of new `ROOT.RDF.Experimental.AddProgressbar` utility
    - see ROOT.RDF.Experimental.AddProgressbar
    - no action for ROOT < 6.30 
    """
    
    if   isinstance ( frame , frame_types ) : pass
    elif isinstance ( frame , ROOT.TTree  ) :
        length = len       ( frame )
        frame  = DataFrame ( frame ) 
    else                                    :
        frame = as_rnode   ( frame )
        
    cnt = frame.Count ()
    
    if    not isatty ()          : return frame , cnt  ## ATTENTION 
    elif  root_info < ( 6 , 30 ) : return frame , cnt  ## ATTENTION! 

    ## add progress bar 
    ROOT.ROOT.RDF.Experimental.AddProgressBar ( frame )    
    return frame, cnt 

# =========================================================================
## Add progress bar to the frame
#  @code
#  frame = ...
#  frame , cnt = frame_progrees ( frame , nevents ) 
#  @endcode 
def frame_progress ( frame , length = -1 ) :
    """ Add progress bar to the frame
    >>> frame = ...
    >>> frame , cnt = frame_progrees ( frame , nevents ) 
    """
    
    if   isinstance ( frame , frame_types ) : pass
    elif isinstance ( frame , ROOT.TTree  ) :
        length = len       ( frame )
        frame  = DataFrame ( frame ) 
    else                                    :
        frame = as_rnode   ( frame )
        
    if 1 <= length : return frame_progress1 ( frame , length )
    else           : return frame_progress2 ( frame , length )
    
# ==============================================================================
## Prescale the frame
#  @code
#  frame = ...
#  prescaled1 = frame_prescale ( frame , 0.15 ) ##  0.0 < prescale < 1.0
#  prescaled2 = frame_prescale ( frame , 20   ) ##  1   < prescale 
#  @endcode 
def frame_prescale ( frame , prescale , name = '' ) :
    """ Prescale the frame
    >>> frame = ...
    >>> prescaled1 = frame_prescale ( frame , 0.15 ) ##  0.0 < prescale < 1.0
    >>> prescaled2 = frame_prescale ( frame , 20   ) ##  1   < prescale 
    """
    
    if isinstance   ( frame  , ROOT.TTree ) : node = DataFrame ( frame )
    else                                    : node = as_rnode  ( frame )
    
    if isinstance ( prescale , integer_types ) and 1 < prescale :

        name = name if name else 'PRESCALE_%d' % prescale        
        code = '( 0 == ( ( rdfentry_ + %d * rdfslot_ ) %% %d ) )'
        ## 16777213 and 16777199 are just large prime numbers 
        code = code % ( 16777213 , prescale )
        return node.Filter ( code , name ) 
        
    elif isinstance ( prescale , float ) and 0 < prescale < 1 :

        name = name if name else 'PRESCALE_%.6g' % prescale        
        code = 'gRandom->Rndm() <= %.12g'        % prescale
        return node.Filter ( code , name )
    
    elif 1 == prescale : return node 

    raise TypeError ( "Invalid type/value for 'prescale' %s/%s" %( prescale , type ( prescale ) ) )

# =============================================================================
## Simplified print out for the frame 
#  @code 
#  frame = ...
#  @endcode 
def frame_print ( frame ) :
    """ Simplified print out for the  frame 
    
    >>> frame = ...
    >>> print frame
    """
    ## 
    if isinstance   ( frame  , ROOT.TTree ) : frame = DataFrame ( frame )
    ## 
    node = as_rnode ( frame ) 
    res  = "DataFrame Enries/#%d" %  len ( frame )  
    ##
    cols = frame_columns ( node ) 
    res += "\nColumns:\n%s" % multicolumn ( cols , indent = 2 , pad = 1 )
    return res

# ==============================================================================
## Data frame as table
#  @code
#  frame = ...
#  table = frame_table ( frame , '.*PT.*' , cuts = ... , more_vars = [ 'x*x/y' , 'y+z'] )
#  print ( table ) 
#  @endcode 
def frame_table ( frame , pattern = None , cuts = '' , more_vars = () , title = '' ,  prefix = '' ) :
    """ Data frame as table
    >>> frame = ...
    >>> table = frame_table ( frame , '.*PT.*' , cuts = ... , more_vars = [ 'x*x/y' , 'y+z'] )
    >>> print ( table )
    """
    
    if isinstance ( frame  , ROOT.TTree ) : frame = DataFrame ( frame  )
    node = as_rnode ( frame )
    
    ## the basic column type 
    def col_type ( var ) :
        """ The basic column type"""
        if var in cols : t = frame.GetColumnType ( var )
        else           : return ''
        
        if   'Double_t'  == t : return 'double'
        elif 'Float_t'   == t : return 'float'
        elif 'Bool_t'    == t : return 'bool'
        elif 'Char_t'    == t : return 'char'
        elif 'Short_t'   == t : return 'short'
        elif 'Int_t'     == t : return 'int'
        elif 'Long64_t'  == t : return 'long'
        elif 'UChar_t'   == t : return 'unsigned char'
        elif 'UShort_t'  == t : return 'unsigned short'
        elif 'UInt_t'    == t : return 'unsigned int'
        elif 'ULong64_t' == t : return 'unisgned long'
        
        return t 
        
    ## get all variables/columns 
    cols  = frame_columns ( node  )
    
    vcols = cols
    
    if pattern :
        import re
        cpat  = re.compile ( pattern )
        vcols = tuple ( [ v for v in vcols if cpat.match ( v ) ] ) 

    svars = vcols + tuple   ( more_vars ) 
    stats = frame_statistic ( node , svars , cuts = cuts , lazy = False ) 

    if   len ( vcols ) < 10    : ifmt = '%1d.'
    elif len ( vcols ) < 100   : ifmt = '%2d.'
    elif len ( vcols ) < 1000  : ifmt = '%3d.'
    else                       : ifmt = '%d6.'

    header =  ( '#' , 'Variable' , 'Type' , 'mean' , 'RMS' , 'min' , 'max' ) 
    table  = [ header ]
    
    for i , e in enumerate ( sorted ( vcols ) , start = 1 ) :

        stat = stats [ e ]
               
        n    = stat.nEntries() 
        mnmx = stat.minmax ()
        mean = stat.mean   () 
        rms  = stat.rms    ()

        row = [ ifmt % i , e , col_type ( e ) ,
                ( '%+.5g' % mean.value() ).strip()                  , ## 4
                ( '%-.5g'  % rms         ).strip()                  , ## 5 
                ( '%+.5g' % mnmx[0]      ).strip()                  , ## 6
                ( '%-+.5g' % mnmx[1]     ).strip()                  , ## 7
                ]

        table.append ( row ) 

    nv = len ( table )
    
    if   len ( more_vars ) < 10    : ifmt = '%1d.'
    elif len ( more_vars ) < 100   : ifmt = '%2d.'
    elif len ( more_vars ) < 1000  : ifmt = '%3d.'
    else                           : ifmt = '%6d.'

    for i , e in enumerate ( sorted ( more_vars ) , start = 1 ) :

        stat = stats [ e ]
               
        n    = stat.nEntries() 
        mnmx = stat.minmax ()
        mean = stat.mean   () 
        rms  = stat.rms    ()

        row = [ ifmt %  i , e , col_type ( e ) ,
                ( '%+.5g' % mean.value() ).strip()                  , ## 4
                ( '%-.5g'  % rms         ).strip()                  , ## 5 
                ( '%+.5g' % mnmx[0]      ).strip()                  , ## 6
                ( '%-+.5g' % mnmx[1]     ).strip()                  , ## 7
                ]

        table.append ( row ) 
        
    import ostap.logger.table as T
    title  = title if title  else 'DataFrame variables'
    if pattern : title = '%s (pattern %s)' % ( title , pattern )
    t  = T.table (  table , title , prefix = prefix , alignment = 'rllrlrl' )
    
    return t 

# =============================================================================
## Get the length/size of the data frame
#  @code
#  frame = ...
#  print ( len(frame) )
#  len   = frame_length ( frame ) 
#  len   = frame_size   ( frame ) ## ditto 
#  @endcode 
def frame_length ( frame            , 
                   cuts     = ''    , * , 
                   progress = False ,
                   report   = False ,
                   lazy     = False ) :
    """ Get the length/size of the data frame
    >>> frame = ...
    >>> print len(frame)
    >>> len   = frame_length ( frame ) 
    >>> len   = frame_size   ( frame ) ## ditto 
    """
    
    current, _ , cut_name , _ = \
        _fr_helper_ ( frame , '1' , cuts , progress = progress )
    
    ## create counter 
    cnt  = current.Count ()

    if lazy :
        return cnt , current                           ## RETURN 
    ##
    return get_values ( cnt                     ,
                        frame  = current        ,
                        report = report         ,
                        title  = 'frame_length' )

# =================================================================================
## number of "good" etnties 
frame_size = frame_length
# =================================================================================

# ==================================================================================
## Generic counters
# ==================================================================================

# ==================================================================================
## get statistics of variable(s) in a form of moments 
#  @code
#  frame = ....
#  stat  = frame_statistic ( 'pt'           )
#  stat  = frame_statistic ( 'pt' , 'eta>0' )
#  @endcode
#  @see Ostap::Math::StatEntity 
#  @see Ostap::Math::WStatEntity 
def frame_statistic ( frame               , 
                      expressions         , 
                      cuts        = ''    , * , 
                      as_weight   = True  , 
                      progress    = False ,
                      report      = False ,
                      lazy        = False ) :
    """ Get statistics of variable(s)
    >>> frame = ....
    >>> stat  = frame_statistic ( frame , 'pt' )
    """
    current, vname , cname , input_string = \
        _fr_helper_ ( frame , expressions , cuts ) 

    ## ATTENTION HERE!!
    if cname and not as_weight :
        logger.warning ( "The cut is treated as boolean: %s" % cuts ) 
        cname = ''

    def screator ( node , var_name , cut_name ) :
        if cut_name : 
            TT = SA1w [ Ostap.WStatEntity ] 
            return node.Book ( std_move ( TT () ) , CNT ( [ var_name , cut_name ] ) ) 
        else :
            TT = SA1  [ Ostap.StatEntity  ] 
            return node.Book ( std_move ( TT () ) , CNT ( 1 , var_name ) ) 
        
    return _fr_helper2_ ( current                         , 
                          screator                        ,
                          expressions                     ,
                          cuts        = cname             ,
                          progress    = progress          ,
                          report      = report            ,
                          lazy        = lazy              ,
                          title       = 'frame_statistic' )

# ==================================================================================
## get min/max for variable 
#  @code
#  frame = ....
#  stat  = frame_minmax  ( frame , 'pt' )
#  @endcode
#  @see Ostap::Math::StatEntity 
#  @see Ostap::Math::WStatEntity 
def frame_minmax ( frame               , 
                   expressions         , 
                   cuts        = ''    , * , 
                   as_weight   = True  , 
                   progress    = False ,
                   report      = False ,
                   lazy        = False ) :
    """ Get min/max variable(s)
    >>> frame = ....
    >>> stat  = frame_minmax ( frame , 'pt' )
    """

    ## get counters 
    results = frame_statistic ( frame                    ,
                                expressions              ,
                                cuts         = cuts      ,
                                as_weight    = as_weight ,
                                progress     = progress  ,
                                report       = report    ,
                                lazy         = lazy      )

    
    ##
    return results if lazy else \
        get_values ( results , transform = lambda s : ( s.min(), s.max() ) )


# ==================================================================================
## get suitable range(s) for variables 
#  @code
#  frame = ....
#  stat  = frame_range ( frame , 'pt' )
#  @endcode
def frame_range ( frame               , 
                  expressions         , 
                  cuts        = ''    , * , 
                  as_weight   = True  , 
                  progress    = False ,
                  report      = False ,
                  delta       = 0.01  ) : 
    """ Get suitable ranges of variable(s)
    >>> frame = ....
    >>> stat  = frame_range ( frame , 'pt' )
    """

    ## get minmax (non LAZY!)
    result = frame_minmax ( frame                    ,
                            expressions              ,
                            cuts         = cuts      ,
                            as_weight    = as_weight ,
                            progress     = progress  ,
                            report       = report    ,
                            lazy         = False     )
    
    if isinstance ( result , dictlike_types ) :
        ranges = {}
        for key , value in loop_items ( result ) :
            mn , mx = value 
            if mx <= mn and cuts :
                ## REDO WITH NO CUTS 
                mn , mx = frame_minmax ( frame , key , lazy = False )
            ## ATTENTION! 
            if mx <= mn : return None            
            ranges [ key ] = axis_range ( mn , mx , delta = delta )            
        return ranges 

    mn , mx = result
    if mx <= mn and cuts :
        ## REDO WITH NO CUTS 
        mn , mx = frame_minmax ( frame , expressions , lazy = False )
    ## ATTENTION 
    if mx <= mn : return None         
    return axis_range ( mn , mx , delta = delta )            

# ==================================================================================
# specfic counters 
# ==================================================================================

# ==================================================================================
## get nEff through the action
#  @code
#  frame = ...
#  stat  = frame_nEff ( dramw , cuts = 'y>0' )  
#  @endcode 
def frame_nEff ( frame              ,  
                 cuts      = ''     , * , 
                 as_weight = True   , 
                 progress  = False  ,
                 report    = False  ,
                 lazy      = False  ) : 
    """ Get nEff through action
    >>> frame = ...
    >>> nEff = frame_the_nEff ( 'x*x' , 'y>0' )
    """
    results = frame_statistic ( frame                 ,
                               "1"                   ,
                               cuts      = cuts      ,
                               as_weight = as_weight ,
                               progress  = progress  ,
                               report    = report    ,
                               lazy      = lazy      )
    
    return results if lazy else \
        get_values ( results , transform = lambda s : s.nEff()  )

# ==================================================================================
## get mean-values through the action
#  @code
#  frame = ...
#  stat  = frame_mean ( frame , 'x, t, z ' , cuts = 'y>0' )  
#  @endcode 
def frame_mean ( frame              ,
                 expressions        , 
                 cuts      = ''     , * ,
                 as_weight = True   , 
                 progress  = False  ,
                 report    = False  ,
                 lazy      = False  ) : 
    """ Get mean through action
    >>> frame = ...
    >>> stat  = frame_mean  ( 'x*x' , 'y>0' )
    """
    results = frame_statistic ( frame                 ,
                                expressions           ,
                                cuts      = cuts      ,
                                as_weight = as_weight ,
                                progress  = progress  ,
                                report    = report    ,
                                lazy      = lazy     )
    
    return results if lazy else \
        get_values ( results , transform = lambda s : s.mean  () )

# ==================================================================================
## get statistics of variable(s) in a form of moments 
#  @code
#  frame = ....
#  stat  = frame_the_moment ( 5 , 'pt'           )
#  stat  = frame_the_moment ( 5 , 'pt' , 'eta>0' )
#  @endcode
#  @see Ostap::Math::Moment_
#  @see Ostap::Math::WMoment_
def frame_the_moment ( frame , N           ,
                       expressions         , 
                       cuts        = ''    , * , 
                       as_weight   = True  , 
                       progress    = False ,
                       report      = False ,
                       lazy        = False ) :
    """ Get statistics of variable(s)
    >>> frame = ....
    >>> stat  = frame_the_moment ( frame , 5 , 'pt'                  )
    >>> stat  = frame_the_moment ( frame , 5 , 'pt' , cuts = 'eta>0' )
    """
    assert isinstance ( N , integer_types ) and 0 <= N , 'Invalid order!'

    current, vname , cname , input_string = _fr_helper_ ( frame , expressions , cuts ) 
    
    ## ATTENTION HERE!!
    if cname and not as_weight :
        logger.warning ( "The cut is treated as boolean: %s" % cuts ) 
        cname = ''

    def mcreator ( node , var_name , cut_name ) :
        
        if cut_name : 
            TT = SA1w [ Ostap.Math.WMoment_[N] ] 
            return node.Book ( std_move ( TT () ) , CNT ( [ var_name , cut_name ] ) ) 
        else :
            TT = SA1  [ Ostap.Math.Moment_[N] ] 
            return node.Book ( std_move ( TT () ) , CNT ( 1 , var_name ) ) 
        
    return _fr_helper2_ ( current                 , 
                          mcreator                ,
                          expressions             ,
                          cuts        = cname     ,
                          progress    = progress  ,
                          report      = report    ,
                          lazy        = lazy      ,
                          title       = 'frame_the_moment' )

# ==================================================================================
## get -values through the action
#  @code
#  frame = ...
#  stat  = frame_variance   ( frame , 'x, t, z ' , cuts = 'y>0' )  
#  stat  = frame_dispersion ( frame , 'x, t, z ' , cuts = 'y>0' )  
#  @endcode 
def frame_variance ( frame              ,
                     expressions        , 
                     cuts      = ''     , * , 
                     as_weight = True   , 
                     progress  = False  ,
                     report    = False  ,
                     lazy      = False  ) : 
    """ Get variance through action
    >>> frame = ...
    >>> stat  = frame_variance  ( 'x*x' , cuts = 'y>0' )
    >>> stat  = frame_dispersion ( 'x*x' , cuts = 'y>0' )
    """
    results = frame_the_moment ( frame                 ,
                                 4                     , ## neeeded to ge the uncertainty
                                 expressions           ,
                                 cuts      = cuts      ,
                                 as_weight = as_weight ,
                                 progress  = progress  ,
                                 report    = report    ,
                                 lazy      = lazy      )
    
    return results if lazy else \
        get_values ( results , transform = lambda s : s.variance () )

# ==================================================================================
## get -values through the action
#  @code
#  frame = ...
#  stat  = frame_variance   ( frame , 'x, t, z ' , cuts = 'y>0' )  
#  stat  = frame_dispersion ( frame , 'x, t, z ' , cuts = 'y>0' )  
#  @endcode 
def frame_rms ( frame              ,
                expressions        , 
                cuts      = ''     , * , 
                as_weight = True   , 
                progress  = False  ,
                report    = False  ,
                lazy      = False  ) : 
    """ Get rms through action
    >>> frame = ...
    >>> stat  = frame_ems   ( 'x*x' , cuts = 'y>0' )
    """
    results = frame_the_moment ( frame                 ,
                                 4                     , ## neeeded to ge the uncertainty
                                 expressions           ,
                                 cuts      = cuts      ,
                                 as_weight = as_weight ,
                                 progress  = progress  ,
                                 report    = report    ,
                                 lazy      = lazy      )
    
    return results if lazy else \
        get_values ( results , transform = lambda s : s.rms () )

# ==================================================================================
## get skewness through the action
#  @code
#  frame = ...
#  stat  = frame_skewness  ( frame , 'x, t, z ' , cuts = 'y>0' )  
#  @endcode 
def frame_skewness ( frame              ,
                     expressions        ,  
                     cuts      = ''     , * , 
                     as_weight = True   , 
                     progress  = False  ,
                     report    = False  ,
                     lazy      = False  ) : 
    """ Get skewness through action
    >>> frame = ...
    >>> stat  = frame_skewnexx   ( 'x*x' , cuts = 'y>0' )
    """
    results = frame_the_moment ( frame                 ,
                                 6                     , ## neeeded to ge the uncertainty
                                 expressions           ,
                                 cuts      = cuts      ,
                                 as_weight = as_weight ,
                                 progress  = progress  ,
                                 report    = report    ,
                                 lazy      = lazy      )
    
    return results if lazy else \
        get_values ( results , transform = lambda s : s.skewness () ) 
    
# ==================================================================================
## get kurtosis through the action
#  @code
#  frame = ...
#  stat  = frame_kurtosis  ( frame , 'x, t, z ' , cuts = 'y>0' )  
#  @endcode 
def frame_kurtosis ( frame              ,
                     expressions        , 
                     cuts      = ''     , * , 
                     as_weight = True   , 
                     progress  = False  ,
                     report    = False  ,
                     lazy      = False  ) : 
    """ Get kurtosis through action
    >>> frame = ...
    >>> stat  = frame_kurtosis ( 'x*x' , cuts = 'y>0' )
    """
    results = frame_the_moment ( frame                 ,
                                 8                     , ## neeeded to ge the uncertainty
                                 expressions           ,
                                 cuts      = cuts      ,
                                 as_weight = as_weight ,
                                 progress  = progress  ,
                                 report    = report    ,
                                 lazy      = lazy      )
    
    return results if lazy else \
        get_values ( results , transform = lambda s : s.kurtosis () )

## dispersion 
frame_dispersion = frame_variance


# =============================================================================
## get the statistic for pair of expressions in DataFrame
#  @code
#  frame  = ...
#  sta = frame_covariance ( 'x' , 'y' )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @see   Ostap::Math::Covariance
#  @see   Ostap::Math::WCovariance
#  @date   2018-06-18
def frame_covariance ( frame              ,
                       expression1        ,
                       expression2        ,  
                       cuts       = ''    , * , 
                       as_weight  = True  ,
                       progress   = False ,
                       lazy       = False ) :
    """ Get the statistic for pair of expressions in DataFrame
    
    >>>  frame  = ...
    >>>  cov = framw.covariance(  frame , 'x' , 'y' )

    """
    ## decode expressions & cuts 
    current , items, cname , _ = \
        _fr_helper_ ( frame , [ expression1 , expression2 ] , cuts , progress = progress )
    
    assert 2 == len ( items ) , "Imvalid expressions!"
   
    ## ATTENTION HERE!!
    if cname and not as_weight :
        logger.warning ( "The cut is treated as boolean: %s" % cuts ) 
        cname = ''

    vars = [ v for v in items.values ]
    if cname : vars.append ( cname ) 
    
    if cname : ACTION = SA2w [ Ostap.Math.WCovarinace ] 
    else     : ACTION = SA2  [ Ostap.Math. Covarinace ] 
        
    results = current.Book ( std_move ( ACTION () ) , CNT ( vars ) )

    if lazy :
        return results , current  ## RETURN 
    
    return get_values ( results ,
                        frame  = current ,
                        report = report  ,
                        title  = 'frame_covariance' ) 

# ==============================================================================
## Get an arithmetic mean
#  @see Ostap::Math::ArithmeticMean
#  @code
#  frame = ...
#  mean = frame_the_arithmetic_mean ( frame , 'x*x' , '0<y' ) 
#  @endcode 
def frame_arithmetic_mean ( frame              ,
                            expressions        , 
                            cuts       = ''    , * , 
                            as_weight  = True  , 
                            progress   = False ,
                            report     = False ,
                            lazy       = False ) :
    """ Get an arithmetic mean
    >>> frame = ...
    >>> mean = frame_the_arithmetic_mean ( frame , 'x*x' , '0<y' ) 
    -  see `Ostap.Math.ArithmeticMean`
    """
    
    ## decode expressions & cuts 
    current , items, cname , _ = \
        _fr_helper_ ( frame , expressions , cuts , progress = progress )

    ## ATTENTION HERE!!
    if cname and not as_weight :
        logger.warning ( "The cut is treated as boolean: %s" % cuts ) 
        cname = ''
   
    def acreator ( node , var_name , cut_name ) : 
        if cut_name: 
            TT = SA1w [ Ostap.Math.WArithmeticMean  ] 
            return node.Book ( std_move ( TT () ) , CNT ( [ var_name , cut_name ] ) )
        else  :
            TT = SA1  [ Ostap.Math. ArithmeticMean  ] 
            return node.Book ( std_move ( TT () ) , CNT ( 1 , var_name ) )     
        
    return _fr_helper2_ ( frame               ,
                          acreator            , 
                          expressions         ,
                          cuts     = cuts     ,
                          progress = progress ,
                          report   = report   ,
                          lazy     = lazy     ,
                          title    = 'frame_arithmeric_mean' , 
                          transform = lambda s : s.value ()  )       

# ==============================================================================
## Get harmonic mean
#  @see Ostap::Math::HarmonicMean
#  @see Ostap::Math::WHarmonicMean
#  @code
#  frame = ...
#  mean = frame_harmoni_mean ( frame , 'x*x' , cut = '0<y' ) 
#  @endcode 
def frame_harmonic_mean ( frame              ,
                          expressions        , 
                          cuts       = ''    , * , 
                          as_weight  = True  , 
                          progress   = False ,
                          report     = False ,
                          lazy       = False ) :
    """ Get harmnonic mean
    >>> frame = ...
    >>> mean = frame_harmonic_mean ( frame , 'x*x' , '0<y' ) 
    -  see `Ostap.Math.HarmonicMean`
    -  see `Ostap.Math.WHarmonicMean`
    """
    
    ## decode expressions & cuts 
    current , items, cname , _ = \
        _fr_helper_ ( frame , expressions , cuts , progress = progress )

    ## ATTENTION HERE!!
    if cname and not as_weight :
        logger.warning ( "The cut is treated as boolean: %s" % cuts ) 
        cname = ''
   
    def acreator ( node , var_name , cut_name ) : 
        if cut_name: 
            TT = SA1w [ Ostap.Math.WHarmonicMean  ] 
            return node.Book ( std_move ( TT () ) , CNT ( [ var_name , cut_name ] ) )
        else  :
            TT = SA1  [ Ostap.Math. HarmonicMean  ] 
            return node.Book ( std_move ( TT () ) , CNT ( 1 , var_name ) )     
        
    return _fr_helper2_ ( frame               ,
                          acreator            , 
                          expressions         ,
                          cuts     = cuts     ,
                          progress = progress ,
                          report   = report   ,
                          lazy     = lazy     ,
                          title    = 'frame_harmonic_mean'  , 
                          transform = lambda s : s.value () )       

# ==============================================================================
## Get geometric mean
#  @see Ostap::Math::GeometricMean
#  @see Ostap::Math::GeometricMean
#  @code
#  frame = ...
#  mean = frame_harmoni_mean ( frame , 'x*x' , cut = '0<y' ) 
#  @endcode 
def frame_geometric_mean ( frame              ,
                           expressions        , 
                           cuts       = ''    , * , 
                           as_weight  = True  , 
                           progress   = False ,
                           report     = False ,
                           lazy       = False ) :
    """ Get an geometric mean
    >>> frame = ...
    >>> mean = frame_geonmetric_mean ( frame , 'x*x' , '0<y' ) 
    -  see `Ostap.Math.GeometricMean`
    -  see `Ostap.Math.WGeometricMean`
    """
    
    ## decode expressions & cuts 
    current , items, cname , _ = \
        _fr_helper_ ( frame , expressions , cuts , progress = progress )

    ## ATTENTION HERE!!
    if cname and not as_weight :
        logger.warning ( "The cut is treated as boolean: %s" % cuts ) 
        cname = ''
   
    def acreator ( node , var_name , cut_name ) : 
        if cut_name: 
            TT = SA1w [ Ostap.Math.WGeometricMean  ] 
            return node.Book ( std_move ( TT () ) , CNT ( [ var_name , cut_name ] ) )
        else  :
            TT = SA1  [ Ostap.Math. GeometricMean  ] 
            return node.Book ( std_move ( TT () ) , CNT ( 1 , var_name ) )     
        
    return _fr_helper2_ ( frame               ,
                          acreator            , 
                          expressions         ,
                          cuts     = cuts     ,
                          progress = progress ,
                          report   = report   ,
                          lazy     = lazy     ,
                          title    = 'frame_geometric_mean' , 
                          transform = lambda s : s.value () ) 

# ==============================================================================
## Get power mean
#  @see Ostap::Math::PowerMean
#  @see Ostap::Math::PowerMean
#  @code
#  frame = ...
#  mean = frame_power_mean ( frame , 0.5 , 'x*x' , cut = '0<y' ) 
#  @endcode 
def frame_power_mean ( frame      ,   p   ,
                       expressions        , 
                       cuts       = ''    , * , 
                       as_weight  = True  , 
                       progress   = False ,
                       report     = False ,
                       lazy       = False ) :
    """ Get power mean
    >>> frame = ...
    >>> mean = frame_power_mean ( frame , 0.5 , 'x*x' , '0<y' ) 
    -  see `Ostap.Math.PowerMean`
    -  see `Ostap.Math.WPowerMean`
    """
    assert isinstance ( p , num_types ) , 'Invalid type of p-parameter: %s' % typename ( p ) 

    args = { 'cuts'      : cuts      ,
             'as_weight' : as_weight ,
             'progress'  : progress  ,
             'report'    : report    ,
             'lazy'      : lazy      }
    if    -1 == p or isequal ( p , -1.0 ) : 
        return frame_harmonic_mean   ( frame , expressions , **args )
    elif   0 == p or iszero  ( p        ) : 
        return frame_geometric_mean  ( frame , expressions , **args )
    elif   1 == p or isequal ( p ,  1.0 ) : 
        return frame_arithmetic_mean ( frame , expressions , **args )
    
    ## decode expressions & cuts 
    current , items, cname , _ = \
        _fr_helper_ ( frame , expressions , cuts , progress = progress )

    ## ATTENTION HERE!!
    if cname and not as_weight :
        logger.warning ( "The cut is treated as boolean: %s" % cuts ) 
        cname = ''
   
    def acreator ( node , var_name , cut_name ) : 
        if cut_name:            
            TT = SA1w [ Ostap.Math.WPowerMean  ]
            PP = Ostap.Math.WPowerMean ( p ) 
            return node.Book ( std_move ( TT ( PP ) ) , CNT ( [ var_name , cut_name ] ) )
        else  :
            TT = SA1  [ Ostap.Math. PowerMean  ]
            PP = Ostap.Math. PowerMean ( p )             
            return node.Book ( std_move ( TT ( PP ) ) , CNT ( 1 , var_name ) )     
        
    return _fr_helper2_ ( frame                ,
                          acreator             , 
                          expressions          ,
                          cuts      = cuts     ,
                          progress  = progress ,
                          report    = report   ,
                          lazy      = lazy     ,
                          title     = 'frame_power_mean'    , 
                          transform = lambda s : s.value () ) 

# ==============================================================================
## Get power mean
#  @see Ostap::Math::LehmerMean
#  @see Ostap::Math::WLehmerMean
#  @code
#  frame = ...
#  mean = frame_lehmer_mean ( frame , 0.5 , 'x*x' , cut = '0<y' ) 
#  @endcode 
def frame_lehmer_mean ( frame      , p     ,
                        expressions        , 
                        cuts       = ''    , * , 
                        as_weight  = True  , 
                        progress   = False ,
                        report     = False ,
                        lazy       = False ) :
    """ Get Lehmer  mean
    >>> frame = ...
    >>> mean = frame_lehmer_mean ( frame , 0.5 , 'x*x' , '0<y' ) 
    -  see `Ostap.Math.LehmerMean`
    -  see `Ostap.Math.WLehmerMean`
    """
    assert isinstance ( p , num_types ) , 'Invalid type of p-parameter: %s' % typename ( p ) 

    args = { 'cuts'      : cuts      ,
             'as_weight' : as_weight ,
             'progress'  : progress  ,
             'report'    : report    ,
             'lazy'      : lazy      }
    
    if   0   == p or iszero  ( p        ) : 
        return frame_harmonic_mean   ( frame , expressions , **args )
    elif 1   == p or isequal ( p ,  1.0 ) : 
        return frame_arithmetic_mean ( frame , expressions , **args )
    elif 0.5 == p or isequal ( p , 0.5 ) : 
        return frame_geometric_mean  ( frame , expressions , **args )
    
    ## decode expressions & cuts 
    current , items, cname , _ = \
        _fr_helper_ ( frame , expressions , cuts , progress = progress )

    ## ATTENTION HERE!!
    if cname and not as_weight :
        logger.warning ( "The cut is treated as boolean: %s" % cuts ) 
        cname = ''
   
    def acreator ( node , var_name , cut_name ) : 
        if cut_name:            
            TT = SA1w [ Ostap.Math.WLehmerMean  ]
            PP =       Ostap.Math.WLehmerMean ( p ) 
            return node.Book ( std_move ( TT ( PP ) ) , CNT ( [ var_name , cut_name ] ) )
        else  :
            TT = SA1  [ Ostap.Math. LehmerMean  ]
            PP =        Ostap.Math. LehmerMean ( p )             
            return node.Book ( std_move ( TT ( PP ) ) , CNT ( 1 , var_name ) )     
        
    return _fr_helper2_ ( frame                            ,
                          acreator                         , 
                          expressions                      ,
                          cuts      = cuts                 ,
                          progress  = progress             ,
                          report    = report               ,
                          lazy      = lazy                 ,
                          title     = 'frame_lehmer_mean'  , 
                          transform = lambda s : s.value () ) 


# =============================================================================
## Get empirical CDF 
#  @code
#  data = ...
#  c1 = frame_ECDF ( data , 'S_sw' ) 
#  c2 = frame_ECDF ( data , S_sw' , 'pt>0'  )
#  @endcode
def frame_ECDF ( frame               ,
                 expressions         , 
                 cuts        = ''    , * , 
                 as_weight   = True  , 
                 progress    = False ,
                 report      = False ,
                 lazy        = False ) :
    """ Get empirical CDF 
    >>> data = ...
    >>> c1 = frame_ECDF ( data , 'S_sw' ) 
    >>> c2 = frame_ECDF ( data , S_sw' , 'pt>0'  )
    """
    current, vname , cname , input_string = _fr_helper_ ( frame , expressions , cuts ) 
    
    ## ATTENTION HERE!!
    if cname and not as_weight :
        logger.warning ( "The cut is treated as boolean: %s" % cuts ) 
        cname = ''
        
    def ecreator ( node , var_name , cut_name ) :
        
        if cut_name : 
            TT = SA1w [ Ostap.Math.WECDF ]
            return node.Book ( std_move ( TT () ) , CNT ( [ var_name , cut_name ] ) ) 
        else :
            TT = SA1  [ Ostap.Math. ECDF ]
            return node.Book ( std_move ( TT () ) , CNT ( 1 , var_name ) ) 
        
    return _fr_helper2_ ( current                 , 
                          ecreator                ,
                          expressions             ,
                          cuts        = cname     ,
                          progress    = progress  ,
                          report      = report    ,
                          lazy        = lazy      , 
                          title       = 'frame_ECDF' )

# ==============================================================================
_types_1D =  Ostap.Math.LegendreSum  , Ostap.Math.Bernstein   , Ostap.Math.ChebyshevSum , 
_types_2D =  Ostap.Math.LegendreSum2 , Ostap.Math.Bernstein2D ,
_types_3D =  Ostap.Math.LegendreSum3 , Ostap.Math.Bernstein3D , 
_types_4D =  Ostap.Math.LegendreSum4 ,
_types_nD = _types_1D + _types_2D + _types_3D + _types_4D 
# ============================================================================
## Project of the frame
#  @code
#  frame    = ...
#  h1_model = ...
#  h1       = frame_project ( frame , h1_model , 'pt' )
#  h2_model = ...
#  h2       = frame_project ( frame , h2_model ,  ( 'pt' , 'x' ) )
#  h3_model = ...
#  h3       = frame_project ( frame , h3_model , ( 'pt' , 'x'  , 'y' ) )
#  ...
#  @endcode
#  Cuts/weigth are also can be specified
#  @code
#  frame    = ...
#  h1_model = ...
#  h1       = frame_project ( frame , h1_model , 'pt' , 'w')
#  h2_model = ...
#  h2       = frame_project ( frame , h2_model , 'pt,x' , 'w')
#  h3_model = ...
#  h3       = frame_project ( frame , h3_model , 'pt,x:y' , 'w' )
#  ...
#  @endcode
# 
#  If model is a real 1D/2D/3D-hstogram, action is *not* lazy, otherwise action *is* lazy
# 
#  Model can be a real histogram or 
#  - <code>ROOT::RDF::TH1DModel</code>
#  - <code>ROOT::RDF::TH2DModel</code>
#  - <code>ROOT::RDF::TH3DModel</code>
#  - anything that can be converted to <code>ROOT::RDF::TH1DModel</code>, 
#    <code>ROOT::RDF::TH2DModel</code> or <code>ROOT::RDF::TH3DModel</code> objects 
def frame_project ( frame               , 
                    model               ,
                    expressions         , 
                    cuts        = ''    , * , 
                    as_weight   = True  , 
                    progress    = False , 
                    report      = False ,
                    lazy        = False ) :
    """ Project of the frame into histogram (or other accumulator) 

    >>> frame    = ...
    >>> h1_model = ...
    >>> h1       = frame_project ( frame , h1_model , 'pt' )
    >>> h2_model = ...
    >>> h2       = frame_project ( frame , h2_model , 'pt,x' )
    >>> h3_model = ...
    >>> h3       = frame_project ( frame , h3_model ,  ( 'pt' , 'x'  , 'y' ) )
    >>> ...
    
    - Cuts/weight are also can be specified
    
    >>> frame    = ...
    >>> h1_model = ...
    >>> h1       = frame_project ( frame , h1_model , 'pt' , 'w')
    >>> h2_model = ...
    >>> h2       = frame_project ( frame , h2_model ,  ( 'pt' , 'x' ) , 'w')
    >>> h3_model = ...
    >>> h3       = frame_project ( frame , h3_model , 'pt' , 'x,y,z' , 'w' )
    >>>...
    
    - If model is a real 1D/2D/3D-hstogram, action is *not* lazy, otherwise action *is* lazy
    
    Model can be a real histogram or 
    -  `ROOT.ROOT.RDF.TH1DModel`
    -  `ROOT.ROOT.RDF.TH2DModel`
    -  `ROOT.ROOT.RDF.TH3DModel`
    -   anything that can be converted to :
    -- `ROOT.ROOT.RDF.TH1DModel`, 
    -- `ROOT.ROOT.RDF.TH2DModel` or
    -- `ROOT.ROOT.RDF.TH3DModel` objects 
    
    """

    if progress and isinstance ( frame , ROOT.TTree ) : progress = len ( frame )

    # ===================================================================================
    ## ROOT.TProfile is not supprted yet (ROOT 6.36)
    if isinstance ( model , ROOT.TProfile3D ) :
        if isinstance ( frame , ROOT.TTree ) : 
            logger.warning ( 'Frame-based projection is not supported for %s' % typename ( model ) )
            from ostap.trees.trees import tree_project as _tree_project_ 
            return _tree_project_ ( frame                ,
                                    model                , 
                                    expressions          ,
                                    cuts                 ,
                                    progress  = progress ,
                                    use_frame = False    ,
                                    parallel  = False    )
        logger.error ( 'Frame-based projection is not supported for %s' % typename ( model ) )
        return

    if isinstance ( model , _types_nD ) : 
        return frame_param ( frame                   ,
                             model                   ,
                             expressions             , 
                             cuts        = cuts      ,
                             as_weight   = as_weight ,  
                             progress    = progress  ,
                             report      = report    ,
                             lazy        = lazy      ) 
    
    ## decode expressions & cuts 
    current , items, cname , _ = _fr_helper_ ( frame , expressions , cuts , progress = progress )
    
    ## add the fiducial cuts
    cvars = tuple ( v for v in items.values () ) 
    if isinstance ( model , ROOT.TH3 ) :
        zvar    = cvars [ 2 ] 
        axis    = model.GetZaxis()
        current = current.Filter ( '%.12g <= %s ' %  ( axis.GetXmin() , zvar ) , 'ZMIN-FILTER' )
        current = current.Filter ( '%.12g >= %s ' %  ( axis.GetXmax() , zvar ) , 'ZMAX-FILTER' )

    if isinstance ( model , ROOT.TH2 ) :
        yvar    = cvars [ 1 ]         
        axis    = model.GetYaxis()
        current = current.Filter ( '%.12g <= %s ' %  ( axis.GetXmin() , yvar ) , 'YMIN-FILTER' )
        current = current.Filter ( '%.12g >= %s ' %  ( axis.GetXmax() , yvar ) , 'YMAX-FILTER' )

    if isinstance ( model , ROOT.TH1 ) :
        xvar    = cvars [ 0 ]                 
        axis    = model.GetXaxis()
        current = current.Filter ( '%.12g <= %s ' %  ( axis.GetXmin() , xvar ) , 'XMIN-FILTER' )
        current = current.Filter ( '%.12g >= %s ' %  ( axis.GetXmax() , xvar ) , 'XMAX-FILTER' )

    ## convert histogram-like objects into 'models'
    
    histo = None
    
    ## if   isinstance ( model , ROOT.TProfile3D ) : histo, model  = model, model.model ()        
    if   isinstance ( model , ROOT.TProfile2D ) : histo, model  = model, model.model ()        
    elif isinstance ( model , ROOT.TProfile   ) : histo, model  = model, model.model ()        
    elif isinstance ( model , ROOT.TH3        ) : histo, model  = model, model.model ()        
    elif isinstance ( model , ROOT.TH2        ) : histo, model  = model, model.model ()        
    elif isinstance ( model , ROOT.TH1        ) : histo, model  = model, model.model ()        

    if histo : histo.Reset()

    ## if true histo is specified, the action is NOT lazy!
    if histo and lazy :
        lazy  = False
        logger.warning ( "True histo is specified, processing is *NOT* lazy!" )
            
    ## ATTENTION HERE!!
    if cname and not as_weight :
        logger.warning ( "The cut is treated as boolean: %s" % cuts ) 
        cname = ''
        
    nvars = len ( items )
    pvars = [ v for v in items.values() ] 
    if cname : pvars.append ( cname )

    ## if   4 == nvars and isinstance ( model , DF_P3Model ) : action = current.Profile3D ( model , *pvars )
    if   3 == nvars and isinstance ( model , DF_P2Model ) : action = current.Profile2D ( model , *pvars )
    elif 2 == nvars and isinstance ( model , DF_P1Model ) : action = current.Profiel1D ( model , *pvars )
    elif 3 == nvars and isinstance ( model , DF_H3Model ) : action = current.Histo3D   ( model , *pvars )
    elif 2 == nvars and isinstance ( model , DF_H2Model ) : action = current.Histo2D   ( model , *pvars )
    elif 1 == nvars and isinstance ( model , DF_H1Model ) : action = current.Histo1D   ( model , *pvars )
    elif 1 <  nvars and isinstance ( model , DF_H1Model ) and histo  :
        
        ## book list of actions
        args    = ( cname , ) if cname else () 
        actions = [ current.Histo1D   ( histo.model () , var , *args ) for var in items.values() ]
        ## run all actions & prepare the final resuts 
        for a in actions : histo  += a.GetValue()
        if report :
            title  = 'frame_project (stack)' 
            report = current.Report()
            logger.info ( '%s\n%s' % ( title , report_print ( report , title = title , prefix = '# ') ) )
        return histo 
        
    else :
        raise TypeError ('Invalid model/what objects %s %s ' % ( type ( model ) , str ( pvars ) ) ) 

    if lazy :
        return action , current  ## RETUTRN

    ## make the actual looping 
    result = get_values ( action ,
                          frame  = current ,
                          report = report  ,
                          title  = 'frame_project' )
    
    if histo :
        histo  += result
        result  = histo 
         
    return result

# =============================================================================
## "project/parameterise" frame into polynomial structures (lazy action)
#  @code
#  frame = ...
#
#  ls    = Ostap.Math.LegendreSum ( ... )
#  res   = frame_param ( frame , ls , 'x' , 'y>0' )
#
#  bs    = Ostap.Math.Bernstein   ( ... )
#  res   = frame_param ( frame , bs , 'x' , 'y>0' )
#
#  cs    = Ostap.Math.ChebyshevSum ( ... )
#  res   = frame_param ( frame , cs , 'x' , 'y>0' )
#
#  ls2   = Ostap.Math.LegendreSum2 ( ... )
#  res   = frame_param ( frame , ls2 , 'y' , 'x' , 'z>0' )
# 
#  bs2   = Ostap.Math.Bernstein2D  ( ... )
#  res   = frame_param ( frame , bs2 ,  ( 'y' , 'x' ) , 'z>0' )
#
#  ls3   = Ostap.Math.LegendreSum3 ( ... )
#  res   = frame_param ( frame , ls3 ,  ( 'z' , 'y' , 'x' ) , 'z>0' )
# 
#  bs2   = Ostap.Math.Bernstein2D  ( ... )
#  res   = frame_param ( frame , bs2 , ( 'y' , 'x' ) , 'z>0' )
#  @endcode
def frame_param ( frame               ,
                  target              ,
                  expressions         , 
                  cuts        = ''    , * , 
                  as_weight   = True  , 
                  progress    = False ,
                  report      = False ,
                  lazy        = False ) :
    """ `project/parameterise` frame into polynomial structures
    
    >>> frame = ...
    
    >>> ls    = Ostap.Math.LegendreSum ( ... )
    >>> res   = frame_the_param ( frame , ls , 'x' , 'y>0' )
    
    >>> bs    = Ostap.Math.Bernstein   ( ... )
    >>> res   = frame_the_param ( frame , bs , 'x' , 'y>0' )

    >>> cs    = Ostap.Math.ChebyshevSum ( ... )
    >>> res   = frame_the_param ( frame , cs , 'x' , 'y>0' )

    >>> ls2   = Ostap.Math.LegendreSum2 ( ... )
    >>> res   = frame_the_param ( frame , ls2 , 'y,x' , 'z>0' )

    >>> bs2   = Ostap.Math.Bernstein2D  ( ... )
    >>> res   = frame_the_param ( frame , bs2 , 'y,x' , 'z>0' )

    >>> ls3   = Ostap.Math.LegendreSum3 ( ... )
    >>> res   = frame_the_param ( frame , ls3 , 'z,y;z' , 'z>0' )

    >>> bs2   = Ostap.Math.Bernstein2D  ( ... )
    >>> res   = frame_the_param ( frame , bs2 , 'y,x' , 'z>0' )
    """
    
    if progress and isinstance ( frame , ROOT.TTree ) : progress = len ( frame )

    ## the histogram ? 
    if isinstance ( target , ROOT.TH1 ) :
        return frame_project ( frame                   ,
                               target                  ,
                               expressions             ,
                               cuts        = cuts      ,
                               as_weight   = as_weight ,  
                               progress    = progress  ,
                               report      = report    ,
                               lazy        = lazy      )
    
    ## 
    current , items , cname , _  = \
        _fr_helper_ ( frame , expressions ,
                      cuts     = cuts     , 
                      progress = progress )

    ## ATTENTION HERE!!
    if cname and not as_weight :
        logger.warning ( "The cut is treated as boolean: %s" % cuts ) 
        cname = ''
    
    nvars = len ( items )

    assert \
        ( 1 == nvars and isinstance ( target , _types_1D ) ) or \
        ( 2 == nvars and isinstance ( target , _types_2D ) ) or \
        ( 3 == nvars and isinstance ( target , _types_3D ) ) or \
        ( 4 == nvars and isinstance ( target , _types_4D ) ) ,  \
        "Invalid structure of  polynomial and variables/cuts!"

    ## variables 
    uvars = [ k for k in items.values() ]

    if 4 <= len ( uvars ) :        
        if hasattr ( target , 'umin' ) :
            current = current.Filter ( '%.10g <= %s ' %  ( target.umin() , uvars[3] ) , 'UMIN-FILTER' )
        if hasattr ( target , 'umax' ) :
            current = current.Filter ( '%.10g >= %s ' %  ( target.umax() , uvars[3] ) , 'UMAX-FILTER' )
        
    if 3 <= len ( uvars ) :
        if hasattr ( target , 'zmin' ) :
            current = current.Filter ( '%.10g <= %s ' %  ( target.zmin() , uvars[2] ) , 'ZMIN-FILTER' )
        if hasattr ( target , 'zmax' ) :
            current = current.Filter ( '%.10g >= %s ' %  ( target.zmax() , uvars[2] ) , 'ZMAX-FILTER' )
        
    if 2 <= len ( uvars ) :
        if hasattr ( target , 'ymin' ) :
            current = current.Filter ( '%.10g <= %s ' %  ( target.ymin() , uvars[1] ) , 'YMIN-FILTER' )
        if hasattr ( target , 'ymax' ) :
            current = current.Filter ( '%.10g >= %s ' %  ( target.ymax() , uvars[1] ) , 'YMAX-FILTER' )
        
    if 1 <= len ( uvars ) :
        if hasattr ( target , 'xmin' ) :
            current = current.Filter ( '%.10g <= %s ' %  ( target.xmin() , uvars[0] ) , 'ZMIN-FILTER' )
        if hasattr ( target , 'xmax' ) :
            current = current.Filter ( '%.10g >= %s ' %  ( target.xmax() , uvars[0] ) , 'ZMAX-FILTER' )

    if cname : uvars.append ( cname )
    uvars = CNT ( uvars )

    ## reset the target 
    target.reset()
    
    TT = type ( target )
    
    if   isinstance ( target , Ostap.Math.LegendreSum4 ) :
        action = SA4w [ TT ] ( target ) if cname else SA4 [ TT ] ( target )
        
    elif isinstance ( target , Ostap.Math.LegendreSum3 ) :
        action = SA3w [ TT ] ( target ) if cname else SA3 [ TT ] ( target )
        
    elif isinstance ( target , Ostap.Math.LegendreSum2 ) :
        action = SA2w [ TT ] ( target ) if cname else SA2 [ TT ] ( target )
        
    elif isinstance ( target , Ostap.Math.LegendreSum  ) :
        action = SA1w [ TT ] ( target ) if cname else SA1 [ TT ] ( target )
        
    elif isinstance ( target , Ostap.Math.ChebyshevSum ) :
        action = SA1w [ TT ] ( target ) if cname else SA1 [ TT ] ( target )
        
    elif isinstance ( target , Ostap.Math.Bernstein3D  ) :
        action = SA3w [ TT ] ( target ) if cname else SA3 [ TT ] ( target )
        
    elif isinstance ( target , Ostap.Math.Bernstein2D  ) :
        action = SA2w [ TT ] ( target ) if cname else SA2 [ TT ] ( target )
        
    elif isinstance ( target , Ostap.Math.Bernstein    ) :
        action = SA1w [ TT ] ( target ) if cname else SA1 [ TT ] ( target )
        
    elif isinstance ( target , Ostap.Math.WStatEntity  ) : action = SA1w [ TT ] ( target )    
    elif isinstance ( target , Ostap.Math.NStatEntity  ) : action = SA1  [ TT ] ( target )
    elif isinstance ( target , Ostap.Math. StatEntity  ) : action = SA1  [ TT ] ( target )
     
    elif isinstance ( target , Ostap.Math.WStatictic4  ) : action = SA4w [ TT ] ( target )
    elif isinstance ( target , Ostap.Math. Statictic4  ) : action = SA4  [ TT ] ( target )
    
    elif isinstance ( target , Ostap.Math.WStatictic3  ) : action = SA3w [ TT ] ( target )
    elif isinstance ( target , Ostap.Math. Statictic3  ) : action = SA3  [ TT ] ( target )
    
    elif isinstance ( target , Ostap.Math.WStatictic2  ) : action = SA2w [ TT ] ( target )
    elif isinstance ( target , Ostap.Math. Statictic2  ) : action = SA2  [ TT ] ( target )
    
    elif isinstance ( target , Ostap.Math.WStatictic2  ) : action = SA1w [ TT ] ( target )
    elif isinstance ( target , Ostap.Math. Statictic2  ) : action = SA1  [ TT ] ( target )
    
    else :
        
        raise TypeError ( "Unknown type of target: %s" % typename ( target  ) )

    ## Book the action:
    result = current.Book ( std_move ( action ) , uvars ) 
    
    ## lazy processing ? 
    if lazy :
        return result , current
    
    ## make the real looping
    result  = get_values ( result ,
                           frame  = current ,
                           report = report  ,
                           title  = 'frame_param' )

    ## update target
    target += result
    
    return target 

# ==========================================================================
## draw the variable(s) from the frame
#  @code
#  frame = ...
#  result = frame_draw ( frame , 'x+12/z' , cut = 'z>1' ) 
## @endcode 
def frame_draw ( frame               ,
                 expressions         ,
                 cuts        = ''    ,
                 opts        = ''    , * , 
                 as_weight   = True  , 
                 delta       = 0.01  , 
                 progress    = False ,
                 report      = False , **kwargs ) :
    """ Draw the variable(s) from the frame
    >>> frame = ...
    >>> result = frame_draw ( frame , 'x+12/z' , cut = 'z>1' ) 
    """

    if progress and isinstance ( frame , ROOT.TTree ) : progress = len ( frame )
 
    ## decode expressions & cuts 
    current , items, cname , _ = \
        _fr_helper_ ( frame , expressions , cuts , progress = progress )

    ## ATTENTION HERE!!
    if cname and not as_weight :
        logger.warning ( "The cut is treated as boolean: %s" % cuts ) 
        cname = ''

    nvars = len ( items )
    assert 1 <= nvars <= 3 , 'Invalid expressions: %s' % str ( items ) 

    cvars  = [ v for v in items.values() ]
    uvars  = CNT ( cvars ) if not cname else CNT ( cvars + [ cname ] )

    ## create the cache ... needed ? 
    ## cache   = current.Cache ( uvars )
    cache   = current 

    ## 11st explicit loop 
    ranges = frame_range ( cache               ,
                           cvars               ,
                           cuts     = cname    ,
                           delta    = delta    ,
                           report   = report   ,
                           progress = progress )

    if not ranges :
        ## remove cuts and recalculate the ranges 
        ranges = frame_range ( current ,
                               cvars   ,
                               cuts     = ''       ,
                               delta    = delta    ,
                               report   = report   ,
                               progress = progress )
        if not ranges :         
            logger.warning ( 'frame_draw: nothing to draw, return None' )
            return None

    kw = cidict ( transform = cidict_fun , **kwargs )

    histos = [] 
    for key, var in loop_items ( items ) :
        mn , mx = ranges [ var ]
        item = key , ( mn , mx ) 
        histos.append ( item )
        
    ## book the histogram
    from   ostap.histos.histos       import histo_book2
    histo = histo_book2 ( histos , kw )

    ## fill the histogram 
    histo = frame_project ( cache ,
                            histo ,
                            cvars ,
                            cuts      = cname     ,
                            as_weight = as_weight , 
                            progress  = progress  ,
                            report    = report    ,
                            lazy      = False     )

    ## draw the histogram
    histo.draw ( opts , **kw )

    return histo 

# =============================================================================
## @class SliceHelper
#  helper class to get the slice from the frame
class SliceHelper(object) :
    """ Helper class to get the slice from the frame
    """
    
    def __init__ ( self              ,
                   result            ,   ## RAW result 
                   vardct            ,   ## name <-> name disctionary
                   cname      = ''   ,   ## name for the cuts-variable 
                   structured = True ,   ## structured ?
                   transpose  = True ) : ## transpose 
        
        self.raw_result = result ## RAW result 
        self.vardct     = vardct
        self.cname      = cname 
        self.structured = True if structured else False 
        self.transpose  = True if transpose  else False 
        
    ## get the value 
    def GetValue ( self ) :
        
        values = self.raw_result.GetValue()
        
        num    = 0 
        result = {} 
        dtypes = []

        import numpy
        
        for key, value in self.vardct.items() :
            
            dtypes.append ( ( key , numpy.float64 ) )
            result [ key ] = values [ value ]
            num = max ( num , len ( result [ key ] ) )            
            
        if self.structured :            
            part = numpy.zeros ( num , dtype = dtypes )        
            for v in result : part [ v ] = result [ v ]
            result = part            
        else :            
            result = [ r for r in result.values() ]
            result = numpy.stack ( result )
            ## 
            if self.transpose : result = numpy.transpose ( result )

        ## get the weights 
        weights = None
        if self.cname :
            weights = values [ cname ]
            if numpy.all ( weights == 1 ) : weight = None

        return result, weights 
                
# ==============================================================================
## get slice form frame as numpy array
#  @code
# 
#  @endcode
# ==============================================================================
def frame_slice ( frame               ,
                  expressions         , 
                  cuts        = ""    ,
                  structured  = True  ,
                  transpose   = True  ,
                  progress    = False , 
                  lazy        = False ) :
    
    """ Get the slice from frame as numpy array 
    """
    
    ## decode expressions & cuts 
    current , items , cname , _ = _fr_helper_ ( frame , expressions , cuts , progress = progress )
    
    columns = [ v for v in items.values () ]
    if cname : columns.append ( cname )
    
    result = current.AsNumpy ( columns = columns , lazy = True )
    
    result  = SliceHelper ( result                  ,
                            items                   ,
                            cname                   ,
                            structured = structured ,
                            transpose  = transpose  ) 
    
    return ( result , current ) if lazy else result.GetValue()
    
# =============================================================================
## @class EffHelper
#  helper class to get the efficiency from the frame
#  @see frame_efficiency 
class EffHelper(object) :
    """ Helper class to get the slice from the frame
    - see frame_efficiency 
    """    
    def __init__ ( self     ,
                   accepted  , ## RAW result for accepted histogram
                   rejected  , ## RAW result for rejected histogram
                   title     ) : 
        
        self.accepted = accepted
        self.rejected = rejected 
        self.__title  = title
        
    ## get the value 
    def GetValue ( self ) :
        
        hA = self.accepted.GetValue()
        hR = self.rejected.GetValue()
        hE = 1 / ( 1 + hR / hA )
        hE.SetTitle ( self.__title )

        return hE , hA , hR 

# ==============================================================================
## Get "Efficincy" histogram for frame
#  Actually it makes two distributions/histograms
# for events that are accepted and rejected by certain
# (boolean) <code>criterion</code>
# and calcualted the efficiency histogram usng the rule: 
# \f[ e = \left( 1 + \frac{r}{a} \right)^{-1}\f]
#
#  @code
#  histo_1D = ...
#  frae     = ...
#  eff, accepted, rejected = frame_efficiency
#  ...   ( frame        ,  
#  ...     'DLL>5'      , 
#  ...     histo_1D     ,
#  ...     'PT'         ,
#  ...     cuts = 'A>2' ,
#          lazy = False ) 
#  @code
#
#  @code
#  histo_1D = ...
#  frae     = ...
#  result , frame = frame_efficiency
#  ...   ( frame        ,  
#  ...     'DLL>5'      , 
#  ...     histo_1D     ,
#  ...     'PT'         ,
#  ...     cuts = 'A>2' ,
#          lazy = True  ) 
#  @code
#
#
#  @param frame       (INPUT) input frame (or TTree)
#  @param criterion   (INPUT) (boolean) criterion, e.g <code>PT>10</code>
#  @param histo       (INPUT) the 1D/2D or 3D (template) histogram
#  @param expressions (INPUT) expressions rtaht defined the axes of the histogram
#  @param cuts        (INPUT) only event that satisfy these (boolean) cuts are procesed
#  @param weight      (INPUT) optional weihgt to be used for filling the historgams
#  @param progress    (INPUT) show the porogress bar for the actual processing?
#  @param report      (INPUT) produce the processing report?
#  @param lazy        (INPUT) lazy processing?
#  @attention  criterion and cuts are treated as boolean! 
#  @see tree_efficiency
#  @see data_efficiency
def frame_efficiency ( frame               ,
                       criterion           ,  
                       histo               , 
                       expressions         ,
                       cuts        = ''    , 
                       weight      = ''    , 
                       progress    = False ,
                       report      = False , 
                       lazy        = False ) :
    """ Get "Efficincy" histogram for frame
    Actually it makes two distributions/histograms
    for events that are accepted and rejected by certain
    (boolean) <code>criterion</code>
    and calculated  the efficiency histogram
    
    >>> histo_1D = ...
    >>> frame     = ...
    >>> eff, accepted, rejected = frame_efficiency
    ...     ( frame        ,  
    ...       'DLL>5'      , 
    ...       histo_1D     ,
    ...       'PT'         ,
    ...       cuts = 'A>2' ,
    ...       lazy = False ) ## attention! 

    >>> histo_1D = ...
    >>> framee   = ...
    >>> result , frame = frame_efficiency
    ...     ( frame        ,  
    ...       'DLL>5'      , 
    ...       histo_1D     ,
    ...       'PT'         ,
    ...       cuts = 'A>2' ,
    ...       lazy = True  ) 
    
    frame       : (INPUT) input frame (or TTree)
    criterion   : (INPUT) (boolean) criterion, e.g <code>PT>10</code>
    histo       : (INPUT) the 1D/2D or 3D (template) histogram
    expressions : (INPUT) expressions rtaht defined the axes of the histogram
    cuts        : (INPUT) only event that satisfy these (boolean) cuts are procesed
    weight      : (INPUT) optional weihgt to be used for filling the historgams
    progress    : (INPUT) show the porogress bar for the actual processing?
    report      : (INPUT) produce the processing report?
    lazy        : (INPUT) lazy processing?

    ATTENTION: criterion and cuts are treated as boolean! 
    - see tree_efficiency
    - see data_efficiency
    """
    
    ## decode expressions & cuts 
    current , items , cut_name , _ = _fr_helper_ ( frame , expressions , cuts , progress = progress )
    nvars = len ( items )

    _ , criterion , _ = SV.vars_and_cuts ( expressions , criterion ) 
    assert criterion , "Valid criterion *MUST* be specified!"
    
    _ , weight    , _ = SV.vars_and_cuts ( expressions , weight    ) 
    
    assert isinstance ( histo , ROOT.TH1 ) , "Invalid `histo` type : %s " % typename ( histo )
    
    hdim = histo.GetDimension()    
    assert 1 <= hdim <= 3 and hdim == nvars , \
        "Mismatch histogram dimension/#vars: %s/%d"% ( hdim, nvars ) 

    ## get all current variables

    ## fiducial cuts
    cvars = tuple ( v for v in items.values () ) 
 
    if isinstance ( histo , ROOT.TH3 ) :
        zvar    = cvars [ 2 ] 
        axis    = histo.GetZaxis()
        current = current.Filter ( '%.12g <= %s ' %  ( axis.GetXmin() , zvar ) , 'ZMIN-FILTER' )
        current = current.Filter ( '%.12g >= %s ' %  ( axis.GetXmax() , zvar ) , 'ZMAX-FILTER' )
        
    if isinstance ( histo , ROOT.TH2 ) :
        yvar    = cvars [ 1 ]         
        axis    = histo.GetYaxis()
        current = current.Filter ( '%.12g <= %s ' %  ( axis.GetXmin() , yvar ) , 'YMIN-FILTER' )
        current = current.Filter ( '%.12g >= %s ' %  ( axis.GetXmax() , yvar ) , 'YMAX-FILTER' )
        
    if isinstance ( histo , ROOT.TH1 ) :
        xvar    = cvars [ 0 ]                 
        axis    = histo.GetXaxis()
        current = current.Filter ( '%.12g <= %s ' %  ( axis.GetXmin() , xvar ) , 'XMIN-FILTER' )
        current = current.Filter ( '%.12g >= %s ' %  ( axis.GetXmax() , xvar ) , 'XMAX-FILTER' )
   

    histo.Reset()
    if not histo.GetSumw2() : histo.Sumw2() 
    
    all_cols = frame_columns ( current ) 
    all_vars = set ( all_cols )

    ##
    wname = ''
    if weight :
        wname = var_name ( 'weight_'  , all_vars , *all_cols )
        all_vars.add ( wname ) 
        wcut  = var_name ( 'weight_cut_' , all_vars , *all_cols )
        all_vars.add ( wcut  ) 
        current = current.Define ( wname , '1.0*(%s)' % weight ) ## weight as variable
        current = current.Define ( wcut  , '!!(%s)'   % wname  ) ## weight as cut 
        current = current.Filter ( wcut  , "FILTER-WEIGTH: %s" % weight)
        
    ## add criterion 
    critname = var_name ( 'criterion_' , all_vars , *all_cols )
    current  = current.Define ( critname , '!!(%s)' % criterion ) ## criterion as boolean 

    ## split into two branches
    accepted = current.Filter (           critname , 'ACCEPTED' )
    rejected = current.Filter ( '!(%s)' % critname , 'REJECTED' )


    ## the models 
    modelA = histo.model()
    modelR = histo.model()

    modelA.fTitle = "Distribution for events accepted by %s" % criterion
    modelR.fTitle = "Distribution for events rejected by %s" % criterion
    etitle        = "Efficiency   for criterion: %s"         % criterion 

    pvars = cvars if not wname else cvars + ( wname ,  ) 

    if   1 == hdim : 
        hA = accepted.Histo1D ( modelA , *pvars ) 
        hR = rejected.Histo1D ( modelR , *pvars )
    elif 2 == hdim : 
        hA = accepted.Histo2D ( modelA , *pvars ) 
        hR = rejected.Histo2D ( modelR , *pvars )
    elif 3 == hdim : 
        hA = accepted.Histo3D ( modelA , *pvars ) 
        hR = rejected.Histo3D ( modelR , *pvars )

    result = EffHelper ( hA , hR , etitle )
    
    if lazy :
        return result, current
    
    ## get the actual triplet of histograms 
    result  = result.GetValue()

    histo.SetTitle ( etitle ) 
    histo  += result [ 0 ]  

    if report :
        the_report = current.Report()
        title      = "frame_efficiency"
        logger.info ( '%s\n%s' % ( title , report_print ( the_report , title = title , prefix = '# ') ) )
   
    return result 


# ==============================================================================
# decorate 
# ==============================================================================
for f in frame_types :
    f.columns         = frame_columns          , ## defined columns/branches
    f.branches        = frame_branches         , ## defined columns/branches
    f.length          = frame_length           , ## length, #godo entries 
    f.size            = frame_size             , ## length, #godo entries 
    ## 
    f.statistic       = frame_statistic        , ## statistci for variable(s)
    f.statVar         = frame_statistic        , ## statistci for variable(s) 
    f.statVars        = frame_statistic        , ## statistci for variable(s)
    ##
    f.minmax          = frame_minmax           , ## min/max for varable(s)
    f.ranges          = frame_range            , ## range   for varable(s)
    f.ranges          = frame_range            , ## range   for varable(s)
    ##
    f.the_moment      = frame_the_moment       , ## get statoistics in form os moments
    ##
    f.nEff            = frame_nEff             , ## number of effective entries 
    f.mean            = frame_mean             , ## mean value(s) for variable(s)
    f.variance        = frame_variance         , ## variance(s)   for variable(s)
    f.dispersion      = frame_dispersion       , ## dispersion(s) for variable(s)
    f.rms             = frame_rms              , ## RMS(s)        for variable(s)
    f.skewness        = frame_skewness         , ## skewness(s)   for variable(s)
    f.kurtosis        = frame_kurtosis         , ## skewness(s)   for variable(s)
    ##
    f.arithmetic_mean = frame_arithmetic_mean  , ## Arithmetic mean 
    f.harmonic_mean   = frame_harmonic_mean    , ## Harmonic   mean 
    f.geometric_mean  = frame_geometric_mean   , ## Harmonic   mean 
    f.power_mean      = frame_power_mean       , ## Power      mean 
    f.lehmer_mean     = frame_lehmer_mean      , ## Lehmer     mean 
    ## 
    f.ECDF            = frame_ECDF             , ## Empirical Cumulative Distribution function
    ##
    f.project         = frame_project          , ## project data frame to the (1D/2D/3D) histogram
    f.param           = frame_param            , ## parameterize data on-flight     
    f.draw            = frame_draw             , ## draw variable from the frame
    ##
    f.slice           = frame_slice            , ## get the (numpy) slice from the frame
    f.efficiency      = frame_efficiency       , ## get the efficiency histos from the frame 
    
_new_methods_       = (
    report_print_table    , ## print the report 
    report_as_table       , ## print the report
    report_print          , ## print the report 
    ##
    frame_columns         , ## all defined columns 
    frame_branches        , ## all defined columns  
    ## 
    frame_progress1       , ## progress bar for frame (OK for ROOT<6.32)
    frame_progress2       , ## progress bar for frame (OK for ROOT>6.29) 
    frame_progress        , ## progress bar for frame
    ## 
    frame_prescale        , ## prescale frame (trivial filter)
    frame_print           , ## over-simplified print frame
    frame_table           , ## Print frame as detailed table 
    ## 
    frame_length          , ## length/size/#entries of the frame 
    frame_size            , ## length/size/#entries of the frame 
    ##
    frame_statistic       , ## statistics for the variable(s) 
    frame_minmax          , ## min/max for the variable(s) 
    frame_range           , ## ranges  for the variable(s)
    ##
    frame_the_moment      , ## get moments for the variable(s)
    ##
    frame_nEff            , ## number of effective entries 
    frame_mean            , ## mean value for variable(s)
    frame_variance        , ## variance   for variable(s)
    frame_dispersion      , ## dispersion for variable(s)
    frame_rms             , ## rms for variable(s)
    frame_skewness        , ## skewness for variable(s)
    frame_kurtosis        , ## kurtosis for variable(s)
    ##
    frame_arithmetic_mean , ## get the arithmetc mean for variables 
    frame_harmonic_mean   , ## get the harmonic  mean for variables 
    frame_geometric_mean  , ## get the geometric mean for variables 
    frame_power_mean      , ## get the power     mean for variables 
    frame_lehmer_mean     , ## get the Lehmer    mean for variables 
    ## 
    frame_ECDF            , ## get the empirical cumulative distribution functions for variable(s)
    ## 
    frame_project         , ## project data frame to the (1D/2D/3D) histogram
    frame_param           , ## parameterize data on-flight     
    frame_draw            , ## draw variable from the frame
    ## 
    frame_slice           , ## the get (numpy) slice form the frame 
    frame_efficiency      , ## the get efficiency histos form the frame 
    ## 
)

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
