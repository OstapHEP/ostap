#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/frames/frames.py
#  Module with decoration of TDataFrame objects for efficient use in python
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2018-06-16
# =============================================================================
"""Decoration of Tree/Chain objects for efficient use in python"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'DataFrame'             , ## RDataFrame object
    ## 
    'frame_columns'         , ## defined columns  
    'frame_branches'        , ## defined columns  
    ## 
    'frame_progress1'       , ## progress bar for frame (OK for ROOT<6.32)
    'frame_progress2'       , ## progress bar for frame (OK for ROOT>6.29) 
    'frame_progress'        , ## progress bar for frame
    'frame_prescale'        , ## prescale frame (trivial filter)
    'frame_project'         , ## project data frame to the (1D/2D/3D) histogram
    'frame_draw'            , ## draw variable from the frame  
    'frame_print'           , ## over-simplified print frame
    'frame_table'           , ## Print frame as detailed table 
    ##
    'frame_length'          , ## length/size of the frame 
    'frame_size'            , ## length/size of the frame 
    ##
    'frame_nEff'            , ## nEff function for frames
    'frame_statVar'         , ## stat var for frame 
    'frame_statCov'         , ## stat cov for frame
    ##
    'frame_get_moment'      , ## get the moment or certain order around defined center 
    'frame_moment'          , ## get the moment or certain order 
    'frame_central_moment'  , ## get the central moment or certain order
    ## 
    'frame_mean'           , ## mean value for the variable 
    'frame_rms'            , ## RMS for the variable 
    'frame_variance'       , ## variance for the variable 
    'frame_dispersion'     , ## dispersion for the variable
    'frame_skewness'       , ## skewness for the variable 
    'frame_kurtosis'       , ## kurtosos for the variable
    ## 
    'frame_quantile'       , ## quantile 
    'frame_interval'       , ## quantile interval 
    'frame_median'         , ## median 
    'frame_quantiles'      , ## quantiles 
    'frame_terciles'       , ## terciles 
    'frame_quartiles'      , ## quartiles 
    'frame_quintiles'      , ## quintiles 
    'frame_deciles'        , ## deciles 
    ##
    'report_print'         , ## print the report 
    'report_print_table'   , ## print the report 
    'report_as_table'      , ## print the report
    ##
    'frame_types'          , ## the basic DataFrame/FeameNode types 
    ) 
# =============================================================================
from   ostap.core.meta_info      import root_info 
from   ostap.core.ostap_types    import ( integer_types , dictlike_types , 
                                          num_types     , ordered_dict   )    
from   ostap.core.core           import cpp, Ostap
from   ostap.math.base           import isequal, iszero, axis_range                             
from   ostap.logger.utils        import multicolumn
from   ostap.utils.cidict        import cidict, cidict_fun      
from   ostap.utils.progress_conf import progress_conf 
from   ostap.utils.basic         import isatty, loop_items
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
Frames_OK    = False 
std_move     = None 
has_std_move = False
# =============================================================================
try : # =======================================================================
    # =========================================================================
    std_move     = ROOT.std.move
    has_std_move = True 
    Frames_OK    = has_std_move                                    ## ATTENTION!
    # =========================================================================
except AttributeError : # =====================================================
    # =========================================================================
    std_move     = None 
    has_std_move = False
    Frames_OK    = False 
# =============================================================================
## type for the column names 
CNT                = DataFrame.ColumnNames_t 
# =============================================================================
## get frame-like stuff  as rnode 
def as_rnode ( frame ) :
    """ Get frame-like stuff  as rnode"""
    return ROOT.RDF.AsRNode ( frame ) 

# ==============================================================================
## Helper method to generate new, unique name for the variable 
def var_name ( prefix , used_names , *garbage ) :
    """ Helper method to generate new, unique name for the variable 
    """
    name =     prefix + '%x' % ( hash ( ( prefix , used_names , garbage ) )               % ( 2**32 ) ) 
    while name in used_names :
        name = prefix + '%x' % ( hash ( ( name , prefix , used_names , name , garbage ) ) % ( 2**32 ) )
    return name

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

    vnames   = ordered_dict()     
    for i , expr in enumerate ( exprs , start = 1 ) :
        vname = expr
        if not expr in all_vars :
            used    = tuple ( all_vars | set ( frame_columns ( current ) ) ) 
            vn      = var_name ( 'var%d_' % i , used , expr , *vars )
            all_vars.add ( vn )
            current = current.Define ( vn , expr )
            vname   = vn
        vnames [ expr ] = vname 

    ## Filter frame, if needed!
    if cuts :
        ## add named filter 
        current = current.Filter ( '(bool) (%s)' % cuts, 'FILTER/WEIGHT: %s' % cuts )
        
    cname = cuts
    if cuts and not cuts in all_vars :
        used    = tuple ( all_vars | set ( frame_columns ( current ) ) ) 
        cn      = var_name ( 'cut_' , used , cuts , *vars )
        all_vars.add ( cn )
        ncuts   = '1.0*(%s)' % cuts 
        current = current.Define ( cn , ncuts )
        cname   = cn
        
    return current, vnames, cname, input_string

# ==================================================================================
## The second helper method to implement various "statistics"-related actions  
def _fr_helper2_ ( frame            ,
                   creator          ,
                   expressions      ,
                   cuts     = ''    ,
                   progress = False ,
                   report   = False ,
                   lazy     = True  ) :
    """ The second helper method to implement various statistic-related actions  
    """

    current, var_names, cut_name, input_string = _fr_helper_ ( frame , expressions , cuts , progress = progress )

    results = {}
    for expr, var_name in loop_items ( var_names ) : 
        results [ expr ] = creator ( current , var_name , cut_name ) 

    lasy = False 
    if not lazy :
        for expr, res  in loop_items ( results ) :
            results [ expr ] = res.GetValue()
            rr = results [ expr ]

    if report and not lazy :
        report = current.Report()
        title  = 'DataFrame processing'
        logger.info ( '%s\n%s' % ( title , report_print ( report , title = title , prefix = '# ') ) )
        
    if input_string and 1 == len ( results ) :
        _ , r = results.popitem()
        return r
    
    return results 

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
    node = as_rnode ( frame )
    
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
    
    elif 1 == prescale :
        ## no action
        return node 
        
    raise TypeError ( "Invalid type/value for 'prescale' %s/%s" %( prescale , type ( prescale ) ) )

# ==============================================================================
## Draw the variables for the frame
#  - old variant 
def frame_draw ( frame , *args , **kwargs ) :
    """ Draw the variable(s) for the frames
    - old variant
    """
    node = as_rnode ( frame )
    from ostap.fitting.dataset import ds_draw as _ds_draw_
    return _ds_draw_ ( node , *args , **kwargs ) 

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
def frame_table ( frame , pattern = None ,  cuts = '' , more_vars = () , title = '' ,  prefix = '' ) :
    """ Data frame as table
    >>> frame = ...
    >>> table = frame_table ( frame , '.*PT.*' , cuts = ... , more_vars = [ 'x*x/y' , 'y+z'] )
    >>> print ( table )
    """
    
    if isinstance ( frame  , ROOT.TTree ) : frame = DataFrame ( frame  )
    node = as_rnode ( frame )
    
    ## the basic column type 
    def col_type ( var ) :
        """The basic column type"""
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

    svars = vcols + tuple  ( more_vars ) 
    stats = frame_statVar  ( node , svars , cuts = cuts ) 

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
# Frame -> histogram prjections 
# =============================================================================
    
DF_H1Model = ROOT.ROOT.RDF.TH1DModel 
DF_H1Type  = ROOT.TH1D
DF_H2Model = ROOT.ROOT.RDF.TH2DModel 
DF_H2Type  = ROOT.TH2D
DF_H3Model = ROOT.ROOT.RDF.TH3DModel 
DF_H3Type  = ROOT.TH3D

DF_P1Model = ROOT.ROOT.RDF.TProfile1DModel 
DF_P1Type  = ROOT.TProfile

DF_P2Model = ROOT.ROOT.RDF.TProfile2DModel
DF_P2Type  = ROOT.TProfile2D

# =============================================================================
## convert 1D-histogram to "model" for usage with DataFrames 
def _h1_model_ ( histo ) :
    """ Convert 1D-histogram to 'model' for usage with DataFrame"""
    model = histo 
    if not isinstance ( model , DF_H1Type ) :
        model = DF_H1Type()
        histo.Copy ( model )
    return DF_H1Model ( model ) 
# =============================================================================
## convert 2D-histogram to "model" for usage with DataFrames 
def _h2_model_ ( histo ) :
    """ Convert 2D-histogram to 'model' for usage with DataFrame"""
    model = histo 
    if not isinstance ( model , DF_H2Type ) :
        model = DF_H2Type()
        histo.Copy ( model )
    return DF_H2Model ( model ) 
# =============================================================================
## convert 3D-histogram to "model" for usage with DataFrames 
def _h3_model_ ( histo ) :
    """ Convert 3D-histogram to 'model' for usage with DataFrame"""
    model = histo 
    if not isinstance ( model , DF_H3Type ) :
        model = DF_H3Type()
        histo.Copy ( model )
    return DF_H3Model ( model ) 
# =============================================================================
## convert 1D-profile to "model" for usage with DataFrames 
def _p1_model_ ( histo ) :
    """ Convert 1D-profile to 'model' for usage with DataFrame"""
    model = histo 
    if not isinstance ( model , DF_P1Type ) :
        model = DF_P1Type()
        histo.Copy ( model )
    return DF_P1Model ( model ) 
# =============================================================================
## convert 2D-profile to "model" for usage with DataFrames 
def _p2_model_ ( histo ) :
    """Convert 2D-profile to 'model' for usage with DataFrame"""
    model = histo 
    if not isinstance ( model , DF_P2Type ) :
        model = DF_P2Type()
        histo.Copy ( model )
    return P2Model ( model ) 
# ==============================================================================
ROOT.TH1.model        = _h1_model_
ROOT.TH2.model        = _h2_model_
ROOT.TH3.model        = _h3_model_
ROOT.TProfile  .model = _p1_model_
ROOT.TProfile2D.model = _p1_model_
# ==============================================================================
_types_1D = Ostap.Math.LegendreSum  , Ostap.Math.Bernstein   , Ostap.Math.ChebyshevSum , 
_types_2D = Ostap.Math.LegendreSum2 , Ostap.Math.Bernstein2D ,
_types_3D = Ostap.Math.LegendreSum3 , Ostap.Math.Bernstein3D , 
_types_4D = Ostap.Math.LegendreSum4 ,
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
#  @attention: variables should be listed in reverse order!  
def frame_project ( frame            , 
                    model            ,
                    expressions      ,
                    cuts      = ''   ,
                    progress = False , 
                    report   = False ,
                    lazy     = True  ) :
    """ Project of the frame into histigram 

    - attention: variables should be listed in reverse order! 
    
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

    if isinstance ( model , _types_nD ) : 
        return _fr_param_ ( frame               ,
                            model               ,
                            expression          ,
                            cuts     = cuts     ,
                            progress = progress ,
                            report   = report   ,
                            lazy     = lazy     ) 

    ## decode expressions & cuts 
    current , items, cname , _ = _fr_helper_ ( frame , expressions , cuts , progress = progress )

    ## convert histogram-like objects into 'models'

    histo = None
    if   isinstance ( model , ROOT.TProfile2D )                : histo, model  = model, model.model ()        
    elif isinstance ( model , ROOT.TProfile   )                : histo, model  = model, model.model ()        
    elif isinstance ( model , ROOT.TH3 ) and 3 == model.dim () : histo, model  = model, model.model ()        
    elif isinstance ( model , ROOT.TH2 ) and 2 == model.dim () : histo, model  = model, model.model ()        
    elif isinstance ( model , ROOT.TH1 ) and 1 == model.dim () : histo, model  = model, model.model ()        

    if histo : histo.Reset()

    nvars = len ( items )
    pvars = [ v for v in items.values() ] 
    if cname : pvars.append ( cname )

    if   3 == nvars and isinstance ( model , DF_P2Model ) : action = current.Profile2D ( model , *pvars )
    elif 2 == nvars and isinstance ( model , DF_P1Model ) : action = current.Profiel1D ( model , *pvars )
    elif 3 == nvars and isinstance ( model , DF_H3Model ) : action = current.Histo3D   ( model , *pvars )
    elif 2 == nvars and isinstance ( model , DF_H2Model ) : action = current.Histo2D   ( model , *pvars )
    elif 1 == nvars and isinstance ( model , DF_H1Model ) : action = current.Histo1D   ( model , *pvars )
    else :
        raise TypeError ('Invalid model/what objects %s %s ' % ( type ( model ) , str ( pvars ) ) ) 

    ## if true histo is specified, the action is NOT lazy!
    if histo :
        histo += action.GetValue() 
        result = histo 
    else :        
        result = action if lazy else action.GetValue()

    if report and not lazy :
         report = current.Report()
         title = 'DataFrame project'
         logger.info ( '%s\n%s' % ( title , report_print ( report , title = title , prefix = '# ') ) )
         
    return result 

# =============================================================================
## "project/parameterise" frame into polynomial structures (lazy action)
#  @code
#  frame = ...
#
#  ls    = Ostap.Math.LegendreSum ( ... )
#  res   = frame_the_param ( frame , ls , 'x' , 'y>0' )
#
#  bs    = Ostap.Math.Bernstein   ( ... )
#  res   = frame_the_param ( frame , bs , 'x' , 'y>0' )
#
#  cs    = Ostap.Math.ChebyshevSum ( ... )
#  res   = frame_the_param ( frame , cs , 'x' , 'y>0' )
#
#  ls2   = Ostap.Math.LegendreSum2 ( ... )
#  res   = frame_the_param ( frame , ls2 , 'y' , 'x' , 'z>0' )
# 
#  bs2   = Ostap.Math.Bernstein2D  ( ... )
#  res   = frame_the_param ( frame , bs2 ,  ( 'y' , 'x' ) , 'z>0' )
#
#  ls3   = Ostap.Math.LegendreSum3 ( ... )
#  res   = frame_the_param ( frame , ls3 ,  ( 'z' , 'y' , 'x' ) , 'z>0' )
# 
#  bs2   = Ostap.Math.Bernstein2D  ( ... )
#  res   = frame_the_param ( frame , bs2 , ( 'y' , 'x' ) , 'z>0' )
#  @endcode
def _fr_param_ ( frame            ,
                 poly             ,
                 expressions      ,
                 cuts     = ''    ,
                 progress = False ,
                 report   = False ,
                 lazy     = True  ) :
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
    if isinstance ( poly , ROOT.TH1 ) :
        return _fr_project_ ( frame , poly , expressions , cuts = cuts , progress = progress , report = report , lazy = lazy )

    ## 
    current , items, cname , _ = _fr_helper_ ( frame , expressions , cuts , progress = progress )
    
    nvars = len ( items )
    
    assert \
        ( 1 == nvars and isinstance ( poly , _types_1D ) ) or \
        ( 2 == nvars and isinstance ( poly , _types_2D ) ) or \
        ( 3 == nvars and isinstance ( poly , _types_3D ) ) or \
        ( 4 == nvars and isinstance ( poly , _types_4D ) ) ,  \
        "Invalid structure of  polynomial and variables/cuts!"

    ## variables 
    uvars = [ k for k in items.values() ]
    if cuts : uvars.append ( cname )
    uvars = CNT ( uvars )
    
    ## finally book the actions!
    if   isinstance ( poly , Ostap.Math.LegendreSum  ) :
        result = current.Book ( ROOT.std.move ( Ostap.Actions.LegendrePoly   ( poly ) ) , uvars )
    elif isinstance ( poly , Ostap.Math.LegendreSum2 ) :
        result = current.Book ( ROOT.std.move ( Ostap.Actions.LegendrePoly2  ( poly ) ) , uvars )
    elif isinstance ( poly , Ostap.Math.LegendreSum3 ) :
        result = current.Book ( ROOT.std.move ( Ostap.Actions.LegendrePoly3  ( poly ) ) , uvars )
    elif isinstance ( poly , Ostap.Math.LegendreSum4 ) :
        result = current.Book ( ROOT.std.move ( Ostap.Actions.LegendrePoly4  ( poly ) ) , uvars )
    elif isinstance ( poly , Ostap.Math.ChebyshevSum ) :
        result = current.Book ( ROOT.std.move ( Ostap.Actions.ChebyshevPoly  ( poly ) ) , uvars )
    elif isinstance ( poly , Ostap.Math.Bernstein    ) :
        result = current.Book ( ROOT.std.move ( Ostap.Actions.BernsteinPoly  ( poly ) ) , uvars )
    elif isinstance ( poly , Ostap.Math.Bernstein2D  ) :
        result = current.Book ( ROOT.std.move ( Ostap.Actions.BernsteinPoly2 ( poly ) ) , uvars )
    elif isinstance ( poly , Ostap.Math.Bernstein3D  ) :
        result = current.Book ( ROOT.std.move ( Ostap.Actions.BernsteinPoly3 ( poly ) ) , uvars )

    if not lazy :
        result = result.GetValue()
        poly  *= 0.0 
        poly  += result
        
    if report and not lazy :
         report = current.Report()
         title = 'DataFrame parameterisation'
         logger.info ( '%s\n%s' % ( title , report_print ( report , title = title , prefix = '# ') ) )

    return result

# ==========================================================================
## draw the variable(s) from the frame
#  @code
#  frame = ...
#  result = frame_draw ( frame , 'x+12/z' , cut = 'z>1' ) 
## @endcode 
def _fr_draw_ ( frame            ,
                expressions      ,
                cuts     = ''    ,
                opts     = ''    ,
                delta    = 0.01  , 
                progress = False ,
                report   = False , **kwargs ) :
    """ Draw the variable(s) from the frame
    >>> frame = ...
    >>> result = frame_draw ( frame , 'x+12/z' , cut = 'z>1' ) 
    """

    if progress and isinstance ( frame , ROOT.TTree ) : progress = len ( frame )

    ## decode expressions & cuts 
    current , items, cname , _ = _fr_helper_ ( frame , expressions , cuts , progress = progress )

    nvars = len ( items )
    assert 1 <= nvars <= 3 , 'Invalid expressions: %s' % str ( items ) 

    cvars  = [ v for v in items.values() ]
    uvars  = CNT ( cvars ) if not cname else CNT ( cvars + [ cname ] )

    ## create the cache ... needed ? 
    ## cache   = current.Cache ( uvars )
    cache   = current 

    ## get the ranges 
    ranges = frame_range ( cache , cvars , cname , delta = delta , report = report )
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
    histo = frame_project ( cache , histo , cvars , cname , progress = progress , report = report , lazy = False )

    ## draw the histogram
    histo.draw ( opts , **kw )

    return histo 

# =============================================================================
## Get the length/size of the data frame
#  @code
#  frame = ...
#  print ( len(frame) )
#  len   = frame_length ( frame ) 
#  len   = frame_size   ( frame ) ## ditto 
#  @endcode 
def frame_length ( frame , lazy = False ) :
    """ Get the length/size of the data frame
    >>> frame = ...
    >>> print len(frame)
    >>> len   = frame_length ( frame ) 
    >>> len   = frame_size   ( frame ) ## ditto 
    """
    node = as_rnode ( frame )
    cnt  = node.Count () 
    return cnt if lazy else cnt.GetValue() 

# ==============================================================================
## Size of the frame, same as frame length 
#  @see frame_length 
frame_size   = frame_length 

# =============================================================================
## Get the effective entries in data frame
#  @code
#  data = ...
#  neff = data.nEff('b1*b1')
#  @endcode
def frame_nEff ( frame , cuts = '' ) :
    """ Get the effective entries in data frame 
    >>> data = ...
    >>> neff = data.nEff('b1*b1')
    """
    if isinstance ( frame  , ROOT.TTree ) : frame = DataFrame ( frame  )
    node  = as_rnode ( frame )    
    return SV.data_nEff ( node , cuts )

# =============================================================================
## Get statistics for the  given expression in data frame
#  @code
#  data = ...
#  c1 = frame_statVar ( 'S_sw' , 'pt>10' ) 
#  c2 = frame_statVar ( 'S_sw' , 'pt>0'  )
#  @endcode
def frame_statVar ( frame , expression ,  cuts = '' ) :
    """ Get statistics for the  given expression in data frame
    >>> data = ...
    >>> c1 = frame_statVar( 'S_sw' , 'pt>10' ) 
    >>> c2 = framestatVar( 'S_sw' )
    """
    if isinstance ( frame  , ROOT.TTree ) : frame = DataFrame ( frame  )
    node = as_rnode ( frame ) 
    return SV.data_statistics ( node , expression,  cuts = cuts )

# =============================================================================
## Get empirical CDF 
#  @code
#  data = ...
#  c1 = frame_ECDF ( data , 'S_sw' ) 
#  c2 = frame_ECDF ( data , S_sw' , 'pt>0'  )
#  @endcode
def frame_ECDF ( frame , expression ,  cuts = '' ) :
    """ Get empirical CDF 
    >>> data = ...
    >>> c1 = frame_ECDF ( data , 'S_sw' ) 
    >>> c2 = frame_ECDF ( data , S_sw' , 'pt>0'  )
    """
    if isinstance ( frame  , ROOT.TTree ) : frame = DataFrame ( frame  )
    node = as_rnode ( frame ) 
    return SV.data_ECDF ( node , expression,  cuts = cuts )

# =============================================================================
## get the statistic for pair of expressions in DataFrame
#  @code
#  frame  = ...
#  stat1 , stat2 , cov2 , len = frame.statCov( 'x' , 'y' )
#  # apply some cuts 
#  stat1 , stat2 , cov2 , len = frame.statCov( 'x' , 'y' , 'z>0' )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2018-06-18
def frame_statCov ( frame       ,
                    expression1 ,
                    expression2 ,
                    cuts = ''   ) :
    """ Get the statistic for pair of expressions in DataFrame
    
    >>>  frame  = ...
    >>>  stat1 , stat2 , cov2 , len = framw.statCov( 'x' , 'y' )
    
    Apply some cuts:
    >>> stat1 , stat2 , cov2 , len = frame.statCov( 'x' , 'y' , 'z>0' )
    
    """
    if isinstance ( frame  , ROOT.TTree ) : frame = DataFrame ( frame  )
    node = as_rnode ( frame ) 
    return SV.data_statCov  ( node , expression1 , expression2 , cuts ) 

# ==============================================================================
## Get the moment for the frame object
#  @code
#  frame = ...
#  value = frame_get_moment ( 3 , 1.234 , 'b1*b2' , 'b3>0' )
#  @endcode
#  @see Ostap::StatVar::get_moment 
def frame_get_moment ( frame , order , center , expression , cuts = '' ) :
    """ Get the moment  for the frame object
    >>> frame = ...
    >>> value = frame_get_moment ( 3 , 1.234 , 'b1*b2' , 'b3>0' )
    - see `Ostap.StatVar.get_moment` 
    """
    if isinstance ( frame  , ROOT.TTree ) : frame = DataFrame ( frame  )
    node = as_rnode ( frame )
    return SV.data_get_moment ( node , order = order , center = center , expression = expression , cuts = cuts )

# ==============================================================================
## Get the moment  for the frame object
#  @code
#  frame = ...
#  value = frame_moment ( 3 , 'b1*b2' , 'b3>0' )
#  @endcode
#  @see Ostap::StatVar::moment 
def frame_moment ( frame , order , expression , cuts = '' ) :
    """ Get the moment  for the frame object
    >>> frame = ...
    >>> value = frame_moment ( 3 , 'b1*b2' , 'b3>0' )
    - see `Ostap.StatVar.moment` 
    """
    node = as_rnode ( frame )
    return SV.data_moment ( node , order = order , expression = expression , cuts = cuts )

# ==============================================================================
## Get the central moment  for the frame object
#  @code
#  frame = ...
#  value = frame_central_moment ( 3 , 'b1*b2' , 'b3>0' )
#  @endcode
#  @see Ostap::StatVar::central_moment 
def frame_central_moment ( frame , order , expression , cuts = '' ) :
    """ Get the central_moment  for the frame object
    >>> frame = ...
    >>> value = frame_central_moment ( 3 , 'b1*b2' , 'b3>0' )
    - see `Ostap.StatVar.central_moment` 
    """
    if isinstance ( frame  , ROOT.TTree ) : frame = DataFrame ( frame  )
    node = as_rnode ( frame )
    return SV.data_central_moment ( node , order = order , expression = expression , cuts = cuts )


# ==============================================================================
## Get the mean for the frame object
#  @code
#  frame = ...
#  value = frame_mean ( 'b1*b2' , 'b3>0' )
#  @endcode
#  @see Ostap::StatVar::moment
def frame_mean ( frame , expression , cuts = '' ) :
    """ Get the mean for the frame object
    >>> frame = ...
    >>> value = frame_mean ( 'b1*b2' , 'b3>0' )
    - see `Ostap.StatVar.moment`
    """    
    return frame_moment ( frame , order = 1 , expression = expression , cuts = cuts )

# ==============================================================================
## Get the variance for the frame object
#  @code
#  frame = ...
#  value = frame_variance   ( 'b1*b2' , 'b3>0' , exact = True )
#  value = frame_dispersion ( 'b1*b2' , 'b3>0' , exact = True ) ## ditto
#  @endcode
#  @see Ostap::StatVar::central_moment 
def frame_variance ( frame , expression , cuts = '' ) :
    """ Get the variance for the frame object
    >>> frame = ...
    >>> value = frame_variance   ( 'b1*b2' , 'b3>0' )
    >>> value = frame_dispersion ( 'b1*b2' , 'b3>0' ) ## ditto
    - see `Ostap.StatVar.moment`
    """    
    return frame_central_moment ( frame , order = 2 , expression = expression , cuts = cuts )

# ==============================================================================
## dispersion, same as variance
frame_dispersion = frame_variance

# ==============================================================================
## Get the RMS for the frame object
#  @code
#  frame = ...
#  value = frame_rms ( 'b1*b2' , 'b3>0' , exact = True )
#  @endcode
#  @see Ostap::StatVar::central_moment 
def frame_rms ( frame , expression , cuts = '' ) :
    """ Get the RMS for the frame object
    >>> frame = ...
    >>> value = frame_rms ( 'b1*b2' , 'b3>0' , exact = True  )
    >>> value = frame_rms ( 'b1*b2' , 'b3>0' , exact = False ) ## use P2 algorithm 
    - see `Ostap.StatVar.moment`
    """    
    return frame_variance ( frame , expression = expression , cuts = cuts ) ** 0.5 

# ==============================================================================
## Get the skewness for the frame object
#  @code
#  frame = ...
#  value = frame_skewness ( 'b1*b2' , 'b3>0' )
#  @endcode
#  @see Ostap::StatVar::skewness 
def frame_skewness ( frame , expression , cuts = '' ) :
    """ Get the skewness for the frame object
    >>> frame = ...
    >>> value = frame_skeness ( 'b1*b2' , 'b3>0' )
    - see `Ostap.StatVar.skewness`
    """
    node = as_rnode ( frame )
    return SV.data_skewness ( node , expression = expression , cuts = cuts )

# ==============================================================================
## Get the kurtosis for the frame object
#  @code
#  frame = ...
#  value = frame_kurtosis ( 'b1*b2' , 'b3>0' )
#  @endcode
#  @see Ostap::StatVar::kurtosis 
def frame_kurtosis ( frame , expression , cuts = '' ) :
    """ Get the kurtosis for the frame object
    >>> frame = ...
    >>> value = frame_kurtosis ( 'b1*b2' , 'b3>0' )
    - see `Ostap.StatVar.kurtosis`
    """
    node = as_rnode ( frame )
    return SV.data_kurtosis ( node , expression = expression , cuts = cuts )

# ==============================================================================
## Get the quantile for the frame object
#  @code
#  frame = ...
#  value = frame_quantile (  0.3 , 'b1*b2' , 'b3>0' , exact = True )
#  value = frame_quantile (  0.3 , 'b1*b2' , 'b3>0' , exact = False  ) ## use P2 algorithm 
#  @endcode
#  @see Ostap::StatVar::quantile
#  @see Ostap::StatVar::p2quantile
def frame_quantile ( frame , q , expression , cuts = '' , exact = True ) :
    """ Get the quantile for the frame object
    >>> frame = ...
    >>> value = frame_quantile ( 0.3 , 'b1*b2' , 'b3>0' , exact = True  )
    >>> value = frame_quantile ( 0.3 , 'b1*b2' , 'b3>0' , exact = False ) ## use P2 algorithm 
    - see `Ostap.StatVar.quantile`
    - see `Ostap.StatVar.p2quantile`
    """
    node = as_rnode ( frame )
    return SV.data_quantile ( node , q = q , expression = expression , cuts = cuts , exact = exact  )

# ==============================================================================
## Get the interval  for the frame object
#  @code
#  frame = ...
#  value = frame_interval (  0.1 , 0.9 , 'b1*b2' , 'b3>0' )
#  @endcode
#  @see Ostap::StatVar::interval
#  @see Ostap::StatVar::p2interval 
def frame_interval ( frame , qmin , qmax  , expression , cuts = '' , exact = True ) :
    """ Get the approximarte quantile for the frame object usnig P2-algorithm
    >>> frame = ...
    >>> value = frame_p2quantile ( 0.3 , 'b1*b2' , 'b3>0' )
    - see `Ostap.StatVar.p2quantile`
    """
    node = as_rnode ( frame )
    return SV.data_interval ( node , qmin = qmin , qmax = qmax  ,
                              expression = expression , cuts = cuts , exact = exact )

# ==============================================================================
## Get the median
#  @code
#  frame = ...
#  value = frame_median ( 'b1*b2' , 'b3>0' , exact = True )
#  value = frame_median ( 'b1*b2' , 'b3>0' , exact = False  ) ## use P2 algorithm 
#  @endcode
#  @see Ostap::StatVar::quantile
#  @see Ostap::StatVar::p2quantile
def frame_median ( frame , expression , cuts = '' , exact = True ) :
    """ Get the quantile for the frame object
    >>> frame = ...
    >>> value = frame_quantile ( 0.3 , 'b1*b2' , 'b3>0' , exact = True  )
    >>> value = frame_quantile ( 0.3 , 'b1*b2' , 'b3>0' , exact = False ) ## use P2 algorithm 
    - see `Ostap.StatVar.quantile`
    - see `Ostap.StatVar.p2quantile`
    """
    return frame_quantile ( frame , q = 0.5 , expression = expression , cuts = cuts , exact = exact  )

# ==============================================================================
## Get the quantiles for the frame object
#  @code
#  frame = ...
#  value = frame_quantiles (  ( 0.1 , 0.2 , 0.3 ) , 'b1*b2' , 'b3>0' , exact = True )
#  value = frame_quantiles (  ( 0.1 , 0.2 , 0.3 ) , 'b1*b2' , 'b3>0' , exact = False  ) ## use P2 algorithm 
#  @endcode
#  @see Ostap::StatVar::quantiles
#  @see Ostap::StatVar::p2quantiles
def frame_quantiles ( frame , quantiles , expression , cuts = '' , exact = True ) :
    """ Get the quantile for the frame object
    >>> frame = ...
    >>> value = frame_quantiles ( (0.1,0.2,0.3) , 'b1*b2' , 'b3>0' , exact = True  )
    >>> value = frame_quantiles ( (0.1,0.2 0.3) , 'b1*b2' , 'b3>0' , exact = False ) ## use P2 algorithm 
    - see `Ostap.StatVar.quantiles`
    - see `Ostap.StatVar.p2quantiles`
    """
    node = as_rnode ( frame )
    return SV.data_quantiles ( node , quantiles = quantiles , expression = expression , cuts = cuts , exact = exact  )

# ==============================================================================
## Get the terciles for the frame object
#  @code
#  frame = ...
#  value = frame_terciles ( 'b1*b2' , 'b3>0' , exact = True )
#  value = frame_terciles ( 'b1*b2' , 'b3>0' , exact = False  ) ## use P2 algorithm 
#  @endcode
#  @see Ostap::StatVar::quantiles
#  @see Ostap::StatVar::p2quantiles
def frame_terciles ( frame , expression , cuts = '' , exact = True ) :
    """ Get the quantile for the frame object
    >>> frame = ...
    >>> value = frame_terciles ( 'b1*b2' , 'b3>0' , exact = True  )
    >>> value = frame_terciles ( 'b1*b2' , 'b3>0' , exact = False ) ## use P2 algorithm 
    - see `Ostap.StatVar.quantiles`
    - see `Ostap.StatVar.p2quantiles`
    """
    return frame_quantiles ( frame , quantiles = 3 , expression = expression , cuts = cuts , exact = exact  )

# ==============================================================================
## Get the quartiles for the frame object
#  @code
#  frame = ...
#  value = frame_quartiles ( 'b1*b2' , 'b3>0' , exact = True )
#  value = frame_quartiles ( 'b1*b2' , 'b3>0' , exact = False  ) ## use P2 algorithm 
#  @endcode
#  @see Ostap::StatVar::quantiles
#  @see Ostap::StatVar::p2quantiles
def frame_quartiles ( frame , expression , cuts = '' , exact = True ) :
    """ Get the quantile for the frame object
    >>> frame = ...
    >>> value = frame_quartiles ( 'b1*b2' , 'b3>0' , exact = True  )
    >>> value = frame_quartiles ( 'b1*b2' , 'b3>0' , exact = False ) ## use P2 algorithm 
    - see `Ostap.StatVar.quantiles`
    - see `Ostap.StatVar.p2quantiles`
    """
    return frame_quantiles ( frame , quantiles = 4 , expression = expression , cuts = cuts , exact = exact  )

# ==============================================================================
## Get the quintiles for the frame object
#  @code
#  frame = ...
#  value = frame_quintiles ( 'b1*b2' , 'b3>0' , exact = True )
#  value = frame_quintiles ( 'b1*b2' , 'b3>0' , exact = False  ) ## use P2 algorithm 
#  @endcode
#  @see Ostap::StatVar::quantiles
#  @see Ostap::StatVar::p2quantiles
def frame_quintiles  ( frame , expression , cuts = '' , exact = True ) :
    """ Get the quintiles for the frame object
    >>> frame = ...
    >>> value = frame_quintiles ( 'b1*b2' , 'b3>0' , exact = True  )
    >>> value = frame_quintiles ( 'b1*b2' , 'b3>0' , exact = False ) ## use P2 algorithm 
    - see `Ostap.StatVar.quantiles`
    - see `Ostap.StatVar.p2quantiles`
    """
    return frame_quantiles ( frame , quantiles = 5 , expression = expression , cuts = cuts , exact = exact  )

# ==============================================================================
## Get the deciles for the frame object
#  @code
#  frame = ...
#  value = frame_deciles ( 'b1*b2' , 'b3>0' , exact = True )
#  value = frame_deciles ( 'b1*b2' , 'b3>0' , exact = False  ) ## use P2 algorithm 
#  @endcode
#  @see Ostap::StatVar::quantiles
#  @see Ostap::StatVar::p2quantiles
def frame_deciles  ( frame , expression , cuts = '' , exact = True ) :
    """ Get the quintiles for the frame object
    >>> frame = ...
    >>> value = frame_deciles ( 'b1*b2' , 'b3>0' , exact = True  )
    >>> value = frame_deciles ( 'b1*b2' , 'b3>0' , exact = False ) ## use P2 algorithm 
    - see `Ostap.StatVar.quantiles`
    - see `Ostap.StatVar.p2quantiles`
    """
    return frame_quantiles ( frame , quantiles = 10 , expression = expression , cuts = cuts , exact = exact  )

# ==================================================================================

# ==================================================================================
## Action-based methods
# ==================================================================================

# ==================================================================================
## get statistics of variable(s)
#  @code
#  frame = ....
#  stat  = frame_the_tatVar  ( 'pt'           , lazy = True )
#  stat  = frame_the_statVar ( 'pt' , 'eta>0' , lazy = True )
#  @endcode
def _fr_the_statVar_ ( frame            ,
                       expressions      ,
                       cuts = ''        ,
                       progress = False ,
                       report   = False , 
                       lazy     = True  ) :
    """ Get statistics of variable(s)
    >>> frame = ....
    >>> stat  = frame_the_statVar ( 'pt'           , lazy = True )
    >>> stat  = frame_the_statVar ( 'pt' , 'eta>0' , lazy = True )
    """
    
    def screator ( node , var_name , cut_name ) : 
        if   cut_name: return node.Book( ROOT.std.move ( Ostap.Actions.WStatVar() ) , CNT ( [ var_name , cut_name ] ) )
        else         : return node.Book( ROOT.std.move ( Ostap.Actions. StatVar() ) , CNT ( 1 , var_name ) )     
        
    return _fr_helper2_ ( frame               ,
                          screator            ,
                          expressions         ,
                          cuts                ,
                          progress = progress ,
                          report   = report   ,  
                          lazy     = lazy     )


# ==================================================================================
## get covariance for variables
#  @code
#  frame = ....
#  stat  = frame_the_statCov ( 'pt1' , 'pt2' ,          lazy = True )
#  stat  = frame_the_statCov ( 'pt1' , 'pt2' , 'z<0' ,  lazy = True )
#  @endcode
def _fr_the_statCov_ ( frame            ,
                       expression1      ,
                       expression2      ,
                       cuts     = ''    ,
                       progress = False ,
                       rreport  = False ,                        
                       lazy     = True  ) :
    """ Get statistics of variable(s)
    >>> frame = ....
    >>> stat  = frame_the_statCov ( 'pt1' , 'pt2' ,          lazy = True )
    >>> stat  = frame_the_statCov ( 'pt1' , 'pt2' , 'z<0' ,  lazy = True )
    """
    
    current , exprs1 , cut_name , input_string1 = _fr_helper_ ( frame   , expression1 , cuts     , progress = progress )
    assert exprs1 and input_string1, 'Invalid expression1: %s' % expression2
    
    current , exprs2 , cut_name , input_string2 = _fr_helper_ ( current , expression2 , cut_name )
    assert exprs2 and input_string2, 'Invalid expression2: %s' % expression2 
    
    ## variables 
    uvars = [ k for k in exprs1.values() ] + [ k for k in exprs2.values() ]  
    if cut_name  : uvars.append ( cut_name )
    uvars = CNT ( uvars )
    
    if   cut_name:
        TT = ROOT.Detail.RDF.Stat3Action ( Ostap.Math.WCovariance ) 
        action = node.Book( ROOT.std.move ( TT () ) , CNT ( uvars ) )
    else         :
        TT = ROOT.Detail.RDF.Stat2Action ( Ostap.Math.Covariance  ) 
        action = node.Book( ROOT.std.move ( TT () ) , CNT ( uvars ) )     

    result = action if lazy else action.GetValue()

    if report and not lazy :
         report = current.Report()
         title = 'DataFrame processing'
         logger.info ( '%s\n%s' % ( title , report_print ( report , title = title , prefix = '# ') ) )
             
    return result 

# ==================================================================================
## get statistics of variable(s)
#  @code
#  frame = ....
#  stat  = frame.statVar ( 'pt'           , lazy = True )
#  stat  = frame.statVar ( 'pt' , 'eta>0' , lazy = True )
#  @endcode
def _fr_the_statVar_ ( frame            ,
                       expressions      ,
                       cuts     = ''    ,
                       progress = False ,
                       report   = False , 
                       lazy     = True  ) :
    """ Get statistics of variable(s)
    >>> frame = ....
    >>> stat  = frame.statVar ( 'pt'           , lazy = True )
    >>> stat  = frame.statVar ( 'pt' , 'eta>0' , lazy = True )
    """
    
    def screator ( node , var_name , cut_name ) : 
        if   cut_name: return node.Book( ROOT.std.move ( Ostap.Actions.WStatVar() ) , CNT ( [ var_name , cut_name ] ) )
        else         : return node.Book( ROOT.std.move ( Ostap.Actions. StatVar() ) , CNT ( 1 , var_name ) )     
        
    return _fr_helper2_ ( frame               ,
                          screator            ,
                          expressions         ,
                          cuts     = cuts     ,
                          progress = progress ,
                          report   = report   ,
                          lazy     = lazy     )

# ==================================================================================
## Get ECDF for the given varibale or expression 
#  @code
#  frame = ....
#  ecdf  = frame.get_ECDF ( 'pt'           , lazy = True )
#  ecdf  = frame.get_ECDF ( 'pt' , 'eta>0' , lazy = True )
#  @endcode
def _fr_the_ECDF_ ( frame            ,
                       expressions      ,
                       cuts     = ''    ,
                       progress = False ,
                       report   = False , 
                       lazy     = True  ) :
    """ Get ECDF for the given varibale or expression 
    >>> frame = ....
    >>> ecdf  = frame.get_ECDF ( 'pt'           , lazy = True )
    >>> ecdf  = frame.get_ECDF ( 'pt' , 'eta>0' , lazy = True )
    """
    
    def screator ( node , var_name , cut_name ) : 
        if   cut_name: return node.Book( ROOT.std.move ( Ostap.Actions.WECDF() ) , CNT ( [ var_name , cut_name ] ) )
        else         : return node.Book( ROOT.std.move ( Ostap.Actions. ECDF() ) , CNT ( 1 , var_name ) )     
        
    return _fr_helper2_ ( frame               ,
                          screator            ,
                          expressions         ,
                          cuts     = cuts     ,
                          progress = progress ,
                          report   = report   ,
                          lazy     = lazy     )

# ==================================================================================
## get nEff through action
#  @code
#  frame = ...
#  nEff = frame_the_nEff ( 'x*x' , 'y>0' )  
#  @endcode 
def _fr_the_nEff_ ( frame            ,
                    cuts     = ''    ,
                    progress = False ,
                    report   = False ,
                    lazy     = True  ) :
    """ Get nEff through action
    >>> frame = ...
    >>> nEff = frame_the_nEff ( 'x*x' , 'y>0' )
    """
    return _fr_the_statVar_ ( frame               ,
                              '1.0'               ,
                              cuts      = cuts    ,
                              progress = progress ,
                              report   = report   ,
                              lazy     = lazy     )

# ==================================================================================
## get statistics of variable(s)
#  @code
#  frame = ....
#  stat  = frame_the_moment ( 5 , 'pt'           )
#  stat  = frame_the_moment ( 5 , 'pt' , 'eta>0' )
#  @endcode
#  @see Ostap::Math::Moment_
#  @see Ostap::Math::WMoment_
def _fr_the_moment_ ( frame , N        ,
                      expressions      ,
                      cuts     = ''    ,
                      progress = False ,
                      report   = False ,
                      lazy     = True  ) :
    """ Get statistics of variable(s)
    >>> frame = ....
    >>> stat  = frame_the_moment ( 5 , 'pt'           )
    >>> stat  = frame_the_moment ( 5 , 'pt' , 'eta>0' )
    """
    assert isinstance ( N , integer_types ) and 0 <= N , 'Invalid order!'

    current, vname , cname , input_string = _fr_helper_ ( frame , expressions , cuts ) 
    
    def mcreator ( node , var_name , cut_name ) : 
        if cname :
            TT = ROOT.Detail.RDF.Stat2Action ( Ostap.Math.WMoment_(N) ) 
            return current.Book ( ROOT.std.move ( TT () ) , CNT ( [ var_name , cut_name ] ) ) 
        else :
            TT = ROOT.Detail.RDF.Stat1Action ( Ostap.Math.Moment_(N) ) 
            return current.Book ( ROOT.std.move ( TT () ) , CNT ( 1 , var_name ) ) 

    return _fr_helper2_ ( frame               , 
                          mcreator            ,
                          expressions         ,
                          cuts     = cuts     ,
                          progress = progress ,
                          report   = report   ,
                          lazy     = lazy     )

# ============================================================================
## Get a mean via moment
#  @code
#  frame = ...
#  stat  = frame.the_mean ( 'x'x' ,  'y<0' )
#  mean  = stat.mean() 
#  @ndcode
def _fr_the_mean_ ( frame            ,
                    expressions      ,
                    cuts     = ''    ,
                    errors   = True  , 
                    progress = False ,
                    report   = False ,
                    lazy     = True  ) :
    """ Get a mean via moment
    >>> frame = ...
    >>> stat  = frame.the_mean ( 'x'x' ,  'y<0' )
    >>> mean  = stat.mean() 
    """
    return _fr_the_moment_ ( frame               ,
                             2 if errors else 1  ,
                             expressions         ,
                             cuts     = cuts     ,
                             progress = progress ,
                             report   = report   ,
                             lazy     = lazy     )

# ============================================================================
## Get a rms via moment
#  @code
#  frame = ...
#  stat  = frame.the_rms ( 'x'x' ,  'y<0' )
#  value = stat.rms() 
#  @ndcode
def _fr_the_rms_ ( frame ,
                   expressions      ,
                   cuts     = ''    ,
                   errors   = True  ,                    
                   progress = False ,
                   report   = False ,
                   lazy     = True  ) :
    """ Get a rms via moment
    >>> frame = ...
    >>> stat  = frame.the_rms ( 'x'x' ,  'y<0' )
    >>> value = stat.rms () 
    """
    return _fr_the_moment_ ( frame               ,
                             4 if errors else 2  , 
                             expressions         ,
                             cuts     = cuts     ,
                             progress = progress ,
                             report   = report   ,
                             lazy     = lazy     )

# ============================================================================
## Get a variance via moment
#  @code
#  frame = ...
#  stat  = frame.the_variance ( 'x'x' ,  'y<0' )
#  value = stat.variance () 
#  @ndcode
def _fr_the_variance_ ( frame            , 
                        expressions      ,
                        cuts     = ''    ,
                        errors   = True  ,                    
                        progress = False ,
                        report   = False ,
                        lazy     = True  ) :
    """ Get a variance via moment
    >>> frame  = ...
    >>> stat   = frame.the_variance( 'x'x' ,  'y<0' )
    >>> value  = stat.variance () 
    """
    return _fr_the_moment_ ( frame               ,
                             4 if errors else 2  , 
                             expressions         ,
                             cuts     = cuts     ,
                             progress = progress ,
                             report   = report   ,
                             lazy     = lazy     )

# ============================================================================
## Get a skewness via moment
#  @code
#  frame = ...
#  stat  = frame.the_skewness ( 'x'x' ,  'y<0' )
#  mean  = stat.skewness  () 
#  @ndcode
def _fr_the_skewness_ ( frame ,
                        expressions      ,
                        cuts     = ''    ,
                        errors   = True  ,                    
                        progress = False ,
                        report   = False ,
                        lazy     = True  ) :
    """ Get a skewness via moment
    >>> frame = ...
    >>> stat  = frame.the_skewness ( 'x'x' ,  'y<0' )
    >>> mean  = stat.skreness  () 
    """
    return _fr_the_moment_ ( frame               ,
                             6 if errors else 3  , 
                             expressions         ,
                             cuts     = cuts     ,
                             progress = progress ,
                             report   = report   ,
                             lazy     = lazy     )

# ============================================================================
## Get a kurtosis via moment
#  @code
#  frame = ...
#  stat  = frame.the_kurtosis ( 'x'x' ,  'y<0' )
#  value = stat.kurtosis   () 
#  @ndcode
def _fr_the_kurtosis_ ( frame            ,
                        expressions      ,
                        cuts     = ''    ,
                        errors   = True  ,                    
                        progress = False ,
                        report   = False ,
                        lazy     = True  ) :
    """ Get a kurtosis  via moment
    >>> frame = ...
    >>> stat  = frame.the_kurtosis ( 'x'x' ,  'y<0' )
    >>> value = stat.kurtosis   () 
    """
    return _fr_the_moment_ ( frame               ,
                             8 if errors else 4  ,
                             expressions         ,
                             cuts     = cuts     ,
                             progress = progress ,
                             report   = report   ,
                             lazy     = lazy     )

# ============================================================================
## Get an arithmetic mean
#  @see Ostap::Math::ArithmeticMean
#  @code
#  frame = ...
#  mean = frame_the_arithmetic_mean ( frame , 'x*x' , '0<y' ) 
#  @endcode 
def _fr_the_arithmetic_mean_ ( frame ,
                               expressions      ,
                               cuts     = ''    ,
                               progress = False ,
                               report   = False ,
                               lazy     = True  ) :
    """ Get an arithmetic mean
    >>> frame = ...
    >>> mean = frame_the_arithmetic_mean ( frame , 'x*x' , '0<y' ) 
    -  see `Ostap.Math.ArithmeticMean`
    """
    
    def acreator ( node , var_name , cut_name ) : 
        if cut_name: 
            TT = ROOT.Detail.RDF.Stat2Action ( Ostap.Math.WArithmeticMean ) 
            return node.Book( ROOT.std.move ( TT () ) , CNT ( [ var_name , cut_name ] ) )
        else  :
            TT = ROOT.Detail.RDF.Stat1Action ( Ostap.Math.ArithmeticMean ) 
            return node.Book( ROOT.std.move ( TT () ) , CNT ( 1 , var_name ) )     
        
    return _fr_helper2_ ( frame               ,
                          acreator            , 
                          expressions         ,
                          cuts     = cuts     ,
                          progress = progress ,
                          report   = report   ,
                          lazy     = lazy     )
    
# ============================================================================
## Get a geometric mean
#  @see Ostap::Math::GeometricMean
#  @code
#  frame = ...
#  mean = frame_the_geometric_mean ( frame , 'x*x' , '0<y' ) 
#  @endcode 
def _fr_the_geometric_mean_ ( frame ,
                              expressions      ,
                              cuts     = ''    ,
                              progress = False ,
                              report   = False ,
                              lazy     = True  ) :
    """ Get an geometric mean
    >>> frame = ...
    >>> mean = frame_the_geometric_mean ( frame , 'x*x' , '0<y' ) 
    -  see `Ostap.Math.GeometricMean`
    """
    def gcreator ( node , var_name , cut_name ) : 
        if cut_name: 
            TT = ROOT.Detail.RDF.Stat2Action ( Ostap.Math.WGeometricMean ) 
            return node.Book( ROOT.std.move ( TT () ) , CNT ( [ var_name , cut_name ] ) )
        else  :
            TT = ROOT.Detail.RDF.Stat1Action ( Ostap.Math.GeometricMean ) 
            return node.Book( ROOT.std.move ( TT () ) , CNT ( 1 , var_name ) )     
        
    return _fr_helper2_ ( frame               ,
                          gcreator            ,
                          expressions         ,
                          cuts     = cuts     ,
                          progress = progress ,
                          report   = report   ,
                          lazy     = lazy     )

# ============================================================================
## Get a harmonic mean
#  @see Ostap::Math::HarmonicMean
#  @code
#  frame = ...
#  mean = frame_the_harmonic_mean ( frame , 'x*x' , '0<y' ) 
#  @endcode 
def _fr_the_harmonic_mean_ ( frame            , 
                             expressions      ,
                             cuts     = ''    ,
                             progress = False ,
                             report   = False ,
                             lazy     = True  ) :
    """ Get a harmonic mean
    >>> frame = ...
    >>> mean = frame_the_harmonic_mean ( frame , 'x*x' , '0<y' ) 
    -  see `Ostap.Math.HarmonicMean`
    """
    def creator ( node , var_name , cut_name ) : 
        if cut_name: 
            TT = ROOT.Detail.RDF.Stat2Action ( Ostap.Math.WHarmonicMean ) 
            return node.Book( ROOT.std.move ( TT () ) , CNT ( [ var_name , cut_name ] ) )
        else  :
            TT = ROOT.Detail.RDF.Stat1Action ( Ostap.Math.HarmonicMean ) 
            return node.Book( ROOT.std.move ( TT () ) , CNT ( 1 , var_name ) )     
        
    return _fr_helper2_ ( frame               ,
                          creator             ,
                          expressions         ,
                          cuts     = cuts     ,
                          progress = progress ,
                          report   = report   ,
                          lazy     = lazy     )
 
# ============================================================================
## Get a power mean
#  @see Ostap::Math::PowerMean
#  - p == -1 : harmonic mean 
#  - p ==  0 : geometric mean
#  - p == +1 : arithmetic mean 
#  @code
#  frame = ...
#  mean = frame_the_power_mean ( frame , 3 , 'x*x' , '0<y' ) 
#  @endcode 
def _fr_the_power_mean_ ( frame , p ,
                          expressions      ,
                          cuts     = ''    ,
                          progress = False ,
                          report   = False ,
                          lazy     = True  ) :
    """ Get a power mean
    >>> frame = ...
    >>> mean = frame_the_power_mean ( frame , 5 ,  'x*x' , '0<y' ) 
    -  see `Ostap.Math.PowerMean`
   - p == -1 : harmonic mean 
   - p ==  0 : geometric mean
   - p == +1 : arithmetic mean 
   """
    assert isinstance ( p , num_types ) , 'Invalid type of p-parameter: %s' % type ( p ) 
    
    if    -1 == p or isequal ( p , -1.0 ) : 
        return _fr_harmonic_mean_   ( frame , expressions , cuts = cuts , lazy = lazy ) 
    elif   0 == p or iszero  ( p        ) : 
        return _fr_geometric_mean_  ( frame , expressions , cuts = cuts , lazy = lazy ) 
    elif   1 == p or isequal ( p ,  1.0 ) : 
        return _fr_arithmetic_mean_ ( frame , expressions , cuts = cuts , lazy = lazy ) 
    
    def pcreator ( node , var_name , cut_name ) : 
        if cut_name: 
            TT = ROOT.Detail.RDF.Stat2Action ( Ostap.Math.WPowerMean ) 
            return node.Book( ROOT.std.move ( TT ( p ) ) , CNT ( [ var_name , cut_name ] ) )
        else  :
            TT = ROOT.Detail.RDF.Stat1Action ( Ostap.Math.PowerMean ) 
            return node.Book( ROOT.std.move ( TT ( p ) ) , CNT ( 1 , var_name ) )
        
    return _fr_helper2_ ( frame               ,
                          pcreator            ,
                          expressions         ,
                          cuts     = cuts     ,
                          progress = progress ,
                          report   = report   ,
                          lazy     = lazy     )

# ============================================================================
## Get Lehmer mean
#  @see Ostap::Math::LehmerMean
#  - p ==  0 : harmonic mean 
#  - p == +1 : arithmetic mean 
#  @code
#  frame = ...
#  mean = frame_the_lehmer_mean ( frame , 3 , 'x*x' , '0<y' ) 
#  @endcode 
def _fr_the_lehmer_mean_ ( frame , p ,
                           expressions      ,
                           cuts     = ''    ,
                           progress = False ,
                           report   = False ,
                           lazy     = True  ) :
    """ Get Lehmer mean
    >>> frame = ...
    >>> mean = frame_the_lehmer_mean ( frame , 5 ,  'x*x' , '0<y' ) 
    -  see `Ostap.Math.LehmerMean`
   - p ==  0 : harmonic mean 
   - p == +1 : arithmetic mean 
   """
    assert isinstance ( p , num_types ) , 'Invalid type of p-parameter: %s' % type ( p ) 
    
    if   0 == p or iszero  ( p        ) : 
        return _fr_harmonic_mean_   ( frame , expressions , cuts = cuts , lazy = lazy ) 
    elif 1 == p or isequal ( p ,  1.0 ) : 
        return _fr_arithmetic_mean_ ( frame , expressions , cuts = cuts , lazy = lazy ) 
    
    def lcreator ( node , var_name , cut_name ) : 
        if cut_name: 
            TT = ROOT.Detail.RDF.Stat2Action ( Ostap.Math.WLehmerMean ) 
            return node.Book( ROOT.std.move ( TT ( p ) ) , CNT ( [ var_name , cut_name ] ) )
        else  :
            TT = ROOT.Detail.RDF.Stat1Action ( Ostap.Math.LehmerMean ) 
            return node.Book( ROOT.std.move ( TT ( p ) ) , CNT ( 1 , var_name ) )     
        
    return _fr_helper2_ ( frame               ,
                          lcreator            ,
                          expressions         ,
                          cuts     = cuts     ,
                          progress = progress ,
                          report   = report   ,
                          lazy     = lazy     )

# ============================================================================
## Get min/max for variables 
#  @see Ostap::Math::MinMaxValue
#  @see Ostap::Math::WMinMaxValue
#  @code
#  frame = ...
#  mean = frame_the_minmax ( frame , 'x*x' , '0<y' ) 
#  @endcode 
def _fr_the_minmax_ ( frame            ,
                      expressions      ,
                      cuts     = ''    ,
                      progress = False ,
                      report   = False ,
                      lazy     = True  ) :
    """ Get min/max values for variables
    >>> frame = ...
    >>> mean = frame_the_minmax ( frame , 5 ,  'x*x' , '0<y' ) 
    -  see `Ostap.Math.MinMaxValue`
    -  see `Ostap.Math.WMinMaxValue`
   """
    def mcreator ( node , var_name , cut_name ) : 
        if cut_name: 
            ## TT = ROOT.Detail.RDF.Stat2Action ( Ostap.Math.WMinMaxValue ) 
            TT = ROOT.Detail.RDF.Stat2Action ( Ostap.Math.WMoment_[1] ) 
            return node.Book( ROOT.std.move ( TT () ) , CNT ( [ var_name , cut_name ] ) )
        else  :
            ## TT = ROOT.Detail.RDF.Stat1Action ( Ostap.Math.MinMaxValue  ) 
            TT = ROOT.Detail.RDF.Stat2Action ( Ostap.Math.Moment_[1] ) 
            return node.Book( ROOT.std.move ( TT () ) , CNT ( 1 , var_name ) )     
        
    return _fr_helper2_ ( frame               ,
                          mcreator            ,
                          expressions         ,
                          cuts     = cuts     ,
                          progress = progress ,
                          report   = report   ,
                          lazy     = lazy     )


# ===============================================================================
if Frames_OK :
    
    ## frame_get_moment 
    ## frame_central_moment
    
    frame_the_statVar         = _fr_the_statVar_
    frame_the_statCov         = _fr_the_statCov_
    frame_the_nEff            = _fr_the_nEff_  
    frame_the_ECDF            = _fr_the_ECDF_  
    frame_the_moment          = _fr_the_moment_
    frame_the_mean            = _fr_the_mean_
    frame_the_minmax          = _fr_the_minmax_
    frame_the_rms             = _fr_the_rms_
    frame_the_variance        = _fr_the_variance_
    frame_the_dispersion      = _fr_the_variance_
    frame_the_skewness        = _fr_the_skewness_
    frame_the_kurtosis        = _fr_the_kurtosis_
    frame_the_arithmetic_mean = _fr_the_arithmetic_mean_ 
    frame_the_geometric_mean  = _fr_the_geometric_mean_ 
    frame_the_harmonic_mean   = _fr_the_harmonic_mean_ 
    frame_the_power_mean      = _fr_the_power_mean_ 
    frame_the_lehmer_mean     = _fr_the_lehmer_mean_ 
    frame_the_param           = _fr_param_

    __all__ += (
        'frame_the_statVar'         ,
        'frame_the_statCov'         ,
        'frame_the_nEff'            ,
        'frame_the_ECDF'            ,
        'frame_the_moment'          ,
        'frame_the_mean'            ,
        'frame_the_minmax'          ,
        'frame_the_rms'             ,
        'frame_the_variance'        ,
        'frame_the_dispersion'      ,
        'frame_the_skewness'        ,
        'frame_the_kurtosis'        ,
        'frame_the_arithmetic_mean' ,
        'frame_the_geometric_mean'  ,
        'frame_the_harmonic_mean'   ,
        'frame_the_power_mean'      ,
        'frame_the_lehmer_mean'     ,
        'frame_the_param'           ,
    )
    
    # ==================================================================================
    ## get statistics of variable(s) (via actions) 
    #  @code
    #  frame = ....
    #  stat  = frame_statVar ( frame , 'pt'           )
    #  stat  = frame_statVar ( frame , 'pt' , 'eta>0' )
    #  @endcode
    def frame_statVar ( frame , expressions , cuts = ''  , progress = False , report = False ) :
        """ Get statistics of variable(s)
        >>> frame = ....
        >>> stat  = frame.statVar ( 'pt'           )
        >>> stat  = frame.statVar ( 'pt' , 'eta>0' )
        """
        return frame_the_statVar ( frame               ,
                                   expressions         ,
                                   cuts     = cuts     ,
                                   progress = progress ,
                                   report   = report   ,
                                   lazy     = False    )

    # ==================================================================================
    ## get ECDF
    #  @code
    #  frame = ....
    #  stat  = frame_ECDF ( frame , 'pt'           )
    #  stat  = frame_ECDF ( frame , 'pt' , 'eta>0' )
    #  @endcode
    def frame_ECDF ( frame , expressions , cuts = ''  , progress = False , report = False ) :
        """ Get statistics of variable(s)
        >>> frame = ....
        >>> ecdf  = frame.ECDF ( 'pt'           )
        >>> ecdf  = frame.ECDF ( 'pt' , 'eta>0' )
        """
        return frame_the_ECDF ( frame               ,
                                expressions         ,
                                cuts     = cuts     ,
                                progress = progress ,
                                report   = report   ,
                                lazy     = False    )
    
    # ==================================================================================
    ## get covariance for variables via action 
    #  @code
    #  frame = ....
    #  stat  = frame_statCov ( 'pt1' , 'pt2' ,        )
    #  stat  = frame_statCov ( 'pt1' , 'pt2' , 'z<0'  )
    #  @endcode
    def frame_statCov ( frame            ,
                        expression1      ,
                        expression2      ,
                        cuts     = ''    ,
                        progress = False ,
                        report   = False ) :
        """ Get statistics of variables via actions 
        >>> frame = ....
        >>> stat  = frame_statCov ( 'pt1' , 'pt2' ,        )
        >>> stat  = frame_statCov ( 'pt1' , 'pt2' , 'z<0'  )
        """
        return frame_the_statCov ( frame               ,
                                   expression1         ,
                                   expression2         ,
                                   cuts     = cuts     , 
                                   progress = progress ,
                                   report   = report   ,
                                   lazy     = False    )

    # =============================================================================
    ## Get the effective entries in data frame (via actions)
    #  @code
    #  data = ...
    #  neff = data.nEff('b1*b1')
    #  @endcode
    def frame_nEff ( frame            ,
                     cuts     = ''    ,
                     progress = False ,
                     report   = False ) :
        """ Get the effective entries in data frame 
        >>> data = ...
        >>> neff = frame_nEff('b1*b1')
        """
        return frame_the_nEff ( frame               ,
                                cuts     = cuts     ,
                                progress = progress ,
                                report   = report   ,
                                lazy     = False    ).nEff()         
    # ==================================================================================
    ## get statistics of variable(s)  (via actions)
    #  @code
    #  frame = ....
    #  stat  = frame_the_moment ( 5 , 'pt'           )
    #  stat  = frame_the_moment ( 5 , 'pt' , 'eta>0' )
    #  @endcode
    #  @see Ostap::Math::Moment_
    #  @see Ostap::Math::WMoment_
    def frame_moment ( frame , N , expressions , cuts = '' , progress = False , report = False ) :
        """ Get statistics of variable(s)  (via actions) 
        >>> frame = ....
        >>> stat  = frame_moment ( 5 , 'pt'           )
        >>> stat  = frame_moment ( 5 , 'pt' , 'eta>0' )
        """
        return frame_the_moment ( frame , N           ,
                                  expressions         ,
                                  cuts     = cuts     ,
                                  progress = progress ,
                                  report   = report   , 
                                  lazy     = False    )
    # ============================================================================
    ## Get a mean via moment&action
    #  @code
    #  frame = ...
    #  mean  = frame_mean ( 'x'x' ,  'y<0' )
    #  @ndcode
    def frame_mean ( frame , expressions , cuts = '' , errors = True , progress = False , report = False ) :
        """ Get a mean via moment
        >>> frame = ...
        >>> mean = frame_mean ( 'x'x' ,  'y<0' )
        """
        return frame_the_mean ( frame               ,
                                expressions         ,
                                cuts     = cuts     ,
                                errors   = errors   ,
                                progress = progress ,
                                report   = report   ,
                                lazy     = False    ).mean() 
    # ============================================================================
    ## Get a rms via moment&action 
    #  @code
    #  frame = ...
    #  rms   = frame_rms ( 'x'x' ,  'y<0' )
    #  @ndcode
    def frame_rms ( frame , expressions , cuts = '' , errors = True , progress = False , report = False ) :
        """ Get a rms via moment&action
        >>> frame = ...
        >>> rms   = frame_rms ( 'x'x' ,  'y<0' )
        """
        return frame_the_rms ( frame               ,
                               expressions         ,
                               cuts     = cuts     ,
                               errors   = errors   ,
                               progress = progress ,
                               report   = report   ,
                               lazy     = False    ).rms() 
    # ============================================================================
    ## Get a variance via moment&actions
    #  @code
    #  frame     = ...
    #  variance  = frame_variance   ( 'x'x' ,  'y<0' )
    #  variance  = frame_dispersion ( 'x'x' ,  'y<0' )
    #  @ndcode
    def frame_variance ( frame , expressions , cuts = '' , errors = True , progress = False , report = False ) :
        """ Get a variance via moment&action
        >>> frame     = ...
        >>> variance  = frame_variance   ( 'x'x' ,  'y<0' )
        >>> variance  = frame_dispersion ( 'x'x' ,  'y<0' )        
        """
        return frame_the_variance ( frame               ,
                                    expressions         ,
                                    cuts     = cuts     ,
                                    errors   = errors   ,
                                    progress = progress ,
                                    report   = report   ,
                                    lazy     = False    ).variance()
    # ============================================================================
    ## ditto 
    frame_dispersion = frame_variance 
    # ============================================================================
    ## Get a skewness via moment&action 
    #  @code
    #  frame = ...
    #  skew  = frame_skewness ( 'x'x' ,  'y<0' )
    #  @ndcode
    def frame_skewness ( frame , expressions , cuts = '' , errors = True , progress = False , report = False ) :
        """ Get a variance via moment
        >>> frame = ...
        >>> skew  = frame_skewness ( 'x'x' ,  'y<0' )
        >>> mean  = stat.skreness  () 
        """
        return frame_the_skewness ( frame               ,
                                    expressions         ,
                                    cuts     = cuts     ,
                                    errors   = errors   ,
                                    progress = progress ,
                                    report   = report   ,                                     
                                    lazy     = False    ).skewness()
    # ============================================================================
    ## Get a kurtosis via moment&actions 
    #  @code
    #  frame = ...
    #  kurt  = frame_kurtosis ( 'x'x' ,  'y<0' )
    #  @ndcode
    def frame_kurtosis ( frame            ,
                         expressions      ,
                         cuts     = ''    ,
                         errors   = True  ,
                         progress = False ,
                         report   = False ) :
        """ Get a kurtosis  via moment
        >>> frame = ...
        >>> stat  = frame.the_kurtosis ( 'x'x' ,  'y<0' )
        >>> value = stat.kurtosis   () 
        """
        return frame_the_kurtosis ( frame               ,
                                    expressions         , 
                                    cuts     = cuts     ,
                                    errors   = True     ,
                                    progress = progress ,
                                    report   = report   , 
                                    lazy     = False    ).kurtosis()     
    # ============================================================================
    ## Get an arithmetic mean
    #  @see Ostap::Math::ArithmeticMean
    #  @code
    #  frame = ...
    #  mean = frame_arithmetic_mean ( frame , 'x*x' , '0<y' ) 
    #  @endcode 
    def frame_arithmetic_mean ( frame            ,
                                expressions      ,
                                cuts     = ''    ,  
                                progress = False ,
                                report   = False ) :
        """ Get an arithmetic mean
        >>> frame = ...
        >>> mean = frame_arithmetic_mean ( frame , 'x*x' , '0<y' ) 
        -  see `Ostap.Math.ArithmeticMean`
        """
        return frame_the_arithmetic_mean ( frame ,
                                           expressions         ,
                                           cuts     = cuts     ,
                                           progress = progress ,
                                           report   = report   , 
                                           lazy     = False    )
    # ============================================================================
    ## Get a geometric mean
    #  @see Ostap::Math::GeometricMean&actions 
    #  @code
    #  frame = ...
    #  mean = frame_geometric_mean ( frame , 'x*x' , '0<y' ) 
    #  @endcode 
    def frame_geometric_mean ( frame            ,
                               expressions      ,
                               cuts     = ''    ,  
                               progress = False ,
                               report   = False ) :
        """ Get an geomentric mean
        >>> frame = ...
        >>> mean = frame_geometric_mean ( frame , 'x*x' , '0<y' ) 
        -  see `Ostap.Math.GeometricMean`
        """
        return frame_the_geometric_mean ( frame               ,
                                          expressions         ,
                                          cuts     = cuts     ,
                                          progress = progress ,
                                          report   = report   , 
                                          lazy     = False    )
    # ============================================================================
    ## Get a harmonic mean
    #  @see Ostap::Math::HarmonicMean
    #  @code
    #  frame = ...
    #  mean = frame_harmonic_mean ( frame , 'x*x' , '0<y' ) 
    #  @endcode 
    def frame_harmonic_mean ( frame            ,
                              expressions      ,
                              cuts     = ''    ,
                              progress = False ,
                              report   = False ) :
        """ Get a harmonic mean
        >>> frame = ...
        >>> mean = frame_harmonic_mean ( frame , 'x*x' , '0<y' ) 
        -  see `Ostap.Math.HarmonicMean`
        """
        return frame_the_harmonic_mean ( frame ,
                                         expressions         ,
                                         cuts     = cuts     ,
                                         progress = progress ,
                                         report   = report   , 
                                         lazy     = False    )
    # ============================================================================
    ## Get a power mean
    #  @see Ostap::Math::PowerMean
    #  - p == -1 : harmonic mean 
    #  - p ==  0 : geometric mean
    #  - p == +1 : arithmetic mean 
    #  @code
    #  frame = ...
    #  mean = frame_power_mean ( frame , 3 , 'x*x' , '0<y' ) 
    #  @endcode 
    def frame_power_mean ( frame , p        ,
                           expressions      ,
                           cuts     = ''    ,
                           progress = False ,
                           report   = False ) :
        """ Get a power mean
        >>> frame = ...
        >>> mean = frame_power_mean ( frame , 5 ,  'x*x' , '0<y' ) 
        -  see `Ostap.Math.PowerMean`
        - p == -1 : harmonic mean 
        - p ==  0 : geometric mean
        - p == +1 : arithmetic mean 
        """
        return frame_the_power_mean ( frame , p ,
                                      expressions         ,
                                      cuts     = cuts     ,
                                      progress = progress ,
                                      report   = report   , 
                                      lazy     = False    )
    # ============================================================================
    ## Get Lehmer mean
    #  @see Ostap::Math::LehmerMean
    #  - p ==  0 : harmonic mean 
    #  - p == +1 : arithmetic mean 
    #  @code
    #  frame = ...
    #  mean = frame_lehmer_mean ( frame , 3 , 'x*x' , '0<y' ) 
    #  @endcode 
    def frame_lehmer_mean ( frame , p , expressions , cuts = '' , progress = False , report = False ) :
        """ Get Lehmer mean
        >>> frame = ...
        >>> mean  = frame_lehmer_mean ( frame , 5 ,  'x*x' , '0<y' ) 
        -  see `Ostap.Math.LehmerMean`
        - p ==  0 : harmonic mean 
        - p == +1 : arithmetic mean 
        """
        return frame_the_lehmer_mean ( frame , p  ,
                                       expressions         ,
                                       cuts     = cuts     ,
                                       progress = progress ,
                                       report   = report   ,
                                       lazy     = False    )
    # ==============================================================================
    ## "project/parameterise" frame into polynomial structures 
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
    def frame_param ( frame            ,
                      poly             ,
                      expressions      ,
                      cuts     = ''    ,
                      progress = False ,
                      report   = False ) :
        """ `project/parameterise` frame into polynomial structures
        
        >>> frame = ...
        
        >>> ls    = Ostap.Math.LegendreSum ( ... )
        >>> res   = frame_param ( frame , ls , 'x' , 'y>0' )
        
        >>> bs    = Ostap.Math.Bernstein   ( ... )
        >>> res   = frame_param ( frame , bs , 'x' , 'y>0' )
        
        >>> cs    = Ostap.Math.ChebyshevSum ( ... )
        >>> res   = frame_param ( frame , cs , 'x' , 'y>0' )
        
        >>> ls2   = Ostap.Math.LegendreSum2 ( ... )
        >>> res   = frame_param ( frame , ls2 , 'y,x' , 'z>0' )
        
        >>> bs2   = Ostap.Math.Bernstein2D  ( ... )
        >>> res   = frame_param ( frame , bs2 , 'y,x' , 'z>0' )
        
        >>> ls3   = Ostap.Math.LegendreSum3 ( ... )
        >>> res   = frame_param ( frame , ls3 , 'z,y;z' , 'z>0' )
        
        >>> bs2   = Ostap.Math.Bernstein2D  ( ... )
        >>> res   = frame_param ( frame , bs2 , 'y,x' , 'z>0' )
        """
        return frame_the_param ( frame            ,
                                 poly             ,
                                 expressions      ,
                                 cuts     = cuts  ,
                                 progress = False ,
                                 report   = False , lazy = False ) 
    
    # ============================================================================
    ## Get the approproate range  for drawing the variables 
    #  In case there is no suitable range None is returned 
    #  @code
    #  frame = ...
    #  mean = frame_range ( frame , 'x*x' , '0<y' ) 
    #  @endcode 
    def frame_range ( frame            ,
                      expressions      ,
                      cuts     = ''    ,
                      delta    = 0.01  ,
                      progress = False ,
                      report   = False ) :
        """ Get the approproavte range  values for variables
         - In case there is no suitable range None is returned 
        >>> frame = ...
        >>> mean = frame_the_minmax ( frame , 5 ,  'x*x' , '0<y' ) 
        """
        ranges = frame_the_minmax ( frame               ,
                                    expressions         ,
                                    cuts     = cuts     ,
                                    progress = progress ,
                                    report   = report   ,
                                    lazy     = False    )
        
        if isinstance ( ranges , dictlike_types ) :
            for k , r in loop_items ( ranges ) :
                mn, mx = r.min () , r.max()
                if mx <= mn : return None                         ## ATTENTION!!
                ranges [ k ] = axis_range ( mn , mx , delta = delta ) 
        else : ranges = axis_range ( ranges.xmin() , ranges.xmax() , delta = delta )
        return ranges


    ## use more efficient drawing function 
    frame_draw = _fr_draw_
    
    __all__ += (
        'frame_arithmetic_mean' ,
        'frame_geometric_mean'  ,
        'frame_harmonic_mean'   ,
        'frame_power_mean'      ,
        'frame_lehmer_mean'     ,
        'frame_param'           ,
        'frame_range'           ,
    )
    

# ===============================================================================
frame_statistics = frame_statVar

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
# decorate 
# ==============================================================================
_new_methods_ = [
    ## 
    frame_columns        , 
    frame_branches       ,
    frame_progress       , 
    frame_prescale       , 
    frame_project        , 
    frame_draw           , 
    frame_print          , 
    frame_table          , 
    ##
    frame_length         , 
    frame_size           , 
    ##
    frame_nEff           , 
    frame_statVar        , 
    frame_statCov        ,
    #
    frame_get_moment     , 
    frame_moment         , 
    frame_central_moment ,
    # 
    frame_mean           , 
    frame_rms            , 
    frame_variance       , 
    frame_dispersion     , 
    frame_skewness       , 
    frame_kurtosis       , 
    ## 
    frame_quantile       ,
    frame_interval       ,
    frame_median         ,
    frame_quantiles      ,
    frame_terciles       ,
    frame_quartiles      ,
    frame_quintiles      ,
    frame_deciles        ,
    ##
    frame_statistics     , 
    ## 
    report_print         ,
    report_print_table   ,
    report_as_table      ,    
    ]     
# ================================================================================
for f in frame_types :

    f.__str__        = frame_print
    f.__repr__       = frame_print
    f.__len__        = frame_length

    f.columns        = frame_columns  
    f.branches       = frame_branches
    
    f.project        = frame_project 
    f.draw           = frame_draw
    
    f.nEff           = frame_nEff
    f.statVar        = frame_statVar 
    f.statCov        = frame_statCov
    
    f.get_moment     = frame_get_moment
    f.moment         = frame_moment
    f.central_moment = frame_moment
    #
    f.mean           = frame_mean
    f.rms            = frame_rms 
    f.variance       = frame_variance 
    f.dispersion     = frame_dispersion
    f.skewness       = frame_skewness
    f.kurtosis       = frame_kurtosis
    #
    f.quantile       = frame_quantile 
    f.interval       = frame_interval 
    f.median         = frame_median
    f.quantiles      = frame_quantiles
    f.terciles       = frame_terciles 
    f.quartiles      = frame_quartiles
    f.quintiles      = frame_quintiles    
    f.deciles        = frame_deciles 
    
    _new_methods_ += [
        f.columns        , 
        f.branches       ,         
        f.project        , 
        f.draw           ,
        #
        f.__str__        , 
        f.__repr__       , 
        f.__len__        ,
        #
        f.nEff           , 
        f.statVar        , 
        f.statCov        , 
        #
        f.get_moment     , 
        f.moment         , 
        f.central_moment , 
        #
        f.mean           , 
        f.rms            , 
        f.variance       , 
        f.dispersion     , 
        f.skewness       , 
        f.kurtosis       , 
        #
        f.quantile       , 
        f.interval       , 
        f.median         , 
        f.quantiles      , 
        f.terciles       , 
        f.quartiles      , 
        f.quintiles      , 
        f.deciles        ,
        # 
        ]
    

if Frames_OK :

    _new_methods_ += [
        frame_the_statVar         , 
        frame_the_statCov         , 
        frame_the_nEff            , 
        frame_the_ECDF            , 
        frame_the_moment          , 
        frame_the_mean            , 
        frame_the_rms             ,
        frame_the_variance        ,
        frame_the_dispersion      ,
        frame_the_skewness        ,
        frame_the_kurtosis        ,
        frame_the_arithmetic_mean ,
        frame_the_geometric_mean  ,
        frame_the_harmonic_mean   ,
        frame_the_power_mean      , 
        frame_the_lehmer_mean     , 
        frame_the_param           , 
        ]
    
    for f in frame_types :

        f.the_statVar          = frame_the_statVar
        f.the_statCov          = frame_the_statCov
        f.the_nEff             = frame_the_nEff 
        f.the_moment           = frame_the_moment 
        f.the_mean             = frame_the_mean 
        f.the_rms              = frame_the_rms        
        f.the_variance         = frame_the_variance 
        f.the_dispersion       = frame_the_dispersion 
        f.the_skewness         = frame_the_skewness 
        f.the_kurtosis         = frame_the_kurtosis
        f.the_ECDF             = frame_the_ECDF
        
        f.the_arithmetic_mean  = frame_the_arithmetic_mean 
        f.the_geometric_mean   = frame_the_geometric_mean 
        f.the_harmonic_mean    = frame_the_harmonic_mean 
        f.the_power_mean       = frame_the_power_mean 
        f.the_lehmer_mean      = frame_the_lehmer_mean 
        f.the_param            = frame_the_param

        f.arithmetic_mean      = frame_arithmetic_mean 
        f.geometric_mean       = frame_geometric_mean 
        f.harmonic_mean        = frame_harmonic_mean 
        f.power_mean           = frame_power_mean 
        f.lehmer_mean          = frame_lehmer_mean 
        f.param                = frame_param
        
        f.ECDF                 = frame_ECDF
        
        _new_methods_ += [ 
            f.arithmetic_mean      , 
            f.geometric_mean       , 
            f.harmonic_mean        , 
            f.power_mean           , 
            f.lehmer_mean          , 
            f.param                , 
            f.the_statVar          , 
            f.the_statCov          , 
            f.the_nEff             , 
            f.the_moment           , 
            f.the_mean             , 
            f.the_rms              , 
            f.the_variance         , 
            f.the_dispersion       , 
            f.the_skewness         , 
            f.the_kurtosis         , 
            
            f.the_arithmetic_mean  , 
            f.the_geometric_mean   , 
            f.the_harmonic_mean    , 
            f.the_power_mean       , 
            f.the_lehmer_mean      , 
            f.the_param            , 
            
            f.arithmetic_mean      , 
            f.geometric_mean       , 
            f.harmonic_mean        , 
            f.power_mean           , 
            f.lehmer_mean          , 
            f.param                ,

            f.ECDF                 , 
            ]
        
    ROOT.TTree. fstatVar = frame_statVar
    ROOT.TTree. fstatCov = frame_statCov
    ROOT.TTree. fmoment  = frame_moment 
    
    _new_methods_ .append ( ROOT.TTree. fstatVar ) 
    _new_methods_ .append ( ROOT.TTree. fstatCov ) 
    _new_methods_ .append ( ROOT.TTree. fmoment  ) 

_decorated_classes_ = frame_types  

_new_methods_ += [ 
    ROOT.TH1.model        , 
    ROOT.TH2.model        , 
    ROOT.TH3.model        , 
    ROOT.TProfile  .model , 
    ROOT.TProfile2D.model ,
    ]

# =============================================================================
## Project the tree to the histogram using DataFrame machinery
#  @code
#  tree  = ...
#  histo = ...
#  tree.fproject ( histo , what , ... ) 
#  @endcode
def _rt_fproject_ ( tree , histo , *args ) :
    """ Project the tree to the histogram using DataFrame machinery
    >>> tree  = ...
    >>> histo = ...
    >>> tree.fproject ( histo , what , ... ) 
    """
    assert isinstance ( histo , ROOT.TH1  ) or \
        isinstance ( histo , _types_nD ) ,  \
        '"histo" must be ROOT.TH1 or polynomial type!'
    ## use frame methods
    frame_project ( tree , histo , *args )
    ## return        
    return histo

_rt_fproject_.__doc__ += '\n' + frame_project.__doc__ 
ROOT.TTree.fproject  = _rt_fproject_

_decorated_classes_ += ( ROOT.TTree , )
_new_methods_.append   ( ROOT.TTree.fproject ) 

# =============================================================================
if Frames_OK : 
    # =========================================================================
    frame_param = _fr_param_
    __all__ = __all__ + ( 'frame_param' , ) 
    
    for f in frame_types : f.param = frame_param

    ## Project/parameterise the tree to the polynomial using DataFrame machinery
    #  @code
    #  tree  = ...
    #  poly  = ...
    #  tree.fparam ( histo , what , ... ) 
    #  @endcode
    def _rt_fparam_ ( tree , poly , *args ) :
        """ Project/parameterise the tree into polynomial using DataFrame machinery
        >>> tree  = ...
        >>> poly  = ...
        >>> tree.fparam ( poly , what , ... ) 
        """
        assert isinstance ( poly , _types_nD  ) , '"poly" must be polynomial!'
        ## use frame methods
        result = frame_param ( tree , poly, *args )
        ## return 
        return result

    _rt_fparam_ .__doc__ += '\n' + frame_param  .__doc__ 

        
    ROOT.TTree.fparam    = _rt_fparam_
    
    _decorated_classes_ += ( ROOT.TTree , )
    _new_methods_.append   ( ROOT.TTree.fparam   ) 

_new_methods_       = tuple ( _new_methods_ ) 

# =============================================================================
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
