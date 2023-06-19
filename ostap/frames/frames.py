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
    'DataFrame'            , ## RDataFrame object
    'report_print'         , ## print the report 
    'report_print_table'   , ## print the report 
    'report_as_table'      , ## print the report
    'frame_progress'       , ## progress bar for frame
    'frame_prescale'       , ## prescale frame (trivial filter)
    ##
    'frame_project'        , ## project data frame to the (1D/2D/3D) histogram
    'frame_print'          , ## print frame
    'frame_draw'           , ## draw variable from the frame  
    'frame_columns'        , ## defined columns  
    ##
    'frame_nEff'           , ## nEff function for frames
    'frame_statVar'        , ## stat var for frame 
    'frame_statCov'        , ## sta tcov for frame
    ## 
    'frame_mean'           , ## mean value for the variable 
    'frame_rms'            , ## RMS for the variable 
    'frame_variance'       , ## variance for the variable 
    'frame_dispersion'     , ## dispersion for the variable
    'frame_skewness'       , ## skewness for the variable 
    'frame_kurtosis'       , ## kurtosos for the variable 
    'frame_get_moment'     , ## get the moment 
    'frame_moment'         , ## moment relative to zero 
    'frame_central_moment' , ## central moment 
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
    ) 
# =============================================================================
from   ostap.core.core           import ( cpp     , Ostap        ,
                                          strings , split_string , var_separators )
from   ostap.core.meta_info      import root_info 
from   ostap.core.ostap_types    import integer_types, string_types, sequence_types   
from   ostap.logger.utils        import multicolumn
from   ostap.utils.progress_conf import progress_conf 
from   ostap.utils.basic         import isatty
import ostap.core.config         as     OCC 
import ostap.stats.statvars      as     SV
import ostap.logger.table        as     T 
import ostap.histos.histos
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

## try : 
##     DataFrame = ROOT.ROOT.RDataFrame
## except AttributeError :
##     DataFrame = ROOT.ROOT.Experimental.TDataFrame
## # =============================================================================
## try : 
##     FrameNode = ROOT.ROOT.RDF.RNode 
## except AttributeError :
##     FrameNode = DataFrame 
    
# =============================================================================
if not hasattr ( Ostap ,'DataFrame' ) :  Ostap.DataFrame = DataFrame
if not hasattr ( Ostap ,'FrameNode' ) :  Ostap.FrameNode = FrameNode 

CNT                = DataFrame.ColumnNames_t 

# =============================================================================
if   ( 6 , 18 ) <= root_info :
    # =========================================================================
    ## get frame-like stuff  as rnode 
    def as_rnode ( frame ) :
        """get frame-like stuff  as rnode"""
        return ROOT.RDF.AsRNode ( frame ) 
else :
    # =========================================================================
    ## get frame-like stuff  as rnode  
    def as_rnode ( frame ) :
        """get frame-like stuff  as rnode"""
        return FrameNode  ( frame ) 
# =============================================================================

# =============================================================================
## get all column names
#  @code
#  cols = colums ( frame ) 
#  @endcode 
def columns ( frame ) :
    """Get all column names
    >>> cols = colums ( frame ) 
    """
    names  = [ str(c) for c in frame.GetColumnNames()        ]
    if ( 6 , 16 ) <= root_info : 
        names += [ str(c) for c in frame.GetDefinedColumnNames() ]            
    return tuple ( sorted ( set ( names ) ) )    
    
frame_columns      = columns

# ==============================================================================
## generate new, unique name for the variable 
def var_name ( prefix , used_names , *garbage ) :
    """Generate new, unique name for the variable 
    """
    name =     prefix + '%x' % ( hash ( ( prefix , used_names , garbage ) )               % ( 2**32 ) ) 
    while name in used_names :
        name = prefix + '%x' % ( hash ( ( name , prefix , used_names , name , garbage ) ) % ( 2**32 ) )
    return name

# ==============================================================================
## Is implicit MC globally enabled? 
mt_global = OCC.general.getboolean ( 'ImplicitMT' , fallback = False )
mt_global = True 
    
# ==============================================================================
## modify constructor for RDataFrame to enable/disable Implicit multithreading
#  @code
#  f1 = DataFrame ( ... , enable = True  ) ## default
#  f2 = DataFrame ( ... , enable = False ) ## default
#  @endcode
#  @see ROOT::EnableImplicitMT 
#  @see ROOT::DisableImplicitMT
#  @see ROOT::IsImplicitMTEnabled
def _fr_new_init_ ( self , name , *args , **kwargs ) :
    """Modify the DataFrame constuctor to allow (semi)automatic
    manipulations wth ROOT.ROOT.EnableImplicitMT/DisableImplicitMT
    - see ROOT.ROOT.EnableImplicitMT 
    - see ROOT.ROOT.DisableImplicitMT
    - see ROOT.ROOT.IsImplicitMTEnabled
    >>> f = DataFrame ( .... , enable = True  ) ## default 
    >>> f = DataFrame ( .... , enable = False )
    """
        
    mt = kwargs.pop ( 'enable' , mt_global ) and mt_global 
    
    if       mt and not ROOT.ROOT.IsImplicitMTEnabled() :
        ROOT.ROOT.EnableImplicitMT  ()
        logger.debug ( 'DataFrame:ImplicitMT is %s' % ( 'Enabled' if ROOT.ROOT.IsImplicitMTEnabled() else 'Disabled' ) )
    elif not mt and     ROOT.ROOT.IsImplicitMTEnabled() :
        ROOT.ROOT.DisableImplicitMT ()
        logger.info  ( 'DataFrame:ImplicitMT is %s' % ( 'Enabled' if ROOT.ROOT.IsImplicitMTEnabled() else 'Disabled' ) ) 
        
    self._fr_old_init_ ( name , *args , **kwargs ) 

if not hasattr ( DataFrame , '_fr_old_init_' ) :
    DataFrame._fr_old_init_ = DataFrame.__init__
    DataFrame.__init__      = _fr_new_init_
    
# =============================================================================
## Get the length/size of the data frame
#  @code
#  frame = ...
#  print len(frame)
#  @endcode 
def _fr_len_ ( frame ) :
    """Get the length/size of the data frame
    >>> frame = ...
    >>> print len(frame)
    """
    node = as_rnode ( frame )
    return node.Count().GetValue() 

# =============================================================================
## Draw (lazy) progress bar for the    DataFrame:
#  @code
#  f = DataFrame ( ... )
#  p = frame_progress ( f , 1000 ) ## number of elements!
#  p = f.ProgressBar  ( 1000 ) ## number of elements!

#  ....
#  @endcode 
def frame_progress ( frame  ,
                     length ) :
    """ Draw (lazy) progress bar for the    DataFrame:
    >>> f = DataFrame ( ... )
    >>> p = f.ProgressBar  ( 1000 ) ## number of elements!
    >>> ....
    """
    
    cnt = frame.Count () 
    if not isatty () : return cnt

    if not length : length = len ( frame )
    
    if   2048 < length  : nchunks = 2000 
    elif 1024 < length  : nchunks = 1000 
    elif  512 < length  : nchunks =  500 
    elif  101 < length  : nchunks =  100 
    else                : return cnt     ## no progress bar for short frames 
    
    csize , rr = divmod ( length , nchunks )
    csize      = max    ( csize  , 1       )
    
    ## if rr : nchunks += 1 
    
    fun = Ostap.Utils.frame_progress ( nchunks , progress_conf () )
    cnt.OnPartialResultSlot  ( csize , fun )
    
    return cnt

# =============================================================================
## Get the effective entries in data frame
#  @code
#  data = ...
#  neff = data.nEff('b1*b1')
#  @endcode
def frame_nEff ( frame , cuts = '' ) :
    """Get the effective entries in data frame 
    >>> data = ...
    >>> neff = data.nEff('b1*b1')
    """
    node  = as_rnode ( frame )    
    return Ostap.StatVar.nEff ( node , cuts )

# =============================================================================
## Get statistics for the  given expression in data frame
#  @code
#  data = ...
#  c1 = data.statVar( 'S_sw' , 'pt>10' ) 
#  c2 = data.statVar( 'S_sw' , 'pt>0'  )
#  @endcode
def frame_statVar ( frame , expression ,  cuts = '' ) :
    """Get statistics for the  given expression in data frame
    >>> data = ...
    >>> c1 = data.statVar( 'S_sw' , 'pt>10' ) 
    >>> c2 = data.statVar( 'S_sw' )
    """
    if isinstance ( frame  , ROOT.TTree ) : frame = DataFrame ( frame  )
    node = as_rnode ( frame ) 
    return Ostap.StatVar.statVar ( node , str ( expression ) , str ( cuts )  )

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
    """Get the statistic for pair of expressions in DataFrame
    
    >>>  frame  = ...
    >>>  stat1 , stat2 , cov2 , len = framw.statCov( 'x' , 'y' )
    
    Apply some cuts:
    >>> stat1 , stat2 , cov2 , len = frame.statCov( 'x' , 'y' , 'z>0' )
    
    """
    import ostap.math.linalg 
    stat1  = Ostap.WStatEntity       ()
    stat2  = Ostap.WStatEntity       ()
    cov2   = Ostap.Math.SymMatrix(2) ()
    
    if isinstance ( frame  , ROOT.TTree ) : frame = DataFrame ( frame  )
    node   = as_rnode ( frame ) 
    length = Ostap.StatVar.statCov ( node        ,
                                     expression1 ,
                                     expression2 ,
                                     cuts        ,    
                                     stat1       ,
                                     stat2       ,
                                     cov2        ) 
        
    return stat1 , stat2 , cov2 , length

# ==================================================================================
## get statistics of variable(s)
#  @code
#  frame = ....
#  stat  = frame.statVar ( 'pt'           , lazy = True )
#  stat  = frame.statVar ( 'pt' , 'eta>0' , lazy = True )
#  @endcode
def _fr_statVar_new_ ( frame , expressions , cuts = '' , lazy = False  ) :
    """Get statistics of variable(s)
    >>> frame = ....
    >>> stat  = frame.statVar ( 'pt'           , lazy = True )
    >>> stat  = frame.statVar ( 'pt' , 'eta>0' , lazy = True )
    """
    
    if isinstance ( frame , ROOT.TTree ) : frame = DataFrame ( frame )
        
    node = as_rnode ( frame ) 

    input_string = False 
    if isinstance ( expressions , string_types ) :
        input_string = True 
        expressions = [ expressions ]
    
    ## get the list of currently known names
    vars     = frame_columns ( node ) 
    all_vars = set ( vars ) 
    
    names    = {}
    current  = node  
    for e in expressions :

        e = str ( e )
        if e in vars :
            names [ e ] = e 
            continue
        
        used    = tuple ( all_vars | set ( frame_columns ( current ) ) ) 
        vn      = var_name ( 'var_' , used , e , *vars )
        all_vars.add ( vn )
        current = current.Define ( vn , e )
        names    [ e  ] = vn

    cuts  = str ( cuts )
    cname = cuts 
    if cuts and not cuts in vars :
        used    = tuple ( all_vars | set ( frame_columns ( current ) ) ) 
        vn      = var_name ( 'cut_' , used , cuts , *vars )
        all_vars.add ( vn )
        current = current.Define ( vn , cuts )
        cname   = vn

    results = {}
    for e in names :

        if cname : 
            results [ e ] = current.Book( ROOT.std.move ( Ostap.Actions.WStatVar() ) , CNT ( [ names [ e ] , cname ] ) ) 
        else :
            results [ e ] = current.Book( ROOT.std.move ( Ostap.Actions. StatVar() ) , CNT ( 1 , names[e] ) ) 

    if not lazy :
        for e in results :
            r = results [ e ]
            results [ e ] = r.GetValue()  

    if input_string and 1 == len ( results ) :
        e , r = results.popitem()
        return r
    
    return results 


# ==================================================================================
## get statistics of variable(s)
#  @code
#  frame = ....
#  stat  = frame.the_moment ( 5 , 'pt'           )
#  stat  = frame.the_moment ( 5 , 'pt' , 'eta>0' )
#  @endcode
def _fr_the_moment_ ( frame , N , expression , cuts = '' ) :
    """Get statistics of variable(s)
    >>> frame = ....
    >>> stat  = frame.the_moment ( 5 , 'pt'           )
    >>> stat  = frame.the_moment ( 5 , 'pt' , 'eta>0' )
    """

    assert isinstance ( N , int ) and 0 <= N , 'Invalid order!'
    assert isinstance ( expression , str )   , 'Invalid expression type!'
    
    if isinstance ( frame , ROOT.TTree ) : frame = DataFrame ( frame )
        
    node = as_rnode ( frame ) 

    ## get the list of currently known names
    vars     = frame_columns ( node ) 
    all_vars = set ( vars ) 

    current  = node
    
    vname = expression 
    if not expression in all_vars :
        
        used    = tuple ( all_vars | set ( frame_columns ( current ) ) ) 
        vn      = var_name ( 'var_' , used , expression , *vars )
        all_vars.add ( vn )
        current = current.Define ( vn , expression )
        vname = expression
        
    cuts  = str ( cuts )
    cname = cuts 
    if cuts and not cuts in vars :
        used    = tuple ( all_vars | set ( frame_columns ( current ) ) ) 
        vn      = var_name ( 'cut_' , used , cuts , *vars )
        all_vars.add ( vn )
        current = current.Define ( vn , cuts )
        cname   = vn

    if cname :
        TT = ROOT.Detail.RDF.WMoment_(N) 
        return current.Book( ROOT.std.move ( TT () ) , CNT ( [ vname , cname ] ) ) 
    else :
        TT = ROOT.Detail.RDF.Moment_(N)             
        return current.Book( ROOT.std.move ( TT () ) , CNT ( 1 , vname ) ) 

    return results 

# ============================================================================
_types_1D = Ostap.Math.LegendreSum  , Ostap.Math.Bernstein   , Ostap.Math.ChebyshevSum , 
_types_2D = Ostap.Math.LegendreSum2 , Ostap.Math.Bernstein2D ,
_types_3D = Ostap.Math.LegendreSum3 , Ostap.Math.Bernstein3D , 
_types_4D = Ostap.Math.LegendreSum4 ,
_types_nD = _types_1D + _types_2D + _types_3D + _types_4D 
# =============================================================================
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
#  res   = frame_param ( frame , bs2 , 'y' , 'x' , 'z>0' )
#
#  ls3   = Ostap.Math.LegendreSum3 ( ... )
#  res   = frame_param ( frame , ls3 , 'z' , 'y' , 'x' , 'z>0' )
# 
#  bs2   = Ostap.Math.Bernstein2D  ( ... )
#  res   = frame_param ( frame , bs2 , 'y' , 'x' , 'z>0' )
# 
#  @endcode
def _fr_param_ ( frame , poly , *expressions ) :
    """ `project/parameterise` frame into polynomial structures
    >>> frame = ...
    
    >>> ls    = Ostap.Math.LegendreSum ( ... )
    >>> res   = frame_param ( frame , ls , 'x' , 'y>0' )
    
    >>> bs    = Ostap.Math.Bernstein   ( ... )
    >>> res   = frame_param ( frame , bs , 'x' , 'y>0' )

    >>> cs    = Ostap.Math.ChebyshevSum ( ... )
    >>> res   = frame_param ( frame , cs , 'x' , 'y>0' )

    >>> ls2   = Ostap.Math.LegendreSum2 ( ... )
    >>> res   = frame_param ( frame , ls2 , 'y' , 'x' , 'z>0' )

    >>> bs2   = Ostap.Math.Bernstein2D  ( ... )
    >>> res   = frame_param ( frame , bs2 , 'y' , 'x' , 'z>0' )

    >>> ls3   = Ostap.Math.LegendreSum3 ( ... )
    >>> res   = frame_param ( frame , ls3 , 'z' , 'y' , 'x' , 'z>0' )

    >>> bs2   = Ostap.Math.Bernstein2D  ( ... )
    >>> res   = frame_param ( frame , bs2 , 'y' , 'x' , 'z>0' )
    """
    
    if isinstance ( frame , ROOT.TTree ) : frame = DataFrame ( frame )
    
    node = as_rnode ( frame )
    
    items = []
    for e in expressions :
        items += split_string ( e , ',:;' , strip = True , respect_groups = True )
        
    ## get the list of currently known names
    vars     = frame_columns ( node ) 
    all_vars = set ( vars ) 


    assert ( isinstance ( poly , _types_1D ) and 1 <= len ( items ) <= 2 ) or \
           ( isinstance ( poly , _types_2D ) and 2 <= len ( items ) <= 3 ) or \
           ( isinstance ( poly , _types_3D ) and 3 <= len ( items ) <= 4 ) or \
           ( isinstance ( poly , _types_4D ) and 4 <= len ( items ) <= 5 ) , \
           "Invalid structure of  polynomial and variables/cuts!"
    
    current = node
    
    ## actual variables 
    uvars   = []   
    for e in items :        
        if e in all_vars :
            uvars.append ( e )
            continue
        else :
            used    = tuple ( all_vars | set ( frame_columns ( current ) ) ) 
            vn      = var_name ( 'var_' , used , e , *vars )
            all_vars.add ( vn )
            current = current.Define ( vn , e )        
            uvars.append ( vn )

    if ( isinstance ( poly , _types_1D ) and 2 == len ( uvars ) ) or \
       ( isinstance ( poly , _types_2D ) and 3 == len ( uvars ) ) or \
       ( isinstance ( poly , _types_3D ) and 4 == len ( uvars ) ) or \
       ( isinstance ( poly , _types_4D ) and 5 == len ( uvars ) ) :
        
        cuts    = uvars [ -1]        
        current = current.Filter ( cuts )        
        
        wn      = var_name ( 'weight_' , used , cuts , *vars )
        all_vars.add ( wn )
        current = current.Define ( wn , '1.0*(%s)' % cuts )
        
        uvars   = [ u for u in reversed ( uvars [:-1] ) ] + [ wn ]    ## ATTENTION !! RESERSED HERE!
        
    else :
        uvars   = [ u for u in reversed ( uvars ) ]             ## ATTENTION !! RESERSED HERE!

    
    ## finally book the actions!
    if   isinstance ( poly , Ostap.Math.LegendreSum  ) :
        result = current.Book ( ROOT.std.move ( Ostap.Actions.LegendrePoly   ( poly ) ) , CNT ( uvars ) )
    elif isinstance ( poly , Ostap.Math.LegendreSum2 ) :
        result = current.Book ( ROOT.std.move ( Ostap.Actions.LegendrePoly2  ( poly ) ) , CNT ( uvars ) )
    elif isinstance ( poly , Ostap.Math.LegendreSum3 ) :
        result = current.Book ( ROOT.std.move ( Ostap.Actions.LegendrePoly3  ( poly ) ) , CNT ( uvars ) )
    elif isinstance ( poly , Ostap.Math.LegendreSum4 ) :
        result = current.Book ( ROOT.std.move ( Ostap.Actions.LegendrePoly4  ( poly ) ) , CNT ( uvars ) )
    elif isinstance ( poly , Ostap.Math.ChebyshevSum ) :
        result = current.Book ( ROOT.std.move ( Ostap.Actions.ChebyshevPoly  ( poly ) ) , CNT ( uvars ) )
    elif isinstance ( poly , Ostap.Math.Bernstein    ) :
        result = current.Book ( ROOT.std.move ( Ostap.Actions.BernsteinPoly  ( poly ) ) , CNT ( uvars ) )
    elif isinstance ( poly , Ostap.Math.Bernstein2D  ) :
        result = current.Book ( ROOT.std.move ( Ostap.Actions.BernsteinPoly2 ( poly ) ) , CNT ( uvars ) )
    elif isinstance ( poly , Ostap.Math.Bernstein3D  ) :
        result = current.Book ( ROOT.std.move ( Ostap.Actions.BernsteinPoly3 ( poly ) ) , CNT ( uvars ) )
        
    return result


# =============================================================================
## Simplified print out for the  frame 
#  @code 
#  frame = ...
#  print frame
#  @endcode 
def frame_print ( frame ) :
    """Simplified print out for the  frame 
    
    >>> frame = ...
    >>> print frame
    """
    ## 
    if isinstance ( frame  , ROOT.TTree ) : frame = DataFrame ( frame )
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
def _fr_table_ ( frame , pattern = None ,  cuts = '' , more_vars = () , prefix = '' ) :
    """Data frame as table
    >>> frame = ...
    >>> table = frame_table ( frame , '.*PT.*' , cuts = ... , more_vars = [ 'x*x/y' , 'y+z'] )
    >>> print ( table )
    """
    
    if isinstance ( frame  , ROOT.TTree ) : frame = DataFrame ( frame  )
    frame = as_rnode ( frame )
    
    def col_type ( var ) :
        
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
        elif 'UShort_t'  == t : return 'unisgned short'
        elif 'UInt_t'    == t : return 'unsigned int'
        elif 'ULong64_t' == t : return 'unisgned long'
        
        return t 
        
    ## get all variables/columns 
    cols  = tuple ( frame.GetColumnNames() )
    
    vcols = cols
    
    if pattern :
        import re
        cpat  = re.compile ( pattern )
        vcols = tuple ( [ v for v in vcols if cpat.match ( v ) ] ) 

    svars = vcols + tuple  ( more_vars ) 
    stats = frame_statVars ( frame , svars , cuts = cuts , lazy = False ) 
    
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
    title  = 'DataFrame variables'
    if pattern : title = '%s (pattern %s)' % ( title , pattern )
    t  = T.table (  table , title , prefix = prefix , alignment = 'rllrlrl' )
    
    return t 

# ===============================================================================
## Print the frame report data
def report_print_table ( report , title  = '' , prefix = '' , more_rows = [] ) :
    """Print a frame report data 
    """
    from ostap.core.core import binomEff
    
    n0     = -1 
    lmax   =  5
    table  = []
    
    for name, passed, all in report :

        n0 = max ( n0 , all , passed )
        
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
    """Print a frame report
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
    """Print the data frame report
    """
    table = report_as_table ( report )    
    return report_print_table ( table , title, prefix , more_rows )


if   ( 6 , 14 ) <= root_info :
    
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

elif ( 6, 12 ) <= root_info :

    DF_H1Model = ROOT.ROOT.Experimental.TDF.TH1DModel 
    DF_H1Type  = ROOT.TH1D
    DF_H2Model = ROOT.ROOT.Experimental.TDF.TH2DModel 
    DF_H2Type  = ROOT.TH2D
    DF_H3Model = ROOT.ROOT.Experimental.TDF.TH3DModel 
    DF_H3Type  = ROOT.TH3D

    DF_P1Model = ROOT.ROOT.Experimental.TDF.TProfile1DModel 
    DF_P1Type  = ROOT.TProfile
    
    DF_P2Model = ROOT.ROOT.Experimental.TDF.TProfile2DModel
    DF_P2Type  = ROOT.TProfile2D

else :

    DF_H1Model = ROOT.TH1F
    DF_H1Type  = ROOT.TH1F
    DF_H2Model = ROOT.TH2F
    DF_H2Type  = ROOT.TH2F
    DF_H3Model = ROOT.TH3F 
    DF_H3Type  = ROOT.TH3F
    DF_P1Model = ROOT.TProfile
    DF_P1Type  = ROOT.TProfile    
    DF_P2Model = ROOT.TProfile
    DF_P2Type  = ROOT.TProfile2D

# =============================================================================
## convert 1D-histogram to "model" for usage with DataFrames 
def _h1_model_ ( histo ) :
    """Convert 1D-histogram to 'model' for usage with DataFrames """
    model = histo 
    if not isinstance ( model , DF_H1Type ) :
        model = DF_H1Type()
        histo.Copy ( model )
    return DF_H1Model ( model ) 
# =============================================================================
## convert 2D-histogram to "model" for usage with DataFrames 
def _h2_model_ ( histo ) :
    """Convert 2D-histogram to 'model' for usage with DataFrames """
    model = histo 
    if not isinstance ( model , DF_H2Type ) :
        model = DF_H2Type()
        histo.Copy ( model )
    return DF_H2Model ( model ) 
# =============================================================================
## convert 3D-histogram to "model" for usage with DataFrames 
def _h3_model_ ( histo ) :
    """Convert 3D-histogram to 'model' for usage with DataFrames """
    model = histo 
    if not isinstance ( model , DF_H3Type ) :
        model = DF_H3Type()
        histo.Copy ( model )
    return DF_H3Model ( model ) 
# =============================================================================
## convert 1D-profile to "model" for usage with DataFrames 
def _p1_model_ ( histo ) :
    """Convert 1D-profile to 'model' for usage with DataFrames """
    model = histo 
    if not isinstance ( model , DF_P1Type ) :
        model = DF_P1Type()
        histo.Copy ( model )
    return DF_P1Model ( model ) 
# =============================================================================
## convert 2D-profile to "model" for usage with DataFrames 
def _p2_model_ ( histo ) :
    """Convert 2D-profile to 'model' for usage with DataFrames """
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
## 'Lazy' project of the frame
#  @code
#  frame    = ...
#  h1_model = ...
#  h1       = frame_project ( frame , h1_model , 'pt' )
#  h2_model = ...
#  h2       = frame_project ( frame , h2_model , 'pt' , 'x' )
#  h3_model = ...
#  h3       = frame_project ( frame , h3_model , 'pt' , 'x'  , 'y' )
#  ...
#  @endcode
#  Cuts/weigth are also can be specified
#  @code
#  frame    = ...
#  h1_model = ...
#  h1       = frame_project ( frame , h1_model , 'pt' , 'w')
#  h2_model = ...
#  h2       = frame_project ( frame , h2_model , 'pt' , 'x' , 'w')
#  h3_model = ...
#  h3       = frame_project ( frame , h3_model , 'pt' , 'x'  , 'y' , 'w' )
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
# 
def frame_project ( frame , model , *what ) :
    """ 'Lazy' project of the frame
    >>> frame    = ...
    >>> h1_model = ...
    >>> h1       = frame_project ( frame , h1_model , 'pt' )
    >>> h2_model = ...
    >>> h2       = frame_project ( frame , h2_model , 'pt' , 'x' )
    >>> h3_model = ...
    >>> h3       = frame_project ( frame , h3_model , 'pt' , 'x'  , 'y' )
    >>> ...
    
    - Cuts/weigth are also can be specified
    
    >>> frame    = ...
    >>> h1_model = ...
    >>> h1       = frame_project ( frame , h1_model , 'pt' , 'w')
    >>> h2_model = ...
    >>> h2       = frame_project ( frame , h2_model , 'pt' , 'x' , 'w')
    >>> h3_model = ...
    >>> h3       = frame_project ( frame , h3_model , 'pt' , 'x'  , 'y' , 'w' )
    >>>...
    
    - If model is a real 1D/2D/3D-hstogram, action is *not* lazy, otherwise action *is* lazy

    Model can be a real histogram or 
    - `ROOT.ROOT.RDF.TH1DModel`
    - `ROOT.ROOT.RDF.TH2DModel`
    - `ROOT.ROOT.RDF.TH3DModel`
    - anything that can be converted to `ROOT.ROOT.RDF.TH1DModel`, 
    `ROOT.ROOT.RDF.TH2DModel` or `ROOT.ROOT.RDF.TH3DModel` objects 
    """

    if isinstance ( frame , ROOT.TTree ) : frame = DataFrame ( frame )
    
    frame = as_rnode  ( frame )
    
    if ( 6 , 16 ) <= root_info and isinstance ( model , _types_nD ) :
        return _fr_param_ ( frame , model , *what ) 
       
    if 1 <= len ( what ) <= 2 :
        if   isinstance  ( what [ 0 ] , string_types ) :
            ww = split_string ( what [ 0 ] , var_separators , strip = True, respect_groups = True  )
            if 1 < len ( ww ) : ww.reverse()                ## ATTENTION HERE: REVERSE!! 
            what = tuple ( ww ) + what [1:]                                              
        elif isinstance ( what [ 0 ] , sequence_types ) :
            what = tuple ( w for w in what [ 0 ] ) + what [1:] 

    ## strip blanks
    what = tuple ( w.strip() for w in what )
    
    ## cuts are empty 
    if 1 < len ( what ) and not what [ -1 ] : what = what [ :-1 ]

    histo = None

    #
    ## convert histogram-like objects into 'models'
    #
    
    if   isinstance ( model , ROOT.TProfile2D ) :        
        histo = model
        model = model.model ()        
    elif isinstance ( model , ROOT.TProfile   ) :        
        histo = model
        model = model.model ()        
    elif isinstance ( model , ROOT.TH3 ) and 3 == model.dim () :        
        histo = model
        model = model.model ()                
    elif isinstance ( model , ROOT.TH2 ) and 2 == model.dim () :
        histo = model
        model = model.model ()
    elif isinstance ( model , ROOT.TH1 ) and 1 == model.dim () :
        histo = model
        model = model.model ()
            
    if histo : histo.Reset()

    ## get the list of currently known names
    vars = frame_columns ( frame )
    
    pvars    = [] 
    current  = frame 
    all_vars = set ( vars )     
    for ww in what :

        w = ww
        if isinstance ( ww , ROOT.TCut ) : w = str ( ww )

        if not w : continue
        
        if   w in  vars : pvars.append ( w )
        elif w in pvars : pvars.append ( w )
        else :
            used    = tuple ( all_vars | set ( frame_columns ( current ) ) ) 
            vname   = var_name ( 'var_' , used , *what )
            current = current.Define ( vname , w )
            all_vars.add ( vname ) 
            pvars.append ( vname )

            
    numvars = len ( pvars ) 
            
    if   isinstance ( model , DF_P2Model ) and 3 <= numvars <= 4 : action = current.Profile2D ( model , *pvars )
    elif isinstance ( model , DF_P1Model ) and 2 <= numvars <= 3 : action = current.Profiel1D ( model , *pvars )
    elif isinstance ( model , DF_H3Model ) and 3 <= numvars <= 4 : action = current.Histo3D   ( model , *pvars )
    elif isinstance ( model , DF_H2Model ) and 2 <= numvars <= 3 : action = current.Histo2D   ( model , *pvars )
    elif isinstance ( model , DF_H1Model ) and 1 <= numvars <= 2 : action = current.Histo1D   ( model , *pvars )
    else :
        raise TypeError ('Invalid model/what objects %s %s ' % ( type ( model ) , what ) ) 

    ## ATTENTION! lazy action!
    if not histo :
        return action

    ## make a loop! 
    histo += action.GetValue() 
    
    return histo 

# ==============================================================================
## Prescale the frame
#  @code
#  frame = ...
#  prescaled1 = frame_prescale ( frame , 0.15 ) ##  0.0 < prescale < 1.0
#  prescaled2 = frame_prescale ( frame , 20   ) ##  1   < prescale 
#  @endcode 
def frame_prescale ( frame , prescale , name = '' ) :
    """Prescale the frame
    >>> frame = ...
    >>> prescaled1 = frame_prescale ( frame , 0.15 ) ##  0.0 < prescale < 1.0
    >>> prescaled2 = frame_prescale ( frame , 20   ) ##  1   < prescale 
    """
    node = as_rnode ( frame )
    
    if isinstance ( prescale , integer_types ) and 1 < prescale :

        name = name if name else 'PRESCALE_%d' % prescale
        
        code = '( 0 == ( ( rdfentry_ + %d * rdfslot_ ) %% %d ) )' \
               if ( 6 , 16 ) < root_info else              \
               '( 0 == ( ( tdfentry_ + %d * tdfslot_ ) %% %d ) ) '
        
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
## Get the moment  for the frame object
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
    node = as_rnode ( frame )
    return SV.data_central_moment ( node , order = order , expression = expression , cuts = cuts )

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
    """ Get the skewness for the frame object
    >>> frame = ...
    >>> value = frame_kurtosis ( 'b1*b2' , 'b3>0' )
    - see `Ostap.StatVar.skewness`
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
    >>> value = frame_mean ( 'b1*b2' , 'b3>0' , exact = True  )
    >>> value = frame_mean ( 'b1*b2' , 'b3>0' , exact = False ) ## use P2 algorithm 
    - see `Ostap.StatVar.moment`
    """    
    return frame_central_moment ( frame , order = 2 , expression = expression , cuts = cuts )

from ostap.fitting.dataset import ds_draw as _ds_draw
# ==============================================================================
## draw for frame 
def frame_draw ( frame , *args , **kwargs ) :
    """Draw function for tehe frames
    """
    node = as_rnode ( frame )
    return _ds_draw ( node , *args , **kwargs ) 

# ==============================================================================
# decorate 
# ==============================================================================
_new_methods_ = [
    frame_nEff           , 
    frame_statVar        , 
    frame_statCov        , 
    frame_progress       , 
    frame_draw           , 
    frame_project        , 
    frame_columns        , 
    ##
    frame_mean           , 
    frame_rms            , 
    frame_variance       , 
    frame_dispersion     , 
    frame_skewness       , 
    frame_kurtosis       , 
    frame_get_moment     , 
    frame_moment         , 
    frame_central_moment , 
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
    ] 

# ==============================================================================
if FrameNode is DataFrame : frames = DataFrame ,
else                      : frames = DataFrame , FrameNode

# ===============================================================================
has_std_move = False 
try :
    mv           = ROOT.std.move
    has_std_move = True and (6,25)<= root_info  ## ATTENTION! 
except AttributeError : 
    has_std_move = False

# ================================================================================

## if ( 6 , 25 ) <= root_info :
## if ( 6 , 16 ) <= root_info and has_std_move :
if ( 6 , 25 ) <= root_info and has_std_move :
    frame_statVar       = _fr_statVar_new_
    frame_statVars      = _fr_statVar_new_
    frame_table         = _fr_table_
    frame_the_moment    = _fr_the_moment_ 
    for f in frames : 
        f.statVars   = frame_statVars  
        _new_methods_ .append ( f.statVars    )
        f.table     = frame_table   
        _new_methods_ .append ( f.table    )
        f.the_moment = frame_the_moment  
        _new_methods_ .append ( f.the_moment  )

    _new_methods_ .append ( frame_statVars   ) 
    _new_methods_ .append ( frame_table      ) 
    _new_methods_ .append ( frame_the_moment ) 
    __all__ = __all__ + ( 'frame_statVars' ,
                          'frame_table'    ,
                          'frame_the_moment' ) 

    ROOT.TTree. fstatVars   = frame_statVars
    ROOT.TTree. fthe_moment = frame_statVars
    
    _new_methods_ .append ( ROOT.TTree. fstatVars   ) 
    _new_methods_ .append ( ROOT.TTree. fthe_moment ) 


for f in frames :
    
    if not hasattr ( f , '__len__') :
        f.__len__ = _fr_len_
        _new_methods_ .append ( f.__len__   )
    
    f.__str__        = frame_print
    f.__repr__       = frame_print

    f.nEff           = frame_nEff
    f.statVar        = frame_statVar
    f.statCov        = frame_statCov
    f.ProgressBar    = frame_progress
    f.progress       = frame_progress
    f.prescale       = frame_prescale 
    
    f.columns        = frame_columns
    f.branches       = frame_columns

    f.mean           = frame_mean    
    f.rms            = frame_rms     
    f.variance       = frame_variance
    f.dispersion     = frame_dispersion 
    f.skewness       = frame_skewness 
    f.kurtosis       = frame_kurtosis
    f.get_moment     = frame_get_moment 
    f.moment         = frame_moment 
    f.central_moment = frame_central_moment 

    f.quantile       = frame_quantile       
    f.interval       = frame_interval       
    f.median         = frame_median         
    f.quantiles      = frame_quantiles      
    f.terciles       = frame_terciles      
    f.quartiles      = frame_quartiles     
    f.quintiles      = frame_quintiles      
    f.deciles        = frame_deciles
    
    _new_methods_ .append ( f.__str__        )
    _new_methods_ .append ( f.__repr__       )
    _new_methods_ .append ( f. nEff          )
    _new_methods_ .append ( f. statVar       )
    _new_methods_ .append ( f. ProgressBar   )
    _new_methods_ .append ( f. progress      )
    _new_methods_ .append ( f. columns       )
    _new_methods_ .append ( f. branches      )
    
    _new_methods_ .append ( f.mean           )
    _new_methods_ .append ( f.rms            ) 
    _new_methods_ .append ( f.variance       ) 
    _new_methods_ .append ( f.dispersion     ) 
    _new_methods_ .append ( f.skewness       ) 
    _new_methods_ .append ( f.kurtosis       ) 
    _new_methods_ .append ( f.get_moment     ) 
    _new_methods_ .append ( f.moment         ) 
    _new_methods_ .append ( f.central_moment )

    _new_methods_ .append ( f.quantile       )
    _new_methods_ .append ( f.interval       )
    _new_methods_ .append ( f.median         )
    _new_methods_ .append ( f.quantiles      ) 
    _new_methods_ .append ( f.terciles       )
    _new_methods_ .append ( f.quartiles      )
    _new_methods_ .append ( f.quintiles      )
    _new_methods_ .append ( f.deciles        )

    f.draw       = frame_draw
    
    _new_methods_ .append ( f. draw       )
    


_decorated_classes_ = frames 

_new_methods_ += [ 
    ROOT.TH1.model        , 
    ROOT.TH2.model        , 
    ROOT.TH3.model        , 
    ROOT.TProfile  .model , 
    ROOT.TProfile2D.model ,
    ]


# =============================================================================
if ( 6 , 18 ) <= root_info :
    # =========================================================================
    ## Project the tree to the histogram using DataFrame machinery
    #  @code
    #  tree  = ...
    #  histo = ...
    #  tree.fproject ( histo , what , ... ) 
    #  @endcode
    def _rt_fproject_ ( tree , histo , *args ) :
        """Project the tree to the histogram using DataFrame machinery
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
##  if ( 6 , 16 ) <= root_info and has_std_move :
if ( 6 , 25 ) <= root_info and has_std_move :
    # =========================================================================
    frame_param = _fr_param_
    __all__ = __all__ + ( 'frame_param' , ) 
    
    for f in frames : f.param = frame_param

    ## Project/parameterise the tree to the polynomial using DataFrame machinery
    #  @code
    #  tree  = ...
    #  poly  = ...
    #  tree.fparam ( histo , what , ... ) 
    #  @endcode
    def _rt_fparam_ ( tree , poly , *args ) :
        """Project/parameterise the tree into polynomial using DataFrame machinery
        >>> tree  = ...
        >>> poly  = ...
        >>> tree.fparam ( poly , what , ... ) 
        """
        assert isinstance ( poly , _types_nD  ) , '"poly" must be polynomial!'
        ## use frame methods
        result = frame_param ( tree , poly, *args )
        ## return 
        return result.GetValue() 

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
