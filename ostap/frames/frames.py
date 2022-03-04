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
    'DataFrame'          , ## RDataFrame object
    'report_print'       , ## print the report 
    'report_print_table' , ## print the report 
    'report_as_table'    , ## print the report
    'frame_progress'     , ## progress bar for frame 
    'frame_table'        , ## print data frame as table 
    'frame_project'      , ## project data frame to the (1D/2D/3D) histogram
    'frame_print'        , ## print frame
    'frame_nEff'         , ## nEff function for frames
    'frame_statVar'      , ## stat var for frame 
    'frame_statCov'      , ## sta tcov for frame 
    ) 
# =============================================================================
import ROOT
# =============================================================================
from   ostap.core.core        import cpp, Ostap, strings, split_string
from   ostap.core.meta_info   import root_info 
from   ostap.core.ostap_types import integer_types, string_types  
from   ostap.logger.utils     import multicolumn
from   ostap.utils.basic      import terminal_size, isatty
from   ostap.logger.colorized import allright
import ostap.histos.histos
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.frames.frames' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Some useful decorations for ROOT::RDataFrame objects')
# =============================================================================
try : 
    DataFrame = ROOT.ROOT.RDataFrame
except AttributeError :
    DataFrame = ROOT.ROOT.Experimental.TDataFrame 
# =============================================================================
Ostap.DataFrame    = DataFrame 
CNT                = DataFrame.ColumnNames_t 

DataFrame.columns  = lambda s : tuple ( [ str(c) for c in s.GetColumnNames() ] ) 
DataFrame.branches = DataFrame.columns 


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
        
    mt = kwargs.pop ( 'enable' , True )
    
    if       mt and not ROOT.ROOT.IsImplicitMTEnabled() :
        ROOT.ROOT.EnableImplicitMT  ()
        logger.debug ( 'DataFrame:ImplicitMT is %s' % ( 'Enabled' if ROOT.ROOT.IsImplicitMTEnabled() else 'Disabled' ) )
    elif not mt and     ROOT.ROOT.IsImplicitMTEnabled() :
        ROOT.ROOT.DisableImplicitMT ()
        logger.info  ( 'DataFrame: ImplicitMT is %s' % ( 'Enabled' if ROOT.ROOT.IsImplicitMTEnabled() else 'Disabled' ) ) 
        
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
def _fr_len_ ( f ) :
    """Get the length/size of the data frame
    >>> frame = ...
    >>> print len(frame)
    """
    cpf = Ostap.DataFrame ( f )   ## make independent loop?
    return cpf.Count().GetValue() 

# =============================================================================
## Draw (lazy) progress bar for the    DataFrame:
#  @code
#  f = DataFrame ( ... )
#  p = frame_progress ( f , 1000 ) ## number of elements!
#  p = f.ProgressBar  ( 1000 ) ## number of elements!

#  ....
#  @endcode 
def frame_progress ( self          , 
                     length        ,
                     width  = None ,
                     symbol = "#"  ) :
    """ Draw (lazy) progress bar for the    DataFrame:
    >>> f = DataFrame ( ... )
    >>> p = f.ProgressBar  ( 1000 ) ## number of elements!
    >>> ....
    """
    
    cnt = self.Count () 
    if not isatty() : return cnt
    
    length = length if isinstance ( length , integer_types ) and  0 < length else len ( self )
    width  = width  if isinstance ( width  , integer_types ) and 10 < width  else terminal_size()[1]
    
    if width < 16 : width = 16 
    nchunks = width - 14
    csize   = max ( int ( length / nchunks ) , 1 ) 

    left   = "[ "
    right  = " ]"
    symbol = allright ( symbol )

    fun = Ostap.Utils.frame_progress ( csize , nchunks , symbol , ' ' , left , right )
    cnt.OnPartialResultSlot  ( csize , fun )

    return cnt

# =============================================================================
## Get the effective entries in data frame
#  @code
#  data = ...
#  neff = data.nEff('b1*b1')
#  @endcode
def frame_nEff ( self , cuts = '' ) :
    """Get the effective entries in data frame 
    >>> data = ...
    >>> neff = data.nEff('b1*b1')
    """
    frame = DataFrame ( self ) 
    return Ostap.StatVar.nEff ( frame , cuts )

# =============================================================================
## Get statistics for the  given expression in data frame
#  @code
#  data = ...
#  c1 = data.statVar( 'S_sw' , 'pt>10' ) 
#  c2 = data.statVar( 'S_sw' , 'pt>0'  )
#  @endcode
def frame_statVar ( self , expression ,  cuts = '' ) :
    """Get statistics for the  given expression in data frame
    >>> data = ...
    >>> c1 = data.statVar( 'S_sw' , 'pt>10' ) 
    >>> c2 = data.statVar( 'S_sw' )
    """
    frame = DataFrame ( self ) 
    return Ostap.StatVar.statVar ( frame , expression , cuts , )


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
def frame_statCov ( self        ,
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

    frame  = DataFrame ( self ) 
    length = Ostap.StatVar.statCov ( frame       ,
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
def _fr_statVar_new_ ( self , expressions , cuts = '' , lazy = False  ) :
    """Get statistics of variable(s)
    """

    frame  = DataFrame ( self ) 

    input_string = False 
    if isinstance ( expressions , string_types ) :
        input_string = True 
        expressions = [ expressions ]
    
    ## get the list of currently known names
    vars = tuple ( [ str( c ) for c in frame.GetColumnNames () ] )

    names   = {}

    current = frame    
    for e in expressions :
        
        if e in vars :
            names [ e ] = e 
            continue

        used    = vars + tuple ( current.GetDefinedColumnNames() )
        vn      = var_name ( 'var_' , vars , e , *vars )
        current = current.Define ( vn , e )
        names [ e ] = vn

    cname = cuts 
    if cuts and not cuts in vars :
        used    = vars + tuple ( [ str(c) for c in current.GetDefinedColumnNames() ] )        
        vn      = var_name ( 'cut_' , used , cuts , *vars )
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


# =============================================================================
## Simplified print out for the  frame 
#  @code 
#  frame = ...
#  print frame
#  @endcode 
def frame_print ( t ) :
    """Simplified print out for the  frame 
    
    >>> frame = ...
    >>> print frame
    """
    ##
    res = "DataFrame Enries/#%d" %  len ( t )  
    ##
    _c          = [ str(c) for c in t.GetColumnNames() ] 
    _c.sort ()  
    res        += "\nColumns:\n%s" % multicolumn ( _c , indent = 2 , pad = 1 )
    return res


# ==============================================================================
## Data frame as table
#  @code
#  frame = ...
#  table = frame_table ( frame , '.*PT.*' , cuts = ... , more_vars = [ 'x*x/y' , 'y+z'] )
#  print ( table ) 
#  @endcode 
def frame_table ( self , pattern = None ,  cuts = '' , more_vars = () , prefix = '' ) :
    """Data frame as table
    >>> frame = ...
    >>> table = frame_table ( frame , '.*PT.*' , cuts = ... , more_vars = [ 'x*x/y' , 'y+z'] )
    >>> print ( table )
    """

    frame = DataFrame ( self )
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

    svars = vcols + tuple ( more_vars ) 
    stats = fr_statVars ( frame , svars , cuts = cuts , lazy = False ) 
    
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
    fmt_eff       =  '%8.3g +- %-8.3g'
    fmt_cumulated =  '%8.3g +- %-8.3g'
        
    header = ( ( '{:^%d}' % lmax ).format ( 'Filter'   ) ,               
               ( '{:>10}'        ).format ( '#input '  ) ,
               ( '{:<10}'        ).format ( ' #passed' ) ,
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
    return T.table ( table_data , title , prefix )

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


# ==============================================================================
## get 1D histogram model from the actual  1D histogram
#  @code
#  histo = ...
#  model = histo.model() 
#  @endcode 
def _th1_model ( histo ) :
    """Get 1D histogram model from the actual  1D histogram
    >>> histo = ...
    >>> model = histo.model() 
    """
    htmp = histo 
    if not isinstance ( histo , ROOT.TH1D ) :
        htmp = ROOT.TH1D ()
        histo.Copy ( htmp ) 
    
    return ROOT.RDF.TH1DModel ( htmp )

# ==============================================================================
## get 2D histogram model from the actual  2D histogram
#  @code
#  histo = ...
#  model = histo.model() 
#  @endcode 
def _th2_model ( histo ) :
    """Get 2D histogram model from the actual 2D histogram
    >>> histo = ...
    >>> model = histo.model() 
    """
    htmp = histo 
    if not isinstance ( histo , ROOT.TH2D ) :
        htmp = ROOT.TH2D ()
        histo.Copy ( htmp ) 
    
    return ROOT.RDF.TH2DModel ( htmp )

# ==============================================================================
## get 3D histogram model from the actual 3D histogram
#  @code
#  histo = ...
#  model = histo.model() 
#  @endcode 
def _th3_model ( histo ) :
    """Get 3D histogram model from the actual 3D histogram
    >>> histo = ...
    >>> model = histo.model() 
    """
    htmp = histo 
    if not isinstance ( histo , ROOT.TH3D ) :
        htmp = ROOT.TH3D ()
        histo.Copy ( htmp ) 
    
    return ROOT.RDF.TH3DModel ( htmp )

ROOT.TH1 .model = _th1_model
ROOT.TH2 .model = _th2_model
ROOT.TH3 .model = _th3_model


## def model_1D ( model ) :

##     if   isinstance ( model , ROOT.RDF.TH1Dmodel )            : return model
##     elif isinstance ( model , ROOT.TH1 ) and 1 == model.dim() : return model.model() 

##     try :
##         if len ( model )
        
    
    

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

    frame  = DataFrame ( frame )
    
    if 1 <= len ( what ) <= 2 :
        ww = split_string ( what[0] , ' ,:;' )
        if 1 < len ( ww ) :
            ww.reverse() 
            what = tuple ( ww ) + tuple ( what [1:] ) 
        
    histo = None

    if   isinstance ( model , ROOT.TH3 ) and 3 == model.dim () : 
        assert 3 <= len ( what ) <= 4, 'Invalid number of parameters %s' % str ( what )
        
        m = model.model () 
        
        histo = model
        if not isinstance ( model , ROOT.TH3D ) :
            tmp  = ROOT.TH3D ()
            model.Copy ( tmp )
            model = tmp 
        model = ROOT.RDF.TH3DModel ( model )
                
    elif isinstance ( model , ROOT.TH2 ) and 2 == model.dim () :
        assert 2 <= len ( what ) <= 3, 'Invalid number of parameters %s' % str ( what )
        histo = model        
        if not isinstance ( model , ROOT.TH2D ) :
            tmp   = ROOT.TH2D ()
            model.Copy ( tmp )
            model = tmp 
        model = ROOT.RDF.TH3DModel ( model )

    elif isinstance ( model , ROOT.TH1 ) and 1 == model.dim () : 
        assert 1 <= len ( what ) <= 2, 'Invalid number of parameters %s' % str ( what )
        histo = model 
        if not isinstance ( model , ROOT.TH1D ) :
            tmp   = ROOT.TH1D ()
            model.Copy ( tmp )
            model = tmp 
        model = ROOT.RDF.TH1DModel ( histo ) 

    elif  isinstance ( model , ROOT.RDF.TH3DModel ) :
        assert 3 <= len ( what ) <= 4 , 'Invalid number of parameters %s' % str ( what )

    elif  isinstance ( model , ROOT.RDF.TH2DModel ) :
        assert 2 <= len ( what ) <= 3 , 'Invalid number of parameters %s' % str ( what )
        
    elif  isinstance ( model , ROOT.RDF.TH1DModel ) :
        assert 1 <= len ( what ) <= 2 , 'Invalid number of parameters %s' % str ( what )

    else :

        logger.error ('ERROR HERE!') 
        
    if histo : histo.Reset()
    
    ## get the list of currently known names
    vars = tuple ( ( str(c) for c in frame.GetColumnNames () ) ) 

    nvars = []

    current = frame 
    added   = False 
    for w in what :
        
        if   w in  vars : nvars.append ( w )
        elif w in nvars : nvars.append ( w )
        else :
            used    = vars + tuple ( nvars ) + tuple ( current.GetDefinedColumnNames() )
            ww      = var_name ( 'var_' , used , *what )
            current = current.Define ( ww , w )
            nvars.append ( ww )
            added = True 

    if   isinstance ( model , ROOT.RDF.TH3DModel ) : action = current.Histo3D ( model , *nvars )
    elif isinstance ( model , ROOT.RDF.TH2DModel ) : action = current.Histo2D ( model , *nvars )
    elif isinstance ( model , ROOT.RDF.TH1DModel ) : action = current.Histo1D ( model , *nvars )


    if not histo :
        return action

    histo += action.GetValue() 
    
    return histo 
    

# ==============================================================================
# decorate 
# ==============================================================================
if not hasattr ( DataFrame , '__len__') : DataFrame.__len__ = _fr_len_
DataFrame .__str__     = frame_print
DataFrame .__repr__    = frame_print

DataFrame .nEff        = frame_nEff
DataFrame .statVar     = frame_statVar
DataFrame .statCov     = frame_statCov
DataFrame .ProgressBar = frame_progress
DataFrame .progress    = frame_progress


from ostap.stats.statvars import  data_decorate 
data_decorate ( DataFrame )
del data_decorate

from ostap.fitting.dataset import ds_draw , ds_project
DataFrame.draw    = ds_draw
DataFrame.project = ds_project


_decorated_classes_ = (
    DataFrame   ,
    )

_new_methods_       = (
    DataFrame.__len__          ,
    DataFrame.__repr__         ,
    DataFrame.__str__          ,
    DataFrame.columns          , 
    DataFrame.branches         ,
    #
    DataFrame.ProgressBar      ,
    DataFrame.progress         ,
    #
    DataFrame.draw             , 
    DataFrame.project          ,
    #
    DataFrame.nEff             ,
    DataFrame.statVar          ,
    DataFrame.statCov          ,
    DataFrame.nEff             ,
    #
    DataFrame.get_moment       , 
    DataFrame.central_moment   , 
    DataFrame.mean             ,
    DataFrame.rms              ,
    DataFrame.skewness         ,
    DataFrame.kurtosis         ,
    DataFrame.quantile         ,
    DataFrame.median           ,
    DataFrame.quantiles        ,
    DataFrame.interval         ,
    DataFrame.terciles         ,
    DataFrame.quartiles        ,
    DataFrame.quintiles        ,
    DataFrame.deciles          ,
    )


# =============================================================================
if ( 6 , 25 ) <= root_info :
    frame_statVar       = _fr_statVar_new_
    frame_statVars      = _fr_statVar_new_
    DataFrame.statVars  = _fr_statVar_new_
    __all__ = __all__ + ( 'frame_statVars' , ) 
    _new_methods_      = _new_methods_ + ( DataFrame.statVars , ) 
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
