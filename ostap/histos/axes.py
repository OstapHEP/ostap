#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/histos/axes.py 
#  Module with decoration of ROOT.TAxis objects for efficient use in python
#  @ser TAxis 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
""" Module with decoration of ROOT.TAxis objects for efficient use in python
- see `ROOT.TAxis`
"""
# =============================================================================
__version__ = "$Revision: 207180 $"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
# =============================================================================
__all__     = (
    'axis_from_edges'   , ## Create the axis from sequence of bins
    'make_axis'         , ## make axis from description 
    'axis_same_binning' , ## same binning for two axes ?
    'h1_axis'           , ## create histogram from axis 
)
#
# =============================================================================
from   ostap.core.meta_info           import root_info 
from   ostap.core.ostap_types         import ( sequence_types  , sized_types   ,
                                               num_types       , integer_types )
from   ostap.utils.basic              import typename 
from   ostap.core.core                import hID , Ostap, valid_pointer, rootException
from   ostap.math.base                import isequal 
from   ostap.logger.pretty            import fmt_pretty_values   
from   ostap.logger.symbols           import times 
import ROOT, array 
# =============================================================================
# logging
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.histos.axes' )
else                       : logger = getLogger( __name__          )
# =============================================================================
logger.debug ( "Decoration of histogram' axes")
# =============================================================================
## Create the axis from sequence of bin edges 
#  @code
#  bins = [1,2,3,4,5,10,100]
#  axis = axis_from_edges ( bins ) 
#  @endcode
#  @see TAxis
def axis_from_edges ( edges , check = True  ) :
    """ Create the axis from sequence of bin edges 
    - see `ROOT.TAxis`
    >>> edges = [1,2,3,4,5,10,100]
    >>> axis  = axis_from_edges ( edges ) 
    """
    assert isinstance ( edges , sequence_types ) , "`edges' must be sequence"  
    
    if isinstance ( edges , sized_types ) :
        assert 2 <= len ( edges ) , "at least two edges are required!"
        
    if check and root_info < ( 6 , 37 ) : 
        edges  = tuple ( e for e in edges )
        nedges = len ( edges ) 
        assert 2 <= nedges , "at least two edges are required!"
        ## the fist elements
        e0   = edges [ 0 ]        
        for ei in edges [ 1 : ] :
            if e0 < ei : raise TypeError ( "sequence of axis edges *must* be increasing!" )
            e0 = ei

    ## convert posible ROOT-Error to python exception 
    with rootException () :        
        ## recent ROOT (>=6.36) has proper constructor
        if ( 6 , 36 ) <= root_info : return ROOT.TAxis ( edges )
        ## 
        new_edges = array.array ( 'd' , edges )
        return ROOT.TAxis ( len ( new_edges ) - 1 , new_edges )

# =============================================================================
## Same binning for two axes?
#  @code
#  axis = ...
#  same = axis.same_binning ( another_axis )
#  @endcode 
#  @see Ostap::Utils::same_binning
def axis_same_binning ( axis , other ) :
    """ Same binning? 
    >>> axis = ...
    >>> same = axis.same_binning ( another_axis )
    - see `Ostap.Utils.same_binning`
    """
    assert isinstance ( other , ROOT.TAxis ) and valid_pointer ( other ) , \
        "same_binning: Invalid axis type %s" % typename ( other )
    return Ostap.Utils.same_binning ( axis , other )

# =============================================================================
## iterator for histogram  axis 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _axis_iterator_ ( axis ) :
    """ Iterator for axis bins 
    >>> axis = ...
    >>> for i in axis : 
    """
    nbins  = axis.GetNbins()
    for i in range ( 1 , nbins + 1 ) : yield i

# =============================================================================
## Reversed iterator for histogram  axis 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _axis_iterator_reversed_ ( axis ) :
    """ Reversed iterator for axis
    >>> axis = ...
    >>> for i in reverse ( axis ) : 
    """
    nbins  = axis.GetNbins()
    for i in range ( nbins , 0 , -1 ) : yield i

# ==============================================================================
## Get the edge by index
#  @code
#  axis   = ...
#  edge1 = axis.edge ( 1 )
#  @endcode
def _axis_edge_ ( axis , index )  :
    """ Get the edge by index
    >>> axis   = ...
    >>> edge1 = axis.edge ( 1 )
    """
    nn = axis.GetNbins() + 1 
    if index < 0 : index += nn 
    if not 0 <= index < nn : raise IndexError() 
    ## 
    if 0  == index     : return axis.GetXmin() 
    if nn == index + 1 : return axis.GetXmax() 
    #
    bins = axis.GetXbins() 
    if bins : return bins[index]
    #
    return axis.GetBinLowEdge ( index + 1 )
    
# =============================================================================
## iterate over items in TAxis
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _axis_items_ ( axis ) :
    """ Iterate over items in axis     
    >>> axis = ...
    >>> for ibin , low , center, up in axis.    items() : ...    
    >>> for ibin , low , center, up in axis.iteritems() : ... ## ditto    
    """
    for index in axis :

        low    = axis.edge         ( index - 1 )
        up     = axis.edge         ( index     )
        center = axis.GetBinCenter ( index     )
        
        yield index , low , center, up  

# =============================================================================
## get bin parameters : low- and up-edges 
def _axis_get_item_ ( axis , index ) :
    """ Get bin parameter: low- and up-edges
    >>> axis = ...
    >>> low , center , high = axis[1]    
    """
    
    nbins = axis.GetNbins () 
    
    if isinstance ( index , slice ) :
        
        start = index.start 
        stop  = index.stop 
        step  = index.step 
        
        int1 = isinstance ( start , integer_types )
        int2 = isinstance ( stop  , integer_types )
        int3 = isinstance ( step  , integer_types ) and 0 < step 
            
        none1 = start is None  
        none2 = stop  is None  
        none3 = step  is None  
            
        ## integer slices 
        if ( int1 or none1 ) and ( int2 or none2 ) and ( int3 or none3 ) :
            indices = index.indices ( nbins ) 
            indices = range  ( *indices )
            assert  2 <= len (  indices ) , "At least two edges are required!" 
            return axis_from_edges ( ( axis.edge ( i ) for i in indices ) ) 
        
        num1 = isinstance ( start , num_types )
        num2 = isinstance ( stop  , num_types )
        
        ## numeric ranges 
        if none3 and ( num1 or num2 ) :
            return _axis_range_ ( axis , start , stop )
        
        raise IndexError ( "Invalid index %s/%s" % ( index , typename ( index ) ) ) 
     
    if not isinstance ( index , integer_types ) : 
        raise IndexError ( "Invalid type of index %s" % typename ( index ) ) 
    
    if   nbins < index : raise IndexError 
    elif index < 0     : index += nbins + 1 
     
    if not 1 <= index <= axis.GetNbins() : raise IndexError
    
    low    = axis.edge         ( index - 1 )
    up     = axis.edge         ( index     )
    center = axis.GetBinCenter ( index     )
    
    item = low , center , up 
    return item 

# =============================================================================
## Iterator over axis bins 
#  @code
#  axis = ...
#  for xlow, xcenter, xhigh in axis.bins () :
#  ...   
#  @endcode 
def _axis_bin_iterator_ ( axis ) :
    """ Iterator over axis bins 
    >>> axis = ...
    >>> for xlow, xcenter, xhigh in axis.bins () :
    ...
    """
    nbins = axis.GetNbins () 
    for i in range ( 1 , nbins + 1 ) :
        yield axis[ i ]

# =============================================================================
## equality for axes
def _axis_equal_ ( axis , another ) :
    """ Equality for two axes
    >>> a1 = ...
    >>> a2 = ...
    >>> print ( a1 == a1 ) 
    """
    if      axis   is       another                        : return True 
    if len( axis ) != len ( another )                      : return False
    if not isequal ( axis.GetXmin()  , another.GetXmin() ) : return False
    if not isequal ( axis.GetXmax()  , another.GetXmax() ) : return False
    ## 
    return axis.same_binning ( another )

# =============================================================================
## Iterator for bin edges
#  @code
#  axis = ...
#  for edge in axis.edges() :
#  ...
#  @endcode 
def _axis_edges_ ( axis ) :
    """ Iterator for bin edges
    >>> axis = ...
    >>> for edge in axis.edges() :
    >>> ...
    """
    nbins = axis.GetNbins ()
    bins  = axis.GetXbins ()

    if bins :
        
        for edge in bins : yield edge 
        
    else :
        
        yield axis.GetXmin()        
        for i in range ( 1 , nbins ) : yield axis.GetBinUpEdge ( i )        
        yield axis.GetXmax()

# ==============================================================================
## Uniform bining ?
#  @code
#  axis        = ...
#  uniform_bins = axis.uniform () 
#  @endcode 
def _axis_uniform_ ( axis ) :
    """ Uniform binning
    >>> axis        = ...
    >>> uniform_bins = axis.uniform () 
    """
    nbins = axis.GetNbins()
    if 1 == nbins  : return True 
    bins  = axis.GetXbins()
    if not bins    : return True 
    ##
    return Ostap.Utils.uniform_bins ( axis , True ) 

# ==============================================================================
## Scale axis edges (equivalent to change units)
#  @code
#  axis = ...
#  scaled = axis.scale ( 1000 )
#  scaled = axis * 1000 
#  @endcode
def _axis_scale_ ( axis , scale ) :
    """ Scale axis edges (equivalent to change units)
    >>> axis = ...
    >>> caled = axis.scale ( 1000 )
    >>> scaled = axis * 1000 
    """
    assert isinstance ( scale , num_types ) and 0 < scale , \
           'Invalid "scale" argument %s' % scale
    
    bins = axis.GetXbins ()
    ## non-uniform bins are here 
    if bins and 2 <= len ( bins ) :
        new_bins = tuple ( v * scale for v in bins )
        return axis_from_edges ( new_bins ) 
    
    return ROOT.TAxis ( axis.GetNbins ()         ,
                        axis.GetXmin  () * scale ,
                        axis.GetXmax  () * scale )

# =============================================================================
## Merge the axis bins into axis with wider bin
#  - each group of n-bin is merged into a single bin
#  @code
#  axis  = ...
#  axis2 = axis.merge ( 2 )
#  axis3 = axis.merge ( 3 )
#  axis2 = axis % 2         ## ditto 
#  axis3 = axis % 3         ## ditto 
#  @endcode 
def _axis_merge_ ( axis , n ) :
    """ Merge the axis bins into axis with wider bin
    - each group of n-bins is merged into single bin 
    >>> axis  = ...
    >>> axis2 = axis.merge ( 2 )
    >>> axis3 = axis.merge ( 3 )
    >>> axis2 = axis % 2          ## the same 
    """
    assert isinstance ( n , integer_types ) and 0 < n , \
        "Invalid `n' %s (must be positive integer)!" % n

    ## no action 
    if 1 == n  : return axis
    
    nbins = axis.GetNbins () 
    bins  = axis.GetXbins ()
    ## 
    if axis.uniform () and 0 == nbins % n : 
        return ROOT.TAxis ( nbins // n , axis.GetXmin() , axis.GetXmax() )
    
    new_bins = [ axis.GetXmin() ]
    for i in range ( 1 , nbins + 1 , n ) :
        if 1 == i : continue  
        new_bins.append ( axis.GetBinLowEdge ( i ) )
    new_bins.append ( axis.GetXmax() )
        
    return axis_from_edges ( new_bins ) 

# =============================================================================
## Join indicated bins
#  - each indicated bins is joined/merged with subsequent bin
#  @code
#  axis  = ...
#  axis2 = axis.join ( 2 ) ## join/merge 2&3rd bins
#  axis3 = axis.join ( 2 , 10 ) ## join/merge 2&3rd, and 10&11 bins 
#  @endcode 
def _axis_join_ ( axis , *bins ) :
    """ Join indicated bins
    - each indicated bins is joined/merged with subsequent bin
    >>> axis  = ...
    >>> axis2 = axis.join ( 2 ) ## join/merge 2&3rd bins
    >>> axis3 = axis.join ( 2 , 10 ) ## join/merge 2&3rd, and 10&11 bins 
    """
    ## no action ? 
    if not bins :  return axis
    
    nbins = axis.GetNbins ()
    
    ## nothing to be merged/joined 
    if 1 >= nbins : return axis 

    assert all ( isinstance ( e , integer_types ) and 1 <= e and e + 1 < nbins for e in bins ) , \
        "Invalid `bins': [%s]" %  ( ','.join ( str ( b ) for b in bins ) ) 

    ## aliminate duplicates, sort and reverse 
    bins  = list ( reversed ( sorted ( set ( b for b in bins ) ) ) )

    ## initial list of edges 
    edges = list ( e for e in axis.edges () )

    ## remove bins from edges 
    while bins :
        del edges [ bins.pop ()  ]

    return axis_from_edges ( edges ) 

# =============================================================================
## Split the axis into axis with narrower bin
#  - each bin is split into n bins
#  @code
#  axis  = ...
#  axis2 = axis.split ( 2 )
#  axis3 = axis.split ( 3 )
#  axis2 = axis / 2 
#  axis3 = axis / 3 
#  @endcode 
def _axis_split_ ( axis , n ) :
    """ Split the axis into axis with narrower bin
    - each bin is split into n bins 
    >>> axis  = ...
    >>> axis2 = axis.split ( 2 )
    >>> axis3 = axis.split ( 3 )
    >>> axis2 = axis / 2 
    >>> axis3 = axis / 3 
    """
    assert isinstance ( n , integer_types ) and 0 < n , \
        "Invalid `n' %s (must be positive integer)!" % n

    ## no action 
    if 1 == n  : return axis
    nbins = axis.GetNbins () 
    bins  = axis.GetXbins ()
    ## 
    ## uniform binning ?
    if not bins :
        return ROOT.TAxis ( nbins * n , axis.GetXmin() , axis.GetXmax() ) 
    ##
    prev     = axis.GetXmin()
    new_bins = [ prev ] 
    for i , current in enumerate ( bins ) :
        if prev < current :
            delta  = current - prev
            delta /= n
            for j in range ( 1 , n  ) :  new_bins.append ( prev + j * delta )
            new_bins.append ( current )
        prev = current

    return axis_from_edges ( new_bins )

# =============================================================================
## get axis as sub-axis with the given range
#  @code
#  axis = ..
#  axis1 = axis.range ( xmax = 3.0 ) 
#  axis2 = axis.range ( xmin = 1.0 ) 
#  axis3 = axis.range ( xmin = 1.0 , xmax = 3 ) 
#  @endcode
def _axis_range_ ( axis , xmin = None , xmax = None ) :
    """ Get axis as sub-axis with given range
    >>> axis = ..
    >>> axis1 = axis.range ( xmax = 3.0 ) 
    >>> axis2 = axis.range ( xmin = 1.0 ) 
    >>> axis3 = axis.range ( xmin = 1.0 , xmax = 3 ) 
    """
    no_min = xmin is None
    no_max = xmax is None

    ## no action
    if no_min and no_max : return axis           ## RETURN 

    amin   = axis.GetXmin ()
    amax   = axis.GetXmax ()
    
    assert no_min or isinstance ( xmin , num_types ) , \
        "Invalid `xmin' %s/%s" % ( xmin , typename ( xmin ) )
    
    assert no_max or isinstance ( xmax , num_types ) , \
        "Invalid `xmax' %s/%s" % ( xmax , typename ( xmax ) )

    ## no action 
    if   no_min and no_max : return axis           ## RETURN 
    elif no_min : return _axis_range_ ( axis , amin , xmax ) 
    elif no_max : return _axis_range_ ( axis , xmin , amax ) 

    assert xmin < xmax , "Invalid xmin/xmax: %s/%s" % ( xmin , xmax )

    new_edges = [] 
    new_edges.append ( xmin )    
    for e in axis.edges () :
        if xmin < e  < xmax : new_edges.append ( e )            
    new_edges.append ( xmax ) 

    return axis_from_edges ( new_edges )

# =============================================================================
## print TAxis 
def _axis_str_ ( axis , precision = 3 , width = 5 ) :
    """ Print TAxis
    """
    
    if axis.uniform () :
        xmin , xmax = axis.GetXmin () , axis.GetXmax() 
        fmt , expo = fmt_pretty_values ( xmin ,
                                         xmax ,  
                                         precision = precision , 
                                         width     = width     )
        if expo : 
            scale = 10 ** expo 
            v1 = '%s%s10^%d' % ( fmt % ( xmin / scale ) , times , expo ) 
            v2 = '%s%s10^%d' % ( fmt % ( xmax / scale ) , times , expo ) 
            return 'TAxix(%d,%s,%s)' % ( axis.GetNbins() , v1 , v2 )
        else : 
            return 'TAxix(%d,%s,%s)' % ( axis.GetNbins() , fmt % xmin , fmt % xmax  )
    
    edges = tuple ( e for e in axis.edges() )
    fmt, expo = fmt_pretty_values ( *edges , precision = precision , width = width )
    if expo : 
        scale = 10 ** expo 
        edges = ( fmt % ( v / scale ) for v in edges )
        edges = ','.join ( edges )
        expo  = '%s10^%d' % ( times , expo )
        return 'TAxis([%s]%s)' % ( edges , expo )
    
    edges = ( fmt %  v for v in edges )
    edges = ','.join ( edges )
    return 'TAxis([%s])' % edges 

    
# =============================================================================
## Create <code>TAxis</code>
#  @code
#  a = make_axis ( 2 , 0.0 , 1.0 )                  ## #bins ,min, max 
#  a = make_axis ( 4 , 0.0 , 1.0 , 2.0 , 5.0 , 10 ) ## #bins , e1 , e2 , e3 , ...    
#  a = make_axis ( 2 , [ 0.0 , 1.0 , 2.0 ]  )       ## #bins , [ e1 , e2 . ... ]  
#  a = make_axis ( [ 0.0 , 1.0 , 2.0 ]  )           ## [ e1 , e2 , ... ] 
#  a = make_axis ( 0.0 , 1.0 , 2.0  )               ##   e1 , e2 , ...  
#  @endcode
def make_axis ( nbins , *bins ) :
    """ Create <code>TAxis</code>    
    >>> a = make_axis ( 2 , 0.0 , 1.0 )                  ## #bins ,min, max 
    >>> a = make_axis ( 4 , 0.0 , 1.0 , 2.0 , 5.0 , 10 ) ## #bins , e1 , e2 , e3 , ...    
    >>> a = make_axis ( 2 , [ 0.0 , 1.0 , 2.0 ]  )       ## #bins , [ e1 , e2 . ... ]  
    >>> a = make_axis ( [ 0.0 , 1.0 , 2.0 ]  )           ## [ e1 , e2 , ... ] 
    >>> a = make_axis (  0.0 , 1.0 , 2.0 , 3.0 , 4.0 )   ##   e1 , e2 , ... 
    """
    if   isinstance ( nbins , ROOT.TAxis ) : return nbins 
    
    if isinstance   ( nbins , integer_types ) and 0 < nbins :
        
        if 1 == len ( bins ) and isinstance ( bins [ 0 ] , sequence_types ) :
            abins = [ float ( e ) for e in bins [ 0 ] ]
            if nbins +1 <= len ( abins ) :
                abins = abins [ : nbins + 1 ] 
                if is_sorted ( abins ) : 
                    return ROOT.TAxis ( nbins , array.array ( 'd' , abins ) )
                
        elif 2 == len ( bins ) and \
               isinstance ( bins [ 0 ] , num_types ) and \
               isinstance ( bins [ 1 ] , num_types ) and bins [ 0 ] < bins [ 1 ] :
            return ROOT.TAxis ( nbins , bins [ 0 ] , bins [ 1 ] )
        
        elif nbins + 1 <= len ( bins ) :                    
            abins = [ float ( e ) for e in bins [ : nbins + 1 ] ]
            if is_sorted ( abins ) : 
                return ROOT.TAxis ( nbins , array.array ( 'd' , abins ) ) 

    elif isinstance ( nbins , sequence_types ) and not bins :
        
        abins = [ float ( e ) for e in nbins ]
        if 2 <= len ( abins ) and is_sorted ( abins ) : 
            return ROOT.TAxis ( nbins , array.array ( 'd' , abins ) )
        
    if isinstance ( nbins , num_types ) and bins : 
        abins = [ float ( nbins ) ] + [ float ( e ) for e in bins ]
        if is_sorted ( abins ) : 
            return ROOT.TAxis ( nbins , array.array ( 'd' , abins ) ) 

    raise ArgumentError('make_axis: invalid arguments %s' % str ( ( nbins , ) + bins ) ) 

# =============================================================================
## make 1D-histogram from the axis
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def h1_axis ( axis               ,
              title  = '1D'      , 
              name   = None      ,
              double = ROOT.TH1D ) :
    """ Make 1D-histogram with binning defined by already created axes    
    >>> axis = ...
    >>> h1 = h1_axes ( axis , title = 'MyHisto' )     
    """
    #
    if not name : name = hID ()
    #
    ## create the axis 
    if not isinstance ( axis , ROOT.TAxis ) : axis = axis_from_edges ( axis )

    ## the actual histogram type: 
    if isinstance ( double , type ) and issubclass ( double , ROOT.TH1 ) : htype = double
    elif double : htype = ROOT.TH1D  
    else        : htype = ROOT.TH1F
    #
    title = title if title else 'Histogram:%s' % name 
    
    if axis.uniform () :
        h1 = htype ( name  ,
                     title ,
                     axis.GetNbins()   , 
                     axis.GetXmin ()   , 
                     axis.GetXmax ()   )
    else : 
        h1 = htype ( name              ,
                     title             ,
                     axis.GetNbins ()  ,
                     array.array ( 'd' , axis.GetXbins() ) )         
    ##
    if not h1.GetSumw2() : h1.Sumw2()
    return h1


# =============================================================================
## make 2D-histogram from axes
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def h2_axes ( x_axis             ,
              y_axis             ,
              title  = '2D'      , 
              name   = None      ,
              double = ROOT.TH2D ) :
    """ Make 2D-histogram with binning deifned by already created axes    
    >>> x_axis = ...
    >>> y_axis = ...
    >>> h2 = h2_axes ( x_axis , y_axis , title = 'MyHisto' )     
    """
    #
    if not name : name = hID() 
    #
    ## create the axes
    if not isinstance ( x_axis , ROOT.TAxis ) : x_axis = axis_from_edges ( x_axis )
    if not isinstance ( y_axis , ROOT.TAxis ) : y_axis = axis_from_edges ( y_axis )

    ## the actual histogram type: 
    if isinstance ( double , type ) and issubclass ( double , ROOT.TH2 ) : htype = double
    elif double : htype = ROOT.TH2D
    else        : htype = ROOT.TH2F

    title = title if title else 'Histogram:%s' % name 
    
    xu = x_axis.uniform()
    yu = y_axis.uniform()
    
    if xu and yu :
        
        h2 = htype ( name  ,
                     title ,
                     x_axis.GetNbins() , x_axis.GetXmin () , x_axis.GetXmax () ,
                     y_axis.GetNbins() , y_axis.GetXmin () , y_axis.GetXmax () )
    
    else :

        xbins = array.array ( 'd' , ( e for e in x_axis.edges() ) )
        ybins = array.array ( 'd' , ( e for e in y_axis.edges() ) )

        h2 = htype ( name  ,
                     title ,
                     x_axis.GetNbins() , xbins , 
                     y_axis.GetNbins() , ybins )

    if not h2.GetSumw2 : h2.Sumw2()
    return h2

# =============================================================================
## make 3D-histogram from axes
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-07-18
def h3_axes ( x_axis             ,
              y_axis             ,
              z_axis             ,
              title  = '3D'      , 
              name   = ''        ,
              double = ROOT.TH3D ) :
    """ Make 3D-histogram with binning deifned by already created axes    
    >>> x_axis = ...
    >>> y_axis = ...
    >>> z_axis = ...
    >>> h3 = h3_axes ( x_axis , y_axis , z_axis , title = 'MyHisto' )     
    """
    #
    if not name : name = hID ()
    else        : name = hID ( prefix = name )
    #

     ## create the axes
    if not isinstance ( x_axis , ROOT.TAxis ) : x_axis = axis_from_edges ( x_axis )
    if not isinstance ( y_axis , ROOT.TAxis ) : y_axis = axis_from_edges ( y_axis )
    if not isinstance ( z_axis , ROOT.TAxis ) : z_axis = axis_from_edges ( z_axis )

    ## the actual histogram type: 
    if isinstance ( double , type ) and issubclass ( double , ROOT.TH3 ) : htype = double
    elif double : htype = ROOT.TH3D
    else        : htype = ROOT.TH3F

    title = title if title else 'Histogram:%s' % name 
    
    xu = x_axis.uniform()
    yu = y_axis.uniform()
    zu = z_axis.uniform()
    
    if xu and yu and zu :
        
        h3 = htype ( name  ,
                     title ,
                     x_axis.GetNbins() , x_axis.GetXmin () , x_axis.GetXmax () ,
                     y_axis.GetNbins() , y_axis.GetXmin () , y_axis.GetXmax () ,
                     z_axis.GetNbins() , z_axis.GetXmin () , z_axis.GetXmax () )
    
    else :

        xbins = array.array ( 'd' , ( e for e in x_axis.edges() ) )
        ybins = array.array ( 'd' , ( e for e in y_axis.edges() ) )
        zbins = array.array ( 'd' , ( e for e in z_axis.edges() ) )

        h3 = htype ( name  ,
                     title ,
                     x_axis.GetNbins() , xbins , 
                     y_axis.GetNbins() , ybins ,
                     z_axis.GetNbins() , zbins )

    if not h3.GetSumw2 : h3.Sumw2()
    return h3

# =============================================================================
## Create 2D or 3D historgam from axes
#  @code
#  xaxis = ...
#  yaxis = ...
#  h2 = xaxis @ yaxis ## 2D histogram 
#
#  zaxis = ...
#  h3 = h2    @ zaxis
#  endcode
def _axis_matmul_ ( axis , another ) :
    """ Create 2D or 3D historgam from axes
    >>> xaxis = ...
    >>> yaxis = ...
    >>> h2 = xaxis @ yaxis ## 2D historam 
    
    >>> zaxis = ...
    >>> h3 = h2    @ xaxis  ## 3D histogram 
    """
    if   isinstance ( another , ROOT.TH2   ) and 2 == another.GetDimension () :
        return h3_axes ( axis , another.GetXaxis ()  , another.GetYaxis () )
    elif isinstance ( another , ROOT.TH1   ) and 1 == another.GetDimension () :
        return h2_axes ( axis , another.GetXaxis () ) 
    elif isinstance ( another , ROOT.TAxis ) : return h2_axes ( axis , another )
    ## 
    return NotImplemented

# =============================================================================
## Create 2D or 3D historgam from axes
#  @code
#  xaxis = ...
#  yaxis = ...
#  h2 = xaxis @ yaxis ## 2D histogram 
#
#  zaxis = ...
#  h3 = h2    @ zaxis
#  endcode
def _axis_rmatmul_ ( axis , another ) :
    """ Create 2D or 3D histogram from axes
    >>> xaxis = ...
    >>> yaxis = ...
    >>> h2 = xaxis @ yaxis ## 2D historam 
    
    >>> zaxis = ...
    >>> h3 = h2    @ xaxis  ## 3D histogram 
    """
    if   isinstance ( another , ROOT.TH2   ) and 2 == another.GetDimension ()  :
        return h3_axes ( another.GetXaxis () , another.GetYaxis () , axis )
    elif isinstance ( another , ROOT.TH1   ) and 1 == another.GetDimension () :
        return h2_axes ( another.GetXaxis () , axis ) 
    elif isinstance ( another , ROOT.TAxis ) : return h2_axes ( another , axis )
    ## 
    return NotImplemented

ROOT.TAxis .from_edges     = staticmethod ( axis_from_edges )
ROOT.TAxis .from_bins      = staticmethod ( axis_from_edges )

ROOT.TAxis . __iter__      = _axis_iterator_
ROOT.TAxis . __reversed__  = _axis_iterator_reversed_
ROOT.TAxis . __contains__  = lambda s , i : 1 <= i <= s.GetNbins()
ROOT.TAxis . __getitem__   = _axis_get_item_
ROOT.TAxis . __eq__        =                  _axis_equal_
ROOT.TAxis . __ne__        = lambda a,o : not _axis_equal_ ( a , o ) 

ROOT.TAxis . __str__       = _axis_str_ 
ROOT.TAxis . __repr__      = _axis_str_ 

ROOT.TAxis . __len__       = lambda s : s.GetNbins ()

ROOT.TAxis . bin_iterator  = _axis_bin_iterator_
ROOT.TAxis . bin_edges     = _axis_bin_iterator_

ROOT.TAxis. edges          = _axis_edges_
ROOT.TAxis. edge_iterator  = _axis_edges_
ROOT.TAxis. edge           = _axis_edge_ 

ROOT.TAxis. uniform        = _axis_uniform_

ROOT.TAxis.     items      = _axis_items_
ROOT.TAxis. iteritems      = _axis_items_

## change units
ROOT.TAxis . scale         = _axis_scale_
ROOT.TAxis . __mul__       = _axis_scale_
ROOT.TAxis . __rmul__      = _axis_scale_

## finer binning
ROOT.TAxis. split          = _axis_split_ 
ROOT.TAxis. __div__        = _axis_split_ 
ROOT.TAxis. __truediv__    = _axis_split_ 
ROOT.TAxis. __floordiv__   = _axis_split_ 

## coarse binnig 
ROOT.TAxis . merge         = _axis_merge_ 
ROOT.TAxis . __mod__       = _axis_merge_ 

## multiply axes to get histograms 
ROOT.TAxis. __matmul__     = _axis_matmul_
ROOT.TAxis. __rmatmul__    = _axis_rmatmul_

## join certain bins
ROOT.TAxis . join          = _axis_join_ 

## change the range 
ROOT.TAxis.range           = _axis_range_

## create the histogram for the axis 
ROOT.TAxis.histo           = h1_axis
ROOT.TAxis.histogram       = h1_axis
ROOT.TAxis.histo1          = h1_axis
ROOT.TAxis.histogram1      = h1_axis
ROOT.TAxis.histo2          = h2_axes
ROOT.TAxis.histogram2      = h2_axes
ROOT.TAxis.histo3          = h3_axes
ROOT.TAxis.histogram3      = h3_axes
ROOT.TAxis.histo1d         = h1_axis
ROOT.TAxis.histogram1d     = h1_axis
ROOT.TAxis.histo2d         = h2_axes
ROOT.TAxis.histogram2d     = h2_axes
ROOT.TAxis.histo3d         = h3_axes
ROOT.TAxis.histogram3d     = h3_axes


##  same bining ? 
ROOT.TAxis.same_binning = axis_same_binning


_decorated_classes_ = (
    ROOT.TAxis  ,
)
_new_methods_  = (
    #
    ## create from bins 
    ROOT.TAxis .from_edges     , 
    ROOT.TAxis .from_bins      , 
    #
    ## basic operations
    ROOT.TAxis . __iter__      , 
    ROOT.TAxis . __reversed__  ,
    ROOT.TAxis . __contains__  , 
    ROOT.TAxis . __getitem__   , 
    ROOT.TAxis . __eq__        , 
    ROOT.TAxis . __ne__        , 
    ## 
    ROOT.TAxis . __str__       , 
    ROOT.TAxis . __repr__      , 
    ## 
    ROOT.TAxis . __len__       , 
    ## 
    ROOT.TAxis . bin_iterator  , 
    ROOT.TAxis . bin_edges     , 
    ## 
    ROOT.TAxis. edges          , 
    ROOT.TAxis. edge_iterator  , 
    ROOT.TAxis. edge           , 
    ## 
    ROOT.TAxis. uniform        , 
    ## 
    ROOT.TAxis.     items      , 
    ROOT.TAxis. iteritems      , 
    # 
    ## change units
    ROOT.TAxis . scale         , 
    ROOT.TAxis . __mul__       , 
    ROOT.TAxis . __rmul__      , 
    #
    ## finer binning
    ROOT.TAxis. split          , 
    ROOT.TAxis. __div__        , 
    ROOT.TAxis. __truediv__    , 
    ROOT.TAxis. __floordiv__   , 
    # 
    ## coarse binnig 
    ROOT.TAxis . merge         , 
    ROOT.TAxis . __mod__       , 
    # 
    ## multiply axes to get histograms 
    ROOT.TAxis. __matmul__     , 
    ROOT.TAxis. __rmatmul__    , 
    # 
    ## join certain bins
    ROOT.TAxis . join          , 
    # 
    ## change the range 
    ROOT.TAxis.range           , 
    #
    ## create the histogram from the axis/axes  
    ROOT.TAxis.histo           , 
    ROOT.TAxis.histogram       ,  
    ROOT.TAxis.histo1          ,
    ROOT.TAxis.histogram1      , 
    ROOT.TAxis.histo2          , 
    ROOT.TAxis.histogram2      , 
    ROOT.TAxis.histo3          , 
    ROOT.TAxis.histogram3      , 
    ROOT.TAxis.histo1d         , 
    ROOT.TAxis.histogram1d     , 
    ROOT.TAxis.histo2d         , 
    ROOT.TAxis.histogram2d     , 
    ROOT.TAxis.histo3d         , 
    ROOT.TAxis.histogram3d     , 
    ## 
    )
# =============================================================================



if '__main__' == __name__ :
            
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger ) 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
