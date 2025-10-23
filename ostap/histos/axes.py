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
    'axis_from_bins'    , ## Create the axis from sequence of bins
    'axis_bins'         , ## ditto
    'make_axis'         , ## make axis from description 
    'axis_same_binning' , ## same binning for two axes ?
    'h1_axis'           , ## create histogram from axis 
)
#
# =============================================================================
from   ostap.core.meta_info           import root_info 
from   ostap.core.ostap_types         import ( sequence_types , sized_types   ,
                                               num_types      , integer_types )
from   ostap.utils.basic              import typename 
from   ostap.core.core                import hID , Ostap, valid_pointer  
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
## Create the axis from sequence of bins
#  @code
#  bins = [1,2,3,4,5,10,100]
#  axis = axis_from_bins ( bins ) 
#  @endcode
#  @see TAxis
def axis_from_bins ( bins ) :
    """ Create the axis from sequence of bins
    - see `ROOT.TAxis`
    >>> bins = [1,2,3,4,5,10,100]
    >>> axis = axis_from_bins ( bins ) 
    """
    assert  isinstance ( bins , sequence_types ) \
        and isinstance ( bins ,    sized_types ) and 2 < len ( bins ) , \
        "`bins' must be sequence&sized&2<=len!"
    
    ## ROOT>6.36 has proper constructor: 
    if ( 6 , 36 ) <= root_info : return ROOT.TAxis ( new_bins )
    else                        : 
        new_bins = array.array ( 'd' , new_bins )
        return ROOT.TAxis ( len ( new_bins ) - 1 , new_bins )

# =============================================================================
## alias     
axis_bins = axis_from_bins

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

# =============================================================================
## iterate over items in TAxis
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _axis_items_ ( axis ) :
    """ Iterate over items in axis     
    >>> axis = ...
    >>> for ibin,low,center,high in axis.    items() : ...    
    >>> for ibin,low,center,high in axis.iteritems() : ... ## ditto    
    """
    for i in axis :

        l = axis.GetBinLowEdge ( i )
        c = axis.GetBinCenter  ( i )
        u = axis.GetBinUpEdge  ( i )
        
        yield i,l,c,u

# =============================================================================
## get bin parameters : low- and up-edges 
def _axis_get_item_ ( axis , i ) :
    """ Get bin parameter: low- and up-edges
    >>> axis = ...
    >>> low,high = axis[1]    
    """
    if not 1 <= i <= axis.GetNbins() : raise IndexError
    
    low = axis.GetBinLowEdge ( i )
    up  = axis.GetBinUpEdge  ( i )
    
    return low, up 

# =============================================================================
## Iterator over axis bins 
#  @code
#  axis = ...
#  for xlow, xhigh in axis.bins () :
#  ...   
#  @endcode 
def _axis_bin_iterator_ ( axis ) :
    """ Iterator over axis bins 
    >>> axis = ...
    >>> for xlow, xhigh in axis.bins () :
    ...
    """
    nbins = axis.GetNbins () 
    bins  = axis.GetXbins ()
    if bins :
        for i in range ( nbins ) :
            yield bins [ i ] , bins [ i + 1 ]
    else :
        xmin  = axis.GetXmin  ()
        xmax  = axis.GetXmax  ()
        binw  = ( xmax - xmin ) / nbins
        for i in range ( 1 , nbins + 1 ) :
            left  = axis.GetBinLowEdge()
            right = axis.GetBinUpEdge ()
            yield lef , right 

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
        for i in range ( 1 , nbins ) : yield axis.GetBinUpEdge()        
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
    bins  = axis.GetXbins()
    ## 
    return  1 == nbins  or not bins 
    
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
        return axis_from_bins ( new_bins ) 
    
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
    if 0 == nbins % n : 
        return ROOT.TAxis ( nbins // n , axis.GetXmin() , axis.GetXmax() )
    
    new_bins = [ axis.GetXmin() ]
    for i in range ( 1 , nbins + 1 , n ) :
        if 1 == i : continue  
        new_bins.append ( axis.GetBinLowEdge ( i ) )
    new_bins.append ( axis.GetXmax() )
        
    return axis_from_bins ( new_bins ) 
        
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

    return axis_from_bins ( new_bins )


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

    amin = axis.GetXmin ()
    amax = axis.GetXmax ()
    
    assert no_min or ( isinstance ( xmin , num_types ) and amax > xmin ) , \
    "Invalid `xmin' %s/%s" % ( xmin , typename ( xmin ) )
    
    assert no_max or ( isinstance ( xmax , num_types ) and amin < xmax ) , \
        "Invalid `xmax' %s/%s" % ( xmax , typename ( xmax ) )
    
    ## no action 
    if   no_min and no_max : return axis           ## RETURN 
    elif no_min :
        xmin , xmax = amin , min ( xmax , amax )
        return _axis_range_ ( axis , xmin , xmax ) ## RETURN
    elif no_max :
        xmin , xmax = max ( xmin , amin ) , amax 
        return _axis_range_ ( axis , xmin , xmax ) ## RETURN

    assert xmin < xmax , "Invalid xmin/xmax: %s/%s" % ( xmin , xmax )

    ## no action!
    if xmin <= amin and amax <= xmax : return axis ## RETURN

    ## adjust ximn/xmax: 
    xmin = max ( xmin , amin )
    xmax = min ( xmax , amax )

    mnbin = axis.FindFixBin ( xmin )
    mxbin = axis.FindFixBin ( xmax )
    
    nbins = axis.GetNbins() 
    
    assert 1     <= mnbin <= nbins , "Invalid min #bin %s" % mnbin 
    assert mnbin <= mxbin <= nbins , "Invalid max #bin %s" % mxbin 

    new_bins = [ xmin ] 
    for i in range ( mnbin , mxbin + 1 ) :
        
        low  = axis.GetBinLowEdge ( i )                
        high = axis.GetBinUpEdge  ( i )
        if xmax <= low  : break 
        if xmin >= high : continue

        if   low  < xmin and high <= xmax : continue 

        last = new_bins [ -1 ] 
        if   last < low  and high <= xmax : new_bins.append ( low )         
        elif last < low  and xmax <= high :
            new_bins .append ( low  )
            new_bins .append ( xmax )
            break

    if not isequal ( xmax , new_bins [ -1 ] ) : new_bins .append ( xmax ) 
    assert 2 <= len ( new_bins  ) , "Invalid edges!"

    return axis_from_bins ( new_bins )


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
    if not isinstance ( axis , ROOT.TAxis ) :
        axis = axis_from_bins ( axis )

    ## the actual histogram type: 
    if isinstance ( double , type ) and issubclass ( double , ROOT.TH1 ) : htype = double
    elif double : htype = ROOT.TH1D  
    else        : htype = ROOT.TH1F
    #
    title = title if title else 'Histogram:%s' % name 
    
    if axis.uniform () :
        h1 = htype ( name  ,
                     title ,
                     axis.GetNbins() , 
                     axis.GetXmin () , 
                     axis.GetXmax () )
    else : 
        h1 = htype ( name           ,
                     title          ,
                     axis.nBins ()  ,
                     array.array ( 'd' , axis.GetXbins() ) )         
    ##
    if not h1.GetSum2() : h1.Sumw2()
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
    if not isinstance ( x_axis , ROOT.TAxis ) : x_axis = axis_from_bins ( x_axis )
    if not isinstance ( y_axis , ROOT.TAxis ) : y_axis = axis_from_bins ( y_axis )

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
def h3_axes ( x_axis            ,
              y_axis            ,
              z_axis            ,
              title  = '3D'     , 
              name   = None     ,
              double = True     ) :
    """ Make 3D-histogram with binning deifned by already created axes    
    >>> x_axis = ...
    >>> y_axis = ...
    >>> z_axis = ...
    >>> h3 = h3_axes ( x_axis , y_axis , z_axis , title = 'MyHisto' )     
    """
    #
    if not name : name = hID() 
    #
    
     ## create the axes
    if not isinstance ( x_axis , ROOT.TAxis ) : x_axis = axis_from_bins ( x_axis )
    if not isinstance ( y_axis , ROOT.TAxis ) : y_axis = axis_from_bins ( y_axis )
    if not isinstance ( z_axis , ROOT.TAxis ) : z_axis = axis_from_bins ( z_axis )

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


ROOT.TAxis . __iter__      = _axis_iterator_
ROOT.TAxis . __reversed__  = _axis_iterator_reversed_
ROOT.TAxis . __contains__  = lambda s , i : 1 <= i <= s.GetNbins()
ROOT.TAxis . __getitem__   = _axis_get_item_
ROOT.TAxis . __eq__        =                  _axis_equal_
ROOT.TAxis . __ne__        = lambda a,o : not _axis_equal_ ( a , o ) 

ROOT.TAxis . bin_iterator  = _axis_bin_iterator_
ROOT.TAxis . bin_edges     = _axis_bin_iterator_

ROOT.TAxis. edges          = _axis_edges_
ROOT.TAxis. edge_iterator  = _axis_edges_

ROOT.TAxis. uniform        = _axis_uniform_

ROOT.TAxis.     items      = _axis_items_
ROOT.TAxis. iteritems      = _axis_items_

## change units
ROOT.TAxis . scale        = _axis_scale_
ROOT.TAxis . __mul__      = _axis_scale_
ROOT.TAxis . __rmul__     = _axis_scale_

## finer binning
ROOT.TAxis. split          = _axis_split_ 
ROOT.TAxis. __div__        = _axis_split_ 
ROOT.TAxis. __truediv__    = _axis_split_ 
ROOT.TAxis. __floordiv__   = _axis_split_ 

## coarse binnig 
ROOT.TAxis . merge         = _axis_merge_ 
ROOT.TAxis . __mod__       = _axis_merge_ 

## change the range 
ROOT.TAxis.range           = _axis_range_


## create the historgam for the axis 
ROOT.TAxis.histo           = h1_axis
ROOT.TAxis.histogram       = h1_axis



##  same bining ? 
ROOT.TAxis.same_binning = axis_same_binning







# =============================================================================



if '__main__' == __name__ :
            
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger ) 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
