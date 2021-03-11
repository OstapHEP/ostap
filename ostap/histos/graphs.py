#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  (T)Graph-related decorations 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07 
# =============================================================================
"""TGraph-related decorations:
 - makeGraph    : make graph from primitive data
 - makeGraph2   : make graph from plain input two-column text 
 - makeGraph3   : make graph with errors  from plain input four-column text 
 - makeGraphs3  : make graphs from plain input multicolumn text 
 - makeGraphs4' : make graphs from plain input multicolumn text 
 - hToGraph'    : convert histogram to graph 
 - hToGraph2'   : convert histogram to graph 
 - hToGraph3'   : convert histogram to graph
 - lw_graph'    : make Laffery-Wyatt's graph 
 - fill_area'   : create a graph for the area between two curves/functions
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'makeGraph'   , # make graph from primitive data
    'makeGraph2'  , # make graph from plain input two-column text 
    'makeGraph3'  , # make graph with errors  from plain input four-column text 
    'makeGraphs3' , # make graphs from plain input multicolumn text 
    'makeGraphs4' , # make graphs from plain input multicolumn text 
    'hToGraph'    , # convert histogram to graph 
    'hToGraph2'   , # convert histogram to graph 
    'hToGraph3'   , # convert histogram to graph
    'lw_graph'    , # make Laffery-Wyatt's graph 
    'fill_area'   , # create a graph for the area between two curves/functions
    ##
    ) 
# =============================================================================
import ROOT, ctypes        
from   ostap.core.core                import cpp, VE, grID
from   ostap.core.ostap_types         import num_types, integer_types  
from   builtins                       import range
from   ostap.plotting.draw_attributes import copy_graph_attributes  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.histos.graphs' )
else                       : logger = getLogger( __name__              )
# =============================================================================
logger.debug ( '(T)Graph-related decorations')
# =============================================================================
pos_infinity = float('+inf')
neg_infinity = float('-inf')
# =============================================================================
## make graph from data 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def makeGraph ( x , y = []  , ex = [] , ey = [] ) :
    """ Make graph using the primitive data
    """
    if  isinstance ( x , dict ) and not y and not ex and not ey : 
        _x = []
        _y = []
        keys = x.keys()
        for k in sorted ( keys ):
            _x += [   k  ]
            _y += [ x[k] ] 
        return makeGraph ( _x , _y )
        
    if  not x                   : raise TypeError ( "X is not a proper vector!" ) 
    if  not y                   : raise TypeError ( "Y is not a proper vector!" )
    if len( x ) != len ( y )    : raise TypeError ( "Mismatch X/Y-lengths" ) 

    if ex and len(ex) != len(x) : raise TypeError ( "Mismatch X/eX-lengths" )
    if ey and len(ey) != len(y) : raise TypeError ( "Mismatch Y/eY-lengths" )

    gr = ROOT.TGraphErrors ( len(x) ) 
        
    for i in range ( 0 , len(x) ) :
        
        if ex : _ex = ex[i]
        else  : _ex = 0.0
        if ey : _ey = ey[i]
        else  : _ey = 0.0

        _x = x[i]
        if hasattr ( x[i] , 'value' ) : _x  = x[i].value ()        
        if hasattr ( x[i] , 'error' ) : _ex = x[i].error ()

        _y = y[i]
        if hasattr ( y[i] , 'value' ) : _y  = y[i].value ()        
        if hasattr ( y[i] , 'error' ) : _ey = y[i].error ()
                    
        gr .SetPoint      ( i ,  _x ,  _y )
        gr .SetPointError ( i , _ex , _ey )
        
    return gr

# =============================================================================
## create TGraph from the plain text two-column format with optional comments.
#  E.g. the files provided by FONLL on-line calculator can be parsed
#  @see http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html
#  @code
#  gr = makeGraph2 ( '''
#  # Job started on: Mon Aug 31 10:02:29 CEST 2015 .
#  # FONLL heavy quark hadroproduction cross section
#  # FONLL version and perturbative order: ## FONLL v1.3.2 fonll [ds/dpt^2dy (pb/GeV^2)]
#  # quark = charm
#  # final state = meson (meson = D0). NP params (cm,lm,hm) = 0.1, 0.06, 0.135
#  # BR(q->meson) = 1
#  # ebeam1 = 4000, ebeam2 = 4000
#  # PDF set = CTEQ6.6
#  # ptmin = 0
#  # ptmax = 20
#  # ymin  = 2
#  # ymax  = 2.5
#  # cross section is ds/dpt (pb/GeV)
#  #  pt      central
#     0.0000 0.0000e+00
#     0.2500 3.4330e+07
#     0.5000 6.0499e+07
#     0.7500 8.2892e+07
#     1.0000 9.8492e+07
#     1.2500 1.0412e+08
#     1.5000 1.0060e+08
#     1.7500 9.1271e+07 '''
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @data 2015-08-31
def makeGraph2 ( text ) :
    """Create TGraph from simple two-column text format,
    e.g. the files provided by FONLL on-line calculator can be parsed
    @see http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html
    >>> graph = makeGraph2 ( '''
    ... # Job started on: Mon Aug 31 10:02:29 CEST 2015 .
    ... # FONLL heavy quark hadroproduction cross section
    ... # FONLL version and perturbative order: ## FONLL v1.3.2 fonll [ds/dpt^2dy (pb/GeV^2)]
    ... # quark = charm
    ... # final state = meson (meson = D0). NP params (cm,lm,hm) = 0.1, 0.06, 0.135
    ... # BR(q->meson) = 1
    ... # ebeam1 = 4000, ebeam2 = 4000
    ... # PDF set = CTEQ6.6
    ... # ptmin = 0
    ... # ptmax = 20
    ... # ymin  = 2
    ... # ymax  = 2.5
    ... # cross section is ds/dpt (pb/GeV)
    ... #  pt      central
    ... 0.0000 0.0000e+00
    ... 0.2500 3.4330e+07
    ... 0.5000 6.0499e+07
    ... 0.7500 8.2892e+07
    ... 1.0000 9.8492e+07
    ... 1.2500 1.0412e+08
    ... 1.5000 1.0060e+08
    ... 1.7500 9.1271e+07 '''
    """
    text = text.split('\n')
    if not text : return ROOT.TGraph ()

    vals = {}
    for line in text :
        line = line.strip()
        if not line  : continue
        if '#' == line[0] :
            logger.info   ('Comment: %s' % line )
            continue
        ## remove  traling comments 
        line, s, tail = line.partition('#')
        if not line : continue 
        values = line.split (' ')
        if 2  > len ( values ) :
            logger.warning('Strange line, skip it: "%s"' % line )
            continue
        x = values[0].strip() 
        y = values[1].strip() 
        x = float(x)
        y = float(y)
        vals[ x ] = y
    ## make a graph 
    return makeGraph ( vals )


# =============================================================================
## create TGraphAsymmErrors from the plain text multicolumn-column format with optional comments.
#  E.g. the files provided by FONLL on-line calculator can be parsed
#  @see http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html
#  @code
#  gr = makeGraph3 ( '''
#  # Job started on: Tue Sep  1 09:17:07 CEST 2015 .
#  # FONLL heavy quark hadroproduction cross section
#  # FONLL version and perturbative order: ## FONLL v1.3.2 fonll [ds/dpt^2dy (pb/GeV^2)]
#  # quark = charm
#  # final state = quark
#  # ebeam1 = 3500, ebeam2 = 3500
  # PDF set = CTEQ6.6
#  # ptmin = 5
#  # ptmax = 20
#  # ymin  = -1
#  # ymax  = 1
#  # Uncertainties from scales
#  # cross section is ds/dpt (pb/GeV)
#  #  pt      central      min       max       min_sc     max_sc
#    5.0000 6.6425e+07 3.6583e+07 1.1376e+08 3.6583e+07 1.1376e+08
#   10.0000 5.5774e+06 4.1829e+06 7.8655e+06 4.1829e+06 7.8655e+06
#   15.0000 9.7257e+05 7.7270e+05 1.2580e+06 7.7270e+05 1.2580e+06
#   20.0000 2.5771e+05 2.1161e+05 3.1773e+05 2.1161e+05 3.1773e+05
#  '''
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @data 2015-08-31
def makeGraph3 ( text ) :
    """Create TGraphAsymmErrors from simple multicolumn-column text format,
    e.g. the files provided by FONLL on-line calculator can be parsed
    @see http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html
    >>> graph = makeGraph3 ( '''
    ... # Job started on: Tue Sep  1 09:17:07 CEST 2015 .
    ... # FONLL heavy quark hadroproduction cross section
    ... # FONLL version and perturbative order: ## FONLL v1.3.2 fonll [ds/dpt^2dy (pb/GeV^2)]
    ... # quark = charm
    ... # final state = quark
    ... # ebeam1 = 3500, ebeam2 = 3500
    ... # PDF set = CTEQ6.6
    ... # ptmin = 5
    ... # ptmax = 20
    ... # ymin  = -1
    ... # ymax  = 1
    ... # Uncertainties from scales
    ... # cross section is ds/dpt (pb/GeV)
    ... #  pt      central      min       max       min_sc     max_sc
    ...  5.0000 6.6425e+07 3.6583e+07 1.1376e+08 3.6583e+07 1.1376e+08
    ... 10.0000 5.5774e+06 4.1829e+06 7.8655e+06 4.1829e+06 7.8655e+06
    ... 15.0000 9.7257e+05 7.7270e+05 1.2580e+06 7.7270e+05 1.2580e+06
    ... 20.0000 2.5771e+05 2.1161e+05 3.1773e+05 2.1161e+05 3.1773e+05'''
    """
    text = text.split('\n')
    if not text : return ROOT.TGraphAsymmErrors()

    points = [] 
    for line in text :
        line = line.strip()
        if not line  : continue
        if '#' == line[0] :
            logger.info   ('Comment: %s' % line )
            continue
        ## remove  traling comments 
        line, s, tail = line.partition('#')
        if not line : continue 
        values = line.split (' ')
        if 4  > len ( values ) :
            logger.warning('Strange line, skip it: "%s"' % line )
            continue
        x  = values[0].strip() 
        y  = values[1].strip() 
        mn = values[2].strip() 
        mx = values[3].strip() 
        x  = float(x)
        y  = float(y)
        mn = float(mn)
        mx = float(mx)
        points.append ( (x , 0 , 0 , y , y - mn , mx - y ) ) 
    #
    ## make a graph
    #
    graph = ROOT.TGraphAsymmErrors ( len ( points ) )
    for p in range( len ( points ) ) : graph[p] = points[p]
    return graph 

# =============================================================================
## create *three* TGraphAsymmErrors from the plain text multicolumn-column
#  format with optional comments.
#  E.g. the files provided by FONLL on-line calculator can be parsed
#  @see http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html
#  @code
#  gtot,gscale,gmass  = makeGraphs3 ( '''
#  # Job started on: Tue Sep  1 09:44:37 CEST 2015 .
#  # FONLL heavy quark hadroproduction cross section
#  # FONLL version and perturbative order: ## FONLL v1.3.2 fonll [ds/dpt^2dy (pb/GeV^2)]
#  # quark = bottom
#  # final state = quark
#  # ebeam1 = 3500, ebeam2 = 3500
#  # PDF set = CTEQ6.6
#  # ptmin = 5
#  # ptmax = 20
#  # ymin  = -1
#  # ymax  = 1
#  # Uncertainties from scales, masses combined quadratically
#  # cross section is ds/dpt (pb/GeV)
#  #  pt      central      min       max       min_sc     max_sc     min_mass   max_mass
#     5.0000 8.6807e+06 5.4820e+06 1.2523e+07 5.6934e+06 1.2285e+07 7.5373e+06 1.0012e+07
#    20.0000 2.2821e+05 1.8208e+05 2.9430e+05 1.8256e+05 2.9398e+05 2.2160e+05 2.3465e+05
#  '''
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @data 2015-08-31
def makeGraphs3 ( text ) :
    """Create *three  TGraphAsymmErrors from simple multicolumn-column text format,
    e.g. the files provided by FONLL on-line calculator can be parsed
    @see http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html
    >>> gtot,gscale,gmass = makeGraphs3 ( '''
    ... # Job started on: Tue Sep  1 09:44:37 CEST 2015 .
    ... # FONLL heavy quark hadroproduction cross section
    ... # FONLL version and perturbative order: ## FONLL v1.3.2 fonll [ds/dpt^2dy (pb/GeV^2)]
    ... # quark = bottom
    ... # final state = quark
    ... # ebeam1 = 3500, ebeam2 = 3500
    ... # PDF set = CTEQ6.6
    ... # ptmin = 5
    ... # ptmax = 20
    ... # ymin  = -1
    ... # ymax  = 1
    ... # Uncertainties from scales, masses combined quadratically
    ... # cross section is ds/dpt (pb/GeV)
    ... #  pt      central      min       max       min_sc     max_sc     min_mass   max_mass
    ...  5.0000 8.6807e+06 5.4820e+06 1.2523e+07 5.6934e+06 1.2285e+07 7.5373e+06 1.0012e+07
    ... 20.0000 2.2821e+05 1.8208e+05 2.9430e+05 1.8256e+05 2.9398e+05 2.2160e+05 2.3465e+05
    '''
    """
    text = text.split('\n')
    if not text : return ( ROOT.TGraphAsymmErrors () ,
                           ROOT.TGraphAsymmErrors () ,
                           ROOT.TGraphAsymmErrors () )
    
    ptot   = [] 
    pscale = [] 
    pmass  = [] 
    for line in text :
        line = line.strip()
        if not line  : continue
        if '#' == line[0] :
            logger.info   ('Comment: %s' % line )
            continue
        ## remove  traling comments 
        line, s, tail = line.partition('#')
        if not line : continue 
        values = line.split (' ')
        if 8  > len ( values ) :
            logger.warning('Strange line, skip it: "%s"' % line )
            continue
        x   = values[0].strip() 
        y   = values[1].strip() 
        mnt = values[2].strip() 
        mxt = values[3].strip()
        mnm = values[4].strip() 
        mxm = values[5].strip()
        mns = values[6].strip() 
        mxs = values[7].strip()
        
        x   = float ( x   ) 
        y   = float ( y   )
        mnt = float ( mnt )
        mxt = float ( mxt )
        mns = float ( mns )
        mxs = float ( mxs )
        mnm = float ( mnm )
        mxm = float ( mxm )
        
        ptot  .append ( ( x , 0 , 0 , y , y - mnt , mxt - y ) ) 
        pscale.append ( ( x , 0 , 0 , y , y - mns , mxs - y ) ) 
        pmass .append ( ( x , 0 , 0 , y , y - mnm , mxm - y ) ) 
    #
    ## make a graph
    #
    np = len( ptot ) 
    gtot   = ROOT.TGraphAsymmErrors ( np )
    gscale = ROOT.TGraphAsymmErrors ( np )
    gmass  = ROOT.TGraphAsymmErrors ( np )
    for p in range ( np ) :
        gtot   [p] = ptot   [p]
        gscale [p] = pscale [p]
        gmass  [p] = pmass  [p]
    return gtot,gscale,gmass


# =============================================================================
## create *four* TGraphAsymmErrors from the plain text multicolumn-column
#  format with optional comments.
#  E.g. the files provided by FONLL on-line calculator can be parsed
#  @see http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html
#  @code
#  gtot,gscale,gmass,gpdf  = makeGraphs4 ( '''
#  # Job started on: Tue Sep  1 09:53:00 CEST 2015 .
#  # FONLL heavy quark hadroproduction cross section
#  # FONLL version and perturbative order: ## FONLL v1.3.2 fonll [ds/dpt^2dy (pb/GeV^2)]
#  # quark = bottom
#  # final state = quark
#  # ebeam1 = 3500, ebeam2 = 3500
#  # PDF set = CTEQ6.6
#  # ptmin = 5
#  # ptmax = 20
#  # ymin  = -1
#  # ymax  = 1
#  # Uncertainties from scales, masses, PDFs combined quadratically
#  # cross section is ds/dpt (pb/GeV)
#  pt      central      min       max       min_sc     max_sc     min_mass   max_mass     min_pdf    max_pdf
#   5.0000 8.6807e+06 5.3324e+06 1.2649e+07 5.6934e+06 1.2285e+07 7.5373e+06 1.0012e+07 7.6909e+06 9.6704e+06
#  20.0000 2.2821e+05 1.7966e+05 2.9601e+05 1.8256e+05 2.9398e+05 2.2160e+05 2.3465e+05 2.1304e+05 2.4337e+05
#  '''
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @data 2015-08-31
def makeGraphs4 ( text ) :
    """Create *four*  TGraphAsymmErrors from simple multicolumn-column text format,
    e.g. the files provided by FONLL on-line calculator can be parsed
    @see http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html
    >>> gtot,gscale,gmass,gpdf = makeGraphs4 ( '''
    ... # Job started on: Tue Sep  1 09:53:00 CEST 2015 .
    ... # FONLL heavy quark hadroproduction cross section
    ... # FONLL version and perturbative order: ## FONLL v1.3.2 fonll [ds/dpt^2dy (pb/GeV^2)]
    ... # quark = bottom
    ... # final state = quark
    ... # ebeam1 = 3500, ebeam2 = 3500
    ... # PDF set = CTEQ6.6
    ... # ptmin = 5
    ... # ptmax = 20
    ... # ymin  = -1
    ... # ymax  = 1
    ... # Uncertainties from scales, masses, PDFs combined quadratically
    ... # cross section is ds/dpt (pb/GeV)
    ... pt      central      min       max       min_sc     max_sc     min_mass   max_mass     min_pdf    max_pdf
    ...  5.0000 8.6807e+06 5.3324e+06 1.2649e+07 5.6934e+06 1.2285e+07 7.5373e+06 1.0012e+07 7.6909e+06 9.6704e+06
    ... 20.0000 2.2821e+05 1.7966e+05 2.9601e+05 1.8256e+05 2.9398e+05 2.2160e+05 2.3465e+05 2.1304e+05 2.4337e+05
    ... '''
    """
    text = text.split('\n')
    if not text : return ( ROOT.TGraphAsymmErrors () ,
                           ROOT.TGraphAsymmErrors () ,
                           ROOT.TGraphAsymmErrors () ,
                           ROOT.TGraphAsymmErrors () )
    
    ptot   = [] 
    pscale = [] 
    pmass  = [] 
    ppdf   = [] 
    for line in text :
        line = line.strip()
        if not line  : continue
        if '#' == line[0] :
            logger.info   ('Comment: %s' % line )
            continue
        ## remove  traling comments 
        line, s, tail = line.partition('#')
        if not line : continue 
        values = line.split (' ')
        if 10  > len ( values ) :
            logger.warning('Strange line, skip it: "%s"' % line )
            continue
        x   = values[0].strip() 
        y   = values[1].strip() 
        mnt = values[2].strip() 
        mxt = values[3].strip()
        mnm = values[4].strip() 
        mxm = values[5].strip()
        mns = values[6].strip() 
        mxs = values[7].strip()
        mnp = values[8].strip() 
        mxp = values[9].strip()
        
        x   = float ( x   ) 
        y   = float ( y   )
        mnt = float ( mnt )
        mxt = float ( mxt )
        mns = float ( mns )
        mxs = float ( mxs )
        mnm = float ( mnm )
        mxm = float ( mxm )
        mnp = float ( mnp )
        mxp = float ( mxp )
        
        ptot  .append ( ( x , 0 , 0 , y , y - mnt , mxt - y ) )
        pscale.append ( ( x , 0 , 0 , y , y - mns , mxs - y ) ) 
        pmass .append ( ( x , 0 , 0 , y , y - mnm , mxm - y ) )
        ppdf  .append ( ( x , 0 , 0 , y , y - mnp , mxp - y ) )
    #
    ## make a graph
    #
    np     = len ( ptot ) 
    gtot   = ROOT.TGraphAsymmErrors ( np )
    gscale = ROOT.TGraphAsymmErrors ( np )
    gmass  = ROOT.TGraphAsymmErrors ( np )
    gpdf   = ROOT.TGraphAsymmErrors ( np )
    for p in range ( np ) :
        gtot   [p] = ptot   [p]
        gscale [p] = pscale [p]
        gmass  [p] = pmass  [p]
        gpdf   [p] = ppdf   [p]
    return gtot,gscale,gmass,gpdf 

# =============================================================================
## convert histogram to graph with optional  transformation
#  @code
#  histo = ...
#  graph = histo.toGraph()
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def hToGraph ( h1                   ,
               funcx = lambda s : s ,
               funcv = lambda s : s ) :
    """ Convert  1D-histogram into graph with optional transformations 
    >>> histo = ...
    >>> graph = histo.toGraph()
    """
    #
    ## book graph
    #
    graph = ROOT.TGraphErrors( len ( h1 )  - 2 )
    
    #
    ## copy attributes
    #
    copy_graph_attributes ( h1 , graph )

    for i in h1.items () :

        x = funcx  ( i[1] ) 
        v = funcv  ( i[2] )

        ## note the different convention 
        graph [ i[0] - 1 ] = x , v 
        
    return graph

# =============================================================================
## Helper function convert histogram to graph, applyong some
#  transformation 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _hToGraph_ ( h1 , funcx , funcy ) :
    """ Convert  1D-histogram into TGraphAsymmError
    applying sime in-flight transformations 
    """
    #
    ## book graph
    #
    graph = ROOT.TGraphAsymmErrors( len ( h1 )  - 2 )
    
    #
    ## copy attributes
    #
    copy_graph_attributes ( h1 , graph )

    for i in h1.items () :
        
        ip = i[0] - 1 ## different convention for TH1 and TGraph
        x  = i[1] 
        y  = i[2]

        x0 , xen , xep = funcx ( x , y )
        y0 , yen , yep = funcy ( x , y )

        graph.SetPoint      ( ip , x0  , y0  ) 
        graph.SetPointError ( ip , xen , xep , yen , yep ) 

    return graph

# =============================================================================
## convert histogram to graph
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def hToGraph2 ( h1 , bias ) :
    """Convert  1D-histogram into graph with small shift in x
    Useful for overlay of very similar plots
    >>> h1 = ....
    >>> g2 = h1.asGraph2 ( 0.1 ) ## shift for 10% of bin width    
    """
    if abs ( bias ) > 1 :
        raise ValueError ( ' Illegal value for "bias" parameter ')
    
    funcx = lambda x,y : ( x.value() + x.error()*bias , x.error()*(1+bias) , x.error()*(1-bias) ) 
    funcy = lambda x,y : ( y.value()                  , y.error()          , y.error()          ) 
        
    return _hToGraph_ ( h1 , funcx , funcy ) 

# =============================================================================
## convert histogram to graph
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def hToGraph3 ( h1 , bias ) :
    """ Convert  1D-histogram into graph with small shift in x
    Useful for overlay of very similar plots
    >>> h1 = ....
    >>> g2 = h1.asGraph2 ( 0.1 ) ## shift for 0.1 (absolute)    
    """
    for p in h1.items() :
        x = p[1]
        if x.error() < abs ( bias ) :
            raise ValueError ( ' Illegal value for "bias" parameter ')
        
    funcx = lambda x,y : ( x.value() + bias , x.error()+bias , x.error()-bias )
    funcy = lambda x,y : ( y.value()        , y.error()      , y.error()      ) 
        
    return _hToGraph_ ( h1 , funcx , funcy ) 


ROOT.TH1F.asGraph  = hToGraph
ROOT.TH1D.asGraph  = hToGraph
ROOT.TH1F.toGraph  = hToGraph
ROOT.TH1D.toGraph  = hToGraph
ROOT.TH1F.asGraph2 = hToGraph2
ROOT.TH1D.asGraph2 = hToGraph2
ROOT.TH1F.toGraph2 = hToGraph2
ROOT.TH1D.toGraph2 = hToGraph2
ROOT.TH1F.asGraph3 = hToGraph3
ROOT.TH1D.asGraph3 = hToGraph3
ROOT.TH1F.toGraph3 = hToGraph3
ROOT.TH1D.toGraph3 = hToGraph3

# =============================================================================
##  get point
#   @code
#   graph = ...
#   x, y = graph.point ( 3  ) 
#   @endcode
def _gr_point_ ( graph , point ) :
    """Get the point fmorm the graph :
    >>> graph = ...
    >>> x, y = graph.point ( 3  )
    >>> x, y = graph.get_point ( 3  ) ## ditto 
    """
    if isinstance ( point , integer_types ) and 0 <= point < len ( graph ) :
        
        x = ctypes.c_double ( 0.0 )
        y = ctypes.c_double ( 1.0 )
        
        graph.GetPoint ( point , x , y )

        return float ( x.value ) , float ( y.value )
    
    raise IndexError ( "Invalid index %s" % p ) 
    

ROOT.TGraph.    point = _gr_point_
ROOT.TGraph.get_point = _gr_point_

# =============================================================================
## use graph as function 
#  @code
#  graph = ...
#  y     = graph ( 0.2 )
#  @endcode
#  @see TGraph.
#  @see TGraph::Eval
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gr_call_ ( graph , x , spline = None , opts = '' ) :
    """ Use graph as function
    >>> graph = ...
    >>> y     = graph ( 0.2 ) 
    """
    if spline is None : spline = ROOT.nullptr 
    return graph.Eval ( float( x ) , spline , opts )

# =============================================================================
## Calculate an integral over the range \f$x_{low} \le x \le x_{high}\f$
#  It is not very efficient, but OK 
#  @code
#  graph = ...
#  i     = graph.integral ( 0 , 1 ) 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-08-31
def _gr_integral_ ( graph , xlow , xhigh , numerical = True ) :
    """Calculate an integral over the range \f$x_{low} \le x \le x_{high}\f$
    It is not very efficient, but OK 
    >>> graph = ...
    >>> i     = graph.integral ( 0 , 1 )
    """
    if numerical :
        from ostap.math.integral import integral 
        return integral ( graph , xlow , xhigh )
    
    tf1 = graph.asTF1()
    return tf1.Integral( xlow , xhigh ) 
        
# =============================================================================
## iterate over points in TGraphErrors
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gr_iter_ ( graph ) :
    """Iterate over graph points 
    >>> gr = ...
    >>> for i in gr : ...    
    """
    for ip in range ( 0 , len ( graph ) ) :
        yield ip
        
# =============================================================================
## get the point in TGraph
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gr_getitem_ ( graph , ipoint )  :
    """Get the point from the Graph
    >>> graph = ...
    >>> x,y   = graph[3]
    """
    
    if isinstance ( ipoint , slice ) :
        if ipoint.step and 1 != ipoint.step :
            raise IndexError ("Step must be +1!")
        return _gr0_getslice_ ( graph , ipoint.start , ipoint.stop ) 
    
    if ipoint < 0 : ipoint += len(graph) 
    if not ipoint in graph : raise IndexError 
    #
    return graph.point ( ipoint )

# =============================================================================
## set the point in TGraph
#  @code
#  graph = ...
#  x,y   = graph[1]
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gr_setitem_ ( graph , ipoint , point )  :
    """ Set the point from the Graph
    >>> graph    = ...
    >>> graph[1] = x,y 
    """
    #
    if not ipoint in graph : raise IndexError 
    #
    
    x = float ( point[0] )
    y = float ( point[1] )
    
    graph.SetPoint      ( ipoint , x , y )


# ==============================================================================
## remove the point from the graph
#  @code
#  graph = ...
#  del graph[1] 
#  @endcode 
def _gr_delitem_ ( graph, ipoint ) :
    """Remove the point from the graph
    >>> graph = ...
    >>> del graph[1] 
    >>> del graph[1:-1:2] 
    """
    #
    if isinstance ( ipoint , slice ) :
        old_len   = len ( graph ) 
        istart , istop , istep = ipoint.indices ( old_len  )
        idel = 0
        np   = 0
        for i in  range ( istart , istop , istep ) :
            np += 1 
            d = graph.RemovePoint ( i - idel )
            if 0 <= d : idel += 1
        new_len = len ( graph ) 
        return newlen + np == old_len
    elif not ipoint in graph : raise IndexError 
    #
    d = graph.RemovePoint ( ipoint )
    return  0 <= d

# =============================================================================
## iterate over the points in TGraph
#  @code 
#  gr = ...
#  for i,x,v in gr.    items(): ...
#  for i,x,v in gr.iteritems(): ... ##  ditto 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gr_iteritems_ ( graph ) :
    """Iterate over graph points 
    >>> graph = ...
    >>> for i,x,v in graph.    items(): ...
    >>> for i,x,v in graph.iteritems(): ... ##  ditto
    """
    for ip in graph :
        x , y = graph[ ip ] 
        yield ip , x , y 
        
        
# =============================================================================
## get the point in TGraphErrors
#  @code
#  graph = ...
#  x,y   = graph[4] 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gre_getitem_ ( graph , ipoint )  :
    """Get the point from the Graph
    >>> graph = ...
    >>> x,y = graph[3]
    """
    
    if isinstance ( ipoint , slice ) :
        if ipoint.step and 1 != ipoint.step :
            raise IndexError ("Step must be +1!")
        return _gr1_getslice_ ( graph , ipoint.start , ipoint.stop ) 

    if ipoint < 0 : ipoint += len(graph) 
    if not ipoint in graph : raise IndexError 
    #

    x_, v_ = graph.get_point ( ipoint ) 
    
    x = VE ( x_ , graph.GetErrorX ( ipoint )**2 )
    v = VE ( v_ , graph.GetErrorY ( ipoint )**2 )
    
    return x , v

# =============================================================================
## set the point in TGraphErrors
#  @code
#  graph    = ...
#  graph[4] = x,y 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gre_setitem_ ( graph , ipoint , point )  :
    """ Set the point in Graph
    >>> graph    = ...
    >>> graph[4] = x,y 
    """
    #
    if not ipoint in graph    : raise IndexError
    if not 2 == len ( point ) :
        raise AttributeError("Invalid dimension of 'point'")
    
    x = VE ( point[0] ) 
    v = VE ( point[1] ) 

    graph.SetPoint      ( ipoint , x . value () , v . value () )
    graph.SetPointError ( ipoint , x . error () , v . error () )


# =============================================================================
## represent TGraph as TF1
#  @code
#  graph = ...
#  fun   = grap.asTF1() 
#  @endcode
def _gr_as_TF1_ ( graph , interpolate = 3  ) :
    """ Represent TGraph as TF1
    >>> graph = ...
    >>> fun   = grap.asTF1() 
    """
    from ostap.fitting.param import _h1_as_fun_ 
    return _h1_as_fun_ ( graph , lambda x : x , interpolate = interpolate ) 

# =============================================================================
## iterate over points in TGraphErrors
#  @code
#  gre = ...
#  for i,x,v in gre.    items(): ...
#  for i,x,v in gre.iteritems(): ... ## ditto
#  @endcode
#  @see TGraphErrors
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gre_iteritems_ ( graph ) :
    """ Iterate over graph points 
    >>> gre = ...
    >>> for i,x,v in gre.    items(): ...
    >>> for i,x,v in gre.iteritems(): ... ## ditto
    """
    for ip in graph :        
        x , y = graph[ip]
        yield ip , x , y 

# =============================================================================
## iterate over points in TGraphAsymmErrors
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _grae_getitem_ ( graph , ipoint ) :
    """ Get the point from the graph 
    >>> grae = ...
    >>> x,xl,xh,y,yl,yh = grae[ 1 ]
    """
    
    if isinstance ( ipoint , slice ) :
        if ipoint.step and 1 != ipoint.step :
            raise IndexError ("Step must be +1!")
        return _gr2_getslice_ ( graph , ipoint.start , ipoint.stop ) 
    

    
    if ipoint < 0 : ipoint += len(graph) 
    if not ipoint in graph : raise IndexError 
    #
    
    x_ , v_ = graph.get_point ( ipoint )
    
    exl = graph.GetErrorXlow  ( ipoint )
    exh = graph.GetErrorXhigh ( ipoint )
    eyl = graph.GetErrorYlow  ( ipoint )
    eyh = graph.GetErrorYhigh ( ipoint )
    
    return float( x_ ) , -exl , exh , float( v_ ) , -eyl , eyh 

# =============================================================================
## iterate over points in TGraphAsymmErrors
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _grae_setitem_ ( graph , ipoint , point ) :
    """Set graph point
    
    >>> grae = ...
    >>> grae[1] = x,xl,xh,y,yl,yh
    
    """
    if not ipoint in graph : raise IndexError
    if 6 != len(point)     : raise AttributeError("Invalid lenght of 'point'")
    # 
    x   =       point[0]
    exl = abs ( point[1] )  ## allow them to be negative, to improve input format 
    exh =       point[2]
    y   =       point[3]
    eyl = abs ( point[4] )  ## allow them to be  negative to improve input format
    eyh =       point[5] 
    
    graph.SetPoint      ( ipoint , x   , y )
    graph.SetPointError ( ipoint , exl , exh , eyl , eyh )


# =============================================================================
## iterate over points in TGraphAsymmErrors
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _grae_iteritems_ ( graph ) :
    """Iterate over graph points 
    
    >>> grae = ...
    >>> for i,x,xl,xh,y,yl,yh in grae.    items(): ...
    >>> for i,x,xl,xh,y,yl,yh in grae.iteritems(): ... ##   ditto
    
    """
    for ip in graph :
        vars = graph[ip]        
        yield (ip,) + vars

# =============================================================================
## get minimal-x 
#  @code
#  graph = ...
#  xmin  = graph.xmin () 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gr_xmin_ ( graph ) :
    """Get minimal x for the points
    >>> graph = ...
    >>> xmin  = graph.xmin () 
    """
    return graph.bb() [ 0 ] 

# =============================================================================
## get maximal-x 
#  @code
#  graph = ...
#  xmax  = graph.xmax () 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gr_xmax_ ( graph ) :
    """Get maximal x for the points
    >>> graph = ...
    >>> xmax  = graph.xmax () 
    """    
    return graph.bb() [ 1 ] 

# =============================================================================
## get minimal-y
#  @code
#  graph = ...
#  ymin  = graph.ymin () 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gr_ymin_ ( graph ) :
    """ Get minimal y for the points
    >>> graph = ...
    >>> ymin  = graph.ymin () 
    """
    return graph.bb() [ 2 ] 

# =============================================================================
## get maximal-y
#  @code
#  graph = ...
#  ymax  = graph.ymax () 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gr_ymax_ ( graph ) :
    """ Get maximal x for the points
    >>> graph = ...
    >>> ymax  = graph.ymax () 
    """    
    return graph.bb() [ 3 ] 

# =============================================================================
## get minimal and maximal x for the points
#  @code
#  graph = ...
#  xmin,xmax = graph.xminmax() 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gr_xminmax_ ( graph ) :
    """Get minimal and maximal x for the points
    >>> graph = ...
    >>> xmin,xmax = graph.xminmax() 
    """
    xmin , xmax ,  _ , _ = graph.bb()
    return xmin , xmax

# =============================================================================
## get minimal and maximal value for the points
#  @code
#  graph = ...
#  mn,mx = graph.yminmax() 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gr_yminmax_ ( graph ) :
    """
    Get minimal and maximal  for the points
    >>> graph = ...
    >>> mn,mx = graph.yminmax() 
    """
    _ , _ , ymin , ymax = graph.bb()
    return ymin , ymax


# =============================================================================
## Get a "bounding box" for the graph
#  @code
#  xmin, xmax , ymin , ymax = graph.bb() 
#  @endcode 
def _gr_bb_ ( graph , more = 0.0 ) :
    """ Get a ``bounding box'' for the graph
    >>> xmin, xmax , ymin , ymax = graph.bb() 
    """
    xmin = pos_infinity
    xmax = neg_infinity
    ymin = pos_infinity
    ymax = neg_infinity

    for i , x , y in graph.iteritems() :
        
        xmin = min ( xmin , x )  
        xmax = max ( xmax , x )
        ymin = min ( ymin , y )  
        ymax = max ( ymax , y )

    if more : 
        xmn  = min ( xmin , xmin - more * ( xmax - xmin ) )
        xmx  = max ( xmax , xmax + more * ( xmax - xmin ) )
        ymn  = min ( ymin , ymin - more * ( ymax - ymin ) )
        ymx  = max ( ymax , ymax + more * ( ymax - ymin ) )
        xmin = xmn
        xmax = xmx
        ymin = ymn
        ymax = ymx
        
    return xmin , xmax , ymin , ymax

# =============================================================================
## Get a "bounding box" for the graph
#  @code
#  xmin, xmax , ymin , ymax = graph.bb() 
#  @endcode 
def _gre_bb_ ( graph , more = 0.0 ) :
    """ Get a ``bounding box'' for the graph
    >>> xmin, xmax , ymin , ymax = graph.bb() 
    """
    xmin = pos_infinity
    xmax = neg_infinity
    ymin = pos_infinity
    ymax = neg_infinity

    for i , x , y in graph.iteritems() :
        
        xv = x.value ()
        ex = x.error ()
        
        yv = y.value ()
        ey = y.error ()

        xmin = min ( xmin , xv , xv + ex , xv - ex )  
        xmax = max ( xmax , xv , xv + ex , xv - ex )
        ymin = min ( ymin , yv , yv + ey , yv - ey )  
        ymax = max ( ymax , yv , yv + ey , yv - ey )
        
    if more :
        
        xmn  = min ( xmin , xmin - more * ( xmax - xmin ) )
        xmx  = max ( xmax , xmax + more * ( xmax - xmin ) )
        ymn  = min ( ymin , ymin - more * ( ymax - ymin ) )
        ymx  = max ( ymax , ymax + more * ( ymax - ymin ) )
        xmin = xmn
        xmax = xmx
        ymin = ymn
        ymax = ymx

    return xmin , xmax , ymin , ymax

# =============================================================================
## Get a "bounding box" for the graph
#  @code
#  xmin, xmax , ymin , ymax = graph.bb() 
#  @endcode 
def _grae_bb_ ( graph , more = 0.0 ) :
    """ Get a ``bounding box'' for the graph
    >>> xmin, xmax , ymin , ymax = graph.bb() 
    """
    xmin = pos_infinity
    xmax = neg_infinity
    ymin = pos_infinity
    ymax = neg_infinity

    for i , xv, enx , epx , yv , eny , epy in graph.iteritems() :

        xmin = min ( xmin , xv , xv + abs ( epx ) , xv - abs ( enx ) )
        xmax = max ( xmax , xv , xv + abs ( epx ) , xv - abs ( enx ) )
        ymin = min ( ymin , yv , yv + abs ( epy ) , xv - abs ( eny ) )
        ymax = max ( ymax , yv , yv + abs ( epy ) , xv - abs ( eny ) )

    if more :
        
        xmn  = min ( xmin , xmin - more * ( xmax - xmin ) )
        xmx  = max ( xmax , xmax + more * ( xmax - xmin ) )
        ymn  = min ( ymin , ymin - more * ( ymax - ymin ) )
        ymx  = max ( ymax , ymax + more * ( ymax - ymin ) )
        xmin = xmn
        xmax = xmx
        ymin = ymn
        ymax = ymx

    return xmin , xmax , ymin , ymax

# =============================================================================
## Get a "bounding box" for the graph
#  @code
#  xmin, xmax , ymin , ymax = graph.bb() 
#  @endcode 
def _mg_bb_ ( graph , more = 0.0 ) :
    """ Get a ``bounding box'' for the graph
    >>> xmin, xmax , ymin , ymax = graph.bb() 
    """
    xmin = pos_infinity
    xmax = neg_infinity
    ymin = pos_infinity
    ymax = neg_infinity

    _gs = graph.GetListOfGraps()
    for gr in _gs :
        
        xmn , xmx , ymn , ymx = gr.bb ( more )
        
        xmin = min ( xmin , xmn )
        xmax = max ( xmax , xmx )
        ymin = min ( ymin , ymn )
        ymax = max ( ymax , ymx )

    return xmin , xmax , ymin , ymax


ROOT.TGraph.bb            =   _gr_bb_ 
ROOT.TGraphErrors.bb      =  _gre_bb_ 
ROOT.TGraphAsymmErrors.bb = _grae_bb_ 
ROOT.TMultiGraph.bb       =   _mg_bb_ 
    
# =============================================================================
## get "slice" for graph 
#  @code     
#    >>> graph = ...
#    >>> gr1   = graph[2:10] 
#  @endcode     
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-03-28 
def _gr0_getslice_ ( graph , i , j ) :
    """Get the ``slice'' for TGraph:
    >>> graph = ...
    >>> gr1   = graph[2:10]
    """
    np = len ( graph ) 

    if i is None : i = 0
    if j is None : j = np 
    
    while i < 0 : i += np
    while j < 0 : j += np

    new_graph = ROOT.TGraph( j - i ) if  i < j else  ROOT.TGraph()
    copy_graph_attributes ( graph , new_graph )

    ii = 0 
    while i < j :
        new_graph[ ii ] = graph[i]
        ii += 1 
        i  += 1 
        
    return new_graph 

# =============================================================================
## get "slice" for graph 
#  @code     
#    >>> graph = ...
#    >>> gr1   = graph[2:10] 
#  @endcode     
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-03-28 
def _gr1_getslice_ ( graph , i , j ) :
    """Get the ``slice'' for TGraphErrors:
    >>> graph = ...
    >>> gr1   = graph[2:10]
    """
    np = len ( graph ) 

    if i is None : i = 0
    if j is None : j = np 
    
    while i < 0 : i += np
    while j < 0 : j += np
    
    new_graph = ROOT.TGraphErrors( j - i ) if  i < j else  ROOT.TGraphErrors ()
    copy_graph_attributes ( graph , new_graph )
    
    ii = 0 
    while i < j :
        new_graph[ ii ] = graph[i]
        ii += 1
        i  += 1 
        
    return new_graph 


# =============================================================================
## get "slice" for graph 
#  @code     
#    >>> graph = ...
#    >>> gr1   = graph[2:10] 
#  @endcode     
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-03-28 
def _gr2_getslice_ ( graph , i , j ) :
    """Get the ``slice'' for TGraphAsymmErrors:
    >>> graph = ...
    >>> gr1   = graph[2:10]
    """
    np = len ( graph ) 
    
    if i is None : i = 0
    if j is None : j = np 
    
    while i < 0 : i += np
    while j < 0 : j += np
    
    new_graph = ROOT.TGraphAsymmErrors( j - i ) if  i < j else  ROOT.TGraphAsymmErrors ()
    copy_graph_attributes ( graph , new_graph )
    
    ii = 0 
    while i < j :
        new_graph[ ii ] = graph[i]
        ii += 1
        i  += 1 
        
    return new_graph 

# ============================================================================
## make sorted graph
#  @code
#  graph = ...
#  s     = graph.sorted() 
#  @endcode
#  @date   2016-03-28 
def _gr0_sorted_ ( graph , reverse = False ) :
    """Make sorted graph
    >>> graph = ...
    >>> s     = graph.sorted() 
    """
    
    oitems =        ( i for i in graph.items() ) 
    sitems = sorted ( oitems , key = lambda s :s[1] , reverse = reverse )
    
    new_graph = ROOT.TGraph ( len( graph ) )
    copy_graph_attributes ( graph , new_graph )

    ip = 0 
    for item in sitems :
        new_graph[ip] = item[1:]
        ip += 1

    return new_graph 

# ============================================================================
## make sorted graph
#  @code
#  graph = ...
#  s     = graph.sorted() 
#  @endcode
#  @date   2016-03-28 
def _gr1_sorted_ ( graph , reverse = False ) :
    """Make sorted graph
    >>> graph = ...
    >>> s     = graph.sorted() 
    """
    
    oitems =        ( i for i in graph.items() ) 
    sitems = sorted ( oitems , key = lambda s :s[1].value() , reverse = reverse )
    
    new_graph = ROOT.TGraphErrors ( len( graph ) )
    copy_graph_attributes ( graph , new_graph )

    ip = 0 
    for item in sitems :
        new_graph[ip] = item[1:]
        ip += 1

    return new_graph 

# ============================================================================
## make sorted graph
#  @code
#  graph = ...
#  s     = graph.sorted() 
#  @endcode
#  @date   2016-03-28 
def _gr2_sorted_ ( graph , reverse = False ) :
    """Make sorted graph
    >>> graph = ...
    >>> s     = graph.sorted() 
    """
    
    oitems =        ( i for i in graph.items() ) 
    sitems = sorted ( oitems , key = lambda s :s[1] , reverse = reverse )
    
    new_graph = ROOT.TGraphAsymmErrors ( len( graph ) )
    copy_graph_attributes ( graph , new_graph )

    ip = 0 
    for item in sitems :
        new_graph[ip] = item[1:]
        ip += 1

    return new_graph 

# =============================================================================
## remove points that do not satisfy the criteria
#  @code
#  graph = ...
#  graph.remove ( lambda s : s[0]<0.0 ) 
#  @endcode 
def _gr_remove_ ( graph , remove ) :
    """Remove points that do not satisfy the criteria
    >> graph = ...
    >>> graph.remove ( lambda s : s[0]<0.0 ) 
    """
    old_len = len ( graph ) 

    removed = [] 
    for point in graph :
        if remove ( *graph [ point ] ) :
            removed.append ( point )
        
    ir = 0 
    for i in  removed :
        d = graph.RemovePoint ( i - ir ) 
        if 0 <= d : ir += 1

    return len ( graph ) - old_len 

# =============================================================================
## create new graph, that contais only "good/filtered" points
#  @code
#  graph = ...
#  new_graph = graph.filter ( lambda s : s[0]<0.0 ) 
#  @endcode 
def _gr_filter_ ( graph , accept , name = '' ) :
    """ Create new graph, that contais only ``good/filtered'' points
    >>> graph = ...
    >>> new_graph = graph.filter ( lambda s : s[0]<0.0 ) 
    """

    if not name : name = graph.GetName() + '_filter'
    new_graph = graph.Clone ( name )
    copy_graph_attributes ( graph , new_graph )
    
    new_graph.remove ( lambda *s : not accept ( *s ) )
    
    return new_graph


# =============================================================================
## transform the graph
#  @code
#  nll = ....
#  fun = lambda x, y : math.exp ( -1 * y )  
#  lh  = nll.transform ( fun ) 
#  @endcode e
def _gr_transform_ ( graph , fun = lambda x , y : y ) :
    """Transform the graph
    
    >>> nll = ....
    >>> fun = lambda x, y : math.exp ( -1 * y )  
    >>> lh  = nll.transform ( fun ) 
    
    """
    
    if   isinstance ( graph , ROOT.TGraphAsymmErrors ) :
        raise TypeError("Transformation for ROOT.TGraphAsymmErrors is not defined ")
    
    if   isinstance ( graph , ROOT.TGraphErrors ) : new_graph =  ROOT.TGraphErrors ( graph )
    elif isinstance ( graph , ROOT.TGraph       ) : new_graph =  ROOT.TGraph       ( graph )
        
    ## make a copy 
    copy_graph_attributes ( graph , new_graph )
    
    for i in graph :
        x , y = graph [ i ]
        v     = fun ( x , y )        
        new_graph[i] = x , v
        
    return new_graph 


# =============================================================================
## transform the graph
#  @code
#  nll = ....
#  fun = lambda x, y : math.exp ( -1 * y )  
#  lh  = nll.transform ( fun ) 
#  @endcode e
def _grae_transform_ ( graph , fun = lambda x , y : y ) :
    """Transform the graph
    
    >>> nll = ....
    >>> fun = lambda x, y : math.exp ( -1 * y )  
    >>> lh  = nll.transform ( fun ) 
    
    """
    
    new_graph = ROOT.TGraphAsymmErrors ( graph ) 
    ## make a copy 
    copy_graph_attributes ( graph , new_graph )
    
    for i in graph :
        
        x , exl , exh , y ,  eyl , eyh  = graph [ i ]

        v  = float ( fun ( x , y               ) )
        v1 = float ( fun ( x , y + abs ( eyh ) ) )
        v2 = float ( fun ( x , y - abs ( eyl ) ) ) 

        vmin = min ( v , v1 , v2 )
        vmax = max ( v , v1 , v2 )
        evl  = vmin - v
        evh  = vmax - v
        
        new_graph[i] = x , exl , exh , v , evl , evh 
        
    return new_graph 

# =============================================================================
ROOT.TGraph       . __len__       = ROOT.TGraph . GetN 
ROOT.TGraph       . __contains__  = lambda s,i : i in range(0,len(s))
ROOT.TGraph       . __iter__      = _gr_iter_ 
ROOT.TGraph       . __call__      = _gr_call_

ROOT.TGraph       . xmin          = _gr_xmin_ 
ROOT.TGraph       . ymin          = _gr_ymin_ 
ROOT.TGraph       . xmax          = _gr_xmax_ 
ROOT.TGraph       . ymax          = _gr_ymax_ 

ROOT.TGraph       . xminmax       = _gr_xminmax_ 
ROOT.TGraph       . yminmax       = _gr_yminmax_ 
ROOT.TGraph       .  minmax       = _gr_yminmax_ 

ROOT.TGraph       . __getitem__   = _gr_getitem_ 
ROOT.TGraph       . __setitem__   = _gr_setitem_
ROOT.TGraph       . __delitem__   = _gr_delitem_
ROOT.TGraph       .     items     = _gr_iteritems_
ROOT.TGraph       . iteritems     = _gr_iteritems_

ROOT.TGraph       . transform     = _gr_transform_

ROOT.TGraphErrors . __getitem__   = _gre_getitem_ 
ROOT.TGraphErrors . __setitem__   = _gre_setitem_ 
ROOT.TGraphErrors .     items     = _gre_iteritems_ 
ROOT.TGraphErrors . iteritems     = _gre_iteritems_ 

ROOT.TH1F.asGraph = hToGraph
ROOT.TH1D.asGraph = hToGraph
ROOT.TH1F.toGraph = hToGraph
ROOT.TH1D.toGraph = hToGraph

ROOT.TGraphAsymmErrors.__len__       = ROOT.TGraphAsymmErrors . GetN 
ROOT.TGraphAsymmErrors.__contains__  = lambda s,i : i in range(0,len(s))
ROOT.TGraphAsymmErrors.__iter__      = _gr_iter_ 
ROOT.TGraphAsymmErrors.     items    = _grae_iteritems_
ROOT.TGraphAsymmErrors. iteritems    = _grae_iteritems_ 
ROOT.TGraphAsymmErrors.__getitem__   = _grae_getitem_ 
ROOT.TGraphAsymmErrors.__setitem__   = _grae_setitem_ 



ROOT.TGraphAsymmErrors . transform   = _grae_transform_

ROOT.TGraph       . integral         = _gr_integral_
ROOT.TGraph       . asTF1            = _gr_as_TF1_


ROOT.TGraph            .__getslice__  = _gr0_getslice_
ROOT.TGraphErrors      .__getslice__  = _gr1_getslice_
ROOT.TGraphAsymmErrors .__getslice__  = _gr2_getslice_

ROOT.TGraph            .sorted        = _gr0_sorted_
ROOT.TGraphErrors      .sorted        = _gr1_sorted_
ROOT.TGraphAsymmErrors .sorted        = _gr2_sorted_ 


ROOT.TGraph            .filter        = _gr_filter_
ROOT.TGraph            .remove        = _gr_remove_


# ==========================================================================
import ostap.math.math_ve as mve

ROOT.TGraph.__exp__    = lambda g : g.transform ( fun = lambda x, y : mve.exp    ( y ) )
ROOT.TGraph.__exp2__   = lambda g : g.transform ( fun = lambda x, y : mve.exp2   ( y ) )
ROOT.TGraph.__expm1__  = lambda g : g.transform ( fun = lambda x, y : mve.expm1  ( y ) )
ROOT.TGraph.__sqrt__   = lambda g : g.transform ( fun = lambda x, y : mve.sqrt   ( y ) )
ROOT.TGraph.__cbrt__   = lambda g : g.transform ( fun = lambda x, y : mve.cbrt   ( y ) )
ROOT.TGraph.__log__    = lambda g : g.transform ( fun = lambda x, y : mve.log    ( y ) )
ROOT.TGraph.__log2__   = lambda g : g.transform ( fun = lambda x, y : mve.log2   ( y ) )
ROOT.TGraph.__log10__  = lambda g : g.transform ( fun = lambda x, y : mve.log10  ( y ) )
ROOT.TGraph.__log1p__  = lambda g : g.transform ( fun = lambda x, y : mve.log1p  ( y ) )
ROOT.TGraph.__sin__    = lambda g : g.transform ( fun = lambda x, y : mve.sin    ( y ) )
ROOT.TGraph.__cos__    = lambda g : g.transform ( fun = lambda x, y : mve.cos    ( y ) )
ROOT.TGraph.__tan__    = lambda g : g.transform ( fun = lambda x, y : mve.tan    ( y ) )
ROOT.TGraph.__sinh__   = lambda g : g.transform ( fun = lambda x, y : mve.sinh   ( y ) )
ROOT.TGraph.__cosh__   = lambda g : g.transform ( fun = lambda x, y : mve.cosh   ( y ) )
ROOT.TGraph.__tanh__   = lambda g : g.transform ( fun = lambda x, y : mve.tanh   ( y ) )
ROOT.TGraph.__asin__   = lambda g : g.transform ( fun = lambda x, y : mve.asin   ( y ) )
ROOT.TGraph.__acos__   = lambda g : g.transform ( fun = lambda x, y : mve.acos   ( y ) )
ROOT.TGraph.__atan__   = lambda g : g.transform ( fun = lambda x, y : mve.atan   ( y ) )
ROOT.TGraph.__asinh__  = lambda g : g.transform ( fun = lambda x, y : mve.asinh  ( y ) )
ROOT.TGraph.__acosh__  = lambda g : g.transform ( fun = lambda x, y : mve.acosh  ( y ) )
ROOT.TGraph.__atanh__  = lambda g : g.transform ( fun = lambda x, y : mve.atanh  ( y ) )
ROOT.TGraph.__erf__    = lambda g : g.transform ( fun = lambda x, y : mve.erf    ( y ) )
ROOT.TGraph.__erfc__   = lambda g : g.transform ( fun = lambda x, y : mve.erfc   ( y ) )
ROOT.TGraph.__erfcx__  = lambda g : g.transform ( fun = lambda x, y : mve.erfcx  ( y ) )
ROOT.TGraph.__erfi__   = lambda g : g.transform ( fun = lambda x, y : mve.erfi   ( y ) )
ROOT.TGraph.__tgamma__ = lambda g : g.transform ( fun = lambda x, y : mve.tgamma ( y ) )
ROOT.TGraph.__lgamma__ = lambda g : g.transform ( fun = lambda x, y : mve.lgamma ( y ) )
ROOT.TGraph.__sech__   = lambda g : g.transform ( fun = lambda x, y : mve.sech   ( y ) )
ROOT.TGraph.__probit__ = lambda g : g.transform ( fun = lambda x, y : mve.probit ( y ) )
ROOT.TGraph.__pow__    = lambda g , *o : g.transform ( fun = lambda x, y : mve.pow ( y , *o ) )

# =============================================================================
## set color attributes  
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-01-21 
def _color_ ( self , color = 2 , marker = 20 , size = -1 ) :
    """Set color attributes

    >>> h.color ( 3 ) 
    """
    #
    if hasattr ( self , 'SetLineColor'   ) : self.SetLineColor   ( color  )
    if hasattr ( self , 'SetMarkerColor' ) : self.SetMarkerColor ( color  )
    if hasattr ( self , 'SetMarkerStyle' ) : self.SetMarkerStyle ( marker )
    ##
    if 0 > size and hasattr ( self , 'GetMarkerSize' ) and not marker in ( 1 , 6 , 7 ) :
        size = self.GetMarkerSize()
        if   marker in ( 22 , 23 , 29 , 33 ) : size *= 2.0 ## small filled objects
        elif marker in ( 24 , 25 , 28 , 4  ) : size *= 1.5 ## large open objects  
        elif marker in ( 26 , 27 , 30 , 32 ) : size *= 2.2 ## small open objects
    ##
    if 0 < size and hasattr ( self , 'SetMarkerSize' ) :
        self.SetMarkerSize ( size )
    #
    return self
# =============================================================================
## set color attributes  
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-01-21 
def _red_     ( self , marker   = 20 ) : return _color_( self , 2 , marker ) 
# =============================================================================
## set color attributes  
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-01-21 
def _blue_    ( self , marker   = 25 ) : return _color_( self , 4 , marker )
# =============================================================================
## set color attributes  
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-01-21 
def _magenta_ ( self , marker   = 22 ) : return _color_( self , 6 , marker ) 
# =============================================================================
## set color attributes  
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-01-21 
def _cyan_    ( self , marker   = 23 ) : return _color_( self , 7 , marker ) 
# =============================================================================
## set color attributes  
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-01-21 
def _green_   ( self , marker = 33 ) : return _color_( self , 8 , marker ) 
# =============================================================================
## set color attributes  
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-01-21 
def _yellow_  ( self , marker = 34 ) : return _color_( self , 92 , marker ) 


for _t in  ( ROOT.TH1D , ROOT.TH1F , ROOT.TGraph , ROOT.TGraphErrors ) :
    
    _t . color   = _color_
    _t . red     = _red_
    _t . blue    = _blue_
    _t . magenta = _magenta_
    _t . cyan    = _cyan_
    _t . green   = _green_
    _t . yellow  = _yellow_

# =============================================================================
## min-value for the graph)
def _gr_xmax_ ( graph ) :
    """Get x-max for the graph    
    >>> xmax = graph.xmax()
    """
    #
    _size = len ( graph )
    if 0 == _size : return 1
    #
    _last = _size - 1
    #
    x_ , y_ = graph.get_point ( _last )
    #
    return x_
# =============================================================================
## min-value for the graph)
def _gr_xmin_ ( graph ) :
    """Get x-min for the graph    
    >>> xmin = graph.xmin()    
    """
    #
    _size = len ( graph )
    if 0 == _size : return 0
    #
    x_ , y_ = graph.get_point ( 0 )
    #
    return x_

# =============================================================================
## minmax-value for the graph)
def _gr_xminmax_ ( graph ) :
    """Get x-minmax for the graph
    >>> xmin,xmax = graph.xminmax() 
    """
    #
    return  graph.xmin() , graph.xmax() 

ROOT.TGraph  . xmin    = _gr_xmin_
ROOT.TGraph  . xmax    = _gr_xmax_
ROOT.TGraph  . xminmax = _gr_xminmax_


# =============================================================================
## a bit of trivial math-like operations  with graphs 
# =============================================================================
## scale the graph
#  @code
#  gr = ...
#  ng = gr * 10 
#  @endcode
def _gr_mul_ ( graph , scale ) :
    """Scale the graph
    graph = ...
    newg  = graph * 10 
    """
    
    if not isinstance ( scale , num_types + (VE,) ) : return NotImplemented
    
    new_graph = ROOT.TGraph (  len ( graph ) )
    for i , x , y in graph.iteritems() :
        new_graph [ i ] = x , y * scale
    return new_graph

# ================================================================================
## scale the graph
#  @code
#  gr = ...
#  ng = gr / 10 
#  @endcode
def _gr_div_ ( graph , scale ) :
    """Scale the graph
    graph = ...
    newg  = graph / 10 
    """
    return graph * ( 1.0 / scale ) 

# =============================================================================
## scale the graph
#  @code
#  gr *= 10
#  @endcode
def _gr_imul_ ( graph , scale ) :
    """Scale the graph
    graph *= 10 
    """    
    if not isinstance ( scale , num_types +  ( VE , ) ) : return NotImplemented
    
    for i , x , y in graph.iteritems() :
        graph [ i ] = x , y * scale
    return graph


# ================================================================================
## scale the graph
#  @code
#  gr /= 10 
#  @endcode
def _gr_idiv_ ( graph , scale ) :
    """Scale the graph
    graph /= 10 
    """
    return _gr_imul_ ( graph ,  1.0 / scale )

# =================================================================================
## scale the graph
#  @code
#  gr = ...
#  ng =  10 / gr 
#  @endcode
def _gr_rdiv_ ( graph , scale ) :
    """Scale the graph
    graph = ...
    newg  = 10 / graph
    """
    if not isinstance ( scale , num_types + (VE,) ) : return NotImplemented
    
    new_graph = ROOT.TGraph (  len ( graph ) )
    for i , x , y in graph.iteritems() :
        new_graph [ i ] = x , 1.0 * scale / y 
    return new_graph

# =============================================================================
## shift the graph
#  @code
#  gr = ...
#  ng = gr + 10 
#  @endcode
def _gr_add_ ( graph , shift ) :
    """Shift the graph
    graph = ...
    newg  = graph + 10 
    """
    
    if not isinstance ( shift , num_types +  ( VE , ) ) : return NotImplemented
    
    new_graph = ROOT.TGraph (  len ( graph ) )
    for i , x , y in graph.iteritems() :
        new_graph [ i ] = x , y + shift
    return new_graph

# =============================================================================
## shift the graph
#  @code
#  gr = ...
#  ng = gr - 10 
#  @endcode
def _gr_sub_ ( graph , shift ) :
    """Shift the graph
    graph = ...
    newg  = graph - 10 
    """
    return _gr_add_ (  graph , -1.0 * shift )

# =============================================================================
## shift the graph
#  @code
#  gr = ...
#  ng = 10 - gr
#  @endcode
def _gr_rsub_ ( graph , shift ) :
    """Shift the graph
    graph = ...
    newg  = 10 - graph 
    """
    
    if not isinstance ( shift , num_types +  ( VE , ) ) : return NotImplemented
    
    new_graph = ROOT.TGraph (  len ( graph ) )
    for i , x , y in graph.iteritems() :
        new_graph [ i ] = x , shift - y
    return new_graph

# =============================================================================
## shift the graph
#  @code
#  gr += 10
#  @endcode
def _gr_iadd_ ( graph , shift ) :
    """Shift the graph
    graph += 10 
    """
    
    if not isinstance ( shift , num_types + ( VE , ) ) : return NotImplemented
    
    for i , x , y in graph.iteritems() :
        graph [ i ] = x , y + shift
    return graph

# =============================================================================
## shift the graph
#  @code
#  gr -= 10
#  @endcode
def _gr_isub_ ( graph , shift ) :
    """Shift the graph
    graph -= 10 
    """
    return _gr_iadd_ (  graph , -1.0 *  shift )

# ==============================================================================
## Left shift of the graph
#  @code
#  graph = ...
#  newg  = graph << 14.5 
#  @endcode 
def _gr_lshift_ ( graph , shift ) :
    """Left shift of the graph
    >>> graph = ...
    >>> newg  = graph << 14.5 
    """
    if not isinstance ( shift , num_types + ( VE , ) ) : return NotImplemented
    
    new_graph = ROOT.TGraph (  len ( graph ) )
    for i , x , y in graph.iteritems() :
        new_graph [ i ] = x + shift  , y
    return new_graph

# ==============================================================================
## Right shift of the graph
#  @code
#  graph = ...
#  newg  = graph >> 14.5 
#  @endcode 
def _gr_rshift_ ( graph , shift ) :
    """Right  shift of the graph
    >>> graph = ...
    >>> newg  = graph >> 14.5 
    """
    return _gr_lshift_ ( self , -1.0 * shift )


# ==============================================================================
## Left shift of the graph
#  @code
#  graph <<= 14.5 ...
#  @endcode 
def _gr_ilshift_ ( graph , shift ) :
    """Left shift of the graph
    >>> graph <<= 14.5 
    """
    if not isinstance ( shift , num_types + ( VE , ) ) : return NotImplemented
    
    for i , x , y in graph.iteritems() :
        graph [ i ] = x + shift  , y
    return graph
    
# ==============================================================================
## Right shift of the graph
#  @code
#  graph >>= 14.5
#  @endcode 
def _gr_irshift_ ( graph , shift ) :
    """Right  shift of the graph
    >>> graph >>= 14.5
    """
    return _gr_ilshift_ ( self , -1.0 * shift )


ROOT.TGraph. __mul__      = _gr_mul_
ROOT.TGraph.__rmul__      = _gr_mul_
ROOT.TGraph.__imul__      = _gr_imul_

ROOT.TGraph. __div__      = _gr_div_
ROOT.TGraph.__idiv__      = _gr_idiv_
ROOT.TGraph.__rdiv__      = _gr_rdiv_

ROOT.TGraph. __truediv__  = _gr_div_
ROOT.TGraph.__itruediv__  = _gr_idiv_
ROOT.TGraph.__rtruediv__  = _gr_rdiv_

ROOT.TGraph. __add__      = _gr_add_
ROOT.TGraph.__radd__      = _gr_add_
ROOT.TGraph.__iadd__      = _gr_iadd_
ROOT.TGraph. __sub__      = _gr_sub_
ROOT.TGraph.__rsub__      = _gr_rsub_
ROOT.TGraph.__isub__      = _gr_isub_

ROOT.TGraph.__lshift__    = _gr_lshift_
ROOT.TGraph.__rshift__    = _gr_rshift_
ROOT.TGraph.__ilshift__   = _gr_ilshift_
ROOT.TGraph.__irshift__   = _gr_irshift_


# =============================================================================
## scale the graph
#  @code
#  gr = ...
#  ng = gr * 10 
#  @endcode
def _gre_mul_ ( graph , scale ) :
    """Scale the graph
    graph = ...
    newg  = graph * 10 
    """
    
    if not isinstance ( scale , num_types + (VE,) ) : return NotImplemented
    
    new_graph = ROOT.TGraphErrors (  len ( graph ) )
    for i , x , y in graph.iteritems() :
        new_graph [ i ] = x , y * scale
    return new_graph

# ================================================================================
## scale the graph
#  @code
#  gr = ...
#  ng = gr / 10 
#  @endcode
def _gre_div_ ( graph , scale ) :
    """Scale the graph
    graph = ...
    newg  = graph / 10 
    """
    return graph * ( 1.0 / scale ) 

# =============================================================================
## scale the graph
#  @code
#  gr *= 10
#  @endcode
def _gre_imul_ ( graph , scale ) :
    """Scale the graph
    graph *= 10 
    """    
    if not isinstance ( scale , num_types +  ( VE , ) ) : return NotImplemented
    
    for i , x , y in graph.iteritems() :
        graph [ i ] = x , y * scale
    return graph


# ================================================================================
## scale the graph
#  @code
#  gr /= 10 
#  @endcode
def _gre_idiv_ ( graph , scale ) :
    """Scale the graph
    graph /= 10 
    """
    return _gre_imul_ ( graph ,  1.0 / scale )

# =================================================================================
## scale the graph
#  @code
#  gr = ...
#  ng =  10 / gr 
#  @endcode
def _gre_rdiv_ ( graph , scale ) :
    """Scale the graph
    graph = ...
    newg  = 10 / graph
    """
    if not isinstance ( scale , num_types + (VE,) ) : return NotImplemented
    
    new_graph = ROOT.TGraphErrors (  len ( graph ) )
    for i , x , y in graph.iteritems() :
        new_graph [ i ] = x , 1.0 * scale / y 
    return new_graph

# =============================================================================
## shift the graph
#  @code
#  gr = ...
#  ng = gr + 10 
#  @endcode
def _gre_add_ ( graph , shift ) :
    """Shift the graph
    graph = ...
    newg  = graph + 10 
    """
    
    if not isinstance ( shift , num_types +  ( VE , ) ) : return NotImplemented
    
    new_graph = ROOT.TGraphErrors (  len ( graph ) )
    for i , x , y in graph.iteritems() :
        new_graph [ i ] = x , y + shift
    return new_graph

# =============================================================================
## shift the graph
#  @code
#  gr = ...
#  ng = gr - 10 
#  @endcode
def _gre_sub_ ( graph , shift ) :
    """Shift the graph
    graph = ...
    newg  = graph - 10 
    """
    return _gre_add_ (  graph , -1.0 * shift )

# =============================================================================
## shift the graph
#  @code
#  gr = ...
#  ng = 10 - gr
#  @endcode
def _gre_rsub_ ( graph , shift ) :
    """Shift the graph
    graph = ...
    newg  = 10 - graph 
    """
    
    if not isinstance ( shift , num_types +  ( VE , ) ) : return NotImplemented
    
    new_graph = ROOT.TGraphErrors (  len ( graph ) )
    for i , x , y in graph.iteritems() :
        new_graph [ i ] = x , shift - y
    return new_graph

# =============================================================================
## shift the graph
#  @code
#  gr += 10
#  @endcode
def _gre_iadd_ ( graph , shift ) :
    """Shift the graph
    graph += 10 
    """
    
    if not isinstance ( shift , num_types + ( VE , ) ) : return NotImplemented
    
    for i , x , y in graph.iteritems() :
        graph [ i ] = x , y + shift
    return graph

# =============================================================================
## shift the graph
#  @code
#  gr -= 10
#  @endcode
def _gre_isub_ ( graph , shift ) :
    """Shift the graph
    graph -= 10 
    """
    return _gre_iadd_ (  graph , -1.0 *  shift )

# ==============================================================================
## Left shift of the graph
#  @code
#  graph = ...
#  newg  = graph << 14.5 
#  @endcode 
def _gre_lshift_ ( graph , shift ) :
    """Left shift of the graph
    >>> graph = ...
    >>> newg  = graph << 14.5 
    """
    if not isinstance ( shift , num_types + ( VE , ) ) : return NotImplemented
    
    new_graph = ROOT.TGraphErrors (  len ( graph ) )
    for i , x , y in graph.iteritems() :
        new_graph [ i ] = x + shift  , y
    return new_graph

# ==============================================================================
## Right shift of the graph
#  @code
#  graph = ...
#  newg  = graph >> 14.5 
#  @endcode 
def _gre_rshift_ ( graph , shift ) :
    """Right  shift of the graph
    >>> graph = ...
    >>> newg  = graph >> 14.5 
    """
    return _gre_lshift_ ( self , -1.0 * shift )


# ==============================================================================
## Left shift of the graph
#  @code
#  graph <<= 14.5 ...
#  @endcode 
def _gre_ilshift_ ( graph , shift ) :
    """Left shift of the graph
    >>> graph <<= 14.5 
    """
    if not isinstance ( shift , num_types + ( VE , ) ) : return NotImplemented
    
    for i , x , y in graph.iteritems() :
        graph [ i ] = x + shift  , y
    return graph
    
# ==============================================================================
## Right shift of the graph
#  @code
#  graph >>= 14.5
#  @endcode 
def _gre_irshift_ ( graph , shift ) :
    """Right  shift of the graph
    >>> graph >>= 14.5
    """
    return _gre_ilshift_ ( self , -1.0 * shift )


ROOT.TGraphErrors. __mul__      = _gre_mul_
ROOT.TGraphErrors.__rmul__      = _gre_mul_
ROOT.TGraphErrors.__imul__      = _gre_imul_

ROOT.TGraphErrors. __div__      = _gre_div_
ROOT.TGraphErrors.__idiv__      = _gre_idiv_
ROOT.TGraphErrors.__rdiv__      = _gre_rdiv_

ROOT.TGraphErrors. __truediv__  = _gre_div_
ROOT.TGraphErrors.__itruediv__  = _gre_idiv_
ROOT.TGraphErrors.__rtruediv__  = _gre_rdiv_

ROOT.TGraphErrors. __add__      = _gre_add_
ROOT.TGraphErrors.__radd__      = _gre_add_
ROOT.TGraphErrors.__iadd__      = _gre_iadd_
ROOT.TGraphErrors. __sub__      = _gre_sub_
ROOT.TGraphErrors.__rsub__      = _gre_rsub_
ROOT.TGraphErrors.__isub__      = _gre_isub_

ROOT.TGraphErrors.__lshift__    = _gre_lshift_
ROOT.TGraphErrors.__rshift__    = _gre_rshift_
ROOT.TGraphErrors.__ilshift__   = _gre_ilshift_
ROOT.TGraphErrors.__irshift__   = _gre_irshift_




# =============================================================================
## scale the graph
#  @code
#  gr = ...
#  ng = gr * 10 
#  @endcode
def _grae_mul_ ( graph , scale ) :
    """Scale the graph
    graph = ...
    newg  = graph * 10 
    """
    
    if not isinstance ( scale , num_types + (VE,) ) : return NotImplemented
    
    new_graph = ROOT.TGraphAsymmErrors (  len ( graph ) )
    for i , x , exl , exh , y , eyl , eyh in graph.iteritems() :
        new_graph [ i ] = x , exl , exh , y * scale , eyl * scale , eyh * scale 
    return new_graph

# ================================================================================
## scale the graph
#  @code
#  gr = ...
#  ng = gr / 10 
#  @endcode
def _grae_div_ ( graph , scale ) :
    """Scale the graph
    graph = ...
    newg  = graph / 10 
    """
    return graph * ( 1.0 / scale ) 

# =============================================================================
## scale the graph
#  @code
#  gr *= 10
#  @endcode
def _grae_imul_ ( graph , scale ) :
    """Scale the graph
    graph *= 10 
    """    
    if not isinstance ( scale , num_types +  ( VE , ) ) : return NotImplemented
    
    for i , x , exl , exh , y , eyl , eyh in graph.iteritems() :
        graph [ i ] = x , exl , exh , y * scale , eyl * scale , eyh * scale 
    return graph


# ================================================================================
## scale the graph
#  @code
#  gr /= 10 
#  @endcode
def _grae_idiv_ ( graph , scale ) :
    """Scale the graph
    graph /= 10 
    """
    return _grae_imul_ ( graph ,  1.0 / scale )

# =================================================================================
## scale the graph
#  @code
#  gr = ...
#  ng =  10 / gr 
#  @endcode
def _grae_rdiv_ ( graph , scale ) :
    """Scale the graph
    graph = ...
    newg  = 10 / graph
    """
    if not isinstance ( scale , num_types + (VE,) ) : return NotImplemented
    
    new_graph = ROOT.TGraphAsymmErrors (  len ( graph ) )
    for i , x , y in graph.iteritems() :
        new_graph [ i ] = x , 1.0 * scale / y 
    return new_graph

# =============================================================================
## shift the graph
#  @code
#  gr = ...
#  ng = gr + 10 
#  @endcode
def _grae_add_ ( graph , shift ) :
    """Shift the graph
    graph = ...
    newg  = graph + 10 
    """
    if not isinstance ( shift , num_types +  ( VE , ) ) : return NotImplemented
    
    new_graph = ROOT.TGraphAsymmErrors (  len ( graph ) )
    for i , x , exl , exh , y , eyl , eyh in graph.iteritems() :
        new_graph [ i ] = x , exl , exh , y + scale , eyl , eyh 
    return new_graph

# =============================================================================
## shift the graph
#  @code
#  gr = ...
#  ng = gr - 10 
#  @endcode
def _grae_sub_ ( graph , shift ) :
    """Shift the graph
    graph = ...
    newg  = graph - 10 
    """
    return _grae_add_ (  graph , -1.0 * shift )

# =============================================================================
## shift the graph
#  @code
#  gr = ...
#  ng = 10 - gr
#  @endcode
def _grae_rsub_ ( graph , shift ) :
    """Shift the graph
    graph = ...
    newg  = 10 - graph 
    """
    
    if not isinstance ( shift , num_types +  ( VE , ) ) : return NotImplemented
    
    new_graph = ROOT.TGraphAsymmErrors (  len ( graph ) )
    for i , x , exl , exh , y , eyl , eyh in graph.iteritems() :
        new_graph [ i ] = x , exl , exh , shift - y , -1.0 * abs ( eyh ) , abs ( eyl )  
    return new_graph

# =============================================================================
## shift the graph
#  @code
#  gr += 10
#  @endcode
def _grae_iadd_ ( graph , shift ) :
    """Shift the graph
    graph += 10 
    """
    
    if not isinstance ( shift , num_types + ( VE , ) ) : return NotImplemented
    
    for i , x , exl , exh , y , eyl , eyh in graph.iteritems() :
        graph [ i ] = x , exl , exh , y + scale , eyl , eyh 
    return graph

# =============================================================================
## shift the graph
#  @code
#  gr -= 10
#  @endcode
def _grae_isub_ ( graph , shift ) :
    """Shift the graph
    graph -= 10 
    """
    return _grae_iadd_ (  graph , -1.0 *  shift )

# ==============================================================================
## Left shift of the graph
#  @code
#  graph = ...
#  newg  = graph << 14.5 
#  @endcode 
def _grae_lshift_ ( graph , shift ) :
    """Left shift of the graph
    >>> graph = ...
    >>> newg  = graph << 14.5 
    """
    if not isinstance ( shift , num_types + ( VE , ) ) : return NotImplemented
    
    new_graph = ROOT.TGraphAsymmErrors (  len ( graph ) )
    for i , x , exl , exh , y , eyl , eyh in graph.iteritems() :
        new_graph [ i ] = x +  shift , exl , exh , y , eyl , eyh 
    return new_graph

# ==============================================================================
## Right shift of the graph
#  @code
#  graph = ...
#  newg  = graph >> 14.5 
#  @endcode 
def _grae_rshift_ ( graph , shift ) :
    """Right  shift of the graph
    >>> graph = ...
    >>> newg  = graph >> 14.5 
    """
    return _grae_lshift_ ( self , -1.0 * shift )

# ==============================================================================
## Left shift of the graph
#  @code
#  graph <<= 14.5 ...
#  @endcode 
def _grae_ilshift_ ( graph , shift ) :
    """Left shift of the graph
    >>> graph <<= 14.5 
    """
    if not isinstance ( shift , num_types + ( VE , ) ) : return NotImplemented
    
    for i , x , exl , exh , y , eyl , eyh in graph.iteritems() :
        graph [ i ] = x +  shift , exl , exh , y , eyl , eyh 
    return graph
    
# ==============================================================================
## Right shift of the graph
#  @code
#  graph >>= 14.5
#  @endcode 
def _grae_irshift_ ( graph , shift ) :
    """Right  shift of the graph
    >>> graph >>= 14.5
    """
    return _grae_ilshift_ ( self , -1.0 * shift )


ROOT.TGraphAsymmErrors. __mul__      = _grae_mul_
ROOT.TGraphAsymmErrors.__rmul__      = _grae_mul_
ROOT.TGraphAsymmErrors.__imul__      = _grae_imul_

ROOT.TGraphAsymmErrors. __div__      = _grae_div_
ROOT.TGraphAsymmErrors.__idiv__      = _grae_idiv_
ROOT.TGraphAsymmErrors.__rdiv__      = _grae_rdiv_

ROOT.TGraphAsymmErrors. __truediv__  = _grae_div_
ROOT.TGraphAsymmErrors.__itruediv__  = _grae_idiv_
ROOT.TGraphAsymmErrors.__rtruediv__  = _grae_rdiv_

ROOT.TGraphAsymmErrors. __add__      = _grae_add_
ROOT.TGraphAsymmErrors.__radd__      = _grae_add_
ROOT.TGraphAsymmErrors.__iadd__      = _grae_iadd_
ROOT.TGraphAsymmErrors. __sub__      = _grae_sub_
ROOT.TGraphAsymmErrors.__rsub__      = _grae_rsub_
ROOT.TGraphAsymmErrors.__isub__      = _grae_isub_

ROOT.TGraphAsymmErrors.__lshift__    = _grae_lshift_
ROOT.TGraphAsymmErrors.__rshift__    = _grae_rshift_
ROOT.TGraphAsymmErrors.__ilshift__   = _grae_ilshift_
ROOT.TGraphAsymmErrors.__irshift__   = _grae_irshift_


# =============================================================================
## append the graph with new point
#  @code
#  point = x , y
#  graph.append ( *point ) 
#  @endcode
def _gr_append_ ( graph , *point ) :
    """ append the graph with new point
    >>> point = x , y
    >>> graph.append ( *point ) 
    """
    last = len ( graph )
    graph.SetPoint ( last , 0 , 0 )
    graph [ last ] = point
    ##
    return len ( graph )

# =============================================================================
## pop the point fro mthe graph
#  @code
#  graph = ...
#  graph.pop ( 3 ) ## pop the point #3
#  graph.pop (   ) ## pop th elast point 
#  @endcode
def _gr_pop_  ( graph , i = None ) :
    """Pop the point fro mthe graph
    >>> graph = ...
    >>> graph.pop ( 3 ) ## pop the point #3
    >>> graph.pop (   ) ## pop th elast point 
    """

    if i is None :
        
        last  = len ( graph )
        if 1 <= last : 
            point = graph [ -1 ] 
            graph.RemovePoint ( last - 1 )
            return point
        
        return None
    
    if i < 0 : i += len ( graph )
    if not i in graph : raise IndexError ( "Point #%s is not in graph!" % i )

    point = graph [ i ]
    graph.RemovePoint ( i )
    return point

# =============================================================================
## swap two points in the graph
#  @code
#  graph = ...
#  graph.swap ( 1 , 6 ) 
#  @endcode 
def _gr_swap_ ( grap , i , j ) :

    if i < 0 : i += len ( graph )
    if not i in graph : raise IndexError ( "Point #%s is not in graph!" % i )
    
    if j < 0 : j += len ( graph )
    if not j in graph : raise IndexError ( "Point #%s is not in graph!" % j )

    pi = graph [ i ]
    pj = graph [ j ]
    
    graph [ i ] = pj 
    graph [ j ] = pi
    
    return graph

ROOT.TGraph.append    =  _gr_append_
ROOT.TGraph.swap      =  _gr_swap_ 
ROOT.TGraph.pop       =  _gr_pop_ 

# =============================================================================
## Transpose the graphs
# =============================================================================

# =============================================================================
## transpose the graph
#  @code
#  graph   = ...
#  graph_T = graph.transpose ()  
#  graph_T = graph.T() ## ditto 
#  @endcode
def _gr_transpose_ ( self ) :
    """Transpose the graph
    >>> graph   = ...
    >>> graph_T = graph.transpose ()  
    >>> graph_T = graph.T() ## ditto 
    """
    new_graph = ROOT.TGraph( len ( self ) )
    for i , x , y in self.iteritems() :
        new_graph[i] = y , x

    copy_graph_attributes ( self , new_graph ) 
    return new_graph 

# =============================================================================
## transpose the graph
#  @code
#  graph   = ...
#  graph_T = graph.transpose ()  
#  graph_T = graph.T() ## ditto 
#  @endcode
def _gre_transpose_ ( self ) :
    """Transpose the graph:
    >>> graph   = ...
    >>> graph_T = graph.transpose ()  
    >>> graph_T = graph.T() ## ditto 
    """
    new_graph = ROOT.TGraphErrors ( len ( self ) )
    for i , x , y in self.iteritems() :
        new_graph[i] = y , x
        
    copy_graph_attributes ( self , new_graph ) 
    return new_graph 

# =============================================================================
## transpose the graph
#  @code
#  graph   = ...
#  graph_T = graph.transpose ()  
#  graph_T = graph.T() ## ditto 
#  @endcode
def _grae_transpose_ ( self ) :
    """Transpose the graph:
    >>> graph   = ...
    >>> graph_T = graph.transpose ()  
    >>> graph_T = graph.T() ## ditto 
    """
    new_graph = ROOT.TGraphAsymmErrors ( len ( self ) )
    for item in self.iteritems() :
        ip, x , exl , exh , y , eyl , eyh =  item 
        new_graph [ ip ] = y , eyl , eyh , x , exl , exh
        
    copy_graph_attributes ( self , new_graph ) 
    return new_graph 


ROOT.TGraph.transpose            =   _gr_transpose_ 
ROOT.TGraph.T                    =   _gr_transpose_ 
ROOT.TGraphErrors.transpose      =  _gre_transpose_ 
ROOT.TGraphErrors.T              =  _gre_transpose_ 
ROOT.TGraphAsymmErrors.transpose = _grae_transpose_ 
ROOT.TGraphAsymmErrors.T         = _grae_transpose_ 


# ===============================================================================

# =============================================================================
## propagate the color for each graph in multigraph 
def _mg_color_ ( mgraph , color = 2 , marker = 20 , size = -1 ) :
    """ Propagate color settion to each graph in multigraph
    >>> mg = .. ## multigraph
    >>> mg.color ( 2 ) 
    """
    _graphs = mgraph.GetListOfGraps()
    for g in _graphs : g.color ( g , color , marker , size ) 

def _mg_red_      ( mgraph , marker = 20 , size = -1 ) :
    return _mg_color_ ( mgraph , 2 , marker , size )
def _mg_blue_     ( mgraph , marker = 25 , size = -1 ) :
    return _mg_color_ ( mgraph , 2 , marker , size )
def _mg_magenta_  ( mgraph , marker = 22 , size = -1 ) :
    return _mg_color_ ( mgraph , 2 , marker , size )
def _mg_cyan_     ( mgraph , marker = 23 , size = -1 ) :
    return _mg_color_ ( mgraph , 2 , marker , size )
def _mg_green_    ( mgraph , marker = 33 , size = -1 ) :
    return _mg_color_ ( mgraph , 2 , marker , size )
def _mg_yellow_   ( mgraph , marker = 34 , size = -1 ) :
    return _mg_color_ ( mgraph , 2 , marker , size )

ROOT.TMultiGraph.color    = _mg_color_
ROOT.TMultiGraph.red      = _mg_red_
ROOT.TMultiGraph.blue     = _mg_blue_
ROOT.TMultiGraph.magenta  = _mg_magenta_
ROOT.TMultiGraph.cyan     = _mg_cyan_
ROOT.TMultiGraph.green    = _mg_green_
ROOT.TMultiGraph.yellow   = _mg_yellow_


# =============================================================================
##  get graph  with certaoniono index
#   @code
#    mg = ....
#    g0 = mg[0] 
#   @endcode  
def _mg_getitem_ ( mgraph , item ) :
    """Get graph with certain index
    >>> mg = ....
    >>> g0 = mg[0] 
    """

    for i , g in enumerate  ( mgraph ) :
        if i == item  : return g
        
    raise IndexError("Invalid index for multigraph %s" % item )

ROOT.TMultiGraph.__getitem__ = _mg_getitem_             

# ==============================================================================
## length  of the multigraph
def _mg_len_  ( mgraph ) :
    """ Len of multigraph
    """
    _graphs = mgraph.GetListOfGraphs()
    return 0 if not _graphs else len ( _graphs ) 

ROOT.TMultiGraph.__len__ = _mg_len_             

# =============================================================================
## transpose the graph
#  @code
#  graph   = ...
#  graph_T = graph.transpose ()  
#  graph_T = graph.T() ## ditto 
#  @endcode
def _mg_transpose_ ( graph ) :
    """Transpose the graph:
    >>> graph   = ...
    >>> graph_T = graph.transpose ()  
    >>> graph_T = graph.T() ## ditto 
    """
    new_graph         = ROOT.TMultiGraph()
    new_graph._graphs = [] 
    
    _graphs = graph.GetListOfGraphs()
    for g in _graphs :
        tg =  g.T() 
        opt = graph.GetGraphDrawOption( g )
        new_graph.Add            ( tg , opt )
        new_graph._graphs.append ( tg )
        
    return new_graph 

ROOT.TMultiGraph.transpose = _mg_transpose_ 
ROOT.TMultiGraph.T         = _mg_transpose_ 


ROOT.TMultiGraph.xminmax = _gr_xminmax_
ROOT.TMultiGraph.yminmax = _gr_yminmax_
ROOT.TMultiGraph. minmax = _gr_yminmax_

ROOT.TMultiGraph.xmin    = _gr_xmin_ 
ROOT.TMultiGraph.ymin    = _gr_ymin_  
ROOT.TMultiGraph.xmax    = _gr_xmax_ 
ROOT.TMultiGraph.ymax    = _gr_ymax_ 


# =============================================================================
## Convert the histogram to into "Laffery-Wyatt" graph
#  See G.D. Lafferty and T.R. Wyatt,
#  ``Where to stick your data points: The treatment of measurements within wide bins,''
#  Nucl. Instrum. Meth. A355, 541 (1995).
#  @param histo  the histogram
#  @param func   the model
#  @attention: the model can be any reasonable model.
#  No need in multiplicative and additive terms:
#  the affine transformations do not affect the result.
#  e.g. following three graphs are equivalent: 
#  @code
#  >>> histo = ...
#  >>> gr1   = histo.lw_graph ( lambda x :     math.exp(-x)    )
#  >>> gr2   = histo.lw_graph ( lambda x : 100*math.exp(-x)    )
#  >>> gr3   = histo.lw_graph ( lambda x : 100*math.exp(-x)-10 )
#  @endcode 
#  If no reasonable model is known, the splines can be used instead: 
#  @code
#  >>> histo  = 
#  >>> spline = histo.(p,i,d)spline( .... )
#  >>> graph  = histo.lw_graph ( spline[2] ) 
#  @endcode 
#  @see https://doi.org/10.1016/0168-9002(94)01112-5
#  @author  Vanya BELYAEV  Ivan.Belyaev@itep.ru
#  @date    2014-12-08
def _lw_graph_ ( histo , func ) :
    """Convert the histogram to into ``Laffery-Wyatt'' graph
    See G.D. Lafferty and T.R. Wyatt,
    ``Where to stick your data points: The treatment of measurements within wide bins,''
    Nucl. Instrum. Meth. A355, 541 (1995).
    >>> histo = ... ## the histogram    
    ## the explicit model:
    >>> graph  = histo.lw_graph ( lambda x : math.exp(-x) )    
    ## use splines:
    >>> spline = histo.(p,i,d)spline( .... ) 
    >>> graph  = histo.lw_graph ( spline[2] )
    >>> histo.Draw('e1')
    >>> graph.Draw('e1p same')    
    """    
    #
    ## book graph
    #
    graph = ROOT.TGraphAsymmErrors( len ( histo )  - 2 )
        
    #
    ## copy attributes
    #
    copy_graph_attributes ( histo , graph )

    #
    ## start actual evaluations
    #
    
    from ostap.math.intergal   import integral as _integral 
    from ostap.math.rootfinder import findroot 
    
    for item in histo.items () :

        ibin  = item[0]
        x     = item[1]
        y     = item[2]

        yv    = y.value()
        ye    = y.error()
        
        xv    = x.value() 
        xe    = x.error()
        
        xmx   = xv + xe 
        xmn   = xv - xe 

        #
        ##  solve the equation f(x) = 1/dx*int(f,xmin,xmax)
        #

        ##  1) calculate int(f,xmin,xmax)
        fint  = _integral ( lambda x : float ( func ( x ) ) , xmn , xmx ) 

        ## bin-width 
        dx    = 2.0 * xe 
        
        fx    = float ( fint[0]/dx  ) 
        
        fxmin = float ( func ( xmn ) ) 
        fxmax = float ( func ( xmx ) ) 
        
        if 0 <= ( fxmin - fx ) * ( fxmax - fx ) : 
            logger.warning('Lafferty-Wyatt graph: invalid point: %s ' % x )
            r0 = x.value()
        else :
            ##  solve the equation f(x) - 1/dx*int(f,xmin,xmax) = 0 
            r0 = findroot ( lambda x :   ( float ( func ( x ) ) - fx ) ,
                            xmn               ,
                            xmx               ,
                            xtol = 0.005 * dx )
            
        ## fill graph
            
        ip   = ibin - 1 ## different conventions for TGraph and  TH1 

        xep  = xmx - r0
        xen  =       r0 - xmn
        
        graph.SetPoint      ( ip , r0  ,  yv           )
        graph.SetPointError ( ip , xen , xep , ye, ye  )

    return graph

ROOT.TH1D.lw_graph = _lw_graph_
ROOT.TH1F.lw_graph = _lw_graph_

# =============================================================================
## Convert the histogram to into "Laffery-Wyatt" graph
#  See G.D. Lafferty and T.R. Wyatt,
#  ``Where to stick your data points: The treatment of measurements within wide bins,''
#  Nucl. Instrum. Meth. A355, 541 (1995).
#  @param histo  the histogram
#  @param func   the model
#  @attention: the model can be any reasonable model.
#  No need in multiplicative and additive terms:
#  the affine transformations do not affect the result.
#  e.g. following three graphs are equivalent: 
#  @code
#  >>> histo = ...
#  >>> gr1   = lw_graph ( histo , lambda x :     math.exp(-x)    )
#  >>> gr2   = lw_graph ( histo , lambda x : 100*math.exp(-x)    )
#  >>> gr3   = lw_graph ( histo , lambda x : 100*math.exp(-x)-10 )
#  @endcode 
#  If no reasonable model is known, the splines can be used instead: 
#  @code
#  >>> histo  = 
#  >>> spline = histo.(p,i,d)spline( .... )
#  >>> graph  = lw_graph ( histo ,  spline[2] ) 
#  @endcode 
#  @see https://doi.org/10.1016/0168-9002(94)01112-5
#  @author  Vanya BELYAEV  Ivan.Belyaev@itep.ru
#  @date    2014-12-08
def lw_graph ( histo , func ) :
    """Convert the histogram to into ``Laffery-Wyatt'' graph
    See G.D. Lafferty and T.R. Wyatt,
    ``Where to stick your data points: The treatment of measurements within wide bins,''
    Nucl. Instrum. Meth. A355, 541 (1995).
    >>> histo = ... ## the histogram    
    ## the explicit model:
    >>> graph  = lw_graph ( histo , lambda x : math.exp(-x) )    
    ## use splines:
    >>> spline = histo.[p,i,d]spline( .... ) 
    >>> graph  = lw_graph ( histo , spline[2] )
    >>> histo.Draw('e1')
    >>> graph.Draw('e1p same')    
    """
    return _lw_graph_ ( histo , func ) 


# =============================================================================
## Create a graph, that represents the area between two curves/functions:
#  @code
#  import math 
#  graph = fill_area ( math.sin , math.cos , xmin = 0 , xmax = 5 )
#  graph.Draw('f')
#  @endcode
#  ``Functions'' could be
#  - plain functions
#  - function objects 
#  - histograms
#  - graphs
#  - ...
#  Inspired by Rene Brun's example
#  @see https://root.cern.ch/phpBB3/viewtopic.php?t=6346
def fill_area ( fun1                     ,
                fun2                     ,
                n         = 100          ,
                xmin      = neg_infinity ,
                xmax      = pos_infinity ,  
                log_scale = False        ) :
    """Create a graph, that represents the area between
    two curves/functions:
    >>> import math 
    >>> graph = fill_area ( math.sin , math.cos , xmin = 0 , xmax = 5 )
    >>> graph.Draw('f')
    ``Functions'' could be
    - plain functions
    - function objects 
    - histograms
    - graphs
    - ...
    Inspired by Rene Brun's example
    - see https://root.cern.ch/phpBB3/viewtopic.php?t=6346
    """

    #
    ## try to define proper x-range for graph..
    #  - from input arguments
    #  - from fun1 and fun2. 

    x1mn , x1mx = neg_infinity , pos_infinity 

    if hasattr   ( fun1 , 'xminmax' ) :
        x1mn,x1mx = fun1.xminmax()
    elif hasattr ( fun1 , 'xmin'    )    and hasattr ( fun1 , 'xmax' ) :
        x1mn,x1mx = fun1.xmin(), fun1.xmax()
    elif hasattr ( fun1 , 'GetXmin'    ) and hasattr ( fun1 , 'GetXmax' ) :
        x1mn,x1mx = fun1.GetXmin(), fun1.GetXmax()
    elif hasattr ( fun1 , 'GetXaxis' ) :
        axis = fun1.GetXaxis() 
        x1mn,x1mx = axis.GetXmin(), axis.GetXmax()

    x1mn = max ( x1mn , xmin )
    x1mx = min ( x1mx , xmax )        


    x2mn , x2mx = neg_infinity , pos_infinity 
    if hasattr   ( fun2 , 'xminmax' ) :
        x2mn,x2mx = fun2.xminmax()
    elif hasattr ( fun2 , 'xmin'    )    and hasattr ( fun2 , 'xmax' ) :
        x2mn,x2mx = fun2.xmin(), fun2.xmax()
    elif hasattr ( fun2 , 'GetXmin'    ) and hasattr ( fun2 , 'GetXmax' ) :
        x2mn,x2mx = fun2.GetXmin(), fun2.GetXmax()
    elif hasattr ( fun2 , 'GetXaxis' ) :
        axis = fun2.GetXaxis() 
        x2mn,x2mx = axis.GetXmin(), axis.GetXmax()

    x2mn = max ( x2mn , xmin )
    x2mx = min ( x2mx , xmax )        

    ## to be replaced with numpy.isfinite 
    if x1mn == neg_infinity and x2mn != neg_infinity : x1mn = x2mn 
    if x1mn != neg_infinity and x2mn == neg_infinity : x2mn = x1mn 
    if x1mx == pos_infinity and x2mn != pos_infinity : x1mx = x2mx 
    if x1mx != pos_infinity and x2mn == pos_infinity : x2mx = x1mx 

    if x1mn == neg_infinity or x2mn == neg_infinity or \
           x1mx == pos_infinity or x2mn == pos_infinity :
        raise ArrtibuteError("Can't define proper xmin/xmax")
            
    ## create the graph 
    graph = ROOT.TGraph ( 2 * n + 3 )

    if log_scale and 0 < x1mn < x1mx and 0 < x2mn < x2mx :

        dx1 = (x1mx/x1mn)**(1.0/n)
        dx2 = (x2mx/x2mn)**(1.0/n)
        
        x1f =  lambda i : x1mn*(dx1**i)
        x2f =  lambda i : x2mx/(dx2**i)

    else :
        
        dx1 = float(x1mx-x1mn)/n
        dx2 = float(x2mx-x2mn)/n
        
        x1f =  lambda i : x1mn + i * dx1 
        x2f =  lambda i : x2mx - i * dx2
        
    
    for i in range ( n + 1 ) :
        xi           = x1f  ( i ) 
        yi           = float( fun1( xi ) ) 
        graph[i]     = xi , yi 
        
    for i in range ( n + 1 ) :
        xi           = x2f  ( i ) 
        yi           = float( fun2 ( xi ) ) 
        graph[i+n+1] = xi , yi 

    ## the last point is the same as the first one 
    graph[-1] = graph[0] 
        
    graph.SetFillStyle(3013) 
    return graph 


# ==============================================================================
##  transpose the arrow
#   @code
#   a = ROOT.TArrow ( ... )
#   aT1 = a.transpose ()
#   aT2 = a.T()   ## ditto
#   @endcode 
def _ar_transpose_ (  arrow ) :
    """Transpose the arrow
    >>>  = ROOT.TArrow ( ... )
    >>> aT1 = a.transpose ()
    >>> aT2 = a.T()   ## ditto
    """
    na = ROOT.TArrow()
    arrow.Copy ( na )

    na.SetX1 ( arrow.GetY1 () )  
    na.SetY1 ( arrow.GetX1 () )  
    na.SetX2 ( arrow.GetY2 () )  
    na.SetY2 ( arrow.GetX2 () )  

    return na

ROOT.TArrow.transpose = _ar_transpose_
ROOT.TArrow.T         = _ar_transpose_

# ==============================================================================
##  transpose the box
#   @code
#   a = ROOT.TBox ( ... )
#   aT1 = a.transpose ()
#   aT2 = a.T()   ## ditto
#   @endcode 
def _box_transpose_ ( box  ) :
    """Transpose the box
    >>>  = ROOT.TBox ( ... )
    >>> aT1 = a.transpose ()
    >>> aT2 = a.T()   ## ditto
    """
    na = ROOT.TBox()
    box.Copy ( na )

    na.SetX1 ( box.GetY1 () )  
    na.SetY1 ( box.GetX1 () )  
    na.SetX2 ( box.GetY2 () )  
    na.SetY2 ( box.GetX2 () )  

    return na

ROOT.TBox.transpose = _box_transpose_
ROOT.TBox.T         = _box_transpose_

# ==============================================================================
##  transpose the line 
#   @code
#   a = ROOT.TLine ( ... )
#   aT1 = a.transpose ()
#   aT2 = a.T()   ## ditto
#   @endcode 
def _line_transpose_ ( line  ) :
    """Transpose the line
    >>>  = ROOT.TBox ( ... )
    >>> aT1 = a.transpose ()
    >>> aT2 = a.T()   ## ditto
    """
    na = ROOT.TLine()
    line.Copy ( na )

    na.SetX1 ( line.GetY1 () )  
    na.SetY1 ( line.GetX1 () )  
    na.SetX2 ( line.GetY2 () )  
    na.SetY2 ( line.GetX2 () )  

    return na

ROOT.TLine.transpose = _line_transpose_
ROOT.TLine.T         = _line_transpose_

# ==============================================================================
##  transpose the text
#   @code
#   a = ROOT.TText ( ... )
#   aT1 = a.transpose ()
#   aT2 = a.T()   ## ditto
#   @endcode 
def _text_transpose_ ( text ) :
    """Transpose the line
    >>>  = ROOT.TText ( ... )
    >>> aT1 = a.transpose ()
    >>> aT2 = a.T()   ## ditto
    """
    TEXT = type ( text ) 
    na   = TEXT ( text )
    ##
    na.SetX ( text.GetY () )
    na.SetY ( text.GetX () )
    na.SetTextAngle ( 90 - text.GetTextAngle () )
    
    ## al = text.GetTextAlign ()
    ## a1 , a2  = divmod ( al , 10 )
    ## na.SetTextAlign ( 10 * ( 4 - a1 ) + 4 - a2  )

    return na 
    
ROOT.TText.transpose = _text_transpose_
ROOT.TText.T         = _text_transpose_

# =============================================================================
_decorated_classes_ = (
    ROOT.TH1F              ,
    ROOT.TH1D              ,
    ROOT.TGraph            , 
    ROOT.TGraphErrors      ,
    ROOT.TGraphAsymmErrors ,
    ROOT.TArrow            ,  
    ROOT.TBox              ,
    ROOT.TLine             , 
    ROOT.TText
    )


_new_methods_      = (
    #
    ROOT.TGraph       . __len__       ,
    ROOT.TGraph       . __contains__  ,
    ROOT.TGraph       . __iter__      ,
    ROOT.TGraph       . __call__      , 
    #
    ROOT.TGraph       . xmin          ,
    ROOT.TGraph       . ymin          ,
    ROOT.TGraph       . xmax          ,
    ROOT.TGraph       . ymax          ,
    #
    ROOT.TMultiGraph  . xmin          , 
    ROOT.TMultiGraph  . ymin          ,
    ROOT.TMultiGraph  . xmax          , 
    ROOT.TMultiGraph  . ymax          ,
    #
    ROOT.TMultiGraph  . xminmax       , 
    ROOT.TMultiGraph  . yminmax       ,
    ROOT.TMultiGraph  .  minmax       ,
    #
    ROOT.TGraph       . xminmax       ,
    ROOT.TGraph       . yminmax       ,
    ROOT.TGraph       .  minmax       ,
    ROOT.TGraph       .  bb           ,
    #
    ROOT.TGraphErrors . xminmax       ,
    ROOT.TGraphErrors . yminmax       ,
    ROOT.TGraphErrors .  minmax       ,
    #
    ROOT.TGraphAsymmErrors . xminmax       ,
    ROOT.TGraphAsymmErrors . yminmax       ,
    ROOT.TGraphAsymmErrors .  minmax       ,
    #
    ROOT.TGraph       . __getitem__   ,
    ROOT.TGraph       . __setitem__   ,
    ROOT.TGraph       .     items     ,
    ROOT.TGraph       . iteritems     ,
    #
    ROOT.TGraph       .  __mul__      , 
    ROOT.TGraph       . __rmul__      , 
    ROOT.TGraph       . __imul__      ,
    #
    ROOT.TGraph       .  __div__      ,
    ROOT.TGraph       . __rdiv__      ,
    ROOT.TGraph       . __idiv__      ,
    # 
    ROOT.TGraph       .  __truediv__  ,
    ROOT.TGraph       . __itruediv__  ,
    ROOT.TGraph       . __rtruediv__  ,
    #
    ROOT.TGraph       .  __add__      ,
    ROOT.TGraph       . __radd__      ,
    ROOT.TGraph       . __iadd__      ,
    ROOT.TGraph       .  __sub__      ,
    ROOT.TGraph       . __rsub__      ,
    ROOT.TGraph       . __isub__      ,
    # 
    ROOT.TGraph       .  __lshift__   , 
    ROOT.TGraph       .  __rshift__   , 
    ROOT.TGraph       . __ilshift__   , 
    ROOT.TGraph       . __irshift__   , 
    #
    ROOT.TGraph       . transform     ,
    #
    ROOT.TGraph . __exp__    , 
    ROOT.TGraph . __exp2__   ,
    ROOT.TGraph . __expm1__  , 
    ROOT.TGraph . __sqrt__   , 
    ROOT.TGraph . __cbrt__   , 
    ROOT.TGraph . __log__    , 
    ROOT.TGraph . __log2__   , 
    ROOT.TGraph . __log10__  , 
    ROOT.TGraph . __log1p__  ,
    ROOT.TGraph . __sin__    ,
    ROOT.TGraph . __cos__    ,
    ROOT.TGraph . __tan__    ,
    ROOT.TGraph . __sinh__   ,
    ROOT.TGraph . __cosh__   ,
    ROOT.TGraph . __tanh__   ,
    ROOT.TGraph . __asin__   ,
    ROOT.TGraph . __acos__   ,
    ROOT.TGraph . __atan__   ,
    ROOT.TGraph . __asinh__  ,
    ROOT.TGraph . __acosh__  ,
    ROOT.TGraph . __atanh__  ,
    ROOT.TGraph . __erf__    ,
    ROOT.TGraph . __erfc__   ,
    ROOT.TGraph . __erfcx__  ,
    ROOT.TGraph . __erfi__   ,
    ROOT.TGraph . __tgamma__ ,
    ROOT.TGraph . __lgamma__ ,
    ROOT.TGraph . __sech__   ,
    ROOT.TGraph . __probit__ ,
    ROOT.TGraph . __pow__    ,
    #
    ROOT.TGraph . append     ,
    ROOT.TGraph . swap       ,    
    #
    ROOT.TGraphErrors . __getitem__   ,
    ROOT.TGraphErrors . __setitem__   ,
    ROOT.TGraphErrors .     items     ,
    ROOT.TGraphErrors . iteritems     ,
    #
    ROOT.TGraphErrors . xmin          ,
    ROOT.TGraphErrors . ymin          ,
    ROOT.TGraphErrors . xmax          ,
    ROOT.TGraphErrors . ymax          , 
    ROOT.TGraphErrors .  bb           ,
    #
    ROOT.TGraphErrors .  __mul__      , 
    ROOT.TGraphErrors . __rmul__      , 
    ROOT.TGraphErrors . __imul__      ,
    #
    ROOT.TGraphErrors .  __div__      ,
    ROOT.TGraphErrors . __rdiv__      ,
    ROOT.TGraphErrors . __idiv__      ,
    # 
    ROOT.TGraphErrors .  __truediv__  ,
    ROOT.TGraphErrors . __itruediv__  ,
    ROOT.TGraphErrors . __rtruediv__  ,
    #
    ROOT.TGraphErrors .  __add__      ,
    ROOT.TGraphErrors . __radd__      ,
    ROOT.TGraphErrors . __iadd__      ,
    ROOT.TGraphErrors .  __sub__      ,
    ROOT.TGraphErrors . __rsub__      ,
    ROOT.TGraphErrors . __isub__      ,
    # 
    ROOT.TGraphErrors .  __lshift__   , 
    ROOT.TGraphErrors .  __rshift__   , 
    ROOT.TGraphErrors . __ilshift__   , 
    ROOT.TGraphErrors . __irshift__   , 
    #
    ROOT.TGraphAsymmErrors.__len__       ,
    ROOT.TGraphAsymmErrors.__contains__  ,
    ROOT.TGraphAsymmErrors.__iter__      ,
    ROOT.TGraphAsymmErrors.     items    ,
    ROOT.TGraphAsymmErrors. iteritems    ,
    ROOT.TGraphAsymmErrors.__getitem__   ,
    ROOT.TGraphAsymmErrors.__setitem__   ,
    # 
    ROOT.TGraphAsymmErrors . xmin        ,
    ROOT.TGraphAsymmErrors . ymin        ,
    ROOT.TGraphAsymmErrors . xmax        ,
    ROOT.TGraphAsymmErrors . ymax        ,
    ROOT.TGraphAsymmErrors . bb          ,
    #
    ROOT.TGraphAsymmErrors .  __mul__      , 
    ROOT.TGraphAsymmErrors . __rmul__      , 
    ROOT.TGraphAsymmErrors . __imul__      ,
    #
    ROOT.TGraphAsymmErrors .  __div__      ,
    ROOT.TGraphAsymmErrors . __rdiv__      ,
    ROOT.TGraphAsymmErrors . __idiv__      ,
    # 
    ROOT.TGraphAsymmErrors .  __truediv__  ,
    ROOT.TGraphAsymmErrors . __itruediv__  ,
    ROOT.TGraphAsymmErrors . __rtruediv__  ,
    #
    ROOT.TGraphAsymmErrors .  __add__      ,
    ROOT.TGraphAsymmErrors . __radd__      ,
    ROOT.TGraphAsymmErrors . __iadd__      ,
    ROOT.TGraphAsymmErrors .  __sub__      ,
    ROOT.TGraphAsymmErrors . __rsub__      ,
    ROOT.TGraphAsymmErrors . __isub__      ,
    # 
    ROOT.TGraphAsymmErrors .  __lshift__   , 
    ROOT.TGraphAsymmErrors .  __rshift__   , 
    ROOT.TGraphAsymmErrors . __ilshift__   , 
    ROOT.TGraphAsymmErrors . __irshift__   , 
    #
    ROOT.TGraph       . integral          ,
    ROOT.TGraph       . asTF1             ,
    # 
    ROOT.TGraph            .__getslice__  ,
    ROOT.TGraphErrors      .__getslice__  ,
    ROOT.TGraphAsymmErrors .__getslice__  ,
    #
    ROOT.TGraph            .sorted        ,
    ROOT.TGraphErrors      .sorted        ,
    ROOT.TGraphAsymmErrors .sorted        ,
    #
    ROOT.TGraph            .filter        ,
    ROOT.TGraph            .remove        ,
    ##
    ROOT.TGraph.transpose                 ,
    ROOT.TGraph.T                         ,     
    ROOT.TGraphErrors.transpose           ,
    ROOT.TGraphErrors.T                   ,
    ROOT.TGraphAsymmErrors.transpose      ,
    ROOT.TGraphAsymmErrors.T              ,
    ##
    _color_   ,
    _red_     ,
    _blue_    ,
    _magenta_ ,
    _cyan_    ,
    _green_   ,
    _yellow_  ,
    #
    ROOT.TGraph  . xmin    ,
    ROOT.TGraph  . xmax    ,
    ROOT.TGraph  . xminmax ,
    #
    ROOT.TMultiGraph.color     ,
    ROOT.TMultiGraph.red       ,
    ROOT.TMultiGraph.blue      , 
    ROOT.TMultiGraph.magenta   , 
    ROOT.TMultiGraph.cyan      ,
    ROOT.TMultiGraph.green     ,
    ROOT.TMultiGraph.yellow    ,
    #
    ROOT.TMultiGraph.transpose ,
    ROOT.TMultiGraph.T         ,
    #
    ROOT.TH1D.lw_graph         , 
    ROOT.TH1F.lw_graph         ,
    ##
    ROOT.TArrow.transpose      ,
    ROOT.TArrow.T              ,
    ROOT.TBox.transpose        ,
    ROOT.TBox.T                ,
    ROOT.TLine.transpose       ,
    ROOT.TLine.T               ,
    ROOT.TText.transpose       ,
    ROOT.TText.T               ,
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
