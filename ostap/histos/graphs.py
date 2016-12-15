#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
# =============================================================================
## @file
#  (T)Graph-related decorations 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07 
# =============================================================================
"""TGraph-related decorations"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'makeGraph'   , ## make graph from primitive data
    'makeGraph2'  , ## make graph from plain input two-column text 
    'makeGraph3'  , ## make graph with errors  from plain input four-column text 
    'makeGraphs3' , ## make graphs from plain input multicolumn text 
    'makeGraphs4' , ## make graphs from plain input multicolumn text 
    'hToGraph'    , ## convert histogram to graph 
    'hToGraph2'   , ## convert histogram to graph 
    'hToGraph3'   , ## convert histogram to graph
    'lw_graph'    , ## make Laffery-Wyatt's graph
    ##
    ) 
# =============================================================================
import ROOT, cppyy              ## attention here!!
from ostap.core.core import cpp, VE
# 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.histos.graphs' )
else                       : logger = getLogger( __name__              )
# =============================================================================
logger.debug ( '(T)Graph-related decorations')
# =============================================================================
## copy graph attributes
#  - LineColor
#  - LineWidth
#  - LineStyle
#  - MarkerColor 
#  - MarkerSize
#  - MarkerStyle 
#  - FillColor
#  - FillStyle  
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07  
def copy_graph_attributes ( o_from , o_to ) :
    """Copy graph attributes
    - LineColor
    - LineWidth
    - LineStyle
    - MarkerColor 
    - MarkerSize
    - MarkerStyle 
    - FillColor
    - FillStyle  
    """
    ##
    ## line attributes:
    ## 
    if hasattr ( o_from , 'GetLineColor' ) and hasattr ( o_to , 'SetLineColor' ) :
        o_to.SetLineColor   ( o_from.GetLineColor () ) 
    if hasattr ( o_from , 'GetLineWidth' ) and hasattr ( o_to , 'SetLineWidth' ) :
        o_to.SetLineWidth   ( o_from.GetLineWidth () ) 
    if hasattr ( o_from , 'GetLineStyle' ) and hasattr ( o_to , 'SetLineStyle' ) :
        o_to.SetLineStyle   ( o_from.GetLineStyle () ) 
    ##
    ## marker attributes:
    ## 
    if hasattr ( o_from , 'GetMarkerColor' ) and hasattr ( o_to , 'SetMarkerColor' ) :
        o_to.SetMarkerColor ( o_from.GetMarkerColor () ) 
    if hasattr ( o_from , 'GetMarkerSize'  ) and hasattr ( o_to , 'SetMarkerSize'  ) :
        o_to.SetMarkerSize  ( o_from.GetMarkerSize  () ) 
    if hasattr ( o_from , 'GetMarkerStyle' ) and hasattr ( o_to , 'SetMarkerStyle' ) :
        o_to.SetMarkerStyle ( o_from.GetMarkerStyle () ) 
    ##
    ## Fill attributes:
    ##
    if hasattr ( o_from , 'GetFillColor' ) and hasattr ( o_to , 'SetFillColor' ) :
        o_to.SetFillColor   ( o_from.GetFillColor () ) 
    if hasattr ( o_from , 'GetFillStyle' ) and hasattr ( o_to , 'SetFillStyle' ) :
        o_to.SetFillStyle   ( o_from.GetFillStyle () ) 
    
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
        keys.sort()
        for k in keys :
            _x += [   k  ]
            _y += [ x[k] ] 
        return makeGraph ( _x , _y )
        
    if  not x : raise TypeError, "X is not a proper vector!"
    if  not y : raise TypeError, "Y is not a proper vector!"
    if len( x ) != len ( y ) :
        raise TypeError, "Mismatch X/Y-lengths"

    if ex and len(ex) != len(x) : raise TypeError, "Mismatch X/eX-lengths"
    if ey and len(ey) != len(y) : raise TypeError, "Mismatch Y/eY-lengths"

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
#  # PDF set = CTEQ6.6
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

    for i in h1.iteritems () :

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

    for i in h1.iteritems () :
        
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
        raise ValueError, ' Illegal value for "bias" parameter '
    
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
    for p in h1.iteritems() :
        x = p[1]
        if x.error() < abs ( bias ) :
            raise ValueError, ' Illegal value for "bias" parameter '
        
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
## use graph as function 
#  @code
#  graph = ...
#  y     = graph ( 0.2 )
#  @endcode
#  @see TGraph.
#  @see TGraph::Eval
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gr_call_ ( graph , x , spline = None , opts = 'S' ) :
    """ Use graph as function
    >>> graph = ...
    >>> y     = graph ( 0.2 ) 
    """
    if not spline : spline = ROOT.MakeNullPointer(ROOT.TSpline)
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
def _gr_integral_ ( graph , xlow , xhigh , scipy = True ) :
    """Calculate an integral over the range \f$x_{low} \le x \le x_{high}\f$
    It is not very efficient, but OK 
    >>> graph = ...
    >>> i     = graph.integral ( 0 , 1 )
    """
    if scipy :
        from LHCbMath.deriv import integral 
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
    if not ipoint in graph : raise IndexError 
    #
    
    x_ = ROOT.Double(0)
    v_ = ROOT.Double(0)
    
    graph.GetPoint ( ipoint , x_ , v_ )
    
    return x_,v_

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


# =============================================================================
## iterate over the points in TGraph
#  @code 
#  gr = ...
#  for i,x,v in gr.iteritems(): ...
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gr_iteritems_ ( graph ) :
    """Iterate over graph points 
    >>> graph = ...
    >>> for i,x,v in graph.iteritems(): ...
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
    if not ipoint in graph : raise IndexError 
    #
    
    x_ = ROOT.Double(0)
    v_ = ROOT.Double(0)
    
    graph.GetPoint ( ipoint , x_ , v_ )
    
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
def _gr_as_TF1_ ( graph ) :
    """ Represent TGraph as TF1
    >>> graph = ...
    >>> fun   = grap.asTF1() 
    """
    from HParamDeco import _h1_as_fun_ 
    return _h1_as_fun_ ( graph , lambda x : x ) 
    
# =============================================================================
## iterate over points in TGraphErrors
#  @code
#  gre = ...
#  for i,x,v in gre.iteritems(): ...
#  @endcode
#  @see TGraphErrors
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gre_iteritems_ ( graph ) :
    """ Iterate over graph points 
    >>> gre = ...
    >>> for i,x,v in gre.iteritems(): ...
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
    if not ipoint in graph : raise IndexError 
    #
    
    x_ = ROOT.Double(0)
    v_ = ROOT.Double(0)
    
    graph.GetPoint ( ipoint , x_ , v_ )
    
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
    >>> for i,x,xl,xh,y,yl,yh in grae.iteritems(): ...
    
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
    xmn  = None
    np   = len(graph) 
    for ip in range( np ) :
        x , y = graph[ip] 
        if None == xmn or x <= xmn : xmn = x
    return xmn

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
    xmx  = None
    np   = len(graph) 
    for ip in range( np ) :
        x , y = graph[ip] 
        if None == xmx or x >= xmx : xmx = x
    return xmx

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
    ymn  = None
    np   = len(graph) 
    for ip in range( np ) :
        x , y = graph[ip] 
        if None == ymn or y <= ymn : ymn = y
    return ymn

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
    yxm  = None
    np   = len(graph) 
    for ip in range( np ) :
        x , y = graph[ip] 
        if None == ymx or y >= ymx : ymx = y
    return xmx

# =============================================================================
## get minimal-x 
#  @code
#  graph = ...
#  xmin  = graph.xmin () 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gre_xmin_ ( graph ) :
    """Get minimal x for the points
    >>> graph = ...
    >>> xmin  = graph.xmin () 
    """    
    xmn  = None
    np   = len(graph) 
    for ip in range( np ) :
        x , y = graph[ip]
        x = x.value() - x.error() 
        if None == xmn or x <= xmn : xmn = x
    return xmn

# =============================================================================
## get maximal-x 
#  @code
#  graph = ...
#  xmax  = graph.xmax () 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gre_xmax_ ( graph ) :
    """Get maximal x for the points
    >>> graph = ...
    >>> xmax  = graph.xmax () 
    """    
    xmx  = None
    np   = len(graph) 
    for ip in range( np ) :
        x , y = graph[ip] 
        x = x.value() + x.error() 
        if None == xmx or x >= xmx : xmx = x
    return xmx

# =============================================================================
## get minimal-y
#  @code
#  graph = ...
#  ymin  = graph.ymin () 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gre_ymin_ ( graph ) :
    """Get minimal y for the points
    >>> graph = ...
    >>> ymin  = graph.ymin () 
    """    
    ymn  = None
    np   = len(graph) 
    for ip in range( np ) :
        x , y = graph[ip] 
        y = y.value() - y.error() 
        if None == ymn or y <= ymn : ymn = y
    return ymn

# =============================================================================
## get maximal-y
#  @code
#  graph = ...
#  ymax  = graph.ymax () 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _gre_ymax_ ( graph ) :
    """Get maximal x for the points
    >>> graph = ...
    >>> ymax  = graph.ymax () 
    """    
    yxm  = None
    np   = len(graph) 
    for ip in range( np ) :
        x , y = graph[ip] 
        y = y.value() + y.error() 
        if None == ymx or y >= ymx : ymx = y
    return xmx


# =============================================================================
## get minimal-x 
#  @code
#  graph = ...
#  xmin  = graph.xmin () 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _grae_xmin_ ( graph ) :
    """Get minimal x for the points
    >>> graph = ...
    >>> xmin  = graph.xmin () 
    """    
    xmn  = None
    np   = len(graph) 
    for ip in range( np ) :
        x , exl , exh , y , eyl , eyh = graph[ip] 
        x = x - abs( exl ) 
        if None == xmn or x <= xmn : xmn = x
    return xmn

# =============================================================================
## get maximal-x 
#  @code
#  graph = ...
#  xmax  = graph.xmax () 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _grae_xmax_ ( graph ) :
    """Get maximal x for the points
    >>> graph = ...
    >>> xmax  = graph.xmax () 
    """    
    xmx  = None
    np   = len(graph) 
    for ip in range( np ) :
        x , exl , exh , y , eyl , eyh = graph[ip] 
        x = x + abs( exh ) 
        if None == xmx or x >= xmx : xmx = x
    return xmx

# =============================================================================
## get minimal-y
#  @code
#  graph = ...
#  ymin  = graph.ymin () 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _grae_ymin_ ( graph ) :
    """Get minimal y for the points
    >>> graph = ...
    >>> ymin  = graph.ymin () 
    """    
    ymn  = None
    np   = len(graph) 
    for ip in range( np ) :
        x , exl , exh , y , eyl , eyh = graph[ip] 
        y = y - abs( eyl ) 
        if None == ymn or y <= ymn : ymn = y
    return ymn

# =============================================================================
## get maximal-y
#  @code
#  graph = ...
#  ymax  = graph.ymax () 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _grae_ymax_ ( graph ) :
    """Get maximal x for the points
    >>> graph = ...
    >>> ymax  = graph.ymax () 
    """    
    yxm  = None
    np   = len(graph) 
    for ip in range( np ) :
        x , exl , exh , y , eyl , eyh = graph[ip] 
        y = y + abs( eyh ) 
        if None == ymx or y >= ymx : ymx = y
    return xmx


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
    xmn  = graph.xmin()
    xmx  = graph.xmax()
    return xmn , xmx 

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
    ymn  = graph.ymin()
    ymx  = graph.ymax()
    return ymn , ymx 

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
    
    while i < 0 : i += nb
    while j < 0 : j += nb

    new_graph = ROOT.TGraph( j - i ) if  i < j else  ROOT.TGraph()
    copy_graph_attributes ( graph , new_graph )

    ii = 0 
    while i < j :
        new_graph[ ii ] = graph[i]
        ii +=1 
        i  +=1 
        
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
    
    while i < 0 : i += nb
    while j < 0 : j += nb
    
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
    
    while i < 0 : i += nb
    while j < 0 : j += nb
    
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
    
    oitems =        ( i for i in graph.iteritems() ) 
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
    
    oitems =        ( i for i in graph.iteritems() ) 
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
    
    oitems =        ( i for i in graph.iteritems() ) 
    sitems = sorted ( oitems , key = lambda s :s[1] , reverse = reverse )
    
    new_graph = ROOT.TGraphAsymmErrors ( len( graph ) )
    copy_graph_attributes ( graph , new_graph )

    ip = 0 
    for item in sitems :
        new_graph[ip] = item[1:]
        ip += 1

    return new_graph 

# =============================================================================
## filter points from the graph
#  @code
#  graph = ...
#  f     = graph.filter( lambda s : s[1]>0 ) 
#  @endcode
#  @date   2016-03-28 
def _gr0_filter_ ( graph , accept ):
    """Filter points from the graph
    >>> graph = ...
    >>> f     = graph.filter( lambda s : s[1]>0 ) 
    """
    oitems =        ( i for i in graph.iteritems() ) 
    fitems = filter ( accept , oitems ) 
    
    new_graph = ROOT.TGraph ( len( fitems ) )
    copy_graph_attributes ( graph , new_graph )
    
    ip = 0 
    for item in fitems :
        new_graph[ip] = item[1:]
        ip += 1

    return new_graph

# =============================================================================
## filter points from the graph
#  @code
#  graph = ...
#  f     = graph.filter( lambda s : s[1]>0 ) 
#  @endcode
#  @date   2016-03-28 
def _gr1_filter_ ( graph , accept ):
    """Filter points from the graph
    >>> graph = ...
    >>> f     = graph.filter( lambda s : s[1]>0 ) 
    """
    oitems =        ( i for i in graph.iteritems() ) 
    fitems = filter ( accept , oitems ) 
    
    new_graph = ROOT.TGraphErrors ( len( fitems ) )
    copy_graph_attributes ( graph , new_graph )
    
    ip = 0 
    for item in fitems :
        new_graph[ip] = item[1:]
        ip += 1

    return new_graph


# =============================================================================
## filter points from the graph
#  @code
#  graph = ...
#  f     = graph.filter( lambda s : s[1]>0 ) 
#  @endcode
#  @date   2016-03-28 
def _gr2_filter_ ( graph , accept ):
    """Filter points from the graph
    >>> graph = ...
    >>> f     = graph.filter( lambda s : s[1]>0 ) 
    """
    
    oitems =        ( i for i in graph.iteritems() ) 
    fitems = filter ( accept , oitems ) 
    
    new_graph = ROOT.TGraphAsymmErrors ( len( fitems ) )
    copy_graph_attributes ( graph , new_graph )
    
    ip = 0 
    for item in fitems :
        new_graph[ip] = item[1:]
        ip += 1
        
    return new_graph

# =============================================================================
ROOT.TGraph       . __len__       = ROOT.TGraphErrors . GetN 
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
ROOT.TGraph       . iteritems     = _gr_iteritems_

ROOT.TGraphErrors . __getitem__   = _gre_getitem_ 
ROOT.TGraphErrors . __setitem__   = _gre_setitem_ 
ROOT.TGraphErrors . iteritems     = _gre_iteritems_ 


ROOT.TGraphErrors . xmin          = _gre_xmin_ 
ROOT.TGraphErrors . ymin          = _gre_ymin_ 
ROOT.TGraphErrors . xmax          = _gre_xmax_ 
ROOT.TGraphErrors . ymax          = _gre_ymax_ 


ROOT.TH1F.asGraph = hToGraph
ROOT.TH1D.asGraph = hToGraph
ROOT.TH1F.toGraph = hToGraph
ROOT.TH1D.toGraph = hToGraph

ROOT.TGraphAsymmErrors.__len__       = ROOT.TGraphAsymmErrors . GetN 
ROOT.TGraphAsymmErrors.__contains__  = lambda s,i : i in range(0,len(s))
ROOT.TGraphAsymmErrors.__iter__      = _gr_iter_ 
ROOT.TGraphAsymmErrors. iteritems    = _grae_iteritems_ 
ROOT.TGraphAsymmErrors.__getitem__   = _grae_getitem_ 
ROOT.TGraphAsymmErrors.__setitem__   = _grae_setitem_ 

ROOT.TGraphAsymmErrors . xmin        = _grae_xmin_ 
ROOT.TGraphAsymmErrors . ymin        = _grae_ymin_ 
ROOT.TGraphAsymmErrors . xmax        = _grae_xmax_ 
ROOT.TGraphAsymmErrors . ymax        = _grae_ymax_ 


ROOT.TGraph       . integral         = _gr_integral_
ROOT.TGraph       . asTF1            = _gr_as_TF1_


ROOT.TGraph            .__getslice__  = _gr0_getslice_
ROOT.TGraphErrors      .__getslice__  = _gr1_getslice_
ROOT.TGraphAsymmErrors .__getslice__  = _gr2_getslice_

ROOT.TGraph            .sorted        = _gr0_sorted_
ROOT.TGraphErrors      .sorted        = _gr1_sorted_
ROOT.TGraphAsymmErrors .sorted        = _gr2_sorted_ 


ROOT.TGraph            .filter        = _gr0_filter_
ROOT.TGraphErrors      .filter        = _gr1_filter_ 
ROOT.TGraphAsymmErrors .filter        = _gr2_filter_ 


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
    if 0 == _sise : return 1
    #
    _last = _size - 1
    x_ = ROOT.Double(0)
    v_ = ROOT.Double(0)    
    g.GetPoint ( _last , x_ , v_ )
    #
    return x_
# =============================================================================
## min-value for the graph)
def _gr_xmin_ ( g ) :
    """Get x-min for the graph    
    >>> xmin = graph.xmin()    
    """
    #
    _size = len ( graph )
    if 0 == _sise : return 0
    #
    x_ = ROOT.Double(0)
    v_ = ROOT.Double(0)    
    graph.GetPoint ( 0 , x_ , v_ )
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
#  @see http://dx.doi.org/10.1016/0168-9002(94)01112-5
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
    
    from scipy import integrate
    from scipy import optimize 
    
    for item in histo.iteritems () :

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
        fint  = integrate.quad ( lambda x : float ( func ( x ) ) , xmn , xmx ) 

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
            r0 = optimize.brentq ( lambda x : (float(func(x))-fx)  ,
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
#  @see http://dx.doi.org/10.1016/0168-9002(94)01112-5
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
_decorated_classes_ = (
    ROOT.TH1F              ,
    ROOT.TH1D              ,
    ROOT.TGraph            , 
    ROOT.TGraphErrors      ,
    ROOT.TGraphAsymmErrors 
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
    ROOT.TGraph       . xminmax       ,
    ROOT.TGraph       . yminmax       ,
    ROOT.TGraph       .  minmax       ,
    #
    ROOT.TGraph       . __getitem__   ,
    ROOT.TGraph       . __setitem__   ,
    ROOT.TGraph       . iteritems     ,
    #
    ROOT.TGraphErrors . __getitem__   ,
    ROOT.TGraphErrors . __setitem__   ,
    ROOT.TGraphErrors . iteritems     ,
    #
    ROOT.TGraphErrors . xmin          ,
    ROOT.TGraphErrors . ymin          ,
    ROOT.TGraphErrors . xmax          ,
    ROOT.TGraphErrors . ymax          , 
    #
    ROOT.TGraphAsymmErrors.__len__       ,
    ROOT.TGraphAsymmErrors.__contains__  ,
    ROOT.TGraphAsymmErrors.__iter__      ,
    ROOT.TGraphAsymmErrors. iteritems    ,
    ROOT.TGraphAsymmErrors.__getitem__   ,
    ROOT.TGraphAsymmErrors.__setitem__   ,
    # 
    ROOT.TGraphAsymmErrors . xmin        ,
    ROOT.TGraphAsymmErrors . ymin        ,
    ROOT.TGraphAsymmErrors . xmax        ,
    ROOT.TGraphAsymmErrors . ymax        ,
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
    ROOT.TGraphErrors      .filter        ,
    ROOT.TGraphAsymmErrors .filter        ,
    #
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
    ROOT.TMultiGraph.color    ,
    ROOT.TMultiGraph.red      ,
    ROOT.TMultiGraph.blue     ,
    ROOT.TMultiGraph.magenta  ,
    ROOT.TMultiGraph.cyan     ,
    ROOT.TMultiGraph.green    ,
    ROOT.TMultiGraph.yellow   ,
    #
    ROOT.TH1D.lw_graph        , 
    ROOT.TH1F.lw_graph        ,
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
