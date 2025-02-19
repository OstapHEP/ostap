#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  redcue (T)Graph-like classes 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07 
# =============================================================================
"""reduce (T)Graph-like classes 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    ) 
# =============================================================================
from   ostap.math.reduce              import root_factory
from   ostap.core.ostap_types         import num_types, integer_types, string_types 
from   ostap.plotting.draw_attributes import copy_graph_attributes  
import ROOT, array, itertools  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.histos.graph_redcue' )
else                       : logger = getLogger( __name__                    )
# =============================================================================
logger.debug ( 'reduce (T)Graph-like classes')
# ===============================================================================
## reconstruct/deserialize/unpickle <code>TGraph</code> object
#  @see TGraph
def graph_factory ( klass     ,
                    name      ,
                    title     ,
                    minmax    ,
                    attline   ,
                    attfill   ,
                    attmarker , 
                    xvalues   ,
                    yvalues   ) :
    """Reconstruct/deserialize/unpickle `ROOT.TGraph` object
    -see `ROOT.TGraph`
    """
    graph = klass ()
    ## TNamed
    graph.SetName   ( name  ) 
    graph.SetTitle  ( title )
    ## min/max 
    graph.SetMinimum ( minmax[0] )
    graph.SetMaximum ( minmax[1] )    
    ## TAttLine
    att = ROOT.TAttLine    ( *attline ) 
    att.Copy (  graph )
    ## TAttFill
    att = ROOT.TAttFill   ( *attfill   ) 
    att.Copy (  graph )
    ## TAttMarker
    att = ROOT.TAttMarker ( *attmarker ) 
    att.Copy (  graph )
    ## reconstruct the content of TGraph
    graph.Set ( len ( xvalues )  )
    for i, xv, yv in zip ( range ( graph.GetN () ) , xvalues , yvalues ) :
        graph.SetPoint ( i , xv , yv )        
    ## 
    return graph 

# =============================================================================
## Reduce/serialize/pickle  simple <code>TGraph</code> object
#  @see ROOT.TGraph
def graph_reduce ( graph ) :
    """Reduce/serialize/pickle simple `ROOT.TGraph` object\
    - see `ROOT.TGraph`
    """
    xvalues = array.array  ( 'd' ,  graph.GetX () ) 
    yvalues = array.array  ( 'd' ,  graph.GetY () )     
    return graph_factory , ( type ( graph ) ,
                             ## TNamed                             
                             str ( graph.GetName    () ) ,
                             str ( graph.GetTitle   () ) ,                             
                             ## minmax
                             ( graph.GetMinimum     () ,
                               graph.GetMaximum     () ) ,  
                             ## TAttLine 
                             ( graph.GetLineColor   () ,
                               graph.GetLineStyle   () ,
                               graph.GetLineWidth   () ) , 
                             ## TAttFill
                             ( graph.GetFillColor   () ,
                               graph.GetFillStyle   () ) ,  
                             ## TAttMarker 
                             ( graph.GetMarkerColor () ,
                               graph.GetMarkerStyle () ,
                               graph.GetMarkerSize  () ) ,                              
                             ## finally: the content 
                             xvalues ,
                             yvalues )

ROOT.TGraph.__reduce__ = graph_reduce

# ===============================================================================
## reconstruct/deserialize/unpickle <code>TGraphErrors</code> object
#  @see TGraphErrors
def graph_errors_factory  ( klass     ,
                            name      ,
                            title     ,
                            minmax    ,
                            attline   ,
                            attfill   ,
                            attmarker , 
                            xvalues   ,
                            yvalues   , 
                            xerrors   ,
                            yerrors   ) : 
    """Reconstruct/deserialize/unpickle `ROOT.TGraphErrors` object
    -see `ROOT.TGraphErrors`
    """
    graph = graph_factory ( klass , name , title , minmax , attline , attfill , attmarker , xvalues , yvalues )
    ##
    for i, ex, ey in zip ( range ( graph.GetN()  ) , xerrors , yerrors ) : 
        graph.SetPointError ( i , ex , ey )
    ##
    return graph

# =============================================================================
## Reduce/serialize/pickle  simple <code>TGraphErrors</code> object
#  @see ROOT.TGraphErorrs
def graph_errors_reduce ( graph ) :
    """Reduce/serialize/pickle simple `ROOT.TGraphErrors` object
    - see `ROOT.TGraphErrors`
    """
    ##
    _ , content = graph_reduce ( graph )
    ##
    N = graph.GetN() 
    xerrors = array.array ( 'd' , ( graph.GetErrorX ( i ) for i in range ( N ) ) )
    yerrors = array.array ( 'd' , ( graph.GetErrorY ( i ) for i in range ( N ) ) )
    ##
    return graph_errors_factory , content + ( xerrors , yerrors ) 

ROOT.TGraphErrors.__reduce__ = graph_errors_reduce

# ===============================================================================
## reconstruct/deserialize/unpickle <code>TGraphAsymmErrors</code> object
#  @see TGraphAsymmErrors
def graph_asymerrors_factory  ( klass     ,
                                name      ,
                                title     ,
                                minmax    ,
                                attline   ,
                                attfill   ,
                                attmarker , 
                                xvalues   ,
                                yvalues   , 
                                xlerrors  ,
                                xherrors  ,
                                ylerrors  ,
                                yherrors  ) :
    """Reconstruct/deserialize/unpickle `ROOT.TGraphAsymmErrors` object
    -see `ROOT.TGraphAsymmErrors`
    """
    graph = graph_factory ( klass , name , title , minmax , attline , attfill , attmarker , xvalues , yvalues )
    ##
    for i, exl, exh , eyl , eyh  in zip ( range ( graph.GetN() ) ,
                                          xlerrors , xherrors ,
                                          ylerrors , yherrors ) :  
        graph.SetPointEXlow  ( i , exl )
        graph.SetPointEXhigh ( i , exh )
        graph.SetPointEYlow  ( i , eyl )
        graph.SetPointEYhigh ( i , eyh )
    ##
    return graph

# =============================================================================
## Reduce/serialize/pickle  simple <code>TGraphAsymmErrors</code> object
#  @see ROOT.TGraphErorrs
def graph_asymerrors_reduce ( graph ) :
    """Reduce/serialize/pickle simple `ROOT.TGraphAsymmErrors` object
    - see `ROOT.TGraphAsymmErrors`
    """
    ##
    _ , content = graph_reduce ( graph )
    ## 
    N = graph.GetN() 
    xlerrors = array.array ( 'd' , ( graph.GetErrorXlow  ( i ) for i in range ( N ) ) )
    xherrors = array.array ( 'd' , ( graph.GetErrorXhigh ( i ) for i in range ( N ) ) )
    ylerrors = array.array ( 'd' , ( graph.GetErrorYlow  ( i ) for i in range ( N ) ) )
    yherrors = array.array ( 'd' , ( graph.GetErrorYhigh ( i ) for i in range ( N ) ) )
    ##
    return graph_asymerrors_factory , content + ( xlerrors ,
                                                  xherrors ,
                                                  ylerrors ,
                                                  yherrors )

ROOT.TGraphAsymmErrors.__reduce__ = graph_asymerrors_reduce

# =============================================================================
## reconstruct/deserialize/unpickle <code>TGraphMultiErrors</code> object
#  @see TGraphMultiErrors
def graph_multierrors_factory  ( klass     ,
                                 name      ,
                                 title     ,
                                 minmax    ,
                                 attline   ,
                                 attfill   ,
                                 attmarker , 
                                 xvalues   ,
                                 yvalues   ,
                                 mode      , 
                                 xlerrors  ,
                                 xherrors  ,
                                 ylerrors  ,
                                 yherrors  ,
                                 attsline  ,
                                 attsfill  ) :
    """ Reconstruct/deserialize/unpickle `ROOT.TGraphAsymmErrors` object
    - see `ROOT.TGraphAsymmErrors`
    """
    graph = graph_factory ( klass , name , title , minmax , attline , attfill , attmarker , xvalues , yvalues )
    ##
    graph.SetSumErrorsMode( mode ) 
    N  = graph.GetN()
    ## adjust number of Y-errors 
    NE = len ( ylerrors )
    za = array.array ( 'd' , N * [ 0.0 ] )
    while graph.GetNYErrors() < NE : graph.AddYError( N , za , za )
    ## 
    for i, exl, exh  in zip ( range ( N ) , xlerrors , xherrors ) :
        graph.SetPointEXlow  ( i , exl )
        graph.SetPointEXhigh ( i , exh )
        for e, yle , yhe in zip ( range ( NE ) , xlerrors , yherrors ) :
            graph.SetPointEYlow  ( i , e , eyl [ i ] )
            graph.SetPointEYhigh ( i , e , eyh [ i ] )
            
    ## line attributes 
    for e , aline in zip ( range ( NE ) , attsline ) :
        la = ROOT.TAttLine ( *aline )
        al = graph.GetAttLine ( e )
        if al : la.Copy ( al )
        
    ## fill attributes 
    for e , afill in zip ( range ( NE ) , attsfill ) :
        fa = ROOT.TAttFill ( *afill )
        af = graph.GetAttFill ( e )
        if af : fa.Copy ( af )  

    return graph


# =============================================================================
## Reduce/serialize/pickle  simple <code>TGraphMultiErrors</code> object
#  @see ROOT.TGraphMultiErorrs
def graph_multierrors_reduce ( graph ) :
    """ Reduce/serialize/pickle simple `ROOT.TGraphMultiErrors` object
    - see `ROOT.TGraphMultiErrors`
    """
    ##
    _ , content = graph_reduce ( graph )
    ## 
    N  = graph.GetN        ()
    xlerrors = array.array ( 'd' , ( graph.GetErrorXlow  ( i ) for i in range ( N ) ) )
    xherrors = array.array ( 'd' , ( graph.GetErrorXhigh ( i ) for i in range ( N ) ) )
    ##
    NE = graph.GetNYErrors ()
    ylerrors = []
    yherrors = []
    attsline = []
    attsfill = []
    for e in range ( NE ) :  
        ylerrors.append ( array.array ( 'd' , ( graph.GetErrorYlow  ( i , e ) for i in range ( N ) ) ) )
        yherrors.append ( array.array ( 'd' , ( graph.GetErrorYhigh ( i , e ) for i in range ( N ) ) ) )
        ## line attributes 
        aline = graph.GetAttLine ( e )
        if aline : attl = aline.GetLineColor () , aline.GetLineStyle (), aline.GetLineWidth ()
        else     : attl = 1 , 1 , 1 
        attsline.append ( attl )
        ## fill attribtes 
        afill = graph.GetAttFill ( e )
        if afill : attf = afill.GetFillColor () , afill.GetFillStyle ()
        else     : attf = 1 , 1000 
        attsfill.append ( attf )
        
    ##
    ylerrors = tuple ( ylerrors )
    yherrors = tuple ( yherrors )
    attsline = tuple ( attsline )
    attsfill = tuple ( attsfill )
    
    ##
    return graph_multierrors_factory , content + ( graph.GetSumErrorsMode() ,
                                                   xlerrors ,
                                                   xherrors ,
                                                   ylerrors ,
                                                   yherrors ,
                                                   attsline ,
                                                   attsfill )


ROOT.TGraphMultiErrors.__reduce__ = graph_multierrors_reduce

# =============================================================================

_new_methods_ = [
    #
    'graph_factory'             ,
    'graph_errors_factory'      ,
    'graph_asymmerrors_factory' ,
    #
    'graph_reduce'              ,
    'graph_errors_reduce'       ,
    'graph_asymmerrors_reduce'  ,    
    ] 

_decorated_classes_  = (
    ROOT.TGraph            ,
    ROOT.TGraphErrors      ,
    ROOT.TGraphAsymmErrors ,
    )

for t in _decorated_classes_ :
    _new_methods_.append ( t.__reduce__  )

_new_methods_ = tuple ( _new_methods_ ) 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
