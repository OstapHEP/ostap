#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/pllotting/draw_attributes.py
#  Draw attrribute for line,fill, marker,... 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2020-04-11
# =============================================================================
"""Draw attrributes for line,fill, marker,... 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2020-04-11"
# =============================================================================
__all__     = (
    'set_line_attributes'   ,  ## set line attributes
    'set_fill_attributes'   ,  ## set fill attributes
    'set_marker_attributes' ,  ## set marker attributes
    'set_text_attributes'   ,  ## set text   attributes
    'copy_graph_attributes' ,  ## copy draw attributed from one object to another
    'copy_draw_attributes'  ,  ## ditto 
  )  
# =============================================================================
from   ostap.utils.cidict import cidict
import ROOT 
# =============================================================================
## set line attributes for the object
#  @code
#  set_line_attribute ( obj , line_color = 5 , LineStyle = 6 , ... ) 
#  @endocode
def set_line_attributes  ( obj , **kwargs ) :
    """Set line attributes for the object
    >>> set_line_attribute ( obj , line_color = 5 , LineStyle = 6 , ... ) 
    """
    
    key_transform = lambda k : k.lower().replace('_','').replace('line','')
    keys = cidict ( transform = key_transform )
    keys.update ( kwargs ) 

    if hasattr ( obj , 'SetLineStyle' ) :
        l = keys.get ( 'style' , None )
        if not l is None : obj.SetLineStyle ( l )  
    
    if hasattr ( obj , 'SetLineWidth' ) :
        l = keys.get ( 'line_width' , None )
        if not l is None : obj.SetLineWidth ( l )  

    if hasattr ( obj , 'SetLineColor' ) :
        l = keys.get ( 'color' , None )
        if not l is None : obj.SetLineColor ( l )  

    if hasattr ( obj , 'SetLineColorAlpha' ) and hasattr ( obj , 'GetLineColor' ) :
        l = keys.get ( 'color_alpha' , None )
        if not l is None : obj.SetLineColorAlpha ( object.GetLineColor() , l )  

# =============================================================================
## set fill attributes for the object
#  @code
#  set_fill_attribute ( obj , fill_color = 5 , FillStyle = 6 , ... ) 
#  @endocode
def set_fill_attributes  ( obj , **kwargs ) :
    """Set fill attributes for the object
    >>> set_fill_attribute ( obj , fill_color = 5 , FillStyle = 6 , ... ) 
    """
    
    key_transform = lambda k : k.lower().replace('_','').replace('fill','').replace('area','')
    keys = cidict ( transform = key_transform )
    keys.update ( kwargs ) 

    if hasattr ( obj , 'SetFillStyle' ) :
        l = keys.get ( 'style' , None )
        if not l is None : obj.SetFillStyle ( l )  

    if hasattr ( obj , 'SetFillColor' ) :
        l = keys.get ( 'color' , None )
        if not l is None : obj.SetFillColor ( l )  

    if hasattr ( obj , 'SetFillColorAlpha' ) and hasattr ( obj , 'GetFillColor' ) :
        l = keys.get ( 'color_alpha' , None )
        if not l is None : obj.SetFillColorAlpha ( object.GetFillColor() , l )  


# =============================================================================
## set marker attributes for the object
#  @code
#  set_marker_attribute ( obj , marker_color = 5 , MarkerStyle = 6 , ... ) 
#  @endocode
def set_marker_attributes  ( obj , **kwargs ) :
    """Set marker attributes for the object
    >>> set_marker_attribute ( obj , marker_color = 5 , markerStyle = 6 , ... ) 
    """
    
    key_transform = lambda k : k.lower().replace('_','').replace('marker','')
    keys = cidict ( transform = key_transform )
    keys.update ( kwargs ) 

    if hasattr ( obj , 'SetMarkerStyle' ) :
        l = keys.get ( 'style' , None )
        if not l is None : obj.SetMarkerStyle ( l )  

    if hasattr ( obj , 'SetMarkerColor' ) :
        l = keys.get ( 'color' , None )
        if not l is None : obj.SetMarkerColor ( l )  

    if hasattr ( obj , 'SetMarkerSize' ) :
        l = keys.get ( 'size' , None )
        if not l is None : obj.SetMarkerSize  ( l )  

    if hasattr ( obj , 'SetMarkerColorAlpha' ) and hasattr ( obj , 'GetMarkerColor' ) :
        l = keys.get ( 'color_alpha' , None )
        if not l is None : obj.SetMarkerColorAlpha ( object.GetMarkerColor() , l )  

# =============================================================================
## set text attributes for the object
#  @code
#  set_text_attribute ( obj , marker_color = 5 , MarkerStyle , ... ) 
#  @endocode
def set_text_attributes  ( obj , **kwargs ) :
    """Set text attributes for the object
    >>> set_text_attribute ( obj , marker_color = 5 , markerStyle = 6 , ... ) 
    """
    
    key_transform = lambda k : k.lower().replace('_','').replace('text','')
    keys = cidict ( transform = key_transform )
    keys.update ( kwargs ) 

    if hasattr ( obj , 'SetTextAlign' ) :
        l = keys.get ( 'align' , None )
        if not l is None : obj.SetTextAlign ( l )  

    if hasattr ( obj , 'SetTextAngle' ) :
        l = keys.get ( 'angle' , None )
        if not l is None : obj.SetTextAngle ( l )  

    if hasattr ( obj , 'SetTextFont' ) :
        l = keys.get ( 'font' , None )
        if not l is None : obj.SetTextFont  ( l )  

    if hasattr ( obj , 'SetTextSize' ) :
        l = keys.get ( 'size' , None )
        if not l is None : obj.SetTextSize  ( l )  

    if hasattr ( obj , 'SetTextColor' ) :
        l = keys.get ( 'color' , None )
        if not l is None : obj.SetTextColor ( l )  


    if hasattr ( obj , 'SetTextColorAlpha' ) and hasattr ( obj , 'GetTextColor' ) :
        l = keys.get ( 'color_alpha' , None )
        if not l is None : obj.SetTextColorAlpha ( object.GetTextColor() , l )  

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

    ##
    ## Copy All 
    ##
    if isinstance ( o_to , ROOT.TAttLine   ) and isinstance ( o_from , ROOT.TAttLine   ) :
        tmp = ROOT.TAttLine   ( o_from )
        tmp.Copy    ( o_to )
        
    if isinstance ( o_to , ROOT.TAttFill   ) and isinstance ( o_from , ROOT.TAttFill   ) :
        tmp = ROOT.TAttFill   ( o_from )
        tmp.Copy    ( o_to )
        
    if isinstance ( o_to , ROOT.TAttMarker ) and isinstance ( o_from , ROOT.TAttMarker ) :
        tmp = ROOT.TAttMarker ( o_from )
        tmp.Copy    ( o_to )

    if isinstance ( o_to , ROOT.TAttText   ) and isinstance ( o_from , ROOT.TAttText   ) :
        tmp = ROOT.TAttText ( o_from )
        tmp.Copy    ( o_to )

# =============================================================================
##  ditto 
copy_draw_attributes = copy_graph_attributes

        
        

ROOT.TAttLine   .set_line_attributes   = set_line_attributes
ROOT.TAttFill   .set_fill_attributes   = set_fill_attributes
ROOT.TAttMarker .set_marker_attributes = set_marker_attributes
ROOT.TAttText   .set_text_attributes   = set_text_attributes

_decorated_classes_ = (
    ROOT.TAttLine   ,
    ROOT.TAttFill   ,
    ROOT.TAttMarker ,
    ROOT.TAttText   ,
    )

# =============================================================================
if '__main__' == __name__ :
    
    from   ostap.logger.logger import getLogger
    logger = getLogger( 'ostap.plotting.draw_attributes' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
