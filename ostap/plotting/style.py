#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file
# LHCb style file for ROOT-plots
#
#                   $Revision$
# Last modification $Date$
#                by $Author$
# =============================================================================
"""LHCb Style for ROOT-plots"""
# =============================================================================
import ROOT
__all__ = (
    'Style'      ,
    'label'      ,
    'latex'      , 
    'font'       , 
    'ostapStyle' ,
    'ostapLabel' ,
    'ostapLatex' ,
    'ostapFont'  ,
    'OstapStyle' 
    )
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.plotting.style' )
else                       : logger = getLogger( __name__ )
# =============================================================================

## font 
font      = 132     ## Times-Roman 
## line thickness
lineWidth =   1

##define style for text
label = ROOT.TText   (    )
label . SetTextFont  ( font )
label . SetTextColor (  1 )
label . SetTextSize  (  0.04 )
label . SetTextAlign ( 12 )

## define style of latex text
latex = ROOT.TLatex   () 
latex . SetTextFont   ( font )
latex . SetTextColor  ( 1    )
latex . SetTextSize   ( 0.04 )
latex . SetTextAlign  ( 12   )

ostapLatex = latex
ostapLabel = label
ostapFont  = font

# =============================================================================
## define LHCb style for plots 
def OstapStyle ( name      = "OstapStyle"           ,
                 desc      = "Standard plots style" ,
                 lineWidth = lineWidth              ,
                 font      = font                   ,
                 makeNew   = False                  ,
                 force     = True                   ) :
    """
    Define Ostap style for plots
    """
    obj = ROOT.gROOT.FindObject ( name )
    if obj and issubclass ( type( obj  ) , ROOT.TStyle ) and not makeNew : 
        logger.info              ('The style %s is reused' % obj.GetName() )
        if force : 
            logger.info          ('The style %s is forced' % obj.GetName() )
            ROOT.gROOT.SetStyle   ( obj.GetName()  )
            ROOT.gROOT.ForceStyle ( )
        return obj
    
    nam = name
    i   = 1
    while obj :
        nam  = name + '_%d' % i
        obj  = ROOT.gROOT.FindObject ( nam )
        i   += 1
        
    style = ROOT.TStyle ( nam , desc )
    logger.debug ('New style %s is created' % style.GetName() )
    
    ## use plain black on white colors
    style . SetFrameBorderMode  ( 0 )
    style . SetCanvasBorderMode ( 0 )
    style . SetPadBorderMode    ( 0 )
    style . SetPadColor         ( 0 )
    style . SetCanvasColor      ( 0 )
    style . SetStatColor        ( 0 )
    style . SetPalette          ( 1 )
    
    ## set the paper & margin sizes
    style . SetPaperSize         ( 20 , 26 )
    style . SetPadTopMargin      ( 0.05    )
    style . SetPadRightMargin    ( 0.05    ) ## increase for colz plots
    style . SetPadBottomMargin   ( 0.16    )
    style . SetPadLeftMargin     ( 0.14    )
        
    ##  use large fonts
    style . SetTextFont          ( font )
    style . SetTextSize          ( 0.08 )
    style . SetLabelFont         ( font , "x" ) 
    style . SetLabelFont         ( font , "y" ) 
    style . SetLabelFont         ( font , "z" ) 
    style . SetLabelSize         ( 0.05 , "x" )
    style . SetLabelSize         ( 0.05 , "y" )
    style . SetLabelSize         ( 0.05 , "z" )
    style . SetTitleFont         ( font )
    style . SetTitleSize         ( 0.06 , "x" )
    style . SetTitleSize         ( 0.06 , "y" )
    style . SetTitleSize         ( 0.06 , "z" ) 
        
    ## use bold lines and markers
    style . SetLineWidth         ( lineWidth )
    style . SetFrameLineWidth    ( lineWidth )
    style . SetHistLineWidth     ( lineWidth )
    style . SetFuncWidth         ( lineWidth )
    style . SetGridWidth         ( lineWidth )
    style . SetLineStyleString   ( 2 , "[12 12]" ) ##  postscript dashes
    style . SetMarkerStyle       ( 20  )
    style . SetMarkerSize        ( 1.2 )

    ## ## colors for 2D
    ## from array import array
    ## NRGBs =   4
    ## NCont = 999
    ## #
    ## stops_ = array ( 'd' , [ 0.00 , 0.33 , 0.66 , 1.0 ] ) 
    ## red_   = array ( 'd' , [ 0.00 , 0.00 , 0.4  , 1.0 ] )  
    ## green_ = array ( 'd' , [ 0.00 , 0.00 , 1.0  , 1.0 ] )
    ## blue_  = array ( 'd' , [ 0.00 , 1.00 , 0.4  , 0.6 ] )
    ## ROOT.TColor.CreateGradientColorTable ( NRGBs  ,
    ##                                        stops_ ,
    ##                                        red_   ,
    ##                                        green_ ,
    ##                                        blue_  ,
    ##                                        NCont  )
    ## style.SetNumberContours ( NCont ) 
    
    # const Int_t NRGBs = 4;
    # const Int_t NCont = 999;
    # Double_t stops[NRGBs] = { 0.00,0.33,0.66,1.0};
    # Double_t red[NRGBs]     = {0.0,0.00,0.4,1.0};
    # Double_t green[NRGBs] = { 0.0,0.00,1.0,1.0};
    # Double_t blue[NRGBs] = { 0.0,1.00,0.4,0.6};
    # TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue , NCont);
    # gStyle->SetNumberContours(NCont);
    
    ## label offsets
    style . SetLabelOffset(0.015);
    
    ## by default, do not display histogram decorations:
    style . SetOptStat    ( 0 )  
    ## lhcbStyle->SetOptStat("emr");  ##  show only nent -e , mean - m , rms -r
    ## full opts at http://root.cern.ch/root/html/TStyle.html#TStyle:SetOptStat
    style . SetStatFormat ("6.3g") ##  specified as c printf options
    style . SetOptTitle   ( 0 )
    style . SetOptFit     ( 0 )
    ## lhcbStyle . SetOptFit(1011); // order is probability, Chi2, errors, parameters

    ## size of small lines at the end of error bars
    style.  SetEndErrorSize ( 5 ) 

    ## look of the statistics box:
    style . SetStatBorderSize         ( 0 )
    style . SetStatFont               ( font )
    style . SetStatFontSize           ( 0.05 )
    style . SetStatX                  ( 0.9  )
    style . SetStatY                  ( 0.9  )
    style . SetStatW                  ( 0.25 )
    style . SetStatH                  ( 0.15 )
    ## put tick marks on top and RHS of plots
    style . SetPadTickX           ( 1 )
    style . SetPadTickY           ( 1 )
    
    ## histogram divisions: only 5 in x to avoid label overlaps
    style . SetNdivisions    ( 505 , "x" )
    style . SetNdivisions    ( 510 , "y" )

    
    ## few useful line styles 
    style.SetLineStyleString(11,"76  24"       );
    style.SetLineStyleString(12,"60  16 8 16"  );
    style.SetLineStyleString(13,"168 32"       );
    style.SetLineStyleString(14,"32  32"       );

    if force : 
        logger.debug ('The style %s is forced' % style.GetName() )
        ROOT.gROOT.SetStyle   ( style.GetName()  )
        ROOT.gROOT.ForceStyle ()
    
    return style 

Style      = OstapStyle() 
ostapStyle = Style


# =============================================================================
if '__main__' == __name__ :
         
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
