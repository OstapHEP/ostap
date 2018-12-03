#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/plotting/style.py
# Ostap style file for ROOT-plots
# =============================================================================
"""Ostap Style for ROOT-plots"""
# =============================================================================
import ROOT
__all__ = (
    'UseStyle'    ,  ## context manager  for the style (class) 
    'useStyle'    ,  ## context manager  for the style (function)
    ##
    'Style'       ,  ## the default Ostap style 
    'label'       ,  ## the style for the text 
    'latex'       ,  ## the style for LaTeX 
    'font'        ,  ## the default  font 
    'ostapStyle'  ,  ## the default Ostap style 
    'ostapLabel'  ,  ## the default Ostap style for labels 
    'ostapLatex'  ,  ## the default Ostap style for LaTeX
    'ostapFont'   ,  ## he default Ostap font 
    'OstapStyle'  ,  ## the function to instantiate the style
    ## 
    'StyleZ'      ,  ## the style for COLZ plots
    'Style2'      ,  ## the style for downscaled 2-in-row plots
    'Style3'      ,  ## the style for downscaled 3-in-row plots
    'Style2Z'     ,  ## the style for downscaled 2-in-row COLZ plots
    'Style3Z'     ,  ## the style for downscaled 3-in-row COLZ plots
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
                 force     = True                   ,
                 scale     = 1.0                    ,
                 colz      = False                  ) :
    """

    """
    obj = ROOT.gROOT.FindObject  ( name )
    if obj and isinstance ( obj , ROOT.TStyle ) and not makeNew : 
        logger.info              ('The style %s is reused' % obj.GetName() )
        if force :
            obj.cd                () 
            logger.info           ('The style %s is forced' % obj.GetName() )
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


       ## set the paper & margin sizes
    style . SetPaperSize         ( 20 , 26 )
    style . SetPadTopMargin      ( 0.05    )
    style . SetPadRightMargin    ( 0.14 if colz else 0.05 )
    style . SetPadBottomMargin   ( 0.14    )
    style . SetPadLeftMargin     ( 0.14    )
        
    ##  use large fonts
    style . SetTextFont          ( font )
    style . SetTextSize          ( 0.08 * scale )
    style . SetLabelFont         ( font , "x" ) 
    style . SetLabelFont         ( font , "y" ) 
    style . SetLabelFont         ( font , "z" ) 
    style . SetLabelSize         ( 0.05 * scale , "x" )
    style . SetLabelSize         ( 0.05 * scale , "y" )
    style . SetLabelSize         ( 0.05 * scale , "z" )
    style . SetTitleFont         ( font )
    style . SetTitleSize         ( 0.06 * scale , "x" )
    style . SetTitleSize         ( 0.06 * scale , "y" )
    style . SetTitleSize         ( 0.06 * scale , "z" ) 
        
    ## use bold lines and markers
    style . SetLineWidth         ( lineWidth )
    style . SetFrameLineWidth    ( lineWidth )
    style . SetHistLineWidth     ( lineWidth )
    style . SetFuncWidth         ( lineWidth )
    style . SetGridWidth         ( lineWidth )
    style . SetLineStyleString   ( 2 , "[12 12]" ) ##  postscript dashes
    style . SetMarkerStyle       ( 20  )
    style . SetMarkerSize        ( 1.2 )

    ## label offsets
    style . SetLabelOffset(0.015);
    
    ## by default, do not display histogram decorations:
    style . SetOptStat    ( 0 )  
    style . SetStatFormat ("6.3g") ##  specified as c printf options
    style . SetOptTitle   ( 0 )
    style . SetOptFit     ( 0 )

    ## size of small lines at the end of error bars
    style.  SetEndErrorSize ( 3 ) 

    ## look of the statistics box:
    style . SetStatBorderSize         ( 0 )
    style . SetStatFont               ( font )
    style . SetStatFontSize           ( 0.05 * scale )
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
    
    ## add several useful line types:
    style . SetLineStyleString ( 11  ,"76  24"       )
    style . SetLineStyleString ( 12 , "60  16 8 16"  )
    style . SetLineStyleString ( 13 , "168 32"       )
    style . SetLineStyleString ( 14 , "32  32"       )

    ## Palettes:
    #  @see  https://root.cern.ch/doc/master/classTColor.html
    ## kDeepSea=51,          kGreyScale=52,    kDarkBodyRadiator=53,
    ## kBlueYellow= 54,      kRainBow=55,      kInvertedDarkBodyRadiator=56,
    ## kBird=57,             kCubehelix=58,    kGreenRedViolet=59,
    ## kBlueRedYellow=60,    kOcean=61,        kColorPrintableOnGrey=62,
    ## kAlpine=63,           kAquamarine=64,   kArmy=65,
    ## kAtlantic=66,         kAurora=67,       kAvocado=68,
    ## kBeach=69,            kBlackBody=70,    kBlueGreenYellow=71,
    ## kBrownCyan=72,        kCMYK=73,         kCandy=74,
    ## kCherry=75,           kCoffee=76,       kDarkRainBow=77,
    ## kDarkTerrain=78,      kFall=79,         kFruitPunch=80,
    ## kFuchsia=81,          kGreyYellow=82,   kGreenBrownTerrain=83,
    ## kGreenPink=84,        kIsland=85,       kLake=86,
    ## kLightTemperature=87, kLightTerrain=88, kMint=89,
    ## kNeon=90,             kPastel=91,       kPearl=92,
    ## kPigeon=93,           kPlum=94,         kRedBlue=95,
    ## kRose=96,             kRust=97,         kSandyTerrain=98,
    ## kSienna=99,           kSolar=100,       kSouthWest=101,
    ## kStarryNight=102,     kSunset=103,      kTemperatureMap=104,
    ## kThermometer=105,     kValentine=106,   kVisibleSpectrum=107,
    ## kWaterMelon=108,      kCool=109,        kCopper=110,
    ## kGistEarth=111,       kViridis=112,     kCividis=113
    
    style . SetPalette         ( ROOT.kDarkBodyRadiator ) ##  dark-body radiator pallete
    style . SetNumberContours  ( 255 )

    if force : 
        style . cd() 
        logger.debug ('The style %s is forced' % style.GetName() )
        ROOT.gROOT.SetStyle   ( style.GetName()  )
        ROOT.gROOT.ForceStyle ()
        
    return style     


# =============================================================================
## style for half-page-width COLZ plots
Style2Z = OstapStyle (
    'OstapStyle2Z' ,
    desc    = "Style, suitable to produce downscaled (2-in-row) plot with large right margin (COLZ)",
    scale   =  1.41 , makeNew = True , colz = True )

## style for one-third-page-width COLZ plots 
Style3Z = OstapStyle (
    'OstapStyle3Z' ,
    desc    = "Style, suitable to produce downscaled (3-in-row) plot with large right margin (COLZ)",
    scale   =  1.71 , makeNew = True , colz = True )

## style for standard COLZ plots 
StyleZ  = OstapStyle (
    'OstapStyleZ'  ,
    desc    = "Style, suitable to produce plot with large right margin (COLZ)",
    makeNew = True , colz = True )

## style for half-page-width plots 
Style2  = OstapStyle (
    'OstapStyle2' , 
    desc    = "Style, suitable to produce downscaled (2-in-row) plots",
    scale   =  1.41 , makeNew = True )

## style for one-third-page-width plots 
Style3  = OstapStyle (
    'LHCbStyle3' ,
    desc    = "Style, suitable to produce downscaled (3-in-row) plots",
    scale   =  1.71 , makeNew = True )


## default style 
Style      = OstapStyle()
## default style 
ostapStyle = Style

# =============================================================================
## Use some (temporary) style   as context manager 
#  @code
#  with UseStyle ( lhcbStyle2 ) :
#  ...     h1 = ...
#  ...     h1.Draw() 
#  @endcode
class UseStyle(object):
    """ Use some (temporary) style   as context manager 
    >>> with UseStyle ( lhcbStyle2 ) :
    ...     h1 = ...
    ...     h1.Draw() 
    """
    def __init__ ( self, style = None ) :
        
        if isinstance ( style , int ):
            if    1  == style : style = Style 
            elif  2  == style : style = Style2
            elif  3  == style : style = Style3
            
        if isinstance ( style , str ):
            if   style.upper() in ( '' , '0' , '1' )    : style = Style 
            elif '2'  == style                          : style = Style2
            elif '3'  == style                          : style = Style3
            elif style.upper() in ( 'Z' , '1Z' , 'Z1' ) : style = StyleZ 
            elif style.upper() in ( 'Z' , '2Z' , 'Z2' ) : style = Style2Z 
            elif style.upper() in ( 'Z' , '3Z' , 'Z3' ) : style = Style3Z 

        ## use the style by name 
        if isinstance   ( style , str ) :
            styles = ROOT.gROOT.GetListOfStyles()
            for s in styles :
                if s.GetName() == style :
                    style = s
                    break
                
        if style is None : pass  
        elif not isinstance  ( style , ROOT.TStyle ) :
            logger.warning ( 'No valid style "%s" is found, use default style' % style )
            style = ostapStyle

        self.__new_style = style
        self.__old_style = None 

    ## context  manager: enter 
    def __enter__ ( self )      :
        
        self.__force_style = ROOT.gROOT.GetForceStyle()
        if self.__new_style : 
            self.__old_style = ROOT.gStyle 
            self.__new_style.cd   ()
            ROOT.gROOT.ForceStyle ( True )
        
    ## context  manager: exit
    def __exit__  ( self , *_ ) :

        if self.__old_style: 
            self.__old_style.cd()
            ROOT.gROOT.ForceStyle ( self.__force_style ) 

        self.__new_style = None
        self.__old_style = None

    @property
    def new_style ( self ) :
        "``new_style'' : new style"
        return self.__new_style
    
    @property
    def old_style ( self ) :
        "``old_style'' : old style"
        return self.__old_style


# =============================================================================
## Use some (temporary) style   as context manager 
#  @code
#  with useStyle ( lhcbStyle2 ) :
#  ...     h1 = ...
#  ...     h1.Draw() 
#  @endcode
def useStyle ( style = Style ) :
    """ Use some (temporary) style   as context manager 
    >>> with useStyle ( Style2 ) :
    ...     h1 = ...
    ...     h1.Draw() 
    """
    return UseStyle (  style )

    
# =============================================================================
if '__main__' == __name__ :
         
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
