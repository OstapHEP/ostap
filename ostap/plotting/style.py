#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/plotting/style.py
# Ostap style file for ROOT-plots
# =============================================================================
"""Ostap Style for ROOT-plots"""
# =============================================================================
__all__ = (
    'UseStyle'         ,  ## context manager  for the style (class) 
    'useStyle'         ,  ## context manager  for the style (function)
    ##
    'Style'            ,  ## the default Ostap style 
    'ostap_label'      ,  ## the style for the text 
    'ostap_latex'      ,  ## the style for LaTeX 
    'ostap_font'       ,  ## the default  font 
    'ostap_line_width' ,  ## the default  font 
    'OstapStyle'       ,  ## the function to instantiate the style
    ## 
    'StyleZ'           ,  ## the style for COLZ plots
    'Style1'           ,  ## the default Ostap style 
    'Style2'           ,  ## the style for downscaled 2-in-row plots
    'Style3'           ,  ## the style for downscaled 3-in-row plots
    'Style1Z'          ,  ## the style for COLZ plots
    'Style2Z'          ,  ## the style for downscaled 2-in-row COLZ plots
    'Style3Z'          ,  ## the style for downscaled 3-in-row COLZ plots
    )
# =============================================================================
import ostap.plotting.color 
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.plotting.style' )
else                       : logger = getLogger( __name__ )
# =============================================================================
from ostap.plotting.makestyles import ( ostap_font       ,
                                        ostap_label      ,
                                        ostap_line_width ,
                                        ostap_latex      ,
                                        set_style        )
# =============================================================================
## the dictionary of known  styles 
styles = {}

# =============================================================================
## Create Ostap-style for the plots
def OstapStyle ( name                           ,
                 description                    ,
                 base        = ''               , 
                 line_width  = ostap_line_width ,
                 font        = ostap_font       ,                 
                 force       = True             ,
                 scale       = 1.0              ,
                 colz        = False            ) :
    """Create Ostap-style for the plots    
    """

    # ================================================================
    ## check the configuration 
    import ostap.core.config         as CONFIG

    config = {} 
    n2     = name.strip()

    base_conf = {}
    if base :
        groot = ROOT.ROOT.GetROOT      ()
        slst  = groot.GetListOfStyles  ()
        for s in slst :
            if base == s.GetName() :
                base_config = dump_style ( s ) 
                break
        else :
            logger.warning ( "OstapStyle: no base style `%s' is found! " % base )
            
    ## Look for "[Style:Nick]" in configuration 
    for key in CONFIG.config :
        lkey = key.lower() 
        if not lkey.startswith('style') : continue
        s , c , n = key.partition(':')
        if not c : continue
        n1 = n.strip() if c else s.strip()
        if n1 == n2 :
            config = CONFIG.config[key]
            logger.debug ( 'Use existing configuration style %s' % name )
            break
        
    import ostap.plotting.makestyles as MS 
    style = MS.make_ostap_style ( name        = name        ,
                                  description = description ,
                                  base        = base_conf   , 
                                  config      = config      ,
                                  colz        = colz        ,
                                  scale       = scale       ,
                                  font        = font        ,
                                  line_width  = line_width  )
    # =======================================================================
    if force : 
        style . cd() 
        logger.debug ('The style %s is forced' % style.GetName() )
        groot = ROOT.ROOT.GetROOT() 
        groot.SetStyle   ( style.GetName()  )
        groot.ForceStyle ()
        
    return style     

# =============================================================================
## style for standard COLZ plots 
Style1Z = OstapStyle (
    'Ostap1Z'   ,
    description = "Ostap style, suitable to produce plots with large right margin (COLZ)",
    colz        = True )
## ditto 
StyleZ  = Style1Z 
# =============================================================================
## style for half-page-width COLZ plots
Style2Z = OstapStyle (
    'Ostap2Z'   ,
    description = "Style, suitable to produce downscaled (2-in-row) plot with large right margin (COLZ)",
    scale       = 1.41 ,
    colz        = True )
# =============================================================================
## style for one-third-page-width COLZ plots 
Style3Z = OstapStyle (
    'Ostap3Z'   ,
    description = "Style, suitable to produce downscaled (3-in-row) plot with large right margin (COLZ)",
    scale       = 1.71 ,
    colz        = True )
# ===========================================================================
## style for half-page-width plots 
Style2  = OstapStyle (
    'Ostap2'    , 
    description = "Ostap style, suitable to produce downscaled (2-in-row) plots",
    scale       =  1.41 )
# ============================================================================
## style for one-third-page-width plots 
Style3  = OstapStyle (
    'Ostap3'    ,
    description = "Ostap style, suitable to produce downscaled (3-in-row) plots",
    scale       =  1.71 )
# ============================================================================
## default Ostap style 
Style = OstapStyle (
    'Ostap'     ,
    description = 'The basic/default Ostap style for plots' ,
    force       = True 
)
## ditto
Style1 = Style 
# ============================================================================
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
    def __init__ ( self, style = None , **config ) :
        
        if   style is None : style = ostapStyle 
        elif isinstance ( style , int ):
            
            if    0  == style : style = ostapStyle
            elif  1  == style : style = Style1
            elif  2  == style : style = Style2
            elif  3  == style : style = Style3
            
        elif isinstance ( style , str ):
            
            import ostap.plotting.makestyles as MS            
            if   style in MS.StyleStore.styles : style = MS.StyleStore.styles [ style ] 
            elif style.upper() in ( '' , '0' , '1' , 'def' , 'default') : style = ostapStyle
            elif '2'  == style                          : style = Style2
            elif '3'  == style                          : style = Style3
            elif style.upper() in ( 'Z' , '1Z' , 'Z1' ) : style = StyleZ 
            elif style.upper() in (       '2Z' , 'Z2' ) : style = Style2Z 
            elif style.upper() in (       '3Z' , 'Z3' ) : style = Style3Z 
            else :
                groot = ROOT.ROOT.GetROOT ()        
                slst  = groot.GetListOfStyles()
                for s in slst :
                    if s and s.GetName () == style :
                        style = s
                        break
                    
        if ( not style ) or ( not isinstance  ( style , ROOT.TStyle ) ) :
            unknown = style 
            style   = ostapStyle
            logger.warning ( "No valid style `%s'is found, use default `%s' style" % ( str ( unknown ) , style.GetName() ) ) 

        self.__new_style = style
        self.__old_style = None 

        self.__config    = config 
        self.__changed   = {}
        
    ## context  manager: enter 
    def __enter__ ( self )      :
        
        groot  = ROOT.ROOT.GetROOT()        
        self.__force_style = groot.GetForceStyle()
        if self.__new_style : 
            self.__old_style = ROOT.gStyle 
            self.__new_style.cd   ()
            if self.__config :
                self.__changed = set_style ( self.__new_style , self.__config ) 
            groot.ForceStyle ( True )
            if ROOT.gPad : ROOT.gPad.UseCurrentStyle()
            logger.debug ( "Swith to `%s' style!" % self.__new_style.GetName() ) 
                           
    ## context  manager: exit
    def __exit__  ( self , *_ ) :

        if self.__changed :
            set_style ( self.__new_style , self.__changed ) 
            
        if self.__old_style : 
            self.__old_style.cd()
            groot = ROOT.ROOT.GetROOT()        
            groot.ForceStyle ( self.__force_style ) 
            logger.debug ( "Swith to `%s' style!" % self.__old_style.GetName() ) 

        self.__new_style = None
        self.__old_style = None

    @property
    def new_style ( self ) :
        "`new_style' : new style"
        return self.__new_style
    
    @property
    def old_style ( self ) :
        "`old_style' : old style"
        return self.__old_style

    @property
    def config   ( self ) :
        """`config' : addtional configuration parameters """
        return self.__config

    @property
    def changed   ( self ) :
        """`changed' : changed configuration parameters """
        return self.__changed
    
# =============================================================================
## Use some (temporary) style   as context manager 
#  @code
#  with useStyle ( lhcbStyle2 ) :
#  ...     h1 = ...
#  ...     h1.Draw() 
#  @endcode
def useStyle ( style = None , **config  ) :
    """ Use some (temporary) style   as context manager 
    >>> with useStyle ( Style2 ) :
    ...     h1 = ...
    ...     h1.Draw() 
    """
    return UseStyle ( style , **config )

    
# =============================================================================
if '__main__' == __name__ :
         
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
