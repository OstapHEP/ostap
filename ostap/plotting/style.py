#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/plotting/style.py
# Ostap style file for ROOT-plots
# =============================================================================
"""Ostap Style for ROOT-plots"""
# =============================================================================
import ROOT
import ostap.plotting.color 
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
    'Style2'           ,  ## the style for downscaled 2-in-row plots
    'Style3'           ,  ## the style for downscaled 3-in-row plots
    'Style2Z'          ,  ## the style for downscaled 2-in-row COLZ plots
    'Style3Z'          ,  ## the style for downscaled 3-in-row COLZ plots
    )
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
                 line_width  = ostap_line_width , 
                 font        = ostap_font       ,
                 makeNew     = False            ,
                 force       = True             ,
                 scale       = 1.0              ,
                 colz        = False            ) :
    """Create Ostap-style for the plots    
    """
    groot = ROOT.ROOT.GetROOT() 
    obj   = groot.FindObject  ( name )
    if obj and isinstance ( obj , ROOT.TStyle ) and not makeNew : 
        logger.info               ('The style %s is reused' % obj.GetName() )
        if force :
            obj.cd                () 
            logger.info           ('The style %s is forced' % obj.GetName() )
            groot.SetStyle   ( obj.GetName()  )
            groot.ForceStyle ( )
        return obj

    nam = name
    i   = 1
    while obj :
        nam  = name + '_%d' % i
        obj  = groot.FindObject ( nam )
        i   += 1

    # ================================================================
    ## check the configuration 
    import ostap.core.config         as CONFIG

    config = {} 
    for key in CONFIG.config :
        if not key.upper().startswith('STYLE') : continue
        s , c , n = key.partition(':')
        if not c : continue
        n1 = n.strip()
        n2 = name.strip()
        if n1 == n2 :
            config = CONFIG.config[key]
            logger.info ( 'Use existing configuration style %s' % name )
            break
        
    import ostap.plotting.makestyles as MS 
    style = MS.make_ostap_style ( name        = nam         ,
                                  description = description , 
                                  config      = config      ,
                                  colz        = colz        ,
                                  scale       = scale       ,
                                  font        = font        ,
                                  line_width  = line_width  )
        
    if force : 
        style . cd() 
        logger.debug ('The style %s is forced' % style.GetName() )
        groot.SetStyle   ( style.GetName()  )
        groot.ForceStyle ()
        
    return style     

# =============================================================================
## style for half-page-width COLZ plots
Style2Z = OstapStyle (
    'Style2Z' ,
    description = "Style, suitable to produce downscaled (2-in-row) plot with large right margin (COLZ)",
    scale       =  1.41 , makeNew = True , colz = True )

## style for one-third-page-width COLZ plots 
Style3Z = OstapStyle (
    'Style3Z' ,
    description = "Style, suitable to produce downscaled (3-in-row) plot with large right margin (COLZ)",
    scale       =  1.71 , makeNew = True , colz = True )

## style for standard COLZ plots 
StyleZ  = OstapStyle (
    'Style1Z' ,
    description = "Style, suitable to produce plot with large right margin (COLZ)",
    makeNew     = True , colz = True )

## style for half-page-width plots 
Style2  = OstapStyle (
    'Style2' , 
    description = "Style, suitable to produce downscaled (2-in-row) plots",
    scale       =  1.41 , makeNew = True )

## style for one-third-page-width plots 
Style3  = OstapStyle (
    'Style3' ,
    description = "Style, suitable to produce downscaled (3-in-row) plots",
    scale       =  1.71 , makeNew = True )

# ============================================================================
## default style 
Style      = OstapStyle (
    'Style' ,
    'The basic Ostap style for plots' )

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
        
        if   isinstance ( style , int ):
            
            if    1  == style : style = Style 
            elif  2  == style : style = Style2
            elif  3  == style : style = Style3
            
        if isinstance ( style , str ):
            
            import ostap.plotting.makestyles as MS
            if   style in MS.StyleStore.styles : style = MC.StyleStore.styles [ style ] 
            elif style.upper() in ( '' , '0' , '1' )    : style = Style 
            elif '2'  == style                          : style = Style2
            elif '3'  == style                          : style = Style3
            elif style.upper() in ( 'Z' , '1Z' , 'Z1' ) : style = StyleZ 
            elif style.upper() in (       '2Z' , 'Z2' ) : style = Style2Z 
            elif style.upper() in (       '3Z' , 'Z3' ) : style = Style3Z 

        ## use the style by name 
        if isinstance   ( style , str ) :
            groot  = ROOT.ROOT.GetROOT()
            styles = groot.GetListOfStyles()
            for s in styles :
                if s.GetName() == style :
                    style = s
                    break
                
        if   style is None : style = ROOT.gStyle 
        elif not isinstance  ( style , ROOT.TStyle ) :
            logger.warning ( 'No valid style "%s" is found, use default style' % style )
            style = ostapStyle

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
            
    ## context  manager: exit
    def __exit__  ( self , *_ ) :

        if self.__changed :
            self.__changed = set_style ( self.__new_style , self.__changed ) 
            
        if self.__old_style: 
            self.__old_style.cd()
            groot = ROOT.ROOT.GetROOT()        
            groot.ForceStyle ( self.__force_style ) 

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

    @property
    def config   ( self ) :
        """``config'' : addtional configuration parameters """
        return self.__config

    @property
    def changed   ( self ) :
        """``changed'' : changed configuration parameters """
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
