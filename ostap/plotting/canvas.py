#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# =============================================================================
## @file canvas.py
#  
#     .oooooo.                .                        
#    d8P'  `Y8b             .o8                        
#   888      888  .oooo.o .o888oo  .oooo.   oo.ooooo.  
#   888      888 d88(  "8   888   `P  )88b   888' `88b 
#   888      888 `"Y88b.    888    .oP"888   888   888 
#   `88b    d88' o.  )88b   888 . d8(  888   888   888 
#    `Y8bood8P'  8""888P'   "888" `Y888""8o  888bod8P' 
#                                            888       
#                                           o888o      
#                                                    
#  Simple helper module to get ROOT TCanvas
# 
#  This file is a part of 
#  <a href="http://cern.ch/lhcb-comp/Analysis/Bender/index.html">Bender project</a>
#  <b>``Python-based Interactive Environment for Smart and Friendly Physics Analysis''</b>
#
#  The package has been designed with the kind help from
#  Pere MATO and Andrey TSAREGORODTSEV. 
#  And it is based on the 
#  <a href="http://cern.ch/lhcb-comp/Analysis/LoKi/index.html">LoKi project:</a>
#  <b>``C++ ToolKit for Smart and Friendly Physics Analysis''</b>
#
#  By usage of this code one clearly states the disagreement 
#  with the smear campaign of Dr.O.Callot et al.: 
#  ``No Vanya's lines are allowed in LHCb/Gaudi software''
#
#  @date   2012-02-15
#  @author Vanya BELYAEV Ivan.Belyaevitep.ru
#
# =============================================================================
""" Simple helper module to get ROOT TCanvas
    
    This file is a part of BENDER project:

  ``Python-based Interactive Environment for Smart and Friendly Physics Analysis''

The project has been designed with the kind help from Pere MATO and Andrey TSAREGORODTSEV. 

And it is based on the LoKi project:
 
   ``C++ ToolKit for Smart and Friendly Physics Analysis''

By usage of this code one clearly states the disagreement with the smear campaign 
of Dr.O.Callot et al.:

   ``No Vanya's lines are allowed in LHCb/Gaudi software''
"""
# =============================================================================
__author__  = 'Vanya BELYAEV Ivan.Belyaev@itep.ru'
__date__    = "2014-10-19"
__version__ = '$Revision$'
__all__     = (
    'getCanvas'   ,
    'getCanvases' ,
    'AutoPlots'   , ## context manager to activate the auto-plotting machinery
    'auto_plots'  , ## ditto, but as function 
    )
# =============================================================================
import ROOT, os, tempfile  
import ostap.core.core    
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.plotting.canvas' )
else                       : logger = getLogger( __name__ )
# =============================================================================
_canvases = {} 
# =============================================================================
## get the canvas
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-10-19
def getCanvas ( name   = 'glCanvas' ,   ## canvas name 
                title  = 'Ostap'    ,   ## canvas title
                width  = 1000       ,   ## canvas width
                height = 800        ) : ## canvas height 
    """Get create canvas/create new canvas
    
    >>> cnv = getCanvas ( 'glnewCanvas' , width = 1200 , height = 1000 )
    """
    cnv   = _canvases.get ( name , None )
    if not cnv :
        ## create new canvas 
        ## cnv  = ROOT.TCanvas ( 'glCanvas', 'Ostap' , width , height )
        cnv  = ROOT.TCanvas ( name , 'Ostap' , width , height )
        ## adjust newly created canvas
        ## @see http://root.cern.ch/root/html/TCanvas.html#TCanvas:TCanvas@4
        if not ROOT.gROOT.IsBatch() :
            dw = width  - cnv.GetWw()
            dh = height - cnv.GetWh()
            cnv.SetWindowSize ( width + dw , height + dh )
            
        ## 
        _canvases [ name ] = cnv
        
    return cnv

# =============================================================================
all_extensions = (
    'pdf'  , 'png'  , 'gif' ,
    'eps'  , 'ps'   ,
    'cxx'  , 'c'    , 
    'jpg'  , 'jpeg' , 'svg' , 
    'root' , 'xml'  , 'xpm' , 
    'tiff' , 'tex' 
    )

# =============================================================================
## define simplified print for TCanvas 
def _cnv_print_ ( cnv , fname , exts = ( 'pdf' , 'png' , 'eps', 'C' ) ) :
    """A bit simplified version for TCanvas print
    >>> canvas.print ( 'fig' )    
    """
    #
    cnv.Update () 
    from ostap.logger.utils import rootWarning 
    n,d,e = fname.rpartition('.')
    if d and e.lower() in all_extensions : 
        with rootWarning () :
            cnv.Update   () 
            cnv.Print    ( fname )
            logger.debug ( 'Canvas --> %s' % fname )
            return cnv
        
    for ext in exts :
        with rootWarning () :
            name = fname + '.' + ext
            cnv.Print   ( name )
            logger.debug('Canvas --> %s' % name )
            
    return cnv 

# =============================================================================
## define streamer for canvas
#  @code
#  canvas >> 'a'    
#  @endcode 
def _cnv_rshift_ ( cnv , fname ) :
    """Very simple print for canvas:
    >>> canvas >> 'a'    
    """
    return _cnv_print_ ( cnv , fname )

ROOT.TVirtualPad.print_     = _cnv_print_
ROOT.TVirtualPad.__rshift__ = _cnv_rshift_

# =============================================================================
# Auto-plotting
# =============================================================================
from collections import defaultdict 
# =============================================================================
## @class AutoPlots
#  Helper structure/context manager to setup "auto-plotting"
#  all produced plots will be saved
class AutoPlots ( object ) :
    """helper structure to setup ``auto-plotting''
    all produced plots will be saved
    """
    
    __auto_plots  = []
    __plot_counts = defaultdict(int) 
    
    @classmethod
    def plot ( cls ) :
        """Get the plot name for auto-plotting:
        >>> plotname = AutoPlot.plot()
        >>> if plotname : canvas >> plotname 
        """
        if not cls.__auto_plots : return False 
        p  = cls.__auto_plot[-1]
        c  = cls._counts[p] + 1
        c %=  1000
        cls._counts[p] = c 
        return p % c 
    
    def __init__ ( self                      ,
                   pattern   = 'ostap_%0.4d' ,
                   directory = ''            ) :

        if directory and not os.path.exists ( directory ) :
            try :
                os.path.makedirs ( _dname )
            except OSError :
                logger.error ( "Can't create directory \"%s\"" %  directory)
                directory = None
                
        if directory and os.path.exists ( directory ) :
            if not os.path.isdir ( directory ) :
                logger.error ( "Invalid directory \"%s\"" %  directory)
                directory = None
                
        if directory and os.path.exists ( directory ) and os.path.isdir ( directory ) : 
            if not os.access ( directory , os.W_OK ) :
                logger.error ( "Non-writeable directory \"%s\"" %  directory)
                directory = None
                
        if not directory :
            directory = tempfile.mkdtemp ( prefix = 'plots_' )
            logger.info ( 'AutoPlots: use directory "%s"' % directory ) 

        ## check the validity of pattern
        try :
            r = pattern % 999
        except TypeError :
            pattern = 'ostap_%0.4d' 
            logger.info ( 'Pattern "%s"' % pattern)

        self.__pattern = os.path.join ( directory , pattern )

    @property
    def pattern ( self ):
        """``pattern'' -  the pattern for the file name"""
        return self.__pattern

    ## context manager  
    def __enter__ ( self ) :        
        AutoPlots.__auto_plots.append ( self.pattern )
        return self 
        
    def __exit__  ( self, *_ ) :
        AutoPlots.__auto_plots.pop  ()  


# =============================================================================
## Helper function /context manager to setup "auto-plotting"
#  all produced plots will be saved
#  @code
#  with auto_plot ( 'all_%d'  , directory  = 'plots' ) :
#  ...     a.draw()
#  ...     b.draw()
#  ...     c.draw()
#  ...     d.draw()
#  @endcode
def auto_plot ( pattern   = 'ostap_%0.4d' ,
                directory = ''            ) :
    """Helper function /context manager to setup "auto-plotting"
    all produced plots will be saved
    with auto_plots ( 'all_%d'  , directory  = 'plots' ) :
    ...     a.draw()
    ...     b.draw()
    ...     c.draw()
    ...     d.draw()
    """
    return AutoPlots ( pattern = pattern , directory = directory )

# =============================================================================
##  new draw method: silent draw
def _TO_draw_ ( obj , *args , **kwargs ) :
    """ (silent) Draw of ROOT object
    >>> obj
    >>> obj.Draw()  ##
    >>> obj.draw()  ## ditto
    """
    from ostap.logger.utils import rootWarning, rooSilent 
    with rootWarning() , rooSilent ( 2 ) :
        result = obj.Draw ( *args , **kwargs )
        if ROOT.gPad : 
            plot = AutoPlots.plot()
            if plot : ROOT.gPad >> plot
        return result

# =============================================================================
## decorate ROOT.TObject
if not hasattr ( ROOT.TObject , 'draw_with_autoplot' ) :
    
    ## add new method  
    ROOT.TObject.draw_with_autoplot = _TO_draw_
    ## save old method 
    if hasattr ( ROOT.TObject ,  'draw' ) :
        ROOT.TObject._draw_backup = ROOT.TObject.draw
        
    ROOT.TObject.draw       = ROOT.TObject.draw_with_autoplot
        
# =============================================================================
## get all known canvases 
def getCanvases () :
    """ Get all known canvases """ 
    return _canvases.keys() 

def _remove_canvases_() :
    keys = _canvases.keys() 
    for k in keys : del _canvases[k]
        
import atexit
atexit.register ( _remove_canvases_ )

# =============================================================================
_decorated_classes_  = (
    ROOT.TVirtualPad , 
    ROOT.TObject     , 
    )

_new_methods_        = (
    ROOT.TVirtualPad.__rshift__ , 
    ROOT.TObject    .draw       , 
    ROOT.TObject    .draw_with_autoplot, 
    )
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
