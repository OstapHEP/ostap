#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# =============================================================================
## @file canvas.py
#  Simple helper module to get ROOT TCanvas
#  @date   2012-02-15
#  @author Vanya BELYAEV Ivan.Belyaevitep.ru
# =============================================================================
""" Simple helper module to get ROOT TCanvas
"""
# =============================================================================
__author__  = 'Vanya BELYAEV Ivan.Belyaev@itep.ru'
__date__    = "2014-10-19"
__version__ = '$Revision$'
__all__     = (
    'getCanvas'        , ## get/create canvas 
    'getCanvases'      , ## get all created canvases 
    'canvas_partition' , ## split canvas into several pads with no space between pads 
    'canvas_pull'      , ## split canvas into two pads with no vertical interspace
    'draw_pads'        , ## plot sequence of object on sequence of pads, adjustinng axis label size
    'AutoPlots'        , ## context manager to activate the auto-plotting machinery
    'auto_plots'       , ## ditto, but as function 
    )
# =============================================================================
import ROOT, os, tempfile  
import ostap.core.core
import ostap.plotting.style
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.plotting.canvas' )
else                       : logger = getLogger( __name__ )
# =============================================================================
_canvases = {} 
# =============================================================================
## The default configuration 
# =============================================================================
import ostap.core.config as OCC
width_  = OCC.canvas.get ( 'Width'  , fallback = '1000' )
try               : width_  = int  ( width_ ) 
except ValueError : width_  = 1000 
height_ = OCC.canvas.get ( 'Height' , fallback =  '800' )
try               : height_ = int  ( height_ ) 
except ValueError : height_ = 800 

# =============================================================================
## get the canvas
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-10-19
def getCanvas ( name   = 'glCanvas' ,   ## canvas name 
                title  = 'Ostap'    ,   ## canvas title
                width  = width_     ,   ## canvas width
                height = height_    ) : ## canvas height 
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
def auto_plots ( pattern   = 'ostap_%0.4d' ,
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
        pad = ROOT.gROOT.GetSelectedPad()
        if pad and not ROOT.gPad :
            c = pad.GetCanvas()
            if c : c.Update()
        if ROOT.gPad :
            c = ROOT.gPad
            if c : c.Update() 
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
    keys = list( _canvases.keys() )
    for k in keys : del _canvases[k]
        
import atexit
atexit.register ( _remove_canvases_ )

# =============================================================================

# =============================================================================
## perform partition of Canvas into
#  @code
#  canvas = ...
#  pads   = canvas.partition ( 3 , 2 )
#  for i in range(3) :
#    for j in range(2) :
#       histo_ij = ...
#       canvas.cd(0)
#       pad = pads[i,j]
#       pad.Draw()
#       pad.cd()
#       histo_ij.Draw()
#  @endcode
#  @see https://root.cern/doc/master/canvas2_8C.html
def canvas_partition ( canvas               , 
                       nx                   ,
                       ny                   ,
                       left_margin   = 0.10 ,
                       right_margin  = 0.03 ,   
                       bottom_margin = 0.10 ,
                       top_margin    = 0.03 ,
                       hSpacing      = 0.0  ,
                       vSpacing      = 0.0  ) :
    """Perform partition of Canvas into pads with no inter-margins

    canvas = ...
    nx = 3 , ny = 2 
    pads   = canvas.partition ( nx  , ny )
    for i in range(nx) :
    ... for j in range(ny) :
    ... ... histo_ij = ...
    ... ... canvas.cd(0)
    ... ... pad = pads[i,j]
    ... ... pad.Draw()
    ... ... pad.cd()
    ... ... histo_ij.Draw()
    
    @see https://root.cern/doc/master/canvas2_8C.html
    
    """

    if not isinstance ( nx , int ) or nx<= 0 :
        raise AttributeError('partition: invalid nx=%s' % nx )
    if not isinstance ( ny , int ) or ny<= 0 :
        raise AttributeError('partition: invalid ny=%s' % ny )

    ## get the window size
    wsx = abs ( canvas.GetWindowWidth  () ) 
    wsy = abs ( canvas.GetWindowHeight () ) 

    #
    ## if parametes given in the absolute units, convert them into relative coordinates
    #
    
    if   left_margin < 0 :   left_margin = abs (   left_margin ) / wsx
    if  right_margin < 0 :  right_margin = abs (  right_margin ) / wsx
    if bottom_margin < 0 : bottom_margin = abs ( bottom_margin ) / wsy
    if    top_margin < 0 :    top_margin = abs (    top_margin ) / wsy
    
    if hSpacing      < 0 : hSpacing = abs ( hSpacing ) / wsx
    if vSpacing      < 0 : vSpacing = abs ( vSpacing ) / wsy

    #
    ## check consistency 
    # 
    if 1 <=   left_margin :
        raise AttributeError('partition: invalid   left margin=%f' %   left_margin )
    if 1 <=  right_margin :
        raise AttributeError('partition: invalid  right margin=%f' %  right_margin )
    if 1 <= bottom_margin :
        raise AttributeError('partition: invalid bottom margin=%f' % bottom_margin )
    if 1 <=    top_margin :
        raise AttributeError('partition: invalid    top margin=%f' %    top_margin )
    
    if   hasattr ( canvas , 'pads' ) and isinstance ( canvas.pads , dict ) :
        while canvas.pads :
            i,p = canvas.pads.popitem()
            if p :
                logger.verbose ( 'delete pas %s' % p.GetName() )
                del p
        del canvas.pads
        
    elif hasattr ( canvas , 'pads' ) and isinstance ( canvas.pads , tuple ) :
        for p in canvas.pads :
            if p :
                logger.verbose ( 'delete pas %s' % p.GetName() )
                del p
        del canvas.pads
        
    ## make new empty dictionary 
    canvas.pads = {} 
            
    vStep    = ( 1.0 - bottom_margin - top_margin   - (ny-1) * vSpacing ) / ny
    if 0 > vStep : raise AttributeError('partition: v-step=%f' % vStep  )
        
    hStep    = ( 1.0 - left_margin   - right_margin - (nx-1) * vSpacing ) / nx 
    if 0 > hStep : raise AttributeError('partition: h-step=%f' % hStep  )

    hposr, hposl, hmarr, hmarl, hfactor = 0.,0.,0.,0.,0.
    vposr, vposd, vmard, vmaru, vfactor = 0.,0.,0.,0.,0.
    
    for ix in range ( nx ) :
        
        if 0 == ix : 
            hposl   = 0
            hposr   = left_margin + hStep
            hfactor = hposr - hposl
            hmarl   = left_margin / hfactor
            hmarr   = 0.0 
        elif nx == ix + 1 :
            hposl   = hposr + hSpacing 
            hposr   = hposl + hStep + right_margin
            hfactor = hposr - hposl 
            hmarl   = 0.0
            hmarr   = right_margin / hfactor 
        else : 
            hposl   = hposr + hSpacing
            hposr   = hposl + hStep
            hfactor = hposr - hposl
            hmarl   = 0.0
            hmarr   = 0.0

        for iy in range(ny) :
            if 0 == iy : 
                vposd   = 0.0
                vposu   = bottom_margin + vStep
                vfactor = vposu - vposd
                vmard   = bottom_margin / vfactor
                vmaru   = 0.0 
            elif ny == iy + 1 : 
                vposd   = vposu + vSpacing
                vposu   = vposd + vStep + top_margin
                vfactor = vposu - vposd;
                vmard   = 0.0
                vmaru   = top_margin    / vfactor 
            else :
                vposd   = vposu + vSpacing
                vposu   = vposd + vStep
                vfactor = vposu - vposd
                vmard   = 0.0
                vmaru   = 0.0

            canvas.cd(0)
            pname = 'pad_%s_%d_%d' % ( canvas.GetName() , ix , iy )
            pad   = ROOT.gROOT.FindObject ( pname )
            if pad : del pad
            pad   = ROOT.TPad ( pname , '' ,  hposl , vposd  , hposr , vposu )

            logger.verbose ( ' Create pad[%d,%d]=(%f,%f,%f,%f),[%f,%f,%f,%f] %s ' % (
                ix    , iy    ,
                hposl , vposd , hposr , vposu , 
                hmarl , hmarr , vmard , vmaru , pad.GetName() ) ) 
                             
            pad.SetLeftMargin      ( hmarl )
            pad.SetRightMargin     ( hmarr )
            pad.SetBottomMargin    ( vmard )
            pad.SetTopMargin       ( vmaru )
            
            pad.SetFrameBorderMode ( 0 )
            pad.SetBorderMode      ( 0 )
            pad.SetBorderSize      ( 0 )

            ROOT.SetOwnership ( pad , True )
            
            if not hasattr ( canvas , 'pads' ) : canvas.pads = {}
            canvas.pads[  ( ix , iy ) ] = pad

    pds  = canvas.pads
    keys = pds.keys()
    import collections as _C
    _p  =_C.OrderedDict ()
    for k in  sorted ( keys ) :
        _p[k] = pds[k]
    canvas.pads = _p 
    
    return canvas.pads 


ROOT.TCanvas.partition = canvas_partition


# ==============================================================================
## Perform partition of Canvas into 1x2 non-equal pads with no inter-margins
#  @code
#  canvas    = ...
#  pad_u, pud_b= canvas.pull_partition ( 0.20 )    
#  canvas.cd(0)
#  pad_u.Draw()
#  pad_u.cd() 
#  object1.Draw()    
#  canvas.cd(0)
#  pad_b.Draw()
#  pad_b.cd() 
#  object2.Draw()
#  @endcode 
def canvas_pull ( canvas               ,
                  ratio         = 0.80 ,
                  left_margin   = 0.14 ,
                  right_margin  = 0.05 ,   
                  bottom_margin = 0.14 ,
                  top_margin    = 0.05 ,
                  hSpacing      = 0.0  ,
                  vSpacing      = 0.0  ) :
    """Perform partition of Canvas into 1x2 non-equal pads with no inter-margins
    
    >>> canvas    = ...
    >>> pad_u, pud_b= canvas.pull_partition ( 0.20 )
    
    >>> canvas.cd(0)
    >>> pad_u.Draw()
    >>> pad_u.cd() 
    >>> object1.Draw()
    
    >>> canvas.cd(0)
    >>> pad_b.Draw()
    >>> pad_b.cd() 
    >>> object2.Draw()
    """

    ## get the window size
    wsx = abs ( canvas.GetWindowWidth  () ) 
    wsy = abs ( canvas.GetWindowHeight () ) 

    #
    ## if parametes given in the absolute units, convert them into relative coordinates
    #
    
    if   left_margin < 0 :   left_margin = abs (   left_margin ) / wsx
    if  right_margin < 0 :  right_margin = abs (  right_margin ) / wsx
    if bottom_margin < 0 : bottom_margin = abs ( bottom_margin ) / wsy
    if    top_margin < 0 :    top_margin = abs ( bottom_margin ) / wsy
    
    if hSpacing      < 0 : hSpacing = abs ( hSpacing ) / wsx
    if vSpacing      < 0 : vSpacing = abs ( vSpacing ) / wsy
    
    hposr, hposl, hmarr, hmarl, hfactor = 0.,0.,0.,0.,0.
    vposr, vposd, vmard, vmaru, vfactor = 0.,0.,0.,0.,0.

    nx = 1
    ny = 2
    
    vStep    = ( 1.0 - bottom_margin - top_margin   - ( ny - 1 ) * vSpacing ) / ny
    if 0 > vStep : raise AttributeError('partition: v-step=%f' % vStep  )
    
    hStep    = ( 1.0 - left_margin   - right_margin - ( nx - 1 ) * vSpacing ) / nx 
    if 0 > hStep : raise AttributeError('partition: h-step=%f' % hStep  )
    
    hposl   = 0
    hposr   = left_margin + hStep
    hfactor = hposr - hposl
    hmarl   = left_margin / hfactor
    hmarr   = 0.0
    
    if hasattr ( canvas ,  'pads' ) :
        del canvas.pads


    vStep0 = 2 * vStep * 1     / ( 1 + ratio )
    vStep1 = 2 * vStep * ratio / ( 1 + ratio )
    
    ix = 0 
    for iy in range(2) :
        
        if 0 == iy : 
            vposd   = 0.0
            vposu   = bottom_margin + vStep0
            vfactor = vposu - vposd 
            vmard   = bottom_margin / vfactor
            vmaru   = 0.0 
        else : 
            vposd   = vposu + vSpacing
            vposu   = vposd + vStep1 + top_margin
            vfactor = vposu - vposd
            vmard   = 0.0
            vmaru   = top_margin    / vfactor
            
        canvas.cd(0)
        pname = 'pad_%s_%d_%d' % ( canvas.GetName() , ix , iy )
        pad   = ROOT.gROOT.FindObject ( pname )
        if pad : del pad
        pad   = ROOT.TPad ( pname , '' ,  hposl , vposd  , hposr , vposu )
        
        logger.verbose ( ' Create pad[%d,%d]=(%f,%f,%f,%f),[%f,%f,%f,%f] %s ' % (
            ix    , iy    ,
            hposl , vposd , hposr , vposu , 
            hmarl , hmarr , vmard , vmaru , pad.GetName() ) ) 
        
        pad.SetLeftMargin      ( hmarl )
        pad.SetRightMargin     ( hmarr )
        pad.SetBottomMargin    ( vmard )
        pad.SetTopMargin       ( vmaru )
        
        pad.SetFrameBorderMode ( 0 )
        pad.SetBorderMode      ( 0 )
        pad.SetBorderSize      ( 0 )
        
        ROOT.SetOwnership ( pad , True )
        
        if not hasattr ( canvas , 'pads' ) : canvas.pads=[]
        canvas.pads.append ( pad ) 
        
    pads =  tuple ( reversed ( canvas.pads ) )
    
    import collections as _C
    _pads =_C.OrderedDict ()
    for k, p in enumerate ( pads ) :
        _pads [k] = p
    canvas.pads = _pads 

    return canvas.pads 


ROOT.TCanvas.pull_partition = canvas_pull

# =============================================================================
## draw sequence of object on sequence of pads,
#  - the label font size is adjusted to be uniform (in pixels) 
#  @code
#  pads   = ...
#  frames = ...
#  draw_pads ( frames , pads , fontsize = 25 ) 
#  @endcode 
def draw_pads ( objects , pads , fontsize = 25 ) :
    """Draw sequence of object on sequence of pads,
    - the label font size is adjusted to be uniform (in pixels)     
    >>> pads   = ...
    >>> frames = ...
    >>> draw_pads ( frames , pads , fontsize = 25 ) 
    """

    assert isinstance  ( fontsize , int ) and 1 <= fontsize , 'Invalid fontsize %s [pixels] ' % fontsize
    
    for obj , pad_ in zip ( objects , pads ) : 
        
        if isinstance ( pad_ , ROOT.TPad ) : pad = pad_
        else                               : pad = pads [ pad_ ] 
        
        c = pad.GetCanvas()
        if c : c.cd(0)
        
        pad.draw ()
        pad.cd   ()

        ## redefne label font and size 
        for attr in ( 'GetXaxis' , 'GetYaxis' , 'GetZaxis' ) :
            if not hasattr ( obj , attr ) : continue
            
            axis = getattr ( obj , attr )()
            if not axis : continue
            
            fnp  = axis.GetLabelFont ()
            fn , prec = divmod  ( fnp , 10 ) 
            if 3 != prec :
                ## redefine  label  font 
                fnp = fn * 10 + 3
                axis.SetLabelFont ( fnp )

            ## redefine label size 
            axis.SetLabelSize ( fontsize  )

        ## draw object on the pad 
        obj.draw ()
        
        if c : c.cd(0)
        
        
# =============================================================================
_decorated_classes_  = (
    ROOT.TVirtualPad , 
    ROOT.TCanvas     , 
    ROOT.TObject     , 
    )

_new_methods_        = (
    ROOT.TVirtualPad.__rshift__ ,
    ROOT.TCanvas.partition      ,
    ROOT.TCanvas.pull_partition ,
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
