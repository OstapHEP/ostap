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
    'getCanvases'      , ## get all existing canvases 
    'canvas_partition' , ## split canvas into several pads with no space between pads 
    'partition_draw'   , ## split canvas into several pads with no space between pads & draw  
    'canvas_pull'      , ## split canvas into two pads with no vertical interspace
    'draw_pads'        , ## plot sequence of object on sequence of pads, adjustinng axis label size
    'AutoPlots'        , ## context manager to activate the auto-plotting machinery
    'auto_plots'       , ## ditto, but as function
    'UseWeb'           , ## constext manmager to use Web Displya 
    'useWeb'           , ## constext manmager to use Web Displya
    'setWebDisplay'    , ## set WebDisplay, see ROOT.TROOT.SetWebDisplay
    'use_pad'          , ## context manager to modifty TPad
    'KeepCanvas'       , ## context manager to keep/preserve currect canvas 
    'keepCanvas'       , ## context manager to keep/preserve currect canvas 
    'Canvas'           , ## context manager to create currect canvas
    'use_canvas'       , ## context manager to create currect canvas
    )
# =============================================================================
from   sys                     import version_info as python_version
from   ostap.core.ostap_types  import ordered_dict 
from   ostap.utils.cidict      import cidict
from   ostap.utils.utils       import KeepCanvas, keepCanvas 
from   ostap.core.core         import cidict_fun
from   ostap.core.core         import rootWarning
from   ostap.utils.utils       import which
import ostap.plotting.style
import ROOT, os, tempfile, math   
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.plotting.canvas' )
else                       : logger = getLogger( __name__ )
# =============================================================================
from ostap.plotting.makestyles import  ( canvas_width , canvas_height ,
                                         margin_left  , margin_right  ,
                                         margin_top   , margin_bottom )


# =============================================================================
#$ define WebDispla
#  @see TROOT::SetWebDisplay
def setWebDisplay ( web ) :
    """Define WebDispllay
    - see `ROOT.TROOT.SetWebDisplay`
    """
    wlow = web.lower()
    if    not web         : web = 'off'
    elif 'root'   == wlow : web = 'off'
    elif 'chrome' == wlow :
        if   which ( web  ) : pass
        elif which ( wlow ) : web = wlow 
        elif which ( 'google-chrome'        ) : web = 'google-chrome'
        elif which ( 'google-chrome-stable' ) : web = 'google-chrome-stable'
    ROOT.gROOT.SetWebDisplay( web )
    return ROOT.gROOT.GetWebDisplay()

# =============================================================================
## @class UseWeb
#  Context manager to redefien the web-display
#  @code
#  with UseWe ('chrome') :
#  ...
#  @endcode
#  @see TROOT::GetWebDisplay
#  @see TROOT::SetWebDisplay
class UseWeb(object) :
    """ Context manager to redefien the web-display
    >>> with UseWe ('chrome') :
    >>> ...
    - see `ROOT.TROOT.GetWebDisplay`
    - see `ROOT.TROOT.SetWebDisplay`
    """
    def __init__ ( self , web = "default" ) :

        self.__web  = web
        self.__prev = None 

    ## ENTER the context 
    def __enter__ ( self ) :
        self.__prev = ROOT.gROOT.GetWebDisplay()
        setWebDisplay ( self.__web )
        return self 
    
    ## EXIT the context 
    def __exit__ ( self , *_ ) :        
        if not self.__prev is None :
            setWebDisplay ( self.__prev )
            
# =============================================================================
## Context manager to redefien the web-display
#  @code
#  with useWeb ('chrome') :
#  ...
#  @endcode
#  @see TROOT::GetWebDisplay
#  @see TROOT::SetWebDisplay
def  useWeb( web  = 'default' ) :
    """ Context manager to redefien the web-display
    >>> with useWeb ('chrome') :
    >>> ...
    - see `ROOT.TROOT.GetWebDisplay`
    - see `ROOT.TROOT.SetWebDisplay`
    """
    return UseWeb ( web )
        
# =============================================================================
import ostap.core.config as cnf
web = cnf.general['WebDisplay']
if web :
    logger.debug ( 'Set WebDisplay to be `%s`' % web ) 
    setWebDisplay ( web )

# =============================================================================
## list of created canvases 
_canvases = [] 
# =============================================================================
## get the canvas
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-10-19
def getCanvas ( name   = 'glCanvas'    ,   ## canvas name 
                title  = 'Ostap'       ,   ## canvas title
                width  = canvas_width  ,   ## canvas width
                height = canvas_height ,   ## canvas height 
                **kwargs               ) : ## other properties 
    """Get create canvas/create new canvas
    
    >>> cnv = getCanvas ( 'glnewCanvas' , width = 1200 , height = 1000 )
    """
    if not name : name = 'glCanvas'

    cnvlst = ROOT.gROOT.GetListOfCanvases()
    cnv    = cnvlst.get ( name , None ) 
    if cnv and isinstance ( cnv , ROOT.TCanvas ) :
        _canvases.append ( cnv )
        set_pad ( cnv , **kwargs )
        return cnv 

    ## create new canvas
    mx    = max ( 100 , int ( math.floor ( 0.85 * width  ) ) + 1 )
    my    = max ( 100 , int ( math.floor ( 0.70 * height ) ) + 1 )
    
    wtopx = int ( ( 30 * len ( cnvlst ) ) % mx )
    wtopy = int ( ( 25 * len ( cnvlst ) ) % my ) 

    ## cnv  = ROOT.TCanvas ( 'glCanvas', 'Ostap' , width , height )
    cnv  = ROOT.TCanvas ( name , title , wtopx , wtopy , width , height )
    ## adjust newly created canvas
    ## @see http://root.cern.ch/root/html/TCanvas.html#TCanvas:TCanvas@4
    groot = ROOT.ROOT.GetROOT() 
    if not groot.IsBatch() :
        dw = width  - cnv.GetWw()
        dh = height - cnv.GetWh()
        cnv.SetWindowSize ( width + dw , height + dh )

    set_pad ( cnv , **kwargs )
    
    _canvases.append ( cnv ) 
    return cnv

# =============================================================================
import atexit
@atexit.register 
def clean_canvases () :
    while _canvases : _canvases.pop()
    
# =============================================================================
all_extensions = (
    'pdf'  , 'png'  , 'gif' ,
    'eps'  , 'ps'   ,
    'cxx'  , 'c'    , 
    'jpg'  , 'jpeg' , 'svg' , 
    'root' , 'xml'  , 'xpm' , 
    'tiff' , 'tex'  , 
    'json'
    )

# =============================================================================
## Define simplified print for TCanvas
#  Alows to create several output file types  at once
#  - if extension is equal to <code>tar</code> or <code>tgz</code>,
#    (gzipped) tar-files is created
#  - if extension is equal to <code>zip</code>, zip-archive is created
#  @code
#  canvas.print_ ( 'A' ) ## 
#  canvas.save   ( 'A' ) ## ditto 
#  canvas >> 'A'         ## ditto
#  @endcode 
def _cnv_print_ ( cnv , fname , exts = ( 'pdf'  , 'png' , 'eps'  , 'C'   ,
                                         'jpg'  , 'gif' , 'json' , 'svg' ) ) :
    """A bit simplified version for TCanvas print
    It Alows to create several output file types  at once
    - if extension is equal to `tar` or `tgz`, single (gzipped) tar-files is created
    - if extension is equal to `zip`, single zip-archive is created 
    >>> canvas.print_ ( 'A' )
    >>> canvas.save   ( 'A' ) ## ditto
    >>> canvas >> 'fig'       ## ditto
    """
    #
    cnv.Update ()

    dirname = os.path.dirname ( fname )
    if dirname and not os.path.exists ( dirname ) :
        from ostap.utils.basic import make_dirs
        dirname = os.path.abspath ( dirname ) 
        logger.debug ( "create directory %s" % os.path.abspath ( dirname ) ) 
        make_dirs ( dirname ) 
    
    n , e  = os.path.splitext ( fname )

    el = e.lower()
    if el.startswith('.') : el = el[1:]
    
    if n and el and ( el in all_extensions ) : 
        with rootWarning () :
            cnv.Update   () 
            cnv.Print    ( fname )
            logger.debug ( 'Canvas --> %s' % fname )
            return cnv
        
    if n and el in ( 'tgz' , 'gztar' , 'targz' , 'tar' ,
                     'zip' ,
                     'tbz' , 'tbz2'  , 'tarbz' , 'tarbz2' , 'bztar' , 'bz2tar' ,                     
                     'txz' , 'tlz'   , 'tarxz' , 'tarlz'  , 'xztar' , 'lztar'  ) :            
        files = [] 
        for ext in exts :
            with rootWarning () :
                name = n + '.' + ext
                cnv.Print    ( name )
                logger.debug ( 'Canvas --> %s' % name )
                if os.path.exists ( name ) and os.path.isfile ( name ) : 
                    files.append ( name )
                    
        if   files and el in  ( 'tgz' , 'targz'  , 'gztar' ) : 
            import tarfile
            with tarfile.open ( fname , "w:gz" ) as output :
                for f in files : output.add   ( f )
            if os.path.exists ( fname ) :
                logger.debug  ( 'tgz-archive created %s' % fname )
                
        elif files and el in  ( 'tar' , ) : 
            import tarfile
            with tarfile.open ( fname , "w" ) as output :
                for f in files : output.add   ( f )
            if os.path.exists ( fname ) :
                logger.debug  ( 'tar-archive created %s' % fname )
                
        elif files and el in  ( 'tbz' , 'tarbz' , 'tarbz2' , 'tbz2' , 'bztar' , 'bz2tar' ) : 
            import tarfile
            with tarfile.open ( fname , "w:bz2" ) as output :
                for f in files : output.add   ( f )
            if os.path.exists ( fname ) :
                logger.debug  ( 'tbz-archive created %s' % fname )

        elif files and el in  ( 'txz'   , 'tlz'   ,
                                'tarxz' , 'tarlz' ,
                                'xztar' , 'lztar' ) and 3 <= python_version.major : 
            import tarfile
            with tarfile.open ( fname , "w:xz" ) as output :
                for f in files : output.add   ( f )
            if os.path.exists ( fname ) :
                logger.debug  ( 'txz-archive created %s' % fname )
                            
        elif files and el in  ( 'zip' , ) : 
            import zipfile
            with zipfile.ZipFile( fname , "w"  ) as output :
                for f in files : output.write ( f )
            if os.path.exists ( fname ) :
                logger.debug  ( 'zip-archive created %s' % fname )

        for f in files :
            try :
                os.remove ( f )
            except OSError :
                pass

        return cnv
    
    for ext in exts :
        with rootWarning () :
            name = fname + '.' + ext
            cnv.Print    (  name )
            logger.debug ( 'Canvas --> %s' % name )

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

ROOT.TVirtualPad.save       = _cnv_print_
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
            directory = tempfile.mkdtemp ( prefix = 'ostap-plots-' )
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
## decorate ROOT.TObject
if not hasattr ( ROOT.TObject , 'draw_with_autoplot' ) :
    
    
    def _TO_draw_with_auto_ ( obj  , option = '' , *args , **kwargs ) :
        """ Draw an object     
        """
        result = obj.draw_ostap ( option, *args , **kwargs ) 

        if ROOT.gPad :
            plot = AutoPlots.plot()
            if plot :
                cnv = ROOT.gPad.GetCanvas ()
                if cnv : 
                    cnv >> plot
            
        return result
    
    _TO_draw_with_auto_.__doc__ += '\n' + ROOT.TObject.draw_ostap.__doc__
    
    ## add new method
    ROOT.TObject.draw_with_autoplot = _TO_draw_with_auto_ 
    
    ROOT.TObject.draw = ROOT.TObject.draw_with_autoplot
        
# =============================================================================
## get all known canvases 
def getCanvases () :
    """ Get all known canvases """
    
    return tuple ( ( c.GetName() for c in ROOT.gROOT.GetListOfCanvases() if  ( c and isintance ( c , ROOT.TCanvas ) ) ) ) 
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
def canvas_partition ( canvas                        , 
                       nx                            ,
                       ny                            ,
                       left_margin   = margin_left   , 
                       right_margin  = margin_right  , 
                       bottom_margin = margin_bottom , 
                       top_margin    = margin_top    ,
                       hSpacing      = 0.0           ,
                       vSpacing      = 0.0           ) :
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
    ## if parameters given in the absolute units, convert them into relative coordinates
    #
    
    if not 0 <  left_margin  < 1 :   left_margin = abs (   left_margin ) / wsx
    if not 0 < right_margin  < 1 :  right_margin = abs (  right_margin ) / wsx
    if not 0 < bottom_margin < 1 : bottom_margin = abs ( bottom_margin ) / wsy
    if not 0 < top_margin    < 1 :    top_margin = abs (    top_margin ) / wsy    
    if not 0 < vSpacing      < 1 :      vSpacing = abs ( vSpacing )      / wsy
    if not 0 < hSpacing      < 1 :      hSpacing = abs ( hSpacing )      / wsx

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

    ## delete the pad dictionary 
    del canvas.pads 
        
    ## make new empty dictionary 
    pads = {} 
            
    vStep    = ( 1.0 - bottom_margin - top_margin   - (ny-1) * vSpacing ) / ny
    if 0 > vStep : raise AttributeError('partition: v-step=%f' % vStep  )
        
    hStep    = ( 1.0 - left_margin   - right_margin - (nx-1) * hSpacing ) / nx 
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

        for iy in range ( ny ) :
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
            pname = 'glPad_%s_%d_%d' % ( canvas.GetName() , ix , iy )
            groot = ROOT.ROOT.GetROOT()
            pad   = groot.FindObject ( pname )
            if pad : del pad

            hposl = max ( 0 , min ( 1.0 , hposl ) )
            hposr = max ( 0 , min ( 1.0 , hposr ) )
            vposu = max ( 0 , min ( 1.0 , vposu ) )
            vposd = max ( 0 , min ( 1.0 , vposd ) )
            
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
            pads[  ( ix , iy ) ] = pad

    ## fill pads structure 
    for iy in reversed ( range ( ny ) ) : 
        for ix in range ( nx ) :
            key = ix , iy 
            canvas.pads [ key ] = pads[ key ]
            
    return canvas.pads 


# =============================================================================
## Make partition of Canvas into (nx)x(ny) adjuncent pads and plot objects
#  @code
#  canvas  = ...
#  objects = ...
#  pads    = canvas.partition_draw ( objects , 3 , 2 )
def partition_draw ( canvas                        ,
                     objects                       , 
                     nx                            ,
                     ny                            ,
                     left_margin   = margin_left   , 
                     right_margin  = margin_right  , 
                     bottom_margin = margin_bottom , 
                     top_margin    = margin_top    ,
                     hSpacing      = 0.0           ,
                     vSpacing      = 0.0           ) :
    """Perform partition of Canvas into pads with no inter-marginspads and draw objects
    >>> canvas  = ...
    >>> objects = ...
    >>> pads    = canvas.partition_draw ( objects , 3 , 2 )
    """
    pads = canvas_partition ( canvas  ,
                              nx = nx ,
                              ny = ny ,
                              left_margin   = left_margin   ,
                              right_margin  = right_margin  ,
                              bottom_margin = bottom_margin ,
                              top_margin    = top_margin    ,
                              hSpacing      = hSpacing      ,
                              vSpacing      = vSpacing      )

    canvas.Clear()
    for key , obj in zip ( pads, objects ) :
        
        canvas.cd(0)
        pad = pads [ key ] 
        pad.draw ()
        pad.cd   ()
        obj.draw ()

    canvas.Update()
    canvas.cd(0)
    return pads
    
    
# =============================================================================
def _cnv_pads_ ( self ) :
    """`pads' : get (an ordered) dict of pads for the given canvas partition (if prepared)
    """
    if not hasattr ( self , '__pads' ) :
        self.__pads  = ordered_dict ()
    return self.__pads
# =============================================================================
def _cnv_del_pads_ ( self ) :
    """ deleter for the created pad structure 
    """
    while self.pads :
        key , pad = self.pads.popitem ()
        if pad :
            logger.verbose ( 'delete pad %s' % pad .GetName() )
            del pad
# ==============================================================================
## property! 
ROOT.TCanvas.pads = property ( _cnv_pads_ , None , _cnv_del_pads_ )

ROOT.TCanvas.partition      = canvas_partition
ROOT.TCanvas.partition_draw = partition_draw

# ==============================================================================
## Split canvas in y-direction into non-equal pads,
#  proportionally to <code>heights</code> 
#  @code
#  canvas = ...
#  pads   = canvas.vsplit ( [1,2,1] )    
#  @endcode 
def canvas_vsplit ( canvas                        ,
                    heights                       ,
                    left_margin   = margin_left   , 
                    right_margin  = margin_right  , 
                    bottom_margin = margin_bottom , 
                    top_margin    = margin_top    ,
                    vSpacing      = 0.0           ) :
    """ Split canvas in y-direction into non-equal pads, proportionally to heights
    >>> canvas = ...
    >>> pads   = canvas.vsplit ( [1,2,1] )    
    """

    ## get the window size
    wsx = abs ( canvas.GetWindowWidth  () ) 
    wsy = abs ( canvas.GetWindowHeight () ) 

    #
    ## if parametes given in the absolute units, convert them into relative coordinates
    #
    
    if not 0 <  left_margin  < 1 :   left_margin = abs (   left_margin ) / wsx
    if not 0 < right_margin  < 1 :  right_margin = abs (  right_margin ) / wsx
    if not 0 < bottom_margin < 1 : bottom_margin = abs ( bottom_margin ) / wsy
    if not 0 < top_margin    < 1 :    top_margin = abs (    top_margin ) / wsy    
    if not 0 < vSpacing      < 1 :      vSpacing = abs (      vSpacing ) / wsy

    hSpacing = 0
    
    hposr, hposl, hmarr, hmarl, hfactor = 0.,0.,0.,0.,0.
    vposr, vposd, vmard, vmaru, vfactor = 0.,0.,0.,0.,0.
    
    nx = 1
    ny = len ( heights ) 

    vSize    = ( 1.0 - bottom_margin - top_margin   - ( ny - 1 ) * vSpacing ) 
    hSize    = ( 1.0 - left_margin   - right_margin - ( nx - 1 ) * hSpacing ) 

    vStep    = ( 1.0 - bottom_margin - top_margin   - ( ny - 1 ) * vSpacing ) / ny
    if 0 > vStep : raise AttributeError('partition: v-step=%f' % vStep  )
    
    hStep    = ( 1.0 - left_margin   - right_margin - ( nx - 1 ) * hSpacing ) / nx 
    if 0 > hStep : raise AttributeError('partition: h-step=%f' % hStep  )

    sumy = sum ( heights ) / vSize 
    hy   = [ h*vSize/sum(heights) for h in reversed ( heights ) ]

    hposl   = 0
    hposr   = left_margin + hStep 
    hfactor = hposr - hposl
    hmarl   = left_margin / hfactor
    hmarr   = 0.0

    del canvas.pads
    pads   = {}
    
    ix     = 0
    
    for iy , height in enumerate ( hy ) : 
        
        if 0 == iy : 
            vposd   = 0.0
            vposu   = bottom_margin + height
            vfactor = vposu - vposd 
            vmard   = bottom_margin / vfactor
            vmaru   = 0.0 
        elif ny == iy + 1 : 
            vposd   = vposu + vSpacing
            vposu   = vposd + height + top_margin
            vfactor = vposu - vposd
            vmard   = 0.0
            vmaru   = top_margin    / vfactor
        else :
            vposd   = vposu + vSpacing
            vposu   = vposd + height
            vfactor = vposu - vposd
            vmard   = 0.0
            vmaru   = 0.0
            
        canvas.cd ( 0 )
        pname = 'glPad_%s_%d_%d' % ( canvas.GetName() , ix , iy )
        groot = ROOT.ROOT.GetROOT()
        pad   = groot.FindObject ( pname )
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

        pads[ (0,iy) ] = pad 

    ## fill pads structure 
    for iy in reversed ( range ( ny ) ) : 
        key = 0 , iy 
        canvas.pads [ key ] = pads [ key ]
            
    return canvas.pads 

ROOT.TCanvas.vsplit = canvas_vsplit 

# ==============================================================================
## Perform partition of Canvas into 1x2 non-equal pads with no inter-margins
#  @code
#  canvas = ...
#  pads   = canvas.pull_partition ( 0.20 ) ## top-pad 4-times larger    
#  @endcode 
def canvas_pull ( canvas                        ,
                  ratio         = 4.0           ,
                  left_margin   = margin_left   , 
                  right_margin  = margin_right  , 
                  bottom_margin = margin_bottom , 
                  top_margin    = margin_top    ,
                  vSpacing      = 0.0           ) :
    """ Perform partition of Canvas into 1x2 non-equal pads with no inter-margins
    
    >>> canvas = ...
    >>> pads   = canvas.pull_partition ( 4.0 ) ## top-pad 4-times larger    
    
    """
    return canvas_vsplit ( canvas                        ,
                           heights       = ( 1 , ratio ) ,
                           left_margin   =   left_margin ,
                           right_margin  =  right_margin ,
                           bottom_margin = bottom_margin ,
                           top_margin    =    top_margin ,
                           vSpacing      = vSpacing      )

ROOT.TCanvas.pull_partition = canvas_pull

# =============================================================================
## draw sequence of object on sequence of pads,
#  - the label font size is adjusted to be uniform (in pixels) 
#  @code
#  pads   = ...
#  frames = ...
#  draw_pads ( frames , pads , fontsize = 25 ) 
#  @endcode 
def draw_pads ( objects            ,
                pads               ,
                fontsize   = 36    ,
                trim_left  = False ,
                trim_right = False ) :
    """ Draw sequence of object on sequence of pads,
    - the label font size is adjusted to be uniform (in pixels)     
    >>> pads   = ...
    >>> frames = ...
    >>> draw_pads ( frames , pads , fontsize = 25 ) 
    """
    assert isinstance  ( fontsize , int ) and 5 < fontsize , 'Invalid fontsize %s [pixels] ' % fontsize
    
    for obj , pad_ in zip ( objects , pads ) : 
        
        if isinstance ( pad_ , ROOT.TPad ) : pad = pad_
        else                               : pad = pads [ pad_ ] 
        
        c = pad.GetCanvas()
        if c : c.cd(0)
        
        pad.draw ()
        pad.cd   ()

        ## redefine the label font and size 
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

        if  ( trim_left or trim_right ) and hasattr ( obj , 'GetXaxis' ) :
            
            axis  = obj.GetXaxis()
            xmin  = axis.GetXmin()
            xmax  = axis.GetXmax()
            delta = xmax - xmin
            
            if   trim_left and isinstance ( trim_left , float ) :
                xmin += trim_left * delta
            elif trim_left :
                xmin += 0.001 * delta
                
            if   trim_right and isinstance ( trim_right , float ) :
                xmax -= trim_right * delta
            elif trim_right :
                xmax -= 0.001 * delta 
                
            axis.SetLimits ( xmin , xmax )

        ## draw object on the pad 
        obj.draw ()
        
        if c : c.cd(0)

# =============================================================================
## change main parametes of TAttPad
def set_pad ( pad , **config ) :
    """ Change main parametes of `TAttPad`"""

    conf = cidict ( transform = cidict_fun )
    conf.update ( config )

    changed = {}

    if 'frame_border_mode' in conf :
        changed[ 'frame_border_mode' ] = pad.GetFrameBorderMode()
        pad.SetFrameBorderMode ( conf.pop ('frame_border_mode') )

    if 'frame_border_size' in conf :
        changed[ 'frame_border_size' ] = pad.GetFrameBorderSize ()
        pad.SetFrameBorderSize ( conf.pop ('frame_border_size') )
        
    if 'frame_fill_color' in conf :
        changed[ 'frame_fill_color' ] = pad.GetFrameFillColor ()
        pad.SetFrameFillColor ( conf.pop ('frame_fill_color') )

    if 'frame_fill_style' in conf :
        changed[ 'frame_fill_style' ] = pad.GetFrameFillStyle ()
        pad.SetFrameFillStyle ( conf.pop ('frame_fill_style') )

    if 'frame_line_color' in conf :
        changed[ 'frame_line_color' ] = pad.GetFrameLineColor ()
        pad.SetFrameLineColor ( conf.pop ('frame_line_color') )

    if 'frame_line_style' in conf :
        changed[ 'frame_line_style' ] = pad.GetFrameLineStyle ()
        pad.SetFrameLineStyle ( conf.pop ('frame_line_style') )

    if 'frame_line_width' in conf :
        changed[ 'frame_line_width' ] = pad.GetFrameLineWidth ()
        pad.SetFrameLineWidth ( conf.pop ('frame_line_width') )

    if 'a_file' in conf :
        changed[ 'a_file' ] = pad.GetAfile ()
        pad.SetAfile ( conf.pop ('a_file') )

    if 'a_stat' in conf :
        changed[ 'a_stat' ] = pad.GetAstat ()
        pad.SetAstat ( conf.pop ('a_stat') )

    if 'x_file' in conf :
        changed[ 'x_file' ] = pad.GetXfile ()
        pad.SetXfile ( conf.pop ('x_file') )

    if 'x_stat' in conf :
        changed[ 'x_stat' ] = pad.GetXstat ()
        pad.SetXstat ( conf.pop ('x_stat') )

    if 'y_file' in conf :
        changed[ 'y_file' ] = pad.GetYfile ()
        pad.SetYfile ( conf.pop ('y_file') )

    if 'y_stat' in conf :
        changed[ 'y_stat' ] = pad.GetYstat ()
        pad.SetYstat ( conf.pop ('y_stat') )
    
    if 'top_margin' in conf or 'margin_top' in conf :                
        changed ['margin_top']  = pad.GetTopMargin()
        if 'top_margin' in conf    : pad.SetTopMargin    ( conf.pop ( 'top_margin'   ) )
        else                       : pad.SetTopMargin    ( conf.pop ( 'margin_top'   ) )
        
    if 'bottom_margin' in conf or 'margin_bottom' in conf :            
        changed ['margin_bottom']  = pad.GetBottomMargin()
        if 'bottom_margin' in conf : pad.SetBottomMargin ( conf.pop ( 'bottom_margin' ) )
        else                       : pad.SetBottomMargin ( conf.pop ( 'margin_bottom' ) ) 
        
    if 'left_margin' in conf or 'margin_left' in conf :                
        changed ['margin_left']  = pad.GetLeftMargin()
        if 'left_margin' in conf   : pad.SetLeftMargin   ( conf.pop ( 'left_margin'   ) )
        else                       : pad.SetLeftMargin   ( conf.pop ( 'margin_left'   ) )
        
    if 'right_margin' in conf or 'margin_right' in conf :                
        changed ['margin_right']  = pad.GetRightMargin()
        if 'right_margin' in conf  : pad.SetRightMargin  ( conf.pop ( 'right_margin'  ) )
        else                       : pad.SetRightMargin  ( conf.pop ( 'margin_right'  ) )

    if 'grid_x' in conf or 'x_grid' in conf :
        changed ['grid_x']  = pad.GetGridx () 
        if 'grid_x' in conf : pad.SetGridx ( conf.pop ( 'grid_x' ) ) 
        else                : pad.SetGridx ( conf.pop ( 'x_grid' ) )

    if 'grid_y' in conf or 'y_grid' in conf :
        changed ['grid_y']  = pad.GetGridy () 
        if 'grid_y' in conf : pad.SetGridy ( conf.pop ( 'grid_y' ) ) 
        else                : pad.SetGridy ( conf.pop ( 'y_grid' ) ) 

    if 'log_x' in conf or 'x_log' in conf :
        changed ['log_x']  =  pad.GetLogx () 
        if 'log_x' in conf  : pad.SetLogx ( conf.pop ( 'log_x' ) ) 
        else                : pad.SetLogx ( conf.pop ( 'x_log' ) )
        
    if 'log_y' in conf or 'y_log' in conf :
        changed ['log_y']  =  pad.GetLogy () 
        if 'log_y' in conf  : pad.SetLogy ( conf.pop ( 'log_y' ) ) 
        else                : pad.SetLogy ( conf.pop ( 'y_log' ) )
        
    if conf :
        logger.warning ("set_pad: unprocessed items: %s" % conf ) 

    return changed 
    
# =============================================================================
## helper context manager for <code>TAttPad</code> objects
#  @see TAttPad 
class UsePad(object) :
    """ Helper context manager for `TAttPad` objects
    - see `TAttPad`
    """

    def  __init__ ( self , pad = None , **config ) :

        self.__pad     = pad if ( pad and isinstance ( pad , ROOT.TAttPad ) ) else None 
        self.__config  = config
        self.__changed = {} 
        
    def __enter__ ( self ) :
            
        if not self.__pad and ROOT.gPad : self.__pad = ROOT.gPad
        
        if self.pad : 
            self.__changed = set_pad ( self.pad , **self.config )

        return self 
        
    def __exit__ ( self , *_ ) :

        if self.pad and self.changed :
            set_pad ( self.pad , **self.changed )

    @property
    def pad ( self ) :
        """``pad'' : pad to be configured"""
        return self.__pad
    
    @property
    def config ( self ) :
        """``config'' : cofiguration pad to be"""
        return self.__config

    @property
    def changed ( self ) :
        """``changed'' : changed parameters"""
        return self.__changed
    
# =============================================================================
## helper context manager to modify <code>TAttPad</code>
#  @see TAttPad 
def use_pad ( pad , **config ) :
    """Helper context manager for `TAttPad` objects
    - see `TAttPad`
    """
    return UsePad ( pad , **config ) 

usePad = use_pad

# =============================================================================
## @class Canvas
#  helper context manager to create and configure a canvas (and pad)
#  @code
#  with Canvas ( title = 'Canvas #2' , width = 1000 ) :
#  ... 
#  @endcode
class Canvas(KeepCanvas) :
    """ Helper context manager to create and configure a canvas (and pad)
    >>> with Canvas ( title = 'Canvas #2' , width = 1000 ) :
    >>> ... 
    """
    def __init__ ( self                   ,
                   name   = ''            ,
                   title  = ''            ,
                   width  = canvas_width  ,   ## canvas width
                   height = canvas_height ,   ## canvas height 
                   wait   = 0             ,   ## pause before exit
                   plot   = ''            ,   ## produce the plot at __exit__
                   **kwargs               ) : ## Pad configuration
        
        self.__name   = name
        self.__title  = title 
        self.__width  = width
        self.__height = height
        self.__kwargs = kwargs
        self.__cnv    = None

        if plot :
            plot = plot.strip() 
            while '  ' in plot : plot = plot.replace ( '  ' , ' ' )                
            plot = plot.strip().replace ( ' ' , '_' )
            
        self.__plot   = plot            
        ## 
        KeepCanvas.__init__ ( self , wait ) 
        
    ## context manager: exit 
    def __enter__ ( self ) :

        ## 1) use context manager 
        KeepCanvas.__enter__ ( self )

        if not self.__name :
            cnvlst = ROOT.gROOT.GetListOfCanvases() 
            self.__name = 'gl_canvas#%d' % len ( cnvlst )
            while self.__name in cnvlst : 
                h = self.__name , title , width , height , len ( cnvlst ) 
                self.__name = 'gl_canvas#%d' % hash ( h ) 
                
        if not self.__title :
            self.__title = self.__name

        ## 2) create/use new canvas 
        self.__cnv = getCanvas ( name   = self.__name   ,
                                 title  = self.__title  ,
                                 width  = self.__width  ,
                                 height = self.__height )
        
        self.__name  = self.__cnv.GetName  () 
        self.__title = self.__cnv.GetTitle () 

        ## 3) make it active 
        self.__cnv.cd() 

        ## 4) apply pad settings
        if self.__kwargs :
            set_pad ( ROOT.gPad , **self.__kwargs ) 
            
        return self.__cnv  ## return current canvas 
    
    ## context manager: exit 
    def __exit__ ( self , *_ ) :

        if self.__cnv :
            self.__cnv.Update()            
            if self.__plot :
                self.__cnv >> self.__plot 
                    
        KeepCanvas.__exit__ ( self , *_ ) 
        self.__cnv = None 
    
# =============================================================================
## helper context manager to create and configure a canvas (and pad)
#  @code
#  with use_canvas ( title = 'Canvas #2' , width = 1000 ) :
#  ... 
#  @endcode
def use_canvas ( name   = ''            ,
                 title  = ''            ,
                 width  = canvas_width  ,   ## canvas width
                 height = canvas_height ,   ## canvas height
                 wait   = 0             ,   ## pause before exit 
                 **kwargs               ) : ## Pad configuration
    """Helper context manager to create and configure a canvas (and pad)
    >>> with use_canvas ( title = 'Canvas #2' , width = 1000 ) :
    >>> ... 
    """
    return Canvas ( name   = name   ,
                    title  = title  ,
                    width  = width  ,   ## canvas width
                    height = height ,   ## canvas height
                    wait   = wait   ,   ## pause before exit 
                    **kwargs        )   ## Pad configuration


    
# =============================================================================
_decorated_classes_  = (
    ROOT.TVirtualPad , 
    ROOT.TCanvas     , 
    ROOT.TObject     , 
    )

_new_methods_        = (
    ROOT.TVirtualPad.__rshift__ ,
    ROOT.TCanvas.partition      ,
    ROOT.TCanvas.partition_draw ,
    ROOT.TCanvas.pull_partition ,
    ROOT.TObject    .draw       , 
    ROOT.TObject    .draw_with_autoplot, 
    )
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
