#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# =============================================================================
# $Id: Canvas$ 
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
#                    $Revision$
#  Last modification $Date$
#                 by $Author$
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
__all__     = ( 'getCanvas' , 'getCanvases' )
# =============================================================================
import ROOT
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
    """
    Get create canvas/create new canvas
    
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
    from ostap.logger.utils import rootWarning 
    n,d,e = fname.rpartition('.')
    if d and e.lower() in all_extensions : 
        with rootWarning () :  
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

ROOT.TCanvas.print_     = _cnv_print_
ROOT.TCanvas.__rshift__ = _cnv_rshift_

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
if '__main__' == __name__ :
    
    from ostap.logger.line import line 
    logger.info ( __file__ + '\n' + line  )
    logger.info ( 80*'*' )
    logger.info ( __doc__  )
    logger.info ( 80*'*' )
    logger.info ( ' Author  : %s' %         __author__    ) 
    logger.info ( ' Version : %s' %         __version__   ) 
    logger.info ( ' Date    : %s' %         __date__      )
    logger.info ( ' Symbols : %s' %  list ( __all__     ) )
    logger.info ( 80*'*' ) 

# =============================================================================
# The END 
# =============================================================================
