#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/histos/histos.py 
#  Module with decoration of TH* objects for efficient use in python
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
""" Module with decoration of TH* objects for efficient use in python
"""
# =============================================================================
__version__ = "$Revision: 207180 $"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
# =============================================================================
__all__     = (
    #
    'binomEff_h1'     , ## calculate binomial efficiency for 1D-histos
    'binomEff_h2'     , ## calculate binomial efficiency for 2D-ihstos
    'binomEff_h3'     , ## calculate binomial efficiency for 3D-ihstos
    'h1_axis'         , ## book 1D-histogram from axis 
    'h2_axes'         , ## book 2D-histogram from axes
    'h3_axes'         , ## book 3D-histogram from axes
    'axis_bins'       , ## convert list of bin edges to axis
    've_adjust'       , ## adjust the efficiency to be in physical range
    'Histo1DFun'      , ## 1D-histogram as function object 
    'Histo2DFun'      , ## 2D-histogram as function object 
    'Histo3DFun'      , ## 3D-histogram as function object
    #
    )
# =============================================================================
from   ostap.core.core                import ( cpp      , Ostap     , cidict_fun      , 
                                               ROOTCWD  , rootID    , typename        , 
                                               funcID   , funID     , fID             ,
                                               histoID  , hID       , dsID            ,
                                               VE       , SE        , WSE             ,
                                               binomEff , binomEff2 ,
                                               zechEff  , wilsonEff , agrestiCoullEff , 
                                               iszero   , isequal   , inrange         , 
                                               isint    , islong    , is_sorted       ,
                                               natural_entry        ,
                                               natural_number       )
from   ostap.math.base                import ( frexp10      , isequalf      ,
                                               pos_infinity , neg_infinity  ) 
from   ostap.math.math_ve             import significance
from   ostap.core.ostap_types         import ( integer_types  , num_types   ,
                                               long_type      , sized_types , 
                                               sequence_types ) 
from   ostap.utils.progress_bar       import progress_bar
from   ostap.core.meta_info           import root_info, python_info
from   ostap.math.random_ext          import poisson
from   ostap.utils.utils              import accumulate
from   ostap.utils.cidict             import cidict
import ostap.logger.table             as     T 
import ostap.stats.moment 
import ostap.plotting.draw_attributes 
import ostap.io.root_file 
import ROOT, sys, math, ctypes, array  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.histos.histos' )
else                       : logger = getLogger( __name__          )
# =============================================================================
logger.debug ( 'Decoration of historams')
# =============================================================================
_new_methods_  = []  
# =============================================================================
## ensure that object/histogram is created in ROOT.gROOT
#  @attention clone is always goes to ROOT main memory!
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-08-14
def _h_new_init_ ( self , *args ) :
    """ Modified TH* constructor:
    - ensure that created object/histogram goes to ROOT main memory
    """
    with ROOTCWD () :
        groot = ROOT.ROOT.GetROOT()
        groot.cd() 
        init = self._old_init_   ( *args )
        self.SetDirectory ( groot )  ## NB! 
        ## optionally:
        if not self.GetSumw2() : self.Sumw2()
        return init 
    
# =============================================================================
## a bit modified 'Clone' function for histograms
#  - it automatically assign unique ID
#  - it ensures that cloned histogram is not going to die with
#    the accidentally opened file/directory
#  - a title can be optionally redefined 
#  @attention clone is always goes to ROOT main memory!
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h_new_clone_ ( self , name = '' , title = '' , prefix = 'h_' ) :
    """ Modifiled Clone-function
    - it automatically assigns unique ID
    - it ensures that cloned histogram is not going to die with
    the accidentally opened file/directory
    - a title can be optionally redefined 
    """
    if not name : name = hID ( prefix = prefix )
    #
    with ROOTCWD() :
        groot = ROOT.ROOT.GetROOT()
        groot.cd() 
        nh = self._old_clone_ ( name )
        nh.SetDirectory ( groot )
        ## optionally 
        if not nh.GetSumw2() : nh.Sumw2()
        
    if title : nh.SetTitle ( title ) 
    return nh

# =============================================================================
## apply decoration
for h in ( ROOT.TH1F , ROOT.TH1D ,
           ROOT.TH2F , ROOT.TH2D ,
           ROOT.TH3F , ROOT.TH3D ) :

    if hasattr ( h , '_new_clone_' ) and hasattr ( h , '_old_clone_' ) : pass
    else : 
        h._old_clone_ =  h.Clone
        h._new_clone_ = _h_new_clone_
        h.Clone       = _h_new_clone_
        h.clone       = _h_new_clone_

    if hasattr ( h , '_new_init_' ) and hasattr ( h , '_old_init_' ) : pass
    else : 
        h._old_init_ =  h.__init__ 
        h._new_init_ = _h_new_init_
        ## h.__init__   = _h_new_init_  

# =============================================================================
# Decorate histogram axis and iterators 
# =============================================================================

# =============================================================================
## iterator for histogram  axis 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _axis_iter_1_ ( a ) :
    """ Iterator for axis
    >>> axis = ...
    >>> for i in axis : 
    """
    i = 1
    s = a.GetNbins()
    while i <= s :
        yield i
        i+=1        

# =============================================================================
## iterator for histogram  axis (reversed order) 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _axis_iter_reversed_ ( a ) :
    """Iterator for axis
    >>> axis = ...
    >>> for i in reverse(axis) : 
    """
    s = a.GetNbins()
    i = int ( s )
    while 1 <= i :
        yield i
        i-=1

# =============================================================================
## get new "scaled" axis 
#  @code
#  axis   = ...
#  scaled = axis.scale ( 5 ) 
#  @endcode
def _axis_scale_ ( axis , scale ) :
    """ Get new `scaled' axis 
    >>> axis   = ...
    >>> scaled = axis.scale ( 5 ) 
    """
    assert isinstance ( scale , num_types ) and 0 < scale , \
           'Invalid "scale" argument %s' % scale
    bins = axis.GetXbins ()
    ## non-uniform bins are here 
    if bins and 2 <= len ( bins ) :
        newbins = array.array ('d' , ( v * scale for v in bins ) )
        return ROOT.TAxis ( len ( newbins ) - 1 , newbins )
    
    return ROOT.TAxis ( axis.GetNbins ()         ,
                        axis.GetXmin  () * scale ,
                        axis.GetXmax  () * scale )


ROOT.TAxis . __iter__     = _axis_iter_1_
ROOT.TAxis . __reversed__ = _axis_iter_reversed_
ROOT.TAxis . __contains__ = lambda s , i : 1 <= abs(i) <= s.GetNbins()
ROOT.TAxis .   scale      = _axis_scale_
        
# =============================================================================
## get item for the 1-D histogram 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_get_item_ ( h1 , ibin ) :
    """``Get-item'' for the 1D-histogram :
    >>> histo = ...
    >>> ve    = histo[ibin]    
    """
    if isinstance ( ibin , slice ) :
        start = ibin.start
        step  = ibin.step
        stop  = ibin.stop
        if step and step != 1: raise IndexError ('Step must be +1!')
        if start is None : start = 0 
        if stop  is None : stop  = len ( h1 ) 
        return _h1_getslice_ ( h1 , start , stop ) 
        
    a  = h1.GetXaxis()
    nb = a.GetNbins ()
    #
    if     1 <=  ibin <= nb : pass
    elif   1 <= -ibin <= nb : ibin += ( nb + 1 )
    else                    : raise IndexError 
    #
    val = h1.GetBinContent ( ibin ) 
    err = h1.GetBinError   ( ibin )
    #
    return VE ( val , err * err )

# ==========================================================================
## get item for the 2D histogram 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_set_item_ ( h1 , ibin , v ) :
    """ ``Set-item'' for the 1D-histogram :    
    >>> histo[ibin] = value     
    """
    #
    ## check the validity of the bin
    #
    a  = h1.GetXaxis()
    nb = a.GetNbins ()
    #
    if     1 <=  ibin <= nb : pass
    elif   1 <= -ibin <= nb : ibin += ( nb + 1 ) 
    else                    : raise IndexError 
    #
    ## treat the value:
    #
    vv = VE ( v ) 
    if   isinstance ( v , integer_types ) :
        
        if   0  < v   : vv = VE ( v , v )  
        elif 0 == v   : vv = VE ( 0 , 1 ) 
        else          : vv = VE ( v , 0 ) 

    elif isinstance ( v , float ) :
        
        if   islong ( v ) : return _h1_set_item_ ( h1 , ibin , long_type ( v ) )
        else              : vv = VE ( v , 0 ) 

    #
    h1.SetBinContent ( ibin , vv.value () )
    h1.SetBinError   ( ibin , vv.error () )
    
ROOT.TH1F. __setitem__ = _h1_set_item_
ROOT.TH1D. __setitem__ = _h1_set_item_


# ==========================================================================
## get item for the 2D histogram
#  @code 
#  histo[ ix , iy ] = value     
#  histo[ ix][ iy ] = value ## ditto
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_set_item_ ( h2 , ibin , v ) :
    """``Set-item'' for the 2D-histogram :    
    >>> histo[ ix , iy ] = value     
    >>> histo[ ix][ iy ] = value ## ditto    
    """
    ##
    vv = VE ( v )
    if isinstance   ( v , integer_types ) :
        
        if   0  < v     : vv = VE ( v , v ) 
        elif 0 == v     : vv = VE ( 0 , 1 ) 
        else            : vv = VE ( v , 0 ) 

    elif isinstance ( v , float ) :
        
        if islong ( v ) : return _h2_set_item_ ( h2 , ibin , long_type ( v )  )
        else            : vv = VE ( v , 0 ) 
        
    ## check and adjust bins 
    ix = ibin[0] 
    ax = h2.GetXaxis()
    nx = ax.GetNbins ()
    #
    if     1 <=  ix <= nx : pass
    elif   1 <= -ix <= nx : ix += ( nx + 1 ) 
    else                  : raise IndexError 
    #
    iy = ibin[1] 
    ay = h2.GetYaxis()
    ny = ay.GetNbins ()
    #
    if     1 <=  iy <= ny : pass
    elif   1 <= -iy <= ny : ix += ( ny + 1 ) 
    else                  : raise IndexError 
    #
    
    h2.SetBinContent ( ix , iy , vv.value () )
    h2.SetBinError   ( ix , iy , vv.error () )
    
ROOT.TH2F. __setitem__ = _h2_set_item_
ROOT.TH2D. __setitem__ = _h2_set_item_

# ==========================================================================
## get item for the 3D histogram
#  @code
#  histo[ ix , iy ,  iz ] = value
#  histo[ ix , iy ][ iz ] = value ## ditto 
#  histo[ ix][ iy ,  iz ] = value ## ditto 
#  histo[ ix][ iy ][ iz ] = value ## ditto 
#  @endcode    
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_set_item_ ( h3 , ibin , v ) :
    """``Set-item'' for the 3D-histogram :
    >>> histo[ ix , iy ,  iz ] = value
    >>> histo[ ix , iy ][ iz ] = value ## ditto 
    >>> histo[ ix][ iy ,  iz ] = value ## ditto 
    >>> histo[ ix][ iy ][ iz ] = value ## ditto     
    """
    #
    vv = VE ( v ) 
    if   isinstance ( v , integer_types ) :
        
        if   0  < v     : vv = VE ( v , v )
        elif 0 == v     : vv = VE ( 0 , 1 ) 
        else            : vv = VE ( v , 0 ) 
        
    elif isinstance ( v , float ) :
        
        if islong ( v ) : return _h3_set_item_ ( h3 , ibin , long_type ( v ) )
        else            : vv = VE ( v , 0 ) 

    ## check and adjust bins 
    ix = ibin[0] 
    ax = h3.GetXaxis()
    nx = ax.GetNbins ()
    #
    if     1 <=  ix <= nx : pass
    elif   1 <= -ix <= nx : ix += ( nx + 1 ) 
    else                  : raise IndexError 
    #
    iy = ibin[1] 
    ay = h3.GetYaxis()
    ny = ay.GetNbins ()
    #
    if     1 <=  iy <= ny : pass
    elif   1 <= -iy <= ny : ix += ( ny + 1 ) 
    else                  : raise IndexError 
    #
    iz = ibin[2] 
    az = h3.GetZaxis()
    nz = az.GetNbins ()
    #
    if     1 <=  iz <= nz : pass
    elif   1 <= -iz <= nz : ix += ( nz + 1 ) 
    else                  : raise IndexError 
    #
    #
    h3.SetBinContent ( ix , iy , iz , vv.value () )
    h3.SetBinError   ( ix , iy , iz , vv.error () )
    
ROOT.TH3F. __setitem__ = _h3_set_item_
ROOT.TH3D. __setitem__ = _h3_set_item_

# =============================================================================
## @class _H2_item_
#  Helper class to extend 'get/set'item
#  @code
#  h2 = ...
#  v  = h3[1][2]              
#  v  = h3[1][2](0.44)    ## as 1D-function, use interpolation
#  h2[1][2] = v2
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
class _H2_item_(object):
    """Helper class to extend 'get/set'item
    >>> h2 = ...
    >>> v  = h3[1][2]              
    >>> v  = h3[1](0.44)    ## as 1D-function, use interpolation
    >>> h2[1][2] = v2
    """
    def __init__ ( self , h2 , xbin ) :
        if not xbin in h2.GetXaxis() : raise IndexError
        self._h2   = h2
        self._xbin = xbin 
    def __getitem__ ( self , ybin         ) :
        return self._h2[ self._xbin , ybin ]
    def __setitem__ ( self , ybin , value ) :
        self._h2[ self._xbin , ybin ] = value
    def __call__    ( self , y , func = lambda s : s , interpolate = True )  :
        x = self._h2.GetXaxis().GetBinCenter ( self._xbin )
        return self._h2( x , y , func , interpolate ) 
                      
# =============================================================================
## get item for the 2D histogram
#  @code
#  histo = ...
#  ve    = histo[ix, iy]
#  ve    = histo[ix][iy] ## ditto
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_get_item_ ( h2 , ibin ) :
    """ ``Get-item'' for the 2D-histogram :
    >>> histo = ...
    >>> ve    = histo[ix, iy]
    >>> ve    = histo[ix][iy] ## ditto 
    """
    #
    if isinstance ( ibin , integer_types ) :
        return _H2_item_ ( h2 , ibin )
    
    ix = ibin[0]
    ax = h2.GetXaxis()
    nx = ax.GetNbins()
    #
    if   1<=  ix <= nx : pass
    elif 1<= -ix <= nx : ix += ( nx + 1 ) 
    else               : raise IndexError
    #
    iy = ibin[1]
    ay = h2.GetYaxis()
    ny = ay.GetNbins()
    #
    if   1<=  iy <= ny : pass
    elif 1<= -iy <= ny : iy += ( ny + 1 ) 
    else               : raise IndexError
    #
    val = h2.GetBinContent ( ix , iy )  
    err = h2.GetBinError   ( ix , iy ) 
    #
    return VE ( val , err * err ) 

# =============================================================================
## @class _H3_item_2_
#  Helper  class to extend 'get/set'item
#  @code
#  h3 = ...
#  v  = h3[1][2 ,3]              
#  v  = h3[1][2][3]
#  v  = h3[1](0.13, 0.44) ## as 2D-function, use interpolation
#  v  = h3[1][2](0.44)    ## as 1D-function, use interpolation
#  h3[1][2 ,3] = v1  
#  h3[1][2][3] = v2
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
class _H3_item_2_(object):
    """Helper class to extend 'get/set'item
    >>> h3 = ...
    >>> v  = h3[1][2 ,3]              
    >>> v  = h3[1][2][3]
    >>> v  = h3[1](0.13, 0.44) ## as 2D-function, use interpolation
    >>> v  = h3[1][2](0.44)    ## as 1D-function, use interpolation
    >>> h3[1][2 ,3] = v1  
    >>> h3[1][2][3] = v2
    """
    def __init__ ( self , h3 , xbin , ybin ) :
        self._h3   = h3 
        if not xbin in h3.GetXaxis () : raise IndexError
        if not ybin in h3.GetYaxis () : raise IndexError
        self._xbin = xbin
        self._ybin = ybin
        
    def __getitem__ ( self , zbin  ) :
        return self._h3 [ self._xbin , self._ybin , zbin ]
    
    def __setitem__ ( self , zbin , value ) :
        self._h3[ self._xbin , self._ybin , zbin ] = value

    def __call__    ( self , z , func = lambda s : s , interpolate = True )  :
        x = self._h3.GetXaxis().GetBinCenter ( self._xbin )
        y = self._h3.GetYaxis().GetBinCenter ( self._ybin )
        return self._h3( x , y , z , func , interpolate ) 
    
# =============================================================================
## @class _H3_item_1_
#  Helper  class to extend 'get/set'item
#  @code
#  h3 = ...
#  v  = h3[1][2 ,3]              
#  v  = h3[1][2][3]
#  v  = h3[1](0.13, 0.44) ## as 2D-function, use interpolation
#  v  = h3[1][2](0.44)    ## as 1D-function, use interpolation
#  h3[1][2 ,3] = v1  
#  h3[1][2][3] = v2
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
class _H3_item_1_(object):
    """Helper class to extend 'get/set'item
    >>> h3 = ...
    >>> v  = h3[1][2 ,3]              
    >>> v  = h3[1][2][3]
    >>> v  = h3[1](0.13, 0.44) ## as 2D-function, use interpolation
    >>> v  = h3[1][2](0.44)    ## as 1D-function, use interpolation
    >>> h3[1][2 ,3] = v1  
    >>> h3[1][2][3] = v2
    """
    def __init__ ( self , h3 , xbin ) :
        self._h3   = h3 
        if not xbin in h3.GetXaxis () : raise IndexError
        self._xbin = xbin
        
    def __getitem__ ( self , ybin  ) :
        if isinstance ( ybin , integer_types ) : 
            return _H3_item_2_ ( self._h3 , self._xbin , ybin )
        elif 2 == len ( ybin ) :
            return self._h3 [ self._xbin , ybin[0] , ybin[1] ] 
        else :
            raise IndexError

    def __setitem__ ( self , bins , value ) :
        if 2 != len ( bins ) : raise IndexError 
        self._h3 [ self._xbin , bins[0] , bins[1] ] = value 

    def __call__    ( self , y , z , func = lambda s : s , interpolate = True )  :
        x = self._h3.GetXaxis().GetBinCenter ( self._xbin )
        return self._h3( x , y , z , func , interpolate ) 
        
# =============================================================================
## get item for the 3D histogram
#  @code
#  histo = ...
#  ve    = histo[ix, iy, iz]
#  ve    = histo[ix][iy][iz] ## ditto 
#  ve    = histo[ix, iy][iz] ## ditto
#  ve    = histo[ix][iy, iz] ## ditto    
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_get_item_ ( h3 , ibin ) :
    """``Get-item'' for the 3D-histogram :    
    >>> histo = ...
    >>> ve    = histo[ix, iy, iz]
    >>> ve    = histo[ix][iy][iz] ## ditto 
    >>> ve    = histo[ix, iy][iz] ## ditto
    >>> ve    = histo[ix][iy, iz] ## ditto    
    """
    #
    if isinstance ( ibin , integer_types ) : return _H3_item_1_ ( h3 ,  ibin )
    elif 2 == len( ibin )                  : return _H3_item_2_ ( h3 , *ibin )

    ix = ibin[0]
    ax = h3.GetXaxis()
    nx = ax.GetNbins()
    #
    if   1<=  ix <= nx : pass
    elif 1<= -ix <= nx : ix += ( nx + 1 ) 
    else               : raise IndexError
    #
    iy = ibin[1]
    ay = h3.GetYaxis()
    ny = ay.GetNbins()
    #
    if   1<=  iy <= ny : pass
    elif 1<= -iy <= ny : iy += ( ny + 1 ) 
    else               : raise IndexError
    #
    iz = ibin[2]
    az = h3.GetZaxis()
    nz = az.GetNbins()
    #
    if   1<=  iz <= nz : pass
    elif 1<= -iz <= nz : iz += ( nz + 1 ) 
    else               : raise IndexError
    #
    val = h3.GetBinContent ( ix , iy , iz )   
    err = h3.GetBinError   ( ix , iy , iz ) 
    #
    return VE ( val , err * err ) 

ROOT.TH1F  . __getitem__  = _h1_get_item_
ROOT.TH1D  . __getitem__  = _h1_get_item_
ROOT.TH2F  . __getitem__  = _h2_get_item_
ROOT.TH2D  . __getitem__  = _h2_get_item_
ROOT.TH3F  . __getitem__  = _h3_get_item_
ROOT.TH3D  . __getitem__  = _h3_get_item_

# =============================================================================
## iterator for 1D-histogram 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_iter_ ( h1 ) :
    """ Iterator over 1D-histogram
    >>> for i in h1 : print (i) 
    """
    ax = h1.GetXaxis () 
    sx = ax.GetNbins ()
    for i in range ( 1 , sx + 1 ) : 
        yield i

# =============================================================================
## iterator for 1D-histogram in reverse order 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-06
def _h1_iter_reversed_ ( h1 ) :
    """ Iterator over 1D-histogram    
    >>> for i in reversed(h1) : print (i)
    """
    ax = h1.GetXaxis () 
    sx = ax.GetNbins ()
    for i in range ( sx , 0 , -1 ) : 
        yield i
        
ROOT.TH1  . __iter__     = _h1_iter_ 
ROOT.TH1F . __iter__     = _h1_iter_ 
ROOT.TH1D . __iter__     = _h1_iter_ 
ROOT.TH1F . __reversed__ = _h1_iter_reversed_ 
ROOT.TH1D . __reversed__ = _h1_iter_reversed_ 

# =============================================================================
## iterator for 2D-histogram 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_iter_ ( h2 ) :
    """ Iterator over 2D-histogram    
    >>> for i in h2 : print (i)
    """
    #
    ax = h2.GetXaxis()
    ay = h2.GetYaxis()
    #
    sx = ax.GetNbins()
    sy = ay.GetNbins()
    #
    for ix in range ( 1 , sx + 1 ) : 
        for iy in range ( 1 , sy + 1 ) : 
            yield ix , iy

# =============================================================================
## iterator for 2D-histogram 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_iter_reversed_ ( h2 ) :
    """Iterator over 2D-histogram    
    >>> for i in reversed(h2) : print (i)
    """
    #
    ax = h2.GetXaxis()
    ay = h2.GetYaxis()
    #
    sx = ax.GetNbins()
    sy = ay.GetNbins()
    #
    for iy in range ( sy , 0 , -1 ) : 
        for ix in range ( sx , 0 , -1  ) : 
            yield ix , iy


ROOT.TH2  . __iter__     = _h2_iter_ 
ROOT.TH2F . __iter__     = _h2_iter_ 
ROOT.TH2D . __iter__     = _h2_iter_ 
ROOT.TH2F . __reversed__ = _h2_iter_reversed_ 
ROOT.TH2D . __reversed__ = _h2_iter_reversed_ 


# =============================================================================
## iterator for 3D-histogram 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_iter_ ( h3 ) :
    """ Iterator over 3D-histogram    
    >>> for i in h3 : print (i)
    """
    #
    ax = h3.GetXaxis()
    ay = h3.GetYaxis()
    az = h3.GetZaxis()
    #
    sx = ax.GetNbins()
    sy = ay.GetNbins()    
    sz = az.GetNbins()
    #
    for ix in range ( 1 , sx + 1 ) : 
        for iy in range ( 1 , sy + 1 ) : 
            for iz in range ( 1 , sz + 1 ) : 
                yield  ix , iy , iz

# =============================================================================
## iterator for 3D-histogram 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_iter_reversed_ ( h3 ) :
    """ Reversed iterator over 3D-histogram
    >>> for i in reversed(h3) : print (i)
    """
    #
    ax = h3.GetXaxis()
    ay = h3.GetYaxis()
    az = h3.GetZaxis()
    #
    sx = ax.GetNbins()
    sy = ay.GetNbins()    
    sz = az.GetNbins()
    #
    for iz in range ( sz , 0 , -1 ) : 
        for iy in range ( sy , 0 , -1  ) : 
            for ix in range ( sx , 0 , -1  ) : 
                yield  ix , iy , iz
                                    
ROOT.TH3  . __iter__     = _h3_iter_ 
ROOT.TH3F . __iter__     = _h3_iter_ 
ROOT.TH3D . __iter__     = _h3_iter_ 
ROOT.TH3F . __reversed__ = _h3_iter_reversed_ 
ROOT.TH3D . __reversed__ = _h3_iter_reversed_ 

# =============================================================================
## C++ function for 1D-histogram interpolation
#  @see Ostap::Math::HistoInterpolation
#  @see Ostap::Math::HistoInterpolation::interpolate_1D
_interpolate_1D_ = Ostap.Math.HistoInterpolation.interpolate_1D    
# =============================================================================
## C++ function for 2D-histogram interpolation
#  @see Ostap::Math::HistoInterpolation
#  @see Ostap::Math::HistoInterpolation::interpolate_2D
_interpolate_2D_ = Ostap.Math.HistoInterpolation.interpolate_2D
# =============================================================================
## C++ function for 3D-histogram interpolation
#  @see Ostap::Math::HistoInterpolation
#  @see Ostap::Math::HistoInterpolation::interpolate_3D
_interpolate_3D_ = Ostap.Math.HistoInterpolation.interpolate_3D
# =============================================================================

# =============================================================================
## histogram as function 
#  @see Ostap::Math::HistoInterpolation
#  @see Ostap::Math::HistoInterpolation::interpolate_1D
#  @param func        (INPUT) apply this function to result of interpolation
#  @param interpolate (INPUT) interpolation rule
#  @param edges       (INPUT) use special treatment of edge bins 
#  @param extrapolate (INPUT) allow extrapolation
#  @param density     (INPUT) interpolate  "density" instead of bin-content, where
#                             "density" is defined as (bin-content)/(bin-volume)
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_call_ ( h1                    ,
                x                     ,
                func = lambda s : s   ,
                interpolate   = True  ,
                edges         = True  ,
                extrapolate   = False ,
                density       = False ) :
    """ Histogram as function:
    >>> histo = ....
    >>> ve    = histo ( x )                       ## default interpolation
    >>> ve    = histo ( x , interpolate = False ) ## no interpolation 
    >>> ve    = histo ( x , interpolate = True  ) ## default interpolation  
    >>> ve    = histo ( x , interpolate = 2     ) ## parabolic interpolation
    >>> ve    = histo ( x , edges = False ) ## no special treatment of edge bins
    >>> ve    = histo ( x , extrapolate = True )  ## allow extrapolation outside histogram range
    >>> ve    = histo ( x , density = True ) ## interpolate for density: (bin-content)/(bin-volume)
    >>> ve    = histo ( x , func = lambda v : v*v ) ## interpolate and square the result 
    """
    #
    x = float ( x )
    #
    tx  = 1 
    if isinstance ( interpolate , int ) and 0 <= interpolate :
        tx = interpolate
    elif not interpolate :
        tx  = 0
        
    ## use C++ function for fast interpolation 
    result = _interpolate_1D_ ( h1 , x , tx , edges , extrapolate , density )
    #
    return func ( result )

# =============================================================================
## histogram as trivial function object
#  @code
#  histo = ...
#  fun   = Histo1DFun( histo )
#  @endcode
#  @author Vanya BELYAEV IvanBelyaev@itep.ru
#  @date 2019-05-14
class Histo1DFun (object) :
    """ Histogram as trivial function object
    >>> histo = ...
    >>> fun   = Histo1DFun( histo )
    """
    def __init__ ( self                     ,
                   h1                       ,
                   func        = lambda : s ,
                   interpolate = True       ,
                   edges       = True       ,
                   extrapolate = False      ,
                   density     = False      ) :
        
        assert isinstance ( h1 , ROOT.TH1 ) and 1 == h1.GetDimension() , 'Invalid histogram type!'
        
        self.__histo  = h1
        self.__config = {
            'func'        : func        ,
            'interpolate' : interpolate ,
            'edges'       : edges       ,
            'extrapolate' : extrapolate ,
            'density'     : density
            }
    ##  the only one important method
    def __call__ ( x ) :
        return self.__histo ( x , **self.__config )
    @property
    def histo ( self )  :
        """`histo': the  histogram itself"""
        return self.__histo 
    @property
    def config ( self ) :
        """`config'  : the configuration (immutable)"""
        cnf = {}
        cnf.update ( self.__config )
        return cnf 

        
ROOT.TH1F  . __call__     = _h1_call_
ROOT.TH1D  . __call__     = _h1_call_

ROOT.TH1   . __len__      = lambda s : s.size() 
ROOT.TH1   .   size       = lambda s : s.GetNbinsX() * s.GetNbinsY() * s.GetNbinsZ() 
ROOT.TH1   . __contains__ = lambda s , i : 1 <= i <= s.size() 

ROOT.TH1F  . __contains__ = lambda s , i : i in s.GetXaxis() 
ROOT.TH1D  . __contains__ = lambda s , i : i in s.GetXaxis() 

ROOT.TH2   . __len__      = lambda s : s.size() 
ROOT.TH2   .   size       = lambda s : s.GetNbinsX() * s.GetNbinsY() * s.GetNbinsZ()

ROOT.TH3   . __len__      = lambda s : s.size() 
ROOT.TH3   .   size       = lambda s : s.GetNbinsX() * s.GetNbinsY() * s.GetNbinsZ()

ROOT.TH1D  . nbins        = lambda s : s.GetNbinsX() 
ROOT.TH1F  . nbins        = lambda s : s.GetNbinsX()
ROOT.TH1D  . nbinsx       = lambda s : s.GetNbinsX() 
ROOT.TH1F  . nbinsx       = lambda s : s.GetNbinsX()
ROOT.TH1D  .  binsx       = lambda s : s.GetNbinsX() 
ROOT.TH1F  .  binsx       = lambda s : s.GetNbinsX()

ROOT.TH1D  .  bins        = lambda s : s.GetNbinsX() 
ROOT.TH1F  .  bins        = lambda s : s.GetNbinsX()
ROOT.TH2D  . nbinsx       = lambda s : s.GetNbinsX()
ROOT.TH2D  . nbinsy       = lambda s : s.GetNbinsY()
ROOT.TH2F  . nbinsx       = lambda s : s.GetNbinsX()
ROOT.TH2F  . nbinsy       = lambda s : s.GetNbinsY()
ROOT.TH2D  .  binsx       = lambda s : s.GetNbinsX()
ROOT.TH2D  .  binsy       = lambda s : s.GetNbinsY()
ROOT.TH2F  .  binsx       = lambda s : s.GetNbinsX()
ROOT.TH2F  .  binsy       = lambda s : s.GetNbinsY()

ROOT.TH3D  . nbinsx       = lambda s : s.GetNbinsX()
ROOT.TH3D  . nbinsy       = lambda s : s.GetNbinsY()
ROOT.TH3D  . nbinsz       = lambda s : s.GetNbinsZ()
ROOT.TH3F  . nbinsx       = lambda s : s.GetNbinsX()
ROOT.TH3F  . nbinsy       = lambda s : s.GetNbinsY()
ROOT.TH3F  . nbinsz       = lambda s : s.GetNbinsZ()
ROOT.TH3D  .  binsx       = lambda s : s.GetNbinsX()
ROOT.TH3D  .  binsy       = lambda s : s.GetNbinsY()
ROOT.TH3D  .  binsz       = lambda s : s.GetNbinsZ()
ROOT.TH3F  .  binsx       = lambda s : s.GetNbinsX()
ROOT.TH3F  .  binsy       = lambda s : s.GetNbinsY()
ROOT.TH3F  .  binsz       = lambda s : s.GetNbinsZ()


# =============================================================================
## check bin in 2D-histo 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_contains_ ( h2 , ibin ) :
    """ Check if the bin contains in 2D-histogram:
    >>> (3,5) in h2    
    """
    return ibin[0] in h2.GetXaxis() and ibin[1] in h2.GetYaxis()

ROOT.TH2   . __contains__ = _h2_contains_
ROOT.TH2F  . __contains__ = _h2_contains_
ROOT.TH2D  . __contains__ = _h2_contains_

# ============================================================================
## check bin in 3D-histo 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_contains_ ( h3 , ibin ) :
    """ Check if the bin contains in 3D-histogram:
    >>> (3,5,10) in h3    
    """
    return ibin[0] in h2.GetXaxis() and \
           ibin[1] in h2.GetYaxis() and \
           ibin[2] in h2.GetZaxis() 

ROOT.TH3   . __contains__ = _h3_contains_
ROOT.TH3F  . __contains__ = _h3_contains_
ROOT.TH3D  . __contains__ = _h3_contains_


# =============================================================================
## number of "empty" bins
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-20
def _num_empty_ ( h ) :
    """Check number of empty bins :
    >>> h = ...
    >>> e = h.numEmpty()
    """
    ne  = 0 
    for i in h.items() : 
        v = i[-1]
        if iszero ( v.value() ) and iszero  ( v.cov2 () ) : ne +=1
    return ne 

ROOT.TH1 . numEmpty = _num_empty_

# =============================================================================
## find bin in 1D-histogram
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_find_ ( h1 , x ) :
    """ Find the bin in 1D-histogram
    >>> ibin = h1.findBin ( x ) 
    """
    x = float ( x )
    #
    ax = h1.GetXaxis()
    #
    return ax.FindBin ( x )
# =============================================================================
## find bin in 2D-histogram
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_find_ ( h2 , x , y ) :
    """Find the bin in 3D-histogram
    >>> ibin = h2.findBin ( x , y ) 
    """
    x = float ( x )
    y = float ( y )
    #
    ax = h2.GetXaxis()
    ay = h2.GetYaxis()
    #
    return  ( ax.FindBin ( x ) ,
              ay.FindBin ( y ) )
# =============================================================================
## find bin in 3D-histogram
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_find_ ( h3 , x , y , z ) :
    """ Find the bin in 3D-histogram
    >>> ibin = h3.findBin ( x , y , z ) 
    """
    x = float ( x )
    y = float ( y )
    z = float ( z )
    #
    ax = h3.GetXaxis()
    ay = h3.GetYaxis()
    az = h3.GetZaxis()
    #
    return  ( ax.FindBin ( x ) ,
              ay.FindBin ( y ) ,
              az.FindBin ( z ) )


ROOT.TH1F . findBin  = _h1_find_
ROOT.TH1D . findBin  = _h1_find_
ROOT.TH2F . findBin  = _h2_find_
ROOT.TH2D . findBin  = _h2_find_
ROOT.TH3F . findBin  = _h3_find_
ROOT.TH3D . findBin  = _h3_find_

# =============================================================================
## get mean for 2D-histogram 
#  @code 
#  >>> histo = ...
#  >>> x , y = histo.mean ()
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-03-20
def _h2_mean_ ( h2 ) :
    """ Get the mean
    >>> histo = ...
    >>> x,y   = histo.mean () 
    """
    x  = h2.GetMean      ( 1 )
    ex = h2.GetMeanError ( 1 )
    y  = h2.GetMean      ( 2 )
    ey = h2.GetMeanError ( 2 )
    #
    return VE ( x , ex * ex ), VE ( y , ey * ey ) 

# ============================================================================
## Get x-mean for 2D-histogram
#  @code
#  h2 = ...
#  xm = h2.xmean() 
#  @endcode
def _h2_xmean_ ( h2 ) :
    """ Get x-mean for 2D-histogram
    >>> h2 = ...
    >>> xm = h2.xmean() 
    """
    x  = h2.GetMean      ( 1 )
    ex = h2.GetMeanError ( 1 )
    #
    return VE ( x , ex * ex )

# ============================================================================
## Get y-mean for 2D-histogram
#  @code
#  h2 = ...
#  xm = h2.ymean() 
#  @endcode
def _h2_ymean_ ( h2 ) :
    """ Get y-mean for 2D-histogram
    >>> h2 = ...
    >>> xm = h2.ymean() 
    """
    y  = h2.GetMean      ( 2 )
    ey = h2.GetMeanError ( 2 )
    #
    return VE ( y , ey * ey )

# ============================================================================
## Get x-rms for 2D-histogram
#  @code
#  h2 = ...
#  xm = h2.xrms() 
#  @endcode
def _h2_xrms_ ( h2 ) :
    """ Get x-rms for 2D-histogram
    >>> h2 = ...
    >>> xm = h2.xrms() 
    """

    order = 2
    
    mu2    = h2.central_moment ( 2 , 0 )
    mu4    = h2.central_moment ( 4 , 0 )
    ##
    cov2   = mu4
    cov2  -= mu2 * mu2 
    cov2  /= h2.nEff()
    ##
    cov2   = max ( cov2 , 0.0 )
    ##
    return VE ( mu2 , cov2 ) ** 0.5 

# ============================================================================
## Get y-rms for 2D-histogram
#  @code
#  h2 = ...
#  xm = h2.yrms() 
#  @endcode
def _h2_yrms_ ( h2 ) :
    """ Get y-rms for 2D-histogram
    >>> h2 = ...
    >>> xm = h2.yrms() 
    """
    order = 2
    
    mu2    = h2.central_moment ( 0 , 2 )
    mu4    = h2.central_moment ( 0 , 4 )
    ##
    cov2   = mu4
    cov2  -= mu2 * mu2 
    cov2  /= h2.nEff()
    ##
    cov2   = max ( cov2 , 0.0 )
    ##
    return VE ( mu2 , cov2 ) ** 0.5 

# ===============================================================================
## Get RMS for x&y-axes
#  @code
#  xrms , yrms = h2.rms() 
#  @endcode
def _h2_rms_ ( h2 ) :
    """ Get RMS for x&y-axes
    >>> h2 = ...
    >>> xrms , yrms = h2.rms() 
    """
    return h2.xrms() , h2.yrms()

# ===============================================================================
## get correlation coefficient between x and y
#  @code
#  h2 = ...
#  xy = h2.xycorr()
#  @endcode 
def _h2_xycorr_  ( h2 ) :
    """ Get correlation coefficient between x and y
    >>> h2 = ...
    >>> xy = h2.xycorr()
    """
    c2 = h2.central_moment ( 1 , 1 )
    x2 = h2.central_moment ( 2 , 0 )
    y2 = h2.central_moment ( 0 , 2 )
    ## 
    return c2 / ( x2 * y2 ) ** 0.5


for h in ( ROOT.TH2F , ROOT.TH2D ) :
    h. mean  =  _h2_mean_
    h.xmean  =  _h2_xmean_
    h.ymean  =  _h2_ymean_
    h.xrms   =  _h2_xrms_
    h.yrms   =  _h2_yrms_
    h. rms   =  _h2_rms_
    h.xycorr =  _h2_xycorr_
    h.xycor  =  _h2_xycorr_


# =============================================================================
## get mean for 3D-histogram 
#  @code 
#  >>> histo = ...
#  >>> x , y , z = histo.mean ()
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-03-20
def _h3_mean_ ( h3 ) :
    """ Get the mean
    >>> histo = ...
    >>> x , y , z = histo.mean () 
    """
    x  = h3.GetMean      ( 1 )
    ex = h3.GetMeanError ( 1 )
    y  = h3.GetMean      ( 2 )
    ey = h3.GetMeanError ( 2 )
    z  = h3.GetMean      ( 3 )
    ez = h3.GetMeanError ( 4 )
    #
    return VE ( x , ex * ex ), VE ( y , ey * ey ) , VE ( z , ez * ez )

# ============================================================================
## Get x-mean for 3D-histogram
#  @code
#  h3 = ...
#  xm = h3.xmean() 
#  @endcode
def _h3_xmean_ ( h3 ) :
    """ Get x-mean for 2D-histogram
    >>> h2 = ...
    >>> xm = h2.xmean() 
    """
    x  = h3.GetMean      ( 1 )
    ex = h3.GetMeanError ( 1 )
    #
    return VE ( x , ex * ex )

# ============================================================================
## Get y-mean for 3D-histogram
#  @code
#  h3 = ...
#  ym = h3.ymean() 
#  @endcode
def _h3_ymean_ ( h3 ) :
    """ Get y-mean for 3D-histogram
    >>> h3 = ...
    >>> xm = h3.ymean() 
    """
    y  = h3.GetMean      ( 2 )
    ey = h3.GetMeanError ( 2 )
    #
    return VE ( y , ey * ey )


# ============================================================================
## Get z-mean for 3D-histogram
#  @code
#  h3 = ...
#  zm = h3.zmean() 
#  @endcode
def _h3_zmean_ ( h3 ) :
    """ Get z-mean for 3D-histogram
    >>> h3 = ...
    >>> zm = h3.zmean() 
    """
    z  = h3.GetMean      ( 3 )
    ez = h3.GetMeanError ( 3 )
    #
    return VE ( z , ez * ez )

# ============================================================================
## Get x-rms for 3D-histogram
#  @code
#  h3 = ...
#  xm = h3.xrms() 
#  @endcode
def _h3_xrms_ ( h2 ) :
    """ Get x-rms for 3D-histogram
    >>> h3 = ...
    >>> xm = h3.xrms() 
    """

    order = 2
    
    mu2    = h2.central_moment ( 2 , 0 , 0 )
    mu4    = h2.central_moment ( 4 , 0 , 0)
    ##
    cov2   = mu4
    cov2  -= mu2 * mu2 
    cov2  /= h2.nEff()
    ##
    cov2   = max ( cov2 , 0.0 )
    ##
    return VE ( mu2 , cov2 ) ** 0.5 

# ============================================================================
## Get y-rms for 3D-histogram
#  @code
#  h3 = ...
#  ym = h3.yrms() 
#  @endcode
def _h3_yrms_ ( h2 ) :
    """ Get y-rms for 3D-histogram
    >>> h3 = ...
    >>> ym = h3.yrms() 
    """
    order = 2
    
    mu2    = h2.central_moment ( 0 , 2 , 0 )
    mu4    = h2.central_moment ( 0 , 4 , 0 )
    ##
    cov2   = mu4
    cov2  -= mu2 * mu2 
    cov2  /= h2.nEff()
    ##
    cov2   = max ( cov2 , 0.0 )
    ##
    return VE ( mu2 , cov2 ) ** 0.5 


# ============================================================================
## Get z-rms for 3D-histogram
#  @code
#  h3 = ...
#  zm = h3.zrms() 
#  @endcode
def _h3_zrms_ ( h3 ) :
    """ Get z-rms for 3D-histogram
    >>> h3 = ...
    >>> zm = h3.zrms() 
    """
    order = 2
    ## 
    mu2    = h3.central_moment ( 0 , 0 , 2 )
    mu4    = h3.central_moment ( 0 , 0 , 4 )
    ##
    cov2   = mu4
    cov2  -= mu2 * mu2 
    cov2  /= h3.nEff()
    ##
    cov2   = max ( cov2 , 0.0 )
    ##
    return VE ( mu2 , cov2 ) ** 0.5 

# ===============================================================================
## Get RMS for x,y&z-axes
#  @code
#  xrms , yrms , zrms = h3.rms() 
#  @endcode
def _h3_rms_ ( h3 ) :
    """ Get RMS for x,y&z-axes
    >>> h2 = ...
    >>> xrms , yrms . zrms = h3.rms() 
    """
    return h3.xrms() , h3.yrms() , h3.zrms()

# ===============================================================================
## get correlation coefficient between x and y
#  @code
#  h3 = ...
#  xy = h3.xycorr()
#  @endcode 
def _h3_xycorr_  ( h3 ) :
    """ Get correlation coefficient between x and y
    >>> h3 = ...
    >>> xy = h3.xycorr()
    """
    c2 = h3.central_moment ( 1 , 1 , 0 )
    x2 = h3.central_moment ( 2 , 0 , 0 )
    y2 = h3.central_moment ( 0 , 2 , 0 )
    ## 
    return c2 / ( x2 * y2 ) ** 0.5

# ===============================================================================
## get correlation coefficient between y and z
#  @code
#  h3 = ...
#  yz = h3.yzcorr()
#  @endcode 
def _h3_yzcorr_  ( h3 ) :
    """ Get correlation coefficient between y and z 
    >>> h3 = ...
    >>> xy = h3.yzcorr()
    """
    c2 = h3.central_moment ( 0 , 1 , 1 )
    x2 = h3.central_moment ( 0 , 2 , 0 )
    y2 = h3.central_moment ( 0 , 0 , 2 )
    ## 
    return c2 / ( x2 * y2 ) ** 0.5

# ===============================================================================
## get correlation coefficient between x and z
#  @code
#  h3 = ...
#  xz = h3.xzcorr()
#  @endcode 
def _h3_xzcorr_  ( h3 ) :
    """ Get correlation coefficient between x and z 
    >>> h3 = ...
    >>> xz = h3.xzcorr()
    """
    c2 = h3.central_moment ( 1 , 0 , 1 )
    x2 = h3.central_moment ( 2 , 0 , 0 )
    y2 = h3.central_moment ( 0 , 0 , 2 )
    ## 
    return c2 / ( x2 * y2 ) ** 0.5

for h in ( ROOT.TH3F , ROOT.TH3D ) :
    
    h. mean  =  _h3_mean_
    h.xmean  =  _h3_xmean_
    h.ymean  =  _h3_ymean_
    h.zmean  =  _h3_zmean_
    h.xrms   =  _h3_xrms_
    h.yrms   =  _h3_yrms_
    h.zrms   =  _h3_zrms_
    h.rms    =  _h3_rms_
    
    h.xycorr = _h3_xycorr_
    h.yxcorr = _h3_xycorr_
    h.xycor  = _h3_xycorr_
    h.yxcor  = _h3_xycorr_
    
    h.yzcorr = _h3_yzcorr_
    h.zycorr = _h3_yzcorr_
    h.yzcor  = _h3_yzcorr_
    h.zycor  = _h3_yzcorr_
    
    h.xzcorr = _h3_xzcorr_
    h.zycorr = _h3_xzcorr_
    h.xzcor  = _h3_xzcorr_
    h.zxcor  = _h3_xzcorr_

# ============================================================================
## find the first X-value for the certain Y-value
#  Actually solve the equation
#  \f$ y = f(x) \f$ using linear interpolation 
#  @code
#  histo = ...
#  x = histo.find_X( 100 ) 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-03-20
def _h1_find_X ( self             ,
                 v                ,
                 forward   = True ) : 
    """ Find the first x-value for certain Y-value
    >>> histo = ...
    >>> x = histo.find_X( 100 ) 
    Actually it solves the equation: \f$ y = f(x) \f$ using linear interpolation 
    """
    ##
    v = float ( v ) 
    ##
    mn , mx  = self.minmax()
    ##
    if v < mn.value() :
        return hist.xmin() if forward else hist.xmax ()  
    if v > mx.value() :
        return hist.xmax() if forward else hist.xmin () 
    
    nb  = len ( hist ) 
    
    ax  = hist.GetXaxis()
    
    for i in hist :
        
        j  = i+1
        if not j in hist : continue
        
        ib = i if forward else nb + 1 - i 
        jb = j if forward else nb + 1 - j 
        
        vi = hist [ ib ].value () 
        vj = hist [ jb ].value ()
        
        if  vi <= v <= vj or vj <= v <= vi :
            
            xi = ax.GetBinCenter ( ib )
            xj = ax.GetBinCenter ( jb )
            
            if   vi == v             : return xi
            elif vj == v             : return xj 
            
            dv = vi - vj
            dx = xi - xj

            if   vi == vj or 0 == dv : return 0.5 * ( xi + xj ) 
            
            return (v*dx+vi*xj-vj*xi)/dv


# =============================================================================
## histogram as function
#  @see Ostap::Math::HistoInterpolation
#  @see Ostap::Math::HistoInterpolation::interpolate_2D
#  @param func        (INPUT) apply this function to result of interpolation
#  @param interpolate (INPUT) interpolation rule
#  @param edges       (INPUT) use special treatment of edge bins 
#  @param extrapolate (INPUT) allow extrapolation
#  @param density     (INPUT) interpolate  "density" instead of bin-content, where
#                             "density" is defined as (bin-content)/(bin-volume)
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_call_ ( h2 ,
                x  ,
                y                          ,
                func        = lambda s : s ,
                interpolate = True         ,
                edges       = True         ,
                extrapolate = False        ,
                density     = False        ) :
    """ 2D-histogram as function:
    >>> histo = ....
    >>> ve    = histo ( x , y ) ## default interpolation
    >>> ve    = histo ( x , y , interpolate = (3,3) )   ## use bi-qubic interpolation
    >>> ve    = histo ( x , y , interpolate = False )   ## no interpolation
    >>> ve    = histo ( x , y , func = lambda v : v*v ) ## interpolate and square the result
    >>> ve    = histo ( x , y , edges = False )         ## no special treatment of edge bins
    >>> ve    = histo ( x , y , extrapolate =  True ) ## allow extrapolation outside histogram range
    >>> ve    = histo ( x , y , density = True ) ## interpolate for density: (content)/(volume)  
    """    
    #
    x = float ( x )
    y = float ( y )
    #
    tx = 1
    ty = 1 
    if   not interpolate :
        tx = 0
        ty = 0
    elif isinstance ( interpolate , (tuple,list) ) and 2<= len ( interpolate ) : 
        tx = int ( interpolate[0] )
        ty = int ( interpolate[1] )
        
    ## use C++ function for fast interpolation 
    result = _interpolate_2D_ ( h2 , x , y , tx , ty , edges , extrapolate , density )

    return func ( result )

# =============================================================================
## histogram as trivial function object
#  @code
#  histo = ...
#  fun   = Histo2DFun( histo )
#  @endcode
#  @author Vanya BELYAEV IvanBelyaev@itep.ru
#  @date 2019-05-14
class Histo2DFun (object) :
    """ Histogram as trivial function object
    >>> histo = ...
    >>> fun   = Histo2DFun( histo )
    """
    def __init__ ( self                     ,
                   h2                       ,
                   func        = lambda : s ,
                   interpolate = True       ,
                   edges       = True       ,
                   extrapolate = False      ,
                   density     = False      ) :
        
        assert isinstance ( h2 , ROOT.TH2 ) and 2 == h2.GetDimension() , 'Invalid histogram type!'
        
        self.__histo  = h2
        self.__config = {
            'func'        : func        ,
            'interpolate' : interpolate ,
            'edges'       : edges       ,
            'extrapolate' : extrapolate ,
            'density'     : density
            }
    ##  the only one imporant method
    def __call__ ( x , y ) :
        return self.__histo ( x , y , **self.__config )
    @property
    def histo ( self )  :
        """``histo'': the  histogram itself"""
        return self.__histo 
    @property
    def config ( self ) :
        """``config'' : the configuration (immutable)"""
        cnf = {}
        cnf.update ( self.__config )
        return cnf 



ROOT.TH2   . __call__     = _h2_call_
ROOT.TH2F  . __getitem__  = _h2_get_item_
ROOT.TH2D  . __getitem__  = _h2_get_item_


# =============================================================================
## histogram as function 
#  @see Ostap::Math::HistoInterpolation
#  @see Ostap::Math::HistoInterpolation::interpolate_3D
#  @param func        (INPUT) apply this function to result of interpolation
#  @param interpolate (INPUT) interpolation rule
#  @param edges       (INPUT) use special treatment of edge bins 
#  @param extrapolate (INPUT) allow extrapolation
#  @param density     (INPUT) interpolate  "density" instead of bin-content, where
#                             "density" is defined as (bin-content)/(bin-volume)
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_call_ ( h3                         ,
                x                          ,
                y                          ,
                z                          ,
                func        = lambda s : s ,
                interpolate = True         ,
                edges       = True         ,
                extrapolate = False        ,
                density     = False        ) :                
    """ 3D-histogram as function:
    >>> histo = ....
    >>> ve = histo ( x , y , z ) ## default interpolation
    >>> ve = histo ( x , y , z , interpolate = (3,3,3) )   ## use tri-qubic interpolation
    >>> ve = histo ( x , y , z , interpolate = False )   ## no interpolation
    >>> ve = histo ( x , y , z , func = lambda v : v*v ) ## interpolate and square the result
    >>> ve = histo ( x , y , z , edges = False )         ## no special treatment of edge bins
    >>> ve = histo ( x , y , z , extrapolate =  True ) ## allow extrapolation outside histogram range
    >>> ve = histo ( x , y , z , density = True ) ## interpolate for density: (content)/(volume)  
    """
    #
    x = float ( x )
    y = float ( y )
    z = float ( z )
    #
    tx = 1
    ty = 1 
    tz = 1 
    if   not interpolate :
        tx = 0
        ty = 0
        tz = 0
    elif isinstance ( interpolate , ( tuple , list ) ) and 3 <= len ( interpolate ) :
        tx = int ( interpolate [ 0 ] ) 
        ty = int ( interpolate [ 1 ] ) 
        tz = int ( interpolate [ 2 ] )

    ## use C++ function for fast interpolation 
    result = _interpolate_3D_ ( h3 , x , y , z , tx , ty , tz , edges , extrapolate , density )

    return func ( result )


# =============================================================================
## histogram as trivial function object
#  @code
#  histo = ...
#  fun   = Histo3DFun( histo )
#  @endcode
#  @author Vanya BELYAEV IvanBelyaev@itep.ru
#  @date 2019-05-14
class Histo3DFun (object) :
    """ Histogram as trivial function object
    >>> histo = ...
    >>> fun   = Histo3DFun( histo )
    """
    def __init__ ( self                     ,
                   h3                       ,
                   func        = lambda : s ,
                   interpolate = True       ,
                   edges       = True       ,
                   extrapolate = False      ,
                   density     = False      ) :
        
        assert isinstance ( h3 , ROOT.TH3 ) and 3 == h3.GetDimension() , 'Invalid histogram type!'
        
        self.__histo  = h3
        self.__config = {
            'func'        : func        ,
            'interpolate' : interpolate ,
            'edges'       : edges       ,
            'extrapolate' : extrapolate ,
            'density'     : density
            }
    ##  the only one imporant method
    def __call__ ( x , y , z ) :
        return self.__histo ( x , y , z , **self.__config )
    @property
    def histo ( self )  :
        """``histo'': the  histogram itself"""
        return self.__histo 
    @property
    def config ( self ) :
        """``config'' : the configuration (immutable)"""
        cnf = {}
        cnf.update ( self.__config )
        return cnf 

ROOT.TH3   . __call__     = _h3_call_
ROOT.TH3F  . __getitem__  = _h3_get_item_
ROOT.TH3D  . __getitem__  = _h3_get_item_

# =============================================================================
# iterate over items
# =============================================================================
## iterate over entries in 1D-histogram 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_iteritems_ ( h1 , low = 1 , high = sys.maxsize ) :
    """ Iterate over histogram items:    
    >>> h1 = ...
    >>> for i,x,y in h1.    items()  : ...    
    >>> for i,x,y in h1.iteritems()  : ...    ## ditto 
    """
    ax = h1.GetXaxis()
    sx = ax.GetNbins()
    if low  < 1   : low = 1 
    
    for ix in range ( max (      1 , low  )  ,
                      min ( sx + 1 , high ) ) : 
        
        x   =       ax.GetBinCenter ( ix )
        xe  = 0.5 * ax.GetBinWidth  ( ix )
        
        y   =       h1.GetBinContent ( ix )
        ye  =       h1.GetBinError   ( ix )
        
        yield ix, VE(x,xe*xe) , VE ( y,ye*ye)
        

ROOT.TH1F  .     items     = _h1_iteritems_
ROOT.TH1D  .     items     = _h1_iteritems_
ROOT.TH1F  . iteritems     = _h1_iteritems_
ROOT.TH1D  . iteritems     = _h1_iteritems_


# =============================================================================
## Iterate over the values
#  @cdoe
#  histo = ...
#  for v in histo.values() : ...
#  @endcode 
def _h_values_ ( histo ) :
    """ Iterate over the values
    >>> histo = ...
    >>> for v in histo.values() : ...
    """
    for item in histo.items() :
        yield item[-1] 

ROOT.TH1F  . itervalues    = _h_values_
ROOT.TH1D  . itervalues    = _h_values_
ROOT.TH1F  .     values    = _h_values_
ROOT.TH1D  .     values    = _h_values_
ROOT.TH2F  . itervalues    = _h_values_
ROOT.TH2D  . itervalues    = _h_values_
ROOT.TH2F  .     values    = _h_values_
ROOT.TH2D  .     values    = _h_values_
ROOT.TH3F  . itervalues    = _h_values_
ROOT.TH3D  . itervalues    = _h_values_
ROOT.TH3F  .     values    = _h_values_
ROOT.TH3D  .     values    = _h_values_

_new_methods_ += [
    ROOT.TH1F  . itervalues  , 
    ROOT.TH1D  . itervalues  , 
    ROOT.TH1F  .     values  , 
    ROOT.TH1D  .     values  , 
    ROOT.TH2F  . itervalues  , 
    ROOT.TH2D  . itervalues  , 
    ROOT.TH2F  .     values  , 
    ROOT.TH2D  .     values  , 
    ROOT.TH3F  . itervalues  , 
    ROOT.TH3D  . itervalues  , 
    ROOT.TH3F  .     values  , 
    ROOT.TH3D  .     values  , 
    ]

# =============================================================================
## return information about the bin center and width
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-07-27
def _h1_bin_ ( h1 , ibin ) :
    """ Get the information about the bin center and width     
    >>> h1 = ... ## the histo
    >>> print ( h1.bin(1) ) 
    """
    if not ibin in h1.GetXaxis() : raise IndexError

    ax  = h1.GetXaxis()
    x   =       ax.GetBinCenter ( ibin )
    xe  = 0.5 * ax.GetBinWidth  ( ibin )

    return VE(x,xe*xe)


# =============================================================================
## return information about the bin center and width
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-07-27
def _h2_bin_ ( h2 , *ibin ) :
    """ Get the information about the bin center and width 
    >>> h2 = ... ## the histo
    >>> print ( h2.bin(1,2) ) 
    """
    if 2 != len ( ibin ) : raise IndexError 
    if not ibin in h2    : raise IndexError
    
    ax  = h2.GetXaxis()
    x   =       ax.GetBinCenter ( ibin[0] )
    xe  = 0.5 * ax.GetBinWidth  ( ibin[0] )

    ay  = h2.GetYaxis()
    y   =       ay.GetBinCenter ( ibin[1] )
    ye  = 0.5 * ay.GetBinWidth  ( ibin[1] )

    return VE(x,xe*xe),VE(y,ye*ye)


# =============================================================================
## return information about the bin center and width
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-07-27
def _h3_bin_ ( h3 , *ibin ) :
    """ Get the information about the bin center and width 
    >>> h3 = ... ## the histo
    >>> print ( h3.bin(1,2,15) ) 
    """
    if 3 != len ( ibin ) : raise IndexError 
    if not ibin in h3    : raise IndexError
    
    ax  = h3.GetXaxis()
    x   =       ax.GetBinCenter ( ibin[0] )
    xe  = 0.5 * ax.GetBinWidth  ( ibin[0] )

    ay  = h3.GetYaxis()
    y   =       ay.GetBinCenter ( ibin[1] )
    ye  = 0.5 * ay.GetBinWidth  ( ibin[1] )

    az  = h3.GetZaxis()
    z   =       az.GetBinCenter ( ibin[2] )
    ze  = 0.5 * az.GetBinWidth  ( ibin[2] )

    return VE(x,xe*xe),VE(y,ye*ye),VE(z,ze*ze) 

ROOT.TH1F.bin = _h1_bin_
ROOT.TH1D.bin = _h1_bin_
ROOT.TH2F.bin = _h2_bin_
ROOT.TH2D.bin = _h2_bin_
ROOT.TH3F.bin = _h3_bin_
ROOT.TH3D.bin = _h3_bin_


# =============================================================================
## iterate over entries in 2D-histogram 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_iteritems_ ( h2 ) :
    """ Iterate over histogram items:
    >>> h2 = ...
    >>> for ix,iy,x,y,z in h2.    items() : ...   
    >>> for ix,iy,x,y,z in h2.iteritems() : ... ## ditto   
    """
    ax = h2.GetXaxis()
    sx = ax.GetNbins()

    ay = h2.GetYaxis()
    sy = ay.GetNbins()

    for ix in range( 1 , sx + 1 ) :  
        x   =       ax.GetBinCenter ( ix )
        xe  = 0.5 * ax.GetBinWidth  ( ix )
        #
        for iy in range ( 1 , sy + 1 ) : 
            #
            y   =       ay.GetBinCenter  ( iy      )
            ye  = 0.5 * ay.GetBinWidth   ( iy      )
            #
            z   =       h2.GetBinContent ( ix , iy )
            ze  =       h2.GetBinError   ( ix , iy )
            #
            yield ix, iy, VE(x,xe*xe) , VE ( y,ye*ye) , VE( z,ze*ze) 
            #
            
        
ROOT.TH2F  .     items     = _h2_iteritems_
ROOT.TH2D  .     items     = _h2_iteritems_
ROOT.TH2F  . iteritems     = _h2_iteritems_
ROOT.TH2D  . iteritems     = _h2_iteritems_

# =============================================================================
## iterate over entries in 3D-histogram 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_iteritems_ ( h3 ) :
    """ Iterate over histogram items:    
    >>> h3 = ...
    >>> for ix,iy,iz,x,y,z,V in h3.    items() : ...    
    >>> for ix,iy,iz,x,y,z,V in h3.iteritems() : ...   ##  ditto
    """
    ax = h3.GetXaxis ()
    sx = ax.GetNbins ()

    ay = h3.GetYaxis ()
    sy = ay.GetNbins ()

    az = h3.GetZaxis ()
    sz = az.GetNbins ()
    
    for ix in range( 1 , sx + 1 ) :
        #
        x   =       ax.GetBinCenter ( ix )
        xe  = 0.5 * ax.GetBinWidth  ( ix )
        #
        for iy in range ( 1 , sy + 1 ) : 
            #
            y   =       ay.GetBinCenter  ( iy      )
            ye  = 0.5 * ay.GetBinWidth   ( iy      )
            #
            for iz in range ( 1 , sz + 1 ) :
                #
                z   =       az.GetBinCenter  ( iz      )
                ze  = 0.5 * az.GetBinWidth   ( iz      )
                #
                
                v   =       h3.GetBinContent ( ix , iy , iz )
                ve  =       h3.GetBinError   ( ix , iy , iz )
                #
                yield ix, iy, iz , VE(x,xe*xe) , VE ( y,ye*ye) , VE( z,ze*ze) , VE( v, ve*ve) 
                #

            
ROOT.TH3F  .     items     = _h3_iteritems_
ROOT.TH3D  .     items     = _h3_iteritems_
ROOT.TH3F  . iteritems     = _h3_iteritems_
ROOT.TH3D  . iteritems     = _h3_iteritems_


# =============================================================================
## iterate over items in TAxis
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _a_iteritems_ ( axis ) :
    """ Iterate over items in axis     
    >>> axis = ...
    >>> for ibin,low,center,high in axis.    items() : ...    
    >>> for ibin,low,center,high in axis.iteritems() : ... ## ditto    
    """
    for i in axis :

        l = axis.GetBinLowEdge ( i )
        c = axis.GetBinCenter  ( i )
        u = axis.GetBinUpEdge  ( i )

        yield i,l,c,u

ROOT.TAxis.     items     = _a_iteritems_
ROOT.TAxis. iteritems     = _a_iteritems_

# =============================================================================
## get bin parameters : low- and up-edges 
def _a_get_item_ ( axis , i ) :
    """ Get bin parameter: low- and up-edges
    >>> axis = ...
    >>> low,high = axis[1]
    
    """
    if not i in axis : raise IndexError
    
    l = axis.GetBinLowEdge ( i )
    u = axis.GetBinUpEdge  ( i )
    
    return l,u

ROOT.TAxis.__getitem__  = _a_get_item_


# =============================================================================
## equality for axes
def _a_equal_ ( axis , another ) :
    """ Equality for two axes
    >>> a1 = ...
    >>> a2 = ...
    >>> print ( a1 == a1 ) 
    """
    if      axis   is       another                        : return True 
    if len( axis ) != len ( another )                      : return False
    if not isequal ( axis.GetXmin()  , another.GetXmin() ) : return False
    if not isequal ( axis.GetXmax()  , another.GetXmax() ) : return False
    for i in axis :
        l1,u1 = axis    [i]
        l2,u2 = another [i]
        if not isequal ( l1 , l2 ) : return False 
        if not isequal ( u1 , u2 ) : return False
    return True

ROOT.TAxis.__eq__ =                  _a_equal_
ROOT.TAxis.__ne__ = lambda a,o : not _a_equal_ ( a , o ) 

# =============================================================================
# some minor decoration for 2D-histos 
# =============================================================================

# =============================================================================
## transpose 2D-ihstogram 
def _h2_transpose_ ( h2 ) :
    """``Transpose'' 2D-histogram
    >>> h2 = ...
    >>> ht = h2.T()
    """
    if hasattr ( h2 , 'GetSumw2' ) and not h2.GetSumw2() : h2.Sumw2()
    #
    xa = h2.GetXaxis()
    ya = h2.GetYaxis()
    #
    hn = h2_axes ( ya , xa )  
    ##
    for i in h2 :
        hn [ (i[1],i[0]) ] = h2[ i ] 

    return hn

ROOT.TH2F.T         = _h2_transpose_
ROOT.TH2D.T         = _h2_transpose_
ROOT.TH2F.Transpose = _h2_transpose_
ROOT.TH2D.Transpone = _h2_transpose_
ROOT.TH2F.transpose = _h2_transpose_
ROOT.TH2D.transpone = _h2_transpose_


# =============================================================================
## Draw 2D-histogram as 'colz'
def _h2_colz_ ( h2 , opts = '' , *args , **kwargs ) :
    """ Draw 2D-histogram as 'colz'    
    >>> h2.colz()    
    """
    return h2.draw ( 'colz ' + opts , *args , **kwargs )


ROOT.TH2F . colz = _h2_colz_
ROOT.TH2D . colz = _h2_colz_
ROOT.TH2F . Colz = _h2_colz_
ROOT.TH2D . Colz = _h2_colz_

## add it also for TF2 
ROOT.TF2  . colz = _h2_colz_

# =============================================================================
## Draw 2D-histogram as 'text'
def _h2_text_ ( h2 , opts = '' , fmt = '' , *args , **kwargs ) :
    """ Draw 2D-histogram as 'text'    
    >>> h2.text( fmt = '5.3f')    
    """
    if fmt : ROOT.gStyle.SetPaintTextFormat ( fmt )        
    return h2.draw ( 'text ' + opts , *args , **kwargs )
    
ROOT.TH2F . text = _h2_text_
ROOT.TH2D . text = _h2_text_

# =============================================================================
## Draw 2D-histogram as 'text'
def _h2_texte_ ( h2 , opts = '' , fmt = '' , *args , **kwargs ) :
    """ Draw 2D-histogram as 'texte'    
    >>> h2.texte( fmt = '5.2f')    
    """
    if fmt : ROOT.gStyle.SetPaintTextFormat ( fmt )        
    return h2.draw ( 'texte ' + opts , *args , **kwargs )
    
ROOT.TH2F . texte = _h2_texte_
ROOT.TH2D . texte = _h2_texte_


# =============================================================================
## calculate the efficiency histogram using the binomial errors
#  @code 
#  >>> accepted   = ...
#  >>> total      = ...
#  >>> efficiency = acepted // total
#  @endcode
#  @see Ostap::Math::binomEff 
#  @see Ostap::Math::wilsonEff 
#  @see Ostap::Math::agrestiCoullEff 
#  @see https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def binomEff_h1 ( h1 , h2 , func = binomEff ) :
    """ Calculate the efficiency histogram using the binomial errors    
    >>> accepted   = ...
    >>> total      = ...
    >>> efficiency = accepted // total                   ## default
    >>> effB       = accepted.       binomEff ( total )  ## default
    >>> effAC      = accepted.agrestiCoullEff ( total ) 
    >>> effW       = accepted.      wilsonEff ( total ) 
    - see Ostap::Math::binomEff 
    - see Ostap::Math::wilsonEff 
    - see Ostap::Math::agrestiCoullEff 
    - see https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    """
    #
    if isinstance ( h1 , ROOT.TProfile ) :
        hh = h1.asH1()
        return binomEff_h1 ( hh , h2 , func )
                    
    if                                 not h1.GetSumw2() : h1.Sumw2()
    if hasattr ( h2 , 'GetSumw2' ) and not h2.GetSumw2() : h2.Sumw2()
    #
    result = h1.Clone( hID () )
    if not result.GetSumw2() : result.Sumw2()
    #
    for i1,x1,y1 in h1.items() :
        #
        if not natural_entry ( y1 ) :
            raise AssertionError ( 'Non-natural entry[%d]=%s' % ( ix1 , y1 ) ) 
        #
        y2 = h2 ( x1.value() ) 
        if not natural_entry ( y2 ) :
            raise AssertionError ( 'Non-natural entry(%s)=%s' % ( x1.value() , y2 ) )
        #
        l1 = long_type ( y1.value () )
        l2 = long_type ( y2.value () )
        #
        if  0 > l1 or 0 > l2 or l1 > l2 :
            raise AssertionError ( 'Incorrect values(%s): %d,%d' % ( x1.value() , l1 , l2 ) )
        #
        v = VE ( func ( l1 , l2 ) ) 
        #
        if not v.isfinite() : continue 
        #
        result.SetBinContent ( i1 , v.value () ) 
        result.SetBinError   ( i1 , v.error () )
        
    return result 


ROOT.TH1F.       binomEff  = binomEff_h1 
ROOT.TH1D.       binomEff  = binomEff_h1 
ROOT.TH1F.       wilsonEff = lambda haccepted,htotal : binomEff_h1 ( haccepted , htotal, func =       wilsonEff )
ROOT.TH1D.       wilsonEff = lambda haccepted,htotal : binomEff_h1 ( haccepted , htotal, func =       wilsonEff )
ROOT.TH1F. agrestiCoullEff = lambda haccepted,htotal : binomEff_h1 ( haccepted , htotal, func = agrestiCoullEff )
ROOT.TH1D. agrestiCoullEff = lambda haccepted,htotal : binomEff_h1 ( haccepted , htotal, func = agrestiCoullEff ) 

# =============================================================================
## @var one_sigma
#  the width of the +-1 sigma confidence interval 
one_sigma   = 0.682689492137086  ## the width of the +-1 sigma confidence interval 
## @var two_sigma
#  the width of the +-2 sigma confidence interval 
two_sigma   = 0.9544997361036415 ## the width of the +-2 sigma confidence interval 
## @var three_sigma
#  the width of the +-3 sigma confidence interval 
three_sigma = 0.9973002039367398 ## the width of the +-3 sigma confidence interval 
# =============================================================================
## calculate the efficiency graph using the binomial intervals 
#  @code 
#  >>> accepted   = ...
#  >>> rejected   = ...
#  >>> efficiency = acepted.eff_wald                    ( rejected )
#  >>> efficiency = acepted.eff_wilson_score            ( rejected )
#  >>> efficiency = acepted.eff_wilson_score_continuity ( rejected )
#  >>> efficiency = acepted.eff_arcsin                  ( rejected )
#  >>> efficiency = acepted.eff_agresti_coull           ( rejected )
#  >>> efficiency = acepted.eff_jeffreys                ( rejected )
#  >>> efficiency = acepted.eff_clopper_pearson         ( rejected )#
#  @endcode
#  @see Ostap::Math::wald_interval
#  @see Ostap::Math::wilson_score_interval
#  @see Ostap::Math::wilson_score_continuity_interval
#  @see Ostap::Math::arcsin_interval
#  @see Ostap::Math::agresti_coull_interval
#  @see Ostap::Math::jeffreys_interval
#  @see Ostap::Math::clopper_pearson_interval
#  @see https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def binom_interval_h1 ( accepted , rejected , func , interval = one_sigma ) :
    """ Calculate the efficiency graph using the binomial errors    
    >>> accepted   = ...
    >>> rejected   = ...
    >>> efficiency = acepted.eff_wald                    ( rejected )
    >>> efficiency = acepted.eff_wilson_score            ( rejected )
    >>> efficiency = acepted.eff_wilson_score_continuity ( rejected )
    >>> efficiency = acepted.eff_arcsin                  ( rejected )
    >>> efficiency = acepted.eff_agresti_coull           ( rejected )
    >>> efficiency = acepted.eff_jeffreys                ( rejected )
    >>> efficiency = acepted.eff_clopper_pearson         ( rejected )
    - see https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    - see Ostap::Math::wald_interval
    - see Ostap::Math::wilson_score_interval
    - see Ostap::Math::wilson_score_continuity_interval
    - see Ostap::Math::arcsin_interval
    - see Ostap::Math::agresti_coull_interval
    - see Ostap::Math::jeffreys_interval
    - see Ostap::Math::clopper_pearson_interval
    """
    #
    h1 = accepted
    h2 = rejected 
    if isinstance ( h1 , ROOT.TProfile ) :
        hh = h1.asH1()
        return binom_interval_h1 ( hh , h2 , func , interval )
    
    if not 0 <= interval <= 1 :
        logger.warning ( 'Invalid confidence interval: 0<=%s<==1' % interval )
        
    if                                 not h1.GetSumw2() : h1.Sumw2()
    if hasattr ( h2 , 'GetSumw2' ) and not h2.GetSumw2() : h2.Sumw2()
    #
    points = [] 
    #
    for i1,x1,y1 in h1.items() :
        #
        if not natural_entry ( y1 ) :
            raise AssertionError ( 'Non-natural entry[%d]=%s' % ( i1 , y1 ) ) 
        #
        y2 = h2 ( x1.value() ) 
        if not natural_entry ( y2 ) :
            raise AssertionError ( 'Non-natural entry(%s)=%s' % ( x1.value() , y2 ) )
        #
        l1 = long_type ( y1.value () )
        l2 = long_type ( y2.value () )
        #
        if  0 > l1 or 0 > l2 :
            raise AssertionError ( 'Incorrect values(%s): %d,%d' % ( x1.value() , l1 , l2 ) )
        #
        
        if   0 == l1 and 0 == l2 : center = 0.5
        elif 0 == l1             : center = 0.0
        elif 0 == l2             : center = 1.0
        else                     : center = 1./(1.+float(l2)/l1)
        
        res    = func ( l1 , l2 , interval )

        minv   = res.first
        maxv   = res.second
        #
        xv = x1.value()
        xe = x1.error()

        points.append ( ( xv , xe , xe , center , abs(center - minv) , abs(maxv - center ) ) ) 

    import ostap.histos.graphs 
    graph = ROOT.TGraphAsymmErrors ( len ( points ) )
    for p in range( len ( points ) ) : graph[p] = points[p]
    return graph 

ROOT.TH1F.eff_wald                    = lambda accepted,rejected,interval=one_sigma : binom_interval_h1 ( accepted , rejected , Ostap.Math.wald_interval                    , interval ) 
ROOT.TH1D.eff_wald                    = lambda accepted,rejected,interval=one_sigma : binom_interval_h1 ( accepted , rejected , Ostap.Math.wald_interval                    , interval ) 
ROOT.TH1F.eff_wilson_score            = lambda accepted,rejected,interval=one_sigma : binom_interval_h1 ( accepted , rejected , Ostap.Math.wilson_score_interval            , interval ) 
ROOT.TH1D.eff_wilson_score            = lambda accepted,rejected,interval=one_sigma : binom_interval_h1 ( accepted , rejected , Ostap.Math.wilson_score_interval            , interval ) 
ROOT.TH1F.eff_wilson_score_continuity = lambda accepted,rejected,interval=one_sigma : binom_interval_h1 ( accepted , rejected , Ostap.Math.wilson_score_continuity_interval , interval ) 
ROOT.TH1D.eff_wilson_score_continuity = lambda accepted,rejected,interval=one_sigma : binom_interval_h1 ( accepted , rejected , Ostap.Math.wilson_score_continuity_interval , interval ) 
ROOT.TH1F.eff_arcsin                  = lambda accepted,rejected,interval=one_sigma : binom_interval_h1 ( accepted , rejected , Ostap.Math.arcsin_interval                  , interval ) 
ROOT.TH1D.eff_arcsin                  = lambda accepted,rejected,interval=one_sigma : binom_interval_h1 ( accepted , rejected , Ostap.Math.arcsin_interval                  , interval ) 
ROOT.TH1F.eff_agresti_coull           = lambda accepted,rejected,interval=one_sigma : binom_interval_h1 ( accepted , rejected , Ostap.Math.agresti_coull_interval           , interval ) 
ROOT.TH1D.eff_agresti_coull           = lambda accepted,rejected,interval=one_sigma : binom_interval_h1 ( accepted , rejected , Ostap.Math.agresti_coull_interval           , interval ) 
ROOT.TH1F.eff_jeffreys                = lambda accepted,rejected,interval=one_sigma : binom_interval_h1 ( accepted , rejected , Ostap.Math.jeffreys_interval                , interval ) 
ROOT.TH1D.eff_jeffreys                = lambda accepted,rejected,interval=one_sigma : binom_interval_h1 ( accepted , rejected , Ostap.Math.jeffreys_interval                , interval ) 
ROOT.TH1F.eff_clopper_pearson         = lambda accepted,rejected,interval=one_sigma : binom_interval_h1 ( accepted , rejected , Ostap.Math.clopper_pearson_interval         , interval ) 
ROOT.TH1D.eff_clopper_pearson         = lambda accepted,rejected,interval=one_sigma : binom_interval_h1 ( accepted , rejected , Ostap.Math.clopper_pearson_interval         , interval ) 

# =============================================================================
## calculate the efficiency histogram using the binomial errors 
#  @code 
#  >>> accepted   = ...
#  >>> total      = ...
#  >>> efficiency = acepted // total
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def binomEff_h2 ( h1 , h2 , func = binomEff ) :
    """ Calculate the efficiency histogram using the binomial errors
    >>> accepted   = ...
    >>> total      = ...
    >>> efficiency = accepted // total    
    """
    #
    if                                 not h1.GetSumw2() : h1.Sumw2()
    if hasattr ( h2 , 'GetSumw2' ) and not h2.GetSumw2() : h2.Sumw2()
    #
    result = h1.Clone( hID () )
    if not result.GetSumw2() : result.Sumw2()
    #
    for ix1,iy1,x1,y1,z1 in h1.items() :
        #
        if not natural_entry ( z1 ) :
            raise AssertionError ( 'Non-natural entry[%d,%d]=%s' % ( ix1 , iy1 , z1 ) ) 
        #
        z2 = h2 ( x1.value() , y1.value() ) 
        if not natural_entry ( z2 ) :
            raise AssertionError ( 'Non-natural entry(%s,%s)=%s' % ( x1.value() , y1.value() , z2 ) )
        #
        l1 = long_type ( z1.value () )
        l2 = long_type ( z2.value () )
        #
        if  0 > l1 or 0 > l2 or l1 > l2 :
            raise AssertionError ( 'Incorrect values(%s,%s): %d,%d' % ( x1.value() , y1.value() , l1 , l2 ) )
        #
        v = VE ( func ( l1 , l2 ) ) 
        #
        if not v.isfinite() : continue 
        #
        result.SetBinContent ( ix1 , iy1 , v.value () ) 
        result.SetBinError   ( ix1 , iy1 , v.error () ) 
        
    return result 

ROOT.TH2F.  binomEff    = binomEff_h2 
ROOT.TH2D.  binomEff    = binomEff_h2 


# =============================================================================
## calculate the efficiency histogram using the binomial errors 
#  @code 
#  >>> accepted   = ...
#  >>> total      = ...
#  >>> efficiency = accepted // total
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def binomEff_h3 ( h1 , h2 , func = binomEff ) :
    """ Calculate the efficiency histogram using the binomial errors
    >>> accepted   = ...
    >>> total      = ...
    >>> efficiency = accepted // total    
    """
    #
    if                                 not h1.GetSumw2() : h1.Sumw2()
    if hasattr ( h2 , 'GetSumw2' ) and not h2.GetSumw2() : h2.Sumw2()
    #
    result = h1.Clone( hID () )
    if not result.GetSumw2() : result.Sumw2()
    #
    for ix1,iy1,iz1,x1,y1,z1,v1 in h1.items() :
        #
        if not natural_entry ( v1 ) :
            raise AssertionError ( 'Non-natural entry[%d,%d,%d]=%s' % ( ix1 , iy1 , iz1, v1 ) ) 
        #
        v2 = h2 ( x1.value() , y1.value() , z1.value() ) 
        if not natural_entry ( v2 ) :
            raise AssertionError ( 'Non-natural entry(%s,%s,%s)=%s' % ( x1.value() , y1.value() , z1.value() , v2 ) )
        
        l1 = long_type ( v1.value () )
        l2 = long_type ( v2.value () )
        #
        if  0 > l1 or 0 > l2 or l1 > l2 :
            raise AssertionError ( 'Incorrect values(%s,%s,%s): %d,%d' % ( x1.value() , y1.value() , z1.value() , l1 , l2 ) )
        #
        v = VE ( func ( l1 , l2 ) ) 
        #
        if not v.isfinite() : continue 
        #
        result.SetBinContent ( ix1 , iy1 , iz1 , v.value () ) 
        result.SetBinError   ( ix1 , iy1 , iz1 , v.error () ) 
        
    return result 

ROOT.TH3F .  binomEff     = binomEff_h3 
ROOT.TH3D .  binomEff     = binomEff_h3 

ROOT.TH1F . __floordiv__  = binomEff_h1 
ROOT.TH1D . __floordiv__  = binomEff_h1 
ROOT.TH2F . __floordiv__  = binomEff_h2
ROOT.TH2D . __floordiv__  = binomEff_h2
ROOT.TH3F . __floordiv__  = binomEff_h3
ROOT.TH3D . __floordiv__  = binomEff_h3
 
# =============================================================================
## calculate the efficiency histogram using the binomial errors
#  @code 
#  >>> accepted   = ...
#  >>> total      = ...
#  >>> efficiency = acepted // total
#  @endcode 
#  @see Ostap::Math::zechEff
#  @see Ostap::Math::wilsonEff
#  @see Ostap::Math::argestiCoullEff
#  @param h1 histogram of "accepted" sample
#  @param h2 histogram of "total"    sample 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def zechEff_h1 ( h1 , h2 , func = zechEff ) :
    """ Calculate the efficiency histogram using the binomial errors
    >>> accepted  = ... ##  histogram for accepted sample 
    >>> total     = ... ##  histogram for total    sample 
    >>> efficiency = accepted % total    
    """
    #
    if isinstance ( h1 , ROOT.TProfile ) :
        hh = h1.asH1()
        return zechEff_h1 ( hh , h2 , func )
    
    if                                 not h1.GetSumw2() : h1.Sumw2()
    if hasattr ( h2 , 'GetSumw2' ) and not h2.GetSumw2() : h2.Sumw2()
    #
    result = h1.Clone( hID () )
    if not result.GetSumw2() : result.Sumw2()
    #
    for i1,x1,y1 in h1.items() :
        #
        y2 = h2 ( x1.value() ) 
        #
        v  = VE ( func ( y1 , y2 ) )
        #
        if v.cov2() < 0 :
            ## if uncertainty is unphysical, make simple divizion
            v = y1/y2 
            #
        if not v.isfinite() : continue 
        #
        result.SetBinContent ( i1 , v.value () ) 
        result.SetBinError   ( i1 , v.error () )
        
    return result 

ROOT.TH1F.  zechEff    = zechEff_h1 
ROOT.TH1D.  zechEff    = zechEff_h1 

# =============================================================================
## calculate the efficiency histogram using the binomial errors
#  @code 
#  >>> accepted   = ...
#  >>> total      = ...
#  >>> efficiency = acepted // total
#  @endcode 
#  @see Ostap::Math::zechEff
#  @param h1 histogram of "accepted" sample
#  @param h2 histogram of "total"    sample 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def zechEff_h2 ( h1 , h2 ) :
    """ Calculate the efficiency histogram using the binomial errors    
    >>> accepted  = ... ##  histogram for accepted sample 
    >>> total     = ... ##  histogram for total    sample 
    >>> efficiency = accepted % total    
    """
    func = zechEff 
    #
    if                                 not h1.GetSumw2() : h1.Sumw2()
    if hasattr ( h2 , 'GetSumw2' ) and not h2.GetSumw2() : h2.Sumw2()
    #
    result = h1.Clone( hID () )
    if not result.GetSumw2() : result.Sumw2()
    #
    for ix1,iy1,x1,y1,z1 in h1.items() :
        #
        z2 = h2 ( x1.value() , y1.value() ) 
        #
        v  = VE ( func ( z1 , z2 ) )
        #
        if v.cov2() < 0 :
            ## if uncertainty is unphysical, make simple divizion
            v = z1/z2 
            #
        if not v.isfinite() : continue 
        #
        result.SetBinContent ( ix1 , iy1 , v.value () ) 
        result.SetBinError   ( ix1 , iy1 , v.error () ) 
        
    return result 

ROOT.TH2F.  zechEff    = zechEff_h2
ROOT.TH2D.  zechEff    = zechEff_h2 


# =============================================================================
## calculate the efficiency histogram using the binomial errors
#  @code 
#  >>> accepted   = ...
#  >>> total      = ...
#  >>> efficiency = acepted // total
#  @endcode 
#  @see Ostap::Math::zechEff
#  @param h1 histogram of "accepted" sample
#  @param h2 histogram of "total"    sample 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def zechEff_h3 ( h1 , h2 ) :
    """ Calculate the efficiency histogram using the binomial errors    
    >>> accepted  = ... ##  histogram for accepted sample 
    >>> total     = ... ##  histogram for total    sample 
    >>> efficiency = accepted % total    
    """
    func = zechEff 
    #
    if                                 not h1.GetSumw2() : h1.Sumw2()
    if hasattr ( h2 , 'GetSumw2' ) and not h2.GetSumw2() : h2.Sumw2()
    #
    result = h1.Clone( hID () )
    if not result.GetSumw2() : result.Sumw2()
    #
    for ix1,iy1,iz1,x1,y1,z1,v1 in h1.items() :
        #
        v2 = h2 ( x1.value() , y1.value() , z1.value() ) 
        #
        v  = VE ( func ( v1 , v2 ) )
        #
        if v.cov2() < 0 :
            ## if uncertainty is unphysical, make simple divizion
            v = v1/v2 
            #
        if not v.isfinite() : continue 
        #
        result.SetBinContent ( ix1 , iy1 , iz1 , v.value () ) 
        result.SetBinError   ( ix1 , iy1 , iz1 , v.error () ) 

        
    return result 

ROOT.TH3F .  zechEff = zechEff_h3
ROOT.TH3D .  zechEff = zechEff_h3

ROOT.TH1F . __mod__  = zechEff_h1 
ROOT.TH1D . __mod__  = zechEff_h1 
ROOT.TH2F . __mod__  = zechEff_h2
ROOT.TH2D . __mod__  = zechEff_h2
ROOT.TH3F . __mod__  = zechEff_h3
ROOT.TH3D . __mod__  = zechEff_h3

# =============================================================================
## get binomial division by ROOT
#  @code 
#  accepted  = ... ##  histogram for accepted sample 
#  total     = ... ##  histogram for total    sample 
#  efficiency = accepted.binomEff2 ( total )
#  @endcode 
#  - <code>ROOT.TH1.Divide ( ... , 1 , 1 , "B")</code> is used 
def _h_binomEff_2_ ( accepted , total ) :
    """ Calculate the efficiency histogram using the binomial errors    
    >>> accepted  = ... ##  histogram for accepted sample 
    >>> total     = ... ##  histogram for total    sample 
    >>> efficiency = accepted.binomEff2 ( total )
    - ROOT.TH1.Divide ( ... , 1 , 1 , 'B') is used 
    """
    if not accepted.same_bins ( total ) :
        total1 =  accepted.clone()
        total1.Reset()
        total1 += total 
        result = _h_binomEff_2_ ( accepted  , total1 )
        del total1
        return result 

    result = accepted.clone()
    result.Reset()
    ok = result.Divide ( accepted , total , 1.0 , 1.0 , 'b' )
    if not ok :
        logger.warning("Can't ``B''-divide histograms, switch to Zech's method")
        return accepted % total 
    return result 
    

ROOT.TH1.binomEff_2 = _h_binomEff_2_


# =============================================================================
## consider object as function:
def objectAsFunction ( obj ) :

    if   isinstance ( obj , num_types ) :
        
        val  = float  ( obj ) 
        func = lambda x,*y : VE ( val , 0 )
        return func                                      ## RETURN
    
    elif isinstance ( obj  ,    VE ) :
        
        val  =     VE ( obj )
        func = lambda x,*y : val
        return func                                      ## RETURN
    
    elif isinstance ( obj ,   ROOT.TF3 ) :
        
        f3   = obj
        func = lambda x,y,z,*t : VE ( f3 ( float ( x ) , float ( y ) , float ( z ) , 0 ) )
        return func                                      ## RETURN 

    elif isinstance ( obj ,   ROOT.TF2 ) :
        
        f2   = obj
        func = lambda x , y , *z   : VE ( f2 ( float ( x ) , float ( y ) , 0 ) )
        return func                                      ## RETURN 

    elif isinstance ( obj ,   ROOT.TF1 ) :
        
        f1   = obj
        func = lambda x,*y     : VE ( f1 ( float ( x ) , 0 ) )
        return func                                      ## RETURN 

    ## the original stuff 
    return obj


# ==============================================================================
## trivial helper base ....
class FUNCX  (object) : pass

# ==============================================================================
## constant 
class FUNC_C0(FUNCX) :
    def __init__ ( self , value       ) : self.value = value 
    def __call__ ( self , x ,      *y ) : return self.value 
    def integral ( self , xmin , xmax ) : return ( xmax - xmin ) * self.value
    
# ==============================================================================
## trivial 1D-function 
class FUNC_F1(FUNCX) :
    def __init__ ( self , func        ) : self.func = func 
    def __call__ ( self , x    ,   *y ) : return self.func ( float ( x ) )  
    def integral ( self , xmin , xmax ) : return self.func.Integral( xmin , xmax )
    
# ==============================================================================
## trivial 2D-function 
class FUNC_F2(FUNCX) :
    def __init__ ( self , func        ) : self.func = func 
    def __call__ ( self , x , y ,  *z ) : return self.func ( float ( x ) , float ( y ) ) 
    def integral ( self ,
                   xmin , xmax        ,
                   ymin , ymax        ) : return self.func.Integral( xmin , xmax ,
                                                                     ymin , ymax )
# ==============================================================================
## trivial 3D-function 
class FUNC_F3(FUNCX) :
    def __init__ ( self , func        ) : self.func = func 
    def __call__ ( self , x , y ,  z  ) : return self.func ( float ( x ) ,
                                                             float ( y ) ,
                                                             float ( z ) ) 
    def integral ( self ,
                   xmin , xmax        ,
                   ymin , ymax        , 
                   zmin , zmax        ) : return self.func.Integral( xmin , xmax ,
                                                                     xmin , xmax ,
                                                                     zmin , zmax )
# ==============================================================================
## everything else 
class FUNC_OTHER(FUNCX) :
    def __init__ ( self , func        ) :
        self.func = func
        if hasattr ( self.func , 'integral' ) :
            self.integral = lambda x,y,*z : self.func.integral ( x , y , *z ) 
        else :
            from ostap.math.integral import integral
            self.integral = lambda x,y,*z : integral ( self.func , x , y , *z ) 
    def __call__ ( self , x , *y      ) : return self.func ( float ( x ) , *y ) 
    def integral ( self , s ,  y , *z ) : return self.integral ( x , y , *z )



# =============================================================================
## consider object as function:
def objectAsFunctionObject ( obj ) :
    
    if   isinstance ( obj , FUNCX     ) : return obj 
    elif isinstance ( obj , num_types ) : return FUNC_CO ( obj ) 
    elif isinstance ( obj , VE        ) : return FUNC_CO ( obj ) 
    elif isinstance ( obj , ROOT.TF1  ) : return FUNC_F1 ( obj ) 
    elif isinstance ( obj , ROOT.TF2  ) : return FUNC_F2 ( obj ) 
    elif isinstance ( obj , ROOT.TF3  ) : return FUNC_F3 ( obj ) 
    #
    ## the original stuff    
    return FUNC_OTHER ( obj ) 

# =============================================================================
## operation with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_oper_ ( h1 , h2 , oper ) :
    """ Operation with the histogram
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1 (oper) h2 
    """
    if isinstance ( h1 , ROOT.TProfile ) :
        hh = h1.asH1()
        return _h1_oper_ ( hh , h2 , oper ) 
    #
    if                                 not h1.GetSumw2() : h1.Sumw2()
    if hasattr ( h2 , 'GetSumw2' ) and not h2.GetSumw2() : h2.Sumw2()
    #
    result = h1.Clone( hID() )
    if not result.GetSumw2() : result.Sumw2()

    ## 
    f2 = objectAsFunction ( h2 )
    
    ##
    for i1,x1,y1 in h1.items() :
        #
        result.SetBinContent ( i1 , 0 ) 
        result.SetBinError   ( i1 , 0 )
        #
        y2 = f2 ( x1.value() ) 
        #
        v  = VE ( oper ( y1 , y2 ) ) 
        #
        if not v.isfinite() : continue 
        #
        result.SetBinContent ( i1 , v.value () ) 
        result.SetBinError   ( i1 , v.error () )
        
    result.ResetStats() 
    return result

# =============================================================================
## operation with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_ioper_ ( h1 , h2 , oper ) :
    """ Operation with the histogram
    >>> obj= ...
    >>> h2     = ...
    >>> h2 (oper) obj    
    """
    ##
    if                                 not h1.GetSumw2() : h1.Sumw2()
    if hasattr ( h2 , 'GetSumw2' ) and not h2.GetSumw2() : h2.Sumw2()
    #
    f2 = objectAsFunction ( h2 ) 
    ##
    for i1,x1,y1 in h1.items() :
        #
        y2 = f2 ( x1.value() ) 
        #
        v  = VE ( oper ( y1 , y2 ) ) 
        #
        if not v.isfinite() : continue 
        #
        h1.SetBinContent ( i1 , v.value () ) 
        h1.SetBinError   ( i1 , v.error () )
        
    h1.ResetStats() 
    return h1 



# =============================================================================
##  Division with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_div_ ( h1 , h2 ) :
    """ Divide the histograms 
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1 / h2  
    """
    #
    return _h1_oper_ ( h1 , h2 , lambda x,y : x/y ) 
# =============================================================================
##  Division with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_mul_ ( h1 , h2 ) :
    """ Multiply the histograms 
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1 * h2  
    """
    return _h1_oper_ ( h1 , h2 , lambda x,y : x*y ) 
# =============================================================================
##  Addition with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_add_ ( h1 , h2 ) :
    """ Add the histograms     
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1 + h2  
    """
    return _h1_oper_ ( h1 , h2 , lambda x,y : x+y ) 
# =============================================================================
## Subtraction of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_sub_ ( h1 , h2 ) :
    """ Subtract the histogram     
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1 - h2  
    """
    return _h1_oper_ ( h1 , h2 , lambda x,y : x-y ) 
# =============================================================================
##  Fraction of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_frac_ ( h1 , h2 ) :
    """ `Fraction' the histogram h1/(h1+h2)    
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1.frac  ( h2 ) 
    """
    return _h1_oper_ ( h1 , h2 , lambda x,y : x.frac(y) ) 
# =============================================================================
##  ``Asymmetry'' of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_asym_ ( h1 , h2 ) :
    """ `Asymmetry' the histogram (h1-h2)/(h1+h2)    
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1.asym ( h2 )     
    """
    return _h1_oper_ ( h1 , h2 , lambda x,y : x.asym(y) ) 


# =============================================================================
## `minimum' of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2023-09-11
def _h1_min_ ( h1 , h2 ) :
    """ `minimum' the histogram min(j1,h2)
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1.min ( h2 )     
    """
    return _h1_oper_ ( h1 , h2 , lambda x,y : min ( x , y ) ) 

# =============================================================================
## `maximum' of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2023-09-11
def _h1_max_ ( h1 , h2 ) :
    """ `maximum' the histogram min(j1,h2)
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1.max ( h2 )     
    """
    return _h1_oper_ ( h1 , h2 , lambda x,y : max ( x , y ) ) 


# =============================================================================
## `Difference' of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_diff_ ( h1 , h2 ) :
    """``Difference'' the histogram 2*(h1-h2)/(h1+h2)    
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1.diff ( h2 )     
    """
    return _h1_oper_ ( h1 , h2 , lambda x,y : 2*x.asym(y) )

# =============================================================================
##  ``Chi2-tension'' of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_chi2_ ( h1 , h2 ) :
    """``Chi2-tension'' the histogram    
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1.chi2  ( h2 )     
    """
    return _h1_oper_ ( h1 , h2 , lambda x,y : VE ( x.chi2 ( y ) , 0 ) ) 
# =============================================================================
##  ``Average'' of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_average_ ( h1 , h2 ) :
    """``Mean'' the histograms    
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1.average  ( h2 )     
    """
    return _h1_oper_ ( h1 , h2 , lambda x,y : x.mean ( y ) ) 

# =============================================================================
## 'pow' the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_pow_ ( h1 , val ) :
    """``pow'' the histogram
    >>> histo = ...
    >>> other = histo**5 
    """
    if isinstance ( h1 , ROOT.TProfile ) :
        hh = h1.asH1()
        return _h1_pow_ ( hh , val )
    #
    if not h1.GetSumw2() : h1.Sumw2()
    #
    result = h1.Clone( hID() )
    if not result.GetSumw2() : result.Sumw2()
    #
    for i1,x1,y1 in h1.items() :
        #
        result.SetBinContent ( i1 , 0 ) 
        result.SetBinError   ( i1 , 0 )
        #
        v = VE ( pow ( y1 , val ) ) 
        #
        if not v.isfinite() : continue 
        #
        result.SetBinContent ( i1 , v.value () ) 
        result.SetBinError   ( i1 , v.error () )

    result.ResetStats() 
    return result 

# =============================================================================
## 'abs' the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-04-29
def _h1_abs_ ( h1 ) :
    """``abs'' the histogram
    >>> h      = ...
    >>> result = abs ( h )
    """
    if not h1.GetSumw2() : h1.Sumw2()
    #
    result = h1.Clone( hID() )
    if not result.GetSumw2() : result.Sumw2()
    #
    for i1,x1,y1 in h1.items() :
        #
        result.SetBinContent ( i1 , 0 ) 
        result.SetBinError   ( i1 , 0 )
        #
        v = abs ( y1 )  
        #
        result.SetBinContent ( i1 , v.value () ) 
        result.SetBinError   ( i1 , v.error () )
        
    result.ResetStats() 
    return result 

# =============================================================================
## Division with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_idiv_ ( h1 , h2 ) :
    """ Divide the histograms     
    >>> h1  = ...
    >>> h2  = ...
    >>> h1 /=  h2     
    """
    return _h1_ioper_ ( h1 , h2 , lambda x,y : x/y ) 

# =============================================================================
## Multiplication with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_imul_ ( h1 , h2 ) :
    """ Multiply the histograms     
    >>> h1  = ...
    >>> h2  = ...
    >>> h1 *=  h2     
    """
    return _h1_ioper_ ( h1 , h2 , lambda x,y : x*y ) 

# =============================================================================
## Addition with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_iadd_ ( h1 , h2 ) :
    """ Add the histograms    
    >>> h1  = ...
    >>> h2  = ...
    >>> h1 +=  h2     
    """
    return _h1_ioper_ ( h1 , h2 , lambda x,y : x+y ) 
# =============================================================================
##  Subtraction of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_isub_ ( h1 , h2 ) :
    """ Subtract the histogram    
    >>> h1  = ...
    >>> h2  = ...
    >>> h1 -=  h2     
    """
    return _h1_ioper_ ( h1 , h2 , lambda x,y : x-y ) 

# =============================================================================
## Division with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_rdiv_ ( h1 , h2 ) :
    """ Divide the histograms
    >>> h1     = ...
    >>> obj    = ...
    >>> result = obj / h1 
    """
    return _h1_oper_ ( h1 , h2 , lambda x,y : y/x ) 
# =============================================================================
## Multiplication with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_rmul_ ( h1 , h2 ) :
    """ Multiply the histograms 
    >>> h1     = ...
    >>> obj    = ...
    >>> result = obj * h1 
    """
    return _h1_oper_ ( h1 , h2 , lambda x,y : y*x ) 

# =============================================================================
## Addition with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_radd_ ( h1 , h2 ) :
    """ Add the histograms    
    >>> h1     = ...
    >>> obj    = ...
    >>> result = obj + h1 
    """
    return _h1_oper_ ( h1 , h2 , lambda x,y : y+x )

# =============================================================================
## Subtraction of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_rsub_ ( h1 , h2 ) :
    """ Subtract the histogram            
    >>> h1     = ...
    >>> obj    = ...
    >>> result = obj - h1 
    """
    return _h1_oper_ ( h1 , h2 , lambda x,y : y-x ) 

# =============================================================================
## Feed the histogram from other object, e.g. function
#  @code
#  h1 = ....         ## create the histogram 
#  math.sin  >> h1   ## feed it using function  
#  @endcode
#  The action is equivalent to
#  @code
#  h1.Reset()        ## reset it 
#  h1 += math.sin    ## update it using function  
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_rrshift_ ( h1 , obj ) :
    """ Feed the histogram from external source, e.g. function    
    >>> h1 = ...        ## the histogram
    >>> math.sin >> h1  ## feed it from the function
    """
    h1.Reset() 
    h1 += obj
    
    return h1

# =============================================================================
## decorate:
for t in ( ROOT.TH1F , ROOT.TH1D ) : 
    
    t . _oper_       = _h1_oper_
    t . _ioper_      = _h1_ioper_
    t . __mul__      = _h1_mul_
    t . __add__      = _h1_add_
    t . __sub__      = _h1_sub_
    t . __pow__      = _h1_pow_
    t . __truediv__  = _h1_div_
    t . __div__      = _h1_div_

    t . __imul__     = _h1_imul_
    t . __iadd__     = _h1_iadd_
    t . __isub__     = _h1_isub_
    t . __itruediv__ = _h1_idiv_
    t . __idiv__     = _h1_idiv_
    
    t . __rtruediv__ = _h1_rdiv_
    t . __rdiv__     = _h1_rdiv_
    t . __rmul__     = _h1_rmul_
    t . __radd__     = _h1_radd_
    t . __rsub__     = _h1_rsub_

    t . __rrshift__  = _h1_rrshift_
    
    t . __abs__      = _h1_abs_
    t .  frac        = _h1_frac_
    t .  asym        = _h1_asym_
    t .  diff        = _h1_diff_
    t .  chi2        = _h1_chi2_
    t .  average     = _h1_average_
    t .  min         = _h1_min_
    t .  max         = _h1_max_

# =============================================================================
## find the first X-value for the given Y-value 
def _h1_xfind_ ( self           ,
                 v              ,
                 forward = True ) :
    """ Find the first X-value for the given Y-value
    >>> h1 = ...
    >>> x = h1.find_X ( 1000 )
    >>> h1 = ...
    >>> x = h1.find_X ( 1000 , False )
    """
    #
    v = float ( v ) 
    #
    mn = self [ self.GetMinimumBin() ]
    mx = self [ self.GetMaximumBin() ]

    #
    if v < mn.value() :
        return self.xmin() if forward else self.xmax ()
    #
    if v > mx.value() :
        return self.xmax() if forward else self.xmin () 
    
    nb = len ( self )
    ax = self.GetXaxis () 
    
    for i in self :
        
        j  = i+1
        if not j in self : continue
        
        ib = i if forward else nb + 1 - i 
        jb = j if forward else nb + 1 - j 
        
        vi = self [ ib ].value () 
        vj = self [ jb ].value ()
        
        if  vi <= v <= vj or vj <= v <= vi :
            
            xi = ax.GetBinCenter ( ib )
            xj = ax.GetBinCenter ( jb )
            
            if   vi == v             : return xi   ## RETURN 
            elif vj == v             : return xj   ## RETURN 
            
            dv = vi - vj
            dx = xi - xj
            
            if   vi == vj or 0 == dv : return 0.5 * ( xi + xj )  ## RETURN
            
            return (v*dx+vi*xj-vj*xi)/dv                         ## RETURN
        

ROOT.TH1F . find_X = _h1_xfind_
ROOT.TH1D . find_X = _h1_xfind_


# =======================================================================
## get "n-sigma" interval assuming the presence of peak
def _h1_CL_interval_ ( self       ,
                       nsigma = 1 ) : 
    """ Get `n-sigma' interval for the distribution.
    >>> h = ...
    >>> xlow,xhigh = h.cl_interval( 1 ) 
    """
    mv   = self.maxv()
    import math 
    
    try :
        #
        f    = math.exp    ( -0.5 * nsigma * nsigma )
        x1   = self.find_X ( mv * f , True  )
        x2   = self.find_X ( mv * f , False )
        #
    except :
        #
        x1,x2 = self.xminmax () 
        
    return x1,x2 
    
ROOT.TH1F . cl_interval = _h1_CL_interval_
ROOT.TH1D . cl_interval = _h1_CL_interval_


# =============================================================================
## get the minumum value for the histogram 
def _h_minv_ ( self ) :
    """ Get the minimum value for the histogram
    >>> h  = ...
    >>> mv = h.minv ()
    """
    mv = VE ( 1.e+100 , -1 ) 
    for ibin in self :
        v = self[ibin]
        if v.value() <= mv.value() or mv.cov2() < 0 : mv = v
    return mv 

# =============================================================================
## get the maximum value for the histogram 
def _h_maxv_ ( self ) :
    """ Get the maximum value for the histogram
    >>> h  = ...
    >>> mv = h.maxv ()
    """
    mv = VE ( -1.e+100 , -1 ) 
    for ibin in self :
        v = self[ibin]
        if v.value() >= mv.value() or mv.cov2() < 0 : mv = v
    return mv 
        
# =============================================================================
## get the minmax values for the histogram 
def _h_minmax_ ( self ) :
    """ Get the minmax pair for the histogram    
    >>> h     = ...
    >>> mn,mx = h.minmax ()
    """
    return self.minv() , self.maxv() 

ROOT.TH1 . minv    = _h_minv_
ROOT.TH1 . maxv    = _h_maxv_
ROOT.TH1 . minmax  = _h_minmax_

# ============================================================================
## Get the minimal positive entry in the historgam
#  - It is useful e.g. to define the proper scale for drawing, especially for log-scale
#  @code
#  histo  = ...
#  minpos = histo.min_positive() 
#  @endcode
#  @return minimal positive entry or negative infinity if no positive entries
def _h_minpos_ ( histo ) :
    """ Get the minimal positive entry in the histogram
    - It is useful e.g. to define the proper scale for drawing, especially for log-scale
    - return minimal positive entry or negative infinity if no positive entries 
    >>> histo  = ...
    >>> minpos = histo.min_positive() 
    """
    result = neg_infinity 
    for items in histo.items() :
        v = items [ -1 ] . value()
        if 0 < v  :
            result = min ( v , result ) if 0 < result else v 
    return result

ROOT.TH1.minpos       = _h_minpos_
ROOT.TH1.min_positive = _h_minpos_


# ============================================================================
## Only positive entries ?
#  @code
#  histo = ...
#  only_positive = histo.all_positive() 
#  @endcode
def _h_all_positive ( histo ) :
    """ Only positive entries ?
    >>> histo = ...
    >>> only_positive = histo.all_positive() 
    """
    return all ( 0 < v.value() for v in histo.values() ) 

# ============================================================================
## Only negative entries ?
#  @code
#  histo = ...
#  only_negative  = histo.all_negative () 
#  @endcode
def _h_all_negative ( histo ) :
    """ Only negative entries ?
    >>> histo = ...
    >>> only_bnegative  = histo.all_negative () 
    """
    return all ( 0 > v.value() for v in histo.values() ) 

# ============================================================================
## Only non-negative entries ?
#  @code
#  histo = ...
#  only_nonnegative  = histo.all_nonnegative () 
#  @endcode
def _h_all_nonnegative ( histo ) :
    """ Only non-negative  entries ?
    >>> histo = ...
    >>> only_positive = histo.all_positive() 
    """
    return all ( 0 <= v.value() or iszero ( v ) for v in histo.values() ) 

# ============================================================================
## Only non-positive entries ?
#  @code
#  histo = ...
#  only_nonpositive  = histo.all_nonpositive  () 
#  @endcode
def _h_all_nonpositive ( histo ) :
    """ Only non-positive  entries ?
    >>> histo = ...
    >>> only_nonpositive = histo.all_nonpositive() 
    """
    return all ( 0 >= v.value() or iszero ( v ) for v in histo.values() ) 

# ============================================================================
## Only zero entries ?
#  @code
#  histo = ...
#  only_zero  = histo.all_zero  () 
#  @endcode
def _h_all_zero ( histo ) :
    """ Only zero entries ?
    >>> histo = ...
    >>> only_zero = histo.all_sero() 
    """
    return all ( iszero ( v ) for v in histo.values() ) 

for h in ( ROOT.TH1F , ROOT.TH1D ,
           ROOT.TH2F , ROOT.TH2D ,
           ROOT.TH3F , ROOT.TH3D ) : 
    h.all_zero         = _h_all_zero 
    h.all_positive     = _h_all_positive
    h.all_negative     = _h_all_negative
    h.all_nonpositive  = _h_all_nonpositive
    h.all_nonnegative  = _h_all_nonnegative
    h.only_zero        = _h_all_zero 
    h.only_positive    = _h_all_positive
    h.only_negative    = _h_all_negative
    h.only_nonpositive = _h_all_nonpositive
    h.only_nonnegative = _h_all_nonnegative
    
    _new_methods_ += [
        h.all_zero         , 
        h.all_positive     , 
        h.all_negative     , 
        h.all_nonpositive  , 
        h.all_nonnegative  , 
        h.only_zero        , 
        h.only_positive    ,
        h.only_negative    , 
        h.only_nonpositive , 
        h.only_nonnegative , 
    ]  

    
           
# ============================================================================
## get the first bin with the minimum content
def _h1_minimum_bin_ ( h1 ) :
    """ Get the first bin with the minimum content
    """
    value = None
    bin   = 0
    for i in h1 :
        v = h1.GetBinContent ( i )
        if value is None or v < value :
            value = v
            bin   = i
    return bin

# ============================================================================
## get the first bin with the maximal content
def _h1_maximum_bin_ ( h1 ) :
    """ Get the first bin with the maximal content
    """
    value = None
    bin   = 0
    for i in h1 :
        v = h1.GetBinContent ( i )
        if value is None or value < v :
            value = v
            bin   = i
    return bin

# =============================================================================
## get the first bin with the minimum content
def _h2_minimum_bin_ ( h2 ) :
    """ Get the first bin with the minimum content
    """
    value = None
    bin   = 0,0
    for i in h2 :
        v = h2.GetBinContent ( *i )
        if value is None or v < value :
            value = v
            bin   = i 
    return bin

# =============================================================================
## get the first bin with the maximal content
def _h2_maximum_bin_ ( h2 ) :
    """ Get the first bin with the maximal content
    """
    value = None
    bin   = 0,0
    for i in h2 :
        v = h2.GetBinContent ( *i )
        if value is None or value < v :
            value = v
            bin   = i 
    return bin

# =============================================================================
## get the first bin with the minimum content
def _h3_minimum_bin_ ( h3 ) :
    """ Get the first bin with the minimum content
    """
    value = None
    bin   = 0,0
    for i in h3 :
        v = h3.GetBinContent ( *i )
        if value is None or v < value :
            value = v
            bin   = i 
    return bin

# =============================================================================
## get the first bin with the maximal content
def _h3_maximum_bin_ ( h3 ) :
    """ Get the first bin with the maximal content
    """
    value = None
    bin   = 0,0
    for i in h3 :
        v = h3.GetBinContent ( *i )
        if value is None or value < v :
            value = v
            bin   = i 
    return bin

ROOT.TH1.minimum_bin = _h1_minimum_bin_
ROOT.TH1.maximum_bin = _h1_maximum_bin_
ROOT.TH2.minimum_bin = _h2_minimum_bin_
ROOT.TH2.maximum_bin = _h2_maximum_bin_
ROOT.TH3.minimum_bin = _h3_minimum_bin_
ROOT.TH3.maximum_bin = _h3_maximum_bin_

_new_methods_ += [ 
    ROOT.TH1.minimum_bin , 
    ROOT.TH1.maximum_bin , 
    ROOT.TH2.minimum_bin , 
    ROOT.TH2.maximum_bin , 
    ROOT.TH3.minimum_bin , 
    ROOT.TH3.maximum_bin , 
    ]

# ============================================================================
## get the minimum value for X-axis 
def _ax_min_ ( self ) :
    """ Get the minimum value for X-axis
    >>> xmin = h.xmin()
    """
    ax = self.GetXaxis () 
    return ax.GetXmin()
# ============================================================================
## get the minimum value for y-axis 
def _ay_min_ ( self ) :
    """ Get the minimum value for Y-axis
    >>> ymin = h.ymin()
    """
    ay = self.GetYaxis () 
    return ay.GetXmin()
# ============================================================================
## get the minimum value for z-axis 
def _az_min_ ( self ) :
    """ Get the minimum value for Z-axis
    >>> zmin = h.zmin()
    """
    az = self.GetZaxis () 
    return az.GetXmin()

# ============================================================================
## get the maximum value for X-axis 
def _ax_max_ ( self ) :
    """ Get the maximum value for X-axis
    >>> xmax = h.xmax()
    """
    ax = self.GetXaxis () 
    return ax.GetXmax()
# ============================================================================
## get the maximum value for y-axis 
def _ay_max_ ( self ) :
    """ Get the maximum value for Y-axis
    >>> ymax = h.ymax()
    """
    ay = self.GetYaxis () 
    return ay.GetXmax()
# ============================================================================
## get the maximum value for z-axis 
def _az_max_ ( self ) :
    """Get the maximum value for Z-axis
    >>> zmax = h.zmax()
    """
    az = self.GetZaxis () 
    return az.GetXmax()


ROOT.TH1D. xmin    = _ax_min_
ROOT.TH1D. xmax    = _ax_max_
ROOT.TH1D. ymin    = _h_minv_
ROOT.TH1D. ymax    = _h_maxv_

ROOT.TH1F. xmin    = _ax_min_
ROOT.TH1F. xmax    = _ax_max_
ROOT.TH1F. ymin    = _h_minv_
ROOT.TH1F. ymax    = _h_maxv_

ROOT.TH2 . xmin    = _ax_min_
ROOT.TH2 . xmax    = _ax_max_
ROOT.TH2 . ymin    = _ay_min_
ROOT.TH2 . ymax    = _ay_max_
ROOT.TH2 . zmin    = _h_minv_
ROOT.TH2 . zmax    = _h_maxv_

ROOT.TH3 . xmin    = _ax_min_
ROOT.TH3 . xmax    = _ax_max_
ROOT.TH3 . ymin    = _ay_min_
ROOT.TH3 . ymax    = _ay_max_
ROOT.TH3 . zmin    = _az_min_
ROOT.TH3 . zmax    = _az_max_

ROOT.TH1D. xminmax = lambda s : ( s.xmin() , s.xmax() )
ROOT.TH1D. yminmax = lambda s : ( s.ymin() , s.ymax() )
ROOT.TH1F. xminmax = lambda s : ( s.xmin() , s.xmax() )
ROOT.TH1F. yminmax = lambda s : ( s.ymin() , s.ymax() )

ROOT.TH2 . xminmax = lambda s : ( s.xmin() , s.xmax() )
ROOT.TH2 . yminmax = lambda s : ( s.ymin() , s.ymax() )
ROOT.TH2 . zminmax = lambda s : ( s.zmin() , s.zmax() )

ROOT.TH3 . xminmax = lambda s : ( s.xmin() , s.xmax() )
ROOT.TH3 . yminmax = lambda s : ( s.ymin() , s.ymax() )
ROOT.TH3 . zminmax = lambda s : ( s.zmin() , s.zmax() )



# =============================================================================
## get the random pair from 2D-histogram
#  @code 
#  >>> h2  = ...
#  >>> x,y = h2.random() 
#  @endcode
#  @see ROOT::TH2::GetRandom2
#  @see ROOT::TH2
def _h2_random_ ( h2 ) :
    """ Get the random value from 2D-histogram    
    >>> h2  = ...
    >>> x,y = h2.random() 
    """
    #
    x = ctypes.c_double ( 0.0 ) 
    y = ctypes.c_double ( 1.0 )
    #
    h2.GetRandom2( x , y )
    #
    return float( x.value ) , float( y.value )

# =============================================================================
## get the random triplet from 3D-histogram
#  @code 
#  >>> h3    = ...
#  >>> x,y,z = h3.random() 
#  @endcode
#  @see ROOT::TH3::GetRandom3
#  @see ROOT::TH3
def _h3_random_ ( h3 ) :
    """ Get the random value from 2D-histogram    
    >>> h3    = ...
    >>> x,y,z = h3.random() 
    """
    #
    x = ctypes.c_double ( 0.0 ) 
    y = ctypes.c_double ( 1.0 )
    z = ctypes.c_double ( 2.0 )
    #
    h3.GetRandom3 ( x , y , z )
    #
    return float( x.value ) , float( y.value ) , float( z.value )


ROOT.TH2F.random = _h2_random_
ROOT.TH2D.random = _h2_random_

ROOT.TH3F.random = _h3_random_
ROOT.TH3D.random = _h3_random_

ROOT.TH1F.random = lambda s : s.GetRandom() 
ROOT.TH1D.random = lambda s : s.GetRandom()

# =============================================================================
## Shoot several random numbers from 1D histogram
#  @code
#  h1 = ...
#  for v in h1.shoot ( 10 ) :
#  .... 
#  @endcode
def _h1_shoot_ ( h1 , N , accept = lambda v : True ) :
    """ Shoot several random numbers from 1D histogram
    >>> h1 = ...
    >>> for v in h1.shoot ( 10 ) :
    >>> .... 
    """
    
    assert isinstance ( N , integer_types ) and 0 <= N ,\
           'Invalid number of random shoots!'
    
    for i in range ( N ) : 
        
        v = h1.random()
        while not accept ( v ) :
            v = h1.random()

        yield v 
        
# =============================================================================
## Shoot several random numbers from 2D histogram
#  @code
#  h2 = ...
#  for x,y in h2.shoot ( 10 ) :
#  .... 
#  @endcode
def _h2_shoot_ ( h2 , N , accept = lambda x,y : True ) :
    """ Shoot several random numbers from 2D histogram
    >>> h2 = ...
    >>> for x,y in h2.shoot ( 10 ) :
    >>> .... 
    """
    
    assert isinstance ( N , integer_types ) and 0 <= N ,\
           'Invalid number of random shoots!'
    
    for i in range ( N ) : 
        
        x, y  = h2.random()
        while not accept ( x , y ) :
            x, y  = h2.random()
            
        yield x , y 

# =============================================================================
## Shoot several random numbers from 3D histogram
#  @code
#  h3 = ...
#  for x,y,z in h3.shoot ( 10 ) :
#  .... 
#  @endcode
def _h3_shoot_ ( h3 , N , accept = lambda x,y,z : True ) :
    """ Shoot several random numbers from 3D histogram
    >>> h3 = ...
    >>> for x,y,z in h3.shoot ( 10 ) :
    >>> .... 
    """
    
    assert isinstance ( N , integer_types ) and 0 <= N ,\
           'Invalid number of random shoots!'
    
    for i in range ( N ) : 
        
        x, y , z  = h3.random()
        while not accept ( x , y , z ) :
            x , y , z  = h3.random()
            
        yield x , y , z 

ROOT.TH1F.shoot = _h1_shoot_
ROOT.TH1D.shoot = _h1_shoot_
ROOT.TH2F.shoot = _h2_shoot_
ROOT.TH2D.shoot = _h3_shoot_
ROOT.TH3F.shoot = _h3_shoot_
ROOT.TH3D.shoot = _h3_shoot_

# =============================================================================
## operation with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_oper_ ( h1 , h2 , oper ) :
    """ Operation with the histogram        
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1 (oper) h2    
    """
    if                                 not h1.GetSumw2() : h1.Sumw2()
    if hasattr ( h2 , 'GetSumw2' ) and not h2.GetSumw2() : h2.Sumw2()
    #
    result = h1.Clone( hID() )
    if not result.GetSumw2() : result.Sumw2()
    #
    f2 = objectAsFunction ( h2 )
    # 
    for ix1,iy1,x1,y1,z1 in h1.items() :
        #
        result.SetBinContent ( ix1 , iy1 , 0 ) 
        result.SetBinError   ( ix1 , iy1 , 0 )
        #
        z2 = f2 ( x1.value() , y1.value() ) 
        #
        v  = VE ( oper ( z1 , z2 ) ) 
        #
        if not v.isfinite() : continue 
        #
        result.SetBinContent ( ix1 , iy1 , v.value () ) 
        result.SetBinError   ( ix1 , iy1 , v.error () )
        
    result.ResetStats() 
    return result

# =============================================================================
## operation with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-06-03
def _h2_ioper_ ( h1 , h2 , oper ) :
    """
    Operation with the histogram 
    """
    if                                 not h1.GetSumw2() : h1.Sumw2()
    if hasattr ( h2 , 'GetSumw2' ) and not h2.GetSumw2() : h2.Sumw2()
    #
    f2 = objectAsFunction ( h2 )
    # 
    for ix1 , iy1 , x1 , y1 , z1 in h1.items() :
        #
        h1.SetBinContent ( ix1 , iy1 , 0 ) 
        h1.SetBinError   ( ix1 , iy1 , 0 )
        #
        z2 = f2 ( x1.value() , y1.value() ) 
        #
        v  = VE ( oper ( z1 , z2 ) )

        #
        if not v.isfinite() : continue 
        #
        h1.SetBinContent ( ix1 , iy1 , v.value () ) 
        h1.SetBinError   ( ix1 , iy1 , v.error () )

    h1.ResetStats() 
    return h1

# =============================================================================
## Division with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_div_ ( h1 , h2 ) :
    """ Divide the histograms    
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1 / h2
    
    """
    return _h2_oper_ ( h1 , h2 , lambda x,y : x/y ) 
# =============================================================================
## Division with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_mul_ ( h1 , h2 ) :
    """ Multiply the histograms     
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1 * h2 
    """
    return _h2_oper_ ( h1 , h2 , lambda x,y : x*y ) 
# =============================================================================
## Addition with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_add_ ( h1 , h2 ) :
    """ Add the histograms     
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1 + h2 
    """
    return _h2_oper_ ( h1 , h2 , lambda x,y : x+y ) 
# =============================================================================
## Subtraction of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_sub_ ( h1 , h2 ) :
    """ Subtract the histogram    
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1 - h2 
    """
    return _h2_oper_ ( h1 , h2 , lambda x,y : x-y ) 


# =============================================================================
## Division with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_rdiv_ ( h1 , h2 ) :
    """ Divide the histograms    
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1 / h2
    
    """
    return _h2_oper_ ( h1 , h2 , lambda x,y : y/x ) 
# =============================================================================
## Division with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_rmul_ ( h1 , h2 ) :
    """ Multiply the histograms     
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1 * h2 
    """
    return _h2_oper_ ( h1 , h2 , lambda x,y : y*x ) 
# =============================================================================
## Addition with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_radd_ ( h1 , h2 ) :
    """ Add the histograms     
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1 + h2 
    """
    return _h2_oper_ ( h1 , h2 , lambda x,y : y+x ) 
# =============================================================================
## Subtraction of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_rsub_ ( h1 , h2 ) :
    """ Subtract the histogram    
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1 - h2 
    """
    return _h2_oper_ ( h1 , h2 , lambda x,y : y-x ) 


# =============================================================================
##  ``Fraction'' of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_frac_ ( h1 , h2 ) :
    """``Fraction'' the histogram h1/(h1+h2)    
    >>> h1     = ...
    >>> h2     = ...
    >>> frac   = h1.frac ( h2 )
    """
    return _h2_oper_ ( h1 , h2 , lambda x,y : x.frac(y) ) 
# =============================================================================
## ``Asymmetry'' of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_asym_ ( h1 , h2 ) :
    """``Asymmetry'' the histogram (h1-h2)/(h1+h2)    
    >>> h1     = ...
    >>> h2     = ...
    >>> asym   = h1.asym ( h2 )
    """
    return _h2_oper_ ( h1 , h2 , lambda x,y : x.asym(y) )


# =============================================================================
## ``minimum'' of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2023-08-11
def _h2_min_ ( h1 , h2 ) :
    """``minimum'' the histogram
    >>> h1     = ...
    >>> h2     = ...
    >>> asym   = h1.min ( h2 )
    """
    return _h2_oper_ ( h1 , h2 , lambda x,y : min ( x , y ) )

# =============================================================================
## ``maximum'' of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2023-08-11
def _h2_max_ ( h1 , h2 ) :
    """``maximum'' the histogram
    >>> h1     = ...
    >>> h2     = ...
    >>> asym   = h1.max ( h2 )
    """
    return _h2_oper_ ( h1 , h2 , lambda x,y : max ( x , y ) )


# =============================================================================
## ``Difference'' of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-06-03
def _h2_diff_ ( h1 , h2 ) :
    """``Difference'' the histogram 2*(h1-h2)/(h1+h2)    
    >>> h1     = ...
    >>> h2     = ...
    >>> diff   = h1.diff ( h2 )
    """
    return _h2_oper_ ( h1 , h2 , lambda x,y : 2*x.asym(y) ) 
# =============================================================================
##  ``Chi2-tension'' the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_chi2_ ( h1 , h2 ) :
    """``Chi2-tension'' for the histograms    
    >>> h1     = ...
    >>> h2     = ...
    >>> chi2   = h1.chi2 ( h2 ) 
    """
    return _h2_oper_ ( h1 , h2 , lambda x,y : VE ( x.chi2 ( y ) , 0 ) ) 

# =============================================================================
##  ``Average'' the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_average_ ( h1 , h2 ) :
    """``Average'' for the histograms    
    >>> h1     = ...
    >>> h2     = ...
    >>> mean   = h1.average ( h2 ) 
    """
    return _h2_oper_ ( h1 , h2 , lambda x,y : x.mean ( y ) )  

# =============================================================================
## 'pow' the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_pow_ ( h1 , val ) :
    """``pow'' the histogram    
    >>> h2     = ...
    >>> result = h2 ** 
    """
    if not h1.GetSumw2() : h1.Sumw2()
    #
    result = h1.Clone( hID() )
    if not result.GetSumw2() : result.Sumw2()
    #
    for ix1,iy1,x1,y1,z1 in h1.items() :
        #
        result.SetBinContent ( ix1 , iy1 , 0 ) 
        result.SetBinError   ( ix1 , iy1 , 0 )
        #
        v = VE ( pow ( z1 , val ) ) 
        #
        if not v.isfinite() : continue 
        #
        result.SetBinContent ( ix1 , iy1 , v.value () ) 
        result.SetBinError   ( ix1 , iy1 , v.error () )
        
    result.ResetStats() 
    return result 

# =============================================================================
## 'abs' for histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_abs_ ( h1 ) :
    """``abs'' the histogram
    >>> h2     = ...
    >>> result = abs ( h2 ) 
    """
    if not h1.GetSumw2() : h1.Sumw2()
    #
    result = h1.Clone( hID() )
    if not result.GetSumw2() : result.Sumw2()
    #
    for ix1,iy1,x1,y1,z1 in h1.items() :
        #
        result.SetBinContent ( ix1 , iy1 , 0 ) 
        result.SetBinError   ( ix1 , iy1 , 0 )
        #
        v = abs  ( z1 ) 
        #
        result.SetBinContent ( ix1 , iy1 , v.value () ) 
        result.SetBinError   ( ix1 , iy1 , v.error () )
        
    result.ResetStats() 
    return result 

# =============================================================================
## Division with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-06-03
def _h2_idiv_ ( h1 , h2 ) :
    """ Divide the histograms     
    >>> h1  = ...
    >>> h2  = ...
    >>> h1 /=  h2     
    """
    return _h2_ioper_ ( h1 , h2 , lambda x,y : x/y ) 
# =============================================================================
## Division with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-06-03
def _h2_imul_ ( h1 , h2 ) :
    """ Multiply the histograms     
    >>> h1  = ...
    >>> h2  = ...
    >>> h1 *=  h2     
    """
    return _h2_ioper_ ( h1 , h2 , lambda x,y : x*y ) 

# =============================================================================
## Addition with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-06-03
def _h2_iadd_ ( h1 , h2 ) :
    """ Add the histograms    
    >>> h1  = ...
    >>> h2  = ...
    >>> h1 +=  h2     
    """
    return _h2_ioper_ ( h1 , h2 , lambda x,y : x+y ) 
# =============================================================================
## Subtraction of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-06-03
def _h2_isub_ ( h1 , h2 ) :
    """ Subtract the histogram    
    >>> h1  = ...
    >>> h2  = ...
    >>> h1 -=  h2     
    """
    return _h2_ioper_ ( h1 , h2 , lambda x,y : x-y ) 
# =============================================================================

def _h2_box_   ( self , opts = '' ) : return self.Draw ( opts + ' box'   )
def _h2_lego_  ( self , opts = '' ) : return self.Draw ( opts + ' lego'  )
def _h2_surf_  ( self , opts = '' ) : return self.Draw ( opts + ' surf'  )
def _h2_surf2_ ( self , opts = '' ) : return self.Draw ( opts + ' surf2' )

# =============================================================================
## decorate 
for t in ( ROOT.TH2F , ROOT.TH2D ) : 
    
    t . _oper_       = _h2_oper_
    t . __truediv__  = _h2_div_
    t . __div__      = _h2_div_
    t . __mul__      = _h2_mul_
    t . __add__      = _h2_add_
    t . __sub__      = _h2_sub_
    t . __pow__      = _h2_pow_
    t . __abs__      = _h2_abs_

    t . __rtruediv__ = _h2_rdiv_
    t . __rdiv__     = _h2_rdiv_
    t . __rmul__     = _h2_rmul_
    t . __radd__     = _h2_radd_
    t . __rsub__     = _h2_rsub_

    t . __itruediv__ = _h2_idiv_
    t . __idiv__     = _h2_idiv_
    t . __imul__     = _h2_imul_
    t . __iadd__     = _h2_iadd_
    t . __isub__     = _h2_isub_
    
    t .  frac        = _h2_frac_
    t .  asym        = _h2_asym_
    t .  diff        = _h2_diff_
    t .  chi2        = _h2_chi2_
    t .  average     = _h2_average_
    t .  min         = _h2_min_
    t .  max         = _h2_max_

    t .  box         = _h2_box_
    t .  lego        = _h2_lego_
    t .  surf        = _h2_surf_
    t .  surf2       = _h2_surf2_



# =============================================================================
## operation with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_oper_ ( h1 , h2 , oper ) :
    """ Operation with the 3D-histogram     
    >>> h1 = ...
    >>> h2 = ...
    >>> h3 = h1 (oper) h2    
    """
    if                                 not h1.GetSumw2() : h1.Sumw2()
    if hasattr ( h2 , 'GetSumw2' ) and not h2.GetSumw2() : h2.Sumw2()
    #

    result = h1.Clone( hID() )
    if not result.GetSumw2() : result.Sumw2()
    #
    f2 = objectAsFunction ( h2 ) 
    #
    
    for ix1,iy1,iz1,x1,y1,z1,v1 in h1.items() :
        #
        result.SetBinContent ( ix1 , iy1 , iz1 , 0 ) 
        result.SetBinError   ( ix1 , iy1 , iz1 , 0 )
        #
        v2 = f2 ( x1.value() , y1.value() , z1.value() ) 
        #
        v  = VE ( oper ( v1 , v2 ) )
        #
        if not v.isfinite() : continue 
        #
        result.SetBinContent ( ix1 , iy1 , iz1 , v.value () ) 
        result.SetBinError   ( ix1 , iy1 , iz1 , v.error () )

    result.ResetStats() 
    return result



# =============================================================================
## operation with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-06-03
def _h3_ioper_ ( h1 , h2 , oper ) :
    """ Operation with the 3D-histogram 
    """
    if                                 not h1.GetSumw2() : h1.Sumw2()
    if hasattr ( h2 , 'GetSumw2' ) and not h2.GetSumw2() : h2.Sumw2()
    #
    f2 = objectAsFunction ( h2 ) 
    # 
    for ix1,iy1,iz1,x1,y1,z1,v1 in h1.items() :
        #
        h1.SetBinContent ( ix1 , iy1 , iz1 , 0 ) 
        h1.SetBinError   ( ix1 , iy1 , iz1 , 0 )
        #
        v2 = f2 ( x1.value() , y1.value() , z1.value() ) 
        #
        v  = VE ( oper ( v1 , v2 ) ) 
        #
        if not v.isfinite() : continue 
        #
        h1.SetBinContent ( ix1 , iy1 , iz1 , v.value () ) 
        h1.SetBinError   ( ix1 , iy1 , iz1 , v.error () )

    h1.ResetStats() 
    return h1

# =============================================================================
##  Division with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_div_ ( h1 , h2 ) :
    """ Divide the histograms     
    >>> h1 = ...
    >>> h2 = ...
    >>> h3 = h1 / h2 
    """
    return _h3_oper_ ( h1 , h2 , lambda x,y : x/y ) 
# =============================================================================
##  Division with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_mul_ ( h1 , h2 ) :
    """ Multiply the histograms    
    >>> h1 = ...
    >>> h2 = ...
    >>> h3 = h1 * h2 
    """
    return _h3_oper_ ( h1 , h2 , lambda x,y : x*y ) 
# =============================================================================
##  Addition with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_add_ ( h1 , h2 ) :
    """ Add the histograms     
    >>> h1 = ...
    >>> h2 = ...
    >>> h3 = h1 + h2 
    """
    return _h3_oper_ ( h1 , h2 , lambda x,y : x+y ) 
# =============================================================================
##  Subtraction of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_sub_ ( h1 , h2 ) :
    """ Subtract the histogram    
    >>> h1 = ...
    >>> h2 = ...
    >>> h3 = h1 - h2 
    """
    return _h3_oper_ ( h1 , h2 , lambda x,y : x-y ) 
# =============================================================================
##  ``Fraction'' of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_frac_ ( h1 , h2 ) :
    """``Fraction'' the histogram h1/(h1+h2)    
    >>> h1 = ...
    >>> h2 = ...
    >>> h3 = h1.frac ( h2 )    
    """
    return _h3_oper_ ( h1 , h2 , lambda x,y : x.frac(y) )

# =============================================================================
##  ``Asymmetry'' of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_asym_ ( h1 , h2 ) :
    """``Asymmetry'' the histogram (h1-h2)/(h1+h2)    
    >>> h1 = ...
    >>> h2 = ...
    >>> h3 = h1.asym ( h2 )    
    """
    return _h3_oper_ ( h1 , h2 , lambda x,y : x.asym(y) ) 


# =============================================================================
## ``minimum'' of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2023-08-11
def _h3_min_ ( h1 , h2 ) :
    """``minimum'' the histogram
    >>> h1     = ...
    >>> h2     = ...
    >>> asym   = h1.min ( h2 )
    """
    return _h23_oper_ ( h1 , h2 , lambda x,y : min ( x , y ) )

# =============================================================================
## ``maximum'' of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2023-08-11
def _h3_max_ ( h1 , h2 ) :
    """``maximum'' the histogram
    >>> h1     = ...
    >>> h2     = ...
    >>> asym   = h1.max ( h2 )
    """
    return _h3_oper_ ( h1 , h2 , lambda x,y : max ( x , y ) )


# =============================================================================
##  ``Chi2-tension'' the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_chi2_ ( h1 , h2 ) :
    """``Chi2-tension'' for the histograms    
    >>> h1 = ...
    >>> h2 = ...
    >>> h3 = h1.chi2 ( h2 )    
    """
    return _h3_oper_ ( h1 , h2 , lambda x,y : VE ( x.chi2 ( y ) , 0 ) ) 
# =============================================================================
##  ``Average'' the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_average_ ( h1 , h2 ) :
    """``Average'' for the histograms
    >>> h1 = ...
    >>> h2 = ...
    >>> h3 = h1.average ( h2 ) 
    """
    return _h3_oper_ ( h1 , h2 , lambda x,y : x.mean ( y ) )  



# =============================================================================
## 'pow' the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_pow_ ( h1 , val ) :
    """``pow'' the histogram     
    >>> h1      = ...
    >>> result  = h1 ** 3 
    """
    if not h1.GetSumw2() : h1.Sumw2()
    #
    result = h1.Clone( hID() )
    if not result.GetSumw2() : result.Sumw2()
    #
    for ix1,iy1,iz1,x1,y1,z1,v1 in h1.items() :
        #
        result.SetBinContent ( ix1 , iy1 , iz1 , 0 ) 
        result.SetBinError   ( ix1 , iy1 , iz1 , 0 )
        #
        v = VE ( pow ( v1 , val ) ) 
        #
        if not v.isfinite() : continue 
        #
        result.SetBinContent ( ix1 , iy1 , iz1 , v.value () ) 
        result.SetBinError   ( ix1 , iy1 , iz1 , v.error () )
        
    result.ResetStats() 
    return result 


# =============================================================================
## Division with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_rdiv_ ( h1 , h2 ) :
    """ Divide the histograms    
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1 / h2    
    """
    return _h3_oper_ ( h1 , h2 , lambda x,y : y/x ) 
# =============================================================================
## Division with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_rmul_ ( h1 , h2 ) :
    """ Multiply the histograms     
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1 * h2 
    """
    return _h3_oper_ ( h1 , h2 , lambda x,y : y*x ) 

# =============================================================================
## Addition with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_radd_ ( h1 , h2 ) :
    """ Add the histograms     
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1 + h2 
    """
    return _h3_oper_ ( h1 , h2 , lambda x,y : y+x ) 
# =============================================================================
## Subtraction of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_rsub_ ( h1 , h2 ) :
    """ Subtract the histogram    
    >>> h1     = ...
    >>> h2     = ...
    >>> result = h1 - h2 
    """
    return _h3_oper_ ( h1 , h2 , lambda x,y : y-x ) 

# =============================================================================
## 'abs' the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_abs_ ( h1 ) :
    """``abs'' the histogram
    >>> h3     = ...
    >>> result = abs ( h3 ) 
    """
    if not h1.GetSumw2() : h1.Sumw2()
    # 
    result = h1.Clone( hID() )
    if not result.GetSumw2() : result.Sumw2()
    #
    for ix1,iy1,iz1,x1,y1,z1,v1 in h1.items() :
        #
        result.SetBinContent ( ix1 , iy1 , iz1 , 0 ) 
        result.SetBinError   ( ix1 , iy1 , iz1 , 0 )
        #
        v = abs  ( v1 ) 
        #
        result.SetBinContent ( ix1 , iy1 , iz1 , v.value () ) 
        result.SetBinError   ( ix1 , iy1 , iz1 , v.error () )
        
    result.ResetStats() 
    return result 

# =============================================================================
## Division with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-06-03
def _h3_idiv_ ( h1 , h2 ) :
    """ Divide the histograms     
    >>> h1  = ...
    >>> h2  = ...
    >>> h1 /=  h2 
    """
    return _h3_ioper_ ( h1 , h2 , lambda x,y : x/y ) 
# =============================================================================
## Division with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-06-03
def _h3_imul_ ( h1 , h2 ) :
    """ Multiply the histograms     
    >>> h1  = ...
    >>> h2  = ...
    >>> h1 *=  h2     
    """
    return _h3_ioper_ ( h1 , h2 , lambda x,y : x*y ) 
# =============================================================================
## Addition with the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-06-03
def _h3_iadd_ ( h1 , h2 ) :
    """ Add the histograms    
    >>> h1  = ...
    >>> h2  = ...
    >>> h1 +=  h2     
    """
    return _h3_ioper_ ( h1 , h2 , lambda x,y : x+y ) 
# =============================================================================
## Subtraction of the histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-06-03
def _h3_isub_ ( h1 , h2 ) :
    """ Subtract the histogram    
    >>> h1  = ...
    >>> h2  = ...
    >>> h1 -=  h2     
    """
    return _h3_ioper_ ( h1 , h2 , lambda x,y : x-y ) 
# =============================================================================

ROOT.TH3._oper_    = _h3_oper_
    
ROOT.TH3.__truediv__  = _h3_div_
ROOT.TH3.__div__      = _h3_div_
ROOT.TH3.__mul__      = _h3_mul_
ROOT.TH3.__add__      = _h3_add_
ROOT.TH3.__sub__      = _h3_sub_

ROOT.TH3.__rtruediv__ = _h3_rdiv_
ROOT.TH3.__rdiv__     = _h3_rdiv_
ROOT.TH3.__rmul__     = _h3_rmul_
ROOT.TH3.__radd__     = _h3_radd_
ROOT.TH3.__rsub__     = _h3_rsub_

ROOT.TH3.__itruediv__ = _h3_idiv_
ROOT.TH3.__idiv__     = _h3_idiv_
ROOT.TH3.__imul__     = _h3_imul_
ROOT.TH3.__iadd__     = _h3_iadd_
ROOT.TH3.__isub__     = _h3_isub_

ROOT.TH3.__abs__   = _h3_abs_
ROOT.TH3.__pow__   = _h3_pow_


ROOT.TH3.  frac    = _h3_frac_
ROOT.TH3.  asym    = _h3_asym_
ROOT.TH3.  chi2    = _h3_chi2_
ROOT.TH3.  average = _h3_average_
ROOT.TH3.  min     = _h3_min_
ROOT.TH3.  max     = _h3_max_


# =============================================================================
## Update value of histogram from the function.
#  For each bin the integral of the function within the bin is used.
#  It is useful for wide or non-uniform bins.
#  @code
#  func     = math.exp
#  histo    = ...
#  histo.addFunctionIntegral( func )
#  histo   += math.exp                ## different!!!! 
#  @endcode
#  @attention It is a differrent from  <c>+=</c>
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-08-31
def _h1_add_function_integral_ ( h1 , func ) :
    """ Update value of histogram from the function.
    For each bin the integral of the function within the bin is used.
    It is useful for wide or non-uniform bins
    >>> func     = math.exp
    >>> histo    = ...
    >>> histo.addFunctionIntegral( func )  
    >>> histo   += math.exp                ## different!!!! 
    """
    ## new stuff 
    ff = objectAsFunctionObject( func )
    if not h1.GetSumw2() : h1.Sumw2()
    
    for i in h1.items() :

        ibin      = i [ 0 ]
        x         = i [ 1 ]
        y         = i [ 2 ]
        
        ## calculate the integral of the function withing the bin 
        ii        = ff.integral ( x.value() - x.error() , x.value()+x.error() )

        ## update 
        h1[ibin]  = y + ii

    h1.ResetStats() 
    return h1

ROOT.TH1F.addFunctionIntegral = _h1_add_function_integral_
ROOT.TH1D.addFunctionIntegral = _h1_add_function_integral_

# =============================================================================
## get the running sum over the histogram bins
#  @code
#  h   = ...
#  ht = h.sumv ( forward = True  )
#  hf = h.sumv ( forward = False )
#  @endcode
#  It creates a new historgma with the same binning
#  such that each bin<code>i</code> contains the sum over
#  all bins `j` such as 
#  - \f$ j\le i \f$ if <code>forward=True</code>
#  - \f$ j\ge i \f$ if <code>forward=False</code>
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def h1_sumv ( histo , forward = True ) :
    """ Create the `running sum' over the histogram bins 
    >>> h   = ...
    >>> h1 = h.sumv( forward = True  )
    >>> h2 = h.sumv( forward = False )        
    It creates a new histogram with the same binning
    such that each bin `i` contains the sum over all bins `j` such as 
    - `j<=i`  if `forward=True`
    - `i<=j`  if `forward=False`
    """
    assert isinstance ( histo , ROOT.TH1 ) and 1 == histo.dim() , \
        "Invalid `histo' type %s" % type ( histo ) 
                        
    result = histo.Clone ( hID () )
    result.Reset() 
    if not result.GetSumw2() : result.Sumw2()

    if forward :
        
        sumi = VE ( 0 , 0 ) 
        for i , x , y  in histo.items()  :
            sumi         += y 
            result [ i ]  = sumi

        return result

    sumi = VE ( 0 , 0 ) 
    for i in reversed ( histo ) :
        sumi += histo [ i ]
        result [ i ] = sumi 
        
    return result 

# ==============================================================================
for t in ( ROOT.TH1F , ROOT.TH1D ) : 
    t . sumv   = h1_sumv 

# =============================================================================
## Calculate the "cut-efficiency from the histogram
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_effic_ ( h , cut_low ) :
    """ Calculate the cut efficiency for the histogram
    >>> h  = ...
    >>> he = h.effic ( 14.2 , cut_low = True )    
    """
    
    result = h.Clone ( hID() )
    result.Reset() 
    if not result.GetSumw2() : result.Sumw2()

    for ibin in h :

        s1 = VE ( 0 , 0 )
        s2 = VE ( 0 , 0 )
        
        for jbin in h :
            
            if jbin < ibin : s1 += h [ jbin ]
            else           : s2 += h [ jbin ]

        result [ibin] = s2.frac ( s1 ) if cut_low else s1.frac ( s2 ) 

    result.ResetStats() 
    return result 

# =============================================================================
## Calculate the "cut-efficiency from the histogram
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_effic2_ ( h , value , cut_low ) :
    """ Calculate the cut efficiency for the histogram
    >>> h  = ...
    >>> he = h.efficiency ( 14.2 , cut_low = True )
    """

    s1 = VE ( 0 , 0 )
    s2 = VE ( 0 , 0 )

    for i,x,y in h.iteritems () :

        if x.value() < value : s1 += y 
        else                 : s2 += y 
        
    return s2.frac ( s1 ) if cut_low else s1.frac ( s2 ) 

# =============================================================================
## Convert historgam into "efficinecy" histogram
#  @code
#  histo = ...
#  effic = histo.eff ( ... )
#  @endcode
#  It adds two extra narrow fake bins!
def _h1_effic3_ ( h1 , cut_low ) :
    """ Convert historgam into "efficinecy" histogram
    >>> histo = ...
    >>> effic = histo.eff ( ... )
    - It adds two extra narrow fake bins!
    """

    xa    = h1.GetXaxis()
    edges = list ( xa.edges() )
    
    if h1.uniform_bins () :
        minbin =  ( xa.GetXmax() - xa.GetXmin() ) / xa.GetNbins()
    else :        
        nn     = len ( edges )
        minbin =  ( xa.GetXmax() - xa.GetXmin() )
        for i in xa :
            bw = xa.GetBinUpEdge() - xa.GetBinLowEdge()
            if abs ( bw ) < minbin : minbin = abs ( bw )
    
    xf = edges [  0 ] + 0.001 * minbin
    xl = edges [ -1 ] - 0.001 * minbin
    
    edges.insert (  1 , xf )
    edges.insert ( -1 , xl )

    result = h1_axis ( edges , title = 'Efficiency histo for %s' % h1.title , double = type ( h1 ) )

    def _my_eff_ ( a , r , c ) :
        
        R  = r + 0.5 * c
        A  = a + 0.5 * c
        
        fA = float ( A )
        fR = float ( R )

        if not fA :
            return 1.0 - _my_eff_ ( r , a , c ) 
        
        d  =   fR  / fA
        
        dr =   1.0 / fA
        da = - d   / fA 
        dc =  0.5  / fA - 0.5 * d / fA  
        
        c2 = dr * dr * r.cov2() + da  * da * a.cov2() + dc * dc * c.cov2()
        d  = VE ( d , c2 )
        
        return 1.0 / ( 1.0 + d )
    
    N    = len ( h1 )
    sumi = VE ( 0 , 0 ) 
    for i in h1 :        

        c = h1 [ i ]
        
        rest     = VE ( 0 , 0 ) 
        for j in reversed ( h1 ) :
            if i < j : rest += h1 [ j ]
            else     : break 

        a  = sumi
        r  = rest

        result [ i + 1 ] = _my_eff_ ( r , a , c ) if cut_low  else _my_eff_ ( a , r , c ) 
        
        sumi += c 
        
    if h1.natural () :
        nn = max ( 1 , int ( math.ceil ( float ( h1.Integral () ) ) ) ) 
        e0 = binomEff ( 0  , nn )
        e1 = binomEff ( nn , nn )
    else :
        e0 = VE ( 0 , 0 )
        e1 = VE ( 1 , 0 ) 
        
    result [  1 ] = e1 if cut_low else e0 
    result [ -1 ] = e0 if cut_low else e1 
        
    return result 

# ===============================================================================
## Get the cut efficiency in graph form 
#  - useful for efficiency visualisation
#  - a bit better treatment of binnig effects fo wide bins
#  @code
#  histo     = ...
#  eff_graph = histo.eff_graph ( cut_low )  
#  @endcode
def _h1_effic4_ ( histo , cut_low  ) :
    """ Get the cut efficiency in graph forms 
    - useful for drawing,
    - a bit better treatment of binnig effects for wide bins 
    >>> histo     = ...
    >>> eff_graph = histo.eff_graph ( cut_low = True  ) 
    >>> eff_graph = histo.eff_graph ( cut_low = False )    
    """

    c1 = [ histo [ i ] for i in histo ]
    c2 = c1 [ : ]
    c2.reverse  ()

    ## for "natural" histograms make better treatment of first/last uncertainties 
    vz = VE ( 0 , 1 ) if histo.natural() else VE ( 0 , 0 )

    s1 = [ vz ] + [ s for s in accumulate ( c1 ) ]
    s2 = [ vz ] + [ s for s in accumulate ( c2 ) ]

    s2.reverse  ()

    import ostap.histos.graphs
    
    np    = len ( s1 )
    graph = ROOT.TGraphErrors ( np  )

    the_eff = lambda a,b : a.frac ( b ) 
    
    ## special treatment of the first point 
    e0       = the_eff ( s2 [ 0 ] , s1 [ 0 ] ) if cut_low else the_eff ( s1 [ 0 ] , s2 [ 0 ] ) 
    xmin , _ = histo.xminmax()
    
    graph.SetPoint      ( 0 , xmin , e0.value () )
    graph.SetPointError ( 0 , 0    , e0.error () )

    for i, x , _  in histo.items() :
        xx = x.value() + x.error()
        ei = the_eff ( s2 [ i ] , s1 [ i ] ) if cut_low else the_eff ( s1 [ i ] , s2 [ i ] ) 
        graph.SetPoint      ( i , xx , ei.value () )
        graph.SetPointError ( i , 0  , ei.error () )
        
    return graph 
    
ROOT.TH1F.effic      = _h1_effic_ 
ROOT.TH1D.effic      = _h1_effic_ 
ROOT.TH1F.efficiency = _h1_effic2_ 
ROOT.TH1D.efficiency = _h1_effic2_ 
ROOT.TH1F.eff        = _h1_effic3_ 
ROOT.TH1D.eff        = _h1_effic3_ 
ROOT.TH1F.eff_graph  = _h1_effic4_
ROOT.TH1D.eff_graph  = _h1_effic4_ 

# ===============================================================================
## Convert the natural histogram into (p-value,significance) pair.
#  E.g. one can start frpm distirbtion of signal yields from
#  ensemble of presudoexperimens produced with `background-only` hypotehsis.
#  And using the obtained `pvalue` and `signif` histograms
#  one can convert the observed signal yeild into `p-value` or `significance'
#  @code
#  yields = ...
#  pvalue, signif = yields.significance ()
#  value = ...
#  print ( 'p-value      : %s ' % pvalue ( value ) )
#  print ( 'signigicance : %s ' % signif ( value ) ) 
#  @endcode
def _h1_signif_ ( h1 ) :
    """ Convert the natural histogram into (p-value,significance) pair: 
    - E.g. one can start from distirbtion of signal yields from
    ensemble of presudoexperimens produced with `background-only` hypotehsis.
    And using the obtained `pvalue` and `signif` histograms
    one can convert the observed signal yeild into `p-value` or `significance'
    >>> yields = ...
    >>> pvalue, signif = yields.significance ( cut_low = True ) 
    >>> print ( 'p-value      : %s ' % pvalue ( value ) ) 
    >>> print ( 'signigicance : %s ' % signif ( value ) ) 
    """
    
    if not h1.natural() : logger.warning ( "The input histogram is not `natural`!")
    ## get-pvalue 
    hp = h1.effic ( cut_low = True )
    hs = hp.clone ()

    the_last = VE ()
    for i, _ , y in hs.items() :
        if 0 < y.value() < 1 :
            s        = significance ( y ) 
            the_last = s
        else : s = the_last
        hs [ i ] = s 
        
    hp.blue ()
    hs.red  () 

    hp.SetTitle ( 'p-value'      ) 
    hs.SetTitle ( 'significance' )
    
    return hp, hs 

ROOT.TH1F.significance = _h1_signif_ 
ROOT.TH1D.significance = _h1_signif_ 

# ================================================================================
if  ( 2 , 7 ) <= python_info : _erfc_ = math.erfc 
else                         : _erfc_ = ROOT.Math.erfc
# ================================================================================
_sqrt_2_ = math.sqrt( 2.0 )
## helper function : convolution of gaussian with the single pulse 
def _cnv_ ( x , x0 , dx , sigma ) :
    """ Simple wrapper over error-function:
    convolution of gaussian with the single pulse 
    """
    #
    s = abs ( float ( sigma ) ) 
    #
    h = ( x - ( x0 + 0.5 * dx ) ) / ( _sqrt_2_ * s )
    l = ( x - ( x0 - 0.5 * dx ) ) / ( _sqrt_2_ * s )
    #
    return 0.5 * ( _erfc_ ( h ) - _erfc_ ( l ) ) / dx 

# =============================================================================
## "Smear" : make a convolution of the histogram with gaussian function
#  @param sigma     the gaussian resolutuon
#  @param addsigmas the parameter to treat the boundary conditions 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-09-20
def _smear_ ( h1 , sigma , addsigmas = 5 , silent = True ) :
    """ Smear the histogram through the convolution with the gaussian function:
    >>> histo    = ...
    >>> smeared1 = h.smear ( sigma = 0.01 )  ## smear using 'constant' sigma
    >>> smeared2 = h.smear ( sigma = lambda x :  0.1*x ) ## smear using 'running' sigma  
    """

    if   isinstance ( sigma , num_types ) :
        if 0 > sigma or iszero ( sigma ) : return h1.clone()  ## no need to smear 
        sig    = sigma  
        sigma  = lambda x : sig
    elif isinstance ( sigma , VE ) :
        sigma  = float ( sigma )  
        if 0 > sigma or iszero ( sigma ) : return h1.clone()  ## no need to smear 
        sig    = sigma  
        sigma  = lambda x : sig 
    elif callable  ( sigma ) :
        sigfun = sigma
        sigma  = lambda x : max ( float ( sigfun ( float ( x ) ) ) , 0.0 ) 
    else :
        raise TypeError("Can't deduce the resolution type %s/%s" % ( sigma , type ( sigma ) ) )
    
    #
    ## clone the input histogram
    #
    h2  = h1.Clone ( hID () ) ; h2.Reset() ;
    if not h2.GetSumw2() : h2.Sumw2()

    first_bin = None
    last_bin  = None

    ## collect all bins
    bins = []

    for ibin in h1.items() : bins.append ( ibin[1:] )

    ## add few artificially replicated bins
    fb    = bins [  0 ]
    lb    = bins [ -1 ]
    xmin  = fb[0].value() - fb[0].error()
    xmax  = lb[0].value() + lb[0].error()
    
    s_xmin = max ( sigma ( fb[0].value() ) , 2*fb[0].error() )  
    s_xmax = max ( sigma ( lb[0].value() ) , 2*lb[0].error() )  
    
    x_min = xmin - abs ( addsigmas * s_xmin ) 
    x_max = xmax + abs ( addsigmas * s_xmax ) 
    
    ## add few fictive bins 
    while  xmin > x_min :
        
        bin0 = bins[0]
        
        xc   = bin0[0]
        val  = bin0[1]
        
        bin  = (xc-2*xc.error(),val)
        bins.insert (  0, bin )
        fb    = bins [  0 ]
        
        xmin  = bin[0].value()-bin[0].error()

    ## add few fictive bins 
    while  xmax < x_max :
        
        bin0 = bins[-1]
        
        xc   = bin0[0]
        val  = bin0[1]
        
        bin  = (xc+2*xc.error(),val)
        bins.append ( bin )
        fb    = bins [  0 ]
        
        xmax  = bin[0].value()+bin[0].error()
    
    for ibin1 in progress_bar ( bins , silent = silent , description = 'Bins:'  ) :

        x1c  =     ibin1 [0].value () 
        x1w  = 2 * ibin1 [0].error () 
        val1 =     ibin1 [1]

        ## nothing to smear 
        if iszero ( val1.value() ) : continue

        ## calculate the resolution  
        sig  = abs  ( float ( sigma ( x1c ) ) )
        
        ## if very narrow, skip it...
        if sig < 0.01 * x1w  :
            
            i2 = h2.findBin  ( x1c )
            if i2 in h2 : h2 [ i2 ] += val1
            
            continue 
        
        for ibin2 in h2.items() :

            i2      =     ibin2 [0]
            x2c     =     ibin2 [1].value () 
            x2w     = 2 * ibin2 [1].error () 

            val2    = VE ( val1 )
            val2   *= x2w
            val2   *= _cnv_ ( x2c , x1c , x1w , sig ) 
            
            h2[i2] += val2

    h2.ResetStats() 
    return h2


ROOT.TH1F. smear = _smear_
ROOT.TH1D. smear = _smear_

# =============================================================================
## Smooth the histogram
#  @see TH1::Smooth
def  _h1_smooth_ ( h1 , ntimes = 1 ) :
    """ Smooth the histogram
    - see `ROOT.TH1.Smooth`
    """
    assert isinstance ( ntimes , integer_types ) and 1 <= ntimes ,\
           "Invalid `ntimes` parameter! "
    result = h1.clone() 
    result .Smooth ( ntimes )
    return result

ROOT.TH1F. smooth = _h1_smooth_
ROOT.TH1D. smooth = _h1_smooth_

# =============================================================================
## Get the integral with error
def _h1_integral_err_ ( h1 ) :
    """ Get the integral with error"""
    
    xaxis = h1.GetXaxis() 
    err   = ctypes.c_double ( 0.0 ) 
    value = h1.IntegralAndError ( xaxis.GetFirst() , xaxis.GetLast() , err )
    error = err.value 
    return VE ( value , error ** 2 )

# =============================================================================
### Get the integral with error
def _h2_integral_err_ ( h2 ) :
    """ Get the integral with error"""
    
    xaxis = h2.GetXaxis() 
    yaxis = h2.GetYaxis() 
    err   = ctypes.c_double ( 0.0 ) 
    value = h2.IntegralAndError ( xaxis.GetFirst() , xaxis.GetLast() , 
                                  yaxis.GetFirst() , yaxis.GetLast() , err )
    error = err.value 
    return VE ( value , error ** 2 )

# =============================================================================
### Get the integral with error
def _h3_integral_err_ ( h3 ) :
    """ Get the integral with error"""
    
    xaxis = h3.GetXaxis() 
    yaxis = h3.GetYaxis() 
    zaxis = h3.GetYaxis() 
    err   = ctypes.c_double ( 0.0 ) 
    value = h23.IntegralAndError ( xaxis.GetFirst() , xaxis.GetLast() , 
                                   yaxis.GetFirst() , yaxis.GetLast() ,
                                   zaxis.GetFirst() , zaxis.GetLast() , err )
    error = err.value 
    return VE ( value , error ** 2 ) 

ROOT.TH1F.the_integral = _h1_integral_err_
ROOT.TH1D.the_integral = _h1_integral_err_

ROOT.TH2F.the_integral = _h2_integral_err_
ROOT.TH2D.the_integral = _h2_integral_err_

ROOT.TH3F.the_integral = _h3_integral_err_
ROOT.TH3D.the_integral = _h3_integral_err_

# =============================================================================
## make transformation of histogram content 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-10-23  
def _h1_transform_ ( h1 , func ) :
    """ Make the transformation of the histogram content     
    >>> func = lambda x,y: y   ## identical transformation/copy
    >>> h1 = ...
    >>> h2 = h1.fransform ( func )     
    """
    #
    if not h1.GetSumw2() : h1.Sumw2()
    h2 = h1.Clone( hID() )
    if not h2.GetSumw2() : h2.Sumw2()
    #
    for i,x,y in h1.items() :
        
        h2 [ i ] = func ( x, y ) 
        
    h2.ResetStats() 
    return h2 

ROOT.TH1F. transform = _h1_transform_ 
ROOT.TH1D. transform = _h1_transform_ 

import ostap.math.math_ve as math_ve
for h in ( ROOT.TH1F , ROOT.TH1D ) : 
    
    h. __exp__    = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.exp    ( y ) ) 
    h. __expm1__  = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.expm1  ( y ) ) 
    h. __log__    = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.log    ( y ) )
    h. __log10__  = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.log10  ( y ) )
    h. __log1p__  = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.log1p  ( y ) )
    h. __sqrt__   = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.sqrt   ( y ) )
    h. __cbrt__   = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.cbrt   ( y ) )
    h. __sin__    = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.sin    ( y ) )
    h. __cos__    = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.cos    ( y ) )
    h. __tan__    = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.tan    ( y ) )
    h. __sinh__   = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.sinh   ( y ) )
    h. __cosh__   = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.cosh   ( y ) )
    h. __tanh__   = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.tanh   ( y ) )
    h. __sech__   = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.sech   ( y ) )
    h. __asin__   = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.asin   ( y ) )
    h. __acos__   = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.acos   ( y ) )
    h. __atan__   = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.atan   ( y ) )
    h. __asinh__  = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.asinh  ( y ) )
    h. __acosh__  = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.acosh  ( y ) )
    h. __atanh__  = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.atanh  ( y ) )
    h. __erf__    = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.erf    ( y ) )
    h. __erfc__   = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.erfc   ( y ) )
    h. __erfi__   = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.erfi   ( y ) )
    h. __erfcx__  = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.erfcx  ( y ) )
    h. __gamma__  = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.gamma  ( y ) )
    h. __tgamma__ = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.tgamma ( y ) )
    h. __lgamma__ = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.lgamma ( y ) )
    h. __igamma__ = lambda s : _h1_transform_ ( s , lambda x,y : math_ve.igamma ( y ) )
    
    
# =============================================================================
## make transformation of histogram content 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-10-23  
def _h2_transform_ ( h2 , func ) :
    """ Make the transformation of the 2D-histogram content 
    >>> func = lambda x,y,z: z   ## identical transformation/copy
    >>> h2   = ...
    >>> h3   = h2.fransform ( func ) 
    """
    #
    if not h2.GetSumw2() : h2.Sumw2()
    h3 = h2.Clone( hID() )
    if not h3.GetSumw2() : h3.Sumw2()
    #
    for ix,iy,x,y,z in h2.items() :
        
        h3 [ ix , iy ] = func ( x, y , z ) 
        
    h3.ResetStats() 
    return h3 

ROOT.TH2F. transform = _h2_transform_ 
ROOT.TH2D. transform = _h2_transform_ 

import ostap.math.math_ve as math_ve
for h in ( ROOT.TH2F , ROOT.TH2D ) : 
    
    h. __exp__    = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.exp    ( z ) ) 
    h. __expm1__  = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.expm1  ( z ) ) 
    h. __log__    = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.log    ( z ) )
    h. __log10__  = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.log10  ( z ) )
    h. __logp1__  = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.logp1  ( z ) )
    h. __sqrt__   = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.sqrt   ( z ) )
    h. __cbrt__   = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.cbrt   ( z ) )
    h. __sin__    = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.sin    ( z ) )
    h. __cos__    = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.cos    ( z ) )
    h. __tan__    = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.tan    ( z ) )
    h. __sinh__   = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.sinh   ( z ) )
    h. __cosh__   = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.cosh   ( z ) )
    h. __tanh__   = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.tanh   ( z ) )
    h. __asin__   = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.asin   ( z ) )
    h. __acos__   = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.acos   ( z ) )
    h. __atan__   = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.atan   ( z ) )
    h. __asinh__  = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.asinh  ( z ) )
    h. __acosh__  = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.acosh  ( z ) )
    h. __atanh__  = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.atanh  ( z ) )
    h. __erf__    = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.erf    ( z ) )
    h. __erfc__   = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.erfc   ( z ) )
    h. __gamma__  = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.gamma  ( z ) )
    h. __tgamma__ = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.tgamma ( z ) )
    h. __lgamma__ = lambda s : _h2_transform_ ( s , lambda x,y,z : math_ve.lgamma ( z ) )

# =============================================================================
# Few "specific" transformations
# =============================================================================

# =============================================================================
## transform the histogram into "precision" histogram 
def _h1_precision_ ( self ) :
    """ Make precision histogram, each bin constains ``precision''
    >>> h =
    >>> p = h.precision()    
    """
    return _h1_transform_ ( self ,  lambda x,y  : y.precision() )
# =============================================================================
## transform the histogram into "precision" histogram 
def _h2_precision_ ( self ) :
    """ Make precision histogram, each bin constains ``precision''
    >>> h =
    >>> p = h.precision()
    """
    return _h2_transform_ ( self ,  lambda x,y,z: z.precision() )
# =============================================================================
## transform the histogram into "B/S" histogram 
def _h1_b2s_       ( self ) :
    """ Make B/S histogram, each bin constains ``B/S''
    >>> h    =
    >>> btos = h.b2s()    
    """
    return _h1_transform_ ( self ,  lambda x,y  : y.b2s () )
# =============================================================================
## transform the histogram into "B/S" histogram 
def _h2_b2s_       ( self ) :
    """ Make B/S histogram, each bin constains ``B/S''
    >>> h    =
    >>> btos = h.b2s()    
    """
    return _h2_transform_ ( self ,  lambda x,y,z: z.b2s () )

ROOT.TH1F. precision = _h1_precision_
ROOT.TH1D. precision = _h1_precision_
ROOT.TH2F. precision = _h2_precision_
ROOT.TH2D. precision = _h2_precision_
ROOT.TH1F. b2s       = _h1_b2s_
ROOT.TH1D. b2s       = _h1_b2s_
ROOT.TH2F. b2s       = _h2_b2s_
ROOT.TH2D. b2s       = _h2_b2s_

# =============================================================================
## rescale the histogram for effective uniform bins
#  new_content = old_content * factor / bin_width
#  @code
#  >>> h1 = ...
#  >>> h2 = h1.rescale_bins ( h1 , 1 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-03-19   
def _h1_rescale_ ( h1 , factor = 1 ) :
    """ Rescale the histogram for effective uniform bins :
    new_bin_content = old_bin_content * factor / bin_width
    >>> h1 = ...
    >>> h2 = h1.rescale_bins ( h1 , 1 )    
    """
    return _h1_transform_ ( h1  ,  lambda x , y : ( 0.5 * factor / x.error() ) * y )

ROOT.TH1F. rescale_bins = _h1_rescale_ 
ROOT.TH1D. rescale_bins = _h1_rescale_ 

# =============================================================================
## rescale the histogram for effective uniform bins
#  new_content = old_content * factor / bin_area
#  @code
#  >>> h2 = ...
#  >>> h3 = h2.rescale_bins ( h2 , 1 )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-03-19   
def _h2_rescale_ ( h2 , factor = 1 ) :
    """ Rescale the histogram for effective uniform bins : 
    new_bin_content = old_bin_content * factor / bin_area
    >>> h1 = ...
    >>> h2 = h1.rescale_bins ( h1 , 1 )
    
    """
    return _h2_transform_ ( h2  ,  lambda x , y , z : ( 0.25 * factor / x.error() / y.error() ) * z )

ROOT.TH2F. rescale_bins = _h2_rescale_ 
ROOT.TH2D. rescale_bins = _h2_rescale_ 

# =============================================================================
## sample the histogram using gaussian hypothesis
#  @code
#  >>> h = ...          ## original histogram
#  >>> s1 = h.sample ()  ## the sampled histogram
#  >>> s2 = h.gauss  ()  ## ditto 
#  >>> s3 = h.gauss  ( accept = lambda s : 0<s ) ## only positive numbers 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
def _h_sample_ ( histo , accept = lambda s : True , nmax = 1000 ) :
    """ Sample the histogram using gaussian hypothesis 
    >>> h  = ...         ## original histogram 
    >>> s1 = h.sample()  ## the sampled hist
    """
    #
    result = histo.Clone ( hID () )     
    result.Reset() 
    if not result.GetSumw2() : result.Sumw2()
    
    for bin in histo :

        ## getbin content
        v1 = histo[bin]
        
        ## sample it! 
        v2 = VE    ( v1.gauss ( accept = accept , nmax = nmax ) )
        
        v2.setCov2 ( v1.cov2() )
        
        result [ bin ] = v2
        
    result.ResetStats() 
    return result

ROOT.TH1 ._sample_ = _h_sample_
ROOT.TH1 .sample   = _h_sample_
ROOT.TH1 .gauss    = _h_sample_

# =============================================================================
## sample the histogram using the Poisson hypothesis
#  @code
#  histo = ...
#  h1 = histo.poisson () ## sample using poisson hypothesis
#  h2 = histo.poisson ( fluctuate =   True ) ## fluctuate Poisson mean 
#  h3 = histo.poisson ( accept = lambda s : s>10 ) ## require all numbers are large
#  @endcode
def _h_poisson_ ( histo , fluctuate = False , accept = lambda s : True ) :
    """ Sample the histogram using the Poisson hypothesis
    >>> histo = ...
    >>> h1 = histo.poisson () ## sample using poisson hypothesis
    >>> h2 = histo.poisson ( fluctuate =   True ) ## fluctuate Poisson mean 
    >>> h3 = histo.poisson ( accept = lambda s : s>10 ) ## require all numbers are large
    """
    result = histo.Clone ( hID () )
    result.Reset() 
    if not result.GetSumw2() : result.Sumw2()
    
    for bin in histo :

        ## getbin content
        v1 = histo [ bin ]
        
        ## sample it!
        p  = v1.poisson ( fluctuate , accept )

        v2 = VE ( p , p  )
        
        result [ bin ] = v2
        
    result.ResetStats() 
    return result
    
ROOT.TH1 .poisson = _h_poisson_

# =============================================================================
## Get the Figure-of-Merit (FoM) for the pure signal distribution,
#  e.g. from sPlot)
#  the FoM is defined from the relative precision of the signal yield
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-10-23  
def _fom_2_ ( h1 , increase = True ) :
    """ Get figure-of-merit (FOM) distribution for signal
    >>> h1 = ...  ## signal distribution
    >>> f1 = h1.FoM2 () 
    """
    #
    h = h1.sumv( forward = increase )
    #
    return _h1_transform_ ( h , func = lambda x,y : y.precision() )

# =============================================================================
## Calculate S/sqrt(S+a*B)
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-10-23   
def _sb_ ( s , b  , a = 1 ) :
    """ Calculate S/sqrt(S+a*B) 
    """
    #
    v     = ( s + a * b ).value()
    if 0 >= v : return VE ( 0 , 0 )
    #
    F     = s.value() / pow ( v , 0.5 )
    #
    # (dF/dS)**2 
    dFdS2 = ( 0.5 * s  + a * b ).value() 
    dFdS2 = dFdS2**2  / v**3
    #
    # (dR/dB)**2 
    dFdB2 = ( -0.5 * a * s  ).value() 
    dFdB2 = dFdB2**2  / v**3
    #
    return VE ( F , dFdS2 * s.cov2() + dFdB2 * b.cov2() ) 

# =============================================================================
## Get the figure-of-merit (FoM) for the signal and background distributions 
#  the FoM is defined as S/sqrt(S+alpha*B) 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-10-23  
def _fom_1_ ( s , b , alpha = 1 , increase = True ) :
    """ Get figure-of-merit FoM = S/sqrt(S+a*B)
    >>> s = ... ## signal distributions
    >>> b = ... ## background distributions
    >>> fom = s.FoM1( b , alpha = 1.0 )      
    """
    #
    if not s.GetSumw2() : s.Sumw2()
    h  = s.Clone( hID() )
    if not h.GetSumw2() : h.Sumw2()
    #
    hs = s.sumv ( forward = increase )
    hb = b.sumv ( forward = increase )
    #
    from math import sqrt, pow 
    #
    for i,x,y in hs.items() :
        
        ## the signal 
        si = y
        
        ## the background 
        bi = hb ( x ) 
        
        h [i] = _sb_ ( si , bi , alpha ) 
        
    h.ResetStats() 
    return h 

ROOT.TH1D . fom_1 = _fom_1_ 
ROOT.TH1D . fom_2 = _fom_2_ 
ROOT.TH1F . fom_1 = _fom_1_ 
ROOT.TH1F . fom_2 = _fom_2_

ROOT.TH1D . FoM_1 = _fom_1_ 
ROOT.TH1D . FoM_2 = _fom_2_ 
ROOT.TH1F . FoM_1 = _fom_1_ 
ROOT.TH1F . FoM_2 = _fom_2_ 


# =============================================================================
# rebin the histograms
# =============================================================================
## get the overlap for 1D-bins 
def _bin_overlap_1D_ ( x1 , x2 ) :
    """
    """
    #
    x1v = x1.value()
    x1e = x1.error()
    #
    xmin_1 = x1v - x1e 
    xmax_1 = x1v + x1e
    if xmin_1 >= xmax_1  : return 0         ## RETURN
    #
    x2v = x2.value()
    x2e = x2.error()
    #
    xmin_2 = x2v - x2e 
    xmax_2 = x2v + x2e
    if xmin_2 >= xmax_2  : return 0         ## RETURN
    #
    xmin = max ( xmin_1 , xmin_2 )
    xmax = min ( xmax_1 , xmax_2 )
    #
    if xmin >= xmax      : return 0         ## RETURN 
    #
    return ( xmax - xmin ) / ( 2.0 * x1e ) 
# =============================================================================
## get the overlap for 2D-bins 
def _bin_overlap_2D_ ( x1 , y1 , x2 , y2 ) :
    """
    """
    #
    x1v = x1.value()
    x1e = x1.error()
    #
    xmin_1 = x1v - x1e
    xmax_1 = x1v + x1e
    if xmin_1 >= xmax_1  : return 0         ## RETURN
    #
    x2v = x2.value()
    x2e = x2.error()    
    #
    xmin_2 = x2v - x2e 
    xmax_2 = x2v + x2e
    if xmin_2 >= xmax_2  : return 0         ## RETURN
    #
    y1v = y1.value()
    y1e = y1.error()    
    #
    ymin_1 = y1v - y1e 
    ymax_1 = y1v + y1e
    if ymin_1 >= ymax_1  : return 0         ## RETURN
    #
    y2v = y2.value()
    y2e = y2.error()        
    #
    ymin_2 = y2v - y2e 
    ymax_2 = y2v + y2e
    if ymin_2 >= ymax_2  : return 0         ## RETURN
    #
    xmin = max ( xmin_1 , xmin_2 )
    xmax = min ( xmax_1 , xmax_2 )
    if xmin >= xmax      : return 0         ## RETURN 
    #
    ymin = max ( ymin_1 , ymin_2 )
    ymax = min ( ymax_1 , ymax_2 )
    if ymin >= ymax      : return 0         ## RETURN 
    #
    #
    return ( xmax - xmin ) * ( ymax - ymin ) / ( 4.0 * x1e * y1e )


# ==============================================================================
## rebin 1D-histogram with NUMBERS 
def _rebin_nums_1D_ ( h1 , template ) :
    """ Rebin 1D-histogram assuming it is a histogram with *NUMBERS*
    >>> horig    = ...  ## the original histogram 
    >>> template = ...  ## the template with binnings
    >>> h = horig.rebinNumbers ( template ) 
    """
    ##
    if   isinstance ( template , ROOT.TH1 ) : pass 
    elif isinstance ( template  , integer_types  ) and 1 < template :
        xaxis    = h1.GetXaxis()
        htype    = ROOT.TH1D if isinstance ( h1  , ROOT.TH1D ) else ROOT.TH1F
        template = htype  ( hID() , 'template' , template , xaxis.GetXmin() , xaxis.GetXmax() )
    elif isinstance ( template  , sized_types      ) and \
             isinstance ( template  , sequence_types   ) and \
             3 == len   ( template  )                    and \
             isinstance ( template [0] , integer_types ) and  1 < template [0] and \
             isinstance ( template [1] , num_types     ) and \
             isinstance ( template [2] , num_types     ) and template [ 1] < template [ 2 ] :
        htype    = ROOT.TH1D if isinstance ( h1  , ROOT.TH1D ) else ROOT.TH1F
        template = htype  ( hID() , 'template' , *template )
        
    # clone it!
    h2 = template.Clone( hID() )
    if not h2.GetSumw2() : h2.Sumw2()
    #
    ## reset the histogram 
    for i2 in h2 : h2[i2] = VE(0,0)
    #
    for i2 in h2.items() :

        xb  = i2[1]
        xbv = xb.value ()
        xbe = xb.error ()
        
        bl = h1.findBin ( xbv - xbe ) - 1
        bh = h1.findBin ( xbv + xbe ) + 1
        
        for i1 in h1.items( bl , bh + 1 ) :
            
            o = _bin_overlap_1D_ ( i1[1] , i2[1] )
            
            h2 [ i2[0] ] +=  o * i1[2] 

    h2.ResetStats() 
    return h2 
# =============================================================================
## rebin 1D-histogram as FUNCTION 
def _rebin_func_1D_ ( h1 , template ) :
    """ Rebin 1D-histogram assuming it is a FUNCTION
    >>> horig    = ...  ## the original histogram 
    >>> template = ...  ## the template with binnings
    >>> h = horig.rebinFunction ( template ) 
    """
    
    if   isinstance ( template , ROOT.TH1 ) : pass 
    elif isinstance ( template  , integer_types  ) and 1 < template :
        xaxis    = h1.GetXaxis()
        htype    = ROOT.TH1D if isinstance ( h1  , ROOT.TH1D ) else ROOT.TH1F
        template = htype  ( hID() , 'template' , template , xaxis.GetXmin() , xaxis.GetXmax() )
    elif isinstance ( template  , sized_types      ) and \
             isinstance ( template  , sequence_types   ) and \
             3 == len   ( template  )                    and \
             isinstance ( template [0] , integer_types ) and  1 < template [0] and \
             isinstance ( template [1] , num_types     ) and \
             isinstance ( template [2] , num_types     ) and template [ 1] < template [ 2 ] :  
        htype    = ROOT.TH1D if isinstance ( h1  , ROOT.TH1D ) else ROOT.TH1F
        template = htype  ( hID() , 'template' , *template )
        
    # clone it!
    h2 = template.Clone( hID() )
    if not h2.GetSumw2() : h2.Sumw2()
    ## reset the histogram 
    for i2 in h2 : h2[i2] = VE(0,0)
    #
    for i2 in h2.items() :
        
        xb  = i2[1]
        xbv = xb.value ()
        xbe = xb.error ()
        
        bl = h1.findBin ( xbv - xbe ) - 1  
        bh = h1.findBin ( xbv + xbe ) + 1
        
        for i1 in h1.items( bl , bh + 1 ) :

            o = _bin_overlap_1D_ ( i2[1] , i1[1] ) ## NOTE THE ORDER!!! 
            
            h2 [ i2[0] ] +=  o * i1[2]
            
    h2.ResetStats() 
    return h2 


# ==============================================================================
## rebin 2D-histogram with NUMBERS 
def _rebin_nums_2D_ ( h1 , template ) :
    """ Rebin 2D-histogram assuming it is a histogram with *NUMBERS*
    >>> horig    = ...  ## the original histogram 
    >>> template = ...  ## the template with binnings
    >>> h = horig.rebinNumbers ( template ) 
    """
    #
    # clone it!
    h2 = template.Clone( hID() )
    if not h2.GetSumw2() : h2.Sumw2()
    ## reset the histogram 
    for i2 in h2 : h2[i2] = VE(0,0)
    #
    for i2 in h2.items() :

        for i1 in h1.items() :
            
            o = _bin_overlap_2D_ ( i1[2] , i1[3] , i2[2] , i2[3] )
            
            h2 [ i2[0] , i2[1] ] +=  o * i1[4] 
             
    h2.ResetStats() 
    return h2

# =============================================================================
## rebin 2D-histogram as FUNCTION 
def _rebin_func_2D_ ( h1 , template ) :
    """ Rebin 2D-histogram assuming it is a FUNCTION

    >>> horig    = ...  ## the original histogram 
    >>> template = ...  ## the template with binnings

    >>> h = horig.rebinFunction ( template ) 
    """
    # clone it!
    h2 = template.Clone( hID() )
    if not h2.GetSumw2() : h2.Sumw2()
    ## reset the histogram 
    for i2 in h2 : h2[i2] = VE(0,0)
    #
    for i2 in h2.items() :
        
        for i1 in h1.items() :
            
            o = _bin_overlap_1D_ ( i2[2] , i2[3] , i1[2] , i2[3] ) ## NOTE THE ORDER!!! 
            
            h2 [ i2[0] , i2[1] ] +=  o * i1[4]
            
    h2.ResetStats() 
    return h2 

for t in ( ROOT.TH1F , ROOT.TH1D ) :
    t.rebinNumbers  = _rebin_nums_1D_
    t.rebinFunction = _rebin_func_1D_
    
for t in ( ROOT.TH2F , ROOT.TH2D ) :
    t.rebinNumbers  = _rebin_nums_2D_
    t.rebinFunction = _rebin_func_2D_

# =============================================================================
## Create NULL-line for the histogram and (optionally) draw it
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-11-01
def _h1_null_ ( h1 , draw = False , style = 1 ) :
    """ Create NULL-line for the histogram and (optionally) draw it    
    """
    axis = h1.GetXaxis()
    line = ROOT.TLine ( axis.GetXmin() , 0 ,
                        axis.GetXmax() , 0 )
    
    line.SetLineStyle ( style )
    
    if draw :
        line.Draw()

    return line 

ROOT.TH1D.null = _h1_null_
ROOT.TH1F.null = _h1_null_


# =============================================================================
## get edges from the axis:
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _edges_ ( axis ) :
    """ Get list of edges from the TAxis
    >>> axis
    >>> edges = axis.edges() 
    """
    #
    bins  = [ axis.GetBinLowEdge ( i ) for i in axis ] + [ axis.GetXmax() ]
    #
    return tuple ( bins )

# =============================================================================
ROOT.TAxis.edges = _edges_

# =============================================================================
## get "slice" for 1D histogram
#  @code     
#    >>> h1 = ...
#    >>> nh = h1[2:10] ## keep only bins from 2nd (inclusive) till 10 (exclusive)
#  @endcode     
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-20 
def _h1_getslice_ ( h1 , i , j ) :
    """ Get the ``slice'' for 1D-histogram:
    >>> h1 = ...
    >>> nh = h1[2:10] ## keep only bins from 2nd (inclusive) till 10 (exclusive)    
    """
    axis = h1  .GetXaxis()
    nb   = axis.GetNbins()

    if i is None : i = 1
    if j is None : j = nb + 1 
    
    while i < 0 : i += nb
    while j < 0 : j += nb
    
    i = max ( 1 , min ( nb + 1 , i ) ) 
    j = max ( 1 , min ( nb + 1 , j ) ) 

    if i >= j :
        raise IndexError 
    
    edges = axis.edges ()
    edges = edges [i-1:j]
    
    typ = h1.__class__
    result = typ ( hID  ()       ,
                   h1.GetTitle() ,
                   len ( edges ) - 1 , array.array ( 'd' , edges ) )
    
    result.Sumw2()
    result += h1
    
    return result 

ROOT.TH1F  . __getslice__  =   _h1_getslice_ 
ROOT.TH1D  . __getslice__  =   _h1_getslice_ 

# =============================================================================
## Create <code>TAxis</code>
#  @code
#  a = make_axis ( 2 , 0.0 , 1.0 )                  ## #bins ,min, max 
#  a = make_axis ( 4 , 0.0 , 1.0 , 2.0 , 5.0 , 10 ) ## #bins , e1 , e2 , e3 , ...    
#  a = make_axis ( 2 , [ 0.0 , 1.0 , 2.0 ]  )       ## #bins , [ e1 , e2 . ... ]  
#  a = make_axis ( [ 0.0 , 1.0 , 2.0 ]  )           ## [ e1 , e2 , ... ] 
#  a = make_axis ( 0.0 , 1.0 , 2.0  )               ##   e1 , e2 , ...  
#  @endcode
def make_axis ( nbins , *bins ) :
    """ Create <code>TAxis</code>    
    >>> a = make_axis ( 2 , 0.0 , 1.0 )                  ## #bins ,min, max 
    >>> a = make_axis ( 4 , 0.0 , 1.0 , 2.0 , 5.0 , 10 ) ## #bins , e1 , e2 , e3 , ...    
    >>> a = make_axis ( 2 , [ 0.0 , 1.0 , 2.0 ]  )       ## #bins , [ e1 , e2 . ... ]  
    >>> a = make_axis ( [ 0.0 , 1.0 , 2.0 ]  )           ## [ e1 , e2 , ... ] 
    >>> a = make_axis (  0.0 , 1.0 , 2.0 , 3.0 , 4.0 )   ##   e1 , e2 , ... 
    """
    
    if   isinstance ( nbins , ROOT.TAxis ) : return nbins 

    if isinstance ( nbins , integer_types ) and 0 < nbins :
        
        if 1 == len ( bins ) and isinstance ( bins [ 0 ] , sequence_types ) :
            abins = [ float ( e ) for e in bins [ 0 ] ]
            if nbins +1 <= len ( abins ) :
                abins = abins [ : nbins + 1 ] 
                if is_sorted ( abins ) : 
                    return ROOT.TAxis ( nbins , array.array ( 'd' , abins ) )
                
        elif 2 == len ( bins ) and \
               isinstance ( bins [ 0 ] , num_types ) and \
               isinstance ( bins [ 1 ] , num_types ) and bins [ 0 ] < bins [ 1 ] :
            return ROOT.TAxis ( nbins , bins [ 0 ] , bins [ 1 ] )
        
        elif nbins + 1 <= len ( bins ) :                    
            abins = [ float ( e ) for e in bins [ : nbins + 1 ] ]
            if is_sorted ( abins ) : 
                return ROOT.TAxis ( nbins , array.array ( 'd' , abins ) ) 

    elif isinstance ( nbins , sequence_types ) and not bins :
        
        abins = [ float ( e ) for e in nbins ]
        if 2 <= len ( abins ) and is_sorted ( abins ) : 
            return ROOT.TAxis ( nbins , array.array ( 'd' , abins ) )
        
    if isinstance ( nbins , num_types ) and bins : 
        abins = [ float ( nbins ) ] + [ float ( e ) for e in bins ]
        if is_sorted ( abins ) : 
            return ROOT.TAxis ( nbins , array.array ( 'd' , abins ) ) 

    raise ArgumentError('make_axis: invalid arguments %s' % str ( ( nbins , ) + bins ) ) 

# =============================================================================
## create TH1D from model 
def h1_from_model ( model , *args ) :
    
    if isinstance ( model , ROOT.ROOT.RDF.TH1DModel ) :
        
        bins = model.fBinXEdges
        if bins : return ROOT.TH1D ( model.fName , model.fNbinsX , array.array  ('d' , model.fBinXEdges ) ) 
        else    : return ROOT.TH1D ( model.fName , model.fNbinsX , model.fXLow , model.fXUp               ) 

    if isinstance ( model , ROOT.TAxis ) :
        return h1_axis ( model , title = '1D histogram' , name = hID () , double = True )

    preps = []    
    aargs = [ model ] + list ( args )
    
    if 1 < len ( aargs ) and isistance ( aargs[0] , string_types ) :
        preps.append ( aargs[0] )
        aargs = aargs [1:]
    if 1 < len ( aargs ) and isistance ( aargs[0] , string_types ) :
        preps.append ( aargs[0] )
        aargs = aargs [1:]

    axis = None 
    if 1 <= len  ( aargs ) :
        try :
            axis = make_axis ( *aargs )
        except :
            pass
        
    if not axis : return None  ## RETURN
    
    if   2 == len ( preps ) : name , title = preps 
    elif 1 == len ( preps ) : name , title = hID () , preps [0]  
    else                    : name , title = hID () , '1D histo'
    
    return h1_axis ( axis , title = title , name = name , double = True )

        
# =============================================================================
## make 1D-histogram from axis
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def h1_axis ( axis           ,
              title  = '1D'  , 
              name   = None  ,
              double = True  ) :
    """ Make 1D-histogram with binning defined by already created axes    
    >>> axis = ...
    >>> h1 = h1_axes ( axis , title = 'MyHisto' )     
    """
    #
    if not name : name = hID()
    #
    if not issubclass ( type ( axis ) , ROOT.TAxis ) : axis = axis_bins   ( axis )
    #
    bins  = axis.edges()
    #
    if isinstance ( double , type ) and issubclass ( double , ROOT.TH1 ) : typ = double
    else : typ = ROOT.TH1D if double else ROOT.TH1F
    #
    h1  = typ ( name  ,
                title ,
                len ( bins ) - 1 , array.array ( 'd' , bins ) )
    ##
    h1.Sumw2()
    return h1

# =============================================================================
## make 2D-histogram from axes
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def h2_axes ( x_axis            ,
              y_axis            ,
              title  = '2D'     , 
              name   = None     ,
              double = True     ) :
    """ Make 2D-histogram with binning deifned by already created axes    
    >>> x_axis = ...
    >>> y_axis = ...
    >>> h2 = h2_axes ( x_axis , y_axis , title = 'MyHisto' )     
    """
    #
    if not name : name = hID() 
    #
    if not issubclass ( type ( x_axis ) , ROOT.TAxis ) : x_axis = axis_bins   ( x_axis )
    if not issubclass ( type ( y_axis ) , ROOT.TAxis ) : y_axis = axis_bins   ( y_axis )
    #
    # 
    x_bins  = x_axis.edges()
    y_bins  = y_axis.edges()
    #
    if isinstance ( double , type ) and issubclass ( double , ROOT.TH2 ) : typ = double
    else : typ = ROOT.TH2D if double else ROOT.TH2F
    #
    h2  =  typ ( name  ,
                 title ,
                 len ( x_bins ) - 1 , array.array ( 'd' , x_bins ) ,
                 len ( y_bins ) - 1 , array.array ( 'd' , y_bins ) )
    ##
    h2.Sumw2()
    return h2

# =============================================================================
## make 3D-histogram from axes
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-07-18
def h3_axes ( x_axis            ,
              y_axis            ,
              z_axis            ,
              title  = '3D'     , 
              name   = None     ,
              double = True     ) :
    """ Make 3D-histogram with binning deifned by already created axes    
    >>> x_axis = ...
    >>> y_axis = ...
    >>> z_axis = ...
    >>> h3 = h3_axes ( x_axis , y_axis , z_axis , title = 'MyHisto' )     
    """
    #
    if not name : name = hID() 
    #
    if not issubclass ( type ( x_axis ) , ROOT.TAxis ) : x_axis = axis_bins   ( x_axis )
    if not issubclass ( type ( y_axis ) , ROOT.TAxis ) : y_axis = axis_bins   ( y_axis )
    if not issubclass ( type ( z_axis ) , ROOT.TAxis ) : z_axis = axis_bins   ( z_axis )
    #
    x_bins  = x_axis.edges ()
    y_bins  = y_axis.edges ()
    z_bins  = z_axis.edges ()
    #
    if isinstance ( double , type ) and issubclass ( double , ROOT.TH3 ) : typ = double
    else : typ = ROOT.TH3D if double else ROOT.TH3F
    #
    return typ ( name  ,
                 title ,
                 len ( x_bins ) - 1 , array.array ( 'd' , x_bins ) ,
                 len ( y_bins ) - 1 , array.array ( 'd' , y_bins ) , 
                 len ( z_bins ) - 1 , array.array ( 'd' , z_bins ) ) 



# =======================================================================
## calculate the `difference' between two histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h_diff_ ( h1 , h2 , func = lambda s1,s2 : (s1/s2).value() ) :
    """ Estimate the `difference' between two histograms
    """
    se = SE()
    
    for bin in h1 :
        
        v1 = h1 [ bin ] 
        v2 = h2 [ bin ] 
        
        se += func ( v1 , v2 )
        
    return se 


ROOT.TH1F.histoDiff = _h_diff_
ROOT.TH1D.histoDiff = _h_diff_
ROOT.TH2F.histoDiff = _h_diff_
ROOT.TH2D.histoDiff = _h_diff_
ROOT.TH3F.histoDiff = _h_diff_
ROOT.TH3D.histoDiff = _h_diff_

# =============================================================================
## perform some accumulation for the histogram
#  @code
#  h =...
#  sum = h.accumulate ()
#  sum = h.accumulate ( cut = lambda s :  0.4<=s[1].value()<0.5 ) 
#  sum = h.accumulate ( low  = 1    , high = 14 ) ## accumulate over    1<= ibin <14
#  sum = h.accumulate ( xmin = 0.14 , xmax = 14 ) ## accumulate over xmin<= x    <xmax
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_accumulate_ ( h                        ,
                     func = lambda s,v : s + v ,
                     cut  = lambda s   : True  , 
                     init = VE ()              ,
                     low  =  1                 ,
                     high = -1                 ,
                     xmin = neg_infinity       ,
                     xmax = pos_infinity       ) : 
    """ Accumulate the function value over the histogram
    >>> h =...
    >>> sum = h.accumulate ()
    >>> sum = h.accumulate ( cut = lambda s :  0.4<=s[1].value()<0.5 ) 
    >>> sum = h.accumulate ( low  = 1    , high = 14 ) ## accumulate over    1<= ibin <14
    >>> sum = h.accumulate ( xmin = 0.14 , xmax = 14 ) ## accumulate over xmin<= x    <xmax
    """
    result = init
    ##
    if high <= 0 : high = h.bins() + 1
    xmin = max ( xmin , h.xmin() )
    xmax = min ( xmax , h.xmax() )
    ## check sums 
    if high <= low  or  xmax <= xmin : return result
    ##
    for i in h.items() :
        ibin = i[0] 
        if low <= ibin < high :
            xval = i[1].value()
            dx   = i[1].error()
            if inrange ( xmin , xval , xmax ) :
                if cut ( i ) : result = func ( result , i[-1] )
    return result


# =============================================================================
## perform some accumulation for the histogram
#  @code
#  h =...
#  sum = h.accumulate ()
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _hh_accumulate_ ( histo ) :    
    """ Accumulate the function value over the histogram
    >>> h =...
    >>> sum = h.accumulate ()
    """
    result = VE(0) 
    for i in  histo :
        result += histo [ i ]
    return result


# =============================================================================
## get the sum of entries 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_sum_ ( h           ,
               xmin = None ,
               xmax = None ) :
    """ Get the sum of histogram entries over the specified range xmin<x<xmax
    >>> h = ....
    >>> h.sum ( 1.0 , 20.0 )    
    """
    if xmin is None  or xmax is None :
        xmn , xmx = h.xminmax()
        if xmin is None : xmin = xmn
        if xmax is None : xmax = xmx
        
    return _h1_accumulate_ ( h , cut = lambda s : xmin<=s[1].value()<=xmax ) 

# =============================================================================
## simple scaling
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h_scale_ ( histo , val = 1.0 ) :
    """ Scale the histogram to certain integral
    >>> h = ...
    >>> h.scale ( 15 )    
    """
    factor  = 0.0
    val     = float(val)
    if 0   != val :

        total = VE()
        for ibin in histo :
            total += histo[ ibin ]
        total = total.value()
        
        if 0 != total :
            factor = val/total

    if not histo.GetSumw2() : histo.Sumw2()

    for ibin in histo :

        value          = histo [ ibin ]
        value         *= factor
        histo[ ibin ]  = value 

    return histo



# =============================================================================
## simple shift of the histogram
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_shift_ ( h , bias ) :
    """ Simple shift of the histogram :    
    >>> h = ... # the histogram
    >>> h2 = h.shift ( -5 * MeV )    
    """
    #
    if not h     .GetSumw2()  : h    .Sumw2()
    result = h.clone() ;
    result.Reset() ;
    if not result.GetSumw2()  : result.Sumw2()
    #
    for i,x,y in result.items() :
        
        x         += bias
        result[i]  = h ( x )
        
    result.ResetStats() 
    return result


# =============================================================================
## simple shift of the histogram
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_irshift_ ( h , ibias ) :
    """ Simple shift of the histogram :    
    >>> h   = ... # the histogram
    >>> h >>= 5   # shift for 5 bins left     
    """
    #
    ##
    if not isinstance ( ibias , integer_types ) : return NotImplemented 
    if not h     .GetSumw2()  : h    .Sumw2()
    ##
    if   0 == ibias : return h ## RETURN 
    elif ibias < 0  :
        h <<= abs(ibias)
        return h               ## RETURN
    ##
    for i in reversed ( h ) : 
        j  = i - ibias
        if j in h :  h [ i ] = h[ j ]
        else      :  h [ i ] = VE() 
        
    h.ResetStats() 
    return h     

# =============================================================================
## simple shift of the histogram
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_ilshift_ ( h , ibias ) :
    """ Simple shift of the histogram :    
    >>> h   = ...  # the histogram
    >>> h <<= 5    # shift for 5 bins left     
    """
    #
    ##
    if not isinstance ( ibias , integer_types ) : return NotImplemented 
    if not h     .GetSumw2()  : h    .Sumw2()
    ##
    if   0 == ibias : return h  ## RETURN 
    elif ibias < 0  :
        h >>= abs(ibias) 
        return h                ## RETURN
    ##
    for i in h :
        j = i + ibias
        if j in h :  h [ i ] = h[ j ]
        else      :  h [ i ] = VE() 
        
    h.ResetStats() 
    return h     

# =============================================================================
## simple shift of the histogram
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_lshift_ ( h , ibias ) :
    """ Simple shift of the histogram :    
    >>> h = ...      # the histogram
    >>> h2 = h << 5  # shift for 5 bins left     
    """
    #
    ##
    if not isinstance ( ibias , integer_types ) : return NotImplemented 
    ## 
    if not h     .GetSumw2()  : h    .Sumw2()
    result = h.Clone( hID() ) ;
    result.Reset() ;
    if not result.GetSumw2()  : result.Sumw2()
    #
    result <<= ibias

    return result 

# =============================================================================
## simple shift of the histogram
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_rshift_ ( h , ibias ) :
    """ Simple shift of the histogram :
    
    >>> h = ...      # the histogram
    >>> h2 = h >> 5  # shift for 5 bins right     
    """
    #
    if isinstance ( ibias , ROOT.TH1 ) : return _h1_rrshift_ ( ibias , h ) 
    ##
    if not isinstance ( ibias , integer_types ) : return NotImplemented 
    ##
    return _h1_lshift_ ( h , -1 * ibias )


# =============================================================================
## perform some integration(taking into acount the bin-width) for the histogram
#  @code
#  h =...
#  sum = h.integrate ()
#  sum = h.integrate ( cut = lambda s :  0.4<=s[1].value()<0.5 ) 
#  sum = h.integrate ( lowx = 1    , highx = 14  ) ## accumulate over    1<= ibin <14
#  sum = h.interrate ( xmin = 0.14 , xmax = 21.2 ) ## accumulate over xmin<= x    <xmax
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_integrate_ ( h                           ,
                     func   = lambda s,v : s + v ,
                     cut    = lambda s   : True  , 
                     init   = VE ()              ,
                     lowx   =  1                 ,
                     highx  = -1                 ,
                     xmin   = neg_infinity       ,
                     xmax   = pos_infinity       ) : 
    """ Perform some integration (taking into acount the bin-width) for the histogram
    >>> h =...
    >>> sum = h.integrate ()
    >>> sum = h.integrate ( cut  = lambda s :  0.4<=s[1].value()<0.5 ) 
    >>> sum = h.integrate ( lowx = 1    , highx = 14   ) ## integrate over    1<= xbin <14
    >>> sum = h.integrate ( xmin = 0.14 , xmax  = 21.1 ) ## integrate over xmin<= x    <xmax
    """
    result = init
    ##
    if highx <= 0 : highx = h.bins() + 1

    xmin = max ( xmin , h.xmin() )
    xmax = min ( xmax , h.xmax() )
    
    ## check
    if highx <= lowx  or xmax <= xmin : return result
    ##
    for i in h.items() :

        ibinx = i[0]

        if lowx <= ibinx < highx : 

            x    = i[1]
            y    = i[2]
            
            xval = x.value()
            dx   = x.error()
            
            yval = y.value()
            dy   = y.error()

            xl   = max ( xmin , xval - dx )
            xh   = min ( xmax , xval + dx )
            if xh <= xl : continue
            
            vol  = ( xh - xl )
            
            if cut ( i ) : result = func ( result , i[-1] * vol )
                
    return result 


# =============================================================================
## perform some integration(taking into acount the bin-width) for the histogram
#  @code
#  h =...
#  sum = h.integrate ()
#  sum = h.integrate ( cut = lambda s :  0.4<=s[1].value()<0.5 ) 
#  sum = h.integrate ( lowx = 1    , highx = 14  ) ## accumulate over    1<= ibin <14
#  sum = h.interrate ( xmin = 0.14 , xmax = 21.2 ) ## accumulate over xmin<= x    <xmax
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_integrate_ ( h                         ,
                     func   = lambda s,v : s + v ,
                     cut    = lambda s   : True  , 
                     init   = VE ()              ,
                     lowx   =  1                 ,
                     highx  = -1                 ,
                     lowy   =  1                 ,
                     highy  = -1                 ,
                     xmin   = neg_infinity            ,
                     xmax   = pos_infinity            ,
                     ymin   = neg_infinity            ,
                     ymax   = pos_infinity       ) : 
    """ Perform some integration (taking into acount the bin-width) for the histogram
    >>> h =...
    >>> sum = h.integrate ()
    >>> sum = h.integrate ( cut  = lambda s :  0.4<=s[1].value()<0.5 ) 
    >>> sum = h.integrate ( lowx = 1    , highx = 14   ) ## integrate over    1<= xbin <14
    >>> sum = h.integrate ( xmin = 0.14 , xmax  = 21.1 ) ## integrate over xmin<= x    <xmax
    """
    result = init
    ##
    if highx <= 0 : highx = h.binsx() + 1
    if highy <= 0 : highy = h.binsy() + 1

    xmin = max ( xmin , h.xmin() )
    xmax = min ( xmax , h.xmax() )
    ymin = max ( ymin , h.ymin() )
    ymax = min ( ymax , h.ymax() )
    
    ## check
    if highx <= lowx  or xmax <= xmin : return result
    if highy <= lowy  or ymax <= ymin : return result
    ##
    for i in h.iteritems() :

        ibinx = i[0]
        ibiny = i[1]
        
        if lowx <= ibinx < highx and lowy <= ibiny < highy : 

            x    = i[2]
            y    = i[3]
            
            xval = x.value()
            dx   = x.error()
            
            yval = y.value()
            dy   = y.error()
            
            xl   = max ( xmin , xval - dx )
            xh   = min ( xmax , xval + dx )
            if xh <= xl : continue

            vol  = xh - xl
            
            yl   = max ( ymin , yval - dy )
            yh   = min ( ymax , yval + dy )
            if yh <= yl : continue

            vol *= yh - yl
            
            if cut ( i ) : result = func ( result , i[-1] * vol )
                
    return result 


# =============================================================================
## perform some integration(taking into acount the bin-width) for the histogram
#  @code
#  h =...
#  sum = h.integrate ()
#  sum = h.integrate ( cut = lambda s :  0.4<=s[1].value()<0.5 ) 
#  sum = h.integrate ( lowx = 1    , highx = 14  ) ## accumulate over    1<= ibin <14
#  sum = h.interrate ( xmin = 0.14 , xmax = 21.2 ) ## accumulate over xmin<= x    <xmax
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_integrate_ ( h                         ,
                     func   = lambda s,v : s + v ,
                     cut    = lambda s   : True  , 
                     init   = VE ()              ,
                     lowx   =  1                 ,
                     highx  = -1                 ,
                     lowy   =  1                 ,
                     highy  = -1                 ,
                     lowz   =  1                 ,
                     highz  = -1                 ,
                     xmin   = neg_infinity       ,
                     xmax   = pos_infinity       ,
                     ymin   = neg_infinity       ,
                     ymax   = pos_infinity       , 
                     zmin   = neg_infinity       ,
                     zmax   = pos_infinity       ) : 
    """ Perform some integration (taking into acount the bin-width) for the histogram
    >>> h =...
    >>> sum = h.integrate ()
    >>> sum = h.integrate ( cut  = lambda s :  0.4<=s[1].value()<0.5 ) 
    >>> sum = h.integrate ( lowx = 1    , highx = 14   ) ## integrate over    1<= xbin <14
    >>> sum = h.integrate ( xmin = 0.14 , xmax  = 21.1 ) ## integrate over xmin<= x    <xmax
    """
    result = init
    ##
    if highx <= 0 : highx = h.binsx() + 1
    if highy <= 0 : highy = h.binsy() + 1
    if highz <= 0 : highz = h.binsz() + 1

    xmin = max ( xmin , h.xmin() )
    xmax = min ( xmax , h.xmax() )
    ymin = max ( ymin , h.ymin() )
    ymax = min ( ymax , h.ymax() )
    zmin = max ( zmin , h.zmin() )
    zmax = min ( zmax , h.zmax() )
    
    ## check
    if highx <= lowx  or xmax <= xmin : return result
    if highy <= lowy  or ymax <= ymin : return result
    if highz <= lowz  or zmax <= zmin : return result
    ##
    for i in h.items() :

        ibinx = i [ 0 ]
        ibiny = i [ 1 ]
        ibinz = i [ 2 ]
        
        if lowx <= ibinx < highx and lowy <= ibiny < highy and lowz <= ibinz < highz :
        
            x    = i [ 3 ]
            y    = i [ 4 ]
            z    = i [ 5 ]
            
            xval = x.value()
            dx   = x.error()
            
            yval = y.value()
            dy   = y.error()

            zval = z.value()
            dz   = z.error()
            
            xl   = max ( xmin , xval - dx )
            xh   = min ( xmax , xval + dx )
            if xh <= xl : continue

            vol  = xh - xl
            
            yl   = max ( ymin , yval - dy )
            yh   = min ( ymax , yval + dy )
            if yh <= yl : continue

            vol *= yh - yl
            
            zl   = max ( zmin , zval - dz )
            zh   = min ( zmax , zval + dz )
            if zh <= zl : continue
            
            vol *= zh - zl

            if cut ( i ) : result = func ( result , i[-1] * vol )
                
    return result 


# =============================================================================
## make transformation of histogram content 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-10-23  
def _h3_transform_ ( h3 , func ) :
    """ Make the transformation of the 2D-histogram content 
    >>> func = lambda x,y,z,v: v   ## identical transformation/copy
    >>> h3   = ...
    >>> h3   = h3.fransform ( func ) 
    """
    #
    if not h3.GetSumw2() : h3.Sumw2()
    hn = h3.Clone ( hID() )
    if not hn.GetSumw2() : hn.Sumw2()
    #
    for ix,iy,iz, x,y,z, v  in h3.items() :        
        hn [ ix , iy , iz ] = func ( x, y , z , v ) 
        
    hn.ResetStats() 
    return hn 

ROOT.TH3F. transform = _h3_transform_ 
ROOT.TH3D. transform = _h3_transform_ 

# =============================================================================
## rescale the histogram for effective uniform bins
#  new_content = old_content * factor / bin_area
#  @code
#  >>> h2 = ...
#  >>> h3 = h2.rescale_bins ( h2 , 1 )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-03-19   
def _h3_rescale_ ( h3 , factor = 1 ) :
    """ Rescale the histogram for effective uniform bins : 
    new_bin_content = old_bin_content * factor / bin_area
    >>> h1 = ...
    >>> h2 = h1.rescale_bins ( h1 , 1 )
    
    """
    return _h3_transform_ ( h3 , lambda x , y , z, v : ( 0.125 * factor / ( x.error() * y.error() * z.error () ) * v ) )

ROOT.TH3F. rescale_bins = _h3_rescale_ 
ROOT.TH3D. rescale_bins = _h3_rescale_ 



ROOT.TH1F .   shift     = _h1_shift_
ROOT.TH1D .   shift     = _h1_shift_
ROOT.TH1D . __rshift__  = _h1_rshift_
ROOT.TH1F . __rshift__  = _h1_rshift_
ROOT.TH1D . __lshift__  = _h1_lshift_
ROOT.TH1F . __lshift__  = _h1_lshift_
ROOT.TH1D . __ilshift__ = _h1_ilshift_
ROOT.TH1F . __ilshift__ = _h1_ilshift_
ROOT.TH1D . __irshift__ = _h1_irshift_
ROOT.TH1F . __irshift__ = _h1_irshift_


# =============================================================================    
for t in ( ROOT.TH1F , ROOT.TH1D ) :    
    t . accumulate  = _h1_accumulate_ 
    t . sum         = _h1_sum_ 
    t . integrate   = _h1_integrate_ 

for t in ( ROOT.TH2F , ROOT.TH2D ) :    
    t . accumulate = _hh_accumulate_ 
    t . sum        = _hh_accumulate_
    t . integrate  = _h2_integrate_ 

for t in ( ROOT.TH3F , ROOT.TH3D ) :    
    t . accumulate = _hh_accumulate_ 
    t . sum        = _hh_accumulate_
    t . integrate  = _h3_integrate_ 

# generic sscaling
ROOT.TH1 . scale      = _h_scale_

HStats   = Ostap.Utils.HistoStat
# =============================================================================
## calculate bin-by-bin momenta 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_moment_ ( h1 , order ) :
    """ Get `bin-by-bin'-moment around the specified value    
    >>> histo = ...
    >>> mom   = histo.moment ( 4 , 0 ) 
    """
    #
    h1.ResetStats()
    #
    m = HStats.moment    ( h1 , order )
    e = HStats.momentErr ( h1 , order )
    #
    return VE ( m , e*e )  if 0<= e  else VE ( m , -e*e )

_h1_moment_ .__doc__ += '\n' + HStats.moment    .__doc__
_h1_moment_ .__doc__ += '\n' + HStats.momentErr .__doc__

# =============================================================================
## calculate bin-by-bin central momenta 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_central_moment_ ( h1 , order ) :
    """ Get ``bin-by-bin'' central moment    
    >>> histo = ...
    >>> cmom  = histo.centralMoment ( 4 ) 
    """
    #
    h1.ResetStats()
    #
    m = HStats.centralMoment    ( h1 , order )
    e = HStats.centralMomentErr ( h1 , order )
    #
    return VE ( m , e*e )  if 0 <= e else VE ( m , -e*e ) 

_h1_central_moment_ .__doc__ += '\n' + HStats.centralMoment    .__doc__
_h1_central_moment_ .__doc__ += '\n' + HStats.centralMomentErr .__doc__

# =============================================================================
## calculate bin-by-bin ``standardized moment''
#  \f$ s_n = \frac{\mu_n}{\sigma^n} = \frac{\int (x-\mu)^n h(x) dx }{\sigma^n} \f$
#  @see https://en.wikipedia.org/wiki/Standardized_moment
#  for <code>exp_moment=True</code> it calculates
#  \f$ \sign \mu_n \\frac{ \left| \mu_n \right|^{1/n} } {\sigma }\f$  
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2019-09-04
def _h1_std_moment_ ( h1 , order ,  exp_moment = False ) :
    """ Get ``bin-by-bin'' ``standardized moment''
    >>> histo = ...
    >>> cmom  = histo.stdMoment ( 4 ) 
     - see https://en.wikipedia.org/wiki/Standardized_moment
     For `exp_moment==True`, it calcualtes the value
      sign(cmom)*(abs(cmom)**1.0/order)
     
    """
    #

    if   1 == order : return VE ( 0 , 0 ) 
    elif 2 == order : return VE ( 1 , 1 )
    #
    mom   = h1.centralMoment ( order )    
    sigma = math.sqrt ( h1.rms().value() )
    
    if exp_moment :
        
        mm = abs ( mom ) ** ( 1. / order )
        if mom.value() < 0 : mm *= -1 
        moment = mm / sigma

    else :

        moment = mom / ( sigma ** order ) 
        

    return moment 
    
# =============================================================================
## get skewness
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_skewness_ ( h1 ) :
    """ Get the skewness
    >>> histo = ...
    >>> skew  = histo.skewness () 
    """
    #
    h1.ResetStats()
    #
    m = HStats.skewness    ( h1 )
    e = HStats.skewnessErr ( h1 )
    #
    return VE ( m , e*e )  if 0<= e  else VE ( m , -e*e ) 

_h1_skewness_ .__doc__ += '\n' + HStats.skewness    .__doc__
_h1_skewness_ .__doc__ += '\n' + HStats.skewnessErr .__doc__

# =============================================================================
## get kurtosis
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_kurtosis_ ( h1 ) :
    """ Get the kurtosis
    >>> histo = ...
    >>> k     = histo.kurtosis () 
    """
    #
    h1.ResetStats()
    #
    m = HStats.kurtosis    ( h1 )
    e = HStats.kurtosisErr ( h1 )
    #
    return VE ( m , e*e )  if 0<= e  else VE ( m , -e*e ) 

_h1_kurtosis_ .__doc__ += '\n' + HStats.kurtosis    .__doc__
_h1_kurtosis_ .__doc__ += '\n' + HStats.kurtosisErr .__doc__

# =============================================================================
## get mean
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_mean_ ( h1 ) :
    """ Get the mean
    >>> histo = ...
    >>> k     = histo.mean () 
    """
    #
    h1.ResetStats()
    #
    m = HStats.mean    ( h1 )
    e = HStats.meanErr ( h1 )
    #
    return VE ( m , e*e )  if 0<= e  else VE ( m , -e*e ) 

_h1_mean_ .__doc__ += '\n' + HStats.mean    .__doc__
_h1_mean_ .__doc__ += '\n' + HStats.meanErr .__doc__

# =============================================================================
## get RMS
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_rms_ ( h1 ) :
    """ Get the rms    
    >>> histo = ...
    >>> s     = histo.rms () 
    """
    #
    h1.ResetStats()
    #
    m = HStats.rms    ( h1 )
    e = HStats.rmsErr ( h1 )
    #
    return VE ( m , e*e )  if 0<= e  else VE ( m , -e*e ) 

_h1_rms_ .__doc__ += '\n' + HStats.rms    .__doc__
_h1_rms_ .__doc__ += '\n' + HStats.rmsErr .__doc__

for h in ( ROOT.TH1F , ROOT.TH1D ) :
    
    h.mean           = _h1_mean_
    h.rms            = _h1_rms_
    h.skewness       = _h1_skewness_
    h.kurtosis       = _h1_kurtosis_
    h.moment         = _h1_moment_
    h.centralMoment  = _h1_central_moment_
    h.central_moment = _h1_central_moment_
    h.stdMoment      = _h1_std_moment_
    h.std_moment     = _h1_std_moment_
    #

ROOT.TH1.nEff = ROOT.TH1.GetEffectiveEntries 


# =============================================================================
## Get the histogram moments
#  @code
#  histo = ....
#  m     = histo.the_moment ( 10 ) 
#  @edncode 
def _h1_the_moment_ ( h1 , N ) :
    """ Get the histogram moments
    >>> histo = ....
    >>> m     = histo.the_moment ( 10 ) 
    """
    assert isinstance ( N , int ) and 0 <= N , 'Invalid order!'

    M = Ostap.Math.WMoment_(N)
    m = M()

    for i , x , y in h1.items() :
        m.add ( x.value() , y.value() )
        
    return m 
                
    
ROOT.TH1F.the_moment = _h1_the_moment_
ROOT.TH1D.the_moment = _h1_the_moment_

# =============================================================================
## calculate bin-by-bin moments
#  \f$ m(k_x,k_y; x , y  ) \equiv 
#   \frac{ \sum_i (x_i - x)^{n_x}
#                 (y_i - y)^{n_y} N_i}
#        { \sum_j N_j } \f$ 
#  @code
#  h2 = ...
#  m  = h2.moment ( 2  , 3 , 0.0 , 0.0 )
#  @endcode
#  @see Ostap::Utils::HistoStat::moment2 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_moment_ ( h2 , orderx , ordery  , x0 = 0.0 , y0 = 0.0 ) :
    """ Get ``bin-by-bin''-moment around the specified value    
    >>> h2   = ...
    >>> mom  = histo.moment  ( 2 , 3 , 0.0 , 0.0 ) 
    """
    #
    h2.ResetStats()
    #
    return HStats.moment2 ( h2 , orderx , ordery , x0 , y0 )

_h2_moment_ . __doc__ +=  '\n' + HStats.moment2.__doc__  

# =============================================================================
## calculate bin-by-bin central moments
#  \f$ m(k_x,k_y; x , y  ) \equiv 
#   \frac{ \sum_i (x_i - \mu_x)^{k_x}
#                 (y_i - \mu_y)^{k_y} N_i }
#        { \sum_j N_j } \f$ 
#  @code
#  h2 = ...
#  m  = h2.central_moment ( 2  , 3 )
#  @endcode
#  @see Ostap::Utils::HistoStat::cental_moment2 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_cmoment_ ( h2 , orderx , ordery ) :
    """ Get ``bin-by-bin''-central moment 
    >>> h2   = ...
    >>> mom  = histo.central_moment  ( 2 , 3 ) 
    """
    #
    h2.ResetStats()
    #
    return HStats.central_moment2 ( h2 , orderx , ordery )

_h2_cmoment_ . __doc__ +=  '\n' + HStats.central_moment2.__doc__  

# =============================================================================
## calculate bin-by-bin standartized moment
#  \f$ m(k_x,k_y; x , y  ) \equiv
#   \frac{1}{\sigma_x^{k_x}\sigma_y^{k_y}}
#   \frac{ \sum_i (x_i - \mu_x)^{k_x}
#                 (y_i - \mu_y)^{k_y} N_i }
#        { \sum_j N_j } \f$ 
#  @code
#  h2 = ...
#  m  = h2.central_moment ( 2  , 3 )
#  @endcode
#  @see Ostap::Utils::HistoStat::cental_moment2 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_smoment_ ( h2 , orderx , ordery ) :
    """ Get ``bin-by-bin''-standartized moment 
    >>> h2   = ...
    >>> mom  = histo.std_moment  ( 2 , 3 ) 
    """
    #
    h2.ResetStats()
    #
    return HStats.std_moment2 ( h2 , orderx , ordery )

_h2_smoment_ . __doc__ +=  '\n' + HStats.std_moment2.__doc__  


for h in ( ROOT.TH2F , ROOT.TH2D ) :
    h.moment         = _h2_moment_
    h.central_moment = _h2_cmoment_
    h.std_moment     = _h2_smoment_

# =============================================================================
## calculate bin-by-bin moment
#  \f$ m(k_x,k_y,k_z; x , y  , z ) \equiv 
#   \frac{ \sum_i (x_i - x)^{k_x}
#                 (y_i - y)^{k_y}
#                 (z_i - z )^{k_z} N_i}
#        { \sum_j N_j } \f$ 
#  @code
#  h3 = ...
#  m  = h3.moment ( 2  , 3 , 4 , 0.0 , 0.0 , 0.0 )
#  @endcode#  @see Ostap::Utils::HistoStat::moment3 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_moment_ ( h3 , orderx , ordery , orderz  , x0 = 0.0 , y0 = 0.0 , z0 = 0.0 ) :
    """ Get ``bin-by-bin''-moment around the specified value    
    >>> h3   = ...
    >>> mom  = histo.moment  ( 2 , 3 , 4 , 0.0 , 0.0 , 0.0 ) 
    """
    #
    h2.ResetStats()
    #
    return HStats.moment3 ( h3 , orderx , ordery , orderz ,  x0 , y0 , z0 )

_h3_moment_ . __doc__ +=  '\n' + HStats.moment3.__doc__ 


# =============================================================================
## calculate bin-by-bin central moment
#  \f$ m(k_x,k_y,k_z; x , y  , z ) \equiv 
#   \frac{ \sum_i (x_i - \mu_x)^{k_x}
#                 (y_i - \mu_y)^{k_y}
#                 (z_i - \mu_z)^{k_z} N_i}
#        { \sum_j N_j } \f$ 
#  @code
#  h3 = ...
#  m  = h3.moment ( 2  , 3 , 4 , 0.0 , 0.0 , 0.0 )
#  @endcode
#  @see Ostap::Utils::HistoStat::moment3 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_cmoment_ ( h3 , orderx , ordery , orderz  ) :
    """ Get ``bin-by-bin''-central moment around the specified value    
    >>> h3   = ...
    >>> mom  = histo.central_moment  ( 2 , 3 , 4 ) 
    """
    #
    h3.ResetStats()
    #
    return HStats.central_moment3 ( h3 , orderx , ordery , orderz  )

_h3_cmoment_ . __doc__ +=  '\n' + HStats.central_moment3.__doc__ 

# =============================================================================
## calculate bin-by-bin standartized  moment
#  \f$ m(k_x,k_y,k_z; x , y  , z ) \equiv 
#   \frac{1}{\sigma_x^{k_x}\sigma_y^{k_y}\sigma_z^{k_z}}
#   \frac{ \sum_i (x_i - \mu_x)^{k_x}
#                 (y_i - \mu_y)^{k_y}
#                 (z_i - \mu_z)^{k_z} N_i}
#        { \sum_j N_j } \f$ 
#  @code
#  h3 = ...
#  m  = h3.moment ( 2 , 3 , 4 , 0.0 , 0.0 , 0.0 )
#  @endcode
#  @see Ostap::Utils::HistoStat::std_moment3 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h3_smoment_ ( h3 , orderx , ordery , orderz  ) :
    """ Get ``bin-by-bin''-standartized moment
    >>> h3   = ...
    >>> mom  = histo.std_moment  ( 2 , 3 , 4 ) 
    """
    #
    h3.ResetStats()
    #
    return HStats.std_moment3 ( h3 , orderx , ordery , orderz  )


_h3_smoment_ . __doc__ +=  '\n' + HStats.std_moment3.__doc__ 

for h in ( ROOT.TH3F , ROOT.TH3D ) :
    h.moment         = _h3_moment_
    h.central_moment = _h3_cmoment_
    h.std_moment     = _h3_smoment_


# =============================================================================
## get some statistic information on the histogram content
#  @code 
#  >>> histo = ... 
#  >>> stat  = histo.stat()
#  >>> print stat
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-03-26
def _h_stat_ ( h ) :
    """ Get some statistic information on the histogram content
    >>> histo = ... 
    >>> stat  = histo.stat()
    >>> print ( stat )
    """
    cnt = SE() 
    for i in h : cnt += h[i].value()  
    return cnt

ROOT.TH1.stat = _h_stat_


# =============================================================================
## get some (weighted) statistic information on the histogram content
#  @code 
#  >>> histo = ... 
#  >>> wstat = histo.wstat()
#  >>> print ( wstat ) 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-03-26
def _h_wstat_ ( h ) :
    """ Get some weighted statistic infomration on the histogram content
    >>> histo = ... 
    >>> wstat  = histo.wstat()
    >>> print wstat 
    """
    cnt = WSE() 
    for i in h :
        v = h[i]
        cnt.add ( v.value() , 1.0/v.cov2() ) 
    return cnt

ROOT.TH1.wstat = _h_wstat_

# =============================================================================
## get some (weighted) statistic information on the histogram content
#  @code 
#  >>> histo = ... 
#  >>> wstat = histo.wstat()
#  >>> print wstat
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-03-26
def _h_wstat_ ( h ) :
    """ Get some weighted statistic information on the histogram content
    >>> histo = ... 
    >>> wstat  = histo.wstat()
    >>> print wstat 
    """
    cnt = WSE() 
    for i in h :
        v = h[i]
        cnt.add ( v.value() , 1.0/v.cov2() ) 
    return cnt

ROOT.TH1.wstat = _h_wstat_

# =============================================================================
## get some (weighted) statistic information on the histogram content
#  @code 
#  >>> histo = ... 
#  >>> wstat = histo.wstat()
#  >>> print ( wstat)
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-03-26
def _h_wstat_ ( h ) :
    """ Get some weighted statistic information on the histogram content
    >>> histo = ... 
    >>> wstat  = histo.wstat()
    >>> print wstat 
    """
    cnt = WSE() 
    for i in h :
        v = h[i]
        cnt.add ( v.value() , 1.0/v.cov2() ) 
    return cnt

ROOT.TH1.wstat = _h_wstat_

# =============================================================================
## get some (weighted) statistic information on the histogram content
#  @code 
#  >>> histo = ... 
#  >>> xstat = histo.xstat()
#  >>> print ( xstat ) 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-03-26
def _h_xstat_ ( h ) :
    """ Get some weighted statistic information on the histogram content
    >>> histo = ... 
    >>> xstat  = histo.xstat()
    >>> print ( wstat ) 
    """
    cnt = WSE() 
    for i,x,y in h.iteritems() : 
        cnt.add ( x.value() , y.value() )
    return cnt

ROOT.TH1F.xstat = _h_xstat_
ROOT.TH1D.xstat = _h_xstat_


# =============================================================================
## get some (weighted) statistic information on the histogram X-content
#  @code 
#  >>> histo2 = ... 
#  >>> xstat = histo2.xstat()
#  >>> print xstat
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-09-20
def _h2_xstat_ ( h2 ) :
    """ Get some (weighted) statistic information on the histogram X-content
    >>> histo2 = ... 
    >>> xstat = histo2.xstat()
    >>> print ( xstat ) 
    """
    _h2 =  h2.ProjX()
    _r  = _h2.xstat() 
    del _h2
    return _r

# =============================================================================
## get some (weighted) statistic information on the histogram X-content
#  @code 
#  >>> histo2 = ... 
#  >>> ystat = histo2.ystat()
#  >>> print ( ystat ) 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-09-20
def _h2_ystat_ ( h2 ) :
    """ Get some (weighted) statistic information on the histogram X-content
    >>> histo2 = ... 
    >>> ystat = histo2.ystat()
    >>> print ( xstat ) 
    """
    _h2 =  h2.ProjY()
    _r  = _h2.xstat() 
    del _h2
    return _r

ROOT.TH2F.xstat = _h2_xstat_
ROOT.TH2D.xstat = _h2_xstat_
ROOT.TH2F.ystat = _h2_ystat_
ROOT.TH2D.ystat = _h2_ystat_

# =============================================================================
## get some (weighted) statistic information on the histogram X-content
#  @code 
#  >>> histo3 = ... 
#  >>> zstat = histo3.zstat()
#  >>> print ( zstat ) 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-09-20
def _h3_zstat_ ( h3 ) :
    """ Get some (weighted) statistic information on the histogram X-content
    >>> histo3 = ... 
    >>> zstat = histo3.zstat()
    >>> print ( zstat ) 
    """
    _h3 =  h3.ProjZ()
    _r  = _h3.xstat() 
    del _h3
    return _r

ROOT.TH3.xstat = _h2_xstat_  ## sic!
ROOT.TH3.ystat = _h2_ystat_  ## sic! 
ROOT.TH3.zstat = _h3_zstat_

# =============================================================================
## adjust the "efficiency"
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def ve_adjust ( ve , mn = 0 , mx = 1.0 ) :
    """ Adjust the efficiency
    """
    if ve.value() < mn : ve.setValue ( mn )
    if ve.value() > mx : ve.setValue ( mx )
    #
    return ve

# =============================================================================
## draw the line for the histogram 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-01-21 
def _level_ ( self , level = 0 , **kwargs ) :
    """ Draw ``NULL''-line for the histogram
    >>> h.level ( 5 )    
    >>> h.hline ( 5 ) ## ditto
    """
    mn,mx = self.xminmax()
    line = ROOT.TLine ( mn , level , mx , level )
    line.draw ( **kwargs )
    if not hasattr ( self , '_lines_' ) : self._lines_ = [] 
    self._lines_.append ( line ) 
    return line

# =============================================================================
## draw the vertical line for the histogram 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-01-21 
def _vline_ ( self , x , ymin = 0 , ymax = None , **kwargs ) :
    """ Draw vertical line for the histogram
    >>> h.vline ( 5 )    
    """
    if ymax is None : ymax = self.GetMaximum() 
    mn , mx = ymin , ymax 
    line = ROOT.TLine ( x , mn , x , mx )
    line.draw( **kwargs ) 
    if not hasattr ( self , '_lines_' ) : self._lines_ = [] 
    self._lines_.append ( line ) 
    return line 

# =============================================================================
## draw null-level for histogram  
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-01-21 
def _null_ ( self , linestyle = 2 ) :
    """ Draw NULL-line for the histogram    
    >>> h.null() 
    """
    return _level_ ( self , 0 , linestyle ) 

ROOT.TH1D. level = _level_
ROOT.TH1F. level = _level_
ROOT.TH1D. null  = _null_
ROOT.TH1F. null  = _null_

ROOT.TH1F. hline = _level_
ROOT.TH1F. vline = _vline_
ROOT.TH1D. hline = _level_
ROOT.TH1D. vline = _vline_


# =============================================================================
## add "fake" bin into the histogram
#  It is useful to control the functional behaviour at edge bins, e.g. f(0)=0...
#  @code
#  histo1 = ...
#  histo2 = histo1.add_fake_bin ( left = True )  
#  @endcode
#  Fake bin can extend the histogram range (for width>0) and
#  can be in the histogram range (for width<0) 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2016-08-15
def _h1_add_fake_bin_ ( h1 , left = True , value = 0 , width = 1.e-6 ) :
    """ Add ``fake'' bin into the histogram
    It is useful to control the functional behaviour at the edges
    >>> histo1 = ...
    >>> histo2 = histo1.add_fake_bin ( left = True )  
    Fake bin can extend the histogram range (for width>0) or
    it can be in the histogram range (for width<0) 
    """
    _a = h1.GetXaxis()
    _e = list ( _a.edges() ) 
    if left :
        _en = _e[ 0] - width * ( _e[ 1] - _e[ 0] )
        if 0 <= width : _e.insert ( 0 , _en )
        else          : _e.insert ( 1 , _en ) 
    else   :
        _en = _e[-1] + width * ( _e[-1] - _e[-2] )
        if 0 <= width : _e.append (      _en ) 
        else          : _e.insert ( -1 , _en )

    ## create new histogram 
    hn = h1_axis ( _e , title = h1.GetTitle() , double = type(h1) )
    
    if left : hn[ 1] = value
    else    : hn[-1] = value
    
    for i in h1:
        j = i + 1 if left else i
        hn [ j ] = h1 [ i ]
        
    #
    hn.ResetStats()
    #
    return hn 

# =============================================================================
## add "fake" bin into the histogram
#  It is useful to control the functional behaviour at edge bins, e.g. f(0)=0...
#  @code
#  histo1 = ...
#  histo2 = histo1.add_fake_bin_left() 
#  @endcode
#  Fake bin can extend the histogram range (for width>0) or 
#  it can be in the histogram range (for width<0) 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2016-08-15
def _h1_add_fake_bin_left_  ( h1 , value = 0 , width = 1.e-6 ) :
    """ Add ``fake'' bin into the histogram
    It is useful to control the functional behaviour at the edges,
    e.g. ensure f(0)=0
    >>> histo1 = ...
    >>> histo2 = histo1.add_fake_bin_left()   
    Fake bin can extend the histogram range (for width>0) or 
    it can be in the histogram range (for width<0) 
    """    
    return _h1_add_fake_bin_ ( h1 , True , value , width )

# =============================================================================
## add "fake" bin into the histogram
#  It is useful to control the functional behaviour at edge bins, e.g. f(1)=0...
#  @code
#  histo1 = ...
#  histo2 = histo1.add_fake_bin_right() 
#  @endcode
#  Fake bin can extend the histogram range (for width>0) or 
#  it can be in the histogram range (for width<0) 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2016-08-15
def _h1_add_fake_bin_right_ ( h1 , value = 0 , width = 1.e-6 ) :
    """ Add ``fake'' bin into the histogram
    It is useful to control the functional behaviour at the edges
    e.g. ensure f(1)=0
    >>> histo1 = ...
    >>> histo2 = histo1.add_fake_bin_right()   
    Fake bin can extend the histogram range (for width>0) or 
    it can be in the histogram range (for width<0) 
    """    
    return _h1_add_fake_bin_ ( h1 , False , value , width ) 


for _h in ( ROOT.TH1F , ROOT.TH1D ) :
    _h.add_fake_bin       = _h1_add_fake_bin_
    _h.add_fake_bin_left  = _h1_add_fake_bin_left_
    _h.add_fake_bin_right = _h1_add_fake_bin_right_


# =============================================================================
## add "fake" side(row of bins) into the histogram
#  It is useful to control the functional behaviour at edge bins, e.g. f(0)=0...
#  @code
#  histo1 = ...
#  histo2 = histo1.add_fake_side ( side = 1 )  
#  @endcode
#  Sides:
#  - 0: min x
#  - 1: max x 
#  - 2: min y
#  - 3: max y 
#  Fake bin can extend the histogram range (for width>0) and
#  can be in the histogram range (for width<0) 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2016-08-15
def _h2_add_fake_side_ ( h2 , side , value = 0 , width = 1.e-6 ) :
    """ Add ``fake'' side (row of bins) into the histogram
    It is useful to control the functional behaviour at the edges
    >>> histo1 = ...
    >>> histo2 = histo1.add_fake_side ( side = 1 )  
    Sides:
    - 0: min x
    - 1: max x 
    - 2: min y
    - 3: max y 
    Fake bin can extend the histogram range (for width>0) or
    it can be in the histogram range (for width<0) 
    """
    
    _ax = h2.GetXaxis()
    _ay = h2.GetYaxis()
    
    _ex = list ( _ax.edges() )
    _ey = list ( _ay.edges() )

    
    if   0 == side % 4  :  
        _en = _ex[ 0] - width * ( _ex[ 1] - _ex[ 0] )
        if 0 <= width : _ex.insert ( 0 , _en )
        else          : _ex.insert ( 1 , _en )
    elif 1 == side % 4  :  
        _en = _ex[-1] + width * ( _ex[-1] - _ex[-2] )
        if 0 <= width : _ex.append (      _en ) 
        else          : _ex.insert ( -1 , _en )
    elif 2 == side % 4  :  
        _en = _ey[ 0] - width * ( _ey[ 1] - _ey[ 0] )
        if 0 <= width : _ey.insert ( 0 , _en )
        else          : _ey.insert ( 1 , _en )
    elif 3 == side % 4  :  
        _en = _ey[-1] + width * ( _ey[-1] - _ey[-2] )
        if 0 <= width : _ey.append (      _en ) 
        else          : _ey.insert ( -1 , _en )

    ## create new histogram 
    hn  = h2_axes ( _ex , _ey , title = h2.GetTitle() , double = type ( h2 ) )

    _nx = hn.nbinsx ()
    _ny = hn.nbinsy ()
    
    if   0 == side % 4 :
        for j in range(1,_ny+1) : hn [  1 ,  j ] = value 
    elif 1 == side % 4 :
        for j in range(1,_ny+1) : hn [ -1 ,  j ] = value
    elif 2 == side % 4 :
        for i in range(1,_nx+1) : hn [  i ,  1 ] = value 
    elif 3 == side % 4 :
        for i in range(1,_nx+1) : hn [  i , -1 ] = value
        
    for i in h2:

        if   0 == side % 4 : j = i[0]+1,i[1]
        elif 1 == side % 4 : j = i
        elif 2 == side % 4 : j = i[0]  ,i[1]+1
        elif 3 == side % 4 : j = i
        
        hn [ j ] = h2 [ i ]
        
    hn.ResetStats()
    return hn 


for _h in ( ROOT.TH2F , ROOT.TH2D ) :
    _h.add_fake_side       = _h2_add_fake_side_

# =============================================================================
## add "fake" side(row of bins) into the histogram
#  It is useful to control the functional behaviour at edge bins, e.g. f(0)=0...
#  @code
#  histo1 = ...
#  histo2 = histo1.add_fake_side ( side = 1 )  
#  @endcode
#  Sides:
#  - 0: min x
#  - 1: max x 
#  - 2: min y
#  - 3: max y 
#  - 4: min z
#  - 5: max z 
#  Fake bin can extend the histogram range (for width>0) and
#  can be in the histogram range (for width<0) 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2016-08-15
def _h3_add_fake_side_ ( h2 , side , value = 0 , width = 1.e-6 ) :
    """ Add ``fake'' side (row of bins) into the histogram
    It is useful to control the functional behaviour at the edges
    >>> histo1 = ...
    >>> histo2 = histo1.add_fake_side ( side = 1 )  
    - 0: min x
    - 1: max x 
    - 2: min y
    - 3: max y 
    - 4: min z
    - 5: max z 
    Fake bin can extend the histogram range (for width>0) or
    it can be in the histogram range (for width<0) 
    """
    
    _ax = h3.GetXaxis()
    _ay = h3.GetYaxis()
    _az = h3.GetZaxis()
    
    _ex = list ( _ax.edges() )
    _ey = list ( _ay.edges() )
    _ez = list ( _ay.edges() )
    
    
    if   0 == side % 6  :   ## left  
        _en = _ex[ 0] - width * ( _ex[ 1] - _ex[ 0] )
        if 0 <= width : _ex.insert ( 0 , _en )
        else          : _ex.insert ( 1 , _en )
    elif 1 == side % 6  :   ## right 
        _en = _ex[-1] + width * ( _ex[-1] - _ex[-2] )
        if 0 <= width : _ex.append (      _en ) 
        else          : _ex.insert ( -1 , _en )
    elif 2 == side % 6  :   ## bottom
        _en = _ey[ 0] - width * ( _ey[ 1] - _ey[ 0] )
        if 0 <= width : _ey.insert ( 0 , _en )
        else          : _ey.insert ( 1 , _en )
    elif 3 == side % 6  :   ## top 
        _en = _ey[-1] + width * ( _ey[-1] - _ey[-2] )
        if 0 <= width : _ey.append (      _en ) 
        else          : _ey.insert ( -1 , _en )
    elif 4 == side % 6  :   ## right 
        _en = _ez[ 0] - width * ( _ez[ 1] - _ez[ 0] )
        if 0 <= width : _ez.insert ( 0 , _en )
        else          : _ez.insert ( 1 , _en )
    elif 5 == side % 6  :   ## bottom
        _en = _ez[-1] + width * ( _ez[-1] - _ez[-2] )
        if 0 <= width : _ez.append (      _en ) 
        else          : _ez.insert ( -1 , _en )

    ## create new histogram 
    hn  = h3_axes ( _ex , _ey , _ez , title = h3.GetTitle() , double = type ( h3 ) )

    _nx = hn.nbinsx ()
    _ny = hn.nbinsy ()
    _ny = hn.nbinsy ()
    
    if   0 == side % 6 :
        for j in range(1,_ny+1) :
            for k in range(1,_nz+1) : hn [  1 ,  j ,  k ] = value 
    elif 1 == side % 6 :
        for j in range(1,_ny+1) :
            for k in range(1,_nz+1) : hn [ -1 ,  j ,  k ] = value 
    elif 2 == side % 6 :
        for i in range(1,_nx+1) :
            for k in range(1,_nz+1) : hn [  i ,  1 ,  k ] = value 
    elif 3 == side % 6 :
        for i in range(1,_nx+1) :
            for k in range(1,_nz+1) : hn [  i , -1 ,  k ] = value 
    elif 4 == side % 6 :
        for i in range(1,_nx+1) :
            for j in range(1,_ny+1) : hn [  i ,  j ,  1 ] = value 
    elif 5 == side % 6 :
        for i in range(1,_nx+1) :
            for j in range(1,_ny+1) : hn [  i ,  j , -1 ] = value 
        
    for i in h3:

        if   0 == side % 6 : j = i[0]+1, i[1]   , i[2]
        elif 1 == side % 6 : j = i
        elif 2 == side % 6 : j = i[0]  , i[1]+1 , i[2] 
        elif 3 == side % 6 : j = i
        elif 4 == side % 6 : j = i[0]  , i[1]   , i[2]+1
        elif 5 == side % 6 : j = i 
        
        hn [ j ] = h3 [ i ]
        
    hn.ResetStats()
    return hn 

for _h in ( ROOT.TH3F , ROOT.TH3D ) :
    _h.add_fake_side       = _h3_add_fake_side_


# =============================================================================
## make the solution for equation   h(x)=v
#
#  @code 
#    >>> histo = ...
#    >>> value = ...
#    >>> solutions = histo.solve ( values )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-05-13
def _solve_ ( h1 , value ) :
    """ Make a solution for equation h(x)=v    
    >>> histo = ...
    >>> value = ...
    >>> solutions = histo.solve ( values )    
    """
    #
    value = float ( value )
    #
    solutions = []
    _size = len ( h1 ) 
    for i in h1.items() :

        ibin = i[0]
        x    = i[1]
        y    = i[2]

        yi   = y.value()
        if value == yi : 
            solutions.append  ( x.value () )
            continue

        di   = value - yi 
        xi   = x.value()

        j = ibin + 1
        if not j <= _size : continue

        yj  = h1.GetBinContent ( j )
        xj  = h1.GetBinCenter  ( j )
        
        dj   = value - yj 

        if 0 <= di * dj : continue  ## the same sign, skip 
        
        dd   = yi - yj
        
        if 0 == dd    : continue
        
        xs   = ( xi * dj  - xj * di ) / dd 

        solutions.append ( xs )


    return tuple  ( solutions ) 


# =============================================================================
## propose edge for "equal" bins 
#
#  @code
#
#    >>> histo = ....
#    >>> edges = histo.equal_edges ( 10 )
#
#  @endcode
#  @author Vanya BELYAVE Ivan.Belyaev@itep.ru
#  @date 2013-05-13
def _equal_edges_ ( h1 , num  , wmax = -1 ) :
    """ Propose the edges for ``equal-bins''
    >>> histo = ....
    >>> edges = histo.equal_edges ( 10 )    
    """
    if not isinstance ( num , integer_types ) :
        return TypeError, "'num' is not integer!"
    elif 1 >  num :
        return TypeError, "'num' is invalid!"
    elif 1 == num :
        _bins = [ h1.xmin() , h1.xmax() ] ## trivial binnig scheme
    else :
        
        ## integrate it! 
        _eff = h1.effic()
        
        _bins = [ h1.xmin()  ]
        d     = 1.0 / num
        
        for i in range ( 1 , num ) :
            vi  = float ( i ) / num
            eqs = _solve_ ( _eff , vi )
            if not eqs : continue
            _bins.append ( eqs[0] )
            
        _bins.append ( h1.xmax() )

    #
    ## split if needed
    #
    if 0 < wmax :
        _lst  = []
        _lst.append ( _bins[0] )
        for i in range ( 1 , len( _bins ) ) :
            b    = _bins[i]
            last = _lst[-1] 
            dist = b - last 
            if dist > wmax :
                n   = 1 + int ( dist // wmax ) 
                dlt = dist / float( n )
                for j in range(1,n) :
                    _lst.append ( last + j * dlt )
            _lst.append ( b )
            
        _bins = _lst 

    return tuple ( _bins ) 
    
ROOT.TH1F . solve       = _solve_                              
ROOT.TH1D . solve       = _solve_                              
ROOT.TH1F . equal_edges = _equal_edges_                              
ROOT.TH1D . equal_edges = _equal_edges_                         


# =============================================================================
# slices
# =============================================================================
## define 2D slice for 2D-histogram
def _h2_get_slice_ ( h2 , axis , ibins ) :
    """ Get 1D-slice for 2D-histogram
    >>> histo2 = ...
    >>> h1 = histo2.sliceX ( 1        ) ## along 1st X-bin
    >>> h2 = histo2.sliceY ( [2,4,19] ) ## along 2nd,4th and 19th Y-bins    
    """
    #
    if   1 == axis :
        if isinstance ( ibins , integer_types ) :
            if not ibins in h2.GetXaxis() :
                raise TypeError ( 'Illegal bin  index %s' % ibins ) 
            ibins = ibins,
        h_slice = h1_axis ( h2.GetYaxis () ,
                            title  = h2.GetTitle() + ":X-slice %s" % ibins , 
                            double = isinstance ( h2 , ROOT.TH2D )         )              
    elif 2 == axis : 
        if isinstance ( ibins , integer_types ) :
            if not ibins in h2.GetYaxis() :
                raise TypeError ( 'Illegal bin  index %s' % ibins )
            ibins = ibins,
        h_slice = h1_axis ( h2.GetXaxis () ,
                            title  = h2.GetTitle() + ":Y-slice %s" % ibins , 
                            double = isinstance ( h2 , ROOT.TH2D )         )  
            
    else :
        raise TypeError ( 'Illegal axis index %s' % axis ) 

    ## fill the slice histogram
    for i in h2 :
        
        ix = i[0]
        iy = i[1]
        
        if   1 == axis and ix in ibins : h_slice [ iy ] += h2 [ i ]
        elif 2 == axis and iy in ibins : h_slice [ ix ] += h2 [ i ]
        
    h_slice.ResetStats()
    return h_slice 

ROOT.TH2F . slice  = _h2_get_slice_ 
ROOT.TH2F . sliceX = lambda s , ibin : _h2_get_slice_ ( s , 1 , ibin ) 
ROOT.TH2F . sliceY = lambda s , ibin : _h2_get_slice_ ( s , 2 , ibin ) 

ROOT.TH2D . slice  = _h2_get_slice_
ROOT.TH2D . sliceX = lambda s , ibin : _h2_get_slice_ ( s , 1 , ibin ) 
ROOT.TH2D . sliceY = lambda s , ibin : _h2_get_slice_ ( s , 2 , ibin ) 


# =============================================================================
## "Slice" iteration over x-bin
#   @code
#   histo = ...
#   for i in histo.islice_X( 3 ) :
#            print i, histo[i] 
#   @endcode
def _h2_islice_X_ (  h2 , xbin ) :
    """``Slice'' iteration over x-bin
    >>> histo = ...
    >>> for i in histo.islice_X( 3 ) :
    ...         print (  i, histo[i] )
    """
    ax = h2.GetXaxis  ()
    nx = ax.GetNbins ()
    #
    if     1 <=  xbin <= nx : pass
    elif   1 <= -xbin <= nx : xbin += ( nx + 1 ) 
    else                    : raise IndexError 
    ##
    ay = h2.GetYaxis  ()
    ny = ay.GetNbins ()
    
    for ybin in range(1,ny+1) :
        yield xbin,ybin

# =============================================================================
## "Slice" iteration over y-bin
#   @code
#   histo = ...
#   for i in histo.islice_Y( 3 ) :
#            print ( i, histo[i] ) 
#   @endcode
def _h2_islice_Y_ (  h2 , ybin ) :
    """``Slice'' iteration over y-bin
    >>> histo = ...
    >>> for i in histo.islice_Y( 3 ) :
    ...         print ( i, histo[i] ) 
    """
    ay = h2.GetYaxis  ()
    ny = ay.GetNbins ()
    #
    if     1 <=  ybin <= ny : pass
    elif   1 <= -ybin <= ny : ybin += ( ny + 1 ) 
    else                    : raise IndexError 
    ##
    ax = h2.GetXaxis ()
    nx = ax.GetNbins ()
    
    for xbin in range(1,nx+1) :
        yield xbin,ybin
        
        
for _h in ( ROOT.TH2F , ROOT.TH2D ) :
    _h.islice_X = _h2_islice_X_ 
    _h.islice_Y = _h2_islice_Y_ 


# =============================================================================
## define 2D slice for 3D-histogram
def _h3_get_slice_ ( h3 , axis , ibins ) :
    """ Get 2D-slice for 3D-histogram    
    >>> histo3 = ...
    >>> h1 = histo3.sliceX ( 1          ) ## along 1st X-bin
    >>> h2 = histo3.sliceY ( [2,4,19]   ) ## along 2nd,4th and 19th Y-bins
    >>> h3 = histo3.sliceZ ( range(3,8) ) ## along 3-8th Z-bins

    """
    #
    if   1 == axis :
        if isinstance ( ibins , integer_types ) : 
            if not ibins in h3.GetXaxis() :
                raise TypeError ( 'Illegal bin  index %s' % ibins ) 
            ibins = ibins ,
        h_slice = h2_axes ( h3.GetYaxis () ,
                            h3.GetZaxis () ,
                            title  = h3.GetTitle()+ ":X-slice %s" % ibins , 
                            double = isinstance ( h3 , ROOT.TH3D )        )
        
    elif 2 == axis : 
        if isinstance ( ibins , integer_types ) : 
            if not ibins in h3.GetYaxis() :
                raise TypeError ( 'Illegal bin  index %s' % ibins ) 
            ibins = ibins ,
        h_slice = h2_axes ( h3.GetXaxis () ,
                            h3.GetZaxis () ,
                            title  = h3.GetTitle()+ ":Y-slice %s" % ibnis , 
                            double = isinstance ( h3 , ROOT.TH2D )        )             
    elif 3 == axis : 
        if isinstance ( ibins , integer_types ) : 
            if not ibins in h3.GetZaxis() :
                raise TypeError ( 'Illegal bin  index %s' % ibins )
            ibins = ibins ,
        h_slice = h2_axes ( h3.GetXaxis () , 
                            h3.GetYaxis () ,
                            title  = h3.GetTitle()+ ":Z-slice %s" % ibins , 
                            double = isinstance ( h3 , ROOT.TH2D ) )     
    else :
        raise TypeError ( 'Illegal axis index %s' % axis )

    ## fill the slice histogram
    for i in h3 :
        
        ix = i[0]
        iy = i[1]
        iz = i[2]
        
        if   1 == axis and ix in ibins : h_slice [ iy , iz ] += h3 [ i ]
        elif 2 == axis and iy in ibins : h_slice [ ix , iz ] += h3 [ i ]
        elif 3 == axis and iz in ibins : h_slice [ ix , iy ] += h3 [ i ] 

    h_slice.ResetStats()
    return h_slice 

ROOT.TH3F . slice  = _h3_get_slice_
ROOT.TH3F . sliceX = lambda s,ibin : _h3_get_slice_ ( s , 1 , ibin ) 
ROOT.TH3F . sliceY = lambda s,ibin : _h3_get_slice_ ( s , 2 , ibin ) 
ROOT.TH3F . sliceZ = lambda s,ibin : _h3_get_slice_ ( s , 3 , ibin ) 

ROOT.TH3D . slice  = _h3_get_slice_ 
ROOT.TH3D . sliceX = lambda s,ibin : _h3_get_slice_ ( s , 1 , ibin ) 
ROOT.TH3D . sliceY = lambda s,ibin : _h3_get_slice_ ( s , 2 , ibin ) 
ROOT.TH3D . sliceZ = lambda s,ibin : _h3_get_slice_ ( s , 3 , ibin ) 



# =============================================================================
## "Slice" iteration over x-bin
#   @code
#   histo = ...
#   for i in histo.islice_X( 3 ) :
#            print ( i, histo[i] ) 
#   @endcode
def _h3_islice_X_ (  h3 , xbin ) :
    """``Slice'' iteration over x-bin
    >>> histo = ...
    >>> for i in histo.islice_X( 3 ) :
    ...         print ( i, histo[i] ) 
    """
    ax = h3.GetXaxis ()
    nx = ax.GetNbins ()
    #
    if     1 <=  xbin <= nx : pass
    elif   1 <= -xbin <= nx : xbin += ( nx + 1 ) 
    else                    : raise IndexError 
    ##
    ay = h3.GetYaxis ()
    ny = ay.GetNbins ()
    
    az = h3.GetZaxis ()
    nz = az.GetNbins ()
    
    for zbin in range(1,nz+1) :
        for ybin in range(1,ny+1) :
            yield xbin,ybin,zbin
            
# =============================================================================
## "Slice" iteration over Y-bin
#   @code
#   histo = ...
#   for i in histo.islice_Y( 3 ) :
#            print ( i, histo[i] ) 
#   @endcode
def _h3_islice_Y_ (  h3 , ybin ) :
    """``Slice'' iteration over y-bin
    >>> histo = ...
    >>> for i in histo.islice_Y( 3 ) :
    ...         print ( i, histo[i] ) 
    """
    ay = h3.GetYaxis ()
    ny = ay.GetNbins ()
    #
    if     1 <=  ybin <= ny : pass
    elif   1 <= -ybin <= ny : ybin += ( ny + 1 ) 
    else                    : raise IndexError 
    ##
    ax = h3.GetXaxis ()
    nx = ax.GetNbins ()
    
    az = h3.GetZaxis ()
    nz = az.GetNbins ()
    
    for zbin in range(1,nz+1) :
        for xbin in range(1,nx+1) :
            yield xbin,ybin,zbin
            

# =============================================================================
## "Slice" iteration over z-bin
#   @code
#   histo = ...
#   for i in histo.islice_Z( 3 ) :
#            print ( i, histo[i] ) 
#   @endcode
def _h3_islice_Z_ (  h3 , zbin ) :
    """``Slice'' iteration over Z-bin
    >>> histo = ...
    >>> for i in histo.islice_Z( 3 ) :
    ...         print ( i, histo[i] ) 
    """
    az = h3.GetZaxis ()
    nz = az.GetNbins ()
    #
    if     1 <=  zbin <= nz : pass
    elif   1 <= -zbin <= nz : zbin += ( nz + 1 ) 
    else                    : raise IndexError 
    ##
    ay = h3.GetYaxis ()
    ny = ay.GetNbins ()
    
    ax = h3.GetXaxis ()
    nx = ax.GetNbins ()
    
    for ybin in range(1,ny+1) :
        for xbin in range(1,nx+1) :
            yield xbin,ybin,zbin


for _h in ( ROOT.TH3F , ROOT.TH3D ) :
    _h.islice_X = _h3_islice_X_ 
    _h.islice_Y = _h3_islice_Y_ 
    _h.islice_Z = _h3_islice_Z_ 

# ============================================================================
## make projections
#  see ROOT::TH2::ProjectionX
#  see ROOT::TH2::ProjectionY
def _h2_proj_ ( h2 , axis , *args ) :
    """ Make histogram projection along the axis
    >>> histo2 = ...
    >>> hx = histo2.projX ()
    >>> hy = histo2.projY ()    
    """
    if not args :
        args = (0,-1,'e')
    if   1 == axis : return h2.ProjectionX ( hID () , *args ) 
    elif 2 == axis : return h2.ProjectionY ( hID () , *args )
    
    raise TypeError ( 'Illegal axis index %s' % axis )

_h2_proj_ . __doc__ += '\n'
_h2_proj_ . __doc__ += ROOT.TH2.ProjectionX.__doc__
_h2_proj_ . __doc__ += ROOT.TH2.ProjectionY.__doc__


ROOT.TH2F . proj   = _h2_proj_
ROOT.TH2F . projX  = lambda s , *args : _h2_proj_ ( s , 1 , *args )
ROOT.TH2F . projY  = lambda s , *args : _h2_proj_ ( s , 2 , *args )

ROOT.TH2D . proj   = _h2_proj_
ROOT.TH2D . projX  = lambda s , *args : _h2_proj_ ( s , 1 , *args )
ROOT.TH2D . projY  = lambda s , *args : _h2_proj_ ( s , 2 , *args )


ROOT. TH3 . projX  = lambda s : s.Project3D ( 'xe'  )
ROOT. TH3 . projY  = lambda s : s.Project3D ( 'ye'  )
ROOT. TH3 . projZ  = lambda s : s.Project3D ( 'ze'  )
ROOT. TH3 . projXY = lambda s : s.Project3D ( 'xye' )
ROOT. TH3 . projYX = lambda s : s.Project3D ( 'yxe' )
ROOT. TH3 . projXZ = lambda s : s.Project3D ( 'xze' )
ROOT. TH3 . projZX = lambda s : s.Project3D ( 'zxe' )
ROOT. TH3 . projYZ = lambda s : s.Project3D ( 'yze' )
ROOT. TH3 . projZY = lambda s : s.Project3D ( 'zye' )



# =============================================================================
## convert TProfile to TH1D  (needed for the proper math) 
#  @see ROOT::TProfile
#  @see ROOT::TH1D 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-08-28
def _prof_asH1_ ( p ) :
    """ Convert TProfile to TH1D (needed for the proper math)
    >>> profile = ...
    >>> h1      = profile.asH1()
    """
    h1 = h1_axis ( p.GetXaxis() , double = True )
    ## copy the content bin-by-bin
    for i in h1 : h1[i] = p[i]
    ##
    return h1 

ROOT.TProfile . asH1     = _prof_asH1_ 
ROOT.TProfile . toH1     = _prof_asH1_
ROOT.TProfile . asHisto  = _prof_asH1_

ROOT.TProfile . __idiv__ = NotImplemented 
ROOT.TProfile . __imul__ = NotImplemented 
ROOT.TProfile . __iadd__ = NotImplemented 
ROOT.TProfile . __isub__ = NotImplemented
ROOT.TProfile . __abs__  = NotImplemented 


# ============================================================================
## get x-range for the function:
def _f_xminmax_ ( self ) :
    """ Get x-range for the function
    >>> f = ...
    >>> xmn,xmx = f.xminmax()
    """
    return self.GetXmin(),self.GetXmax()

# ============================================================================
## get minmax-range for the function:
def _f_minmax_ ( self ) :
    """ Get minmax-range for the function
    >>> f = ...
    >>> mn,mx = f.minmax()
    """
    return self.GetMinimum(),self.GetMaximum()

ROOT.TF1 .  minv    =  lambda s : s.GetMinimum()
ROOT.TF1 .  maxv    =  lambda s : s.GetMaximum()
ROOT.TF1 .  minmax  = _f_minmax_
ROOT.TF1 . xminmax  = _f_xminmax_

import ostap.stats.moments as moments
_tf1 = ROOT.TF1
if not hasattr ( _tf1 , 'mean'     ) : _tf1.mean     = moments.mean 
if not hasattr ( _tf1 , 'variance' ) : _tf1.variance = moments.variance 
if not hasattr ( _tf1 , 'rms'      ) : _tf1.rms      = moments.rms  
if not hasattr ( _tf1 , 'median'   ) : _tf1.median   = moments.median
if not hasattr ( _tf1 , 'mode'     ) : _tf1.mode     = moments.mode 
if not hasattr ( _tf1 , 'moment'   ) : _tf1.moment   = moments.moment
if not hasattr ( _tf1 , 'quantile' ) : _tf1.quantile = moments.quantile            



# =============================================================================
## make axis from bin-edges 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def axis_bins ( bins         ) :
    """ Make axis according to the binning 
    >>> bins = [ ... ] 
    >>> axis = axis_bins ( bins )      
    """
    #
    bins = set  ( bins )
    bins = [ i for i in bins ]
    bins.sort()
    #
    if 2 > len ( bins ) :
        raise AttributeError("axis_bins: insufficient length of bins: %s" % bins )
    #
    return ROOT.TAxis ( len ( bins ) - 1 , array.array ( 'd' , bins ) )
    #

# =============================================================================
## prepare "slice" for the axis
#  @code 
#    >>> axis  = ...
#    >>> naxis = axis[2:10] ## keep only bins from 2nd (inclusive) till 10 (exclusive)
#  @endcode     
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-20
def _axis_getslice_ ( self , i , j ) :
    """ Make a ``slice'' for the axis:

    >>> axis  = ...
    >>> naxis = axis[2:10] ## keep only bins from 2nd (inclusive) till 10 (exclusive)
    
    """
    nb = self.GetNbins()
    
    while i < 1 : i += nb
    while j < 1 : j += nb
    
    i = min ( nb , i )
    j = min ( nb , j )
    
    if i >= j : raise IndexError 
    
    edges = self.edges()
    
    return axis_bins ( edges [i-1:j] ) 


ROOT.TAxis. __getslice__ = _axis_getslice_

# =============================================================================
## get the guess for three major parameters of the histogram:
#    - number of signal events
#    - background level under the signal (per bin)
#    - background slope
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def histoGuess ( histo , mass , sigma ) :
    """ Get the guess for three major parameters of the histogram:
    - number of signal events
    - background level under the signal (per bin)
    - background slope
    - (minimal and maximal slopes to ensure 'positive' background)

    >>> histo = ...
    >>> signal, background, slope, slope_min, slope_max = histoGuess ( histo , mass , sigma )
    
    """
    tot0  = 0
    bin0  = 0
    tot3r = 0
    bin3r = 0
    tot3l = 0
    bin3l = 0
    tot4  = 0
    bin4  = 0

    axis  = histo.GetXaxis()   
    for ibin in axis:

        xbin = axis.GetBinCenter   ( ibin )
        dx   = ( xbin - mass ) / sigma
        val  = histo.GetBinContent ( ibin )
        
        if   abs ( dx ) < 2 :  ## +-2sigma
            tot0  += val 
            bin0  += 1
        elif   2 < dx  <  4 :  ## 'near' right sideband: 2 4 sigma 
            tot3r += val 
            bin3r += 1
        elif  -4 < dx  < -2 :  ## 'near' left sideband: -4 -2 sigma 
            tot3l += val 
            bin3l += 1
        else :
            tot4  += val 
            bin4  += 1

    bin3 = bin3r + bin3l + bin4
    tot3 = tot3r + tot3l + tot4 

    p00 = 0
    p03 = 0
    p04 = 0
    
    if bin3          : p03 = float(tot3)/bin3 
    if bin0 and bin3 : p00 = max ( float(tot0)-bin0*p03 , 0 ) 

    if bin3r and bin3l and tot3r and tot3l : 
        p04 = ( tot3r - tot3l ) / ( tot3r + tot3l ) / 3 / sigma 
    
    _xmin = axis.GetXmin() - 0.5* axis.GetBinWidth ( 1               )
    _xmax = axis.GetXmax() + 0.5* axis.GetBinWidth ( axis.GetNbins() )
    
    s1   = -1.0 / ( _xmin - mass )
    s2   = -1.0 / ( _xmax - mass )
    
    smin = min ( s1 , s2 )
    smax = max ( s1 , s2 )
    
    # if   p04 < smin : p04 = smin
    # elif p04 > smax : p04 = smax 
    
    return p00, p03, p04 , smin , smax  


ROOT.TH1F . histoGuess = histoGuess
ROOT.TH1D . histoGuess = histoGuess


# =============================================================================
## use likelihood in histogram fit ? 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def useLL ( histo         ,
            minc  = 10    ,
            diff  = 1.e-5 ) :
    """ Use Likelihood in histogram fit?
    """    
    minv = 1.e+9    
    for ibin in histo : 

        v = histo [ ibin ]

        if not natural_entry ( v ) : return False 
        minv = min ( minv , v.value() )
        
    return  minv < abs ( minc ) 


ROOT.TH1.useLL = useLL

# =============================================================================
## Natural histogram with all integer entries ? (generic version)
#  @code
#  h = ...
#  if h.natural() : print ( 'OK!') 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def allInts ( histo         ,
              diff  = 1.e-4 ) :
    """ Natural histogram with all integer entries ?
    >>> h = ...
    >>> if h.natural() : print ( 'OK!' ) 
    """
    
    for ibin in histo : 
        v = histo [ ibin ]
        if not natural_entry ( v ) : return False   ## RETURN 
        
    return True

ROOT.TH1.allInts = allInts
ROOT.TH1.natural = allInts

# =============================================================================
## Natural histogram with all integer entries ? (a bit faster vection for TH1)
#  @code
#  h1 = ...
#  if h1.natural() : print ( 'OK!') 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_allInts ( histo         ,
                  diff  = 1.e-4 ) :
    """ Natural histogram with all integer entries ?
    >>> h1 = ...
    >>> if h1.natural() : print(  'OK!') 
    """    
    for i in histo :  
        if not natural_entry ( histo[i] ) : return False   ## RETURN 
        
    return True

# =============================================================================
## check if  histogram is "natural" one, that represent event count
#  @see Ostap::Math::round 
#  @see Ostap::Math::islong
#  @see Ostap::Math::isint
#  @author Vanya BELYAEV Ivan.Belayev@itep.ru
#  @date   2015-03-31
def _natural_ ( histo ) :
    """ Check is the histogram is the ``natural'' one:
    - entries are non-negative integers
    - errors  are sqrt(entry)
    [ for null entries, the errors of  0 or 1 are both allowed] 
    
    >>> histo = ...
    >>> print ( 'natural? ', histo.natural()    ) 
    """
    ## loop over all histogram bins 
    for i in histo :
        v  = histo [ i ]
        if not natural_entry ( v ) : return False

    return True

ROOT.TH1.natural = _natural_

for t in ( ROOT.TH1F , ROOT.TH1D ) : 
    t.allInts = _h1_allInts
    t.natural = _h1_allInts


# =============================================================================
## Axis with uniform bins?
#  @code
#  axis = ...
#  print ( axis.uniform () ) 
#  @endcode
def _a_uniform_bins_ ( axis ) :
    """ Axis with unuform bins? 
    >>> axis = ...
    >>> print ( axis.uniform () ) 
    """
    edges = axis.GetXbins()
    nbins = len ( edges ) 
    if nbins <= 2  : return True
    ## 
    return all ( isequalf ( edges [ i     ] - edges [ i - 1 ] ,
                            edges [ i - 1 ] - edges [ i - 2 ] ) for i in range ( 2 , nbins ) )

ROOT.TAxis.uniform       = _a_uniform_bins_
ROOT.TAxis.uniform_bins  = _a_uniform_bins_

# =============================================================================
## check if histogram has uniform binnings
def _h_uniform_bins_ ( histo ) :
    """ Check if histogram has uniform binnings
    >>> histo = ...
    >>> uni   = histo.uniform_bins() 
    """
    xaxis = histo.GetXaxis()
    if not xaxis.uniform() : return False
    yaxis = histo.GetYaxis()
    if not yaxis.uniform() : return False
    zaxis = histo.GetZaxis()
    if not zaxis.uniform() : return False
    #
    return True

ROOT.TH1.uniform_bins = _h_uniform_bins_
ROOT.TH1.uniform      = _h_uniform_bins_

# =============================================================================
## same dims?
#  @see TH1::GetDimension 
def _h_same_dims_ ( histo , another ) :
    """ Two histograms with the same dimension?
    - see `ROOT.TH1.GetDimension`
    """
    
    return histo.GetDimension() == another.GetDimension()

    ## if   isinstance ( histo , ROOT.TH3 ) :
    ##     return isinstance ( another , ROOT.TH3 )
    ## elif isinstance ( histo , ROOT.TH2 ) :
    ##     return isinstance ( another , ROOT.TH2 ) \
    ##            and not isinstance ( another , ROOT.TH3 )
    ## elif isinstance ( histo , ROOT.TH1 ) :
    ##     return isinstance ( another , ROOT.TH1 ) \
    ##            and not isinstance ( another , ( ROOT.TH2 , ROOT.TH3 ) ) 
    ## return False 

ROOT.TH1. same_dims = _h_same_dims_

# =============================================================================
## same binning?
def _h_same_bins_ ( histo , another ) :
    """ Same binning for two histograms?
    >>> h1 = ...
    >>> h2 = ...
    >>> print ( h1.same_bins( h2 ) ) 
    """
    ## same dimensions ? 
    if not histo.same_dims ( another ) : return False
    ## 
    if   isinstance ( histo , ROOT.TH3 ) and isinstance ( another , ROOT.TH3 ) : 
        a1 = histo  .GetXaxis()
        a2 = another.GetXaxis()
        if a1 != a2 : return False
        a1 = histo  .GetYaxis()
        a2 = another.GetYaxis()
        if a1 != a2 : return False
        a1 = histo  .GetZaxis()
        a2 = another.GetZaxis()
        if a1 != a2 : return False
        return True
    elif isinstance ( histo , ROOT.TH2 ) and isinstance ( another , ROOT.TH2 ) : 
        a1 = histo  .GetXaxis()
        a2 = another.GetXaxis()
        if a1 != a2 : return False
        a1 = histo  .GetYaxis()
        a2 = another.GetYaxis()
        return True
    elif isinstance ( histo , ROOT.TH1 ) and isinstance ( another , ROOT.TH1 ) : 
        a1 = histo  .GetXaxis()
        a2 = another.GetXaxis()
        return True
    
    return False 

ROOT.TH1. same_bins = _h_same_bins_
    
# =============================================================================
## histogram dimension
#  Get dimension of the histogram
#  @code
#  histo = ...
#  d     =  histo.dim() 
#  @endcode
def _h_dim_ ( self  ) :
    """ Get dimension of the histogram
    >>> histo = ...
    >>> d     =  histo.dim() 
    """
    ## if   isinstance ( self , ROOT.TH3 ) : return 3
    ## elif isinstance ( self , ROOT.TH2 ) : return 2
    ## return 1
    #
    return self.GetDimension()
    
ROOT.TH1 .dim = _h_dim_

# =============================================================================
## Get tuple of bins-dimensions 
#  @code
#  histo = ...
#  nb    = histo.n_bins () 
#  @endcode
def _h_n_bins_ ( self ) :
    """ Get tuple of bins-dimensions 
    >>> histo = ...
    >>> nb    = histo.n_bins () 
    """
    if   isinstance ( self , ROOT.TH3 ) :
        ax = self.SetXaxis()
        ay = self.SetYaxis()
        az = self.SetZaxis()
        return ax.GetNbins() , ay.GetNbins(), az.xGetNbins() 
    elif isinstance ( self , ROOT.TH2 ) :
        ax = self.SetXaxis()
        ay = self.SetYaxis()
        return ax.GetNbins() , ay.GetNbins()
    
    ax = self.SetXaxis()
    return ax.GetNbins() ,

ROOT.TH1 .n_bins = _h_n_bins_

# =============================================================================
## transfrom the x-axis for the 1D-histogram
def _h1_transform_x_ ( h1 , fun , numbers = False , deriv = None ) :
    """ Transfrom the x-axis for the 1D-histogram
    """
    ax =  h1.GetXaxis()
    e  =  ax.edges()

    _fun = lambda x : float ( fun ( float(x) ) ) 
    
    ## new histogram 
    nh =  h1_axis ( ( _fun ( x ) for x in e  ) , title  = h1.GetTitle() , double = type(h1) )  

    if numbers :
        deriv = lambda x,h : 1
    elif deriv :
        deriv = lambda x,h : deriv ( x ) 
    else :
        from ostap.math.derivative import derivative as _deriv_
        deriv = lambda x,h : _deriv_ ( _fun , x , h )

    ## loop over the content of original histogram 
    for i,x,y in h1.items() :

        ## center bin 
        xc     = x.value()
        nc     = _fun ( xc )
        
        ## value
        dd     = deriv ( xc , 0.2 * x.error() )
        if iszero ( dd ) :
            logger.warning('derivative is zero, skip this bin')
            continue
        
        value  = VE(y) / dd 
        if value.isinf() or value.isnan() : continue
        
        ib = nh.findBin( nc )        
        if ib in nh : nh[ib] += value
        
    nh.ResetStats()
    return nh 

# =============================================================================
## make a histogram transfromation:
#  H(x)  ->  H'(y(x))
#  where transformation is defined  y=fun(x)
#  @code
#  h  = ...
#  from math import tanh 
#  ht = h.transform_X_numbers ( lambda x : tanh(x) ) 
#  @endcode
#  Histogram is transformed as raw histigram (no Jacobian correction is applied)
def _h1_transform_x_numbers_ ( h1 , fun ) :
    """ Make a histogram transfromation:
    H(x)  ->  H'(y(x))
    where transformation is defined  y=fun(x)
    >>> h  = ...
    >>> from math import tanh 
    >>>ht = h.transform_X_numbers ( lambda x : tanh(x) ) 
    Histogram is transformed as raw histigram (no Jacobian correction is applied)
    """
    return _h1_transform_x_ ( h1 , fun , numbers = True )

# =============================================================================
## make a histogram transformation:
#  H(x)  ->  H'(y(x)) / (dy/dx)
#  where transformation is defined  y=fun(x)
#  @code
#  h  = ...
#  from math import tanh 
#  ht = h.transform_X_function ( lambda x : tanh(x) ) 
#  @endcode    
#  - Histogram is transformed as function (Jacobian correction is applied)
#  - Derivative for Jacobian calculation can be explicitely specified, otherwise
#  they are calculated numerically 
def _h1_transform_x_function_ ( h1 , fun , deriv = None ) :
    """ Make a histogram transformation:
    H(x)  ->  H'(y(x)) / (dy/dx)
    where transformation is defined  y=fun(x)
    >>> h  = ...
    >>> from math import tanh 
    >>> ht = h.transform_X_function ( lambda x : tanh(x) ) 
    - Histogram is transformed as function (Jacobian correction is applied)
    - Derivative for Jacobian calculation can be explicitely specified, otherwise
    they are calculated numerically 
    """    
    return _h1_transform_x_ ( h1 , fun , numbers = False , deriv = deriv )
    

for t in ( ROOT.TH1F , ROOT.TH1D ) :
    t.transform_X_numbers  =  _h1_transform_x_numbers_
    t.transform_X_function =  _h1_transform_x_function_

# =============================================================================
## the actual function for text dump of the histogram
#  @code 
#  histo = ...
#  print ( dumpHisto  ( histo ) ) 
#  print ( histo.dump () ) 
#  print ( histo.dump ( 20 , 20 ) ) 
#  print ( histo.dump ( 20 , 20 , True ) ) 
#  @endcode
#  @see Ostap::Utils::Histo::.histoDump 
def dumpHisto ( histo , *args ) :
    """ Dump the histogram/profile in text format (a'la HBOOK)
    >>> histo = ...
    >>> print ( dumpHisto ( histo ) )     
    >>> print ( histo.dump () ) 
    >>> print ( histo.dump ( 20 , 20 ) ) 
    >>> print ( histo.dump ( 20 , 20 , True ) ) 
    
    Uses:
    - see Ostap.Utils.Histos.histoDump 
    
    """
    if not args :
        from ostap.utils.basic import terminal_size 
        w , h = terminal_size()
        if 0 < h and 0 < w :
            w = max ( 20 ,       w - 25     ) 
            h = max ( 20 , int ( h * 0.55 ) )
            return Ostap.Utils.Histos.histoDump ( histo , w , h )
                        
    return Ostap.Utils.Histos.histoDump ( histo , *args )


dumpHisto .__doc__ = '\n'  + Ostap.Utils.Histos.histoDump . __doc__
for t in  ( ROOT.TH1D             ,
            ROOT.TH1F             ,
            ROOT.TProfile         ) :
    for method in ( 'dump'       ,
                    'dumpHisto'  ,
                    'dumpASCII'  ,
                    'dumpAsText' ) :
        if not hasattr ( t , method ) :
            setattr ( t , method , dumpHisto )

# =============================================================================
## Hashing function for the histograms
#  @code
#  histo = ...
#  h     = hash ( histo ) 
#  @endcode
def _h_hash_ ( histo) :
    """ Hash for the histogram 
    >>> histo = ...
    >>> h     = hash ( histo ) 
    """
    hv = hash ( ( histo.GetName  () , histo.GetTitle() , type ( histo ) ) ) 
    for item in histo.iteritems() :
        hv = hash ( item + ( hv, ) )
    return hv

ROOT.TH1.__hash__ = _h_hash_ 
# =============================================================================

# =============================================================================
## Get "Riemann sum" of the histogram 
#  \f$ \sum_{i} f_i \Delta x_i\f$ 
def _h1_riemann_sum_ ( h1 ) :
    """ Get 'Riemann sum' of the historgam
    """
    ax     = h1.GetXaxis()
    result = 0.0
    for ix in h1 :
        volume  = ax.GetBinWidth   ( ix )        
        content = h1.GetBinContent ( ix )  
        result += content * volume 
    return result

# =============================================================================
## Get "Riemann sum" of the histogram 
#  \f$ \sum_{i} f_i \Delta x_i \Delta y_i\f$ 
def _h2_riemann_sum_ ( h2 ) :
    """ Get 'Riemann sum' of the historgam
    """
    ax     = h2.GetXaxis()
    ay     = h2.GetYaxis()
    result = 0.0
    for ix, iy in h2 :
        volume  = ax.GetBinWidth   ( ix      ) * ay.GetBinWidth ( iy )        
        content = h2.GetBinContent ( ix , iy ) 
        result += content * volume 
    return result

# =============================================================================
## Get "Riemann sum" of the histogram 
#  \f$ \sum_{i} f_i \Delta x_i \Delta y_i \Delta z_i \f$ 
def _h3_riemann_sum_ ( h3 ) :
    """ Get 'Riemann sum' of the historgam
    """
    ax     = h3.GetXaxis()
    ay     = h3.GetYaxis()
    az     = h3.GetZaxis()
    result = 0.0
    for ix, iy, iz  in h3 :
        volume  = ax.GetBinWidth   ( ix ) * ay.GetBinWidth ( iy ) * az.GetBinWidth ( iz )        
        content = h3.GetBinContent ( ix , iy , iz ) 
        result += content * volume 
    return result

ROOT.TH1.riemann_sum = _h1_riemann_sum_
ROOT.TH2.riemann_sum = _h2_riemann_sum_
ROOT.TH3.riemann_sum = _h3_riemann_sum_

# =============================================================================
## does this histogram reprsent the density?
def _h_isDensity_ ( h1 ) :
    """ Does this histogram reprsent the density?
    """
    return isequalf ( h1.riemann_sum () , 1.0 ) 
    
# =============================================================================
## represent histogram as ``density''
#  - the function with unit integral over the range
def _h_density_ ( h1 ) :
    """ Represent histogram as  `density'
    - the function with unit integral over the range
    """
    if h1.is_density () :
        ## it is already density!, return the clone
        return h1.clone()
    
    ## take into account bin width
    hh  = h1.rescale_bins ( 1 )
    
    ##  rescale to unit integral in range
    rs  = hh.riemann_sum  () 
    hh /= rs

    ## check if properly rescaled: 
    rs  = hh.riemann_sum  () 
    if not isequal ( rs , 1.0 ) : hh /= rs
        
    return hh

for t in ( ROOT.TH1F , ROOT.TH1D ,
           ROOT.TH2F , ROOT.TH2D ,
           ROOT.TH3F , ROOT.TH3D ) :
    t.density    = _h_density_
    t.isDensity  = _h_isDensity_
    t.is_density = _h_isDensity_
    
# =============================================================================
## Fill the histogram from iterable/generator/...
#  @code
#
#  h1 = ... ## 1D histogram 
#  h1.fill_loop (  ( random.gauss(5,1) , 1 ) for i in range ( 1000 ) ) ## use generator 
#  
#  h2  = ... ## 2D histogram#
#  h2.fill_loop (  ( random.gauss(5,1) , random.gauss(4,1)     ) for i in range ( 1000 ) ) ## use generator 
#  h2.fill_loop (  ( random.gauss(5,1) , random.gauss(4,1) , 1 ) for i in range ( 1000 ) ) ## use generator 
#
#  h3  = ... ## 2D histogram#
#  h3.fill_loop (  ( random.gauss(5,1) , random.gauss(4,1) , random.gauss (3,0.5)     ) for i in range ( 1000 ) ) ## use generator 
#  h3.fill_loop (  ( random.gauss(5,1) , random.gauss(4,1) , random.gauss (3,0.5) , 1 ) for i in range ( 1000 ) ) ## use generator 
#
#  @endcode 
def _h_fill_loop_ ( histo , filler ) :
    """ Fill the histogram from iterable/generator/...
    >>> 
    >>> h1 = ... ## 1Dhistogram 
    >>> h1.fill_loop (  ( random.gauss(5,1) , 1 ) for i in range ( 1000 ) ) ## use generator 
    
    >>> h2  = ... ## 2D-histogram
    >>> h2.fill_loop (  ( random.gauss(5,1) , random.gauss(4,1)     ) for i in range ( 1000 ) ) ## use generator 
    >>> h2.fill_loop (  ( random.gauss(5,1) , random.gauss(4,1) , 1 ) for i in range ( 1000 ) ) ## use generator 
    
    >>> h3  = ... ## 3D-histogram
    >>> h3.fill_loop (  ( random.gauss(5,1) , random.gauss(4,1) , random.gauss (3,0.5)     ) for i in range ( 1000 ) ) ## use generator 
    >>> h3.fill_loop (  ( random.gauss(5,1) , random.gauss(4,1) , random.gauss (3,0.5) , 1 ) for i in range ( 1000 ) ) ## use generator 
    """

    num = 0 
    for item in filler :
        num += 1
        histo.Fill ( *item ) 

    return num  

ROOT.TH1.fill_loop = _h_fill_loop_

# ==============================================================================
## split the histogram into <code>n</code> sub-histograms 
#  @code
#  h1 = ...
#  histos = h1.split(10) 
#  @endcode 
def _h1_split1_ ( h1 , n ) :
    """ Split the histogram into `n` sub-histograms
    >>> h1 = ...
    >>> histos = h1.split(10) 
    """
    assert isinstance ( n , integer_types ) and 0 < n ,\
           'Invalid number of sub-histograms %s' % n
    
    histos = []
    N      = len ( h1 ) + 1 
    from ostap.utils.utils import divide as _divide 
    bins   = _divide ( n , range ( 1 , N ) )
    for item in bins  :
        lst = list ( item )
        if not lst : return
        h   = h1[ lst [ 0 ] : lst [ -1 ] + 1 ]
        histos.append ( h )
    return tuple ( histos ) 

# ==============================================================================
## split the histogram into sub-histograms with at most <code>n</code> bins 
#  @code
#  h1 = ...
#  histos = h1.split_bins(10) 
#  @endcode 
def _h1_split2_ ( h1 , n ) :
    """ Split the histogram into sub-histograms with at most <code>n</code> bins 
    >>> h1 = ...
    >>> histos = h1.split_bins ( 10 ) 
    """
    assert isinstance ( n , integer_types ) and 0 < n ,\
           'Invalid number of bins  %s' % n

    histos = []
    N      = len ( h1 ) + 1 
    from ostap.utils.utils import chunked as _chunked
    bins   = _chunked ( range ( 1 , N ) , n )
    for item in bins  :
        lst = list ( item )
        if not lst : return
        h   = h1[ lst [ 0 ] : lst [ -1 ] +1 ]
        histos.append ( h )
    return tuple ( histos ) 

for h in ( ROOT.TH1F , ROOT.TH1D ) :
    h.split      = _h1_split1_
    h.split_bins = _h1_split2_


# =============================================================================
if root_info < ( 3 , 0 ) : 
    ## Number of entries (as int)
    #  @code
    #  histo
    #  n = histo.nEntries() 
    #  @endcode
    #  @see TH1::GetEntries 
    def _h1_nEntries_ ( histo ) :
        """ Number of entries (as int)
        >>> histo = ..
        >>> n = histo.nEntries() 
        """
        return long ( histo.GetEntries() )
else :
    ## Number of enbtries (as int)
    #  @code
    #  histo
    #  n = histo.nEntries() 
    #  @endcode
    #  @see TH1::GetEntries 
    def _h1_nEntries_ ( histo ) :
        """ Number of enbtries (as int)
        >>> histo = ..
        >>> n = histo.nEntries() 
        """
        return int ( histo.GetEntries() )
    
ROOT.TH1.nEntries = _h1_nEntries_

# =============================================================================
## Scale axis (change units) for 1D histogram
#  @code
#  histo1 = ...
#  histo2 = histo1.scale_axis ( 1000 ) 
#  @endcode 
def _h1_scale_axis_ ( h1 , scale ) :
    """ Scale axis (change units) for 1D histogram
    >>> histo1 = ...
    >>> histo2 = histo1.scale_axis ( 1000  ) 
    """
    xaxis   = h2.GetXaxis() 
    newx    = xaxis . scale ( scale )
    ## create new histogram 
    result  = h1_axis ( newx , title = h1.title , double = type ( h1 ) )
    ## copy content of the histogram
    for i in h1 : result [ i ] = h1 [ i ]
    ##
    return result

# =============================================================================
## Scale axes (change units) for 2D histogram
#  @code
#  histo1 = ...
#  histo2 = histo1.scale_axes ( 1000 , 5000 ) 
#  @endcode 
def _h2_scale_axes_ ( h2 , scalex , scaley = 1.0 ) :
    """ Scale axes (change units) for 2D histogram
    >>> histo1 = ...
    >>> histo2 = histo1.scale_axes ( 1000  , 5000 ) 
    """
    xaxis   = h2.GetXaxis() 
    newx    = xaxis . scale ( scalex )
    yaxis   = h2.GetYaxis() 
    newy    = yaxis . scale ( scaley )
    ## create new histogram 
    result  = h2_axes ( newx , newy  , title = h2.title , double = type ( h2 ) )
    ## copy content of the histogram
    for i in h2 : result  [ i ] = h2 [ i ]
    ##
    return result 

# =============================================================================
## Scale axes (change units) for 3D histogram
#  @code
#  histo1 = ...
#  histo2 = histo1.scale_axes ( 1 , 100 , 3 ) 
#  @endcode 
def _h3_scale_axes_ ( h3 , scalex , scaley = 1.0 , scalez = 1.0 ) :
    """ Scale axes (change units) for 3D histogram
    >>> histo1 = ...
    >>> histo2 = histo1.scale_axes ( 1000  , 5000 , 1000 ) 
    """
    xaxis   = h3.GetXaxis() 
    newx    = xaxis . scale ( scalex )
    yaxis   = h3.GetYaxis() 
    newy    = yaxis . scale ( scaley )
    zaxis   = h3.GetZaxis() 
    newz    = zaxis . scale ( scalez )
    ## create new histogram 
    result  = h3_axes ( newx , newy  , newz , title = h3.title , double = type ( h3 ) )
    ## copy content of the histogram
    for i in h3 : result [ i ] = h3 [ i ]
    ##
    return result 

ROOT.TH1F. scale_axis = _h1_scale_axis_
ROOT.TH1D. scale_axis = _h1_scale_axis_

ROOT.TH2F. scale_axes = _h2_scale_axes_
ROOT.TH2D. scale_axes = _h2_scale_axes_

ROOT.TH3F. scale_axes = _h3_scale_axes_
ROOT.TH3D. scale_axes = _h3_scale_axes_


# =============================================================================
## keys used for booing 1/2/3-dimensiona histograms via the `histo_book` method 
histo_keys = ( 'xbins' , 'nbinsx' , 'binsx' , 'nbins' ,
               'ybins' , 'nbinsy' , 'binsy' , 
               'zbins' , 'nbinsz' , 'binsz' ,
               'xmin'  , 'xmax'   ,
               'ymin'  , 'ymax'   ,
               'zmin'  , 'zmax'   ,
               'title' ,
               'xtitle', 'titlex' ,
               'ytitle', 'titley' ,
               'ztitle', 'titlez' )
               
# =============================================================================
## helper method to book/create 1,2&3-dimension histograms
#  Valid keys: 
#  @see ostap.histos.histos.histo_keys 
def histo_book2 ( ranges , kwargs ) :
    """ Helper method to book/create 1,2&3-dimension histograms
    Valid keys: 
    - see `ostap.histos.histos.histo_keys`
    """
    
    nvars = len ( ranges  ) 
    assert 1 <= nvars <= 3 , 'Invalid histogram dimension %s' % str ( items )

    if 1 == nvars :
        xbins = kwargs.pop ( 'xbins' , kwargs.pop ( 'nbinsx' , kwargs.pop ( 'binsx' , kwargs.pop ( 'nbins' , 100 ) ) ) )
    elif 2 == nvars :
        xbins = kwargs.pop ( 'xbins' , kwargs.pop ( 'nbinsx' , kwargs.pop ( 'binsx' , 50 ) ) ) 
        ybins = kwargs.pop ( 'ybins' , kwargs.pop ( 'nbinsy' , kwargs.pop ( 'binsy' , 50 ) ) ) 
    elif 3 == nvars :
        xbins = kwargs.pop ( 'xbins' , kwargs.pop ( 'nbinsx' , kwargs.pop ( 'binsx' , 20 ) ) ) 
        ybins = kwargs.pop ( 'ybins' , kwargs.pop ( 'nbinsy' , kwargs.pop ( 'binsy' , 20 ) ) ) 
        zbins = kwargs.pop ( 'zbins' , kwargs.pop ( 'nbinsz' , kwargs.pop ( 'binsz' , 20 ) ) )
        
    assert isinstance ( xbins , integer_types ) and 0 < xbins , "Invalid xbins setting!"
    if 2 <= nvars : assert isinstance ( ybins , integer_types ) and 0 < ybins , "Invalid ybins setting!"
    if 3 <= nvars : assert isinstance ( zbins , integer_types ) and 0 < zbins , "Invalid zbins setting!"
        
    xvar , xrng = ranges [ 0 ]
    xmin , xmax = xrng if xrng else ( None , None ) 
    xmin        = kwargs.pop ( 'xmin' , xmin )
    xmax        = kwargs.pop ( 'xmax' , xmax )    
    assert  isinstance ( xmin , num_types ) \
        and isinstance ( xmax , num_types ) and xmin < xmax , "Invalid xmin/xmax setting!"
    
    if 2 <= nvars:
        yvar , yrng = ranges [ 1 ]
        ymin , ymax = yrng if yrng else ( None , None ) 
        ymin        = kwargs.pop ( 'ymin' , ymin )
        ymax        = kwargs.pop ( 'ymax' , ymax )
        assert  isinstance ( ymin , num_types ) \
            and isinstance ( ymax , num_types ) and ymin < ymax , "Invalid ymin/ymax setting!"

    if 3 <= nvars:
        zvar , zrng = ranges [ 2 ]
        zmin , zmax = zrng if zrng else ( None , None ) 
        zmin        = kwargs.pop ( 'zmin' , zmin )
        zmax        = kwargs.pop ( 'zmax' , zmax )
        assert  isinstance ( zmin , num_types ) \
            and isinstance ( zmax , num_types ) and zmin < zmax , "Invalid zmin/zmax setting!"

    histo = None 
    if   1 == nvars :
        
        title  = kwargs.pop ( 'title' , "x=`%s'" % xvar )
        histo  = ROOT.TH1D  ( hID()   , title , xbins  , xmin , xmax )  ; histo.Sumw2()
        
        xtitle = kwargs.pop ( 'xtitle' , kwargs.pop ( 'titlex' , xvar ) )        
        histo.GetXaxis().SetTitle ( xtitle )
        
    elif 2 == nvars :
        
        title  = kwargs.pop ( 'title' , "x=`%s' y=`%s'" % ( xvar, yvar ) ) 
        histo  = ROOT.TH2F  ( hID()   , title ,
                              xbins   , xmin , xmax ,
                              ybins   , ymin , ymax ) ; histo.Sumw2()
        
        xtitle = kwargs.pop ( 'xtitle' , kwargs.pop ( 'titlex' , xvar ) )
        ytitle = kwargs.pop ( 'ytitle' , kwargs.pop ( 'titley' , yvar ) )
        
        histo.GetXaxis().SetTitle ( xtitle )
        histo.GetYaxis().SetTitle ( ytitle )
        
    elif 3 == nvars :
        
        title  = kwargs.pop ( 'title' , "x=`%s' y=`%s' z=`%s'" % ( xvar, yvar , zvar ) ) 
        histo  = ROOT.TH3F  ( hID()   , title ,
                              xbins   , xmin , xmax ,
                              ybins   , ymin , ymax , 
                              zbins   , zmin , zmax ) ; histo.Sumw2()
        
        xtitle = kwargs.pop ( 'xtitle' , kwargs.pop ( 'titlex' , xvar ) )
        ytitle = kwargs.pop ( 'ytitle' , kwargs.pop ( 'titley' , yvar ) )
        ztitle = kwargs.pop ( 'ztitle' , kwargs.pop ( 'titlez' , zvar ) )
        
        histo.GetXaxis().SetTitle ( xtitle )
        histo.GetYaxis().SetTitle ( ytitle )        
        histo.GetZaxis().SetTitle ( ztitle )        

    return histo 
# =============================================================================

# =============================================================================
## helper method to book/create 1,2&3-dimension histograms
#  @code
#  r1 = [  ('xvar', (0,1) } ]
#  h1 = histo_book ( r1 , xmin = -1 , xbins =  100 )
#
#  r2 = [  ('xvar', (0,1) } , ( 'yvar' , (10,20) ) ] 
#  h2 = histo_book ( r2 , xmin = -1 , xbins =  100 , ymax = 30 )
# 
#  r3 = [  ('xvar', (0,1) } , ( 'yvar' , (10,20) ) , ( 'z' , (0,100) ) ] 
#  h3 = histo_book ( r2 , xmin = -1 , xbins =  100 , ymax = 30 , zbins = 10 )
#  @endocode
#  Valid keys: 
#  @see ostap.histos.histos.histo_keys 
def histo_book ( ranges , **kwargs ) :
    """ Helper method to book/create 1,2&3-dimension histograms
    Valid keys: 
    - see ostap.histos.histos.histo_keys 
    >>> r1 = [  ('xvar', ( 0 , 1 ) } ]
    >>> h1 = histo_book ( r1 , xmin = -1 , xbins =  100 )
    
    >>> r2 = [  ('xvar', (0,1) } , ( 'yvar' , (10,20) ) ] 
    >>> h2 = histo_book ( r2 , xmin = -1 , xbins =  100 , ymax = 30 )
    
    >>> r3 = [  ('xvar', (0,1) } , ( 'yvar' , (10,20) ) , ( 'z' , (0,100) ) ] 
    >>>h3 = histo_book ( r2 , xmin = -1 , xbins =  100 , ymax = 30 , zbins = 10 )
    """
    ## 
    kw = cidict ( transform = cidict_fun , **kwargs )
    ## 
    histo = histo_book2 ( ranges , kw )
    if kw : 
        rows = [ ( 'Argument' , 'Value' ) ]
        for k , v  in loop_items ( kw ) :
            row = k , str ( v )
            rows.append ( row )
        title = 'histo_book: %d unused arguments' % len ( kw ) 
        table = T.table ( rows , title = 'Unused argumens' , prefix = '# ' , alignment = 'll' )    
        logger.warning ( '%s\n%s' % ( title , table ) )
    ## 
    return histo
    
# =============================================================================
## Get a kind of ROC-like curve/graph from two histograms
def _h1_roc_ ( h1 , h2 ) :
    """ Get a kind of ROC-like curve/graph from two histograms
    """
    h1sum = h1.sumv ()
    h2sum = h2.sumv ()

    import ostap.histos.graphs
    graph = ROOT.TGraphErrors ( len ( h1 ) + 2 )
    
    graph [  0 ] = h1sum [  0 ] , h2sum [ 0  ] 
    graph [ -1 ] = h1sum [ -1 ] , h2sum [ -1 ] 
    
    for i , x , y1 in h1sum.items() :        
        y2 = h2sum ( x.value () )        
        graph [ i ] = y1 , y2 
        
    return graph 

# =============================================================================
## print histogram summary as table
#  @code
#  histo = ...
#  print ( histo.table() )   
#  @endcode
def _h1_table_ ( h1 , title = '' , prefix = '' , width = 5 , precision = 3 ) :
    """ Print histogram summary as table 
    >>> histo = ...
    >>> print ( histo.table() )   
    """
    from ostap.logger.pretty import pretty_float , pretty_ve 
    rows = [ ( '' , 'value' , 'scale' ) ]

    row = 'Name' , h1.GetName() , '' 
    rows.append ( row )

    htitle = h1.GetTitle ()
    if htitle : 
        row = 'Title' , htitle , '' 
        rows.append ( row )
    
    rdir  = h1.GetDirectory()
    if rdir : 
        row = 'directory' , rdir.GetName() , '' 
        rows.append ( row )
    groot = ROOT.ROOT.GetROOT()
    row   = 'gROOT?' , 'True' if groot is rdir else 'False' , '' 
    rows.append ( row )
                
    row = 'SumW2' , 'True' if h1.GetSumw2() else 'False' , '' 
    rows.append ( row )

    xmin , expo = pretty_float ( h1.xmin () , width = width , precision = precision )    
    row = 'x-min' , xmin , '10^%+d' % expo if expo else '' 
    rows.append ( row )
    
    xmax , expo = pretty_float ( h1.xmax () , width = width , precision = precision )
    row = 'x-max' , xmax , '10^%+d' % expo if expo else '' 
    rows.append ( row )

    nbins = h1.GetXaxis().GetNbins() 
    row = '#bins' , '%d/%s' %  ( nbins , 'uniform' if h1.GetXaxis().uniform() else 'no-uniform' ) , '' 
    rows.append ( row )
    
    row = '#entries' , '%d' % h1.GetEntries() , '' 
    rows.append ( row )

    unfl = h1.GetBinContent ( 0 )
    if unfl :
        unfl , expo = pretty_float ( unfl , width = width , precision = precision )
        row = 'underflow', '%s' % unfl  ,  '10^%+d' % expo if expo else '' 
        rows.append ( row )

    ovfl = h1.GetBinContent ( nbins + 1 )
    if ovfl :
        ovfl , expo = pretty_float ( ovfl , width = width , precision = precision )        
        row = 'overflow' , ovfl , '10^%+d' % expo if expo else '' 
        rows.append ( row )
            
    integral  = h1.Integral()
    nEff      = h1.GetEffectiveEntries()
    
    if integral != h1.GetEntries () or integral != nEff : 
        integral , expo  = pretty_float ( integral , width = width , precision = precision ) 
        row = 'in-range' , integral ,  '10^%+d' % expo if expo else '' 
        rows.append ( row )

    if nEff != h1.GetEntries() or nEff != h1.Integral() :
        nEff , expo  = pretty_float ( nEff , width = width , precision = precision ) 
        row = '#eff' , nEff ,  '10^%+d' % expo if expo else '' 
        rows.append ( row )

    row = 'natural?'   , '%s' % h1.natural () , '' 
    rows.append ( row )

    row = 'all ints?' , '%s' % h1.allInts ()  , '' 
    rows.append ( row )

    row = 'all >0/>=0/<0/<=0/=0?', '%s/%s/%s/%s/%s' % ( 'True' if h1.all_positive    () else 'False' ,
                                                        'True' if h1.all_nonnegative () else 'False' ,
                                                        'True' if h1.all_negative    () else 'False' ,
                                                        'True' if h1.all_nonpositive () else 'False' ,
                                                        'True' if h1.all_zero        () else 'False' ) , ''
    rows.append ( row )

    
    rsum , expo  = pretty_float ( h1.riemann_sum ()   , width = width , precision = precision ) 
    row = 'Riemann sum' , rsum  ,  '10^%+d' % expo if expo else '' 
    rows.append ( row )

    if h1.nEff() and not h1.all_zero () :
        
        row = 'density?'   , '%s' % h1.isDensity() , '' 
        rows.append ( row )
        
        mean , expo = pretty_ve ( h1.mean() , width = width , precision = precision , parentheses = False )
        row = 'x-mean' , mean, '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        rms , expo = pretty_ve ( h1.rms  () , width = width , precision = precision , parentheses = False )
        row = 'x-rms' , rms, '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        skew , expo = pretty_ve ( h1.skewness  () , width = width , precision = precision , parentheses = False )
        row = 'skewness' , skew, '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        kurt , expo = pretty_ve ( h1.kurtosis  () , width = width , precision = precision , parentheses = False )
        row = 'kurtosis' , kurt, '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        stat = h1.stat()
        
        ymin , expo = pretty_float ( stat.min () , width = width , precision = precision )    
        row = 'y-min' , ymin , '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        ymax , expo = pretty_float ( stat.max () , width = width , precision = precision )    
        row = 'y-max' , ymax , '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        ymean , expo = pretty_ve ( stat.mean  () , width = width , precision = precision , parentheses = False )
        row = 'y-mean' , ymean , '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        yrms , expo = pretty_float ( stat.rms () , width = width , precision = precision )    
        row = 'y-rms' , yrms , '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        
    has3 = False
    for row in rows [1:] :
        if row [ 2 ] :
            has3 = True
            break
    
    if not has3 :
        new_rows = []
        for row in rows : new_rows.append ( row [ :-1] ) 
        rows = new_rows 
            

    title = title if title else 'Summary of histogram %s(%s)' % ( typename ( h1 ) , h1.GetName() ) 
    
    table = T.table ( rows ,title = title , prefix = prefix , alignment = 'lc' ) 
    return table 

# =============================================================================
## print histogram summary as table
#  @code
#  histo = ...
#  print ( histo.table() )   
#  @endcode
def _h2_table_ ( h2 , title = '' , prefix = '' , width = 5 , precision = 3 ) :
    """ Print histogram summary as table 
    >>> histo = ...
    >>> print ( histo.table() )   
    """
    from ostap.logger.pretty import pretty_float , pretty_ve 
    rows = [ ( '' , 'value' , 'scale' ) ]

    row = 'Name' , h2.GetName() , '' 
    rows.append ( row )

    htitle = h2.GetTitle ()
    if htitle : 
        row = 'title' , htitle , ''
        rows.append ( row )
    
    rdir  = h2.GetDirectory()
    if rdir : 
        row = 'directory' , rdir.GetName() , '' 
        rows.append ( row )
    groot = ROOT.ROOT.GetROOT()
    row   = 'gROOT?' , 'True' if groot is rdir else 'False' , '' 
    rows.append ( row )
                
    row = 'SumW2' , 'True' if h2.GetSumw2() else 'False' , '' 
    rows.append ( row )

    xmin , expo = pretty_float ( h2.xmin () , width = width , precision = precision )    
    row = 'x-min' , xmin , '10^%+d' % expo if expo else '' 
    rows.append ( row )
    
    xmax , expo = pretty_float ( h2.xmax () , width = width , precision = precision )
    row = 'x-max' , xmax , '10^%+d' % expo if expo else '' 
    rows.append ( row )

    ymin , expo = pretty_float ( h2.ymin () , width = width , precision = precision )    
    row = 'y-min' , ymin , '10^%+d' % expo if expo else '' 
    rows.append ( row )
    
    ymax , expo = pretty_float ( h2.ymax () , width = width , precision = precision )
    row = 'y-max' , ymax , '10^%+d' % expo if expo else '' 
    rows.append ( row )

    xaxis = h2.GetXaxis() 
    row = '#xbins' , '%d/%s' %  ( xaxis.GetNbins()  , 'uniform' if xaxis.uniform() else 'no-uniform' ) , '' 
    rows.append ( row )
    yaxis = h2.GetYaxis() 
    row = '#ybins' , '%d/%s' %  ( yaxis.GetNbins()  , 'uniform' if yaxis.uniform() else 'no-uniform' ) , '' 
    rows.append ( row )
    
    row = '#entries' , '%d' % h2.GetEntries()  , '' 
    rows.append ( row )

    integral  = h2.Integral()
    nEff      = h2.GetEffectiveEntries()
    
    if integral != h2.GetEntries () or integral != nEff : 
        integral , expo  = pretty_float ( integral , width = width , precision = precision ) 
        row = 'in-range' , integral ,  '10^%+d' % expo if expo else '' 
        rows.append ( row )

    if nEff != h2.GetEntries() or nEff != h2.Integral() :
        nEff , expo  = pretty_float ( nEff , width = width , precision = precision ) 
        row = '#eff' , nEff ,  '10^%+d' % expo if expo else '' 
        rows.append ( row )

    row = 'natural?'   , '%s' % h2.natural () , '' 
    rows.append ( row )

    row = 'all ints?' , '%s' % h2.allInts ()  , '' 
    rows.append ( row )

    row = 'all >0/>=0/<0/<=0/=0?', '%s/%s/%s/%s/%s' % ( 'True' if h2.all_positive    () else 'False' ,
                                                        'True' if h2.all_nonnegative () else 'False' ,
                                                        'True' if h2.all_negative    () else 'False' ,
                                                        'True' if h2.all_nonpositive () else 'False' ,
                                                        'True' if h2.all_zero        () else 'False' ) , '' 
    rows.append ( row )

    rsum , expo  = pretty_float ( h2.riemann_sum ()   , width = width , precision = precision ) 
    row = 'Riemann sum' , rsum  ,  '10^%+d' % expo if expo else '' 
    rows.append ( row )

    if h2.nEff() and not h2.all_zero () : 
        
        row = 'density?'   , '%s' % h2.isDensity() , '' 
        rows.append ( row )
        
        xmean, ymean = h2.mean () 
        xrms , yrms  = h2.rms  () 
        
        xmean , expo = pretty_ve ( xmean , width = width , precision = precision , parentheses = False )
        row = 'x-mean' , xmean, '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        xrms , expo = pretty_ve ( xrms  , width = width , precision = precision , parentheses = False )
        row = 'x-rms' , xrms, '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        ymean , expo = pretty_ve ( ymean , width = width , precision = precision , parentheses = False )
        row = 'y-mean' , ymean, '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        yrms , expo = pretty_ve ( yrms  , width = width , precision = precision , parentheses = False )
        row = 'y-rms' , yrms, '10^%+d' % expo if expo else '' 
        rows.append ( row )

        xyc , expo = pretty_float ( h2.xycorr() , width = width , precision = precision )
        row = 'xy-corr' , xyc , '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        stat = h2.stat()
        
        zmin , expo = pretty_float ( stat.min () , width = width , precision = precision )    
        row = 'z-min' , zmin , '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        zmax , expo = pretty_float ( stat.max () , width = width , precision = precision )    
        row = 'z-max' , zmax , '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        zmean , expo = pretty_ve ( stat.mean  () , width = width , precision = precision , parentheses = False )
        row = 'z-mean' , zmean , '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        zrms , expo = pretty_float ( stat.rms () , width = width , precision = precision )    
        row = 'z-rms' , zrms , '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
    has3 = False
    for row in rows [ 1 : ] :
        if row [ 2 ] :
            has3 = True
            break
        
    if not has3 :
        new_rows = []
        for row in rows : new_rows.append ( row [ : -1 ] )
        rows = new_rows 
                
    title = title if title else 'Summary of histogram %s(%s)' % ( typename ( h2 ) , h2.GetName() ) 
    
    table = T.table ( rows ,title = title , prefix = prefix , alignment = 'lc' ) 
    return table 


# =============================================================================
## print histogram summary as table
#  @code
#  histo = ...
#  print ( histo.table() )   
#  @endcode
def _h3_table_ ( h3 , title = '' , prefix = '' , width = 5 , precision = 3 ) :
    """ Print histogram summary as table 
    >>> histo = ...
    >>> print ( histo.table() )   
    """
    from ostap.logger.pretty import pretty_float , pretty_ve 
    rows = [ ( '' , 'value' , 'scale' ) ]

    row = 'Name' , h3.GetName() , '' 
    rows.append ( row )

    htitle = h3.GetTitle ()
    if htitle : 
        row = 'title' , htitle , ''
        rows.append ( row )
    
    rdir  = h3.GetDirectory()
    if rdir : 
        row = 'directory' , rdir.GetName() , '' 
        rows.append ( row )
    groot = ROOT.ROOT.GetROOT()
    row   = 'gROOT?' , 'True' if groot is rdir else 'False' , '' 
    rows.append ( row )
                
    row = 'SumW2' , 'True' if h3.GetSumw2() else 'False' , '' 
    rows.append ( row )

    xmin , expo = pretty_float ( h3.xmin () , width = width , precision = precision )    
    row = 'x-min' , xmin , '10^%+d' % expo if expo else '' 
    rows.append ( row )
    
    xmax , expo = pretty_float ( h3.xmax () , width = width , precision = precision )
    row = 'x-max' , xmax , '10^%+d' % expo if expo else '' 
    rows.append ( row )

    ymin , expo = pretty_float ( h3.ymin () , width = width , precision = precision )    
    row = 'y-min' , ymin , '10^%+d' % expo if expo else '' 
    rows.append ( row )
    
    ymax , expo = pretty_float ( h3.ymax () , width = width , precision = precision )
    row = 'y-max' , ymax , '10^%+d' % expo if expo else '' 
    rows.append ( row )

    zmin , expo = pretty_float ( h3.zmin () , width = width , precision = precision )    
    row = 'z-min' , zmin , '10^%+d' % expo if expo else '' 
    rows.append ( row )
    
    zmax , expo = pretty_float ( h3.zmax () , width = width , precision = precision )
    row = 'z-max' , zmax , '10^%+d' % expo if expo else '' 
    rows.append ( row )


    xaxis = h3.GetXaxis() 
    row = '#xbins' , '%d/%s' %  ( xaxis.GetNbins()  , 'uniform' if xaxis.uniform() else 'no-uniform' ) , '' 
    rows.append ( row )
    yaxis = h3.GetYaxis() 
    row = '#ybins' , '%d/%s' %  ( yaxis.GetNbins()  , 'uniform' if yaxis.uniform() else 'no-uniform' ) , '' 
    rows.append ( row )
    zaxis = h3.GetZaxis() 
    row = '#zbins' , '%d/%s' %  ( zaxis.GetNbins()  , 'uniform' if zaxis.uniform() else 'no-uniform' ) , '' 
    rows.append ( row )
    
    row = '#entries' , '%d' % h3.GetEntries()  , '' 
    rows.append ( row )

    integral  = h3.Integral()
    nEff      = h3.GetEffectiveEntries()
    
    if integral != h3.GetEntries () or integral != nEff : 
        integral , expo  = pretty_float ( integral , width = width , precision = precision ) 
        row = 'in-range' , integral ,  '10^%+d' % expo if expo else '' 
        rows.append ( row )

    if nEff != h3.GetEntries() or nEff != h3.Integral() :
        nEff , expo  = pretty_float ( nEff , width = width , precision = precision ) 
        row = '#eff' , nEff ,  '10^%+d' % expo if expo else '' 
        rows.append ( row )

    row = 'natural?'   , '%s' % h3.natural () , '' 
    rows.append ( row )

    row = 'all ints?' , '%s' % h3.allInts ()  , '' 
    rows.append ( row )

    row = 'all >0/>=0/<0/<=0/=0?', '%s/%s/%s/%s/%s' % ( 'True' if h3.all_positive    () else 'False' ,
                                                        'True' if h3.all_nonnegative () else 'False' ,
                                                        'True' if h3.all_negative    () else 'False' ,
                                                        'True' if h3.all_nonpositive () else 'False' ,
                                                        'True' if h3.all_zero        () else 'False' ) , '' 
    rows.append ( row )

    rsum , expo  = pretty_float ( h3.riemann_sum ()   , width = width , precision = precision ) 
    row = 'Riemann sum' , rsum  ,  '10^%+d' % expo if expo else '' 
    rows.append ( row )

    if h3.nEff() and not h3.all_zero () : 
        
        row = 'density?'   , '%s' % h3.isDensity() , '' 
        rows.append ( row )
        
        xmean, ymean, zmean = h3.mean () 
        xrms , yrms , zrms  = h3.rms  () 
        
        xmean , expo = pretty_ve ( xmean , width = width , precision = precision , parentheses = False )
        row = 'x-mean' , xmean, '10^%+d' % expo if expo else '' 
        rows.append ( row )

        xrms , expo = pretty_ve ( xrms  , width = width , precision = precision , parentheses = False )
        row = 'x-rms' , xrms, '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        ymean , expo = pretty_ve ( ymean , width = width , precision = precision , parentheses = False )
        row = 'y-mean' , ymean, '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        yrms , expo = pretty_ve ( yrms  , width = width , precision = precision , parentheses = False )
        row = 'y-rms' , yrms, '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        zmean , expo = pretty_ve ( zmean , width = width , precision = precision , parentheses = False )
        row = 'z-mean' , zmean, '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        zrms , expo = pretty_ve ( zrms  , width = width , precision = precision , parentheses = False )
        row = 'z-rms' , zrms, '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        xyc , expo = pretty_float ( h3.xycorr() , width = width , precision = precision )
        row = 'xy-corr' , xyc , '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        xzc , expo = pretty_float ( h3.xzcorr() , width = width , precision = precision )
        row = 'xz-corr' , xzc , '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        yzc , expo = pretty_float ( h3.yzcorr() , width = width , precision = precision )
        row = 'yz-corr' , yzc , '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        stat = h3.stat()
        
        vmin , expo = pretty_float ( stat.min () , width = width , precision = precision )    
        row = 'v-min' , vmin , '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        vmax , expo = pretty_float ( stat.max () , width = width , precision = precision )    
        row = 'v-max' , vmax , '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        vmean , expo = pretty_ve ( stat.mean  () , width = width , precision = precision , parentheses = False )
        row = 'v-mean' , vmean , '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
        vrms , expo = pretty_float ( stat.rms () , width = width , precision = precision )    
        row = 'v-rms' , vrms , '10^%+d' % expo if expo else '' 
        rows.append ( row )
        
    has3 = False
    for row in rows [ 1 : ] :
        if row [ 2 ] :
            has3 = True
            break
        
    if not has3 :
        new_rows = []
        for row in rows : new_rows.append ( row [ : -1 ] )
        rows = new_rows 
                
    title = title if title else 'Summary of histogram %s(%s)' % ( typename ( h3 ) , h3.GetName() ) 
    
    table = T.table ( rows ,title = title , prefix = prefix , alignment = 'lc' ) 
    return table 


ROOT.TH1F.table   = _h1_table_
ROOT.TH1D.table   = _h1_table_
ROOT.TH1F.summary = _h1_table_
ROOT.TH1D.summary = _h1_table_

ROOT.TH2F.table   = _h2_table_
ROOT.TH2D.table   = _h2_table_
ROOT.TH2F.summary = _h2_table_
ROOT.TH2D.summary = _h2_table_

ROOT.TH3F.table    = _h3_table_
ROOT.TH3D.table    = _h3_table_
ROOT.TH3F.summary  = _h3_table_
ROOT.TH3D.summary  = _h3_table_


ROOT.TH1F.__repr__ = _h1_table_
ROOT.TH1D.__repr__ = _h1_table_
ROOT.TH2F.__repr__ = _h2_table_
ROOT.TH2D.__repr__ = _h2_table_
ROOT.TH3F.__repr__ = _h3_table_
ROOT.TH3D.__repr__ = _h3_table_


# =============================================================================
_decorated_classes_ = (
    ROOT.TH1   ,
    #
    ROOT.TH1F  ,
    ROOT.TH1D  ,
    #
    ROOT.TH2F  ,
    ROOT.TH2D  ,
    #
    ROOT.TH3F  ,
    ROOT.TH3D  ,
    #
    ROOT.TProfile ,
    #
    ROOT.TAxis , 
    ROOT.TF1   , 
    )

_new_methods_  = tuple ( _new_methods_ ) 
_new_methods_  += (
    #
    ROOT.TH1F.__init__ , 
    ROOT.TH2F.__init__ , 
    ROOT.TH3F.__init__ ,
    #
    ROOT.TH3F.Clone    ,
    ROOT.TH1F.Clone    ,
    ROOT.TH2F.Clone    ,
    #
    ROOT.TH3F.clone    ,
    ROOT.TH1F.clone    ,
    ROOT.TH2F.clone    ,
    #
    ROOT.TAxis . __iter__     ,
    ROOT.TAxis . __reversed__ ,
    ROOT.TAxis . __contains__ ,
    ROOT.TAxis .   scale      ,
    #
    ROOT.TAxis.uniform        , 
    ROOT.TAxis.uniform_bins   , 
    #
    ROOT.TH1   . same_dims    , 
    ROOT.TH1   . same_bins    , 
    #
    ROOT.TH1F. __setitem__    ,
    ROOT.TH1D. __setitem__    ,
    #
    ROOT.TH2F. __setitem__    ,
    ROOT.TH2D. __setitem__    ,
    #
    ROOT.TH3F. __setitem__    ,
    ROOT.TH3D. __setitem__    ,
    #
    ROOT.TH1F  . __getitem__  ,
    ROOT.TH1D  . __getitem__  ,
    ROOT.TH2F  . __getitem__  ,
    ROOT.TH2D  . __getitem__  ,
    ROOT.TH3F  . __getitem__  ,
    ROOT.TH3D  . __getitem__  ,
    #
    ROOT.TH1  . __iter__      ,
    ROOT.TH1F . __iter__      ,
    ROOT.TH1D . __iter__      ,
    ROOT.TH1F . __reversed__  ,
    ROOT.TH1D . __reversed__  ,
    #
    ROOT.TH2  . __iter__      , 
    ROOT.TH2F . __iter__      ,
    ROOT.TH2D . __iter__      ,
    ROOT.TH2F . __reversed__  ,
    ROOT.TH2D . __reversed__  ,
    #
    ROOT.TH3  . __iter__      ,
    ROOT.TH3F . __iter__      ,
    ROOT.TH3D . __iter__      ,
    ROOT.TH3F . __reversed__  ,
    ROOT.TH3D . __reversed__  ,
    #
    ROOT.TH1F  . __call__     ,
    ROOT.TH1D  . __call__     ,
    #
    ROOT.TH1   . __len__      ,
    ROOT.TH1   .   size       ,
    ROOT.TH1   . __contains__ ,
    #
    ROOT.TH1F  . __contains__ ,
    ROOT.TH1D  . __contains__ ,
    #
    ROOT.TH2   . __len__      ,
    ROOT.TH2   .   size       ,
    #
    ROOT.TH3   . __len__      ,
    ROOT.TH3   .   size       ,
    #
    ROOT.TH1D  . nbins        ,
    ROOT.TH1F  . nbins        ,
    ROOT.TH1D  .  bins        ,
    ROOT.TH1F  .  bins        ,
    #
    ROOT.TH2D  . nbinsx       ,
    ROOT.TH2D  . nbinsy       ,
    ROOT.TH2F  . nbinsx       ,
    ROOT.TH2F  . nbinsy       ,
    ROOT.TH2D  .  binsx       ,
    ROOT.TH2D  .  binsy       ,
    ROOT.TH2F  .  binsx       ,
    ROOT.TH2F  .  binsy       ,
    #
    ROOT.TH3D  . nbinsx       ,
    ROOT.TH3D  . nbinsy       ,
    ROOT.TH3D  . nbinsz       ,
    ROOT.TH3F  . nbinsx       ,
    ROOT.TH3F  . nbinsy       ,
    ROOT.TH3F  . nbinsz       ,
    ROOT.TH3D  .  binsx       ,
    ROOT.TH3D  .  binsy       ,
    ROOT.TH3D  .  binsz       ,
    ROOT.TH3F  .  binsx       ,
    ROOT.TH3F  .  binsy       ,
    ROOT.TH3F  .  binsz       ,
    #
    ROOT.TH3   . __contains__ ,
    ROOT.TH3F  . __contains__ ,
    ROOT.TH3D  . __contains__ ,
    #
    ROOT.TH1 . numEmpty       ,
    #
    ROOT.TH1F . findBin       ,
    ROOT.TH1D . findBin       ,
    ROOT.TH2F . findBin       ,
    ROOT.TH2D . findBin       ,
    ROOT.TH3F . findBin       ,
    ROOT.TH3D . findBin       ,
    #
    ROOT.TH2F . mean          , 
    ROOT.TH2D . mean          ,
    #
    ROOT.TH2   . __call__     ,
    ROOT.TH2F  . __getitem__  ,
    ROOT.TH2D  . __getitem__  ,
    #
    ROOT.TH3   . __call__     ,
    ROOT.TH3F  . __getitem__  ,
    ROOT.TH3D  . __getitem__  ,
    #
    ROOT.TH1F  . iteritems    ,
    ROOT.TH1D  . iteritems    ,
    #
    ROOT.TH2F  . iteritems    ,
    ROOT.TH2D  . iteritems    ,
    #
    ROOT.TH3F  . iteritems    ,
    ROOT.TH3D  . iteritems    ,
    #
    ROOT.TH1F  .     items    ,
    ROOT.TH1D  .     items    ,
    #
    ROOT.TH2F  .     items    ,
    ROOT.TH2D  .     items    ,
    #
    ROOT.TH3F  .     items    ,
    ROOT.TH3D  .     items    ,
    #
    ROOT.TAxis . iteritems    ,
    ROOT.TAxis .     items    ,
    ROOT.TAxis . __eq__       ,
    ROOT.TAxis . __ne__       ,
    #
    ROOT.TH1F.bin , 
    ROOT.TH1D.bin ,
    ROOT.TH2F.bin ,
    ROOT.TH2D.bin ,
    ROOT.TH3F.bin ,
    ROOT.TH3D.bin ,
    #
    ROOT.TH2F.T         ,
    ROOT.TH2D.T         ,
    ROOT.TH2F.Transpose ,
    ROOT.TH2D.Transpone ,
    ROOT.TH2F.transpose ,
    ROOT.TH2D.transpone ,
    #
    ROOT.TH2F . colz    ,
    ROOT.TH2D . colz    ,
    ROOT.TH2F . Colz    ,
    ROOT.TH2D . Colz    ,
    #
    ROOT.TH2F . texte   ,
    ROOT.TH2D . texte   ,
    #
    ROOT.TH1F.       binomEff   ,
    ROOT.TH1D.       binomEff   ,
    ROOT.TH1F.       wilsonEff  ,
    ROOT.TH1D.       wilsonEff  ,
    ROOT.TH1F. agrestiCoullEff  ,
    ROOT.TH1D. agrestiCoullEff  ,
    #
    ROOT.TH1 .       binomEff_2 ,
    #
    ROOT.TH1F.eff_wald                   ,
    ROOT.TH1D.eff_wald                   , 
    ROOT.TH1F.eff_wilson_score           , 
    ROOT.TH1D.eff_wilson_score           , 
    ROOT.TH1F.eff_wilson_score_continuity ,
    ROOT.TH1D.eff_wilson_score_continuity ,
    ROOT.TH1F.eff_arcsin                  ,
    ROOT.TH1D.eff_arcsin                  ,
    ROOT.TH1F.eff_agresti_coull           ,
    ROOT.TH1D.eff_agresti_coull           ,
    ROOT.TH1F.eff_jeffreys                ,
    ROOT.TH1D.eff_jeffreys                ,
    ROOT.TH1F.eff_clopper_pearson         ,
    ROOT.TH1D.eff_clopper_pearson         ,
    #
    ROOT.TH2F.  binomEff    ,
    ROOT.TH2D.  binomEff    ,
    #
    ROOT.TH3F.  binomEff    ,
    ROOT.TH3D.  binomEff    ,
    #
    ROOT.TH1F . __floordiv__  ,
    ROOT.TH1D . __floordiv__  ,
    ROOT.TH2F . __floordiv__  ,
    ROOT.TH2D . __floordiv__  ,
    ROOT.TH3F . __floordiv__  ,
    ROOT.TH3D . __floordiv__  ,
    ##
    ROOT.TH1F .  zechEff ,
    ROOT.TH1D .  zechEff ,
    ROOT.TH2F .  zechEff ,
    ROOT.TH2D .  zechEff ,
    ROOT.TH3F .  zechEff ,
    ROOT.TH3D .  zechEff ,
    ##
    ROOT.TH1F . __mod__  ,
    ROOT.TH1D . __mod__  ,
    ROOT.TH2F . __mod__  ,
    ROOT.TH2D . __mod__  ,
    ROOT.TH3F . __mod__  ,
    ROOT.TH3D . __mod__  ,
    ##
    _h1_oper_    ,
    _h1_ioper_   ,
    _h1_div_     ,
    _h1_mul_     ,
    _h1_add_     ,
    _h1_sub_     ,
    _h1_pow_     ,
    _h1_idiv_    ,
    _h1_idiv_    ,
    _h1_imul_    ,
    _h1_iadd_    ,
    _h1_isub_    ,    
    _h1_rdiv_    ,   
    _h1_rmul_    ,
    _h1_radd_    ,
    _h1_rsub_    ,
    _h1_rrshift_ , 
    _h1_abs_     , 
    _h1_frac_    , 
    _h1_asym_    , 
    _h1_diff_    , 
    _h1_chi2_    , 
    _h1_mean_    ,
    #
    ROOT.TH1F . find_X ,
    ROOT.TH1D . find_X , 
    #
    ROOT.TH1F . cl_interval , 
    ROOT.TH1D . cl_interval , 
    #
    ROOT.TH1 . minv    ,
    ROOT.TH1 . maxv    ,
    ROOT.TH1 . minmax  ,
    #
    ROOT.TH1D. xmin    ,
    ROOT.TH1D. xmax    ,
    ROOT.TH1D. ymin    ,
    ROOT.TH1D. ymax    ,
    #
    ROOT.TH2 . xmin    ,
    ROOT.TH2 . xmax    ,
    ROOT.TH2 . ymin    ,
    ROOT.TH2 . ymax    ,
    ROOT.TH2 . zmin    ,
    ROOT.TH2 . zmax    ,
    #
    ROOT.TH3 . xmin    ,
    ROOT.TH3 . xmax    ,
    ROOT.TH3 . ymin    ,
    ROOT.TH3 . ymax    ,
    ROOT.TH3 . zmin    ,
    ROOT.TH3 . zmax    ,
    #
    ROOT.TH1D. xminmax ,
    ROOT.TH1D. yminmax ,
    ROOT.TH1F. xminmax ,
    ROOT.TH1F. yminmax ,
    #
    ROOT.TH2 . xminmax ,
    ROOT.TH2 . yminmax ,
    ROOT.TH2 . zminmax ,
    #
    ROOT.TH3 . xminmax ,
    ROOT.TH3 . yminmax ,
    ROOT.TH3 . zminmax ,
    #
    ROOT.TH1 . min_positive  ,
    ROOT.TH1 . minpos        ,
    #
    ROOT.TH2F.random   ,
    ROOT.TH2D.random   ,
    #
    ROOT.TH3F.random   ,
    ROOT.TH3D.random   ,
    #
    ROOT.TH1F.random   ,
    ROOT.TH1D.random   ,
    #
    ROOT.TH2F.shoot    ,
    ROOT.TH2D.shoot    ,
    #
    ROOT.TH3F.shoot    ,
    ROOT.TH3D.shoot    ,
    #
    ROOT.TH1F.shoot    ,
    ROOT.TH1D.shoot    ,
    #
    _h2_oper_    ,
    _h2_ioper_   ,
    _h2_div_     ,
    _h2_mul_     ,
    _h2_add_     ,
    _h2_sub_     ,
    _h2_pow_     ,
    _h2_idiv_    ,
    _h2_imul_    ,
    _h2_iadd_    ,
    _h2_isub_    ,    
    _h2_rdiv_    ,   
    _h2_rmul_    ,
    _h2_radd_    ,
    _h2_rsub_    ,
    _h2_abs_     , 
    _h2_frac_    , 
    _h2_asym_    , 
    _h2_diff_    , 
    _h2_chi2_    , 
    _h2_mean_    ,
    #    
    _h2_box_     ,
    _h2_lego_    , 
    _h2_surf_    , 
    _h2_surf2_   ,
    #
    _h3_oper_    ,
    _h3_ioper_   ,
    _h3_div_     ,
    _h3_mul_     ,
    _h3_add_     ,
    _h3_sub_     ,
    _h3_pow_     ,
    _h3_idiv_    ,
    _h3_imul_    ,
    _h3_iadd_    ,
    _h3_isub_    ,    
    _h3_rdiv_    ,   
    _h3_rmul_    ,
    _h3_radd_    ,
    _h3_rsub_    ,
    _h3_abs_     , 
    _h3_frac_    , 
    _h3_asym_    , 
    _h3_chi2_    , 
    _h3_mean_    ,
    #    
    ROOT.TH1F.addFunctionIntegral ,
    ROOT.TH1D.addFunctionIntegral ,
    ###
    h1_sumv                ,
    ## 
    ROOT.TH1F.eff          ,  
    ROOT.TH1D.eff          ,
    ROOT.TH1F.effic        , 
    ROOT.TH1D.effic        ,
    ROOT.TH1F.efficiency   ,
    ROOT.TH1D.efficiency   ,
    ROOT.TH1F.significance ,
    ROOT.TH1D.significance ,
    #
    ROOT.TH1F. smear     ,
    ROOT.TH1D. smear     , 
    ROOT.TH1F. smooth    ,
    ROOT.TH1D. smooth    , 
    #
    ROOT.TH1F. transform ,
    ROOT.TH1D. transform , 
    #
    ROOT.TH2F. transform ,
    ROOT.TH2D. transform , 
    #
    ROOT.TH1F. precision ,
    ROOT.TH1D. precision ,
    ROOT.TH2F. precision ,
    ROOT.TH2D. precision ,
    ROOT.TH1F. b2s       ,
    ROOT.TH1D. b2s       ,
    ROOT.TH2F. b2s       ,
    ROOT.TH2D. b2s       ,
    #
    ROOT.TH1F. rescale_bins ,
    ROOT.TH1D. rescale_bins ,
    #
    ROOT.TH2F. rescale_bins ,
    ROOT.TH2D. rescale_bins ,
    #
    ROOT.TH1  . sample  ,
    ROOT.TH1  . gauss   ,
    ROOT.TH1  . poisson ,
    #
    ROOT.TH1D . fom_1 ,
    ROOT.TH1D . fom_2 ,
    ROOT.TH1F . fom_1 ,
    ROOT.TH1F . fom_2 ,
    #
    ROOT.TH1D . FoM_1 ,
    ROOT.TH1D . FoM_2 ,
    ROOT.TH1F . FoM_1 ,
    ROOT.TH1F . FoM_2 ,
    #
    ROOT.TH1F.rebinNumbers  , 
    ROOT.TH1D.rebinNumbers  , 
    ROOT.TH1F.rebinFunction , 
    ROOT.TH1D.rebinFunction ,
    ROOT.TH2F.rebinNumbers  , 
    ROOT.TH2D.rebinNumbers  , 
    ROOT.TH2F.rebinFunction , 
    ROOT.TH2D.rebinFunction ,
    #
    ROOT.TH1D.null   ,
    ROOT.TH1F.null   ,
    ROOT.TH1D.level  ,
    ROOT.TH1F.level  ,
    ROOT.TH1D.vline  ,
    ROOT.TH1F.vline  ,
    ROOT.TH1D.hline  ,
    ROOT.TH1F.hline  ,
    #
    ROOT.TAxis.edges ,
    #
    ROOT.TH1F  . __getslice__  ,
    ROOT.TH1D  . __getslice__  ,
    #
    ROOT.TH1F.histoDiff ,
    ROOT.TH1D.histoDiff ,
    ROOT.TH2F.histoDiff ,
    ROOT.TH2D.histoDiff ,
    ROOT.TH3F.histoDiff ,
    ROOT.TH3D.histoDiff ,
    #
    ROOT.TH1F .   shift     ,
    ROOT.TH1D .   shift     ,
    ROOT.TH1D . __rshift__  ,
    ROOT.TH1F . __rshift__  ,
    ROOT.TH1D . __lshift__  ,
    ROOT.TH1F . __lshift__  ,
    ROOT.TH1D . __ilshift__ ,
    ROOT.TH1F . __ilshift__ ,
    ROOT.TH1D . __irshift__ ,
    ROOT.TH1F . __irshift__ ,
    #
    ROOT.TH1D.  accumulate  , 
    ROOT.TH1D.  sum         , 
    ROOT.TH1D.  integrate   , 
    #
    ROOT.TH2D.  accumulate  , 
    ROOT.TH2D.  sum         , 
    ROOT.TH2D.  integrate   , 
    #
    ROOT.TH3D.  accumulate  , 
    ROOT.TH3D.  sum         , 
    ROOT.TH3D.  integrate   , 
    #
    ROOT.TH1D.mean          ,
    ROOT.TH1D.rms           ,
    ROOT.TH1D.skewness      ,
    ROOT.TH1D.kurtosis      ,
    ROOT.TH1D.moment        ,
    ROOT.TH1D.centralMoment ,
    ROOT.TH1D.nEff          ,
    # 
    ROOT.TH1.stat           ,
    ROOT.TH1.wstat          ,
    #
    ROOT.TH1F.xstat         ,
    ROOT.TH2F.xstat         ,
    ROOT.TH2F.ystat         ,
    ROOT.TH3F.xstat         ,
    ROOT.TH3F.ystat         ,
    ROOT.TH3F.zstat         ,
    #
    ROOT.TH1F.add_fake_bin       ,
    ROOT.TH1F.add_fake_bin_left  ,
    ROOT.TH1F.add_fake_bin_right ,
    #
    ROOT.TH2F.add_fake_side      ,
    ROOT.TH3F.add_fake_side      ,
    #
    ROOT.TH1F.solve              , 
    ROOT.TH1F.equal_edges        ,
    #
    ROOT.TH2F.slice              , 
    ROOT.TH2F.sliceX             , 
    ROOT.TH2F.sliceY             ,
    #
    ROOT.TH2F.islice_X           , 
    ROOT.TH2F.islice_Y           ,
    #
    ROOT.TH3F.slice              , 
    ROOT.TH3F.sliceX             , 
    ROOT.TH3F.sliceY             ,
    ROOT.TH3F.sliceZ             ,
    #
    ROOT.TH3F.islice_X           , 
    ROOT.TH3F.islice_Y           ,
    ROOT.TH3F.islice_Z           ,
    #
    ROOT.TH2F . proj   ,
    ROOT.TH2F . projX  ,
    ROOT.TH2F . projY  ,
    #
    ROOT.TH2D . proj   ,
    ROOT.TH2D . projX  ,
    ROOT.TH2D . projY  ,
    #
    ROOT. TH3 . projX  ,
    ROOT. TH3 . projY  ,
    ROOT. TH3 . projZ  ,
    ROOT. TH3 . projXY ,
    ROOT. TH3 . projYX ,
    ROOT. TH3 . projXZ ,
    ROOT. TH3 . projZX ,
    ROOT. TH3 . projYZ ,
    ROOT. TH3 . projZY ,
    #
    ROOT.TProfile . asH1     ,
    ROOT.TProfile . toH1     ,
    ROOT.TProfile . asHisto  ,
    #
    ROOT.TF1 .  minv    ,
    ROOT.TF1 .  maxv    ,
    ROOT.TF1 .  minmax  ,
    ROOT.TF1 . xminmax  ,
    ROOT.TF1 . mean     ,
    ROOT.TF1 . variance ,
    ROOT.TF1 . rms      ,
    ROOT.TF1 . median   ,
    ROOT.TF1 . mode     ,
    ROOT.TF1 . moment   ,
    ROOT.TF1 . quantile ,
    #
    ROOT.TAxis. __getslice__ ,
    #
    ROOT.TH1F . histoGuess   ,
    ROOT.TH1D . histoGuess   ,
    ROOT.TH1  . useLL        ,
    ROOT.TH1  . allInts      ,
    ROOT.TH1  . natural      ,
    #
    ROOT.TH1  . dim          ,
    ROOT.TH1  . n_bins       ,
    #
    ROOT.TH1.uniform_bins    ,
    ROOT.TH1.uniform         ,
    #
    ROOT.TH1F.dump           ,
    ROOT.TH1F.dumpHisto      ,
    ROOT.TH1F.dumpASCII      ,
    ROOT.TH1F.dumpAsText     ,
    #
    ROOT.TH1F.transform_X_numbers  ,
    ROOT.TH1F.transform_X_function ,
    ##
    ROOT.TH1   . fill_loop   ,
    ## 
    ROOT.TH1F  . split       ,
    ROOT.TH1F  . split_bins  ,
    ##
    ROOT.TH1   . nEntries    , 
    ## 
    ROOT.TH1   . minimum_bin , 
    ROOT.TH1   . maximum_bin , 
    ROOT.TH2   . minimum_bin , 
    ROOT.TH2   . maximum_bin , 
    ROOT.TH3   . minimum_bin , 
    ROOT.TH3   . maximum_bin ,
    ##
    ROOT.TH1F.the_moment     , 
    ROOT.TH1D.the_moment     , 
    ##
    ROOT.TH1F. scale_axis    , 
    ROOT.TH1D. scale_axis    , 
    ##
    ROOT.TH2F. scale_axes    ,
    ROOT.TH2D. scale_axes    , 
    ##
    ROOT.TH3F. scale_axes    , 
    ROOT.TH3D. scale_axes    ,
    ##
    ROOT.TH1F. the_integral  , 
    ROOT.TH1D. the_integral  ,
    ##
    ROOT.TH2F. the_integral  , 
    ROOT.TH2D. the_integral  ,
    ##
    ROOT.TH3F. the_integral  , 
    ROOT.TH3D. the_integral  ,
    ##
    histo_book               ,
    ##
    ROOT.TH1F.table          , 
    ROOT.TH1D.table          , 
    ROOT.TH1F.summary        , 
    ROOT.TH1D.summary        ,
    ##
    ROOT.TH2F.table          , 
    ROOT.TH2D.table          , 
    ROOT.TH2F.summary        , 
    ROOT.TH2D.summary        ,
    ##
    ROOT.TH3F.table          , 
    ROOT.TH3D.table          , 
    ROOT.TH3F.summary        , 
    ROOT.TH3D.summary        ,
    ##
    )

# =============================================================================
if '__main__' == __name__ :
    
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger ) 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
