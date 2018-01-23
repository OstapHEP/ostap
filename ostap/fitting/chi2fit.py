#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  @see Ostap::Math::Chi2Fit
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-09-12
# =============================================================================
"""Simple chi2-fit
- see TH1(F/D).hFit 
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@nikhef.nl"
__date__    = "2009-09-12"
__version__ = ""
# =============================================================================
__all__     = (
    'C2FIT' , ## simple chi2-fit 
    ) 
# =============================================================================
import ROOT, cppyy
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.chi2fit' )
else                       : logger = getLogger ( __name__                )
# =============================================================================
from ostap.math.ve import VE,cpp 
cpp   = cppyy.gbl
Ostap = cpp.Ostap
C2FIT = Ostap.Math.Chi2Fit

C2FIT . __str__  = lambda s : s.toString ()
C2FIT . __repr__ = lambda s : s.toString ()
# =============================================================================
## chi2-probabilty
def _c2_prob_  ( s ) :
    """Chi2 probabiilty
    >>> r = h.hfit ( ... )
    >>> r.Prob()
    """
    dofs = s.points() - s.size()
    return ROOT.TMath.Prob ( s.chi2() , dofs )

C2FIT . Prob        = _c2_prob_
C2FIT . Probability = _c2_prob_
C2FIT . prob        = _c2_prob_
C2FIT . probability = _c2_prob_
C2FIT . __len__     = lambda s     : s.size  (   )
C2FIT . __getitem__ = lambda s , i : s.param ( i )

# =============================================================================
## (Chi2)fit the histogram by sum of components
#  The ``components'' could be histograms, functions and other
#  callable object :
#  @code
#
#  h0 = ...
#  h .hFit ( h0 )
# 
#  h0 = ...
#  h1 = ...
#  h .hFit ( [ h0 , h1 ] )
#  @endcode 
def _h_Fit_ ( self                              ,
              components                        ,
              draw = False                      ,
              interpolate = True                ,
              selector    = lambda i,x,y : True ) :
    """(Chi_2)-fit the histogram with the set of ``components''
    
    The ``components'' could be histograms, functions and other
    callable object :
    
    >>> h0 = ...
    >>> h .hFit ( h0 )
    
    >>> h0 = ...
    >>> h1 = ...
    >>> h .hFit ( [ h0 , h1 ] )
    
    """
    DATA =   VE.Vector
    CMPS = DATA.Vector

    if   isinstance ( components , ROOT.TH1D ) : components = [ components ]
    elif isinstance ( components , ROOT.TH1F ) : components = [ components ]
    elif not isinstance ( components , ( tuple , list ) ) :
        components = [ components ]
        
    data = DATA ()
    cmps = CMPS ()
    
    while len ( cmps ) < len ( components )  :
        cmps.push_back( DATA() )

    for i,x,y in self.iteritems () :

        if not selector ( i , x , y ) : continue
        
        dp = VE ( y )
        data.push_back ( dp )
        
        for j in range ( 0 , len ( components ) ) :

            cmp = components[j]
            if isinstance ( cmp , ( ROOT.TH1F , ROOT.TH1D )) :
                cp = cmp ( x , interpolate   )
            else :
                cp =  VE ( cmp ( x.value() ) )  ## CALLABLE !!! 
            
            cmps[ j ].push_back ( cp ) 
            

    _c2Fit = Ostap.Math.Chi2Fit ( data , cmps )

    if draw :
        
        if not hasattr ( self , '_histos_' ) :
            self._histos_ = []
            
        self.Draw( 'e1' ) 
        for j in range ( 0 , len ( components ) ) :
            cmpj = components [j]
            if not isinstance ( cmpj , ( ROOT.TH1F , ROOT.TH1D ) ) : continue
            sc = _c2Fit[j].value ()
            nh = cmpj * sc 
            nh.Draw ( 'same' )
            self._histos_.append ( nh ) 

    return _c2Fit


ROOT.TH1F. hFit = _h_Fit_ 
ROOT.TH1D. hFit = _h_Fit_ 

# =============================================================================
if '__main__' == __name__  :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
