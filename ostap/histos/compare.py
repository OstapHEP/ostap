#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## Copyright (c) Ostap developpers.
# =============================================================================  
## @file
#  Module with utilities for specific comparison of histograms/functions/shapes 
#  @date   2014-05-10
# =============================================================================
"""Module with utilities for specific comparison of histograms/functions/shapes
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = () 
# =============================================================================
import ROOT, math 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.histos.compare' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Some specific comparison of histo-objects')
# =============================================================================
from   ostap.core.core        import hID , VE 
import ostap.histos.histos 
import ostap.histos.param
import ostap.fitting.param 
from   ostap.logger.colorized import allright  
# =============================================================================
## Can 1D-histogram can be considered as ``constant'' ?
#  @code
#  histo = ...
#  print 'Is constant? %s ' % histo.is_constant( prob = 0.01 )
#  @endcode 
def _h1_constant_ ( h1 , prob = 0.10 , opts = '0Q' , density = False ) :
    """Can  1D-histogram be considered as constant ?
    >>> histo = ...
    >>> print 'Is constant? %s ' % histo.is_constant( prob = 0.01 ) 
    """
    #
    if density :
        h1_ = h1.density() 
        cmp = _h1_constant_ ( h1_ , prob , opts , density = False )
        del h1_
        return cmp
    
    if not isinstance ( h1 , ( ROOT.TH1D , ROOT.TH1F ) ) : return False 
    #
    r  = h1.Fit ( 'pol0', 'S' + opts )
    if 0 != r.Status() : 
        logger.warning("Can't fit with constant function %s" % r.Status() )
        r  = h1.Fit ( 'pol0', 'S' + opts )
        if 0 != r.Status() : return False
    #
    return prob <= r.Prob()

ROOT.TH1D.is_constant = _h1_constant_
ROOT.TH1F.is_constant = _h1_constant_

# =============================================================================
## compare the 1D-histograms trying to fit one with other
def _h1_cmp_fit_ ( h1              ,
                   h2              ,
                   density = False ,  
                   opts    = '0Q'  ) :
    """Compare histograms by refit of the first with functions,
    extracted from the second one

    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo
    >>> r  = h1.cmp_fit ( h2 )
    >>> if r : print r.Prob()    
    """
    
    assert isinstance ( h1 , ROOT.TH1 ) and 1 == h1.dim () , \
           "cmp_dist: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 1 == h2.dim () , "cmp_dist: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )

    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1 
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2 
        cmp = _h1_cmp_fit_ ( h1_ , h2_ ,  density = False , opts = opts )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp

    from ostap.fitting.param import C1Fun 
    f2 = C1Fun  ( h2 , *h1.xminmax() ) 
    f2.release  ( 0     )
    f2.set      ( 0 , 1 )
    f2.fix      ( 1 , 0 )
    f2.fix      ( 2 , 1 )

    rf = f2.Fit ( h1 , 'S' + opts )    
    if 0 != rf.Status() :
        logger.error ( "Can't fit with function " % rf.Status() )
        return None

    return rf

ROOT.TH1D.cmp_fit = _h1_cmp_fit_
ROOT.TH1F.cmp_fit = _h1_cmp_fit_ 


# =============================================================================
## compare the 1D-histograms trying to fit one with another 
def _h1_cmp_pdf_ ( h1               ,
                   h2               ,
                   density  = False , 
                   draw     = False ,
                   silent   = True  ) :
    """Compare histograms by refit of the first with functions,
    extracted from the second one

    >>> h1      = ... ## the first histo
    >>> h2      = ... ## the second histo
    >>> r, chi2 = h1.cmp_pdf ( h2 )
    """
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2
        cmp = _h1_cmp_pdf_ ( h1_ , h2_ , density = False ,draw =  draw , silent = silent )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp
    
    from ostap.fitting.basic      import H1D_pdf, Fit1D
    from ostap.fitting.background import Bkg_pdf 

    pdf2   = H1D_pdf    ( '_H1', h2 , density   = True , silent = True )
    model2 = Fit1D ( signal = pdf2 , background = Bkg_pdf( '_B2' , mass = pdf2.mass , power = 0 , tau = 0 ) ) 
    model2.b.fix(0)

    from ostap.utils.utils import invisibleCanvas 
    with invisibleCanvas() : 
        r2 , f2 = model2.chi2fitTo ( h1  , silent = silent , draw = draw , density = False , sumw2 = True )

    return r2 , f2.chiSquare()

ROOT.TH1D.cmp_pdf = _h1_cmp_pdf_
ROOT.TH1F.cmp_pdf = _h1_cmp_pdf_ 

# =============================================================================
## compare the 1D-histograms by chi2
#  @code
#  h1 = ... ## the first histo
#  h2 = ... ## the second histo (or function or anything else) 
#  chi2ndf , probability  = h1.cmp_chi2 ( h2 )
#  @endcode 
def _h1_cmp_chi2_ ( h1              ,
                    h2              ,
                    density = False ) :
    """Compare histograms by chi2
    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo (or function or anything else) 
    >>> chi2ndf , probability  = h1.cmp_chi2 ( h2 )    
    """
    
    assert isinstance ( h1 , ROOT.TH1 ) and 1 == h1.dim () , \
           "cmp_dist: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 1 == h2.dim () , "cmp_dist: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )
        
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2
        cmp = _h1_cmp_chi2_ ( h1_ , h2_ , density = False )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp

    chi2  = 0.0
    ndf   = 0  
    for i , x , v1  in h1.items() :        
        v2    = h2 ( x.value() )
        chi2 += v1.chi2 ( v2 )
        ndf  += 1 

    c2ndf = chi2/ndf 
    return c2ndf, ROOT.TMath.Prob ( chi2 , ndf ) 


ROOT.TH1D.cmp_chi2 = _h1_cmp_chi2_
ROOT.TH1F.cmp_chi2 = _h1_cmp_chi2_ 

# =============================================================================
## compare the 2D-histograms by chi2
#  @code
#  h1 = ... ## the first histo
#  h2 = ... ## the second histo (or function or anything else) 
#  chi2ndf , probability  = h1.cmp_chi2 ( h2 )
#  @endcode 
def _h2_cmp_chi2_ ( h1              ,
                    h2              ,
                    density = False ) :
    """Compare 2D-histograms by chi2
    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo (or function or anything else) 
    >>> chi2ndf , probability  = h1.cmp_chi2 ( h2 )    
    """
    
    assert isinstance ( h1 , ROOT.TH2 ) and 2 == h1.dim () , \
           "cmp_dist: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 2 == h2.dim () , "cmp_dist: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )
        
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2
        cmp = _h2_cmp_chi2_ ( h1_ , h2_ , density = False )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp

    chi2  = 0.0
    ndf   = 0 
    for ix , iy , x , y , v1  in h1.items() :        
        v2    = h2 ( x.value () , y.value () )
        chi2 += v1.chi2 ( v2 )
        ndf  += 1 

    c2ndf = chi2/ndf 
    return c2ndf, ROOT.TMath.Prob ( chi2 , ndf ) 

ROOT.TH2F.cmp_chi2 = _h2_cmp_chi2_
ROOT.TH2D.cmp_chi2 = _h2_cmp_chi2_

# =============================================================================
## compare the 3D-histograms by chi2
#  @code
#  h1 = ... ## the first histo
#  h2 = ... ## the second histo (or function or anything else) 
#  chi2ndf , probability  = h1.cmp_chi2 ( h2 )
#  @endcode 
def _h3_cmp_chi2_ ( h1              ,
                    h2              ,
                    density = False ) :
    """Compare 3D-histograms by chi2
    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo (or function or anything else) 
    >>> chi2ndf , probability  = h1.cmp_chi2 ( h2 )    
    """
    
    assert isinstance ( h1 , ROOT.TH3 ) and 3 == h1.dim () , \
           "cmp_dist: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 3 == h2.dim () , "cmp_dist: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )
        
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2
        cmp = _h3_cmp_chi2_ ( h1_ , h2_ , density = False )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp

    chi2  = 0.0
    ndf   = 0   
    for ix , iy , iz , x , y , z , v1  in h1.items() :        
        v2    = h2 ( x.value () , y.value () , z.value () )
        chi2 += v1.chi2 ( v2 )
        ndf  += 1 

    c2ndf = chi2/ndf 
    return c2ndf, ROOT.TMath.Prob ( chi2 , ndf ) 


ROOT.TH3F.cmp_chi2 = _h3_cmp_chi2_
ROOT.TH3D.cmp_chi2 = _h3_cmp_chi2_

# =============================================================================
## Calculate chi2 for histogram and ``function''
def _h1_chi2_cmp_ ( h1                                    ,
                    func                                  ,
                    integral = False                      ,
                    select   = lambda x,y,v : True        ,
                    chi2     = lambda v1,v2 : v1.chi2(v2) ) :
    
    """Calculate chi2 for histogram and ``function''
    >>> h1   = ... ## the first histo
    >>> func = ... ## the the function 
    >>> chi2ndf,prob  = h1.chi2_cmp ( func , integral = False )
    """
    c2  = 0
    ndf = 0

    _func_  = lambda  x , xl , xr : func ( x )
    if   integral and hasattr ( func , 'integral' ) :
        _func_  = lambda x,xl,xr : func.integral ( xl , xr ) / ( xr - xl ) 
    elif integral and hasattr ( func , 'Integral' ) :  
        _func_  = lambda x,xl,xr : func.Integral ( xl , xr ) / ( xr - xl ) 
    elif integral :
        ## use numerical integration 
        from ostap.math.intergal import integral as _integral_
        _func_  = lambda x , xl , xr : _integral_ ( func , xl , xr ) / ( xr - xl )


    ## helper function
    def _chi2_ ( c , histo , func , accept , funchi2 )  :

        c2   = 0.0
        ndf  = 1

        for entry in histo.items() :
        
            x    = entry [ 1 ]
            y1   = entry [ 2 ]
            
            xv   = x.value()
            xe   = x.error()
            xl   = xv - xe
            xr   = xv + xe
            
            y2   = func ( x , xl , xr )    
            if not accept ( x, y1 , y2 ) : continue

            c2  += funchi2 ( y1 , c * y2 )
            ndf += 1

        return c2 , ndf 

    if not scale : 
        c2 , ndf = _chi2_ ( 1.0 , h1 , _func_ , select , chi2 )
        c2ndf = c2/ndf 
        return c2ndf, ROOT.TMath.Prob( c2 , ndf )
    
    fun = lambda c : _chi2_ ( 1.0 , h1 , _func_ , select , chi2 )[0]

    from ostap.math.minimize import minimize_scalar 
    r = minimize_scalar ( fun )

    c2 , ndf =  _chi2_ ( r.x , h1 , _func_ , select , chi2 )
    
    c2ndf = c2/ndf 
    return c2ndf, ROOT.TMath.Prob( c2 , ndf ) , r.x 
    

ROOT.TH1D.chi2_cmp = _h1_chi2_cmp_
ROOT.TH1F.chi2_cmp = _h1_chi2_cmp_ 

# =============================================================================
## compare the 1D-historgams (as functions)
#  - calculate the scalar product and get cos(theta) from it 
# \f$cos \theta = \frac{ f_1 \cdot f_2 } { \left|f_1\right|\left|f_2\right| }\f$
def _h1_cmp_costheta_ ( h1              ,
                        h2              ,
                        density = False ) :  
    """Compare the 1D-historgams (as functions)
    Calculate the scalar product and get ``cos(theta)'' from it
    
    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo
    >>> cos_theta  = h1.cmp_cos ( h2 )
    
    """
    assert isinstance ( h1 , ROOT.TH1 ) and 1 == h1.dim () , \
           "cmp_cos: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 1 == h2.dim () , "cmp_cos: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )
    
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2
        cmp = _h1_cmp_costheta_ ( h1_ , h2_ , density  = False )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp
        
    f1 = lambda x : float ( h1 ( x ) ) 
    f2 = lambda x : float ( h2 ( x ) ) 

    lims = h1.xminmax()
    
    params = lims [ 0 ] , lims [ 1 ]     
    
    from ostap.math.integral import integral as _integral_
    r1   = _integral_ ( lambda x : f1 ( x ) ** 2    , *params )
    r2   = _integral_ ( lambda x : f2 ( x ) ** 2    , *params )
    r12  = _integral_ ( lambda x : f1 ( x ) * f2(x) , *params )
 
    return r12 / ( r1 * r2 ) ** 0.5 

ROOT.TH1D.cmp_cos = _h1_cmp_costheta_
ROOT.TH1F.cmp_cos = _h1_cmp_costheta_ 


# =============================================================================
## compare the 2D-historgams (as functions)
#  - calculate the scalar product and get cos(theta) from it 
# \f$cos \theta = \frac{ f_1 \cdot f_2 } { \left|f_1\right|\left|f_2\right| }\f$
# @attention It could be rather slow    
def _h2_cmp_costheta_ ( h1              ,
                        h2              ,
                        density = False ) :  
    """Compare the 2D-historgams (as functions)
    Calculate the scalar product and get ``cos(theta)'' from it
    
    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo
    >>> cos_theta  = h1.cmp_cos ( h2 )
    - it could be rather slow    
    
    """
    assert isinstance ( h1 , ROOT.TH2 ) and 2 == h1.dim () , \
           "cmp_cos: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 2 == h2.dim () , "cmp_cos: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )
    
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2
        cmp = _h2_cmp_costheta_ ( h1_ , h2_ , density  = False )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp
        
    f1 = lambda x , y : float ( h1 ( x , y ) ) 
    f2 = lambda x , y : float ( h2 ( x , y ) )
    
    xlims  = h1.xminmax()
    ylims  = h1.yminmax()    
    params = xlims [ 0 ] , xlims [ 1 ] , ylims [ 0 ] , ylims [ 1 ] 
    
    from ostap.math.integral import integral2 as _integral2_
    r1   = _integral2_ ( lambda x , y : f1 ( x , y ) ** 2           , *params ) 
    r2   = _integral2_ ( lambda x , y : f2 ( x , y ) ** 2           , *params ) 
    r12  = _integral2_ ( lambda x , y : f1 ( x , y ) * f2 ( x , y ) , *params )
    
    return r12 / ( r1 * r2 ) ** 0.5 

ROOT.TH2F.cmp_cos = _h2_cmp_costheta_
ROOT.TH2D.cmp_cos = _h2_cmp_costheta_

# =============================================================================
## compare the 3D-historgams (as functions)
#  - calculate the scalar product and get cos(theta) from it 
# \f$cos \theta = \frac{ f_1 \cdot f_2 } { \left|f_1\right|\left|f_2\right| }\f$
# @attention It could be rather slow    
def _h3_cmp_costheta_ ( h1              ,
                        h2              ,
                        density = False ) :  
    """Compare the 3D-historgams (as functions)
    Calculate the scalar product and get ``cos(theta)'' from it
    
    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo
    >>> cos_theta  = h1.cmp_cos ( h2 )
    - it could be rather slow    
    """
    assert isinstance ( h1 , ROOT.TH3 ) and 3 == h1.dim () , \
           "cmp_cos: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 3 == h2.dim () , "cmp_cos: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )
    
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2
        cmp = _h3_cmp_costheta_ ( h1_ , h2_ , density  = False )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp
        
    f1 = lambda x , y , z  : float ( h1 ( x , y , z ) ) 
    f2 = lambda x , y , z  : float ( h2 ( x , y , z ) )
    
    xlims  = h1.xminmax()
    ylims  = h1.yminmax()    
    zlims  = h1.zminmax()
    params = xlims [ 0 ] , xlims [ 1 ] , ylims [ 0 ] , ylims [ 1 ] , zlims [ 0 ] , zlims [ 1 ] 
        
    from ostap.math.integral import integral3 as _integral3_
    r1   = _integral3_ ( lambda x , y , z : f1 ( x , y , z ) ** 2               , *params )
    r2   = _integral3_ ( lambda x , y , z : f2 ( x , y , z ) ** 2               , *params )
    r12  = _integral3_ ( lambda x , y , z : f1 ( x , y , z ) * f2 ( x , y , z ) , *params ) 
    
    return r12 / ( r1 * r2 ) ** 0.5 

ROOT.TH3F.cmp_cos = _h3_cmp_costheta_
ROOT.TH3D.cmp_cos = _h3_cmp_costheta_

# =============================================================================
## compare the 1D-histograms 
#  - calculate "discrete" scalar products and take cos(theta) 
# \f$cos \theta = \frac{ f_1 \cdot f_2 } { \left|f_1\right|\left|f_2\right| }\f$
def _h1_cmp_dcostheta_ ( h1              ,
                         h2              ,
                         density = False ) :  
    """Compare the 1D-historgams (as functions)
    Calculate discrete scalar product and get ``the angle'' from it
    
    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo
    >>> cos_theta  = h1.cmp_dcos ( h2 )
    
    """
    assert isinstance ( h1 , ROOT.TH1 ) and 1 == h1.dim () , \
           "cmp_dcos: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 1 == h2.dim () , "cmp_dcos: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )

    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2 
        cmp = _h1_cmp_dcostheta_ ( h1_ , h2_ , density  = False )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp

    r1 , r2 , r12 = 0.0 , 0.0 , VE () 
    for i ,  x , v1  in h1.items () :
        xv   = x.value ()
        xe   = x.error () 
        v2   = VE ( h2 ( xv ) )
        r1  += 2 * xe * ( float ( v1 ) ** 2 )
        r2  += 2 * xe * ( float ( v2 ) ** 2 ) 
        r12 += 2 * xe * ( v1 * v2 ) 

    return r12 / ( r1 * r2 ) ** 0.5 

ROOT.TH1D.cmp_dcos = _h1_cmp_dcostheta_
ROOT.TH1F.cmp_dcos = _h1_cmp_dcostheta_ 

# =============================================================================
## compare the 2D-histograms 
#  - calculate "discrete" scalar products and take cos(theta) 
# \f$cos \theta = \frac{ f_1 \cdot f_2 } { \left|f_1\right|\left|f_2\right| }\f$
def _h2_cmp_dcostheta_ ( h1              ,
                         h2              ,
                         density = False ) :  
    """Compare the 2D-histograms (as functions)
    Calculate discrete scalar product and get ``the angle'' from it
    
    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo
    >>> cos_theta  = h1.cmp_dcos ( h2 )
    
    """
    assert isinstance ( h1 , ROOT.TH2 ) and 2 == h1.dim () , \
           "cmp_dcos: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 2 == h2.dim () , "cmp_dcos: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )

    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2 
        cmp = _h2_cmp_dcostheta_ ( h1_ , h2_ , density  = False )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp

    r1 , r2 , r12 = 0.0 , 0.0 , VE () 
    for ix , iy , x , y , v1 in h1.items () :
        xv   = x.value ()
        xe   = x.error () 
        yv   = y.value ()
        ye   = y.error () 
        v2   = VE ( h2 ( xv , yv ) )
        r1  += 4 * xe * ye * ( float ( v1 ) ** 2 )
        r2  += 4 * xe * ye * ( float ( v2 ) ** 2 ) 
        r12 += 4 * xe * ye * ( v1 * v2 ) 

    return r12 / ( r1 * r2 ) ** 0.5 

ROOT.TH2F.cmp_dcos = _h2_cmp_dcostheta_
ROOT.TH2D.cmp_dcos = _h2_cmp_dcostheta_

# =============================================================================
## compare the 3D-histograms 
#  - calculate "discrete" scalar products and take cos(theta) 
# \f$cos \theta = \frac{ f_1 \cdot f_2 } { \left|f_1\right|\left|f_2\right| }\f$
def _h3_cmp_dcostheta_ ( h1              ,
                         h2              ,
                         density = False ) :  
    """Compare the 3D-histograms (as functions)
    Calculate discrete scalar product and get ``the angle'' from it
    
    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo
    >>> cos_theta  = h1.cmp_dcos ( h2 )
    
    """
    assert isinstance ( h1 , ROOT.TH3 ) and 3 == h1.dim () , \
           "cmp_dcos: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 3 == h2.dim () , "cmp_dcos: invalid type of h2 %s/%s" % ( h2 , type ( h2 ) )

    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2 
        cmp = _h3_cmp_dcostheta_ ( h1_ , h2_ , density  = False )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp

    r1 , r2 , r12 = 0.0 , 0.0 , VE () 
    for ix , iy , iz , x , y , z , v1  in h1.items () :
        xv   = x.value ()
        xe   = x.error () 
        yv   = y.value ()
        ye   = y.error () 
        zv   = z.value ()
        ze   = z.error () 
        v2   = VE ( h2 ( xv , yv , zv ) )
        r1  += 8 * xe * ye * ze * ( float ( v1 ) ** 2 )
        r2  += 8 * xe * ye * ze * ( float ( v2 ) ** 2 ) 
        r12 += 8 * xe * ye * ze * ( v1 * v2 ) 

    return r12 / ( r1 * r2 ) ** 0.5 

ROOT.TH3F.cmp_dcos = _h3_cmp_dcostheta_
ROOT.TH3D.cmp_dcos = _h3_cmp_dcostheta_


# =============================================================================
## calculate the "distance" for two histograms or function-like objects
# Distance is defined as  
#  \f$ d = \left| f_1^{*} - f_2^{*} \right|^{1/2} \f$,
#  where:
#  - \f$ f^* \f$-are scaled functions, such \f$ \left| f^*\right| = 1 \f$
#  -  and the norm is   defined as
#   \f$ |left| f(x) \right| \ equiv = \frac{1}{b-a}\int^{b}_{a}  f^2 ( s ) dx \f$ 
def _h1_cmp_dist_ ( h1              ,
                    h2              ,
                    density = False ) : 
    """Calculate the norm of difference of scaled histograms/functions 
    |f1-f2|, such |f1|=1 and |f2|=1
    
    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo
    >>> diff = h1.cmp_dist ( h2 )
    
    """
    assert isinstance ( h1 , ROOT.TH1 ) and 1 == h1.dim () , \
           "cmp_dist: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 1 == h2.dim () , "cmp_dist: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )
    
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1 
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2
        cmp = _h1_cmp_dist_ ( h1_ , h2_ , density = False )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp

    f1   = lambda x : float ( h1 ( x ) )
    f2   = lambda x : float ( h2 ( x ) )

    lims   = h1.xminmax ()
    volume = lims [ 1 ] - lims [ 0 ]
    params = lims [ 0 ] , lims [ 1 ]
    
    from ostap.math.integral import integral as _integral_    
    r1 = _integral_ ( lambda x : f1 ( x ) ** 2 , *params ) / volume 
    r2 = _integral_ ( lambda x : f2 ( x ) ** 2 , *params ) / volume 
        
    sf1  = 1.0 / r1 ** 0.5 
    sf2  = 1.0 / r2 ** 0.5 

    df   = lambda x : ( sf1 * f1 ( x ) - sf2 * f2 ( x ) ) ** 2 
    d12  = _integral_ ( df  , *params ) / volume 
    
    return d12 ** 0.5 

ROOT.TH1D.cmp_dist = _h1_cmp_dist_
ROOT.TH1F.cmp_dist = _h1_cmp_dist_ 


# =============================================================================
## calculate the "distance" for two histograms or function-like objects
# Distance is defined as  
#  \f$ d = \left| f_1^{*} - f_2^{*} \right|^{1/2} \f$,
#  where:
#  - \f$ f^* \f$-are scaled functions, such \f$ \left| f^*\right| = 1 \f$
#  -  and the norm is   defined as
#   \f$ |left| f(x) \right| \ equiv = \frac{1}{b-a}\int^{b}_{a}  f^2 ( s ) dx \f$
#  @attention It could be rather slow
def _h2_cmp_dist_ ( h1              ,
                    h2              ,
                    density = False ) : 
    """Calculate the norm of difference of scaled histograms/functions 
    |f1-f2|, such |f1|=1 and |f2|=1
    
    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo
    >>> diff = h1.cmp_dist ( h2 )
    - it could be rather slow
    """
    assert isinstance ( h1 , ROOT.TH2 ) and 2 == h1.dim () , \
           "cmp_dist: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 2 == h2.dim () , "cmp_dist: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )
    
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2
        cmp = _h2_cmp_dist_ ( h1_ , h2_ , density = False )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp

    f1   = lambda x , y  : float ( h1 ( x , y ) )
    f2   = lambda x , y  : float ( h2 ( x , y ) )

    xlims = h1.xminmax ()
    ylims = h1.yminmax ()
    
    volume = ( xlims [ 1 ] - xlims [ 0 ] ) * ( ylims [ 1 ] - ylims [ 0 ] )
    params = xlims [ 0 ] , xlims [ 1 ] , ylims [ 0 ] , ylims [ 1 ]

    from ostap.math.integral import integral2 as _integral2_    
    r1 = _integral2_ ( lambda x , y : f1 ( x , y ) ** 2 , *params ) / volume 
    r2 = _integral2_ ( lambda x , y : f2 ( x , y ) ** 2 , *params ) / volume 
        
    sf1  = 1.0 / r1 ** 0.5 
    sf2  = 1.0 / r2 ** 0.5 

    df   = lambda x , y : ( sf1 * f1 ( x , y ) - sf2 * f2 ( x , y ) ) ** 2 
    d12  = _integral2_ ( df , *params ) / volume 
    
    return d12 ** 0.5 

ROOT.TH2F.cmp_dist = _h2_cmp_dist_
ROOT.TH2D.cmp_dist = _h2_cmp_dist_

# =============================================================================
## calculate the "distance" for two histograms or function-like objects
# Distance is defined as  
#  \f$ d = \left| f_1^{*} - f_2^{*} \right|^{1/2} \f$,
#  where:
#  - \f$ f^* \f$-are scaled functions, such \f$ \left| f^*\right| = 1 \f$
#  -  and the norm is   defined as
#   \f$ |left| f(x) \right| \ equiv = \frac{1}{b-a}\int^{b}_{a}  f^2 ( s ) dx \f$ 
#  @attention It could be rather slow
def _h3_cmp_dist_ ( h1              ,
                    h2              ,
                    density = False ) : 
    """Calculate the norm of difference of scaled histograms/functions 
    |f1-f2|, such |f1|=1 and |f2|=1
    
    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo
    >>> diff = h1.cmp_dist ( h2 )
    - it could be rather slow    
    """
    assert isinstance ( h1 , ROOT.TH3 ) and 3 == h1.dim () , \
           "cmp_dist: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 3 == h2.dim () , "cmp_dist: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )
    
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2
        cmp = _h3_cmp_dist_ ( h1_ , h2_ , density = False )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp

    f1   = lambda x , y , z : float ( h1 ( x , y , z ) )
    f2   = lambda x , y , z : float ( h2 ( x , y , z ) )

    xlims = h1.xminmax ()
    ylims = h1.yminmax ()
    zlims = h1.zminmax ()
    
    volume = ( xlims [ 1 ] - xlims [ 0 ] ) * \
             ( ylims [ 1 ] - ylims [ 0 ] ) * \
             ( zlims [ 1 ] - zlims [ 0 ] )
    params = xlims [ 0 ] , xlims [ 1 ] , ylims [ 0 ] , ylims [ 1 ] , zlims [ 0 ] , zlims [ 1 ]

    from ostap.math.integral import integral3 as _integral3_    
    r1 = _integral3_ ( lambda x , y , z : f1 ( x , y , z ) ** 2 , *params ) / volume 
    r2 = _integral3_ ( lambda x , y , z : f2 ( x , y , z ) ** 2 , *params ) / volume 
        
    sf1  = 1.0 / r1 ** 0.5 
    sf2  = 1.0 / r2 ** 0.5 

    df   = lambda x , y , z : ( sf1 * f1 ( x , y , z ) - sf2 * f2 ( x , y , z ) ) ** 2 
    d12  = _integral3_ ( df , *params ) / volume 
    
    return d12 ** 0.5 

ROOT.TH3F.cmp_dist = _h3_cmp_dist_
ROOT.TH3D.cmp_dist = _h3_cmp_dist_



# =============================================================================
## calculate the "discrete distance" for two scaled histograms or function-like objects
# Distance is defined as  
#  \f$ d = \left| f_1^{*} - f_2^{*} \right|^{1/2} \f$,
#  where:
#  - \f$ f^* \f$-are scaled functions, such \f$ \left| f^*\right| = 1 \f$
#  -  and the norm is   defined as
#   \f$ |left| f(x) \right| \ equiv = \frac{1}{b-a}\int^{b}_{a}  f^2 ( s ) dx \f$ 
def _h1_cmp_ddist_ ( h1              ,
                     h2              ,
                     density = False ) : 
    """Calculate the norm of difference of scaled histograms/functions 
    |f1-f2|, such |f1|=1 and |f2|=1
    
    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo
    >>> diff = h1.cmp_dist ( h2 )
    
    """
    assert isinstance ( h1 , ROOT.TH1 ) and 1 == h1.dim () , \
           "cmp_dist: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 1 == h2.dim () , "cmp_dist: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )
    
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2 
        cmp = _h1_cmp_ddist_ ( h1_ , h2_ , density = False )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp

    r1 , r2 = 0.0 , 0.0
    for i , x , v1  in h1.items () :
        xv  = x.value ()
        xe  = x.error () 
        v2  = VE ( h2 ( xv ) )
        r1 += 2 * xe * ( float ( v1 ) ** 2 )
        r2 += 2 * xe * ( float ( v2 ) ** 2 ) 

    lims   = h1.xminmax ()    
    volume = lims [ 1 ] - lims [ 0 ] 

    r1 /= volume
    r2 /= volume

    sf1  = 1.0 / r1 ** 0.5 
    sf2  = 1.0 / r2 ** 0.5
    
    d12 = VE() 
    for i , x , v1  in h1.items () :
        xv   = x.value ()
        xe   = x.error () 
        v2   = VE ( h2 ( xv ) )
        d12 += 2 * xe * ( ( sf1 * v1 - sf2 * v2 ) ** 2 )

    d12 /= volume  

    return d12 ** 0.5

ROOT.TH1D.cmp_ddist = _h1_cmp_ddist_
ROOT.TH1F.cmp_ddist = _h1_cmp_ddist_ 



# =============================================================================
## calculate the "discrete distance" for two scaled histograms or function-like objects
# Distance is defined as  
#  \f$ d = \left| f_1^{*} - f_2^{*} \right|^{1/2} \f$,
#  where:
#  - \f$ f^* \f$-are scaled functions, such \f$ \left| f^*\right| = 1 \f$
#  -  and the norm is   defined as
#   \f$ |left| f(x) \right| \ equiv = \frac{1}{b-a}\int^{b}_{a}  f^2 ( s ) dx \f$ 
def _h2_cmp_ddist_ ( h1              ,
                     h2              ,
                     density = False ) : 
    """Calculate the norm of difference of scaled histograms/functions 
    |f1-f2|, such |f1|=1 and |f2|=1
    
    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo
    >>> diff = h1.cmp_dist ( h2 )
    
    """
    assert isinstance ( h1 , ROOT.TH2 ) and 2 == h1.dim () , \
           "cmp_dist: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 2 == h2.dim () , "cmp_dist: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )
    
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2
        cmp = _h2_cmp_ddist_ ( h1_ , h2_ , density = False )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp

    r1 , r2 = 0.0 , 0.0
    for ix , iy , x , y , v1  in h1.items () :
        xv  = x.value ()
        xe  = x.error () 
        yv  = y.value ()
        ye  = y.error () 
        v2  = VE ( h2 ( xv , yv ) )
        r1 += 4 * xe * ye * ( float ( v1 ) ** 2 )
        r2 += 4 * xe * ye * ( float ( v2 ) ** 2 ) 

    xlims  = h1.xminmax ()    
    ylims  = h1.yminmax ()    
    volume = ( xlims [ 1 ] - xlims [ 0 ] ) * ( ylims [ 1 ] - ylims [ 0 ] ) 

    r1 /= volume
    r2 /= volume

    sf1  = 1.0 / r1 ** 0.5 
    sf2  = 1.0 / r2 ** 0.5
    
    d12 = VE() 
    for ix , iy , x , y , v1  in h1.items () :
        xv   = x.value ()
        xe   = x.error () 
        yv   = y.value ()
        ye   = y.error () 
        v2   = VE ( h2 ( xv , yv ) )
        d12 += 4 * xe * ye * ( ( sf1 * v1 - sf2 * v2 ) ** 2 )

    d12 /= volume  

    return d12 ** 0.5

ROOT.TH2F.cmp_ddist = _h2_cmp_ddist_
ROOT.TH2D.cmp_ddist = _h2_cmp_ddist_


# =============================================================================
## calculate the "discrete distance" for two scaled histograms or function-like objects
# Distance is defined as  
#  \f$ d = \left| f_1^{*} - f_2^{*} \right|^{1/2} \f$,
#  where:
#  - \f$ f^* \f$-are scaled functions, such \f$ \left| f^*\right| = 1 \f$
#  -  and the norm is   defined as
#   \f$ |left| f(x) \right| \ equiv = \frac{1}{b-a}\int^{b}_{a}  f^2 ( s ) dx \f$ 
def _h3_cmp_ddist_ ( h1              ,
                     h2              ,
                     density = False ) : 
    """Calculate the norm of difference of scaled histograms/functions 
    |f1-f2|, such |f1|=1 and |f2|=1
    
    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo
    >>> diff = h1.cmp_dist ( h2 )
    
    """
    assert isinstance ( h1 , ROOT.TH3 ) and 3 == h1.dim () , \
           "cmp_dist: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 3 == h2.dim () , "cmp_dist: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )
    
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2
        cmp = _h3_cmp_ddist_ ( h1_ , h2_ , density = False )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp

    r1 , r2 = 0.0 , 0.0
    for ix , iy , iz , x , y , z , v1  in h1.items () :
        xv  = x.value ()
        xe  = x.error () 
        yv  = y.value ()
        ye  = y.error () 
        zv  = z.value ()
        ze  = z.error () 
        v2  = VE ( h2 ( xv , yv , zv ) )
        r1 += 8 * xe * ye * ze * ( float ( v1 ) ** 2 )
        r2 += 8 * xe * ye * ze * ( float ( v2 ) ** 2 ) 

    xlims  = h1.xminmax ()    
    ylims  = h1.yminmax ()    
    zlims  = h1.zminmax ()    
    volume = ( xlims [ 1 ] - xlims [ 0 ] ) * \
             ( ylims [ 1 ] - ylims [ 0 ] ) * \
             ( zlims [ 1 ] - zlims [ 0 ] ) 

    r1 /= volume
    r2 /= volume

    sf1  = 1.0 / r1 ** 0.5 
    sf2  = 1.0 / r2 ** 0.5
    
    d12 = VE() 
    for ix, iy , iz , x , y , z , v1 in h1.items () :
        xv  = x.value ()
        xe  = x.error () 
        yv  = y.value ()
        ye  = y.error () 
        zv  = z.value ()
        ze  = z.error () 
        v2  = VE ( h2 ( xv , yv , zv ) )
        d12 += 8 * xe * ye * ze * ( ( sf1 * v1 - sf2 * v2 ) ** 2 )

    d12 /= volume  

    return d12 ** 0.5

ROOT.TH3F.cmp_ddist = _h3_cmp_ddist_
ROOT.TH3D.cmp_ddist = _h3_cmp_ddist_



# =============================================================================
## compare two histograms and find the largest difference
#  @code
#  h1 = ...
#  h2 = ...
#  ( x1_min , dy1_min) , ( x1_max , dy1_max ) = h1.cmp_minmax ( h2 )
#  @endcode 
def _h1_cmp_minmax_ ( h1                         ,
                      h2                         ,                      
                      density = False            ,
                      diff    = lambda a,b : b-a , **kwargs ) :
    """Compare two histograms and find the largest difference
    >>> h1 = ...
    >>> h2 = ...
    >>>  ( x1_min , dy1_min ) , ( x1_max , dy1_max ) = h1.cmp_minmax ( h2 )
    """
    
    assert isinstance ( h1 , ROOT.TH1 ) and 1 == h1.dim () , \
           "cmp_minmax: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 1 == h2.dim () , "cmp_minmax: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )

    if density :
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2
        cmp = _h1_cmp_minmax_ ( h1_ , h2_ , density = False , diff = diff , **kwargs )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp

    mn_x   = None
    mx_x   = None 
    mn_val = None
    mx_val = None

    ## loop over  bnis in the first histo 
    for i , x , v in h1.items() :
        
        dv = diff ( v , h2 ( x , **kwargs ) ) ## NOTE  THE ARGUMENTS! 
        if mn_val is None or dv < mn_val :
            mn_val = dv
            mn_x   = x.value()
        if mx_val is None or dv > mx_val :
            mx_val = dv
            mx_x   = x.value() 
            
    if isinstance ( h2 , ROOT.TH1 ) and 1 == h2.dim () : 

        ## loop over  bins in the second histo 
        for i , x , v in h2.items() :            
            dv = diff ( h1 ( x , **kwargs ) , v ) ## NOTE  THE ARGUMENTS! 
            if mn_val is None or dv < mn_val :
                mn_val = dv
                mn_x   = x.value()
            if mx_val is None or dv > mx_val :
                mx_val = dv
                mx_x   = x.value() 


    return ( mn_x , mn_val ) , ( mx_x , mx_val )
                
ROOT.TH1D.cmp_minmax = _h1_cmp_minmax_
ROOT.TH1F.cmp_minmax = _h1_cmp_minmax_ 


# =============================================================================
## compare two historgams and find the largest difference
#  @code
#  h1 = ...
#  h2 = ...
#  (x1_min,y1_min,dz1_min),(x1_max,y1_min,dz1_max) = h1.cmp_minmax ( h2 )
#  @endcode 
def _h2_cmp_minmax_ ( h1                         ,
                      h2                         ,
                      density = False            ,
                      diff    = lambda a,b : b-a , **kwargs ) :
    """Compare two historgams and find the largest difference
    >>> h1 = ...
    >>> h2 = ...
    >>> (x_min,y_min,v_min),(x_max,y_max,v_max) = h1.cmp_minmax ( h2 )
    """

    assert isinstance ( h1 , ROOT.TH2 ) and 2 == h1.dim () , \
           "cmp_minmax: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 2 == h2.dim () , "cmp_minmax: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )

    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2
        r   = _h2_cmp_minmax_  ( h1_  , h2_ , density = False , diff = diff , **kwargs )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return r 

    mn_x   = None
    mx_x   = None
    mn_y   = None 
    mx_y   = None 
    mn_val = None
    mx_val = None
    
    ## loop over the 1st histo bins 
    for ix , iy , x , y , v in h1.items() :
        
        dv = diff ( v , h2 ( x , y , **kwargs ) ) ## NOTE ORDER OF ARGUMENTS 
        if mn_val is None or dv < mn_val :
            mn_val = dv
            mn_x   = x.value() 
            mn_y   = y.value() 
        if mx_val is None or dv > mx_val :
            mx_val = dv
            mx_x   = x.value()
            mx_y   = y.value() 


    if isinstance ( h2 , ROOT.TH2 ) and 2 == h2.dim () : 
        
        ## loop over the 2nd histo bins 
        for ix , iy , x , y , v in h2.items() :
            
            dv = diff ( h1 ( x , y , **kwargs ) , v  ) ## NOTE ORDER OF ARGUMENTS 
            if mn_val is None or dv < mn_val :
                mn_val = dv
                mn_x   = x.value() 
                mn_y   = y.value() 
            if mx_val is None or dv > mx_val :
                mx_val = dv
                mx_x   = x.value()
                mx_y   = y.value() 
                

    return ( mn_x , mn_y , mn_val ) , ( mx_x , mx_y , mx_val )
                
ROOT.TH2F.cmp_minmax = _h2_cmp_minmax_
ROOT.TH2D.cmp_minmax = _h2_cmp_minmax_



# =============================================================================
## compare two historgams and find the largest difference
#  @code
#  h1 = ...
#  h2 = ...
#  (x1_min,y1_min,dz1_min),(x1_max,y1_min,dz1_max) = h1.cmp_minmax ( h2 )
#  @endcode 
def _h3_cmp_minmax_ ( h1                         ,
                      h2                         ,
                      density = False            ,
                      diff    = lambda a,b : b-a , **kwargs ) :
    """Compare two historgams and find the largest difference
    >>> h1 = ...
    >>> h2 = ...
    >>> (x_min,y_min,z_min,v_min),(x_max,y_max,z_max,v_max) = h1.cmp_minmax ( h2 )
    """
    assert isinstance ( h1 , ROOT.TH3 ) and 3 == h1.dim () , \
           "cmp_minmax: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 3 == h2.dim () , "cmp_minmax: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )

    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2
        r   = _h3_cmp_minmax_  ( h1_  , h2_ , density = False , diff = diff , **kwargs )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return r 

    mn_x   = None
    mx_x   = None
    mn_y   = None 
    mx_y   = None 
    mn_z   = None 
    mx_z   = None 
    mn_val = None
    mx_val = None
    
    ## loop over the 1st histo bins 
    for ix , iy , x , y , z , v in h1.items() :
        
        dv = diff ( v , h2 ( x , y , z , **kwargs ) ) ## NOTE ORDER OF ARGUMENTS 
        if mn_val is None or dv < mn_val :
            mn_val = dv
            mn_x   = x.value() 
            mn_y   = y.value() 
            mn_z   = z.value() 
        if mx_val is None or dv > mx_val :
            mx_val = dv
            mx_x   = x.value()
            mx_y   = y.value() 
            mx_z   = z.value() 


    if isinstance ( h2 , ROOT.TH3 ) and 3 == h2.dim () : 
        
        ## loop over the 2nd histo bins 
        for ix , iy , x , y , z , v in h2.items() :
            
            dv = diff ( h1 ( x , y , z , **kwargs ) , v  ) ## NOTE ORDER OF ARGUMENTS 
            if mn_val is None or dv < mn_val :
                mn_val = dv
                mn_x   = x.value() 
                mn_y   = y.value() 
                mn_z   = y.value() 
            if mx_val is None or dv > mx_val :
                mx_val = dv
                mx_x   = x.value()
                mx_y   = y.value() 
                mz_y   = y.value() 
                

    return ( mn_x , mn_y , mn_z , mn_val ) , ( mx_x , mx_y , mx_z , mx_val )
                
ROOT.TH3F.cmp_minmax = _h3_cmp_minmax_
ROOT.TH3D.cmp_minmax = _h3_cmp_minmax_






# =============================================================================
## calculate and print some statistic for comparison
#  @code
#  h1 , h2 = ...
#  h1.cmp_prnt ( h2 )
#  @endcode 
def _h1_cmp_prnt_ ( h1                  ,
                    h2                  ,
                    head1       = ''    ,
                    head2       = ''    ,
                    title       = ''    ,
                    density     = False ,
                    max_moment  = 10    ,
                    exp_moment  = True  , 
                    prefix      = ''    ) : 
    """ Calculate and print some statistic information for two histos
    >>> h1 , h2 = ...
    >>> h1.cmp_prnt ( h2 ) 
    """
    assert isinstance ( h1 , ROOT.TH1 ) and 1 == h1.dim () , \
           "cmp_prnt: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 1 == h2.dim () , "cmp_prnt: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )

    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2
        cmp = _h1_cmp_prnt_ ( h1_ , h2_ , head1 = head1 , head2 = head2 , title = title , density = False , max_moment = max_moment , exp_moment = exp_moment )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp

    if not head1 : head1 = h1.GetName() 
    if not head2 : head2 = h2.GetName()

    fmt    = '%+11.4g +- %-10.4g'
    wid0   = 25

    values = [ 'Mean'     ,
               'Rms'      ,
               'Skewness' , 
               'Kurtosis' ] 

    numbers = []

    mean =  h1.mean     ()  , h2.mean     ()
    rms  =  h1.rms      ()  , h2.rms      ()   
    skew =  h1.skewness ()  , h2.skewness ()   
    kurt =  h1.kurtosis ()  , h2.kurtosis ()
    numbers.append ( mean )
    numbers.append ( rms  )
    numbers.append ( skew )
    numbers.append ( kurt )

    if 4 < max_moment  : 
        for i in range ( 5 , max_moment + 1 ) :
            v1   = h1.stdMoment ( i , exp_moment )
            v2   = h2.stdMoment ( i , exp_moment )
            item = v1 , v2
            numbers.append ( item )
            if   exp_moment : values .append ( 'ExpMom/%d' % i )
            else            : values .append ( 'StdMom/%d' % i )

    numbers = tuple ( numbers )
    values  = tuple ( values  ) 
        
    wid1 = max ( len ( v )  for v in values    )
    wid1 = max ( wid1  , len ( 'Quantity' ) )
    wid2 = max ( wid0  , len ( head1   ) )
    wid3 = max ( wid0  , len ( head2   ) )
    wid4 = max ( wid0  , len ( 'Delta' ) )
    
    header = (  ( '{:^%d}' % wid1 ).format ( 'Quantity' ) ,
                ( '{:^%d}' % wid2 ).format ( head1      ) ,
                ( '{:^%d}' % wid3 ).format ( head2      ) ,
                ( '{:^%d}' % wid4 ).format ( 'Delta'    ) )
    
    table_data = [ header ]

    for v , item in zip ( values , numbers ) :
        v1 , v2 =  item 
        dv  = v1 - v2 
        row = allright ( v ) , v1.toString ( fmt ) , v2.toString( fmt ) , dv.toString ( fmt )         
        table_data.append ( row ) 

    title = title if title else '%s vs %s' % ( head1 , head2 ) 
    import ostap.logger.table as T
    return T.table (  table_data , title , prefix )

    
ROOT.TH1D.cmp_prnt = _h1_cmp_prnt_
ROOT.TH1F.cmp_prnt = _h1_cmp_prnt_ 



# =============================================================================
## calculate and print some statistic for comparison
#  @code
#  h1 , h2 = ...
#  h1.cmp_prnt ( h2 )
#  @endcode 
def _h2_cmp_prnt_ ( h1              ,
                    h2              ,
                    head1   = ''    ,
                    head2   = ''    ,
                    title   = ''    ,
                    density = False ,
                    prefix  = ''    ) : 
    """ Calculate and print some statistic information for two histos
    >>> h1 , h2 = ...
    >>> h1.cmp_prnt ( h2 ) 
    """
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2
        cmp = _h2_cmp_prnt_ ( h1_ , h2_ , head1 = head1 , head2 = head2 , title = title , density = False )
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp

    if not head1 : head1 = h1.GetName() 
    if not head2 : head2 = h2.GetName()

    fmt1   = '%+11.4g +- %-10.4g'
    fmt    = '%+12.5g'
    wid0   = 25
    
    values = ( 'x-mean'      ,
               'y-mean'      ,
               'x-rms'       ,
               'y-rms'       ,
               'xy-corr'     ,
               ##
               'StdMom[0,3]' ,
               'StdMom[1,2]' ,
               'StdMom[2,1]' ,
               'StdMom[3,0]' ,
               ##
               'StdMom[0,4]' ,
               'StdMom[1,3]' ,
               'StdMom[2,2]' ,
               'StdMom[3,1]' ,
               'StdMom[4,0]' ,
               ##
               'StdMom[0,5]' ,
               'StdMom[1,4]' ,
               'StdMom[2,3]' ,
               'StdMom[3,2]' ,
               'StdMom[4,1]' ,
               'StdMom[5,0]' )
    
    functions  =  (
        lambda h : h.xmean      (       ) ,
        lambda h : h.ymean      (       ) ,
        lambda h : h.xrms       (       ) ,
        lambda h : h.yrms       (       ) ,
        lambda h : h.xycorr     (       ) ,
        lambda h : h.std_moment ( 0 , 3 ) , 
        lambda h : h.std_moment ( 1 , 2 ) , 
        lambda h : h.std_moment ( 2 , 1 ) , 
        lambda h : h.std_moment ( 3 , 0 ) ,
        lambda h : h.std_moment ( 0 , 4 ) , 
        lambda h : h.std_moment ( 1 , 3 ) , 
        lambda h : h.std_moment ( 2 , 2 ) , 
        lambda h : h.std_moment ( 3 , 1 ) ,
        lambda h : h.std_moment ( 4 , 0 ) ,       
        lambda h : h.std_moment ( 0 , 5 ) , 
        lambda h : h.std_moment ( 1 , 4 ) , 
        lambda h : h.std_moment ( 2 , 3 ) , 
        lambda h : h.std_moment ( 3 , 2 ) ,
        lambda h : h.std_moment ( 4 , 1 ) ,       
        lambda h : h.std_moment ( 5 , 0 ) ,       
        )
    
    wid1 = max ( len(v)  for v in values    )
    wid1 = max ( wid1  , len ( 'Quantity' ) )
    wid2 = max ( wid0  , len ( head1   ) )
    wid3 = max ( wid0  , len ( head2   ) )
    wid4 = max ( wid0  , len ( 'Delta' ) )
    
    header = (  ( '{:^%d}' % wid1 ).format ( 'Quantity' ) ,
                ( '{:^%d}' % wid2 ).format ( head1      ) ,
                ( '{:^%d}' % wid3 ).format ( head2      ) ,
                ( '{:^%d}' % wid4 ).format ( 'Delta'    ) )
    
    table_data = [ header ]


    fmt2 = '{:^%d}' % wid2
    fmt3 = '{:^%d}' % wid3
    fmt4 = '{:^%d}' % wid4
    
    import itertools 
    for i , v , f in zip ( itertools.count() , values , functions ) :
        v1 = f ( h1 )
        v2 = f ( h2 )
        dv = v1 - v2
        if i <= 3 :
            row = allright ( v ) , v1.toString( fmt1 ) , v2.toString( fmt1 ) , dv.toString  ( fmt1 )
        else :
            row = allright ( v ) , fmt2. format ( fmt % v1 ) , fmt3. format ( fmt % v2 ) , fmt4. format ( fmt % dv  ) 
            
        table_data.append ( row ) 

    title = title if title else '%s vs %s' % ( head1 , head2 ) 
    import ostap.logger.table as T
    return T.table (  table_data , title , prefix )


ROOT.TH2F.cmp_prnt = _h2_cmp_prnt_
ROOT.TH2D.cmp_prnt = _h2_cmp_prnt_


# =============================================================================
## calculate and print some statistic for comparison   of 1D-histograms 
#  @code
#  h1 , h2 = ...
#  h1.cmp_diff_prnt ( h2 )
#  @endcode 
def _h1_cmp_diff_prnt_ ( h1                             ,
                         h2                             ,
                         title           = ''           ,
                         density         = False        ,
                         ##
                         distance        = 'distance'   ,
                         ddistance       = 'distance/D' ,
                         diffneg         = 'diff(neg)'  ,
                         diffpos         = 'diff(pos)'  ,
                         angle           = 'angle'      ,
                         dangle          = 'angle/D'    ,
                         chi2ndf         = 'chi2/ndf'   ,
                         probchi2        = 'prob(chi2)' ,
                         probfit         = 'prob(fit)'  ,
                         prefix          = ''           ) : 
    """ Calculate and print some statistic information for two 1D-histos
    >>> h1 , h2 = ...
    >>> h1.cmp_diff_prnt ( h2 ) 
    """
    assert isinstance ( h1 , ROOT.TH1 ) and 1 == h1.dim () , \
           "cmp_diff_prnt: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 1 == h2.dim () , "cmp_diff_prnt: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )

    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2
        cmp = _h1_cmp_diff_prnt_ ( h1_                   ,
                                   h2_                   ,
                                   title     = title     ,
                                   density   = False     ,
                                   distance  = distance  ,
                                   ddistance = ddistance ,
                                   diffneg   = diffneg   ,
                                   diffpos   = diffpos   ,                                   
                                   angle     = angle     ,
                                   dangle    = dangle    ,
                                   chi2ndf   = chi2ndf   ,
                                   probchi2  = probchi2  ,
                                   probfit   = probfit   ,
                                   prefix    = prefix    ) 
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp

    rows = [ ('Quantity' , '' , 'value' ) ]

    histo2 = isinstance ( h2 , ROOT.TH1 ) and 1 == h2.dim() 

    from   ostap.logger.utils import pretty_float, pretty_ve
    import ostap.math.math_ve as     ME
    
    if distance :
        value = h1.cmp_dist ( h2 , density = density )
        v , n = pretty_float ( value  )
        row   = distance , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )
        
    if ddistance :
        value = h1.cmp_ddist ( h2 , density = density )
        v , n = pretty_ve ( value  , parentheses = False )
        row   = ddistance , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )

    if ddistance and histo2 :
        value = h2.cmp_ddist ( h1 , density = density )
        v , n = pretty_ve ( value  , parentheses = False )
        row   = ddistance , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )

    def fmt (  v  ) :
        vv  = abs ( v.value() )
        if   100 <= vv : return '%+8.4f +/- %-08f'
        elif  10 <= vv : return '%+8.5f +/- %-08f'
        return                  '%+8.6f +/- %-08f'
        
    if diffneg :        
        dmn , _ = h1.cmp_minmax ( h2                     ,
                                  density = density      ,
                                  diff    = lambda a , b : 2 * a.asym ( b ) )
        value = dmn [ -1 ] * 100        
        row   = diffneg , '[%]' , value.toString( fmt ( value ) ) 
        rows.append ( row  )
        
    if diffpos :        
        _ , dmx  = h1.cmp_minmax ( h2                     ,
                                   density = density      ,
                                   diff    = lambda a , b : 2 * a.asym ( b ) )         
        value = dmx [ -1 ]  * 100       
        row   = diffpos , '[%]' , value.toString( fmt ( value ) ) 
        rows.append ( row  )
        
    if angle :
        value = h1.cmp_cos ( h2 , density = density )
        value = math.acos    ( value ) 
        v , n = pretty_float ( value )
        row   = angle , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )

    if dangle :
        value = h1.cmp_dcos ( h2 , density = density )
        value = ME.acos     ( value ) 
        v , n = pretty_ve ( value  , parentheses = False )
        row   = dangle , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )

    if dangle and histo2 :
        value = h2.cmp_dcos ( h1 , density = density )
        value = ME.acos     ( value ) 
        v , n = pretty_ve ( value  , parentheses = False )
        row   = dangle , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )

    if chi2ndf :
        chi2, prob = h1.cmp_chi2 ( h2 , density = density )

        value = chi2
        v , n = pretty_float ( value )
        row   = chi2ndf , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )

    if chi2ndf and histo2 :
        
        chi2, prob = h2.cmp_chi2 ( h1 , density = density )
        
        value = chi2
        v , n = pretty_float ( value )
        row   = chi2ndf , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )        
            
    if probchi2 : 
        chi2, prob = h1.cmp_chi2 ( h2 , density = density )

        value = prob
        v , n = pretty_float ( value )
        row   = probchi2 , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )

    if probchi2 and histo2 :
        
        chi2, prob = h2.cmp_chi2 ( h1 , density = density )

        value = prob
        v , n = pretty_float ( value )
        row   = probchi2 , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )

    if probfit and histo2 :

        rf1   = h1.cmp_fit ( h2 , density = density ) 
        prob  = rf1.Prob()
        
        value = prob 
        v , n = pretty_float ( value )
        row   = probfit , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )

        rf2   = h2.cmp_fit ( h1 , density = density ) 
        prob  = rf2.Prob()
        
        value = prob 
        v , n = pretty_float ( value )
        row   = probfit , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )
        
        
    import ostap.logger.table as T
    return T.table ( rows , title =  title , prefix = prefix , alignment = 'lcl' )

ROOT.TH1D.cmp_diff_prnt = _h1_cmp_diff_prnt_
ROOT.TH1F.cmp_diff_prnt = _h1_cmp_diff_prnt_ 
    

# =============================================================================
## calculate and print some statistic for comparison of 2D-histograms 
#  @code
#  h1 , h2 = ...
#  h1.cmp_diff_prnt ( h2 )
#  @endcode 
def _h2_cmp_diff_prnt_ ( h1                             ,
                         h2                             ,
                         title           = ''           ,
                         density         = False        ,
                         ##
                         distance        = 'distance'   ,
                         ddistance       = 'distance/D' ,
                         diffneg         = 'diff(neg)'  ,
                         diffpos         = 'diff(pos)'  ,
                         angle           = 'angle'      ,
                         dangle          = 'angle/D'    ,
                         chi2ndf         = 'chi2/ndf'   ,
                         probchi2        = 'prob(chi2)' ,
                         prefix          = ''           ) : 
    """ Calculate and print some statistic information for two histos
    >>> h1 , h2 = ...
    >>> h1.cmp_prnt ( h2 ) 
    """
    assert isinstance ( h1 , ROOT.TH2 ) and 2 == h1.dim () , \
           "cmp_diff_prnt: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    
    if isinstance ( h2 , ROOT.TH1 ) :
        assert 2 == h2.dim () , "cmp_diff_prnt: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )

    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1 
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2 
        cmp = _h2_cmp_diff_prnt_ ( h1_                   ,
                                   h2_                   ,
                                   title     = title     ,
                                   density   = False     ,
                                   distance  = distance  ,
                                   ddistance = ddistance ,
                                   angle     = angle     ,
                                   dangle    = dangle    ,
                                   chi2ndf   = chi2ndf   ,
                                   probchi2  = probchi2  ,
                                   prefix    = prefix    ) 
        if h1_ is not h1 : del h1_
        if h2_ is not h2 : del h2_
        return cmp

    rows = [ ('Quantity' , '' , 'value' ) ]

    histo2 = isinstance ( h2 , ROOT.TH2 ) and 2 == h2.dim() 

    from   ostap.logger.utils import pretty_float, pretty_ve
    import ostap.math.math_ve as     ME

    ## if distance :
    ##     value = h1.cmp_dist ( h2 , density = density )
    ##     v , n = pretty_float ( value  )
    ##     row   = distance , '[10^%d]' %n if  n  else  '' , v 
    ##     rows.append ( row  )
        
    if ddistance :
        value = h1.cmp_ddist ( h2 , density = density )
        v , n = pretty_ve ( value  , parentheses = False )
        row   = ddistance , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )

    if ddistance and histo2 :
        value = h2.cmp_ddist ( h1 , density = density )
        v , n = pretty_ve ( value  , parentheses = False )
        row   = ddistance , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )

    def fmt (  v  ) :
        vv  = abs ( v.value() )
        if   100 <= vv : return '%+8.4f +/- %-08f'
        elif  10 <= vv : return '%+8.5f +/- %-08f'
        return                  '%+8.6f +/- %-08f'
        
    if diffneg :        
        dmn , _ = h1.cmp_minmax ( h2                     ,
                                  density = density      ,
                                  diff    = lambda a , b : 2 * a.asym ( b ) )
        value = dmn [ -1 ] * 100        
        row   = diffneg , '[%]' , value.toString( fmt ( value ) ) 
        rows.append ( row  )
        
    if diffpos :        
        _ , dmx  = h1.cmp_minmax ( h2                     ,
                                   density = density      ,
                                   diff    = lambda a , b : 2 * a.asym ( b ) )         
        value = dmx [ -1 ]  * 100       
        row   = diffpos , '[%]' , value.toString( fmt ( value ) ) 
        rows.append ( row  )
        
    ## if angle :
    ##     value = h1.cmp_cos ( h2 , density = density )
    ##     value = math.acos    ( value ) 
    ##     v , n = pretty_float ( value )
    ##     row   = angle , '[10^%d]' %n if  n  else  '' , v 
    ##     rows.append ( row  )

    if dangle :
        value = h1.cmp_dcos ( h2 , density = density )
        value = ME.acos     ( value ) 
        v , n = pretty_ve ( value  , parentheses = False )
        row   = dangle , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )

    if dangle and histo2 :
        value = h2.cmp_dcos ( h1 , density = density )
        value = ME.acos     ( value ) 
        v , n = pretty_ve ( value  , parentheses = False )
        row   = dangle , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )

    if chi2ndf :
        chi2, prob = h1.cmp_chi2 ( h2 , density = density )

        value = chi2
        v , n = pretty_float ( value )
        row   = chi2ndf , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )

    if chi2ndf and histo2 :
        
        chi2, prob = h2.cmp_chi2 ( h1 , density = density )
        
        value = chi2
        v , n = pretty_float ( value )
        row   = chi2ndf , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )        
            
    if probchi2 : 
        chi2, prob = h1.cmp_chi2 ( h2 , density = density )

        value = prob
        v , n = pretty_float ( value )
        row   = probchi2 , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )

    if probchi2 and histo2 :
        
        chi2, prob = h2.cmp_chi2 ( h1 , density = density )

        value = prob
        v , n = pretty_float ( value )
        row   = probchi2 , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )

    import ostap.logger.table as T
    return T.table ( rows , title =  title , prefix = prefix , alignment = 'lcl' )

ROOT.TH2F.cmp_diff_prnt = _h2_cmp_diff_prnt_
ROOT.TH2D.cmp_diff_prnt = _h2_cmp_diff_prnt_
    

        

# =============================================================================
_decorated_classes_ = (
    ROOT.TH1  , 
    ROOT.TH2  , 
    ROOT.TH3  , 
    ROOT.TH1F , 
    ROOT.TH1D ,
    ROOT.TH2F , 
    ROOT.TH2D ,
    ROOT.TH3F , 
    ROOT.TH3D ,
    )
_new_methods_       = (
    #
    ROOT.TH1D.is_constant   ,
    ROOT.TH1F.is_constant   ,
    #
    ROOT.TH1D.cmp_fit       , 
    ROOT.TH1F.cmp_fit       ,
    #
    ROOT.TH1D.cmp_pdf       , 
    ROOT.TH1F.cmp_pdf       ,
    #
    ROOT.TH1D.cmp_chi2      , 
    ROOT.TH1F.cmp_chi2      ,
    #
    ROOT.TH2F.cmp_chi2       , 
    ROOT.TH2D.cmp_chi2       , 
    ROOT.TH3F.cmp_chi2       ,
    ROOT.TH3D.cmp_chi2       ,
    #
    ROOT.TH1D.chi2_cmp      , 
    ROOT.TH1F.chi2_cmp      ,
    #
    ROOT.TH1D.cmp_cos       , 
    ROOT.TH1F.cmp_cos       ,
    #
    ROOT.TH2F.cmp_cos        , 
    ROOT.TH2D.cmp_cos        , 
    ROOT.TH3F.cmp_cos        ,
    ROOT.TH3D.cmp_cos        ,
    #
    ROOT.TH1D.cmp_dcos      , 
    ROOT.TH1F.cmp_dcos      ,
    #
    ROOT.TH2F.cmp_dcos       , 
    ROOT.TH2D.cmp_dcos       , 
    ROOT.TH3F.cmp_dcos       ,
    ROOT.TH3D.cmp_dcos       ,
    #
    ROOT.TH1D.cmp_dist      , 
    ROOT.TH1F.cmp_dist      ,
    ROOT.TH2D.cmp_dist      , 
    ROOT.TH2F.cmp_dist      , 
    ROOT.TH3D.cmp_dist      , 
    ROOT.TH3F.cmp_dist      ,
    #
    ROOT.TH1D.cmp_ddist     , 
    ROOT.TH1F.cmp_ddist     ,
    #
    ROOT.TH1D.cmp_prnt      ,  
    ROOT.TH1F.cmp_prnt      ,
    ##
    ROOT.TH2F.cmp_prnt      , 
    ROOT.TH2D.cmp_prnt      , 
    ##
    ROOT.TH1D.cmp_minmax    ,
    ROOT.TH1F.cmp_minmax    ,
    ## 
    ROOT.TH2F.cmp_minmax    ,
    ROOT.TH2D.cmp_minmax    ,
    ##
    ROOT.TH3F.cmp_minmax    ,
    ROOT.TH3D.cmp_minmax    ,
    ##
    ROOT.TH1D.cmp_diff_prnt , 
    ROOT.TH1F.cmp_diff_prnt ,
    ## 
    ROOT.TH2F.cmp_diff_prnt , 
    ROOT.TH2D.cmp_diff_prnt , 
    )
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
