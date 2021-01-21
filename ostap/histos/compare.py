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
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1.clone() 
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2.clone() 
        cmp = _h1_cmp_fit_ ( h1_ , h2_ ,  density = False , opts = opts )
        del h1_
        del h2_
        return cmp
    
    f2 = h2.asTF () 
    f2.ReleaseParameter ( 0 ) 

    rf = h1.Fit ( f2 , 'S' + opts ) 
    if 0 != rf.Status() :
        logger.warning("Can't fit with function " % rf.Status() )
        return None

    return rf

ROOT.TH1D.cmp_fit = _h1_cmp_fit_
ROOT.TH1F.cmp_fit = _h1_cmp_fit_ 


# =============================================================================
## compare the 1D-histograms trying to fit one 
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
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1.clone() 
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2.clone() 
        cmp = _h1_cmp_pdf_ ( h1_ , h2_ , density = False ,draw =  draw , silent = silent )
        del h1_
        del h2_
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
def _h1_cmp_chi2_ ( h1              ,
                    h2              ,
                    density = False ) :
    """Compare histograms by chi2
    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo (or function or anything else) 
    >>> chi2ndf,prob  = h1.cmp_chi2 ( h2 )    
    """
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1.clone() 
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2.clone() 
        cmp = _h1_cmp_chi2_ ( h1_ , h2_ , density = False )
        del h1_
        del h2_
        return cmp

    c2  = 0
    ndf = 0  
    for entry in h1.items() :
        
        x     = entry[1]
        y1    = entry[2]
        
        y2    = h2 ( x.value() )

        c2   += y1.chi2 ( y2 )
        ndf  += 1 

    c2ndf = c2/ndf 
    return c2ndf, ROOT.TMath.Prob ( c2 , ndf ) 

ROOT.TH1D.cmp_chi2 = _h1_cmp_chi2_
ROOT.TH1F.cmp_chi2 = _h1_cmp_chi2_ 

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
        _func_  = lambda x,xl,xr : func.integral(xl,xr)/(xr-xl) 
    elif integral and hasattr ( func , 'Integral' ) :  
        _func_  = lambda x,xl,xr : func.Integral(xl,xr)/(xr-xl) 
    elif integral :
        ## use numerical integration 
        from ostap.math.intergal import integral as _integral_
        _func_  = lambda x , xl , xr : _integral_ ( func , xl , xr ) / ( xr - xl )
        
    for entry in h1.items() :
        
        x    = entry[1]
        y1   = entry[2]
        
        xv   = x.value()
        xe   = x.error()
        xl   = xv - xe
        xr   = xv + xe
        
        y2   = _func_ ( x , xl , xr )        
        if not select ( x, y1 , y2 ) : continue

        c2  += chi2 ( y1 , y2 )
        ndf += 1
        
    c2ndf = c2/ndf 
    return c2ndf, ROOT.TMath.Prob( c2 , ndf ) 

ROOT.TH1D.chi2_cmp = _h1_chi2_cmp_
ROOT.TH1F.chi2_cmp = _h1_chi2_cmp_ 

# =============================================================================
## compare the 1D-historgams (as functions)
#  calculate
# \f$cos \theta = \frac{ f_1 \cdot f_2 } { \left|f_1\right|\left|f_2\right| }\f$
def _h1_cmp_costheta_ ( h1              ,
                        h2              ,
                        density = False ) :  
    """Compare the 1D-historgams (as functions)
    Calculate scalar product and get ``the angle'' from it
    
    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo
    >>> cos_theta  = h1.cmp_cos ( h2 )
    
    """
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1.clone() 
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2.clone() 
        cmp = _h1_cmp_costheta_ ( h1_ , h2_ , density  = False )
        del h1_
        del h2_
        return cmp
        
    f1 = h1.asFunc   ()
    f2 = h2.asFunc   ()

    lims = h1.xminmax()
    
    from ostap.math.integral import integral as _integral_
    vr1   = _integral_ ( lambda x : f1( x )**2    , lims[0] , lims[1] , limit = 200 , err = False )
    vr2   = _integral_ ( lambda x : f2( x )**2    , lims[0] , lims[1] , limit = 200 , err = False )
    vr12  = _integral_ ( lambda x : f1( x )*f2(x) , lims[0] , lims[1] , limit = 200 , err = False )

    return vr12 / ( vr1 * vr2 ) ** 0.5 

ROOT.TH1D.cmp_cos = _h1_cmp_costheta_
ROOT.TH1F.cmp_cos = _h1_cmp_costheta_ 

# =============================================================================
## compare the 1D-historgams (as functions)
#  calculate "discrete"
# \f$cos \theta = \frac{ f_1 \cdot f_2 } { \left|f_1\right|\left|f_2\right| }\f$
def _h1_cmp_dcostheta_ ( h1              ,
                        h2              ,
                        density = False ) :  
    """Compare the 1D-historgams (as functions)
    Calculate scalar product and get ``the angle'' from it
    
    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo
    >>> cos_theta  = h1.cmp_dcos ( h2 )
    
    """
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1.clone() 
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2.clone() 
        cmp = _h1_cmp_dcostheta_ ( h1_ , h2_ , density  = False )
        del h1_
        del h2_
        return cmp

    r1 , r2 , r12 = 0.0 , 0.0 , VE() 
    for i ,  x , y1  in h1.items () :
        xv   = x.value ()
        xe   = x.error () 
        y2   = VE ( h2 ( xv ) )
        r1  += 2 * xe * ( y1.value () ** 2 )
        r2  += 2 * xe * ( y2.value () ** 2 ) 
        r12 += 2 * xe * ( y1 * y2 ) 

    return r12 / ( r1 * r2 ) ** 0.5 

ROOT.TH1D.cmp_dcos = _h1_cmp_dcostheta_
ROOT.TH1F.cmp_dcos = _h1_cmp_dcostheta_ 


# =============================================================================
## calculate the "distance" for two scaled histograms or function-like objects
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
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1.clone() 
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2.clone() 
        cmp = _h1_cmp_dist_ ( h1_ , h2_ , density = False )
        del h1_
        del h2_
        return cmp
    
    f1   = h1.asFunc  ()
    f2   = h2.asFunc  () if isinstance ( h2 , ROOT.TH1 ) else h2 

    lims = h1.xminmax ()
    
    volume = lims [ 1 ] - lims [ 0 ]

    from ostap.math.integral import integral as _integral_
    
    r1 = _integral_ ( lambda x : f1 ( x ) ** 2 , lims [ 0 ] , lims [ 1 ] , limit = 300 , err = False ) / volume 
    r2 = _integral_ ( lambda x : f2 ( x ) ** 2 , lims [ 0 ] , lims [ 1 ] , limit = 300 , err = False ) / volume 
    
    import math 
    
    sf1  = 1.0 / r1 ** 0.5 
    sf2  = 1.0 / r2 ** 0.5 

    dfn  = lambda x : ( sf1 * f1 ( x ) - sf2 * f2 ( x ) ) ** 2 
    d12  = _integral_ ( dfn , lims [ 0 ] , lims [ 1 ] , limit = 300 , err = False ) / volume 
    
    return d12 ** 0.5 

ROOT.TH1D.cmp_dist = _h1_cmp_dist_
ROOT.TH1F.cmp_dist = _h1_cmp_dist_ 

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
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1.clone() 
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2.clone() 
        cmp = _h1_cmp_ddist_ ( h1_ , h2_ , density = False )
        del h1_
        del h2_
        return cmp

    r1 , r2 = 0.0 , 0.0
    for i , x , y1  in h1.items () :
        xv  = x.value ()
        xe  = x.error () 
        y2  = VE ( h2 ( xv ) )
        r1 += 2 * xe * ( y1.value () ** 2 )
        r2 += 2 * xe * ( y2.value () ** 2 ) 

    lims = h1.xminmax ()    
    volume = lims [ 1 ] - lims [ 0 ] 

    r1 /= volume
    r2 /= volume

    sf1  = 1.0 / r1 ** 0.5 
    sf2  = 1.0 / r2 ** 0.5
    
    d12 = VE() 
    for i , x , y1  in h1.items () :
        xv  = x.value ()
        xe  = x.error () 
        y2  = VE ( h2 ( xv ) )

        d12 += 2 * xe * ( ( sf1 * y1 - sf2 * y2 ) ** 2 )

    d12 /= volume  

    return d12 ** 0.5

ROOT.TH1D.cmp_ddist = _h1_cmp_ddist_
ROOT.TH1F.cmp_ddist = _h1_cmp_ddist_ 


# =============================================================================
## calculate the norm of difference of scaled histograms/functions 
#  \f$ d = \left| f_1^{*} - f_2^{*}\right| \f$,
#  where \f$ f^* \f$-are scaled functions, such \f$ \left| f^*\right| = 1 \f$ 
def _h2_cmp_dist_ ( h1              ,
                    h2              ,
                    density = False ) : 
    """Calculate the norm of difference of scaled histograms/functions 
    |f1-f2|, such |f1|=1 and |f2|=1
    
    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo
    >>> diff = h1.cmp_dist ( h2 )
    
    """
    assert isinstance ( h1 , ROOT.TH2 ) and not isinstance ( h1 , ROOT.TH3 ) , \
           "cmp_dist: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    assert isinstance ( h2 , ROOT.TH2 ) and not isinstance ( h2 , ROOT.TH3 ) , \
           "cmp_dist: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )
    
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1.clone() 
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2.clone() 
        cmp = _h2_cmp_dist_ ( h1_ , h2_ , density = False )
        del h1_
        del h2_
        return cmp

    

    xlims = h1.xminmax()
    ylims = h1.yminmax()
    
    f1 = lambda x,y : h1(x,y).value() 
    f2 = lambda x,y : h2(x,y).value()
    
    from ostap.math.integral import integral2 as _integral2_
    r1   = _integral2_ ( lambda x,y : f1 ( x , y )**2    ,
                         xlims[0] , xlims[1] ,
                         ylims[0] , ylims[1] , err = True )
    r2   = _integral2_ ( lambda x,y : f2 ( x , y )**2    ,
                         xlims[0] , xlims[1] ,
                         ylims[0] , ylims[1] , err = True )
    
    import math 
    
    sf1  = 1.0 / math.sqrt ( r1.value() ) 
    sf2  = 1.0 / math.sqrt ( r2.value() ) 
    
    d12  = _integral2_ ( lambda x,y : (sf1*f1(x,y)-sf2*f2(x,y))**2 ,
                         xlims[0] , xlims[1] ,
                         ylims[0] , ylims[1] , err = True )
    
    volume = ( xlims[1] - xlims[0] ) * ( ylims[1] - ylims[0] ) 
        
    return d12 / volume 


ROOT.TH2D.cmp_dist = _h2_cmp_dist_
ROOT.TH2F.cmp_dist = _h2_cmp_dist_ 


# =============================================================================
## calculate the norm of difference of scaled histograms/functions 
#  \f$ d = \left| f_1^{*} - f_2^{*}\right| \f$,
#  where \f$ f^* \f$-are scaled functions, such \f$ \left| f^*\right| = 1 \f$ 
def _h3_cmp_dist_ ( h1              ,
                    h2              ,
                    density = False ) : 
    """Calculate the norm of difference of scaled histograms/functions 
    |f1-f2|, such |f1|=1 and |f2|=1
    
    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo
    >>> diff = h1.cmp_dist ( h2 )
    
    """
    assert isinstance ( h3 , ROOT.TH3 ) , \
           "cmp_dist: invalid type of h1  %s/%s" % ( h1 , type ( h1 ) )
    assert isinstance ( h2 , ROOT.TH3 ) , \
           "cmp_dist: invalid type of h2  %s/%s" % ( h2 , type ( h2 ) )
    
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1.clone() 
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2.clone() 
        cmp = _h3_cmp_dist_ ( h1_ , h2_ , density = False )
        del h1_
        del h2_
        return cmp

    
    xlims = h1.xminmax()
    ylims = h1.yminmax()
    zlims = h1.zminmax()
    
    f1 = lambda x,y,z : h1(x,y,z).value() 
    f2 = lambda x,y,z : h2(x,y,z).value()
    
    from ostap.math.integral import integral3 as _integral3_
    r1   = _integral3_ ( lambda x,y,z : f1 ( x , y , z )**2    ,
                         xlims[0] , xlims[1] ,
                         ylims[0] , ylims[1] ,
                         zlims[0] , zlims[1] , err = True )
    r2   = _integral3_ ( lambda x,y,z : f2 ( x , y , z )**2    ,
                         xlims[0] , xlims[1] ,
                         ylims[0] , ylims[1] ,
                         zlims[0] , zlims[1] , err = True )
    import math 
    
    sf1  = 1.0 / math.sqrt ( r1.value() ) 
    sf2  = 1.0 / math.sqrt ( r2.value() ) 
    
    d12  = _integral3_ ( lambda x,y,z : (sf1*f1(x,y,z)-sf2*f2(x,y,z))**2 ,
                         xlims[0] , xlims[1] ,
                         ylims[0] , ylims[1] ,
                         zlims[0] , zlims[1] , err = True )
    
    volume = ( xlims[1] - xlims[0] ) * ( ylims[1] - ylims[0] ) * ( zlims[1] - zlims[0] ) 
    
    return d12 / volume 


ROOT.TH3D.cmp_dist = _h3_cmp_dist_
ROOT.TH3F.cmp_dist = _h3_cmp_dist_ 


# =============================================================================



# =============================================================================
## calculate the norm of difference of scaled histograms/functions 
#  \f$ d = \left| (f_1^{*}-f_2^{*})^2/(f_1^{*}f_2^*(x))\right| \f$,
#  where \f$ f^* \f$-are scaled functions, such \f$ \left| f^*\right| = 1 \f$ 
def _h1_cmp_dist2_ ( h1              ,
                     h2              ,
                     density = False ) :   
    """Calculate the norm of difference of scaled histograms/functions 
    |(f1-f2)**2/(f1*f2)|, such |f1|=1 and |f2|=1

    >>> h1 = ... ## the first histo
    >>> h2 = ... ## the second histo
    >>> diff = h1.cmp_dist2 ( h2 )
    
    """
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1.clone() 
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2.clone() 
        cmp = _h1_cmp_dist2_ ( h1_ , h2_ , density = False )
        del h1_
        del h2_
        return cmp

    f1 = h1.asFunc   ()
    f2 = h2.asFunc   ()

    lims = h1.xminmax()
         
    from ostap.math.integral import integral as _integral_
    r1   = _integral_( lambda x : f1 ( x )**2    , lims[0] , lims[1] , limit = 200 , err = True )
    r2   = _integral_( lambda x : f2 ( x )**2    , lims[0] , lims[1] , limit = 200 , err = True )
    
    import math 
    
    sf1 = 1.0 / math.sqrt ( r1.value() ) 
    sf2 = 1.0 / math.sqrt ( r2.value() ) 
    
    def  _func_   ( x ) :
        v1 =  sf1 * f1 ( x )
        v2 =  sf2 * f2 ( x )
        v  = (v1-v2)*(v1-v2)/abs(v1*v2)
        return v*v 
    
    d12  = _integral_ ( _func_ , lims[0] , lims[1] , limit = 200 , err = True )

    return d12 

ROOT.TH1D.cmp_dist2 = _h1_cmp_dist2_
ROOT.TH1F.cmp_dist2 = _h1_cmp_dist2_ 

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
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1.clone() 
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2.clone() 
        cmp = _h1_cmp_prnt_ ( h1_ , h2_ , head1 = head1 , head2 = head2 , title = title , density = False , max_moment = max_moment , exp_moment = exp_moment )
        del h1_
        del h2_
        return cmp

    if not head1 : head1 = h1.GetName() 
    if not head2 : head2 = h1.GetName()

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
#  h1.cmp_diff_prnt ( h2 )
#  @endcode 
def _h1_cmp_diff_prnt_ ( h1                             ,
                         h2                             ,
                         title           = ''           ,
                         density         = False        ,
                         ##
                         distance        = 'distance'   ,
                         ddistance       = 'distance/D' ,
                         angle           = 'angle'      ,
                         dangle          = 'angle/D'    ,
                         chi2ndf         = 'chi2/ndf'   ,
                         probchi2        = 'prob(chi2)' ,
                         probfit         = 'prob(fit)'  ,
                         prefix          = ''           ) : 
    """ Calculate and print some statistic information for two histos
    >>> h1 , h2 = ...
    >>> h1.cmp_prnt ( h2 ) 
    """
    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1.clone() 
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2.clone() 
        cmp = _h1_cmp_diff_prnt_ ( h1_                   ,
                                   h2_                   ,
                                   title     = title     ,
                                   density   = False     ,
                                   distance  = distance  ,
                                   ddistance = ddistance ,
                                   angle     = angle     ,
                                   dangle    = dangle    ,
                                   chi2ndf   = chi2ndf   ,
                                   probchi2  = probchi2  ,
                                   probfit   = probfit   ,
                                   prefix    = prefix    ) 
        del h1_
        del h2_
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

    if angle :
        value = h1.cmp_cos ( h2 , density = density )
        value = math.acos    ( value ) 
        v , n = pretty_float ( value )
        row   = angle , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )

    if dangle :
        value = h1.cmp_dcos ( h2 , density = density )
        import ostap.math.math_ve as ME 
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

        rf1   = h1.cmp_fit ( h2 , density = density , opts = '' ) 
        prob  = rf1.Prob()
        
        value = prob 
        v , n = pretty_float ( value )
        row   = probfit , '[10^%d]' %n if  n  else  '' , v 
        rows.append ( row  )

        rf2   = h2.cmp_fit ( h1 , density = density , opts = '' ) 
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
    
    if density :
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1.clone() 
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2.clone() 
        cmp = _h1_cmp_minmax_ ( h1_ , h2_ , density = False , diff = diff , **kwargs )
        del h1_
        del h2_
        return cmp

    mn_x   = None
    mx_x   = None 
    mn_val = None
    mx_val = None

    ## loop over  bnis in the first histo 
    for i , x , y in h1.items() :
        
        dy = diff ( y , h2 ( x , **kwargs ) ) ## NOTE  THE ARGUMENTS! 
        if mn_val is None or dy.value() < mn_val.value() :
            mn_val = dy
            mn_x   = x.value()
        if mx_val is None or dy.value() > mx_val.value() :
            mx_val = dy
            mx_x   = x.value() 

    ## loop over  bnis in the second histo 
    for i , x , y in h2.items() :
        
        dy = diff ( h1 ( x , **kwargs ) , y ) ## NOTE  THE ARGUMENTS! 
        if mn_val is None or dy.value() < mn_val.value() :
            mn_val = dy
            mn_x   = x.value()
        if mx_val is None or dy.value() > mx_val.value() :
            mx_val = dy
            mx_x   = x.value() 

    return ( mn_x , mn_val) , ( mx_x , mx_val )
                
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
    (x1_min,y1_min,dz1_min),(x1_max,y1_min,dz1_max) = h1.cmp_minmax ( h2 )
    (x2_min,y2_min,dz2_min),(x2_max,y2_min,dz2_max) = h2.cmp_minmax ( h1 )
    """

    if density : 
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1.clone() 
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2.clone() 
        r   = _h2_cmp_minmax_  ( h1_  , h2_ , density = False , diff = diff , **kwargs )
        del h1_
        del h2_
        return r 

    mn_x   = None
    mx_x   = None
    mn_y   = None 
    mx_y   = None 
    mn_val = None
    mx_val = None
    
    ## loop over the 1st histo bins 
    for ix , iy , x , y , z in h1.items() :
        
        dz = diff ( z , h2 ( x , y , **kwargs ) ) ## NOTE ORDER OF ARGUMENTS 
        if mn_val is None or dz.value() < mn_val.value() :
            mn_val = dz
            mn_x   = x.value() 
            mn_y   = y.value() 
        if mx_val is None or dz.value() > mx_val.value() :
            mx_val = dz
            mx_x   = x.value()
            mx_y   = y.value() 

    ## loop over the 2nd histo bins 
    for ix , iy , x , y , z in h2.items() :
        
        dz = diff ( h1 ( x , y , **kwargs ) , z  ) ## NOTE ORDER OF ARGUMENTS 
        if mn_val is None or dz.value() < mn_val.value() :
            mn_val = dz
            mn_x   = x.value() 
            mn_y   = y.value() 
        if mx_val is None or dz.value() > mx_val.value() :
            mx_val = dz
            mx_x   = x.value()
            mx_y   = y.value() 


    return ( mn_x , mn_y , mn_val) , ( mx_x , mx_y , mx_val )
                
ROOT.TH2D.cmp_minmax = _h2_cmp_minmax_
ROOT.TH2F.cmp_minmax = _h2_cmp_minmax_ 
        

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
        h1_ = h1.density() if hasattr ( h1 , 'density' ) else h1.clone() 
        h2_ = h2.density() if hasattr ( h2 , 'density' ) else h2.clone() 
        cmp = _h2_cmp_prnt_ ( h1_ , h2_ , head1 = head1 , head2 = head2 , title = title , density = False )
        del h1_
        del h2_
        return cmp

    if not head1 : head1 = h1.GetName() 
    if not head2 : head2 = h1.GetName()

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


ROOT.TH2D.cmp_prnt = _h2_cmp_prnt_
ROOT.TH2F.cmp_prnt = _h2_cmp_prnt_ 

# =============================================================================
_decorated_classes_ = (
    ROOT.TH1  , 
    ROOT.TH1F , 
    ROOT.TH1D ,
    ROOT.TH2F , 
    ROOT.TH2D ,
    ROOT.TH3F , 
    ROOT.TH3D ,
    )
_new_methods_       = (
    ROOT.TH1D.is_constant ,
    ROOT.TH1F.is_constant , 
    ROOT.TH1D.cmp_fit     , 
    ROOT.TH1F.cmp_fit     , 
    ROOT.TH1D.cmp_pdf     , 
    ROOT.TH1F.cmp_pdf     , 
    ROOT.TH1D.cmp_chi2    , 
    ROOT.TH1F.cmp_chi2    , 
    ROOT.TH1D.chi2_cmp    , 
    ROOT.TH1F.chi2_cmp    , 
    ROOT.TH1D.cmp_cos     , 
    ROOT.TH1F.cmp_cos     ,
    #
    ROOT.TH1D.cmp_dist2   , 
    ROOT.TH1F.cmp_dist2   ,
    #
    ROOT.TH1D.cmp_dist    , 
    ROOT.TH1F.cmp_dist    ,
    ROOT.TH2D.cmp_dist    , 
    ROOT.TH2F.cmp_dist    , 
    ROOT.TH3D.cmp_dist    , 
    ROOT.TH3F.cmp_dist    , 
    #
    ROOT.TH1D.cmp_prnt   , 
    ROOT.TH1F.cmp_prnt   ,
    ROOT.TH2D.cmp_prnt   , 
    ROOT.TH2F.cmp_prnt   ,
    ##
    ROOT.TH1D.cmp_minmax ,
    ROOT.TH1F.cmp_minmax ,
    ROOT.TH2D.cmp_minmax ,
    ROOT.TH2F.cmp_minmax ,
    ##
    )
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
