#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with utilities for parameterization of historgams 
#
## (1) using fit of histograms (ROOT::TF1::Fit)
# 
# - as      Bernstein/Bezier sum 
# - as even Bernstein/Bezier sum 
# - as      Chebyshev        sum  
# - as      Legendre         sum 
# - as      Fourier          sum 
# - as      Cosine/Fourier   sum 
# - as      plain monomal    sum
# - as positive                           Bernstein/Bezier sum
# - as positive even                      Bernstein/Bezier sum
# - as positive monotonic                Bernstein/Bezier sum
# - as positive            convex/concave Bernstein/Bezier sum
# - as positive monotonic convex/concave Bernstein/Bezier sum
# - as      generic                            spline  (b-spline)
# - as      positive                           spline  (p-spline)
# - as      positive monotonic                spline  (m-spline) 
# - as      positive monotonic convex/concave spline  (c-spline) 
# - as      positive            convex/concave spline  (convex/concave-spline) 
#
# where possible, fit starts from reasonable approximatuion, taken from (1)
# 
# Typical usage:
# @code
# histo = ...                ## the historgam
# b = histo.bernstein ( 5 )  ## make a fit... 
# 
# tf1        = b[0]    ## TF1 object
# obj        = b[1]    ## helper object 
# fun        = b[2]    ## underlying normalzed C++ object 
# fit_result = b[3]    ## fit result & status
# norm       = b[4]    ## normalization coefficient 
# x = ...
# print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
# print 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) 
# @endcode
#
## (2) using RooFit pdfs:
#
# - as positive                           Bernstein/Bezier sum
# - as positive even                      Bernstein/Bezier sum
# - as positive monotonic                Bernstein/Bezier sum
# - as positive            convex/concave Bernstein/Bezier sum
# - as positive monotonic convex/concave Bernstein/Bezier sum
# - as      positive                           spline  (p-spline)
# - as      positive monotonic                spline  (m-spline) 
# - as      positive monotonic convex/concave spline  (c-spline) 
# - as      positive            convex/concave spline  (convex/concave-spline)
#
# Typical usage:
# @code
# histo = ...                                 ## the histogram
# b = histo.pdf_positive ( 5 , draw = True )  ## make a fit... 
# result = b[0]  ## RooFit result
# pdf    = b[1]  ## PDF used in the fit
# fun    = b[2]  ## the actual function from PDF
# frame  = b[3]  ## frame/RooPlot object if option 'draw' was activated
# @endcode 
#  
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date   2011-06-07  
# =============================================================================
"""Module with utilities for parameterization of historgams

## (1) using fit of histograms (ROOT::TF1::Fit)

- as      Bernstein/Bezier sum 
- as even Bernstein/Bezier sum 
- as      Chebyshev        sum  
- as      Legendre         sum 
- as      Fourier          sum 
- as      Cosine/Fourier   sum 
- as      plain monomal    sum
- as positive                           Bernstein/Bezier sum
- as positive even                      Bernstein/Bezier sum
- as positive monotonic                Bernstein/Bezier sum
- as positive            convex/concave Bernstein/Bezier sum
- as positive monotonic convex/concave Bernstein/Bezier sum
- as      generic                            spline  (b-spline)
- as      positive                           spline  (p-spline)
- as      positive monotonic                spline  (m-spline) 
- as      positive monotonic convex/concave spline  (c-spline) 
- as      positive            convex/concave spline  (convex/concave-spline) 

where possible, fit starts from reasonable approximatuion, taken from (1)

Typical usage:

>>> histo = ...                ## the historgam
>>> b = histo.bernstein ( 5 )  ## make a fit... 

>>> tf1        = b[0]    ## TF1 object
>>> obj        = b[1]    ## helper object 
>>> fun        = b[2]    ## underlying normalzed C++ object 
>>> fit_result = b[3]    ## fit result & status
>>> norm       = b[4]    ## normalization coefficient 

>>> x = ...
>>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
>>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) 

## (2) using RooFit pdfs:

- as positive                           Bernstein/Bezier sum
- as positive even                      Bernstein/Bezier sum
- as positive monotonic                Bernstein/Bezier sum
- as positive            convex/concave Bernstein/Bezier sum
- as positive monotonic convex/concave Bernstein/Bezier sum
- as      positive                           spline  (p-spline)
- as      positive monotonic                spline  (m-spline) 
- as      positive monotonic convex/concave spline  (c-spline) 
- as      positive            convex/concave spline  (convex/concave-spline)

Typical usage:

>>> histo = ...                   ## the historgam
>>> b = histo.pdf_positive ( 5 )  ## make a fit... 

>>> result = b[0]  ## RooFit result
>>> pdf    = b[1]  ## PDF used in the fit
>>> fun    = b[2]  ## the actual function from PDF
>>> frame  = b[3]  ## frame/RooPlot object if option 'draw' was activated
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = () ## nothing to import 
# =============================================================================
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.histos.param' )
else                       : logger = getLogger( __name__             )
# =============================================================================
logger.debug ( 'Some parameterization utilities for Histo objects')
# =============================================================================
from ostap.core.core     import cpp, VE, funID
from ostap.math.param    import ( legendre_sum      ,
                                  chebyshev_sum     ,
                                  fourier_sum       ,
                                  cosine_sum        ,
                                  bezier_sum        ,
                                  bernstein_sum     , 
                                  beziereven_sum    ,
                                  bernsteineven_sum )
from ostap.fitting.param import H_fit, H_Nfit 
# =============================================================================
inf_pos =  float('inf') ## positive infinity
inf_neg = -float('inf') ## negative infinity
# =============================================================================

# =============================================================================
## represent 1D-histo as polynomial sum 
def _h1_param_sum_ ( h1              ,
                     fun_obj         ,
                     fit_type        ,  
                     opts  = 'SQ0I'  ,
                     xmin  = inf_neg ,
                     xmax  = inf_pos ,
                     fixes = ()      ) :
    """ Represent histo as polynomial sum    
    """
    ## 
    b     = fun_obj  
    #
    bfit  = fit_type ( b )    
    bfit.fun.SetNpx  ( max ( 100 , 3 * h1.bins() ) )  
    
    bfit.histo     = h1
    
    xmin  = max ( xmin , h1.xmin() )
    xmax  = min ( xmax , h1.xmax() )

    ## calculate the integral in range 
    _integral_ = h1.integrate ( xmin = xmin , xmax = xmax )
    
    fun = bfit.fun

    if hasattr ( bfit , 'norm' ) and bfit.norm() : 
        fun.SetParameter ( 0 , _integral_ ) 
        from math import pi 
        for i in range( 0, b.npars() ) :
            fun.SetParameter ( i + 1 , 0  )
            fun.SetParLimits ( i + 1 , -1.75 * pi  , 1.75 * pi ) 
            
    if not opts                   : opts  = 'S'
    if 0 > opts.upper().find('S') : opts += 'S'
    
    ## fitting options:
    fopts = opts,'',xmin,xmax 

    ## fix parameters 
    for i,v in fixes :
        fun.FixParameter ( i , v )
        
    if hasattr ( bfit , 'norm' ) and bfit.norm() :
        fun.FixParameter    ( 0 , _integral_ )
        ## specific options here: 
        r  = fun.Fit( h1 , opts+'0Q', '', xmin , xmax )
        fun.ReleaseParameter(0)
        
    r = fun.Fit( h1 , *fopts )
        
    if 0 != r.Status() :
        logger.info ( 'Fit status is  %d [%s]' % ( r.Status() , type(b).__name__ ) )
        if hasattr ( bfit , 'norm' ) and bfit.norm() : 
            fun.FixParameter ( 0 , _integral_ )
            for i in range ( 0 , b.npars() ) : fun.SetParameter ( i + 1 , 0  )
            r = fun.Fit( h1 , *fopts )
            fun.ReleaseParameter ( 0 )
        r = fun.Fit( h1 , *fopts )
            
    if 0 != r.Status() :
        logger.warning ( 'Fit status is  %d [%s]' % ( r.Status() , type(b).__name__ ) )
        r = fun.Fit( h1, *fopts )        
        if hasattr ( bfit , 'norm' ) and bfit.norm() : 
            fun.FixParameter ( 0 , _integral_ )
            for i in range( 0 , b.npars() ) : fun.SetParameter ( i + 1 , 0  )
            r = fun.Fit(h1, *fopts)
            fun.ReleaseParameter( 0 )
        r = fun.Fit ( h1 , *fopts )
        if 0 != r.Status() :
            logger.error ('Fit status is  %d [%s]' % ( r.Status() , type(b).__name__ ) )

    bfit.fitresult = r
    
    norm = VE(1,0)
    if hasattr ( bfit , 'norm' ) and bfit.norm() : 
        bfit.fitnorm = r[0]
        norm         = r[0]
        
    return bfit.fun , bfit , b , bfit.fitresult, norm 

# =============================================================================
## represent 1D-histo as Bernstein polynomial
#  @code
#  h = ...                   ## the historgam
#  b = h.bernstein ( 5 )     ## make a fit... 
# 
#  tf1        = b[0]         ## TF1 object
#  obj        = b[1]         ## helper object 
#  fun        = b[2]         ## underlying normalzed C++ object 
#  fit_result = b[3]         ## fit result & status
#  norm       = b[4]         ## normalization
# 
#  x = ...
#  print 'TF1(%s) = %s' % ( x , tf1 ( x )        ) 
#  print 'fun(%s) = %s' % ( x , fun ( x ) * norm )
#  @endcode 
def _h1_bernstein_ ( h1 , degree , opts = 'SQ0' , xmin = inf_neg , xmax = inf_pos , fixes = () ) :
    """Represent histo as Bernstein polynomial
    
    >>> h = ...                # the historgam
    >>> b = h.bernstein ( 5 )  ## make a fit... 
    
    >>> tf1        = b[0]    ## TF1 object
    >>> obj        = b[1]    ## helper object 
    >>> fun        = b[2]    ## underlying normalzed C++ object 
    >>> fit_result = b[3]    ## fit result & status
    >>> norm       = b[4]    ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) 
    """
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    # make reasonable approximation
    func  = bezier_sum ( h1   , degree , xmin , xmax )
    #
    return _h1_param_sum_ ( h1 , func , H_fit , opts , xmin , xmax , fixes )  


# =============================================================================
## represent 1D-histo as even Bernstein polynomial
#  @code
#  h = ...                   ## the historgam
#  b = h.bernsteineven ( 3 ) ## make a fit... 
# 
#  tf1        = b[0]         ## TF1 object
#  obj        = b[1]         ## helper object 
#  fun        = b[2]         ## underlying normalzed C++ object 
#  fit_result = b[3]         ## fit result & status
#  norm       = b[4]         ## normalization
# 
#  x = ...
#  print 'TF1(%s) = %s' % ( x , tf1 ( x )        ) 
#  print 'fun(%s) = %s' % ( x , fun ( x ) * norm )
#  @endcode 
def _h1_bernsteineven_ ( h1 , halfdegree , opts = 'SQ0' , xmin = inf_neg , xmax = inf_pos , fixes = () ) :
    """Represent histo as even Bernstein polynomial
    
    >>> h = ...                    ## the historgam
    >>> b = h.bernsteineven ( 3 )  ## make a fit... 
    
    >>> tf1        = b[0]    ## TF1 object
    >>> obj        = b[1]    ## helper object 
    >>> fun        = b[2]    ## underlying normalzed C++ object 
    >>> fit_result = b[3]    ## fit result & status
    >>> norm       = b[4]    ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) 
    """
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    # make reasonable approximation
    func  = beziereven_sum ( h1   , halfdegree , xmin , xmax )
    # make a fit 
    return _h1_param_sum_ ( h1 , func , H_fit , opts , xmin , xmax , fixes )  


# =============================================================================
## represent 1D-histo as Chebyshev polynomial
#  @code
#  h = ...                  ## the historgam
#  b = h.chebyshev ( 5 )     ## make a fit... 
# 
#  tf1        = b[0]        ## TF1 object
#  obj        = b[1]        ## helper object 
#  fun        = b[2]        ## underlying normalzed C++ object 
#  fit_result = b[3]        ## fit result & status
#  norm       = b[4]        ## normalization
# 
#  x = ...
#  print 'TF1(%s) = %s' % ( x , tf1 ( x )        ) 
#  print 'fun(%s) = %s' % ( x , fun ( x ) * norm )
#  @endcode 
def _h1_chebyshev_ ( h1 , degree , opts = 'SQ0I' , xmin = inf_neg , xmax = inf_pos ) :
    """Represent histo as Chebyshev sum 
    
    >>> h = ... # the historgam
    >>> b = h.chebyshev ( 5 )  ## make a fit... 

    >>> tf1        = b[0]    ## TF1 object
    >>> obj        = b[1]    ## helper object 
    >>> fun        = b[2]    ## underlying normalzed C++ object 
    >>> fit_result = b[3]    ## fit result & status
    >>> norm       = b[4]    ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) 
    """
    #
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    ## make reasonable approximation: 
    func = chebyshev_sum ( h1 , degree ,  xmin , xmax )
    ## fit it!
    return _h1_param_sum_ ( h1 , func , H_fit , opts , xmin , xmax )  

# =============================================================================
## represent 1D-histo as Legendre polynomial
#  @code
#  h = ...                  ## the historgam
#  b = h.legendre ( 5 )     ## make a fit... 
# 
#  tf1        = b[0]        ## TF1 object
#  obj        = b[1]        ## helper object 
#  fun        = b[2]        ## underlying normalzed C++ object 
#  fit_result = b[3]        ## fit result & status
#  norm       = b[4]        ## normalization
# 
#  x = ...
#  print 'TF1(%s) = %s' % ( x , tf1 ( x )        ) 
#  print 'fun(%s) = %s' % ( x , fun ( x ) * norm )
#  @endcode 
def _h1_legendre_ ( h1 , degree , opts = 'SQ0' , xmin = inf_neg , xmax = inf_pos ) :
    """Represent histo as Legendre sum 
    
    >>> h = ... # the historgam
    >>> b = h.legendre ( 5 )  ## make a fit...

    >>> tf1        = b[0]    ## TF1 object
    >>> obj        = b[1]    ## helper object 
    >>> fun        = b[2]    ## underlying normalzed C++ object 
    >>> fit_result = b[3]    ## fit result & status
    >>> norm       = b[4]    ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) 
    """
    ## make reasonable approximation:
    mn,mx  = h1.xminmax()
    ##
    if 1 > h1.GetEntries() : sw = 1
    else :
        vmn,vmx = h1.minmax()     
        sw      = max ( abs ( h1.GetSumOfWeights() ) ,
                        abs ( h1.Integral()        ) , abs ( vmn ) , abs ( vmx ) ) 
    ##
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    func  = legendre_sum   ( h1 , degree , xmin , xmax , 
                             epsabs = 1.e-4 * sw    ,
                             epsrel = 1.e-3         ,
                             limit  = 3 * h1.bins() ) 
    ## fit it!
    return _h1_param_sum_ ( h1 , func , H_fit , opts , xmin , xmax )  


# =============================================================================
## represent 1D-histo as Fourier polynomial
#  @code
#  h = ...                  ## the historgam
#  b = h.fourier ( 3 )      ## make a fit... 
# 
#  tf1        = b[0]        ## TF1 object
#  obj        = b[1]        ## helper object 
#  fun        = b[2]        ## underlying normalzed C++ object 
#  fit_result = b[3]        ## fit result & status
#  norm       = b[4]        ## normalization
# 
#  x = ...
#  print 'TF1(%s) = %s' % ( x , tf1 ( x )        ) 
#  print 'fun(%s) = %s' % ( x , fun ( x ) * norm )
#  @endcode 
def _h1_fourier_ ( h1 , degree , fejer = False , opts = 'SQ0I' , xmin = inf_neg , xmax = inf_pos ) :
    """Represent histo as Fourier sum 
    
    >>> h = ... # the historgam
    >>> b = h.fourier ( 3 )  ## make a fit... 
        
    >>> tf1        = b[0]    ## TF1 object
    >>> obj        = b[1]    ## helper object 
    >>> fun        = b[2]    ## underlying normalzed C++ object 
    >>> fit_result = b[3]    ## fit result & status
    >>> norm       = b[4]    ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) 
    """
    #

    ## make reasonable approximation:
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    func = fourier_sum ( h1 , degree , xmin , xmax , fejer )

    ## fit it!
    return _h1_param_sum_ ( h1 , func , H_fit , opts , xmin , xmax )  

# =============================================================================
## represent 1D-histo as cosine Fourier polynomial
#  @code
#  h = ...                  ## the historgam
#  b = h.cosine ( 3 )       ## make a fit... 
# 
#  tf1        = b[0]        ## TF1 object
#  obj        = b[1]        ## helper object 
#  fun        = b[2]        ## underlying normalzed C++ object 
#  fit_result = b[3]        ## fit result & status
#  norm       = b[4]        ## normalization
# 
#  x = ...
#  print 'TF1(%s) = %s' % ( x , tf1 ( x )        ) 
#  print 'fun(%s) = %s' % ( x , fun ( x ) * norm )
#  @endcode 
def _h1_cosine_ ( h1 , degree , fejer = False , opts = 'SQ0I' , xmin = inf_neg , xmax = inf_pos ) :
    """Represent histo as Cosine Fourier sum 
    
    >>> h = ... # the historgam
    >>> b = h.cosine ( 3 )  ## make a fit... 
    
    >>> tf1        = b[0]    ## TF1 object
    >>> obj        = b[1]    ## helper object 
    >>> fun        = b[2]    ## underlying normalzed C++ object 
    >>> fit_result = b[3]    ## fit result & status
    >>> norm       = b[4]    ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) 
    """
    #

    ## make reasonable approximation:
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    func = cosine_sum ( h1 , degree , xmin , xmax , fejer )

    ## fit it!
    return _h1_param_sum_ ( h1 , func , H_fit , opts , xmin , xmax )  

# =============================================================================
## represent 1D-histo as plain vanilla polynomial
#  @code
#  h = ...                  ## the historgam
#  b = h.polinomial ( 5 )  ## make a fit... 
# 
#  tf1        = b[0]        ## TF1 object
#  obj        = b[1]        ## helper object 
#  fun        = b[2]        ## underlying normalzed C++ object 
#  fit_result = b[3]        ## fit result & status
#  norm       = b[4]        ## normalization
# 
#  x = ...
#  print 'TF1(%s) = %s' % ( x , tf1 ( x )        ) 
#  print 'fun(%s) = %s' % ( x , fun ( x ) * norm )
#  @endcode 
def _h1_polinomial_ ( h1 , degree , opts = 'SQ0' , xmin = inf_neg , xmax = inf_pos ) :
    """Represent histo as plain vanilla polynomial    
    >>> h = ... # the historgam    
    >>> b = h.polinomial ( 5 )  ## make a fit... 
    
    >>> tf1        = b[0]    ## TF1 object
    >>> obj        = b[1]    ## helper object 
    >>> fun        = b[2]    ## underlying normalzed C++ object 
    >>> fit_result = b[3]    ## fit result & status
    >>> norm       = b[4]    ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) 
    """
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    func  = cpp.Ostap.Math.Polynomial ( degree , xmin , xmax ) 
    #
    my = h1.accumulate().value()/h1.bins()
    func.setPar( 0, my )
    ##
    return _h1_param_sum_ ( h1 , func , H_fit , opts , xmin , xmax )  


# =============================================================================
## represent 1D-histo as B-spline
#  @code
#  h = ...                  ## the historgam
#  b = h.bSpline ( degree = 3 , innerknots = 3  )
#  b = h.bSpline ( degree = 3 , innerknots = [ 0.1 , 0.2, 0.8, 0.9 ]  )
# 
#  tf1        = b[0]        ## TF1 object
#  obj        = b[1]        ## helper object 
#  fun        = b[2]        ## underlying normalzed C++ object 
#  fit_result = b[3]        ## fit result & status
#  norm       = b[4]        ## normalization
# 
#  x = ...
#  print 'TF1(%s) = %s' % ( x , tf1 ( x )        ) 
#  print 'fun(%s) = %s' % ( x , fun ( x ) * norm )
#  @endcode 
def _h1_bspline_ ( h1 , degree = 3 , knots = 3 , opts = 'SQ0' , xmin = inf_neg , xmax = inf_pos ) :
    """Represent histo as B-spline polynomial    
    >>> h = ... # the historgam
    >>> b = h.bSpline ( degree = 3 , innerknots = 3  )
    >>> b = h.bSpline ( degree = 3 , innerknots = [ 0.1 , 0.2, 0.8, 0.9 ]  )
    
    >>> tf1        = b[0]    ## TF1 object
    >>> obj        = b[1]    ## helper object 
    >>> fun        = b[2]    ## underlying normalzed C++ object 
    >>> fit_result = b[3]    ## fit result & status
    >>> norm       = b[4]    ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) 
    
    """
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    #
    if isinstance ( knots , ( int , long ) ) :
        func = cpp.Ostap.Math.BSpline ( xmin , xmax , knots , degree )
    else :
        from ostap.math.base import doubles
        _knots = doubles ( xmin , xmax ) 
        for k in knots : _knots.push_back( k )
        func = cpp.Ostap.Math.BSpline ( _knots , degree )
        
    ##
    return _h1_param_sum_ ( h1 , func , H_fit , opts , xmin , xmax )  


# =============================================================================
## represent 1D-histo as POSITIVE bernstein polynomial
#  @code
#  h = ...                  ## the historgam
#  b = h.positive     ( 5 ) ## 5 is degree
# 
#  tf1        = b[0]        ## TF1 object
#  obj        = b[1]        ## helper object 
#  fun        = b[2]        ## underlying normalzed C++ object 
#  fit_result = b[3]        ## fit result & status
#  norm       = b[4]        ## normalization
# 
#  x = ...
#  print 'TF1(%s) = %s' % ( x , tf1 ( x )        ) 
#  print 'fun(%s) = %s' % ( x , fun ( x ) * norm )
#  @endcode 
def _h1_positive_ ( h1 , N , opts = 'SQ0' , xmin = inf_neg , xmax = inf_pos  ) :
    """Represent histo as Positive Bernstein polynomial
    >>> h = ...              ## the historgam
    >>> b = h.positive ( 5 ) ## 5 is degree
    
    >>> tf1        = b[0]    ## TF1 object
    >>> obj        = b[1]    ## helper object 
    >>> fun        = b[2]    ## underlying normalzed C++ object 
    >>> fit_result = b[3]    ## fit result & status
    >>> norm       = b[4]    ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) 
    """
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    func = cpp.Ostap.Math.Positive ( N , xmin , xmax )
    # 
    return _h1_param_sum_ ( h1 , func , H_Nfit , opts , xmin , xmax ) 

# =============================================================================
## represent 1D-histo as POSITIVE EVEN bernstein polynomial
#  @code
#  h = ...                  ## the historgam
#  b = h.positiveeven ( 2 ) ## 2 is half-degree
# 
#  tf1        = b[0]        ## TF1 object
#  obj        = b[1]        ## helper object 
#  fun        = b[2]        ## uderlying normalzed C++ object 
#  fit_result = b[3]        ## fit result & status
#  norm       = b[4]        ## normalization
# 
#  x = ...
#  print 'TF1(%s) = %s' % ( x , tf1 ( x )        ) 
#  print 'fun(%s) = %s' % ( x , fun ( x ) * norm )
#  @endcode 
def _h1_positiveeven_ ( h1 , N , opts = 'SQ0' , xmin = inf_neg , xmax = inf_pos  ) :
    """Represent histo as Positive Even Bernstein polynomial
    
    >>> h = ... # the historgam
    >>> b = h.positiveeven ( 2 ) ## 2 is half-degree
    
    >>> tf1        = b[0] ## TF1 object
    >>> obj        = b[1] ## helper object 
    >>> fun        = b[2] ## underlying normalzed C++ object 
    >>> fit_result = b[3] ## fit result & status
    >>> norm       = b[4] ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) 
    """
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    func = cpp.Ostap.Math.PositiveEven ( N , xmin , xmax )
    # 
    return _h1_param_sum_ ( h1 , func , H_Nfit , opts , xmin , xmax ) 


# =============================================================================
## represent 1D-histo as MONOTONIC bernstein polynomial
#  @code
#  h = ...                  ## the historgam
#  b = h.monotonic ( 5 , increasing = True )
# 
#  tf1        = b[0]        ## TF1 object
#  obj        = b[1]        ## helper object 
#  fun        = b[2]        ## uderlying normalzed C++ object 
#  fit_result = b[3]        ## fit result & status
#  norm       = b[4]        ## normalization
# 
#  x = ...
#  print 'TF1(%s) = %s' % ( x , tf1 ( x )        ) 
#  print 'fun(%s) = %s' % ( x , fun ( x ) * norm )
#  @endcode 
def _h1_monotonic_ ( h1 , N , increasing = True , opts = 'SQ0' , xmin = inf_neg , xmax = inf_pos  ) :
    """Represent histo as Monotonic Bernstein polynomial
    
    >>> h = ...           ## the historgam
    >>> b = h.monotonic ( 5 , increasing = True )
    
    >>> tf1        = b[0] ## TF1 object
    >>> obj        = b[1] ## helper object 
    >>> fun        = b[2] ## underlying normalzed C++ object 
    >>> fit_result = b[3] ## fit result & status
    >>> norm       = b[4] ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) 
    """
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    func  = cpp.Ostap.Math.Monotonic ( N , xmin , xmax , increasing )
    # 
    return _h1_param_sum_ ( h1 , func , H_Nfit , opts , xmin , xmax ) 


# =============================================================================
## represent 1D-histo as MONOTONIC CONVEX/CONCAVE bernstein polynomial
#  @code
#  h = ...                  ## the historgam
#  b = h.convex ( 5 , increasing = True , convex = False )    
# 
#  tf1        = b[0]        ## TF1 object
#  obj        = b[1]        ## helper object 
#  fun        = b[2]        ## uderlying normalzed C++ object 
#  fit_result = b[3]        ## fit result & status
#  norm       = b[4]        ## normalization
# 
#  x = ...
#  print 'TF1(%s) = %s' % ( x , tf1 ( x )        ) 
#  print 'fun(%s) = %s' % ( x , fun ( x ) * norm )
#  @endcode 
def _h1_convex_ ( h1 , N , increasing = True , convex = True , opts = 'SQ0' ,  xmin = inf_neg , xmax = inf_pos ) :
    """Represent histo as Monotonic Convex/Concave  Bernstein polynomial    
    >>> h = ...           ## the historgam
    >>> b = h.convex ( 5 , increasing = True , convex = False )    
    
    >>> tf1        = b[0] ## TF1 object
    >>> obj        = b[1] ## helper object 
    >>> fun        = b[2] ## underlying normalzed C++ object 
    >>> fit_result = b[3] ## fit result & status
    >>> norm       = b[4] ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) 
    >>> h = ... # the historgam
    """
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    func = cpp.Ostap.Math.Convex ( N , xmin , xmax , increasing , convex )
    # 
    return _h1_param_sum_ ( h1 , func , H_Nfit , opts , xmin , xmax ) 


# =============================================================================
## represent 1D-histo as CONVEX bernstein polynomial
#  @code
#  h = ...                  ## the historgam
#  b = h.convexpoly ( 5 )    
# 
#  tf1        = b[0]        ## TF1 object
#  obj        = b[1]        ## helper object 
#  fun        = b[2]        ## uderlying normalzed C++ object 
#  fit_result = b[3]        ## fit result & status
#  norm       = b[4]        ## normalization
# 
#  x = ...
#  print 'TF1(%s) = %s' % ( x , tf1 ( x )        ) 
#  print 'fun(%s) = %s' % ( x , fun ( x ) * norm )
#  @endcode 
def _h1_convexpoly_ ( h1 , N , opts = 'SQ0' , xmin = inf_neg , xmax = inf_pos ) :
    """Represent histo as Convex Bernstein polynomial
    >>> h = ...           ## the historgam
    >>> b = h.convexpoly ( 5 )    
    
    >>> tf1        = b[0] ## TF1 object
    >>> obj        = b[1] ## helper object 
    >>> fun        = b[2] ## underlying normalzed C++ object 
    >>> fit_result = b[3] ## fit result & status
    >>> norm       = b[4] ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) 
    """
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    func = cpp.Ostap.Math.ConvexOnly ( N , xmin , xmax , True )
    # 
    return _h1_param_sum_ ( h1 , func , H_Nfit , opts , xmin , xmax ) 

# =============================================================================
## represent 1D-histo as CONCAVE bernstein polynomial
#  @code
#  h = ...                  ## the historgam
#  b = h.concavepoly ( 5 )    
# 
#  tf1        = b[0]        ## TF1 object
#  obj        = b[1]        ## helper object 
#  fun        = b[2]        ## uderlying normalzed C++ object 
#  fit_result = b[3]        ## fit result & status
#  norm       = b[4]        ## normalization
# 
#  x = ...
#  print 'TF1(%s) = %s' % ( x , tf1 ( x )        ) 
#  print 'fun(%s) = %s' % ( x , fun ( x ) * norm )
#  @endcode 
def _h1_concavepoly_ ( h1 , N , opts = 'SQ0' , xmin = inf_neg , xmax = inf_pos ) :
    """Represent histo as Concave  Bernstein polynomial

    >>> h = ...           ## the historgam
    >>> b = h.concavepoly ( 5 )    
    
    >>> tf1        = b[0] ## TF1 object
    >>> obj        = b[1] ## helper object 
    >>> fun        = b[2] ## underlying normalzed C++ object 
    >>> fit_result = b[3] ## fit result & status
    >>> norm       = b[4] ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) 
    """
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    func = cpp.Ostap.Math.ConvexOnly ( N , xmin , xmax , False )
    # 
    return _h1_param_sum_ ( h1 , func , H_Nfit , opts , xmin , xmax ) 

# =============================================================================
## represent 1D-histo as positive B-spline
#  @code
#  h = ...                  ## the historgam
#  b  = h.pSpline ( degree = 3 , knots = 3  )
#  b  = h.pSpline ( degree = 3 , knots = [ 0.1 , 0.2, 0.8, 0.9 ]  )
# 
#  tf1        = b[0]        ## TF1 object
#  obj        = b[1]        ## helper object 
#  fun        = b[2]        ## uderlying normalzed C++ object 
#  fit_result = b[3]        ## fit result & status
#  norm       = b[4]        ## normalization
# 
#  x = ...
#  print 'TF1(%s) = %s' % ( x , tf1 ( x )        ) 
#  print 'fun(%s) = %s' % ( x , fun ( x ) * norm )
#  @endcode 
def _h1_pspline_ ( h1 , degree = 3 , knots = 3 , opts = 'SQ0I' , xmin = inf_neg , xmax = inf_pos ) :
    """Represent histo as positive B-spline 
    
    >>> h  = ... # the historgam
    >>> b  = h.pSpline ( degree = 3 , knots = 3  )
    >>> b  = h.pSpline ( degree = 3 , knots = [ 0.1 , 0.2, 0.8, 0.9 ]  )

    >>> tf1        = b[0] ## TF1 object
    >>> obj        = b[1] ## helper object 
    >>> fun        = b[2] ## underlying normalzed C++ object 
    >>> fit_result = b[3] ## fit result & status
    >>> norm       = b[4] ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) )     
    """
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    #
    if isinstance ( knots , ( int , long ) ) :
        func = cpp.Ostap.Math.PositiveSpline ( xmin , xmax , knots , degree  )
    else :
        from ostap.math.base import doubles
        _knots = doubles ( xmin , xmax ) 
        for k in knots : _knots.push_back( k )
        func = cpp.Ostap.Math.PositiveSpline ( _knots , degree )
    #
    return _h1_param_sum_ ( h1 , func , H_Nfit , opts , xmin , xmax ) 

# =============================================================================
## represent 1D-histo as positive monotonic spline
#  @code
#  h = ...                  ## the historgam
#  b = h.mSpline ( degree = 3 , knots = 3  , increasing = True  )
#  b = h.mSpline ( degree = 3 , knots = [ 0.1 , 0.2, 0.8, 0.9 ] )
#  b = h.mSpline ( degree = 3 , knots = 3  , increasing = False )
# 
#  tf1        = b[0]        ## TF1 object
#  obj        = b[1]        ## helper object 
#  fun        = b[2]        ## uderlying normalzed C++ object 
#  fit_result = b[3]        ## fit result & status
#  norm       = b[4]        ## normalization
# 
#  x = ...
#  print 'TF1(%s) = %s' % ( x , tf1 ( x )        ) 
#  print 'fun(%s) = %s' % ( x , fun ( x ) * norm )
#  @endcode 
def _h1_mspline_ ( h1 , degree = 3 , knots = 3 , increasing = True , opts = 'SQ0I' , xmin = inf_neg , xmax = inf_pos ) :
    """Represent histo as positive monotonic  spline 
    
    >>> h  = ... # the historgam
    
    >>> b  = h.mSpline ( degree = 3 , knots = 3  , increasing = True  )
    >>> b  = h.mSpline ( degree = 3 , knots = [ 0.1 , 0.2, 0.8, 0.9 ] )
    >>> b  = h.mSpline ( degree = 3 , knots = 3  , increasing = False )
    
    >>> tf1        = b[0] ## TF1 object
    >>> obj        = b[1] ## helper object 
    >>> fun        = b[2] ## underlying normalzed C++ object 
    >>> fit_result = b[3] ## fit result & status
    >>> norm       = b[4] ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) )     
    
    """
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    #
    if isinstance ( knots , ( int , long ) ) :
        func = cpp.Ostap.Math.MonotonicSpline ( xmin , xmax , knots , degree , increasing )
    else :
        from ostap.math.base import doubles
        _knots = doubles ( xmin , xmax ) 
        for k in knots : _knots.push_back( k )
        func = cpp.Ostap.Math.MonotonicSpline ( knots , degree , increasing )
    #
    return _h1_param_sum_ ( h1 , func , H_Nfit , opts , xmin , xmax ) 

# =============================================================================
## represent 1D-histo as monotonic convex/concave spline
#  @code
#  h = ...                  ## the historgam
#  b = h.cSpline ( degree = 3 , knots = 3  , increasing = True , convex = False )
#  b = h.cSpline ( degree = 3 , knots = [ 0.1 , 0.2, 0.8, 0.9 ]  )
# 
#  tf1        = b[0]        ## TF1 object
#  obj        = b[1]        ## helper object 
#  fun        = b[2]        ## uderlying normalzed C++ object 
#  fit_result = b[3]        ## fit result & status
#  norm       = b[4]        ## normalization
# 
#  x = ...
#  print 'TF1(%s) = %s' % ( x , tf1 ( x )        ) 
#  print 'fun(%s) = %s' % ( x , fun ( x ) * norm )
#  @endcode 
def _h1_cspline_ ( h1                   ,
                   degree     = 3       ,
                   knots      = 3       ,
                   increasing = True    ,
                   convex     = True    , 
                   opts       = 'SQ0I'  ,
                   xmin       = inf_neg ,
                   xmax       = inf_pos ) :
    """Represent histo as positive monotonic convex/concave spline  
    
    >>> h = ... # the historgam
    >>> b = h.cSpline ( degree = 3 , knots = 3  , increasing = True , convex = False )
    >>> b = h.cSpline ( degree = 3 , knots = [ 0.1 , 0.2, 0.8, 0.9 ]  )

    >>> tf1        = b[0] ## TF1 object
    >>> obj        = b[1] ## helper object 
    >>> fun        = b[2] ## underlying normalzed C++ object 
    >>> fit_result = b[3] ## fit result & status
    >>> norm       = b[4] ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) )     
    """
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    #
    if isinstance ( knots , ( int , long ) ) :
        func = cpp.Ostap.Math.ConvexSpline ( xmin , xmax , knots , degree , increasing , convex  )
    else :
        from ostap.math.base import doubles
        _knots = doubles ( xmin , xmax ) 
        for k in knots : _knots.push_back( k )
        func   = cpp.Ostap.Math.ConvexSpline ( _knots , order , increasing , convex )
        
    return _h1_param_sum_ ( h1 , func , H_Nfit , opts ) 


# =============================================================================
## represent 1D-histo as positive convex spline
#  @code
#  h = ...                  ## the historgam
#  b = h.convexSpline ( degree = 3 , knots = 3 )
#  b = h.convexSpline ( degree = 3 , knots = [ 0.1 , 0.2, 0.8, 0.9 ]  )
# 
#  tf1        = b[0]        ## TF1 object
#  obj        = b[1]        ## helper object 
#  fun        = b[2]        ## uderlying normalzed C++ object 
#  fit_result = b[3]        ## fit result & status
#  norm       = b[4]        ## normalization
# 
#  x = 10
#  print 'TF1(%s) = %s' % ( x , tf1 ( x )        ) 
#  print 'fun(%s) = %s' % ( x , fun ( x ) * norm )
#  @endcode 
def _h1_convexspline_ ( h1                 ,
                        degree   = 3       ,
                        knots    = 3       ,
                        opts     = 'SQ0I'  ,
                        xmin     = inf_neg ,
                        xmax     = inf_pos ) :
    """Represent histo as positive convex spline  
    
    >>> h = ... # the historgam
    >>> b = h.convexSpline ( degree = 3 , knots = 3 )
    >>> b = h.convexSpline ( degree = 3 , knots = [ 0.1 , 0.2, 0.8, 0.9 ]  )
    
    >>> tf1        = b[0] ## TF1 object
    >>> obj        = b[1] ## helper object 
    >>> fun        = b[2] ## underlying normalzed C++ object 
    >>> fit_result = b[3] ## fit result & status
    >>> norm       = b[4] ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) )     
    """
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    #
    if isinstance ( knots , ( int , long ) ) :
        func = cpp.Ostap.Math.ConvexOnlySpline ( xmin , xmax , knots , degree , True )
    else :
        from ostap.math.base import doubles
        _knots = doubles ( mn , mx ) 
        for k in knots : _knots.push_back( k )
        func   = cpp.Ostap.Math.ConvexOnlySpline ( _knots , order , True )
        
    return _h1_param_sum_ ( h1 , func , H_Nfit , opts , xmin , xmax ) 

# =============================================================================
## represent 1D-histo as positive concave spline
#  @code
#  h = ...                  ## the historgam
#  b = h.concaveSpline ( degree = 3 , knots = 3 )
#  b = h.concaveSpline ( degree = 3 , knots = [ 0.1 , 0.2, 0.8, 0.9 ]  )
#
#  tf1        = b[0]        ## TF1 object
#  obj        = b[1]        ## helper object 
#  fun        = b[2]        ## uderlying normalzed C++ object 
#  fit_result = b[3]        ## fit result & status
#  norm       = b[4]        ## normalization
# 
#  x = ...
#  print 'TF1(%s) = %s' % ( x , tf1 ( x )        ) 
#  print 'fun(%s) = %s' % ( x , fun ( x ) * norm )
#  @endcode 
def _h1_concavespline_ ( h1               ,
                         degree = 3       ,
                         knots  = 3       ,
                         opts   = 'SQ0I'  ,
                         xmin   = inf_neg ,
                         xmax   = inf_pos ) :
    """Represent histo as positive convcave spline  
    
    >>> h = ... # the historgam
    >>> b = h.concaveSpline ( degree = 3 , knots = 3 )
    >>> b = h.concaveSpline ( degree = 3 , knots = [ 0.1 , 0.2, 0.8, 0.9 ]  )
    
    >>> tf1        = b[0] ## TF1 object
    >>> obj        = b[1] ## helper object 
    >>> fun        = b[2] ## underlying normalzed C++ object 
    >>> fit_result = b[3] ## fit result & status
    >>> norm       = b[4] ## normalization 
    
    >>> x = ... 
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) )         
    """
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    #
    #
    if isinstance ( knots , ( int , long ) ) :
        func = cpp.Ostap.Math.ConvexOnlySpline ( xmin , xmax , knots , degree , False )
    else :
        from ostap.math.base import doubles
        _knots = doubles ( xmin , xmax ) 
        for k in knots : _knots.push_back( k )
        func   = cpp.Ostap.Math.ConvexOnlySpline ( _knots , order , False )
        
    return _h1_param_sum_ ( h1 , func , H_Nfit , opts ) 

# =============================================================================

_new_methods_ = []
# =============================================================================
## decorate histograms 
for t in ( ROOT.TH1D , ROOT.TH1F ) :
    t.bernstein      = _h1_bernstein_
    t.bezier         = _h1_bernstein_ ## ditto 
    t.bernsteineven  = _h1_bernsteineven_
    t.beziereven     = _h1_bernsteineven_ ## ditto 
    t.chebyshev      = _h1_chebyshev_
    t.legendre       = _h1_legendre_
    t.fourier        = _h1_fourier_
    t.cosine         = _h1_cosine_
    t.polynomial     = _h1_polinomial_
    t.positive       = _h1_positive_
    t.positiveeven   = _h1_positiveeven_
    t.monotonic     = _h1_monotonic_
    t.convex         = _h1_convex_
    t.convexpoly     = _h1_convexpoly_
    t.concavepoly    = _h1_concavepoly_
    t.bSpline        = _h1_bspline_
    t.pSpline        = _h1_pspline_
    t.mSpline        = _h1_mspline_
    t.cSpline        = _h1_cspline_
    t.convexSpline   = _h1_convexspline_ 
    t.concaveSpline  = _h1_concavespline_ 
    t.convexspline   = _h1_convexspline_ 
    t.concavespline  = _h1_concavespline_ 

    _new_methods_ += [
        _h1_bernstein_     ,
        _h1_bernsteineven_ ,
        _h1_chebyshev_     ,
        _h1_legendre_      ,
        _h1_fourier_       ,
        _h1_polinomial_    ,
        _h1_positive_      ,
        _h1_positiveeven_  ,
        _h1_monotonic_    ,
        _h1_convex_        ,
        _h1_convexpoly_    ,
        _h1_concavepoly_   ,
        _h1_bspline_       ,
        _h1_pspline_       ,
        _h1_mspline_       ,
        _h1_cspline_       ,
        _h1_convexspline_  ,
        _h1_concavespline_ ,
        ]
    
# =============================================================================
## parameterize positive histogram with certain PDF
def _h1_pdf_ ( h1 , pdf_type , pars , *args, **kwargs ) :
    """Parameterize positive histogram with certain PDF
    """
    ##
    mn,mx = h1.minmax()
    if mn.value() < 0 or mx.value() <= 0 :
        raise AttributeError("Histo goes to negative %s/%s" % ( mn , mx ) )
    ##
    if not hasattr ( h1 , 'xvar' ) :
        h1.xvar = ROOT.RooRealVar ( 'x' + h1.GetName() , 'xvar(%s)' % h1.GetName() , *h1.xminmax() )

    ## create pdf 
    pdf  = pdf_type     ( 'pdf_' + h1.GetName() , h1.xvar , *pars )
    ## fit the histogram 
    r,f  = pdf.fitHisto ( h1 , *args, **kwargs )
    ##
    func = pdf.pdf.function()
    ##
    from ostap.fitting.roofit import PDF_fun
    pdf_fun = PDF_fun( pdf.pdf , h1.xvar , *h1.xminmax() )
    ##
    norm = VE(  h1.integrate().value() , 0 ) 
    return r , pdf , func , norm , pdf_fun, f  

# =============================================================================
## parameterize/fit histogram with the positive polynomial
#  @code
#  h1 = ...
#  results = h1.pdf_positive ( 3 )
#  results = h1.pdf_positive ( 3 , draw = True , silent = True )
#  print     results[ 0]
#  pdf     = results[ 1]
#  func    = results[ 2]
#  norm    = results[ 3]
#  pdf_fun = results[ 4] 
#  frame   = results[-1] 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-26
def _h1_pdf_positive_ ( h1 , degree , *args , **kwargs ) :
    """Parameterize/fit histogram with the positive polynomial
    >>> h1 = ...
    >>> results = h1.pdf_positive ( 3 )
    >>> results = h1.pdf_positive ( 3 , draw = 3 , silent = True )
    >>> print     results[ 0] ## fit results 
    >>> pdf     = results[ 1] ## get PDF
    >>> func    = results[ 2] ## normalized function  
    >>> norm    = results[ 3] ## normalization
    >>> pdf_fun = results[ 4] ## pdf-function (normalized)
    >>> frame   = results[-1] ## frame/RooPlot 
    """
    from ostap.fitting.background import PolyPos_pdf
    return _h1_pdf_ ( h1 , PolyPos_pdf , (degree,) , *args , **kwargs )


# =============================================================================
## parameterize/fit histogram with the positive polynomial
#  @code
#  h1 = ...
#  results = h1.pdf_even         ( 3 )
#  results = h1.pdf_positiveeven ( 3 )                               ## ditto 
#  results = h1.pdf_even         ( 3 , draw = True , silent = True )
#  results = h1.pdf_positiveeven ( 3 , draw = True , silent = True ) ## ditto 
#  print     results[ 0]
#  pdf     = results[ 1]
#  func    = results[ 2]
#  norm    = results[ 3]
#  pdf_fun = results[ 4] 
#  frame   = results[-1] 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-26
def _h1_pdf_even_ ( h1 , halfdegree , *args , **kwargs ) :
    """Parameterize/fit histogram with the positive even polynomial
    >>> h1 = ...
    >>> results = h1.pdf_even ( 2 )
    >>> results = h1.pdf_positiveeven ( 2 )                            ## ditto 
    >>> results = h1.pdf_even ( 2 , draw = 3 , silent = True )
    >>> results = h1.pdf_positiveeven ( 2 , draw = 3 , silent = True ) ## ditto
    >>> print     results[ 0] ## fit results 
    >>> pdf     = results[ 1] ## get PDF
    >>> func    = results[ 2] ## normalized function  
    >>> norm    = results[ 3] ## normalization
    >>> pdf_fun = results[ 4] ## pdf-function (normalized)
    >>> frame   = results[-1] ## frame/RooPlot 
    """
    from ostap.fitting.background import PolyEven_pdf
    return _h1_pdf_ ( h1 , PolyEven_pdf , (halfdegree,) , *args , **kwargs )


# =============================================================================
## parameterize/fit histogram with the monotonic positive polynomial
#  @code
#  h1 = ...
#  results = h1.pdf_monotonic ( 3 , increasing = True )
#  results = h1.pdf_monotonic ( 3 , increasing = True , draw = True , silent = True )
#  print     results[ 0]
#  pdf     = results[ 1]
#  func    = results[ 2]
#  norm    = results[ 3]
#  pdf_fun = results[ 4] 
#  frame   = results[-1] 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-26
def _h1_pdf_monotonic_ ( h1 , degree , increasing , *args , **kwargs ) :
    """Parameterize/fit histogram with the monotonic positive polynomial
    >>> h1 = ...
    >>> results = h1.pdf_monotonic ( 3 , increasing = True )
    >>> results = h1.pdf_monotonic ( 3 , increasing = True , draw = 3 , silent = True )
    >>> print     results[ 0] ## fit results 
    >>> pdf     = results[ 1] ## get PDF
    >>> func    = results[ 2] ## normalized function  
    >>> norm    = results[ 3] ## normalization
    >>> pdf_fun = results[ 4] ## pdf-function (normalized)
    >>> frame   = results[-1] ## frame/RooPlot 
    """
    from ostap.fitting.background import Monotonic_pdf
    return _h1_pdf_ ( h1 , Monotonic_pdf , (degree,increasing) , *args , **kwargs )

# =============================================================================
## parameterize/fit histogram with the increasing positive polynomial
#  @code
#  h1 = ...
#  results = h1.pdf_increasing ( 3 )
#  results = h1.pdf_increasing ( 3 , draw = 3 , silent = True )
#  print     results[ 0]
#  pdf     = results[ 1]
#  func    = results[ 2]
#  norm    = results[ 3]
#  pdf_fun = results[ 4] 
#  frame   = results[-1] 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-26
def _h1_pdf_increasing_ ( h1 , degree , *args , **kwargs ) :
    """Parameterize/fit histogram with the monotonic positive polynomial
    >>> h1 = ...
    >>> results = h1.pdf_increasing ( 3 )
    >>> results = h1.pdf_increasing ( 3 , draw = True , silent = True )
    >>> print     results[ 0] ## fit results 
    >>> pdf     = results[ 1] ## get PDF
    >>> func    = results[ 2] ## normalized function  
    >>> norm    = results[ 3] ## normalization
    >>> pdf_fun = results[ 4] ## pdf-function (normalized)
    >>> frame   = results[-1] ## frame/RooPlot 
    """
    return _h1_pdf_monotonic_ ( h1 , degree , True , *args , **kwargs ) 

# =============================================================================
## parameterize/fit histogram with the decreasing positive polynomial
#  @code
#  h1 = ...
#  results = h1.pdf_decreasing ( 3 )
#  results = h1.pdf_decreasing ( 3 , draw = True , silent = True )
#  print     results[ 0]
#  pdf     = results[ 1]
#  func    = results[ 2]
#  norm    = results[ 3]
#  pdf_fun = results[ 4] 
#  frame   = results[-1] 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-26
def _h1_pdf_decreasing_ ( h1 , degree , *args , **kwargs ) :
    """Parameterize/fit histogram with the monotonic positive polynomial
    >>> h1 = ...
    >>> results = h1.pdf_decreasing ( 3 )
    >>> results = h1.pdf_decreasing ( 3 , draw = True , silent = True )
    >>> print     results[ 0] ## fit results 
    >>> pdf     = results[ 1] ## get PDF
    >>> func    = results[ 2] ## normalized function  
    >>> norm    = results[ 3] ## normalization
    >>> pdf_fun = results[ 4] ## pdf-function (normalized)
    >>> frame   = results[-1] ## frame/RooPlot 
    """
    return _h1_pdf_monotonic_ ( h1 , degree , False , *args , **kwargs ) 

# =============================================================================
## parameterize/fit histogram with the convex polynomial
#  @code
#  h1 = ...
#  results = h1.pdf_convex ( 3 , increasing = True )
#  results = h1.pdf_convex ( 3 , increasing = True , draw = True , silent = True )
#  print     results[ 0]
#  pdf     = results[ 1]
#  func    = results[ 2]
#  norm    = results[ 3]
#  pdf_fun = results[ 4] 
#  frame   = results[-1] 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-26
def _h1_pdf_convex_ ( h1 , degree , increasing , *args , **kwargs ) :
    """Parameterize/fit histogram with convex polynomial
    >>> h1 = ...
    >>> results = h1.pdf_convex ( 3 , increasing = True )
    >>> results = h1.pdf_convex ( 3 , increasing = True , draw = True , silent = True )
    >>> print     results[ 0] ## fit results 
    >>> pdf     = results[ 1] ## get PDF
    >>> func    = results[ 2] ## normalized function  
    >>> norm    = results[ 3] ## normalization
    >>> pdf_fun = results[ 4] ## pdf-function (normalized)
    >>> frame   = results[-1] ## frame/RooPlot 
    """
    from ostap.fitting.background import Convex_pdf
    return _h1_pdf_ ( h1 , Convex_pdf , (degree,increasing,True) , *args , **kwargs )

# =============================================================================
## parameterize/fit histogram with the convex increasing polynomial
#  @code
#  h1 = ...
#  results = h1.pdf_convex_increasing ( 3 ,)
#  results = h1.pdf_convex_increasing ( 3 , draw = True , silent = True )
#  print     results[ 0]
#  pdf     = results[ 1]
#  func    = results[ 2]
#  norm    = results[ 3]
#  pdf_fun = results[ 4] 
#  frame   = results[-1] 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-26
def _h1_pdf_convex_increasing_ ( h1 , degree , *args , **kwargs ) :
    """Parameterize/fit histogram with convex increasing polynomial
    >>> h1 = ...
    >>> results = h1.pdf_convex_increasing ( 3 ,)
    >>> results = h1.pdf_convex_increasing ( 3 , draw = True , silent = True )
    >>> print     results[ 0] ## fit results 
    >>> pdf     = results[ 1] ## get PDF
    >>> func    = results[ 2] ## normalized function  
    >>> norm    = results[ 3] ## normalization
    >>> pdf_fun = results[ 4] ## pdf-function (normalized)
    >>> frame   = results[-1] ## frame/RooPlot 
    """
    return _h1_pdf_convex_ ( h1 , degree , True , *args , **kwargs )

# =============================================================================
## parameterize/fit histogram with the convex decreasing polynomial
#  @code
#  h1 = ...
#  results = h1.pdf_convex_decreasing ( 3 ,)
#  results = h1.pdf_convex_decreasing  ( 3 , draw = True , silent = True )
#  print     results[ 0]
#  pdf     = results[ 1]
#  func    = results[ 2]
#  norm    = results[ 3]
#  pdf_fun = results[ 4] 
#  frame   = results[-1] 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-26
def _h1_pdf_convex_decreasing_ ( h1 , degree , *args , **kwargs ) :
    """Parameterize/fit histogram with convex decreasing polynomial
    >>> h1 = ...
    >>> results = h1.pdf_convex_decreasing ( 3 ,)
    >>> results = h1.pdf_convex_decreasing ( 3 , draw = True , silent = True )
    >>> print     results[ 0] ## fit results 
    >>> pdf     = results[ 1] ## get PDF
    >>> func    = results[ 2] ## normalized function  
    >>> norm    = results[ 3] ## normalization
    >>> pdf_fun = results[ 4] ## pdf-function (normalized)
    >>> frame   = results[-1] ## frame/RooPlot 
    """
    return _h1_pdf_convex_ ( h1 , degree , False , *args , **kwargs )

# =============================================================================
## parameterize/fit histogram with the concave polynomial
#  @code
#  h1 = ...
#  results = h1.pdf_concave ( 3 , increasing = True )
#  results = h1.pdf_concave ( 3 , increasing = True , draw = True , silent = True )
#  print     results[ 0]
#  pdf     = results[ 1]
#  func    = results[ 2]
#  norm    = results[ 3]
#  pdf_fun = results[ 4] 
#  frame   = results[-1] 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-26
def _h1_pdf_concave_ ( h1 , degree , increasing , *args , **kwargs ) :
    """Parameterize/fit histogram with concave polynomial
    >>> h1 = ...
    >>> results = h1.pdf_concave ( 3 , increasing = True )
    >>> results = h1.pdf_concave ( 3 , increasing = True , draw = True , silent = True )
    >>> print     results[ 0] ## fit results 
    >>> pdf     = results[ 1] ## get PDF
    >>> func    = results[ 2] ## normalized function  
    >>> norm    = results[ 3] ## normalization
    >>> pdf_fun = results[ 4] ## pdf-function (normalized)
    >>> frame   = results[-1] ## frame/RooPlot 
    """
    from ostap.fitting.background import Convex_pdf
    return _h1_pdf_ ( h1 , Convex_pdf , (degree,increasing,False) , *args , **kwargs )

# =============================================================================
## parameterize/fit histogram with the concave increasing polynomial
#  @code
#  h1 = ...
#  results = h1.pdf_concave_increasing ( 3 )
#  results = h1.pdf_concave_increasing ( 3 , draw = True , silent = True )
#  print     results[ 0]
#  pdf     = results[ 1]
#  func    = results[ 2]
#  norm    = results[ 3]
#  pdf_fun = results[ 4] 
#  frame   = results[-1] 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-26
def _h1_pdf_concave_increasing_ ( h1 , degree , *args , **kwargs ) :
    """Parameterize/fit histogram with concave increasing polynomial
    >>> h1 = ...
    >>> results = h1.pdf_concave_increasing ( 3 )
    >>> results = h1.pdf_concave_increasing ( 3 , draw = True , silent = True )
    >>> print     results[ 0] ## fit results 
    >>> pdf     = results[ 1] ## get PDF
    >>> func    = results[ 2] ## normalized function  
    >>> norm    = results[ 3] ## normalization
    >>> pdf_fun = results[ 4] ## pdf-function (normalized)
    >>> frame   = results[-1] ## frame/RooPlot 
    """
    return _h1_pdf_concave_ ( h1 , degree , True , *args , **kwargs )

# =============================================================================
## parameterize/fit histogram with the concave decreasing polynomial
#  @code
#  h1 = ...
#  results = h1.pdf_concave_decreasing ( 3 )
#  results = h1.pdf_concave_decreasing ( 3 , draw = True , silent = True )
#  print     results[ 0]
#  pdf     = results[ 1]
#  func    = results[ 2]
#  norm    = results[ 3]
#  pdf_fun = results[ 4] 
#  frame   = results[-1] 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-26
def _h1_pdf_concave_decreasing_ ( h1 , degree , *args , **kwargs ) :
    """Parameterize/fit histogram with concave decreasing polynomial
    >>> h1 = ...
    >>> results = h1.pdf_concave_decreasing ( 3 )
    >>> results = h1.pdf_concave_decreasing ( 3 , draw = True , silent = True )
    >>> print     results[ 0] ## fit results 
    >>> pdf     = results[ 1] ## get PDF
    >>> func    = results[ 2] ## normalized function  
    >>> norm    = results[ 3] ## normalization
    >>> pdf_fun = results[ 4] ## pdf-function (normalized)
    >>> frame   = results[-1] ## frame/RooPlot 
    """
    return _h1_pdf_concave_ ( h1 , degree , False , *args , **kwargs )

# =============================================================================
## parameterize/fit histogram with the convex polynomial
#  @code
#  h1 = ...
#  results = h1.pdf_convexpoly ( 3 )
#  results = h1.pdf_convexpoly ( 3 , draw = True , silent = True )
#  print     results[ 0]
#  pdf     = results[ 1]
#  func    = results[ 2]
#  norm    = results[ 3]
#  pdf_fun = results[ 4] 
#  frame   = results[-1] 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-26
def _h1_pdf_convexpoly_ ( h1 , degree , *args , **kwargs ) :
    """Parameterize/fit histogram with convex polynomial
    >>> h1 = ...
    >>> results = h1.pdf_convexpoly ( 3 )
    >>> results = h1.pdf_convexpoly ( 3 , draw = True , silent = True )
    >>> print     results[ 0] ## fit results 
    >>> pdf     = results[ 1] ## get PDF
    >>> func    = results[ 2] ## normalized function  
    >>> norm    = results[ 3] ## normalization
    >>> pdf_fun = results[ 4] ## pdf-function (normalized)
    >>> frame   = results[-1] ## frame/RooPlot 
    """
    from ostap.fitting.background import ConvexOnly_pdf
    return _h1_pdf_ ( h1 , ConvexOnly_pdf , (degree,True) , *args , **kwargs )

# =============================================================================
## parameterize/fit histogram with the concave polynomial
#  @code
#  h1 = ...
#  results = h1.pdf_concavepoly ( 3 )
#  results = h1.pdf_concavepoly ( 3 , draw = True , silent = True )
#  print     results[ 0]
#  pdf     = results[ 1]
#  func    = results[ 2]
#  norm    = results[ 3]
#  pdf_fun = results[ 4] 
#  frame   = results[-1] 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-26
def _h1_pdf_concavepoly_ ( h1 , degree , *args , **kwargs ) :
    """Parameterize/fit histogram with convex polynomial
    >>> h1 = ...
    >>> results = h1.pdf_concavepoly ( 3 )
    >>> results = h1.pdf_concavepoly ( 3 , draw = True , silent = True )
    >>> print     results[ 0] ## fit results 
    >>> pdf     = results[ 1] ## get PDF
    >>> func    = results[ 2] ## normalized function  
    >>> norm    = results[ 3] ## normalization
    >>> pdf_fun = results[ 4] ## pdf-function (normalized)
    >>> frame   = results[-1] ## frame/RooPlot 
    """
    from ostap.fitting.background import ConvexOnly_pdf
    return _h1_pdf_ ( h1 , ConvexOnly_pdf , (degree,False) , *args , **kwargs )

# =============================================================================

for t in ( ROOT.TH1D , ROOT.TH1F ) :
    t.pdf_positive           = _h1_pdf_positive_
    t.pdf_positiveeven       = _h1_pdf_even_
    t.pdf_even               = _h1_pdf_even_
    t.pdf_monotonic         = _h1_pdf_monotonic_
    t.pdf_increasing         = _h1_pdf_increasing_
    t.pdf_decreasing         = _h1_pdf_decreasing_
    t.pdf_convex             = _h1_pdf_convex_
    t.pdf_convex_increasing  = _h1_pdf_convex_increasing_
    t.pdf_convex_decreasing  = _h1_pdf_convex_decreasing_
    t.pdf_concave            = _h1_pdf_concave_
    t.pdf_concave_increasing = _h1_pdf_concave_increasing_
    t.pdf_concave_decreasing = _h1_pdf_concave_decreasing_
    t.pdf_convexpoly         = _h1_pdf_convexpoly_
    t.pdf_concavepoly        = _h1_pdf_concavepoly_

    _new_methods_ += [
        _h1_pdf_positive_    ,
        _h1_pdf_even_        ,
        _h1_pdf_monotonic_  ,
        _h1_pdf_increasing_  ,
        _h1_pdf_decreasing_  ,
        _h1_pdf_convex_      ,
        _h1_pdf_convex_increasing_ ,
        _h1_pdf_convex_decreasing_ ,
        _h1_pdf_concave_      ,
        _h1_pdf_concave_increasing_ ,
        _h1_pdf_concave_decreasing_ ,
        _h1_pdf_convexpoly_   ,
        _h1_pdf_concavepoly_  ,        
        ]
# =============================================================================
## parameterize/fit histogram with the positive  b-spline 
#  @code
#  h1 = ...
#  results = h1.pdf_pSpline ( spline = ( 3 ,2 )  ) ## order=3, inner knots=2
#  results = h1.pdf_pSpline ( ( 3 , 2 ) ,  draw = True , silent = True )
#  print results[0]
#  pdf = results[2]
#  print results[3]
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-26
def _h1_pdf_pspline_ ( h1 , spline , *args , **kwargs ) :
    """Parameterize/fit histogram with positive b-spline 
    >>> h1 = ...
    >>> results = h1.pdf_pSpline ( spline = (3,2) )
    >>> results = h1.pdf_pSpline ( (3,2) , draw = True , silent = True )
    >>> print results[0] ## fit results 
    >>> pdf = results[1] ## get PDF 
    >>> print results[2] ## underlying parameterization 
    """
    #
    if isinstance ( spline , (tuple,list) ) : 
        ## create the spline with uniform binning 
        PS     = cpp.Ostap.Math.PositiveSpline
        spline = PS ( h1.xmin() , h1.xmax() , spline[1] , spline[0] )
    #
    from ostap.fitting.background import PSpline_pdf
    return _h1_pdf_ ( h1 , PSpline_pdf , ( spline , ) , *args , **kwargs )


# =============================================================================
## parameterize/fit histogram with the monotonic positive  b-spline 
#  @code
#  h1 = ...
#  results = h1.pdf_mSpline ( spline = ( 3 , 2 , True )  ) ## order=3, inner knots=2
#  results = h1.pdf_mSpline ( ( 3 , 2, False ) ,  draw = True , silent = True )
#  print results[0]
#  pdf = results[2]
#  print results[3]
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-26
def _h1_pdf_mspline_ ( h1 , spline , *args , **kwargs ) :
    """Parameterize/fit histogram with monotonic positive b-spline 
    >>> h1 = ...
    >>> results = h1.pdf_mSpline ( spline = (3,2,True) )
    >>> results = h1.pdf_mSpline ( (3,2,True) , draw = True , silent = True )
    >>> print results[0] ## fit results 
    >>> pdf = results[1] ## get PDF 
    >>> print results[2] ## underlying parameterization 
    """
    #
    if isinstance ( spline , ( tuple , list ) ) :
        ## create the spline with uniform binning 
        MS     = cpp.Ostap.Math.MonotonicSpline
        spline = MS ( h1.xmin() , h1.xmax() , spline[1] , spline[0] , spline[2] )
        #
    from ostap.fitting.background import MSpline_pdf
    return _h1_pdf_ ( h1 , MSpline_pdf , ( spline , ) , *args , **kwargs )

# =============================================================================
## parameterize/fit histogram with the convex/concave monotonic positive  b-spline 
#  @code
#  h1 = ...
#  results = h1.pdf_cSpline ( spline = ( 3 , 2 , True )  ) ## order=3, inner knots=2
#  results = h1.pdf_cSpline ( ( 3 , 2, False ) ,  draw = True , silent = True )
#  print results[0]
#  pdf = results[2]
#  print results[3]
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-26
def _h1_pdf_cspline_ ( h1 , spline , *args , **kwargs ) :
    """Parameterize/fit histogram with convex/concave montonic positive b-spline 
    >>> h1 = ...
    >>> results = h1.pdf_cSpline ( spline = (3,2,True,True) )
    >>> results = h1.pdf_cSpline ( (3,2,True,True), draw = True , silent = True )
    >>> print results[0] ## fit results 
    >>> pdf = results[1] ## get PDF 
    >>> print results[2] ## underlying parameterization 
    """
    #
    if isinstance ( spline , (tuple,list) ) : 
        ## create the spline with uniform binning 
        CS     = cpp.Ostap.Math.ConvexSpline
        spline = CS ( h1.xmin() , h1.xmax() , spline[1] , spline[0] , spline[2] , spline[3] )
    #
    from ostap.fitting.background import CSpline_pdf
    return _h1_pdf_ ( h1 , CSpline_pdf , ( spline , ) , *args , **kwargs )


# =============================================================================
## parameterize/fit histogram with the positive  convex b-spline 
#  @code
#  h1 = ...
#  results = h1.pdf_convexSpline ( spline = ( 3 , 2 )  ) ## order=3, inner knots=2
#  results = h1.pdf_convexSpline ( ( 3 , 2 ) ,  draw = True , silent = True )
#  print results[0]
#  pdf = results[2]
#  print results[3]
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-26
def _h1_pdf_convexSpline_ ( h1 , spline , *args , **kwargs ) :
    """Parameterize/fit histogram with convex positive b-spline 
    >>> h1 = ...
    >>> results = h1.pdf_convexSpline ( spline = (3,2) )
    >>> results = h1.pdf_convexSpline ( ( 3 , 2 ), draw = True , silent = True )
    >>> print results[0] ## fit results 
    >>> pdf = results[1] ## get PDF 
    >>> print results[2] ## underlying parameterization 
    """
    #
    if isinstance ( spline , ( tuple , list ) ) : 
        ## create the spline with uniform binning 
        CS     = cpp.Ostap.Math.ConvexOnlySpline
        spline = CS ( h1.xmin() , h1.xmax() , spline[1] , spline[0] , True )
        #
    from ostap.fitting.background import CPSpline_pdf
    return _h1_pdf_ ( h1 , CPSpline_pdf , ( spline , ) , *args , **kwargs )


# =============================================================================
## parameterize/fit histogram with the positive  concave b-spline 
#  @code
#  h1 = ...
#  results = h1.pdf_concaveSpline ( spline = ( 3 , 2 )  ) ## order=3, inner knots=2
#  results = h1.pdf_concaveSpline ( ( 3 , 2 ) ,  draw = True , silent = True )
#  print results[0]
#  pdf = results[2]
#  print results[3]
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-26
def _h1_pdf_concaveSpline_ ( h1 , spline , *args , **kwargs ) :
    """Parameterize/fit histogram with positive concave b-spline 
    >>> h1 = ...
    >>> results = h1.pdf_concaveSpline ( spline =  ( 3 , 2 ) )
    >>> results = h1.pdf_concaveSpline ( ( 3 , 2 ), draw = True , silent = True )
    >>> print results[0] ## fit results 
    >>> pdf = results[1] ## get PDF 
    >>> print results[2] ## underlying parameterization 
    """
    #
    if isinstance ( spline , ( tuple , list ) ) : 
        ## create the spline with uniform binning 
        CS     = cpp.Ostap.Math.ConvexOnlySpline
        spline = CS ( h1.xmin() , h1.xmax() , spline[1] , spline[0] , False )
        #
    from ostap.fitting.background import CPSpline_pdf
    return _h1_pdf_ ( h1 , CPSpline_pdf , ( spline , ) , *args , **kwargs )



## decorate !
for t in ( ROOT.TH1D , ROOT.TH1F ) :
    t.pdf_pSpline       = _h1_pdf_pspline_
    t.pdf_mSpline       = _h1_pdf_mspline_
    t.pdf_cSpline       = _h1_pdf_cspline_
    t.pdf_convexSpline  = _h1_pdf_convexSpline_
    t.pdf_concaveSpline = _h1_pdf_concaveSpline_

    _new_methods_ += [
        _h1_pdf_pspline_        ,
        _h1_pdf_mspline_        ,
        _h1_pdf_cspline_        ,
        _h1_pdf_convexSpline_   ,
        _h1_pdf_concaveSpline_  ,
        ]
# =============================================================================
## make a histogram representation in terms of Legendre polynomials
#  @code 
#  histo  = ...
#  fsum   = histo.legendre_sum ( 4 )
#  print fsum
#  x = ...
#  print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
#  @endcode 
#  @see Ostap::Math::LegendreSum
#  @author Vanya Belyaev Ivan.Belyaev@iter.ru
#  It is not very CPU efficient (scipy is used for integration), but stable enough...
#  @date 2015-07-26
def _h1_legendre_sum_ ( h1 , N , **kwargs ) :
    """Make a histogram representation in terms of Legendre polynomials
    >>> histo  = ...
    >>> fsum   = histo.legendre_sum ( 4 )
    >>> print fsum
    >>> x = ...
    >>> print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
    It is not very CPU efficient (scipy is used for integration), but stable enough...
    """
    ##
    xmin = max ( kwargs.get( 'xmin' , h1.xmin() ) , h1.xmin () ) 
    xmax = min ( kwargs.get( 'xmax' , h1.xmax() ) , h1.xmax () ) 
    ##
    return legendre_sum ( h1 , N , xmin , xmax , *kwargs )

# =============================================================================
## make a histogram representation in terms of Chebyshev polynomials
#  @code 
#  histo  = ...
#  fsum   = histo.chebyshev_sum ( 4 )
#  print fsum
#  x = ...
#  print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
#  @endcode 
#  @see Ostap::Math::ChebyshevSum
#  @author Vanya Belyaev Ivan.Belyaev@iter.ru
#  It is not very CPU efficient (scipy is used for integration), but stable enough...
#  @date 2015-07-26
def _h1_chebyshev_sum_ ( h1 , N , **kwargs ) :
    """Make a histogram representation in terms of Chebyshev polynomials
    >>> histo  = ...
    >>> fsum   = histo.chebyshev_sum ( 4 )
    >>> print fsum
    >>> x = ...
    >>> print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
    """
    ##
    xmin = max ( kwargs.get( 'xmin' , h1.xmin() ) , h1.xmin () ) 
    xmax = min ( kwargs.get( 'xmax' , h1.xmax() ) , h1.xmax () ) 
    ##
    return chebyshev_sum ( h1 , N , xmin , xmax )

# =============================================================================
## make a histogram representation in terms of Fourier serie
#  @code 
#  histo  = ...
#  fsum   = histo.fourier_sum ( 4 )
#  print fsum
#  x = ...
#  print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
#  @endcode 
#  @see Ostap::Math::FourierSum
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2015-07-26
def _h1_fourier_sum_ ( h1 , N , fejer = False , **kwargs ) :
    """Make a histogram representation in terms of Fourier serie
    >>> histo  = ...
    >>> fsum   = histo.fourier_sum ( 4 )
    >>> print fsum
    >>> x = ...
    >>> print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
    """
    ##
    xmin = max ( kwargs.get( 'xmin' , h1.xmin() ) , h1.xmin () ) 
    xmax = min ( kwargs.get( 'xmax' , h1.xmax() ) , h1.xmax () ) 
    ##
    return fourier_sum ( h1 , N , xmin , xmax , fejer )

# =============================================================================
## make a histogram representation in terms of cosine Fourier serie
#  @code 
#  histo  = ...
#  fsum   = histo.cosine_sum ( 4 )
#  print fsum
#  x = ...
#  print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
#  @endcode 
#  @see Ostap::Math::CosineSum
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2015-07-26
def _h1_cosine_sum_ ( h1 , N , fejer = False , **kwargs ) :
    """Make a histogram representation in terms of cosine Fourier serie
    >>> histo  = ...
    >>> fsum   = histo.cosine_sum ( 4 )
    >>> print fsum
    >>> x = ...
    >>> print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
    """
    ##
    xmin = max ( kwargs.get( 'xmin' , h1.xmin() ) , h1.xmin () ) 
    xmax = min ( kwargs.get( 'xmax' , h1.xmax() ) , h1.xmax () ) 
    ##
    return cosine_sum ( h1 , N , xmin , xmax , fejer )

# =============================================================================
## make a histogram representation in terms of Bezier(Bernstein) sum
#  (sum over Bernstein polynomials)
#  @code 
#  histo  = ...
#  fsum   = histo.bezier_sum    ( 4 )
#  fsum   = histo.bernstein_sum ( 4 ) ## distto
#  print fsum
#  x = ...
#  print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
#  @endcode 
#  @see Ostap::Math::Bernstein
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2015-07-26
def _h1_bezier_sum_ ( h1 , N , **kwargs ) :
    """Make a histogram representation in terms of Bezier/Bernstein sum
    (sum over Bernstein polynomials)
    >>> histo  = ...
    >>> fsum   = histo.bezier_sum    ( 4 )
    >>> fsum   = histo.bernstein_sum ( 4 ) ## ditto 
    >>> print fsum
    >>> x = ...
    >>> print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
    """
    ##
    xmin = max ( kwargs.get( 'xmin' , h1.xmin() ) , h1.xmin () ) 
    xmax = min ( kwargs.get( 'xmax' , h1.xmax() ) , h1.xmax () ) 
    ##
    return bezier_sum ( h1 , N , xmin , xmax )


# =============================================================================
## make a histogram representation in terms of Bezier(Bernstein) sum
#  (sum over Bernstein polynomials) using even polynomials:
#  \f$ f( \frac{x_{min}+x_{max}}{2} - x ) = f( \frac{x_{min}+x_{max}}{2} + x ) \f$ 
#  @code 
#  histo  = ...
#  fsum   = histo.beziereven_sum    ( 2 )
#  fsum   = histo.bernsteineven_sum ( 2 ) ## distto
#  print fsum
#  x = ...
#  print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
#  @endcode 
#  @see Ostap::Math::BernsteinEven
#  @param h1 the historgram
#  @param N  the half-degree actual degree of polynomial is 2*N
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2015-07-26
def _h1_beziereven_sum_ ( h1 , N , **kwargs ) :
    """Make a histogram representation in terms of Bezier/Bernstein sum
    (sum over Bernstein polynomials) using even polynomials 
    >>> histo  = ...
    >>> fsum   = histo.beziereven_sum    ( 2 ) ##  2 is a *half-degree*
    >>> fsum   = histo.bernsteineven_sum ( 2 ) ## ditto 
    >>> print fsum
    >>> x = ...
    >>> print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
    """
    ##
    xmin = max ( kwargs.get( 'xmin' , h1.xmin() ) , h1.xmin () ) 
    xmax = min ( kwargs.get( 'xmax' , h1.xmax() ) , h1.xmax () ) 
    ##
    return beziereven_sum ( h1 , N , xmin , xmax )


for h in ( ROOT.TH1F , ROOT.TH1D ) :

    h.legendre_sum      = _h1_legendre_sum_
    h.chebyshev_sum     = _h1_chebyshev_sum_
    h.fourier_sum       = _h1_fourier_sum_
    h.cosine_sum        = _h1_cosine_sum_
    h.bezier_sum        = _h1_bezier_sum_
    h.bernstein_sum     = _h1_bezier_sum_
    h.beziereven_sum    = _h1_beziereven_sum_
    h.bernsteineven_sum = _h1_beziereven_sum_

    _new_methods_ .append ( legendre_sum   )
    _new_methods_ .append ( chebyshev_sum  )
    _new_methods_ .append ( fourier_sum    )
    _new_methods_ .append ( cosine_sum     )
    _new_methods_ .append ( bezier_sum     )
    _new_methods_ .append ( beziereven_sum )

# =============================================================================
_decorated_classes = (
    ROOT.TH1D ,
    ROOT.TH1F ,
    )

_new_methods_ = tuple ( _new_methods_ )
    

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
