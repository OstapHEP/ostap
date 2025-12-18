#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with utilities for parameterization of histograms 
#
## (1) using fit of histograms (ROOT::TF1::Fit)
# 
# - as      Bernstein/Bezier sum 
# - as even Bernstein/Bezier sum 
# - as      Chebyshev        sum  
# - as      Legendre         sum 
# - as      Fourier          sum (if numpy/scipy available)
# - as      Cosine/Fourier   sum (if numpy/scipy available) 
# - as      plain monomal    sum
# - as positive                          Bernstein/Bezier sum
# - as positive even                     Bernstein/Bezier sum
# - as positive monotonic                Bernstein/Bezier sum
# - as positive           convex/concave Bernstein/Bezier sum
# - as positive monotonic convex/concave Bernstein/Bezier sum
# - as      generic                           spline  (b-spline)
# - as      positive                          spline  (p-spline)
# - as      positive monotonic                spline  (m-spline) 
# - as      positive monotonic convex/concave spline  (c-spline) 
# - as      positive           convex/concave spline  (convex/concave-spline) 
#
# where possible, fit starts from reasonable approximatuion, taken from (1)
# 
# Typical usage:
# @code
# histo = ...                ## the histogram
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
# - as positive                          Bernstein/Bezier sum
# - as positive even                     Bernstein/Bezier sum
# - as positive monotonic                Bernstein/Bezier sum
# - as positive           convex/concave Bernstein/Bezier sum
# - as positive monotonic convex/concave Bernstein/Bezier sum
# - as      positive                          spline  (p-spline)
# - as      positive monotonic                spline  (m-spline) 
# - as      positive monotonic convex/concave spline  (c-spline) 
# - as      positive           convex/concave spline  (convex/concave-spline)
# - as      positive rational function
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
""" Module with utilities for parameterization of histograms

## (1) using fit of histograms (ROOT::TF1::Fit)

- as      Bernstein/Bezier sum 
- as even Bernstein/Bezier sum 
- as      Chebyshev        sum  
- as      Legendre         sum 
- as      Fourier          sum (if numpy/scipy available)
- as      Cosine/Fourier   sum (is numpy/scipy available)
- as      plain monomal    sum
- as positive                           Bernstein/Bezier sum
- as positive even                      Bernstein/Bezier sum
- as positive monotonic                 Bernstein/Bezier sum
- as positive            convex/concave Bernstein/Bezier sum
- as positive monotonic  convex/concave Bernstein/Bezier sum
- as      generic                           spline  (b-spline)
- as      positive                          spline  (p-spline)
- as      positive monotonic                spline  (m-spline) 
- as      positive monotonic convex/concave spline  (c-spline) 
- as      positive           convex/concave spline  (convex/concave-spline) 

where possible, fit starts from reasonable approximatuion, taken from (1)

Typical usage:

>>> histo = ...                ## the histogram
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
- as positive monotonic                 Bernstein/Bezier sum
- as positive            convex/concave Bernstein/Bezier sum
- as positive monotonic convex/concave  Bernstein/Bezier sum
- as      positive                          spline  (p-spline)
- as      positive monotonic                spline  (m-spline) 
- as      positive monotonic convex/concave spline  (c-spline) 
- as      positive           convex/concave spline  (convex/concave-spline)

Typical usage:

>>> histo = ...                   ## the histogram
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
from   ostap.core.ostap_types import integer_types, num_types
from   ostap.core.core        import cpp, VE, funID, Ostap
from   ostap.math.param       import ( legendre_sum      ,
                                       chebyshev_sum     ,
                                       bezier_sum        ,
                                       bernstein_sum     , 
                                       beziereven_sum    ,
                                       bernsteineven_sum ,
                                       rational_fun      ) 
from   ostap.utils.ranges     import vrange 
from   ostap.utils.basic      import typename
from   collections            import namedtuple
import ROOT, math 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger    import getLogger 
# =============================================================================
if '__main__' ==  __name__ : logger = getLogger( 'ostap.histos.param' )
else                       : logger = getLogger( __name__             )
# =============================================================================
logger.debug ( 'Some parameterization utilities for Histo objects')
# =============================================================================
_new_methods_ = []
# =============================================================================
inf_pos =  float('inf') ## positive infinity
inf_neg = -float('inf') ## negative infinity
# =============================================================================
## result of the histogram  parameterisation based on TH1::Fit 
ParamFITInfo = namedtuple ( 'ParamFITInfo' , ( 'tf1'       ,   ## ROOT TF1 object 
                                               'fitobject' ,   ## Fittting object 
                                               'funobject' ,   ## original function object 
                                               'fitresult' ,   ## result of TH1::Fit
                                               'norm'      ,   ## normalization factor
                                               'nfits'     ) ) ## number of refits 
# =============================================================================
## result of the histogram  parameterisation based on RooFit
ParamPDFInfo = namedtuple ( 'ParamPDFInfo' , ( 'fitresult' ,   ## RooFitResult
                                               'pdf'       ,   ## Fit PDF 
                                               'funobject' ,   ## the function object 
                                               'norm'      ,   ## normalization factor 
                                               'pdffun'    ,   ## PDF_fun object
                                               'plot'      ) ) ## RooPlot object 
# =============================================================================
## represent 1D-histo as polynomial-like sum 
def _h1_param_sum_ ( h1               ,
                     fun_obj          ,
                     fit_type         ,  
                     opts   = 'SQ0'   ,
                     xmin   = inf_neg ,
                     xmax   = inf_pos ,
                     fixes  = ()      ,   ## List [ (i1,value1) , .... , (i_n,value_n) ] 
                     params = ()      ,   ## List [ value1 ,value2 , ... , value_n     ]
                     limits = ()      ,   ## Triplets [ ( i , min , max ) , ... ]
                     refit  = False   ) : ## refit ? 
    """ Represent histo as polynomial-like  sum    
    """
    
    import ostap.fitting.funcs 

    xmin  = max ( xmin , h1.xmin() )
    xmax  = min ( xmax , h1.xmax() )

    assert h1.xmin() <= xmin and xmin <  h1.xmax() , "Invalid xmin!"
    assert h1.xmin() <  xmax and xmax <= h1.xmax() , "Invalid xmax!"
    assert xmin < xmax                             , "Invalid xmin/xmax!"
    
    ## 
    b     = fun_obj  
    #
    bfit  = fit_type ( b , xmin = xmin , xmax = xmax )
    
    bfit.fun.SetNpx  ( max ( 100 , 3 * h1.bins() , bfit.fun.GetNpx() ) )   
    
    bfit.histo     = h1
    
    fun = bfit.fun

    normalized = hasattr ( bfit , 'norm' ) and  bfit.norm()
    
    if normalized :
        
        ## calculate the integral in range 
        ## _integral_ = h1.integrate ( xmin = xmin , xmax = xmax ) 
        _integral_ = h1.integral ( xmin = xmin , xmax = xmax )
        fun.SetParameter ( 0 , _integral_ )
        
        for i in range ( 0, b.npars() ) :
            fun.SetParameter ( i + 1 , 0  )
            fun.SetParLimits ( i + 1 , -10 * math.pi  , 10 * math.pi )
            
    else :
        
        for i in range ( 0, b.npars() ) :
            fun.SetParameter ( i  , 0  )
                
    if not opts                : opts  = 'S'
    if not 'S' in opts.upper() : opts += 'S'
    
    if ( len ( h1 ) < 100 or h1.GetXaxis().IsVariableBinSize() ) and not 'I' in opts.upper() :
        logger.info ("param_sum: add fitting option 'I'") 
        opts += 'I'
        
    ## fitting options:
    fopts = opts,'',xmin,xmax 

    np = bfit.fun.GetNpar ()
    ## set limits for parameter (if specified)
    for i, l, h in limits :
        if l < np and l <= h :
            logger.verbose ( 'param_sum: set limits    %d as %s/%s' % ( i , l , h  ) ) 
            fun.SetParLimits ( i , l , h )
            
    ## set parameters (if specified) 
    for i , v in enumerate ( params ) :
        if i < np : 
            logger.verbose ( 'param_sum: set parameter %d at %s' % ( i , v ) ) 
            fun.SetParameter ( i , float ( v ) )
        
    ## fix parameters (if specified) 
    for i , v in fixes :
        if i < np : 
            logger.verbose ( 'param_sum: fix parameter %d at %s' % ( i , v ) ) 
            fun.FixParameter ( i , float ( v ) )
            
    if normalized :
        fun.FixParameter    ( 0 , _integral_ )
        r  = fun.Fit ( h1 , opts+'0Q', '', xmin , xmax )        
        fun.ReleaseParameter ( 0 )
        
    import ostap.fitting.fitresult
    from   ostap.fitting.utils     import fit_status
    from   ostap.logger.colorized  import attention

    ## the fit itself 
    r = fun.Fit( h1 , *fopts )
    
    if isinstance ( refit , integer_types ) : refit = max ( 0 , refit )
    else                                    : refit = 1 if refit else 0 

    nfits = 1 

    while r.Status() and 0 < refit :
        
        status = attention ( fit_status ( r.Status() ) ) 
        b_name = typename ( b ) 
        rtable = r.table ( title = 'Fit result for %s' % b_name , prefix = '# ' )
        logger.debug   ( 'Fit status is %s [%s]\n%s' % ( status , b_name ,  rtable ) )

        if normalized : 
            fun.FixParameter ( 0 , _integral_ )
            r  = fun.Fit( h1 , opts+'0Q', '', xmin , xmax )
            fun.ReleaseParameter ( 0 )
            
        r = fun.Fit ( h1 , *fopts )    
        refit -= 1 
        nfits += 1
        
    if r.Status() :
        status = attention ( fit_status ( r.Status() ) ) 
        b_name = typename ( b ) 
        rtable = r.table ( title = 'Fit result for %s' % b_name , prefix = '# ' )
        logger.error ( 'Fit result [%s]\n%s'   % ( b_name , rtable ) )

    bfit.fitresult = r
    
    scale = VE ( 1 , 0 )
    if normalized : 
        bfit.fitnorm = r [ 0 ]
        scale        = r [ 0 ]
        
    results            = ParamFITInfo ( bfit.fun , bfit , b , bfit.fitresult , scale , nfits  ) 
    h1._param_FIT_info = results 
    return results

# =============================================================================
## represent 1D-histo as Bernstein polynomial
#  @code
#  h = ...                   ## the histogram
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
def _h1_bernstein_ ( h1               ,
                     degree           ,
                     opts   = 'SQ0'   ,
                     xmin   = inf_neg ,
                     xmax   = inf_pos ,
                     fixes  = ()      ,
                     params = ()      ,
                     limits = ()      ,
                     refit  = 1       ) :
    """ Represent histo as Bernstein polynomial
    
    >>> h = ...                # the histogram
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
    ## make a fit 
    if not params : params = tuple ( [ p for p in func.pars() ] ) 
    from ostap.fitting.param import H_fit
    return _h1_param_sum_ ( h1              ,
                            func            ,
                            H_fit           ,
                            opts   = opts   ,
                            xmin   = xmin   ,
                            xmax   = xmax   ,
                            fixes  = fixes  ,
                            params = params ,
                            limits = limits , 
                            refit  = refit  )

# =============================================================================
## represent 1D-histo as even Bernstein polynomial
#  @code
#  h = ...                   ## the histogram
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
def _h1_bernsteineven_ ( h1               ,
                         degree           ,
                         opts   = 'SQ0'   ,
                         xmin   = inf_neg ,
                         xmax   = inf_pos ,
                         fixes  = ()      ,
                         params = ()      ,
                         limits = ()      ,
                         refit  = 1       ) :
    """ Represent histo as even Bernstein polynomial
    
    >>> h = ...                    ## the histogram
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
    func  = beziereven_sum ( h1   , degree , xmin , xmax )
    ## 
    ## make a fit 
    if not params : params = tuple ( [ p for p in func.pars() ] ) 
    from ostap.fitting.param import H_fit
    return _h1_param_sum_ ( h1              ,
                            func            ,
                            H_fit           ,
                            opts   = opts   ,
                            xmin   = xmin   ,
                            xmax   = xmax   ,
                            fixes  = fixes  ,
                            params = params ,
                            limits = limits ,
                            refit  = refit  )


# =============================================================================
## represent 1D-histo as Chebyshev polynomial
#  @code
#  h = ...                  ## the histogram
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
def _h1_chebyshev_ ( h1               ,
                     degree           ,
                     opts   = 'SQ0'   ,
                     xmin   = inf_neg ,
                     xmax   = inf_pos ,
                     fixes  = ()      ,
                     params = ()      ,
                     limits = ()      ,
                     refit  = 1       ) :
    """ Represent histo as Chebyshev sum 
    
    >>> h = ... # the histogram
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
    ## make a fit 
    if not params : params = tuple ( [ p for p in func.pars() ] ) 
    from ostap.fitting.param import H_fit
    return _h1_param_sum_ ( h1              ,
                            func            ,
                            H_fit           ,
                            opts   = opts   ,
                            xmin   = xmin   ,
                            xmax   = xmax   ,
                            fixes  = fixes  ,
                            params = params ,
                            limits = limits ,
                            refit  = refit  )

# =============================================================================
## represent 1D-histo as Legendre polynomial
#  @code
#  h = ...                  ## the histogram
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
def _h1_legendre_ ( h1               ,
                    degree           ,
                    opts   = 'SQ0'   ,
                    xmin   = inf_neg ,
                    xmax   = inf_pos ,
                    fixes  = ()      ,
                    params = ()      ,
                    limits = ()      ,
                    refit  = 1       ) :
    """ Represent histo as Legendre sum 
    
    >>> h = ... # the histogram
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
    ## make a fit 
    if not params : params = tuple ( [ p for p in func.pars() ] ) 
    from ostap.fitting.param import H_fit
    return _h1_param_sum_ ( h1              ,
                            func            ,
                            H_fit           ,
                            opts   = opts   ,
                            xmin   = xmin   ,
                            xmax   = xmax   ,
                            fixes  = fixes  ,
                            params = params ,
                            limits = limits ,
                            refit  = refit  )

# ==============================================================================
## (relatively) fast parameterization of 1D histogram as sum of
#  legendre polynomials
#  @code
#  histo = ...
#  func  = histo.legendre_fast ( 5 ) 
#  @endcode
#  @see Ostap::LegendreSum
#  @see Ostap::LegendreSum::fill 
def _h1_legendre_fast_ ( h1 , degree , xmin = inf_neg , xmax = inf_pos ) :
    """ (relatively) fast parameterization of 1D histogram as sum of Legendre polynomials
    >>> histo = ...
    >>> func  = histo.legendre_fast ( 5 )
    - see Ostap.Math.LegendreSum
    - see Ostap.Math.LegendreSum.fill
    """
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )
    assert isinstance ( degree , int ) and 0 <= degree , \
           "Invalid degree %s" % degree
    assert xmin < xmax                                 , \
           'Invalid xmin/xmax: %s/%s' % ( xmin , xmax )
    
    func = Ostap.Math.LegendreSum ( degree ,  xmin , xmax )
    
    wsum  = 0
    
    for i, x, y in h1.iteritems()  :
        xv = x.value ()
        if  not xmin <= xv <= xmax : continue 
        yv = y.value ()
        w  = 2*yv*x.error() 
        ok = func.fill ( xv , w )
        if ok : wsum +=       w 
    
    return func 
        
# =============================================================================
## represent 1D-histo as rational function. based on Floater-Hormann's
#  rational barycentric interpolant
# 
#  @code
#  h = ...                   ## the histogram
#  b = h.rational ( 5, 2  )  ## make a fit... 
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
def _h1_rational_  ( h1               ,
                     p                , ## degree of numerator 
                     d                ,
                     opts   = 'SQ0'   ,
                     xmin   = inf_neg ,
                     xmax   = inf_pos ,
                     fixes  = ()      ,
                     params = ()      ,
                     limits = ()      ,
                     refit  = 1       ) :
    """ Represent histo as ratioal function, based on Floater-Hormann's
    rational barycentric interpolant 
    
    >>> h = ...                    # the histogram
    >>> b = h.rational ( 5 , 2  )  ## make a fit... 
    
    >>> tf1        = b[0]    ## TF1 object
    >>> obj        = b[1]    ## helper object 
    >>> fun        = b[2]    ## underlying normalzed C++ object 
    >>> fit_result = b[3]    ## fit result & status
    >>> norm       = b[4]    ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) 
    """
    ## 
    assert isinstance ( p , integer_types ) and 0 <= p , "Invalid p: %s" % p
    assert isinstance ( d , integer_types ) and 0 <= d , "Invalid d: %s" % d
    ##
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )
    ##
    # make reasonable approximation
    func  = rational_fun ( h1 , p + d , d , xmin , xmax )
    ## make a fit 
    if not params : params = tuple ( [ v for v in func.pars() ] )
    ## 
    from ostap.fitting.param import H_fit
    return _h1_param_sum_ ( h1              ,
                            func            ,
                            H_fit           ,
                            opts   = opts   ,
                            xmin   = xmin   ,
                            xmax   = xmax   ,
                            fixes  = fixes  ,
                            params = params ,
                            limits = limits , 
                            refit  = refit  )


# =============================================================================
## represent 1D-histo as rational function - ration of Bernstein and positive Bernstein polynomials 
# 
#  @code
#  h = ...                    ## the histogram
#  b = h.brational ( 3 , 3 )  ## make a fit... 
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
def _h1_brational_  ( h1               ,
                      n                , ## degree of numerator 
                      d                , ## degree of denominator 
                      opts   = 'SQ0'   ,
                      xmin   = inf_neg ,
                      xmax   = inf_pos ,
                      fixes  = ()      ,
                      params = ()      ,
                      limits = ()      ,
                      refit  = 1       ) :
    """ Represent 1D-histo as rational function - ratio
    of Bernstein and positive Bernstein polynomials 
    
    >>> h = ...                    # the histogram
    >>> b = h.brational ( 3 , 3  ) ## make a fit... 
    
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

    ii = 1.0 / ( xmax - xmin ) 

    # make reasonable approximation
    func  = Ostap.Math.RationalBernstein ( n, d , xmin , xmax )

    first = not params and not fixes and not limits 
    
    if not params and not fixes : 
        np    = func.pnpars()
        for i , x in enumerate ( vrange ( xmin , xmax , n ) ) :
            func.setPar ( i , float ( h1 ( x ) ) * ii / max ( 1 , np ) )

    from ostap.fitting.param import H_fit

    if not fixes :
        
        nopts = opts.lower()
        if not 'q' in nopts : nopts += 'q'
        if not '0' in nopts : nopts += '0'
        
        ## fix denominator
        p0 = tuple ( [ func.par  ( i )  for p in range ( func.npars() ) ] ) 
        f0 = [ ( i , func.par ( i ) )   for i in range ( np , func.npars() ) ] 
        r0 = _h1_param_sum_ ( h1 ,
                              func            ,
                              H_fit           ,
                              opts   = nopts  ,
                              xmin   = xmin   ,
                              xmax   = xmax   ,
                              fixes  = f0     ,
                              params = p0     )
        
    ## make a fit 
    if not params : params = tuple ( [ p for p in func.pars() ] ) 
    return _h1_param_sum_ ( h1              ,
                            func            ,
                            H_fit           ,
                            opts   = opts   ,
                            xmin   = xmin   ,
                            xmax   = xmax   ,
                            fixes  = fixes  ,
                            params = params ,
                            limits = limits , 
                            refit  = refit  )

# =================================================================================
try : # ===========================================================================
    # =============================================================================
    from ostap.math.param       import fourier_sum    
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
    def _h1_fourier_sum_ ( h1 , N , **kwargs ) :
        """ Make a histogram representation in terms of Fourier serie
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
        return fourier_sum ( h1 , N , xmin , xmax )

    # =============================================================================
    ## represent 1D-histo as Fourier polynomial
    #  @code
    #  h = ...                  ## the histogram
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
    def _h1_fourier_ ( h1 ,  
                       degree           ,
                       opts   = 'SQ0'   ,
                       xmin   = inf_neg ,
                       xmax   = inf_pos ,
                       fixes  = ()      ,
                       params = ()      ,
                       limits = ()      , 
                       refit  = 1       ) :
        
        """ Represent histo as Fourier sum 
        
        >>> h = ... # the histogram
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
        func = fourier_sum ( h1 , degree , xmin , xmax )
    
        ## make a fit 
        if not params : params = tuple ( [ p for p in func.pars() ] ) 
        from ostap.fitting.param import H_fit
        return _h1_param_sum_ ( h1              ,
                                func            ,
                                H_fit           ,
                                opts   = opts   ,
                                xmin   = xmin   ,
                                xmax   = xmax   ,
                                fixes  = fixes  ,
                                params = params ,
                                limits = limits ,
                                refit  = refit  )

    
    for t in ( ROOT.TH1F , ROOT.TH1D ) :
        
        t.fourier_sum    = _h1_fourier_sum_
        t.fourier        = _h1_fourier_

    _new_methods_ .append ( _h1_fourier_sum_  )
    _new_methods_ .append ( _h1_fourier_      )

    # ==============================================================================
except ImportError : # =============================================================
    # ==============================================================================
    pass

# =================================================================================
try : # ===========================================================================
    # =============================================================================
    from ostap.math.param       import cosine_sum    
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
    def _h1_cosine_sum_ ( h1 , N , **kwargs ) :
        """ Make a histogram representation in terms of cosine Fourier serie
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
        return cosine_sum ( h1 , N , xmin , xmax )
        
    # =============================================================================
    ## represent 1D-histo as cosine Fourier polynomial
    #  @code
    #  h = ...                  ## the histogram
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
    def _h1_cosine_ ( h1 ,
                      degree           ,
                      opts   = 'SQ0'   ,
                      xmin   = inf_neg ,
                      xmax   = inf_pos ,
                      fixes  = ()      ,
                      params = ()      ,
                      limits = ()      ,
                      refit  = 1       ) :        
        """ Represent histo as Cosine Fourier sum         
        >>> h = ... # the histogram
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
        func = cosine_sum ( h1 , degree , xmin , xmax )
        
        ## make a fit 
        if not params : params = tuple ( [ p for p in func.pars() ] ) 
        from ostap.fitting.param import H_fit
        return _h1_param_sum_ ( h1              ,
                                func            ,
                                H_fit           ,
                                opts   = opts   ,
                                xmin   = xmin   ,
                                xmax   = xmax   ,
                                fixes  = fixes  ,
                                params = params ,
                                limits = limits , 
                                refit  = refit  )

    
    for t in ( ROOT.TH1F , ROOT.TH1D ) :

        t.cosine_sum     = _h1_cosine_sum_
        t.cosine         = _h1_cosine_

    _new_methods_ .append ( _h1_cosine_sum_   )
    _new_methods_ .append ( _h1_cosine_       )

    # ==========================================================================
except ImportError :     # =====================================================
    # ==========================================================================
    pass
    # ==========================================================================

# =============================================================================
## represent 1D-histo as plain vanilla polynomial
#  @code
#  h = ...                  ## the histogram
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
def _h1_polinomial_ ( h1               ,
                      degree           ,
                      opts   = 'SQ0'   ,
                      xmin   = inf_neg ,
                      xmax   = inf_pos ,
                      fixes  = ()      ,
                      params = ()      ,
                      limits = ()      ,
                      refit  = 1       ) :
    """ Represent histo as plain vanilla polynomial    
    >>> h = ... # the histogram    
    >>> b = h.polinomial ( 5 )  ## make a fit... 
    
    >>> tf1        = b [ 0 ]    ## TF1 object
    >>> obj        = b [ 1 ]    ## helper object 
    >>> fun        = b [ 2 ]    ## underlying normalzed C++ object 
    >>> fit_result = b [ 3 ]    ## fit result & status
    >>> norm       = b [ 4 ]    ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) 
    """
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    func  = Ostap.Math.Polynomial ( degree , xmin , xmax ) 
    #
    my = h1.accumulate().value()/h1.bins()
    func.setPar ( 0, my )
    #
    ## make a fit 
    if not params : params = tuple ( [ p for p in func.pars() ] ) 
    from ostap.fitting.param import H_fit
    return _h1_param_sum_ ( h1              ,
                            func            ,
                            H_fit           ,
                            opts   = opts   ,
                            xmin   = xmin   ,
                            xmax   = xmax   ,
                            fixes  = fixes  ,
                            params = params ,
                            limits = limits ,
                            refit  = refit  )

# =============================================================================
## represent 1D-histo as B-spline
#  @code
#  h = ...                  ## the histogram
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
def _h1_bspline_ ( h1               ,
                   degree = 3       ,
                   knots  = 3       ,
                   opts   = 'SQ0'   ,
                   xmin   = inf_neg ,
                   xmax   = inf_pos ,
                   fixes  = ()      ,
                   params = ()      ,
                   limits = ()      ,
                   refit  = 1       ) :
    """ Represent histo as B-spline polynomial    
    >>> h = ... # the histogram
    >>> b = h.bSpline ( degree = 3 , knots = 3  )
    >>> b = h.bSpline ( degree = 3 , knots = [ 0.1 , 0.2, 0.8, 0.9 ]  )
    
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
    if isinstance ( knots , integer_types ) and 0 <= knots :
        func = Ostap.Math.BSpline ( xmin , xmax , knots , degree )
    else :
        from ostap.math.base import doubles
        _knots = doubles ( xmin , xmax ) 
        for k in knots : _knots.push_back( k )
        func = Ostap.Math.BSpline ( _knots , degree )
        
    ## make a fit 
    if not params : params = tuple ( [ p for p in func.pars() ] ) 
    from ostap.fitting.param import H_fit
    return _h1_param_sum_ ( h1              ,
                            func            ,
                            H_fit           ,
                            opts   = opts   ,
                            xmin   = xmin   ,
                            xmax   = xmax   ,
                            fixes  = fixes  ,
                            params = params ,
                            limits = limits ,
                            refit  = refit  )

# =============================================================================
## represent 1D-histo as POSITIVE bernstein polynomial
#  @code
#  h = ...                  ## the histogram
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
def _h1_positive_ ( h1               ,
                    degree           ,
                    opts   = 'SQ0'   ,
                    xmin   = inf_neg ,
                    xmax   = inf_pos ,
                    fixes  = ()      ,
                    params = ()      ,
                    limits = ()      , 
                    refit  = 1       ) :
    """ Represent histo as Positive Bernstein polynomial
    >>> h = ...              ## the histogram
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
    func = Ostap.Math.Positive ( degree , xmin , xmax )
    #
    ## make a fit 
    from ostap.fitting.param import H_Nfit
    return _h1_param_sum_ ( h1              ,
                            func            ,
                            H_Nfit          ,
                            opts   = opts   ,
                            xmin   = xmin   ,
                            xmax   = xmax   ,
                            fixes  = fixes  ,
                            params = params ,
                            limits = limits ,
                            refit  = refit  )

# =============================================================================
## represent 1D-histo as POSITIVE EVEN bernstein polynomial
#  @code
#  h = ...                  ## the histogram
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
def _h1_positiveeven_ ( h1 ,
                        degree           ,
                        opts   = 'SQ0'   ,
                        xmin   = inf_neg ,
                        xmax   = inf_pos ,
                        fixes  = ()      ,
                        params = ()      ,
                        limits = ()      ,
                        refit  = 1       ) :
    """ Represent histo as Positive Even Bernstein polynomial
    
    >>> h = ... # the histogram
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
    func = Ostap.Math.PositiveEven ( degree , xmin , xmax )
    # 
    ## make a fit 
    from ostap.fitting.param import H_Nfit
    return _h1_param_sum_ ( h1              ,
                            func            ,
                            H_Nfit          ,
                            opts   = opts   ,
                            xmin   = xmin   ,
                            xmax   = xmax   ,
                            fixes  = fixes  ,
                            params = params ,
                            limits = limits ,  
                            refit  = refit  )

# =============================================================================
## represent 1D-histo as MONOTONIC bernstein polynomial
#  @code
#  h = ...                  ## the histogram
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
def _h1_monotonic_ ( h1                   ,
                     degree               ,
                     increasing = True    , 
                     opts       = 'SQ0'   ,
                     xmin       = inf_neg ,
                     xmax       = inf_pos ,
                     fixes      = ()      ,
                     params     = ()      ,
                     limits     = ()      ,
                     refit      = 1       ) :
    """ Represent histo as Monotonic Bernstein polynomial
    
    >>> h = ...           ## the histogram
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
    func = Ostap.Math.Monotonic ( degree , xmin , xmax , increasing )
    # 
    ## make a fit 
    from ostap.fitting.param import H_Nfit
    return _h1_param_sum_ ( h1              ,
                            func            ,
                            H_Nfit          ,
                            opts   = opts   ,
                            xmin   = xmin   ,
                            xmax   = xmax   ,
                            fixes  = fixes  ,
                            params = params ,
                            limits = limits ,
                            refit  = refit  )


# =============================================================================
## represent 1D-histo as MONOTONIC CONVEX/CONCAVE bernstein polynomial
#  @code
#  h = ...                  ## the histogram
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
def _h1_convex_ ( h1                   , 
                  degree               ,
                  increasing = True    ,
                  convex     = True    ,
                  opts       = 'SQ0'   ,
                  xmin       = inf_neg ,
                  xmax       = inf_pos ,
                  fixes      = ()      ,
                  params     = ()      ,
                  limits     = ()      , 
                  refit      = 1       ) :
    """ Represent histo as Monotonic Convex/Concave  Bernstein polynomial    
    >>> h = ...           ## the histogram
    >>> b = h.convex ( 5 , increasing = True , convex = False )    
    
    >>> tf1        = b[0] ## TF1 object
    >>> obj        = b[1] ## helper object 
    >>> fun        = b[2] ## underlying normalzed C++ object 
    >>> fit_result = b[3] ## fit result & status
    >>> norm       = b[4] ## normalization 
    
    >>> x = ...
    >>> print 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) 
    >>> print 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) 
    >>> h = ... # the histogram
    """
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    func = Ostap.Math.Convex ( degree , xmin , xmax , increasing , convex )
    # 
    ## make a fit 
    from ostap.fitting.param import H_Nfit
    return _h1_param_sum_ ( h1              ,
                            func            ,
                            H_Nfit          ,
                            opts   = opts   ,
                            xmin   = xmin   ,
                            xmax   = xmax   ,
                            fixes  = fixes  ,
                            params = params ,
                            limits = limits ,
                            refit  = refit  )

# =============================================================================
## represent 1D-histo as CONVEX bernstein polynomial
#  @code
#  h = ...                  ## the histogram
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
def _h1_convexpoly_ ( h1                   ,
                      degree               ,
                      opts       = 'SQ0'   ,
                      xmin       = inf_neg ,
                      xmax       = inf_pos ,
                      fixes      = ()      ,
                      params     = ()      ,
                      limits     = ()      ,
                      refit      = 1       ) :

    """ Represent histo as Convex Bernstein polynomial
    >>> h = ...           ## the histogram
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
    func = Ostap.Math.ConvexOnly ( degree , xmin , xmax , True )
    # 
    ## make a fit 
    from ostap.fitting.param import H_Nfit
    return _h1_param_sum_ ( h1              ,
                            func            ,
                            H_Nfit          ,
                            opts   = opts   ,
                            xmin   = xmin   ,
                            xmax   = xmax   ,
                            fixes  = fixes  ,
                            params = params ,
                            limits = limits ,  
                            refit  = refit  )

# =============================================================================
## represent 1D-histo as CONCAVE bernstein polynomial
#  @code
#  h = ...                  ## the histogram
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
def _h1_concavepoly_ ( h1                   ,
                       degree               ,
                       opts       = 'SQ0'   ,
                       xmin       = inf_neg ,
                       xmax       = inf_pos ,
                       fixes      = ()      ,
                       params     = ()      ,
                       limits     = ()      ,
                       refit      = 1       ) :
    
    """ Represent histo as Concave  Bernstein polynomial

    >>> h = ...           ## the histogram
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
    func = Ostap.Math.ConvexOnly ( degree , xmin , xmax , False )
    # 
    ## make a fit 
    from ostap.fitting.param import H_Nfit
    return _h1_param_sum_ ( h1              ,
                            func            ,
                            H_Nfit          ,
                            opts   = opts   ,
                            xmin   = xmin   ,
                            xmax   = xmax   ,
                            fixes  = fixes  ,
                            params = params ,
                            limits = limits ,
                            refit  = refit  )

# =============================================================================
## represent 1D-histo as positive B-spline
#  @code
#  h = ...                  ## the histogram
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
def _h1_pspline_ ( h1               ,
                   degree = 3       ,
                   knots  = 3       ,
                   opts   = 'SQ0'   ,
                   xmin   = inf_neg ,
                   xmax   = inf_pos , 
                   fixes  = ()      ,
                   params = ()      ,
                   limits = ()      ,
                   refit  = 1       ) :
    """ Represent histo as positive B-spline 
    
    >>> h  = ... # the histogram
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
    if isinstance ( knots , integer_types ) and 0 <= knots :
        func = Ostap.Math.PositiveSpline ( xmin , xmax , knots , degree  )
    else :
        from ostap.math.base import doubles
        _knots = doubles ( xmin , xmax ) 
        for k in knots : _knots.push_back( k )
        func = Ostap.Math.PositiveSpline ( _knots , degree )
    #
    ## make a fit 
    from ostap.fitting.param import H_Nfit
    return _h1_param_sum_ ( h1              ,
                            func            ,
                            H_Nfit          ,
                            opts   = opts   ,
                            xmin   = xmin   ,
                            xmax   = xmax   ,
                            fixes  = fixes  ,
                            params = params ,
                            limits = limits ,
                            refit  = refit  )

# =============================================================================
## represent 1D-histo as positive monotonic spline
#  @code
#  h = ...                  ## the histogram
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
def _h1_mspline_ ( h1                   ,
                   degree     = 3       ,
                   knots      = 3       ,
                   increasing = True    ,                
                   opts       = 'SQ0'   ,
                   xmin       = inf_neg ,
                   xmax       = inf_pos , 
                   fixes      = ()      ,
                   params     = ()      ,
                   limits     = ()      ,
                   refit      = 1       ) :
    """ Represent histo as positive monotonic  spline 
    
    >>> h  = ... # the histogram
    
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
    if isinstance ( knots , integer_types ) and 0 <= knots :
        func = Ostap.Math.MonotonicSpline ( xmin , xmax , knots , degree , increasing )
    else :
        from ostap.math.base import doubles
        _knots = doubles ( xmin , xmax ) 
        for k in knots : _knots.push_back( k )
        func = Ostap.Math.MonotonicSpline ( knots , degree , increasing )
    #
    ## make a fit 
    from ostap.fitting.param import H_Nfit
    return _h1_param_sum_ ( h1              ,
                            func            ,
                            H_Nfit          ,
                            opts   = opts   ,
                            xmin   = xmin   ,
                            xmax   = xmax   ,
                            fixes  = fixes  ,
                            params = params ,
                            limits = limits ,
                            refit  = refit  )

# =============================================================================
## represent 1D-histo as monotonic convex/concave spline
#  @code
#  h = ...                  ## the histogram
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
                   opts       = 'SQ0'   ,
                   xmin       = inf_neg ,
                   xmax       = inf_pos , 
                   fixes      = ()      ,
                   params     = ()      ,
                   limits     = ()      ,
                   refit      = 1       ) :
    """ Represent histo as positive monotonic convex/concave spline  
    
    >>> h = ... # the histogram
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
    if isinstance ( knots , integer_types ) and 0 <= knots :
        func = Ostap.Math.ConvexSpline ( xmin , xmax , knots , degree , increasing , convex  )
    else :
        from ostap.math.base import doubles
        _knots = doubles ( xmin , xmax ) 
        for k in knots : _knots.push_back( k )
        func   = Ostap.Math.ConvexSpline ( _knots , order , increasing , convex )
        
    ## make a fit 
    from ostap.fitting.param import H_Nfit
    return _h1_param_sum_ ( h1              ,
                            func            ,
                            H_Nfit          ,
                            opts   = opts   ,
                            xmin   = xmin   ,
                            xmax   = xmax   ,
                            fixes  = fixes  ,
                            params = params ,
                            limits = limits ,
                            refit  = refit  )

# =============================================================================
## represent 1D-histo as positive convex spline
#  @code
#  h = ...                  ## the histogram
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
def _h1_convexspline_ ( h1                  ,
                        degree     = 3       ,
                        knots      = 3       ,
                        opts       = 'SQ0'   ,
                        xmin       = inf_neg ,
                        xmax       = inf_pos , 
                        fixes      = ()      ,
                        params     = ()      ,
                        limits     = ()      ,
                        refit      = 1       ) :
    """ Represent histo as positive convex spline  
    
    >>> h = ... # the histogram
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
    if isinstance ( knots , integer_types ) and 0 <= knots :
        func = Ostap.Math.ConvexOnlySpline ( xmin , xmax , knots , degree , True )
    else :
        from ostap.math.base import doubles
        _knots = doubles ( mn , mx ) 
        for k in knots : _knots.push_back( k )
        func   = Ostap.Math.ConvexOnlySpline ( _knots , order , True )
    ##
    ## make a fit 
    from ostap.fitting.param import H_Nfit
    return _h1_param_sum_ ( h1              ,
                            func            ,
                            H_Nfit          ,
                            opts   = opts   ,
                            xmin   = xmin   ,
                            xmax   = xmax   ,
                            fixes  = fixes  ,
                            params = params ,
                            limits = limits ,
                            refit  = refit  )

# =============================================================================
## represent 1D-histo as positive concave spline
#  @code
#  h = ...                  ## the histogram
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
                         degree     = 3       ,
                         knots      = 3       ,
                         opts       = 'SQ0'   ,
                         xmin       = inf_neg ,
                         xmax       = inf_pos , 
                         fixes      = ()      ,
                         params     = ()      ,
                         limits     = ()      ,
                         refit      = 1       ) :
    """ Represent histo as positive convcave spline  
    
    >>> h = ... # the histogram
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
    if isinstance ( knots , integer_types ) :
        func = Ostap.Math.ConvexOnlySpline ( xmin , xmax , knots , degree , False )
    else :
        from ostap.math.base import doubles
        _knots = doubles ( xmin , xmax ) 
        for k in knots : _knots.push_back( k )
        func   = Ostap.Math.ConvexOnlySpline ( _knots , order , False )
        
    ##
    ## make a fit 
    from ostap.fitting.param import H_Nfit
    return _h1_param_sum_ ( h1              ,
                            func            ,
                            H_Nfit          ,
                            opts   = opts   ,
                            xmin   = xmin   ,
                            xmax   = xmax   ,
                            fixes  = fixes  ,
                            params = params ,
                            limits = limits ,
                            refit  = refit  )

# =============================================================================

# =============================================================================
## represent 1D-histo as positive Karlin-Shapley polymonial 
#  @code
#  h = ...                  ## the histogram
#  b = h.karlin_shapley ( degree = 3 )
#
#  tf1        = b[0]        ## TF1 object
#  obj        = b[1]        ## helper object 
#  fun        = b[2]        ## uderlying normalzed C++ object 
#  fit_result = b[3]        ## fit result & status
#  norm       = b[4]        ## normalization
# 
#  x = ...
#  print ( 'TF1(%s) = %s' % ( x , tf1 ( x )        ) ) 
#  print ( 'fun(%s) = %s' % ( x , fun ( x ) * norm ) ) 
#  @endcode 
def _h1_karlinshapley_ ( h1                   ,
                         degree     = 3       ,
                         opts       = 'SQ0'   ,
                         xmin       = inf_neg ,
                         xmax       = inf_pos , 
                         fixes      = ()      ,
                         params     = ()      ,
                         limits     = ()      ,
                         refit      = 1       ) :
    """ Represent histo as positive Karlin-Shapley polynomial
    
    >>> h = ... # the histogram
    >>> b = h.karlin_shapley ( degree = 3 , knots = 3 )
    
    >>> tf1        = b[0] ## TF1 object
    >>> obj        = b[1] ## helper object 
    >>> fun        = b[2] ## underlying normalzed C++ object 
    >>> fit_result = b[3] ## fit result & status
    >>> norm       = b[4] ## normalization 
    
    >>> x = ... 
    >>> print ( 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) )
    >>> print ( 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) )        
    """
    ##
    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    #
    func   = Ostap.Math.KarlinShapley ( degree , xmin , xmax )
    
    my = h1.accumulate().value()/h1.bins()
    func.setPar ( 0, my )   

    newlims = tuple ( l for l in limits ) 
    if degree :
        lims = {}
        for l in newlims : lims [ l[0] ] = l
        for i in range ( 1 , degree + 1 ) :
            if not i in lims : lims [i] = i , -3 * math.pi , 4 * math.pi

        newlims = [] 
        for k in sorted ( lims ) : newlims.append ( lims[k] )
        newlims = tuple ( newlims )
        
    ## make a fit 
    from ostap.fitting.param import H_fit
    return _h1_param_sum_ ( h1               ,
                            func             ,
                            H_fit            ,
                            opts   = opts    ,
                            xmin   = xmin    ,
                            xmax   = xmax    ,
                            fixes  = fixes   ,
                            params = params  ,
                            limits = newlims ,
                            refit  = refit   )

# =============================================================================
## represent 1D-histo as positive Karlin-Studden polymonial 
#  @code
#  h = ...                  ## the histogram
#  b = h.karlin_studden( degree = 3 )
#
#  tf1        = b[0]        ## TF1 object
#  obj        = b[1]        ## helper object 
#  fun        = b[2]        ## uderlying normalzed C++ object 
#  fit_result = b[3]        ## fit result & status
#  norm       = b[4]        ## normalization
# 
#  x = ...
#  print ( 'TF1(%s) = %s' % ( x , tf1 ( x )        ) ) 
#  print ( 'fun(%s) = %s' % ( x , fun ( x ) * norm ) ) 
#  @endcode 
def _h1_karlinstudden_ ( h1               ,
                         degree     = 3       ,
                         opts       = 'SQ0'   ,
                         xmin       = inf_neg ,
                         xmax       = inf_pos ,
                         scale      = None    ,
                         fixes      = ()      ,
                         params     = ()      ,
                         limits     = ()      ,
                         refit      = 1       ) :
    """ Represent histo as positive Karlin-Studden polynomial
    
    >>> h = ... # the histogram
    >>> b = h.karlin_shapley ( degree = 3 , knots = 3 )
    
    >>> tf1        = b[0] ## TF1 object
    >>> obj        = b[1] ## helper object 
    >>> fun        = b[2] ## underlying normalzed C++ object 
    >>> fit_result = b[3] ## fit result & status
    >>> norm       = b[4] ## normalization 
    
    >>> x = ... 
    >>> print ( 'TF1(%s) = %s' % ( x ,        tf1 ( x ) ) )
    >>> print ( 'FUN(%s) = %s' % ( x , norm * fun ( x ) ) )        
    """

    xmin = max ( xmin , h1.xmin() ) 
    xmax = min ( xmax , h1.xmax() )  
    
    if scale is None : scale = xmax - xmin 
    assert isinstance ( scale , num_types ) and 0 < scale < inf_pos , \
           'Invalid scale parameter'

    #
    func   = Ostap.Math.KarlinStudden ( degree , xmin , scale )
    
    my = h1.accumulate().value()/h1.bins()
    func.setPar ( 0, my )   

    newlims = tuple ( l for l in limits ) 
    if degree :
        lims = {}
        for l in newlims : lims [ l[0] ] = l
        for i in range ( 1 , degree + 1 ) :
            if not i in lims : lims [i] = i , -3 * math.pi , 4 * math.pi

        newlims = [] 
        for k in sorted ( lims ) : newlims.append ( lims[k] )
        newlims = tuple ( newlims )
        
    ## make a fit 
    from ostap.fitting.param import H_fit
    return _h1_param_sum_ ( h1               ,
                            func             ,
                            H_fit            ,
                            opts   = opts    ,
                            xmin   = xmin    ,
                            xmax   = xmax    ,
                            fixes  = fixes   ,
                            params = params  ,
                            limits = newlims ,
                            refit  = refit   )


# =============================================================================
## decorate histograms 
for t in ( ROOT.TH1D , ROOT.TH1F ) :
    t.bernstein      = _h1_bernstein_
    t.bezier         = _h1_bernstein_ ## ditto 
    t.bernsteineven  = _h1_bernsteineven_
    t.beziereven     = _h1_bernsteineven_ ## ditto 
    t.chebyshev      = _h1_chebyshev_
    t.legendre       = _h1_legendre_
    t.polynomial     = _h1_polinomial_
    t.rational       = _h1_rational_
    t.brational      = _h1_brational_
    t.positive       = _h1_positive_
    t.positiveeven   = _h1_positiveeven_
    t.monotonic      = _h1_monotonic_
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
    t.legendre_fast  = _h1_legendre_fast_
    t.karlin_shapley = _h1_karlinshapley_ 
    t.karlin_studden = _h1_karlinstudden_ 

_new_methods_ += [
    _h1_bernstein_     ,
    _h1_bernsteineven_ ,
    _h1_chebyshev_     ,
    _h1_legendre_      ,
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
    _h1_legendre_fast_ ,
    _h1_karlinshapley_ , 
    _h1_karlinstudden_ ,  
    _h1_rational_      ,  
    _h1_brational_     ,  
    ]

# =============================================================================
## parameterize positive histogram with certain PDF
#   @param h1 1D-histogram
#   @param pdf_type Either the type of fit-PDF  or ready-to-use PDF
#   @param pars     tuple of parameters for PDF creation in case `pdf_type` is `type` 
#   @param args     addtitional positional argumetts for `PDF.fitHisto ( histo , *args , **kwargs )`
#   @param pdfpars  additional keyword arguments for PDF creation in case `pdf_type` is `type`
#   @param kwargs   addtitional keyword argument for `PDF.fitHisto ( histo , *args , **kwargs )`
#   @return ParamPDFInfo 
def _h1_pdf_ ( h1 , pdf_type , pars , *args , pdfpars = {} , **kwargs ) :
    """ Parameterize the positive histogram with certain PDF
    
    - `h1`       : 1D-histogram
    - `pdf_type` : either the type of the fit-PDF or ready-to-use PDF1 instance 
    - `pars`     : tuple of parameters for PDF creation in case `pdf_type` is `type` 
    - `args`     : addtitional positional argumetts for `PDF.fitHisto ( histo , *args , **kwargs )`
    - `pdfpars`  : additional  keyword arguments for PDF creation in case `pdf_type` is `type`
    - `kwargs`   : addtitional keyword argument for `PDF.fitHisto ( histo , *args , **kwargs )`

    return `ParamPDFInfo` object 

    """
    assert isinstance ( h1 , ROOT.TH1 ) and 1 == h1.dim() ,\
        "Invalid 1D-histogram type: %s" % typename ( h1 ) 
    ##
    mn , mx = h1.minmax()
    assert 0 <= mn.value () and 0 < mx.value() , "Histo goes to negative: %s/%s" % ( mn , mx )
    ##
    if not hasattr ( h1 , 'xvar' ) :
        vname   = 'x' + h1.GetName()
        vname   = vname.strip().replace(' ', '_').replace('*', '_')
        h1.xvar = ROOT.RooRealVar ( vname , 'xvar(%s)' % h1.GetName() , *h1.xminmax() )
        
    ## create pdf
    from ostap.fitting.fithelpers import VarMaker
    from ostap.fitting.pdfbasic   import PDF1
    ## 
    name   = VarMaker.generate_name (  'pdf' , '' , h1.GetName() )
    ##
    ## ready to use PDF is supplied 
    if   isinstance ( pdf_type , PDF1 )  :
        ## use the supplied PDF 
        if pars    : logger.warning ( "Ignore specified `pars` arguments: %s"   % str ( pars    ) )
        if pdfpars : logger.warning ( "Ignore specified `pdfpars` argument: %s" % str ( pdfpars ) )
        pdf = pdt_type
    elif isinstance ( pdf_type , type ) and issubclass ( pdf_type , PDF1 ) :
        ## create new pdf 
        pdf    = pdf_type     ( name , h1.xvar , *pars , **pdfpars )
    else :
        raise TypeError( "Invalid `pdf_type/pars/pdfpars` setting!" ) 

    ## refit ? 
    refit = kwargs.pop ( 'refit' , 8 ) 
    # ===========================================================================
    ## fit the histogram
    # ===========================================================================
    try : # =====================================================================
        # =======================================================================
        r , f  = pdf.fitHisto ( h1 , *args, refit = refit , **kwargs )
        # =======================================================================
    except Exception : # ========================================================
        # =======================================================================
        logger.error ( 'Exception from %s/%s fitHisto' % ( typename ( pdf     ) ,
                                                           typename ( pdf.pdf ) ) , exc_info = True )
        raise

    ## 
    if hasattr ( pdf.pdf , 'setPars'  ) : pdf.pdf.setPars()
    ## 
    if hasattr ( pdf.pdf , 'function' ) : func = pdf.pdf.function()
    else                                : func = lambda x : pdf ( x )
    ##
    from ostap.fitting.roofit import PDF_fun
    pdf_fun = PDF_fun( pdf.pdf , h1.xvar , *h1.xminmax() )
    ##
    norm   = h1.integral() 
    ##
    params = ParamPDFInfo ( r , pdf , func , norm , pdf_fun , f )
    h1._param_PDF_info = params
    return params 
    
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
def _h1_pdf_positive_ ( h1 , degree , *args , pdfpars = {} , **kwargs ) :
    """ Parameterize/fit histogram with the positive polynomial
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
    return _h1_pdf_ ( h1 , PolyPos_pdf , (degree,) , *args , pdfpars = pdfpars , **kwargs )

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
def _h1_pdf_even_ ( h1 , halfdegree , *args , pdfpars = {} , **kwargs ) :
    """ Parameterize/fit histogram with the positive even polynomial
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
    return _h1_pdf_ ( h1 , PolyEven_pdf , (halfdegree,) , *args , pdfpars = pdfpars , **kwargs )

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
def _h1_pdf_monotonic_ ( h1 , degree , increasing , *args , pdfpars = {} , **kwargs ) :
    """ Parameterize/fit histogram with the monotonic positive polynomial
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
    return _h1_pdf_ ( h1 , Monotonic_pdf , (degree,increasing) , *args , pdfpars = {} , **kwargs )

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
def _h1_pdf_increasing_ ( h1 , degree , *args , pdfpars = {} , **kwargs ) :
    """ Parameterize/fit histogram with the monotonic positive polynomial
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
    return _h1_pdf_monotonic_ ( h1 , degree , True , *args , pdfpars = pdfpars , **kwargs ) 

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
def _h1_pdf_decreasing_ ( h1 , degree , *args , pdfpars = {} , **kwargs ) :
    """ Parameterize/fit histogram with the monotonic positive polynomial
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
    return _h1_pdf_monotonic_ ( h1 , degree , False , *args , pdfpars = pdfpars , **kwargs ) 

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
def _h1_pdf_convex_ ( h1 , degree , increasing , *args , pdfpars = {} , **kwargs ) :
    """ Parameterize/fit histogram with convex polynomial
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
    return _h1_pdf_ ( h1 , Convex_pdf , (degree,increasing,True) , *args , pdfpars = {} , **kwargs )

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
def _h1_pdf_convex_increasing_ ( h1 , degree , *args , pdfpars = {} , **kwargs ) :
    """ Parameterize/fit histogram with convex increasing polynomial
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
    return _h1_pdf_convex_ ( h1 , degree , True , *args , pdfpars = pdfpars , **kwargs )

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
def _h1_pdf_convex_decreasing_ ( h1 , degree , *args , pdfpars = {} , **kwargs ) :
    """ Parameterize/fit histogram with convex decreasing polynomial
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
    return _h1_pdf_convex_ ( h1 , degree , False , *args , pdfpars = pdfpars , **kwargs )

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
def _h1_pdf_concave_ ( h1 , degree , increasing , *args , pdfpars = {} , **kwargs ) :
    """ Parameterize/fit histogram with concave polynomial
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
    return _h1_pdf_ ( h1 , Convex_pdf , (degree,increasing,False) , *args , pdfpars = pdfpars , **kwargs )

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
def _h1_pdf_concave_increasing_ ( h1 , degree , *args , pdfpars = {} , **kwargs ) :
    """ Parameterize/fit histogram with concave increasing polynomial
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
    return _h1_pdf_concave_ ( h1 , degree , True , *args , pdfpars = {} , **kwargs )

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
def _h1_pdf_concave_decreasing_ ( h1 , degree , *args , pdfpars = {} , **kwargs ) :
    """ Parameterize/fit histogram with concave decreasing polynomial
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
    return _h1_pdf_concave_ ( h1 , degree , False , *args , pdfpars = pdfpars , **kwargs )

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
def _h1_pdf_convexpoly_ ( h1 , degree , *args , pdfpars = {} , **kwargs ) :
    """ Parameterize/fit histogram with convex polynomial
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
    return _h1_pdf_ ( h1 , ConvexOnly_pdf , (degree,True) , *args , pdfpars = pdfpars , **kwargs )

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
def _h1_pdf_concavepoly_ ( h1 , degree , *args , pdfpars = {} , **kwargs ) :
    """ Parameterize/fit histogram with convex polynomial
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
    return _h1_pdf_ ( h1 , ConvexOnly_pdf , (degree,False) , *args , pdfpars = pdfpars , **kwargs )

# =============================================================================
## parameterize/fit histogram with the rational function 
#  @code
#  h1 = ...
#  results = h1.pdf_rational ( 3 , 3 )
#  results = h1.pdf_rational ( 3 , 3 , draw = True , silent = True )
#  print     results[ 0]
#  pdf     = results[ 1]
#  func    = results[ 2]
#  norm    = results[ 3]
#  pdf_fun = results[ 4] 
#  frame   = results[-1] 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2023-09-19
def _h1_pdf_rational_ ( h1 , p , d  , *args , pdfpars = {} , **kwargs ) :
    """ Parameterize/fit histogram with the rational function 
    >>> h1 = ...
    >>> results = h1.pdf_rational ( 3 , 3 )
    >>> results = h1.pdf_rational ( 3 , 3 , draw = 3 , silent = True )
    >>> print     results[ 0] ## fit results 
    >>> pdf     = results[ 1] ## get PDF
    >>> func    = results[ 2] ## normalized function  
    >>> norm    = results[ 3] ## normalization
    >>> pdf_fun = results[ 4] ## pdf-function (normalized)
    >>> frame   = results[-1] ## frame/RooPlot 
    """

    assert isinstance ( p , integer_types ) and 0 <= p , "Invalid p: %s" % p
    assert isinstance ( d , integer_types ) and 0 <= d , "Invalid d: %s" % d

    from ostap.fitting.background import Rational_pdf
    from ostap.fitting.variables  import FIXVAR 
    
    mn , mx = h1.minmax()
    if mn.value() < 0 or mx.value() <= 0 :
        raise AttributeError("Histo goes to negative %s/%s" % ( mn , mx ) )
    ##
    if not hasattr ( h1 , 'xvar' ) :
        h1.xvar = ROOT.RooRealVar ( 'x' + h1.GetName() , 'xvar(%s)' % h1.GetName() , *h1.xminmax() )

    ## create pdf
    from ostap.fitting.fithelpers import VarMaker 
    name   = VarMaker.generate_name (  'pdf' , '' , h1.GetName() )
    ##
    ## pdf    = pdf_type     ( name , h1.xvar , *pars )
    pdf    = Rational_pdf ( name , h1.xvar , p , d , *pdfpars )

    # ===========================================================================
    ## fit the histogram
    # ===========================================================================
    try : # =====================================================================
        # =======================================================================
        ## all pars 
        apars = tuple ( v for v in pdf.pars )
        pvars = tuple ( v for v in apars [   : p ] ) ## numerator pars 
        qvars = tuple ( v for v in apars [ p :   ] ) ## denominator parameters 
        
        for qv in qvars :
            if hasattr ( qv , 'setValue' ) and not qv.isConstant() : qv.setValue ( 0 ) 
                
        kw = {}
        kw.update ( kwargs )
        kw [ 'draw'   ] = False
        kw [ 'silent' ] = True
        kw [ 'refit'  ] = 5 

        if pvars and qvars :
            
            ## (1) fix denominator and fit with numerator-only
            with FIXVAR ( qvars ) : r , _  = pdf.fitHisto ( h1 , *args , **kw )        
            ## (2) fix numerator and re-fit with denominator-only
            with FIXVAR ( pvars ) : r , _  = pdf.fitHisto ( h1 , *args , **kw )
            ## (3) fix denominator and re-fit with numerator-only
            with FIXVAR ( qvars ) : r , _  = pdf.fitHisto ( h1 , *args , **kw  )            
                
        ## (4) the final fit
        r , f  = pdf.fitHisto ( h1 , *args, **kwargs )
            
        # =======================================================================
    except Exception : # ========================================================
        # =======================================================================
        logger.error ( 'Exception from %s/%s fitHisto' % ( typename ( pdf     ) ,
                                                           typename ( pdf.pdf ) ) , exc_info = True )

        pdf.pdf.setPars ()

        fun  = pdf.pdf.function()
        fun1 = fun.numerator   ()
        fun2 = fun.denominator () 

        from ostap.utils.ranges import vrange

        print ('histo' , h1 )
        
        print ('FUN/1' , fun1 )
        print ('FUN/2' , fun2 )
        print ('FUN  ' , fun  )

        xmin, xmax = h1.xminmax()
        
        print ('FUN/1' , list ( fun1 ( x ) for x in vrange ( xmin , xmax , 20 ) ) )
        print ('FUN/2' , list ( fun2 ( x ) for x in vrange ( xmin , xmax , 20 ) ) )
        print ('FUN  ' , list ( fun  ( x ) for x in vrange ( xmin , xmax , 20 ) ) )
        
        pdf.numerator  .draw ()
        pdf.denominator.draw ()
        
        raise
    
    ## 
    if hasattr ( pdf.pdf , 'setPars'  ) : pdf.pdf.setPars()
    ## 
    if hasattr ( pdf.pdf , 'function' ) : func = pdf.pdf.function()
    else                                : func = lambda x : pdf ( x )
        
    ##
    from ostap.fitting.roofit import PDF_fun
    pdf_fun = PDF_fun( pdf.pdf , h1.xvar , *h1.xminmax() )
    ##
    norm   = VE ( h1.integrate() ).value() 
    ##
    params = ParamPDFInfo ( r , pdf , func , norm , pdf_fun , f )
    h1._param_PDF_info = params
    return params 

    
    ## return _h1_pdf_ ( h1 , Rational_pdf , ( p , q ) , *args , **kwargs )

# =============================================================================

for t in ( ROOT.TH1D , ROOT.TH1F ) :
    t.pdf_positive           = _h1_pdf_positive_
    t.pdf_positiveeven       = _h1_pdf_even_
    t.pdf_even               = _h1_pdf_even_
    t.pdf_monotonic          = _h1_pdf_monotonic_
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
    t.pdf_rational           = _h1_pdf_rational_

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
        _h1_pdf_rational_     ,        
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
def _h1_pdf_pspline_ ( h1 , spline , *args , pdfpars = {} , **kwargs ) :
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
        PS     = Ostap.Math.PositiveSpline
        spline = PS ( h1.xmin() , h1.xmax() , spline[1] , spline[0] )
    #
    from ostap.fitting.background import PSpline_pdf
    return _h1_pdf_ ( h1 , PSpline_pdf , ( spline , ) , *args , pdfpars = pdfpars , **kwargs )


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
def _h1_pdf_mspline_ ( h1 , spline , *args , pdfpars = {} , **kwargs ) :
    """ Parameterize/fit histogram with monotonic positive b-spline 
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
        MS     = Ostap.Math.MonotonicSpline
        spline = MS ( h1.xmin() , h1.xmax() , spline[1] , spline[0] , spline[2] )
        #
    from ostap.fitting.background import MSpline_pdf
    return _h1_pdf_ ( h1 , MSpline_pdf , ( spline , ) , *args , pdfpars = pdfpars , **kwargs )

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
def _h1_pdf_cspline_ ( h1 , spline , *args , pdfpars = {} , **kwargs ) :
    """ Parameterize/fit histogram with convex/concave montonic positive b-spline 
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
        CS     = Ostap.Math.ConvexSpline
        spline = CS ( h1.xmin() , h1.xmax() , spline[1] , spline[0] , spline[2] , spline[3] )
    #
    from ostap.fitting.background import CSpline_pdf
    return _h1_pdf_ ( h1 , CSpline_pdf , ( spline , ) , *args , pdfpars = pdfpars , **kwargs )


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
def _h1_pdf_convexSpline_ ( h1 , spline , *args , pdfpars = {} , **kwargs ) :
    """ Parameterize/fit histogram with convex positive b-spline 
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
        CS     = Ostap.Math.ConvexOnlySpline
        spline = CS ( h1.xmin() , h1.xmax() , spline[1] , spline[0] , True )
        #
    from ostap.fitting.background import CPSpline_pdf
    return _h1_pdf_ ( h1 , CPSpline_pdf , ( spline , ) , *args , pdfpars = pdfpars , **kwargs )


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
def _h1_pdf_concaveSpline_ ( h1 , spline , *args , pdfpars = {} , **kwargs ) :
    """ Parameterize/fit histogram with positive concave b-spline 
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
        CS     = Ostap.Math.ConvexOnlySpline
        spline = CS ( h1.xmin() , h1.xmax() , spline[1] , spline[0] , False )
        #
    from ostap.fitting.background import CPSpline_pdf
    return _h1_pdf_ ( h1 , CPSpline_pdf , ( spline , ) , *args , pdfpars = pdfpars , **kwargs )

# =============================================================================
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
#  fsum   = histo.legendre_sum_orig ( 4 )
#  print fsum
#  x = ...
#  print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
#  @endcode 
#  @see Ostap::Math::LegendreSum
#  @author Vanya Belyaev Ivan.Belyaev@iter.ru
#  It is not very CPU efficient, but stable enough...
#  @date 2015-07-26
def _h1_legendre_sum_orig_ ( h1 , N , **kwargs ) :
    """ Make a histogram representation in terms of Legendre polynomials
    >>> histo  = ...
    >>> fsum   = histo.legendre_sum_orig ( 4 )
    >>> print fsum
    >>> x = ...
    >>> print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
    It is not very CPU efficient, but stable enough...
    """
    ##
    xmin = max ( kwargs.get( 'xmin' , h1.xmin() ) , h1.xmin () ) 
    xmax = min ( kwargs.get( 'xmax' , h1.xmax() ) , h1.xmax () ) 
    ##
    return legendre_sum ( h1 , N , xmin , xmax , *kwargs )

# =============================================================================
## make a histogram representation in terms of Legendre polynomials
#  @code 
#  histo  = ...
#  fsum   = histo.legendre_sum_fill ( 4 )
#  print fsum
#  x = ...
#  print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
#  @endcode 
#  @see Ostap::Math::LegendreSum
#  @see Ostap::Math::LegendreSum::fill
#  @author Vanya Belyaev Ivan.Belyaev@iter.ru
#  @date 2015-07-26
def _h1_legendre_sum_fill_ ( h1 , N , **kwargs ) :
    """ Make a histogram representation in terms of Legendre polynomials
    >>> histo  = ...
    >>> fsum   = histo.legendre_sum_fill ( 4 )
    >>> print fsum
    >>> x = ...
    >>> print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
    - see `Ostap::Math::LegendreSum`
    - see `Ostap::Math::LegendreSum::fill`
    """
    ##
    xmin = max ( kwargs.get( 'xmin' , h1.xmin() ) , h1.xmin () ) 
    xmax = min ( kwargs.get( 'xmax' , h1.xmax() ) , h1.xmax () ) 
    ##
    if h1.nbinsx () <= N + 1 :
        logger.warning ( "Legendre-sum parameterisation: N=%d is too large for %d-bins" % ( N , h1.nbinss () ) )
        
    lsum = Ostap.Math.LegendreSum ( N , xmin , xmax ) 
    for i, x , y in h1.items() :
        xx     = x.value()
        if  not xmin <= xx <= xmax : continue                 
        yy     = y.value() 
        volume = 2 * x.error() 
        lsum.fill ( xx , yy * volume )
        
    return lsum 

# =============================================================================
## make a histogram representation in terms of Chebyshev polynomials
#  @code 
#  histo  = ...
#  fsum   = histo.chebyshev_sum_orig ( 4 )
#  print fsum
#  x = ...
#  print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
#  @endcode 
#  @see Ostap::Math::ChebyshevSum
#  @author Vanya Belyaev Ivan.Belyaev@iter.ru
#  It is not very CPU efficient, but stable enough...
#  @date 2015-07-26
def _h1_chebyshev_sum_orig_ ( h1 , N , **kwargs ) :
    """ Make a histogram representation in terms of Chebyshev polynomials
    >>> histo  = ...
    >>> fsum   = histo.chebyshev_sum_orig ( 4 )
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
## make a histogram representation in terms of Chebyshev polynomials
#  @code 
#  histo  = ...
#  fsum   = histo.chebyshev_sum_fill ( 4 )
#  print fsum
#  x = ...
#  print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
#  @endcode 
#  @see Ostap::Math::ChebyshevSum
#  @see Ostap::Math::ChebyshevSum::fill
#  @author Vanya Belyaev Ivan.Belyaev@iter.ru
#  It is not very CPU efficient, but stable enough...
#  @date 2015-07-26
def _h1_chebyshev_sum_fill_ ( h1 , N , **kwargs ) :
    """ Make a histogram representation in terms of Chebyshev polynomials
    >>> histo  = ...
    >>> fsum   = histo.chebyshev_sum_fill ( 4 )
    >>> print fsum
    >>> x = ...
    >>> print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
    """
    ##
    xmin = max ( kwargs.get( 'xmin' , h1.xmin() ) , h1.xmin () ) 
    xmax = min ( kwargs.get( 'xmax' , h1.xmax() ) , h1.xmax () ) 
    ##
    if h1.nbinsx () <= N + 1 :
        logger.warning ( "Chebyshev-sum parameterisation: N=%d is too large for %d-bins" % ( N , h1.nbinsx () ) )

    csum = Ostap.Math.ChebyshevSum ( N , xmin , xmax ) 
    for i, x , y in h1.items() :
        xx     = x.value()
        if  not xmin <= xx <= xmax : continue                 
        yy     = y.value() 
        volume = 2 * x.error() 
        csum.fill ( xx , yy * volume ) 
    return csum 

# =============================================================================
## make a histogram representation in terms of Bezier(Bernstein) sum
#  (sum over Bernstein polynomials)
#  @code 
#  histo  = ...
#  fsum   = histo.bezier_sum_orig    ( 4 )
#  fsum   = histo.bernstein_sum_orig ( 4 ) ## ditto
#  print fsum
#  x = ...
#  print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
#  @endcode 
#  @see Ostap::Math::Bernstein
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2015-07-26
def _h1_bezier_sum_orig_ ( h1 , N , **kwargs ) :
    """ Make a histogram representation in terms of Bezier/Bernstein sum
    (sum over Bernstein polynomials)
    >>> histo  = ...
    >>> fsum   = histo.bezier_sum_orig    ( 4 )
    >>> fsum   = histo.bernstein_sum_orig ( 4 ) ## ditto 
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
#  (sum over Bernstein polynomials)
#  @code 
#  histo  = ...
#  fsum   = histo.bezier_sum_fill    ( 4 )
#  fsum   = histo.bernstein_sum_fill ( 4 ) ## ditto
#  print fsum
#  x = ...
#  print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
#  @endcode 
#  @see Ostap::Math::Bernstein
#  @see Ostap::Math::Bernstein::fill
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2015-07-26
def _h1_bezier_sum_fill_ ( h1 , N , **kwargs ) :
    """ Make a histogram representation in terms of Bezier/Bernstein sum
    (sum over Bernstein polynomials)
    >>> histo  = ...
    >>> fsum   = histo.bezier_sum_fill    ( 4 )
    >>> fsum   = histo.bernstein_sum_fill ( 4 ) ## ditto 
    >>> print fsum
    >>> x = ...
    >>> print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
    - see `Ostap/Math.Bernstein`
    - see `Ostap.Math.Bernstein.fill`
    """
    ##
    xmin = max ( kwargs.get( 'xmin' , h1.xmin() ) , h1.xmin () ) 
    xmax = min ( kwargs.get( 'xmax' , h1.xmax() ) , h1.xmax () ) 
    ##
    if h1.nbinsx () <= N + 1 :
        logger.warning ( "Bernstein-sum parameterisation: N=%d is too large for %d-bins" % ( N , h1.nbinss () ) )

    bsum = Ostap.Math.Bernstein ( N , xmin , xmax ) 
    for i, x , y in h1.items() :
        xx     = x.value()
        if not xmin <= xx <= xmax : continue                 
        yy     = y.value() 
        volume = 2 * x.error() 
        bsum.fill ( xx , yy * volume ) 
    return bsum 

# =============================================================================
## make a histogram representation in terms of Bezier(Bernstein) sum
#  (sum over Bernstein polynomials) using even polynomials:
#  \f$ f( \frac{x_{min}+x_{max}}{2} - x ) = f( \frac{x_{min}+x_{max}}{2} + x ) \f$ 
#  @code 
#  histo  = ...
#  fsum   = histo.beziereven_sum    ( 2 )
#  fsum   = histo.bernsteineven_sum ( 2 ) ## ditto
#  print fsum
#  x = ...
#  print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
#  @endcode 
#  @see Ostap::Math::BernsteinEven
#  @param h1 the histogram
#  @param N  the degree of polynomial 
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2015-07-26
def _h1_beziereven_sum_ ( h1 , N , **kwargs ) :
    """ Make a histogram representation in terms of Bezier/Bernstein sum
    (sum over Bernstein polynomials) using even polynomials 
    >>> histo  = ...
    >>> fsum   = histo.beziereven_sum    ( 2 ) ##  2 is a *half-degree*
    >>> fsum   = histo.bernsteineven_sum ( 2 ) ## ditto 
    >>> print fsum
    >>> x = ...
    >>> print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
    """
    ##
    xmin = max ( kwargs.get ( 'xmin' , h1.xmin() ) , h1.xmin () ) 
    xmax = min ( kwargs.get ( 'xmax' , h1.xmax() ) , h1.xmax () ) 
    ##
    return beziereven_sum ( h1 , N , xmin , xmax )

# =============================================================================
## make a histogram representation in terms of Bezier(Bernstein) sum
#  (sum over Bernstein polynomials) using even polynomials:
#  \f$ f( \frac{x_{min}+x_{max}}{2} - x ) = f( \frac{x_{min}+x_{max}}{2} + x ) \f$ 
#  @code 
#  histo  = ...
#  fsum   = histo.beziereven_sum    ( 2 )
#  fsum   = histo.bernsteineven_sum ( 2 ) ## ditto
#  print fsum
#  x = ...
#  print 'FUN(%s) = %s ' % ( x , fsum ( x ) ) 
#  @endcode 
#  @see Ostap::Math::BernsteinEven
#  @param h1 the histogram
#  @param N  the degree of polynomial 
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2015-07-26
def _h1_rational_fun_ ( h1 , n , d  , **kwargs ) :
    """ Make a histogram representation in terms of Rational function
    inspired by Floated-Hormann's rational bary centric interpolant 
    >>> histo  = ...
    >>> rat    = rational_fun ( 2 ) ## ditto 
    >>> print ( rat )
    >>> x = ...
    >>> print 'FUN(%s) = %s ' % ( x , rat ( x ) ) 
    """
    ##
    xmin = max ( kwargs.get ( 'xmin' , h1.xmin() ) , h1.xmin () ) 
    xmax = min ( kwargs.get ( 'xmax' , h1.xmax() ) , h1.xmax () ) 
    ##
    return rational_fun  ( h1 , n , d  , xmin , xmax )

# =============================================================================
for h in ( ROOT.TH1F , ROOT.TH1D ) :

    h.beziereven_sum     = _h1_beziereven_sum_
    h.bernsteineven_sum  = _h1_beziereven_sum_
    
    h.chebyshev_sum_orig = _h1_chebyshev_sum_orig_
    h.chebyshev_sum_fill = _h1_chebyshev_sum_fill_

    h.legendre_sum_orig  = _h1_legendre_sum_orig_
    h.legendre_sum_fill  = _h1_legendre_sum_fill_


    h.bezier_sum_fill    = _h1_bezier_sum_fill_
    h.bernstein_sum_fill = _h1_bezier_sum_fill_

    h.bezier_sum_orig    = _h1_bezier_sum_orig_
    h.bernstein_sum_orig = _h1_bezier_sum_orig_

    h.rational_fun       = _h1_rational_fun_ 

    h.legendre_sum       = _h1_legendre_sum_orig_  ## ATTENTION 
    h.chebyshev_sum      = _h1_chebyshev_sum_orig_ ## ATTENTION 
    h.bezier_sum         = _h1_bezier_sum_orig_    ## ATTENTION
    h.bernstein_sum      = _h1_bezier_sum_orig_    ## ATTENTION
    
    _new_methods_ .append ( h.legendre_sum_orig  )
    _new_methods_ .append ( h.legendre_sum_fill  )
    _new_methods_ .append ( h.legendre_sum       )

    _new_methods_ .append ( h.chebyshev_sum_orig )
    _new_methods_ .append ( h.chebyshev_sum_fill )
    _new_methods_ .append ( h.chebyshev_sum      )

    _new_methods_ .append ( h.bezier_sum_fill    )
    _new_methods_ .append ( h.bernstein_sum_fill )
    _new_methods_ .append ( h.bezier_sum_orig    )
    _new_methods_ .append ( h.bernstein_sum_orig )
    _new_methods_ .append ( h.bezier_sum         )
    _new_methods_ .append ( h.bernstein_sum      )
    
    _new_methods_ .append ( h.beziereven_sum     )
    _new_methods_ .append ( h.bernsteineven_sum  )

    _new_methods_ .append ( h.rational_fun       )
    
# ==============================================================================
## (relatively) fast parameterization of 2D histogram as sum of
#  2D-Legendre polynomials
#  @code
#  histo = ...
#  func  = histo.legendre      ( 5 , 3 ) 
#  func  = histo.legendre_fast ( 5 , 3 ) ## ditto
#  func  = histo.legendre_fill ( 5 , 3 ) ## ditto
#  @endcode
#  @see Ostap::LegendreSum2
#  @see Ostap::LegendreSum2::fill 
def _h2_legendre_fast_ ( h2  , NX , NY , 
                         xmin = inf_neg , xmax = inf_pos , 
                         ymin = inf_neg , ymax = inf_pos ) :
    """ (relatively) fast parameterization of 2D histogram as sum of 2D-Legendre polynomials
    >>> histo = ...
    >>> func  = histo.legendre      ( 5 , 3  )
    >>> func  = histo.legendre_fast ( 5 , 3  )
    >>> func  = histo.legendre_fill ( 5 , 3  )
    - see Ostap.Math.LegendreSum2
    - see Ostap.Math.LegendreSum2.fill
    """
    
    assert isinstance ( NX , int ) and 0 <= NX , "Invalid NX=%s" % NX     
    assert isinstance ( NY , int ) and 0 <= NY , "Invalid NY=%s" % NY

    xmin = max ( xmin , h2.xmin() ) 
    xmax = min ( xmax , h2.xmax() )
    assert xmin < xmax , 'Invalid xmin/xmax: %s/%s' % ( xmin , xmax )

    ymin = max ( ymin , h2.ymin() ) 
    ymax = min ( ymax , h2.ymax() )
    assert ymin < ymax , 'Invalid ymin/ymax: %s/%s' % ( ymin , ymax )
    
    if h2.nbinsx() <= NX + 1 :
        logger.warning ( "Legendre-sum parameterisation: NX=%d is too large for %d-bins" % ( NX , h2.nbinsx () ) )
    if h2.nbinsy() <= NY + 1 :
        logger.warning ( "Legendre-sum parameterisation: NY=%d is too large for %d-bins" % ( NY , h2.nbinsy () ) )

    func = Ostap.Math.LegendreSum2 ( NX , NY ,  xmin , xmax , ymin , ymax )
    
    sum  = 0
    for ix , iy , x , y , v in h2.iteritems()  :
        
        xv = x.value()
        if not xmin <= xv <= xmax : continue        
        yv = y.value()
        if not ymin <= yv <= ymax : continue        
        vv = v.value()        
        w  = 4 * vv * x.error() * y.error()        
        ok = func.fill ( xv , yv , w )        
        if ok : sum +=             w

    return func 

# ==============================================================================
## (relatively) fast parameterization of 2D histogram as sum of
#  2D-Bernstein polynomials
#  @code
#  histo = ...
#  func  = histo.bernstein      ( 5 , 3 ) 
#  func  = histo.bernstein_fast ( 5 , 3 ) ## ditto
#  func  = histo.bernstein_fill ( 5 , 3 ) ## ditto
#  @endcode
#  @see Ostap::Bernstein2D
#  @see Ostap::Bernstein2D::fill 
def _h2_bernstein_fast_ ( h2  , NX , NY , 
                          xmin = inf_neg , xmax = inf_pos , 
                          ymin = inf_neg , ymax = inf_pos ) :
    """ (relatively) fast parameterization of 2D histogram as sum of 2D-Bernstein polynomials
    >>> histo = ...
    >>> func  = histo.bernstein      ( 5 , 3  )
    >>> func  = histo.bernstein_fast ( 5 , 3  )
    >>> func  = histo.bernstein_fill ( 5 , 3  )
    - see Ostap.Math.Bernstein2D
    - see Ostap.Math.Bernstein2D.fill
    """
    
    assert isinstance ( NX , int ) and 0 <= NX , "Invalid NX=%s" % NX     
    assert isinstance ( NY , int ) and 0 <= NY , "Invalid NY=%s" % NY
    
    xmin = max ( xmin , h2.xmin() ) 
    xmax = min ( xmax , h2.xmax() )
    assert xmin < xmax , 'Invalid xmin/xmax: %s/%s' % ( xmin , xmax )

    ymin = max ( ymin , h2.ymin() ) 
    ymax = min ( ymax , h2.ymax() )
    assert ymin < ymax , 'Invalid ymin/ymax: %s/%s' % ( ymin , ymax )
    
    if h2.nbinsx() <= NX + 1 :
        logger.warning ( "Bernstein-sum parameterisation: NX=%d is too large for %d-bins" % ( NX , h2.nbinsx () ) )
    if h2.nbinsy() <= NY + 1 :
        logger.warning ( "Bernstein-sum parameterisation: NY=%d is too large for %d-bins" % ( NY , h2.nbinsy () ) )

    func = Ostap.Math.Bernstein2D ( NX , NY , xmin , xmax , ymin , ymax )
    sum  = 0
    for ix , iy , x , y , v in h2.iteritems()  :        
        xv = x.value()
        if not xmin <= xv <= xmax : continue        
        yv = y.value()
        if not ymin <= yv <= ymax : continue
        vv = v.value()        
        w  = 4 * vv * x.error() * y.error()        
        ok = func.fill ( xv , yv , w )        
        if ok : sum +=             w

    return func 

# =============================================================================
for h in ( ROOT.TH2F , ROOT.TH2D ) :

    h.legendre       = _h2_legendre_fast_
    h.legendre_fast  = _h2_legendre_fast_
    h.legendre_fill  = _h2_legendre_fast_

    h.bernstein      = _h2_bernstein_fast_
    h.bernstein_fast = _h2_bernstein_fast_
    h.bernstein_fill = _h2_bernstein_fast_

    h.bezier         = _h2_bernstein_fast_
    h.bezier_fast    = _h2_bernstein_fast_
    h.bezier_fill    = _h2_bernstein_fast_

    _new_methods_ .append ( h.legendre       )
    _new_methods_ .append ( h.legendre_fast  )
    _new_methods_ .append ( h.legendre_fill  )
    
    _new_methods_ .append ( h.bernstein      )
    _new_methods_ .append ( h.bernstein_fast )
    _new_methods_ .append ( h.bernstein_fill )

    _new_methods_ .append ( h.bezier         )
    _new_methods_ .append ( h.bezier_fast    )
    _new_methods_ .append ( h.bezier_fill    )

# ==============================================================================
## (relatively) fast parameterization of 3D histogram as sum of
#  3D-Legendre polynomials
#  @code
#  histo = ...
#  func  = histo.legendre      ( 5 , 3 , 2 ) 
#  func  = histo.legendre_fast ( 5 , 3 , 2 ) ## ditto
#  func  = histo.legendre_fill ( 5 , 3 , 2 ) ## ditto
#  @endcode
#  @see Ostap::LegendreSum3
#  @see Ostap::LegendreSum3::fill 
def _h3_legendre_fast_ ( h3 , NX , NY , NZ , 
                         xmin = inf_neg , xmax = inf_pos  ,
                         ymin = inf_neg , ymax = inf_pos  ,
                         zmin = inf_neg , zmax = inf_pos  ) :
    """(relatively) fast parameterization of 2D histogram as sum of 3D-Legendre polynomials
    >>> histo = ...
    >>> func  = histo.legendre      ( 5 , 3  )
    >>> func  = histo.legendre_fast ( 5 , 3  )
    >>> func  = histo.legendre_fill  ( 5 , 3  )
    - see Ostap.Math.LegendreSum3
    - see Ostap.Math.LegendreSum3.fill
    """
    
    assert isinstance ( NX , int ) and 0 <= NX , "Invalid NX=%s" % NX     
    assert isinstance ( NY , int ) and 0 <= NY , "Invalid NY=%s" % NY
    assert isinstance ( NZ , int ) and 0 <= NZ , "Invalid NZ=%s" % NZ

    xmin = max ( xmin , h3.xmin() ) 
    xmax = min ( xmax , h3.xmax() )
    assert xmin < xmax , 'Invalid xmin/xmax: %s/%s' % ( xmin , xmax )

    ymin = max ( ymin , h3.ymin() ) 
    ymax = min ( ymax , h3.ymax() )
    assert ymin < ymax , 'Invalid ymin/ymax: %s/%s' % ( ymin , ymax )

    zmin = max ( zmin , h3.zmin() ) 
    zmax = min ( zmax , h3.zmax() )
    assert zmin < zmax , 'Invalid zmin/zmax: %s/%s' % ( zmin , zmax )
    
    if h3.nbinsx () <= NX + 1 :
        logger.warning ( "Legendre-sum parameterisation: NX=%d is too large for %d-bins" % ( NX , h3.nbinsx () ) )
    if h3.nbinsy () <= NY + 1 :
        logger.warning ( "Legendre-sum parameterisation: NY=%d is too large for %d-bins" % ( NY , h3.nbinsy () ) )
    if h3.nbinsz () <= NZ + 1 :
        logger.warning ( "Legendre-sum parameterisation: NZ=%d is too large for %d-bins" % ( NZ , h3.nbinsz () ) )

    func = Ostap.Math.LegendreSum3 ( NX , NY , NZ , 
                                     xmin , xmax ,
                                     ymin , ymax ,
                                     zmin , zmax )
    
    sum  = 0    
    for ix , iy, iz , x , y , z , v in h3.iteritems()  :
        
        xv = x.value()
        if not xmin <= xv <= xmax : continue        
        yv = y.value()
        if not ymin <= yv <= ymax : continue
        zv = z.value()
        if not zmin <= zv <= zmax : continue

        vv = v.value()
        w  = 8 * vv * x.error() * y.error() * z.error() 
        ok = func.fill ( xv , yv , zv , w )
        if ok : sum +=                  w

    return func 


# ==============================================================================
## (relatively) fast parameterization of 3D histogram as sum of
#  3D-Bernstein polynomials
#  @code
#  histo = ...
#  func  = histo.bernstein      ( 5 , 3 , 2 ) 
#  func  = histo.bernstein_fast ( 5 , 3 , 2 ) ## ditto
#  func  = histo.bernstein_fill ( 5 , 3 , 2 ) ## ditto
#  @endcode
#  @see Ostap::Bernstein3D
#  @see Ostap::Bernstein3D::fill 
def _h3_bernstein_fast_ ( h3 , NX , NY , NZ ,  
                          xmin = inf_neg , xmax = inf_pos  ,
                          ymin = inf_neg , ymax = inf_pos  ,
                          zmin = inf_neg , zmax = inf_pos  ) :
    """(relatively) fast parameterization of 2D histogram as sum of 3D-Bernstein polynomials
    >>> histo = ...
    >>> func  = histo.bernstein      ( 5 , 3  )
    >>> func  = histo.bernstein_fast ( 5 , 3  )
    >>> func  = histo.bernstein_fill ( 5 , 3  )
    - see Ostap.Math.Bernstein3D
    - see Ostap.Math.Bernstein3D.fill
    """
    
    assert isinstance ( NX , int ) and 0 <= NX , "Invalid NX=%s" % NX     
    assert isinstance ( NY , int ) and 0 <= NY , "Invalid NY=%s" % NY
    assert isinstance ( NZ , int ) and 0 <= NZ , "Invalid NZ=%s" % NZ

    xmin = max ( xmin , h3.xmin() ) 
    xmax = min ( xmax , h3.xmax() )
    assert xmin < xmax , 'Invalid xmin/xmax: %s/%s' % ( xmin , xmax )

    ymin = max ( ymin , h3.ymin() ) 
    ymax = min ( ymax , h3.ymax() )
    assert ymin < ymax , 'Invalid ymin/ymax: %s/%s' % ( ymin , ymax )

    zmin = max ( zmin , h3.zmin() ) 
    zmax = min ( zmax , h3.zmax() )
    assert zmin < zmax , 'Invalid zmin/zmax: %s/%s' % ( zmin , zmax )
    
    if h3.nbinsx () <= NX + 1 :
        logger.warning ( "Bernstein-sum parameterisation: NX=%d is too large for %d-bins" % ( NX , h3.nbinsx () ) )
    if h3.nbinsy () <= NY + 1 :
        logger.warning ( "Bernstein-sum parameterisation: NY=%d is too large for %d-bins" % ( NY , h3.nbinsy () ) )
    if h3.nbinsz () <= NZ + 1 :
        logger.warning ( "Bernstein-sum parameterisation: NZ=%d is too large for %d-bins" % ( NZ , h3.nbinsz () ) )

    func = Ostap.Math.Bernstein3D ( NX , NY , NZ , 
                                    xmin , xmax ,
                                    ymin , ymax ,
                                    zmin , zmax )
    
    sum  = 0
    
    for ix , iy, iz , x , y , z , v in h3.iteritems()  :
        
        xv = x.value()
        if not xmin <= xv <= xmax : continue        
        yv = y.value()
        if not ymin <= yv <= ymax : continue
        zv = z.value()
        if not zmin <= zv <= zmax : continue

        vv = v.value()

        w  = 8 * vv * x.error() * y.error() * z.error() 
        ok = func.fill ( xv , yv , zv , w )
        if ok : sum +=                  w

    return func 


for h in ( ROOT.TH3F , ROOT.TH3D ) :

    h.legendre       = _h3_legendre_fast_
    h.legendre_fast  = _h3_legendre_fast_
    h.legendre_fill  = _h3_legendre_fast_

    h.bernstein      = _h3_bernstein_fast_
    h.bernstein_fast = _h3_bernstein_fast_
    h.bernstein_fill = _h3_bernstein_fast_

    h.bezier         = _h3_bernstein_fast_
    h.bezier_fast    = _h3_bernstein_fast_
    h.bezier_fill    = _h3_bernstein_fast_

    _new_methods_ .append ( h.legendre       )
    _new_methods_ .append ( h.legendre_fast  )
    _new_methods_ .append ( h.legendre_fill  )
    
    _new_methods_ .append ( h.bernstein      )
    _new_methods_ .append ( h.bernstein_fast )
    _new_methods_ .append ( h.bernstein_fill )

    _new_methods_ .append ( h.bezier         )
    _new_methods_ .append ( h.bezier_fast    )
    _new_methods_ .append ( h.bezier_fill    )

_new_methods_ .append ( _h3_legendre_fast_   )

# =============================================================================
_decorated_classes = (
    ROOT.TH1D ,
    ROOT.TH1F ,
    ROOT.TH2D ,
    ROOT.TH2F ,
    ROOT.TH3D ,
    ROOT.TH3F ,
    )

_new_methods_ = tuple ( _new_methods_ )
    

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
