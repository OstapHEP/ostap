#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/funcs.py
#  Tiny decoration for ROOT.TF objects 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Tiny decoration for ROOT.TF objects 
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
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.funcs' )
else                       : logger = getLogger( __name__              )
# =============================================================================
logger.debug ( 'Tiny decoration for ROOT.TF objects')
# =============================================================================
from ostap.core.core   import cpp, Ostap, VE, funID
from ostap.core.types  import num_types , integer_types
# =============================================================================

# =============================================================================
## Generate random tuple according to TF2
#  @code
#  f2 = ...          ## ROOT TF2
#  x,y = f2.random() ## get random number 
#  @endcode
def _f2_random_ ( f2 ) :
    """Generate random tuple according to TF2
    >>> f2  = ...         ## ROOT TF2
    >>> x,y = f2.random() ## get random number 
    """
    _x = ROOT.Double(0.0)
    _y = ROOT.Double(0.0)
    f2.GetRandom2( _x , _y )
    return float(_x) , float(_y)

ROOT.TF2.random = _f2_random_ 

_new_methods_ = [
    ROOT.TF2.random ,
    ]

# =============================================================================
## fit histo
#  @see TH1::Fit
#  @code
#  func  = ...
#  histo = ...
#  func.Fit ( histo , .... )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-09-28
def _f_fit_ ( func , histo , *args ) :
    """Fit histogram (Actually delegate to TH1::Fit method)
    >>> func  = ...
    >>> histo = ...
    >>> func.Fit ( histo , .... )
    """
    return histo.Fit( func , *args )

ROOT.TF1 . Fit      = _f_fit_ 
ROOT.TF1 . fitHisto = _f_fit_ 
ROOT.TF1 . fit      = _f_fit_ 
ROOT.TH1 . fit      = ROOT.TH1.Fit

_new_methods_ += [
    ROOT.TF1 . Fit      ,
    ROOT.TF1 . fitHisto ,
    ROOT.TF1 . fit      ,
    ROOT.TH1 . fit      ,
    ]
# =============================================================================
## check existence parameter for the function
#  @code 
#  fun = ...         ## function
#    >>> if i   in fun : ... ## check if i   is valid parameter number 
#    >>> if 'm' in fun : ... ## check if 'm' is valid parameter name
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-13
def _tf1_contains_ ( func , par ) :
    """Check existence parameter for the function
    >>> fun = ...         ## function
    >>> if i   in fun : ... ## check if i   is valid parameter number 
    >>> if 'm' in fun : ... ## check if 'm' is valid parameter name  
    """
    ## check name 
    if   isinstance ( par , str           ) : return 0<=func.GetParNumber ( par ) 
    elif isinstance ( par , integer_types ) : return 0<= par<func.GetNpar (     )
    #
    return False 

# =============================================================================
## Fix parameter for TF1
#  @code
#  fun =  ...     ## function
#  fun.fix(1,0)   ## fix parameter #1  at 0 
#  fun.fix(2)     ## fix parameter #2  at current value
#  fun.fix('m',1) ## fix parameter 'm' at 1 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-13
def _tf1_fix_ ( func , par , value = None ) :
    """Fix parameter for TF1
    >>> fun =  ...   ## function 
    >>> fun.fix(1,0) ## fix parameter #1 at 0 
    >>> fun.fix(2)   ## fix parameter #2 at current value 
    """
    if not par in func : raise IndexError("Invalid parameter index %s" % par )
    if     isinstance ( par , str  ) : par = func.GetParNumber( par )
    ##
    if not isinstance ( value , num_types )  :
        value = func.GetParameter(par)
    #
    func.FixParameter( par , value ) 

# =============================================================================
## Release parameter for TF1
#  @code
#  fun =  ...       ## function
#  fun.release(1)   ## release parameter #1 
#  fun.release('m') ## release parameter 'm'
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _tf1_release_ ( func , par ) :
    """Release parameter for TF1
    >>> fun =  ...       ## function
    >>> fun.release(1)   ## release parameter #1 
    >>> fun.release('m') ## release parameter 'm'
    """
    #
    if not par in func : raise IndexError("Invalid parameter index %s" % par )
    #
    if     isinstance ( par , str  ) : par = func.GetParNumber( par )
    func.ReleaseParameter( par ) 


# =============================================================================
## get the parameter from TF1
#  @code
#  fun =  ...   ## function
#  p = fun.par(1) 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _tf1_par_ ( func , par ) :
    """Get parameter from TF1
    >>> fun =  ...        ## function 
    >>> p2 = fun.par(2)   ## get parameter #2 
    >>> pm = fun.par('m') ## get parameter 'm'
    """
    if not par in func : raise IndexError("Invalid parameter index %s" % par )
    #
    if isinstance ( par , str  ) : par = func.GetParNumber( par )    
    v = func.GetParameter ( par )
    e = func.GetParError  ( par )
    #
    return VE ( v , e * e )


# =============================================================================
## set parameter of TF1
#  @code
#  fun =  ...   ## function
#  fun.setPar(1,1) 
#  fun.setPar('m',2) 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _tf1_setpar_ ( func , par , value ) :
    """Set parameter of TF1 
    >>> fun =  ...          ## function 
    >>> fun.setPar(1,1)     ## set parameter #1 to be 1 
    >>> fun.setPar('m',2)   ## set parameter 'm' to be 2
    """
    if not par in func : raise IndexError("Invalid parameter index %s" % par )
    #
    if isinstance ( par , str  ) : par = func.GetParNumber( par )
    #
    func.SetParameter ( par , float ( value ) )

# =============================================================================
## primitive iteration over parameters:
#  @code
#  fun =  ...        ## function
#  for p in fun: print fun(p)
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _tf1_iter_ ( func ) :
    """Primitive iteration over parameters 
    >>> fun =  ...        ## function 
    >>> for p in fun: print fun(p)
    """
    s = func.GetNpar()
    i = 0
    while i < s :
        yield i
        i += 1
        
# =============================================================================
## get parameter as attribute
#  @code
#  fun =  ...   ## function
#  pm  = fun.m  ## get parameter 'm'
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _tf1_getattr_ ( func , par ) :
    """Get parameter as attribute
    >>> fun =  ...   ## function
    >>> pm  = fun.m  ## get parameter 'm'
    """
    if not isinstance ( par , str ) : raise AttributeError('TF1:Invalid attribute %s' % par )
    if not par in func              : raise AttributeError('TF1:Invalid attribute %s' % par )
    return _tf1_par_  ( func , par )

    
ROOT.TF1.__contains__ = _tf1_contains_
ROOT.TF1.__len__      = lambda s : s.GetNpar() 
    
ROOT.TF1.par          = _tf1_par_
ROOT.TF1.param        = _tf1_par_
ROOT.TF1.parameter    = _tf1_par_

ROOT.TF1.setPar       = _tf1_setpar_
ROOT.TF1.__setitem__  = _tf1_setpar_

ROOT.TF1.fix          = _tf1_fix_
ROOT.TF1.rel          = _tf1_release_
ROOT.TF1.release      = _tf1_release_

ROOT.TF1.__iter__     = _tf1_iter_
ROOT.TF1.__getitem__  = _tf1_par_
ROOT.TF1.__getattr__  = _tf1_getattr_



ROOT.TF1.integral     = ROOT.TF1.Integral 
ROOT.TF2.integral     = ROOT.TF2.Integral 


ROOT.TF2.xminmax = lambda s : ( s.GetXmin() , s.GetXmax() )
ROOT.TF2.yminmax = lambda s : ( s.GetYmin() , s.GetYmax() )


_new_methods_ += [
    ROOT.TF1.__contains__ ,
    ROOT.TF1.__len__      ,
    #
    ROOT.TF1.par          ,
    ROOT.TF1.param        ,
    ROOT.TF1.parameter    ,
    #
    ROOT.TF1.setPar       ,
    ROOT.TF1.__setitem__  ,
    #
    ROOT.TF1.fix          ,
    ROOT.TF1.rel          ,
    ROOT.TF1.release      ,
    #
    ROOT.TF1.__iter__     ,
    ROOT.TF1.__getitem__  ,
    ROOT.TF1.__getattr__  ,
    ]
# =============================================================================
## Integrate TF2 over range in Y
#  \f$ f = \int_{y_{low}}^{y_{high}} f(x,y) dy \f$
#  @code
#  func = ROOT.TF2( ... )
#  a = func.integrate_Y ( x , ylow , yhigh ) 
#  @endcode
#  @see ostap.math.integral
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2016-05-03
def _tf2_integrate_Y_ ( tf2 , x , ylow = None , yhigh = None ) :
    r"""Integrate TF2 over range in Y
    f = \int_{y_{low}}^{y_{high}} f(x,y) dy
    
    >>> func = ROOT.TF2( ... )
    >>> a = func.integrate_Y ( x , ylow , yhigh )
    """
    ## check X 
    xmin,xmax = tf2.xminmax()
    if not xmin <= x <= xmax : return 0
    ## check Y
    ymin,ymax = tf2.yminmax() 
    if None is ylow  : ylow  = ymin
    if None is yhigh : yhigh = ymax
    ymin = max ( ymin , ylow  )
    ymax = min ( ymax , yhigh )
    ##
    func_y = lambda y : tf2 ( x , y )
    from ostap.math.integral import integral 
    return integral ( func_y , ymin , ymax )

# =============================================================================
## Integrate TF2 over range in Y
#  \f$ f = \int_{x_min}^{x_max} f(x,y) dx \f$
#  @code
#  func = ROOT.TF2( ... )
#  a = func.integrate_X ( y , xlow , xhigh ) 
#  @endcode
#  @see ostap.math.integral
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2016-05-03
def _tf2_integrate_X_ ( tf2 , y , xlow = None , xhigh = None ) :
    r"""Integrate TF2 over range in X
    f = \int_{x_{low}}^{x_{high}} f(x,y) dx 
    
    >>> func = ROOT.TF2( ... )
    >>> a = func.integrate_X ( y , xlow , xhigh )
    """
    ## check Y 
    ymin,ymax = tf2.yminmax()
    if not ymin <= y <= ymax : return 0
    ## check X
    xmin,xmax = tf2.yminmax() 
    if None is xlow  : xlow  = xmin
    if None is xhigh : xhigh = xmax
    xmin = max ( xmin , xlow  )
    xmax = min ( xmax , xhigh )
    ##
    func_x = lambda x : tf2 ( x , y )
    from ostap.math.integral import integral 
    return integral ( func_x , xmin , xmax )

    
ROOT.TF2.integrate_X = _tf2_integrate_X_
ROOT.TF2.integrate_Y = _tf2_integrate_Y_

_new_methods_ += [
    ROOT.TF2.integrate_X ,
    ROOT.TF2.integrate_Y ,
    ]

# =============================================================================
_decorated_classes_ = (
    ROOT.TF1 ,
    ROOT.TF2 ,
    ROOT.TH1 ,
    )
_new_methods_ = tuple ( _new_methods_ )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
