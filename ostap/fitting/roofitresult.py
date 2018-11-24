#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/roofitresult.py
#  Module with decoration of some RooFit objects for efficient use in python
#  @see RooFitResult
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Decoration of some RooFit objects for efficient use in python
- see RooFitResult
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    ) 
# =============================================================================
import ROOT
from   ostap.core.core          import Ostap, VE, valid_pointer 
import ostap.fitting.variables     
import ostap.fitting.printable     
# =============================================================================
_new_methods_ = []

# =============================================================================
## ``easy'' print of RooFitResult
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _rfr_print_ ( self , opts = 'v' ) :
    """Easy print of RooFitResult
    >>> result = ...
    >>> print result    
    """
    if not valid_pointer ( self ) : return 'Invalid RooFitResult'
    return self.print_multiline ( content = 1 , verbose = True )

# =============================================================================
## get parameters from RooFitResult
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _rfr_params_ (self , float_only = True ) :
    """GetParameters from RooFitResult:
    >>> result = ...
    >>> params = results
    >>> p0     = params()['A'][0]  ## get the value
    >>> p0s    = params()['A'][1]  ## get the parameter itself     
    """
    pars  = self.floatParsFinal()
    pars_ = {}
    for p in pars :
        pars_ [ p.GetName() ] = p.as_VE(), p

    ## also fixed parameters? 
    if not float_only :
        fixed = self.constPars()
        for p in fixed :
            pars_ [ p.GetName() ] = p.as_VE(), p
            
    return pars_

# =============================================================================
## get parameter by name  from RooFitResult
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _rfr_param_  ( self , pname , float_only = False ) :
    """Get Parameter from RooFitResult by name 
    >>> result = ...
    >>> signal = results.param('Signal')
    >>> print signal
    """
    if not isinstance ( pname , str ) :
        if   hasattr ( pname , 'GetName' ) : pname = pname.GetName ()
        elif hasattr ( pname , 'getName' ) : pname = pname.getName ()
        elif hasattr ( pname , 'name'    ) : pname = pname.   name () 
    p = self.parameters ( float_only )[ pname ] 
    return p 

# =============================================================================
## iterator over fit results 
def _rfr_iter_ ( self ) :
    """Iterator over fit results :
    >>> fit_result = ...
    >>> for i in fit_results : print i 
    """
    pars  = self.floatParsFinal()
    for p in pars  : yield p
    fixed = self.constPars     ()
    for f in fixed : yield f

# =============================================================================
## iterator over fit items  
def _rfr_iteritems_ ( self , float_only = False ) :
    """Iterator over fit items:
    >>> fit_result = ...
    >>> for name,var in fit_results.iteritems() :
    ...                   print name,var.as_VE()  
    """
    pars  = self.floatParsFinal()
    for p in pars  :
        yield p.GetName() , p
        
    if not float_only :  
        fixed = self.constPars ()
        for f in fixed :
            yield f.GetName() , f

# =============================================================================
## get the correlation coefficient
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _rfr_corr_  ( self , var1 , var2 ) :
    """Get correlation coefficient for two parameters 
    >>> result = ...
    >>> corr = results.corr('Signal', 'Background')
    >>> print corr
    """
    if isinstance ( var1 ,  str ) : var1 = self.param ( var1 )[1]
    if isinstance ( var2 ,  str ) : var2 = self.param ( var2 )[1]
    #
    if var1 in self.constPars() : return 0.0
    if var2 in self.constPars() : return 0.0
    #
    return self.correlation ( var1 , var2 ) 

# =============================================================================
## get the covariance (sub) matrix 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _rfr_cov_matrix_  ( self , var1 , var2 , *vars ) :
    """Get covariance (sub) matrix 
    >>> result = ...
    >>> cov = results.cov_matrix('Signal', 'Background')
    >>> print corr
    """
    if isinstance ( var1 , str ) : var1 = self.param (   var1 ) [1] 
    if isinstance ( var2 , str ) : var2 = self.param (   var2 ) [1]
    
    args = ROOT.RooArgList ( var1 , var2 )
    for v in vars :
        if isinstance ( v , str ) : v = self.param ( v ) [1] 
        args.add ( v ) 
        
    cm = self.reducedCovarianceMatrix (  args )
    N  = cm.GetNrows()

    import ostap.math.linalg 
    m  = Ostap.Math.SymMatrix ( N )()

    for i in range ( N ) :
        for j in  range ( i , N ) :
            m [i,j] = cm(i,j)
            
    return m  

# =============================================================================
## get the covariance (sub) matrix 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _rfr_cov_  ( self , name1 , name2 ) :
    """Get covariance (sub) matrix 
    >>> result = ...
    >>> cov = results.cov('Signal', 'Background')
    >>> print corr
    """
    if isinstance ( var1 ,  str ) : var1 = self.param ( var1 )[1]
    if isinstance ( var2 ,  str ) : var2 = self.param ( var2 )[1]
    #
    if var1 in self.constPars() : return 0.0
    if var2 in self.constPars() : return 0.0
    #
    r  = self.correlation ( var1 , var2 )
    #
    v1 = var1.error
    v2 = var2.error
    # 
    return v1 * v2 * r 

# ===============================================================================
## get fit-parameter as attribute
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-08-03
def _rfr_getattr_ ( self , att ) :
    """Get fit-parameter as attribute
    >>> r = ....
    >>> print r.sigma 
    """
    ##
    pars = self.floatParsFinal()
    for p in pars :
        if att == p.GetName() : return p      
    #
    pars = self.constPars()
    for p in pars :
        if att == p.GetName() : return p
        
    raise AttributeError ( 'RooFitResult: invalid attribute %s ' % att )

# ===========================================================================
## get correct estimate of sum of two (or more) variables,
#  taking into account correlations
#  @code
#  >>> r = ...
#  >>> print r.sum( 'S' , 'B' )  ## S+B
#  @endcode
def _rfr_sum_ ( self , var1 , var2 , *vars ) :
    """Get correct estimate of sum of two or more variables,
    taking into account correlations
    >>> r = ...
    >>> print r.sum( 'S' , 'B' ) ## S+B
    """
    allvars = ( var1 , var2 ) + vars 
    n       = len ( allvars ) 
    s  = 0
    c2 = 0
    for i in range ( n ) :
        vi = allvars [ i ]
        if isinstance ( vi , str ) : vi = self.param ( vi ) [1]        
        v   = vi.value
        v   = VE ( v ) 
        s  += v . value ()
        vc  = v.cov2() 
        if 0 >= vc or vi in self.constPars() : continue        
        c2 += vc 
        for j in range ( i + 1 , n ) :
            vj  = allvars [ j ]
            if isinstance ( vj , str ) : vj = self.param ( vj ) [1]
            if vj in self.constPars()  : continue        
            c2 += 2 * self.correlation ( vi , vj ) 
            
    return VE ( s , c2 ) 
 
# ===========================================================================
## get correct estimate of product of two variables,
#  taking into account correlations
#  @code
#  >>> r = ...
#  >>> print r.multiply( 'S' , 'B' ) ## S*B
#  @endcode
#  @see Gaudi:Math::multiply 
def _rfr_multiply_ ( self , var1 ,  var2 , *vars ) :
    """Get correct estimate of product of two variables,
    taking into account correlations
    >>> r = ...
    >>> print r.multiply( 'S' , 'B' ) ## S*B
    """
    allvars = ( var1 , var2 ) + vars 
    n       = len ( allvars ) 
    
    m  = 1.0
    c2 = 0
    for i in range ( n ) :
        vi = allvars[i]
        if isinstance ( vi , str ) : vi = self.param ( vi )[1]
        v  = vi.value
        v  = VE ( v ) 
        vv = v.value ()
        if iszero ( vv ) or iszero ( m ) : return  VE ( 0.0 , 0.0 )   ## RETURN HERE
        m  *= vv
        vc  = v.cov2()
        if 0 >= vc or vi in self.constPars() : continue        
        c2 += vc / ( vv * vv )        
        for j in range ( i + 1 , n ) :            
            vj  = allvars [ j ]
            if isinstance ( vj , str ) : vj = self.param( vj )[1]
            if vj in self.constPars()  : continue        
            w   = vj . value
            w   = VE ( w ) 
            ww  = w.value() 
            c2 += 2 * self.correlation ( vi , vj ) / ( vv * ww ) 
            
    return  VE ( m , c2 * m * m ) 
    
# ===========================================================================
## get correct estimate of division  of two variables,
#  taking into account correlations
#  @code
#  >>> r = ...
#  >>> print r.divide( 'S' , 'B' ) ## S/B
#  @endcode
#  @see Ostap:Math::divide
def _rfr_divide_ ( self , var1 , var2 ) :
    """Get correct estimate of division of two variables,
    taking into account correlations
    >>> r = ...
    >>> print r.divide( 'S' , 'B' ) ## S/B
    """
    if isinstance ( var1 , str ) : var1 = self.param ( var1 ) [1]
    if isinstance ( var2 , str ) : var2 = self.param ( var2 ) [1]    
    _v1  = var1.value
    _v2  = var2.value 
    _cor = self.corr  ( var1 , var2 ) 
    return Ostap.Math.divide ( _v1 , _v2 , _cor ) 

# ===========================================================================
## get correct estimate of subtraction of two variables,
#  taking into account correlations
#  @code
#  >>> r = ...
#  >>> print r.subtract( 'S' , 'B' ) ## S-B
#  @endcode
#  @see Ostap:Math::subtract
def _rfr_subtract_ ( self , var1 , var2 ) :
    """Get correct estimate of subtraction of two variables,
    taking into account correlations
    >>> r = ...
    >>> print r.subtract( 'S' , 'B' ) ## S-B
    """
    if isinstance ( var1 , str ) : var1 = self.param ( var1 ) [1]
    if isinstance ( var2 , str ) : var2 = self.param ( var2 ) [1]    
    _v1  = var1.value 
    _v2  = var2.value 
    _cor = self.corr  ( var1 , var2 ) 
    return Ostap.Math.subtract ( _v1 , _v2 , _cor ) 

# ===========================================================================
## get correct estimate of fraction  of two variables,
#  taking into account correlations
#  @code
#  >>> r = ...
#  >>> print r.fraction( 'S' , 'B' ) ## S/(S+B)
#  @endcode
#  @see Gaudi:Math::divide
def _rfr_fraction_ ( self , var1 , var2 ) :
    """Get correct estimate of fraction of two variables,
    taking into account correlations
    >>> r = ...
    >>> print r.fraction( 'S' , 'B' ) ##   S/(S+B)
    """
    if isinstance ( var1 , str ) : var1 = self.param ( var1 ) [1]
    if isinstance ( var2 , str ) : var2 = self.param ( var2 ) [1]    
    _av1  = abs ( var1.value.value() ) 
    _av2  = abs ( var2.value.value() ) 
    if _av1 > _av2 : return 1 / ( 1 + self.ratio ( var2 , var1 ) )
    return 1.0 - self.fraction ( var2 , var1 ) 

# ============================================================================
## get the required results in form of SVectorWithError object
#  @code
#  fit_resuts = ...
#  res   = fit_results.results( 'A', 'B' , 'C' )
#  print res, res.cov2() 
#  @endcode
#  @see Ostap::Math::SVectorWithError
def _rfr_results_( self , *vars ) :
    """Get the required results in form of SVectorWithError object
    >>> fit_resuts = ...
    >>> res   = fit_results.results( 'A', 'B' , 'C' )
    >>> print res, res.cov2() 
    """
    _n = len ( vars )
    _r = Ostap.Math.SVectorWithError(_n,'double')()
    _i = 0 
    for _i1 in range( 0 , _n ) :
        _v1                = vars[_i1]
        _vv                = self.param ( _v1 ) [0]
        _r       [ _i1   ] = _vv
        _r.cov2()[_i1,_i1] = _vv.cov2() 
        for _i2 in range ( _i1 + 1 , _n ) :
            _v2  = vars[_i2]
            _c12 = self.cov ( _v1 , _v2 )(0,1) 
            _r.cov2()[_i1,_i2] = _c12 
    return _r 
        
# =============================================================================
## some decoration over RooFitResult
ROOT.RooFitResult . __repr__    = _rfr_print_
ROOT.RooFitResult . __str__     = _rfr_print_
ROOT.RooFitResult . __call__    = _rfr_param_
ROOT.RooFitResult . __getattr__ = _rfr_getattr_ 
ROOT.RooFitResult . __iter__    = _rfr_iter_
ROOT.RooFitResult . iteritems   = _rfr_iteritems_
ROOT.RooFitResult . parameters  = _rfr_params_
ROOT.RooFitResult . params      = _rfr_params_
ROOT.RooFitResult . param       = _rfr_param_
ROOT.RooFitResult . parameter   = _rfr_param_
ROOT.RooFitResult . corr        = _rfr_corr_
ROOT.RooFitResult . cor         = _rfr_corr_
ROOT.RooFitResult . cov         = _rfr_cov_
ROOT.RooFitResult . covariance  = _rfr_cov_
ROOT.RooFitResult . cov_matrix  = _rfr_cov_matrix_
ROOT.RooFitResult . parValue    = lambda s,n : s.parameter(n)[0]
ROOT.RooFitResult . sum         = _rfr_sum_
ROOT.RooFitResult . plus        = _rfr_sum_
ROOT.RooFitResult . multiply    = _rfr_multiply_
ROOT.RooFitResult . product     = _rfr_multiply_
ROOT.RooFitResult . subtract    = _rfr_subtract_
ROOT.RooFitResult . minus       = _rfr_subtract_
ROOT.RooFitResult . divide      = _rfr_divide_
ROOT.RooFitResult . ratio       = _rfr_divide_
ROOT.RooFitResult . fraction    = _rfr_fraction_
ROOT.RooFitResult . results     = _rfr_results_

_new_methods_ += [
    ROOT.RooFitResult . __repr__    ,
    ROOT.RooFitResult . __str__     ,
    ROOT.RooFitResult . __call__    ,
    ROOT.RooFitResult . __getattr__ ,
    ROOT.RooFitResult . __iter__    ,
    ROOT.RooFitResult . iteritems   ,
    ROOT.RooFitResult . parameters  ,
    ROOT.RooFitResult . params      ,
    ROOT.RooFitResult . param       ,
    ROOT.RooFitResult . parameter   ,
    ROOT.RooFitResult . corr        ,
    ROOT.RooFitResult . cor         ,
    ROOT.RooFitResult . cov         ,
    ROOT.RooFitResult . covariance  ,
    ROOT.RooFitResult . parValue    ,
    ROOT.RooFitResult . sum         ,
    ROOT.RooFitResult . plus        ,
    ROOT.RooFitResult . multiply    ,
    ROOT.RooFitResult . product     ,
    ROOT.RooFitResult . subtract    ,
    ROOT.RooFitResult . minus       ,
    ROOT.RooFitResult . divide      ,
    ROOT.RooFitResult . ratio       ,
    ROOT.RooFitResult . fraction    ,
    ROOT.RooFitResult . results     ,
    ]

# =============================================================================
_decorated_classes_ = (
    ROOT.RooFitResult , 
    )

_new_methods_ = tuple ( _new_methods_ ) 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
