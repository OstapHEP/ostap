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
from   ostap.core.core          import Ostap, VE, valid_pointer, iszero, isequal
from   ostap.core.ostap_types   import string_types , integer_types
from   ostap.utils.valerrors    import ValWithErrors, AsymErrors   
from   ostap.logger.colorized   import allright, attention, attstr 
from   ostap.logger.pretty      import pretty_float, pretty_error2
from   ostap.logger.symbols     import show 
import ostap.math.linalg        as     LA 
import ostap.fitting.variables     
import ostap.fitting.printable
import ROOT, math, sys, ctypes  
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.roofitresult' )
else                       : logger = getLogger ( __name__                     )
# =============================================================================        
_new_methods_ = []

# =============================================================================
## get parameters from RooFitResult
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _rfr_params_ (self , float_only = True ) :
    """ GetParameters from RooFitResult:
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

# ==============================================================================
## get numerical values of all parameters as simple dictionary
#  @code
#  result = ...
#  dct = results.dct_params() 
#  @endcode
def _rfr_dct_params_ ( self ) :
    """ Get numerical values of all parameters as simple dictionary
    >>> result = ...
    >>> dct = results.dct_params() 
    """
    pars_ = {}

    for p in self.floatParsFinal(): 
        pars_ [ p.GetName () ] = p.as_VE()
        
    for p in self.constPars() : 
        pars_ [ p.GetName () ] = p.getValue() 

    return pars_ 
    
# =============================================================================
## get parameter by name  from RooFitResult
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _rfr_param_  ( self , pname , float_only = False ) :
    """ Get Parameter from RooFitResult by name 
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
## check if certain parameter (or index) is in <code>RooFitResult</code> object
#  @code
#  fit_result = ...
#  print ( 1      in fit_results ) 
#  print ( 'mean' in fit_results ) 
#  print ( 'A'    in fit_results ) 
#  @endcode 
def _rfr_contains_ ( self , label ) :
    """ Check if certain parameter (or index) is in <code>RooFitResult</code> object
    >>> fit_result = ...
    >>> print ( 1      in fit_results ) 
    >>> print ( 'mean' in fit_results ) 
    >>> print ( 'A'    in fit_results ) 
    """
    
    if   isinstance ( label , integer_types ) :
        return 0 <= label < len ( self.floatParsFinal() ) + len ( self.constPars () )
    elif isinstance ( label , string_types  ) :
        return label in  self.floatParsFinal() or label in self.constPars ()
    elif isinstance ( label , ROOT.RooAbsArg ) :
        return label.GetName() in self
    
    return False 

# =============================================================================
## iterator over fit results 
def _rfr_iter_ ( self ) :
    """ Iterator over fit results :
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
    """ Iterator over fit items:
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
    """ Get correlation coefficient for two parameters 
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
    """ Get covariance (sub) matrix 
    >>> result = ...
    >>> cov = results.cov_matrix('Signal', 'Background')
    >>> print ( cov ) 
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
## get the covariance matrix 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _rfr_covmatrix_  ( self , tmatrix = False ) :
    """ Get covariance ( matrix 
    >>> result = ...
    >>> cov = results.covmatrix()[
    >>> print ( corr ) 
    """
        
    cm = self.covarianceMatrix ()
    
    if tmatrix : return cm ## RETURN NATIVE ROOT TMAtrixSym<double

    N  = cm.GetNrows()
    
    import ostap.math.linalg 
    m  = Ostap.Math.SymMatrix ( N )()

    for i in range ( N ) :
        for j in  range ( i , N ) :
            m [i,j] = cm(i,j)
            
    return m  

# ============================================================================
## get the required results in form of SVectorWithError object
#  @code
#  fit_resuts = ...
#  res   = fit_results.results( 'A', 'B' , 'C' )
#  print res, res.cov2() 
#  @endcode
#  @see Ostap::Math::SVectorWithError
def _rfr_results_( self , *vars ) :
    """ Get the required results in form of SVectorWithError object
    >>> fit_resuts = ...
    >>> res   = fit_results.results( 'A', 'B' , 'C' )
    >>> print res, res.cov2() 
    """

    if not vars :
        pars = self.floatParsFinal()
        vars = [ p.name  for p in pars ]
        logger.debug ( 'Use parameters %s' % vars )

    if not vars :
        return None
    

    elif 1 == len ( vars ) :
        return self.param ( vars [0]) [0] 
        
    args = ROOT.RooArgList ()
    for v in vars :
        if isinstance ( v , str ) : v = self.param ( v ) [1] 
        args.add ( v ) 
        
    cm = self.reducedCovarianceMatrix (  args )
    N  = cm.GetNrows()

    import ostap.math.linalg    
    m  = Ostap.Math.SymMatrix ( N )()
    v  = Ostap.Math.Vector    ( N )()
    
    S  = Ostap.Math.SVectorWithError  ( N )
    
    v  = S()
    c2 = v.cov2()
    
    for i , a  in enumerate ( args ) :
        v [ i ] = float ( a )
        for j in range ( i , N ) :
            c2 [ i , j ] = cm ( i , j )

    return v

    
# ==============================================================================
## Get vector of eigenvalues for the covariance matrix
#  @code
#  result = ...
#  eigenvalues =  result.cov_eigenvalues() 
#  @endcode 
def _rfr_eigenvalues_ ( self , sorted = True ) :
    """ Get vector of eigenvalues for the covariance matrix
    >>> result = ...
    >>> eigenvalues =  result.cov_eigenvalues() 
    """
    cm = _rfr_covmatrix_ ( self )
    return Ostap.Math.EigenSystems.eigenValues ( cm , sorted )

# =============================================================================
## get the covariance (sub) matrix 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _rfr_cov_  ( self , var1 , var2 ) :
    """ Get covariance (sub) matrix 
    >>> result = ...
    >>> cov = results.cov('Signal', 'Background')
    >>> print (cov)
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
## find a maximal correlation coefficient for the given variable
#  @code
#  result = ...
#  coefficient , variable = result.max_cor ( 'x' ) 
#  @endcode
def _rfr_max_cor_ ( self , v ) :
    """ Find a maximal correlation coefficient for the given variable
    >>> result = ...
    >>> coefficient , variable = result.max_cor ( 'x' ) 
    """
    
    if isinstance ( v ,  str ) : v = self.param ( v )[1]
    if v in self.constPars() : return 0 , ''
    
    pars = self.floatParsFinal()
    assert v in pars, 'Unknown variable %s' % v

    if 1 == len ( pars ) : return 0 , '' 
    
    rmax = None
    pmax = None 
    for p in pars :
        if v is p : continue
        r =  self.correlation ( v , p )
        assert -1 <= r <= 1 or isequal ( abs ( r ) , 1 ) ,\
               'Invalid correlation coefficient for (%s,%s):%s' % ( v.name , p.name , r )
        
        if ( rmax is None ) or ( pmax is None ) or abs ( r ) > abs ( rmax ) :
            rmax = r
            pmax = p.name
            
    return rmax , pmax 

    
# ===============================================================================
## get fit-parameter as attribute
#  @code
#  fit_result = ...
#  sigma = fit_resul.sigma 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-08-03
def _rfr_getattr_ ( self , att ) :
    """ Get fit-parameter as attribute
    >>> r = ....
    >>> print( 'sigma is', r.sigma)
    """
    ##
    pars = self.floatParsFinal()
    for p in pars :
        if att == p.GetName() : return p      
    #
    pars = self.constPars()
    for p in pars :
        if att == p.GetName() : return p
        
    raise AttributeError ( "RooFitResult: invalid attribute `%s'" % att )

# =============================================================================
## get fit-parameter through the key/name 
#  @code
#  fit_result = ...
#  sigma = fit_result['sigma']
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-08-03
def _rfr_getitem_ ( self , key  ) :
    """ Get fit-parameter through the key/name 
    >>> fit_result = ...
    >>> sigma = fit_result['sigma']
    """
    
    ##
    pars = self.floatParsFinal()
    for p in pars :
        if key == p.GetName() : return p      
    #
    pars = self.constPars()
    for p in pars :
        if key == p.GetName() : return p
        
    raise KeyError ( 'RooFitResult: invalid key %s ' % key  )
    
# ===========================================================================
## get correct estimate of sum of two (or more) variables,
#  taking into account correlations
#  @code
#  >>> r = ...
#  >>> print r.sum( 'S' , 'B' )  ## S+B
#  @endcode
def _rfr_sum_ ( self , var1 , var2 , *vars ) :
    """ Get correct estimate of sum of two or more variables,
    taking into account correlations
    >>> r = ...
    >>> print r.sum( 'S' , 'B' ) ## S+B
    """
    allvars = ( var1 , var2 ) + vars 
    n       = len ( allvars ) 
    s  = 0.0
    c2 = 0.0
    for i in range ( n ) :
        vi = allvars [ i ]
        if isinstance ( vi , str ) : vi = self.param ( vi ) [1]        
        v   = vi.value
        v   = VE ( v ) 
        s  += v . value ()
        vc  = v.cov2() 
        if 0 >= vc or vi in self.constPars() : continue        
        c2 += vc 
        for j in range ( i ) :
            vj  = allvars [ j ]
            if isinstance ( vj , str ) : vj = self.param ( vj ) [1]
            if vj in self.constPars() : continue        
            c2 += 2 * self.covariance ( vi , vj ) 
            
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
    """ Get correct estimate of product of two variables,
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
        for j in range ( i  ) :            
            vj  = allvars [ j ]
            if isinstance ( vj , str ) : vj = self.param( vj )[1]
            if vj in self.constPars()  : continue        
            w   = vj . value
            w   = VE ( w ) 
            ww  = w.value() 
            c2 += 2 * self.covariance ( vi , vj ) / ( vv * ww ) 
            
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
    """ Get correct estimate of division of two variables,
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
    """ Get correct estimate of subtraction of two variables,
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
    """ Get correct estimate of fraction of two variables,
    taking into account correlations
    >>> r = ...
    >>> print r.fraction( 'S' , 'B' ) ##   S/(S+B)
    """
    if isinstance ( var1 , str ) : var1 = self.param ( var1 ) [1]
    if isinstance ( var2 , str ) : var2 = self.param ( var2 ) [1]    
    _av1  = abs ( var1.value.value() ) 
    _av2  = abs ( var2.value.value() ) 
    if _av1 >= _av2 : return 1 / ( 1 + self.ratio ( var2 , var1 ) )
    return 1.0 - self.fraction ( var2 , var1 ) 

# ===========================================================================
## get correct estimate of asymmetry of two variables (v1-v2)/(v1+v2)
#  taking into account correlations
#  @code
#  >>> r = ...
#  >>> print r.asymmetry( 'A' , 'B' ) ## (A-B)/(A+B)
#  @endcode
#  @see Gaudi:Math::divide
def _rfr_asymmetry_ ( self , var1 , var2 ) :
    """ Get correct estimate of asymmetry of two variables, (v1-v2)/(v1+v2)
    taking into account correlations
    >>> r = ...
    >>> print r.asymmetry( 'A' , 'B' ) ##   (A-B)/(A+B)
    """
    if isinstance ( var1 , str ) : var1 = self.param ( var1 ) [1]
    if isinstance ( var2 , str ) : var2 = self.param ( var2 ) [1]
    ##
    _av1 = abs ( var1.value.value() ) 
    _av2 = abs ( var2.value.value() )
    _one = VE  ( 1 , 0 )

    if _av1 <= _av2 :
        return      self.ratio ( var1 , var2 ) . asym ( _one )
    else : 
        return -1 * self.ratio ( var2 , var1 ) . asym ( _one )

# ===========================================================================
## get correct estimate of power:  \f$ a^b \f$ 
#  taking into account correlations
#  @code
#  >>> r = ...
#  >>> print r.pow( 'A' , 'B' ) ## A**B
#  @endcode
#  @see Gaudi:Math::pow
def _rfr_pow_ ( self , var1 , var2 ) :
    """ Get correct estimate of pow-function of two variables, v1**v2, 
    taking into account correlations
    >>> r = ...
    >>> print r.pow( 'A' , 'B' ) ##  A**B 
    """
    if isinstance ( var1 , str ) : var1 = self.param ( var1 ) [1]
    if isinstance ( var2 , str ) : var2 = self.param ( var2 ) [1]
    ##
    _v1  = VE ( var1.value ) 
    _v2  = VE ( var2.value ) 
    _cor = self.corr  ( var1 , var2 )
    return Ostap.Math.pow ( _v1 , _v2 , _cor ) 

    
    
# =============================================================================
## evaluate the certain  function/expression for the fit   result
# @code
# func = lambda x , y , z : x*x+y*y+z*z
# res = ... # RooFitResult object
# val = re.evaluate ( func , ( 'x' , 'y' , 'z' ) ) 
# @endcode
def _rfr_evaluate_ ( self , func , args , partial = () ) :
    """ Evaluate the certain  function/expression for the fit   result
    >>> func = lambda x , y , z : x*x+y*y+z*z
    >>> res = ... # RooFitResult object
    >>> val = re.evaluate ( func , ( 'x' , 'y' , 'z' ) ) 
    """
    res   = _rfr_results_ ( self , *args )
    val   = res.value()
    cov2  = res.cov2 ()
    N     = len ( val ) 
    from ostap.math.derivative import EvalNVEcov
    func2 = EvalNVEcov ( func , N = N , partial = partial )
    return func2 ( args = val , cov2 = cov2 )
    
# =============================================================================
## print <code>RooFitResult</code> as a table
#  @code
#  result = ...
#  result.table() 
#  @endcode 
def _rfr_table_ ( rr , title = '' , prefix = '' , more_vars = {} ) :
    """ Print RooFitResult  as a table
    >>> result = ...
    >>> result.table() 
    """

    if not rr or not valid_pointer ( rr ) :
        return "Invalid RooFitResul object!"
    
    r = Ostap.Utils.FitResults ( rr )

    from  ostap.fitting.utils    import fit_status, cov_qual
    rows = []

    ##  1. fit status
    status = r.status() 
    if status :
        row = attention ( 'Status' )  , '' , attention ( fit_status ( status ) ) , '' 
        rows.append ( row )
    else :
        row = 'Status'                , '' , allright  ( fit_status ( status ) ) , '' 
        rows.append ( row )

    ## 2. minumum NLL
    s , n = pretty_float ( r.minNll() )
    if n : n = '[10^%+d]' % n
    else : n = '' 

    rows.append ( ( "Minimized FCN/NLL value"    , n , '  ' + s , '' ) )

    s , n = pretty_float ( r.edm () )
    if n : n = '[10^%+d]' % n
    else : n = '' 

    rows.append ( ( 'Estimated distance to minimum' , n , '  ' + s , '' ) )

    cq = r.covQual()
    cn = '' 
    if  -1 == cq :
        cn = cov_qual  ( cq ) 
    elif 3 == cq :
        cn = allright  ( cov_qual ( cq ) )
    elif cq in (  0 , 1 , 2 ) :
        cn = attention ( cov_qual ( cq ) )
    else :
        cn = cov_qual  ( cq ) 
        
    rows.append ( ( 'Covariance matrix quality'     , '' , '  ' + cn , '' ) )
    
    for i in  range ( r.numStatusHistory() ) :
        label =  r.statusLabelHistory ( i )
        code  =  r.statusCodeHistory  ( i )
        row   =  'Status: %s '% label   , '' ,             '%d' % code             
        if not code in ( 0 , -1 ) :
            row = attention  ( row [ 0 ] ) , row [ 1 ] , '   ' + attention ( row [ 2 ] ) , '' 
        else    :
            row =              row [ 0 ]   , row [ 1 ] , '   ' + allright  ( row [ 2 ] ) , ''
        rows.append ( row )

    nbadnll = r.numInvalidNLL()
    if 0 < nbadnll :
        rows.append ( ( 'Invalid FCN/NLL evaluations' , '' , '  %d' % nbadnll , '' ) )


    with_globcorr = True 

    pars_all   = r.params ( float_only = False )
    pars_float = r.params ( float_only = True  )

    #
    ## parameter close to the limits? 
    # 
    limits     = {}
    threshold  = 10
    for p in pars_float :
        
        v , a = pars_float [ p ]

        has_min = hasattr ( a , 'hasMin' ) and a.hasMin()
        has_max = hasattr ( a , 'hasMax' ) and a.hasMax()
        
        if has_min or has_max :
            
            value = v.value()
            
            if a.hasAsymError () : 
                v = ValWithErrors ( a.getVal() , AsymErrors ( negative = a.getAsymErrorLo() , positive = a.getAsymErrorHi() ) )            

            dmin , dmax  = -1, -1
            if has_min and hasattr ( a , 'getMin' ) :
                vmin = a.getMin()
                verr = abs ( v.neg_error )
                if   isequal ( vmin , value )  : dmin = 0   
                elif vmin < value and 0 < verr : dmin = abs ( value - vmin ) / verr 
                    
            if has_max and hasattr ( a , 'getMax' ) :
                vmax = a.getMax()
                verr = abs ( v.pos_error )
                if   isequal ( vmax , value )  : dmax = 0   
                elif value < vmax and 0 < verr : dmax = abs ( value - vmax ) / verr 
                
            dmn = ( 0 <= dmin < threshold ) 
            dmx = ( 0 <= dmax < threshold ) 
            
            if   dmn and dmx : limits [ a.name ] = min ( dmin , dmax ) 
            elif dmn         : limits [ a.name ] = dmin
            elif dmx         : limits [ a.name ] = dmax


    if with_globcorr : header = ( '', 'Unit', 'Value' , 'Global/max correlation [%]') 
    else             : header = ( '', 'Unit', 'Value' , 'Max correlation [%]') 
    
    if limits : header   +=  ( '@limit?', ) 

    rows = [ header ] + rows
    
    ## constant/fix parameters 
    crows = [] 
    for p in pars_all :

        if p in pars_float : continue 
        v , a = pars_all [ p ]

        s , n = pretty_float  ( v.value()  ) 

        if n : n = '[10^%+d]' % n
        else : n = '' 
        row = p , n , '  ' + s + ' (fix)' , ''  
        crows.append ( row ) 

    ## floating parameters
    max_corr = False
    frows    = []
    for p in pars_float :

        v , a = pars_float [ p ]

        if not a.hasAsymError () :
            s , n = v.pretty_print ( parentheses = False ) 
        else :
            s , n = pretty_error2  ( a.getVal() , a.getAsymErrorLo() , a.getAsymErrorHi() , parentheses = False )

        if n : n = '[10^%+d]' % n
        else : n = '' 

        if 0 <= cq and 1 < len ( pars_float ) :

            mxr , mxv = r.max_cor    ( p )

            if with_globcorr :
                
                gc = -1.0 
                ## gc    = r.globalCorr ( p ) if 3 == cq else -1.00
                gc    = r.global_corr ( p ) if 3 == cq else -1.00
                
                if 0 <= gc :  cc = '% +5.1f/(% +5.1f,%s)' % ( gc*100 , mxr*100   , mxv )
                else       :  cc = '% +5.1f : %-s'        % (          mxr * 100 , mxv )                
                if 0.95 < abs ( mxr ) or 0.95 < gc : cc = attention ( cc )

            else :
                
                cc = '% +5.1f : %-s'  % ( mxr * 100 , mxv )                
                if 0.95 < abs ( mxr ) : cc = attention ( cc )

                max_corr = True

            row = p , n , s , cc
            
        else :
            
            row = p , n , s

        dist = limits.get ( a.name , None )
        if not dist is None :
            item = '%.1f%s' % ( dist , '\U000003C3' if show else 's' ) 
            if   dist <= 3 : item = attention ( item )
            elif dist <= 5 : item = attstr    ( item )
            row += ( item , )
            
        frows.append ( row ) 

    ## more parameters
    mrows = []
    for p in sorted ( more_vars ) :

        func  = more_vars [ p ]
        
        v     = func ( r )
        
        s , n = v.pretty_print()  
        
        if n : n = '[10^%+d]' % n
        else : n = '' 

        cc = 'derived'
        row = p , n , s , cc
        mrows.append ( row ) 

    crows.sort()
    frows.sort()

    all = rows + crows + frows + mrows
    
    import ostap.logger.table as T
    if all : all = T.remove_empty_columns ( all ) 
    return T.table ( all , title = title if title else r.GetTitle() , prefix = prefix , alignment = 'lcclcccc' )

# =============================================================================
## 'easy' print of RooFitResult
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _rfr_print_ ( self , opts = 'v' ) :
    """ Easy print of RooFitResult
    >>> result = ...
    >>> print ( result ) 
    """
    if not valid_pointer ( self ) : return 'Invalid RooFitResult'

    ## 1. try to use table print 
    table = _rfr_table_ ( self )
    lmax  = -1
    from ostap.utils.basic      import terminal_size 
    width , _ = terminal_size()
    
    from ostap.logger.colorized import decolorize 
    for row in table :
        lmax = max ( lmax , len ( decolorize ( row ) ) )
        if width < lmax : break

    ## if the table is narrow enough, print it 
    if 10 < lmax and lmax < width : return table 

    ## otherwise use naive RooFit print 
    result = self.print_multiline ( content = 1 , verbose = True )
    if sys.version_info < ( 3 , 0 ) :
        if isinstance ( result , unicode ) :
            result = result.encode ('utf-8')
    return result 
    
# =============================================================================
## Get global correlation coeffcient
#  @code
#  fit_reuslt = ///
#  coeff      = fit_results.global_cc ( 3 ) 
#  @endcode 
#  @see Ostap::Utils::global_cc 
def _rfr_global_cc_ ( rfr , par ) :
    """ Get global correlation coeffcient
    
    - see `Ostap.Utils.global_cc `
    
    >>> fit_reuslt = ///
    >>> coeff      = fit_results.global_cc ( 3 ) 
    """
    if par is None :
        
        ## get all values at once 
        values = Ostap.Utils.global_cc ( rfr ) 
        if len ( values ) != len ( rfr ) :
            raise ValueError( "Error from Ostap::Utils::global_cc" )
        return tuple ( v for v in values )
    
    elif  isinstance  ( par , string_types    ) :
        
        index = rfr.floatParsFinal().index ( par )        
        if index < 0 : raise KeyError ( 'Invalid parameter %s' % par )
        return _rfr_global_cc_ ( rfr , index )
    
    elif isinstance ( par , ROOT.RooAbsArg ) :
        return _rfr_global_cc_ ( rfr , par.name )

    index = par 
    assert isinstance ( index , integer_types ) and 0 <= index  , \
           'Invald index %s/%s' % ( par , index )
    
    if not 0 <= index < len ( rfr ) :
        raise IndexError("Invalid index %s/%s" % ( par , index ) ) 

    return Ostap.Utils.global_cc ( rfr , index )

# =============================================================================
## Get global correlation coefficient for the parameter
#  \f$ \rho_k = \sqrt{    1 - \left[ C_{kk} V_{kk}\right]^{-1} } \f$
#  where \f$ C \f$ is covarinace matrix and \f$ V = C^{-1}\$ is inverse
#  @code
#  resut = ...
#  result.global_corr ( 'S' ) 
#  @endcode 
#  It should be accessible via <code>RooFitResult::globalCorr</code> method
#  but often it results in segfault
#  @see RooFitResult::globalCorr
def _rfr_global_corr_ ( rfr , par ) :
    """ Get global correlation coefficient for the parameter
    rho_k = sqrt ( 1 - 1/( C_{kk} V_{kk}) )
    - where C is covariance  matrix and V is inverse
    
    It should be accessible via `RooFitResult.globalCorr`method
    but often it results in segfauls
    - see `RooFitResult.globalCorr`
    """
    if  isinstance  ( par , string_types    ) :
        index = rfr.floatParsFinal().index ( par ) 
        if index < 0 :                                    return -1
    elif isinstance ( par , integer_types   ) :
        if not 0 <= par < len ( rfr.floatParsFinal () ) : return -1
        index = par 
    elif isinstance ( par , ROOT.RooAbsArg ) :
        return _rfr_global_corr_ ( rfr , par.name )
    else : return -1 
    
    np    = len ( rfr ) 

    return _rfr_global_cc_ ( rfr , par ) 

    ## ## ostap.math.linalg2 machinery fails for large matrices
    ## if 10 >= np or ( np , 'double' ) in LA.LinAlgT.known_ssymmatrices : 

    ##     ## get the covariance matrix and invert it 
    ##     mc = self.covmatrix () 
    ##     if not mc.InvertChol () : return -1   
        
    ##     ## get the covariance matrix 
    ##     cm = self.covmatrix()

    ## else : ## use TMatrixTSym
        

    ##     cm  = self.covarianceMatrix()
    ##     mc  = type ( cm ) ( cm )
    ##     det = ctypes.c_double(0)
    ##     mc.Invert( det )
    ##     det = det.value
    ##     if iszero ( det ) : return -1

                
    ## cv = cm ( index , index ) * mc ( index , index ) 
    ## rho2 = 1.0 - 1.0 / cv
    ## if rho2 < 0 : return -1
    ## return  math.sqrt ( rho2 )


# =============================================================================
## The size of the problem == number of floating parameters
def _rfr_len_ ( rfr ) :
    """ The size of the problem == number  of floating parameters
    """
    return len ( rfr.floatParsFinal() )

# =============================================================================
## Symmetrized Kullback-Leibler divergency between two RooFitResult objects
#  @code
#  r1 = ...
#  r2 = ...
#  r1.kullback( r2 )
#  @endcode
def _rfr_kullback_ ( r1 , r2 ) :
    """ Symmetrized Kullback-Leibler divergency between two RooFitResult objects
    >>> r1 = ...
    >>> r2 = ...
    >>> r1.kullback( r2 )
    """    
    result = -999.0  
    
    ## 1. check structure: constant parameters  
    s1 = set ( p.name for p in r1.constPars() )
    s2 = set ( p.name for p in r2.constPars() )
    if s1 != s2 :
        logger.error ("Constant parameters are different!")
        return result 

    ## 2. check structure: floating parameters  
    s1 = set ( p.name for p in r1.floatParsFinal() )
    s2 = set ( p.name for p in r2.floatParsFinal() )
    if s1 != s2 :
        logger.error ("Floating parameters are different!") 
        return result 

    ## list of parameter names 
    vars = tuple ( sorted ( s1 ) ) 

    ## 3. get the results
    results1 = r1.results ( *vars )
    results2 = r2.results ( *vars )

    ## 4. calculate the divergency
    result = results1.kullback_leibler ( results2 )
    if result <= -999  :
        logger.error ("Error from SVectorWithError::kullack_leibler")

    return result 
    
# =============================================================================
## Run MIGRAD for RooMinimizer object
#  @code
#  pdf = ...
#  minuit = pdf.minuit()
#  minuit.migrad() 
#  @endcode
#  @see RooMinimizer 
def _rm_migrad_  ( self , refit = 0 , tag = "" , minos = () ) :
    """ Run MIGRAD for RooMinimizer
    - see ROOT.RooMinimizer 
    >>> pdf = ...
    >>> minuit = pdf.minuit()
    >>> minuit.migrad() 
    """
    status = self._old_migrad_()
    if 0 != status :
        from   ostap.fitting.utils import fit_status
        message = "MIGRAD %s status %s"  % ( tag , fit_status ( status ) )
        if   isinstance ( refit , integer_types ) and 0 < refit :
            logger.error ( message +  ', refit... ')
            return _rm_migrad_ ( self ,  refit - 1 , tag ) 
        elif isinstance ( refit , bool          ) and     refit :
            logger.error ( message +  ', refit... ')
            return _rm_migrad_ ( self ,  False     , tag )
        logger.error ( message )        
    elif minos :        
        return self.minos ( *minos ) 
        
        
    return status 

# =============================================================================
## run MINOS for set of parameters
#  @code
#  pdf = ...
#  minuit = pdf.minuit( ... )
#  minuit.migrad() 
#  minuit.minos ( 'sigma' , 'mean' )  
#  @endcode
def _rm_minos_ ( self , *variables ) :
    """ Run MINOS for set of parameters
    >>> pdf = ...
    >>> minuit = pdf.minuit( ... )
    >>> minuit.migrad() 
    >>> minuit.minos ( 'sigma' , 'mean' )  
    """
    
    if not variables :
        status = self._old_minos_ ()
        if 0 != status :
            from   ostap.fitting.utils import fit_status 
            logger.error ( "MINOS status %s"  % fit_status ( status ) )
        return status

    aset = ROOT.RooArgSet()
    res  = self.save()

    for v in variables :

        if   isinstance   ( v , string_types ) or isinstance ( v , ROOT.RooAbsReal ) :
            par = res.param ( v ) [ 1 ]
            if par : aset.add ( par )
        elif isinstance ( v , ROOT.RooAbsCollection ) or \
             isinstance ( v , list_types            ) :
            for a in v :
                if isinstance ( a , ROOT.RooAbsReal ) :
                    par = res.param ( a ) [ 1 ]
                    if par : aset.add ( par )
                    
    del res
    
    if aset : status = self._old_minos_ ( aset ) 
    else    : status = self._old_minos_ (      )
    
    if 0 != status :
        from   ostap.fitting.utils import fit_status 
        logger.error ( "MINOS status %s"  % fit_status ( status ) )
        
    return status
    
    
if not hasattr ( ROOT.RooMinimizer , '_old_migrad_' ) :
    ROOT.RooMinimizer._old_migrad_ = ROOT.RooMinimizer.migrad
    _rm_migrad_.__doc__ += '\n' + ROOT.RooMinimizer.migrad.__doc__
    ROOT.RooMinimizer.     migrad  = _rm_migrad_ 

if not hasattr ( ROOT.RooMinimizer , '_old_minos_' ) :
    ROOT.RooMinimizer._old_minos_  = ROOT.RooMinimizer.minos
    _rm_minos_.__doc__ += '\n' + ROOT.RooMinimizer.minos.__doc__
    ROOT.RooMinimizer.     minos   = _rm_minos_ 

# =============================================================================
## make 2D contours  in units of 'sigma'
#  @code
#  pdf    = ...
#  minuit = pdf.minuit( ... )
#  minuit.migrad() 
#  contour = minuit.contour ( 'A' , 'B' ) 
#  @endcode 
def _rm_contour_ ( self                   ,
                   var1                   ,
                   var2                   ,
                   npoints = 100          ,
                   *levels                ) : 

    """ Make 2D contours  in uniys  of 'sigma'
    >>> pdf    = ...
    >>> minuit = pdf.minuit( ... )
    >>> minuit.migrad() 
    >>> contour = minuit.contour ( 'A' , 'B' ) 
    """
    
    assert isinstance ( npoints , integer_types ) and 2 < npoints ,\
           'Invalid number of points %s' % npoints
    
    res  = self.save () 
    if not isinstance ( var1 , ROOT.RooAbsReal ) : 
        var1 = res.param ( var1 ) [ 1 ]
    if not isinstance ( var2 , ROOT.RooAbsReal ) :         
        var2 = res.param ( var2 ) [ 1 ] 
            
    n = 6 * [ 0.0 ]
    for i , l in enumerate ( levels ) :
        if len ( n ) <= i : break
        ll = float ( l )
        assert 0 <= l , 'Invalid level %s' % ll 
        n [ i ] = ll  

    return self._old_contour_ ( var1    , var2    ,
                                n [ 0 ] , n [ 1 ] ,
                                n [ 2 ] , n [ 3 ] ,
                                n [ 4 ] , n [ 5 ] , npoints )


if not hasattr ( ROOT.RooMinimizer , '_old_contour_' ) :
    ROOT.RooMinimizer._old_contour_ = ROOT.RooMinimizer.contour 
    _rm_contour_.__doc__ += '\n' + ROOT.RooMinimizer.contour.__doc__
    ROOT.RooMinimizer. new_contour  = _rm_contour_ 
    ROOT.RooMinimizer. contour      = _rm_contour_ 

# =============================================================================
## some decoration over RooFitResult
## ROOT.RooFitResult . __repr__        = _rfr_print_
## ROOT.RooFitResult . __str__         = _rfr_print_
ROOT.RooFitResult . __repr__         = _rfr_table_
ROOT.RooFitResult . __str__          = _rfr_table_
ROOT.RooFitResult . __call__         = _rfr_param_
ROOT.RooFitResult . __getattr__      = _rfr_getattr_ 
ROOT.RooFitResult . __getitem__      = _rfr_getitem_ 
ROOT.RooFitResult . __iter__         = _rfr_iter_
ROOT.RooFitResult . __contains__     = _rfr_contains_
ROOT.RooFitResult . iteritems        = _rfr_iteritems_
ROOT.RooFitResult . dct_params       = _rfr_dct_params_
ROOT.RooFitResult . parameters       = _rfr_params_
ROOT.RooFitResult . params           = _rfr_params_
ROOT.RooFitResult . param            = _rfr_param_
ROOT.RooFitResult . parameter        = _rfr_param_
ROOT.RooFitResult . results            = _rfr_results_
ROOT.RooFitResult . corr               = _rfr_corr_
ROOT.RooFitResult . cor                = _rfr_corr_
ROOT.RooFitResult . max_cor            = _rfr_max_cor_
ROOT.RooFitResult . max_corr           = _rfr_max_cor_
ROOT.RooFitResult . max_correlation    = _rfr_max_cor_
ROOT.RooFitResult . global_cor         = _rfr_global_corr_
ROOT.RooFitResult . global_corr        = _rfr_global_corr_
ROOT.RooFitResult . global_correlation = _rfr_global_corr_
ROOT.RooFitResult . global_cc          = _rfr_global_cc_
ROOT.RooFitResult . cov                = _rfr_cov_
ROOT.RooFitResult . covariance         = _rfr_cov_
ROOT.RooFitResult . cov_matrix         = _rfr_cov_matrix_
ROOT.RooFitResult . covmatrix          = _rfr_covmatrix_ 
ROOT.RooFitResult . parValue           = lambda s,n : s.parameter(n)[0]
ROOT.RooFitResult . sum                = _rfr_sum_
ROOT.RooFitResult . plus               = _rfr_sum_
ROOT.RooFitResult . multiply           = _rfr_multiply_
ROOT.RooFitResult . product            = _rfr_multiply_
ROOT.RooFitResult . subtract           = _rfr_subtract_
ROOT.RooFitResult . minus              = _rfr_subtract_
ROOT.RooFitResult . divide             = _rfr_divide_
ROOT.RooFitResult . ratio              = _rfr_divide_
ROOT.RooFitResult . fraction           = _rfr_fraction_
ROOT.RooFitResult . asymmetry          = _rfr_asymmetry_
ROOT.RooFitResult . pow                = _rfr_pow_
ROOT.RooFitResult . power              = _rfr_pow_
ROOT.RooFitResult . evaluate           = _rfr_evaluate_

ROOT.RooFitResult . kullback           = _rfr_kullback_
ROOT.RooFitResult . kullback_leibler   = _rfr_kullback_
ROOT.RooFitResult . cov_eigenvalues    = _rfr_eigenvalues_

ROOT.RooFitResult . table              = _rfr_table_
ROOT.RooFitResult . __len__            = _rfr_len_ 

_new_methods_ += [
    ROOT.RooFitResult . __repr__            ,
    ROOT.RooFitResult . __str__             ,
    ROOT.RooFitResult . __call__            ,
    ROOT.RooFitResult . __getattr__         ,
    ROOT.RooFitResult . __iter__            ,
    ROOT.RooFitResult . __contains__        ,
    ROOT.RooFitResult . iteritems           ,
    ROOT.RooFitResult . parameters          ,
    ROOT.RooFitResult . params              ,
    ROOT.RooFitResult . dct_params          ,
    ROOT.RooFitResult . results             ,
    ROOT.RooFitResult . results             ,
    ROOT.RooFitResult . param               ,
    ROOT.RooFitResult . parameter           ,
    #
    ROOT.RooFitResult . corr                ,
    ROOT.RooFitResult . cor                 ,
    #
    ROOT.RooFitResult . max_cor             ,
    ROOT.RooFitResult . max_corr            ,
    ROOT.RooFitResult . max_correlation     ,
    #
    ROOT.RooFitResult . global_cor          ,
    ROOT.RooFitResult . global_corr         ,
    ROOT.RooFitResult . global_correlation  ,
    ROOT.RooFitResult . global_cc           ,
    # 
    ROOT.RooFitResult . cov              ,
    ROOT.RooFitResult . covariance       ,
    ROOT.RooFitResult . parValue         ,
    ROOT.RooFitResult . sum              ,
    ROOT.RooFitResult . plus             ,
    ROOT.RooFitResult . multiply         ,
    ROOT.RooFitResult . product          ,
    ROOT.RooFitResult . subtract         ,
    ROOT.RooFitResult . minus            ,
    ROOT.RooFitResult . divide           ,
    ROOT.RooFitResult . ratio            ,
    ROOT.RooFitResult . fraction         ,
    ROOT.RooFitResult . asymmetry        ,
    ROOT.RooFitResult . pow              ,
    ROOT.RooFitResult . power            ,
    ROOT.RooFitResult . evaluate         ,
    ROOT.RooFitResult . table            ,
    ROOT.RooFitResult . cov_eigenvalues  ,
    ROOT.RooFitResult . covmatrix        ,
    ROOT.RooFitResult . kullback         , 
    ROOT.RooFitResult . kullback_leibler ,
    ]

# =============================================================================
_decorated_classes_ = (
    ROOT.RooFitResult , 
    ROOT.RooMinimizer , 
    )

_new_methods_ = tuple ( _new_methods_ ) 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                     The END 
# =============================================================================
