#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Tiny decoration for ROOT.FitResult object
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Tiny decoration for ROOT.FitResult object
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
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.fitresult' )
else                       : logger = getLogger( __name__                  )
# =============================================================================
logger.debug ( 'Tiny decoration for ROOT.FitResult object')
# =============================================================================
from ostap.math.ve import VE 
# =============================================================================
## representation of TFitResult object 
#  @code 
#  fit_result = hiisto.Fit( func , 'S' , ... )
#  print fit_result
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _fit_repr_ ( self ) :
    """Representaion of TFitResult object
    >>> fit_result = hiisto.Fit( func , 'S' , ... )
    >>> print fit_result
    """
    _r  = ''
    _r += "\n Status      = %s "    %   self.Status ()
    _r += "\n Chi2/nDoF   = %s/%s " % ( self.Chi2   () , self.Ndf() ) 
    _r += "\n Probability = %s "    %   self.Prob   () 
    _p = self.Parameters ()
    _e = self.Errors     ()
    for i in range( 0 , len(_p) ) :
        v = _p[i]
        e = _e[i]
        a = VE ( v ,e*e )
        _r  += " \n %s " % a 
    return _r

# =============================================================================
## get number of parameters
#  @code 
#  fit_result = hiisto.Fit( func , 'S' , ... )
#  print len(fit_result)
#  @endcode 
def _fit_len_ ( r ) :
    """Get number of parameters
    >>> fit_result = hiisto.Fit( func , 'S' , ... )
    >>> print len(fit_result)
    """
    return len ( r.Parameters() ) 

# =============================================================================
## iterator over fit-result object 
#  @code 
#  fit_result = hiisto.Fit( func , 'S' , ... )
#  for i in fit_results : print i,fit_results[i] 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _fit_iter_ ( r ) :
    """Iterator over fit-result object
    >>> fit_result = hiisto.Fit( func , 'S' , ... )
    >>> for i in fit_results : print i,fit_results[i] 
    """
    l = len( r )
    i = 0
    while i < l :
        yield i
        i += 1

# =============================================================================
## get parameter number
#  @code
#  r    = h1.Fit( ... )
#  name = r.GetParNumber ( 'mass' ) 
#  @endcode
def _fit_parnum_ ( self , par ) : 
    """Get parameter number:
    >>> r    = h1.Fit( ... )
    >>> name = r.GetParNumber ( 'mass' ) 
    """ 
    if isinstance ( par , ( int , long ) ) :
        if 0<= par< len ( self ) : return int( par )   ## RETURN 
        else                     : return       -1     ## RETURN 
    #
    if isinstance   ( par , str )  :
        ll = len ( self )
        for i in range ( 0 , ll ) :
            if self.ParName(i) == par : return i       ## RETURN 
            
    ## nothing is found 
    return -1                                          ## RETURN 

# =============================================================================
## check parameter
#  @code
#  r = h1.Fit(....) ##
#  if  i  in r :   ...  ## check parameter by index  
#  if 'm' in r :   ...  ## check parameter by name  
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-12
def _fit_contains_ ( self , par ) :
    """Check parameter
    >>> r = h1.Fit(....) ##
    >>> if i   in r :   ...  ## check parameter by index  
    >>> if 'm' in r :   ...  ## check parameter by name  
    """
    return  0 <= _fit_parnum_ ( self , par )

    
# =============================================================================
## getitem for fit-result-object
#  @code
#  r = h1.Fit(....) ##
#  print r[0]  ## print 0th parameter 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _fit_getitem_ ( self , par ) :
    """Getitem for fit-result-object            
    >>> r = h1.Fit(....) ##
    >>> print r[0]  ## print 0th parameter 
    """
    ## convert parameter into integer 
    ipar = _fit_parnum_ ( self , par )
    if not 0<= ipar : raise IndexError("TFitResult:illegal index %s" % par)
    #
    _p = self.Parameter ( ipar )
    _e = self.Error     ( ipar )
    #
    return VE( _p , _e * _e )

# =============================================================================
## Get correlation coefficient for parameters 'i' and 'j'
#  @code
#  r = ...
#  print r.cor(1,2)
#  @endcode 
def _fit_cor_ ( self , i , j ) :
    """Get correlation coefficient for parameters 'i' and 'j'
    >>> r = ...
    >>> print r.cor(1,2)
    """
    ipar = _fit_parnum_ ( self , i )
    jpar = _fit_parnum_ ( self , j )
    #
    if  0 > ipar : raise IndexError( "TFitResult:invalid index %s" % i )
    if  0 > jpar : raise IndexError( "TFitResult:invalid index %s" % j )
    #
    _cij = self.CovMatrix ( ipar , jpar )
    _ei  = self.Errors    ( ipar )
    _ej  = self.Errors    ( jpar )
    ##
    if 0 == _ei or 0 == _ej : return 0   ## RETURN 
    #
    return _cij / ( _ei * _ej ) 

# =============================================================================
## Get correlation matrix 
#  @code
#  r = ...
#  print r.corMatrix()
#  @endcode 
def _fit_corm_ ( self , root = False ) :
    """Get correlation matrix 
    >>> r = ...
    >>> print r.corMtrx ()
    """
    _l = len (self) 
    matrix = None

    import ostap.math.linalg
    try :
        matrix = Ostap.Math.SymMatrix(_l)
    except :
        pass    

    ## fill matrix 
    for i in range (0,_l):
        for j in range (i, _l):
            _cij = self.CovMatrix( i , j ) 
            _eij = self.Error( i ) * self.Error( j )
            if 0 != _eij : _vij = _cij / _eij
            else         : _vij = 0
            matrix [ i , j ] = _vij 
            matrix [ j , i ] = _vij
            
    return matrix


ROOT.TFitResultPtr.__contains__ = _fit_contains_ 
ROOT.TFitResultPtr.__repr__     = _fit_repr_ 
ROOT.TFitResultPtr.__str__      = _fit_repr_ 
ROOT.TFitResultPtr.__iter__     = _fit_iter_ 
ROOT.TFitResultPtr.__getitem__  = _fit_getitem_ 
ROOT.TFitResultPtr.__call__     = _fit_getitem_ 
ROOT.TFitResultPtr.__len__      = _fit_len_ 
ROOT.TFitResultPtr.cor          = _fit_cor_ 
ROOT.TFitResultPtr.corMtrx      = _fit_corm_ 
ROOT.TFitResultPtr.GetParNumber = _fit_parnum_ 
ROOT.TFitResultPtr.parnum       = _fit_parnum_ 


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
