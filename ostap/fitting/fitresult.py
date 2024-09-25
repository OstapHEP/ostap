#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/fitresult.py
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
from   ostap.math.ve    import VE 
from   ostap.core.ostap_types import integer_types
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
## representation of TFitResult object 
#  @code 
#  fit_result = histo.Fit( func , 'S' , ... )
#  print fit_result
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _fit_repr_ ( self ) :
    """ Representaion of TFitResult object
    >>> fit_result = histo.Fit( func , 'S' , ... )
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

# ==============================================================================
## print <code>TFitResult</code> as a table 
def _fit_table_ ( rfit , title = '' , prefix = '' ) :
    """ Print <code>TFitResult</code> as a table
    """
    from  ostap.fitting.utils    import fit_status, cov_qual
    from ostap.logger.colorized  import attention, allright
    from ostap.logger.pretty     import pretty_float, pretty_ve, fmt_pretty_ve 
    
    header = ( '', 'Unit' , 'Value' )

    rows = []

    ##  0. minimized type
    row = "Minimizer Type" , '' , rfit.MinimizerType() 
    rows.append ( row )

    ##  0. minimized type
    v = rfit.IsValid()
    if v : row = "Valid"   , '' , 'True'
    else : row = "Valid"   , '' , attention ( 'False') 
    rows.append ( row )
    
    ##  1. fit status
    status = rfit.Status() 
    if status :
        row = attention ( 'Status' )  , '' , attention ( fit_status ( status ) ) 
        rows.append ( row )
    else :
        row =             'Status' , '' , allright ( fit_status ( status ) )   
        rows.append ( row )

    ## 4. covariance status
    cq = rfit.CovMatrixStatus() 
    cn = '' 
    if  -1 == cq              : cn = cov_qual  ( cq ) 
    elif 3 == cq              : cn = allright  ( cov_qual ( cq ) )
    elif cq in (  0 , 1 , 2 ) : cn = attention ( cov_qual ( cq ) )
    else                      : cn = cov_qual  ( cq )         
    rows.append ( ( 'Covariance matrix quality'     , '' , '  ' + cn  ) )


    ## 3-6. chi2,nDoF,chi2/nDoF,minFCN
    chi2 = rfit.Chi2 () 
    s , n = pretty_float ( chi2 )
    if n : n = '[10^%+d]' % n
    else : n = '' 
    rows.append ( ( "Chi2"       , n , '  ' + s ) )
    ##
    ndf = rfit.Ndf()
    rows.append ( ( "nDoF"       , '' , '  ' + '%d' % ndf   ) )
    ##
    c2ndf = rfit.Chi2 () /  ndf  
    s , n = pretty_float ( c2ndf  )
    if n : n = '[10^%+d]' % n
    else : n = '' 
    rows.append ( ( "Chi2/nDoF"  , n , '  ' + s   ) )
    ##
    minfcn = rfit.MinFcnValue() 
    s , n = pretty_float ( minfcn  )
    if n : n = '[10^%+d]' % n
    else : n = '' 
    rows.append ( ( "Minimal FCN"  , n , '  ' + s   ) )

    ## 7.Probability in %[%]
    prob = rfit.Prob() / 100
    rows.append ( ( "Probability"  , '[%]' , '  %5.3e' % prob ) )

    ## 8. distrance to minimum 
    edm  = rfit.Edm()
    s , n = pretty_float ( edm  )
    if n : n = '[10^%+d]' % n
    else : n = '' 
    rows.append ( ( "Estimated distance to minimum" , n , '  ' + s ) )

    ncalls = rfit.NCalls()
    rows.append ( ( "FCN calls" , '' , '  ' + '%d' % ncalls  ) )
    ##
    
    has_minos = False
    for i in rfit :
        if not rfit.HasMinosError( i ) : continue 
        has_minos = True
        break

    if has_minos :
        rows   = [ row + ('','','') for row in rows ]
        header = header + ( 'neg-minos' , 'pos-minos' , 'Global corr.' )   
    else         :
        rows = [ row + ('',)      for row in rows ] 
        header = header + ( 'Global corr.' , )   

    for i in rfit :
        
        pname  = rfit.GetParameterName ( i )
        value  = rfit.Value            ( i )  
        
        fixed  = rfit.IsParameterFixed ( i )
        fmte = '' 
        if fixed :
            v = value
            s  , n = pretty_float ( v )
            s = s + '(fixed)'
            nv = n
        else :
            error = rfit.Error ( i )
            v     = VE ( value , error * error )
            ##
            fmt , fmtv , fmte , n = fmt_pretty_ve  ( v , parentheses = False )
            s  = fmt % ( value / 10**n , error / 10**n )
            nv = n 
        if n : n = '[10^%+d]' % n
        else : n = '' 
        pname = "%-2d: %s"% ( i , pname )

        row = pname , n , '  ' + s             
        if not fixed and rfit.HasMinosError( i ) :
            if fmte : 
                error_low  = fmte   % ( rfit.LowerError ( i ) / 10**nv ) 
                error_up   = fmte   % ( rfit.UpperError ( i ) / 10**nv )
            else : 
                error_low  = "%+8g" % ( rfit.LowerError ( i ) / 10**nv ) 
                error_up   = "%+8g" % ( rfit.UpperError ( i ) / 10**nv ) 
        else :
            error_low  = '' 
            error_up   = '' 

        if has_minos  :
            row  = row  + ( error_low , error_up )

            
        gcc    = rfit.GlobalCC         ( i ) * 100 
        gcc  = '%+5.1f%%' % gcc
        row = row + ( gcc, )
        
        rows.append ( row ) 
        
    if not title :  title = rfit.GetTitle()
    
    import ostap.logger.table as T

    rows  = [ header ] + rows  

    return T.table ( rows , title = title , prefix = prefix )


# =============================================================================
## get number of parameters
#  @code 
#  fit_result = histo.Fit( func , 'S' , ... )
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
    """ Get parameter number:
    >>> r    = h1.Fit( ... )
    >>> name = r.GetParNumber ( 'mass' ) 
    """ 
    if isinstance ( par , integer_types ) :
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
    """ Check parameter
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
    """ Get correlation coefficient for parameters 'i' and 'j'
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
    """ Get correlation matrix 
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


for klass in ( ROOT.TFitResult    ,
               ROOT.TFitResultPtr ) :  
    
    klass.__contains__ = _fit_contains_ 
    klass.__repr__     = _fit_repr_ 
    klass.__str__      = _fit_repr_ 
    klass.__iter__     = _fit_iter_ 
    klass.__getitem__  = _fit_getitem_ 
    klass.__call__     = _fit_getitem_ 
    klass.__len__      = _fit_len_ 
    klass.cor          = _fit_cor_ 
    klass.corMtrx      = _fit_corm_ 
    klass.GetParNumber = _fit_parnum_ 
    klass.parnum       = _fit_parnum_ 
    klass.table           = _fit_table_ 

# =============================================================================
_decorated_classes_ = (
    ROOT.TFitResult    , 
    ROOT.TFitResultPtr , 
    )

_new_methods_       = (
    ##
    ROOT.TFitResult   .__contains__  ,
    ROOT.TFitResult   .__repr__      ,
    ROOT.TFitResult   .__str__       ,
    ROOT.TFitResult   .__iter__      ,
    ROOT.TFitResult   .__getitem__   , 
    ROOT.TFitResult   .__call__      , 
    ROOT.TFitResult   .__len__       ,
    ROOT.TFitResult   .cor           ,
    ROOT.TFitResult   .corMtrx       ,
    ROOT.TFitResult   .GetParNumber  ,
    ROOT.TFitResult   .parnum        ,
    ROOT.TFitResult   .table         ,
    ##
    ROOT.TFitResultPtr.__contains__  ,
    ROOT.TFitResultPtr.__repr__      ,
    ROOT.TFitResultPtr.__str__       ,
    ROOT.TFitResultPtr.__iter__      ,
    ROOT.TFitResultPtr.__getitem__   , 
    ROOT.TFitResultPtr.__call__      , 
    ROOT.TFitResultPtr.__len__       ,
    ROOT.TFitResultPtr.cor           ,
    ROOT.TFitResultPtr.corMtrx       ,
    ROOT.TFitResultPtr.GetParNumber  ,
    ROOT.TFitResultPtr.parnum        ,
    ROOT.TFitResultPtr.table         ,
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
