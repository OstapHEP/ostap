#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
# =============================================================================
## @file
#  Module with decoration of asome RooFit objects for efficient use in python
#  - iterators  for RooArgList
#  - iterators  for RooArgSet
#  - iterators  for RooAbsData
#  - decorators for RooRealVar
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
#
# =============================================================================
"""Decoration of some RooFit objects for efficient use in python
- iterators  for RooArgList
- iterators  for RooArgSet
- iterators  for RooAbsData
- decorators for RooRealVar
- and a lot of other stuff 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'setStorage' , ## define the default storage for  RooDataStore 
    'PDF_fun'    , ## wrapper of PDF to ``simple'' function 
    'SETVAR'     , ## context manager to preserev the current value for RooRealVar
    ) 
# =============================================================================
import ROOT
from   ostap.core.core import cpp, VE, hID, dsID   
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.rootfit' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug( 'Some useful decorations for RooFit objects')
# =============================================================================
## iterator for RooArgList 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2011-06-07
def _ral_iter_ ( self ) :
    """Iterator for RooArgList:
    >>> arg_list = ...
    >>> for p in arg_list : print p    
    """
    l = len ( self )
    for i in range ( 0 , l )  :
        yield self[i]

# =============================================================================
## some decoration over RooArgList 
ROOT.RooArgList . __len__       = lambda s   : s.getSize()
ROOT.RooArgList . __contains__  = lambda s,i :  0<= i < len(s)
ROOT.RooArgList . __iter__      = _ral_iter_
ROOT.RooArgList . __nonzero__   = lambda s   : 0 != len ( s ) 
#
# =============================================================================
# helper function 
def _rs_list_ ( self ) :
    _l = []
    for i in self :
        
        if   hasattr  ( i , 'GetName' ) and hasattr ( i , 'getVal' ) :
            _l.append ( i.GetName() + ":%s" % i.getVal() )
        elif hasattr  ( i , 'GetName' ) :
            _l.append ( i.GetName()   )
        elif hasattr  ( i , 'getVal'  ) :
            _l.append ( "%s" % i.getVal ()  )
        else :
            _l.append (  str ( i )    )
            
    return _l ;


# =============================================================================
## printout for   RooArgList 
ROOT.RooArgList . __str__       = lambda s : str ( _rs_list_ ( s ) )  
ROOT.RooArgList . __repr__      = lambda s : str ( _rs_list_ ( s ) )  

# =============================================================================
## iterator for RooArgSet
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _ras_iter_ ( self ) :
    """Simple iterator for RootArgSet:
    >>> arg_set = ...
    >>> for i in arg_set : print i    
    """    
    it  = cpp.Ostap.Utils.Iterator ( self )
    val = it.Next()
    while val :
        yield val 
        val = it.Next()
        
    del it

# =============================================================================
## get the attibute for RooArgSet 
def _ras_getattr_ ( self , aname ) :
    """Get the attibute from RooArgSet
    >>> aset = ...
    >>> print aset.pt    
    """
    _v = self.find ( aname )
    if not _v : raise  AttributeError
    return _v 

# =============================================================================
## get the attibute for RooArgSet 
def _ras_getitem_ ( self , aname ) :
    """Get the attibute from RooArgSet
    >>> aset = ...
    >>> print aset['pt']    
    """
    _v = self.find ( aname )
    if not _v : raise  IndexError
    return _v 

# =============================================================================
## check the presence of variable in set 
def _ras_contains_ ( self , aname ) :
    """Check the presence of variable in set 
    """
    _v = self.find ( aname )
    if not _v : return False 
    return             True 

# =============================================================================
## some decoration over RooArgSet 
ROOT.RooArgSet . __len__       = lambda s   : s.getSize()
ROOT.RooArgSet . __iter__      = _ras_iter_ 
ROOT.RooArgSet . __getattr__   = _ras_getattr_ 
ROOT.RooArgSet . __getitem__   = _ras_getitem_ 
ROOT.RooArgSet . __contains__  = _ras_contains_ 
ROOT.RooArgSet . __nonzero__   = lambda s   : 0 != len ( s ) 
        
ROOT.RooArgSet . __str__   = lambda s : str ( tuple ( _rs_list_ ( s ) ) )  
ROOT.RooArgSet . __repr__  = lambda s : str ( tuple ( _rs_list_ ( s ) ) )  

# =============================================================================
## iterator for RooAbsData
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _rad_iter_ ( self ) :
    """Iterator for RooAbsData
    >>> dataset = ...
    >>> for i in dataset : ... 
    """
    _l = len ( self )
    for i in xrange ( 0 , _l ) :
        yield self.get ( i )

# =============================================================================
## access to the entries in  RooAbsData
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-03-31
def _rad_getitem_ ( self , i ) :
    """Get the entry from RooDataSet
    >>> dataset = ...
    >>> event = dataset[4]
    """
    if 0<= i < len ( self ) :
        return self.get ( i )
    raise IndexError 

# =============================================================================
## Get variables in form of RooArgList 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-03-31
def _rad_vlist_ ( self ) :
    """Get variables in form of RooArgList 
    """
    vlst     = ROOT.RooArgList()
    vset     = self.get()
    for v in vset : vlst.add ( v )
    #
    return vlst

# =============================================================================
## check the presence of variable with given name in dataset 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-03-31
def _rad_contains_ ( self , aname ) :
    """Check the presence of variable in dataset    
    >>> if 'mass' in dataset : print 'ok!'
    """
    vset = self.get()
    return aname in vset 

# =============================================================================
## some decoration over RooDataSet 
ROOT.RooAbsData . varlist       = _rad_vlist_
ROOT.RooAbsData . varlst        = _rad_vlist_
ROOT.RooAbsData . vlist         = _rad_vlist_
ROOT.RooAbsData . vlst          = _rad_vlist_
ROOT.RooAbsData . varset        = lambda s : s.get()

ROOT.RooAbsData . __len__       = lambda s   : s.numEntries()
ROOT.RooAbsData . __nonzero__   = lambda s   : 0 != len ( s ) 
ROOT.RooAbsData . __contains__  = _rad_contains_
ROOT.RooAbsData . __iter__      = _rad_iter_ 
ROOT.RooAbsData . __getitem__   = _rad_getitem_ 

from ostap.trees.trees import _stat_var_, _stat_cov_ , _sum_var_, _sum_var_old_
ROOT.RooAbsData . statVar       = _stat_var_ 
ROOT.RooAbsData . sumVar        = _sum_var_ 
ROOT.RooAbsData . sumVar_       = _sum_var_old_ 

_new_methods_ = [
   ROOT.RooAbsData . varlist       ,
   ROOT.RooAbsData . varlst        ,
   ROOT.RooAbsData . vlist         ,
   ROOT.RooAbsData . vlst          ,
   ROOT.RooAbsData . varset        ,
   #
   ROOT.RooAbsData . __len__       ,
   ROOT.RooAbsData . __nonzero__   ,
   ROOT.RooAbsData . __contains__  ,
   ROOT.RooAbsData . __iter__      ,
   ROOT.RooAbsData . __getitem__   ,
   #
   ROOT.RooAbsData . statVar       ,
   ROOT.RooAbsData . sumVar        ,
   ROOT.RooAbsData . sumVar_       ,
   ]

# =============================================================================
## ``easy'' print of RooFitResult
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _rfr_print_ ( self , opts = 'v' ) :
    """Easy print of RooFitResult
    >>> result = ...
    >>> print result    
    """
    self.Print( opts )
    return 'RooFitResult'

# =============================================================================
## get parameters from RooFitResult
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _rfr_params_ ( self , float_only = True ) :
    """GetParameters from RooFitResult:
    >>> result = ...
    >>> params = results
    >>> p0     = params['A'][0]  ## get the value
    >>> p0s    = params['A'][1]  ## get the parameter itself     
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
def _rfr_corr_  ( self , name1 , name2 ) :
    """Get correlation coefficient for two parameter 
    >>> result = ...
    >>> corr = results.corr('Signal', 'Background')
    >>> print corr
    """
    p1 = self.param ( name1 )
    p2 = self.param ( name2 )
    #
    return self.correlation ( p1[1] , p2[1] ) 

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
    p1   = self.param ( name1 ) 
    p2   = self.param ( name2 ) 
    args = ROOT.RooArgList ( p1[1] , p2[1] ) 
    return self.reducedCovarianceMatrix (  args )


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
## get correct estimate of sum of two variables,
#  taking into account correlations
#  @code
#  >>> r = ...
#  >>> print r.sum( 'S' , 'B' )  ## S+B
#  @endcode
#  @see Gaudi:Math::sum 
def _rfr_sum_ ( self , var1 , var2 ) :
    """Get correct estimate of sum of two variables,
    taking into account correlations
    >>> r = ...
    >>> print r.sum( 'S' , 'B' ) ## S+B
    """
    _v1  = self.param ( var1 )[0]
    _v2  = self.param ( var2 )[0]
    _cor = self.corr  ( var1 , var2 ) 
    return cpp.Ostap.Math.sum ( _v1 , _v2 , _cor ) 
   
# ===========================================================================
## get correct estimate of product of two variables,
#  taking into account correlations
#  @code
#  >>> r = ...
#  >>> print r.multiply( 'S' , 'B' ) ## S*B
#  @endcode
#  @see Gaudi:Math::multiply 
def _rfr_multiply_ ( self , var1 , var2 ) :
    """Get correct estimate of product of two variables,
    taking into account correlations
    >>> r = ...
    >>> print r.multiply( 'S' , 'B' ) ## S*B
    """
    _v1  = self.param ( var1 )[0]
    _v2  = self.param ( var2 )[0]
    _cor = self.corr  ( var1 , var2 ) 
    return cpp.Ostap.Math.multiply ( _v1 , _v2 , _cor ) 
    
# ===========================================================================
## get correct estimate of division  of two variables,
#  taking into account correlations
#  @code
#  >>> r = ...
#  >>> print r.divide( 'S' , 'B' ) ## S/B
#  @endcode
#  @see Gaudi:Math::divide
def _rfr_divide_ ( self , var1 , var2 ) :
    """Get correct estimate of division of two variables,
    taking into account correlations
    >>> r = ...
    >>> print r.divide( 'S' , 'B' ) ## S/B
    """
    _v1  = self.param ( var1 )[0]
    _v2  = self.param ( var2 )[0]
    _cor = self.corr  ( var1 , var2 ) 
    return cpp.Ostap.Math.divide ( _v1 , _v2 , _cor ) 

# ===========================================================================
## get correct estimate of subtraction of two variables,
#  taking into account correlations
#  @code
#  >>> r = ...
#  >>> print r.subtract( 'S' , 'B' ) ## S-B
#  @endcode
#  @see Gaudi:Math::subtract
def _rfr_subtract_ ( self , var1 , var2 ) :
    """Get correct estimate of subtraction of two variables,
    taking into account correlations
    >>> r = ...
    >>> print r.subtract( 'S' , 'B' ) ## S-B
    """
    _v1  = self.param ( var1 )[0]
    _v2  = self.param ( var2 )[0]
    _cor = self.corr  ( var1 , var2 ) 
    return cpp.Ostap.Math.subtract ( _v1 , _v2 , _cor ) 

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
    _av1  = abs ( self.param ( var1 )[0].value() ) 
    _av2  = abs ( self.param ( var2 )[0].value() ) 
    if _av1 > _av2 :
        return 1 / ( 1 + _rfr_divide ( self , var2 , var1  ) )
    return 1 - _rfr_fraction_ ( self , var2 , var1 ) 

# ============================================================================
## get the required results in form of SVectorWithError object
#  @code
#  fit_resuts = ...
#  res   = fit_results.results( 'A', 'B' , 'C' )
#  print res, res.cov2() 
#  @endcode
#  @see Gaudi::Math::SVectorWithError
def _rfr_results_( self , *vars ) :
    """Get the required results in form of SVectorWithError object
    >>> fit_resuts = ...
    >>> res   = fit_results.results( 'A', 'B' , 'C' )
    >>> print res, res.cov2() 
    """
    _n = len ( vars )
    _r = cpp.Gaudi.Math.SVectorWithError(_n,'double')()
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
## fix parameter at some value
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-09-20
def _fix_par_ ( var , value  = None ) :
    """Fix parameter at some value :

    >>> par = ...
    >>> par.fix ( 10 )     
    """
    #
    if None is value :
        var.setConstant( True )
        return var.ve()
    
    if hasattr ( value , 'value' ) : value = value.value()
    #
    var.setVal      ( value )
    var.setConstant ( True  )
    #
    return var.ve() 

# =============================================================================
## release the parameter
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-09-20
def _rel_par_ ( var )  :
    """Release the parameters

    >>> par = ...
    >>> par.release ()     
    """
    var.setConstant ( False )
    #
    return var.ve()

# ==============================================================================
## Convert RooRealVar into ValueWithError 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-23
def _rrv_ve_ ( var ) :
    """Convert RooRealVar into ValueWithError
    
    >>> par = ...
    >>> ve  = par.ve()    
    """
    v  =      var.getVal()
    e2 = 0 if var.isConstant() else var.getError()**2
    #
    return VE ( v , e2 )

# =============================================================================
## decorate RooRealVar:
ROOT.RooRealVar   . as_VE     = _rrv_ve_ 
ROOT.RooRealVar   . asVE      = _rrv_ve_ 
ROOT.RooRealVar   . ve        = _rrv_ve_
ROOT.RooRealVar   . fix       = _fix_par_
ROOT.RooRealVar   . Fix       = _fix_par_
ROOT.RooRealVar   . release   = _rel_par_
ROOT.RooRealVar   . Release   = _rel_par_
## convert to float 
ROOT.RooRealVar   . __float__ = lambda s : s.getVal()
## print it in more suitable form 
ROOT.RooRealVar   . __repr__  = lambda s : "'%s' : %s " % ( s.GetName() , s.ve() )

ROOT.RooRealVar   . xmin      = lambda s : s.getMin()
ROOT.RooRealVar   . xmax      = lambda s : s.getMax()
ROOT.RooRealVar   . minmax    = lambda s : (s.xmin(),s.xmax()) 

ROOT.RooConstVar    .as_VE    = lambda s : VE( s.getVal() , 0 )
ROOT.RooFormulaVar  .as_VE    = lambda s : VE( s.getVal() , 0 )
ROOT.RooConstVar    .asVE     = lambda s : VE( s.getVal() , 0 )
ROOT.RooFormulaVar  .asVE     = lambda s : VE( s.getVal() , 0 )

_new_methods_ += [
    ROOT.RooRealVar   . as_VE     ,
    ROOT.RooRealVar   . asVE      ,
    ROOT.RooRealVar   . ve        ,
    ROOT.RooRealVar   . fix       ,
    ROOT.RooRealVar   . Fix       ,
    ROOT.RooRealVar   . release   ,
    ROOT.RooRealVar   . Release   ,
    ## convert to float 
    ROOT.RooRealVar   . __float__ ,
    ## print it in more suitable form 
    ROOT.RooRealVar   . __repr__  ,
    #
    ROOT.RooRealVar   . xmin      ,
    ROOT.RooRealVar   . xmax      ,
    ROOT.RooRealVar   . minmax    ,
    #
    ROOT.RooConstVar    .as_VE    ,
    ROOT.RooFormulaVar  .as_VE    ,
    ROOT.RooConstVar    .asVE     ,
    ROOT.RooFormulaVar  .asVE     ,
    #
    ]
# =============================================================================
## Prepare ``soft'' gaussian constraint for the given variable
#  @code 
#    >>> var     = ...                            ## the variable 
#    >>> extcntr = var.constaint( VE(1,0.1**2 ) ) ## create constrains 
#    >>> model.fitTo ( ... , extcntr )            ## use it in the fit 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-23
def _rar_make_constraint_ ( v , value ) :
    """Prepare ``soft'' gaussian constraint for the variable
    
    >>> var     = ...                            ## the variable 
    >>> extcntr = var.constaint( VE(1,0.1**2 ) ) ## create constrains 
    >>> model.fitTo ( ... , extcntr )            ## use it in the fit 
    """
    #
    #
    ## create gaussian constrains
    #
    vn       = 'Constr(%s)' % v.GetName()
    vt       = 'Gauissian constraint(%s) at %s' % ( v.GetName() , value )
    #
    v._cvv   = ROOT.RooFit.RooConst ( value.value () )  ## NB! 
    v._cve   = ROOT.RooFit.RooConst ( value.error () )  ## NB! 
    v._cntr  = ROOT.RooGaussian     ( vn , vt , v , v._cvv , v._cve )
    #
    ## keep it 
    v._cntrs = ROOT.RooArgSet       ( v._cntr )
    #
    return ROOT.RooFit.ExternalConstraints ( v._cntrs ) 

ROOT.RooAbsReal. constraint = _rar_make_constraint_

_new_methods_ += [
    ROOT.RooAbsReal. constraint 
    ]
# ============================================================================
## make a histogram for RooRealVar
#  @see RooRealVar
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-07-14
def _rrv_as_H1_ ( v , bins = 100 , double = True ) :
    """Make TH1 histogram from RooRealVar
    
    >>> variable = ...
    >>> histo = variable.histo ( 100 )    
    """
    _hT = ROOT.TH1D if double else ROOT.TH1F
    _h  = _hT ( hID() , v.GetTitle() , bins , v.getMin()  , v.getMax() )
    _h.Sumw2()
    
    return _h 

ROOT.RooRealVar   . histo = _rrv_as_H1_
ROOT.RooRealVar   . asH1  = _rrv_as_H1_

_RRV_ = ROOT.RooRealVar

_new_methods_ += [
    ROOT.RooRealVar.histo , 
    ROOT.RooRealVar.asH1  
    ]

# ============================================================================
## Addition of RooRealVar and ``number''
def _rrv_add_ ( s , o ) :
    """Addition of RooRealVar and ``number''

    >>> var = ...
    >>> num = ...
    >>> res = var + num
    
    """
    if   isinstance ( o , _RRV_    ) and not o.isConstant() : o = o.ve     () 
    elif hasattr    ( o , 'getVal' )                        : o = o.getVal ()
    #
    v = s.getVal() if s.isConstant() else s.ve()
    #
    return v + o

# ============================================================================
## Subtraction  of RooRealVar and ``number''
def _rrv_sub_ ( s , o ) :
    """Subtraction of RooRealVar and ``number''
    
    >>> var = ...
    >>> num = ...
    >>> res = var - num
    
    """
    if   isinstance ( o , _RRV_    ) and not o.isConstant() : o = o.ve     () 
    elif hasattr    ( o , 'getVal' )                        : o = o.getVal ()
    #
    v = s.getVal() if s.isConstant() else s.ve()
    #
    return v - o

# ============================================================================
## Multiplication of RooRealVar and ``number''
def _rrv_mul_ ( s , o ) :
    """Multiplication  of RooRealVar and ``number''
    
    >>> var = ...
    >>> num = ...
    >>> res = var * num
    
    """
    if   isinstance ( o , _RRV_    ) and not o.isConstant() : o = o.ve     () 
    elif hasattr    ( o , 'getVal' )                        : o = o.getVal ()
    #
    v = s.getVal() if s.isConstant() else s.ve()
    #
    return v * o

# ============================================================================
## Division of RooRealVar and ``number''
def _rrv_div_ ( s , o ) :
    """Division of RooRealVar and ``number''
    
    >>> var = ...
    >>> num = ...
    >>> res = var / num
    
    """
    if   isinstance ( o , _RRV_    ) and not o.isConstant() : o = o.ve     () 
    elif hasattr    ( o , 'getVal' )                        : o = o.getVal ()
    #
    v = s.getVal() if s.isConstant() else s.ve()
    #
    return v / o

# ============================================================================
## (right) Addition of RooRealVar and ``number''
def _rrv_radd_ ( s , o ) :
    """(right) Addition of RooRealVar and ``number''
    
    >>> var = ...
    >>> num = ...
    >>> res = num + var 
    
    """
    if   isinstance ( o , _RRV_    ) and not o.isConstant() : o = o.ve     () 
    elif hasattr    ( o , 'getVal' )                        : o = o.getVal ()
    #
    v = s.getVal() if s.isConstant() else s.ve()
    #
    return o + v 

# ============================================================================
## (right) subtraction  of RooRealVar and ``number''
def _rrv_rsub_ ( s , o ) :
    """(right) subtraction of RooRealVar and ``number''
    
    >>> var = ...
    >>> num = ...
    >>> res = num - var 
    
    """
    if   isinstance ( o , _RRV_    ) and not o.isConstant() : o = o.ve     () 
    elif hasattr    ( o , 'getVal' )                        : o = o.getVal ()
    #
    v = s.getVal() if s.isConstant() else s.ve()
    #
    return o - v 

# ============================================================================
## (right) multiplication of RooRealVar and ``number''
def _rrv_rmul_ ( s , o ) :
    """(right) Multiplication  of RooRealVar and ``number''
    
    >>> var = ...
    >>> num = ...
    >>> res = num * var     
    """
    if   isinstance ( o , _RRV_    ) and not o.isConstant() : o = o.ve     () 
    elif hasattr    ( o , 'getVal' )                        : o = o.getVal ()
    #
    v = s.getVal() if s.isConstant() else s.ve()
    #
    return o * v 

# ============================================================================
## (right) Division of RooRealVar and ``number''
def _rrv_rdiv_ ( s , o ) :
    """(right) Division of RooRealVar and ``number''
    
    >>> var = ...
    >>> num = ...
    >>> res = num / var     
    """
    if   isinstance ( o , _RRV_    ) and not o.isConstant() : o = o.ve     () 
    elif hasattr    ( o , 'getVal' )                        : o = o.getVal ()
    #
    v = s.getVal() if s.isConstant() else s.ve()
    #
    return o / v 

# ============================================================================
## pow of RooRealVar and ``number''
def _rrv_pow_ ( s , o ) :
    """pow of RooRealVar and ``number''
    
    >>> var = ...
    >>> num = ...
    >>> res = var ** num     
    """
    if   isinstance ( o , _RRV_    ) and not o.isConstant() : o = o.ve     () 
    elif hasattr    ( o , 'getVal' )                        : o = o.getVal ()
    #
    v = s.getVal() if s.isConstant() else s.ve()
    #
    return v**o  

# ============================================================================
## (right) pow of RooRealVar and ``number''
def _rrv_rpow_ ( s , o ) :
    """pow of RooRealVar and ``number''
    
    >>> var = ...
    >>> num = ...
    >>> res = num ** var 
    
    """
    if   isinstance ( o , _RRV_    ) and not o.isConstant() : o = o.ve     () 
    elif hasattr    ( o , 'getVal' )                        : o = o.getVal ()
    #
    v = s.getVal() if s.isConstant() else s.ve()
    #
    return o**v   


# ============================================================================
ROOT.RooRealVar . __add__   = _rrv_add_
ROOT.RooRealVar . __sub__   = _rrv_sub_
ROOT.RooRealVar . __div__   = _rrv_div_
ROOT.RooRealVar . __mul__   = _rrv_mul_
ROOT.RooRealVar . __pow__   = _rrv_pow_

ROOT.RooRealVar . __radd__  = _rrv_radd_
ROOT.RooRealVar . __rsub__  = _rrv_rsub_
ROOT.RooRealVar . __rdiv__  = _rrv_rdiv_
ROOT.RooRealVar . __rmul__  = _rrv_rmul_
ROOT.RooRealVar . __rpow__  = _rrv_rpow_


_new_methods_ += [
    ROOT.RooRealVar.__add__  , 
    ROOT.RooRealVar.__sub__  , 
    ROOT.RooRealVar.__div__  , 
    ROOT.RooRealVar.__mul__  , 
    ROOT.RooRealVar.__pow__  , 
    ROOT.RooRealVar.__radd__ , 
    ROOT.RooRealVar.__rsub__ , 
    ROOT.RooRealVar.__rdiv__ , 
    ROOT.RooRealVar.__rmul__ , 
    ROOT.RooRealVar.__rpow__ , 
    ]

# =============================================================================
## (compare RooRealVar and "number"
def _rrv_le_ ( s , o ) :
    """compare RooRealVal and ``number''
    
    >>> var = ...
    >>> num = ...
    >>> iv var <= num : print ' ok! '
    """
    return o >= s.getVal()

# ============================================================================
## (compare RooRealVar and "number"
def _rrv_lt_ ( s , o ) :
    """compare RooRealVal and ``number''
    
    >>> var = ...
    >>> num = ...
    >>> iv var < num : print ' ok! '
    """
    return o > s.getVal()

# ============================================================================
## (compare RooRealVar and "number"
def _rrv_ge_ ( s , o ) :
    """compare RooRealVal and ``number''
    
    >>> var = ...
    >>> num = ...
    >>> iv var >= num : print ' ok! '
    """
    return o <= s.getVal()

# ============================================================================
## (compare RooRealVar and "number"
def _rrv_gt_ ( s , o ) :
    """compare RooRealVal and ``number''
    
    >>> var = ...
    >>> num = ...
    >>> iv var > num : print ' ok! '
    """
    return o < s.getVal()

# ============================================================================
ROOT.RooRealVar . __lt__   = _rrv_lt_
ROOT.RooRealVar . __gt__   = _rrv_gt_
ROOT.RooRealVar . __le__   = _rrv_le_
ROOT.RooRealVar . __ge__   = _rrv_ge_

_new_methods_ += [
    ROOT.RooRealVar.__lt__  ,
    ROOT.RooRealVar.__gt__  ,
    ROOT.RooRealVar.__le__  ,
    ROOT.RooRealVar.__ge__  ,
    ]

# =============================================================================
## get min/max in one go 
#  @see RooRealVar
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-07-14
def _rrv_minmax_ ( s ) :
    """Get min/max in one go

    >>> var = ...
    >>> mn,mx = var.minmax()
    """
    return s.getMin(),s.getMax()

ROOT.RooRealVar   . minmax  = _rrv_minmax_

_new_methods_ += [
    ROOT.RooRealVar.minmax  ,
    ]

## # ==========================================
## _is_equal_ = cpp.LHCb.Math.equal_to_double 
## ## (compare RooRealVar and "number"
## def _rrv_eq_ ( s , o ) :
##     """
##     compare RooRealVal and ``number''
    
##     >>> var = ...
##     >>> num = ...
##     >>> iv var == num : print ' ok! '
##     """
##     return    _is_equal_ ( o , s.getVal() ) 

## ## (compare RooRealVar and "number"
## def _rrv_ne_ ( s , o ) :
##     """
##     compare RooRealVal and ``number''
    
##     >>> var = ...
##     >>> num = ...
##     >>> iv var != num : print ' ok! '
##     """
##     return not _is_equal_ ( o , s.getVal() ) 


## ROOT.RooRealVar . __eq__   = _rrv_eq_
## ROOT.RooRealVar . __ne__   = _rrv_ne_

# ============================================================================
## product of two PDFs 
def _pdf_mul_ ( pdf1 , pdf2 ) :
    """Easy contruct for the product of two PDFs:
    
    >>> pdf1 = ...
    >>> pdf2 = ...
    
    >>> product = pdf1 * pdf2 
    """
    return cpp.Ostap.Models.Product ( '%s*%s'             % ( pdf1.GetName  () ,
                                                              pdf2.GetName  () ) ,
                                      'Product: %s & %s ' % ( pdf1.GetTitle () ,
                                                              pdf2.GetTitle () ) ,
                                      pdf1 , pdf2 )

# ============================================================================
ROOT.RooAbsPdf . __mul__  = _pdf_mul_ 


_new_methods_ += [
    ROOT.RooAbsPdf.__mul__  , 
    ]

# ============================================================================
## "convert" name/expression into variable/formula
def var_from_name ( w , varset ) :
    """ Convert name/expression into variable/formula
    """
    w = w.strip() 
    if    0 < w.find('(') < what.find(')') :
        #print ' function ' , w 
        pass
    elif  0 < w.find('*')                  :
        #print ' multiply ' , w 
        pass
    elif  0 < w.find('/')                  :
        #print ' divide ' , w 
        pass
    elif  0 < w.find('+')                  :
        #print ' add  ' , w 
        pass
    elif  0 < w.find('-')                  :
        #print ' minus  ' , w 
        pass
    else :
        #print ' primitive ' , w 
        v = varset[w]
        return v
    ##
    
    vlst = ROOT.RooArgList()
    for s in varset : vlst.add ( s )
    #
    f = ROOT.RooFormulaVar( w , w , vlst )
    return f 


# =============================================================================
## Helper project method for RooDataSet
#
#  @code 
#    
#    >>> h1   = ROOT.TH1D(... )
#    >>> dataset.project ( h1.GetName() , 'm', 'chi2<10' ) ## project variable into histo
#    
#    >>> h1   = ROOT.TH1D(... )
#    >>> dataset.project ( h1           , 'm', 'chi2<10' ) ## use histo
#
#  @endcode
#
#  @see RooDataSet 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-07-06
def _ds_project_  ( dataset , histo , what , cuts = '' , *args ) :
    """Helper project method for RooDataSet
    
    >>> h1   = ROOT.TH1D(... )
    >>> dataset.project ( h1.GetName() , 'm', 'chi2<10' ) ## project variable into histo
    
    >>> h1   = ROOT.TH1D(... )
    >>> dataset.project ( h1           , 'm', 'chi2<10' ) ## use histo
    
    """
    if isinstance ( cuts , ROOT.TCut ) : cuts = str ( cuts ).strip()  
    if isinstance ( what , str       ) : what = what.strip()
    if isinstance ( cuts , str       ) : cuts = cuts.strip()
    
    ## native RooFit...  I have some suspicion that it does not work properly
    if isinstance ( what , ROOT.RooArgList ) and isinstance ( histo , ROOT.TH1 ) :
        histo.Reset() 
        return dataset.fillHistogram  ( histo , what , cuts , *args )
    
    ## delegate to TTree (only for non-weighted dataset with TTree-based storage type) 
    if not dataset.isWeighted() \
       and isinstance ( what , str ) \
       and isinstance ( cuts , str ) : 
        store = dataset.store()
        if store :
            tree = store.tree()
            if tree : return tree.project ( histo , what , cuts , *args ) 

            
    if   isinstance ( what , ROOT.RooFormulaVar ) : 
        return _ds_project_ ( dataset , histo , what.GetTitle () , cuts , *args )
    
    if   isinstance ( what , ROOT.RooAbsReal ) : 
        return _ds_project_ ( dataset , histo , what.GetName  () , cuts , *args ) 
    
    if isinstance ( what , str ) : 
        vars  = [ v.strip() for v in what.split(':') ]
        return _ds_project_ ( dataset , histo , vars , cuts , *args ) 
    
    if isinstance ( what , (tuple,list) ) :
        vars = []
        for w in what :
            if isinstance ( w , str ) : vars.append ( w.strip() )
            else                      : vars.append ( w ) 
        ##return _ds_project_ ( dataset , histo , vars , cuts , *args ) 

    if isinstance ( what , ROOT.RooArgList ) :
        vars  = [ w for w in what ]
        cuts0 = cuts 
        if ''   == cuts : cuts0 = 0
        elif isinstance ( cuts , str ) :
            cuts0 = ROOT.RooFormulaVar( cuts , cuts , dataset.varlist() )
        return _ds_project_ ( dataset , histo , vars , cuts0 , *args ) 
            
    if isinstance ( histo , str ) :
    
        obj = ROOT.gROOT     .FindObject    ( histo )
        if instance ( obj  , ROOT.TH1 ) :
            return _ds_project_ ( dataset , obj , what , cuts , *args )
        obj = ROOT.gROOT     .FindObjectAny ( histo )
        if instance ( obj  , ROOT.TH1 ) :
            return _ds_project_ ( dataset , obj , what , cuts , *args )
        obj = ROOT.gDirectory.FindObject    ( histo )
        if instance ( obj  , ROOT.TH1 ) :
            return _ds_project_ ( dataset , obj , what , cuts , *args )
        obj = ROOT.gDirectory.FindObjectAny ( histo )
        if instance ( obj  , ROOT.TH1 ) :
            return _ds_project_ ( dataset , obj , what , cuts , *args )

    if  1 <= len(what) and isinstance ( what[0] , ROOT.RooAbsReal ) and isinstance ( cuts , str ) : 
        if '' == cuts : cuts0 = 0 
        elif isinstance ( cuts , str ) :
            cuts0 = ROOT.RooFormulaVar( cuts , cuts , dataset.varlist() )
        return _ds_project_ ( dataset , histo , what , cuts0 , *args )

    if   isinstance ( histo , ROOT.TH3 ) and 3 == len(what)  :
        return cpp.Ostap.HistoProject.project3 ( dataset ,
                                                histo   , 
                                                what[0] ,
                                                what[1] ,
                                                what[2] , cuts , *args) 
    elif isinstance ( histo , ROOT.TH2 ) and 2 == len(what)  :
        return cpp.Ostap.HistoProject.project2 ( dataset ,
                                                 histo   , 
                                                 what[0] ,
                                                 what[1] , cuts , *args )
    elif isinstance ( histo , ROOT.TH1 ) and 1 == len(what)  :
        return cpp.Ostap.HistoProject.project  ( dataset ,
                                                 histo   , 
                                                 what[0] , cuts , *args )
    
    raise AttributeError ( 'DataSet::project, invalid case' )

# =============================================================================
## Helper draw method for RooDataSet
#
#  @code 
#    
#    >>> dataset.draw ( 'm', 'chi2<10' ) ## use histo
#
#  @endcode
#
#  @see RooDataSet 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-07-06
def _ds_draw_ ( dataset , what , cuts = '' , opts = '' , *args ) :
    """Helper draw method for drawing of RooDataSet
    >>> dataset.draw ( 'm', 'chi2<10'                 )
    ## cuts & weight 
    >>> dataset.draw ( 'm', '(chi2<10)*weight'        )
    ## use drawing options 
    >>> dataset.draw ( 'm', '(chi2<10)*weight' , 'e1' )
    ## start form event #1000 
    >>> dataset.draw ( 'm', '(chi2<10)*weight' , 'e1' , 1000 ) 
    ## for event in range 1000< i <10000
    >>> dataset.draw ( 'm', '(chi2<10)*weight' , 'e1' , 1000 , 100000 )
    """
    if isinstance ( cuts , ROOT.TCut ) : cuts = str ( cuts ).strip()  
    if isinstance ( what , str       ) : what = what.strip()
    if isinstance ( cuts , str       ) : cuts = cuts.strip()
    if isinstance ( opts , str       ) : opts = opts.strip()

    ## delegate to TTree for non-weighted datasets with TTree-based storage type 
    if not dataset.isWeighted() \
       and isinstance ( what , str ) \
       and isinstance ( cuts , str ) \
       and isinstance ( opts , str ) : 
        store = dataset.store()
        if store : 
            tree = store.tree()
            if tree : return tree.Draw( what , cuts , opts  , *args )
        
    if   isinstance ( what , str ) : 
        vars  = [ v.strip() for v in what.split(':') ]
        return _ds_draw_ ( dataset , vars , cuts , opts , *args ) 
    
    if   isinstance ( what , ROOT.RooFormulaVar ) : 
        return _ds_draw_ ( dataset , what.GetTitle () , cuts , opts , *args )
    
    if   isinstance ( what , ROOT.RooAbsReal ) : 
        return _ds_draw_ ( dataset , what.GetName  () , cuts , opts , *args ) 
    
    if not 1 <= len ( what ) <= 3 :
        raise AttributeError ( 'DataSet::draw, invalid length %s' % what  )
    
    if 1 == len ( what )  :
        w1        = what[0] 
        mn1 , mx1 = _ds_var_minmax_  ( dataset , w1 , cuts )
        histo = ROOT.TH1F ( hID() , w1 , 200 , mn1 , mx1 )  ; histo.Sumw2()
        _ds_project_ ( dataset , histo , what , cuts , *args  )
        histo.Draw( opts )
        return histo

    if 2 == len ( what )  :
        w1        = what[0] 
        mn1 , mx1 = _ds_var_minmax_  ( dataset , w1 , cuts )
        w2        = what[1] 
        mn2 , mx2 = _ds_var_minmax_  ( dataset , w2 , cuts )
        histo = ROOT.TH2F ( hID() , "%s:%s" % ( w1 , w2 ) ,
                            50 , mn1 , mx1 ,
                            50 , mn2 , mx2 )  ; histo.Sumw2()
        _ds_project_ ( dataset , histo , what , cuts , *args  )
        histo.Draw( opts )
        return histo

    if 3 == len ( what )  :
        w1        = what[0] 
        mn1 , mx1 = _ds_var_minmax_  ( dataset , w1 , cuts )
        w2        = what[1] 
        mn2 , mx2 = _ds_var_minmax_  ( dataset , w2 , cuts )
        w3        = what[2] 
        mn3 , mx3 = _ds_var_minmax_  ( dataset , w3 , cuts )
        histo = ROOT.TH3F ( hID() , "%s:%s:%s" % ( w1 , w2 , w3 ) ,
                            20 , mn1 , mx1 ,
                            20 , mn2 , mx2 ,
                            20 , mn2 , mx2 )  ; histo.Sumw2()
        _ds_project_ ( dataset , histo , what , cuts , *args  )
        histo.Draw( opts )
        return histo

    raise AttributeError ( 'DataSet::draw, invalid case' )

# =============================================================================
## get the attibute for RooDataSet
def _ds_getattr_ ( dataset , aname ) :
    """Get the attibute from RooDataSet 

    >>> dset = ...
    >>> print dset.pt
    
    """
    _vars = dataset.get()
    return getattr ( _vars , aname )  

# =============================================================================
## Get min/max for the certain variable in dataset
#  @code  
#  data = ...
#  mn,mx = data.vminmax('pt')
#  mn,mx = data.vminmax('pt','y>3')
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-09-19
def _ds_var_minmax_ ( dataset , var , cuts = '' , delta = 0.0 )  :
    """Get min/max for the certain variable in dataset
    >>> data = ...
    >>> mn,mx = data.vminmax('pt')
    >>> mn,mx = data.vminmax('pt','y>3')
    """
    if isinstance  (  var , ROOT.RooAbsReal ) : var = var.GetName() 
    if cuts : s = dataset.statVar ( var , cuts )
    else    : s = dataset.statVar ( var )
    mn,mx = s.minmax()
    if mn < mn and 0.0 < delta :
        dx   = delta * 1.0 * ( mx - mn )  
        mx  += dx   
        mn  -= dx   
    return mn , mx


ROOT.RooDataSet .vminmax  = _ds_var_minmax_ 

_new_methods_ += [
    ROOT.RooDataSet .vminmax ,
    ]
# =============================================================================
## print method for RooDataSet
#  @code
#
#   >>> print dataset
#
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-07-06
def _ds_print_ ( dataset , opts = 'v' ) :
    """Helper print method:    
    >>> print dataset 
    """
    #
    dataset.Print( opts )
    #
    return dataset.GetName() 

# =============================================================================
ROOT.RooDataSet.draw        = _ds_draw_
ROOT.RooDataSet.project     = _ds_project_
ROOT.RooDataSet.__repr__    = _ds_print_
ROOT.RooDataSet.__getattr__ = _ds_getattr_

ROOT.RooDataHist.__repr__   = _ds_print_
ROOT.RooDataHist.__len__    = lambda s : s.numEntries() 


_new_methods_ += [
    ROOT.RooDataSet .draw         ,
    ROOT.RooDataSet .project      ,
    ROOT.RooDataSet .__repr__     ,
    ROOT.RooDataSet .__getattr__  ,
    ROOT.RooDataHist.__repr__     ,
    ROOT.RooDataHist.__len__      ,
    ]

# =============================================================================
## add variable to dataset 
def _rds_addVar_ ( dataset , vname , formula ) : 
    """Add/calculate variable to RooDataSet

    >>> dataset.addVar ( 'ratio' , 'pt/pz' )
    """
    vlst     = ROOT.RooArgList()
    vset     = dataset.get()
    for   v     in vset : vlst.add ( v )
    #
    vcol     = ROOT.RooFormulaVar ( vname , formula , formula , vlst )
    dataset.addColumn ( vcol )
    #
    return dataset 

# =============================================================================
ROOT.RooDataSet.addVar = _rds_addVar_

_new_methods_ += [
    ROOT.RooDataSet .addVar       ,
    ]

# =============================================================================
## make weighted data set from unweighted dataset
#  @code
#  >>> dataset = ...
#  >>> wdata   = dataset.makeWeighted ( 'S_sw' ) 
#  @endcode
#  @param wvarname name of weighting variable
#  @param varset   variables to be used in new dataset
#  @param cuts     optional cuts to be applied 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-07-06
def _rds_makeWeighted_ ( dataset , wvarname , varset = None , cuts = '' , vname = '' ) :
    """Make weighted data set from unweighted dataset
    
    >>> dataset = ...
    >>> wdata   = dataset.makeWeighted ( 'S_sw' )    
    """
    if dataset.isWeighted () : 
        logger.warning ("Dataset '%s/%s' is already weighted!" % ( dataset.GetName  () ,
                                                                   dataset.GetTitle () ) ) 

    ##
    formula =  0 <= wvarname.find ( '(' ) and wvarname.find( '(' ) < wvarname.find ( ')' )
    formula = formula or 0 <  wvarname.find ( '*' ) 
    formula = formula or 0 <  wvarname.find ( '/' )     
    formula = formula or 0 <= wvarname.find ( '+' ) 
    formula = formula or 0 <= wvarname.find ( '-' )     
    formula = formula or 0 <  wvarname.find ( '&' )     
    formula = formula or 0 <  wvarname.find ( '|' )     

    if formula :
        wname    = 'W' or vname 
        while wname in dataset : wname += 'W'
        dataset.addVar ( wname , wvarname ) 
        wvarname = wname  
        
    if not varset :
        varset = dataset.get()  
   
    ## make weighted dataset 
    return ROOT.RooDataSet ( dsID()             ,
                             dataset.GetTitle() ,
                             dataset            ,
                             varset             , 
                             cuts               ,
                             wvarname           )

# =============================================================================
ROOT.RooDataSet.makeWeighted = _rds_makeWeighted_

_new_methods_ += [
    ROOT.RooDataSet .makeWeighted ,
    ]

RAD = ROOT.RooAbsData
# =============================================================================
## change the default storage for RooDataSet 
def setStorage ( new_type = RAD.Tree ) :
    """ Redefine the default storage 
    """
    if not new_type in ( RAD.Tree , RAD.Vector ) :
        logger.error ('RooAbsData: Invalid storage type %s, replace with Tree ' % new_type )
        new_type = RAD.Tree
        
    if RAD.getDefaultStorageType() != new_type :
        logger.info  ( 'RooAbsData: DEFINE default storage type to be %d' % new_type ) 
        RAD.setDefaultStorageType ( new_type  ) 

    the_type = RAD.getDefaultStorageType()
    if   RAD.Tree   == the_type : logger.debug ( 'RooAbsData: Default storage type is Tree'   )
    elif RAD.Vector == the_type : logger.debug ( 'RooAbsData: Default storage type is Vector' )
    else : logger.debug ( 'RooAbsData: Default storage type is %s' % the_type  )

# =============================================================================
## make easy print for RooPrintable 
def _rp_print_ ( o , opts = 'v' ) :
    """Make easy print for RooPrintable
    >>> o = ...
    >>> print o 
    """
    return cpp.Ostap.Utils.print_printable ( o , opts )

ROOT.RooPrintable.print_printable = _rp_print_ 
ROOT.RooPrintable.__str__         = _rp_print_
ROOT.RooPrintable.__repr__        = _rp_print_

_new_methods_ += [
    ROOT.RooPrintable.print_printable ,
    ROOT.RooPrintable.__str__         ,
    ROOT.RooPrintable.__repr__        ,
    ]

# =============================================================================
## @class SETVAR
#  Simple context manager to preserve current value for RooAbsVar
#  @code
#  var = ...
#  var.setVal(1) 
#  print '1) value %s ' % var.getVal() 
#  with SETVAR(var) :
#        print '2) value %s ' % var.getVal() 
#        var.setVal(10)
#        print '3) value %s ' % var.getVal() 
#  print '4) value %s ' % var.getVal() 
#  @endcode
class SETVAR(object):
    """ Simple context manager to preserve current value for RooAbsVar
    >>> var = ...
    >>> var.setVal(1) 
    >>> print '1) value %s ' % var.getVal() 
    >>> with SETVAR(var) :
    ...    print '2) value %s ' % var.getVal() 
    ...    var.setVal(10)
    ...    print '3) value %s ' % var.getVal() 
    >>> print '4) value %s ' % var.getVal() 
    """
    def __init__  ( self , xvar ) :
        self.xvar = xvar
    def __enter__ ( self        ) :
        self._old = float ( self.xvar.getVal() ) 
        return self 
    def __exit__  ( self , *_   ) :
        self.xvar.setVal  ( self._old ) 

# =============================================================================
## @class PDF_fun
#  Helper class to wrap PDF as 'function'
#  can be helpful for some pure math-operations
#  @code
#  pdf,var = ....
#  fun     = PDF( fun , var , xmin=0 , xmax=1 )
#  from ostap.stats.moments import mean, mode, median, CL   
#  print 'MEAN    : %s' % mean    ( fun , 0 , 1 )
#  print 'MODE    : %s' % mode    ( fun , 0 , 1 )
#  print 'MEDIAN  : %s' % median  ( fun , 0 , 1 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-03-29
class PDF_fun(object):
    """ Helper class to wrap PDF as 'function'
    >>> pdf,var = ....
    >>> fun     = PDF( pdf , var , xmin=0 , xmax=1 )
    >>> print fun(0.1),fun(0.5) 
    >>> from LHCbMath.deriv import mean, mode, median  
    >>> print 'MEAN    : %s' % mean    ( fun , 0 , 1 )
    >>> print 'MODE    : %s' % mode    ( fun , 0 , 1 )
    >>> print 'MEDIAN  : %s' % median  ( fun , 0 , 1 )
    """
    ##
    def __init__ ( self , pdf , xvar , xmin = None , xmax = None ) :
        
        self.pdf     = pdf

        ## ostap stuff: 
        if not isinstance ( pdf , ROOT.RooAbsPdf ) :
            if hasattr ( self.pdf , 'pdf' ) :
                self.pdf_ = pdf 
                self.pdf  = pdf.pdf

        self.xvar    = xvar

        self._xmin   = None 
        self._xmax   = None
        
        if not xmin is None : self._xmin = xmin 
        if not xmax is None : self._xmax = xmax

        if hasattr ( xvar , 'getMin' ) :
            if self._xmin is None : self._xmin = xvar.getMin()
            else                  : self._xmin = max ( self._xmin , xvar.getMin() )
            
        if hasattr ( xvar , 'getMax' ) :
            if self._xmax is None : self._xmax = xvar.getMax()
            else                  : self._xmax = min ( self._xmax , xvar.getMax() )
            
        if self._xmin is None :
            raise AttributeError, "xmin can't be deduced from  input arguments"
        if self._xmax is None :
            raise AttributeError, "xmax can't be deduced from  input arguments"
        
        if self._xmin > self._xmax :
            self._xmin , self._xmax = self._xmax , self._xmin
            
    def xmin     ( self ) : return self._xmin
    def xmax     ( self ) : return self._xmax
    
    ## the main method 
    def __call__ ( self , x , pars = [] ) :

        ## for ROOT.TF1
        if   hasattr ( x , '__len__' ) and hasattr( x , 'getitem' ) and not len( x )   : x = x[0]
        elif hasattr ( x , '__len__' ) and hasattr( x , 'getitem' ) and 0 != x.size()  : x = x[0]
        
        ## try to be efficient 
        if not self._xmin <= x <= self._xmax : return 0 
        
        with SETVAR( self.xvar ) :
            self.xvar.setVal ( x )
            return self.pdf.getVal()
    
# =============================================================================
_decorated_classes_ = (
    ROOT.RooAbsData    ,
    ROOT.RooFitResult  ,
    ROOT.RooRealVar    ,
    ROOT.RooConstVar   ,
    ROOT.RooFormulaVar ,
    ROOT.RooAbsReal    ,
    ROOT.RooDataSet    ,
    ROOT.RooDataHist   ,
    ROOT.RooPrintable  ,
    )
_new_methods_ = tuple ( _new_methods_ ) 
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
