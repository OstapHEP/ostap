#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/variables.py
#  Module with decoration of some RooFit variables for efficient usage in python
#  @see RooAbsReal
#  @see RooRealVar
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Module with decoration of some RooFit variables for efficient usage in python
- see RooAbsReal
- see RooRealVar
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'SETVAR'          , ## context manager to preserve the current value for RooRealVar
    'FIXVAR'          , ## context manager to fix/unfix variables
    'KeepBinning'     , ## context manager to preserve the binning scheme
    ## converters
    'total_ratio'     , ## converter : A,B ->  (T,R) == ( A+B , A/B     )
    'total_ratio'     , ## converter : A,B ->  (T,F) == ( A+B , A/(A+B) ) 
    'two_yields'      , ## converter : T,F ->  (A,B) == ( R*F , T*(1-F) )
    'depends_on'      , ## Is this "RooFit" function depends on the variable?
    'binning'         , ## create RooBinning object
    ##
    'var_sum'        , ## sum                          for RooAbsReal objects           
    'var_mul'        , ## product                      for RooAbsReal objects           
    'var_sub'        , ## subtraction                  for RooAbsReal objects           
    'var_div'        , ## division                     for RooAbsReal objects           
    'var_fraction'   , ## fraction                     for RooAbsReal objects           
    'var_asymmetry'  , ## asymmetry                    for RooAbsReal objects
    'var_pow'        , ## pow                 function for RooAbsReal objects           
    'var_abs'        , ## absolutevalue       function for RooAbsReal objects           
    'var_exp'        , ## exponent            function for RooAbsReal objects           
    'var_log'        , ## logarithm           function for RooAbsReal objects           
    'var_log10'      , ## logarithm           function for RooAbsReal objects           
    'var_erf'        , ## error               function for RooAbsReal objects           
    'var_sin'        , ## sine                function for RooAbsReal objects           
    'var_cos'        , ## cosine              function for RooAbsReal objects           
    'var_tan'        , ## tangent             function for RooAbsReal objects           
    'var_tanh'       , ## hyperbolic tangent  function for RooAbsReal objects           
    'var_sech'       , ## hyperbolic secant   function for RooAbsReal objects           
    'var_atan2'      , ## inverse tangent     function for RooAbsReal objects
    'var_bessel_J'   , ## regular   Bessel function 
    'var_bessel_Y'   , ## irregular Bessel function 
    'var_bessel_I'   , ## modified  Bessel function 
    'var_bessel_K'   , ## modified  Bessel function 
    'var_min'        , ## minimal             function for RooAbsReal objects           
    'var_max'        , ## minimal             function for RooAbsReal objects           
    'var_gamma'      , ## gamma               function for RooAbsReal objects           
    'var_lgamma'     , ## logarithm of gamma  function for RooAbsReal objects           
    'var_igamma'     , ## 1/gamma             function for RooAbsReal objects
    ##
    'scale_var'      , ## var_mul
    'add_var'        , ## var_sum
    'sum_var'        , ## var_sum
    'ratio_var'      , ## var_div
    'fraction_var'   , ## var_fraction
    'asymmetry_var'  , ## var_asymmetry
    ##
    ) 
# =============================================================================
from   builtins                 import range
from   ostap.math.base          import doubles, iszero, isequal 

from   ostap.core.core          import VE, hID, Ostap
from   ostap.math.reduce        import root_factory 
from   ostap.core.meta_info     import root_info 
from   ostap.core.ostap_types   import ( num_types      , list_types   ,
                                         integer_types  , string_types ,
                                         dictlike_types )
import ostap.math.math_ve       as       mve 
import ostap.fitting.rooreduce 
import ROOT, random, array, ctypes
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.variables' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug( 'Some useful decorations for RooFit variables')
# =============================================================================
## is value equal to 1?
isone = lambda x : isequal ( float ( x ) , 1 )
# =============================================================================
_new_methods_ = []
# =============================================================================
## Factory for deserialization of generic objects
#  @attention it stores the constructor parameters as local attributes
def root_store_factory ( klass , *params ) :
    """Factory for deserialization of generic object
    - attention: it stores the constructor parameters as local attributes
    """
    ## create the objects 
    obj = root_factory ( klass , *params )
    ## keep argumets with the newly created obnject  
    obj.__store = params    ## Attention - keep argumetns with newly crfeated object!
    return obj 

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
        if var.isConstant() : return var.ve() 
        var.setConstant( True )
        return var.ve()
    
    if hasattr ( value , 'value' ) : value = value.value()
    #
    var.setVal      ( value )
    if not var.isConstant() : var.setConstant ( True  )
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
    if var.isConstant() : var.setConstant ( False )
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

# ==============================================================================
## check if the given value is in the range of RooRealVar
#  @code 
#  mass_range = ...
#  if v in mass_range : ...
#  @endcode 
def _rrv_contains_ ( var , value ) :
    """check if the given value is in the range of RooRealVar
    >>> mass_range = ...
    >>> if v in mass_range : ... 
    """
    if var.hasMin() and value < var.getMin() : return False 
    if var.hasMax() and value > var.getMax() : return False
    return True 
    
# =============================================================================
## decorate RooRealVar:
ROOT.RooRealVar     . as_VE           = _rrv_ve_ 
ROOT.RooRealVar     . asVE            = _rrv_ve_ 
ROOT.RooRealVar     . ve              = _rrv_ve_
ROOT.RooRealVar     . fix             = _fix_par_
ROOT.RooRealVar     . Fix             = _fix_par_
ROOT.RooRealVar     . release         = _rel_par_
ROOT.RooRealVar     . Release         = _rel_par_
## convert to float 
ROOT.RooRealVar     . __float__       = lambda s : s.getVal()
## print it in more suitable form 
ROOT.RooRealVar     . __repr__        = lambda s : "'%s' : %s " % ( s.GetName() , s.ve() )
ROOT.RooRealVar     . __str__         = lambda s : "'%s' : %s " % ( s.GetName() , s.ve() )


ROOT.RooConstVar    . as_VE          = lambda s : VE ( s.getVal() , 0 )
ROOT.RooFormulaVar  . as_VE          = lambda s : VE ( s.getVal() , 0 )
ROOT.RooConstVar    . asVE           = lambda s : VE ( s.getVal() , 0 )
ROOT.RooFormulaVar  . asVE           = lambda s : VE ( s.getVal() , 0 )

ROOT.RooRealVar     . __float__      = lambda s : s.getVal()
ROOT.RooConstVar    . __float__      = lambda s : s.getVal()
ROOT.RooAbsReal     . __float__      = lambda s : s.getVal() ## NB!!!


ROOT.RooAbsReal       .__contains__ = lambda s,v : False ## ??? do we need it???
ROOT.RooAbsRealLValue .__contains__ = _rrv_contains_ 

# =====================================================================

ROOT.RooAbsReal. minmax  = lambda s : ()
ROOT.RooAbsReal.xminmax  = lambda s : ()
ROOT.RooAbsRealLValue  . xmin            = lambda s : s.getMin()
ROOT.RooAbsRealLValue  . xmax            = lambda s : s.getMax()
ROOT.RooAbsRealLValue  . minmax          = lambda s : ( s.xmin() , s.xmax() ) 
ROOT.RooAbsRealLValue  .xminmax          = lambda s : ( s.xmin() , s.xmax() ) 


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
    ROOT.RooAbsRealLValue .__contains__ , 
    ROOT.RooRealVar   . xmin      ,
    ROOT.RooRealVar   . xmax      ,
    ROOT.RooRealVar   . minmax    ,
    #
    ROOT.RooConstVar    .as_VE    ,
    ROOT.RooFormulaVar  .as_VE    ,
    ROOT.RooConstVar    .asVE     ,
    ROOT.RooFormulaVar  .asVE     ,
    #
    ROOT.RooAbsReal        . minmax  , 
    ROOT.RooAbsReal        .xminmax  , 
    ROOT.RooAbsRealLValue  . xmin    ,
    ROOT.RooAbsRealLValue  . xmax    ,
    ROOT.RooAbsRealLValue  . minmax  ,
    ROOT.RooAbsRealLValue  .xminmax  ,    
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

# =============================================================================
## Math operations 
# =============================================================================
val_types = num_types + ( VE , ROOT.RooAbsReal )

# ============================================================================
## Addition of RooRealVar and "the number"
def _rrv_add_ ( s , o ) :
    """Addition of RooRealVar and 'the number' 

    >>> var = ...
    >>> num = ...
    >>> res = var + num
    
    """
    if not isinstance ( o , val_types ) : return NotImplemented
    
    if   isinstance ( o , _RRV_    ) and not o.isConstant() : o = o.ve     () 
    elif hasattr    ( o , 'getVal' )                        : o = o.getVal ()
    #
    v = s.getVal() if s.isConstant() else s.ve()
    #
    return v + o

# ============================================================================
## Subtraction  of RooRealVar and "the number"
def _rrv_sub_ ( s , o ) :
    """Subtraction of RooRealVar and 'the 'number'
    
    >>> var = ...
    >>> num = ...
    >>> res = var - num
    
    """    
    if not isinstance ( o , val_types ) : return NotImplemented

    if   isinstance ( o , _RRV_    ) and not o.isConstant() : o = o.ve     () 
    elif hasattr    ( o , 'getVal' )                        : o = o.getVal ()
    #
    v = s.getVal() if s.isConstant() else s.ve()
    #
    return v - o

# ============================================================================
## Multiplication of RooRealVar and "the number"
def _rrv_mul_ ( s , o ) :
    """Multiplication  of RooRealVar and 'the number' 
    
    >>> var = ...
    >>> num = ...
    >>> res = var * num
    
    """    
    if not isinstance ( o , val_types ) : return NotImplemented
    
    if   isinstance ( o , _RRV_    ) and not o.isConstant() : o = o.ve     () 
    elif hasattr    ( o , 'getVal' )                        : o = o.getVal ()
    #
    v = s.getVal() if s.isConstant() else s.ve()
    #
    return v * o

# ============================================================================
## Division of RooRealVar and "the number"
def _rrv_div_ ( s , o ) :
    """Division of RooRealVar and 'the number'
    
    >>> var = ...
    >>> num = ...
    >>> res = var / num
    
    """
    if not isinstance ( o , val_types ) : return NotImplemented

    if   isinstance ( o , _RRV_    ) and not o.isConstant() : o = o.ve     () 
    elif hasattr    ( o , 'getVal' )                        : o = o.getVal ()
    #
    v = s.getVal() if s.isConstant() else s.ve()
    #
    return v / o

# ============================================================================
## (right) Addition of RooRealVar and "the number"
def _rrv_radd_ ( s , o ) :
    """(right) Addition of RooRealVar and 'the number'
    
    >>> var = ...
    >>> num = ...
    >>> res = num + var 
    
    """
    if not isinstance ( o , val_types ) : return NotImplemented

    if   isinstance ( o , _RRV_    ) and not o.isConstant() : o = o.ve     () 
    elif hasattr    ( o , 'getVal' )                        : o = o.getVal ()
    #
    v = s.getVal() if s.isConstant() else s.ve()
    #
    return o + v 

# ============================================================================
## (right) subtraction  of RooRealVar and "number"
def _rrv_rsub_ ( s , o ) :
    """(right) subtraction of RooRealVar and 'the number' 
    
    >>> var = ...
    >>> num = ...
    >>> res = num - var 
    
    """
    if not isinstance ( o , val_types ) : return NotImplemented

    if   isinstance ( o , _RRV_    ) and not o.isConstant() : o = o.ve     () 
    elif hasattr    ( o , 'getVal' )                        : o = o.getVal ()
    #
    v = s.getVal() if s.isConstant() else s.ve()
    #
    return o - v 

# ============================================================================
## (right) multiplication of RooRealVar and "the number"
def _rrv_rmul_ ( s , o ) :
    """(right) Multiplication  of RooRealVar and 'the number'
    
    >>> var = ...
    >>> num = ...
    >>> res = num * var     
    """
    if not isinstance ( o , val_types ) : return NotImplemented

    if   isinstance ( o , _RRV_    ) and not o.isConstant() : o = o.ve     () 
    elif hasattr    ( o , 'getVal' )                        : o = o.getVal ()
    #
    v = s.getVal() if s.isConstant() else s.ve()
    #
    return o * v 

# ============================================================================
## (right) Division of RooRealVar and "the number"
def _rrv_rdiv_ ( s , o ) :
    """(right) Division of RooRealVar and 'the number'
    
    >>> var = ...
    >>> num = ...
    >>> res = num / var     
    """
    if not isinstance ( o , val_types ) : return NotImplemented

    if   isinstance ( o , _RRV_    ) and not o.isConstant() : o = o.ve     () 
    elif hasattr    ( o , 'getVal' )                        : o = o.getVal ()
    #
    v = s.getVal() if s.isConstant() else s.ve()
    #
    return o / v 

# ============================================================================
## pow of RooRealVar and "the number"
def _rrv_pow_ ( s , o ) :
    """pow of RooRealVar and 'the number'
    
    >>> var = ...
    >>> num = ...
    >>> res = var ** num     
    """
    if not isinstance ( o , val_types ) : return NotImplemented

    if   isinstance ( o , _RRV_    ) and not o.isConstant() : o = o.ve     () 
    elif hasattr    ( o , 'getVal' )                        : o = o.getVal ()
    #
    v = s.getVal() if s.isConstant() else s.ve()
    #
    return v**o  

# ============================================================================
## (right) pow of RooRealVar and "the number"
def _rrv_rpow_ ( s , o ) :
    """pow of RooRealVar and 'the number'
    
    >>> var = ...
    >>> num = ...
    >>> res = num ** var 
    
    """
    if not isinstance ( o , val_types ) : return NotImplemented

    if   isinstance ( o , _RRV_    ) and not o.isConstant() : o = o.ve     () 
    elif hasattr    ( o , 'getVal' )                        : o = o.getVal ()
    #
    v = s.getVal() if s.isConstant() else s.ve()
    #
    return o**v   

# ============================================================================
def _rrv_iadd_ ( s , o ) :
    s.setVal ( s.getVal() + float ( o ) )
    return s
def _rrv_imul_ ( s , o ) :
    s.setVal ( s.getVal() * float ( o ) )
    return s
def _rrv_isub_ ( s , o ) :
    s.setVal ( s.getVal() - float ( o ) )
    return s
def _rrv_idiv_ ( s , o ) :
    s.setVal ( s.getVal() / float ( o ) )
    return s

# ============================================================================
ROOT.RooRealVar . __add__     = _rrv_add_
ROOT.RooRealVar . __sub__     = _rrv_sub_
ROOT.RooRealVar . __div__     = _rrv_div_
ROOT.RooRealVar . __mul__     = _rrv_mul_
ROOT.RooRealVar . __pow__     = _rrv_pow_

ROOT.RooRealVar . __radd__    = _rrv_radd_
ROOT.RooRealVar . __rsub__    = _rrv_rsub_
ROOT.RooRealVar . __rdiv__    = _rrv_rdiv_
ROOT.RooRealVar . __rmul__    = _rrv_rmul_
ROOT.RooRealVar . __rpow__    = _rrv_rpow_

ROOT.RooRealVar . __iadd__    = _rrv_iadd_
ROOT.RooRealVar . __isub__    = _rrv_isub_
ROOT.RooRealVar . __idiv__    = _rrv_idiv_

ROOT.RooRealVar . __truediv__ = ROOT.RooRealVar . __div__
ROOT.RooRealVar .__rtruediv__ = ROOT.RooRealVar .__rdiv__
ROOT.RooRealVar .__itruediv__ = ROOT.RooRealVar .__idiv__


_new_methods_ += [
    ROOT.RooRealVar.__add__   , 
    ROOT.RooRealVar.__sub__   , 
    ROOT.RooRealVar.__div__   , 
    ROOT.RooRealVar.__mul__   , 
    ROOT.RooRealVar.__pow__   , 
    ROOT.RooRealVar.__radd__  , 
    ROOT.RooRealVar.__rsub__  , 
    ROOT.RooRealVar.__rdiv__  , 
    ROOT.RooRealVar.__rmul__  , 
    ROOT.RooRealVar.__rpow__  ,
    #
    ROOT.RooRealVar.__iadd__  , 
    ROOT.RooRealVar.__isub__  , 
    ROOT.RooRealVar.__idiv__  ,
    #
    ROOT.RooRealVar . __truediv__ ,
    ROOT.RooRealVar .__rtruediv__ , 
    ROOT.RooRealVar .__itruediv__ ,
    ]

# =============================================================================
## (compare RooRealVar and "the number"
def _rrv_le_ ( s , o ) :
    """compare RooRealVal and 'the number'
    
    >>> var = ...
    >>> num = ...
    >>> iv var <= num : print ( ' ok! ') 
    """
    return o >= s.getVal()

# ============================================================================
## (compare RooRealVar and "the number"
def _rrv_lt_ ( s , o ) :
    """compare RooRealVal and 'the number'
    
    >>> var = ...
    >>> num = ...
    >>> iv var < num : print ( ' ok! ' ) 
    """
    return o > s.getVal()

# ============================================================================
## (compare RooRealVar and "the number"
def _rrv_ge_ ( s , o ) :
    """compare RooRealVal and 'the number'
    
    >>> var = ...
    >>> num = ...
    >>> iv var >= num : print ( ' ok! ' )
    """
    return o <= s.getVal()

# ============================================================================
## (compare RooRealVar and "the number"
def _rrv_gt_ ( s , o ) :
    """compare RooRealVal and 'the number'
    
    >>> var = ...
    >>> num = ...
    >>> iv var > num : print ( ' ok! ') 
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
## Properties 
# =============================================================================

# =============================================================================
def _rav_getval_  ( self ) :
    """Get the value, associated with the variable
    >>> var = ...
    >>> print ( var.value )
    """
    return self.getVal()

# =============================================================================
def _rav_getvale_ ( self ) :
    """Get the value(and the error), associated with the variable
    >>> var = ...
    >>> print  ( var.value )
    """
    v = self.getVal()
    e = self.getError() 
    return VE ( v , e * e ) if e > 0 else v

# =============================================================================
def _rav_setval_  ( self , value ) :
    """Assign the value for the variable 
    >>> var = ...
    >>> var.value = 10 
    """
    value = float ( value )
    self.setVal ( value ) 
    return self.getVal()

# =============================================================================
def _rav_setvalc_  ( self , value ) :
    """Assign the valeu for the variable 
    >>> var = ...
    >>> var.value = 10 
    """
    value = float ( value )
    mn , mx  = self.getMin(), self.getMax() 
    if not mn <= value <= mx :
        logger.warning('Value %s is out the range [%s,%s]' %  ( value  , mn , mx ) ) 
    self.setVal ( value ) 
    return self.getVal()

# =============================================================================
def _rav_geterr_  ( self ) :
    """Get the error
    >>> var = ...
    >>> print(var.error)
    """
    return self.getError()

# =============================================================================
def _rav_seterr_  ( self , value ) :
    """Set the error
    >>> var = ...
    >>> var.error = 10 
    """
    value = float ( value )
    if not 0<= value :
        logger.warning('Error %s is not non-negative' % value  ) 
    self.setError ( value )
    return self.getError()

# =============================================================================
## decorate classes 
for t in ( ROOT.RooAbsReal       , 
           ROOT.RooAbsLValue     , 
           ROOT.RooAbsRealLValue , 
           ROOT.RooRealVar       ) :

    _getter_ = None
    _setter_ = None

    if hasattr  ( t , 'getVal' ) and hasattr ( t , 'getError' ) :
        _getter_ = _rav_getvale_
    elif hasattr  ( t , 'getVal' ) :
        _getter_ = _rav_getval_

    if hasattr  ( t , 'setVal' ) and hasattr ( t , 'getMin' ) and hasattr ( t , 'getMax' ) :
        _setter_ = _rav_setvalc_
    elif hasattr  ( t , 'setVal' ) :
        _setter_ = _rav_setval_

    doc1 = """The current value, associated with the variable,
    
    >>> var = ...
    
    get value:
    =========
    
    >>> print (var.value) ## getter
    
    """
    doc2 = """The current value, associated with the variable,
    
    >>> var = ...
    
    get value:
    =========
    
    >>> print (var.value) ## getter
    
    Set value:
    ==========

    >>> var.value = 15 
    
    """
    if   _setter_  : t.value = property ( _getter_ , _setter_ , None  , doc2 )
    elif _getter_  : t.value = property ( _getter_ , _setter_ , None  , doc1 )


    doce1 = """The current error, associated with the variable,
    
    >>> var = ...
    
    Get error:
    =========
    
    >>> print (var.error) ## getter
    
    """
    doce2 = """The current error, associated with the variable,
    
    >>> var = ...
    
    Get error:
    =========
    
    >>> print (var.error) ## getter
    
    Set error:
    ==========

    >>> var.error = 15 
    
    """
    
    _gettere_ = None
    _settere_ = None

    if hasattr  ( t , 'getError' ) and hasattr ( t , 'setError' ) :
        _gettere_ = _rav_geterr_
        _settere_ = _rav_seterr_
    elif hasattr  ( t , 'getError' ) :
        _gettere_ = _rav_geterr_

    if   _settere_  : t.error = property ( _gettere_ , _settere_ , None  , doce2 )
    elif _gettere_  : t.error = property ( _gettere_ , _settere_ , None  , doce1 )

    if hasattr ( t , 'getVal' ) and not hasattr ( t , '__float__' ) :
        t.__float__ = lambda s : s.getVal()


# =============================================================================
## 
def _rar_name_ ( vname ) :
    #
    vname = vname.replace('(','Open'  )
    vname = vname.replace(')','Close' )
    vname = vname.replace('.','stop'  )
    #
    return vname 

# ==============================================================================
##   Categories
# ==============================================================================

# =============================================================================
## redefine ROOT.RooCategory constructor to define the categories
#  @code
#  sample = ROOT.RooCategory('sample','fitting sample', 'signal' , 'control' )
#  @endcode
def _rc_init_ ( self , name , title , *categories , **kwargs) :
    """Modified ROOT.RooCategory constructor to define the categories
    >>> sample = ROOT.RooCategory('sample','fitting sample', 'signal' , 'control' )
    """
    ROOT.RooCategory._old_init_ ( self , name ,  title )
    for item in categories :
        if   isinstance ( item , string_types   ) : self.defineType ( item )
        elif isinstance ( item , dictlike_types ) :
            for k in item : self.defineType ( k ,  item [ k ] )
        else :
            raise TypeError("Invalid type %s/%s" % ( item , type ( item ) ) ) 
            
    for k in kwargs :
        self.defineType ( k , kwargs[k] )
        
if not hasattr ( ROOT.RooCategory , '_old_init_' ) :
    ROOT.RooCategory._old_init_ = ROOT.RooCategory.__init__
    ROOT.RooCategory.__init__   = _rc_init_     

# =============================================================================
if root_info < ( 6, 22 ) :

    # =========================================================================
    ## iterator over label/index pairs
    #  @code
    #  category = ...
    #  for label, index  in category : ... 
    #  for label, index  in category.items()     : ...  ## ditto 
    #  for label, index  in category.iteritems() : ...  ## ditto 
    #  @endcode
    #  @see RooAbsCategory 
    def _racat_items_ ( cat ) :
        """Iterator over labelindex pairs
        >>> category = ...
        >>> for label , index in category : ... 
        >>> for label , index in category.items () : ...     ## ditto
        >>> for label , index in category.iteritems()  : ... ## ditto
        - see RooAbsCategory 
        """
        
        cit  = Ostap.Utils.Iterator ( cat.typeIterator() )
        rcat = cit .Next()
        while rcat :            
            yield rcat.GetName() , rcat.getVal()
            rcat = cit.Next() 
            
        del cit

    ROOT.RooAbsCategory.__iter__  = _racat_items_
    
    # =========================================================================
    ## Is given label (or index) defined fron this categroy?
    #  @see RooAbsCategory.isValidIndex
    #  @see RooAbsCategory.isValidLabel     
    def _racat_contains_ ( cat , item ) :
        """Is given label (or index) defined fron this category?
        - see `ROOT.RooAbsCategory.isValidIndex`
        - see `ROOT.RooAbsCategory.isValidLabel`
        """
        return ( isinstance ( item , string_types  ) and cat.isValidLabel ( item ) ) or \
               ( isinstance ( item , integer_types ) and cat.isValidIndex ( item ) )
    
    # =========================================================================
    ## get the index for the label (or label for the  index)
    #  @code
    #  category = ...
    #  index = category ['MyType']
    #  label = category [ 15 ]
    #  @endcode
    def _racat_getitem_ ( cat , item ) :
        """Get the index for the label (or label for the  index)
        >>> category = ...
        >>> index = category ['MyType']
        >>> label = category [ 15 ]
        """
        if  isinstance ( item , string_types ) :
            result = cat.lookupType ( item , False )
            if not result : raise KeyError("No '%s' label is defined!" % item )
            return result.getVal() 
        elif isinstance ( item , integer_types ) :  
            result = cat.lookupType ( item , False )
            if not result : raise  IndexError("No '%s' index is defined!" % item )
            return result.GetName()
        elif isinstance ( item , ROOT.RooCatType ) :  
            result = cat.lookupType ( item , False )
            if not result : raise  IndexError("No '%s' categroy is defined!" % item )
            return result
        
        raise TypeError("No '%s' label/index is defined!" % item )
    
else :

        
    # =========================================================================
    ## iterator over label/index pairs
    #  @code
    #  category = ...
    #  for label, index  in category : ... 
    #  @endcode
    #  @see RooAbsCategory 
    def _racat_items_ ( cat ) :
        """Iterator over labelindex pairs
        >>> category = ...
        >>> for label , index in category : ... 
        - see RooAbsCategory 
        """

        for label , index  in cat._old_iter_ () :
            yield label , index

    if not hasattr ( ROOT.RooAbsCategory , '_old_iter_' ) :
        ROOT.RooAbsCategory._old_iter_ =  ROOT.RooAbsCategory.__iter__
        
    ROOT.RooAbsCategory.__iter__   =  _racat_items_ 

    # =========================================================================
    ## Is given label (or index) defined fron this categroy?
    #  @see RooAbsCategory.hasIndex
    #  @see RooAbsCategory.hasLabel     
    def _racat_contains_ ( cat , item ) :
        """Is given label (or index) defined fron this category?
        - see `ROOT.RooAbsCategory.isValidIndex`
        - see `ROOT.RooAbsCategory.isValidLabel`
        """
        return ( isinstance ( item , string_types  ) and cat.hasLabel ( item ) ) or \
               ( isinstance ( item , integer_types ) and cat.hasIndex ( item ) )
        
    # =========================================================================
    ## get the index for the label (or label for theindex)
    #  @code
    #  category = ...
    #  index = category ['MyType']
    #  label = category [ 15 ]
    #  @endcode 
    def _racat_getitem_ ( cat , item ) :
        """Get the index for the label (or label for the  index)
        >>> category = ...
        >>> index = category ['MyType']
        >>> label = category [ 15 ]
        """
        if    isinstance ( item , string_types ) :
            if not cat.hasLabel ( item ) : raise KeyError("No '%s' label is defined!" % item )
            return cat.lookupIndex  ( item )
        elif  isinstance ( item , integer_types ) :
            if not cat.hasIndex ( item ) : raise IndexError("No '%s' index is defined!" % item )
            return cat.lookupName ( item )
        
        raise TypeError("No '%s' label/index is defined!" % item )
    
# =========================================================================
## define new category entry 
#  @code
#  category = ...
#  category[ 'MyType' ] = 16
#  category[ 25 ] = 'AnotherType
#  @endcode
def _rcat_setitem_ ( cat , label , value ) :
    """ define new category entry 
    >>> category = ...
    >>> category[ 'MyType' ] = 16
    >>> category[ 25 ] = 'AnotherType
    """
    if   isinstance ( label , string_types ) and isinstance ( value , integer_types ) :
        r = cat.defineType ( label , value )
        if r : raise KeyError('Cannot set [%s] to %s' % ( label , value ) )
    elif isinstance ( value , string_types ) and isinstance ( label , integer_types ) :
        r = cat.defineType ( value , label  )
        if r : raise IndexError('Cannot set [%s] to %s' % ( label , value ) )
        
    raise TypeError ("Invalid type for label/index %s/%s" %  ( label , value ) )


# =============================================================================
## Get the list/tuple of categories 
#  @code
#  cat = ....
#  labels = cat.labels ()
#  labels = cat.names  () ## ditto 
#  labels = cat.keys   () ## ditto 
#  @endcode
#  @see RooCategory
def _racat_labels_ ( cat  ) :
    """Get the list/tuple of categories
    >>> cat = ....
    >>> labels = cat.labels()
    >>> labels = cat.names () ## ditto 
    >>> labels = cat.keys  () ## ditto 
    """
    return tuple ( l  for (l,_) in cat.items() )

# =============================================================================
## print RooCategory instance
def _rcat_str_ ( cat ) : 
    """Print RooCategory instance"""
    return "'%s' : '%s'/%d" % ( cat.name , cat.getLabel() , cat.getIndex() )

ROOT.RooAbsCategory.items         = _racat_items_
ROOT.RooAbsCategory.iteritems     = _racat_items_
ROOT.RooAbsCategory.__contains__  = _racat_contains_
ROOT.RooAbsCategory.__getitem__   = _racat_getitem_
ROOT.RooCategory   .__setitem__   = _rcat_setitem_
ROOT.RooAbsCategory.labels        = _racat_labels_
ROOT.RooAbsCategory.names         = _racat_labels_
ROOT.RooAbsCategory.keys          = _racat_labels_
ROOT.RooCategory   .__str__       = _rcat_str_ 
ROOT.RooCategory   .__repr__      = _rcat_str_ 

_new_methods_       += [
    ROOT.RooAbsCategory.__iter__      , 
    ROOT.RooAbsCategory.items         , 
    ROOT.RooAbsCategory.iteritems     , 
    ROOT.RooAbsCategory.__contains__  , 
    ROOT.RooAbsCategory.__getitem__   , 
    ROOT.RooAbsCategory.labels        ,
    ROOT.RooAbsCategory.names         ,
    ROOT.RooAbsCategory.keys          ,
    ROOT.RooCategory   .__setitem__   ,
    ]


# ==============================================================================
## convert two yields into "total yield" and "ratio"
#  @code
#  yield1 = ROOT.RooRelaVar( ... )
#  yield2 = ROOT.RooRelaVar( ... )
#  total , ratio = total_ratio (  yield1 , yield2 ) 
#  @endcode
def  total_ratio ( var1 , var2 ) :
    """Convert two yields into 'total yield' and 'ratio'
    
    >>> yield1 = ROOT.RooRelaVar( ... )
    >>> yield2 = ROOT.RooRelaVar( ... )
    >>> total , ratio = total_ratio (  yield1 , yield2 ) 
    """
    
    assert isinstance ( var1 , ROOT.RooAbsReal,\
                        "Invalid type of 'var1': %s/%s" % ( var1 , type ( var1 ) ) )
    assert isinstance ( var2 , ROOT.RooAbsReal,\
                        "Invalid type of 'var2': %s/%s" % ( var2 , type ( var2 ) ) )
    
    name     = var1.name ,  var2.name 
    total    = add_var   ( var1 , var2 ,
                           name  = 'Total_%s_%s'               % names ,
                           title = 'Total yields of %s and %s' % names )
    ratio    = ratio_var ( var1 , var2 ,
                           name  = 'Ratio_%s_%s'               % names ,
                           title = 'Ratio of %s and %s yields' % names )
    return total , ratio 

# ==============================================================================
## convert two yields into "total yield" and "fraction"
#  @code
#  yield1 = ROOT.RooRelaVar( ... )
#  yield2 = ROOT.RooRelaVar( ... )
#  total , fraction = total_fraction (yield1 , yield2 ) 
#  @endcode
def total_fraction ( var1 , var2 ) :
    """Convert two yields into 'total yield'' and 'fraction'
    >>> yield1 = ROOT.RooRelaVar( ... )
    >>> yield2 = ROOT.RooRelaVar( ... )
    >>> total , fraction = total_fraction (  yield1 , yield2 ) 
    """
    
    assert isinstance ( var1 , ROOT.RooAbsReal,\
                        "Invalid type of 'var1': %s/%s" % ( var1 , type ( var1 ) ) )
    assert isinstance ( var2 , ROOT.RooAbsReal,\
                        "Invalid type of 'var2': %s/%s" % ( var2 , type ( var2 ) ) )
    
    names    = var1.name ,  var2.name 
    total    = add_var      ( var1 , var2 ,
                              name  = 'Total_%s_%s'                  % names ,
                              title = 'Total yields of %s and %s'    % names )
    fraction = fraction_var ( var1 , var2 ,
                              name  = 'Fraction_%s_%s'               % names ,
                              title = 'Fraction of %s and %s yields' % names )
    return total , fraction 

# ================================================================================
## construct two yields from the total yeilds and fraction
#  @code
#  total = ...
#  fraction = ...
#  y1 , y2 = two_yields ( total , fraction )
#  @endcode
def two_yields ( total , fraction ) :
    """Construct two yields from the total yield and fraction
    >>> total = ...
    >>> fraction = ...
    >>> y1 , y2 = two_yields ( total , fraction )
    """
    var1     = total
    var2     = fraction
    
    vnames   = var1.name , var2.name 
    
    name     = 'Yield1_%s_%s'          % vnames 
    title    = 'Yield1 (%s) and (%s)'  % vnames 
    
    yield1   = scale_var ( total , fraction , name , title )
    
    name     = 'Yield2_%s_%s'          % vnames 
    title    = 'Yield2 (%s) and (%s)'  % vnames 
    
    formula  = '(%s)*(1-(%s))'         % vnames  
    varlist  = ROOT.RooArgList    ( var1 , var2                      )
    yield2   = Ostap.FormulaVar   ( name , title , formula , varlist )
    
    yield2._varlist = [ var1 , var2 , varlist ]
    
    return yield1 , yield2 
    

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
## @class FIXVAR
#  Simple context manager to fix/unfix the variable 
#  @code
#  a   = ...
#  var = ...
#  with FIXVAR( var, a ) :
#      ## variable is fixed here 
#  @endcode
class FIXVAR(object):
    """ Simple context manager to fix/unfix value for RooAbsVar
    >>> a   = ...
    >>> var = ...
    >>> with FIXVAR(var,a ) :
    ...     ## variable is fixed here 
    """
    def __init__  ( self , vars , value = None ) :
        
        if isinstance ( vars , ROOT.RooAbsReal ) : vars = [ vars ]

        variables = [] 
        for v in vars :
            assert isinstance ( v , ROOT.RooAbsReal ),\
                   'FIXVAR: Invalid variable type %s/%s' % ( v   , type ( v ) )
            variables.append ( v )

        self.__variables = tuple ( variables )
        self.__fixed     = ()
        
    def __enter__ ( self        ) :

        self.__fixed = tuple ( [ c.isConstant() for c in self.variables ] )
        for v in self.variables : v.fix ()

        return self
    
    def __exit__  ( self , *_   ) :
        for v , c in zip ( self.variables , self.fixed ) :
            if not c : v.release()

    @property
    def variables ( self ) :
        """'variables' : list of variables"""
        return self.__variables

    @property
    def fixed ( self ) :
        """'fixed' : fixed variable?"""
        return self.__fixed
        
# =============================================================================
## Simple context manager to keep the binning scheme
#  @code
#  var = ...
#  var.setBins      ( 100 )
#  with KeepBinning ( var ) :
#     vars.setBins  ( 10 )
#  @endcode
class KeepBinning(object):
    """Simple context manager to keep the binning scheme
    >>> var = ...
    >>> var.setBins   ( 100 )
    >>> with KeepBins ( var ) :
    ...     vars.setBins ( 10 )
    """
    def __init__  ( self , var ) :
        self.__var  = var
        self.__bins = None
        
    def __enter__ ( self ) :
        self.__bins = self.__var.getBinning()
        return self
    
    def __exit__  ( self , *_ ) :
        ## self.__var.setBinning ( self.__bins )
        if self.__bins : self.__var.bins = self.__bins
    
    @property
    def bins ( self ) :
        """'bins' : the binning scheme"""
        return self.__bins
    
# =============================================================================
## set the binning scheme 
def _rrv_setbins_ ( self , bins ) :
    """Set the binnning scheme"""

    ## prepended/appended with name?
    name = None
    
    if bins and isinstance  ( bins , list_types ) :
        if   2 == len ( bins ) and isinstance ( bins[-1] , string_types ) :
            bins, name = bins 
        elif 2 == len ( bins ) and isinstance ( bins[ 0] , string_types ) :
            name , bins = bins 
        elif isinstance ( bins[-1] , string_types ) :
            bins = bins [:-1] 
            name = bins [ -1]
        elif isinstance ( bins[ 0] , string_types ) :
            bins = bins [1: ] 
            name = bins [0  ]
        
    if   bins and isinstance ( bins , ROOT.RooAbsBinning ) :
        if name : self.setBinning ( bins , name )
        else    : self.setBinning ( bins )
        return

    elif isinstance  ( bins , integer_types ) and 0 < bins :
        if name : self.setBins ( bins , name )
        else    : self.setBins ( bins )
        return

    if isinstance  ( bins , list_types ) :
        if   4 == len ( bins ) and isinstance ( bins[-1] , string_types ) : bins = bins [:-1]
        elif 3 == len ( bins ) and isinstance ( bins[-1] , string_types ) : bins = bins [:-1]
    else :
        logger.error ('bins: invalid binning scheme/1 %s/%s' % ( bins , type ( bins  ) ) ) 
        return 

    if   3 == len ( bins ) : ## low,high,bins or bins,low,high
        
        if   isinstance  ( bins[0] , ROOT.RooAbsReal )                       and \
             isinstance  ( bins[1] , ROOT.RooAbsReal )                       and \
             isinstance  ( bins[2] , integer_types   ) and 0 < nbins[2]      :
            
            low , high , n = bins 
            bs = ROOT.RooParamBinning ( low , high , n )
            if name : self.setBinning ( bs  , name )
            return 
            
        elif isinstance  ( bins[0] , integer_types   ) and 0 < nbins[0]      and \
             isinstance  ( bins[1] , ROOT.RooAbsReal )                       and \
             isinstance  ( bins[2] , ROOT.RooAbsReal ) : 
            
            n , low , high = bins 
            bs = ROOT.RooParamBinning ( low , high , n )
            if name : self.setBinning ( bs  , name )
            return 
            
        elif isinstance  ( bins[0] , num_types       )                       and \
             isinstance  ( bins[1] , num_types       ) and bins[0] < bins[1] and \
             isinstance  ( bins[2] , integer_types   ) and 0 < nbins[2]      :

            low , high , n = bins 
            bs = ROOT.RooUniformBinning ( low , high , n )
            if name : self.setBinning ( bs  , name )
            self.setBinning ( bs )
            return 

        elif isinstance  ( bins[0] , integer_types   ) and 0 < nbins[0]      and \
             isinstance  ( bins[1] , num_types       )                       and \
             isinstance  ( bins[2] , num_types       ) and bins[1] < bins[2] : 
            
            n , low , high = bins 
            bs = ROOT.RooUniformBinning ( low , high , n )
            if name : self.setBinning ( bs  , name )
            else    : self.setBinning ( bs         )
            return 

    elif 2 == len ( bins )  :  ## nbins, bins  or bins, nbins 

        if   isinstance  ( bins[0] , integer_types ) and 0 < bins[0]       and \
             isinstance  ( bins[1] , list_types    ) and \
             len( bins[1] ) == bins[0] + 1 :
        
            from array import array
            a  = array  ( 'd' , bins[1] )
            bs = ROOT.RooBinning ( bins[0] , a )
            if name : self.setBinning ( bs  , name )
            else    : self.setBinning ( bs         )
            return 

        elif isinstance  ( bins[1] , integer_types ) and 0 < bins[1]       and \
             isinstance  ( bins[0] , list_types    ) and \
             len( bins[0] ) == bins[1] + 1 :

            from array import array
            a  = array  ( 'd' , bins[1] )
            bs = ROOT.RooBinning ( bins[0] , a )
            if name : self.setBinning ( bs  , name )
            else    : self.setBinning ( bs         )
            return
        
    elif 4 <= len ( bins ) :  ## bins

        from array import array
        a  = array  ( 'd' , bins )
        bs = ROOT.RooBinning ( len  ( bins ) - 1 , a )
        if name : self.setBinning ( bs  , name )
        else    : self.setBinning ( bs         )
        return
                    
    logger.error ('bins: invalid binning scheme/3 %s/%s' % ( bins , type ( bins  ) ) ) 

_bins_doc_ = """Get bining scheme :

>>> vars =. ..
>>> bins = vars.bins

Set Binning scheme : 

>>> vars.bins = ROOT.RooBinning ( ... )

>>> vars.bins = 100 , [0, 2, 3, 4, 5, 6]

>>> vars.bins = 100 , 1.0 , 2.0

>>> vars.bins = 1.0 , 2.0 , 100 
"""

# =============================================================================
ROOT.RooRealVar.bins = property (  ROOT.RooRealVar.getBinning ,
                                   _rrv_setbins_  , None  , _bins_doc_ )
        
_new_methods_ += [
    ROOT.RooRealVar. bins ,
    ]


ROOT.RooRealVar.getValue = ROOT.RooRealVar.getVal
ROOT.RooRealVar.setValue = ROOT.RooRealVar.setVal


_new_methods_ += [
    ROOT.RooRealVar.getValue , 
    ROOT.RooRealVar.setValue , 
    ]


# =============================================================================
## printout for RooUniformBinning
def _rub_str_ ( bins ) :
    """Printout for RooUniformBinning"""
    l = bins. lowBound ()
    h = bins.highBound ()
    n = bins.numBoundaries () - 1
    x = bins.GetName()
    if not x : return "RooUniformBinning(%s,%s,%d)" % ( l , h , n )
    return "RooUniformBinning(%s,%s,%d,'%s')" % ( l , h , n , x )
    
ROOT.RooUniformBinning.__str__  = _rub_str_
ROOT.RooUniformBinning.__repr__ = _rub_str_

_new_methods_ += [
    ROOT.RooUniformBinning.__str__  , 
    ROOT.RooUniformBinning.__repr__ , 
    ]
# =============================================================================
## printout for RooRangeBinning
def _rrb_str_ ( bins ) :
    """Printout for RooRangeBinning"""
    l = bins. lowBound ()
    h = bins.highBound ()
    return "RooRangeBinning(%s,%s,%d')" % ( l , h , bins.GetName() )
    
ROOT.RooRangeBinning.__str__  = _rrb_str_
ROOT.RooRangeBinning.__repr__ = _rrb_str_

_new_methods_ += [
    ROOT.RooRangeBinning.__str__  ,
    ROOT.RooRangeBinning.__repr__ ,
    ]


# =============================================================================
## Does  this variable depends on another one?
#  @code
#  fun = ...
#  var = ...
#  if fun.depends_on ( var ) :
#     ...
#  @endcode
def depends_on ( fun , var ) :
    """Does  this variable depends on another one?
    
    >>> fun = ...
    >>> var = ...
    >>> if fun.depends_on ( var ) :
    ...
    
    """
    if isinstance ( var , ROOT.RooAbsCollection ) :
        for v in var :
            if depends_on ( fun , v ) : return True
        return False

    fpars = fun.getParameters ( 0 )
    
    ## direct dependency?
    if var in fpars : return True

    ## check indirect dependency
    if var and isinstance ( var , ROOT.RooAbsArg ) : 
        vvars = var.getParameters ( 0 )
        for v in vvars :
            if v in fpars : return True
            
    ##
    return False 

ROOT.RooAbsReal.depends_on  = depends_on

_new_methods_ += [
    ROOT.RooAbsReal.depends_on ,
    ]

# =============================================================================
## Create <code>RooBinnig</code> object
#  @param edges vector of bin edges 
#  @param nbins number of bins
#  @param name  binning name
#  @see RooBinning
def binning ( edges , nbins = 0 , name = '' ) :
    """Create `RooBinnig` object
    - see ROOT.RooBinning
    """
    assert isinstance ( nbins , integer_types ) and 0 <= nbins, \
           "Invalid 'nbins' parameter %s/%s" % ( nbins , type ( nbins ) )

    nb = len ( edges ) 
    assert 2 <= nb , "Invalid length of 'edges' array!"

    if 2 == nb :
        return ROOT.RooBinning ( max ( 1 , nbins ) , edges[0] , edges[1] , name ) 

    buffer = array.array ( 'd', edges )
    return ROOT.RooBinning ( nb - 1 , buffer , name ) 


# =============================================================================
## get the actual expression from <code>RooFormualVar</code>
#  @code
#  fomular = ...
#  expression = formular.expression()  
#  @endcode
def _rfv_expr_ ( var ) :
    """Get the actual expression from `RooFormualVar`
    >>> fomular = ...
    >>> expression = formular.expression()  
    """
    return var.formula().GetTitle() 

# ==============================================================================
## string representaion of the RooFormulaVar
def _rfv_str_ ( var ) :
    """String representaion of the RooFormulaVar
    """
    return '%s : %s' % ( var.expression() , var.getVal() ) 

# ==============================================================================
## string representaion of the RooFormulaVar
def _rfv_repr_ ( var ) :
    """String representation of the RooFormulaVar
    """
    return '%s : %s' % ( var.expression() , var.getVal() ) 

if (6,22) <= root_info : 
    ROOT.RooFormulaVar. expression = _rfv_expr_
    ROOT.RooFormulaVar. __str__    = _rfv_str_
    ROOT.RooFormulaVar. __repr__   = _rfv_repr_
    _new_methods_ += [
        ROOT.RooFormulaVar. expression , 
        ROOT.RooFormulaVar. __str__    , 
        ROOT.RooFormulaVar. __repr__   , 
        ]
else :
    Ostap.FormulaVar. __str__      = _rfv_str_
    Ostap.FormulaVar. __repr__     = _rfv_repr_
    _new_methods_ += [
        Ostap.FormulaVar. __str__    , 
        Ostap.FormulaVar. __repr__   , 
        ]

# =============================================================================

# ==============================================================================
# primitive functions for RooAbsReal objects 
# ==============================================================================

# ==============================================================================
## absolute value   \f$ f = abs{ab} \f$
#  @code
#  var = ...
#  e   = var_abs ( var ) 
#  @endcode 
def var_abs ( a , b = 1 , name = '' , title = '' ) :
    """Absolute value: f(x) = abs(ab)
    >>> var = ...
    >>> e   = var_abs ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.abs ( float ( a ) * float ( b ) )
        return ROOT.RooFit.RooConst ( ab )     
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 0 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 0 ) 
    #
    return Ostap.MoreRooFit.Abs( a, b , name , title )
    #
# ==============================================================================
## exponent  \f$ f = \mathrm{e}^{ab}\f$
#  @code
#  var = ...
#  e   = var_exp ( var ) 
#  @endcode 
def var_exp ( a , b = 1 , name = '' , title = '' ) :
    """Exponent: f(x) = exp(ab)
    >>> var = ...
    >>> e   = var_exp ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.exp ( float ( a ) * float ( b ) )
        return ROOT.RooFit.RooConst ( ab )          ## RETURN
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 1 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 1 ) 
    #
    return Ostap.MoreRooFit.Exp ( a, b , name , title )

# ==============================================================================
## logarithm  \f$ f = \log ab \f$
#  @code
#  var = ...
#  e   = var_log ( var ) 
#  @endcode 
def var_log ( a , b = 1 , name = '' , title = '' ) :
    """logarithm f(x) = log(ab)
    >>> var = ...
    >>> e   = var_log ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.log ( float ( a ) * float ( b ) )       ## RETURN
        return ROOT.RooFit.RooConst ( ab )
    #
    return Ostap.MoreRooFit.Log ( a, b , name , title ) 

# ==============================================================================
## logarithm  \f$ f = \log10 ab \f$
#  @code
#  var = ...
#  e   = var_log10 ( var ) 
#  @endcode 
def var_log10 ( a , b = 1 , name = '' , title = '' ) :
    """logarithm f(x) = log10(ab)
    >>> var = ...
    >>> e   = var_log10 ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.log10 ( float ( a ) * float ( b ) )       ## RETURN
        return ROOT.RooFit.RooConst ( ab )
    #
    return Ostap.MoreRooFit.Log10 ( a, b , name , title ) 


# ==============================================================================
## error function \f$ f = erf ( ab) \f$
#  @code
#  var = ...
#  e   = var_erf ( var ) 
#  @endcode 
def var_erf ( a , b = 1 , name = '' , title = '' ) :
    """Error function f(x) = erf(ab)
    >>> var = ...
    >>> e   = var_erf ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.erf ( float ( a ) * float ( b ) )
        return ROOT.RooFit.RooConst ( ab )          ## RETURN
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 0 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 0 ) 
    #
    return Ostap.MoreRooFit.Erf ( a, b , name , title ) 

# ==============================================================================
## complementary error function \f$ f = erfc ( ab) \f$
#  @code
#  var = ...
#  e   = var_erfc ( var ) 
#  @endcode 
def var_erfc ( a , b = 1 , name = '' , title = '' ) :
    """Error function f(x) = erfc(ab)
    >>> var = ...
    >>> e   = var_erfc ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.erfc ( float ( a ) * float ( b ) )
        return ROOT.RooFit.RooConst ( ab )          ## RETURN
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 1 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 1 ) 
    #
    return Ostap.MoreRooFit.Erf ( a, b , name , title ) 


# ==============================================================================
## Sine \f$ f = \sin ab\f$
#  @code
#  var = ...
#  e   = var_sin ( var ) 
#  @endcode 
def var_sin ( a , b = 1 , name = '' , title = '' ) :
    """Sine  f(x) = sin(ab)
    >>> var = ...
    >>> e   = var_sin ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.sin ( float ( a ) * float ( b ) )
        return ROOT.RooFit.RooConst ( ab )          ## RETURN 
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 0 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 0 ) 
    #
    return Ostap.MoreRooFit.Sin ( a, b , name , title ) 
    
# ==============================================================================
## Cosine \f$ f = \cos ab\f$
#  @code
#  var = ...
#  e   = var_cos ( var ) 
#  @endcode 
def var_cos ( a , b = 1 , name = '' , title = '' ) :
    """Cosine  f(x) = cos(ab)
    >>> var = ...
    >>> e   = var_cos ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.cos ( float ( a ) * float ( b ) )
        return ROOT.RooFit.RooConst ( ab )       ## RETURN
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 1 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 1 ) 
    #
    return Ostap.MoreRooFit.Cos ( a, b , name , title )


# ==============================================================================
## Tangent\f$ f = \tan ab\f$
#  @code
#  var = ...
#  e   = var_tan ( var ) 
#  @endcode 
def var_tan ( a , b = 1 , name = '' , title = '' ) :
    """Tangent  f(x) = tan(ab)
    >>> var = ...
    >>> e   = var_tan ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.tan ( float ( a ) * float ( b ) )   ##  RETURN
        return ROOT.RooFit.RooConst ( ab )
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 0 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 0 ) 
    #
    return Ostap.MoreRooFit.Tan ( a, b , name , title )


# ==============================================================================
## Hyprbolic sine \f$ f = \sinh ab\f$
#  @code
#  var = ...
#  e   = var_sinh ( var ) 
#  @endcode 
def var_sinh ( a , b = 1 , name = '' , title = '' ) :
    """Hyperbolic sine  f(x) = sinh(ab)
    >>> var = ...
    >>> e   = var_sinh ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.sinh ( float ( a ) * float ( b ) )
        return ROOT.RooFit.RooConst ( ab )          ## RETURN 
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 0 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 0 ) 
    #
    return Ostap.MoreRooFit.Sinh ( a, b , name , title ) 
    

# ==============================================================================
## Hyperbolic cosine \f$ f = \cosh ab\f$
#  @code
#  var = ...
#  e   = var_cosh ( var ) 
#  @endcode 
def var_cosh ( a , b = 1 , name = '' , title = '' ) :
    """Hyperbolic cosine  f(x) = cos(ab)
    >>> var = ...
    >>> e   = var_cosh ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.cosh ( float ( a ) * float ( b ) )
        return ROOT.RooFit.RooConst ( ab )       ## RETURN
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 1 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 1 ) 
    #
    return Ostap.MoreRooFit.Cosh ( a, b , name , title )

# ==============================================================================
## Hyperboilic tangent\f$ f = \tanh ab\f$
#  @code
#  var = ...
#  e   = var_tanh ( var ) 
#  @endcode 
def var_tanh ( a , b = 1 , name = '' , title = '' ) :
    """Hyperbolic tangent  f(x) = tanh(ab)
    >>> var = ...
    >>> e   = var_tanh ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :        
        ab = math.tanh ( float ( a ) * float ( b ) ) 
        return ROOT.RooFit.RooConst ( ab )      ## RETURN
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 0 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 0 ) 
    #
    return Ostap.MoreRooFit.Tanh ( a, b , name , title )

# ==============================================================================
## Hyperboilic secant \f$ f = \sech ab\f$
#  @code
#  var = ...
#  e   = var_sech ( var ) 
#  @endcode 
def var_sech ( a , b = 1 , name = '' , title = '' ) :
    """Hyperbolic tangent  f(x) = tanh(ab)
    >>> var = ...
    >>> e   = var_tanh ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :        
        ab = Ostap.Math.sech ( float ( a ) * float ( b ) ) 
        return ROOT.RooFit.RooConst ( ab )      ## RETURN
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 1 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 1 ) 
    #
    return Ostap.MoreRooFit.Sech ( a, b , name , title )

# ==============================================================================
## arctangent \f$ f = atan2 (a,b)\f$
#  @code
#  var = ...
#  e   = var_atan2 ( var ) 
#  @endcode 
def var_atan2 ( a , b = 1 , name = '' , title = '' ) :
    """Inverse tangent  f(x) = atan2(a,b)
    >>> var = ...
    >>> e   = var_atan2 ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.atan2 ( float ( a ) , float ( b ) ) 
        return ROOT.RooFit.RooConst ( ab )      ## RETURN
    return Ostap.MoreRooFit.Atan2 ( a, b , name , title )

# ==============================================================================
## Bessel function \f$ f = J_{\nu}(x)\f$
#  @code
#  x   = ...
#  nu  = ... 
#  e   = var_bessel_J ( x , nu ) 
#  @endcode
def var_bessel_J ( x , nu = 0 , name = '' , title = '' ) :
    """Bessel function f(x) = J_{nu}(x)
    >>> x  = ...
    >>> nu = ... 
    >>> e  = var_bessel_J ( x , nu  ) 
    """
    fx  = isinstance ( x  , num_types )
    fnu = isinstance ( nu , num_types )
    if fa and fb :
        val = mve.bessel_J ( float ( fnu ) , float ( fx ) ) 
        return ROOT.RooFit.RooConst ( val )      ## RETURN
    return Ostap.MoreRooFit.BesselJ ( x , nu , name , title )

# ==============================================================================
## Bessel function \f$ f = Y_{\nu}(x)\f$
#  @code
#  x   = ...
#  nu  = ... 
#  e   = var_bessel_Y ( x , nu ) 
#  @endcode
def var_bessel_Y ( x , nu = 0 , name = '' , title = '' ) :
    """Bessel function f(x) = Y_{nu}(x)
    >>> x  = ...
    >>> nu = ... 
    >>> e  = var_bessel_Y ( x , nu  ) 
    """
    fx  = isinstance ( x  , num_types )
    fnu = isinstance ( nu , num_types )
    if fa and fb :
        val = mve.bessel_Y ( float ( fnu ) , float ( fx ) ) 
        return ROOT.RooFit.RooConst ( val )      ## RETURN
    return Ostap.MoreRooFit.BesselY ( x , nu , name , title )

# ==============================================================================
## Bessel function \f$ f = I_{\nu}(x)\f$
#  @code
#  x   = ...
#  nu  = ... 
#  e   = var_bessel_I ( x , nu ) 
#  @endcode
def var_bessel_I ( x , nu = 0 , name = '' , title = '' ) :
    """Bessel function f(x) = I_{nu}(x)
    >>> x  = ...
    >>> nu = ... 
    >>> e  = var_bessel_I ( x , nu  ) 
    """
    fx  = isinstance ( x  , num_types )
    fnu = isinstance ( nu , num_types )
    if fa and fb :
        val = mve.bessel_I ( float ( fnu ) , float ( fx ) ) 
        return ROOT.RooFit.RooConst ( val )      ## RETURN
    return Ostap.MoreRooFit.BesselI ( x , nu , name , title )

# ==============================================================================
## Bessel function \f$ f = K_{\nu}(x)\f$
#  @code
#  x   = ...
#  nu  = ... 
#  e   = var_bessel_K ( x , nu ) 
#  @endcode
def var_bessel_K ( x , nu = 0 , name = '' , title = '' ) :
    """Bessel function f(x) = K_{nu}(x)
    >>> x  = ...
    >>> nu = ... 
    >>> e  = var_bessel_K ( x , nu  ) 
    """
    fx  = isinstance ( x  , num_types )
    fnu = isinstance ( nu , num_types )
    if fa and fb :
        val = mve.bessel_K ( float ( fnu ) , float ( fx ) ) 
        return ROOT.RooFit.RooConst ( val )      ## RETURN
    return Ostap.MoreRooFit.BesselK ( x , nu , name , title )



# ==============================================================================
## maximal \f$ f = max (a,b)\f$
#  @code
#  var1 = ...
#  var2 = ...
#  var  = var_max ( var1 , var2 ) 
#  @endcode 
def var_max ( a , b = 1 , name = '' , title = '' ) :
    """Maximal from two fnuctions f(x) = max(a,b)
    >>> var = ...
    >>> e   = var_max ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = max ( float ( a ) , float ( b ) ) 
        return ROOT.RooFit.RooConst ( ab )      ## RETURN
    return Ostap.MoreRooFit.MaxV ( a, b , name , title )

# ==============================================================================
## minimal \f$ f = min (a,b)\f$
#  @code
#  var1 = ...
#  var2 = ...
#  var  = var_min ( var1 , var2 ) 
#  @endcode 
def var_min ( a , b = 1 , name = '' , title = '' ) :
    """Minimal from two fnuctions f(x) = min(a,b)
    >>> var = ...
    >>> e   = var_min ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = min ( float ( a ) , float ( b ) ) 
        return ROOT.RooFit.RooConst ( ab )      ## RETURN
    return Ostap.MoreRooFit.MinV ( a, b , name , title )

# ==============================================================================
## Gamma function \f$ f =    \Gamma(ab) \f$
#  @code
#  a = ...
#  e = var_gamma ( a ) 
#  @endcode 
def var_gamma ( a , b = 1 , name = '' , title = '' ) :
    """Gamma function  f = Gamma(ab)
    >>> a = ...
    >>> e = var_gamma ( a ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.gamma ( float ( a ) * float ( b ) )      
        return ROOT.RooFit.RooConst ( ab )           ## RETURN 
    #
    return Ostap.MoreRooFit.Gamma ( a, b , name , title ) 

# ==============================================================================
## logarithm of Gamma function \f$ f = \log \Gamma(ab) \f$
#  @code
#  a = ...
#  e = var_lgamma ( a ) 
#  @endcode 
def var_lgamma ( a , b = 1 , name = '' , title = '' ) :
    """logarithm of Gamma function  f = log Gamma(ab)
    >>> a = ...
    >>> e = var_lgamma ( a ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.lgamma ( float ( a ) * float ( b ) )      
        return ROOT.RooFit.RooConst ( ab )           ## RETURN 
    #
    return Ostap.MoreRooFit.LGamma ( a, b , name , title ) 

# ==============================================================================
## 1/Gamma function \f$ f = \frac{1}{\Gamma(ab)} \f$
#  @code
#  a = ...
#  e = var_igamma ( a ) 
#  @endcode 
def var_igamma ( a , b = 1 , name = '' , title = '' ) :
    """1/Gamma  f = 1/Gamma(ab)
    >>> a = ...
    >>> e = var_igamma ( a ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = Ostap.Math.igamma ( float ( a ) * float ( b ) )
        return ROOT.RooFit.RooConst ( ab )              ## RETURN 
    #
    return Ostap.MoreRooFit.IGamma ( a, b , name , title ) 


# =============================================================================
## Sum of two RooAbsReal objects
#  @code
#  v1 = ...
#  v2 = ...
#  v  = var_sum ( v1 ,  v2 )
#  @endcode
def var_sum ( v1 , v2 , name = '' , title = '' ) :
    """Sum of two RooAbsReal objects
    >>> v1 = ...
    >>> v2 = ...
    >>> v  = var_sum ( v1 , v2 )
    """
    f1 = isinstance ( v1 , num_types )
    f2 = isinstance ( v2 , num_types )
    
    if f1 and f2 :
        r = float ( v1 ) + float ( v2 ) 
        return ROOT.RooFit.RooConst ( r )                 ## RETURN    
    elif f1 and iszero ( v1 ) : return v2
    elif f2 and iszero ( v2 ) : return v1
    #
    return Ostap.MoreRooFit.Addition ( v1 , v2 , name , title )


# =============================================================================
## Subtraction of two RooAbsReal objects
#  @code
#  v1 = ...
#  v2 = ...
#  v  = var_sub ( v1 , v2 )
#  @endcode
def var_sub ( v1 , v2 , name = '' , title = '' ) :
    """Subraction of two RooAbsReal objects
    >>> v1 = ...
    >>> v2 = ...
    >>> v  = var_sub ( v1 , v2 )  
    """
    f1 = isinstance ( v1 , num_types )
    f2 = isinstance ( v2 , num_types )    
    if f1 and f2 :
        r = float ( v1 ) - float ( v2 ) 
        return ROOT.RooFit.RooConst ( r )                 ## RETURN    
    elif f1 and iszero ( v1 ) : return var_mul ( v2 , -1 , name , title ) 
    elif f2 and iszero ( v2 ) : return v1 
    # 
    return Ostap.MoreRooFit.Subtraction ( v1 , v2 , name , title ) 


# =============================================================================
## Product of two RooAbsReal objects
#  @code
#  v1 = ...
#  v2 = ...
#  v  = var_mul ( v1 ,  v2 )
#  @endcode
def var_mul ( v1 , v2 , name = '' , title = '' ) :
    """Product of two RooAbsReal objects
    >>> v1 = ...
    >>> v2 = ...
    >>> v  = var_mul ( v1 ,  v2 )
    """
    f1 = isinstance ( v1 , num_types )
    f2 = isinstance ( v2 , num_types )    
    if f1 and f2 :
        r = float ( v1 ) * float ( v2 ) 
        return ROOT.RooFit.RooConst ( r )                 ## RETURN
    elif f1 and iszero ( v1 ) : return ROOT.RooFit.RooConst ( 0 ) 
    elif f2 and iszero ( v2 ) : return ROOT.RooFit.RooConst ( 0 ) 
    elif f1 and isone  ( v1 ) : return v2 
    elif f2 and iseone ( v2 ) : return v1 
    #
    return Ostap.MoreRooFit.Product ( v1 , v2 , name , title ) 

# =============================================================================
## Division of two RooAbsReal objects
#  @code
#  v1 = ...
#  v2 = ...
#  v  = var_div ( v1 , v2 )
#  @endcode
def var_div ( v1 , v2 , name = '' , title = '' ) :
    """Division of two RooAbsReal objects
    >>> v1 = ...
    >>> v2 = ...
    >>> v  = var_div ( v1 , v2 )  
    """
    f1 = isinstance ( v1 , num_types )
    f2 = isinstance ( v2 , num_types )    
    if f1 and f2 :
        r = float ( v1 ) / float ( v2 ) 
        return ROOT.RooFit.RooConst ( r )                 ## RETURN
    elif f1 and iszero ( v1 ) : return ROOT.RooFit.RooConst ( 0 ) 
    elif f2 and isone  ( v2 ) : return v1
    #
    return Ostap.MoreRooFit.Division ( v1 , v2 , name , title ) 

# =============================================================================
## pow for two RooAbsReal objects
#  @code
#  v1 = ...
#  v2 = ...
#  v  = var_pow ( v1 , v2 )
#  @endcode
def var_pow ( v1 , v2 , name = '' , title = '' ) :
    """pow for two RooAbsReal objects
    >>> v1 = ...
    >>> v2 = ...
    >>> v  = var_pow ( v1 ,  v2 ) 
    """
    f1 = isinstance ( v1 , num_types )
    f2 = isinstance ( v2 , num_types )    
    if f1 and f2 :
        r = float ( v1 ) ** float ( v2 ) 
        return ROOT.RooFit.RooConst ( r )                 ## RETURN    
    elif f2 and iszero  ( v2 ) : return ROOT.RooFit.RooConst ( 1 )
    elif f2 and isone   ( v2 ) : return v1 
    elif f1 and iszero  ( v1 ) : return ROOT.RooFit.RooConst ( 0 )
    elif f1 and isone   ( v1 ) : return ROOT.RooFit.RooConst ( 1 )
    #
    return Ostap.MoreRooFit.Power ( v1 , v2 , name , title ) 

# ==============================================================================
## "Fraction" of two RooAbsReal objects: f = a/(a+b)
#  @code
#  a = ...
#  b = ...
#  e   = var_fraction ( a , b ) 
#  @endcode 
def var_fraction ( a , b , name = '' , title = '' ) :
    """'Fraction'  f(x) = a/(a+b)
    >>> a = ...
    >>> b = ...
    >>> e = var_fraction ( a , b  ) 
    """
    f1 = isinstance ( v1 , num_types )
    f2 = isinstance ( v2 , num_types )    
    if f1 and f2 :
        r = float ( v1) / ( float ( v1 ) + float ( v1 ) )  
        return ROOT.RooFit.RooConst ( r )                 ## RETURN    
    elif f1 and iszero  ( v1 ) : return ROOT.RooFit.RooConst ( 0 ) 
    elif f2 and iszero  ( v2 ) : return ROOT.RooFit.RooConst ( 1 ) 
    #
    return Ostap.MoreRooFit.Fraction ( v1 , v2 , name , title ) 


# ==============================================================================
## "Asymmetry" of two RooAbsReal objects: f = (a-b)/(a+b)
#  @code
#  a = ...
#  b = ...
#  e = var_asymmetry ( a , b ) 
#  @endcode 
def var_asymmetry ( a , b , name = '' , title = '' ) :
    """'Asymmetry'  f(x) = (a-b)/(a+b)
    >>> a = ...
    >>> b = ...
    >>> e = var_asymmetry ( a , b  ) 
    """
    f1 = isinstance ( v1 , num_types )
    f2 = isinstance ( v2 , num_types )    
    if f1 and f2 :
        r = ( float ( v1 ) - float ( v2 ) ) / ( float ( v1 ) + float ( v2 ) )  ## 
        return ROOT.RooFit.RooConst ( r )                                ## RETURN    
    elif f1 and iszero  ( v1 ) : return ROOT.RooFit.RooConst ( -1 ) 
    elif f2 and iszero  ( v2 ) : return ROOT.RooFit.RooConst (  1 ) 
    #
    return Ostap.MoreRooFit.Asymmetry ( v1 , v2 , name , title ) 

scale_var     = var_mul
add_var       = var_sum
sum_var       = var_sum
ratio_var     = var_div
fraction_var  = var_fraction
asymmetry_var = var_asymmetry

_decorated_classes_ = (
    ##
    ROOT.RooRealVar               ,
    ROOT.RooConstVar              ,
    ROOT.RooFormulaVar            ,
    ROOT.RooAbsReal               ,
    ROOT.RooAbsRealLValue         ,
    ROOT.RooUniformBinning        ,
    ROOT.RooBinning               ,
    ROOT.RooRangeBinning          ,
    ##
    ROOT.RooCategory              ,
)

_new_methods_ = tuple ( _new_methods_ ) 


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
