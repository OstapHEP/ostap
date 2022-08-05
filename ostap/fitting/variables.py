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
    ) 
# =============================================================================
from   builtins               import range
from   ostap.math.base        import doubles
from   ostap.core.core        import VE, hID, Ostap
from   ostap.math.reduce      import root_factory 
from   ostap.core.meta_info   import root_info 
from   ostap.core.ostap_types import ( num_types     , list_types   ,
                                       integer_types , string_types )
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
_new_methods_ = []
# =============================================================================
## Factory for deserialization of generic objects
#  @attention it stores the constructor kaprameters as local attributes
def root_store_factory ( klass , *params ) :
    """Factory for deserialization of generic object
    - attention: it stores the constructor kaprameters as local attributes
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
    if var  in fpars : return True
        
    ## check indirect dependency
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
## Dedicated unpickling factory for RooRealVar
def _rrv_factory ( args , errors , binnings , fixed , *attrs ) :
    """ Dedicated unpickling factory for `ROOT.RooRealVar`
    """

    ## create it
    rrv = ROOT.RooRealVar ( *args )

    ## set errors if needed 
    if   errors and  2 == len ( errors ) : rrv.setAsymError ( *errors ) 
    elif errors and  1 == len ( errors ) : rrv.setError     ( *errors )

    for b in binnings :
        if b : rrv.setBinning ( b , b.GetName ()) 

    rrv.setConstant ( fixed )
    
    if attrs :
        battrs = attrs[0] 
        for n , a in battrs : rrv.setAttribute          ( n , a )
        if 1 < len ( attrs ) :
            sattrs = attrs[1] 
            for n , a in sattrs : rrv.setStringAttribute    ( n , a )
            if 2 < len ( attrs ) :
                tattrs = attrs[2]             
                for n , a in tattrs : rrv.setTransientAttribute ( n , a ) 
    
    return rrv

# =============================================================================
## Reducing of <code>RooRealVar</code> for pickling/unpickling 
#  @see RooRooRealVar 
def _rrv_reduce ( rrv ) :
    """ Reducing of `ROOT.RooRealVar` for pickling 
    - see ROOT.RooRooRealVar 
    """
    name    = rrv.name 
    title   = rrv.title
    value   = rrv.getVal () 

    has_min = rrv.hasMin ()
    has_max = rrv.hasMax ()

    ## constructor arguments 
    if has_min and has_max :
        args = name , title , value ,  rrv.getMin () , rrv.getMax ()              , rrv.getUnit ()
    elif has_min : 
        args = name , title , value ,  rrv.getMin () , ROOT.RooNumber.infinity () , rrv.getUnit () 
    elif has_max : 
        args = name , title , value , -ROOT.RooNumber.infinity () , rrv.getMax () , rrv.getUnit () 
    else :
        args = name , title , value ,  rrv.getUnit () 

    ## errors 
    if   rrv.hasAsymError () :
        errors = rrv.getAsymErrorLo () , rrv.getAsymErrorHi ()
    elif rrv.hasError   () :
        errors = rrv.getError() ,
    else :
        errors = () 

    ## binings 
        
    binnings = tuple (   rrv.getBinning ( n , False )       for n in rrv.getBinningNames     () )

    ## attributes:
    
    battrs   = tuple ( ( n , rrv.getAttribute          ( n ) ) for   n       in rrv.attributes          () ) 
    sattrs   = tuple ( ( n , a                               ) for ( n , a ) in rrv.stringAttributes    () ) 
    tattrs   = tuple ( ( n , rrv.getTransientAttribute ( n ) ) for   n       in rrv.transientAttributes () ) 

    ## fixed ? 
    fixed = True if rrv.isConstant() else False 

    if tattrs : 
        content = args , errors , binnings , fixed , battrs , sattrs , tattrs
    elif sattrs :
        content = args , errors , binnings , fixed , battrs , sattrs
    elif battrs : 
        content = args , errors , binnings , fixed , battrs 
    else :        
        content = args , errors , binnings , fixed 
    
    return _rrv_factory , content 


ROOT.RooRealVar.__reduce__ = _rrv_reduce

_new_methods_ += [
    ROOT.RooRealVar.__reduce__ ,
    ]


# =============================================================================
## factory for unpickling of <code>RooFormulaVar</code> and
#  <code>Ostap::FormulaVar</code>
def _rfv_factory ( klass , args , vars ) :
    """Factory for unpickling of `RooFormulaVar` and `Ostap.FormulaVar`
    """

    lst = ROOT.RooArgList ()
    for v in vars : lst.add ( v ) 

    margs = list  ( args  )
    margs.append  ( lst   )
    margs = tuple ( margs ) 
    
    rfv = klass   ( *margs )
    
    rfv.__args = args, vars, lst

    return rfv
    
# =============================================================================
## Reduce <code>RooFormulaVar</code> and <code>Ostap::FormulaVar</code> for pickling
#  @see RooFormulaVar 
#  @see Ostap::FormulaVar 
def _rfv_reduce ( rfv ) : 
    """Reduce `RooFormulaVar` and `Ostap::FormulaVar` for pickling
    - see RooFormulaVar 
    - see Ostap.FormulaVar 
    """

    name       = rfv.GetName     ()
    title      = rfv.GetTitle    ()
    
    rform      = rfv.formula     ()    
    expression = rform.GetTitle  () 

    vars       = tuple ( d for d in rform.actualDependents() ) 
    args       = name , title , expression
    
    return _rfv_factory , ( type ( rfv ) , args , vars ) 


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
    """String representaion of the RooFormulaVar
    """
    return '%s : %s' % ( var.expression() , var.getVal() ) 

if (6,22) <= root_info : 
    ROOT.RooFormulaVar. expression = _rfv_expr_
    ROOT.RooFormulaVar. __str__    = _rfv_str_
    ROOT.RooFormulaVar. __repr__   = _rfv_repr_
    ROOT.RooFormulaVar.__reduce__  = _rfv_reduce
    _new_methods_ += [
        ROOT.RooFormulaVar. expression , 
        ROOT.RooFormulaVar. __str__    , 
        ROOT.RooFormulaVar. __repr__   , 
        ROOT.RooFormulaVar.__reduce__  ,
        ]
else :
    Ostap.FormulaVar.__reduce__    = _rfv_reduce
    Ostap.FormulaVar. __str__      = _rfv_str_
    Ostap.FormulaVar. __repr__     = _rfv_repr_
    Ostap.FormulaVar.__reduce__    = _rfv_reduce
    _new_methods_ += [
        Ostap.FormulaVar.__reduce__  ,
        Ostap.FormulaVar. __str__    , 
        Ostap.FormulaVar. __repr__   , 
        Ostap.FormulaVar.__reduce__  , 
        ]


# =============================================================================
## Reduce <code>Ostap::MoreFooFit::TwoVars</code> objects
#  @see Ostap::MoreRooFit.TwoVars
def _r2v_reduce ( var ) :
    """Reduce `Ostap::MoreFooFit::TwoVars` objects
    - see Ostap.MoreRooFit.TwoVars
    """
    return root_store_factory , ( type ( var ) , var.name , var.title , var.x() , var.y() )

Ostap.MoreRooFit.TwoVars.__reduce__  = _r2v_reduce

_new_methods_ += [
    Ostap.MoreRooFit.TwoVars.__reduce__
    ]

# ===================================================================
## Reduce <code>Ostap::MoreFooFit::Addition</code> objects
#  @see Ostap::MoreRooFit.Addition
def _radd1_reduce ( var ) :
    """Reduce `Ostap.MoreFooFit.Addition` objects
    - see Ostap.MoreRooFit.Addition
    """
    content = type ( var ) , var.name , var.title , var.x () , var.y ()
    return root_store_factory , content 

# ===================================================================
## Reduce <code>Ostap::MoreFooFit::Addition3</code> objects
#  @see Ostap::MoreRooFit.Addition3
def _radd2_reduce ( var ) :
    """Reduce `Ostap.MoreFooFit.Addition` objects
    - see Ostap.MoreRooFit.Addition2
    """
    content = type ( var ) , var.name , var.title , var.x () , var.y () , var.c1 () , var.c2 () 
    return root_store_factory , content

# =============================================================================
## reduce <code>Ostap::MoreRooFit::Id<code> object
#  @see Ostap.MoreRooFit::Id
def _rid_reduce ( var ) :
    """Reduce `Ostap.MoreRooFit.Id` object
    - see Ostap.MoreRooFit.Id
    """
    return root_store_factory , ( type ( var ) ,
                            var.name     ,
                            var.title    ,
                            var.x ()     )


Ostap.MoreRooFit.Addition    .__reduce__  = _radd1_reduce
Ostap.MoreRooFit.Addition2   .__reduce__  = _radd2_reduce

Ostap.MoreRooFit.Subtraction .__reduce__  = _r2v_reduce
Ostap.MoreRooFit.Product     .__reduce__  = _r2v_reduce
Ostap.MoreRooFit.ProductPdf  .__reduce__  = _r2v_reduce

Ostap.MoreRooFit.Id          .__reduce__  = _rid_reduce

_new_methods_ += [
    Ostap.MoreRooFit.Addition    .__reduce__  , 
    Ostap.MoreRooFit.Addition2   .__reduce__  , 
    Ostap.MoreRooFit.Subtraction .__reduce__  , 
    Ostap.MoreRooFit.Product     .__reduce__  , 
    Ostap.MoreRooFit.ProductPdf  .__reduce__  , 
    Ostap.MoreRooFit.Id          .__reduce__  ,
]

# ===================================================================
## Reduce <code>Ostap::MoreFooFit::Combination</code> objects
#  @see Ostap::MoreRooFit.Combination
def _rcomb_reduce ( var ) :
    """Reduce `Ostap.MoreFooFit.Combination` objects
    - see Ostap.MoreRooFit.Combination
    """
    return root_store_factory , ( type ( var ) ,
                            var.name     ,
                            var.title    ,
                            var.x()      ,
                            var.y()      ,
                            var.alpha () ,
                            var.beta  () ,
                            var.gamma () )

# ===================================================================
## Reduce <code>Ostap::MoreFooFit::Asymmetry</code> objects
#  @see Ostap::MoreRooFit.Asymmetry
def _rasym_reduce ( var ) :
    """Reduce  `Ostap::MoreFooFit::Asymmetry` objects
    - see Ostap.MoreRooFit.Asymmetry
    """
    return root_store_factory , ( type ( var ) ,
                            var.name , var.title ,
                            var.x()  , var.y()   ,
                            var.scale () )

Ostap.MoreRooFit.Combination.__reduce__  = _rcomb_reduce
Ostap.MoreRooFit.Asymmetry.  __reduce__  = _rasym_reduce

_new_methods_ += [
    Ostap.MoreRooFit.Combination.__reduce__  , 
    Ostap.MoreRooFit.Asymmetry.  __reduce__  ,
    ]

# =============================================================================
## Reduce <code>RooConstVar</code>
#  @see RooConstVar 
def _rconst_reduce ( var ) :
    """ Reduce `ROOT.RooConstVar`
    - see ROOT.RooConstVar
    """
    return root_store_factory , ( type ( var )  ,
                            var.name      ,
                            var.title     ,
                            float ( var ) )  

ROOT.RooConstVar.__reduce__ = _rconst_reduce 

_new_methods_ += [
    ROOT.RooConstVar.__reduce__ ,
    ]


# =============================================================================
## Reduce <code>Ostap::MoreRooFit::Constant</code>
#  @see RooConstVar 
def _rconst2_reduce ( var ) :
    """ Reduce `ROOT.RooConstVar`
    - see ROOT.RooConstVar
    """
    return root_store_factory , ( type ( var )  ,
                            var.name      ,
                            var.title     ,
                            float ( var ) ,
                            var.vlst()    )

Ostap.MoreRooFit.Constant.  __reduce__  = _rconst2_reduce

_new_methods_ += [
    Ostap.MoreRooFit.Constant.  __reduce__  ,
    ]

# =============================================================================
## Reduce <code>Ostap::MoreRooFit::Bernstein</code> object 
#  @see Ostap::MoreRooDit::Bernstein
def _rbern_reduce ( var ) :
    """Reduce `Ostap.MoreRooFit.Bernstein` object 
    - see Ostap.MoreRooDit.Bernstein
    """
    return root_store_factory , ( type ( var )  ,
                            var.name      ,
                            var.title     ,
                            var.xvar   () ,
                            var.pars   () ,                              
                            var.xmin   () ,
                            var.xmax   () )

Ostap.MoreRooFit.Bernstein.  __reduce__  = _rbern_reduce
_new_methods_ += [
    Ostap.MoreRooFit.Bernstein.  __reduce__ ,
    ]

# =============================================================================
## Reduce <code>Ostap::MoreRooFit::Monotonic</code> object 
#  @see Ostap::MoreRooDit::Monotonic
def _rmono_reduce ( var ) :
    """Reduce `Ostap.MoreRooFit.Monotonic` object 
    - see Ostap.MoreRooDit.Monotonic
    """
    return root_store_factory , ( type ( var )  ,
                            var.name      ,
                            var.title     ,
                            var.xvar ()   ,
                            var.pars ()   ,                               
                            var.monotonic().increasing() , 
                            var.xmin ()   ,
                            var.xmax ()   ,
                            var.a    ()   ,
                            var.b    ()   )

Ostap.MoreRooFit.Monotonic.  __reduce__  = _rmono_reduce

_new_methods_ += [
    Ostap.MoreRooFit.Monotonic.  __reduce__ , 
    ]

# =============================================================================
## Reduce <code>Ostap::MoreRooFit::Convex</code> object 
#  @see Ostap::MoreRooDit::Monotonic
def _rconv1_reduce ( var ) :
    """Reduce `Ostap.MoreRooFit.Monotonic` object 
    - see Ostap.MoreRooDit.Monotonic
    """
    return root_store_factory , ( type ( var )  ,
                            var.name      ,
                            var.title     ,
                            var.xvar ()   ,
                            var.pars ()   ,                               
                            var.convex ().increasing () , 
                            var.convex ().convex     () , 
                            var.xmin ()   ,
                            var.xmax ()   ,
                            var.a    ()   ,
                            var.b    ()   )

Ostap.MoreRooFit.Convex.  __reduce__  = _rconv1_reduce

_new_methods_ += [
    Ostap.MoreRooFit.Convex.  __reduce__ ,
    ]

# =============================================================================
## Reduce <code>Ostap::MoreRooFit::ConvexOnly</code> object 
#  @see Ostap::MoreRooDit::ConvexOnly
def _rconv2_reduce ( var ) :
    """Reduce `Ostap.MoreRooFit.ConvexOnly` object 
    - see Ostap.MoreRooDit.ConvexOnly
    """
    return root_store_factory , ( type ( var )  ,
                            var.name      ,
                            var.title     ,
                            var.xvar ()   ,
                            var.pars ()   ,                               
                            var.convex ().convex () , 
                            var.xmin ()   ,
                            var.xmax ()   ,
                            var.a    ()   ,
                            var.b    ()   )

Ostap.MoreRooFit.ConvexOnly.  __reduce__  = _rconv2_reduce

_new_methods_ += [
    Ostap.MoreRooFit.ConvexOnly.  __reduce__  ,
    ]

# =============================================================================
## unpickle <code>Ostap::MoreFooFit::BSpline</code> objects
#  @see Ostap::MoreRooFit::BSpline
def _rbspl_factory ( klass , name , title , xvar , knots , pars ) :
    """unpickle `Ostap.MoreFooFit.BSpline` objects
    - see Ostap.MoreRooFit.BSpline
    """
    plst  = ROOT.RooArgList()
    for v in pars : plst.add ( v )
    knots = doubles ( knots )
    #
    obj        = klass ( name , title , xvar , knots , pars )
    obj.__args = pars
    return obj 

# =============================================================================
## Reduce <code>Ostap::MoreRooFit::BSpline</code> object 
#  @see Ostap::MoreRooDit::BSpline
def _rbspl_reduce ( spl ) :
    """Reduce `Ostap.MoreRooFit.BSpline` object 
    - see Ostap.MoreRooDit.BSpline
    """
    return _rbspl_factory, ( type ( spl ) ,
                             spl.name     ,
                             spl.title    ,
                             spl.xvar  () ,
                             array.array( 'd' , spl.knots() ) ,
                             spl.pars()   )

Ostap.MoreRooFit.BSpline.  __reduce__  = _rbspl_reduce

_new_methods_ += [
    Ostap.MoreRooFit.BSpline.  __reduce__  ,
    ]

# ============================
## reduce <code>Ostap::Models::Uniform<code> object
#  @see Ostap::Models.Uniform 
def _runi_reduce_ ( uni ) :
    """reduce <code>Ostap::Models::Uniform<code> object
    - see Ostap::Models.Uniform 
    """
    tail = ()
    ##
    d = uni.dim() 
    if   1 == d : tail = uni.x() ,
    elif 2 == d : tail = uni.x() , uni.y()
    elif 3 == d : tail = uni.x() , uni.y() , uni.z()
    ##
    return root_store_factory , ( type ( uni ) ,
                            uni.name     ,
                            uni.title    ) + tail

Ostap.Models.Uniform.__reduce__ = _runi_reduce_

_new_methods_ += [
    Ostap.Models.Uniform.__reduce__ ,
    ]

# =============================================================================
## unpickle RooUniformBinning object
#  @see RooUniformBinnig  
def _rub_factory ( *args ) :
    """unpickle RooUniformBinning object
    -see ROOT.RooUniformBinnig
    """
    return ROOT.RooUniformBinning ( *args )

# =============================================================================
## reduce uniform binning scheme
#  @see RoUniformBinnig 
def _rub_reduce_ ( rub ) :
    """Reduce RooUniformBininkg Object
    - see ROOT.RooUniformBinning
    """
    nbins = rub.numBoundaries()
    if nbins : nbins -= 1
    content = rub.lowBound () , rub.highBound(), nbins , rub.GetName()
    
    return _rub_factory, content

# =============================================================================
## unpickle RooBinning object
#  @see RooBinnig  
def _rb_factory ( data , name  ) :
    """unpickle RooBinning object
    -see ROOT.RooUniformBinnig
    """
    return ROOT.RooBinning ( len ( data ) - 1 , data [ 0 ] , name )
# =============================================================================
## reduce RooBinning object
#  @see RooBinning 
def _rb_reduce_ ( rb  )  :
    """Reduce RooBinning object
    - see ROOT.RooBinning 
    """
    if rb.isUniform() : return _rub_reduce_ ( rb )

    nb    = rb.numBoundaries() 
    ab    = rb.array ()
    data  = array.array ( 'd' , [ 1.0 * ab[i] for i in range ( nb ) ] )

    content = data, rb.GetName()  
    return _rb_factory, content

# ==========================================================================
## unpickle RooRangeBinning object
#  @see RooRangeBinnig 
def _rrb_factory ( low , high , name ) :
    """unpickle RooRangeBinning object"""
    return ROOT.RooRangeBinning( low , high , name )
# ============================================================================
## reduce RooRangeBinnig object
#  @see RooRangeBinnig 
def _rrb_reduce_ ( rrb ) :
    """Reduce RooRangeBinnig object"""
    return _rrb_factory ,  ( rrb.lowBound() , rrb.highBound() , rrb.GetName() ) 

ROOT.RooBinning       .__reduce__ = _rb_reduce_
ROOT.RooUniformBinning.__reduce__ = _rub_reduce_
ROOT.RooRangeBinning  .__reduce__ = _rrb_reduce_

_new_methods_ += [
    ROOT.RooBinning       .__reduce__ ,
    ROOT.RooUniformBinning.__reduce__ , 
    ROOT.RooRangeBinning  .__reduce__ , 
    ]
    
# =============================================================================
__item_store = set() 
# =============================================================================
## Unpickle RooArgList/RooArgSet
#  @attention it stores the constructore arguments  in the local store 
#  @see RooArgList
def _rac_factory ( klass , *args ) :
    """ Unpickle RooArgList/RooArgSet  instance
    0 attention: it stores the constructore arguments  in the local store 
    """
    c = klass () 
    for a in args : c.add ( a )
    __item_store.add ( args ) 
    return c
# =============================================================================
## reduce `RooArgList`
def _ral_reduce_ ( rac ) :
    """Reduce `RooArgList` instances"""    
    return _rac_factory, ( ROOT.RooArgList, ) + tuple ( a for a in rac ) 

# =============================================================================
## reduce `RooArgSet`
def _ras_reduce_ ( rac ) :
    """Reduce `RooArgSet` instances"""    
    return _rac_factory, ( ROOT.RooArgSet, ) + tuple ( a for a in rac )  

ROOT.RooArgSet  .__reduce__ = _ras_reduce_
ROOT.RooArgList .__reduce__ = _ral_reduce_

_new_methods_ += [
    ROOT.RooArgSet  .__reduce__ , 
    ROOT.RooArgList .__reduce__ ,
    ]

# ==========================================================================--
## reduce  RooGaussian object 
def _rgau_reduce_ ( pdf ) :
    """Reduce `RooGaussian` object"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.getX     () ,
                            pdf.getMean  () , 
                            pdf.getSigma () )

ROOT.RooGaussian.__reduce__ = _rgau_reduce_ 

# ==========================================================================--
## Get the original fractions from the <code>RooAddPdf</code>
#  @code
#  addpdf = ...
#  fractions , recursive = addpdf.orig_fracs() 
#  @endcode 
#  @see RooAddPdf 
def _raddpdf_fractions ( radd ) :
    """ Get the original fractions from the <code>RooAddPdf</code>
    >>> addpdf = ...
    >>> fractions , recursive = addpdf.orig_fracs() 
    """
    rec   = ctypes.c_bool ()
    fracs = Ostap.MoreRooFit.fractions ( radd , rec )
    return fracs, rec.value
    
# ==========================================================================--
ROOT.RooAddPdf  .orig_fracs = _raddpdf_fractions

# ==========================================================================--
## reduce  RooAddPdf object 
def _raddpdf_reduce_ ( pdf ) :
    """Reduce `RooAddPdf` object"""
    content = type ( pdf ) , pdf.name , pdf.title , pdf.pdfList()
    pars    = pdf.coefList ()
    if 1 <= len ( pars ) :
        content   = content + pdf.orig_fracs () 
    return root_store_factory , content

ROOT.RooAddPdf.__reduce__ = _raddpdf_reduce_ 


# =============================================================================
## Is this <code>RooProdPdf</code> object conditional ?
#  @see RooProdPdf 
def _rprodpdf_cond_ ( pdf ) :
    """Is this `RooProdPdf` object conditional ?
    - see `ROOT.RooProdPdf` 
    """
    for p in pdf.pdfList() :
        vset = pdf.findPdfNSet ( p )
        if not set             : continue
        elif 1 <= len ( vset ) : return True 
    return False

ROOT.RooProdPdf.conditional = _rprodpdf_cond_ 

# =============================================================================
## reduce RooProdPdf 
#  @see RooProdPdf
def _rprodpdf_reduce_ ( pdf ) :
    """Reduce `ROOT.RooProdPdf`
    - see `ROOT.RooProdPdf`
    """
    if pdf.conditional () :
        import pickle
        raise pickle.PicklingError("RooProdPdf is conditional (cannot be picked)")
    
    content = type ( pdf ) , pdf.name , pdf.title , pdf.pdfList()
    return root_store_factory , content

ROOT.RooProdPdf.__reduce__ = _rprodpdf_reduce_ 

_new_methods_ += [
    ROOT.RooGaussian .__reduce__   ,
    ROOT.RooAddPdf   .__reduce__   ,
    ROOT.RooAddPdf   . orig_fracs  ,
    ROOT.RooProdPdf  .__reduce__   ,
    ROOT.RooProdPdf  . conditional ,
    ]

## # =============================================================================
## ## reduce BreitWigner
## def _rbw_reduce_ ( pdf ):
##     """Reduce BreitWigner"""
##     return root_store_factory , ( type ( pdf )       ,
##                             pdf.name           ,
##                             pdf.title          ,
##                             pdf.x      ()      , 
##                             pdf.mass   ()      ,
##                             pdf.widths ()[0]   ,
##                             pdf.breit_wigner() ) 

## Ostap.Models.BreitWigner.__reduce__ = _rbw_reduce_ 

## # =============================================================================
## ## reduce BreitWignerMC
## def _rbwmc_reduce_ ( pdf ):
##     """Reduce BreitWignerMC"""
##     return root_store_factory , ( type ( pdf )          ,
##                             pdf.name              ,
##                             pdf.title             ,
##                             pdf.x      ()         , 
##                             pdf.mass   ()         ,
##                             pdf.widths ()         ,
##                             pdf.breit_wigner_MC() ) 

## Ostap.Models.BreitWignerMC.__reduce__ = _rbwmc_reduce_ 

# =============================================================================
## reduce Flatte
def _rflatte_reduce_ ( pdf ):
    """Reduce Flatte"""
    return root_store_factory , ( type ( pdf )     ,
                            pdf.name         ,
                            pdf.title        ,
                            pdf.x      ()    , 
                            pdf.mass   ()    ,
                            pdf.widths ()[0] ,
                            pdf.widths ()[1] ,
                            pdf.widths ()[2] ,
                            pdf.flatte  ()   ) 

Ostap.Models.Flatte.__reduce__ = _rflatte_reduce_ 

# =============================================================================
## reduce Voigt
def _rvoigt_reduce_ ( pdf ):
    """Reduce Voigt"""
    return root_store_factory , ( type ( pdf )  ,
                            pdf.name      ,
                            pdf.title     ,
                            pdf.x      () , 
                            pdf.m0     () ,
                            pdf.gamma  () ,
                            pdf.sigma  () )

Ostap.Models.Voigt      .__reduce__ = _rvoigt_reduce_ 
Ostap.Models.PseudoVoigt.__reduce__ = _rvoigt_reduce_ 

# =============================================================================
## reduce CrystalBall
def _rcb_reduce_ ( pdf ):
    """Reduce CristalBall"""
    return root_store_factory , ( type ( pdf )  ,
                            pdf.name      ,
                            pdf.title     ,
                            pdf.x      () , 
                            pdf.m0     () ,
                            pdf.sigma  () ,
                            pdf.alpha  () ,
                            pdf.n      () )

Ostap.Models.CrystalBall   .__reduce__ = _rcb_reduce_ 
Ostap.Models.CrystalBallRS .__reduce__ = _rcb_reduce_ 

# =============================================================================
## reduce CrystalBallDS
def _rcb2_reduce_ ( pdf ):
    """Reduce CristalBallDS"""
    return root_store_factory , ( type ( pdf )  ,
                            pdf.name      ,
                            pdf.title     ,
                            pdf.x      () , 
                            pdf.m0     () ,
                            pdf.sigma  () ,
                            pdf.alphaL () ,
                            pdf.nL     () ,
                            pdf.alphaR () ,
                            pdf.nR     () )

Ostap.Models.CrystalBallDS .__reduce__ = _rcb2_reduce_ 

# =============================================================================
## reduce Needham
def _rneedham_reduce_ ( pdf ):
    """Reduce Needham"""
    return root_store_factory , ( type ( pdf )  ,
                            pdf.name      ,
                            pdf.title     ,
                            pdf.x      () , 
                            pdf.m0     () ,
                            pdf.sigma  () ,                            
                            pdf.a0     () ,
                            pdf.a1     () ,
                            pdf.a2     () )

Ostap.Models.Needham .__reduce__ = _rneedham_reduce_ 

# =============================================================================
## reduce Apollonious
def _rapo_reduce_ ( pdf ):
    """Reduce Apollonios"""
    return root_store_factory , ( type ( pdf )  ,
                            pdf.name      ,
                            pdf.title     ,
                            pdf.x      () , 
                            pdf.m0     () ,
                            pdf.sigma  () ,                            
                            pdf.alpha  () ,                            
                            pdf.n      () ,
                            pdf.b      () )

Ostap.Models.Apollonios.__reduce__ = _rapo_reduce_ 

# =============================================================================
## reduce Apollonious2
def _rapo2_reduce_ ( pdf ):
    """Reduce Apollonios2"""
    return root_store_factory , ( type ( pdf )  ,
                            pdf.name      ,
                            pdf.title     ,
                            pdf.x      () , 
                            pdf.m0     () ,
                            pdf.sigmaL () ,                            
                            pdf.sigmaR () ,                            
                            pdf.beta   () )

Ostap.Models.Apollonios2.__reduce__ = _rapo2_reduce_ 

# =============================================================================
## reduce BifurcatedGauss
def _rgbf_reduce_ ( pdf ):
    """Reduce BifurcatedGauss"""
    return root_store_factory , ( type ( pdf )  ,
                            pdf.name      ,
                            pdf.title     ,
                            pdf.x      () , 
                            pdf.peak   () ,
                            pdf.sigmaL () ,                            
                            pdf.sigmaR () )

Ostap.Models.BifurcatedGauss.__reduce__ = _rgbf_reduce_ 


# =============================================================================
## reduce GenGaussV1
def _rggv1_reduce_ ( pdf ):
    """Reduce GenGaussV1"""
    return root_store_factory , ( type ( pdf )  ,
                            pdf.name      ,
                            pdf.title     ,
                            pdf.x      () , 
                            pdf.mu     () ,
                            pdf.alpha  () ,                            
                            pdf.beta   () )

Ostap.Models.GenGaussV1.__reduce__ = _rggv1_reduce_ 


# =============================================================================
## reduce GenGaussV2
def _rggv2_reduce_ ( pdf ):
    """Reduce GenGaussV2"""
    return root_store_factory , ( type ( pdf )  ,
                            pdf.name      ,
                            pdf.title     ,
                            pdf.x      () , 
                            pdf.xi     () ,
                            pdf.alpha  () ,                            
                            pdf.kappa  () )

Ostap.Models.GenGaussV2.__reduce__ = _rggv2_reduce_ 

# =============================================================================
## reduce SkewGauss
def _rskg_reduce_ ( pdf ):
    """Reduce SkewGauss"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () , 
                            pdf.xi       () ,
                            pdf.omega    () ,                            
                            pdf.alpha    () )

Ostap.Models.SkewGauss.__reduce__ = _rskg_reduce_ 

# =============================================================================
## reduce ExGauss
def _rexg_reduce_ ( pdf ):
    """Reduce ExGauss"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () , 
                            pdf.mu       () ,
                            pdf.varsigma () ,                            
                            pdf.k        () )

Ostap.Models.ExGauss.__reduce__ = _rexg_reduce_ 

# =============================================================================
## reduce NormalLaplace
def _rnl_reduce_ ( pdf ):
    """Reduce NormalLaplace"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () , 
                            pdf.mu       () ,
                            pdf.varsigma () ,                            
                            pdf.kL       () ,
                            pdf.kL       () )

Ostap.Models.NormalLaplace.__reduce__ = _rnl_reduce_ 

# =============================================================================
## reduce Novosibirsk
def _rnovo_reduce_ ( pdf ):
    """Reduce Novisibirsk"""
    return root_store_factory , ( type ( pdf )  ,
                            pdf.name      ,
                            pdf.title     ,
                            pdf.x      () , 
                            pdf.peak   () ,
                            pdf.sigma  () ,                            
                            pdf.tau    () )

Ostap.Models.Novosibirsk.__reduce__ = _rnovo_reduce_ 

# =============================================================================
## reduce Bukin
def _rbukin_reduce_ ( pdf ):
    """Reduce Bukin"""
    return root_store_factory , ( type ( pdf )  ,
                            pdf.name      ,
                            pdf.title     ,
                            pdf.x      () , 
                            pdf.peak   () ,
                            pdf.sigma  () ,                            
                            pdf.xi     () ,
                            pdf.rhoL   () ,
                            pdf.rhoR   () )

Ostap.Models.Bukin.__reduce__ = _rbukin_reduce_ 

# =============================================================================
## reduce StudentT
def _rstt_reduce_ ( pdf ):
    """Reduce StudentT"""
    return root_store_factory , ( type ( pdf )  ,
                            pdf.name      ,
                            pdf.title     ,
                            pdf.x      () , 
                            pdf.mu     () ,
                            pdf.sigma  () ,                            
                            pdf.n      () )

Ostap.Models.StudentT.__reduce__ = _rstt_reduce_ 


# =============================================================================
## reduce BifurcatedStudentT
def _rbstt_reduce_ ( pdf ):
    """Reduce BifurcatedStudentT"""
    return root_store_factory , ( type ( pdf )  ,
                            pdf.name      ,
                            pdf.title     ,
                            pdf.x      () , 
                            pdf.mu     () ,
                            pdf.sigmaL () ,                            
                            pdf.sigmaR () ,                            
                            pdf.nL     () ,
                            pdf.nR     () )

Ostap.Models.BifurcatedStudentT.__reduce__ = _rbstt_reduce_ 

# =============================================================================
## reduce PearsonIV
def _rp4_reduce_ ( pdf ):
    """Reduce PearsonIV"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () , 
                            pdf.mu       () ,
                            pdf.varsigma () ,                            
                            pdf.n        () ,                            
                            pdf.kappa    () )

Ostap.Models.PearsonIV.__reduce__ = _rp4_reduce_ 


# =============================================================================
## reduce GramCharlierA
def _rgca_reduce_ ( pdf ):
    """Reduce GramCharlierA"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () , 
                            pdf.mu       () ,
                            pdf.sigma    () ,                            
                            pdf.kappa3   () ,
                            pdf.kappa4   () )

Ostap.Models.GramCharlierA.__reduce__ = _rgca_reduce_ 

# =============================================================================
## reduce PhaseSpace2
def _rps2_reduce_ ( pdf ):
    """Reduce PhaseSpace2"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () , 
                            pdf.m1       () ,
                            pdf.m2       () )

Ostap.Models.PhaseSpace2.__reduce__ = _rps2_reduce_ 


# =============================================================================
## reduce PhaseSpaceLeft
def _rpsl_reduce_ ( pdf ):
    """Reduce PhaseSpaceLeft"""
    return root_store_factory , ( type ( pdf )     ,
                            pdf.name         ,
                            pdf.title        ,
                            pdf.x         () , 
                            pdf.threshold () ,
                            pdf.scale     () ,
                            pdf.left      () )


Ostap.Models.PhaseSpaceLeft.__reduce__ = _rpsl_reduce_ 

# =============================================================================
## reduce PhaseSpaceRight
def _rpsr_reduce_ ( pdf ):
    """Reduce PhaseSpaceRight"""
    return root_store_factory , ( type ( pdf )     ,
                            pdf.name         ,
                            pdf.title        ,
                            pdf.x         () , 
                            pdf.threshold () ,
                            pdf.L         () ,
                            pdf.N         () )


Ostap.Models.PhaseSpaceRight.__reduce__ = _rpsr_reduce_ 

# =============================================================================
## reduce PhaseSpaceNL
def _rpsnl_reduce_ ( pdf ):
    """Reduce PhaseSpaceNL"""
    return root_store_factory , ( type ( pdf )     ,
                            pdf.name         ,
                            pdf.title        ,
                            pdf.x         () , 
                            pdf.low       () ,
                            pdf.high      () ,
                            pdf.N         () ,
                            pdf.L         () )


Ostap.Models.PhaseSpaceNL.__reduce__ = _rpsnl_reduce_ 


# =============================================================================
## reduce PhaseSpacePol
def _rpspol_reduce_ ( pdf ):
    """Reduce PhaseSpacePol"""
    return root_store_factory , ( type ( pdf )     ,
                            pdf.name         ,
                            pdf.title        ,
                            pdf.x         () , 
                            pdf.psNL      () ,
                            pdf.phis      () )

Ostap.Models.PhaseSpacePol.__reduce__ = _rpspol_reduce_ 


# =============================================================================
## reduce PhaseSpaceLeftExpoPol
def _rpslepol_reduce_ ( pdf ):
    """Reduce PhaseSpaceLeftExpoPol"""
    return root_store_factory , ( type ( pdf )     ,
                            pdf.name         ,
                            pdf.title        ,
                            pdf.x         () , 
                            pdf.psleft    () ,
                            pdf.tau       () ,
                            pdf.scale     () ,
                            pdf.phis      () )

Ostap.Models.PhaseSpaceLeftExpoPol.__reduce__ = _rpslepol_reduce_ 


## PhaseSpace23L - not done .... why?


# =============================================================================
## reduce PolyPositive 
def _rpolpos_reduce_ ( pdf ):
    """Reduce PolyPositive"""
    return root_store_factory , ( type ( pdf )     ,
                            pdf.name         ,
                            pdf.title        ,
                            pdf.x         () , 
                            pdf.phis      () ,
                            pdf.xmin      () ,
                            pdf.xmax      () )

Ostap.Models.PolyPositive.__reduce__ = _rpolpos_reduce_

# =============================================================================
## reduce PolyPositiveEven 
def _rpolpose_reduce_ ( pdf ):
    """Reduce PolyPositiveEven"""
    return root_store_factory , ( type ( pdf )     ,
                            pdf.name         ,
                            pdf.title        ,
                            pdf.x         () , 
                            pdf.phis      () ,
                            pdf.xmin      () ,
                            pdf.xmax      () )

Ostap.Models.PolyPositiveEven.__reduce__ = _rpolpose_reduce_ 

# =============================================================================
## reduce PolyMonotonic
def _rpolmon_reduce_ ( pdf ):
    """Reduce PolyMonotonic"""
    return root_store_factory , ( type ( pdf )      ,
                            pdf.name          ,
                            pdf.title         ,
                            pdf.x          () , 
                            pdf.phis       () ,
                            pdf.xmin       () ,
                            pdf.xmax       () ,
                            pdf.increasing () )
                            
Ostap.Models.PolyMonotonic.__reduce__ = _rpolmon_reduce_ 

# =============================================================================
## reduce PolyConvex
def _rpolcon_reduce_ ( pdf ):
    """Reduce PolyConvex"""
    return root_store_factory , ( type ( pdf )      ,
                            pdf.name          ,
                            pdf.title         ,
                            pdf.x          () , 
                            pdf.phis       () ,
                            pdf.xmin       () ,
                            pdf.xmax       () ,
                            pdf.increasing () ,
                            pdf.convex     () )
                            
Ostap.Models.PolyConvex.__reduce__ = _rpolcon_reduce_ 

# =============================================================================
## reduce PolyConvexOnly
def _rpolcono_reduce_ ( pdf ):
    """Reduce PolyConvexOnly"""
    return root_store_factory , ( type ( pdf )      ,
                            pdf.name          ,
                            pdf.title         ,
                            pdf.x          () , 
                            pdf.phis       () ,
                            pdf.xmin       () ,
                            pdf.xmax       () ,
                            pdf.convex     () )
                            
Ostap.Models.PolyConvexOnly.__reduce__ = _rpolcono_reduce_ 

# =============================================================================
## reduce ExpoPositive
def _rexppos_reduce_ ( pdf ):
    """Reduce ExpoPoisitive"""
    return root_store_factory , ( type ( pdf )      ,
                            pdf.name          ,
                            pdf.title         ,
                            pdf.x          () , 
                            pdf.tau        () ,
                            pdf.phis       () ,
                            pdf.xmin       () ,
                            pdf.xmax       () )
                            
Ostap.Models.ExpoPositive.__reduce__ = _rexppos_reduce_ 


# =============================================================================
## reduce PolySigmoid
def _rpolsigm_reduce_ ( pdf ):
    """Reduce PoLySigmoid"""
    return root_store_factory , ( type ( pdf )      ,
                            pdf.name          ,
                            pdf.title         ,
                            pdf.x          () , 
                            pdf.phis       () ,
                            pdf.xmin       () ,
                            pdf.xmax       () ,
                            pdf.alpha      () ,
                            pdf.x0         () )

Ostap.Models.PolySigmoid.__reduce__ = _rpolsigm_reduce_ 

# =============================================================================
## reduce TwoExpoPositive
def _r2exppos_reduce_ ( pdf ):
    """Reduce TwoExpoPositive"""
    return root_store_factory , ( type ( pdf )      ,
                            pdf.name          ,
                            pdf.title         ,
                            pdf.x          () ,                            
                            pdf.alpha      () ,
                            pdf.delta      () ,
                            pdf.x0         () ,
                            pdf.phis       () ,
                            pdf.xmin       () ,
                            pdf.xmax       () )

Ostap.Models.TwoExpoPositive.__reduce__ = _r2exppos_reduce_ 

# =============================================================================
## reduce GammaDist
def _rgamdist_reduce_ ( pdf ):
    """Reduce GammaDist"""
    return root_store_factory , ( type ( pdf )      ,
                            pdf.name          ,
                            pdf.title         ,
                            pdf.x          () ,                            
                            pdf.k          () ,
                            pdf.theta      () )

Ostap.Models.GammaDist     .__reduce__ = _rgamdist_reduce_ 
Ostap.Models.LogGammaDist  .__reduce__ = _rgamdist_reduce_ 
Ostap.Models.Log10GammaDist.__reduce__ = _rgamdist_reduce_ 

# =============================================================================
## reduce GenGammaDist
def _rggamdist_reduce_ ( pdf ):
    """Reduce GenGammaDist"""
    return root_store_factory , ( type ( pdf )      ,
                            pdf.name          ,
                            pdf.title         ,
                            pdf.x          () ,                            
                            pdf.k          () ,
                            pdf.theta      () ,
                            pdf.p          () ,
                            pdf.low        () )

Ostap.Models.GenGammaDist.__reduce__ = _rggamdist_reduce_ 


# =============================================================================
## reduce Amoroso
def _ramoroso_reduce_ ( pdf ):
    """Reduce Amoroso"""
    return root_store_factory , ( type ( pdf )      ,
                            pdf.name          ,
                            pdf.title         ,
                            pdf.x          () ,                            
                            pdf.theta      () ,
                            pdf.alpha      () ,
                            pdf.beta       () ,
                            pdf.a          () )

Ostap.Models.Amoroso.__reduce__ = _ramoroso_reduce_ 


# =============================================================================
## reduce LogGamma
def _rloggam_reduce_ ( pdf ):
    """Reduce LogGamma"""
    return root_store_factory , ( type ( pdf )      ,
                            pdf.name          ,
                            pdf.title         ,
                            pdf.x          () ,                            
                            pdf.nu         () ,
                            pdf.lambd      () ,
                            pdf.alpha      () )

Ostap.Models.LogGamma   .__reduce__ = _rloggam_reduce_ 

# =============================================================================
## reduce BetaPrime
def _rbetap_reduce_ ( pdf ):
    """Reduce BetaPrime"""
    return root_store_factory , ( type ( pdf )      ,
                            pdf.name          ,
                            pdf.title         ,
                            pdf.x          () ,                            
                            pdf.alpha      () ,
                            pdf.beta       () ,
                            pdf.scale      () ,
                            pdf.shift      () )

Ostap.Models.BetaPrime.__reduce__ = _rbetap_reduce_ 


# =============================================================================
## reduce Landau
def _rlandau_reduce_ ( pdf ):
    """Reduce Landau"""
    return root_store_factory , ( type ( pdf )      ,
                            pdf.name          ,
                            pdf.title         ,
                            pdf.x          () ,                            
                            pdf.scale      () ,
                            pdf.shift      () )

Ostap.Models.Landau.__reduce__ = _rlandau_reduce_ 

# =============================================================================
## reduce SinhAsinh
def _rshash_reduce_ ( pdf ):
    """Reduce SinhAsinh"""
    return root_store_factory , ( type ( pdf )   ,
                            pdf.name       ,
                            pdf.title      ,
                            pdf.x       () ,                            
                            pdf.mu      () ,
                            pdf.sigma   () ,
                            pdf.epsilon () ,
                            pdf.delta   () )

Ostap.Models.SinhAsinh.__reduce__ = _rshash_reduce_ 

# =============================================================================
## reduce JohnsonSU
def _rjsu_reduce_ ( pdf ):
    """Reduce JohnsonSU"""
    return root_store_factory , ( type ( pdf )   ,
                            pdf.name       ,
                            pdf.title      ,
                            pdf.x       () ,                            
                            pdf.xi      () ,
                            pdf.lambd   () ,
                            pdf.delta   () ,
                            pdf.gamma   () )

Ostap.Models.JohnsonSU.__reduce__ = _rjsu_reduce_ 


# =============================================================================
## reduce ATLAS
def _ratlas_reduce_ ( pdf ):
    """Reduce Atlas"""
    return root_store_factory , ( type ( pdf )   ,
                            pdf.name       ,
                            pdf.title      ,
                            pdf.x       () ,                            
                            pdf.mu      () ,
                            pdf.sigma   () )

Ostap.Models.Atlas    .__reduce__ = _ratlas_reduce_ 
Ostap.Models.Sech     .__reduce__ = _ratlas_reduce_ 
Ostap.Models.Logistic .__reduce__ = _ratlas_reduce_ 
Ostap.Models.Hat      .__reduce__ = _ratlas_reduce_ 
Ostap.Models.Up       .__reduce__ = _ratlas_reduce_ 



# =============================================================================
## reduce Slash
def _rslash_reduce_ ( pdf ):
    """Reduce Slasj"""
    return root_store_factory , ( type ( pdf )   ,
                            pdf.name       ,
                            pdf.title      ,
                            pdf.x       () ,                            
                            pdf.mu      () ,
                            pdf.scale   () )
Ostap.Models.Slash    .__reduce__ = _rslash_reduce_ 


# =============================================================================
## reduce ARGUS
def _rargus_reduce_ ( pdf ):
    """Reduce ARGUS"""
    return root_store_factory , ( type ( pdf )   ,
                            pdf.name       ,
                            pdf.title      ,
                            pdf.x       () ,                            
                            pdf.mu      () ,
                            pdf.c       () ,
                            pdf.chi     () )

Ostap.Models.Argus  .__reduce__ = _rargus_reduce_ 

# =============================================================================
## reduce GenArgus
def _rgargus_reduce_ ( pdf ):
    """Reduce GenArgus"""
    return root_store_factory , ( type ( pdf )   ,
                            pdf.name       ,
                            pdf.title      ,
                            pdf.x       () ,                            
                            pdf.mu      () ,
                            pdf.c       () ,
                            pdf.chi     () ,
                            pdf.dp     () )

Ostap.Models.GenArgus .__reduce__ = _rgargus_reduce_ 


# =============================================================================
## reduce Losev
def _rlosev_reduce_ ( pdf ):
    """Reduce Losev"""
    return root_store_factory , ( type ( pdf )   ,
                            pdf.name       ,
                            pdf.title      ,
                            pdf.x       () ,                            
                            pdf.mu      () ,
                            pdf.alpha   () ,
                            pdf.beta    () )

Ostap.Models.Losev.__reduce__ = _rlosev_reduce_ 


# =============================================================================
## reduce AsymmetricLaplace
def _ralap_reduce_ ( pdf ):
    """Reduce AsymmetricLaplace"""
    return root_store_factory , ( type ( pdf )   ,
                            pdf.name       ,
                            pdf.title      ,
                            pdf.x       () ,                            
                            pdf.mu      () ,
                            pdf.lambdaL () ,
                            pdf.lambdaR () )

Ostap.Models.AsymmetricLaplace.__reduce__ = _ralap_reduce_ 


# =============================================================================
## reduce FupN
def _rfupn_reduce_ ( pdf ):
    """Reduce FupN"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.N        () ,                            
                            pdf.mu       () ,
                            pdf.varsigma () )

Ostap.Models.FupN.__reduce__ = _rfupn_reduce_ 

# =============================================================================
## reduce Tsallis
def _rtsal_reduce_ ( pdf ):
    """Reduce Tsallis"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.N        () ,                            
                            pdf.T        () ,
                            pdf.mass     () )

Ostap.Models.Tsallis.__reduce__ = _rtsal_reduce_ 


# =============================================================================
## reduce QGSM
def _rqgsm_reduce_ ( pdf ):
    """Reduce QGSM"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.b        () ,                            
                            pdf.mass     () )

Ostap.Models.QGSM.__reduce__ = _rqgsm_reduce_ 


# =============================================================================
## reduce TwoExpos
def _r2expo_reduce_ ( pdf ):
    """Reduce TwoExpos"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.alpha    () ,                            
                            pdf.delta    () ,
                            pdf.x0       () )

Ostap.Models.TwoExpos.__reduce__ = _r2expo_reduce_ 


# =============================================================================
## reduce DoubleGauss
def _r2gau_reduce_ ( pdf ):
    """Reduce DoubleGauss"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.sigma    () ,                            
                            pdf.fraction () ,
                            pdf.scale    () ,
                            pdf.mean     () )

Ostap.Models.DoubleGauss.__reduce__ = _r2gau_reduce_ 


# =============================================================================
## reduce Gumbel
def _rgumbel_reduce_ ( pdf ):
    """Reduce Gumbel"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.mu       () ,                            
                            pdf.beta     () )

Ostap.Models.Gumbel.__reduce__ = _rgumbel_reduce_ 


# =============================================================================
## reduce Weibull
def _rweibull_reduce_ ( pdf ):
    """Reduce Weibull"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.scale    () ,                            
                            pdf.shape    () ,
                            pdf.shift    () )

Ostap.Models.Weibull.__reduce__ = _rweibull_reduce_ 


# =============================================================================
## reduce Rice
def _rrice_reduce_ ( pdf ):
    """Reduce Rice"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name         ,
                            pdf.title        ,
                            pdf.x         () ,                            
                            pdf.nu        () ,                            
                            pdf.varshigma () ,
                            pdf.shift     () )

Ostap.Models.Rice.__reduce__ = _rrice_reduce_ 


# =============================================================================
## reduce RasingCosine 
def _rrcos_reduce_ ( pdf ):
    """Reduce RaisingCosine"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.mean     () ,                            
                            pdf.scale    () )

Ostap.Models.RaisingCosine.__reduce__ = _rrcos_reduce_ 


# =============================================================================
## reduce q-Gaussian
def _rqgau_reduce_ ( pdf ):
    """Reduce QGaussian"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.mean     () ,                            
                            pdf.q        () ,
                            pdf.scale    () )

Ostap.Models.QGaussian.__reduce__ = _rqgau_reduce_ 

# =============================================================================
## reduce Hyperbolic
def _rhyp_reduce_ ( pdf ):
    """Reduce Hyperbolic"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.mu       () ,                            
                            pdf.sigma    () ,
                            pdf.zeta     () ,
                            pdf.kappa    () )

Ostap.Models.Hyperbolic.__reduce__ = _rhyp_reduce_ 


# =============================================================================
## reduce GenHyperbolic
def _rghyp_reduce_ ( pdf ):
    """Reduce GenHyperbolic"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.mu       () ,                            
                            pdf.sigma    () ,
                            pdf.zeta     () ,
                            pdf.kappa    () ,
                            pdf.lambd    () )

Ostap.Models.GenHyperbolic.__reduce__ = _rghyp_reduce_ 

# =============================================================================
## reduce Das
def _rdas_reduce_ ( pdf ):
    """Reduce Das"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.mu       () ,                            
                            pdf.sigma    () ,
                            pdf.kL       () ,
                            pdf.kR       () )

Ostap.Models.Das.__reduce__ = _rdas_reduce_ 


# =============================================================================
## reduce GenInvGauss
def _rgig_reduce_ ( pdf ):
    """Reduce GenInvGauss"""
    return root_store_factory , ( type ( pdf )  ,
                            pdf.name      ,
                            pdf.title     ,
                            pdf.x      () ,                            
                            pdf.theta  () ,                            
                            pdf.eta    () ,
                            pdf.p      () ,
                            pdf.shift  () )

Ostap.Models.GenInvGauss.__reduce__ = _rgig_reduce_ 


# =============================================================================
## reduce PositiveSpline 
def _rsplpos_reduce_ ( pdf ):
    """Reduce PositiveSpline"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.spline   () ,                            
                            pdf.phis     () )

Ostap.Models.PositiveSpline   .__reduce__ = _rsplpos_reduce_ 
Ostap.Models.MonotonicSpline  .__reduce__ = _rsplpos_reduce_ 
Ostap.Models.ConvexSpline     .__reduce__ = _rsplpos_reduce_ 
Ostap.Models.ConvexOnlySpline .__reduce__ = _rsplpos_reduce_ 

# =============================================================================
## reduce HORNSdini
def _rdini_reduce_ ( pdf ):
    """Reduce HORNSdini"""
    return root_store_factory , ( type ( pdf )  ,
                            pdf.name      ,
                            pdf.title     ,
                            pdf.x      () ,                            
                            pdf.a      () ,                            
                            pdf.delta  () ,
                            pdf.phi    () )

Ostap.Models.HORNSdini.__reduce__ = _rdini_reduce_ 
Ostap.Models.HILLdini .__reduce__ = _rdini_reduce_ 


# =============================================================================
## reduce Histo1D
def _rh1d_reduce_ ( pdf ) :
    """reduce Histo1D"""
    return root_store_factory , ( type ( pdf )  ,
                            pdf.name      ,
                            pdf.title     ,
                            pdf.x      () ,                            
                            pdf.histo  () )

Ostap.Models.Histo1D.__reduce__ = _rh1d_reduce_ 

# =============================================================================
## reduce Histo2D
def _rh2d_reduce_ ( pdf ) :
    """reduce Histo2D"""
    return root_store_factory , ( type ( pdf )  ,
                            pdf.name      ,
                            pdf.title     ,
                            pdf.x      () ,                            
                            pdf.y      () ,                            
                            pdf.histo  () )

Ostap.Models.Histo2D.__reduce__ = _rh2d_reduce_ 

# =============================================================================
## reduce Histo3D
def _rh3d_reduce_ ( pdf ) :
    """reduce Histo3D"""
    return root_store_factory , ( type ( pdf )  ,
                            pdf.name      ,
                            pdf.title     ,
                            pdf.x      () ,                            
                            pdf.y      () ,                            
                            pdf.z      () ,                            
                            pdf.histo  () )

Ostap.Models.Histo3D.__reduce__ = _rh3d_reduce_ 


# =============================================================================
## reduce CutoffGauss
def _rcutgau_reduce_ ( pdf ):
    """Reduce CutOffGauss"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.x0       () ,                            
                            pdf.sigma    () ,
                            pdf.cutoff   () )

Ostap.Models.CutOffGauss.__reduce__ = _rcutgau_reduce_ 

# =============================================================================
## reduce CutoffStudent
def _rcutstt_reduce_ ( pdf ):
    """Reduce CutOffStudent"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.x0       () ,                            
                            pdf.nu       () ,
                            pdf.sigma    () ,
                            pdf.cutoff   () )

Ostap.Models.CutOffStudent.__reduce__ = _rcutstt_reduce_ 


# =============================================================================
##  2D PDFs
# =============================================================================

# =============================================================================
## reduce Poly2DPositive
def _rpol2d_reduce_ ( pdf ):
    """Reduce Poly2DPositive"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.y        () ,                            
                            pdf.nX       () ,                            
                            pdf.nY       () ,
                            pdf.phis     () )

Ostap.Models.Poly2DPositive.__reduce__ = _rpol2d_reduce_ 


# =============================================================================
## reduce Poly2DSymPositive
def _rpol2ds_reduce_ ( pdf ):
    """Reduce Poly2DSymPositive"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.y        () ,                            
                            pdf.n        () ,                            
                            pdf.phis     () )

Ostap.Models.Poly2DSymPositive.__reduce__ = _rpol2ds_reduce_ 

# =============================================================================
## reduce PS2DPol
def _rps2dpol_reduce_ ( pdf ):
    """Reduce PS2DPol"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.y        () ,
                            pdf.psX      () ,                            
                            pdf.psY      () ,                            
                            pdf.nX       () ,                            
                            pdf.nX       () ,                            
                            pdf.phis     () )

Ostap.Models.PS2DPol.__reduce__ = _rps2dpol_reduce_ 

# =============================================================================
## reduce PS2DPolSym
def _rps2dpols_reduce_ ( pdf ):
    """Reduce PS2DPolSym"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.y        () ,
                            pdf.psX      () ,                            
                            pdf.n        () ,                            
                            pdf.phis     () )

Ostap.Models.PS2DPolSym.__reduce__ = _rps2dpols_reduce_ 


# =============================================================================
## reduce PS2DPol2
def _rps2dpol2_reduce_ ( pdf ):
    """Reduce PS2DPol2"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.y        () ,
                            pdf.psX      () ,                            
                            pdf.psY      () ,                            
                            pdf.mmax     () ,                            
                            pdf.nX       () ,                            
                            pdf.nY       () ,                            
                            pdf.phis     () )

Ostap.Models.PS2DPol2.__reduce__ = _rps2dpol2_reduce_ 
Ostap.Models.PS2DPol3.__reduce__ = _rps2dpol2_reduce_ 

# =============================================================================
## reduce PS2DPol2Sym
def _rps2dpol2s_reduce_ ( pdf ):
    """Reduce PS2DPol2Sym"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.y        () ,
                            pdf.psX      () ,                            
                            pdf.mmax     () ,                            
                            pdf.nX       () ,                            
                            pdf.phis     () )

Ostap.Models.PS2DPol2Sym.__reduce__ = _rps2dpol2s_reduce_ 
Ostap.Models.PS2DPol3Sym.__reduce__ = _rps2dpol2s_reduce_ 

# =============================================================================
## Reduce Expo2DPol
def _rexp2d_reduce_ ( pdf ):
    """Reduce Expo2DPol"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.y        () ,
                            pdf.taux     () ,                            
                            pdf.tauy     () ,                            
                            pdf.nX       () ,                            
                            pdf.nY       () ,                            
                            pdf.phis     () )

Ostap.Models.Expo2DPol.__reduce__ = _rexp2d_reduce_ 


# =============================================================================
## Reduce Expo2DPolSym
def _rexp2ds_reduce_ ( pdf ):
    """Reduce Expo2DPolSym"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.y        () ,
                            pdf.tau      () ,                            
                            pdf.nX       () ,                            
                            pdf.phis     () )

Ostap.Models.Expo2DPolSym.__reduce__ = _rexp2ds_reduce_ 


# =============================================================================
## Reduce ExpoPS2DPol
def _rexpps2d_reduce_ ( pdf ):
    """Reduce ExpoPS2DPol"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.y        () ,
                            pdf.tau      () ,                            
                            pdf.psY      () ,                            
                            pdf.nX       () ,                            
                            pdf.nY       () ,                            
                            pdf.phis     () )

Ostap.Models.ExpoPS2DPol.__reduce__ = _rexpps2d_reduce_ 


# =============================================================================
## reduce Spline2D
def _rspl2d_reduce_ ( pdf ):
    """Reduce Spline2D"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.y        () ,                            
                            pdf.spline   () ,                            
                            pdf.phis     () )

Ostap.Models.Spline2D .__reduce__ = _rspl2d_reduce_ 

# =============================================================================
## reduce Spline2DSym
def _rspl2ds_reduce_ ( pdf ):
    """Reduce Spline2DSym"""
    return root_store_factory , ( type ( pdf )    ,
                            pdf.name        ,
                            pdf.title       ,
                            pdf.x        () ,                            
                            pdf.y        () ,                            
                            pdf.spline   () ,                            
                            pdf.phis     () )

Ostap.Models.Spline2DSym.__reduce__ = _rspl2ds_reduce_ 

# =============================================================================
## reduce Gauss2D
def _rgauss2d_reduce_ ( pdf ):
    """Reduce Gauss2D"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,                            
                                  pdf.muX      () ,                            
                                  pdf.muY      () ,                            
                                  pdf.sigmaX   () ,                            
                                  pdf.sigmaY   () ,                            
                                  pdf.theta    () )

Ostap.Models.Gauss2D.__reduce__ = _rgauss2d_reduce_ 


# =============================================================================
## reduce Poly3DPositive
def _rpol3d_reduce_ ( pdf ):
    """Reduce Poly3DPositive"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,                            
                                  pdf.z        () ,                            
                                  pdf.nX       () ,                            
                                  pdf.nX       () ,                            
                                  pdf.nY       () ,                            
                                  pdf.nZ       () ,                            
                                  pdf.phis     () )

Ostap.Models.Poly3DPositive.__reduce__ = _rpol3d_reduce_ 


# =============================================================================
## reduce Poly3DSymPositive
def _rpol3ds_reduce_ ( pdf ):
    """Reduce Poly3DSymPosiitive"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,                            
                                  pdf.z        () ,                            
                                  pdf.nX       () ,                            
                                  pdf.phis     () )

Ostap.Models.Poly3DSymPositive.__reduce__ = _rpol3ds_reduce_ 


# =============================================================================
## reduce Poly3DMixPositive
def _rpol3dm_reduce_ ( pdf ):
    """Reduce Poly3DMixPosiitive"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,                            
                                  pdf.z        () ,                            
                                  pdf.nX       () ,                            
                                  pdf.nZ       () ,                            
                                  pdf.phis     () )

Ostap.Models.Poly3DMixPositive.__reduce__ = _rpol3dm_reduce_ 


# =============================================================================
## reduce Gauss3D
def _rgauss3d_reduce_ ( pdf ):
    """Reduce Gauss3D"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,                            
                                  pdf.z        () ,                            
                                  pdf.muX      () ,                            
                                  pdf.muY      () ,                            
                                  pdf.muZ      () ,                            
                                  pdf.sigmaX   () ,                            
                                  pdf.sigmaY   () ,                            
                                  pdf.sigmaZ   () ,                            
                                  pdf.phi      () ,
                                  pdf.theta    () ,
                                  pdf.psi      () )

Ostap.Models.Gauss3D.__reduce__ = _rgauss3d_reduce_ 





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
    Ostap.MoreRooFit.Id           ,
    Ostap.MoreRooFit.TwoVars      ,
    Ostap.MoreRooFit.Addition     ,
    Ostap.MoreRooFit.Addition2    ,
    Ostap.MoreRooFit.Subtraction  ,
    Ostap.MoreRooFit.Product      ,
    Ostap.MoreRooFit.ProductPdf   ,
    Ostap.MoreRooFit.Combination  ,
    Ostap.MoreRooFit.Asymmetry    ,
    Ostap.MoreRooFit.Constant     ,
    Ostap.MoreRooFit.Bernstein    ,
    Ostap.MoreRooFit.Monotonic    ,
    Ostap.MoreRooFit.Convex       ,
    Ostap.MoreRooFit.ConvexOnly   ,
    ##
    Ostap.Models.Uniform          ,
    ##
    ROOT.RooArgSet                ,
    ROOT.RooArgList               ,
)

_new_methods_ = tuple ( _new_methods_ ) 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
