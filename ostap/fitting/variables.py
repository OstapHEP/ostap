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
    ## ``converters''
    'total_ratio'     , ## ``converter'': A,B ->  (T,R) == ( A+B , A/B     )
    'total_ratio'     , ## ``converter'': A,B ->  (T,F) == ( A+B , A/(A+B) ) 
    'two_yields'      , ## ``converter'': T,F ->  (A,B) == ( R*F , T*(1-F) )
    'depends_on'      , ## Is this "RooFit" function depends on the variable?
    'binning'         , ## create RooBinning object 
    ) 
# =============================================================================
import ROOT, random, array 
from   builtins               import range
from   ostap.core.core        import VE, hID, Ostap
from   ostap.core.meta_info   import root_info 
from   ostap.core.ostap_types import ( num_types     , list_types   ,
                                       integer_types , string_types )

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
## Addition of RooRealVar and ``number''
def _rrv_add_ ( s , o ) :
    """Addition of RooRealVar and ``number''

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
## Subtraction  of RooRealVar and ``number''
def _rrv_sub_ ( s , o ) :
    """Subtraction of RooRealVar and ``number''
    
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
## Multiplication of RooRealVar and ``number''
def _rrv_mul_ ( s , o ) :
    """Multiplication  of RooRealVar and ``number''
    
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
## Division of RooRealVar and ``number''
def _rrv_div_ ( s , o ) :
    """Division of RooRealVar and ``number''
    
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
## (right) Addition of RooRealVar and ``number''
def _rrv_radd_ ( s , o ) :
    """(right) Addition of RooRealVar and ``number''
    
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
## (right) subtraction  of RooRealVar and ``number''
def _rrv_rsub_ ( s , o ) :
    """(right) subtraction of RooRealVar and ``number''
    
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
## (right) multiplication of RooRealVar and ``number''
def _rrv_rmul_ ( s , o ) :
    """(right) Multiplication  of RooRealVar and ``number''
    
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
## (right) Division of RooRealVar and ``number''
def _rrv_rdiv_ ( s , o ) :
    """(right) Division of RooRealVar and ``number''
    
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
## pow of RooRealVar and ``number''
def _rrv_pow_ ( s , o ) :
    """pow of RooRealVar and ``number''
    
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
## (right) pow of RooRealVar and ``number''
def _rrv_rpow_ ( s , o ) :
    """pow of RooRealVar and ``number''
    
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

ROOT.RooRealVar . __iadd__  = _rrv_iadd_
ROOT.RooRealVar . __isub__  = _rrv_isub_
ROOT.RooRealVar . __idiv__  = _rrv_idiv_

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


# =============================================================================
## Properties 
# =============================================================================

# =============================================================================
def _rav_getval_  ( self ) :
    """Get the value, associated with the variable
    >>> var = ...
    >>> print var.value 
    """
    return self.getVal()

# =============================================================================
def _rav_getvale_ ( self ) :
    """Get the value(and the error), associated with the variable
    >>> var = ...
    >>> print  var.value 
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
    """Convert two yields into ``total yield'' and ``ratio''
    
    >>> yield1 = ROOT.RooRelaVar( ... )
    >>> yield2 = ROOT.RooRelaVar( ... )
    >>> total , ratio = total_ratio (  yield1 , yield2 ) 
    """
    
    assert isinstance ( var1 , ROOT.RooAbsReal,\
                        "Invalid type of ``var1'' %s/%s" % ( var1 , type ( var1 ) ) )
    assert isinstance ( var2 , ROOT.RooAbsReal,\
                        "Invalid type of ``var2'' %s/%s" % ( var2 , type ( var2 ) ) )
    
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
    """Convert two yields into ``total yield'' and ``fraction''
    >>> yield1 = ROOT.RooRelaVar( ... )
    >>> yield2 = ROOT.RooRelaVar( ... )
    >>> total , fraction = total_fraction (  yield1 , yield2 ) 
    """
    
    assert isinstance ( var1 , ROOT.RooAbsReal,\
                        "Invalid type of ``var1'' %s/%s" % ( var1 , type ( var1 ) ) )
    assert isinstance ( var2 , ROOT.RooAbsReal,\
                        "Invalid type of ``var2'' %s/%s" % ( var2 , type ( var2 ) ) )
    
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
        """``variables'' : list of variables"""
        return self.__variables

    @property
    def fixed ( self ) :
        """``fixed'' : fixed variable?"""
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
        """``bins'' : the binning scheme"""
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
           "Invalid ``nbins'' parameter %s/%s" % ( nbins , type ( nbins ) )

    nb = len ( edges ) 
    assert 2 <= nb , "Invalid length of ``edges'' array!"

    if 2 == nb :
        return ROOT.RooBinning ( max ( 1 , nbins ) , edges[0] , edges[1] , name ) 

    buffer = array.array ( 'd', edges )
    return ROOT.RooBinning ( nb - 1 , buffer , name ) 


    
# =============================================================================
## Dedicated unpickling factory for RooRealVar
def rrv_factory ( args , errors , binnings , fixed ) :
    """ Dedicated unpickling factory for `ROOT.RooRealVar`
    """
    ## create it
    rrv = ROOT.RooRealVar ( *args )

    ## set errors if needed 
    if errors :
        if   2 == len ( errors ) : rrv.setAsymError ( *errors ) 
        elif 1 == len ( errors ) : rrv.setError     ( *errors )

    for b in binnings :
        nb = b.GetName() 
        if nb : rrv.setBinning ( b , nb ) 
        
    rrv.setConstant ( fixed )
    
    return rrv

# =============================================================================
## Reducing of <code>RooRealVar</code> for pickling/unpickling 
#  @see RooRooRealVar 
def rrv_reduce ( rrv ) :
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

    ## bining schemes 
    bnames   = rrv.getBinningNames()
    binnings = [] 
    for nb in bnames:
        if not nb : continue
        bb = rrv.getBinning ( nb , False ) 
        if bb.GetName () : binnings.append ( bb )
        
    binnings = tuple ( binnings ) 

    ## fixed ? 
    fixed = True if rrv.isConstant() else False 
    
    content = args , errors , binnings , fixed
    
    return rrv_factory , content 


ROOT.RooRealVar.__reduce__ = rrv_reduce


# =============================================================================
## factory for unpickling of <code>RooFormulaVar</code> and
#  <code>Ostap::FormulaVar</code>
def rfv_factory ( klass , args , vars ) :
    """Factory for unpickling of `RooFormulaVar` and `Ostap.FormulaVar`
    """

    lst = ROOT.RooArgList ()
    for v in vars : lst.add ( v ) 

    margs = list  ( args  )
    margs.append  ( lst   )
    margs = tuple ( margs ) 
    
    rfv = klass ( *margs )
    
    rfv.__vlst = lst
    
    return rfv
    
# =============================================================================
## Reduce <code>RooFormulaVar</code> and <code>Ostap::FormulaVar</code> for pickling
#  @see RooFormulaVar 
#  @see Ostap::FormulaVar 
def rfv_reduce ( rfv ) : 
    """Reduce `RooFormulaVar` and `Ostap::FormulaVar` for pickling
    - see RooFormulaVar 
    - see Ostap.FormulaVar 
    """

    name       = rfv.GetName  ()
    title      = rfv.GetTitle ()
    
    rform      = rfv.formula  ()
    
    expression = rform.GetTitle() 

    vars       = []
    deps       = rform.actualDependents()
    for d in deps : vars.append ( d )
    vars = tuple ( vars ) 

    args  = name , title , expression 
    
    return rfv_factory , ( type ( rfv ) , args , vars ) 


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
    ROOT.RooFormulaVar.__reduce__  = rfv_reduce
else :
    Ostap.FormulaVar.__reduce__    = rfv_reduce
    Ostap.FormulaVar. __str__      = _rfv_str_
    Ostap.FormulaVar. __repr__     = _rfv_repr_
    Ostap.FormulaVar.__reduce__    = rfv_reduce


# =============================================================================
## unpickle <code>Ostap::MoreFooFit::TwoVars</code> objects
#  @see Ostap::MoreRooFit.TwoVars
def r2v_factory ( klass , *args ) :
    """unpickle `Ostap::MoreFooFit::TwoVars` objects
    - see Ostap.MoreRooFit.TwoVars
    """
    return klass ( *args )

# =============================================================================
## Reduce <code>Ostap::MoreFooFit::TwoVars</code> objects
#  @see Ostap::MoreRooFit.TwoVars
def r2v_reduce ( vars ) :
    """Reduce `Ostap::MoreFooFit::TwoVars` objects
    - see Ostap.MoreRooFit.TwoVars
    """
    return r2v_factory , ( type ( var ) , var.name , var.title , var.x() , var.y() )


Ostap.MoreRooFit.TwoVars.__reduce__  = r2v_reduce


# ===================================================================
## Reduce <code>Ostap::MoreFooFit::Combination</code> objects
#  @see Ostap::MoreRooFit.Combination
def _rc_reduce ( vars ) :
    """Reduce `Ostap.MoreFooFit.Combination` objects
    - see Ostap.MoreRooFit.Combination
    """
    return r2v_factory , ( type ( var ) ,
                           var.name , var.title ,
                           var.x()  , var.y()   ,
                           var.alpha() , var.beta() , var.gamma () )

# ===================================================================
## Reduce <code>Ostap::MoreFooFit::Asymmetry</code> objects
#  @see Ostap::MoreRooFit.Asymmetry
def _ra_reduce ( vars ) :
    """Reduce  `Ostap::MoreFooFit::Asymmetry` objects
    - see Ostap.MoreRooFit.Asymmetry
    """
    return r2v_factory , ( type ( var ) ,
                           var.name , var.title ,
                           var.x()  , var.y()   ,
                           var.scale () )


Ostap.MoreRooFit.Combination.__reduce__  = _rc_reduce
Ostap.MoreRooFit.Asymmetry.  __reduce__  = _ra_reduce



# =============================================================================
## Reduce <code>RooConstVar</code>
#  @see RooConstVar 
def rcv_reduce ( var ) :
    """ Reduce `ROOT.RooConstVar`
    - see ROOT.RooConstVar
    """
    return r2v_factory , ( type ( var )  ,
                           var.name      ,
                           var.title     ,
                           float ( var ) )  

ROOT.RooConstVar.__reduce__ = rcv_reduce 


# =============================================================================
## unpickle RooUniformBinning object
#  @see RooUniformBinnig  
def rub_factory ( *args ) :
    """unpickle RooUniformBinning object
    -see ROOT.RooUniformBinnig
    """
    print ( 'reconstruct U-bining', args ) 
    return ROOT.RooUniformBining ( *args )

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
    return rub_factory, content


# =============================================================================
## unpickle RooBinning object
#  @see RooBinnig  
def rb_factory ( data , name  ) :
    """unpickle RooBinning object
    -see ROOT.RooUniformBinnig
    """
    return ROOT.RooBining ( len(data) -1 , data[0] , name )
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
    return rb_factory, content


ROOT.RooBinning       .__reduce__ = _rb_reduce_
ROOT.RooUniformBinning.__reduce__ = _rub_reduce_
        
# =============================================================================


_decorated_classes_ = (
    ROOT.RooRealVar        ,
    ROOT.RooConstVar       ,
    ROOT.RooFormulaVar     ,
    ROOT.RooAbsReal        ,
    ROOT.RooAbsRealLValue  ,
    ROOT.RooUniformBinning ,
)

_new_methods_ = tuple ( _new_methods_ ) 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
