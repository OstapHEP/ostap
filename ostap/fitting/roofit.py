#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/roofit.py
#  Module with decoration of some RooFit objects for efficient use in python
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
    'useStorage' , ## define (as context) the default storage for  RooDataStore 
    'PDF_fun'    , ## wrapper of PDF to ``simple'' function 
    'SETVAR'     , ## context manager to preserev the current value for RooRealVar
    ) 
# =============================================================================
import ROOT, random
from   ostap.core.core import cpp, Ostap, VE, hID, dsID , valid_pointer  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger , allright,  attention
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.rootfit' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug( 'Some useful decorations for RooFit objects')
# =============================================================================\
_new_methods_ = []
# =============================================================================\
##  Use RooPrintable::printMultiline function
def print_multiline ( o , content = 1 , verbose = False , indent = '' ) :
    """ Use RooPrintable::printMultiline function
    """
    if not valid_pointer ( o ) : return 'Invalid object'
    return Ostap.Utils.print_printable1 ( o , content , verbose , indent )
# =============================================================================
##  Use RooPrintable::printStream function
def print_stream  ( o , content = 1 , style = 3 , indent = '' ) :
    """ Use RooPrintable::printStream function
    """
    if not valid_pointer ( o ) : return 'Invalid object'
    return Ostap.Utils.print_printable2 ( o , content , style , indent )
# =============================================================================
##  Use RooPrintable::printTree function
def print_tree ( o , indent = '' ) :
    """ Use RooPrintable::printTree function
    """
    if not valid_pointer ( o ) : return 'Invalid object'
    return Ostap.Utils.print_printable_tree ( o , indent )

# =============================================================================
## make easy print for RooPrintable 
def _rp_print_ ( obj , opts = 'vv' , *style ) :
    """Make easy print for RooPrintable
    >>> o = ...
    >>> print o 
    """
    return Ostap.Utils.print_printable ( obj , opts , *style )

ROOT.RooPrintable.print_multiline = print_multiline
ROOT.RooPrintable.print_stream    = print_stream 
ROOT.RooPrintable.print_tree      = print_tree 
ROOT.RooPrintable.print_printable = _rp_print_ 
ROOT.RooPrintable.__str__         = _rp_print_
ROOT.RooPrintable.__repr__        = _rp_print_

_new_methods_ += [
    ROOT.RooPrintable.print_printable ,
    ROOT.RooPrintable.print_multiline ,
    ROOT.RooPrintable.print_stream    ,
    ROOT.RooPrintable.print_tree      ,
    ROOT.RooPrintable.__str__         ,
    ROOT.RooPrintable.__repr__        ,
    ]

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
## check presence of eleemnt or inndex in nthe list
def _ral_contains_ ( self , i ) :
    """Check the presence of element or index in the list
    """
    if isinstance ( i , int ) : return 0<= i < len(self)
    return  0 <= self.index ( i )

# =============================================================================
def _ral_getitem_ ( self , index ) :
    if not isinstance ( index , int ) or not index in self :
        raise IndexError('List Index %s is out of the range [%d,%d)' % ( index , 0 , len(self) ) )
    return self.at ( index  )

# =============================================================================
## some decoration over RooArgList 
ROOT.RooArgList . __len__       = lambda s   : s.getSize()
ROOT.RooArgList . __contains__  = _ral_contains_ 
ROOT.RooArgList . __iter__      = _ral_iter_
ROOT.RooArgList . __nonzero__   = lambda s   : 0 != len ( s ) 
ROOT.RooArgList . __getitem__   = _ral_getitem_
ROOT.RooArgList . __setitem__   = lambda s,*_ : NotImplemented 

_new_methods_ += [
    ROOT.RooArgList. __len__       ,
    ROOT.RooArgList. __contains__  ,
    ROOT.RooArgList. __iter__      ,
    ROOT.RooArgList. __nonzero__   ,
    ROOT.RooArgList. __getitem__   ,
    ROOT.RooArgList. __setitem__   ,
]    

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
    it  = Ostap.Utils.Iterator ( self )
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
        
ROOT.RooArgSet     . __str__   = lambda s : str ( tuple ( _rs_list_ ( s ) ) )  
ROOT.RooArgSet     . __repr__  = lambda s : str ( tuple ( _rs_list_ ( s ) ) )  
ROOT.RooLinkedList . __repr__  = lambda s : str (  _rs_list_ ( s ) )
ROOT.RooLinkedList . __iter__  = _ras_iter_ 

_new_methods_ += [
    ROOT.RooArgSet . __len__      ,
    ROOT.RooArgSet . __iter__     ,
    ROOT.RooArgSet . __getattr__  , 
    ROOT.RooArgSet . __setattr__  ,
    ROOT.RooArgSet . __contains__ ,
    ROOT.RooArgSet . __nonzero__  ,
    ROOT.RooArgSet . __str__      ,
    ROOT.RooArgSet . __repr__     ,
    ]

ROOT.RooLinkedList.add = ROOT.RooLinkedList.Add

_new_methods_ += [
    ROOT.RooLinkedList . __len__  ,
    ROOT.RooLinkedList . __repr__ ,
    ROOT.RooLinkedList . __iter__ ,
    ROOT.RooLinkedList . add 
    ]
    

# =============================================================================
## add more data into list/set
def _ral_iadd_ ( self , other ) :
    """Update/increment collections
    >>> lst = ....
    >>> lst += another_lst
    """
    from collections import Container as _CNT
    _RAC = ROOT.RooAbsCollection
    _RAC = ROOT.RooAbsCollection
    _RAA = ROOT.RooAbsArg      
    if not isinstance ( other , ( _CNT, _RAC , _RAA ) ) : return NotImplemented
    if     isinstance ( other , str )                   : return NotImplemented

    ##
    if isinstance ( other , _RAA ) and not isinstance ( other , _RAC ) : other = [ other ]

    for o in other : self.add ( o )
    return self

# =============================================================================
## add more data into list/set
def _ral_add_ ( self , other ) :
    """Make a sum of two lists/sets/collections
    >>> lst1 = ...
    >>> set2 = ...
    >>> lst2 = lst1 + set2 
    """
    from collections import Container as _CNT
    _RAC = ROOT.RooAbsCollection
    _RAA = ROOT.RooAbsArg  
    if not isinstance ( other , ( _CNT, _RAC , _RAA ) ) : return NotImplemented
    if     isinstance ( other , str )                   : return NotImplemented
    _clone = self.clone('')
    _clone += other
    return _clone

# =============================================================================
## add two list/sets 
def _ral_radd_ ( self , other ) : 
    """Make a sum of two lists/sets/collections
    >>> lst1 = ...
    >>> set2 = ...
    >>> lst2 = lst1 + set2 
    """
    return self + other

# ============================================================================
def _ral_clone_  ( self , name = '' ) :
    return self.Clone(name)


# =============================================================================
for t in ( ROOT.RooArgList , ROOT.RooArgSet , ROOT.RooLinkedList ) :
    t. clone    =  _ral_clone_ 
    t. __add__  =  _ral_add_
    t.__iadd__  = _ral_iadd_
    t.__radd__  = _ral_radd_
    t.append    = _ral_iadd_

_new_methods_ += [
    ROOT.RooArgList    . clone    ,
    ROOT.RooArgSet     . clone    ,
    ROOT.RooLinkedList . clone    ,
    ROOT.RooArgList    . __add__  ,
    ROOT.RooArgSet     . __add__  ,
    ROOT.RooLinkedList . __add__  ,
    ROOT.RooArgList    . __iadd__ ,
    ROOT.RooArgSet     . __iadd__ ,
    ROOT.RooLinkedList . __iadd__ ,
    ROOT.RooArgList    . __radd__ ,
    ROOT.RooArgSet     . __radd__ ,
    ROOT.RooLinkedList . __radd__ ,
    ROOT.RooArgList    . append   ,
    ROOT.RooArgSet     . append   ,
    ROOT.RooLinkedList . append   ,
    ]

    
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
#  @code
#  dataset = ...
#  event   = dataset[4]
#  events  = dataset[0:1000]
#  events  = dataset[0:-1:10]
#  @eendcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-03-31
def _rad_getitem_ ( self , i ) :
    """Get the entry from RooDataSet
    >>> dataset = ...
    >>> event  = dataset[4]
    >>> events = dataset[0:1000]
    >>> events = dataset[0:-1:10]
    """
    if   isinstance ( i , slice ) :
        
        start , stop , step = i.indices ( len ( self ) )
                              
        if 1 == step : return self.reduce ( ROOT.RooFit.EventRange ( start , stop ) )
        
        result = self.emptyClone( dsID() )
        for j in xrange ( start , stop , step ) : result.add ( self [j] ) 
        return result
    
    elif isinstance ( i , ( int , long ) ) and 0<= i < len ( self ) :
        return self.get ( i )
    
    raise IndexError ( 'Invalid index %s'% i )

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
## merge/append two datasets into a single one
# @code
# dset1  = ...
# dset2  = ...
# dset1 += dset2
# @endcode 
def _rad_iadd_ ( self , another ) :
    """ Merge/append two datasets into a single one
    - two datasets must have identical structure 
    >>> dset1  = ...
    >>> dset2  = ...
    >>> dset1 += dset2
    """
    if isinstance ( self , ROOT.RooDataSet ) :
        if isinstance ( another , ROOT.RooDataSet ) :
            self.append ( another )
            return self
        
    return NotImplemented

# =============================================================================
## merge/append two datasets into a single one
#  @code
#  dset1  = ...
#  dset2  = ...
#  dset   = dset1 + dset2 
#  @endcode 
def _rad_add_ ( self , another ) :
    """ Merge/append two datasets into a single one
    - two datasets must have identical structure 
    >>> dset1  = ...
    >>> dset2  = ...
    >>> dset   = dset1 + dset2 
    """
    result = self.emptyClone( dsID() ) 
    result.append ( self    )
    result.append ( another )
    #
    return result 


# =============================================================================
# merge/append two datasets into a single one 
def _rad_imul_ ( self , another ) :
    """ Merge/append two datasets into a single one
    - two datasets must have the  same number of entries!
    >>> dset1  = ...
    >>> dset2  = ...
    >>> dset1 *= dset2
    """
    if  isinstance ( another , ROOT.RooAbsData ) :
        if len ( self ) == len ( another ) :
            self.merge ( another )
            return self
        
    return NotImplemented 

# =============================================================================
## merge two dataset (of same  length) OR get small (random) fraction of  dataset
#  @code
#  ## get smaller dataset:
#  dataset = ....
#  small   = dataset * 0.1
#  ## merge two dataset of the same lenth
#  merged  = dataset1 * dataset2 
#  @endcode
def _rad_mul_ ( self , another ) :
    """
    - (1) Get small (random) fraction of  dataset:
    >>> dataset = ....
    >>> small   = 0.1 * dataset
    - (2) Merge two dataset (of the same length)
    >>> dataset3 = dataset1 * dataset2 
    """

    if isinstance ( another , ROOT.RooAbsData ) :
        
        if len ( self ) == len ( another ) :
            
            result  = self.emptyClone( dsID() )
            result.append ( self    )
            result.merge  ( another )
            return result

        return NotImplemented 
    
    fraction = another    
    if  isinstance ( fraction , float ) and 0 < fraction < 1 :

        res = self.emptyClone()
        l    = len ( self )
        for i in xrange ( l ) :
            if random.uniform(0,1) < fraction : res.add ( self[i] ) 
        return res
    
    elif 1 == fraction : return self.clone      ()
    elif 0 == fraction : return self.emptyClone () 

    return NotImplemented


# =============================================================================
## get small (random) fraction of  dataset
#  @code
#  dataset = ....
#  small   = dataset / 10  
#  @endcode
def  _rad_div_ ( self , fraction ) :
    """ Get small (random) fraction
    >>> dataset = ....
    >>> small   = dataset / 10 
    """
    if  isinstance ( fraction , ( int , long ) ) and 1 < fraction :
        return _rad_mul_ ( self , 1.0 / fraction )
    elif 1 == fraction : return self.clone      ()
    
    return NotImplemented


# =============================================================================
## get small (fixed) fraction of  dataset
#  @code
#  dataset = ....
#  small   = dataset % 10  
#  @endcode
def  _rad_mod_ ( self , fraction ) :
    """ Get small (fixed) fraction of  dataset
    >>> dataset = ....
    >>> small   = dataset % 10 
    """
    if  isinstance ( fraction , ( int , long ) ) and 1 < fraction :

        res = self.emptyClone()
        s    = slice ( 0 , -1 , fraction )
        for i in xrange ( *s.indices ( len ( self ) ) ) : 
            res.add ( self[i] ) 
        return res 
        
    elif 1 == fraction : return self.clone      ()

    return NotImplemented

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

ROOT.RooAbsData . __add__       = _rad_add_
ROOT.RooDataSet . __iadd__      = _rad_iadd_

ROOT.RooAbsData . __mul__       = _rad_mul_
ROOT.RooAbsData . __rmul__      = _rad_mul_
ROOT.RooAbsData . __imul__      = _rad_imul_
ROOT.RooAbsData . __div__       = _rad_div_
ROOT.RooAbsData . __mod__       = _rad_mod_

from ostap.trees.trees import _stat_var_, _stat_cov_ , _stat_covs_ , _sum_var_, _sum_var_old_
ROOT.RooAbsData . sumVar        = _sum_var_ 
ROOT.RooAbsData . sumVar_       = _sum_var_old_ 
ROOT.RooAbsData . statVar       = _stat_var_ 
ROOT.RooAbsData . statCov       = _stat_cov_ 
ROOT.RooAbsData . statCovs      = _stat_covs_ 


_new_methods_ += [
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
   ROOT.RooAbsData . __add__       ,
   ROOT.RooDataSet . __iadd__      ,
   #
   ROOT.RooAbsData . __mul__       ,
   ROOT.RooAbsData . __rmul__      ,
   ROOT.RooAbsData . __imul__      ,
   ROOT.RooAbsData . __div__       ,
   ROOT.RooAbsData . __mod__       ,
   #
   ROOT.RooAbsData . statVar       ,
   ROOT.RooAbsData . sumVar        ,
   ROOT.RooAbsData . sumVar_       ,
   #
   ROOT.RooAbsData . statCov       ,
   ROOT.RooAbsData . statCovs      ,
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


ROOT.RooConstVar    . as_VE          = lambda s : VE( s.getVal() , 0 )
ROOT.RooFormulaVar  . as_VE          = lambda s : VE( s.getVal() , 0 )
ROOT.RooConstVar    . asVE           = lambda s : VE( s.getVal() , 0 )
ROOT.RooFormulaVar  . asVE           = lambda s : VE( s.getVal() , 0 )


ROOT.RooAbsReal       .__contains__ = lambda s,v : False ## ??? do we need it???
ROOT.RooAbsRealLValue .__contains__ = _rrv_contains_ 

# =====================================================================
ROOT.RooAbsReal. minmax  = lambda s : ()
ROOT.RooAbsReal.xminmax  = lambda s : ()
ROOT.RooAbsRealLValue  . xmin            = lambda s : s.getMin()
ROOT.RooAbsRealLValue  . xmax            = lambda s : s.getMax()
ROOT.RooAbsRealLValue  . minmax          = lambda s : (s.xmin(),s.xmax()) 
ROOT.RooAbsRealLValue  .xminmax          = lambda s : (s.xmin(),s.xmax()) 


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
    return Ostap.Models.Product (
        '%s*%s'             % ( pdf1.GetName  () , pdf2.GetName  () ) ,
        'Product: %s & %s ' % ( pdf1.GetTitle () , pdf2.GetTitle () ) ,
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
    if    0 <= w.find('(') < what.find(')') : pass
    elif  0 <  w.find('*')                  : pass
    elif  0 <  w.find('/')                  : pass
    elif  0 <  w.find('%')                  : pass 
    elif  0 <  w.find('+')                  : pass
    elif  0 <  w.find('-')                  : pass
    else :
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
    if isinstance ( what , ROOT.RooArgList ) and \
       isinstance ( histo , ROOT.TH1       ) and \
       hasattr ( dataset , 'fillHistogram' ) :
        histo.Reset() 
        return dataset.fillHistogram  ( histo , what , cuts , *args )
    
    ## delegate to TTree (only for non-weighted dataset with TTree-based storage type) 
    if hasattr ( dataset , 'isWeighted') and not dataset.isWeighted() \
       and isinstance ( what , str ) \
       and isinstance ( cuts , str ) :
        if hasattr ( dataset , 'store' ) : 
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
        return Ostap.HistoProject.project3 ( dataset ,
                                                histo   , 
                                                what[2] ,
                                                what[1] ,
                                                what[0] , cuts , *args) 
    elif isinstance ( histo , ROOT.TH2 ) and 2 == len(what)  :
        return Ostap.HistoProject.project2 ( dataset ,
                                                 histo   , 
                                                 what[1] ,
                                                 what[0] , cuts , *args )
    elif isinstance ( histo , ROOT.TH1 ) and 1 == len(what)  :
        return Ostap.HistoProject.project  ( dataset ,
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
    if hasattr ( dataset , 'isWeighted') and not dataset.isWeighted() \
       and isinstance ( what , str ) \
       and isinstance ( cuts , str ) \
       and isinstance ( opts , str ) :
        if hasattr ( dataset , 'store' ) : 
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
## clear dataset storage
if not hasattr ( ROOT.RooDataSet , '_old_reset_' ) :
    ROOT.RooDataSet._old_reset_ = ROOT.RooDataSet.reset
    def _ds_new_reset_ ( self ) :
        """Clear dataset storage
        >>> print ds
        >>> ds.clear()
        >>> ds.erase() ## ditto
        >>> ds.reset() ## ditto
        >>> ds.Reset() ## ditto
        >>> print ds
        """
        s = self.store()
        if s : s.reset()
        self._old_reset_()
        return len(self)
    ROOT.RooDataSet.reset = _ds_new_reset_

ROOT.RooDataSet.clear = ROOT.RooDataSet.reset
ROOT.RooDataSet.erase = ROOT.RooDataSet.reset
ROOT.RooDataSet.Reset = ROOT.RooDataSet.reset

# =============================================================================
ROOT.RooDataSet.draw        = _ds_draw_
ROOT.RooDataSet.project     = _ds_project_
ROOT.RooDataSet.__getattr__ = _ds_getattr_

ROOT.RooDataHist.__len__    = lambda s : s.numEntries() 


# =============================================================================
## print method for RooDataSet
#  @code
#
#   >>> print dataset
#
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-07-06
def _ds_print_ ( dataset ) :
    """Helper print method:    
    >>> print dataset 
    """
    if not  valid_pointer ( dataset ) : return 'Invalid dataset'
    return dataset.print_multiline ( verbose = True ) 

ROOT.RooDataSet.draw        = _ds_draw_
ROOT.RooDataSet.project     = _ds_project_
ROOT.RooDataSet.__getattr__ = _ds_getattr_

for d in ( ROOT.RooAbsData  ,
           ROOT.RooDataSet  ,
           ROOT.RooDataHist ) :
    d.__repr__    = _ds_print_
    d.__str__     = _ds_print_
    d.__len__     = lambda s : s.numEntries() 

_new_methods_ += [
    ROOT.RooDataSet .draw         ,
    ROOT.RooDataSet .project      ,
    ROOT.RooDataSet .__getattr__  ,
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
## @class UseStorage
#  Context manager to change the storage type
class UseStorage(object) :
    """Context manager to change the storage type
    >>> with UseStorage() :
    ...
    """
    def __init__  ( self , new_storage = RAD.Tree ) :
        if not new_storage in ( RAD.Tree , RAD.Vector )  :
            raise AttributeError( 'Invalid storage type %s' % new_storage )
        self.new_storage = new_storage
        self.old_storage = RAD.getDefaultStorageType()
    def __enter__ ( self ) :
        self.old_storage = RAD.getDefaultStorageType()
        setStorage (  self.new_storage )
    def __exit__ (  self , *_ ) :
        setStorage (  self.old_storage )

# =============================================================================
## context manager to change the storage type
def useStorage ( storage = RAD.Tree ) :
    """Context manager to change the storage type
    >>> with useStorage() :
    ...
    """
    return UseStorage ( storage )


# =============================================================================
## fitting
#  @code
#  model = ...
#  data  = ...
#  data.fitTo ( model , ... )
#  data.Fit   ( model , ... ) ## ditto
#  @endcode
def _rad_fit_ ( data , model , *args , **kwargs ) :
    """ fitting
    >>> model = ...
    >>> data  = ...
    >>> data.fitTo ( model , ... )
    >>> data.Fit   ( model , ... ) ## ditto 
    """
    return model.fitTo ( data , *args , **kwargs )

RAD.Fit   = _rad_fit_
RAD.fitTo = _rad_fit_

_new_methods_ += [
    RAD.Fit ,
    RAD.fitTo
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
    >>> from ostap.stats.moments import mean, mode, median  
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
    return VE ( v , e*e ) if e>0 else v

# =============================================================================
def _rav_setval_  ( self , value ) :
    """Assign the valeu for the variable 
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
    mn,mx  = self.getMin(), self.getMax() 
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
def _rad_moment_ ( data , var , order , value = 0 , error = True , *args ) :
    """ Get n-th moment of the distribution
    >>> data = ...
    >>> print data.moment ( 'mass' , 3 ) 
    """
    assert isinstance ( order , int ) and 0 <= order, 'Invalid "order"  %s' % order
    
    if isintance  ( var  , str ) :
        varset =  data.get()
        assert  var in varset, 'Invalid varibale %s' % var 
        var = getarrt ( varset , var ) 
        return _rad_moment_  ( data , var  , order , value , error , *args )

    m  = data._old_moment_ ( var , order , value , *args )
    if not error : return m
    
    n     = data.sumEntries( *args ) 
    sigma = data.sigma ( var , *args )
    
    if  abs  ( value - 0 ) < 0.01 * sigma :
        
        m2  = data._old_moment ( var , 2 * order , *args )
        c2  = ( m2  - m * m )
        c2 /= n
        
    elif  abs  ( value - data._old_moment_ ( var  , 1 , 0  , *args  ) ) < 0.01 * sigma :

        m2  = data._old_moment_ ( var , 2             , value , *args )
        m2o = data._old_moment_ ( var , 2 * order     , value , *args )
        mum = data._old_moment_ ( var ,     order - 1 , value , *args )
        mup = data._old_moment_ ( var ,     order + 1 , value , *args )

        c2  = m2o
        c2 -= 2.0 * order * mum * mup
        c2 -= m * m
        c2 += order  * order * m2 * mum * mup
        c2 /= n
        
    else  :
        logger.error ("Uncertainnry can be calcualted onlyfro moment/central moment") 
        
    return VE ( m , c2 )


# ===============================================================================
## Get n-th central moment of the distribution
#  @code
#  data = ...
#  print data.central_moment ( 'mass' , 3 )
#  @endcode 
def _rad_central_moment_ ( data , var  , order , error = True , *args  ) :
    """ Get n-th central moment of the distribution
    >>> data = ...
    >>> print data.central_moment ( 'mass' , 3 ) 
    """
    ##  get the men-value:
    mu = _rad_moment_  ( data , var  , 1 , error = False , *args )
    ##  calcualte moments  
    return _rad_moment_ ( data , var , order , mu  , error , *args )


# =============================================================================
def _rad_skewness_ ( data , var , error = True , *args ) :
    
    if isintance  ( var  , str ) :
        varset =  data.get()
        assert var in varset, 'Invalid variable %s' % var  
        var = getarrt ( varset , var ) 
        return _rad_skewness_  ( data , var , error , *args )
    
    s  = data._old_skewness_ ( var , *args )
    if not error : return s
    
    n = dat.sumEntries( *args ) 

    if 2 > n : return VE ( s , 0 )

    c2  = 6.0
    c2 *= ( n - 2 )
    c2 /= ( n + 1 ) * (  n + 3 )

    return VE ( s , c2  )

# =============================================================================
def _rad_kurtosis_ ( data , var , error = True , *args ) :
    
    if isintance  ( var  , str ) :
        varset =  data.get()
        assert var in varset, 'Invalid variable %s' % var 
        var = getarrt ( varset , var ) 
        return _rad_kurtisis_  ( data , var , error , *args )

    k  = data._old_kurtosis_ ( var , *args )
    if not error : return k
    
    n = dat.sumEntries( *args ) 
    
    if 3 > n : return VE ( k , 0 )

    c2 = 24.0 * n 
    c2 *= ( n - 2 ) * ( n - 3 )
    c2 /= ( n + 1 ) * ( n + 1 )
    c2 /= ( n + 3 ) * ( n + 5 )
    
    return VE  ( k , c2 )


RAD  = ROOT.RooAbsData
if  not hasattr ( RAD , '_new_moment_' ) :
    RAD.__old_moment_   = RAD.moment
    RAD.__new_moment_   = _rad_moment_
    RAD.moment          = _rad_moment_
    
if  not hasattr ( RAD , '_new_skewness_' ) :
    RAD.__old_skewness_ = RAD.skewness
    RAD.__new_skewness_ = _rad_skewness_
    RAD.skewness        = _rad_skewness_

if  not hasattr ( RAD , '_new_kurtosis_' ) :
    RAD.__old_kurtosis_ = RAD.kurtosis
    RAD.__new_kurtosis_ = _rad_kurtosis_
    RAD.kurtosis        = _rad_kurtosis_

RAD.central_moment = _rad_central_moment_

# ==============================================================================
## get the list/tuple of variable names 
#  @code
#  data = ...
#  br1 = data.branches() 
#  br2 = data.branches('.*(Muon).*'   , re.I ) 
#  br3 = data.branches('.*(Probnn).*' , re.I ) 
#  @endcode
def _rad_branches_ (  self , pattern = '' , *args ) :
    """Get the list/tuple of variable names 
    >>> data = ...
    >>> br1 = data.branches() 
    >>> br2 = data.branches('.*(Muon).*'   , re.I ) 
    >>> br3 = data.branches('.*(Probnn).*' , re.I )
    >>> br1 = data.leaves  () 
    >>> br2 = data.leaves  ('.*(Muon).*'   , re.I ) 
    >>> br3 = data.leaves  ('.*(Probnn).*' , re.I )
    """
    
    vlst = self.varset()
    if not vlst : return tuple()

    if pattern :        
        try : 
            import re
            c  =  re.compile ( pattern , *args )
            lst  = [ v.GetName() for v in vlst if c.match ( v.GetName () ) ]
            lst.sort()
            return tuple ( lst ) 
        except :
            logger.error ('branches: exception is caught, skip it' , exc_info = True ) 
            
    lst  = [ v.GetName() for v in vlst  ]
    lst.sort()
    return tuple ( lst ) 


RAD.branches = _rad_branches_
RAD.leaves   = _rad_branches_

# ==============================================================================
def _ds_table_0_ ( dataset , variables = [] , cuts = '' , first = 0 , last = 2**62 ) :
    """Print data set as table
    """
    varset = dataset.get()
    if not valid_pointer ( varset ) :
        logger.error('Invalid dataset')
        return ''

    if isinstance ( variables ,  str ) :
        variables = variables.strip ()
        variables = variables.replace ( ',' , ' ' ) 
        variables = variables.replace ( ';' , ' ' )
        variables = variables.split ()
        
    if 1 == len ( variables ) : variables = variables [0]

    if isinstance ( variables ,  str ) :
        
        if variables in varset :
            vars = [ variables ]
        else :
            vars = list ( dataset.branches ( variables ) ) 
            
    elif variables : vars = [ i.GetName() for i in varset if i in variables ]        
    else           : vars = [ i.GetName() for i in varset                   ]
        
    #
    _vars = []
    for v in vars :
        vv   = getattr ( varset , v ) 
        s    = dataset.statVar( v , cuts , first , last )  
        mnmx = s.minmax ()
        mean = s.mean   ()
        rms  = s.rms    ()
        r    = ( vv.GetName  () ,                      ## 0 
                 vv.GetTitle () ,                      ## 1 
                 ('%+.5g' % mean.value() ).strip() ,   ## 2
                 ('%.5g'  % rms          ).strip() ,   ## 3 
                 ('%+.5g' % mnmx[0]      ).strip() ,   ## 4
                 ('%+.5g' % mnmx[1]      ).strip() )   ## 5
            
        _vars.append ( r )
        
    _vars.sort()

    report  = '# %s("%s","%s"):' % ( dataset.__class__.__name__ ,
                                     dataset.GetName  () ,
                                     dataset.GetTitle () )
    report += allright ( '%d entries, %d variables' %  ( len ( dataset )   ,
                                                         len ( varset  ) ) )

    if not _vars :
        return report , 120 


    weight = None 
    if dataset.isWeighted() :
        report += attention ( ' Weighted' )

        dstmp = None 
        wvar  = None
        
        ## 1) try to get the name of the weight variable
        store = dataset.store()
        if not valid_pointer ( store ) : store = None
        if store and not isinstance ( store , ROOT.RooTreeDataStore ) :
            dstmp = dataset.emptyClone ()
            dstmp.convertToTreeStore   ()
            store = dstmp.store        ()
            if not valid_pointer ( store ) : store = None
        if hasattr ( store , 'tree' ) and valid_pointer ( store.tree() ) : 
            tree = store.tree() 
            branches = set ( tree.branches() )
            vvars    = set ( [ i.GetName() for i in  varset ] )
            wvars    = branches - vvars
            if 1 == len ( wvars ):
                wvar = wvars.pop() 
                report += attention ( ' with "%s"' % wvar )
                
        store = None 
        if not dstmp is None :            
            dstmp.reset()            
            del dstmp
            dstmp = None 

        ## 2) if weight name is known, try to get information about the weight
        if wvar :
            store = dataset.store()
            if not valid_pointer ( store ) : store = None
            if store and not isinstance ( store , ROOT.RooTreeDataStore ) :

                rargs = ROOT.RooFit.EventRange ( first , last ) , 
                if cuts :
                    ## need all variables 
                    dstmp = dataset.reduce ( ROOT.RooFit.Cut  ( cuts ) , *rargs ) 
                else    :
                    ## enough to keep only 1 variable
                    vvs   = ROOT.RooArgSet ( varset[vars[0]] )
                    dstmp = dataset.reduce ( ROOT.RooFit.SelectVars ( vvs ) , *rargs )

                dstmp.convertToTreeStore ()
                store = dstmp.store()
                cuts , first , last = '' , 0 , 2**62
                
            if hasattr ( store , 'tree' ) and valid_pointer ( store.tree() ) : 
                tree =  store.tree() 
                s = tree.statVar ( wvar , cuts , first , last ) ## no cuts here... 
                mnmx = s.minmax ()
                mean = s.mean   ()
                rms  = s.rms    ()
                weight = '*%s*' % wvar
                r    = (  weight                           ,   ## 0 
                         'Weight variable'                 ,   ## 1 
                         ('%+.5g' % mean.value() ).strip() ,   ## 2
                         ('%.5g'  % rms          ).strip() ,   ## 3 
                         ('%+.5g' % mnmx[0]      ).strip() ,   ## 4
                         ('%+.5g' % mnmx[1]      ).strip() )   ## 5
                _vars.append ( r ) 
                with_weight = True
                
            store = None 
            if not dstmp is None :
                dstmp.reset ()                
                del dstmp
                dstmp = None 

    # ==============================================================================================
    # build the actual table 
    # ==============================================================================================
    
    name_l  = len ( 'Variable'    ) + 2 
    desc_l  = len ( 'Description' ) + 2 
    mean_l  = len ( 'mean' ) + 2 
    rms_l   = len ( 'rms'  ) + 2
    min_l   = len ( 'min'  ) + 2 
    max_l   = len ( 'max'  ) + 2 
    for v in _vars :
        name_l = max ( name_l , len ( v[0] ) )
        desc_l = max ( desc_l , len ( v[1] ) )
        mean_l = max ( mean_l , len ( v[2] ) )
        rms_l  = max ( rms_l  , len ( v[3] ) )
        min_l  = max ( min_l  , len ( v[4] ) )
        max_l  = max ( max_l  , len ( v[5] ) )
        
    sep      = '# +%s+%s+%s+%s+' % ( ( name_l       + 2 ) * '-' ,
                                     ( desc_l       + 2 ) * '-' ,
                                     ( mean_l+rms_l + 5 ) * '-' ,
                                     ( min_l +max_l + 5 ) * '-' )
    fmt = '# | %%-%ds | %%-%ds | %%%ds / %%-%ds | %%%ds / %%-%ds |'  % (
        name_l ,
        desc_l ,
        mean_l ,
        rms_l  ,
        min_l  ,
        max_l  )
    
                
    header  = fmt % ( 'Variable'    ,
                      'Description' ,
                      'mean'        ,
                      'rms'         ,
                      'min'         ,
                      'max'         )
    
    report += '\n' + sep
    report += '\n' + header
    report += '\n' + sep

    vlst   = _vars
    
    if weight : vlst = _vars[:-1]
    
    for v in vlst :
        line    =  fmt % ( v[0] , v[1] , v[2] , v[3] , v[4] , v[5]  )
        report += '\n' + line  
    report += '\n' + sep
    
    if weight :
        v = _vars[-1]
        line    =  fmt % ( v[0] , v[1] , v[2] , v[3] , v[4] , v[5]  )
        report += '\n' + line.replace ( weight , attention ( weight ) ) 
        report += '\n' + sep
        
    return report , len ( sep ) 


# ==============================================================================
## print dataset in  a form of the table
#  @code
#  dataset = ...
#  print dataset.table() 
#  @endcode
def _ds_table_ (  dataset ,  variables = [] ) :
    """print dataset in a form of the table
    >>> dataset = ...
    >>> print dataset.table()
    """
    return _ds_table_0_ ( dataset ,  variables )[0]

# =============================================================================
##  print DataSet
def _ds_print2_ ( dataset ) :
    """Print dataset"""
    if dataset.isWeighted() :
        store = dataset.store()
        if valid_pointer ( store ) and isinstance ( store , ROOT.RooTreeDataStore ) : pass
        else : return _ds_print_ ( dataset )        
    from ostap.utils.basic import terminal_size, isatty 
    if not isatty() : return _ds_table_ ( dataset )
    th  , tw  = terminal_size()
    rep , wid = _ds_table_0_ ( dataset ) 
    if wid < tw     : return rep
    return _ds_print_ ( dataset )

for t in ( ROOT.RooDataSet , ) :
    t.__repr__    = _ds_print2_
    t.__str__     = _ds_print2_
    t.table       = _ds_table_
    t.pprint      = _ds_print_ 
    


# =============================================================================
from  ostap.stats.statvars import data_decorate as _dd
_dd ( ROOT.RooAbsData )

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
