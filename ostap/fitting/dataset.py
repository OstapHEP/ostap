#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/dataset.py
#  Module with decoration for RooAbsData and related RooFit classes
#  @see RooAbsData 
#  @see RooDataSet 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Module with decoration for RooAbsData and related RooFit classes
- see RooAbsData 
- see RooDataSet 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'setStorage' , ## define the default storage for  RooDataStore 
    'useStorage' , ## define (as context) the default storage for  RooDataStore
    'ds_draw'    , ## draw varibales from RooDataSet 
    'ds_project' , ## project variables from RooDataSet to histogram 
    )
# =============================================================================
import ROOT, random, math, sys 
from   builtins               import range
from   ostap.core.core        import Ostap, VE, hID, dsID , valid_pointer
from   ostap.core.ostap_types import integer_types, string_types  
import ostap.fitting.variables 
import ostap.fitting.roocollections
import ostap.fitting.printable
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.dataset' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug( 'Some useful decorations for RooAbsData object')
# =============================================================================
from ostap.logger.colorized import allright,  attention
_new_methods_ = []
# =============================================================================
_maxv =  0.99 * sys.float_info.max
_minv = -0.99 * sys.float_info.max
# =============================================================================

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
    for i in range ( 0 , _l ) :
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
        for j in range ( start , stop , step ) : result.add ( self [j] ) 
        return result
    
    elif isinstance ( i , integer_types ) and 0<= i < len ( self ) :
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
    if isinstance ( self , ROOT.RooDataSet ) :
        if isinstance ( another , ROOT.RooDataSet ) :    
            result = self.emptyClone( dsID() ) 
            result.append ( self    )
            result.append ( another )
            return result
    
    return NotImplemented


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
        for i in range ( l ) :
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
    if  isinstance ( fraction , integer_types ) and 1 < fraction :
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
    if  isinstance ( fraction , integer_types ) and 1 < fraction :

        res = self.emptyClone()
        s    = slice ( 0 , -1 , fraction )
        for i in range ( *s.indices ( len ( self ) ) ) : 
            res.add ( self[i] ) 
        return res 
        
    elif 1 == fraction : return self.clone      ()

    return NotImplemented


# =============================================================================
## get (random) sub-sample from the dataset
#  @code
#  data   = ...
#  subset =  data.sample ( 100  )  ## get 100   events 
#  subset =  data.sample ( 0.01 )  ## get 1% of events 
#  @endcode 
def _rad_sample_ ( self , num ) :
    """Get (random) sub-sample from the dataset
    >>> data   = ...
    >>> subset =  data.sample ( 100  )  ## get 100   events 
    >>> subset =  data.sample ( 0.01 )  ## get 1% of events 
    """
    if   0 == num : return self.emptyClone ( dsID () ) 
    elif isinstance ( num , integer_types ) and 0 < num :
        num = min ( num , len ( self ) )
    elif isinstance ( num , float ) and 0 < num < 1 :
        from ostap.math.random_ext import poisson 
        num = poisson ( num * len ( self ) )
        return _rad_sample_ ( self , num )
    else :
        raise TypeError("Unknown ``num''=%s" % num )
    
    result  = self.emptyClone ( dsID () )
    indices = random.sample ( range ( len ( self ) ) , num )
    
    while indices :
        i = indices.pop()
        result.add ( self[i] )
        
    return result 

# =============================================================================
## get the shuffled sample
#  @code
#  data = ....
#  shuffled = data.shuffle()
#  @endcode 
def _rad_shuffle_ ( self ) :
    """Get the shuffled sample
    >>> data = ....
    >>> shuffled = data.shuffle()
    """
    result  = self.emptyClone ( dsID () )
    
    indices = [ i for i in range( len ( self ) ) ]  
    random.shuffle ( indices )

    while indices :
        i = indices.pop()
        result.add ( self[i] )
        
    return result 

# ==============================================================================
## Imporved reduce
#  @code
#  data = ...
#  data =  data.subset ( RooArgSet( ... ) , 'a>0' )
#  data =  data.subset ( [a,b,c]       , 'a>0' )
#  data =  data.subset ( ['a','b','c'] , 'a>0' )
#  @endcode
#  @see RooAbsData::reduce
def _rad_subset_ ( dset , vars = [] , cuts = '' ) :
    """ Improved reduce
    >>> data = ...
    >>> data =  data.subset( RooArgSet( ... ) , 'a>0' )
    >>> data =  data.subset ( [a,b,c]       , 'a>0' )
    >>> data =  data.subset ( ['a','b','c'] , 'a>0' )
    - see ROOT.RooAbsData.reduce
    """

    if   isinstance ( cuts , ROOT.TCut ) :   cuts = str(cuts)
    
    if   isinstance ( vars , ROOT.RooArgSet  ) :
        return dset.reduce ( vars , cuts )
    elif isinstance ( vars , string_types    ) : vars = [ vars ]
    
    aset   = ROOT.RooArgSet()
    varset = dset.get()
    for v in vars :
        if   isinstance ( v , ROOT.RooAbsArg ) : aset.add ( v )
        elif isinstance ( v , string_types   ) and v in varset : aset.add ( varset[v] )
        else :
            raise TypeError("``subset'': unknown type %s/%s" % ( v , type(v) ) )

    return dset.reduce ( aset , cuts ) if aset else dset.reduce( cuts )


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

ROOT.RooAbsData . subset        = _rad_subset_ 

ROOT.RooAbsData . __add__       = _rad_add_
ROOT.RooDataSet . __iadd__      = _rad_iadd_

ROOT.RooAbsData . __mul__       = _rad_mul_
ROOT.RooAbsData . __rmul__      = _rad_mul_
ROOT.RooAbsData . __imul__      = _rad_imul_
ROOT.RooAbsData . __div__       = _rad_div_
ROOT.RooAbsData . __mod__       = _rad_mod_


ROOT.RooAbsData . sample        = _rad_sample_
ROOT.RooAbsData . shuffle       = _rad_shuffle_

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
   ROOT.RooDataSet . __add__       ,
   ROOT.RooDataSet . __iadd__      ,
   #
   ROOT.RooAbsData . __mul__       ,
   ROOT.RooAbsData . __rmul__      ,
   ROOT.RooAbsData . __imul__      ,
   ROOT.RooAbsData . __div__       ,
   ROOT.RooAbsData . __mod__       ,
   #
   ROOT.RooAbsData . sample        ,
   ROOT.RooAbsData . shuffle       ,
   #
   ROOT.RooAbsData . statVar       ,
   ROOT.RooAbsData . sumVar        ,
   ROOT.RooAbsData . sumVar_       ,
   #
   ROOT.RooAbsData . statCov       ,
   ROOT.RooAbsData . statCovs      ,
   ]



# =============================================================================
## Helper project method for RooDataSet/DataFrame/... and similar objects 
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
def ds_project  ( dataset , histo , what , cuts = '' , *args ) :
    """Helper project method for RooDataSet/DataFrame/... and similar objects 
    
    >>> h1   = ROOT.TH1D(... )
    >>> dataset.project ( h1.GetName() , 'm', 'chi2<10' ) ## project variable into histo
    
    >>> h1   = ROOT.TH1D(... )
    >>> dataset.project ( h1           , 'm', 'chi2<10' ) ## use histo
    
    """
    if isinstance ( cuts , ROOT.TCut ) : cuts = str ( cuts ).strip()  
    if isinstance ( what , str       ) : what = what.strip()
    if isinstance ( cuts , str       ) : cuts = cuts.strip()
    
    ## ## native RooFit...  I have some suspicion that it does not work properly
    ## if isinstance ( what  , ROOT.RooArgList ) and \
    ##    isinstance ( histo , ROOT.TH1        ) and \
    ##    hasattr ( dataset , 'fillHistogram'  ) :
    ##     histo.Reset() 
    ##     return dataset.fillHistogram  ( histo , what , cuts , *args )
    
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
        return ds_project ( dataset , histo , what.GetTitle () , cuts , *args )
    
    if   isinstance ( what , ROOT.RooAbsReal ) : 
        return ds_project ( dataset , histo , what.GetName  () , cuts , *args ) 
    
    if isinstance ( what , str ) : 
        vars  = [ v.strip() for v in what.split(':') ]
        return ds_project ( dataset , histo , vars , cuts , *args ) 
    
    if isinstance ( what , ( tuple , list ) ) :
        vars = []
        for w in what :
            if isinstance ( w , str ) : vars.append ( w.strip() )
            else                      : vars.append ( w ) 
        ##return ds_project ( dataset , histo , vars , cuts , *args ) 

    ## what is it ? 
    if isinstance ( what , ROOT.RooArgList ) :
        vars  = [ w for w in what ]
        cuts0 = cuts 
        if ''   == cuts : cuts0 = 0
        elif isinstance ( cuts , str ) :
            cuts0 = ROOT.RooFormulaVar( cuts , cuts , dataset.varlist() )
        return ds_project ( dataset , histo , vars , cuts0 , *args ) 
            
    if isinstance ( histo , str ) :
    
        obj = ROOT.gROOT     .FindObject    ( histo )
        if instance ( obj  , ROOT.TH1 ) :
            return ds_project ( dataset , obj , what , cuts , *args )
        obj = ROOT.gROOT     .FindObjectAny ( histo )
        if instance ( obj  , ROOT.TH1 ) :
            return ds_project ( dataset , obj , what , cuts , *args )
        obj = ROOT.gDirectory.FindObject    ( histo )
        if instance ( obj  , ROOT.TH1 ) :
            return ds_project ( dataset , obj , what , cuts , *args )
        obj = ROOT.gDirectory.FindObjectAny ( histo )
        if instance ( obj  , ROOT.TH1 ) :
            return ds_project ( dataset , obj , what , cuts , *args )

    ## what it is ????
    if  1 <= len ( what ) \
           and isinstance ( what[0] , ROOT.RooAbsReal ) \
           and isinstance ( cuts , str ) :
        
        if   '' == cuts : cuts0 = 0 
        elif isinstance ( cuts , str ) :
            cuts0 = ROOT.RooFormulaVar( cuts , cuts , dataset.varlist() )
        return ds_project ( dataset , histo , what , cuts0 , *args )

    
    if   isinstance ( histo , ROOT.TH3 ) and 3 == len ( what )  :
        sc = Ostap.HistoProject.project3 ( dataset ,
                                           histo   , 
                                           what[2] ,
                                           what[1] ,
                                           what[0] , cuts , *args)
        if not sc.isSuccess() :
            logger.error ( "Error from Ostap.HistoProject.project3 %s" % sc )
            return None
        return histo
    elif isinstance ( histo , ROOT.TH2 ) and 2 == len ( what )  :
        sc = Ostap.HistoProject.project2 ( dataset ,
                                           histo   , 
                                           what[1] ,
                                           what[0] , cuts , *args )
        if not sc.isSuccess() :
            logger.error ( "Error from Ostap.HistoProject.project2 %s" % sc )
            return None
        return histo
    elif isinstance ( histo , ROOT.TH1 ) and 1 == len ( what )  :
        sc = Ostap.HistoProject.project  ( dataset ,
                                           histo   , 
                                           what[0] , cuts , *args )
        if not sc.isSuccess() :
            logger.error ( "Error from Ostap.HistoProject.project  %s" % sc )
            return None
        return histo
    
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
def ds_draw ( dataset , what , cuts = '' , opts = '' , *args ) :
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
        return ds_draw ( dataset , vars , cuts , opts , *args ) 
    
    if   isinstance ( what , ROOT.RooFormulaVar ) : 
        return ds_draw ( dataset , what.GetTitle () , cuts , opts , *args )
    
    if   isinstance ( what , ROOT.RooAbsReal ) : 
        return ds_draw ( dataset , what.GetName  () , cuts , opts , *args ) 
    
    if not 1 <= len ( what ) <= 3 :
        raise AttributeError ( 'DataSet::draw, invalid length %s' % what  )
    
    if 1 == len ( what )  :
        w1        = what[0] 
        mn1 , mx1 = ds_var_minmax  ( dataset , w1 , cuts )
        histo = ROOT.TH1F ( hID() , w1 , 200 , mn1 , mx1 )  ; histo.Sumw2()
        ds_project ( dataset , histo , what , cuts , *args  )
        histo.Draw( opts )
        return histo

    if 2 == len ( what )  :
        w1        = what[0] 
        mn1 , mx1 = ds_var_minmax  ( dataset , w1 , cuts )
        w2        = what[1] 
        mn2 , mx2 = ds_var_minmax  ( dataset , w2 , cuts )
        histo = ROOT.TH2F ( hID() , "%s:%s" % ( w1 , w2 ) ,
                            50 , mn1 , mx1 ,
                            50 , mn2 , mx2 )  ; histo.Sumw2()
        ds_project ( dataset , histo , what , cuts , *args  )
        histo.Draw( opts )
        return histo

    if 3 == len ( what )  :
        w1        = what[0] 
        mn1 , mx1 = ds_var_minmax ( dataset , w1 , cuts )
        w2        = what[1] 
        mn2 , mx2 = ds_var_minmax ( dataset , w2 , cuts )
        w3        = what[2] 
        mn3 , mx3 = ds_var_minmax ( dataset , w3 , cuts )
        histo = ROOT.TH3F ( hID() , "%s:%s:%s" % ( w1 , w2 , w3 ) ,
                            20 , mn1 , mx1 ,
                            20 , mn2 , mx2 ,
                            20 , mn2 , mx2 )  ; histo.Sumw2()
        ds_project ( dataset , histo , what , cuts , *args  )
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
def ds_var_minmax ( dataset , var , cuts = '' , delta = 0.0 )  :
    """Get min/max for the certain variable in dataset
    >>> data = ...
    >>> mn,mx = data.vminmax('pt')
    >>> mn,mx = data.vminmax('pt','y>3')
    """
    if isinstance ( var , ROOT.RooAbsReal ) : var = var.GetName() 
    if cuts : s = dataset.statVar ( var , cuts )
    else    : s = dataset.statVar ( var )
    mn,mx = s.minmax()
    if mn < mn and 0.0 < delta :
        dx   = delta * 1.0 * ( mx - mn )  
        mx  += dx   
        mn  -= dx   
    return mn , mx


ROOT.RooDataSet .vminmax  = ds_var_minmax 

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

_new_methods_ += [
    ROOT.RooDataSet .clear ,
    ROOT.RooDataSet .erase ,
    ROOT.RooDataSet .Reset ,
    ]

# =============================================================================
ROOT.RooDataSet.draw         = ds_draw
ROOT.RooAbsData.draw         = ds_draw
ROOT.RooDataSet.project      = ds_project
ROOT.RooDataSet .__getattr__ = _ds_getattr_
ROOT.RooDataHist.__getattr__ = _ds_getattr_

ROOT.RooDataHist.__len__    = lambda s : s.numEntries() 

_new_methods_ += [
    ROOT.RooDataSet.draw    ,
    ROOT.RooAbsData.draw    ,
    ROOT.RooDataSet.project ,
    ]

# =============================================================================
## get the s-factor for   (weighted) dataset, where
#  s-factor is defined as
#  \f$ s_{w} \equiv \frac{\sum w_i}{\sum w_i^2} \f$
#  @see Ostap::SFactor::sFactor 
#  @code
#  dataset = ...
#  sf = dataset.sFactor() 
#  @endcode
#  when the weigths comes from sPlot, the factor effectively accounts
#  statitical fluctuations in background subtraction
#  @see W. T. Eadie et al., Statistical methods in experimental physics,
#       North Holland, Amsterdam, 1971.
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2019-05-30
# =============================================================================
def _rad_sFactor_ ( data ) :
    """Get the s-factor for   (weighted) dataset, where
    s-factor is defined as
     s_{w} equiv frac{ sum w_i}{ sum w_i^2}
     
     - see Ostap::SFactor::sFactor 
     - see W. T. Eadie et al., Statistical methods in experimental physics,
     ...   North Holland, Amsterdam, 1971.
     
     >>> dataset = ...
     >>> sf = dataset.sFactor() 
     """
    if 0 == data.numEntries() :
        logger.warning ("RooAbsData.sFactor: dataset is empty, return 1.0")
        return 1.0
    
    sf = Ostap.SFactor.sFactor ( data )
    if    0 >  sf.cov2() :
        logger.error   ('Ostap::SFactor::sFactor %s, return 1.0' % sf )
        return 1.0 
    elif  0 == sf.cov2() :
        logger.warning ('Ostap::SFactor::sFactor %s, return 1.0' % sf )
        return 1.0 
    
    return  sf.value() / sf.cov2()


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

ROOT.RooDataSet.draw        = ds_draw
ROOT.RooDataSet.project     = ds_project
ROOT.RooDataSet.__getattr__ = _ds_getattr_
ROOT.RooAbsData.sFactor     = _rad_sFactor_


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
    ROOT.RooDataHist.__getattr__  ,
    ROOT.RooDataHist.__len__      ,
    ROOT.RooAbsData .sFactor      
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
## Add/calculate/sample variable to RooDataSet
#  - Use formula expression 
#  @code
#  dataset.add_new_var ( 'ratio' , 'pt/pz' )  ## use RooFormulaVar
#  @endcode 
#  - Use  function:
#  @code
#  func = ...  ## Ostap.IFuncData object 
#  dataset.add_new_var ( 'value' ,   func )
#  @endcode
#  - Sample from 1D-historgam
#  @code
#  h1 = ...## 1D histogram
#  dataset.add_new_var ( 'nTracks' , h1 ) ## sample from 1D histogram
#  @endcode 
#  - Sample from 2D histogram
#  @code
#  h2 = ...## 2D histogram
#  dataset.add_new_var ( 'Pt' , 'eta' , h2 ) ## sample from 2D histogram
#  @encode
#  - Sample  from 3D-histogram
#  @code
#  h3 = ...## 3D histogram
#  dataset.add_new_var ( 'Pt' , 'eta' , 'A' , h3 ) ## sample from 3D histogram
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
def add_new_var ( dataset , *args ) : 
    """Add/calculate/sample variable to RooDataSet

    >>> dataset.add_new_var ( 'ratio' , 'pt/pz' )  ## use RooFormulaVar
    
    >>> func = ...  ## Ostap.IFuncData object 
    >>> dataste.add_new_var ( 'value' ,   func )

    >>> h1 = ...## 1D histogram
    >>> dataset.add_new_var ( 'nTracks' , h1 ) ## sample from 1D histogram
    
    >>> h2 = ...## 2D histogram
    >>> dataset.add_new_var ( 'Pt' , 'eta' , h2 ) ## sample from 2D histogram
    
    >>> h3 = ...## 3D histogram
    >>> dataset.add_new_var ( 'Pt' , 'eta' , 'A' , h3 ) ## sample from 3D histogram
    
    """
    
    vv = Ostap.Functions.add_var ( dataset , *args )
    if not  vv : logger.error('add_new_var: NULLPTR from Ostap.Functions.add_var')
    #
    return dataset 

# =============================================================================
ROOT.RooDataSet.addVar      = _rds_addVar_
ROOT.RooDataSet.add_new_var = add_new_var 
ROOT.RooDataSet.add_var     = add_new_var 

_new_methods_ += [
    ROOT.RooDataSet .addVar      ,
    ROOT.RooDataSet .add_new_var ,
    ROOT.RooDataSet .add_var     ,
    ]

# =============================================================================
## Add specific re-weighting information into dataset
#  @see ostap.tools.reweight
#  @see ostap.tools.reweight.Weight 
#  @see ostap.tools.reweight.W2Data 
#  @code
#  w    = Weight ( ... ) ## weighting object ostap.tools.reweight.Weight 
#  data = ...
#  data.add_reweighting ( w ) 
#  @endcode 
def add_reweighting ( data , weighter , name = 'weight' ) :
    """Add specific re-weighting information into dataset
    
    >>> w    = Weight ( ... ) ## weighting object ostap.tools.reweight.Weight 
    >>> data = ...
    >>> data.add_reweighting ( w )
    - see ostap.tools.reweight
    - see ostap.tools.reweight.Weight 
    - see ostap.tools.reweight.W2Data 
    """
    
    import ostap.tools.reweight as W
    
    assert isinstance ( weighter , W.Weight ), "Invalid type of ``weighting''!"
    
    ## create the weigthting function 
    wfun = W.W2Data ( weighter  )
    
    return data.add_new_var ( name , wfun ) 


ROOT.RooDataSet.add_reweighting = add_reweighting
_new_methods_ += [
    ROOT.RooDataSet .add_reweighting
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

ROOT.RooDataSet.makeWeighted = _rds_makeWeighted_

# =============================================================================
## ``Unweight'' weighted  dataset
#  @code
#  wdata = ...
#  data  = wdata.unweight() 
#  @endcode
def _rds_unWeighted_ ( dataset , weight = '' ) :
    """``Unweight'' weighted  dataset
    >>> wdata = ...
    >>> data  = wdata.unweight() 
    """
    if not dataset.isWeighted() :
        logger.error ("unweight: dataset is not weighted!") 
        return None, ''
    
    ds , w = Ostap.Utils.unweight (  dataset , weight )  
    return ds , w 

# =============================================================================
ROOT.RooDataSet.unWeighted = _rds_unWeighted_

_new_methods_ += [
    ROOT.RooDataSet .makeWeighted ,
    ROOT.RooDataSet .unWeighted   ,
    ]
# =============================================================================


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
## get nth moment of the distribution
def _rad_moment_ ( data , var , order , value = 0 , error = True , *args ) :
    """ Get n-th moment of the distribution
    >>> data = ...
    >>> print data.moment ( 'mass' , 3 ) 
    """
    assert isinstance ( order , int ) and 0 <= order, 'Invalid "order"  %s' % order
    
    if isintance  ( var  , str ) :
        varset =  data.get()
        assert  var in varset, 'Invalid variable %s' % var 
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

_new_methods_ += [
    RAD.moment           , 
    RAD.central_moment   , 
    RAD.skewness         , 
    RAD.kurtosis         , 
    ]


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

_new_methods_ += [
    RAD.branches , 
    RAD.leaves
    ]


# ==============================================================================
def _ds_table_0_ ( dataset        ,
                   variables = [] ,
                   cuts      = '' ,
                   first     = 0  ,
                   last      = 2**62 ,
                   prefix    =  '' ) :
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

    tt = dataset.GetTitle()
    if tt and tt != dataset.GetName() : 
        title = '%s("%s","%s"):' % ( dataset.__class__.__name__ , dataset.GetName () , tt ) 
    else :
        title = '%s("%s"):'      % ( dataset.__class__.__name__ , dataset.GetName () ) 
        
    title +=  '%d entries, %d variables' %  ( len ( dataset ) , len ( varset ) )

    if not _vars :
        return report , 120 


    weight = None
    if   isinstance ( dataset , ROOT.RooDataHist ) :
        if dataset.isNonPoissonWeighted() : title += ' Binned/Weighted' 
        else                              : title += ' Binned'
    elif dataset.isWeighted () :
        
        if dataset.isNonPoissonWeighted() : title += ' Weighted' 
        else :                              title += ' Weighted/Poisson'

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

        if store and hasattr ( store , 'tree' ) and valid_pointer ( store.tree() ) :

            tree = store.tree() 
            branches = set ( tree.branches() )
            vvars    = set ( [ i.GetName() for i in  varset ] )
            wvars    = branches - vvars
            
            if 1 == len ( wvars ):
                wvar = wvars.pop()
                
        if not wvar : wvar   = Ostap.Utils.getWeight ( dataset )
        if     wvar : title += ' with "%s"' % wvar
                
        store = None 
        if dstmp :            
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
                if wvar in tree.branches () : 
                    s = tree.statVar ( wvar , cuts , first , last ) ## no cuts here... 
                    mnmx = s.minmax ()
                    mean = s.mean   ()
                    rms  = s.rms    ()
                    weight = '*%s*' % wvar
                    r    = (  weight                            ,   ## 0 
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
        
    index_l =   int ( math.ceil ( math.log10( len ( _vars ) + 1 ) ) )
    
    fmt_name = '%%%ds. %%-%ds' % ( index_l , name_l )
    fmt_desc = '%%-%ds' % desc_l
    fmt_mean = '%%%ds'  % mean_l
    fmt_rms  = '%%-%ds' % rms_l
    fmt_min  = '%%%ds'  % min_l
    fmt_max  = '%%-%ds' % max_l

    title_l = index_l + 2 + name_l  
    header = [ ( '{:^%d}' % title_l ).format ( 'Variable'    ) ,
               ( '{:^%d}' % desc_l  ).format ( 'Description' ) ,
               ( '{:^%d}' % mean_l  ).format ( 'mean'        ) ,
               ( '{:^%d}' % rms_l   ).format ( 'rms'         ) ,
               ( '{:^%d}' % min_l   ).format ( 'min'         ) ,
               ( '{:^%d}' % max_l   ).format ( 'max'         ) ]

    if weight : header.append ( 'W' )
        
    table_data = [ tuple  ( header ) ]

    vlst = vars

    for i , v in enumerate ( _vars ) :
                
        cols = [ fmt_name %  ( i + 1 , v [ 0 ] ) ,
                 fmt_desc %            v [ 1 ] ,
                 fmt_mean %            v [ 2 ] ,
                 fmt_rms  %            v [ 3 ] ,
                 fmt_min  %            v [ 4 ] ,
                 fmt_max  %            v [ 5 ] ]
        
        if   weight and i + 1 == len ( _vars ) :
            cols.append ( 'W' )
            cols = [ allright (  c ) for c in cols ] 
        elif weight                            : cols.append ( ' ' )
        
        table_data.append ( tuple ( cols ) ) 

    import ostap.logger.table as T
    t  = T.table ( table_data , title , prefix =  prefix )
    w  = T.table_width ( t ) 
    return t , w 

# ==============================================================================
## print dataset in a form of the table
#  @code
#  dataset = ...
#  print dataset.table() 
#  @endcode
def _ds_table_ (  dataset ,  variables = [] , prefix = '' ) :
    """print dataset in a form of the table
    >>> dataset = ...
    >>> print dataset.table()
    """
    return _ds_table_0_ ( dataset ,  variables , prefix = prefix )[0]

# =============================================================================
##  print DataSet
def _ds_print2_ ( dataset ) :
    """Print dataset"""
    if dataset.isWeighted() and not isinstance ( dataset , ROOT.RooDataHist ) :
        store = dataset.store()
        if valid_pointer ( store ) and isinstance ( store , ROOT.RooTreeDataStore ) : pass
        else : return _ds_print_ ( dataset )         
    from ostap.utils.basic import terminal_size, isatty 
    if not isatty() : return _ds_table_ ( dataset )
    th  , tw  = terminal_size()
    rep , wid = _ds_table_0_ ( dataset ) 
    if wid < tw  : return rep
    return _ds_print_ ( dataset )


for t in ( ROOT.RooDataSet , ROOT.RooDataHist ) :
    t.__repr__    = _ds_print2_
    t.__str__     = _ds_print2_
    t.table       = _ds_table_
    t.pprint      = _ds_print_ 

    
_new_methods_ += [
    ROOT.RooDataSet.table    , 
    ROOT.RooDataSet.pprint   , 
    ROOT.RooDataSet.__repr__ ,
    ROOT.RooDataSet.__str__  ,
    ]

# =============================================================================
## make symmetrization/randomization of the dataset
#  @code
#  ds     = ...
#  ds_sym = ds.symmetrize ( 'var1' , 'var2' )
#  ds_sym = ds.symmetrize ( 'var1' , 'var2' , 'var3')
#  @endcode
def _ds_symmetrize_ ( ds , var1 , var2 , *vars ) :
    """Make symmetrization/randomization of the dataset
    >>> ds     = ...
    >>> ds_sym = ds.symmetrize ( 'var1' , 'var2' )
    >>> ds_sym = ds.symmetrize ( 'var1' , 'var2' , 'var3')
    """
    
    varset = ds.varset() 
    lvars  = [ var1 , var2 ] + list ( vars )
    nvars  = [] 
    for v in lvars :
        if not v in varset : raise NameError ( "Variable %s not in dataset" % v )
        if not isinstance ( v , ROOT.RooAbsReal ) : v = varset[ v ]
        nvars.append ( v )

    mnv = min ( [ v.getMin () for v in nvars if hasattr ( v , 'getMin' ) ] ) 
    mxv = max ( [ v.getMax () for v in nvars if hasattr ( v , 'getMax' ) ] ) 

    names   = [ v.name for v in nvars ]

    nds     = ds.emptyClone ()
    nvarset = nds.varset    ()
    
    for v in nvarset :
        if v.name in names : 
            if hasattr ( v ,  'setMin' ) : v.setMin ( mnv )
            if hasattr ( v ,  'setMax' ) : v.setMax ( mxv )        
    
    ## loop over the data set 
    for entry in ds :

        values = [ v.getVal() for v in entry if v in varset ]
        random.shuffle ( values )

        for v in nvarset :
            n = v.name 
            if not n in names : v.setVal ( entry[n].value )                
            else              : v.setVal ( values.pop()   )

        nds.add ( nvarset )
        
    return nds


ROOT.RooDataSet.symmetrize = _ds_symmetrize_

_new_methods_ += [
    ROOT.RooDataSet.symmetrize , 
    ]


# =============================================================================

from  ostap.stats.statvars import data_decorate as _dd
_dd ( ROOT.RooAbsData )

_decorated_classes_ = (
    ROOT.RooAbsData ,
    ROOT.RooDataSet ,
    )

_new_methods_ = tuple ( _new_methods_ ) 
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
