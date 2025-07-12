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
""" Module with decoration for RooAbsData and related RooFit classes
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
    'ds_draw'    , ## draw variables from RooDataSet 
    'ds_project' , ## project variables from RooDataSet to histogram
    'ds_combine' , ## combine two datasets with weights 
    )
# =============================================================================
from   collections               import defaultdict
from   ostap.core.meta_info      import root_info
from   ostap.core.core           import ( Ostap         ,
                                          VE , SE , hID , dsID ,
                                          strings       , 
                                          valid_pointer ,
                                          ROOTCWD       )
from   ostap.core.ostap_types    import ( integer_types , string_types   ,
                                          num_types     , dictlike_types , 
                                          list_types    , sequence_types )
from   ostap.utils.strings       import ( split_string         , 
                                          split_string_respect ,
                                          var_separators       )
from   ostap.utils.basic         import loop_items  , typename            
from   ostap.math.base           import islong, evt_range, FIRST_ENTRY, LAST_ENTRY 
from   ostap.utils.random_seed   import random_seed
from   ostap.fitting.variables   import valid_formula, make_formula 
from   ostap.trees.cuts          import expression_types, vars_and_cuts, order_warning
from   ostap.stats.statvars      import data_decorate, data_range 
from   ostap.utils.valerrors     import VAE
from   ostap.logger.symbols      import cabinet, weight_lifter 
from   ostap.utils.progress_conf import progress_conf
from   ostap.utils.progress_bar  import progress_bar 
import ostap.fitting.roocollections
import ostap.fitting.printable
import ostap.io.root_file 
import ROOT, random, math, os, sys, ctypes  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.dataset' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Some useful decorations for RooAbsData object')
# =============================================================================
from ostap.logger.colorized import allright,  attention
_new_methods_ = []
# =============================================================================
_maxv =  0.99 * sys.float_info.max
_minv = -0.99 * sys.float_info.max
# =============================================================================
## iterator for RooAbsData entries 
#  @cdoe
#  dataset = ...
#  for entry, weight  in dataset : ...
#  @endcode 
#  - For unweighted datatsets, `weight` is `None`
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _rad_iter_ ( self ) :
    """ Iterator for RooAbsData
    >>> dataset = ...
    >>> for entry,weight  in dataset : ... 
    - for unweighted datatsets, `weight` is `None`
    """
    _l = len ( self )
    for i in range ( 0 , _l ) :
        yield self [ i ]

# ===========================================================================
## Iterator over "good" events in dataset 
#  @code
#  dataset = ...
#  for index , entry, weight in dataset.loop ( 'pt>1' ) :
#     print (index, entry, weight) 
#  @endcode 
def _rad_loop_ ( dataset                 ,
                 cuts      = ''          ,
                 cut_range = ''          ,
                 first     = FIRST_ENTRY ,
                 last      = LAST_ENTRY  ,
                 progress = False        ) :
    """ Iterator for `good' events  in dataset
    >>> dataset = ...
    >>> for index, entry, weight in dataset.loop ( ''pt>1' ) :
    >>>    print (index, entry, weight) 
    """

    first, last = evt_range ( dataset , first , last ) 
    if last <= first : return

    assert isinstance ( cuts    , expression_types  ) or not cuts, \
        "Invalid type of cuts: %s" % type ( cuts )
    assert isinstance ( cut_range , expression_types ) or not cut_range, \
        "Invalid type of cut_range: %s" % type ( cut_range )
    
    cuts      = str(cuts).strip()      if cuts      else ''
    cut_range = str(cut_range).strip() if cut_range else ''

    fcuts = None 
    if cuts : fcuts = make_formula ( cuts , cuts , dataset.varlist() )

    weighted          = dataset.isWeighted                     ()
    store_errors      = weighted and dataset.store_errors      ()
    store_asym_errors = weighted and dataset.store_asym_errors () 
    simple_weight     = weighted and ( not store_errors ) and ( not store_asym_errors )
    
    ## loop over dataset
    source = range ( first , last )
    if progress : source = progress_bar ( source , description = 'Entries:' )
    
    nevents = 0 
    for event in source : 

        entry, weight = dataset [ event ]        
        if not entry: break

        if cut_range and not entry.allInRange ( cut_range ) : continue

        wc = fcuts.getVal() if fcuts else 1.0 
        if not wc  : continue
        
        if weight is None : weight = wc if cuts else weight 
        else              : weight = weight * wc

        nevents += 1 
        yield event , entry , weight 
        
    del fcuts
    ## report summary 
    if progress : logger.info ( 'loop: %d from %d entries' % ( nevents , last - first ) ) 


ROOT.RooAbsData.loop = _rad_loop_         


# =============================================================================
## access to the entries in  RooAbsData
#  @code
#  dataset        = ...
#  event , weight = dataset [4]           ## index 
#  events         = dataset[0:1000]       ## slice 
#  events         = dataset[0:-1:10]      ## slice 
#  events         = dataset[ (1,2,3,10) ] ## sequence of indices  
#  @endcode
#  - For unweighted ddatasets `weight` is `None`
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-03-31
def _rad_getitem_ ( data , index ) :
    """ Get the entry from RooDataSet

    >>> dataset  = ... 
    >>> event, weight  = dataset [4]                 ## index 

    - For unweighted ddatsets `weight` is `None`    

    >>> events = dataset[0:1000]            ## slice
    >>> events = dataset[0:-1:10]           ## slice 
    >>> events = dataset[ (1,2,3,4,10) ]    ## sequnce of indices
 
    """
    
    N = len ( data )
    
    if isinstance ( index , integer_types ) and index < 0 : index += N

    ## simple entry index 
    if isinstance ( index , integer_types ) and 0 <= index  < N :
        
        entry = data.get ( index )
        if not data.isWeighted() : return entry, None       ## RETURN 

        weight = data.weight()        
        if   data.store_asym_error () :
            wel , weh = data.weight_errors ()
            weight    = VAE ( weight , wel , weh ) 
        elif data.store_error      () :
            we = data.weightError  ()
            if 0 <= we : weight = VE ( weight , we * we ) 
            
        return entry, weight                               ## RETUR N

    ## range -> range 
    elif isinstance ( index , range ) :

        ## simple case 
        start , stop , step = index.start , index.stop , index.step
        if 1 == step and start <= stop : return data.reduce ( ROOT.RooFit.EventRange ( start , stop ) )  ## RETURN 
        index = range ( start , stop , step ) 

    ## slice -> range 
    elif isinstance ( index , slice ) :
        
        start , stop , step = index.indices ( N )                              
        if 1 == step and start <= stop : return data.reduce ( ROOT.RooFit.EventRange ( start , stop ) ) ## RETURN 
        index = range ( start , stop , step ) 

    ## require sequence of indices here 
    if not isinstance ( index , sequence_types ) :
        raise IndexError ( "Invalid type of `index':%s" % type ( index ) )

    weighted          = data.isWeighted                     ()
    store_errors      = weighted and data.store_errors      ()
    store_asym_errors = weighted and data.store_asym_errors ()

    ## preare the result 
    result = data.emptyClone ( dsID () )
    
    ## the actual loop over set of entries 
    for i , j in enumerate ( index ) :

        ## if not isinstance ( j , integer_types ) :
        ##    raise IndexError ( 'Invalid index [%s]=%s/%s' % ( i , j , type ( j ) ) )
        
        j = int ( j )            ## the content must be convertible to integer 
        if j < 0 : j += N        ## allow `slightly-negative' indices             
        if not 0 <= j < N :      ## adjusted integer in the proper range ?
            raise IndexError ( 'Invalid index [%s]=%s,' % ( i , j ) ) 

        vars = data.get ( j )
        if not vars : raise IndexError ( 'Invalid index %s' % j )  ## 
        
        if   store_asym_errors :  
            wel , weh = data.weight_errors ()
            result.add ( vars , data.weight () , wel , weh ) 
        elif store_errors      :
            we        = data.weightError   ()
            result.add ( vars , data.weight () , we  ) 
        elif weighted         :
            result.add ( vars , data.weight () ) 
        else : 
            result.add ( vars ) 
            
    return result

# ==============================================================================


# ==============================================================================
## Get (asymmetric) weigth errors for the current entry in dataset
#  @code
#  dataset = ...
#  weight_error_low, weight_error_high = dataset.weight_errors () 
#  @endcode
#  @see RooAbsData::weightError
def _rad_weight_errors ( data , *etype ) :
    """ Get (asymmetric) weigth errors for the current entry in dataset
    >>> dataset = ...
    >>> weight_error_low, weigth_error_high = dataset.weight_errors () 
    - see ROOT.RooAbsData.weightError
    """
    ##
    if not w.isWeighted () : return 0.0, 0.0
    ##
    wel = ctypes.c_double ( 0.0 )
    weh = ctypes.c_double ( 0.0 )
    data.weightError ( wel , weh )
    #
    return float ( wel.value ) , float ( weh.value )

# =============================================================================
## Get variables in form of RooArgList 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-03-31
def _rad_vlist_ ( self ) :
    """ Get variables in form of RooArgList 
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
    """ Check the presence of variable in dataset    
    >>> if 'mass' in dataset : print 'ok!'
    """
    if isinstance ( aname , integer_types ) :
        return 0 <= aname < len ( self )
    
    vset = self.get()
    return aname in vset 

# =============================================================================
## merge/append two datasets into a single one
# @code
# dset1  = ...
# dset2  = ...
# dset1 += dset2
# @endcode 
def _rad_iadd_ ( ds1 , ds2 ) :
    """ Merge/append two datasets into a single one
    - two datasets must have identical structure 
    >>> dset1  = ...
    >>> dset2  = ...
    >>> dset1 += dset2
    """
    if isinstance ( ds1 , ROOT.RooDataSet ) and isinstance ( ds2 , ROOT.RooDataSet ) :
        ##
        w1 = ds1.isWeighted()
        w2 = ds2.isWeighted()
        ##
        if   w1 and w2 : pass
        elif w1        : return NotImplemented 
        elif w2        : return NotImplemented 
        ## 
        npw1 = ds1.IsNonPoissonWeighted()
        npw2 = ds2.IsNonPoissonWeighted()
        ##
        if   npw1 and npw2 : pass
        elif npw1          : return NotImplemented 
        elif npw2          : return NotImplemented 
        ## 
        vs1 = set ( v.name for v in ds1.get() )
        vs2 = set ( v.name for v in ds2.get() )
        ## 
        if vs1 != vs2  : return NotImplemented
        ## 
        ds1.append ( ds2 )
        return ds1 
    
    return NotImplemented

# =============================================================================
## merge/append two datasets into a single one
#  @code
#  dset1  = ...
#  dset2  = ...
#  dset   = dset1 + dset2 
#  @endcode 
def _rad_add_ ( ds1 , ds2 ) :
    """ Merge/append two datasets into a single one
    - two datasets must have identical structure 
    >>> dset1  = ...
    >>> dset2  = ...
    >>> dset   = dset1 + dset2 
    """
    if isinstance ( ds1 , ROOT.RooDataSet ) and isinstance ( ds2 , ROOT.RooDataSet ) :
        ##
        w1 = ds1.isWeighted()
        w2 = ds2.isWeighted()
        ##
        if   w1 and w2    : pass
        elif w1           : return NotImplemented 
        elif w2           : return NotImplemented 
        ##
        npw1 = ds1.isNonPoissonWeighted()
        npw2 = ds2.isNonPoissonWeighted()
        ##
        if   npw1 and npw2 : pass
        elif npw1          : return NotImplemented 
        elif npw2          : return NotImplemented 
        ## 
        vs1 = set ( v.name for v in ds1.get() )
        vs2 = set ( v.name for v in ds2.get() )
        ##
        if vs1 != vs2  : return NotImplemented
        ##        
        result = ds1.emptyClone( dsID() ) 
        result.append ( ds1 )
        result.append ( ds2 )
        return result
    
    return NotImplemented


# =============================================================================
# merge/append two datasets into a single one 
def _rad_imul_ ( ds1 , ds2 ) :
    """ Merge/append two datasets into a single one
    - two datasets must have the  same number of entries!
    >>> dset1  = ...
    >>> dset2  = ...
    >>> dset1 *= dset2
    """
    if isinstance ( ds1 , ROOT.RooDataSet ) and isinstance ( ds2 , ROOT.RooDataSet ) :
        if len ( ds1 ) != len ( ds2 ) : return NotImplemented
        ##
        w1 = ds1.isWeighted()
        w2 = ds2.isWeighted()
        ##
        if   w1 and w2    : pass
        elif w1           : return NotImplemented 
        elif w2           : return NotImplemented 
        ##
        npw1 = ds1.isNonPoissonWeighted()
        npw2 = ds2.isNonPoissonWeighted()
        ##
        if   npw1 and npw2 : pass
        elif npw1          : return NotImplemented 
        elif npw2          : return NotImplemented 
        ## 
        ds1.merge ( ds2 )
        return ds1
    
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
def _rad_mul_ ( ds1 , ds2 ) :
    """
    - (1) Get small (random) fraction of  dataset:
    >>> dataset = ....
    >>> small   = 0.1 * dataset
    - (2) Merge two dataset (of the same length)
    >>> dataset3 = dataset1 * dataset2 
    """

    if isinstance ( ds1 , ROOT.RooDataSet ) and isinstance ( ds2 , ROOT.RooDataSet ) :
        if len ( ds1 ) != len ( ds2 ) : return NotImplemented 
        ## 
        w1 = ds1.isWeighted()
        w2 = ds2.isWeighted()
        ##
        if   w1 and w2    : pass
        elif w1           : return NotImplemented 
        elif w2           : return NotImplemented 
        ##
        npw1 = ds1.isNonPoissonWeighted()
        npw2 = ds2.isNonPoissonWeighted()
        ##
        if   npw1 and npw2 : pass
        elif npw1          : return NotImplemented 
        elif npw2          : return NotImplemented 
        ##         
        result = ds1.emptyClone( dsID() )
        result.append ( ds1 )
        result.merge  ( ds2 )
        return ds1 

    # =======================================================================
    fraction = ds2 
    if  isinstance ( fraction , float ) and 0 < fraction < 1 :

        weighted          = ds1.isWeighted ()
        store_errors      = weighted and ds1.store_error      () 
        store_asym_errors = weighted and ds1.store_asym_error ()
        
        res  = ds1.emptyClone()
        l    = len ( ds1 )
        for i in range ( l ) :            
            if random.uniform ( 0 , 1 ) <= fraction :
                if   store_asym_errors :
                    wel , weh = ds1.weight_errors () 
                    entry = ds1.get ( i ) , ds1.weight () , wel , weh 
                elif store_error :
                    entry = ds1.get ( i ) , ds1.weight () , ds1.weightError ()
                elif weighted :
                    entry = ds1.get ( i ) , ds1.weight ()
                else: 
                    entry = ds1.get ( i ) ,
                ##    
                res.append ( *entry ) 
                
        return res
    
    elif 1 == fraction : return ds1.clone      ()
    elif 0 == fraction : return ds1.emptyClone () 

    return NotImplemented

# =============================================================================
## merge two dataset (of same  length) OR get small (random) fraction of  dataset
#  @code
#  ## get smaller dataset:
#  dataset = ....
#  small   = 0.1 * dataset 
#  ## merge two dataset of the same lenth
#  merged  = dataset1 * dataset2 
#  @endcode
def _rad_rmul_ ( ds1 , ds2 ) :
    """
    - (1) Get small (random) fraction of  dataset:
    >>> dataset = ....
    >>> small   = 0.1 * dataset
    - (2) Merge two dataset (of the same length)
    >>> dataset3 = dataset1 * dataset2 
    """
    return _rad_mul_ ( ds1 , ds2 )

    
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
    if  isinstance ( fraction , num_types ) :
        if    1.0 <  fraction : return _rad_mul_  ( self , 1.0 / fraction )
        elif  1   == fraction : return self.clone ()
    
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
            res.add ( self [ i ] ) 
        return res 
        
    elif 1 == fraction : return self.clone      ()

    return NotImplemented

# =============================================================================
## Make dataset with removed i-th element  (for Jackknife/bootstrapping)
#  @code
#  dataset = ...
#  N = len ( dataset)
#  for index in range ( N ) :
#    ds_i = dataset - index
#  @endcode
def _rds_remevt_ ( dataset , index ) :
    """ Make dataset with removed i-th element  (for Jackknife/bootstrapping)
    >>> dataset = ...
    >>> N = len ( dataset)
    >>> for index in range ( N ) :
    >>> ... ds_i = dataset - index
    """
    N = len ( dataset )
    
    if 1 < N and isinstance ( index , integer_types ) :

        if index < 0 : index += N  ## allow negative indices 

        if not 0 <= index < N : return NotImplemented
        
        if   0 == index     : return dataset.reduce ( ROOT.RooFit.EventRange ( 1 , N     ) )
        elif N == index + 1 : return dataset.reduce ( ROOT.RooFit.EventRange ( 0 , N - 1 ) )

        ds1 = dataset.reduce ( ROOT.RooFit.EventRange ( 0         , index ) )
        ds2 = dataset.reduce ( ROOT.RooFit.EventRange ( index + 1 , N     ) )

        result = ds1 + ds2

        assert len ( result ) + 1 == N , 'Invalid length of the resulting dataset!'

        ds1 = Ostap.MoreRooFit.delete_data ( ds1 )
        del ds1
        
        ds2 = Ostap.MoreRooFit.delete_data ( ds2 )
        del ds2
        
        return result
    
    return NotImplemented

# ============================================================================
## Remove certain columns from dataset
#  @code
#  dataset  = ...
#  another1 = dataset - 'pt'
#  another2 = dataset - ( 'pt', 'mass' )
#  another3 = dataset - ( 'pt,mass,z' )
#  @endcode
def _rds_remvar_ ( dataset , variables ) :
    """ Remove certain columns from dataset
    >>> dataset  = ...
    >>> another1 = dataset - 'pt'
    >>> another2 = dataset - ( 'pt', 'mass' )
    >>> another3 = dataset - ( 'pt,mass,z' )
    """
    ## decode list of variables
    try :
        vars , _ , _ = vars_and_cuts ( variables , '' )
        if not vars : return NotImplemented             ## RETURN 
    except AssertionError :
        return NotImplemented                           ## RETURN 
    ##
    if not all ( v in dataset for v in vars ) : return NotImplemented   
    ## 
    varset  = dataset.get() 
    newvars = ROOT.RooArgSet()
    for v in varset :
        if not v.name in vars : newvars.add ( v )
    ## 
    return dataset.reduce (  ROOT.RooFit.SelectVars ( newvars ) )  

# =============================================================================
## remove entry or variable(s) from dataset
#  @code
#  dataset  = ...
#  another- = dataset - 100   ## remove entry #100 
#  another1 = dataset - 'pt'
#  another2 = dataset - ( 'pt', 'mass' )
#  another3 = dataset - ( 'pt,mass,z' )
#  @endcode
def _rds_sub_ ( dataset , what ) : 
    """Remove entry or variable(s) from dataset
    >>> dataset  = ...
    >>> another- = dataset - 100   ## remove entry #100 
    >>> another1 = dataset - 'pt'
    >>> another2 = dataset - ( 'pt', 'mass' )
    >>> another3 = dataset - ( 'pt,mass,z' )
    """
    if isinstance ( what , integer_types ) :
        return _rds_remevt_ ( dataset , what )
    return _rds_remvar_ ( dataset , what )
    
# ============================================================================
## Jackknife generator: generates data sets with removed i-th element
#  @code
#  dataset = ...
#  for ds in ds.jackknife() :
#      ...
#  @endcode
#  - Dataset need to be deleted explicitely:
#  @code
#  dataset = ...
#  for ds in ds.jackknife() :
#      ...
#      ds = Ostap.MoreRooFit.delete_data ( ds )
#      del ds 
#  @endcode
#  - Alternatively one can use `delete=True` :
#  dataset = ...
#  for ds in ds.jackknife( delete = True ) :
#      ...
#  @endcode
#  @see Ostap::MoreRooFit::delete_data
def _rds_jackknife_ ( dataset ,
                      first    = FIRST_ENTRY  ,
                      last     = LAST_ENTRY   ,
                      delete   = False        ,
                      progress = False        ) :
    """ Jackknife generator

    >>> dataset = ...
    >>> for ds in ds.jackknife() :
    >>> ...
    
    - Dataset need to be deleted explicitely:
    >>> dataset = ...
    >>> for ds in ds.jackknife() :
    >>>     ...
    >>>     ds = Ostap.MoreRooFit.delete_data ( ds )
    >>>     del ds 
    
    - Alternatively one can use `delete=True` :

    >>> dataset = ...
    >>> for ds in ds.jackknife( delete = True ) :
    >>>     ...

    - see `Ostap.MoreRooFit.delete_data`

    """
    first , last = evt_range ( dataset , first , last )    
    for i in progress_bar ( range ( first , last ) , silent = not progress ) :
        ds = dataset - i                    ## this is the line! 
        yield ds
        if delete :
            ds = Ostap.MoreRooFit.delete_data ( ds )
            del ds
            
# =============================================================================
## Boostrap generator
#  @code
#  dataset = ...
#  for ds in dataset.bootstrap ( 100 ) :
#  ...
#  @endcode
#  The dataset must be remove explicitely:
#  @code
#  dataset = ...
#  for ds in dataset.bootstrap ( 100 ) :
#      ...
#      ds = Ostap.MoreRooFit.delete_data ( ds )
#      del ds 
#  @endcode
#  Alternatively one can add `delete=True`
#  @code
#  for ds in dataset.bootstrap ( 100 , delete = True ) :
#      ...
#  @endcode 
def _rds_bootstrap_ ( dataset , size = 100 , extended = False , delete = False , progress = False ) :
    """ Boostrap generator:

    >>> dataset = ...
    >>> for ds in dataset.bootstrap ( 100 ) :
    >>>     ...

    - The dataset must be remove explicitely:

    >>> dataset = ...
    >>> for ds in dataset.bootstrap ( 100 ) :
    >>>    ...
    >>>    ds = Ostap.MoreRooFit.delete_data ( ds )
    >>>    del ds 

     - Alternatively one can add `delete=True`

    >>> dataset = ... 
    >>> for ds in dataset.bootstrap ( 100 , delete = True ) :
    >>>     ... 
    """
    from   ostap.stats.bootstrap  import bootstrap_indices, extended_bootstrap_indices 
    N    = len ( dataset )
    bgen = bootstrap_indices ( N , size = size ) if not extended else extended_bootstrap_indices ( N , size = size )    
    for indices in progress_bar ( bgen , silent = not progress , max_value = N ) :
        ds = dataset [ indices ] 
        yield ds
        if delete :
            ds = Ostap.MoreRooFit.delete_data ( ds )
            del ds
            

# =============================================================================
## get (random) unique sub-sample from the dataset
#  @code
#  data   = ...
#  subset =  data.sample ( 100  )  ## get 100   events 
#  subset =  data.sample ( 0.01 )  ## get 1% of events 
#  @endcode 
def _rad_sample_ ( self , num ) :
    """ Get (random) unique sub-sample from the dataset
    >>> data   = ...
    >>> subset =  data.sample ( 100  )  ## get 100   events 
    >>> subset =  data.sample ( 0.01 )  ## get 1% of events 
    """
    N = len ( self ) 
    if   0 == num : return self.emptyClone ( dsID () )
    elif num == N : return _rad_shuffle_ ( self )
    elif isinstance ( num , integer_types ) and 0 < num < N : pass 
    elif isinstance ( num , float ) and 0 < num < 1 :
        from ostap.math.random_ext import poisson 
        num = poisson ( num * N )
        return _rad_sample_ ( self , num )
    else :
        raise TypeError("Unknown `num':%s" % num )
    ## 
    indices = random.sample ( range ( N )  , num )
    return self [ indices ]

# =============================================================================
## get (random) sub-sample from the dataset with replacement 
#  @code
#  data   = ...
#  subset =  data.choice ( 100  )  ## get 100   events 
#  subset =  data.choice ( 0.01 )  ## get 1% of events 
#  @endcode 
def _rad_choice_ ( self , num ) :
    """ Get (random) sub-sample from the dataset with replacement 
    >>> data   = ...
    >>> subset =  data.choice ( 100  )  ## get 100   events 
    >>> subset =  data.choice ( 0.01 )  ## get 1% of events 
    """
    N = len ( self )     
    if   0 == num  or 0 == N : return self.emptyClone ( dsID () )
    elif isinstance ( num , integer_types ) and 0 < num <= N : pass 
    elif isinstance ( num , float ) and 0 < num < 1 :
        from ostap.math.random_ext import poisson 
        num = poisson ( num * N  )
        return _rad_choice_ ( self , num )
    else :
        raise TypeError("Unknown `num':%s" % num )
    ##
    indices =   random.choices   ( range ( N )  , k = num )
    ## 
    return self [ indices ]

# =============================================================================
## get the shuffled sample
#  @code
#  data = ....
#  shuffled = data.shuffle()
#  @endcode 
def _rad_shuffle_ ( self ) :
    """ Get the shuffled sample
    >>> data = ....
    >>> shuffled = data.shuffle()
    """
    
    indices = [ i for i in range ( len ( self ) ) ]  
    random.shuffle ( indices )
    return self [ indices ]

# ==============================================================================
## Imporved reduce
#  @code
#  data = ...
#  data =  data.subset ( RooArgSet( ... ) , 'a>0' )
#  data =  data.subset ( [a,b,c]       , 'a>0' )
#  data =  data.subset ( ['a','b','c'] , 'a>0' )
#  @endcode
#  @see RooAbsData::reduce
def _rad_subset_ ( dataset                  ,
                   variables  = []          ,
                   cuts       = ''          ,
                   cut_range  = ''          ,
                   first      = FIRST_ENTRY ,
                   last       = LAST_ENTRY  ) :

    """ Improved reduce
    >>> data = ...
    >>> data =  data.subset( RooArgSet( ... ) , 'a>0' )
    >>> data =  data.subset ( [a,b,c]       , 'a>0' )
    >>> data =  data.subset ( ['a','b','c'] , 'a>0' )
    - see ROOT.RooAbsData.reduce
    """
    
    assert isinstance ( cut_range , expression_types ) or not cut_range, \
        "Invalid type of cut_range: %s" % type ( cut_range )
    
    cut_range   = str ( cut_range ).strip() if cut_range else ''

    ## variables in this dataset 
    vset = dataset.get()

    var_types   = ROOT.RooAbsReal , ROOT.RooAbsCategory 
    vars        = variables 
    if   isinstance ( variables , expression_types                    ) : vars = [ str ( variables ) ]
    elif isinstance ( variables , ROOT.RooAbsCollection               ) : vars = [ str ( v.GetName() ) for v in variables ]
    elif isinstance ( variables , var_types                           ) : vars = [ str ( variables.GetName() ) ]
    elif all ( isinstance ( v , var_types        ) for v in variables ) : vars = [ str ( v.GetName() ) for v in variables ]
    elif all ( isinstance ( v , expression_types ) for v in variables ) : vars = [ str ( v           ) for v in variables ]
    elif not vars                                                       : vars = [ str ( v.GetName() ) for v in vset      ] 
    ## 
    ## decode list of variables 
    varlst, cuts, _ = vars_and_cuts ( vars , cuts )
    ## 
    ## extra items? 
    extra  = [ v for v in varlst if not v in dataset ]
    assert varlst and not extra , 'Variables are not in dataset: %s' % str ( extra ) 
    ## 
    ## create cuts as RooFormulaVar 
    if cuts : cuts = make_formula ( cuts , cuts , dataset.varlist() )
    ## 
    aset = ROOT.RooArgSet()
    for v in varlst : aset.add ( vset [ v ] )
    ##
    args = []
    ## 
    if aset      : args.append ( ROOT.RooFit.SelectVars ( aset      ) ) 
    if cuts      : args.append ( ROOT.RooFit.Cut        ( cuts      ) )
    if cut_range : args.append ( ROOT.RooFit.CutRange   ( cut_range ) )
    ## 
    ## check the range:
    first, last = evt_range ( dataset , first , last )
    if 0 < first or last < len ( dataset ) : 
        args.append ( ROOT.RooFit.EventRange ( first , last ) )
    ## 
    return dataset.reduce ( *args ) 

# =============================================================================
## helper technical method to seek for unique and/or duplicated entries
def _rds_seek_for_duplicates_ ( dataset           , 
                                entrytag          , 
                                criterium = ''    ,
                                progress  = False ) : 
    """ Helper technical method to seek for unique and/or duplicated entries
    """
    
    if   isinstance ( entrytag , ROOT.RooAbsArg ) : entrytag = [ entrytag ]
    elif isinstance ( entrytag , string_types   ) : entrytag = [ entrytag ]
    
    assert isinstance ( entrytag , sequence_types ) , 'Invalid "entrytag" %s' % str ( entrytag )

    fields = [] 
    for e in entrytag :
        if   isinstance ( e , ROOT.RooAbsArg  ) and e in dataset : fields.append ( e )
        elif isinstance ( e , string_types    ) and e in dataset :
            fields.append ( getattr ( dataset , e ) )
        elif isinstance ( e , string_types    ) and valid_formula ( e , dataset ) :
            fields.append ( make_formula ( '' , e , dataset ) )            
        elif isinstance ( e , ROOT.RooAbsReal ) : fields.append ( e )
        else :
            logger.error ( 'Unknown entry %s, skip!' % e )
            
    tag = tuple ( fields )

    if  criterium :        
        if   isinstance ( criterium , ROOT.RooAbsArg ) and criterium in dataset : crit_var = criterum
        elif isinstance ( criterium , string_types   ) and criterium in dataset :
            crit_var = getattr ( dataset , criterium )
        elif isinstance ( criterium  , string_types  ) and valid_formula ( criterium , dataset ) :
            crit_var = make_formula ( '' , criterium , dataset ) 
        elif isinstance ( e , ROOT.RooAbsReal ) :
            crit_var = criterium
        else :
            raise TypeError ( 'Invalid criterium type!' ) 
    else :
        crit_var = None
    
    snapshot = defaultdict(list)

    if crit_var : content = lambda i : ( float ( crit_var ) , i )
    else        : content = lambda i : ( 0                  , i )
    
    ## make a loop over dataset
    for i, e in progress_bar ( enumerate ( dataset ) ,
                               max_value   = len ( dataset ) ,
                               description = '!st loop:'     , 
                               silent      = not progress    ) : 

        entry = tuple (  float ( v ) for v in tag ) 
        snapshot [ entry ].append ( content ( i )  ) 

    return snapshot

# =============================================================================
## Iterator over the duplicated groups 
#  @code
#  for group in dataset.duplicates ( ( 'evt' , 'run' ) ) :
#  ... for entry in group :
#  ... ...
def _rds_duplicates_ ( dataset          ,
                       entrytag         ,
                       progress = False ) :
    """ Iterator over the duplicated groups 
    >>> for group in dataset.duplicates ( ( 'evt' , 'run' ) ) :
    >>> ... for entry in group :
    >>> ... ...
    """
    snapshot = _rds_seek_for_duplicates_ ( dataset              ,
                                           entrytag  = entrytag ,
                                           criterium = ''       ,
                                           progress  = progress )
    
    for e, lst in progress_bar ( loop_items ( snapshot ) ,
                                 max_value   = len ( snapshot ) ,
                                 description = '2nd loop:'      ,
                                 silent      = not progress     ) : 
        if 2 <= len ( lst ) :
            yield tuple ( sorted ( l [ 1 ] for l in lst ) ) 
            
# =============================================================================        
## Iterator over the unique entries in dataset
#  @code
#  dataset = ...
#  for ientry in dataset.unique_entries ( ( 'evt' , 'run' ) , choice = 'random' ) :
#  ...
#  for ientry in dataset.unique_entries ( ( 'evt' , 'run' ) , choice = 'first'  ) :
#  ...
#  for ientry in dataset.unique_entries ( ( 'evt' , 'run' ) , choice = 'last'    ) :
#  ...
#  for ientry in dataset.unique_entries ( ( 'evt' , 'run' ) , choice = 'max'  , criterium = 'PT') :
#  ...
#  for ientry in dataset.unique_entries ( ( 'evt' , 'run' ) , choice = 'min'  , criterium = 'PT') :
#  ...
#  @endcode
def _rds_unique_entries_ ( dataset           ,
                           entrytag          ,
                           choice            , 
                           criterium = ''    , 
                           seed      = None  ,
                           progress  = False , 
                           report    = True  ,
                           style     = ''    ) :
    
    if criterium  :
        assert choice in ( 'min' , 'max' , 'minimal' , 'maximal' , 'minimum' , 'maximum' ) , \
               "Invalid 'choice' for criterium!"
    else : 
        assert choice in ( 'first' , 'last' , 'random' , 'rndm'  , 'rand' ) , \
               "Invalid 'choice'"
 
    snapshot = _rds_seek_for_duplicates_ ( dataset               ,
                                           entrytag  = entrytag  ,
                                           criterium = criterium ,
                                           progress  = progress  ) 

    choice = choice.lower()
    
    first  = 'first'  == choice
    last   = 'last'   == choice
    rand   = choice in ( 'random' , 'rndm'  , 'rand' ) 
    minv   = criterium and choice in ( 'min' , 'minimal' , 'minimum' )
    maxv   = criterium and choice in ( 'max' , 'maximal' , 'maximum' )

    with random_seed ( seed ) :

        cnt    = SE ()
        unique = 0 
        for e, lst in progress_bar ( loop_items ( snapshot )        ,
                                     max_value   = len ( snapshot ) ,
                                     description = '2nc loop:'      ,
                                     silent      = not progress     ) : 

            unique += 1
            
            num  = len ( lst )
            if 1 < num : cnt += num 
            
            if   1 == num  : yield lst [   0 ] [ 1 ]        
            elif first     : yield lst [   0 ] [ 1 ]
            elif last      : yield lst [  -1 ] [ 1 ]
            elif rand      : yield random.choice ( lst ) [ 1 ] ## seed is needed here! 
            elif minv      : yield min ( lst ) [ 1 ]
            elif maxv      : yield max ( lst ) [ 1 ]
                
    if report or progress :
        title = 'Unique'
        rows  = [ ( '' , 'value' ) ]
        row   =  'Total       entries'      , '%d' % len ( dataset )
        rows.append ( row )
        row   =  'Unique      events'       , '%d' % unique 
        rows.append ( row )
        row   =  'Events with duplicates'   , '%d' % cnt.nEntries() 
        rows.append ( row )
        row   =  'Duplicated  mean +/- rms' , '%.3f +/- %-.3f' % ( cnt.mean() , cnt.rms() )        
        rows.append ( row )
        row   =  'Duplicated  max'          , '%d' % cnt.max()         
        rows.append ( row )
        import ostap.logger.table as T 
        logger.info ( '%s:\n%s' % ( title , T.table ( rows , title = title , prefix = '# ' , alignment = 'lc' , style = style ) ) )
        
# =============================================================================        
## Make a copy of dataset only with unique  entries 
#  @code
#  dataset = ...
#  unique = dataset.make_unique ( ( 'evt' , 'run' ) , choice = 'random' )
#  unique = dataset.make_unique ( ( 'evt' , 'run' ) , choice = 'first'  )
#  unique = dataset.make_unique ( ( 'evt' , 'run' ) , choice = 'last'    )
#  unique = dataset.make_unique ( ( 'evt' , 'run' ) , choice = 'min' , criterium = 'PT' )
#  unique = dataset.make_unique ( ( 'evt' , 'run' ) , choice = 'max' , criterium = 'PT' )
#  @endcode
#  - CPU performance is more or less reasonable up to dataset with 10^7 entries 
def _rds_make_unique_ ( dataset           ,
                        entrytag          ,
                        choice            , 
                        criterium = ''    , 
                        seed      = None  ,
                        progress  = False , 
                        report    = True  ) :
    """ Make a copy of dataset only with unique  entries 
    >>> dataset = ...
    >>> unique = dataset.make_unique ( ( 'evt' , 'run' ) , choice = 'random' )
    >>> unique = dataset.make_unique ( ( 'evt' , 'run' ) , choice = 'first'  )
    >>> unique = dataset.make_unique ( ( 'evt' , 'run' ) , choice = 'last'    )
    >>> unique = dataset.make_unique ( ( 'evt' , 'run' ) , choice = 'min' , criterium = 'PT' )
    >>> unique = dataset.make_unique ( ( 'evt' , 'run' ) , choice = 'max' , criterium = 'PT' )
    - CPU performance is more or less reasonable up to dataset with 10^7 entries 
    """

    weighted          = dataset.isWeighted                     () 
    store_errors      = weighted and dataset.store_errors      ()
    store_asym_errors = weighted and dataset.store_asym_errors ()

    ds       = dataset.emptyClone()
    for index in dataset.unique_entries ( entrytag  = entrytag  ,
                                          choice    = choice    ,
                                          criterium = criterium ,
                                          seed      = seed      ,
                                          progress  = progress  , 
                                          report    = report    ) :

        entry, weight = dataset [ index ]
        
        if   store_asym_errors and isinstance ( weight , VAE ) :
            ds.add ( entry , weight.value , weight.neg_error , weight.pos_erroe )
        elif store_errors      and isinstance ( weight , VE  ) :
            ds.add ( entry , weight.value () ) 
        elif weighted          and isinstance ( weight , num_types ) : 
            ds.add ( entry , float ( weight ) )
        elif weighted : 
            raise TypeError ( 'Inconsistent sae/se/w %s/%s/%s ' % ( store_asym_errors  ,
                                                                    store_errors       ,
                                                                    type ( weight )    ) ) 
        else  :
            ds.add ( entry  )

    return ds


ROOT.RooAbsData.make_unique     = _rds_make_unique_
ROOT.RooAbsData.unique_entries  = _rds_unique_entries_
ROOT.RooAbsData.duplicates      = _rds_duplicates_ 

_new_methods_ += [
   ROOT.RooAbsData . make_unique    , 
   ROOT.RooAbsData . unique_entries , 
   ROOT.RooAbsData . duplicates     ,
   ]


# =============================================================================
## some decoration over RooDataSet 
ROOT.RooAbsData . varlist       = _rad_vlist_
ROOT.RooAbsData . varlst        = _rad_vlist_
ROOT.RooAbsData . vlist         = _rad_vlist_
ROOT.RooAbsData . vlst          = _rad_vlist_

ROOT.RooAbsData . varset        = lambda s : s.get()
ROOT.RooAbsData . __len__       = lambda s : s.numEntries()
ROOT.RooAbsData . __nonzero__   = lambda s : 0 != len ( s )


ROOT.RooAbsData .vars       = property ( lambda s : s.get() ,  None , None , """Variables (as ROOT.RooArgSet)""" ) 
ROOT.RooAbsData .variables  = property ( lambda s : s.get() ,  None , None , """Variables (as ROOT.RooArgSet)""" ) 
ROOT.RooAbsData .content    = property ( lambda s : s.get() ,  None , None , """Variables (as ROOT.RooArgSet)""" ) 
                                    
ROOT.RooAbsData . __contains__  = _rad_contains_

ROOT.RooAbsData . __iter__      = _rad_iter_ 
ROOT.RooAbsData . __getitem__   = _rad_getitem_

ROOT.RooAbsData . subset        = _rad_subset_ 

ROOT.RooDataSet . __add__       = _rad_add_
ROOT.RooDataSet . __iadd__      = _rad_iadd_

ROOT.RooDataSet . __mul__       = _rad_mul_
ROOT.RooDataSet . __rmul__      = _rad_rmul_
ROOT.RooDataSet . __imul__      = _rad_imul_
ROOT.RooDataSet . __mod__       = _rad_mod_
ROOT.RooDataSet . __div__       = _rad_div_
ROOT.RooDataSet . __truediv__   = ROOT.RooDataSet . __div__


ROOT.RooDataSet . sample        = _rad_sample_
ROOT.RooDataSet . choice        = _rad_choice_
ROOT.RooDataSet . shuffle       = _rad_shuffle_



ROOT.RooDataSet . __sub__       = _rds_sub_
ROOT.RooDataSet . jackknife     = _rds_jackknife_
ROOT.RooDataSet . bootstrap     = _rds_bootstrap_


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
   ROOT.RooDataSet . __mul__       ,
   ROOT.RooDataSet . __rmul__      ,
   ROOT.RooDataSet . __imul__      ,
   ROOT.RooDataSet . __div__       ,
   ROOT.RooDataSet . __mod__       ,
   ROOT.RooDataSet . __truediv__   ,
   #
   ROOT.RooDataSet . sample        ,
   ROOT.RooDataSet . shuffle       ,
   #
   ]


ROOT.RooAbsDataStore . __len__      = lambda s : s.numEntries()
ROOT.RooAbsDataStore . __iter__     = _rad_iter_ 
ROOT.RooAbsDataStore . __contains__ = _rad_contains_ 


_new_methods_ += [
   ROOT.RooAbsDataStore . __len__       ,
   ROOT.RooAbsDataStore . __iter__      , 
   ROOT.RooAbsDataStore . __contains__
   ]
# =============================================================================\
## number of warninf prints 
_printed = 10 
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
#  @see RooDataSet 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-07-06
def ds_project  ( dataset                 ,
                  histo                   ,
                  what                    , 
                  cuts      = ''          , * , 
                  cut_range = ''          , 
                  first     = FIRST_ENTRY ,
                  last      = LAST_ENTRY  ,
                  progress  = False       ) :
    """ Helper project method for RooDataSet/DataFrame/... and similar objects 

    >>> h1   = ROOT.TH1D(... )
    >>> dataset.project ( h1.GetName() , 'm', 'chi2<10' ) ## project variable into histo
    
    >>> h1   = ROOT.TH1D(... )
    >>> dataset.project ( h1           , 'm', 'chi2<10' ) ## use histo
    """
    from ostap.stats.statvars import data_project
    first , last = evt_range ( dataset , first , last )
    
    ## 2) if the histogram is specified by the name, try to locate it in ROOT memory  
    if isinstance ( histo , string_types ) :
        groot = ROOT.ROOT.GetROOT()
        h     = groot.FindObject ( histo )
        assert h , "Cannot get locate histogram by name:`%s'" % histo       
        assert isinstance ( h , ROOT.TH1 ) , "Object `%s' exists, but not ROOT.TH1" % typename ( h ) 
        histo = h
    
    return data_project ( dataset  ,
                          histo     ,
                          what      ,
                          cuts      , first , last  , 
                          cut_range = cut_range ,
                          progress  = progress  ,
                          use_frame = False     ,
                          parallel  = False     ) 

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
def ds_draw ( dataset ,
              what                    , * , 
              cuts      = ''          ,
              opts      = ''          ,
              cut_range = ''          ,
              first     = FIRST_ENTRY ,
              last      = LAST_ENTRY  ,
              delta     = 0.01        ,
              progress  = False       ,
              use_frame = False       ,
              paralell  = False       , **kwargs ) :
    """ Helper draw method for drawing of RooDataSet
    >>> dataset.draw ( 'm', 'chi2<10'                 )
    ## cuts & weight 
    >>> dataset.draw ( 'm', '(chi2<10)*weight'        )
    ## use drawing options 
    >>> dataset.draw ( 'm', '(chi2<10)*weight' , 'e1' )
    ## start form event #1000     >>> dataset.draw ( 'm', '(chi2<10)*weight' , 'e1' , 1000 ) 
    ## for event in range 1000< i <10000
    >>> dataset.draw ( 'm', '(chi2<10)*weight' , 'e1' , 1000 , 100000 )
    """
    first , last = evt_range ( dataset , first , last )
    
    ## check type of opts 
    assert isinstance ( opts , string_types ) , "Invalid type of `opts' : %s" % type ( opts ) 

    ## decode variables/cuts 
    varlst, cuts, input_string = vars_and_cuts  ( what , cuts )
    if input_string and 2 <= len ( varlst ) and order_warning :
        vv = ' ; '.join  ( varlst  ) 
        logger.attention ("draw: from v1.10.1.9 variables are in natural order [x;y;..]=[ %s ]" % vv  )

    nvars = len ( varlst ) 
    assert 1 <= nvars <= 3 , "Invalid number of variables: %s" % str ( varlst )
    
    ## get the suitable ranges for the variables
    if isinstance ( dataset , ROOT.TTree ) :
        assert not cut_range , "ds_draw(ROOT.TTree): cut_range `%s'is not allowed!" % cut_range 
        from ostap.trees.trees import tree_draw 
        return tree_draw ( dataset               ,
                           what                  ,
                           cuts      = cuts      ,
                           opts      = opts      ,
                           first     = first     ,
                           last      = last      ,
                           delta     = delta     ,
                           progress  = progress  ,
                           use_frame = use_frame , 
                           parallel  = parallel  , **kwargs )
    
    elif isinstance ( dataset , ROOT.RooAbsData ) :
        
        ranges = ds_range ( dataset     ,
                            varlst      ,
                            cuts      = cuts      ,
                            cut_range = cut_range ,
                            first     = first     ,
                            last      = last      ,
                            delta     = delta     )
    else :
        
        ## something else ? e.g. DataFrame 
        assert not cut_range                   , "ds_draw: `cut_range' is not allowed!"
        assert ( first , last ) == ALL_ENTRIES , "ds_draw: `first'/`last' are not allowed!"
        ranges = data_range ( dataset , varlst , cuts = cuts , delta = delta )

    if not ranges :
        logger.warning ("ds_draw: nothing to draw, return None" ) 
        return None 
        
    assert len ( ranges ) == nvars , 'Invalid ranges: %s' % str ( ranges )

    from ostap.utils.cidict import cidict, cidict_fun
    kw = cidict ( transform = cidict_fun , **kwargs )

    histos = []
    for var in varlst :
        mn, mx = ranges [ var ]
        item   = var, ( mn, mx) 
        histos.append ( item ) 

    ## book the histogram
    from   ostap.histos.histos       import histo_book2
    histo = histo_book2 ( histos , kw )
    ## fill the histogram
    histo = ds_project ( dataset , histo , varlst , cuts = cuts , cut_range = cut_range , first = first , last = last )
    ## draw the histogram 
    histo.draw ( opts , **kw )
    return histo

# =============================================================================
## get the attibute for RooDataSet
def _ds_getattr_ ( dataset , attname ) :
    """ Get the attibute from RooDataSet 

    >>> dset = ...
    >>> print dset.pt
    
    """
    _vars = dataset.get()
    return getattr ( _vars , attname )  

# =============================================================================
## get the attibute for RooDataSet
# =============================================================================
def get_var ( self, aname ) :
    _vars = self.get()
    return getattr ( _vars , aname )  



# =============================================================================\
## Get suitable ranges for drawing expressions/variables
## @code
#  dataset = ...
#  result  = ds_range ( dataset , 'sin(x)*100*y' , 'x<0' )
#  results = ds_range ( dataset , 'x,y,z,t,u,v'  , 'x<0' )
#  @endcode 
def ds_range  ( dataset                 ,
                expressions             ,
                cuts      = ''          , * , 
                cut_range = ''          ,
                first     = FIRST_ENTRY , 
                last      = LAST_ENTRY  ,
                delta     = 0.05        ,
                progress  = False       ) :
    """ Get suitable ranges for drawing expressions/variables
    >>> dataset = ...
    >>> result  = ds_range ( dataset , 'sin(x)*100*y' , 'x<0' )
    >>> results = ds_range ( dataset , 'x,y,z,t,u,v'  , 'x<0' )
    """
    first , last = evt_range ( dataset , first , last ) 
    return data_range ( dataset     ,
                        expressions ,
                        cuts        , first , last , 
                        cut_range   = cut_range    ,
                        delta       = delta        ,
                        progress    = progress     )

# =============================================================================
## clear dataset storage
if not hasattr ( ROOT.RooDataSet , '_old_reset_' ) :
    ROOT.RooDataSet._old_reset_ = ROOT.RooDataSet.reset
    def _ds_new_reset_ ( self ) :
        """ Clear dataset storage
        >>> print ds
        >>> ds.clear()
        >>> ds.erase() ## ditto
        >>> ds.reset() ## ditto
        >>> ds.Reset() ## ditto
        >>> print ds
        """        
        store = self.store()
        if store :
            store.reset        ()
            store.resetBuffers ()
            store.resetCache   ()
        
        ## self.resetCache    ()
        self.resetBuffers  ()
        self._old_reset_   ()
        return len ( self )
    
    ## ROOT.RooDataSet.reset = _ds_new_reset_
    ROOT.RooDataSet.clear = _ds_new_reset_
    ROOT.RooDataSet.erase = _ds_new_reset_
    ## ROOT.RooDataSet.Reset = _ds_new_reset_


# =============================================================================
## clear dataset
#  @see Ostap::MoreRooFit::reset_data
#  @code
#  data = ...
#  data.clear ()
#  data.clean () ## ditto
#  data.erase () ## ditto
#  @endcode 
def _ds_clear_ ( data ) :
    """ Clear dataset
    - see `Ostap.MoreRooFit.reset_data`
    >>> data = ...
    >>> data.clear ()
    >>> data.clean () ## ditto
    >>> data.erase () ## ditto
    """
    Ostap.MoreRooFit.reset_data ( data )
        
ROOT.RooDataSet.clear = _ds_clear_ 
ROOT.RooDataSet.clean = _ds_clear_ 
ROOT.RooDataSet.erase = _ds_clear_ 

ROOT.RooAbsData.get_var= get_var

_new_methods_ += [
    ROOT.RooDataSet .clear   ,
    ROOT.RooDataSet .erase   ,    
    ROOT.RooAbsData .get_var ,
    ]

# =============================================================================
ROOT.RooDataSet.draw         =  ds_draw
ROOT.RooAbsData.draw         =  ds_draw
ROOT.RooDataSet.project      =  ds_project
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
    """ Get the s-factor for (weighted) dataset, where
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

    if not data.isWeighted() :
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
    result = dataset.print_multiline ( verbose = True )
    return result 
    

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
## clone dataset
#  @code
#  dataset = ...
#  cloned  = datatset.clone ( 'new_name') 
#  @endcode
def _rds_clone_ ( dataset , name = '' ) :
    """ Clone dataset
    >>> dataset = ...
    >>> cloned  = datatset.clone ( 'new_name') 
    """
    name = name if name else dsID ()     
    return ROOT.RooDataSet ( dataset , name ) 

# =============================================================================
## clone dataset
#  @code
#  dataset = ...
#  cloned  = datatset.clone ( 'new_name') 
#  @endcode
def _rdh_clone_ ( dataset , name = '' ) :
    """ Clone dataset
    >>> dataset = ...
    >>> cloned  = datatset.clone ( 'new_name') 
    """
    name = name if name else dsID ()     
    return ROOT.RooDataHist ( dataset , name ) 

if not hasattr ( ROOT.RooDataSet  , 'clone' ) :
    ROOT.RooDataSet .clone = _rds_clone_

if not hasattr ( ROOT.RooDataHist , 'clone' ) :
    ROOT.RooDataHist.clone = _rdh_clone_
    
# =============================================================================
## add variable to dataset
#  @code
#  dataset = ...
#  dataset.addVar ( 'NewVar' , 'A+B/3' )
#  @endcode 
def _rds_addVar_ ( dataset , vname , formula ) : 
    """ Add/calculate variable to RooDataSet

    >>> dataset.addVar ( 'ratio' , 'pt/pz' )
    """
    vcol = make_formula ( vname , formula , dataset.get() )
    dataset.addColumn ( vcol )
    del vcol 
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
#  - Sample from 1D-histogram
#  @code
#  h1 = ...## 1D histogram
#  dataset.add_new_var ( 'nTracks' , h1 ) ## sample from 1D histogram
#  @endcode 
#  - Sample from 2D histogram
#  @code
#  h2 = ...## 2D histogram
#  dataset.add_new_var ( 'Pt' , 'eta' , h2 ) ## sample from 2D histogram
#  @endcode
#  - Sample  from 3D-histogram
#  @code
#  h3 = ...## 3D histogram
#  dataset.add_new_var ( 'Pt' , 'eta' , 'A' , h3 ) ## sample from 3D histogram
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
def add_new_var ( dataset , varname , what , *args , progress = False ) : 
    """ Add/calculate/sample variable to RooDataSet

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
    if isinstance ( varname , string_types ) :
        if    isinstance ( what , string_types   ) : pass
        elif  isinstance ( what , sequence_types )  and \
                 len ( what ) == len ( dataset    ) and \
                 all ( isinstance ( v , num_types ) for v in what ) :
            
            vvar  = ROOT.RooRealVar ( varname , 'variable %s' % varname , -999.999 )
            vset  = ROOT.RooArgSet  ( vvar )
            dset  = ROOT.RooDataSet ( dsID() , 'dataset with %s' % varname , vset )
            for v in what :
                vvar.setVal ( float ( v ) )
                dset.add ( vset )
                
            if len ( dataset ) == len ( dset ) :
                dataset *= dset
                return dataset
            
            dset.clear() 
            del dset
            del vset
            del vvar 
            raise TypeError("Invalid type/length of 'what' argument")

    progress = progress_conf ( progress )
    adder    = Ostap.AddVars ( progress ) 
    vv = adder.add_var ( dataset , varname , what , *args )
    if not  vv : logger.error('add_new_var: NULLPTR from Ostap.AddVars.add_var')
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
def _rds_makeWeighted_ ( dataset           ,
                         weightvar         ,
                         cuts   = ''       ,
                         wname  = 'Weight' ) :
    """ Make weighted data set from unweighted dataset
    
    >>> dataset = ...
    >>> wdata   = dataset.makeWeighted ( 'S_sw' )    
    """
    assert not dataset.isWeighted () , "Dataset '%s/%s' is already weighted!" % ( dataset.GetName  () ,
                                                                                  dataset.GetTitle () )
    
    assert isinstance ( weightvar , expression_types ) , \
        "Invalid type of `weigthvar':%s" % type ( weightvar )
    assert isinstance ( cuts      , expression_types ) or not cuts , \
        "Invalid type of `cuts':%s" % type ( cuts )

    ## 
    weightvar = str ( weightvar ).strip()
    cuts      = str ( cuts      ).strip()
    
    if not weightvar in dataset :
        ## is it a formula ?
        wname = wname or 'Weight'
        while wname in dataset :  wname += 'W'
        dataset.addVar ( wname , weightvar )
        weightvar = wname

    ## content
    varset = dataset.get()

    ## make weighted dataset
    if root_info < ( 6, 36 ) : 
        args = ( dsID()             ,
                 dataset.GetTitle() ,
                 dataset            ,
                 varset             , 
                 cuts               ,
                 weightvar          )
    else :
        args = ( dsID()                              ,
                 dataset.GetTitle()                  ,
                 varset                              , 
                 ROOT.RooFit.Import    ( dataset   ) ,
                 ROOT.RooFit.Cut       ( cuts      ) ,
                 ROOT.RooFit.WeightVar ( weightvar ) ) 
        
    return ROOT.RooDataSet ( *args )

    
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
        return dataset , ''
    
    ds , w = Ostap.Utils.unweight ( dataset , weight )
    
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
def setStorage ( new_type = RAD.Tree , silent = False ) :
    """ Redefine the default storage 
    """
    if not new_type in ( RAD.Tree , RAD.Vector ) :
        logger.error ('RooAbsData: Invalid storage type %s, replace with Tree ' % new_type )
        new_type = RAD.Tree
        
    if RAD.getDefaultStorageType() != new_type :
        if   new_type == RAD.Tree and not silent : 
            logger.info  ( 'RooAbsData: DEFINE default storage type to be RooAbsData.Tree'      )
        elif new_type == RAD.Vector    and not silent : 
            logger.info  ( 'RooAbsData: DEFINE default storage type to be RooAbsData.Vector'    )
        elif new_type == RAD.Composite and not silent : 
            logger.info  ( 'RooAbsData: DEFINE default storage type to be RooAbsData.Composite' )
        elif not silent :
            logger.info  ( 'RooAbsData: DEFINE default storage type to be %s' % new_type        )
            
        RAD.setDefaultStorageType ( new_type  ) 

    the_type = RAD.getDefaultStorageType()
    if   RAD.Tree   == the_type : logger.debug ( 'RooAbsData: Default storage type is Tree'   )
    elif RAD.Vector == the_type : logger.debug ( 'RooAbsData: Default storage type is Vector' )
    else : logger.debug ( 'RooAbsData: Default storage type is %s' % the_type  )

# =============================================================================
## @class UseStorage
#  Context manager to change the storage type
class UseStorage(object) :
    """ Context manager to change the storage type
    >>> with UseStorage() :
    ...
    """
    def __init__  ( self , new_storage = RAD.Tree , silent = True ) :
        if not new_storage in ( RAD.Tree , RAD.Vector )  :
            raise AttributeError( 'Invalid storage type %s' % new_storage )
        self.new_storage = new_storage
        self.old_storage = RAD.getDefaultStorageType()
        self.silent      = silent
    def __enter__ ( self ) :
        self.old_storage = RAD.getDefaultStorageType()
        setStorage (  self.new_storage , silent = self.silent )
    def __exit__ (  self , *_ ) :
        setStorage (  self.old_storage , silent = self.silent )

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


# =============================================================================
## Print data set as table
def _ds_table_0_ ( dataset                 ,
                   variables = []          ,
                   cuts      = ''          ,
                   cut_range = ''          , 
                   first     = FIRST_ENTRY ,
                   last      = LAST_ENTRY  ,
                   prefix    = ''          ,
                   title     = ''          ,
                   style     = ''          ) :
    """ Print data set as table
    """
    varset = dataset.get()
    if not valid_pointer ( varset ) :
        logger.error('Invalid dataset')
        return ''

    assert isinstance ( cuts    , expression_types  ) or not cuts, \
        "Invalid type of cuts: %s" % type ( cuts )
    assert isinstance ( cut_range , expression_types ) or not cut_range, \
        "Invalid type of cut_range: %s" % type ( cut_range )
    
    cuts      = str(cuts).strip()      if cuts      else ''
    cut_range = str(cut_range).strip() if cut_range else ''
    
    if variables : vars , cuts , _ = vars_and_cuts ( variables , cuts ) 
    else         : vars , cuts     = [ v.name for v in varset ] , str ( cuts ).strip() 

    ## adjust first/last 
    first, last = evt_range ( dataset , first , last )
    
    vars = [ i.GetName() for i in varset if i.GetName() in vars ]
        
    _vars = []
    
    stat = dataset.statVars ( vars , cuts , first , last , cut_range = cut_range )
    for v in  stat :
        vv  = getattr ( varset , v )
        s   = stat [ v ] 
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
    if not title :    
            
        if  tt and tt != dataset.GetName()  : 
            title = '%s("%s","%s"):' % ( typename ( dataset ) , dataset.GetName () , tt ) 
        else :
            title = '%s("%s"):'      % ( typename ( dataset ) , dataset.GetName () )
        title =  '%s %d#' %  ( title , len ( dataset ) )
        if cabinet : title = '%s %s' % ( cabinet , title )
        
    if not _vars :
        return title , 120 
        ## return report , 120 

    weight = None
    if   isinstance ( dataset , ROOT.RooDataHist ) :
        
        if dataset.isNonPoissonWeighted() : title += ' Binned/Weighted' 
        else                              : title += ' Binned'
        
    elif dataset.isWeighted () :
        
        if dataset.isNonPoissonWeighted() : title += ' Weighted' 
        else :                              title += ' Weighted/Poisson'

        ## name of weight variabe
        weight = dataset.wname ()
        wcnt   = dataset.statVar ( '1' , cuts , first , last , cut_range = cut_range )
        wcnt   = wcnt.weights ()            
        r    = (  weight                            ,   ## 0 
                  'Weight variable'                 ,   ## 1 
                  ('%+.5g' % wcnt.mean().value() ).strip() ,   ## 2
                  ('%.5g'  % wcnt.rms ()         ).strip() ,   ## 3 
                  ('%+.5g' % wcnt.min ()         ).strip() ,   ## 4
                  ('%+.5g' % wcnt.max ()         ).strip() )   ## 5
        _vars.append ( r ) 
        with_weight = True

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

    if weight : header.append ( weight_lifter+' ' if weight_lifter else 'W' )
        
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
            cols.append ( weight_lifter+' ' if weight_lifter else 'W' )
            cols = [ allright ( c ) for c in cols ] 
        elif weight      : cols.append ( ' ' )
        
        table_data.append ( tuple ( cols ) ) 

    import ostap.logger.table as T
    t  = T.table ( table_data , title = title , prefix =  prefix , style = style )
    w  = T.table_width ( t ) 
    return t , w 


# =============================================================================
## Print data set as table
def _ds_table_1_ ( dataset                 ,
                   variables = []          ,
                   cuts      = ''          ,
                   cut_range = ''          , 
                   first     = FIRST_ENTRY ,
                   last      = LAST_ENTRY  ,
                   prefix    = ''          ,
                   title     = ''          ,
                   style     = ''          ) :
    """ Print data set as table
    """
    
    first, last = evt_range ( dataset , first , last ) 
    
    assert isinstance ( cuts    , expression_types  ) or not cuts, \
        "Invalid type of cuts: %s" % type ( cuts )
    assert isinstance ( cut_range , expression_types ) or not cut_range, \
        "Invalid type of cut_range: %s" % type ( cut_range )
    
    cuts      = str(cuts).strip()      if cuts      else ''
    cut_range = str(cut_range).strip() if cut_range else ''

    varset = dataset.get() 
    if variables : vars , cuts , _ = vars_and_cuts ( variables  , cuts ) 
    else         : vars , cuts     = [ v.name for v in varset ] , cuts 

    _vars = []    
    vvars = tuple ( sorted ( vars ) )
    stat = dataset.statVars ( vvars  , cuts , first , last , cut_range = cut_range ) 
    for v in  stat :
        s   = stat [ v ] 
        mnmx = s.minmax ()
        mean = s.mean   ()
        rms  = s.rms    ()

        r    = ( v                                 ,   ## 0 
                 ('%+.5g' % mean.value() ).strip() ,   ## 1
                 ('%.5g'  % rms          ).strip() ,   ## 2 
                 ('%+.5g' % mnmx[0]      ).strip() ,   ## 3
                 ('%+.5g' % mnmx[1]      ).strip() )   ## 4            
        _vars.append ( r )
    _vars.sort() 
        
    tt = dataset.GetTitle()
    if not title :        
        title = '%s("%s"):'     % ( typename ( dataset ) , dataset.GetName () )
        title = '%s %d entries' % ( title , len ( dataset ) )

    if not _vars :
        return title , 120 
        ## return report , 120 

    weight = None

    # ==============================================================================================
    # build the actual table 
    # ==============================================================================================
    
    name_l  = len ( 'Variable'    ) + 2 
    mean_l  = len ( 'mean' ) + 2 
    rms_l   = len ( 'rms'  ) + 2
    min_l   = len ( 'min'  ) + 2 
    max_l   = len ( 'max'  ) + 2 
    for v in _vars :
        name_l = max ( name_l , len ( v[0] ) )
        mean_l = max ( mean_l , len ( v[1] ) )
        rms_l  = max ( rms_l  , len ( v[2] ) )
        min_l  = max ( min_l  , len ( v[3] ) )
        max_l  = max ( max_l  , len ( v[4] ) )
        
    index_l =   int ( math.ceil ( math.log10( len ( _vars ) + 1 ) ) )
    
    fmt_name = '%%%ds. %%-%ds' % ( index_l , name_l )
    fmt_mean = '%%%ds'  % mean_l
    fmt_rms  = '%%-%ds' % rms_l
    fmt_min  = '%%%ds'  % min_l
    fmt_max  = '%%-%ds' % max_l

    title_l = index_l + 2 + name_l  
    header = [ ( '{:^%d}' % title_l ).format ( 'Variable'    ) ,
               ( '{:^%d}' % mean_l  ).format ( 'mean'        ) ,
               ( '{:^%d}' % rms_l   ).format ( 'rms'         ) ,
               ( '{:^%d}' % min_l   ).format ( 'min'         ) ,
               ( '{:^%d}' % max_l   ).format ( 'max'         ) ]

    table_data = [ tuple  ( header ) ]

    vlst = vars

    for i , v in enumerate ( _vars ) :
                
        cols = [ fmt_name %  ( i + 1 , v [ 0 ] ) ,
                 fmt_mean %            v [ 1 ] ,
                 fmt_rms  %            v [ 2 ] ,
                 fmt_min  %            v [ 3 ] ,
                 fmt_max  %            v [ 4 ] ]
        
        table_data.append ( tuple ( cols ) ) 

    title = title    
    import ostap.logger.table as T
    t  = T.table ( table_data , title = title , prefix =  prefix , style = style )
    w  = T.table_width ( t ) 
    return t , w 

# ==============================================================================
## print dataset in a form of the table
#  @code
#  dataset = ...
#  print dataset.table() 
#  @endcode
def _ds_table_ (  dataset                 ,
                  variables = []          ,
                  cuts      = ''          ,
                  cut_range = ''          ,
                  first     = FIRST_ENTRY ,
                  last      = LAST_ENTRY  ,                 
                  prefix    = ''          ,
                  title     = ''          ,
                  style     = ''          ) :
    """ Print dataset in a form of the table
    >>> dataset = ...
    >>> print dataset.table()
    """
    return _ds_table_0_ ( dataset               ,
                          variables             ,
                          cuts      = cuts      ,
                          cut_range = cut_range ,
                          first     = first     ,
                          last      = last      , 
                          prefix    = prefix    ,
                          title     = title     ,
                          style     = style     ) [ 0 ]

# ==============================================================================
## print dataset in a form of the table
#  @code
#  dataset = ...
#  print dataset.table2() 
#  @endcode
def _ds_table2_ (  dataset                 ,
                   variables = []          ,
                   cuts      = ''          ,
                   cut_range = ''          ,
                   first     = FIRST_ENTRY ,
                   last      = LAST_ENTRY  ,                 
                   prefix    = ''          ,
                   title     = ''          ,
                   style     = ''          ) :
    """ Print dataset in a form of the table
    >>> dataset = ...
    >>> print dataset.table()
    """
    return _ds_table_1_ ( dataset               ,
                          variables             ,
                          cuts      = cuts      ,
                          cut_range = cut_range ,
                          first     = first     ,
                          last      = last      ,
                          prefix    = prefix    ,
                          title     = title     ,
                          style     = style     ) [ 0 ]

# =============================================================================
##  print DataSet
def _ds_print2_ ( dataset ) :
    """ Print dataset"""
    if dataset.isWeighted() and not isinstance ( dataset , ROOT.RooDataHist ) :
        store = dataset.store()
        if valid_pointer ( store ) and isinstance ( store , ROOT.RooTreeDataStore ) : pass
        else : return _ds_print_ ( dataset )         
    from ostap.utils.basic import terminal_size, isatty 
    if not isatty() : return _ds_table_ ( dataset )
    tw  , th  = terminal_size()
    rep , wid = _ds_table_0_ ( dataset ) 
    if wid < tw  : return rep
    return _ds_print_ ( dataset )


for t in ( ROOT.RooDataSet , ROOT.RooDataHist ) :
    t.__repr__    = _ds_print2_
    t.__str__     = _ds_print2_
    t.table       = _ds_table_
    t.table2      = _ds_table2_
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
#  ds_sym = ds.symmetrize ( [ 'var1' , 'var2'  ] )
#  ds_sym = ds.symmetrize ( [ 'var1' , 'var2' , 'var3' ] )
#  @endcode
def _ds_symmetrize_ ( dataset , variables ) :
    """ Make symmetrization/randomization of the dataset
    >>> ds     = ...
    >>> ds_sym = ds.symmetrize ( [ 'var1' , 'var2' ] )
    >>> ds_sym = ds.symmetrize ( [ 'var1' , 'var2' , 'var3' [ )
    """

    varlst , _ , _ = vars_and_cuts ( variables , '' )

    varset = dataset.get() 
    extra  = [ v for v in varlst if not v in varset ]
    assert 2 <= len ( varlst ) and not extra , 'Variables are not in dataset: %s' % str ( extra ) 

    
    nvars  = [ varset [ v ] for v in varlst ]
    
    mnv    = min ( [ v.getMin () for v in nvars if hasattr ( v , 'getMin' ) ] ) 
    mxv    = max ( [ v.getMax () for v in nvars if hasattr ( v , 'getMax' ) ] ) 

    names   = tuple ( v.name for v in nvars )

    weighted          = dataset.isWeighted                     () 
    store_errors      = weighted and dataset.store_errors      ()
    store_asym_errors = weighted and dataset.store_asym_errors ()
    
    newds   = dataset.emptyClone ()
    nvarset = newds.varset    ()    
    for v in nvarset :
        if v.name in varlst : 
            if hasattr ( v ,  'setMin' ) : v.setMin ( mnv )
            if hasattr ( v ,  'setMax' ) : v.setMax ( mxv )        

    ## loop over the data set 
    for entry , weight in dataset :

        values = [ v.getVal() for v in entry if v.name in names ]        
        random.shuffle ( values )
        
        for v in nvarset :
            n = v.name 
            if not n in names : v.setVal ( entry [ n ] .value )                
            else              : v.setVal ( values.pop()   )
            
        if   store_asym_errors and isinstance ( weight , VAE ) :
            newds.add ( nvarset , weight.value , weight.neg_error , weight.pos_erroe )
        elif store_errors      and isinstance ( weight , VE  ) :
            newds.add ( nvarset , weight.value () ) 
        elif weighted          and isinstance ( weight , num_types ) : 
            newds.add ( nvarset , float ( weight ) )
        elif weighted : 
            raise TypeError ( 'Inconsistent sae/se/w %s/%s/%s ' % ( store_asym_errors  ,
                                                                    store_errors       ,
                                                                    type ( weight )    ) ) 
        else  :
            newds.add ( nvarset )
        
    return newds


ROOT.RooDataSet.symmetrize = _ds_symmetrize_

_new_methods_ += [
    ROOT.RooDataSet.symmetrize , 
    ]


# =============================================================================
## get the name of weight variable in dataset
#  @code
#  dataset = ...
#  wname   = dataset.wname() 
#  @endcode 
#  @see Ostap::Utils::getWeight
def _ds_wname_ ( dataset ) :
    """ Get the name of weigth variable in dataset
    >>> dataset = ...
    >>> wname   = dataset.wname() 
    """
    if not dataset.isWeighted() : return '' ## UNWEIGHTED!
    ## 
    attr = '_weight_var_name'
    if not hasattr ( dataset , attr ) :
        wn = Ostap.Utils.getWeight ( dataset )    
        setattr ( dataset , attr , wn ) 
    ## 
    return getattr ( dataset , attr , '' )
# =============================================================================

# =============================================================================
## Get the weight variable from the dataset
#  @code
#  dataset = ...
#  wvar    = dataset.weightVar() 
#  @endcode 
def _ds_weight_var_ ( dataset ) :
    """ Get the weight variable from the dataset

    >>> dataset = ...
    >>> wvar    = dataset.weightVar()
    
    """
    wvar = Ostap.Utils.getWeightVar ( dataset )
    return wvar 

# =============================================================================
## Are weight errors stored in dataset?
#  @code
#  dataset     = ...
#  store_error = dataset.store_error () 
#  @endcode
#  The function checks the <code>StoreError</code> and 
#   <code>StoreAsymError</code> attributes for the weight variable 
#  @see Ostap::Utils::storeError
def _ds_store_error_ ( dataset ) :
    """ Are weight errors stored in dataset?
    >>> dataset     = ...
    >>> store_error = dataset.store_error () 
    The function checks the `StoreError` and 
    `StoreAsymError` attributes for the weight variable 
    - see Ostap::Utils::storeError
    """
    ## 
    if not dataset.isWeighted() : return False ## UNWEIGHTED!
    ## 
    attr = '_store_weight_error_'
    if not hasattr ( dataset , attr ) :        
        wn = Ostap.Utils.storeError  (  dataset )
        wn = True if wn else False        
        setattr ( dataset , attr , wn ) 
        
    return getattr ( dataset , attr , '' )
# =============================================================================

# =============================================================================
## Are asymmetric weight errors stored in dataset?
#  @code
#  dataset     = ...
#  store_error = dataset.store_asym_error () 
#  @endcode
#  The function checks the <code>StoreAsymError</code> attribute for the weight variable 
#  @see Ostap::Utils::storeError
def _ds_store_asym_error_ ( dataset ) :
    """ Are weight errors stored in dataset?
    >>> dataset     = ...
    >>> store_error = dataset.store_asym_error () 
    The function checks the `StoreAsymError` attributes for the weight variable 
    - see Ostap::Utils::storeAsymError
    """
    
    if not dataset.isWeighted() : return False ## UNWEIGHTED!

    attr = '_store_asym_weight_error_'
    if not hasattr ( dataset , attr ) :
        ## 
        wn = Ostap.Utils.storeAsymError  (  dataset )
        wn = True if wn else False
        ## 
        setattr ( dataset , attr , wn ) 
                
    return getattr ( dataset , attr , '' )

# =============================================================================

ROOT.RooDataSet.wname             = _ds_wname_
ROOT.RooDataSet.weight_name       = _ds_wname_
ROOT.RooAbsData.store_error       = _ds_store_error_
ROOT.RooAbsData.store_errors      = _ds_store_error_
ROOT.RooAbsData.store_asym_error  = _ds_store_asym_error_
ROOT.RooAbsData.store_asym_errors = _ds_store_asym_error_

if not hasattr ( ROOT.RooDataSet , 'weightVar' ) :
    ROOT.RooDataSet.weightVar = _ds_weight_var_ 

_new_methods_ += [
    ROOT.RooDataSet.wname             , 
    ROOT.RooDataSet.weight_name       , 
    ROOT.RooAbsData.store_error       , 
    ROOT.RooAbsData.store_errors      , 
    ROOT.RooAbsData.store_asym_error  ,
    ROOT.RooAbsData.store_asym_errors ,
    ]

def f_open ( name , mode , **kwargs ) :
    return open ( name , mode , **kwargs )

# =============================================================================
## Convert dataset to CSV format
#  @code
#  data = ...
#  data.cvs ( 'data.csv' )
#  data.cvs ( 'data.csv' , dialect = 'unix'      )
#  data.cvs ( 'data.csv' , dialect = 'excel'     )
#  data.cvs ( 'data.csv' , dialect = 'excel-tab' )
#  data.cvs ( 'data.csv' , vars= ( 'a' , 'b' ) )    ## only subset of variables 
#  data.cvs ( 'data.csv' , more_vars = ( 'a+b/c' , 'sin(a)/b' ) ) ## more variables 
#  @endcode 
def ds_to_csv ( dataset , fname , vars = () , more_vars = () , weight_var = '' , progress = False , mode = 'w' , **kwargs ) :
    """ Convert dataset to CSV format
    >>> data = ...
    >>> data.cvs ( 'data.csv' )
    >>> data.cvs ( 'data.csv' , dialect = 'unix'      )
    >>> data.cvs ( 'data.csv' , dialect = 'excel'     )
    >>> data.cvs ( 'data.csv' , dialect = 'excel-tab' )
    >>> data.cvs ( 'data.csv' , vars= ( 'a' , 'b' ) )    ## only subset of variables 
    >>> data.cvs ( 'data.csv' , more_vars = ( 'a+b/c' , 'sin(a)/b' ) ) ## add more derived variables 
    """

    vvars = vars if vars else [ v for v in dataset.varlist() ]
    
    for v in vvars :
        assert v in dataset, 'ds_to_csv: variable %s is not in dataset' % v

    # more variables 
    mvars = []
    for v in more_vars :
        if   isinstance ( v , ROOT.RooAbsReal ) : mvars.append ( v )
        elif isinstance ( v , string_types    ) :            
            vvar = make_formula ( vv , vv , dataset.get() )
            mvars.append ( vvar )
        else :
            raise TypeError('ds_to_csv: invalid variable %s/%s' % ( v , type(v) ) ) 

    vnames1 = []
    for v in vvars :
        if isinstance ( v , ROOT.RooAbsReal ) : vnames1.append ( v.name )
        else                                  : vnames1.append ( v      )

    vnames = vnames1 + [ v.name for v in mvars ]

    weighted = dataset.isWeighted  ()
    se       = weighted and dataset.store_error      ()
    sae      = weighted and dataset.store_asym_error ()

    if weighted and not weight_var :
        weight_var = dataset.wname()
                
    if   weighted and sae : vnames += [ weight_var , '%sErrorLow' % weight_var , '%sErrorHigh' % weigth_var  ]
    elif weighted and se  : vnames += [ weight_var , '%sError'    % weight_var ]
    elif weighted         : vnames += [ weight_var ]

    import csv
    from   ostap.utils.progress_bar  import progress_bar 

    with f_open ( fname , mode , newline = '' ) as csv_file :
        writer = csv.writer ( csv_file , **kwargs )
        ## write header row 
        writer.writerow ( vnames  )

        ## loop over entries in the dataset
        for entry, _  in progress_bar ( dataset , max_value = len ( dataset ) , silent = not progress , description = 'Entries:' ) :
            
            values =  [ entry [ a ].getVal() for a in vnames1 ]
            values += [ v.getVal()           for v in mvars   ]

            if  sae :
                e1 , e2 = dataset.weight_errors () 
                values += [ dataset.weight() , e1 , e2 ]
            elif se  :
                values += [ dataset.weight() , dataset.weightError () ]
            elif weighted :
                values += [ dataset.weight() ]
                                
            writer.writerow ( values )

# =============================================================================
ROOT.RooDataSet.as_csv            = ds_to_csv
ROOT.RooDataSet.to_csv            = ds_to_csv
ROOT.RooDataSet.toCsv             = ds_to_csv
ROOT.RooDataSet.asCsv             = ds_to_csv

_new_methods_ += [
    ROOT.RooDataSet.as_csv ,
    ROOT.RooDataSet.to_csv ,
    ROOT.RooDataSet.toCsv  , 
    ROOT.RooDataSet.asCsv
    ]

# ============================================================================
## Get dataset as <code>TTree</code>
#  - for Tree-based datasets gets the internal tree
#  - otherwise tree will be created
#  @code
#  dataset = ...
#  tree    = dataset.asTree() 
#  @endcode
def ds_to_tree ( dataset , filename = '' , silent = True ) :
    """ Get dataset as `ROOT.TTree`
    - for Tree-based datasets gets the internal tree
    - otherwise tree will be created
    
    >>> dataset = ...
    >>> tree    = dataset.asTree() 
    """
    
    import ostap.io.root_file
    import ostap.trees.trees

    ## the simplest case
    dataset.convertToTreeStore();

    tree = dataset.tree()
    if valid_pointer ( tree ) :  return tree
    
    return dataset.GetClonedTree()
        
    
    print ( 'T-TREE/1' )
    store = dataset.store()
    if store and isinstance ( store , ROOT.RooTreeDataStore ) :
        print ( 'T-TREE/2.1' )
        tree = store.tree()
        if valid_pointer ( tree ) :
            print ( 'found-TREE/1', type ( tree ) )                 
            return tree
        print ( 'T-TREE/2.2' )
        
    if not filename :
        import ostap.utils.cleanup as CU 
        filename = CU.CleanUp.tempfile ( suffix = '.root' )
        if not silent : logger.info ( "Temporary ROOT file is created: %s" % filename ) 
        
    print ( 'T-TREE/3', filename )
    with useStorage ( RAD.Tree ) :
        dataset.convertToTreeStore ()
        store = dataset.store()        
        print ( 'T-TREE/4', filename )
        if store and isinstance ( store , ROOT.RooTreeDataStore ) :
            print ( 'T-TREE/5', filename )
            tree = store.tree()
            print ( 'T-TREE/6', filename )
            if valid_pointer ( tree ) :
                print ( 'T-TREE/7', filename )
                print ( 'found-TREE/1', type ( tree ) )                 
                return tree
            
    print ( 'T-TREE/8', filename )
    return dataset.GetClonedTree()
            
# ============================================================================

ROOT.RooDataSet.as_tree            = ds_to_tree
ROOT.RooDataSet.to_tree            = ds_to_tree
ROOT.RooDataSet.toTree             = ds_to_tree
ROOT.RooDataSet.asTree             = ds_to_tree


_new_methods_ += [
    ROOT.RooDataSet.as_tree ,
    ROOT.RooDataSet.to_tree ,
    ROOT.RooDataSet.toTree  ,
    ROOT.RooDataSet.asTree
    ]

    
# =============================================================================d
## Get "slice" from <code>RooAbsData</code> in form of numpy array
#  @code
#  data = ...
#  varr , weights = data.slice ( 'a : b : c' , 'd>0' )
#  varr , weights = data.slice ( 'a , b , c' , 'd>0' )
#  varr , weights = data.slice ( 'a ; b ; c' , 'd>0' )
#  @endcode
def ds_slice ( data                     ,
               expressions              ,
               cuts       = ''          ,
               cut_range  = ''          ,
               structured = True        ,
               transpose  = True        , 
               first      = FIRST_ENTRY ,
               last       = LAST_ENTRY  ) :
    
    """ Get "slice" from <code>RooAbsData</code> in form of numpy array
    >>> data = ...
    >>> varr , weights = data.slice ( 'a : b : c' , 'd>0' )
    >>> varr , weights = data.slice ( 'a , b , c' , 'd>0' )
    >>> varr , weights = data.slice ( 'a ; b ; c' , 'd>0' )
    """
    
    ## adjust first/last indices 
    first , last = evt_range ( data , first , last ) 
    if last <= first :
        return () , None 
    
    ## decode cuts & the expressions 
    varlst, cuts , _ = vars_and_cuts  ( expressions , cuts )
    

    ## (4) display progress ? 
    progress = progress_conf ( progress )
    
    ## (5) create the driver 
    sv = StatVar ( progress )

    ## 
    table = Ostap.StatVar.Table  ()
    for v in varlst : table [ v ] 
    
    ## get data
    with rootException() :        
        sc = Ostap.StatVar.get_table ( dataset   ,
                                       table     ,
                                       cuts      ,
                                       cut_range ,
                                       first     ,
                                       last      )        
        assert sc.isSuccess () , "Error code from Ostap::StatVar::get_table %s" % s
        
    result = []
    for var in varlst :
        column = table [ var ]
        result.append  ( numpy.asarray ( column , dtype = numpy.float64 ) )
        column.clear   () 
        table.erase    ( var ) 

    if not result :
        return () , None 
    
    weights = None 
    if cuts or data.isWeighted() :
        assert 1 == table.size()  , "Here table must have size equal to 1!"
        for key in table :
            column  = table [ key ]        
            weights = numpy.asarray ( column , dtype = numpy.float64 ) 
            result.append ( weights ) 
            column.clear () 
        del table
        
    if structured :
        
        dt    = numpy.dtype (  ( v , numpy.float64 ) for v in varlst )
        dt    = numpy.dtype ( dt )
        part  = numpy.zeros ( num , dtype = dt )        
        for v in result : part [ v ] = result [ v ]
        result = part
        
    else :

        if weights :
            weights = result [ -1]
            result  = result [:-1]
            if numpy.all ( weights == 1 ) : weights = None 
   
        result = numpy.stack ( result )
        if transpose : result = numpy.transpose ( result )

    return result, weights 


ROOT.RooAbsData.slice = ds_slice 

_new_methods_ += [
    ROOT.RooAbsData.slice , 
    ds_slice             
    ]

# =============================================================================
## Combine two datasets with some weights
#  @code
#  dataset1 = ... 
#  dataset2 = ...
#  dataset  = ds_combine ( dataset1 , dataset2 , 1.0 , -0.1 )
#  @endcode
#  - Input datasets may be weighted
#  - Output dataset  is weighted 
def ds_combine ( ds1 , ds2 , r1 , r2 , weight = '' , silent = False , title = '' , logger = logger ) :
    """ Combine two datasets with some weights
    >>> dataset1 = ... 
    >>> dataset2 = ...
    >>> dataset  = ds_combine ( dataset1 , dataset2 , 1.0 , -0.1 ) 
    - Input datasets may be weighted
    - Output dataset  is always weighted 
    """

    r1 = float ( r1 )
    r2 = float ( r2 )

    w1 , w2  = '', ''

    ## number of entries 
    n1  , n2  = len      ( ds1 ) , len ( ds2 )
    ## statistics of weights
    st1 , st2 = ds1.statVar('1') , ds2.statVar('1')
    ## sum of weights
    sw1 , sw2 = VE ( st1.sum()   , st1.sumw2() ) , VE ( st2.sum() , st2.sumw2() )
    ## s-factor
    sf1 , sf2 = st1.sum() / st1.sumw2() , st2.sum() / st2.sumw2()
    
    ## 
    if ds1.isWeighted() : ds1 , w1 = ds1.unWeighted ()
    if ds2.isWeighted() : ds2 , w2 = ds2.unWeighted ()
    ##
    if w2 and not w2 in ds1 : ds1.addVar ( w2 , '1' )
    if w1 and not w1 in ds2 : ds2.addVar ( w1 , '1' )
    ## 
    v1 = set ( ds1.branches () )
    v2 = set ( ds2.branches () )
    ##
    if v1 != v2 :

        v12 = v1 - v2
        if v12 : logger.warning ("ds_combine: first  dataset contains %d extra columns: %s" % ( len ( v12 ) , list ( v1 - v2 ) ) )
        
        v21 = v2 - v1
        if v21 : logger.warning ("ds_combine: second dataset contains %d extra columns: %s" % ( len ( v21 ) , list ( v2 - v1 ) ) ) 

        cv = v1.intersection ( v2 )

        assert cv, 'ds_combine: datasets have no common columns!'
        
        ds1s = ds1.subset ( cv )
        ds2s = ds2.subset ( cv )
        
        if w1 : ds1.clear () 
        if w2 : ds2.clear ()
        
        ds1 = ds1s        
        ds2 = ds2s

    ## construct the name for new common weight variable
    new_weight = weight if weight else 'weight'
    i = 0 
    while ( new_weight in ds1 ) or ( new_weight in ds2 )  :
        new_weight = "%s_%d" % ( new_weight , i )
        i += 1

    ## define new weigths for datasets
    if 1 == r1 : weight1 =       '%s' %        w1   if w1 else '1'
    else       : weight1 = '%.16g*%s' % ( r1 , w1 ) if w1 else '%.16g' % r1
    if 1 == r2 : weight2 =       '%s' %        w2   if w2 else '1'
    else       : weight2 = '%.16g*%s' % ( r2 , w2 ) if w2 else '%.16g' % r2
    
    ds1.addVar ( new_weight , weight1 )
    ds2.addVar ( new_weight , weight2 )
    
    ## add two datasets together 
    ds = ds1 + ds2

    ## apply common weight 
    dsw = ds.makeWeighted ( new_weight )

    ## statistics 
    n , sw , sf = len ( dsw ) , dsw.sumVar('1') , dsw.sFactor() 

    if not silent :

        rows = [  ( ''        , 'A' , 'B' , '%+.6g*A%+.6g*B' % ( r1 , r2 ) ) ] 

        row  =  'Size'        , '%d' % n1 , '%s' % n2  , '%d' % n
        rows.append ( row )

        if w1 or w2 :
            row  = 'Weight      (original)' , w1 , w2 , ''
            rows.append ( row )
            
        row = 'Weight      (updated)' , new_weight , new_weight , new_weight
        rows.append ( row )
                
        row  =  'Sum weights (original) ' ,  \
               ( "%+13.6g +- %-13.6g" % ( sw1.value()  , sw1.error() ) ) , \
               ( "%+13.6g +- %-13.6g" % ( sw2.value()  , sw2.error() ) ) , '' 
        rows.append ( row )

        st1n , st2n = ds1.statVar ( new_weight ) , ds2.statVar( new_weight )
        
        sw1n , sw2n = VE ( st1n.sum()   , st1n.sumw2() ) , VE ( st2n.sum() , st2n.sumw2() )
        row  =  'Sum weights (updated) ' , \
               ( "%+13.6g +- %-13.6g" % ( sw1n.value() , sw1n.error() ) ) , \
               ( "%+13.6g +- %-13.6g" % ( sw2n.value() , sw2n.error() ) ) , \
               ( "%+13.6g +- %-13.6g" % ( sw  .value() , sw  .error() ) ) ,
        
        rows.append ( row )
        
        mw1 , mw2 , mw = st1.mean() , st2.mean() ,  sw / n
        row  =  'Mean weight (original)' , \
               ( "%+13.6g +- %-13.6g" % ( mw1.value() , mw1.error() ) ) , \
               ( "%+13.6g +- %-13.6g" % ( mw2.value() , mw2.error() ) ) , ''               
        rows.append ( row )

        mw1n , mw2n = st1n.mean () , st2n.mean ()  
        row  =  'Mean weight (updated)' , \
               ( "%+13.6g +- %-13.6g" % ( mw1n.value() , mw1n.error() ) ) , \
               ( "%+13.6g +- %-13.6g" % ( mw2n.value() , mw2n.error() ) ) , \
               ( "%+13.6g +- %-13.6g" % ( mw  .value() , mw  .error() ) ) 
        rows.append ( row )

        row  =  's-factor    (original)' , "%+13.6g" % sf1 , "%+13.6g" % sf2 , ''
        rows.append ( row )
        
        ## s-factor
        sf1n , sf2n = st1n.sum() / st1n.sumw2() , st2n.sum() / st2n.sumw2()
        
        row  =  's-factor    (updated)' , "%+13.6g" % sf1n , "%+13.6g" % sf2n , "%+13.6g" % sf
        rows.append ( row )

        row  =  'R'                  , "%+13.6g" % r1  , "%+13.6g" % r2  , '' 
        rows.append ( row )
        
        import ostap.logger.table as Table
        title = title if title else 'Combine two datasets: %+.6g*A%+.6g*B' % ( r1 , r2 )
        table = Table.table ( rows               ,
                              title     = title  ,
                              alignment = 'lccc' ,
                              prefix    = "# "   )
        logger.info ( '%s:\n%s' % ( title  , table ) )
    
    ## cleanup
    
    if w1 or v1 != v2 : ds1.clear()
    if w2 or v1 != v2 : ds2.clear()
    
    ds.clear()

    
    return dsw

# ============================================================================
## binned dataset ?
#  @code
#  ds = ...
#  ds.binned ()     ## ditto 
#  @endcode 
def _ds_binned_ ( ds ) :
    """ Binned dataset ?
    >>> ds = ...
    >>> ds.binned ()     ## ditto 
    """
    return isinstance ( ds , ROOT.RooDataHist ) 
# ============================================================================

ROOT.RooAbsData.binned  = _ds_binned_

_new_methods_ += [
    ROOT.RooAbsData.binned 
    ]

# ============================================================================
## Are two datastes equal by content?
#  @code
#  ds1 = ...
#  ds2 = ...
#  ds_equal ( ds1 , ds2 )
#  ds1 == ds2 
#  @endcode 
def ds_equal ( ds1 , ds2 ) : 
    """ Are two datasets equal by content?
    ds1 = ...
    ds2 = ...
    ds_equal ( ds1 , ds2 ) 
    ds1 == ds2 
    """

    if ds1 is ds2 : return True    
    ## same length ? 
    if len ( ds1 ) != len ( ds2 ) :
        logger.debug ("compare datasets: different lengths") 
        return False
    
    ## same variables ? 
    vars1 = set ( ( v.name for v in ds1.varset () ) ) 
    vars2 = set ( ( v.name for v in ds2.varset () ) )
    if vars1 != vars2 :
        logger.debug ("compare datasets: different variables") 
        return False

    ## both binned or unbinned?
    b1 = ds1.binned ()
    b2 = ds2.binned ()
    if b1 and not b2 :
        logger.debug ("compare datasets: binned vs unbinned") 
        return False
    if b2 and not b1 :
        logger.debug ("compare datasets: unbinned vs binned") 
        return False

    ## both weighted or non-weighted?
    w1 = ds1.isWeighted ()
    w2 = ds2.isWeighted ()
    if w1 and not w2 :
        logger.debug ("compare datasets: weighted vs non-weighted") 
        return False
    if w2 and not w1 :
        logger.debug ("compare datasets: non-weighted vs weighted") 
        return False

    ## both non-poissoned 
    n1 = ds1.isNonPoissonWeighted ()
    n2 = ds2.isNonPoissonWeighted ()
    if n1 and not n2 :
        logger.debug ("compare datasets: NP-weighted vs non-weighted") 
        return False
    if n2 and not n1 :
        logger.debug ("compare datasets: non-weighted vs NP-weighted") 
        return False

    if w1 and w2 :
        
        if ds1.store_error      () and not ds2.store_error      () :
            logger.debug ("compare datasets: weight-errors vs no errors") 
            return False
        if ds2.store_error      () and not ds1.store_error      () :
            logger.debug ("compare datasets: no errors vs weight-errors") 
            return False
        if ds1.store_asym_error () and not ds2.store_asym_error () :
            logger.debug ("compare datasets: asym-errors vs no errors") 
            return False
        if ds2.store_asym_error () and not ds1.store_asym_error () :
            logger.debug ("compare datasets: no errors vs asym-errors") 
            return False

        wv1 = Ostap.Utils.getWeight ( ds1 )
        wv2 = Ostap.Utils.getWeight ( ds2 )

        if wv1 != wv2 :
            logger.debug ("compare datasets: different weight names") 
            return False
        
    ## check content
    if 0 == len ( ds1 ) : return True 
    
    st1   = ds1.statVars ( list ( vars1 ) )
    st2   = ds2.statVars ( list ( vars2 ) )
    keys1 = set ( st1.keys() )
    keys2 = set ( st2.keys() )
    if keys1 != keys2 :
        logger.debug ("compare datasets: different vars") 
        return False

    from ostap.math.base import isequal, isequalf  
    for k in keys1 :

        s1 = st1 [ k ]
        s2 = st2 [ k ]
        if s1 != s2 :
            logger.debug ("compare datasets: different %s" % k ) 
            return False 

    return True

# ============================================================================
## Are two datasets non-equal by content?
#  @code
#  ds1 = ...
#  ds2 = ...
#  ds_nonequal ( ds1 , ds2 )
#  ds1 != ds2 
#  @endcode 
def ds_nonequal ( ds1 , ds2 ) : 
    """ Are two datastes non-equal by content?
    ds1 = ...
    ds2 = ...
    ds_nonequal ( ds1 , ds2 ) 
    ds1 != ds2 
    """
    return not ds_equal ( ds1 , ds2 ) 

# =============================================================================
## Are two datasets equal by content?
def _ds_eq_ ( ds1 , ds2 ) :
    """ Are two datasets equal by content?"""
    if not isinstance  ( ds2 , ROOT.RooAbsData ) : return NotImplemented 
    return ds_equal ( ds1 , ds2 ) 

# =============================================================================
## Are two datasets non-equal by content?
def _ds_ne_ ( ds1 , ds2 ) :
    """ Are two datasets non-equal by content?"""
    if not isinstance  ( ds2 , ROOT.RooAbsData ) : return NotImplemented 
    return ds_nonequal ( ds1 , ds2 ) 

# ==========================================================================
ROOT.RooAbsData.__eq__ = _ds_eq_ 
ROOT.RooAbsData.__ne__ = _ds_ne_ 
_new_methods_ += [
    ROOT.RooAbsData.__eq__ , 
    ROOT.RooAbsData.__ne__ 
    ]

# ============================================================================
try : # ======================================================================
    # ========================================================================
    from numpy import array as _array 
    def get_result ( data ) : return _array (  data , dtype = float )
    # =======================================================================
except ImportError : # ======================================================
    # =======================================================================
    from array import array as _array 
    def get_result ( data ) : return _array ( 'd' , data )
    # = =====================================================================
# ===========================================================================
## Iterator for rows in dataset
#  @code
#  dataset = ...
#  for index, row , weight in dataset.rows ( 'pt, pt/p, mass ' , 'pt>1' ) :
#     print (index, row, weight) 
#  @endcode 
def _rad_rows_ ( dataset                  ,
                 variables  = []          ,
                 cuts       = ''          ,
                 cut_range  = ''          ,
                 first      = FIRST_ENTRY ,
                 last       = LAST_ENTRY  ,
                 progress   = False       ) :
    """ Iterator for rows in dataset
    >>> dataset = ...
    >>> for index, row , weight in dataset.rows ( 'pt, pt/p, mass ' , 'pt>1' ) :
    >>>    print (index, row, weight) 
    """

    first  , last    = evt_range      ( dataset   , first , last ) 
    varlst , cuts, _ = vars_and_cuts  ( variables , cuts )
    
    formulas = []
    varlist  = dataset.varlist () 
    for v in varlst :
        f0 = make_formula     ( v , v , varlist  ) 
        formulas.append ( f0 ) 

    ## loop over dataset
    for index, entry, weight in _rad_loop_ ( dataset               ,
                                             cuts      = cuts      ,
                                             cut_range = cut_range ,
                                             first     = first     ,
                                             last      = last      ,
                                             progress  = progress  ) :
        
        result = tuple ( float ( f ) for f in formulas )
        yield index , get_result ( result ) , weight 

    del formulas
    
ROOT.RooAbsData.rows = _rad_rows_         
    
_new_methods_ += [
    ROOT.RooAbsData.rows 
    ]

# ============================================================================
## Convert RooDataSet to TTree
#  @code
#  dataset = ...
#  data    = dataset.ds2tree ( name = 'my_tree' , filename= 'aa.root'
#  @endcode
#  - result of the type <code>ostap.tree.data_utils.Data</code>
#  @param name tree name (if not specified,  dataset name will be used)
#  @param filename tile name (if not specified, temporary file will be used) 
def _ds_2tree_ ( dataset , name = '' , filename = '' , cuts = '' , vars = () , cut_range = '' ) :
    """ Convert `ROOT.RooDataSet` to `ROOT.TTree`
    >>> dataset = ...
    >>> data    = dataset.ds2tree ( name = 'my_tree' , filename= 'aa.root'
    - name     : tree name (if not specified, dataset name will be used)
    - filename : file name (if not specified, temporary file will be used) 
    - result of the type `ostap.tree.data_utils.Data`
    """
    
    if not name : name = dataset.GetName()  
    if not filename :
        import ostap.utils.cleanup as CU 
        filename = CU.CleanUp.tempfile ( suffix = '.root' )

    cuts      = str ( cuts      ).strip() if cuts      else '' 
    cut_range = str ( cut_range ).strip() if cut_range else '' 
    if cuts or cut_range or vars : 
        dsaux = dataset.subset ( variables = vars      ,
                                 cuts      = cuts      ,
                                 cut_range = cut_range )
        result = _ds_2tree_ ( dsaux , name = name , filename = filename )
        if dsaux and isinstance ( dsaux , ROOT.RooDataSet ) :
            dsaux = Ostap.MoreRooFit.delete_data ( dsaux)
            del dsaux
        return result 
                    
    import ostap.io.root_file
    from   ostap.trees.data_utils import Data

    store = dataset.store()
    dstmp = None 
    with ROOTCWD() :  
        with ROOT.TFile ( filename , 'c' ) as rfile :
            rfile.cd() 
            if not isinstance ( store , ROOT.RooTreeDataStore ) :
                with useStorage ( ROOT.RooAbsData.Tree ) : 
                    dstmp = ROOT.RooDataSet ( dataset , dsID() )
                    store = dstmp.store()
                    if not isinstance ( store , ROOT.RooTreeDataStore ) :                        
                        dstmp.convertToTreeStore()
                        store = dstmp.store() 
            assert isinstance ( store , ROOT.RooTreeDataStore ) , \
                'Store type %s is not RooTreeDataStore!' % ( typename ( store ) ) 
            rfile [ name ] = store.tree()

    with ROOT.TFile ( filename , 'r' ) as rfile : rfile.ls()
            
    if dstmp and isinstance ( dstmp , ROOT.RooDataSet ) :
        dstmp = Ostap.MoreRooFit.delete_data ( dstmp )
        del dstmp
        
    return Data ( chain       = name         ,
                  files       = [ filename ] ,
                  description = "TTree from dataset  %s/%s " % ( dataset.name , dataset.title ) ,
                  silent      = True        ) 

ROOT.RooDataSet.ds2tree = _ds_2tree_

_new_methods_ += [
    ROOT.RooDataSet.ds2tree 
    ]

# ============================================================================

_new_methods_ += list ( data_decorate ( ROOT.RooAbsData ) )
del data_decorate

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
##                                                                      The END 
# =============================================================================
