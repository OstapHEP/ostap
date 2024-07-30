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
    'ds_draw'    , ## draw variables from RooDataSet 
    'ds_project' , ## project variables from RooDataSet to histogram
    'ds_combine' , ## combine two datasets with weights 
    )
# =============================================================================
from   builtins                  import range
from   collections               import defaultdict
from   ostap.core.meta_info      import root_info, ostap_version 
from   ostap.core.core           import ( Ostap, VE, SE ,
                                          hID  , dsID , strings , 
                                          valid_pointer , split_string   ,
                                          ROOTCWD       , var_separators ,
                                          loop_items    , cidict_fun     )
from   ostap.core.ostap_types    import ( integer_types , string_types   ,
                                          num_types     , dictlike_types , 
                                          list_types    , sequence_types )
from   ostap.math.base           import islong
from   ostap.fitting.variables   import valid_formula, make_formula 
from   ostap.trees.cuts          import expression_types, vars_and_cuts
from   ostap.utils.utils         import evt_range, LAST_ENTRY

from   ostap.stats.statvars      import data_decorate, data_range 

import ostap.fitting.roocollections
import ostap.fitting.printable
import ROOT, random, math, sys, ctypes  
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
#  event   = dataset[4]            ## index 
#  events  = dataset[0:1000]       ## slice 
#  events  = dataset[0:-1:10]      ## slice 
#  events  = dataset[ (1,2,3,10) ] ## seqeucne of indices  
#  @eendcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-03-31
def _rad_getitem_ ( data , index ) :
    """Get the entry from RooDataSet
    >>> dataset = ...
    >>> event  = dataset[4]                 ## index 
    >>> events = dataset[0:1000]            ## slice
    >>> events = dataset[0:-1:10]           ## slice 
    >>> events = dataset[ (1,2,3,4,10) ]    ## sequnce of indices 
    """

    N = len ( data )
    
    if isinstance ( index , integer_types ) and index < 0 :
        index += N

    if isinstance ( index , integer_types ) and 0 <= index  < N :
        
        return data.get ( index )  ## should we add weight here? 
   
    elif isinstance ( index , range ) :

        ## simpel case 
        start , stop , step = index.start , index.stop , index.step
        if 1 == step : return data.reduce ( ROOT.RooFit.EventRange ( start , stop ) )
                
    elif isinstance ( index , slice ) :
        
        start , stop , step = index.indices ( N )
                              
        if 1 == step : return data.reduce ( ROOT.RooFit.EventRange ( start , stop ) )

        index = range ( start , stop , step ) 


    ## the actual loop over entries 
    if isinstance ( index , sequence_types ) :

        weighted = data.isWeighted                    ()
        se       = weighted and data.store_error      ()
        sae      = weighted and data.store_asym_error ()

        result = data.emptyClone ( dsID () )
        for i in index :

            j = int ( i )                 ## the content must be convertible to integers 

            if j < 0 : j += N             ## allow negative indices 
            
            if not 0 <= j < N :           ## is adjusted integer in the proper range ? 
                logger.error ( 'Invalid index [%s]=%s, skip it' % ( i , j ) ) 
                continue  
            
            vars = data.get ( j )
            
            if   weighted and sae :
                wel , weh = data.weight_errors()
                result.add ( vars , data.weight () , wel , weh ) 
            elif weighted and se  :
                we        = data.weightError()
                result.add ( vars , data.weight () , we  ) 
            elif weighted         :
                result.add ( vars , data.weight () ) 
            else : 
                result.add ( vars ) 
            
        return result

    raise IndexError ( 'Invalid index %s'% index )

# ==============================================================================
## Get (asymmetric) weigth errors for the current entry in dataset
#  @code
#  dataset = ...
#  weight_error_low, weight_error_high = dataset.weightErrors() 
#  @endcode
#  @see RooAbsData::weightError
def _rad_weight_errors( data , *etype ) :
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
    if isinstance ( aname , integer_types ) :
        return 0<= aname < len ( self )
    
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
## Make dataset with removed i-th element  (for Jackknife/bootstrapping)
#  @code
#  dataset = ...
#  N = len ( dataset)
#  for index in range ( N ) :
#    ds_i = dataset - index
#  @endcode 
def _rds_sub_ ( dataset , index ) :
    """Make dataset with removed i-th element  (for Jackknife/bootstrapping)
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

        ds1.clear()
        ds2.clear()

        assert len ( result ) + 1 == N , 'Invalid length of the resulting dataset!'
        
        return result
    
    return NotImplemented
        
# ============================================================================
## Jackknife generator: generates data sets with removed i-th element
#  @code
#  dataset = ...
#  for ds in ds.jackknife() :
#  ...
#  @endcode 
def _rds_jackknife_ ( dataset , low = 0 , high = None ) :
    """Jacknife generator
    >>> dataset = ...
    >>> for ds in ds.jackknife() :
    >>> ...
    """
    N = len ( dataset )
    
    if high == None : high = N
    
    if 1 < N : 
        for i in range ( low , high ) :
            ds_i = dataset - i        ## this is the line 
            yield ds_i               
            ds_i.clear()             
            del ds_i

# =============================================================================
## Boostrap generator
#  @code
#  dataset = ...
#  for ds in dataset.bootstrap ( 100 ) :
#  ...
#  @endcode
def _rds_bootstrap_ ( dataset , size = 100 , extended = False ) :
    """ Boostrap generator
    >>> dataset = ...
    >>> for ds in dataset.bootstrap ( 100 ) :
    >>> ...
    """

    from   ostap.stats.bootstrap  import bootstrap_indices, extended_bootstrap_indices 

    N    = len ( dataset )
    bgen = bootstrap_indices ( N , size = size ) if not extended else extended_bootstrap_indices ( N , size = size ) 
    
    for indices in bgen : 
        ds = dataset [ indices ] 
        yield ds
        ds.clear()
        del ds

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

    if ( not vars ) and ( not cuts ) : return dset.__class__ ( dset )  ##  return copy 
    elif not vars : vars = dset.varset()
    
    if   isinstance ( cuts , ROOT.TCut ) :   cuts = str ( cuts )

    if cuts and isinstance ( cuts , string_types ) :
        cuts  = make_formula ( cuts , cuts , dset.varlist() )

    if   isinstance ( vars , ROOT.RooArgSet  ) : return dset.reduce ( vars , cuts )
    elif isinstance ( vars , string_types    ) : vars = [ vars ]

    
    aset   = ROOT.RooArgSet()
    varset = dset.get()
    for v in vars :
        if   isinstance ( v , ROOT.RooAbsArg ) : aset.add ( v )
        elif isinstance ( v , string_types   ) and v in varset : aset.add ( varset[v] )
        else :
            raise TypeError("``subset'': unknown type %s/%s" % ( v , type(v) ) )

    return dset.reduce ( aset , cuts ) if aset else dset.reduce ( cuts )


# =============================================================================
## helper technical method to seek for unique or duplicated entries
def _rds_seek_for_duplicates_ ( dataset          , 
                                entrytag         , 
                                criterium = ''   ) : 
    """Helper technical method to seek for unique or duplicated entries
    """
    
    if   isinstance ( entrytag , ROOT.RooAbsArg ) : entrytag = [ entrytag ]
    elif isinstance ( entrytag , string_types   ) : entrytag = [ entrytag ]
    
    assert isinstance ( entrytag , sequence_types ) , 'Invalid "enttytag" %s' % str ( entrytag )

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
    for i, e in enumerate ( dataset ) :

        entry = tuple (  float ( v ) for v in tag ) 
        snapshot [ entry ].append ( content ( i )  ) 

    return snapshot

# =============================================================================
## Iterator over the duplicated groups 
#  @code
#  for group in dataset.duplicates ( ( 'evt' , 'run' ) ) :
#  ... for entry in group :
#  ... ...
def _rds_duplicates_ ( dataset  ,
                       entrytag ) :
    """Iterator over the duplicated groups 
    >>> for group in dataset.duplicats ( ( 'evt' , 'run' ) ) :
    >>> ... for entry in group :
    >>> ... ...
    """
    snapshot = _rds_seek_for_duplicates_ ( dataset              ,
                                           entrytag  = entrytag ,
                                           criterium = ''       ) 
    for e, lst in loop_items ( snapshot ) :
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
def _rds_unique_entries_ ( dataset          ,
                           entrytag         ,
                           choice           , 
                           criterium = ''   , 
                           seed      = None ,
                           report    = True ) :
    
    if criterium  :
        assert choice in ( 'min' , 'max' , 'minimal' , 'maximal' , 'minimum' , 'maximum' ) , \
               "Invalid 'choice' for criterium!"
    else : 
        assert choice in ( 'first' , 'last' , 'random' , 'rndm'  , 'rand' ) , \
               "Invalid 'choice'"
 
    snapshot = _rds_seek_for_duplicates_ ( dataset               ,
                                           entrytag  = entrytag  ,
                                           criterium = criterium ) 

    choice = choice.lower()
    
    first  = 'first'  == choice
    last   = 'last'   == choice
    rand   = choice in ( 'random' , 'rndm'  , 'rand' ) 
    minv   = criterium and choice in ( 'min' , 'minimal' , 'minimum' )
    maxv   = criterium and choice in ( 'max' , 'maximal' , 'maximum' )

    from ostap.utils.utils import random_seed

    with random_seed ( seed ) :

        cnt    = SE ()
        unique = 0 
        for e, lst in loop_items ( snapshot ):

            unique += 1
            
            num  = len ( lst )
            if 1 < num : cnt += num 
            
            if   1 == num  : yield lst [   0 ] [ 1 ]        
            elif first     : yield lst [   0 ] [ 1 ]
            elif last      : yield lst [  -1 ] [ 1 ]
            elif rand      : yield random.choice ( lst ) [ 1 ] ## seed is needed here! 
            elif minv      : yield min ( lst ) [ 1 ]
            elif maxv      : yield max ( lst ) [ 1 ]
                
    if report :
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
        logger.info ( '%s:\n%s' % ( title , T.table ( rows , title = title , prefix = '# ' , alignment = 'lc' ) ) )
        
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
def _rds_make_unique_ ( dataset          ,
                        entrytag         ,
                        choice           , 
                        criterium = ''   , 
                        seed      = None ,
                        report    = True ) :
    """Make a copy of dataset only with unique  entries 
    >>> dataset = ...
    >>> unique = dataset.make_unique ( ( 'evt' , 'run' ) , choice = 'random' )
    >>> unique = dataset.make_unique ( ( 'evt' , 'run' ) , choice = 'first'  )
    >>> unique = dataset.make_unique ( ( 'evt' , 'run' ) , choice = 'last'    )
    >>> unique = dataset.make_unique ( ( 'evt' , 'run' ) , choice = 'min' , criterium = 'PT' )
    >>> unique = dataset.make_unique ( ( 'evt' , 'run' ) , choice = 'max' , criterium = 'PT' )
    - CPU performance is more or less reasonable up to dataset with 10^7 entries 
    """
    
    ds = dataset.emptyClone()
    for i in dataset.unique_entries ( entrytag  = entrytag  ,
                                      choice    = choice    ,
                                      criterium = criterium ,
                                      seed      = seed      ,
                                      report    = report    ) :
        ds.add ( dataset [ i ] )

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

ROOT.RooAbsData . __add__       = _rad_add_
ROOT.RooDataSet . __iadd__      = _rad_iadd_

ROOT.RooAbsData . __mul__       = _rad_mul_
ROOT.RooAbsData . __rmul__      = _rad_mul_
ROOT.RooAbsData . __imul__      = _rad_imul_
ROOT.RooAbsData . __mod__       = _rad_mod_
ROOT.RooAbsData . __div__       = _rad_div_
ROOT.RooAbsData . __truediv__   = ROOT.RooAbsData . __div__


ROOT.RooAbsData . sample        = _rad_sample_
ROOT.RooAbsData . shuffle       = _rad_shuffle_



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
   ROOT.RooAbsData . __mul__       ,
   ROOT.RooAbsData . __rmul__      ,
   ROOT.RooAbsData . __imul__      ,
   ROOT.RooAbsData . __div__       ,
   ROOT.RooAbsData . __mod__       ,
   ROOT.RooAbsData . __truediv__   ,
   #
   ROOT.RooAbsData . sample        ,
   ROOT.RooAbsData . shuffle       ,
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
def ds_project  ( dataset                ,
                  histo                   ,
                  what                    ,
                  cuts      = ''          ,
                  cut_range = ''          , 
                  first     =  0          ,
                  last      = LAST_ENTRY  ,
                  progress  = False       ) :
    """Helper project method for RooDataSet/DataFrame/... and similar objects 

    >>> h1   = ROOT.TH1D(... )
    >>> dataset.project ( h1.GetName() , 'm', 'chi2<10' ) ## project variable into histo
    
    >>> h1   = ROOT.TH1D(... )
    >>> dataset.project ( h1           , 'm', 'chi2<10' ) ## use histo
    """

    ## 1) redirect to approproate method
    ## ROOT.Tree 
    if isinstance ( dataset , ROOT.TTree ) :
        assert not cut_range, "ds_project(tree,...) cannot be used with cut_range %s" % cut_range 
        from ostap.trees.trees import tree_project
        return tree_project ( datatet , histo ,
                              what    = what   ,
                              cuts    = cuts   ,
                              first   = first  ,
                              last    = last   )
    target = histo

    ## 2) if the histogram is specified by the name, try to locate it ROOT memory  
    if isinstance ( target , string_types ) :
        hname  = str ( target ) 
        groot  = ROOT.ROOT.GetROOT()
        h      = groot.FindObject ( hname )
        assert h , 'Cannot get locate histo by name %s' % histos        
        assert isinstance ( h , ROOT.TH1 ) , 'the object %s i snot ROOT.TH1' % type ( h ) 
        target = h
        

    ## 3) input histogram?
    input_histo = isinstance ( target , ROOT.TH1 )

    ## 4) vali dtarget? 
    from ostap.trees.param import param_types_nD    
    assert input_histo or isinstance ( target , param_types_nD ) , 'Invalid target/histo type %s' % type ( target ) 

    if   isinstance ( histo , ROOT.TProfile2D ) :
        histo.Reset() 
        logger.error ('ds_project: TProfile2D is not (yet) supported')
        return histo 
    elif isinstance ( histo , ROOT.TProfile  ) :
        logger.error ('ds_project: TProfile   is not (yet) supported')
        return histo

    ## copy/clone the target 
    def target_copy  ( t ) : return t.Clone() if isinstance ( t , ROOT.TH1 ) else type ( t ) ( t )
    ## reset the target 
    def target_reset ( t ) :
        if isinstance ( t , ROOT.TH1 ) : t.Reset()
        else                           : t *= 0.0
        return t
    
    ## 5) reset the target 
    target = target_reset ( target ) 

    
    ## (5) dimension of the target 
    dim = target.dim ()
    assert 1 <= dim <= 4 , 'Invalid dimension of target: %s' % dim  

    ## (6) parse input expressions
    varlst, cuts, input_string = vars_and_cuts  ( what , cuts )
    if input_string and 2 <= len ( varlst ) and ostap_info < (1,11) : 
        if 0 < _to_print  : 
            logger.attention ("From ostap v1.10.1.9 variables are treted in natural order (no reverse!)")
            _to_print -= 1 

    print ( 'VARIABLES/1', varlst, cuts ) 
    nvars = len ( varlst ) 
    assert ( 1 == dim and dim <= nvars ) or dim == nvars , \
        'Mismatch between the target/histo dimension %d and input variables %s' % ( dim , varlst )

    
    print ( 'PROJECT-7' )


    ## RooAbsData 
    if isinstance ( dataset , ROOT.RooAbsData ) :
        
        ## 3) adjust the first/last
        first , last = evt_range ( len ( dataset ) , first , last )
        
        if not first < last : return target                             ## RETURN
        
        tail = cuts , cut_range , first , last

    ## soemthinug convertibel to DataFrame/Frame/Node
    else :
        
        assert not cut_range , 'cut-range is not allowed!'
        if  ( firts, last ) != ( 0 , LAST_RNTRY ) : 
            logger.warning  ( "ds_project(%s): ignored  `first'/`last' %s/%s " % ( type ( dataset ) , first , last ) ) 
        from ostap.frames.frames import as_rnode
        frame   = as_rnode ( dataset )
        dataset = frame 
        tail    = ()
        
    from ostap.utils.progress_conf import progress_conf
    
    the_args = ( target , ) + varlst + tail 
    print ( 'DS_PROJECT' , dim , nvars, varlst , the_args ) 
    ## very special case of projection several expressions into the same target 
    if 1 == dim and dim < nvars : 
        print ( 'DS_PROJECT/1' , dim , nvars, varlst , tail ) 
        ## very special case of projection several expressions into the same target 
        htmp  = target_copy ( target )  ## prepare temporary object 
        for var in varlst :
            htmp = target_reset ( htmp ) ## rest the temporary object 
            args = ( htmp , var ) + tail 
            if progress : sc = Ostap.HistoProject.project  ( dataset , progress_conf () , *args )
            else        : sc = Ostap.HistoProject.project  ( dataset ,                    *args )
            if not sc.isSuccess() : logger.error ( "Error from Ostap.HistoProject.project %s" % sc )
            ## update results 
            target += htmp
            del htmp 
    elif 1 == dim :
        print ( 'DS_PROJECT/2' , dim , nvars, varlst , tail ) 
        if progress : sc = Ostap.HistoProject.project  ( dataset , progress_conf () , *the_args )
        else        : sc = Ostap.HistoProject.project  ( dataset ,                    *the_args )
        if not sc.isSuccess() : logger.error ( "Error from Ostap.HistoProject.project  %s" % sc )
    elif 2 == dim : 
        print ( 'DS_PROJECT/3' , dim , nvars, varlst , tail ) 
        if progress : sc = Ostap.HistoProject.project2 ( dataset , progress_conf () , *the_args )
        else        : sc = Ostap.HistoProject.project2 ( dataset ,                    *the_args )
        if not sc.isSuccess() : logger.error ( "Error from Ostap.HistoProject.project2 %s" % sc )
        print ( 'DS_PROJECT/4' , dim , nvars, varlst , tail ) 
    elif 3 == dim : 
        if progress : sc = Ostap.HistoProject.project3 ( dataset , progress_conf () , *the_args )
        else        : sc = Ostap.HistoProject.project3 ( dataset ,                    *the_args )
        if not sc.isSuccess() : logger.error ( "Error from Ostap.HistoProject.project2 %s" % sc )
    elif 4 == dim : 
        if progress : sc = Ostap.HistoProject.project4 ( dataset , progress_conf () , *the_args )
        else        : sc = Ostap.HistoProject.project4 ( dataset ,                    *the_args )
        if not sc.isSuccess() : logger.error ( "Error from Ostap.HistoProject.project2 %s" % sc )
        
    return target if sc.isSuccess() else None 

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
              what                    ,
              cuts      = ''          ,
              opts      = ''          ,
              cut_range = ''          ,
              first     =  0          ,
              last      =  LAST_ENTRY ,
              delta     = 0.05        , **kwargs ) :
    """Helper draw method for drawing of RooDataSet
    >>> dataset.draw ( 'm', 'chi2<10'                 )
    ## cuts & weight 
    >>> dataset.draw ( 'm', '(chi2<10)*weight'        )
    ## use drawing options 
    >>> dataset.draw ( 'm', '(chi2<10)*weight' , 'e1' )
    ## start form event #1000     >>> dataset.draw ( 'm', '(chi2<10)*weight' , 'e1' , 1000 ) 
    ## for event in range 1000< i <10000
    >>> dataset.draw ( 'm', '(chi2<10)*weight' , 'e1' , 1000 , 100000 )
    """

    ## decode variables/cuts 
    varlst, cuts, input_string = vars_and_cuts  ( what , cuts )
    
    nvars = len ( varlst ) 
    assert 1 <= nvars <= 3 , "Invalid number of variables: %s" % str ( varlst )
    
    ## get the suitable ranges for the variables 
    ranges = ds_range ( dataset , varlst , cuts = cuts , cut_range = cut_range , first = first , last  = last , delta = delta )
    for _, r in loop_items ( ranges ) :
        mn , mx = r
        ## no useful entries 
        if mx < mn :
            logger.warning ("No entries, no draw, return None" ) 
            return None 
        
    assert len ( ranges ) == nvars , 'Invalid ranges: %s' % str ( ranges )

    print ('RANGES' , ranges )
    
    
    from ostap.utils.cidict        import cidict
    kw = cidict ( transform = cidict_fun , **kwargs )
    
    if 1 == nvars :
        xvar       = varlst [ 0    ] 
        xmin, xmax = ranges [ xvar ]
        #
        xmin       = kw.pop ( 'xmin' , xmin )
        xmax       = kw.pop ( 'xmax' , xmax )
        assert xmin < xmax , "Invalid xmin/xmax setting!"
        #
        xbins      = kw.pop ( 'xbins' , kw.pop ( 'bins' , kw.pop ( 'nbinsx' , kw.pop ( 'nbins' , kw.pop ( 'binsx' , 100 ) ) ) ) )
        #
        # book the histogram 
        histo      = ROOT.TH1F  ( hID()  , "%s" % ( xvar ) , xbins  , xmin , xmax )  ; histo.Sumw2()
        # 
        ds_project ( dataset , histo , varlst , cuts = cuts , cut_range = cut_range , first = first , last = last )
        histo.draw ( opts , **kw )
        return histo
    
    elif 2 == nvars :
        
        xvar       = varlst [ 0 ] 
        yvar       = varlst [ 1 ] 
        xmin, xmax = ranges [ xvar ] 
        ymin, ymax = ranges [ yvar ]
        #
        xmin       = kw.pop ( 'xmin' , xmin )
        xmax       = kw.pop ( 'xmax' , xmax )
        assert xmin < xmax , "Invalid xmin/xmax setting!"
        #
        ymin       = kw.pop ( 'ymin' , ymin )
        ymax       = kw.pop ( 'ymax' , ymax )
        assert ymin < xmax , "Invalid ymin/ymax setting!"
        #
        xbins      = kw.pop ( 'xbins' , kw.pop ( 'nbinsx' ,kw.pop ( 'binsx' , 50 ) ) ) 
        ybins      = kw.pop ( 'ybins' , kw.pop ( 'nbinsy' ,kw.pop ( 'binsy' , 50 ) ) )
        #
        histo      = ROOT.TH12  ( hID()  , "x = %s , y = %s" % ( xvar , yvar ) ,
                                  xbins , xmin , xmax , 
                                  ybins , ymin , ymax ) ; histo.Sumw2() 
        ds_project ( dataset , histo , varlst , cuts = cuts , cut_range = cut_range , first = first , last = last )
        histo.draw ( opts , **kw  )
        return histo

    elif 3 == nvars :
        
        xvar       = varlst [ 0 ] 
        yvar       = varlst [ 1 ] 
        zvar       = varlst [ 2 ]
        #
        xmin, xmax = ranges [ xvar ]
        ymin, ymax = ranges [ yvar ]
        zmin, zmax = ranges [ zvar ]
        #
        xmin       = kw.pop ( 'xmin' , xmin )
        xmax       = kw.pop ( 'xmax' , xmax )
        assert xmin < xmax , "Invalid xmin/xmax setting!"
        #
        ymin       = kw.pop ( 'ymin' , ymin )
        ymax       = kw.pop ( 'ymax' , ymax )
        assert ymin < xmax , "Invalid ymin/ymax setting!"
        #
        zmin       = kw.pop ( 'zmin' , zmin )
        zmax       = kw.pop ( 'zmax' , zmax )
        assert zmin < zmax , "Invalid zmin/zmax setting!"
        #
        #
        xbins      = kw.pop ( 'xbins' , kw.pop ( 'nbinsx' , kw.pop ( 'binsx' , 20 ) ) ) 
        ybins      = kw.pop ( 'ybins' , kw.pop ( 'nbinsy' , kw.pop ( 'binsy' , 20 ) ) ) 
        zbins      = kw.pop ( 'zbins' , kw.pop ( 'nbinsz' , kw.pop ( 'binsz' , 20 ) ) )
        #
        histo      = ROOT.TH13  ( hID()  , "x = %s , y = %s , z = %s" % ( xvar , yvar , zvar ) ,
                                  xbins , xmin , xmax , 
                                  ybins , ymin , ymax ,
                                  zbins , zmin , zmax ) ; histo.Sumw2() 
        ds_project ( dataset , histo , varlst , cuts = cuts , cut_range = cut_range , first = first , last = last )
        histo.draw ( opts , **kw )
        return histo

# =============================================================================
## get the attibute for RooDataSet
def _ds_getattr_ ( dataset , aname ) :
    """Get the attibute from RooDataSet 

    >>> dset = ...
    >>> print dset.pt
    
    """
    _vars = dataset.get()
    return getattr ( _vars , aname )  
## get the attibute for RooDataSet
# =============================================================================
def get_var ( self, aname ) :
    _vars = self.get()
    return getattr ( _vars , aname )  



# =============================================================================
## Get min/max for the certain variable/expression in dataset
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

    assert isinstance ( var  , expression_types ) , 'Invalid expression!'
    assert isinstance ( cuts , expression_types ) , 'Invalid expression!'
    
    var  = str ( var  )
    cuts = str ( cuts ).strip()
    
    
    if isinstance ( var , ROOT.RooAbsReal ) : var = var.GetName() 
    if cuts : s = dataset.statVar ( var , cuts )
    else    : s = dataset.statVar ( var )
    mn , mx = s.minmax()

    ## if range is invalid, recalculate the range without cuts 
    if mx < mn and cuts :
        mn , mx = ds_var_minmax ( dataset , var )

    ## if range is invalid and variable is in dataset, try to use the native range 
    if mx < mn and var in dataset :
        vv = getattr ( dataset , var ) 
        mn, mx = vv.minmax ()
        
    if mn < mx and 0.0 < delta :
        dx   = delta * 1.0 * ( mx - mn )  
        mx  += dx   
        mn  -= dx
        
    return mn , mx


ROOT.RooDataSet .vminmax  = ds_var_minmax 

_new_methods_ += [
    ROOT.RooDataSet .vminmax ,
    ]




# =============================================================================
## Is there at least one entry that satisfy selection criteria?
#  @code
#  dataset = ...
#  dataset.hasEntry ( 'pt>100' )
#  dataset.hasEntry ( 'pt>100' , 0, 1000 ) ## 
#  dataset.hasEntry ( 'pt>100' , 'fit_range' ) ## 
#  dataset.hasEntry ( 'pt>100' , 'fit_range' , 0 , 1000 ) ## 
#  @endcode
#  @see Ostap::StatVar::hasEntry
def _ds_has_entry_ ( dataset , selection , *args ) : 
    """Is there at leats one entry that satoisfies selection criteria?
    >>> dataset = ...
    >>> dataset.hasEntry ( 'pt>100' )
    >>> dataset.hasEntry ( 'pt>100' , 0, 1000 ) ## 
    >>> dataset.hasEntry ( 'pt>100' , 'fit_range' ) ## 
    >>> dataset.hasEntry ( 'pt>100' , 'fit_range' , 0 , 1000 ) ## 
    - see Ostap.StatVar.hasEntry
    """
    
    assert isinstance ( selection , expression_types ) , 'Invalid expression!'
    selection  = str  ( selection ).strip() 

    result = Ostap.StatVar.hasEntry ( dataset , selection , *args ) 
    return True if result else False 

ROOT.RooAbsData.hasEntry  = _ds_has_entry_
ROOT.RooAbsData.has_entry = _ds_has_entry_

_new_methods_ += [
    ROOT.RooAbsData.hasEntry  , 
    ROOT.RooAbsData.has_entry , 
    ]

# =============================================================================
## find sutable range for drawing a variable
#  @code
#  data      = ...
#  min . max = ds_var_range ( data , 'variname' , cuts = ... ) 
#  @endcode
def ds_var_range ( dataset , var , cuts = '' ) :
    """Find suitable range for drawing a variable
    >>> data      = ...
    >>> min , max = ds_var_range ( data , 'variname' , cuts = ... ) 
    """
    ## min/max values
    mn , mx = ds_var_minmax ( dataset , var , cuts )
    return axis_range ( mn , mx , delta = 0.05 )


# =============================================================================\
## Get suitable ranges for drawing expressions/variables
## @code
#  dataset = ...
#  result  = ds_range ( dataset , 'sin(x)*100*y' , 'x<0' )
#  results = ds_range ( dataset , 'x,y,z,t,u,v'  , 'x<0' )
#  @endcode 
def ds_range  ( dataset                ,
                expressions            ,
                cuts      = ''         ,
                cut_range = ''         ,
                first     = 0          , 
                last      = LAST_ENTRY ,
                delta     = 0.05       ) :
    """ Get suitable ranges for drawing expressions/variables
    >>> dataset = ...
    >>> result  = ds_range ( dataset , 'sin(x)*100*y' , 'x<0' )
    >>> results = ds_range ( dataset , 'x,y,z,t,u,v'  , 'x<0' )
    """
    return data_range ( dataset     ,
                        expressions ,
                        cuts        ,
                        delta       ,
                        cut_range   ,
                        first       ,
                        last        )

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
        return len ( self )
    
    ## ROOT.RooDataSet.reset = _ds_new_reset_
    ROOT.RooDataSet.clear = _ds_new_reset_
    ROOT.RooDataSet.erase = _ds_new_reset_
    ## ROOT.RooDataSet.Reset = _ds_new_reset_

# ROOT.RooDataSet.clear = ROOT.RooDataSet.reset
# ROOT.RooDataSet.erase = ROOT.RooDataSet.reset
# ROOT.RooDataSet.Reset = ROOT.RooDataSet.reset

ROOT.RooDataSet.get_var       = get_var

_new_methods_ += [
    ROOT.RooDataSet .clear ,
    ROOT.RooDataSet .erase ,
    
    ROOT.RooDataSet .get_var ,
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
    """Get the s-factor for (weighted) dataset, where
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
    if sys.version_info < ( 3 , 0 ) :
        if isinstance ( result , unicode ) :
            result = result.encode ('utf-8')
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
    """Clone dataset
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
    """Clone dataset
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
    """Add/calculate variable to RooDataSet

    >>> dataset.addVar ( 'ratio' , 'pt/pz' )
    """
    vcol = make_formula    ( vname , formula , dataset.get() )
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
def add_new_var ( dataset , varname , what , *args ) : 
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
        
    vv = Ostap.Functions.add_var ( dataset , varname , what , *args )
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
    formula = formula or 0 <= wvarname.find ( '*' ) 
    formula = formula or 0 <= wvarname.find ( '/' )     
    formula = formula or 0 <= wvarname.find ( '+' ) 
    formula = formula or 0 <= wvarname.find ( '-' )     
    formula = formula or 0 <= wvarname.find ( '&' )     
    formula = formula or 0 <= wvarname.find ( '|' )     
    formula = formula or 0 <= wvarname.find ( '^' )     
    formula = formula or 0 <= wvarname.find ( '%' )
    
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
    """Context manager to change the storage type
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
def _ds_table_0_ ( dataset           ,
                   variables = []    ,
                   cuts      = ''    ,
                   first     = 0     ,
                   last      = 2**62 ,
                   prefix    = ''    ,
                   title     = ''    ) :
    """Print data set as table
    """
    varset = dataset.get()
    if not valid_pointer ( varset ) :
        logger.error('Invalid dataset')
        return ''
    
    if variables : vars , cuts , _ = vars_and_cuts ( variables , cuts ) 
    else         : vars , cuts     = [ v.name for v in varset ] , str ( cuts ).strip() 
    
    vars = [ i.GetName() for i in varset if i.GetName() in vars ]
        
    _vars = []
    
    stat = dataset.statVars ( vars , cuts , first , last ) 
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
            title = '%s("%s","%s"):' % ( dataset.__class__.__name__ , dataset.GetName () , tt ) 
        else :
            title = '%s("%s"):'      % ( dataset.__class__.__name__ , dataset.GetName () )

        title =  '%s %d entries, %d variables' %  ( title , len ( dataset ) , len ( varset ) )
        
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

                last  = min ( last , len ( dataset ) + 1 ) 
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


# =============================================================================
## Print data set as table
def _ds_table_1_ ( dataset           ,
                   variables = []    ,
                   cuts      = ''    ,
                   first     = 0     ,
                   last      = 2**62 ,
                   prefix    = ''    ,
                   title     = ''    ) :
    """Print data set as table
    """

    
    if variables : vars , cuts , _ = vars_and_cuts ( variables , cuts ) 
    else         : vars , cuts     = [ v.name for v in varset ] , str ( cuts ).strip() 
    
    _vars = []    
    vvars = tuple ( sorted ( vars ) )
    stat = dataset.statVars ( vvars  , cuts , first , last ) 
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
        
        if  tt and tt != dataset.GetName()  : 
            title = '%s("%s","%s"):' % ( dataset.__class__.__name__ , dataset.GetName () , tt ) 
        else :
            title = '%s("%s"):'      % ( dataset.__class__.__name__ , dataset.GetName () )

        title =  '%s %d entries, %d variables' %  ( title , len ( dataset ) , len ( varset ) )
        
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
def _ds_table_ (  dataset ,  variables = [] , cuts = '' , prefix = '' , title = '' ) :
    """print dataset in a form of the table
    >>> dataset = ...
    >>> print dataset.table()
    """
    return _ds_table_0_ ( dataset ,  variables , cuts = cuts , prefix = prefix , title = title )[0]

# ==============================================================================
## print dataset in a form of the table
#  @code
#  dataset = ...
#  print dataset.table2() 
#  @endcode
def _ds_table2_ (  dataset ,  variables , cuts = '' , prefix = '' , title = '' ) :
    """print dataset in a form of the table
    >>> dataset = ...
    >>> print dataset.table()
    """
    return _ds_table_1_ ( dataset ,  variables , cuts = cuts , prefix = prefix , title = title )[0]

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

    names   = tuple ( v.name for v in nvars )

    nds     = ds.emptyClone ()
    nvarset = nds.varset    ()
    
    for v in nvarset :
        if v.name in names : 
            if hasattr ( v ,  'setMin' ) : v.setMin ( mnv )
            if hasattr ( v ,  'setMax' ) : v.setMax ( mxv )        
    
    ## loop over the data set 
    for entry in ds :

        values = [ v.getVal() for v in entry if v.name in names ]        
        random.shuffle ( values )
        
        for v in nvarset :
            n = v.name 
            if not n in names : v.setVal ( entry [ n ] .value )                
            else              : v.setVal ( values.pop()   )

        nds.add ( nvarset )
        
    return nds


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
    """Get the name of weigth variable in dataset
    >>> dataset = ...
    >>> wname   = dataset.wname() 
    """
    
    if not dataset.isWeighted() : return '' ## UNWEIGHTED!

    attr = '_weight_var_name'
    if not hasattr ( dataset , attr ) :
        wn = Ostap.Utils.getWeight ( dataset )    
        setattr ( dataset , attr , wn ) 
        
    return getattr ( dataset , attr , '' )
# =============================================================================


# =============================================================================
## Get the weight variable from the dataset
#  @code
#  dataset = ...
#  wvar    = dataset.weightVar() 
#  @endcode 
def _ds_weight_var_ ( dataset ) :
    """Get the weight variable from the dataset
    
    >>> dataset = ...
    >>> wvar    = dataset.weightVar()
    
    """
    wvar = Ostap.Utils.getWeightVar ( dataset )
    return wvar 
    ## return wvar if valid_pointer ( wvar ) else None

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
    """Are weight errors stored in dataset?
    >>> dataset     = ...
    >>> store_error = dataset.store_error () 
    The function checks the `StoreError` and 
    `StoreAsymError` attributes for the weight variable 
    - see Ostap::Utils::storeError
    """
    
    if not dataset.isWeighted() : return False ## UNWEIGHTED!
    
    attr = '_store_weight_error'
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
    """Are weight errors stored in dataset?
    >>> dataset     = ...
    >>> store_error = dataset.store_asym_error () 
    The function checks the `StoreAsymError` attributes for the weight variable 
    - see Ostap::Utils::storeAsymError
    """
    
    if not dataset.isWeighted() : return False ## UNWEIGHTED!
    
    attr = '_store_asym_weight_error'
    if not hasattr ( dataset , attr ) :
        
        wn = Ostap.Utils.storeAsymError  (  dataset )
        wn = True if wn else False
        
        setattr ( dataset , attr , wn ) 
        
    return getattr ( dataset , attr , '' )

# =============================================================================

ROOT.RooDataSet.wname            = _ds_wname_
ROOT.RooDataSet.weight_name      = _ds_wname_
ROOT.RooDataSet.store_error      = _ds_store_error_
ROOT.RooDataSet.store_asym_error = _ds_store_asym_error_

if not hasattr ( ROOT.RooDataSet , 'weightVar' ) :
    ROOT.RooDataSet.weightVar = _ds_weight_var_ 


_new_methods_ += [
    ROOT.RooDataSet.wname            , 
    ROOT.RooDataSet.weight_name      , 
    ROOT.RooDataSet.store_error      , 
    ROOT.RooDataSet.store_asym_error ,
    ]

if (3,0) <= sys.version_info :
    def f_open ( name , mode , **kwargs ) :
        return open ( name , mode , **kwargs )
else : 
    def f_open ( name , mode , **kwargs ) :
        return open ( name , mode )
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
    """Convert dataset to CSV format
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
        for entry in progress_bar ( dataset , max_value = len ( dataset ) , silent = not progress ) :
            
            values =  [ entry [ a ].getVal() for a in vnames1 ]
            values += [ v.getVal()           for v in mvars   ]

            if   weighted and sae :
                e1 , e2 = dataset.weight_errors () 
                values += [ dataset.weight() , e1 , e2 ]
            elif weighted and se  :
                e1 , e2 = dataset.weight_errors () 
                values += [ dataset.weight() , 0.5*(e1+e2) ]
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
def ds_to_tree ( dataset , weight = '' , filename = '' ) :
    """Get dataset as `ROOT.TTree`
    - for Tree-based datasets gets the internal tree
    - otherwise tree will be created
    
    >>> dataset = ...
    >>> tree    = dataset.asTree() 
    """
    
    import ostap.io.root_file
    import ostap.trees.trees
    
    ## first, unweight it, if weighted 
    if dataset.isWeighted () :
        
        with useStorage ( RAD.Tree ) , ROOTCWD() :
            
            if not filename :
                import ostap.utils.cleanup as CU 
                filename = CU.CleanUp.tempfile ( suffix = '.root' )
                logger.info ( "Temporary ROOT file is created: %s" % filename ) 

            with ROOT.TFile ( filename , 'u' ) as rfile :
                
                rfile.cd ()
                uwds , wname = dataset.unweight ( weight )
                
                tree = uwds.GetClonedTree()
                
                tree.SetName ( 'tree_%s' % dataset.GetName() )
                tree .Write()
                ## rfile.Write() 
                tname = tree.GetName()
            
                del uwds
                logger.debug ( 'ROOT file %s\n%s' % ( filename , rfile.as_table ( prefix = '# ' ) ) ) 
                
            chain = ROOT.TChain ( tname )
            chain.AddFile ( filename )
            return chain 

    
    store = dataset.store()
    if hasattr ( store , 'tree' ) :
        tree = store.tree()
        if valid_pointer ( tree ) and len ( tree ) == len ( dataset ) :
            return dataset.GetClonedTree()
        
    with useStorage ( RAD.Tree ) , ROOTCWD() :
        
        if not filename :
            import ostap.utils.cleanup as CU 
            filename = CU.CleanUp.tempfile ( suffix = '.root' )
            logger.info ( "Temporary ROOT file is created: %s" % filename ) 
            
        with ROOT.TFile ( filename , 'u' ) as rfile :

            rfile.cd ()
            tree = dataset.GetClonedTree()
            tree.SetName ( 'tree_%s' % dataset.GetName() )
            tree .Write()
            ## rfile.Write()
            tname = tree.GetName()
            logger.debug ( 'ROOT file %s\n%s' % ( filename , rfile.as_table ( prefix = '# ' ) ) ) 
            
        chain = ROOT.TChain ( tname )
        chain.AddFile ( filename )
        return chain 
            
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
def _rda_slice_ ( dataset , variables , cuts = '' , transpose = False , cut_range = '' , *args ) :
    """Get "slice" from <code>RooAbsData</code> in form of numpy array
    >>> data = ...
    >>> varr , weights = data.slice ( 'a : b : c' , 'd>0' )
    >>> varr , weights = data.slice ( 'a , b , c' , 'd>0' )
    >>> varr , weights = data.slice ( 'a ; b ; c' , 'd>0' )
    """
    
    if isinstance ( variables , string_types ) :
        variables = split_string ( variables , var_separators , strip = True , respect_groups = True )
    
    names = []
    for v in variables :
        names += split_string ( v , var_separators , strip = True , respect_groups = True )

    names = strings ( names )
    
    
    tab = Ostap.StatVar.Table  ()
    col = Ostap.StatVar.Column ()

    ## get data 
    n   = Ostap.StatVar.get_table ( dataset   ,
                                    names     ,
                                    cuts      ,
                                    tab       ,
                                    col       , 
                                    cut_range ,
                                    *args     )
    nc = len ( col )    
    assert ( dataset.isWeighted() and nc == n ) or ( 0 == nc and not dataset.isWeighted() ), \
           'slide: invalid size of ``weights'' column! %s/%s/%s' % ( n , nc , dataset.isWeighted() ) 
    
    if 0 == n :
        return () , () 

    import numpy

    result = []
    for column in tab :
        
        nc = len ( column )
        assert 0 <= nc and nc == n, 'slice: invalid column size! %s/%s' % ( nc , n ) 

        l = column.begin().__follow__()
        
        result.append ( numpy.array ( numpy.frombuffer ( l , count = n ) , copy = True ) )
        
        column.clear()
        
    del tab 

    
    if not result :
        return None, None 
    
    result = numpy.stack ( result )

    if transpose : result = result.transpose()

    if not dataset.isWeighted() :
        return result, None 
        
    nn = len ( col ) 
    l  = col.begin().__follow__()        
    weights = numpy.array ( numpy.frombuffer ( l , count = nn ) , copy = True )
    col.clear()
    
    del col
    
    return result , weights 


ROOT.RooAbsData.slice = _rda_slice_

_new_methods_ += [
    ROOT.RooAbsData.slice 
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
    """Ninned dataset ?
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
    """Are two datasets equal by content?
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
## Are two datastes non-equal by content?
#  @code
#  ds1 = ...
#  ds2 = ...
#  ds_nonequal ( ds1 , ds2 )
#  ds1 != ds2 
#  @endcode 
def ds_nonequal ( ds1 , ds2 ) : 
    """Are two datastes non-equal by content?
    ds1 = ...
    ds2 = ...
    ds_nonequal ( ds1 , ds2 ) 
    ds1 != ds2 
    """
    return not ds_equal ( ds1 , ds2 ) 

# =============================================================================
## Are two datasets equal by content?
def _ds_eq_ ( ds1 , ds2 ) :
    """Are two datasets equal by content?"""
    if not isinstance  ( ds2 , ROOT.RooAbsData ) : return NotImplemented 
    return ds_equal ( ds1 , ds2 ) 

# =============================================================================
## Are two datasets non-equal by content?
def _ds_ne_ ( ds1 , ds2 ) :
    """Are two datasets non-equal by content?"""
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



# ============================================================================
try : 
    from numpy import array as _array 
    def get_result ( data ) :
        return _array (  data , dtype = float ) 
except ImportError :
    from array import array as _array 
    def get_result ( data ) : return _array ( 'd' , data )


# ===========================================================================
## Iterator for rows in dataset
#  @code
#  dataset = ...
#  for row , weight in dataset.rows ( 'pt, pt/p, mass ' , 'pt>1' ) :
#     print (row, weight) 
#  @endcode 
def _rad_rows_ ( dataset , variables = [] , cuts = '' , cutrange = '' , first = 0 , last = -1 ) :
    """Iterator for rows in dataset
    >>> dataset = ...
    >>> for row , weight in dataset.rows ( 'pt, pt/p, mass ' , 'pt>1' ) :
    >>>    print (row, weight) 
    """
    
    if last < 0 : last = ROOT.TTree.kMaxEntries    
    last  = min ( last , len ( dataset ) )
    first = max ( 0    , first           ) 

    if isinstance ( variables , string_types ) :
        variables = split_string ( variables , var_separators , strip = True , respect_groups = True )
    vars = []
    for v in variables :
        vars += split_string ( v , var_separators , strip = True , respect_groups = True )
    vars = strings ( vars ) 

    formulas = []
    varlist  = dataset.varlist () 
    for v in vars :
        f0 = make_formula     ( v , v , varlist  ) 
        formulas.append ( f0 ) 

    fcuts = None 
    if cuts : fcuts = make_formula ( cuts , cuts , varlist  )

    weighted = dataset.isWeighted()
    
    ## loop over dataset 
    for event in range ( first , last ) :
        
        vars = dataset.get ( event ) 
        if not vars : break

        if cutrange and not vars.allInRange ( cutrange ) : continue

        wc = fcuts.getVal() if fcuts else 1.0 

        if not wc   : continue

        wd = dataset.weight() if weighted else 1.0 
        
        w  = wc * wd 
        if not w    : continue

        weight = w if weighted else None 

        result = tuple ( tuple ( float(f) for f in formulas ) ) 
        yield  get_result ( result ) , weight 
        

    del fcuts
    del formulas

    
ROOT.RooAbsData.rows = _rad_rows_         
    
_new_methods_ += [
    ROOT.RooAbsData.rows 
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
