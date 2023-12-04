#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/ds2numpy.py
#  Helper module to convert RooDtaSet to numpy array
#  @see RooAbsData 
#  @see RooDataSet 
#  @author Artem Egorychev Artem.Egorychev@cern.ch 
#  @date   2023-12-04
# =============================================================================
""" Helper module to convert RooDtaSet to numpy array
Module with decoration for RooAbsData and related RooFit classes
- see RooAbsData 
- see RooDataSet 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Artem Egorychev Artem.Egorychev@cern.ch"
__date__    = "2011-06-07"
__all__     = (
    )
# =============================================================================
from   ostap.core.meta_info         import root_info
from   ostap.core.ostap_types       import string_types 
import ostap.fitting.roocollections
import ostap.fitting.dataset 
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.ds2numpy' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug( 'Convert RooDataSet uito numpy array')

# =============================================================================
try :
    import numpy as np
except ImportError :
    np = None
# =============================================================================
_new_methods_ = [] 



# input: 
# dataset - initial dataset
# var_lst - name of variables to add in numpy array 
# weight - Bool value, which work with weights vars in dataset
#ds = DS_to_Numpy(data, ['evt', 'run'], weight)
#ds = DS_to_Numpy_for_old_version(data, ['evt', 'run']) - for old ROOT package version
# output: 
# data - numpy array with values of the required variables

#Check the list of variables for duplicates
def find_dublicates_in_var_list(var_lst):
    return len(var_lst) != len(set(var_lst))

#add weight variable in numpy array
def add_weight ( ds , data ):
    
    if not ds.isWeighted() : return data
    
    weight  = ds.weightVar().GetName()
    
    ## creathe the weight array 
    weights = np.zeros( len ( ds ) , dtype=np.float64)

    ## fill it 
    for i in ds : weight_array[i] = ds.weight()

    data [ weight ] = weights 

    return data

# =============================================================================
if   np and ( 6 , 26 ) <= root_info :  ## 6.26 <= ROOT 
# =============================================================================

    # =========================================================================
    ## Convert dataset into numpy array using <code>ROOT.RooDataSet.to_numpy</code>
    #  @see `ROOT.RooDataSet.to_numpy`
    # 
    def ds2numpy ( dataset , var_lst ) :
        """ Convert dataset into numpy array using `ROOT.RooDataSet.to_numpy` methdod from new ROOT
        - see `ROOT.RooDataSet.to_numpy` 
        """
        
        ## 1) check that all variables are present in dataset 
        if   all ( isinstance ( v , string_types   ) for v in var_lst ) :
            vnames = [ v          for  v in var_lst ]
        elif all ( isinstance ( v , ROOT.RooAbsArg ) for v in var_lst ) :
            vnames = [ v.GetName() for v in var_lst ]
        else :
            raise TypeError ( "Invalid type of `var_list`!" ) 

        ## 2) check that all variables are present in the dataset 
        assert all ( ( v in dataset ) for v in var_lst ) , 'Not all variables are in dataset!'

        ## remove duplicated
        new_names = []
        for v in vnames :
            if not v in new_names : new_names.append ( v )
        vnames = new_names
        
        ## 3) reduce dataset if only a small subset of variables is requested 
        nvars = len( dataset.get() )
        if 2 * len ( vnames )  <= nvars :
            dstmp  = dataset.subset ( vnames )
            result = ds2numpy ( dstmp , vnames )
            dstmp.erase()
            del dstmp 
            return result  

        ## 4) convert to numpy 
        data = dataset.to_numpy()
        
        ## 5) convert to named/structured array 
        
        dtypes = [ ( name , data [ name ].dtype ) for name in vnames if name in data ]
        lst    = [          data [ name ]         for name in vnames if name in data ]

        ## 6) add the weight
        if dataset.isWeighted() : 
            weight = dataset.weightVar().GetName()
            if not weight in vnames : 
                dtypes.append ( ( weight , data [ weight ] .dtype ) )
                lst   .append (            data [ weight ]          )

        ## is there a better way to avoid a creation of lists ??? 
        data  = np.array ( list ( zip ( *lst ) ) , dtype = dtypes )
        
        return data
    

    __all__  = __all__ + ( 'ds2numpy' , )
    ROOT.RooDataSet.tonumpy = ds2numpy 
    ROOT.RooDataSet.tonp    = ds2numpy
    ROOT.RooDataSet.to_np   = ds2numpy
    _new_methods_ += [ ROOT.RooDataSet.tonumpy ,
                       ROOT.RooDataSet.tonp    ,
                       ROOT.RooDataSet.to_np   ]
    
# =============================================================================
elif np and ( 6 , 24 ) <= root_info : ## 6.24 <= ROOT < 6.26 
# =============================================================================

    # =========================================================================
    ## Convert dataset into numpy array using <code>ROOT.RooVectorDataStore.getArrays</code>
    #  @see `ROOT.RooVectorDataStore.getArrays`
    # 
    def ds2numpy ( dataset , var_lst ) :
        """ Convert dataset into numpy array using `ROOT.RooVectorDataStore.getArrays`
        - see `ROOT.RooVectorDataStore.getArrays` 
        """

        ## 1) check that all variables are present in dataset 
        if   all ( isinstance ( v , string_types   ) for v in var_lst ) :
            vnames = [ v          for  v in var_lst ]
        elif all ( isinstance ( v , ROOT.RooAbsArg ) for v in var_lst ) :
            vnames = [ v.GetName() for v in var_lst ]
        else :
            raise TypeError ( "Invalid type of `var_list`!" ) 

        ## 2) check that all variables are present in the dataset 
        assert all ( ( v in dataset ) for v in var_lst ) , 'Not all variables are in dataset!'

        ## remove duplicated
        new_names = []
        for v in vnames :
            if not v in new_names : new_names.append ( v )
        vnames = new_names
        
        ## 3) reduce dataset if only a small subset of variables is requested 
        nvars = len( dataset.get() )
        if 2 * len ( vnames )  <= nvars :
            dstmp  = dataset.subset ( vnames )
            result = ds2numpy ( dstmp , vnames )
            dstmp.erase()
            del dstmp 
            return result  

        ## 4) here we need RooVectorDataStore 
        store = dataset.store()
        if not isinstance ( store , ROOT.RooVectorDataStore ) : 
            dataset.ConvertToVectorStore()
            store = dataset.store()
            
        new_store = False 
        if not isinstance ( store , ROOT.RooVectorDataStore ) : 
            variables  = store.get()
            store      = ROOT.RooVectorDataStore ( store, variables , store.GetName() )
            new_store  = True

        ## 5) get arrays from the store 

        array_info = store.getArrays()
        n          = array_info.size

        ## 6) using numpy structured array
        dtypes = [ ( name , 'f8') for name in vnames ]
        
        ## 7) weight?
        if dataset.isWeighted() : 
            weight = dataset.weightVar().GetName()
            if not weight in vnames : 
                dtypes.append ( ( weight , 'f8' ) )

        ## 8) create the structured array 
        data   = np.zeros ( len ( dtypes ) , dtype = dtypes )
        
        for x in array_info.reals:
            if x.name in vnames :
                data [ x.name ] = np.frombuffer ( x.data , dtype = np.float64 , count = n )
                
        for x in array_info.cats:
            if x.name in vnames :
                data [ x.name ] = np.frombuffer ( x.data , dtype = np.int32   , count = n )

        if new_store : ## delete newly created store 
            store.reset() 
            del store

        ## check here!!! 
        return add_weight ( dataset , data )
    

    __all__  = __all__ + ( 'ds2numpy' , )
    ROOT.RooDataSet.tonumpy = ds2numpy 
    ROOT.RooDataSet.tonp    = ds2numpy
    ROOT.RooDataSet.to_np   = ds2numpy
    _new_methods_ += [ ROOT.RooDataSet.tonumpy ,
                       ROOT.RooDataSet.tonp    ,
                       ROOT.RooDataSet.to_np   ]
    
    
# =============================================================================
elif np and ( 6, 20 ) <= root_info :  ## 6.20 <= ROOT < 6.24 
# =============================================================================

    # =========================================================================
    ## Convert dataset into numpy array using <code>ROOT.RooVectorDataStore.getBatches</code>
    #  @see `ROOT.RooVectorDataStore.getBatches`
    # 
    def ds2numpy ( dataset , var_lst ) :
        """ Convert dataset into numpy array using `ROOT.RooVectorDataStore.getBatches`
        - see `ROOT.RooVectorDataStore.getBatches` 
        """
        

        ## 1) check that all variables are present in dataset 
        if   all ( isinstance ( v , string_types   ) for v in var_lst ) :
            vnames = [ v          for  v in var_lst ]
        elif all ( isinstance ( v , ROOT.RooAbsArg ) for v in var_lst ) :
            vnames = [ v.GetName() for v in var_lst ]
        else :
            raise TypeError ( "Invalid type of `var_list`!" ) 

        ## 2) check that all variables are present in the dataset 
        assert all ( ( v in dataset ) for v in var_lst ) , 'Not all variables are in dataset!'

        ## remove duplicated
        new_names = []
        for v in vnames :
            if not v in new_names : new_names.append ( v )
        vnames = new_names
        
        ## 3) reduce dataset if only a small subset of variables is requested 
        nvars = len( dataset.get() )
        if 2 * len ( vnames )  <= nvars :
            dstmp  = dataset.subset ( vnames )
            result = ds2numpy ( dstmp , vnames )
            dstmp.erase()
            del dstmp 
            return result  

        ## 4) here we need RooVectorDataStore 
        store = dataset.store()
        if not isinstance ( store , ROOT.RooVectorDataStore ) : 
            dataset.ConvertToVectorStore()
            store = dataset.store()
            
        new_store = False 
        if not isinstance ( store , ROOT.RooVectorDataStore ) : 
            variables  = store.get()
            store      = ROOT.RooVectorDataStore ( store, variables, store.GetName() )
            new_store  = True

            
        #$ 5) using numpy structed array
        dtypes = [ ( name , 'f8') for name in vnames ]

        ## 6) weight?
        weight = None 
        if dataset.isWeighted() : 
            weight = dataset.weightVar().GetName()
            if not weight in vnames : 
                dtypes.append ( ( weight , 'f8' ) )

        ## 7) book the array 

        # for large datasets
        # check batch size * var size < 10^6 
        num_entries = len ( dataset ) 
        data_limit  = num_entries * nvars 
        num_limit   = 110000000
        nb , r      = divmod ( n , num_limit )

        ##
        ##
        ## REWRITE: should be RunContext here!!! 
        ##
        
        if data_limit < num_limit:
            data    = np.zeros ( len ( dtypes )  , dtype = dtypes )            
            batches = store.getBatches ( 0 , n)
            count   = 0
            for name in vnames :
                for x in batches :
                    if name == x.first.__follow__().GetName() :
                        data [ name ] = x.second
                        break
            if weight :
                data [ weight ] = store.getWeightBatch ( 0 , n )
            
        else:
            
            rargs = [ ( i * num_limit , num_limit ) for i in range ( nb ) ] + [ ( nb * num_limit , r ) ]

            data  = None 
            for first , num in rargs :

                part = np.zeros ( num , dtype = dtypes )                 
                batches = store.getBatches  ( first, num)
                for x in vnames :
                    for x in batches :
                        if name == x.first.__follow__().GetName() :
                            part [ name ] = x.second
                            break
                if weight : part [ weight ] = store.getWeightBatch ( 0 , n )

                if data : data = np.concatenate ( [ data , part ] )
                else    : data = part 
                
        if new_store : ## delete newly created store 
            store.reset() 
            del store
            
        return data

    __all__  = __all__ + ( 'ds2numpy' , )
    ROOT.RooDataSet.tonumpy = ds2numpy 
    ROOT.RooDataSet.tonp    = ds2numpy
    ROOT.RooDataSet.to_np   = ds2numpy
    _new_methods_ += [ ROOT.RooDataSet.tonumpy ,
                       ROOT.RooDataSet.tonp    ,
                       ROOT.RooDataSet.to_np   ]
    
# =============================================================================
else    :
# =============================================================================

    logger.warning ( "Numpy is not available: no action" )


# =============================================================================
_decorated_classes_ = (
    ROOT.RooDataSet , 
    )
_new_methdos_ = tuple ( _new_methods_ ) 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
