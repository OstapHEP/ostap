#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/ds2numpy.py
#  Helper module to convert RooDtaSet to numpy array
#  @see RooAbsData 
#  @see RooDataSet 
#  @author Artem Egorychev Artem.Egorychev@cern.ch 
#  @date   2023-12-12
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
from   ostap.core.ostap_types       import string_types, dictlike_types
from   ostap.core.core              import loop_items 
from   ostap.utils.utils            import split_range 
from   ostap.fitting.dataset        import useStorage
from   ostap.fitting.funbasic       import AFUN1 
from   ostap.utils.progress_bar     import progress_bar 
import ostap.fitting.roocollections
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



# =============================================================================
if   np and ( 6 , 28 ) <= root_info  :  ## 6.26 <= ROOT 
# =============================================================================

    # =========================================================================
    ## Convert dataset into numpy array using <code>ROOT.RooAbsData</code> interface 
    #  @see ROOT.RooAbsData.getBatches
    #  @see ROOT.RooAbsData.getCategoryBatches
    #  @see ROOT.RooAbsData.getWeightBatche    
    #  @see ROOT.RooAbsDataStore.getBatches
    #  @see ROOT.RooAbsDataStore.getCategoryBatches
    #  @see ROOT.RooAbsDataStore.getWeightBatche
    #  @attention conversion to ROOT.RooVectorDataStore is used! 
    def ds2numpy ( dataset , var_lst , silent = True , more_vars = {} ) :
        """ Convert dataset into numpy array using `ROOT.RooAbsData` iterface 
        - see ROOT.RooAbsData.getBatches
        - see ROOT.RooAbsData.getCategoryBatches
        - see ROOT.RooAbsData.getWeightBatche    
        - see ROOT.RooAbsDataStore.getBatches
        - see ROOT.RooAbsDataStore.getCategoryBatches
        - see ROOT.RooAbsDataStore.getWeightBatche    
        - attention: Conversion to `ROOT.RooVectorDataStore` is used! 
        """
        
        ## 1) get names of all requested variables
        if   all ( isinstance ( v , string_types   ) for v in var_lst ) :
            vnames = [ v          for  v in var_lst ]
        elif all ( isinstance ( v , ROOT.RooAbsArg ) for v in var_lst ) :
            vnames = [ v.GetName() for v in var_lst ]
        else :
            raise TypeError ( "Invalid type of `var_list`!" ) 

        ## 2) check that all variables are present in the dataset 
        assert all ( ( v in dataset ) for v in var_lst ) , 'Not all variables are in dataset!'

        funcs = [] 
        if more_vars and isinstance ( more_vars , dictlike_types ) :
            for name , fun in loop_items ( more_vars ) :
                if isinstance ( fun , AFUN1 ) :
                    absreal = fun.fun
                elif isinstance( fun , ROOT.RooAbsPdf  ) : absreal = fun
                elif isinstance( fun , ROOT.RooAbsReal ) : absreal = fun
                else :
                    raise TypeError ( "Invald type ofun/pdf" )
                obsvars = absreal.getObservables ( dataset )
                item    = name , absreal , obsvars
                funcs.append ( item )
        elif more_vars and all ( ( isinstance ( v , ( ROOT.RooAbsReal , AFUN1 ) ) for v in more_vars ) ) :
            for v in more_vars :
                if isinstance ( v  , AFUN1 ) : absreal = v.fun
                else                         : absreal = v  
                obsvars = absreal.getObservables ( dataset )
                item    = v.name , absreal , obsvars
                funcs.append ( item )
        elif more_vars :
            for name, var in more_vars :
                if   isinstance ( var , AFUN1           ) : absreal = var.fun
                elif isinstance ( var , ROOT.RooAbsReal ) : absreal = var   
                else :
                    raise TypeError ( "Invalid content of 'more_vars'!" )
                obsvars = absreal.getObservables ( dataset )
                item    = name , absreal , obsvars
                funcs.append ( item )
                
        ## 3) reduce dataset if only a small subset of variables is requested 
        nvars = len ( dataset.get() )
        if 2 * len ( vnames )  <= nvars :            
            with useStorage ( ROOT.RooAbsData.Vector ) : 
                dstmp  = dataset.subset ( vnames )
            result = ds2numpy ( dstmp , vnames , more_vars = more_vars )
            dstmp.erase()
            del dstmp 
            return result

        ## 4) convert to VectorDataStore
        #  batches are not (yet) implemented for Tree & Composite stores 
        dataset.convertToVectorStore()
        ## dataset.convertToTreeStore()

        ## 5) convert to VectorStore again...
        #  batches are not (yet) implemented for Tree & Composite stores         
        store     = dataset.store()
        source    = dataset
        twoargs   = False
        delsource = False 
        if not isinstance ( store , ROOT.RooVectorDataStore ) :
            source    = ROOT.RooVectorDataStore ( store , dataset.get() , store.name + '_vct' )
            twoargs   = True 
            delsource = True 
        
        vars       = source.get()
        vars       = [ v for v in vars if v.name in vnames ]
        doubles    = [ v.name for v in vars if isinstance ( v , ROOT.RooAbsReal     ) ] 
        categories = [ v.name for v in vars if isinstance ( v , ROOT.RooAbsCategory ) ]

        ## name of weight variable 
        weight = '' if not dataset.isWeighted() else dataset.wname () 

        dtypes = [] 
        for v in vnames :
            if   v in doubles    : dtypes.append ( ( v      , np.float64 ) ) 
            elif v in categories : dtypes.append ( ( v      , np.int64   ) )

        for f in funcs           : dtypes.append ( ( f[0]   , np.float64 ) ) 
        if weight                : dtypes.append ( ( weight , np.float64 ) ) 
            
        ## get data in batches
        nevts  = len ( dataset ) 

        data   = None

        ## maximal size of data chunk 
        nmax   = max ( nevts // 6 , 30000 // nvars )
        
        ## get data is chunks/batches 
        for first, last in progress_bar ( split_range ( 0 , nevts , nmax ) , silent = silent ) :
            
            num   = last - first            
            wget  = False 
            part  = np.zeros ( num , dtype = dtypes )

            if doubles :
                dpart   = source.getBatches ( first , num )
                for d in dpart :
                    dname = d.first.name
                    if   dname in doubles :
                        part [ dname ] = d.second

                        
                    elif dname == weight  :
                        part [ dname ] = d.second
                del dpart
                
            if categories :
                cpart   = source.getCategoryBatches ( first , num )
                for c in cpart :
                    cname = c.first.name
                    if cname in categories : 
                        part [ cname ] = c.second
                del cpart
                
            if weight and not wget :
                if twoargs : weights = source.getWeightBatch ( first , num         )
                else       : weights = source.getWeightBatch ( first , num , False )
                if weights : part [ weight ] = weights 
                else       : part [ weight ] = np.full ( num , source.weight() , dtype = np.float64 )

            if data is None : data = part
            else            :  
                data = np.concatenate ( [ data , part ] )
                del part

        ## add PDF values
        if funcs :  
            for i, entry in enumerate ( source ) :
                for vname , func , obsvars in funcs :
                    obsvars.assign ( entry )
                    data [ vname ] [ i ] = func.getVal()   
                
        if delsource : 
            source.reset()
            del source
            
        return data
    

    __all__  = __all__ + ( 'ds2numpy' , )
    ROOT.RooDataSet.tonumpy = ds2numpy 
    ROOT.RooDataSet.tonp    = ds2numpy
    ROOT.RooDataSet.to_np   = ds2numpy
    _new_methods_ += [ ROOT.RooDataSet.tonumpy ,
                       ROOT.RooDataSet.tonp    ,
                       ROOT.RooDataSet.to_np   ]
    
# =============================================================================
elif   np  :  ## ROOT < 6.26 
# =============================================================================

    # =========================================================================
    ## Convert dataset into numpy array using (slow) explicit loops 
    def ds2numpy ( dataset , var_lst , silent = False , more_vars = {} ) :
        """ Convert dataset into numpy array using (slow) explicit loops
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
 
        funcs = [] 
        if more_vars and isinstance ( more_vars , dictlike_types ) :
            for name , fun in loop_items ( more_vars ) :
                if isinstance ( fun , AFUN1 ) :
                    absreal = fun.fun
                elif isinstance( fun , ROOT.RooAbsPdf  ) : absreal = fun
                elif isinstance( fun , ROOT.RooAbsReal ) : absreal = fun
                else :
                    raise TypeError ( "Invald type ofun/pdf" )
                obsvars = absreal.getObservables ( dataset )
                item    = name , absreal , obsvars
                funcs.append ( item )
        elif more_vars and all ( ( isinstance ( v , ( ROOT.RooAbsReal , AFUN1 ) ) for v in more_vars ) ) :
            for v in more_vars :
                if isinstance ( v  , AFUN1 ) : absreal = v.fun
                else                         : absreal = v  
                obsvars = absreal.getObservables ( dataset )
                item    = v.name , absreal , obsvars
                funcs.append ( item )
        elif more_vars :
            for name, var in more_vars :
                if   isinstance ( var , AFUN1           ) : absreal = var.fun
                elif isinstance ( var , ROOT.RooAbsReal ) : absreal = var   
                else :
                    raise TypeError ( "Invalid content of 'more_vars'!" )
                obsvars = absreal.getObservables ( dataset )
                item    = name , absreal , obsvars
                funcs.append ( item )
                
        ## 3) reduce dataset if only a small subset of variables is requested 
        nvars = len ( dataset.get() )
        if 2 * len ( vnames )  <= nvars :
            dstmp  = dataset.subset ( vnames )
            result = ds2numpy ( dstmp , vnames , more_vars = more_vars )
            dstmp.erase()
            del dstmp 
            return result  
        
        vars       = dataset.get()
        vars       = [ v for v in vars if v.name in vnames ]
        doubles    = [ v.name for v in vars if isinstance ( v , ROOT.RooAbsReal     ) ] 
        categories = [ v.name for v in vars if isinstance ( v , ROOT.RooAbsCategory ) ]

        ## name of weight variable 
        weight = '' if not dataset.isWeighted() else dataset.wname () 

        dtypes = [] 
        for v in vnames :
            if   v in doubles    : dtypes.append ( ( v      , np.float64 ) ) 
            elif v in categories : dtypes.append ( ( v      , np.int64   ) )
        for f in funcs           : dtypes.append ( ( f[0]   , np.float64 ) ) 
        if weight                : dtypes.append ( ( weight , np.float64 ) ) 
            
        
        ## create data 
        data = np.zeros ( len ( dataset )  , dtype = dtypes )

        ## make an explict loop:
        for i , evt in enumerate ( progress_bar ( dataset , silent = silent ) ) :

            for v in evt :
                vname = v.name
                if   vname in doubles    : data [ vname  ] [ i ] = float ( v  )
                elif vname in categories : data [ vname  ] [ i ] = int   ( v  )
                
            if weight                    : data [ weight ] [ i ] = dataset.weight()  

            ## add PDF values
            for vname , func , obsvars in funcs :
                obsvars.assign ( evt )
                data [ vname ] [ i ] = func.getVal()   
        
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
