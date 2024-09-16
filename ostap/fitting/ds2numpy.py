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
from   ostap.core.ostap_types       import string_types, dictlike_types, sized_types
from   ostap.core.core              import Ostap, loop_items 
from   ostap.utils.utils            import split_range 
from   ostap.fitting.dataset        import useStorage
from   ostap.fitting.funbasic       import AFUN1 
from   ostap.utils.progress_bar     import progress_bar
from   ostap.trees.cuts             import vars_and_cuts 
import ostap.fitting.roocollections
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.ds2numpy' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug( 'Convert RooDataSet into numpy array')

# =============================================================================
try :
    import numpy as np
except ImportError :
    np = None
# =============================================================================
_new_methods_ = []

# =============================================================================
if  np and ( 6 , 28 ) <= root_info  :  ## 6.26 <= ROOT     
    # =========================================================================
    ## Convert dataset into numpy array using <code>ROOT.RooAbsData</code> interface 
    #  @see ROOT.RooAbsData.getBatches
    #  @see ROOT.RooAbsData.getCategoryBatches
    #  @see ROOT.RooAbsData.getWeightBatche    
    #  @see ROOT.RooAbsDataStore.getBatches
    #  @see ROOT.RooAbsDataStore.getCategoryBatches
    #  @see ROOT.RooAbsDataStore.getWeightBatche
    #  @attention conversion to ROOT.RooVectorDataStore is used! 
    def ds2numpy ( dataset          , 
                   var_lst          ,
                   cuts      = ''   ,
                   cut_range = ''   , 
                   silent    = True ,
                   more_vars = {}   ) :
        """ Convert dataset into numpy array using `ROOT.RooAbsData` iterface 
        - see ROOT.RooAbsData.getBatches
        - see ROOT.RooAbsData.getCategoryBatches
        - see ROOT.RooAbsData.getWeightBatche    
        - see ROOT.RooAbsDataStore.getBatches
        - see ROOT.RooAbsDataStore.getCategoryBatches
        - see ROOT.RooAbsDataStore.getWeightBatche    
        - attention: Conversion to `ROOT.RooVectorDataStore` is used! 
        """
        
        if isinstance ( var_lst , string_types ) : var_lst = [ var_lst ]
        
        # =====================================================================
        ## 1) get names of all requested variables
        if   all ( isinstance ( v , string_types   ) for v in var_lst ) :
            vnames , cuts , _ = vars_and_cuts ( var_lst , cuts ) 
        elif all ( isinstance ( v , ROOT.RooAbsArg ) for v in var_lst ) :
            vnames = [ v.GetName() for v in var_lst ]
        else :
            raise TypeError ( "Invalid type of `var_list`!" ) 

        cuts      = str ( cuts      ).strip () 
        cut_range = str ( cut_range ).strip () if cut_range else ''
        
        # =====================================================================
        ## 2) check that all variables are present in the dataset 
        assert all ( ( v in dataset ) for v in vnames ) , 'Not all variables are in dataset!'

        # =====================================================================
        ## 3) reduce dataset if only a small subset of variables is requested 
        nvars = len ( dataset.get() )
        if 2 * len ( vnames )  <= nvars and not more_vars  :            
            with useStorage ( ROOT.RooAbsData.Vector ) : 
                dstmp  = dataset.subset ( vnames , cuts = cuts , cut_range = cut_range )
            result = ds2numpy ( dstmp , vnames , more_vars = more_vars )
            ## dstmp.erase()
            dstmp = Ostap.MoreRooFit.delete_data ( dstmp ) 
            del dstmp 
            return result

        # =========================================================================
        ## 4) if cuts or cut-range is specified, assume cuts are hash and make a filtering 
        if cuts or cut_range :
            with useStorage ( ROOT.RooAbsData.Vector ) :
                dstmp = dataset.subset ( vnames if not more_vars else [] , cuts = cuts , cut_range = cut_range )
            result = ds2numpy ( dstmp , vnames , more_vars = more_vars )
            ## dstmp.erase()
            dstmp = Ostap.MoreRooFit.delete_data ( dstmp ) 
            del dstmp 
            return result

        funcs    = []
        formulas = []
        varlst   = dataset.varlst () 
        if more_vars and isinstance ( more_vars , dictlike_types ) :
            for name , fun in loop_items ( more_vars ) :
                assert not name in dataset, 'da2numpy: no way to redefine variable `%s`!' % name 
                if   isinstance ( fun , AFUN1           ) : absreal = fun.fun
                elif isinstance ( fun , ROOT.RooAbsPdf  ) : absreal = fun
                elif isinstance ( fun , ROOT.RooAbsReal ) : absreal = fun
                elif isinstance ( fun , string_types    ) :
                    formula = Ostap.FormulaVar ( name , fun , varlst )
                    formulas.append ( formula ) 
                    absreal = formula
                else :
                    raise TypeError ( "Invalid type of fun/pdf" )
                obsvars = absreal.getObservables ( dataset )
                item    = name , absreal , obsvars
                funcs.append ( item )
        elif more_vars :
            for item in more_vars :                
                len2 = isinstance ( item , sized_types ) and 2 == len ( item )
                name = '' 
                if   isinstance ( item , AFUN1           ) : absreal = item.fun
                elif isinstance ( item , ROOT.RooAbsPdf  ) : absreal = item
                elif isinstance ( item , ROOT.RooAbsReal ) : absreal = item
                elif len2 and isinstance ( item [ 0 ] , string_types ) and isinstance ( item [1] , AFUN1          ) :
                    absreal, name = item[1].fun , item[0]
                elif len2 and isinstance ( item [ 0 ] , string_types ) and isinstance ( item [1] , ROOT.RooAbsPdf ) :
                    absreal, name = item[1]     , item[0]
                elif len2 and isinstance ( item [ 0 ] , string_types ) and isinstance ( item [1] , ROOT.RooAbsPdf ) :
                    absreal, name = item[1]     , item[0]
                elif len2 and isinstance ( item [ 0 ] , string_types ) and isinstance ( item [1] , string_types   ) :
                    formula = Ostap.FormulaVar ( item [ 0 ] , item [ 1 ] , varlst )
                    formulas.append ( formula ) 
                    absreal, name = formula , item [ 0 ]                    
                else :
                    raise TypeError ( "Invalid type of fun/pdf!" )                
                obsvars = absreal.getObservables ( dataset )
                item    = name if name else absreal.name , absreal , obsvars
                funcs.append ( item )

        # =====================================================================
        for item in funcs :
            name , _ , _ = item 
            assert not name in dataset, 'da2numpy: no way to redefine variable `%s`!' % name 

        # =========================================================================
        ## 5) convert to VectorDataStore
        #  batches are not (yet) implemented for Tree & Composite stores 
        dataset.convertToVectorStore()
        ## dataset.convertToTreeStore()

        # =====================================================================
        ## 6) convert to VectorStore again...
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
            for i, item  in enumerate ( source ) :
                entry , weight = item 
                for vname , func , obsvars in funcs :
                    obsvars.assign ( entry )
                    data [ vname ] [ i ] = func.getVal()   
                
        if delsource :
            source.reset()
            del source
            
        del funcs
        del formulas 
        return data

    # =========================================================================
    __all__  = __all__ + ( 'ds2numpy' , )
    ROOT.RooDataSet.tonumpy = ds2numpy 
    ROOT.RooDataSet.tonp    = ds2numpy
    ROOT.RooDataSet.to_np   = ds2numpy
    _new_methods_ += [ ROOT.RooDataSet.tonumpy ,
                       ROOT.RooDataSet.tonp    ,
                       ROOT.RooDataSet.to_np   ]
    
# =============================================================================
elif   np  :  ## ROOT < 6.26 
    # =========================================================================
    ## Convert dataset into numpy array using (slow) explicit loops 
    def ds2numpy ( dataset           ,
                   var_lst           ,
                   cuts      = ''    ,
                   cut_range = ''    , 
                   silent    = False ,
                   more_vars = {}    ) :
        """ Convert dataset into numpy array using (slow) explicit loops
        """
        

        if isinstance ( var_lst , string_types ) : var_lst = [ var_lst ]
        
        # =====================================================================
        ## 1) check that all variables are present in dataset 
        if   all ( isinstance ( v , string_types   ) for v in var_lst ) :
            vnames , cuts , _  = vars_and_cuts ( vnames , cuts )
        elif all ( isinstance ( v , ROOT.RooAbsArg ) for v in var_lst ) :
            vnames = [ v.GetName() for v in var_lst ]
        else :
            raise TypeError ( "Invalid type of `var_list`!" ) 

        cuts      = str ( cuts      ).strip () 
        cut_range = str ( cut_range ).strip () if cut_range else ''
        
        # =====================================================================
        ## 2) check that all variables are present in the dataset 
        assert all ( ( v in dataset ) for v in vnames ) , 'Not all variables are in dataset!'

        # =====================================================================
        ## 3) reduce dataset if only a small subset of variables is requested 
        nvars = len ( dataset.get() )
        if 2 * len ( vnames )  <= nvars and not more_vars  :
            dstmp  = dataset.subset ( vnames )
            result = ds2numpy ( dstmp , vnames , cuts = cuts , cut_range = cut_range )
            ## dstmp.erase()
            dstmp = Ostap.MoreRooFit.delete_data ( dstmp ) 
            del dstmp 
            return result
        
        # =========================================================================
        ## 4) if cuts or cut-range is specified, assume cuts are hash and make a filtering 
        if cuts or cut_range :
            with useStorage ( ROOT.RooAbsData.Vector ) :
                dstmp = dataset.subset ( vnames if not more_vars else [] ,  cuts = cuts , cut_range = cut_range )
            result = ds2numpy ( dstmp , vnames , more_vars = more_vars )
            ## dstmp.erase()
            dstmp = Ostap.MoreRooFit.delete_data ( dstmp ) 
            del dstmp 
            return result
        
        # =======================================================
        funcs    = []
        formulas = []
        varlst   = dataset.varlst () 
        if more_vars and isinstance ( more_vars , dictlike_types ) :
            for name , fun in loop_items ( more_vars ) :
                if   isinstance ( fun , AFUN1           ) : absreal = fun.fun
                elif isinstance ( fun , ROOT.RooAbsPdf  ) : absreal = fun
                elif isinstance ( fun , ROOT.RooAbsReal ) : absreal = fun
                elif isinstance ( fun , string_types    ) :
                    formula = Ostap.FormulaVar ( name , fun , varlst )
                    formulas.append ( formula ) 
                    absreal = formula
                else :
                    raise TypeError ( "Invalid type of fun/pdf" )
                obsvars = absreal.getObservables ( dataset )
                item    = name , absreal , obsvars
                funcs.append ( item )
        elif more_vars :
            for item in more_vars :
                
                len2 = isinstance ( item , sized_types ) and 2 == len ( item )
                name = '' 
                if   isinstance ( item , AFUN1           ) : absreal = item.fun
                elif isinstance ( item , ROOT.RooAbsPdf  ) : absreal = item
                elif isinstance ( item , ROOT.RooAbsReal ) : absreal = item
                elif len2 and isinstance ( item [ 0 ] , string_types ) and isinstance ( item [1] , AFUN1          ) :
                    absreal, name = item[1].fun , item[0]
                elif len2 and isinstance ( item [ 0 ] , string_types ) and isinstance ( item [1] , ROOT.RooAbsPdf ) :
                    absreal, name = item[1]     , item[0]
                elif len2 and isinstance ( item [ 0 ] , string_types ) and isinstance ( item [1] , ROOT.RooAbsPdf ) :
                    absreal, name = item[1]     , item[0]
                elif len2 and isinstance ( item [ 0 ] , string_types ) and isinstance ( item [1] , string_types   ) :
                    formula = Ostap.FormulaVar ( item [ 0 ] , item [ 1 ] , varlst )
                    formulas.append ( formula ) 
                    absreal, name = formula , item [ 0 ]                    
                else :
                    raise TypeError ( "Invalid type of fun/pdf!" )                
                obsvars = absreal.getObservables ( dataset )
                item    = name if name else absreal.name , absreal , obsvars
                funcs.append ( item )
        
        # =====================================================================
        for item in funcs :
            name , _ , _ = item 
            assert not name in dataset, 'da2numpy: no way to redefine variable `%s`!' % name 

        # =====================================================================
        vars       = dataset.get()
        vars       = [ v for v in vars if v.name in vnames ]
        doubles    = [ v.name for v in vars if isinstance ( v , ROOT.RooAbsReal     ) ] 
        categories = [ v.name for v in vars if isinstance ( v , ROOT.RooAbsCategory ) ]

        ## name of weight variable 
        weighted = dataset.isWeighted ()
        weight   = '' if not weighted else dataset.wname () 

        dtypes = [] 
        for v in vnames :
            if   v in doubles    : dtypes.append ( ( v      , np.float64 ) ) 
            elif v in categories : dtypes.append ( ( v      , np.int64   ) )
        for f in funcs           : dtypes.append ( ( f[0]   , np.float64 ) ) 
        if weight                : dtypes.append ( ( weight , np.float64 ) ) 
                    
        ## create data 
        data = np.zeros ( len ( dataset )  , dtype = dtypes )
    
        ## make an explict loop:
        for i , item in enumerate ( progress_bar ( dataset , silent = silent ) ) :

            evt, the_weight = item
            
            for v in evt :
                vname = v.name
                if   vname in doubles    : data [ vname  ] [ i ] = float ( v  )
                elif vname in categories : data [ vname  ] [ i ] = int   ( v  )

            if weighted and weight  and not ( the_weight is None ) :
                data [ weight ] [ i ] = float ( the_weight ) 

            ## add PDF values
            for vname , func , obsvars in funcs :
                obsvars.assign ( evt )
                data [ vname ] [ i ] = func.getVal()   

        del funcs
        del formulas 
        return data
    
    # =========================================================================
    __all__  = __all__ + ( 'ds2numpy' , )
    ROOT.RooDataSet.tonumpy = ds2numpy 
    ROOT.RooDataSet.tonp    = ds2numpy
    ROOT.RooDataSet.to_np   = ds2numpy
    _new_methods_ += [ ROOT.RooDataSet.tonumpy ,
                       ROOT.RooDataSet.tonp    ,
                       ROOT.RooDataSet.to_np   ]
    

# =============================================================================
else    : # ===================================================================
# =============================================================================

    logger.warning ( "Numpy is not available: no action" )


# ==============================================================================
if np : # ======================================================================
    # ==========================================================================
    ## Get the dict of empirical cumulative distribution functions from dataset
    #  @code
    #  dataset = ...
    #  ecdfs   = dataset.ecdfs ( 'a,b,c' , 'pt>10' ) 
    #  @endcode
    #  @see Ostap::Math::ECDF
    def ds2cdfs ( dataset           ,
                  variables         ,
                  cuts       = ''   ,
                  cut_range  = ''   ,
                  more_vars  = {}   , 
                  silent     = True ) :
        """ Get the dict of empirical cumulative distribution functions from dataset
        - see `Ostap.Math.ECDF`
        >>> dataset = ...
        >>> ecdfs   = dataset.ecdfs ( 'a,b,c' , 'pt>10' ) 
        """

        assert not more_vars or isinstance ( more_vars , dictlike_types ) , \
            "ds2cdfs: invalid type of `more_vars`" % type ( more_vars )
        
        ## decode the list of variables 
        varlst, cuts ,  _ = vars_and_cuts ( variables , cuts )
        
        extra  = [ v for v in varlst if not v in dataset ]
        assert varlst and not extra , 'Variables are not in dataset: %s' % str ( extra ) 

        if dataset.isWeighted() :
            logger.warning ( 'ds2cdfs: dataset is weighted! Weight will be ignored!')

        ## 1) get data as numpy-array 
        data = ds2numpy ( dataset               ,
                          varlst                ,
                          cuts      = cuts      ,
                          cut_range = cut_range ,
                          more_vars = more_vars , 
                          silent    = silent    )
        
        result = {}
        for vname in varlst    : result [ vname ] =  Ostap.Math.ECDF ( data [ vname ] )
        for vname in more_vars : result [ vname ] =  Ostap.Math.ECDF ( data [ vname ] )
            
        del data
        return result 

    __all__  = __all__ + ( 'ds2cdfs' , )
    ROOT.RooDataSet.cdfs  = ds2cdfs
    ROOT.RooDataSet.ecdfs = ds2cdfs
    _new_methods_ += [ ROOT.RooDataSet.cdfs  ,
                       ROOT.RooDataSet.ecdfs ]
        
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
