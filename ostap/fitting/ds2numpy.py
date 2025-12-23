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
    'ds2numpy' , ## conversion of RoDataSet into numpy array
    'ds2cdfs'  , ## get empirical cumulative distributions from RooDataSet
)
# =============================================================================
from   ostap.core.meta_info         import root_info
from   ostap.core.ostap_types       import string_types, dictlike_types, sized_types
from   ostap.core.core              import Ostap
from   ostap.utils.basic            import loop_items, typename  
from   ostap.utils.utils            import split_range
from   ostap.math.base              import doubles
from   ostap.fitting.dataset        import useStorage
from   ostap.fitting.funbasic       import AFUN1 
from   ostap.utils.progress_bar     import progress_bar
from   ostap.trees.cuts             import vars_and_cuts 
import ostap.fitting.roocollections
import ROOT, numpy 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.ds2numpy' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug( 'Convert RooDataSet into numpy array')
# =============================================================================
_new_methods_ = []
# =============================================================================
if ( 6 , 28 ) <= root_info  :  ## 6.26 <= ROOT     
    # =========================================================================
    ## Convert dataset into numpy array using <code>ROOT.RooAbsData</code> interface 
    #  @see ROOT.RooAbsData.getBatches
    #  @see ROOT.RooAbsData.getCategoryBatches
    #  @see ROOT.RooAbsData.getWeightBatche    
    #  @see ROOT.RooAbsDataStore.getBatches
    #  @see ROOT.RooAbsDataStore.getCategoryBatches
    #  @see ROOT.RooAbsDataStore.getWeightBatche
    #  @attention conversion to ROOT.RooVectorDataStore is used!
    #
    #  Unlike <code>to_numpy</code> method it allows more flexible outp
    #   - structured array (default) vs unstructured array
    #   - split data and weight columns     
    def ds2numpy ( dataset                , 
                   var_lst                , * , 
                   cuts        = ''       ,
                   cut_range   = ''       , 
                   silent      = True     ,
                   more_vars   = {}       ,
                   structured  = True     ,
                   weight_name = 'weight' , 
                   weight_split = False   ) :
        """ Convert dataset into numpy array using `ROOT.RooAbsData` iterface 
        - see ROOT.RooAbsData.getBatches
        - see ROOT.RooAbsData.getCategoryBatches
        - see ROOT.RooAbsData.getWeightBatche    
        - see ROOT.RooAbsDataStore.getBatches
        - see ROOT.RooAbsDataStore.getCategoryBatches
        - see ROOT.RooAbsDataStore.getWeightBatche    
        - attention: Conversion to `ROOT.RooVectorDataStore` is used! 

        Unlike `ROOT.RooDataSet.to_numpy` method it allows more flexible outp
        - structured array (default) vs unstructured array
        - optional split data and weight columns     

        """
        
        if isinstance ( var_lst , string_types ) : var_lst = [ var_lst ]
        
        assert not more_vars or isinstance ( more_vars , dictlike_types ) , \
            "ds2numpy: invalid type of `more_vars`" % typename ( more_vars )
        
        # =====================================================================
        ## (1) get names of all requested variables
        if   all ( isinstance ( v , string_types   ) for v in var_lst ) :
            vnames , cuts , _ = vars_and_cuts ( var_lst , cuts , allow_empty = more_vars )
        elif all ( isinstance ( v , ROOT.RooAbsArg ) for v in var_lst ) :
            vnames = [ v.GetName() for v in var_lst ]
        else :
            raise TypeError ( "Invalid type of `var_list`!" ) 

        cuts      = str ( cuts      ).strip () 
        cut_range = str ( cut_range ).strip () if cut_range else ''

        assert vnames or more_vars , "No variables are specified! "
        
        # =====================================================================
        ## (2) check that all variables are present in the dataset 
        assert all ( ( v in dataset ) for v in vnames ) , 'Not all variables are in dataset!'

        # =====================================================================
        ## (3) reduce dataset if only a small subset of variables is requested 
        nvars = len ( dataset.get() )
        if 2 * len ( vnames )  <= nvars and not more_vars  :            
            with useStorage ( ROOT.RooAbsData.Vector ) : 
                dstmp  = dataset.subset ( vnames , cuts = cuts , cut_range = cut_range )
            result = ds2numpy ( dstmp     ,
                                vnames    ,
                                more_vars    = more_vars    ,
                                silent       = silent       ,
                                structured   = structured   , 
                                weight_name  = weight_name  ,
                                weight_split = weight_split )
            ## dstmp.erase()
            dstmp = Ostap.MoreRooFit.delete_data ( dstmp ) 
            del dstmp 
            return result

        # =========================================================================
        ## (4) if cuts or cut-range is specified, assume cuts are harsh and make a filtering 
        if cuts or cut_range :
            with useStorage ( ROOT.RooAbsData.Vector ) :
                dstmp = dataset.subset ( vnames if not more_vars else [] , cuts = cuts , cut_range = cut_range )
            result = ds2numpy ( dstmp        ,
                                vnames       ,
                                more_vars    = more_vars    ,
                                silent       = silent       ,
                                structured   = structured   ,
                                weight_name  = weight_name  ,
                                weight_split = weight_split )
            ## dstmp.erase()
            dstmp = Ostap.MoreRooFit.delete_data ( dstmp ) 
            del dstmp 
            return result

        weighted          = dataset.isWeighted ()
        store_errors      = weighted and dataset.store_errors      ()
        store_asym_errors = weighted and dataset.store_asym_errors ()         
        if store_errors or store_asym_errors :
            logger.warning ("Weight uncertainties are defined, but will be ignored!") 
  
        funcs    = []
        formulas = []
        varlst   = dataset.varlst () 
        if more_vars and isinstance ( more_vars , dictlike_types ) :
            for name , fun in loop_items ( more_vars ) :
                assert not name in dataset, 'ds2numpy: no way to redefine variable `%s`!' % name 
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
            assert not name in dataset, 'ds2numpy: no way to redefine variable `%s`!' % name 

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
        weight = '' if not dataset.isWeighted() else str ( dataset.wname () ) 

        dtypes = [] 
        for v in vnames :
            if   v in doubles    : dtypes.append ( ( str ( v      ) , numpy.float64 ) ) 
            elif v in categories : dtypes.append ( ( str ( v      ) , numpy.int64   ) )

        for f in funcs           : dtypes.append ( ( str ( f[0]   ) , numpy.float64 ) ) 
        if weight                : dtypes.append ( ( str ( weight ) , numpy.float64 ) ) 
            
        ## get data in batches
        nevts  = len ( dataset ) 

        data   = None

        ## maximal size of data chunk 
        nmax   = max ( nevts // 6 , 30000 // nvars )
        
        ## get data in chunks/batches 
        for first, last in progress_bar ( split_range ( 0 , nevts , nmax ) , silent = silent , description = 'Chunks:' ) :
            
            num   = last - first            
            wget  = False 
            part  = numpy.zeros ( num , dtype = dtypes )

            if doubles :
                dpart   = source.getBatches ( first , num )
                for d in dpart :
                    dname   = d.first.name
                    buffer  = d.second
                    payload = numpy.frombuffer ( buffer.data() , dtype = numpy.float64 , count = buffer.size() )
                    if   dname in doubles : part [ dname ] = payload
                    elif dname == weight  : part [ dname ] = payload
                del dpart
                
            if categories :
                cpart   = source.getCategoryBatches ( first , num )
                for c in cpart :
                    cname   = c.first.name
                    buffer  = c.second 
                    payload = numpy.frombuffer ( buffer.data() , dtype = numpy.int64 , count = buffer.size() )
                    if cname in categories : part [ cname ] = payload 
                del cpart
                
            if weight and not wget :
                if twoargs : buffer = source.getWeightBatch ( first , num         )
                else       : buffer = source.getWeightBatch ( first , num , False )
                weights = numpy.frombuffer ( buffer.data(), dtype = numpy.float64, count = buffer.size() )
                if buffer  : part [ weight ] = weights 
                else       : part [ weight ] = numpy.full ( num , source.weight() , dtype = numpy.float64 )

            if data is None : data = part
            else            :  
                data = numpy.concatenate ( [ data , part ] )
                del part

        ## add PDF values (explicit loop)
        if funcs :  
            for i, item  in progress_bar ( enumerate ( source ) , silent = silent , description = 'Entries:' ) : 
                entry , _ = item 
                for vname , func , obsvars in funcs :
                    obsvars.assign ( entry )
                    data [ vname ] [ i ] = func.getVal()   
                    
        if delsource :
            source.reset()
            del source

        if not cuts and not cut_range :
            assert len ( data ) == len ( dataset ) , "Mismatch in input/output lengths!"
        else :
            assert len ( data ) <= len ( dataset ) , "Mismatch in input/output lengths!"

        del    funcs
        del    formulas
        
        # =======================================================================
        ## The final actions:
        # =======================================================================

        ## (1) rename weight column 
        if weight and weight_name and weight_name != weight and not weight_name in data.dtype.names : 
            from numpy.lib import recfunctions as rfn
            data   = rfn.rename_fields ( data , { weight : weight_name } ) 
            weight = weight_name

        ## (2) split `weight` column from `data` columns 
        weights  = None        
        if weight and weight_split :
            weights = data [ weight ]
            data    = data [ [ v for v in data.dtype.names if v != weight ] ]

        ## (3) convert to unstructured array, if requested 
        if not structured :
            from numpy.lib.recfunctions import structured_to_unstructured as s2u
            data = s2u ( data , copy = False )
        
        ## 
        return ( data , weights ) if weight_split else data          ## RETURN

# =============================================================================
else : 
    # =========================================================================
    ## Convert dataset into numpy array using (slow) explicit loops 
    #  Unlike <code>to_numpy</code> method it allows more flexible outp
    #   - structured array (default) vs unstructured array
    #   - split data and weight columns     
    def ds2numpy ( dataset                 ,
                   var_lst                 ,
                   cuts         = ''       , * , 
                   cut_range    = ''       ,  
                   silent       = False    ,
                   more_vars    = {}       ,
                   structured   = True     ,
                   weight_name  = 'weight' , 
                   weight_split = False    ) :
        """ Convert dataset into numpy array using (slow) explicit loops
        Unlike `ROOT.RooDataSet.to_numpy` method it allows more flexible outp
        - structured array (default) vs unstructured array
        - optional split data and weight columns             
        """
        
        if isinstance ( var_lst , string_types ) : var_lst = [ var_lst ]
        
        assert not more_vars or isinstance ( more_vars , dictlike_types ) , \
            "ds2numpy: invalid type of `more_vars`" % typename ( more_vars )

        # =====================================================================
        ## 1) check that all variables are present in dataset 
        if   all ( isinstance ( v , string_types   ) for v in var_lst ) :
            vnames , cuts , _  = vars_and_cuts ( var_lst , cuts , allow_empty = more_vars )
        elif all ( isinstance ( v , ROOT.RooAbsArg ) for v in var_lst ) :
            vnames = [ v.GetName() for v in var_lst ]
        else :
            raise TypeError ( "Invalid type of `var_list`!" ) 

        cuts      = str ( cuts      ).strip () 
        cut_range = str ( cut_range ).strip () if cut_range else ''
        
        assert vnames or more_vars , "No variables are specified! "
        
        # =====================================================================
        ## 2) check that all variables are present in the dataset 
        assert all ( ( v in dataset ) for v in vnames ) , 'Not all variables are in dataset!'

        # =====================================================================
        ## 3) reduce dataset if only a small subset of variables is requested 
        nvars = len ( dataset.get() )
        if 2 * len ( vnames )  <= nvars and not more_vars  :
            dstmp  = dataset.subset ( vnames )
            result = ds2numpy ( dstmp                       ,
                                vnames                      ,
                                cuts         = cuts         ,
                                cut_range    = cut_range    ,
                                structured   = structured   ,
                                weight_name  = weight_name  ,
                                weight_split = weight_split )
            ## dstmp.erase()
            ## dstmp = Ostap.MoreRooFit.delete_data ( dstmp )
            dstmp.clear() 
            del dstmp 
            return result
        
        # =========================================================================
        ## 4) if cuts or cut-range is specified, assume cuts are hash and make a filtering 
        if cuts or cut_range :
            with useStorage ( ROOT.RooAbsData.Vector ) :
                dstmp = dataset.subset ( vnames if not more_vars else [] ,  cuts = cuts , cut_range = cut_range )
            result = ds2numpy ( dstmp        ,
                                vnames       ,
                                more_vars    = more_vars    ,
                                structured   = structured   ,
                                weight_name  = weight_name  ,
                                weight_split = weight_split )
            ## dstmp.erase()
            ## dstmp = Ostap.MoreRooFit.delete_data ( dstmp ) 
            dstmp.clear() 
            del dstmp 
            return result
        
        # =======================================================
        weighted          = dataset.isWeighted ()
        store_errors      = weighted and dataset.store_errors      ()
        store_asym_errors = weighted and dataset.store_asym_errors ()         
        if store_errors or store_asym_errors :
            logger.warning ( "Weight uncertainties are defined, but will be ignored!" ) 
  
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
        weight   = '' if not weighted else str ( dataset.wname () ) 

        dtypes = [] 
        for v in vnames :
            if   v in doubles    : dtypes.append ( ( str ( v      ) , numpy.float64 ) ) 
            elif v in categories : dtypes.append ( ( str ( v      ) , numpy.int64   ) )
            
        for f in funcs           : dtypes.append ( ( str ( f[0]   ) , numpy.float64 ) ) 
        if weight                : dtypes.append ( ( str ( weight ) , numpy.float64 ) ) 
                    
        ## create data 
        data = numpy.zeros ( len ( dataset )  , dtype = dtypes )
    
        ## make an explict loop:
        for i , item in enumerate ( progress_bar ( dataset , silent = silent , description = 'Entries:' ) ) :

            evt, the_weight = item
            
            for v in evt :
                vname = str ( v.name ) 
                if   vname in doubles    : data [ vname  ] [ i ] = float ( v  )
                elif vname in categories : data [ vname  ] [ i ] = int   ( v  )

            if weighted and weight  and not ( the_weight is None ) :
                data [ weight ] [ i ] = float ( the_weight ) 

            ## add PDF values
            for vname , func , obsvars in funcs :
                obsvars.assign ( evt )
                data [ vname ] [ i ] = func.getVal()   

        if not cuts and not cut_range :
            assert len ( data ) == len ( dataset ) , "Mismatch in input/output lengths!"
        else :
            assert len ( data ) <= len ( dataset ) , "Mismatch in input/output lengths!"
            
        del funcs
        del formulas

       # =======================================================================
        ## The final actions:
        # =======================================================================

        ## (1) rename weight column 
        if weight and weight_name and weight_name != weight and not weight_name in data.dtype.names : 
            from numpy.lib import recfunctions as rfn
            data   = rfn.rename_fields ( data , { weight : weight_name } ) 
            weight = weight_name

        ## (2) split `weight` column from `data` columns 
        weights  = None        
        if weight and weight_split :
            weights = data [ weight ]
            data    = data [ [ v for v in data.dtype.names if v != weight ] ]

        ## (3) convert to unstructured array, if requested 
        if not structured :
            from numpy.lib.recfunctions import structured_to_unstructured as s2u
            data = s2u ( data , copy = False )
        
        ## 
        return ( data , weights ) if weight_split else data          ## RETURN

# =========================================================================
ROOT.RooDataSet.tonumpy = ds2numpy 
ROOT.RooDataSet.tonp    = ds2numpy
ROOT.RooDataSet.to_np   = ds2numpy
_new_methods_ += [ ROOT.RooDataSet.tonumpy ,
                   ROOT.RooDataSet.tonp    ,
                   ROOT.RooDataSet.to_np   ]
    
# ==============================================================================
if  ( 6 , 32 ) <= root_info : data2vct = lambda s : s
else                        : data2vct = lambda s : doubles ( s ) 
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
        "ds2cdfs: invalid type of `more_vars`" % typename ( more_vars )
    
    ## decode the list of variables 
    varlst, cuts ,  _ = vars_and_cuts ( variables , cuts , allow_empty = more_vars )
        
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
    for vname in varlst    : result [ vname ] = Ostap.Math.ECDF ( data2vct ( data [ vname ] ) ) 
    for vname in more_vars : result [ vname ] = Ostap.Math.ECDF ( data2cvt ( data [ vname ] ) ) 
    
    del data
    return result 



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
