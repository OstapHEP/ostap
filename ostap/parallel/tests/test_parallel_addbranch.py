#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/trees/tests/test_parallel_addbranch.py
# Four ways to add branch into TTree/Tchain
# - using string formula (TTreeFormula-based)
# - using pure python function
# - using histogram/function
# - using histogram sampling
#
# There are more specific similar actions
# - adding to TTree/TChain TMVA decision
# - adding to TTree/TChain TMVA/Choppnig decision
# - adding to TTree/TChain reweigting information
# Copyright (c) Ostap developers.
# ============================================================================= 
""" Test module for ostap/trees/cuts.py.
Four ways to add branch into TTree/Tchain
- using string formula (TTreeFormula-based)
- using pure python function
- using histogram/function
- using histogram sampling
"""
# ============================================================================= 
from   ostap.core.meta_info               import root_info
import ostap.trees.trees
import ostap.histos.histos
from   ostap.core.pyrouts                 import hID , Ostap 
from   ostap.trees.data                   import Data
from   ostap.math.make_fun                import make_fun1, make_fun2, make_fun3 
from   ostap.utils.timing                 import timing 
from   ostap.utils.progress_bar           import progress_bar
from   ostap.parallel.parallel            import pickles 
import ostap.parallel.parallel_add_branch 
import ROOT, math, random, array  
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_parallel_addbranch' )
else                       : logger = getLogger ( __name__                  )
# =============================================================================

# =============================================================================
## create a file with tree 
def create_tree ( fname , nentries = 1000 ) :
    """Create a file with a tree
    >>> create_tree ( 'file.root' ,  1000 ) 
    """
    
    import ROOT, random 
    import ostap.io.root_file
    
    from array import array 
    var1 = array ( 'd', [ 0 ] )
    var2 = array ( 'd', [ 0 ] )
    var3 = array ( 'd', [ 0 ] )
    
    from ostap.core.core import ROOTCWD

    with ROOTCWD() , ROOT.TFile.Open( fname , 'new' ) as root_file:
        root_file.cd () 
        tree = ROOT.TTree ( 'S','tree' )
        tree.SetDirectory ( root_file  ) 
        tree.Branch ( 'mass'  , var1 , 'mass/D'  )
        tree.Branch ( 'c2dtf' , var2 , 'c2dtf/D' )
        tree.Branch ( 'pt'    , var3 , 'pt/D'    )
        
        for i in range ( nentries ) : 
            
            m  = random.gauss        ( 3.1 ,  0.015 )
            c2 = random.gammavariate ( 2.5 , 0.5    ) / 5 
            pt = random.uniform      ( 0   , 10     )
            
            var1[0] = m
            var2[0] = c2 
            var3[0] = pt
            
            tree.Fill()
            
        root_file.Write()
        
# =============================================================================
def prepare_data ( nfiles = 50 ,  nentries = 500  ) :

    from ostap.utils.cleanup import CleanUp    
    files = [ CleanUp.tempfile ( prefix = 'ostap-test-trees-addbranch-%d-' % i ,
                                 suffix = '.root' ) for i in range ( nfiles)  ]
    
    for f in progress_bar ( files ) : create_tree ( f , nentries )
    return files

# =============================================================================

## top-level function 
def fun_ftwo ( x ) : return 2 * x

            
# =============================================================================
# Many ways to add branch into TTree/Tchain
# - using string formula (TTreeFormula-based)
# - using pure python function
# - using histogram/function
# - using histogram sampling
def test_addbranch() :
    """Many ways to add branch into TTree/Tchain
    - using string formula (TTreeFormula-based)
    - using pure python function
    - using histogram/function
       - using histogram sampling
    """
    
    ## files = prepare_data ( 100 , 1000 )
    files = prepare_data ( 10 , 1000 )
    
    logger.info ( '#files:    %s'  % len ( files ) )  
    data = Data ( 'S' , files )
    logger.info ( 'Initial Tree/Chain:\n%s' % data.chain.table ( prefix = '# ' ) )

    # =========================================================================
    ## 1) add new branch as TTree-formula:
    # =========================================================================
    with timing ('expression' , logger = logger ) :          
        chain = data.chain 
        chain.padd_new_branch ( 'et','sqrt(pt*pt+mass*mass)' )        
    ## reload the chain and check: 
    logger.info ( 'With formula:\n%s' % data.chain.table ( prefix = '# ' ) )
    assert 'et' in data.chain , "Branch ``et'' is  not here!"

    # =========================================================================
    ## 2) add several new branches as TTree-formula:
    # =========================================================================
    with timing ('simultaneous' , logger = logger ) :          
        chain = data.chain  
        chain.padd_new_branch ( { 'Et1' : 'sqrt(pt*pt+mass*mass)'   ,
                                  'Et2' : 'sqrt(pt*pt+mass*mass)*2' ,
                                  'Et3' : 'sqrt(pt*pt+mass*mass)*3' } , None )        
    ## reload the chain and check: 
    logger.info ( 'With formula:\n%s' % data.chain.table ( prefix = '# ' ) )
    assert 'Et1' in data.chain , "Branch `Et1' is  not here!"
    assert 'Et2' in data.chain , "Branch `Et2' is  not here!"
    assert 'Et3' in data.chain , "Branch `Et3' is  not here!"

    # =========================================================================
    ## 2) add new branch as pure python function 
    # =========================================================================
    with timing ( 'pyfunc' , logger = logger ) :
        et2 = lambda tree : tree.pt**2 + tree.mass**2        
        chain = data.chain
        if pickles ( et2 ) :
            chain.padd_new_branch ( 'et2', et2 )
        else               :
            logger.warning ( "pyfunc: switch to sequential processing(lambda cannot be pickled" )
            chain. add_new_branch ( 'et2', et2 )            
    ## reload the chain and check: 
    logger.info ( 'With python:\n%s' % data.chain.table ( prefix = '# ' ) )
    assert 'et2' in data.chain , "Branch `et2' is  not here!"

    # =========================================================================
    ## 3) add new branch as histogram-function 
    # =========================================================================
    with timing ('histo-1' , logger = logger ) :          
        h1  = ROOT.TH1D ( hID () , 'some pt-correction' , 100 , 0 , 10 )
        h1 += lambda x :  1.0 + math.tanh( 0.2* ( x - 5 ) )         
        from   ostap.trees.funcs  import FuncTH1
        ptw   = FuncTH1 ( h1 , 'pt' )
        chain = data.chain 
        if pickles ( ptw ) : 
            chain.padd_new_branch ( 'ptw', ptw )
        else :
            logger.warning ( "histo-1: switch to sequential processing (object cannot be pickled)" )
            chain.add_new_branch ( 'ptw', ptw )            
        ## reload the chain and check: 
        logger.info ( 'With histogram:\n%s' % data.chain.table ( prefix = '# ' ) )
        assert 'ptw' in data.chain , "Branch `ptw' is  not here!"
            
    # =========================================================================
    ## 4) add several functions simultanepusly 
    # =========================================================================
    with timing ('histo-3' , logger = logger ) :          
        from   ostap.trees.funcs  import FuncTH1
        hh    = ROOT.TH1D ( hID() , 'some pt-correction' , 100 , 0 , 10 )
        h1    = hh + ( lambda x :  1.0 + math.tanh ( 0.1 * ( x - 5 ) ) )
        ptw1  = FuncTH1 ( h1 , 'pt' )
        h2    = h1 + ( lambda x :  1.0 + math.tanh ( 0.2 * ( x - 5 ) ) )
        ptw2  = FuncTH1 ( h2 , 'pt' )
        h3    = h1 + ( lambda x :  1.0 + math.tanh ( 0.3 * ( x - 5 ) ) ) 
        ptw3  = FuncTH1 ( h3 , 'pt' )
        chain = data.chain
        brs = { 'ptw1' : ptw1 , 'ptw2' : ptw2 , 'ptw3' : ptw1 }
        if pickles ( brs ) :             
            chain.padd_new_branch ( None , brs )
        else :
            logger.warning ( "histo-3: switch to sequential processing (object cannot be pickled)" )
            chain. add_new_branch ( None , brs )            
        ## reload the chain and check: 
        logger.info ( 'With histogram:\n%s' % data.chain.table ( prefix = '# ' ) )
        assert 'ptw' in data.chain , "Branch `ptw1' is  not here!"
        assert 'ptw' in data.chain , "Branch `ptw2' is  not here!"
        assert 'ptw' in data.chain , "Branch `ptw3' is  not here!"
            
    # =========================================================================
    ## 5) add the variable sampled from the histogram
    # =========================================================================
    with timing ('histo-2' , logger = logger ) :          
        h2 = ROOT.TH1D( hID() , 'Gauss' , 120 , -6 , 6 )
        for i in range ( 100000 ) : h2.Fill ( random.gauss ( 0 , 1 ) )
        chain = data.chain 
        if pickles ( h2 ) : 
            chain.padd_new_branch ( 'hg', h2 )
        else :
            logger.warning ( "histo-2: switch to sequential processing (object cannot be pickled)" )
            chain. add_new_branch ( 'hg', h2 )            
        ## reload the chain and check: 
        logger.info ( 'With sampled:\n%s' % data.chain.table ( prefix = '# ' ) )
        assert 'hg' in data.chain , "Branch `hg' is  not here!"
            
    # =========================================================================
    ## 6) python function again 
    # =========================================================================
    with timing ('gauss' , logger = logger ) :          
        def gauss ( *_ ) : return random.gauss(0,1)    
        chain = data.chain
        if pickles ( gauss ) : 
            chain.padd_new_branch ( 'gauss', gauss )
        else :
            logger.warning ( "gauss: switch to sequential processing (object cannot be pickled)" )
            chain. add_new_branch ( 'gauss', gauss )            
        ## reload the chain and check: 
        logger.info ( 'With gauss:\n%s' % data.chain.table ( prefix = '# ' ) )
        assert 'gauss' in data.chain , "Branch `gauss' is  not here!"
        
    # =========================================================================
    ## 7) add numpy array 
    # =========================================================================
    try : 
        import numpy
    except ImportError :
        numpy  = None

    if numpy : ## ATTENTION! 

        with timing ('numpy float16' , logger = logger ) :
           adata  = numpy.full ( 10000 , +0.1 , dtype = numpy.float16 )
           chain  = data.chain
           chain.padd_new_branch ( 'np_f16' , adata )
        ## reload the chain and check: 
        logger.info ( 'With numpy.float16:\n%s' % data.chain.table ( prefix = '# ' ) )
        assert 'np_f16' in data.chain , "Branch `np_f16' is  not here!"

        with timing ('numpy float32' , logger = logger ) :
           adata = numpy.full ( 10000 , -0.2 , dtype = numpy.float32 )
           chain = data.chain
           chain.padd_new_branch ( 'np_f32' , adata )            
        ## reload the chain and check: 
        logger.info ( 'With numpy.float32:\n%s' % data.chain.table ( prefix = '# ' ) )
        assert 'np_f32' in data.chain , "Branch `np_f32' is  not here!"
        
        with timing ('numpy float64' , logger = logger ) :
            adata  = numpy.full ( 10000 , +0.3 , dtype = numpy.float64 )
            chain  = data.chain            
            chain.padd_new_branch ( 'np_f64' , adata )
        ## reload the chain and check: 
        logger.info ( 'With numpy.float64:\n%s' % data.chain.table ( prefix = '# ' ) )
        assert 'np_f64' in data.chain , "Branch `np_f64' is  not here!"

        ## with timing ('numpy int8 ' , logger = logger ) :  
        ##     adata  = numpy.full ( 10000 , -1 , dtype = numpy.int8 )
        ##     chain  = data.chain            
        ##     chain.padd_new_branch ( 'np_i8' , adata )
        ## ## reload the chain and pickles: 
        ## logger.info ( 'With numpy.int8:\n%s' % data.chain.table ( prefix = '# ' ) )
        ## assert 'np_i8' in data.chain , "Branch `np_i8' is  not here!"

        ## with timing ('numpy uint8 ' , logger = logger ) :  
        ##   adata  = numpy.full ( 10000 , +2 , dtype = numpy.uint8 )
        ##   chain  = data.chain            
        ##   chain.padd_new_branch ( 'np_ui8' , adata )
        ## ## reload the chain and pickles: 
        ## logger.info ( 'With numpy.uint8:\n%s' % data.chain.table ( prefix = '# ' ) )
        ## assert 'np_ui8' in data.chain , "Branch `np_ui8' is  not here!"

        with timing ('numpy int16 ' , logger = logger ) :  
           adata  = numpy.full ( 10000 , -3 , dtype = numpy.int16 )
           chain  = data.chain            
           chain.padd_new_branch ( 'np_i16' , adata )
        ## reload the chain and check: 
        logger.info ( 'With numpy.int16:\n%s' % data.chain.table ( prefix = '# ' ) )
        assert 'np_i16' in data.chain , "Branch `np_i16' is  not here!"

        with timing ('numpy uint16 ' , logger = logger ) :  
           adata  = numpy.full ( 10000 , +4 , dtype = numpy.uint16 )
           chain  = data.chain            
           chain.padd_new_branch ( 'np_ui16' , adata )
        ## reload the chain and check: 
        logger.info ( 'With numpy.uint16:\n%s' % data.chain.table ( prefix = '# ' ) )
        assert 'np_ui16' in data.chain , "Branch `np_ui16' is  not here!"

        with timing ('numpy int32 ' , logger = logger ) :  
           adata  = numpy.full( 10000 , -5 , dtype = numpy.int32 )
           chain  = data.chain            
           chain.padd_new_branch ( 'np_i32' , adata )
        ## reload the chain and pickles: 
        logger.info ( 'With numpy.int32:\n%s' % data.chain.table ( prefix = '# ' ) )
        assert 'np_i32' in data.chain , "Branch `np_i32' is  not here!"

        with timing ('numpy uint32 ' , logger = logger ) :  
            adata  = numpy.full ( 10000 , +6 , dtype = numpy.uint32 )
            chain  = data.chain            
            chain.padd_new_branch ( 'np_ui32' , adata )
        ## reload the chain and check: 
        logger.info ( 'With numpy.uint32:\n%s' % data.chain.table ( prefix = '# ' ) )
        assert 'np_ui32' in data.chain , "Branch `np_ui32' is  not here!"
        
        with timing ('numpy int64 ' , logger = logger ) :  
            adata  = numpy.full ( 10000 , -7 , dtype = numpy.int64 )
            chain  = data.chain
            chain.padd_new_branch ( 'np_i64' , adata )
        ## reload the chain and check: 
        logger.info ( 'With numpy.int64:\n%s' % data.chain.table ( prefix = '# ' ) )
        assert 'np_i64' in data.chain , "Branch `np_i64' is  not here!"

        with timing ('numpy uint64 ' , logger = logger ) :  
            adata  = numpy.full ( 10000 , +8 , dtype = numpy.uint64 )
            chain  = data.chain
            chain.padd_new_branch ( 'np_ui64' , adata )
        ## reload the chain and check: 
        logger.info ( 'With numpy.uint64:\n%s' % data.chain.table ( prefix = '# ' ) )
        assert 'np_ui64' in data.chain , "Branch `np_ui64' is  not here!"


    for l,v in ( ('f', +100.1 )  ,
                 ('d', -200.2 ) ,
                 ('i',-3) , ('l',-4) ,
                 ('I',5)  ,  ('L',6) ,
                 ('h',7)  ,  ('H',8) ) :

        with timing ('array %s'% l , logger = logger ) :  
            adata  = array.array ( l ,  10000*[ v ] ) 
            chain  = data.chain
            vname  = 'arr_%s' % l 
            chain.padd_new_branch ( vname , adata )
            ## reload the chain : 
            logger.info ( "With array '%s':\n%s" % ( l ,  data.chain.table ( prefix = '# ' ) ) ) 
            assert vname in data.chain , "Branch `%s' is  not here!" % vname 

    
    ## add function
    with timing ('1D-function ' , logger = logger ) :

        logger.warning  ('1D-function is not allowed for parallel add_branch' ) 
        
        ## use top-level function 
        fun         =  ( make_fun1 ( fun_ftwo , forcepc = True ) , 'pt' )
        
        chain    = data.chain
        vname    = 'doubled_pt1'
        chain. add_new_branch ( vname , fun  )
        
        logger.info ( "With doubled pt:\n%s" % data.chain.table ( prefix = '# ' ) )
        assert vname in data.chain , "Branch `%s' is  not here!" % vname 
        
    ## add lambda 
    with timing ('1D-lambda' , logger = logger ) :
        
        logger.warning  ('1D-lambda  is not allowed for parallel add_branch' )
        
        ftwo     = lambda x : 2 * x 
        fun      =  ( make_fun1 ( ftwo , forcepc = True ) , 'pt' )
        
        chain    = data.chain
        vname    = 'doubled_pt2'
        chain. add_new_branch ( vname , fun  )
        
        logger.info ( "With doubled pt:\n%s" % data.chain.table ( prefix = '# ' ) )
        assert vname in data.chain , "Branch `%s' is  not here!" % vname 
    
    ## add callable
    with timing ('1D-callable' , logger = logger ) :            
        
        logger.warning  ('1D-callable is not allowed for parallel add_branch' )
        
        ## top-level callabke 
        class CALL(object):
            def __call__ ( self , x ) : return 2.0 * x
            
        ftwo = CALL()
        fun      =  ( make_fun1 ( ftwo , forcepc = True ) , 'pt' )
        
        chain    = data.chain
        vname    = 'doubled_pt3'
        chain. add_new_branch ( vname , fun  )
        
        logger.info ( "With doubled pt:\n%s" % data.chain.table ( prefix = '# ' ) )
        assert vname in data.chain , "Branch `%s' is  not here!" % vname 
            
    ## add 2D-function
    with timing ('2D-function ' , logger = logger ) :
        
        logger.warning  ('2D-function is not allowed for parallel add_branch' ) 

        def fff ( x , y  ) : return x * y 
        fun         =  ( make_fun2 ( fff ) , 'pt' , 'et')
        
        chain    = data.chain
        vname    = 'pt_mult_et'
        chain. add_new_branch ( vname , fun  )
        
        logger.info ( "With pt*et:\n%s" % data.chain.table ( prefix = '# ' ) )
        assert vname in data.chain , "Branch `%s' is  not here!" % vname 
            
    ## add 3D-function
    with timing ('3D-function ' , logger = logger ) :
        
        logger.warning  ('3d-function is not allowed for parallel add_branch' ) 

        def fff ( x , y  , z ) : return x * y * z  
        fun         =  ( make_fun3 ( fff ) , 'pt' , 'et' , 'et2' )
        
        chain    = data.chain
        vname    = 'pt_mult_et_e2'
        chain. add_new_branch ( vname , fun  )
        
        logger.info ( "With pt*et*et2:\n%s" % data.chain.table ( prefix = '# ' ) )
        assert vname in data.chain , "Branch `%s' is  not here!" % vname 
    
# =============================================================================
if '__main__' ==  __name__  :

    test_addbranch()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
