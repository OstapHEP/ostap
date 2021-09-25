#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/trees/tests/test_trees_addbranch.py
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
from   __future__               import print_function
import ROOT, math, random, array  
import ostap.trees.trees
import ostap.histos.histos
from   ostap.trees.data         import Data
from   ostap.utils.timing       import timing 
from   ostap.utils.progress_bar import progress_bar
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_trees_cuts' )
else                       : logger = getLogger ( __name__           )
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

# =============================================================================
# Four ways to add branch into TTree/Tchain
# - using string formula (TTreeFormula-based)
# - using pure python function
# - using histogram/function
# - using histogram sampling
def test_addbranch() : 
    """Four ways to add branch into TTree/Tchain
    - using string formula (TTreeFormula-based)
    - using pure python function
    - using histogram/function
    - using histogram sampling
    """
    
    files = prepare_data ( 10 , 1000 )
    ## files = prepare_data ( 2 , 10 )
    
    logger.info ( '#files:    %s'  % len ( files ) )  
    data = Data ( 'S' , files )
    logger.info ( 'Initial Tree/Chain:\n%s' % data.chain.table ( prefix = '# ' ) )
    
    # =========================================================================
    ## 1) add new branch as TTree-formula:
    # =========================================================================
    chain = data.chain 
    chain.add_new_branch ( 'et','sqrt(pt*pt+mass*mass)' )

    ## reload the chain and check: 
    logger.info ( 'With formula:\n%s' % data.chain.table ( prefix = '# ' ) )
    assert 'et' in data.chain , "Branch ``et'' is  not here!"


    # =========================================================================
    ## 2) add several new branches as TTree-formula:
    # =========================================================================
    chain = data.chain  
    chain.add_new_branch ( { 'Et1' : 'sqrt(pt*pt+mass*mass)'   ,
                             'Et2' : 'sqrt(pt*pt+mass*mass)*2' ,
                             'Et3' : 'sqrt(pt*pt+mass*mass)*3' } , None )
    

    ## reload the chain and check: 
    logger.info ( 'With formula:\n%s' % data.chain.table ( prefix = '# ' ) )
    assert 'Et1' in data.chain , "Branch ``Et1'' is  not here!"
    assert 'Et2' in data.chain , "Branch ``Et2'' is  not here!"
    assert 'Et3' in data.chain , "Branch ``Et3'' is  not here!"

    # =========================================================================
    ## 2) add new branch as pure python function 
    # =========================================================================
    et2 = lambda tree : tree.pt**2 + tree.mass**2

    chain = data.chain
    chain.add_new_branch ( 'et2', et2 )

    ## reload the chain and check: 
    logger.info ( 'With python:\n%s' % data.chain.table ( prefix = '# ' ) )
    assert 'et2' in data.chain , "Branch ``et2'' is  not here!"

    # =========================================================================
    ## 3) add new branch as histogram-function 
    # =========================================================================
    h1  = ROOT.TH1D ( 'h1' , 'some pt-correction' , 100 , 0 , 10 )
    h1 += lambda x :  1.0 + math.tanh( 0.2* ( x - 5 ) ) 
    
    from   ostap.trees.funcs  import FuncTH1
    ptw = FuncTH1 ( h1 , 'pt' )
    
    chain = data.chain 
    chain.add_new_branch ( 'ptw', ptw ) 
    
    ## reload the chain and check: 
    logger.info ( 'With histogram:\n%s' % data.chain.table ( prefix = '# ' ) )
    assert 'ptw' in data.chain , "Branch ``ptw'' is  not here!"

    # =========================================================================
    ## 4) add the variable sampled from the histogram
    # =========================================================================
    h2 = ROOT.TH1D('h2', 'Gauss' , 120 , -6 , 6 )
    for i in range ( 100000 ) :
        h2.Fill ( random.gauss ( 0 , 1 ) ) 

    chain = data.chain 
    chain.add_new_branch ( 'hg', h2 ) 
    
    ## reload the chain and check: 
    logger.info ( 'With sampled:\n%s' % data.chain.table ( prefix = '# ' ) )
    assert 'hg' in data.chain , "Branch ``hg'' is  not here!"

    # =========================================================================
    ## 5) python function again 
    # =========================================================================
    def gauss ( *_ ) : return random.gauss(0,1)
    
    chain = data.chain 
    chain.add_new_branch ( 'gauss', gauss ) 
    
    ## reload the chain and check: 
    logger.info ( 'With gauss:\n%s' % data.chain.table ( prefix = '# ' ) )
    assert 'gauss' in data.chain , "Branch ``gauss'' is  not here!"

    # =========================================================================
    ## 6) add numpy array 
    # =========================================================================
    try : 
        import numpy
    except ImportError :
        numpy  = None

        
    if numpy :

        with timing ('numpy float64' , logger = logger ) :
            adata  = numpy.ones ( 1000 , dtype = numpy.float64 )
            chain  = data.chain            
            chain.add_new_branch ( 'np_f64' , adata )
        ## reload the chain and check: 
        logger.info ( 'With numpy.float64:\n%s' % data.chain.table ( prefix = '# ' ) )
        assert 'np_f64' in data.chain , "Branch ``np_f64'' is  not here!"

        with timing ('numpy float32' , logger = logger ) :
            adata = numpy.ones ( 1000 , dtype = numpy.float32 )
            chain = data.chain
            chain.add_new_branch ( 'np_f32' , adata )            
        ## reload the chain and check: 
        logger.info ( 'With numpy.float32:\n%s' % data.chain.table ( prefix = '# ' ) )
        assert 'np_f32' in data.chain , "Branch ``np_f32'' is  not here!"
        
        with timing ('numpy float16' , logger = logger ) :
            adata  = numpy.ones ( 1000 , dtype = numpy.float16 )
            chain  = data.chain
            chain.add_new_branch ( 'np_f16' , adata )
        ## reload the chain and check: 
        logger.info ( 'With numpy.float16:\n%s' % data.chain.table ( prefix = '# ' ) )
        assert 'np_f16' in data.chain , "Branch ``np_f16'' is  not here!"

        with timing ('numpy int62 ' , logger = logger ) :
            adata  = numpy.ones ( 10000 , dtype = numpy.int32)
            chain  = data.chain            
            chain.add_new_branch ( 'np_i32' , adata )
        ## reload the chain and check: 
        logger.info ( 'With numpy.int32:\n%s' % data.chain.table ( prefix = '# ' ) )
        assert 'np_i32' in data.chain , "Branch ``np_i32'' is  not here!"

        with timing ('numpy int64 ' , logger = logger ) :  
            adata  = numpy.ones ( 10000 , dtype = numpy.int64 )
            chain  = data.chain
            chain.add_new_branch ( 'np_i64' , adata )
        ## reload the chain and check: 
        logger.info ( 'With numpy.int64:\n%s' % data.chain.table ( prefix = '# ' ) )
        assert 'np_i64' in data.chain , "Branch ``np_i64'' is  not here!"

        with timing ('numpy int16 ' , logger = logger ) :  
            adata  = numpy.ones ( 10000 , dtype = numpy.int16 )
            chain  = data.chain            
            chain.add_new_branch ( 'np_i16' , adata )
        ## reload the chain and check: 
        logger.info ( 'With numpy.int16:\n%s' % data.chain.table ( prefix = '# ' ) )
        assert 'np_i16' in data.chain , "Branch ``np_i16'' is  not here!"

        with timing ('numpy int8 ' , logger = logger ) :  
            adata  = numpy.ones ( 10000 , dtype = numpy.int8 )
            chain  = data.chain            
            chain.add_new_branch ( 'np_i8' , adata )
        ## reload the chain and check: 
        logger.info ( 'With numpy.int8:\n%s' % data.chain.table ( prefix = '# ' ) )
        assert 'np_i8' in data.chain , "Branch ``np_i8'' is  not here!"


    for l in ( 'f' , 'd' , 'i' , 'l' , 'I' , 'L' ) :

        with timing ('array %s'% l , logger = logger ) :  
            adata  = array.array ( l ,  10000*[0] ) 
            chain  = data.chain
            vname  = 'arr_%s' % l 
            chain.add_new_branch ( vname , adata )
            ## reload the chain and check: 
        logger.info ( "With array '%s':\n%s" % ( l ,  data.chain.table ( prefix = '# ' ) ) ) 
        assert vname in data.chain , "Branch ``%s'' is  not here!" % vname 
            


    
        
# =============================================================================
if '__main__' ==  __name__  :

    test_addbranch()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
