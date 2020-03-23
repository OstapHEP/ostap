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
import ROOT, math, random 
import ostap.trees.trees
import ostap.histos.histos
from   ostap.trees.data         import Data 
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
    
    files = prepare_data ( 100 , 1000 )
    ## files = prepare_data ( 2 , 10 )
        
    logger.info ( '#files:    %s'  % len ( files ) )  
    data = Data ( 'S' , files )
    logger.info ( 'Initial Tree/Chain:\n%s' % data.chain.table ( prefix = '# ' ) )
    
    # =========================================================================
    ## 1) add new branch as TTree-formula:
    # =========================================================================
    data.chain.add_new_branch ( 'et','sqrt(pt*pt+mass*mass)' )

    ## reload the chain and check: 
    data.reload ()
    logger.info ( 'With formula:\n%s' % data.chain.table ( prefix = '# ' ) )
    assert 'et' in data.chain , "Branch ``et'' is  not here!"


    # =========================================================================
    ## 2) add several new branches as TTree-formula:
    # =========================================================================
    data.chain.add_new_branch ( { 'Et1' : 'sqrt(pt*pt+mass*mass)'   ,
                                  'Et2' : 'sqrt(pt*pt+mass*mass)*2' ,
                                  'Et3' : 'sqrt(pt*pt+mass*mass)*3' } , None )
    
    ## reload the chain and check: 
    data.reload ()
    logger.info ( 'With formula:\n%s' % data.chain.table ( prefix = '# ' ) )
    assert 'Et1' in data.chain , "Branch ``Et1'' is  not here!"
    assert 'Et2' in data.chain , "Branch ``Et2'' is  not here!"
    assert 'Et3' in data.chain , "Branch ``Et3'' is  not here!"




    # =========================================================================
    ## 2) add new branch as pure python function 
    # =========================================================================
    et2 = lambda tree : tree.pt**2 + tree.mass**2

    data.chain.add_new_branch ( 'et2', et2 )

    ## reload the chain and check: 
    data.reload ()
    logger.info ( 'With python:\n%s' % data.chain.table ( prefix = '# ' ) )
    assert 'et2' in data.chain , "Branch ``et2'' is  not here!"

    # =========================================================================
    ## 3) add new branch as histogram-function 
    # =========================================================================
    h1  = ROOT.TH1D ( 'h1' , 'some pt-correction' , 100 , 0 , 10 )
    h1 += lambda x :  1.0 + math.tanh( 0.2* ( x - 5 ) ) 
    
    from   ostap.trees.funcs  import FuncTH1
    ptw = FuncTH1 ( h1 , 'pt' )
    data.chain.add_new_branch ( 'ptw', ptw ) 
    
    ## reload the chain and check: 
    data.reload ()
    logger.info ( 'With histogram:\n%s' % data.chain.table ( prefix = '# ' ) )
    assert 'ptw' in data.chain , "Branch ``ptw'' is  not here!"

    # =========================================================================
    ## 4) add the variable sampled from the histogram
    # =========================================================================
    h2 = ROOT.TH1D('2', 'Gauss' , 120 , -6 , 6 )
    for i in range ( 100000 ) :
        h2.Fill ( random.gauss ( 0 , 1 ) ) 
        
    data.chain.add_new_branch ( 'hg', h2 ) 
    
    ## reload the chain and check: 
    data.reload ()
    logger.info ( 'With sampled:\n%s' % data.chain.table ( prefix = '# ' ) )
    assert 'hg' in data.chain , "Branch ``g'' is  not here!"
    
# =============================================================================
if '__main__' ==  __name__  :

    test_addbranch()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
