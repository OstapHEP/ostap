#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/trees/tests/test_trees_addpexps.py
# Add pseudoexperiments into TTree/TChain/RooDataSet
# Copyright (c) Ostap developers.
# ============================================================================= 
""" Add pseudoexperiments into TTree/TChain/RooDataSet
"""
# ============================================================================= 
from   __future__               import print_function
import ROOT, math, random 
import ostap.trees.trees
import ostap.histos.histos
from   ostap.core.pyrouts       import Ostap, VE,   SE  
from   ostap.trees.data         import Data 
from   ostap.utils.progress_bar import progress_bar
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_trees_addpexps' )
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
        tree.Branch ( 'pt'    , var2 , 'pt/D'    )
        tree.Branch ( 'eta'   , var3 , 'eta/D'   )
        
        for i in range ( nentries ) : 
            
            m   = random.gauss        ( 3.1 ,  0.015 )
            pt  = random.uniform      ( 0   , 10     )
            eta = random.uniform      ( 2   ,  5     )
            
            var1[0] = m
            var2[0] = pt
            var3[0] = eta
            
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
## let h2 be a weigth histogram 
h2  = ROOT.TH2D ( 'h2' , '', 20 , 0 , 10 , 15 , 2 , 5 )
h2 += lambda x,y : VE( 10 + x + y , ( 0.25 * ( 10 + x + y ) ) **2 ) 

cut  = ROOT.TCut ( 'abs(mass-3.1)<1*0.015' )

# =============================================================================
## Add pseudoexepriments into TTree/TChain
def test_modify_initial_tree ( NEXP =  10 ) :
    """Add pseudoexepriments into TTree/TChain
    """
    
    files = prepare_data ( 1 , 100000 )
    
    logger.info ( 'Add %s pseudoexepriments into TTree/TChain'  %   NEXP ) 

    logger.info ( '#files:    %s'  % len ( files ) )  
    data = Data ( 'S' , files )
    logger.info ( 'Initial Tree/Chain:\n%s' % data.chain.table ( prefix = '# ' ) )
    
    ## pseudo experiments
    for e in progress_bar ( range ( NEXP ) ) :
        h2_new = h2.sample()
        func   = Ostap.Functions.FuncTH2 ( h2_new , 'pt' , 'eta' ) 
        data.chain.add_new_branch ( 'w%d' % e , func )
        
    data = Data ( 'S' , files )
    logger.info ( 'Tree/Chain after:\n%s' % data.chain.table ( prefix = '# ' ) )

    counter = SE () 
    for e in range ( NEXP ) :
        weight     = 'w%d' % e 
        accepted   = data.chain.sumVar ( '1' , weight *  cut ) 
        rejected   = data.chain.sumVar ( '1' , weight * ~cut ) 
        efficiency = 1 / ( 1 + rejected / accepted )
        logger.info ( "Experiment %3d, accepted/rejected %s/%s , eff = %s "  % ( e          ,
                                                                                 accepted   ,
                                                                                 rejected   , 
                                                                                 efficiency ) )
        counter    += efficiency          
    logger.info ( 'Statistics of pseudoexperiments %s' % counter ) 
    logger.info ( 'Mean/rms: %s[%%]/%.4f]%%]' % ( counter.mean() * 100 ,
                                                  counter.rms () * 100 ) ) 
# =============================================================================
## Add pseudoexepriments into RooDataSet
def test_add_to_dataset ( NEXP =  10 ) :
    """Add pseudoexepriments into RooDataSet
    """
    
    logger.info ( 'Add %s pseudoexepriments into RooDataSet'  %   NEXP ) 
        
    files = prepare_data ( 1 , 100000 )


    logger.info ( '#files:    %s'  % len ( files ) )  
    data = Data ( 'S' , files )
    logger.info ( 'Initial Tree/Chain:\n%s' % data.chain.table ( prefix = '# ' ) )

    import ostap.fitting.selectors 
    dataset , _ = data.chain.fill_dataset ( ['mass','pt','eta'] )

    logger.info ( 'Initial dataset:\n%s' % dataset.table ( prefix = '# ' ) )
      
    ## pseudo experiments
    for e in progress_bar ( range ( NEXP ) ) :
        h2_new = h2.sample()
        func   = Ostap.Functions.FuncRooTH2 ( h2_new , 'pt' , 'eta' ) 
        dataset.add_new_var ( 'w%d' % e , func )

    logger.info ( 'Dataset after:\n%s' % dataset.table ( prefix = '# ' ) )
        
    counter = SE () 
    for e in range ( NEXP ) :
        weight     = 'w%d' % e 
        accepted   = dataset.sumVar ( '1' , weight *  cut ) 
        rejected   = dataset.sumVar ( '1' , weight * ~cut ) 
        efficiency = 1 / ( 1 + rejected / accepted )
        logger.info ( "Experiment %3d, accepted/rejected %s/%s , eff = %s "  % ( e          ,
                                                                                 accepted   ,
                                                                                 rejected   , 
                                                                                 efficiency ) )

        
        counter    += efficiency          
    logger.info ( 'Statistics of pseudoexperiments %s' % counter ) 
    logger.info ( 'Mean/rms: %s[%%]/%.4f[%%]' % ( counter.mean() * 100 ,
                                                  counter.rms () * 100 ) ) 
    
# =============================================================================
if '__main__' == __name__ :

    test_modify_initial_tree ( 49 )
    test_add_to_dataset      ( 49 )

# =============================================================================
##                                                                      The END 
# =============================================================================

