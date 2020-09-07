#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/parallel/tests/test_parallel_kisa.py
#  Test for parallel data processing 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2015-05-17
# =============================================================================
"""Test for parallel data processing 
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-06-08"
__all__     = ()  ## nothing to be imported 
# =============================================================================
import ROOT,os,  random  
from   builtins                 import range
from   ostap.core.pyrouts       import dsID,   Ostap 
from   ostap.utils.timing       import timing
from   ostap.trees.data         import Data
from   ostap.utils.progress_bar import progress_bar
import ostap.trees.trees       
import ostap.fitting.roofit
import ostap.fitting.dataset
from   ostap.core.meta_info     import root_version_int
from   ostap.utils.timing   import timing 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' == __name__ or '__builtin__' == __name__ : 
    logger = getLogger( 'ostap.test_selectors' )
else : 
    logger = getLogger( __name__ )
# =============================================================================
## create a file with tree 
def create_tree ( fname , nentries ) :
    """Create a file with a tree
    >>> create_tree ( ('1.root' ,  1000 ) ) 
    """
    
    import ROOT, random 
    import ostap.io.root_file
    
    from array import array 
    var1 = array ( 'd', [0])
    var2 = array ( 'd', [0])
    var3 = array ( 'd', [0])
    
    with ROOT.TFile.Open( fname , 'new' ) as root_file:
        
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
        
    return  fname

# ==============================================================================================
def prepare_data ( nfiles   = 10  ,
                   nentries = 100 ) :

    files = [] 
    for i in progress_bar  ( range ( nfiles ) ) :
        
        from ostap.utils.cleanup import CleanUp
        tmpfile = CleanUp.tempfile ( prefix = 'test_kisa_' , suffix = '.root' )        
        files.append ( create_tree ( tmpfile , nentries ) ) 

    files.sort() 
    return Data ( 'S' ,  files ) 


with timing ("Prepare data ", logger ) :
    data = prepare_data ( nfiles = 10 , nentries = 500 ) 

## variables 
mass   = ROOT.RooRealVar ( "mass"  , '' , 3.0 , 3.2  )
c2dtf  = ROOT.RooRealVar ( "c2dtf" , '' , 0   , 1000 )
pt     = ROOT.RooRealVar ( "pt"    , '' , 0   , 15   )

cuts  = ROOT.TCut( "3.015<mass" ) & "mass<3.150"
cuts &= "0 < c2dtf "
cuts &= "c2dtf < 3.0"
cuts &= "1.5 < pt"
cuts &= "pt< 10"

# =============================================================================
def test_no_selector()  :

    with timing ("No selector ", logger ) :
        
        varset  = ROOT.RooArgSet   ( mass , c2dtf , pt )
        dataset = ROOT.RooDataSet  ( dsID() , 'Test Data set-1' , varset )
        
        
        for i in data.chain :
            
            v_mass  = i.mass
            v_c2dtf = i.c2dtf
            v_pt    = i.pt
            
            ## apply cuts 
            if not 3.050 < v_mass  < 3.150 : continue
            if not 0.0   < v_c2dtf < 3.0   : continue
            if not 1.5   < v_pt    < 10    : continue

            ## fill dataset 
            mass.value = i.mass
            mass.c2dtf = i.c2dtf
            mass.pt    = i.pt 
            dataset.add ( varset ) 
            
    logger.info ("Data set (no selector):\n%s"  % dataset.table () )


# =============================================================================
def test_with_cuts ()  :

    with timing ("With cuts", logger ) :
        
        varset = ROOT.RooArgSet   ( mass , c2dtf , pt )
        dataset = ROOT.RooDataSet  ( dsID() , 'Test Data set-1' , varset )
        
        
        for i in data.chain.withCuts ( cuts  ) :
            
            v_mass  = i.mass
            v_c2dtf = i.c2dtf
            v_pt    = i.pt
            
            ## fill dataset 
            mass.value = i.mass
            mass.c2dtf = i.c2dtf
            mass.pt    = i.pt 
            dataset.add ( varset ) 
            
    logger.info ("Data set (with cuts):\n%s"  % dataset.table () )


# =============================================================================
def test_simple_selector ()  :

    from ostap.core.meta_info import root_version_int 
    if root_version_int >= 62200 :
        logger.warning("test_simple_selector: test is disabled for ROOT verison %s" % root_version_int )
        return 
        
    from ostap.fitting.pyselectors import Selector
    
    class MySel(Selector) :
        
        def __init__ ( self , tree , dataset ) :
            Selector.__init__ ( self , tree )
            self.__dataset = dataset

        @property
        def dataset ( self ) : return self.__dataset
        
        def process_entry ( self ) :
            
            t = self.tree() 
            
            v_mass  = t.mass
            v_c2dtf = t.c2dtf
            v_pt    = t.pt
            
            ## apply cuts 
            if not 3.050 < v_mass  < 3.150 : return True
            if not 0.0   < v_c2dtf < 3.0   : return True 
            if not 1.5   < v_pt    < 10    : return True 
            
            ## fill dataset 
            mass.value = v_mass
            mass.c2dtf = v_c2dtf
            mass.pt    = v_pt 
            self.__dataset.add ( varset ) 
            
            return True 
    
        
    with timing ("Simple selector", logger ) :

        varset  = ROOT.RooArgSet   ( mass , c2dtf , pt )
        dataset = ROOT.RooDataSet  ( dsID() , 'Test Data set-1' , varset )


        mySel = MySel  ( data.chain , dataset )
        mySel.set_tree ( data.chain )
        
        data.chain.process ( mySel )

            
    logger.info ("Data set (simple selector):\n%s"  % dataset.table () )


# =============================================================================
def test_selector_with_cuts ()  :

    from ostap.core.meta_info import root_version_int 
    if root_version_int >= 62200 :
        logger.warning("test_selector_with_cuts: test is disabled for ROOT verison %s" % root_version_int )
        return 

    from ostap.fitting.pyselectors import SelectorWithCuts 
    
    class MySel2 (SelectorWithCuts) :
        
        def __init__ ( self , tree , cuts , dataset ) :
            SelectorWithCuts.__init__ ( self , selection = cuts , tree = tree )
            self.__dataset = dataset

        @property
        def dataset ( self ) : return self.__dataset

        def process_entry ( self ) :
            
            t = self.tree() 
            
            v_mass  = t.mass
            v_c2dtf = t.c2dtf
            v_pt    = t.pt
            
            ## fill dataset 
            mass.value = v_mass
            mass.c2dtf = v_c2dtf
            mass.pt    = v_pt 
            self.__dataset.add ( varset ) 
            
            return True 
    
        
    with timing ("Simple selector", logger ) :

        varset = ROOT.RooArgSet   ( mass , c2dtf , pt )
        dataset = ROOT.RooDataSet  ( dsID() , 'Test Data set-1' , varset )


        mySel = MySel2 ( data.chain , cuts , dataset )
        mySel.set_tree ( data.chain )
        mySel.set_cuts ( cuts       )

        data.chain.process ( mySel )

            
    logger.info ("Data set (selector with cuts):\n%s"  % dataset.table () )

    
    
# ==============================================================================================
if '__main__' == __name__ :

    with timing("No selector"        , logger ) :
        test_no_selector         ()
        
    with timing("With cuts"          , logger ) :
        test_with_cuts           ()
        
    with timing("Simple selector"    , logger ) :
        test_simple_selector     ()
        
    with timing("Selector with cuts" , logger ) :
        test_selector_with_cuts  ()
    
# ==============================================================================================
##                                                                                       The END
# ==============================================================================================
