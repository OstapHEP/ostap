#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/parallel/tests/test_parallel_kisa.py
#  Test for parallel data processing 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2015-05-17
# =============================================================================
from __future__                 import  print_function 
"""Test for parallel data processing 
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-06-08"
__all__     = ()  ## nothing to be imported 
# =============================================================================
from   builtins                     import range
from   ostap.core.pyrouts           import dsID,   Ostap 
from   ostap.utils.timing           import timing
from   ostap.trees.data             import Data
from   ostap.utils.progress_bar     import progress_bar
from   ostap.core.core              import ROOTCWD 
import ostap.trees.trees       
import ostap.fitting.roofit
import ostap.fitting.dataset
from   ostap.core.meta_info         import root_version_int
from   ostap.utils.timing           import timing
import ostap.parallel.parallel_fill 
import ROOT,os,sys, random
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
    """ Create a file with a tree
    >>> create_tree ( ('1.root' ,  1000 ) ) 
    """
    
    import ROOT, random 
    import ostap.io.root_file
    
    from array import array 
    var1 = array ( 'd', [0.])
    var2 = array ( 'd', [0.])
    var3 = array ( 'd', [0.])
    
    with ROOTCWD() , ROOT.TFile.Open ( fname , 'new' ) as root_file:
        
        tree = ROOT.TTree ( 'S','tree' )
        tree.SetDirectory ( root_file  ) 
        tree.Branch ( 'mass'  , var1 , 'mass/D'  )
        tree.Branch ( 'c2dtf' , var2 , 'c2dtf/D' )
        tree.Branch ( 'pt'    , var3 , 'pt/D'    )
        
        for i in range ( nentries ) : 
            
            m  = random.gauss        ( 3.1 ,  0.015 )
            c2 = random.gammavariate ( 2.5 , 0.5    ) / 5 
            pt = random.uniform      ( 0   , 20     )
            
            var1[0] = m
            var2[0] = c2 
            var3[0] = pt
            
            tree.Fill()
            
        root_file.Write()
        tree = None
        
    return  fname

# ==============================================================================================
##  prepare data for the  test 
def prepare_data ( nfiles   = 10  , nentries = 100 ) :
    """ Prepare data for the test
    """

    files = [] 
    for i in progress_bar  ( range ( nfiles ) ) :
        
        from ostap.utils.cleanup import CleanUp
        tmpfile = CleanUp.tempfile ( prefix = 'ostap-test-selectors-' , suffix = '.root' )        
        files.append ( create_tree ( tmpfile , nentries ) ) 

    files.sort() 
    return Data ( 'S' ,  files ) 

with timing ("Prepare data ", logger ) :
    data = prepare_data ( nfiles = 10 , nentries = 20000 ) 

## variables 
mass   = ROOT.RooRealVar ( "mass"  , '' , 3.0 , 3.2  )
c2dtf  = ROOT.RooRealVar ( "c2dtf" , '' , 0   , 1000 )
pt     = ROOT.RooRealVar ( "pt"    , '' , 0   , 15   )

cuts  = ROOT.TCut( "3.015<mass" ) & "mass<3.150"
cuts &= "0 < c2dtf "
cuts &= "c2dtf < 3.0"
cuts &= "1.5 < pt"
cuts &= "pt< 10"


# ============================================================================
##  Simple test
#   - loop over the entries in the chain
#   - select good entries
#   - fill dataset 
def test_simple_loop ()  :
    """ Simple test
    - loop over the entries in the chain
    - select good entries
    - fill dataset
    """
    
    logger = getLogger("test_simple_loop")

    with timing ( "Simple loop" , logger ) :
        
        varset  = ROOT.RooArgSet   ( mass , c2dtf , pt )
        dataset = ROOT.RooDataSet  ( dsID() , 'Test Data set-1' , varset )
        
        for i in data.chain :
            
            v_mass  = float ( i.mass  )
            v_c2dtf = float ( i.c2dtf )
            v_pt    = float ( i.pt    )
            
            ## apply cuts 
            if not 3.015 < v_mass  < 3.150 : continue
            if not 0.0   < v_c2dtf < 3.0   : continue
            if not 1.5   < v_pt    < 10    : continue
            
            ## fill dataset 
            mass .value = v_mass
            c2dtf.value = v_c2dtf
            pt   .value = v_pt 
            dataset.add ( varset ) 
        
    logger.info ("Data set (simple-loop):\n%s"  % dataset.table ( prefix = "# " ) )


# =============================================================================
## More advanced test
#   - loop over the good entries in the chain,
#   - fill dataset 
def test_loop_with_cuts ()  :
    """ More advanced test
    - loop over the good entries in the chain,
    - fill dataset    
    """
    
    logger = getLogger("test_loop_with_cuts")

    with timing ( "Loop with cuts" , logger ) :
        
        varset = ROOT.RooArgSet   ( mass , c2dtf , pt )
        dataset = ROOT.RooDataSet  ( dsID() , 'Test Data set-1' , varset )
        
        for i , _ in data.chain.withCuts ( cuts  ) :
            
            v_mass  = i.mass
            v_c2dtf = i.c2dtf
            v_pt    = i.pt
            
            ## fill dataset 
            mass .value = v_mass
            c2dtf.value = v_c2dtf
            pt   .value = v_pt
            
            dataset.add ( varset ) 
                        
    logger.info ("Data set (loop-with-cuts):\n%s"  % dataset.table ( prefix = "# " ) )


# =============================================================================
## Use trivial selector to loop over the  chain 
#   - loop over the entries in the chain using selector
#   - filter good entries 
#   - fill dataset 
def test_simple_selector ()  :
    """ Use trivial selector to loop over the  chain 
    - loop over the entries in the chain using selector
    - filter good entries 
    - fill dataset 
    """
    
    logger = getLogger("test_simple_selector")
    
    from ostap.fitting.pyselectors import Selector
    
    class MySel(Selector) :
        
        def __init__ ( self , tree , dataset ) :
            super ( MySel , self ) .__init__ ( tree )
            self.__dataset = dataset

        @property
        def dataset ( self ) : return self.__dataset
        
        def process_entry ( self ) :

            tree = self.tree
            
            v_mass  = tree.mass
            v_c2dtf = tree.c2dtf
            v_pt    = tree.pt
            
            ## apply cuts 
            if not 3.015 < v_mass  < 3.150 : return True
            if not 0.0   < v_c2dtf < 3.0   : return True 
            if not 1.5   < v_pt    < 10    : return True 
            
            ## fill dataset 
            mass .value = v_mass
            c2dtf.value = v_c2dtf
            pt   .value = v_pt 
            self.__dataset.add ( varset ) 
            
            return True 
            
    with timing ( "Simple selector" , logger ) :
        
        varset  = ROOT.RooArgSet   ( mass , c2dtf , pt )
        dataset = ROOT.RooDataSet  ( dsID() , 'Test Data set-1' , varset )
        
        mySel = MySel ( data.chain , dataset )
        
        data.chain.process ( mySel )
        
    logger.info ("Data set (simple selector):\n%s"  % mySel.dataset.table ( prefix = "# " ) )


# =============================================================================
## Use selector-with-cuts to loop over good entries in  the  chain 
#   - loop over the good entries in the chain using selector
#   - fill dataset 
def test_selector_with_cuts ()  :
    """ Use selector-with-cuts to loop over good entries in  the  chain 
    - loop over the good entries in the chain using selector
    - fill dataset 
    """
    
    logger = getLogger("test_selector_with_cuts")

    from ostap.fitting.pyselectors import SelectorWithCuts 
    
    class MySel2 (SelectorWithCuts) :
        
        def __init__ ( self , tree , cuts , dataset ) :
            super(MySel2,self).__init__ (  selection = cuts , tree = tree , logger = logger )
            self.__dataset = dataset

        @property
        def dataset ( self ) : return self.__dataset

        def process_entry ( self ) :

            t = self.tree
            
            v_mass  = t.mass
            v_c2dtf = t.c2dtf
            v_pt    = t.pt
            
            ## fill dataset 
            mass .value = v_mass
            c2dtf.value = v_c2dtf
            pt   .value = v_pt 
            self.__dataset.add ( varset ) 
            
            return True 
            
    with timing ( "Selector with cuts" , logger ) :
        
        varset = ROOT.RooArgSet   ( mass , c2dtf , pt )
        dataset = ROOT.RooDataSet  ( dsID() , 'Test Data set-1' , varset )
        
        
        mySel = MySel2 ( data.chain , cuts , dataset )
        
        data.chain.process ( mySel )
        
        logger.info ("Data set (selector-with-cuts):\n%s"  % dataset.table ( prefix = "# " ) )
        

# =============================================================================
## Use dedicated selector-with-vars (1) to loop
#  over good entries in  the  chain and fill   dataset 
def test_selector_with_vars1 ()  :
    """ Use dedicated selector-with-vars (1) to loop
    over good entries in  the  chain and fill   dataset 
    """
    
    logger = getLogger("test_selector_with_vars1")

    from ostap.fitting.pyselectors import SelectorWithVars

    mySel = SelectorWithVars ( variables = [ mass , c2dtf , pt ] ,
                               selection = cuts  ,
                               logger    = logger )

    
    with timing ( "Selector with vars" , logger ) :
        Ostap.Utils.process ( data.chain , mySel )
        
    dataset = mySel.data
    
    logger.info ("Data set (selector-with-vars):\n%s"  % dataset.table ( prefix = "# " ) )


# =============================================================================
## Use dedicated selector-with-vars (2) to loop
#  over good entries in  the  chain and fill   dataset 
def test_selector_with_vars2 ()  :
    """Use dedicated selector-with-vars (2) to loop
    over good entries in  the  chain and fill   dataset 
    """
    
    logger = getLogger("test_selector_with_vars2")

    from ostap.fitting.pyselectors import SelectorWithVars
    
    
    mySel = SelectorWithVars ( variables = [ mass , c2dtf , pt ] ,
                               selection = cuts  ,
                               logger    = logger )
    
    with timing ( "Selector with vars&logic" , logger ) :
        data.chain.process ( mySel , shortcut = False )
        
    dataset = mySel.data
        
    logger.info ("Data set (selector-with-vars):\n%s"  % dataset.table ( prefix = "# " ) )

# =============================================================================
## Use dedicated selector-with-vars (3) to loop
#  over good entries in  the  chain and fill   dataset 
def test_selector_with_vars3 ()  :
    """ Use dedicated selector-with-vars (3) to loop
    over good entries in  the  chain and fill   dataset 
    """
    
    logger = getLogger("test_selector_with_vars3")

    from ostap.fitting.pyselectors import SelectorWithVars
    
    
    mySel = SelectorWithVars ( variables = [ mass , c2dtf , pt ] ,
                               selection = cuts  ,
                               logger    = logger )
    
    with timing ( "Selector with vars&logic&kisa" , logger ) :
        data.chain.pprocess ( mySel , shortcut = False )
        
    dataset = mySel.data
        
    logger.info ("Data set (selector-with-vars):\n%s"  % dataset.table ( prefix = "# " ) )

    
# ==============================================================================================
if '__main__' == __name__ :

    
    test_simple_loop         ()
    test_loop_with_cuts      ()    
    test_simple_selector     ()    
    test_selector_with_cuts  ()    
    test_selector_with_vars1 ()    
    test_selector_with_vars2 ()
    test_selector_with_vars3 ()
    
# ==============================================================================================
##                                                                                       The END
# ==============================================================================================
