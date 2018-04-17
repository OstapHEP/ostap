#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file prepare_test_kisa.py
#  Prepare input data for parallel data processing 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2015-05-17
# =============================================================================
"""Prepare input data for parallel data processing 
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-06-08"
__all__     = (
    'prepare_data'  ## prepare input data for parallel data processing 
    ) 
# =============================================================================
import ROOT,os,  random  
import ostap.core.pyrouts 
from   ostap.trees.data        import Data
from   ostap.utils.timing      import timing 
import ostap.parallel.parallel as     Parallel  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' == __name__ or '__builtin__' == __name__ : 
    logger = getLogger( 'prepare_test_kisa' )
else : 
    logger = getLogger( __name__ )
# =============================================================================e
## create a file with tree 
def create_tree ( item ) :
    """Create a file with a tree
    >>> create_tree ( ('1.root' ,  1000 ) ) 
    """
    
    fname , nentries = item 
    
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
        
        for i in xrange ( nentries ) : 
            
            m  = random.gauss        ( 3.1 ,  0.015 )
            c2 = random.gammavariate ( 2.5 , 0.5    ) / 5 
            pt = random.uniform      ( 0   , 10     )
            
            var1[0] = m
            var2[0] = c2 
            var3[0] = pt
            
            tree.Fill()
            
        root_file.Write()
        
    return  fname,

# ==============================================================================================
def prepare_data ( tmpdir , nfiles =  100 ,  nentries = 100 , ppservers = () , silent = True ) :
    
    ## Use generic Task from Kisa 
    from ostap.parallel.kisa import GenericTask as Task  
    task  = Task ( processor = create_tree ) 
    
    ## task  = PrepareTask () 
    wmgr  = Parallel.WorkManager( ppservers = ppservers , silent = silent )

    from ostap.utils.utils import CleanUp
    tmpfile = CleanUp.tempfile ( prefix = 'test_kisa_' , suffix = '.root' , dir = tmpdir )
    
    fname = '%s/test_kisa_%d.root'
    
    files = [ CleanUp.tempfile ( prefix = 'test_kisa_' , suffix = '.root' , dir = tmpdir ) for i in range(nfiles) ]
    
    wmgr.process (  task , [ (f,nentries) for f in files  ] )
    
    the_files = set() 
    for f in task.output :
        if os.path.exists ( f ) :
            the_files.add ( f )
    
    from ostap.trees.data   import Data
    the_files = list( the_files )
    the_files.sort() 
    return Data ( 'S' , list ( the_files ) ) 


# =============================================================================
if '__main__' == __name__ :

    data = prepare_data ( '.' , silent = False ) 
    logger.info ( 'Data: %s' % data )
    
# =============================================================================
# The END 
# =============================================================================
