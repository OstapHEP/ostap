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
from   ostap.utils.timing             import timing
from   ostap.trees.data               import Data
from   ostap.utils.root_utils         import batch_env 
import ostap.parallel.parallel        as     Parallel  
import ostap.parallel.kisa              ## ATTENTION!
import ostap.parallel.parallel_project  ## ATTENTION!
import ostap.parallel.parallel_fill     ## ATTENTION!
import ostap.core.pyrouts 
import ostap.io.root_file 
import ROOT, array, os,  random  
##
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' == __name__ or '__builtin__' == __name__ : 
    logger = getLogger( 'test_parallel_kisa' )
else : 
    logger = getLogger( __name__ )
# =============================================================================
batch_env ( logger ) 
# =============================================================================
## create a file with tree 
def create_tree ( jobid , item ) :
    """ Create a file with a tree
    >>> create_tree ( ('1.root' ,  1000 ) ) 
    """
    
    fname , nentries = item 
    
    var1 = array.array ( 'd', [0])
    var2 = array.array ( 'd', [0])
    var3 = array.array ( 'd', [0])

    with ROOT.TFile( fname , 'new' ) as root_file:
        
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
        tree = None
        
    return  fname,

# ==============================================================================================
def prepare_data ( tmpdir , nfiles =  100 ,  nentries = 100 , ppservers = () , silent = True ) :
    
    ## Use generic Task from Kisa 
    from ostap.parallel.parallel import GenericTask as Task  
    task  = Task ( processor = create_tree ) 
    
    ## task  = PrepareTask () 
    wmgr  = Parallel.WorkManager( ppservers = ppservers , silent = silent )

    from ostap.utils.cleanup import CleanUp
    tmpfile = CleanUp.tempfile ( prefix = 'ostap-test-kisa-' , suffix = '.root' , dir = tmpdir )
    
    fname = '%s/test_kisa_%d.root'
    
    files = [ CleanUp.tempfile ( prefix = 'ostap-test-kisa-' , suffix = '.root' , dir = tmpdir ) for i in range(nfiles) ]
    
    wmgr.process (  task , [  ( f , nentries ) for f in files  ] )
    
    the_files = set() 
    for f in task.results () :
        if os.path.exists ( f ) :
            the_files.add ( f )
    
    from ostap.trees.data   import Data
    the_files = list( the_files )
    the_files.sort() 
    return Data ( the_files , 'S' ) 

with timing('Prepare data') :
    
    logger.info('Prepare data, it could take some time')
    from ostap.utils.cleanup import CleanUp
    tmpdir  = CleanUp().tmpdir

    ## data = prepare_data ( tmpdir , nfiles = 250 , nentries = 200000 , silent = False ) 
    ## data = prepare_data ( tmpdir , nfiles = 100 , nentries = 50000  , silent = False ) 
    ## data = prepare_data ( tmpdir , nfiles = 500  , nentries = 1000    , silent = False )
    ## data = prepare_data ( tmpdir , nfiles = 50  , nentries = 50000  , silent = False ) 
    data = prepare_data ( tmpdir , nfiles = 20  , nentries = 2000  , silent = False ) 
    
    logger.info    ( 'DATA %s' % data  )

# =============================================================================
def test_kisa () : 

    logger = getLogger ( 'test_parallel_kisa' )
    
    h1 = ROOT.TH1D( 'h1' , '' , 200 , 3 , 3.2 )
    h2 = h1.clone()
    
    chain = data.chain
    

    print ( chain.table () )
    
    with timing('SEQUENTIAL(%s):' % len(chain) , logger ) :
        chain. project ( h1              ,
                         'mass'          ,
                         '3<=mass && mass<=3.2 && 0<=c2dtf && c2dtf<5' ,
                         progress = True ) 
                        

    logger.info ( h1.dump(100,30) ) 


    with timing('PARALLEL(%s):' % len(chain) , logger ) :
        chain.project ( h2 ,
                        'emass' ,
                        '3<=mass && mass<=3.2 && 0<=c2dtf && c2dtf<5' ,
                        progress = True  , 
                        parallel = True  )
        
    logger.info ( h2.dump(100,30) ) 

# ============================================================================
class MASS (object):
    def __call__ (  self , s ) :
        return s.mass
def MASS1  ( s ) : return s.mass

# =============================================================================
def test_kisa2 () :
       
    logger = getLogger ( 'test_parallel_kisa2' )

    if 62400 <= ROOT.gROOT.GetVersionInt() < 62406 :
        logger.warning ('Test can fail for %s' % ROOT.gROOT.GetVersion() )

    from ostap.fitting.pyselectors import SelectorWithVars, Variable  
    variables = [
        Variable   ( 'mass'  , 'mass(mu+mu-)' ,  3.09 , 3.11 ) , 
        Variable   ( 'c2dtf' , 'chi2(dtf)'    , -1    , 10   ) , 
        Variable   ( 'mass2' , 'mass(mu+mu-)' , 1 , 5 , MASS()  ) , 
        Variable   ( 'mass3' , 'mass(mu+mu-)' , 1 , 5 , 'mass'  ) 
        ]
    
    ppservers = () ## 'lxplus051' , )
    ## ppservers = 'auto'
    
    nf   = len ( data.files )
    nf //= 40
    nf  += 1 
    nf   = min ( nf , 25 )
    
    with timing('%d files in sequence %s' % ( nf , len( data.chain )  ) ) :
        selector = SelectorWithVars  (
            variables = variables ,
            selection = '2<=mass && mass<4 && 0<=c2dtf && c2dtf<5' ,
            silence   = False
            )
        chain =  data.chain[:nf]
        st = chain.process ( selector , silent = False , shortcut = True )
        ds = selector.data
        del selector 
    logger.info ( 'Dataset (sequential):\n%s' % ds.table()  )
    
    with timing('%s files in parallel %s' % ( len ( data.files ) , len( data.chain ) ) ) :
        selector = SelectorWithVars  (
            variables = variables ,
            selection = '2<=mass && mass<4 && 0<=c2dtf && c2dtf<5' ,
            silence   = True 
            )
        st = data.chain.parallel_fill ( selector               ,
                                        silent     = False     ,
                                        chunk_size = -1        ,
                                        max_files  =  1        ,
                                        ppservers  = ppservers )
        ds = selector.data 
        del selector 
    logger.info ( 'Dataset (paralell):\n%s' % ds.table ( prefix = '# ' ) )

# =============================================================================
def test_kisa3 () :

    logger = getLogger ( 'test_parallel_kisa3' )

    h1  = ROOT.TH1D('h1','',100,0,20)
    h1 += lambda x : x
    
    ## from  ostap.trees.funcs import H1DFunc
    ## xh1 = H1DFunc ( histo = h1 , xvar = 'pt' ) 
    ## from  ostap.trees.funcs import FuncTH1 
    ## xh1 = FuncTH1 ( histo = h1 , xvar = 'pt' ) 
    
    from ostap.fitting.pyselectors import SelectorWithVars, Variable  
    variables = [
        ## Variable ( 'mass1' , 'mass(mu+mu-)' , 2 , 4 , lambda s : s.mass ) , 
        Variable ( 'mass'  , 'mass(mu+mu-)' ,  3.09 , 3.11 ) , 
        Variable ( 'c2dtf' , 'chi2(dtf)'    , -1    , 10   ) , 
        ## Variable ( 'xpt'   , 'xpt'          , -1    , 30   , xh1 ) , 
        Variable ( 'mass2' , 'mass(mu+mu-)' , 1 , 5 , MASS()  ) , 
        Variable ( 'mass3' , 'mass(mu+mu-)' , 1 , 5 , 'mass'  ) 
        ]

    
    with timing('fill it!' ) :
        selector = SelectorWithVars  (
            variables = variables ,
            selection =  '2<=mass && mass<4 && 0<=c2dtf && c2dtf<5' ,
            silence   = True 
            )
        chain =  data.chain
        st = chain.parallel_fill ( selector , silent = False , shortcut = True )
        ds = selector.data
        del selector 
    logger.info ( 'Dataset (paralell):\n%s' % ds.table ( prefix = '# ' ) )
        
# =============================================================================
if '__main__' == __name__ :

    test_kisa  ()
    ## test_kisa2 ()    
    ## test_kisa3 ()
    
# =============================================================================
##                                                                      The END 
# =============================================================================
