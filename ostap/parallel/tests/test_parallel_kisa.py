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
import ostap.core.pyrouts 
from   ostap.utils.timing import timing 
##  
import ostap.parallel.kisa              ## ATTENTION!
import ostap.parallel.parallel_project  ## ATTENTION!
import ostap.parallel.parallel_fill     ## ATTENTION!

# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' == __name__ or '__builtin__' == __name__ : 
    logger = getLogger( 'ostap.test_parallel_kisa' )
else : 
    logger = getLogger( __name__ )
# =============================================================================

with timing('Prepare data') :
    logger.info('Prepare data, it could take some time')
    from ostap.utils.utils import CleanUp
    tmpdir  = CleanUp().tmpdir 
    import prepare_test_kisa as PTK
    ## data = PTK.prepare_data ( tmpdir , nfiles = 250 , nentries = 200000 , silent = False ) 
    data    = PTK.prepare_data ( tmpdir , nfiles = 100 , nentries = 50000  , silent = False ) 
    ## data = PTK.prepare_data ( tmpdir , nfiles = 10  , nentries = 100    , silent = False ) 
    logger.info    ( 'DATA %s' % data  )

# =============================================================================
def test_kisa () : 
    
    h1 = ROOT.TH1D( 'h1' , '' , 200 , 3 , 3.2 )
    h2 = h1.clone()
    
    chain = data.chain
    
    
    with timing('SEQUENTIAL(%s):' % len(chain) , logger ) :
        chain. project ( h1 , 'mass' , '3<=mass && mass<=3.2 && 0<=c2dtf && c2dtf<5' )
        
    logger.info ( h1.dump(100,30) ) 
    
    with timing('PARALLEL(%s):' % len(chain) , logger ) :
        chain.pproject ( h2 , 'mass' , '3<=mass && mass<=3.2 && 0<=c2dtf && c2dtf<5' , silent = False )
        
    logger.info ( h2.dump(100,30) ) 


class MASS (object):
    def __call__ (  self , s ) :
        return s.mass
def MASS1  ( s ) : return s.mass

# =============================================================================
def test_kisa2 () :
       
    from ostap.fitting.selectors import SelectorWithVars, Variable  
    variables = [
        ##Variable ( 'mass1' , 'mass(mu+mu-)' , 2 , 4 , lambda s : s.mass ) , 
        Variable   ( 'mass'  , 'mass(mu+mu-)' ,  3.09 , 3.11 ) , 
        Variable   ( 'c2dtf' , 'chi2(dtf)'    , -1    , 10   ) , 
        ##Variable ( 'mass2' , 'mass(mu+mu-)' , 1 , 5 , MASS()  ) , 
        ##Variable ( 'mass3' , 'mass(mu+mu-)' , 1 , 5 , 'mass'  ) 
        ]
    
    ppservers = () ## 'lxplus051' , )
    
    nf   = len ( data.files )
    nf //= 40
    nf  += 1 
    nf   = min ( nf , 25 )
    
    with timing('%d files in sequence %s' % ( nf , len( data.chain )  ) ) :
        selector = SelectorWithVars  (
            variables = variables ,
            selection =  '2<=mass && mass<4 && 0<=c2dtf && c2dtf<5' ,
            silence   = False
            )
        chain =  data.chain[:nf]
        st = chain.process ( selector , silent = False , shortcut = True )
        ds = selector.data
        del selector 
    ##logger.info ( 'Dataset: %s' % ds )
    
    with timing('%s files in parallel %s' % ( len ( data.files ) , len( data.chain ) ) ) :
        selector = SelectorWithVars  (
            variables = variables ,
            selection =  '2<=mass && mass<4 && 0<=c2dtf && c2dtf<5' ,
            silence   = True 
            )
        st = data.chain.pprocess ( selector , silent = False , chunk_size = -1 , max_files = 1 )
        ds = selector.data 
        del selector 
    ##logger.info ( 'Dataset: %s' % ds )


## # =============================================================================
## def test_kisa3 () :

##     h1  = ROOT.TH1D('h1','',100,0,20)
##     h1 += lambda x : x
    
##     from  ostap.trees.funcs import H1DFunc
##     xh1 = H1DFunc ( histo = h1 , xvar = 'pt' ) 
    
##     from ostap.fitting.selectors import SelectorWithVars, Variable  
##     variables = [
##         ##Variable ( 'mass1' , 'mass(mu+mu-)' , 2 , 4 , lambda s : s.mass ) , 
##         Variable   ( 'mass'  , 'mass(mu+mu-)' ,  3.09 , 3.11 ) , 
##         Variable   ( 'c2dtf' , 'chi2(dtf)'    , -1    , 10   ) , 
##         Variable   ( 'xpt'   , 'chi2(dtf)'    , -1    , 30   , xh1 ) , 
##         ##Variable ( 'mass2' , 'mass(mu+mu-)' , 1 , 5 , MASS()  ) , 
##         ##Variable ( 'mass3' , 'mass(mu+mu-)' , 1 , 5 , 'mass'  ) 
##         ]

    
##     with timing('fill it!' ) :
##         selector = SelectorWithVars  (
##             variables = variables ,
##             selection =  '2<=mass && mass<4 && 0<=c2dtf && c2dtf<5' ,
##             silence   = False
##             )
##         chain =  data.chain
##         st = chain.process ( selector , silent = False , shortcut = True )
##         ds = selector.data
##         del selector 

# =============================================================================
if '__main__' == __name__ :

    test_kisa  ()
    test_kisa2 ()
    #test_kisa3 ()
    
    pass


    
    
# =============================================================================
# The END 
# =============================================================================
