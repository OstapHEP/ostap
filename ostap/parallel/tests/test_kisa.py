#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file test_kisa.py
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
from   ostap.trees.data   import Data
from   ostap.utils.timing import timing 
import ostap.parallel.kisa  ## ATTENTION! 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' == __name__ or '__builtin__' == __name__ : 
    logger = getLogger( 'test_kisa' )
else : 
    logger = getLogger( __name__ )
# =============================================================================
patterns = [] 
def prepare_data ( nfiles =  100 ,  nentries = 100000 ) :

    from array import array 
    var1 = array ( 'd', [0])
    var2 = array ( 'd', [0])
    var3 = array ( 'd', [0])
    
    import  tempfile 
    tmpdir =   tempfile.mkdtemp()
    patterns.append ( os.path.join ( tmpdir , '*.root' ) )
    for i  in xrange( nfiles ) :
        fname = os.path.join ( tmpdir  , '%s.root' %  i )
 
        with ROOT.TFile.Open( fname , 'new' ) as root_file:
            
            tree = ROOT.TTree('S','tree')
            tree.SetDirectory ( root_file ) 
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
            
import atexit
@atexit.register
def cleanup () :
    import   glob
    while patterns:        
        p = patterns.pop()
        for f in glob.iglob ( p ) :
            try    : os.remove (  f )
            except : pass
            

while not patterns :
    with timing('Prepare data') :
        logger.info('Prepare data, it could take some time') 
        prepare_data ( 1000 , 20000 )
    data = Data  ( 'S' , patterns    )
    logger.info  ( 'DATA %s' % data  )
    
def test_kisa () : 
    
    h1 = ROOT.TH1D( 'h1' , '' , 200 , 3 , 3.2 )
    h2 = h1.clone()
    
    chain = data.chain
    
    
    with timing('SEQUENTIAL(%s):' % len(chain) , logger ) :
        chain. project ( h1 , 'mass' , '3<=mass && mass<=3.2 && 0<=c2dtf && c2dtf<5' )
        
    logger.info ( h1.dump(100,30) ) 
    
    with timing('PARALLEL(%s):' % len(chain) , logger ) :
        chain.pproject ( h2 , 'mass' , '3<=mass && mass<=3.2 && 0<=c2dtf && c2dtf<5' , silent = True )
        
    logger.info ( h2.dump(100,30) ) 


class MASS (object):
    def __call__ (  self , s ) :
        return s.mass
def MASS1  ( s ) : return s.mass

def test_kisa2 () :

       
    from ostap.fitting.selectors import SelectorWithVars, Variable  
    variables = [
        ## Variable ( 'mass1' , 'mass(mu+mu-)' , 2 , 4 , lambda s : s.mass ) , 
        Variable ( 'mass2' , 'mass(mu+mu-)' , 1 , 5 , MASS1  ) , 
        Variable ( 'mass3' , 'mass(mu+mu-)' , 1 , 5 , MASS() ) 
        ]
    
    
    ppservers = () ## 'lxplus051' , )

    nf   = len ( data.files )
    nf //= 40
    nf  += 1 
    
    with timing('%d files in sequence %s' % ( nf , len( data.chain )  ) ) :
        selector = SelectorWithVars  (
            variables = variables ,
            selection =  '2<=mass && mass<4 && 0<=c2dtf && c2dtf<5' ,
            silence   = True 
            )
        chain =  data.chain[:nf]
        st = chain.process ( selector , silent = True )
        ds = selector.data 
    logger.info ( 'Dataset: %s' % ds )
    
    with timing('%s files in parallel %s' % ( len ( data.files ) , len( data.chain ) ) ) :
        selector = SelectorWithVars  (
            variables = variables ,
            selection =  '2<=mass && mass<4 && 0<=c2dtf && c2dtf<5' ,
            silence   = True 
            )
        st = data.chain.pprocess ( selector , silent = True )
        ds = selector.data 
    logger.info ( 'Dataset: %s' % ds )


# =============================================================================
if '__main__' == __name__ :

    test_kisa  ()
    test_kisa2 ()

    
# =============================================================================
# The END 
# =============================================================================
