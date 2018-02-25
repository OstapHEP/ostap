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
import ROOT,os 
import ostap.core.pyrouts 
from   ostap.trees.data   import Data
from   ostap.utils.timing import timing 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' == __name__ or '__builtin__' == __name__ : 
    logger = getLogger( 'test_kisa' )
else : 
    logger = getLogger( __name__ )
# =============================================================================
logger.info('Test  Data')
# =============================================================================
# data patterns:
ganga = '/afs/cern.ch/work/i/ibelyaev/public/GANGA/workspace/ibelyaev/LocalXML'
# =============================================================================
patterns = [
    ganga + '/690/*/output/ZC.root' , ## 2k+11,down
    ganga + '/691/*/output/ZC.root' , ## 2k+11,up
    ganga + '/692/*/output/ZC.root' , ## 2k+12,down
    ganga + '/693/*/output/ZC.root' , ## 2k+12,up
    ganga + '/708/*/output/ZC.root' , ## 2k+15,down
    ganga + '/709/*/output/ZC.root' , ## 2k+15,up
    ]
  
data = Data ( 'aZ0/Z0'   , patterns )
logger.info  ( 'DATA %s' % data     )


h1 = ROOT.TH1D( 'h1' , '' , 140 , 50 , 120 )
h2 = h1.clone()

chain = data.chain


with timing('SEQUENTIAL(%s):' % len(chain) , logger ) :
    chain. project ( h1 , 'mass' , '50<=mass && mass<=120 && 0<=c2dtf && c2dtf<5' )

logger.info ( h1.dump(100,30) ) 

import ostap.parallel.kisa

with timing('PARALLEL(%s):' % len(chain) , logger ) :
    chain.pproject ( h2 , 'mass' , '50<=mass && mass<=120 && 0<=c2dtf && c2dtf<5' , silent = True )

logger.info ( h2.dump(100,30) ) 



# =============================================================================
# The END 
# =============================================================================
