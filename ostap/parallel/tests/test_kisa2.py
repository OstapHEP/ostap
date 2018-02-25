#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file test_kisa2.py
#  Test for parallel data processing 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2016-02-07
# =============================================================================
"""Test for parallel data processing 
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2016-02-07"
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
    logger = getLogger( 'test_kisa2' )
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


class MASS (object):
    def __call__ (  self , s ) :
        return s.mass
def MASS1  ( s ) : return s.mass

from ostap.fitting.selectors import SelectorWithVars, Variable  
variables = [
    ## Variable ( 'mass1' , 'mass(mu+mu-)' , 50 , 120 , lambda s : s.mass ) , 
    Variable ( 'mass2' , 'mass(mu+mu-)' , 50 , 120 , MASS1  ) , 
    Variable ( 'mass3' , 'mass(mu+mu-)' , 50 , 120 , MASS() ) 
    ]

import ostap.parallel.kisa as kisa

ppservers = () ## 'lxplus051' , )

with timing('All files in sequence %s' % len( data.chain ) ) :
    selector = SelectorWithVars  (
        variables = variables ,
        selection =  '50<=mass && mass<120 &&  0<c2dtf && c2dtf<5' ,
        silence   = True 
        )
    st = data.chain.process ( selector , silent = True )
    ds = selector.data 
    logger.info ( 'Dataset: %s' % ds )
    
with timing('All files in parallel %s' % len( data.chain ) ) :
    selector = SelectorWithVars  (
        variables = variables ,
        selection =  '50<=mass && mass<120 &&  0<c2dtf && c2dtf<5' ,
        silence   = True 
        )
    st = data.chain.pprocess ( selector , silent = True )
    ds = selector.data 
    logger.info ( 'Dataset: %s' % ds )

# =============================================================================
# The END 
# =============================================================================
