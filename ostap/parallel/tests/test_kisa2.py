#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id:$ 
# =============================================================================
## @file TestKisa2.py
#
#  Test for parallel data processing 
#
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2016-02-07
# 
#                    $Revision: 195426 $
#  Last modification $Date: 2015-10-01 15:34:29 +0200 (Thu, 01 Oct 2015) $
#                 by $Author: ibelyaev $
# =============================================================================
"""
Test for parallel data processing 
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2016-02-07"
__all__     = ()  ## nothing to be imported 
# =============================================================================
import ROOT,os 
import Ostap.PyRoUts 
from   Ostap.Data    import Data, DataAndLumi
from   Ostap.Utils   import timing 
# =============================================================================
# logging 
# =============================================================================
from AnalysisPython.Logger import getLogger
if '__main__' == __name__ or '__builtin__' == __name__ : 
    logger = getLogger( 'Ostap.TestKisa2' )
else : 
    logger = getLogger( __name__ )
# =============================================================================
logger.info('Test  Data')
# =============================================================================
# data patterns:
ganga = '/afs/cern.ch/work/i/ibelyaev/public/GANGA/workspace/ibelyaev/LocalXML'
if os.path.exists ( '/mnt/shared/VMDATA/LocalXML' ) :
    ganga = '/mnt/shared/VMDATA/LocalXML'
# =============================================================================
patterns = [
    ganga + '/358/*/output/CY.root' , ## 2k+11,down
    ganga + '/359/*/output/CY.root' , ## 2k+11,up
    ganga + '/356/*/output/CY.root' , ## 2k+12,down
    ganga + '/357/*/output/CY.root' , ## 2k+12,up
    ]

data1 = Data ( 'YD0/CY'   , patterns )
logger.info  ( 'DATA %s' % data1       )

# =============================================================================
logger.info('Test  Data&Lumi')
# =============================================================================

data7 = DataAndLumi ( 'Y/Y' , patterns[ :2])
data8 = DataAndLumi ( 'Y/Y' , patterns[2:] )
dataY = DataAndLumi ( 'Y/Y' , patterns     )
logger.info ( 'DATA@  7TeV %s' % data7    )
logger.info ( 'DATA@  8TeV %s' % data8    )
logger.info ( 'DATA@7+8TeV %s' % dataY    )

logger.info ( 'TChain  %s' % data7.chain )

variables = {
    ( 'mass' , 'mass(mu+mu-)' , 8.5 , 11.5 , lambda s : s.mass ) 
    }

import Ostap.Kisa as Kisa

ppservers = () ## 'lxplus051' , )


one_file =  dataY.chain[:1]
with timing('One file %s' % len ( one_file ) ) : 
    ds0 = Kisa.fillDataSet (
        one_file        ,
        variables       ,
        '8.5<=mass && mass<11.5 && -0.1<c2dtf && c2dtf<5' , ppservers = ppservers ) 
    logger.info ( 'Dataset: %s' % ds0 )


with timing('All files %s' % len( dataY.chain ) ) :  
    dsa = Kisa.fillDataSet (
        dataY.chain ,
        variables   ,
        '8.5<=mass && mass<11.5 && -0.1<c2dtf && c2dtf<5' , ppservers = ppservers ) 
    
    logger.info ( 'Dataset: %s' % dsa )

# =============================================================================
# The END 
# =============================================================================
