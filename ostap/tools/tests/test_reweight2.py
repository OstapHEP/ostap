#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file test_reweigt2.py
#
#  Test for 2D-reweighting machinery
# 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-05-11
# 
# =============================================================================
"""
Test for 2D-reweighting machinery
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-05-10"
__all__     = ()  ## nothing to be imported 
# =============================================================================
import ROOT, random, math, os, time 
from   ostap.core.pyrouts import *
import ostap.io.zipshelve as     DBASE
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__'  == __name__ : 
    logger = getLogger ( 'ostap.test_reweight2' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================    
logger.info ( 'Test for 2D-Reweighting machinery')
# =============================================================================
testdatadb = 'testdata2.db'
tag_data   = 'DATA2-histogram'
tag_datax  = 'DATAX-histogram'
tag_datay  = 'DATAY-histogram'
tag_mc     = 'MC2-dataset'

# =============================================================================
testdatadb = 'testdata2.db'
tag_data   = 'DATA2-histogram'
tag_datax  = 'DATAX-histogram'
tag_datay  = 'DATAY-histogram'
tag_mc     = 'MC2-dataset'
if not os.path.exists( testdatadb ) :
    #
    logger.info ( 'Test RANDOM data will be generated' )   
    ## prepare "data" histograms:
    # 1) 2D hstograms 
    hdata  = h2_axes ( [ i for i in range(0,21) ] ,
                       [ i for i in range(0,21) ] )
    # 2) non-equal bining 1D histogram for x-component    
    hxdata = h1_axis ( [    i     for i in  range (5 ) ] +
                       [  5+i*0.2 for i in  range (50) ] +
                       [ 15+i     for i in  range (6 ) ] )
    # 2) equal bining 1D histogram for y-component    
    hydata = h1_axis ( [ i*0.5 for i in  range (41) ]   )
    for i in range( 0,5000000 ) :
        v1 = random.gauss(10,3)
        v2 = random.gauss( 0,1)
        x  = v1 + 3 * v2
        y  = v1 - 3 * v2
        hdata .Fill(x,y)
        hxdata.Fill(x)
        hydata.Fill(y)
             
    #
    ## prepare MC-dataset: use some other distribution
    #
    x  = ROOT.RooRealVar ( 'x' , 'x-variable'   , 0 , 20 )
    y  = ROOT.RooRealVar ( 'y' , 'y-variable'   , 0 , 20 )
    ## book very simple data set
    import ostap.fitting.models as     Models
    
    ## MC: product of two exponentials 
    f_mc = Models.ExpoPol2D_pdf( 'B2' , x , y , nx = 0 , ny = 0 , taux =  0.1 , tauy = -0.1 )

    from ostap.fitting.roofit import useStorage
    with useStorage() : 

        varset   = ROOT.RooArgSet    ( x , y )
        dataset  = f_mc.pdf.generate ( varset , 1000000 )
        mctree   = dataset .store().tree()
        
    ## store DATA in DBASE 
    with DBASE.open( testdatadb , 'c' ) as db : 
        logger.info ( 'Test data is stored  in   DBASE "%s"' % testdatadb  )   
        db[ tag_data  ] = hdata 
        db[ tag_datax ] = hxdata 
        db[ tag_datay ] = hydata 
        db[ tag_mc    ] = mctree 
        db.ls() 

#
## Read data from DB
with DBASE.open ( testdatadb , 'r' ) as db :
    logger.info ( 'Test data is fetched from DBASE "%s"' % testdatadb )   
    db.ls()
    hdata  = db[ tag_data  ]
    hxdata = db[ tag_datax ]
    hydata = db[ tag_datay ]
    mctree = db[ tag_mc    ]
    
## prepare template histogram for MC 
hmc  = h2_axes ( [ 20.0/45*i for i in range(0,46) ] ,
                 [ 20.0/30*i for i in range(0,31) ] )
hmcx = h1_axis ( [ 20.0/50*i for i in  range (51) ] )
hmcy = h1_axis ( [ 20.0/45*i for i in  range (46) ] )

## prepare re-weighting machinery 
maxIter = 20  

## check database 
dbname = 'reweighting2.db'
import os
if not os.path.exists( dbname ) :
    logger.info('Create new weights DBASE') 
    db = DBASE.open ( dbname , 'c' ) ##  create new empty db 
    db.close()
else :
    logger.info('Existing weights DBASE will be used') 
    
#
## make reweigthing iterations
# 
from ostap.tools.reweight     import Weight, makeWeights 
from ostap.fitting.selectors  import SelectorWithVars 
from ostap.utils.memory       import memory
from ostap.utils.timing       import timing

def test_reweight2() :

    ## start iterations:
    for iter in range ( 0 , maxIter ) :

        with timing() , memory ( 'Iteration %d' % iter ) :
            
            weighting = [
                ## variable          address in DB    
                ( lambda s : s.x       , 'x-reweighting'  ) , 
                ( lambda s : s.y       , 'y-reweighting'  ) , 
                ( lambda s : (s.x,s.y) , '2D-reweighting' ) , 
                ]
    
            weighter   = Weight( dbname , weighting )
            ## variables to be used in MC-dataset 
            variables  = [
                ( 'x'      , 'vyx'   , 0  , 20 , lambda s : s.x ) ,  
                ( 'y'      , 'vxy'   , 0  , 20 , lambda s : s.y ) , 
                ( 'weight' , 'weight' ,            weighter      )  
                ]
    
            #
            ## create new "weighted" mcdataset
            # 
            selector = SelectorWithVars (
                variables ,
                '0<x && x<20 && 0<y && y<20'
                )
            
            mctree.process ( selector )
            mcds = selector.data             ## new reweighted dataset

            #
            ## update weights
            #
            
            plots    = [
                ( 'x'   , 'weight' , 'x-reweighting'  , hxdata , hmcx ) ,  
                ( 'y'   , 'weight' , 'y-reweighting'  , hydata , hmcy ) , 
                ( 'x:y' , 'weight' , '2D-reweighting' , hdata  , hmc  ) , 
                ]
            
            more = makeWeights ( mcds , plots , dbname , delta = 0.001 )
            
            ## make MC-histogram 
            mcds .project  ( hmcx , 'x'   , 'weight'  )
            mcds .project  ( hmcy , 'y'   , 'weight'  )
            mcds .project  ( hmc  , 'y:x' , 'weight'  )
            
            logger.info    ( 'Compare DATA and MC for iteration #%d' % iter )
            #
            ## compare the basic properties: mean, rms, skewness and kurtosis
            # 
            hxdata.cmp_prnt ( hmcx , 'DATA' , 'MC' , 'DATA(x) vs MC(x)' )
            hydata.cmp_prnt ( hmcy , 'DATA' , 'MC' , 'DATA(y) vs MC(y)' )
            #
            ## calculate the distances
            #
            dist = hxdata.cmp_dist ( hmcx , rescale = True )
            logger.info ('DATA(x)-MC(x) "distance"      %s' % dist )
            dist = hydata.cmp_dist ( hmcy , rescale = True )
            logger.info ('DATA(y)-MC(y) "distance"      %s' % dist )
            #
            ## calculate the 'orthogonality'
            #  
            cost = hxdata.cmp_cos  ( hmcx , rescale = True )
            logger.info ('DATA(x)-MC(x) "orthogonality" %s' % cost )
            cost = hydata.cmp_cos  ( hmcy , rescale = True )
            logger.info ('DATA(y)-MC(y) "orthogonality" %s' % cost )
            #
            
            ## final density on data 
            datax_density = hxdata.density()
            ## final density on mc 
            mcx_density   = hmcx.density()
            
            datax_density.red  ()
            mcx_density  .blue ()
            datax_density.draw ('e1')
            mcx_density  .draw ('e1 same')
            time.sleep(5)
            
            if not more : 
                logger.info    ( 'No more iterations are needed #%d' % iter )
                break
            
            if iter + 1 != maxIter :
                mcds.clear() 
                del mcds , selector
            else :
                del selector 
                
    datax_density.draw ('e1')
    mcx_density  .draw ('e1 same')
    time.sleep(60)

# =============================================================================
if '__main__' == __name__ :

    test_reweight2 () 

# =============================================================================
# The END 
# =============================================================================
