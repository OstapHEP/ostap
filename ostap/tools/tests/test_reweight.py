#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file test_reweight.py
#
#  Test for reweighting machinery
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-05-11
# =============================================================================
"""Test for reweighting machinery in  Ostap
"""
# =============================================================================
import ROOT, random, math, os, time  
from   ostap.core.pyrouts    import * 
import ostap.io.zipshelve as DBASE
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__'  == __name__ : 
    logger = getLogger ( 'test_reweight' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================    
logger.info ( 'Test for Reweighting machinery')
# =============================================================================
testdatadb = 'testdata.db'
tag_data   = 'DATA-histogram'
tag_mc     = 'MC-dataset'
# =============================================================================
## prepare data for tests 
def prepare ( ) :
    """Prepare data for tests
    """
    if not os.path.exists( testdatadb ) :
        #
        logger.info ( 'Test RANDOM data will be generated' )   
        ## prepare"data" histogram: lets use simple exponential and non-equidistance bins 
        hdata = h1_axis( [    i   for i in range(0,20) ] +
                         [ 20+i*2 for i in range(0,20) ] +
                         [ 60+i*4 for i in range(0,10) ] + [ 100 ] ) 
        for i in range( 0,500000 ) :
            v = random.random()
            v = -20*math.log(v)
            hdata.Fill(v)
        #
        ## prepare MC-dataset: use Gamma-distribution 
        #
        x        = ROOT.RooRealVar ( 'x' , 'variable'   ,    0 , 100 )
        ## book very simple data set
        import ostap.fitting.models as  Models 
        m_gamma0 = Models.GammaDist_pdf('GD0'  , x )
        mean     =  35.0
        variance =  30**2 
        m_gamma0.k    .setVal( mean**2 /variance )
        m_gamma0.theta.setVal( variance/mean     )
        
        from ostap.fitting.roofit import useStorage
        with useStorage() : 
            
            varset   = ROOT.RooArgSet  ( x )
            dataset  = m_gamma0.pdf.generate ( varset , 500000 )
            mctree   = dataset .store().tree()
        
        ## store DATA in DBASE 
        with DBASE.open( testdatadb , 'c' ) as db : 
            logger.info ( 'Test data is stored  in   DBASE "%s"' % testdatadb  )   
            db[ tag_data ] = hdata 
            db[ tag_mc   ] = mctree 
            db.ls() 
            
    #
    ## Read data from DB
    with DBASE.open ( testdatadb , 'r' ) as db :
        logger.info ( 'Test data is fetched from DBASE "%s"' % testdatadb )   
        db.ls()
        hdata  = db[ tag_data ]
        mctree = db[ tag_mc   ]

    return hdata, mctree 

## prepare template histogram for MC 
hmc = ROOT.TH1D('hMC','histo-template for MC', 60,0,100 )

## prepare re-weighting machinery 
maxIter = 20  

## check database 
dbname = 'reweighting.db'
import os
if not os.path.exists( dbname ) :
    logger.info('Create new weights DBASE') 
    db = DBASE.open ( dbname , 'c' ) ##  create new empty db 
    db.close()
else :
    logger.info('Existing weights DBASE will be used') 

    
# =============================================================================
##  test reweightinng machinery 
def test_reweight ( ) :
    """test reweightinng machinery
    """

    from ostap.tools.reweight    import Weight, makeWeights 
    from ostap.fitting.selectors import SelectorWithVars 
    
    ## read data histogram and MC-tree 
    hdata, mctree = prepare()  
    
    
    #
    ## make reweigthing iterations
    # 
    
    ## start iterations:
    for iter in range ( 0 , maxIter ) :
        
        weighting = [
            ## variable          address in DB    
            ( lambda s : s.x , 'x-reweighting'  ) , 
            ]
    
        weighter   = Weight( dbname , weighting )
        ## variables to be used in MC-dataset 
        variables  = [
            ( 'pt_x'   , 'pt_x'   , 0  , 100 , lambda s : s.x ) , 
            ( 'weight' , 'weight' ,            weighter       )  
            ]
        
        #
        ## create new "weighted" mcdataset
        # 
        selector = SelectorWithVars (
            variables ,
            '0<x && x<100 '
            )
        
        mctree.process ( selector )
        mcds = selector.data             ## new reweighted dataset
        
        print 'MC-dataset', mcds 
        
        #
        ## update weights
        #
        
        plots    = [
            ( 'pt_x'   , 'weight' , 'x-reweighting'  , hdata , hmc )  
            ]
        
        more = makeWeights ( mcds , plots , dbname , delta = 0.001 )

        ## make MC-histogram 
        mcds .project  ( hmc , 'pt_x' , 'weight'  )
        
        if 0 == iter % 2 or not more : 
            logger.info    ( 'Compare DATA and MC for iteration #%d' % iter )
            
            #
            ## compare the basic properties: mean, rms, skewness and kurtosis
            # 
            hdata.cmp_prnt ( hmc , 'DATA' , 'MC' , 'DATA vs MC' )
            #
            ## calculate the distances
            #
            dist = hdata.cmp_dist ( hmc , rescale = True )
            logger.info ('DATA-MC "distance"      %s' % dist )
            #
            ## calculate the 'orthogonality'
            #  
            cost = hdata.cmp_cos  ( hmc , rescale = True )
            logger.info ('DATA-MC "orthogonality" %s' % cost )
            #
            ## try to fit it DATA with MC and vice versa 
            #
            fit1 = hdata.cmp_fit ( hmc   , rescale = True )
            if fit1 and 0 == fit1.Status() :
                logger.info ( 'Fit DATA with MC   Prob=%.3g[%%] ' % ( fit1.Prob() * 100 ) )
            fit2 = hmc  .cmp_fit ( hdata , rescale = True )
            if fit2 and 0 == fit2.Status() :
                logger.info ( 'Fit MC   with DATA Prob=%.3g[%%] ' % ( fit2.Prob() * 100 ) )
            #
            ## make chi2-comparison between data and MC
            #
            c2ndf,prob = hdata.cmp_chi2 ( hmc   , rescale = True )
            logger.info ( 'DATA/MC: chi2/ndf (%.4g) and Prob %.5g%% ' % ( c2ndf , prob*100 ) )
            c2ndf,prob = hmc  .cmp_chi2 ( hdata , rescale = True )
            logger.info ( 'MC/DATA: chi2/ndf (%.4g) and Prob %.5g%% ' % ( c2ndf , prob*100 ) )
            
        if not more : 
            logger.info    ( 'No more iterations are needed #%d' % iter )
            break

        del    mcds , selector

    #
    ## final density on data 
    data_density = hdata.density()
    ## final density on mc 
    mc_density   = hmc.density()

    data_density.red  ()
    mc_density  .blue ()
    data_density.draw ('e1')
    mc_density  .draw ('e1 same')
    time.sleep ( 60 ) 

# =============================================================================
if '__main__' == __name__ :

    test_reweight() 

# =============================================================================
# The END 
# =============================================================================
