#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_components2.py
# Test module 
# ============================================================================= 
""" Test module 
- It tests sPlotting
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.core       import cpp, VE, dsID, rooSilent 
from   ostap.utils.timing    import timing 
from   ostap.plotting.canvas import use_canvas 
from   ostap.utils.utils     import batch_env 
import ostap.fitting.models  as     Models 
import ostap.fitting.roofit 
import ROOT, random
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_components2' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch form environment 
batch_env ( logger )
# =============================================================================
#
models = set() 

def test_components_2 () :

    logger = getLogger ( 'test_components_2' )
    
    ## make simple test mass 
    mass    = ROOT.RooRealVar ( 'test_mass' , 'Some test mass' , 0 , 10 )
    
    ## book very simple data set
    varset  = ROOT.RooArgSet  ( mass )
    dataset = ROOT.RooDataSet ( dsID() , 'Test Data set' , varset )  

    mmin, mmax = mass.minmax()

    ##  two gaussinan signal cmppnent
    N1   = 5000
    m1   = VE(5,0.2**2)
    
    N2   = 5000
    m2   = VE(3,0.3**2)

    ## background exponential components 
    taub = 5.0
    NB   = 5000
    
    ## generate 1st signal 
    for i in range(0,N1) :
        m = m1.gauss ()
        while not mmin < m < mmax : m =  m1.gauss() 
        mass.value = m
        dataset.add  ( varset )

    ## generate 2nd signal 
    for i in range(0,N2) :
        m = m2.gauss ()
        while not mmin < m < mmax : m =  m2.gauss() 
        mass.value = m
        dataset.add  ( varset )
        
    ##  generate background  
    for i in range ( 0,NB) : 
        m = random.expovariate ( 1./ taub )
        while not mmin < m < mmax : m = random.expovariate ( 1./ taub )
        mass.value = m
        dataset.add ( varset )
        
    logger.info ( "Original dataset\n%s" % dataset.table ( prefix = '# ' ) )
        
    signal1 = Models.Gauss_pdf ( 'G1'  , xvar = mass  , mean = ( 5 , 4, 6 ) , sigma = ( 0.2 , 0.1 , 0.3 ) )
    signal2 = Models.Gauss_pdf ( 'G2'  , xvar = mass  , mean = ( 3 , 2, 4 ) , sigma = ( 0.3 , 0.2 , 0.4 ) )

    S1    = ROOT.RooRealVar       ( 'S1' , '1st signal yeild'  , N1          ,  1   , 10000 )  
    ratio = ROOT.RooRealVar       ( 'R'  , 'ratio of S2 to S1' , 1.0*N2 / N1 , 0.05 , 100   )
    S2    = signal1.vars_multiply ( S1 , ratio , name = 'S2' , title = '2ns signal yield'   )

    model = Models.Fit1D( signals = [ signal1 , signal2 ] , background = 1 , S = ( S1 , S2 ) ) 
    
    model.B = NB

    signal1.mean .fix ( m1.value () )
    signal1.sigma.fix ( m1.error () )
    signal2.mean .fix ( m2.value () )
    signal2.sigma.fix ( m2.error () )
    model.fitTo ( dataset , silent = True )
    signal1.mean .release() 
    signal1.sigma.release() 
    signal2.mean .release() 
    signal2.sigma.release() 

    with use_canvas ( 'test_components_2' ) :
        
        model.fitTo ( dataset , silent = True )
        r , f = model.fitTo ( dataset , silent = True , draw = True , nbins = 50 )
        logger.info ( "Mass fit : fit results\n%s" % r.table ( title = 'Mass fit' , prefix = '# ' ) )
        
        r , f = model.fitTo ( dataset , silent = True , draw = True , nbins = 50 , minos = ('R','S1') )
        logger.info ( "Mass fit : fit results\n%s" % r.table ( title = 'Mass fit' , prefix = '# ' ) )

    models.add ( signal1 )
    models.add ( signal2 )
    models.add ( model   )
    
# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    logger = getLogger ( 'test_db' ) 
    logger.info ( 'Saving all objects into DBASE' )
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db :
        for m in models :
            db['model:' + m.name ] = m
            db['roo_tot:%s' % m.name ] = m.pdf
            for i,s in enumerate ( m.signals ) :
                db['roo_sig%d:%s' % ( i , m.name ) ] = s
            for i, b in enumerate ( m.backgrounds ) : 
                db['roo_bkg%d:%s' % ( i , m.name ) ] = s
            for a in m.alist1 : 
                db['cmp:%s/%s' % ( m.name , a.name ) ] = a
        db['models'   ] = models
        db.ls() 

# =============================================================================

if '__main__' == __name__ :

    with timing ("Components2" , logger ) : 
        test_components_2   () 

    ## check finally that everything is serializeable:
    with timing ('test_db'             , logger ) :
        test_db ()

# =============================================================================
##                                                                      The END 
# =============================================================================
