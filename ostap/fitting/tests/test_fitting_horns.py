#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_horns.py
# test HORSNdini_pdf, HILLdini_pdf and HHdini_pdf
# ============================================================================= 
"""Test for HORSNdini_pdf, HILLdini_pdf and HHdini_pdf
"""
# ============================================================================= 
from   __future__              import print_function
# ============================================================================= 
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models 
from   ostap.core.core          import Ostap, VE, dsID, rooSilent,rootError 
from   builtins                 import range
from   ostap.utils.timing       import timing 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait 
import ROOT, random, warnings 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_horns' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## make simple test mass
x = ROOT.RooRealVar( 'x' , 'x-observable' , 0 , 12 )

models = set()

def test_horns_1 () :
    
    logger = getLogger( "test_horns_1" ) 
    
    horns  = Models.Fit1D (
        signal = Models.HORNSdini_pdf  ( 'Horns' ,
                                         xvar = x  ,
                                         a    =  ( 7 , 6 , 8 ) ,
                                         delta = ( 1 , 1 , 2 ) ,
                                         resolution = 0.1 ) ,  
        S = 1000 ,
        B = 10   ,
        suffix = 'horns'
        )
    
    horns.signal.a    .fix(7)
    horns.signal.delta.fix(1) 
    dataset = horns.generate ( 5000 ) 
    with wait ( 1 ) , use_canvas ( 'Test HORNSdini' ) :
        horns.S = 5000
        horns.B =   50
        result , _ = horns.fitTo ( dataset , draw = False , silent = True , minimizer = 'minuit' ) 
        result , _ = horns.fitTo ( dataset , draw = True  , silent = True , minimiser = 'minuit' ) 


    table  = result.table ( prefix = '# ' ) 

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning ('Fit result\n%s' % table ) 
    else :
        logger.info    ('Fit result\n%s' % table ) 

    models.add ( horns  )

# =============================================================================
def test_hill_1 () :

    logger = getLogger( "test_hill_1" ) 
    
    horns  = Models.Fit1D (
        signal = Models.HILLdini_pdf  ( 'Hill' ,
                                        xvar = x  ,
                                         a    =  ( 7 , 6 , 8 ) ,
                                         delta = ( 1 , 1 , 2 ) ,
                                         resolution = 0.1 ) ,  
        S = 1000 ,
        B = 10   ,
        suffix = 'hill'
        )

    horns.signal.a    .fix(7)
    horns.signal.delta.fix(1) 
    dataset = horns.generate ( 5000 ) 
    with wait ( 1 ) , use_canvas ( 'Test HILLdini' ) : 
        horns.S = 5000
        horns.B =   50
        result , _ = horns.fitTo ( dataset , draw = False , silent = True , minimiser = 'minuit' ) 
        result , _ = horns.fitTo ( dataset , draw = True  , silent = True , minimiser = 'minuit' ) 


    table  = result.table ( prefix = '# ' ) 

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning ('Fit result\n%s' % table ) 
    else :
        logger.info    ('Fit result\n%s' % table ) 

    models.add ( horns  )


# =============================================================================
def test_hh_1 () :
    
    logger = getLogger( "test_hh_1" ) 
    
    horns  = Models.Fit1D (
        signal = Models.HHdini_pdf  ( 'HH' ,
                                      xvar = x  ,
                                      a    =  ( 7   , 6 , 8 ) ,
                                      fL   =  ( 0.5 , 0 , 1 ) , 
                                      delta = ( 1 , 1 , 2 ) ,
                                      resolution = 0.1 ) ,  
        S = 1000 ,
        B = 10   ,
        suffix = 'hill'
        )

    horns.signal.a.    fix(7)
    horns.signal.delta.fix(1) 
    horns.signal.fL.fix(0.4)
    dataset = horns.generate ( 5000 ) 
    with wait ( 1 ) , use_canvas ( 'Test HILLdini' ) :
        horns.S = 5000
        horns.B =   50
        result , _ = horns.fitTo ( dataset , draw = False , silent = True , minimiser = 'minuit' )
        horns.signal.fL.release()         
        result , _ = horns.fitTo ( dataset , draw = True  , silent = True , minimiser = 'minuit' ) 

    table  = result.table ( prefix = '# ' ) 

    if 0 != result.status() or 3 != result.covQual() :
        logger.warning ('Fit result\n%s' % table ) 
    else :
        logger.info    ('Fit result\n%s' % table ) 

    models.add ( horns  )


# =============================================================================
## check that everything is serializable
# =============================================================================
def test_db() :

    logger = getLogger ( 'test_db' )
    
    logger.info('Saving all objects into DBASE')
    import ostap.io.zipshelve   as     DBASE
    from ostap.utils.timing     import timing 
    with timing( 'Save everything to DBASE', logger ), DBASE.tmpdb() as db : 
        db['x'       ] = x
        db['models'  ] = models
        for m in models : db['model:%s' % m.name ] = m 
        db.ls() 
 
# =============================================================================
if '__main__' == __name__ :
    
    with timing (  "HORNS-1"  , logger ) : 
        test_horns_1           ()          

    with timing (  "HILL-1"  , logger ) : 
        test_hill_1           ()          
        
    with timing (  "HH-1"  , logger ) : 
        test_hh_1           ()          
            
    ## check finally that everything is serializeable:
    with timing (  "Save to DB"     , logger ) : 
        test_db           ()          
        
    pass

# =============================================================================
##                                                                      The END
# =============================================================================
