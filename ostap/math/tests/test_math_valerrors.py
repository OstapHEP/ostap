#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
## @file ostap/math/tests/test_math_valerrors.py
#  Test module for the file ostap/utils/valerrors.py
#
#  - These objects are *NOT* for math!
#  - The objects are used for graphs 
#
# ============================================================================= 
""" Test module for ostap/utils/valerrors.py
- These objects are *NOT* for math!
"""
# ============================================================================= 
from   ostap.utils.valerrors import VE , VAE , VME, AE  
import ostap.io.zipshelve    as     DBASE
import pickle 
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_math_valerros' ) 
else                       : logger = getLogger ( __name__             )
# ============================================================================= 

# =============================================================================
## test asymmetric errors:
def test_asymerr () :

    logger = getLogger( 'test_asymerr')
    logger.info ( 'Test AsymErrors' )

    
    ae1 = AE ( -0.5 ,  1.0 )
    ae2 = AE (  1.0 , -0.5 )
    ae3 = AE (  0.5 ,  1.0 )

    logger.info ( '%s %s %s' % ( ae1 , ae2 , ae3 ) )
    
    assert ae1 == ae2 , ' %s != %s ' % ( ae1 , ae2 )
    assert ae1 == ae3 , ' %s != %s ' % ( ae1 , ae3 )
    assert ae2 == ae3 , ' %s != %s ' % ( ae2 , ae3 )

    ae1_  = pickle.loads ( pickle.dumps ( ae1 ) )
    ae2_  = pickle.loads ( pickle.dumps ( ae2 ) )
    ae3_  = pickle.loads ( pickle.dumps ( ae3 ) )

    assert ae1_ == ae2 , ' %s != %s ' % ( ae1_ , ae2 )
    assert ae1_ == ae3 , ' %s != %s ' % ( ae1_ , ae3 )
    assert ae2_ == ae3 , ' %s != %s ' % ( ae2_ , ae3 )

    logger.info ( 'Scale /1 : %s %s %s' % ( ae1/1 , ae2/1 , ae3/1 ) )
    logger.info ( 'Scale *1 : %s %s %s' % ( ae1*1 , ae2*1 , ae3*1 ) )
    logger.info ( 'Scale /2 : %s %s %s' % ( ae1/2 , ae2/2 , ae3/2 ) )
    logger.info ( 'Scale *2 : %s %s %s' % ( ae1*2 , ae2*2 , ae3*2 ) )
    logger.info ( 'Scale 2* : %s %s %s' % ( 2*ae1 , 2*ae2 , 2*ae3 ) )
    
# =============================================================================
## test asymmetric errors:
def test_valasymerr () :

    logger = getLogger( 'test_valasymerr')
    logger.info ( 'Test ValWithErrors' )

    ae1 = VAE ( 10.0              , -0.5 ,  1.0 )
    ae2 = VAE ( VE ( 10.0 , 0.0 ) , -0.5 ,  1.0 )
    ae3 = VAE ( VE ( 10.0 , 0.0 ) , AE ( -0.5 , 1.0 ) )
    
    logger.info ( '%s %s %s' % ( ae1 , ae2 , ae3 ) )

    assert ae1 == ae2 , ' %s != %s ' % ( ae1 , ae2 )
    assert ae1 == ae3 , ' %s != %s ' % ( ae1 , ae3 )
    assert ae2 == ae3 , ' %s != %s ' % ( ae2 , ae3 )

    ae1_  = pickle.loads ( pickle.dumps ( ae1 ) )
    ae2_  = pickle.loads ( pickle.dumps ( ae2 ) )
    ae3_  = pickle.loads ( pickle.dumps ( ae3 ) )

    assert ae1_ == ae2 , ' %s != %s ' % ( ae1_ , ae2 )
    assert ae1_ == ae3 , ' %s != %s ' % ( ae1_ , ae3 )
    assert ae2_ == ae3 , ' %s != %s ' % ( ae2_ , ae3 )
    
    logger.info ( 'Scale /1 : %s %s %s' % ( ae1/1 , ae2/1 , ae3/1 ) )
    logger.info ( 'Scale *! : %s %s %s' % ( ae1*1 , ae2*1 , ae3*1 ) )
    logger.info ( 'Scale /2 : %s %s %s' % ( ae1/2 , ae2/2 , ae3/2 ) )
    logger.info ( 'Scale *2 : %s %s %s' % ( ae1*2 , ae2*2 , ae3*2 ) )
    logger.info ( 'Scale 2* : %s %s %s' % ( 2*ae1 , 2*ae2 , 2*ae3 ) )

    ## conversion to VE without bias 
    ve1 = ae1.asVE ( bias = False ) 
    ve2 = ae2.asVE ( bias = False ) 
    ve3 = ae3.asVE ( bias = False ) 

    logger.info ( '%s %s %s' % ( ae1 , ae2 , ae3 ) )
    logger.info ( 'asVE ( bias = False ) : %s %s %s' % ( ve1.toString ( '( %+-.3f +/- %-.3f )' ) ,
                                                         ve2.toString ( '( %+-.3f +/- %-.3f )' ) ,
                                                         ve3.toString ( '( %+-.3f +/- %-.3f )' ) ) )
    
    ## conversion to VE without bias 
    ve1 = ae1.asVE ( bias = True ) 
    ve2 = ae2.asVE ( bias = True ) 
    ve3 = ae3.asVE ( bias = True ) 

    logger.info ( 'asVE ( bias =  True ) : %s %s %s' % ( ve1.toString ( '( %+-.3f +/- %-.3f )' ) ,
                                                         ve2.toString ( '( %+-.3f +/- %-.3f )' ) ,
                                                         ve3.toString ( '( %+-.3f +/- %-.3f )' ) ) )
    
# =============================================================================
## test multiple asymmetric errors:
def test_valmulterr () :
    
    logger = getLogger( 'test_valmulterr')
    logger.info ( 'Test ValWithMultiErrors' )

    me1 = VME ( 10.0 ,     0.5 , 0.5   ,      -0.6 , 1.0   ,      0.7 , -0.3   ) 
    me2 = VME ( 10.0 , AE( 0.5 , 0.5 ) , AE ( -0.6 , 1.0 ) , AE ( 0.7 , -0.3 ) ) 
    me3 = VME ( VE ( 10 , 0.5 **2 )    , AE ( -0.6 , 1.0 ) , AE ( 0.7 , -0.3 ) ) 

    ae1 = me1.asVAE()

    logger.info ( '%s %s %s' % ( me1 , me2 , me3 ) )
    
    assert me1 == me2 , ' %s != %s ' % ( me1 , me2 )
    assert me1 == me3 , ' %s != %s ' % ( me1 , me3 )
    assert me2 == me3 , ' %s != %s ' % ( me2 , me3 )
    
    me1_  = pickle.loads ( pickle.dumps ( me1 ) )
    me2_  = pickle.loads ( pickle.dumps ( me2 ) )
    me3_  = pickle.loads ( pickle.dumps ( me3 ) )

    assert me1_ == me2 , ' %s != %s ' % ( me1_ , me2 )
    assert me1_ == me3 , ' %s != %s ' % ( me1_ , me3 )
    assert me2_ == me3 , ' %s != %s ' % ( me2_ , me3 )

    logger.info ( 'Scale /1 : %s %s %s' % ( me1/1 , me2/1 , me3/1 ) )
    logger.info ( 'Scale *1 : %s %s %s' % ( me1*1 , me2*1 , me3*1 ) )
    logger.info ( 'Scale /2 : %s %s %s' % ( me1/2 , me2/2 , me3/2 ) )
    logger.info ( 'Scale *2 : %s %s %s' % ( me1*2 , me2*2 , me3*2 ) )
    logger.info ( 'Scale 2* : %s %s %s' % ( 2*me1 , 2*me2 , 2*me3 ) )

    ae1 = me1.asVAE()
    ae2 = me2.asVAE()
    ae3 = me3.asVAE()
    
    logger.info ( '%s %s %s' % ( me1 , me2 , me3 ) )
    logger.info ( 'as VAE                : %s %s %s' % ( ae1 , ae2 , ae3 ) )
    
    ## conversion to VE without bias 
    ve1 = me1.asVE ( bias = False ) 
    ve2 = me2.asVE ( bias = False ) 
    ve3 = me3.asVE ( bias = False ) 

    logger.info ( 'asVE ( bias = False ) : %s %s %s' % ( ve1.toString ( '( %+-.3f +/- %-.3f )' ) ,
                                                         ve2.toString ( '( %+-.3f +/- %-.3f )' ) ,
                                                         ve3.toString ( '( %+-.3f +/- %-.3f )' ) ) )
    
    ## conversion to VE without bias 
    ve1 = me1.asVE ( bias = True ) 
    ve2 = me2.asVE ( bias = True ) 
    ve3 = me3.asVE ( bias = True ) 

    logger.info ( 'asVE ( bias =  True ) : %s %s %s' % ( ve1.toString ( '( %+-.3f +/- %-.3f )' ) ,
                                                         ve2.toString ( '( %+-.3f +/- %-.3f )' ) ,
                                                         ve3.toString ( '( %+-.3f +/- %-.3f )' ) ) )

# =============================================================================
if '__main__' == __name__ :

    test_asymerr    ()
    test_valasymerr ()
    test_valmulterr ()

# =============================================================================
##                                                                      The END 
# =============================================================================
