#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/core/tests/test_core_core.py
# Copyright (c) Ostap developpers.
# =============================================================================
""" Test module or soem core utils 
"""
# =============================================================================
from   ostap.core.core   import *
from   ostap.utils.utils import num_fds
from   ostap.logger.mute import mute
import ROOT, warnings  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_core_core'  )
else                       : logger = getLogger ( __name__          )
# =============================================================================

def test_core_core1() :

    logger = getLogger ( 'test_core_core1' ) 

    with ROOTCWD () :

        for f in ( fID  , funID    ,
                   hID  , histoID  ,
                   dsID , rootID   ) :    
            logger.info  ('ID/%20s : %s' % ( f.__name__ , f() ) )

        for i in ( 0.0 , 1.e-125 , 100 , 2**36 )  :

            logger.info ( '%14s/%-6s is-zero  %s' % ( i , type(i).__name__ , iszero         ( i ) ) )
            logger.info ( '%14s/%-6s is-int   %s' % ( i , type(i).__name__ , isint          ( i ) ) )
            logger.info ( '%14s/%-6s is-long  %s' % ( i , type(i).__name__ , islong         ( i ) ) )

            v = VE ( i , i ) 
            logger.info ( '%14s/%-6s n-entry  %s' % ( i , type(i).__name__ , natural_entry  ( v ) ) ) 
            logger.info ( '%14s/%-6s n-number %s' % ( i , type(i).__name__ , natural_number ( v ) ) ) 
            
        total    = VE ( 100 , 100 )
        accepted = VE (  10 ,  10 )
        rejected = VE (  90 ,  90 )

        logger.info ( 'zechEff         %s%%' % ( zechEff         ( accepted , total    ) * 100 ) ) 
        logger.info ( 'binomEff        %s%%' % ( binomEff        (       10 , 100      ) * 100 ) ) 
        logger.info ( 'binomEff2       %s%%' % ( binomEff2       ( accepted , rejected ) * 100 ) ) 
        logger.info ( 'wilsonEff       %s%%' % ( wilsonEff       (       10 , 100      ) * 100 ) ) 
        logger.info ( 'agrestiCoullEff %s%%' % ( agrestiCoullEff (       10 , 100      ) * 100 ) ) 
        
        a = strings ( 'a','b','c' )
        logger.info ( 'Strings  :%s' % a )
        a = split_string ( "a,b;c,b:e" )        
        logger.info ( 'Split    :%s' % str ( a ) ) 
        a = split_string ( "atan2(a,b);atan2(c,b):e" )        
        logger.info ( 'Split    :%s' % str ( a ) ) 
        

def test_core_core2() :
    
    logger = getLogger ( 'test_core_core2' )
    
    before = num_fds() 
    logger.info ( '#open file descriptors before is %d' % before )
    dummy  = 0 
    for i in range ( 10000  ) :
        with mute ( True , True ) :
            dummy += 1 
        
    after = num_fds() 
    logger.info ( '#open file descriptors after  is %d' % after ) 


def test_core_core3() :
    
    logger = getLogger ( 'test_core_core3' )

    logger.error ("A ROOT error occurs just after:" )  
    ROOT.Error   ( 'test_core_core3' , 'fictive ROOT Error/1'    )

    try : 
        with rootException()  : 
            logger.error ("A ROOT error occurs just after:" )  
            ROOT.Error   ( 'test_core_core3' , 'fictive ROOT Error/2'    )
    except :
        logger.exception ( "An exception is here (here it is perfectly OK)" )
        
    logger.error   ("A ROOT error occurs just after:" )  
    ROOT.Error     ( 'test_core_core3' , 'fictive ROOT Error/3'    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")                

        logger.warning ("A ROOT warning occurs just after:" )  
        ROOT.Warning   ( 'test_core_core3' , 'fictive ROOT Warning/1'  )
        
        with rootException()  : 
            logger.warning ("A ROOT warning occurs just after:" )  
            ROOT.Warning   ( 'test_core_core3' , 'fictive ROOT Warning/2'  )
            
        logger.warning ("A ROOT warning occurs just after:" )  
        ROOT.Warning ( 'test_core_core3' , 'fictive ROOT Warning/3'  )
            
# =============================================================================
if '__main__' == __name__ :

    test_core_core1 () 
    test_core_core2 () 
    test_core_core3 () 

# =============================================================================
##                                                                      The END 
# =============================================================================
