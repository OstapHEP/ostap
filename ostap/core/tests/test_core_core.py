#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/core/tests/test_core_core.py
# Copyright (c) Ostap developpers.
# =============================================================================
""" Test module or soem core utils 
"""
# =============================================================================
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_core_core'  )
else                       : logger = getLogger ( __name__          )
# =============================================================================
from   ostap.core.core       import *



def test_core_core1() :


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
        logger.info ( 'Split    :%s' % a )
        
        
# =============================================================================
if '__main__' == __name__ :

    test_core_core1() 


# =============================================================================
##                                                                      The END 
# =============================================================================
