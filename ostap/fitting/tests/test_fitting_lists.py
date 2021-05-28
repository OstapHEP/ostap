#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
# @file test_lists.py
# Test module for ostap/fitting/roofit.py
# - It tests the operations with lists
# ============================================================================= 
"""# Test module for ostap/fitting/roofit.py
# - It tests the operations with lists
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import 
# ============================================================================= 
import ROOT, random
import ostap.fitting.roofit 
from   ostap.logger.utils   import rooSilent 
from   ostap.utils.timing   import timing 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'ostap/fitting/tests/test_lists' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
logger.info ( 'Test for operations with lists')
# =============================================================================
a  = ROOT.RooRealVar  ( 'a' , 'a' , -10 , 10 )
b  = ROOT.RooRealVar  ( 'b' , 'b' , -10      )
c  = ROOT.RooConstVar ( 'c' , 'c' ,  1       )
d  = ROOT.RooConstVar ( 'd' , 'd' , -1       ) 

l1 = ROOT.RooArgList    ( a , b )
l2 = ROOT.RooArgSet     ( b , c )
l3 = ROOT.RooLinkedList ()
l3.add ( a )
l3.add ( c )

def test_len () :

    logger.info ('Len l1 %d' % len ( l1 ) )
    logger.info ('Len l2 %d' % len ( l2 ) )
    logger.info ('Len l3 %d' % len ( l3 ) )

def test_loop () :

    logger.info ( 'L1: %s ' % [ a  for a in l1 ] )
    logger.info ( 'L2: %s ' % [ a  for a in l2 ] )
    logger.info ( 'L3: %s ' % [ a  for a in l3 ] )
    
def test_contains () :

    for l in ( l1 , l2 , l3 ) :
        logger.info ( ' a in l/"a" in l? %s/%s' % ( a in l , 'a' in l ) )
        logger.info ( ' b in l/"b" in l? %s/%s' % ( b in l , 'b' in l ) )
        logger.info ( ' c in l/"c" in l? %s/%s' % ( c in l , 'c' in l ) )

def test_getitem () :

    logger.info( 'l1[0]  /l1[1]    %s/%s ' % ( l1[0]   , l1[1]   ) )
    logger.info( 'l2["b"]/l2["c"]  %s/%s ' % ( l2["b"] , l2["b"] ) )

def test_sums () :

    logger.info ( 'l1+l1 : %s'  % ( l1 + l1 ) )
    logger.info ( 'l1+l2 : %s'  % ( l1 + l2 ) )
    logger.info ( 'l1+l3 : %s'  % ( l1 + l3 ) )

    logger.info ( 'l2+l1 : %s'  % ( l2 + l1 ) )
    logger.info ( 'l2+l2 : %s'  % ( l2 + l2 ) )
    logger.info ( 'l2+l3 : %s'  % ( l2 + l3 ) )
    
    logger.info ( 'l3+l1 : %s'  % ( l3 + l1 ) )
    logger.info ( 'l3+l2 : %s'  % ( l3 + l2 ) )
    logger.info ( 'l3+l3 : %s'  % ( l3 + l3 ) )

def test_add () :

    logger.info ( 'l1+c  : %s'  % ( l1 + c  ) )
    logger.info ( 'l2+c  : %s'  % ( l2 + c  ) )
    logger.info ( 'l3+c  : %s'  % ( l3 + c  ) )
    
    logger.info ( 'l1+d  : %s'  % ( l1 + d  ) )
    logger.info ( 'l2+d  : %s'  % ( l2 + d  ) )
    logger.info ( 'l3+d  : %s'  % ( l3 + d  ) )

    ## logger.info ( 'c+l1  : %s'  % ( c  + l1 ) )
    ## logger.info ( 'c+l2  : %s'  % ( c  + l2 ) )
    ## logger.info ( 'c+l3  : %s'  % ( c  + l3 ) )

    ## logger.info ( 'd+l1  : %s'  % ( d  + l1 ) )
    ## logger.info ( 'd+l2  : %s'  % ( d  + l2 ) )
    ## logger.info ( 'd+l3  : %s'  % ( d  + l3 ) )
    
    
# =============================================================================
if '__main__' == __name__ :

    with timing ( "len"       , logger ) : 
        test_len      ()
        
    with timing ( "contains"  , logger ) : 
        test_contains ()
        
    with timing ( "getitem"   , logger ) : 
        test_getitem  ()
        
    with timing ( "sums"      , logger ) : 
        test_sums     ()
        
    with timing ( "add"       , logger ) : 
        test_add      ()
    
# =============================================================================
##                                                                      The END 
# ============================================================================= 
