#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ============================================================================= 
# Copyright (c) Ostap developpers.
# ============================================================================= 
""" Test module for ostap/trees/cuts.py.
"""
# ============================================================================= 
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_cuts' )
else                       : logger = getLogger ( __name__        )
# ============================================================================= 
import ROOT
import ostap.trees.cuts 

def test_cuts() :

    ## inversion: 
    a = ROOT.TCut()
    assert ~a, "Empty cut must be ``zero'': ``%s''" % a
    
    a = ROOT.TCut('pt/p')
    
    logger.info ( 'a    = %s' %   a       )
    logger.info ( 'a+1  = %s' % ( a + 1 ) )
    logger.info ( 'a-1  = %s' % ( a - 1 ) )
    logger.info ( 'a*1  = %s' % ( a * 1 ) )
    logger.info ( 'a/1  = %s' % ( a / 1 ) )
    
    logger.info ( '1+a  = %s' % ( 1 + a ) )
    logger.info ( '1-a  = %s' % ( 1 - a ) )
    logger.info ( '1*a  = %s' % ( 1 * a ) ) 
    logger.info ( '1/a  = %s' % ( 1 / a ) )

    q  = ROOT.TCut(a)
    q += 1 
    logger.info ( 'a += 1 %s' % q )
    q  = ROOT.TCut(a)
    q -= 1 
    logger.info ( 'a -= 1 %s' % q )
    q  = ROOT.TCut(a)
    q *= 1 
    logger.info ( 'a *= 1 %s' % q )
    q  = ROOT.TCut(a)
    q /= 1 
    logger.info ( 'a /= 1 %s' % q )
    
    
    b = ROOT.TCut('pt>p')
    c = ROOT.TCut('x>4' )
    
    logger.info ( 'b,c  = %s, %s' %  ( b , c ) )
    logger.info ( 'b&c  = %s '    %  ( b & c ) )
    logger.info ( 'b|c  = %s '    %  ( b | c ) )

    q  = ROOT.TCut(b)
    q &= c 
    logger.info ( 'b&=c   %s '    %  q )
    q  = ROOT.TCut(b)
    q |= c 
    logger.info ( 'b|=c   %s '    %  q )

    q = 'z<y'

    logger.info ( 'b,q  = %s, %s' %  ( b , q ) )
    
    logger.info ( 'b&q  = %s '    %  ( b & q ) )
    logger.info ( 'b|q  = %s '    %  ( b | q ) )
    
    logger.info ( 'q&b  = %s '    %  ( q & b ) )
    logger.info ( 'q|b  = %s '    %  ( q | b ) )
    
    w  = ROOT.TCut(b)
    w &= q 
    logger.info ( 'b&=q   %s '    %  w )
    w  = ROOT.TCut(b)
    w |= q 
    logger.info ( 'b|=q   %s '    %  w )

    
    
    

    

    
# =============================================================================
if '__main__' == __name__ :

    test_cuts    ()
    
# =============================================================================
# The END 
# =============================================================================
