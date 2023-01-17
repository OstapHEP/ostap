#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_roostats.py
# Test module for RooStats 
# ============================================================================= 
""" Test module for RooStats 
"""
# ============================================================================= 
from   __future__               import print_function
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ostap.fitting.roofit 
import ostap.fitting.models     as     Models 
from   ostap.core.core          import cpp, VE, dsID, hID , rooSilent, Ostap 
from   ostap.utils.timing       import timing
from   builtins                 import range
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait
import ostap.logger.table       as     T
from   ostap.fitting.roostats   import ( interval_PL, interval_FC   ,
                                         interval_BC, interval_MCMC )    
import ROOT

# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_roostats' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================


x       = ROOT.RooRealVar('x','x-variable', -10 , 10 )
sigma   = ROOT.RooFit.RooConst( 1 ) 
pdf     = Models.Gauss_pdf    ( 'Gauss' , x , mean = ( 0 , -1 , 1 ) , sigma  = sigma  ) 

dataset = pdf.generate ( 100 )

with use_canvas ( "test_fitting_roostats" ) : 
    result , frame = pdf.fitTo ( dataset , draw = True , nbins = 20 )

# ============================================================================
def test_limits () :

    logger = getLogger("test_limits")

    logger.info ( "Test Limits with RooStats" )

    
    with use_canvas ( "test_limits" ) : 
        result , frame = pdf.fitTo ( dataset , draw = True , nbins = 20 )

    rows = [  ( 'Method' , 'Interval [90%]' , 'Lower [95%]' , 'Upper[95%]' , 'CPU' ) ]
    
    with timing ( "Profile Likelihood" , logger ) as t :
        low, high = interval_PL ( pdf  , pdf.mean , 0.90 , dataset )
        lower     = interval_PL ( pdf  , pdf.mean , 0.95 , dataset , limit = -1 )
        upper     = interval_PL ( pdf  , pdf.mean , 0.95 , dataset , limit = +1 )

    row = 'Profile likelihood' , \
          '[%-+.3f,%+.3f]' %  ( low , high ) , \
          '%+.3f' %  lower , '%+.3f' %  upper , '%.1f' % t.delta
    rows.append ( row )
    
    with use_canvas ( "test_limits" ) : 
        result , frame = pdf.fitTo ( dataset , draw = True , nbins = 20 )
        
    with timing ( "Feldman-Cousins" , logger ) as t :
        low, high = interval_FC ( pdf  , pdf.mean , 0.90 , dataset )
        lower     = interval_FC ( pdf  , pdf.mean , 0.95 , dataset , limit = -1 )
        upper     = interval_FC ( pdf  , pdf.mean , 0.95 , dataset , limit = +1 )
        
    row = 'Feldman-Cousin' , \
          '[%-+.3f,%+.3f]' %  ( low , high ) , \
          '%+.3f' %  lower , '%+.3f' %  upper , '%.1f' % t.delta
    rows.append ( row )
    
    with use_canvas ( "test_limits" ) : 
        result , frame = pdf.fitTo ( dataset , draw = True , nbins = 20 )
        
    with timing ( "Bayesian" , logger ) as t :
        low, high = interval_BC ( pdf  , pdf.mean , 0.90 , dataset )
        lower     = interval_BC ( pdf  , pdf.mean , 0.95 , dataset , limit = -1 )
        upper     = interval_BC ( pdf  , pdf.mean , 0.95 , dataset , limit = +1 )
        
    row = 'Bayesian' , \
          '[%-+.3f,%+.3f]' %  ( low , high ) , \
          '%+.3f' %  lower , '%+.3f' %  upper , '%.1f' % t.delta
    rows.append ( row )
    
    with use_canvas ( "test_limits" ) : 
        result , frame = pdf.fitTo ( dataset , draw = True , nbins = 20 )
        
    with timing ( "MCMC" , logger ) as t :
        low, high = interval_MCMC ( pdf  , pdf.mean , 0.90 , dataset )
        lower     = interval_MCMC ( pdf  , pdf.mean , 0.95 , dataset , limit = -1 )
        upper     = interval_MCMC ( pdf  , pdf.mean , 0.95 , dataset , limit = +1 )
        
    row = 'MCMC' , \
          '[%-+.3f,%+.3f]' %  ( low , high ) , \
          '%+.3f' %  lower , '%+.3f' %  upper , '%.1f' % t.delta
    rows.append ( row )

    title = 'Intervals & Limits'
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'rcccr' )
    logger.info ( '%s\n%s' % ( title , table ) )
    
    
          
    
# =============================================================================
if '__main__' == __name__ :

    test_limits () 

# =============================================================================
##                                                                      The END 
# =============================================================================
