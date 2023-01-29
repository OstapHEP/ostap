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
from   ostap.fitting.roostats   import ( ProfileLikelihoodInterval ,
                                         FeldmanCousinsInterval    ,
                                         BayesianInterval          ,
                                         MCMCInterval              )
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

    
    with timing ( "Profile Likelihood", logger ) as t :
        pli = ProfileLikelihoodInterval ( pdf     = pdf      ,
                                          poi     = pdf.mean , 
                                          dataset = dataset  )
        low, high = pli.interval    ( 0.90 , dataset )
        upper     = pli.upper_limit ( 0.95 , dataset )
        lower     = pli.lower_limit ( 0.95 , dataset )
        
    row = 'Profile likelihood' , \
          '[%-+.3f,%+.3f]' %  ( low , high ) , \
          '%+.3f' %  lower , '%+.3f' %  upper , '%.1f' % t.delta
    rows.append ( row )

    with use_canvas  ( 'Profile Likelihood plot ' , wait = 1 ) :
        pli_plot = pli.plot ()
        if pli_plot : pli_plot.draw () 

    with use_canvas ( "test_limits" ) : 
        result , frame = pdf.fitTo ( dataset , draw = True , nbins = 20 )
        
    with timing ( 'Feldman-Cousins [~30",silent]'   , logger ) as t :
        with rooSilent ( ROOT.RooFit.WARNING ) : 
            fci = FeldmanCousinsInterval ( pdf     = pdf      ,
                                           poi     = pdf.mean , 
                                           dataset = dataset  )  
            low, high = fci.interval    ( 0.90 , dataset )
            upper     = fci.upper_limit ( 0.95 , dataset )
            lower     = fci.lower_limit ( 0.95 , dataset )
            
    row = 'Feldman-Cousins' , \
          '[%-+.3f,%+.3f]' %  ( low , high ) , \
          '%+.3f' %  lower , '%+.3f' %  upper , '%.1f' % t.delta
    rows.append ( row )

    with use_canvas  ( 'Felfman-Cousins plot ' , wait = 1 ) :
        fci_plot = fci.plot ()
        if fci_plot : fci_plot.draw ('ap') 

    with use_canvas ( "test_limits" ) : 
        result , frame = pdf.fitTo ( dataset , draw = True , nbins = 20 )
        
    with timing ( 'Bayesian [~60"]' , logger ) as t :        
        bci = BayesianInterval ( pdf     = pdf      ,
                                 poi     = pdf.mean , 
                                 dataset = dataset  ) 
        low, high = bci.interval    ( 0.90 , dataset )
        upper     = bci.upper_limit ( 0.95 , dataset )
        lower     = bci.lower_limit ( 0.95 , dataset )
        
    row = 'Bayesian' , \
          '[%-+.3f,%+.3f]' %  ( low , high ) , \
          '%+.3f' %  lower , '%+.3f' %  upper , '%.1f' % t.delta
    rows.append ( row )
    
    with use_canvas  ( 'Bayesian  plot ' , wait = 1 ) :
        bci_plot = bci.plot ()
        if bci_plot : bci_plot.draw () 

    with timing ( "Markov Chain MC (~6')" , logger ) as t :
        mci = MCMCInterval ( pdf     = pdf      ,
                             poi     = pdf.mean , 
                             dataset = dataset  )          
        low, high = mci.interval    ( 0.90 , dataset )
        upper     = mci.upper_limit ( 0.95 , dataset )
        lower     = mci.lower_limit ( 0.95 , dataset )


    row = 'Markov Chain MC' , \
          '[%-+.3f,%+.3f]' %  ( low , high ) , \
          '%+.3f' %  lower , '%+.3f' %  upper , '%.1f' % t.delta
    rows.append ( row )

    with use_canvas  ( 'Markov Chain MC plot ' , wait = 1 ) :
        mci_plot = mci.plot ()
        if mci_plot : mci_plot.draw () 

    with use_canvas  ( 'Visualise CL intervals' , wait = 2 ) as cnv :
        
        cnv.Divide(2,2)
        
        cnv.cd(1) 
        pli_plot = pli.plot ()
        if pli_plot : pli_plot.draw () 
        
        cnv.cd(2) 
        fci_plot = fci.plot ()
        if fci_plot : fci_plot.draw ('ap') 
        
        cnv.cd(3) 
        bci_plot = bci.plot ()
        if bci_plot : bci_plot.draw () 
        
        cnv.cd(4) 
        mci_plot = mci.plot ()
        if mci_plot : mci_plot.draw () 
        
        title = 'Intervals & Limits'
        table = T.table ( rows , title = title , prefix = '# ' , alignment = 'rcccr' )
        logger.info ( '%s\n%s' % ( title , table ) )

    del mci


# =============================================================================
if '__main__' == __name__ :

    test_limits () 

# =============================================================================
##                                                                      The END 
# =============================================================================
