# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_minuit.py
# Test module for ostap/fitting/minuit.py
# - It tests decorations for <code>ROOT.TMinuit</code>
# @see TMinuit 
# ============================================================================= 
""" Test module for ostap/fitting/minuit.py
- It tests decorations for <code>ROOT.TMinuit</code>
"""
# ============================================================================= 
from   __future__        import print_function
# ============================================================================= 
import ROOT, random 
import ostap.fitting.minuit
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_minuit' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## trivial FCN 
def fcn ( npar ,  gin , f , par , iflag ) :
    """Trivial FCN
    """
    
    s = 0
    n = npar[0]
    
    for i in range ( n ) :
        p  = par[i]        
        s += ( p - 10*i ) ** 2 
    
    f[0] = s 

# =============================================================================
def test_minuit ( ) :
    
    minuit = ROOT.TMinuit( 5 )
    minuit.SetFCN ( fcn  )
    
    minuit.execute ( "SET ERR" , 1    )

    ## define the papameters: name , init , step , low and high edges 
    minuit.addpar  ( 'p1' ,  1 , 0.01 )
    minuit.addpar  ( 'p2' ,  0 , 0.01 )
    minuit.addpar  ( 'p3' , -1 , 0.01 )
    minuit.addpar  ( 'p4' ,  4 , 0.01 )

    minuit.execute ( "SHOW PAR" )
    ## minuit.show    ( "PAR" )  ## ditto 
    ## minuit.show    ( )        ## ditto 
    
    logger.info    ( "MIGRAD: %s"              % minuit.migrad () ) 
    logger.info    ( "HESSE : %s"              % minuit.hesse  () ) 
    logger.info    ( "Covariance  matrix:\n%s" % minuit.cov    () )
    logger.info    ( "Correlation matrix:\n%s" % minuit.cor    () )

    ## set new value for parameter
    minuit[3] = 8

    minuit.fix ( 0 , 0 )  ##   fix parameter 
    minuit.migrad  ()
    
    minuit.release ( 0 )  ##   release it
    
    minuit.migrad  () 
    minuit.minos   () 
    minuit.show    ()

    ## loop over all parameters
    for i in minuit :
        
        p = minuit[i]  ## get parameter
        
        line = "Parameter %d is  %s" % ( i , p ) 
        ml , mh = minuit.minErr ( i )
        if ml < mh  : line += logger.info ( "; minos errors are (%+g,%+g)" % (-ml,mh) )
        logger.info ( line )

        
# =============================================================================
if  '__main__' == __name__ :

    test_minuit () 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
