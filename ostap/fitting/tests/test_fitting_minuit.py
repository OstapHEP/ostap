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
from   ostap.utils.timing     import timing
from   ostap.utils.root_utils import batch_env 
import ostap.fitting.minuit
import ROOT, random 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_minuit' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================

# ==============================================================================
## *NEW* PyROOT  uses <code>ctypes.c_int</code>, <code>ctypes.c_double</code>
def fcn ( npar ,  gin , f , par , iflag ) :
    """ Trivial FCN
    """
    n = int(npar.value)          ## ATTENTION! 
    
    s = 0   
    for i in range ( n ) :
        p  = par [ i ]        
        s += ( p - 10*i ) ** 2 / ( i + 1 ) ** 2 
        
    f.value = s                  ## ATTENTION! 

# =============================================================================
def test_minuit ( ) :
    
    logger = getLogger("test_minuit")
    
    minuit = ROOT.TMinuit( 5 )
    minuit.SetFCN ( fcn  )
    
    minuit.execute ( "SET ERR" , 1    )

    ## define the papameters: name , init , step , low and high edges 
    minuit.addpar  ( 'p1' ,  1  , 0.01 )
    minuit.addpar  ( 'p2' ,  0  , 0.01 )
    minuit.addpar  ( 'p3' , -1  , 0.01 )
    minuit.addpar  ( 'p4' ,  35 , 0.01 , low = -20 , high = +100 )

    minuit.execute ( "SHOW PAR" )
    minuit.show    ( "PAR" )  ## ditto 
    minuit.show    ( )        ## ditto 
    minuit.show    ( )        ## ditto 
    
    logger.info    ( "MIGRAD: %s"              % minuit.migrad () ) 
    logger.info    ( "HESSE : %s"              % minuit.hesse  () ) 
    logger.info    ( "Covariance  matrix:\n%s" % minuit.cov    () )
    logger.info    ( "Correlation matrix:\n%s" % minuit.cor    () )

    ## set new value for parameter
    minuit['p1'] = 8

    minuit.fix     ( 'p2' , 0 )  ##   fix parameter 
    minuit.migrad  ()
    
    minuit.release ( 'p2')        ##   release it
    
    minuit.migrad  () 
    minuit.migrad  () 
    minuit.minos   ( 'p1' , 'p2' , 'p3' , 'p4' ) 

    ## loop over all parameters
    for i in minuit :
        
        p    = minuit [ i ]  ## get parameter value 
        pname = minuit.par_name ( i )
        
        line = "Parameter %d/%s is %s" % ( i , pname , p ) 
        ml , mh = minuit.minErr ( i )
        if ml < mh  : line += "; minos errors are (%+g,%+g)" % ( -ml , mh )
        logger.info ( line )
        
    logger.info ('Fit results:\n%s' % minuit.table ( prefix = "# " ,
                                                     title  = "TMinuit fit" ) )
# =============================================================================
if  '__main__' == __name__ :

    with timing ( "minuit" , logger ) : 
        test_minuit () 


# =============================================================================
##                                                                      The END 
# =============================================================================
