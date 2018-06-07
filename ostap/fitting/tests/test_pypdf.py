#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_pydf.py
# Test module for ostap/fitting/pypdf.py
# - It tests PyPDF construction
# ============================================================================= 
""" Test module for ostap/fitting/pypdf.py
- It tests ``pure-python'' PDF 
"""
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ROOT, random, math
import ostap.fitting.roofit 
from   ostap.core.core      import VE, dsID
from   ostap.fitting.basic  import MASS
from   ostap.fitting.pypdf  import PyPDF
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_models' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## make simple test mass 
mass     = ROOT.RooRealVar ( 'test_mass' , 'Some test mass' , 3.0 , 3.2 )

## book very simple data set
varset0   = ROOT.RooArgSet  ( mass )
dataset   = ROOT.RooDataSet ( dsID() , 'Test Data set-0' , varset0 )  

mmin,mmax = mass.minmax()

## fill it 
m = VE(3.100,0.015**2)
for i in xrange(0,500) :
    mass.value = m.gauss () 
    dataset.add ( varset0 )

logger.info ('DATASET %s' % dataset )

# =============================================================================
## @class PyGauss
#  local ``pure-python'' PDF 
class PyGauss(MASS,PyPDF) :
    """Local ``pure-python'' PDF """    
    norm = 1.0 / math.sqrt ( 2 * math.pi )
    def __init__ ( self         ,
                   name         ,
                   xvar         ,
                   mean  = (         3.05  , 5.15  ) ,
                   sigma = ( 0.010 , 0.005 , 0.025 ) ,
                   pdf   = None ) :
        
        MASS .__init__ ( self , name      , xvar , mean , sigma ) 
        PyPDF.__init__ ( self , self.name , ( self.xvar  ,
                                              self.mean  ,
                                              self.sigma ) , pdf = pdf )
        self.config = {
            'name'  : self.name ,
            'xvar'  : self.xvar ,
            'mean'  : self.mean ,
            'sigma' : self.mean ,
            'pdf'   : None        ## attention! 
            }        
    ## the  main method 
    def evaluate ( self ) :
        varlist = self.varlist        
        x = float ( varlist [ 0 ] ) 
        m = float ( varlist [ 1 ] ) 
        s = float ( varlist [ 2 ] )        
        dx = ( x - m ) / s        
        return math.exp ( -0.5 * dx * dx ) * self.norm / s 

                
# =============================================================================
## pygauss PDF
# =============================================================================
def test_pygauss() :
    
    logger.info ('Test PyGauss:  simple Gaussian signal' )
    
    gauss = PyGauss( 'PyGauss' , xvar = mass )
    
    result = gauss.fitTo ( dataset , draw = True )

    print result 
    
# =============================================================================
if '__main__' == __name__ :

    test_pygauss          () ## simple Gaussian PDF
    pass

# =============================================================================
# The END 
# ============================================================================= 
