#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright (c) Ostap developers.
# ============================================================================= 
# @file test_fitting_pydf.py
# Test module for ostap/fitting/pypdf.py
# - It tests PyPDF construction
# ============================================================================= 
""" Test module for ostap/fitting/pypdf.py
- It tests ``pure-python'' PDF 
"""
# ============================================================================= 
from   __future__        import print_function
# ============================================================================= 
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
import ROOT, random, math
import ostap.fitting.roofit 
from   ostap.core.core      import VE, dsID, Ostap
from   ostap.fitting.basic  import MASS,     Fit1D , Generic1D_pdf 
from   ostap.fitting.pypdf  import PyPDF, PyPDF2 
from   builtins             import range
from   ostap.utils.utils    import timing
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_pypdf' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## make simple test mass 
mass     = ROOT.RooRealVar ( 'test_mass' , 'Some test mass' , 3.0 , 3.2 )

## book very simple data set
varset0   = ROOT.RooArgSet  ( mass )
dataset   = ROOT.RooDataSet ( dsID() , 'Test Data set-0' , varset0 )  

mmin , mmax = mass.minmax()

## fill it 
m = VE(3.100,0.015**2)
for i in range(0,1000) :
    mass.value = m.gauss () 
    dataset.add ( varset0 )
for i in range(0,100) :
    mass.value = random.uniform ( *mass.minmax() )
    dataset.add ( varset0 )

logger.info ('DATASET %s' % dataset )

NORM = 1.0 / math.sqrt ( 2.0 * math.pi )

# =============================================================================
## @class PyGauss
#  local ``pure-python'' PDF 
## class PyGauss(MASS,PyPDF) :
class PyGauss(MASS,PyPDF) :
    """Local ``pure-python'' PDF """    
    def __init__ ( self         ,
                   name         ,
                   xvar         ,
                   mean  = ( 3.080 , 3.05  , 3.15  ) ,
                   sigma = ( 0.010 , 0.005 , 0.020 ) ,
                   pypdf = None ) :

        
        MASS .__init__ ( self , name      , xvar , mean , sigma ) 
        PyPDF.__init__ ( self , self.name , ( self.xvar  ,
                                              self.mean  ,
                                              self.sigma ) , pypdf = pypdf )

        self.pdf = self.pypdf
        
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'mean'  : self.mean  ,
            'sigma' : self.sigma ,
            'pypdf' : None         ## attention! 
            }
        
    ## the  main method 
    def evaluate ( self ) :
        
        vlist = self.varlist
        
        x = self.variable ( 0 ) 
        m = self.variable ( 1 )  
        s = self.variable ( 2 )
        
        dx = ( x - m ) / s        
        return math.exp ( -0.5 * dx * dx ) * NORM / s 

CDF  = Ostap.Math.gauss_cdf
# =============================================================================
## @class PyGaussAI
#  local ``pure-python'' PDF with analytical integrals 
class PyGaussAI(PyGauss) :
    """Local ``pure-python'' PDF with analytical integrals """    
    def __init__ ( self         ,
                   name         ,
                   xvar         ,
                   mean  = (         3.05  , 3.15  ) ,
                   sigma = ( 0.010 , 0.005 , 0.025 ) ,
                   pypdf = None ) :

        PyGauss.__init__ (  self , name , xvar , mean , sigma , pypdf = pypdf )

    ## declare analytical integral 
    def get_analytical_integral ( self ) :
        """Declare the analytical integral"""
        
        x  = self.varlist[0]
        
        if self.matchArgs ( x ) : return 1 ## get the integration code
        
        return 0
    
    ## calculate the analytical integral 
    def analytical_integral ( self ) :
        """Calculate the analytical integral"""

        assert 1 == self.intCode , 'Invalid integration code!'
        
        vlist = self.varlist

        rn    = self.rangeName        
        xv    = vlist [ 0 ]        
        xmax  = xv.getMax ( rn )
        xmin  = xv.getMin ( rn )
        
        m     = float ( vlist [ 1 ] ) 
        s     = float ( vlist [ 2 ] )        
        
        return CDF ( xmax , m , s  ) - CDF ( xmin , m , s  )

# =============================================================================
## @class PyGauss2
#  local ``pure-python'' PDF 
class PyGauss2(MASS,PyPDF2) :
    """Local ``pure-python'' PDF """    
    def __init__ ( self         ,
                   name         ,
                   function     , 
                   xvar         ,                   
                   mean  = ( 3.080 , 3.05  , 3.15  ) ,
                   sigma = ( 0.010 , 0.005 , 0.020 ) ,
                   title = '' 
                   ) :
        
        MASS  .__init__ ( self , name      , xvar , mean , sigma ) 
        PyPDF2.__init__ ( self      ,
                          name     = self.name ,
                          function = function  ,
                          vars     = ( self.xvar , self.mean , self.sigma ) )
        
        self.config = {
            'name'     : self.name     ,
            'function' : self.function , 
            'xvar'     : self.xvar     ,
            'mean'     : self.mean     ,
            'sigma'    : self.sigma    ,
            }
        
# =============================================================================
## pygauss PDF
# =============================================================================
def test_pygauss() :

    from ostap.core.meta_info import root_version_int 
    if root_version_int >= 62200 :
        logger.warning("test_pygauss: test is disabled for ROOT version %s" % root_version_int )
        return 

    logger.info ('Test PyGauss:  simple Gaussian signal' )
    
    gauss   = PyGauss( 'PyGauss'   , xvar = mass )
    
    ## model 
    model = Fit1D ( signal = gauss , background = None  , name = 'M0' )

    ## r1, f1  = gauss  .fitTo ( dataset , draw = True , silent = True , ncpu=1 )
    ## print (r1) 

    r2, f2  = model  .fitTo ( dataset , draw = True , silent = True , ncpu=1 )
    print (r2) 


# =============================================================================
## pygauss PDF
# =============================================================================
def test_pygauss_AI() :
    
    from ostap.core.meta_info import root_version_int 
    if root_version_int >= 62200 :
        logger.warning("test_pygauss_AI: test is disabled for ROOT version %s" % root_version_int )
        return 

    logger.info ('Test PyGaussAI:  simple Gaussian signal with  analytical integral' )
    
    gauss   = PyGauss( 'PyGaussAI' , xvar = mass )
    
    ## model 
    model   = Fit1D ( signal = gauss , background = None  , name = 'M1' )

    ## r1, f1  = gauss.fitTo ( dataset , draw = True , silent = True , ncpu=1 )
    ## print(r1)
    
    
    r2, f2  = model  .fitTo ( dataset , draw = True , silent = True , ncpu=1 )
    print (r2) 

def test_pygauss2 () :
    
    from ostap.core.meta_info import root_version_int 
    if root_version_int >= 62200 :
        logger.warning("test_pygauss2: test is disabled for ROOT version %s" % root_version_int )
        return 

    ## the function
    def function ( x , m , s ) :
        dx = ( x - m ) / s        
        return math.exp ( -0.5 * dx * dx ) * NORM / s 
    
    ## construct PDF
    gauss = PyGauss2 ( name  ='G2' , function = function , xvar = mass  )
    
    ## model 
    model = Fit1D ( signal = gauss , background = 'p0' , name = 'Q2' )

    ## r1, f1  = gauss  .fitTo ( dataset , draw = True , silent = True , ncpu=1 )
    ## print (r1) 

    r2, f2  = model  .fitTo ( dataset , draw = True , silent = True , ncpu=1 )
    print (r2) 

# =============================================================================
if '__main__' == __name__ :


    
    ## simple Gaussian PDF
    with timing ("PyGauss    ", logger ) : test_pygauss        ()
    
    ## simple Gaussian PDF with analytical integral 
    with timing ("PyGaussAI  ", logger ) : test_pygauss_AI  ()
    
    ## simple Gaussian PDF
    with timing ("PyGauss2   ", logger ) : test_pygauss2       ()
    
    
# =============================================================================
##                                                                      The END 
# ============================================================================= 
