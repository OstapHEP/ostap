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
__author__ = "Ostap developers"
__all__    = () ## nothing to import
# ============================================================================= 
from   ostap.core.core          import VE, dsID, Ostap
from   ostap.fitting.pdfbasic   import Generic1D_pdf 
from   ostap.fitting.fit1d      import PEAK ,  Fit1D 
from   ostap.utils.utils        import timing
from   ostap.core.meta_info     import old_PyROOT, root_info  
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait 
import ostap.fitting.roofit 
import ROOT, random, math, time 
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

logger.info ('DATASET\n%s' % dataset )

NORM  = 1.0 / math.sqrt ( 2.0 * math.pi )
CDF   = Ostap.Math.gauss_cdf

from ostap.fitting.models import Gauss_pdf
pdf   = Gauss_pdf ( "G"                               ,
                    xvar  = mass                      , 
                    mean  = ( 3.080 , 3.05  , 3.15  ) ,
                    sigma = ( 0.010 , 0.005 , 0.020 ) ) 

## build fit model 
model0   = Fit1D ( signal = pdf , background = None  )
model0.S = 900
model0.B = 100

S          = model0.S
B          = model0.B
background = model0.background_components[0]

# =============================================================================
## Test pure python PDF: <code>PyPDF2</code> 
#  @see Ostap::Models::PyPdf2 
def test_PyPDF2 () :
    """Test pure python PDF: PyPDF2
    - see Ostap.Models.PyPdf2 
    """

    logger = getLogger("test_PyPDF2")
    
    logger.info  ("Test pure python PDF: PyPDF2 with python function")
        
    from   ostap.fitting.pypdf  import PyPDF2
    # =============================================================================
    ## @class PyGauss2
    #  local ``pure-python'' PDF 
    class PyGauss2(PEAK,PyPDF2) :
        """Local ``pure-python'' PDF """    
        def __init__ ( self         ,
                       name         ,
                       function     , 
                       xvar         ,                   
                       mean  = ( 3.080 , 3.05  , 3.15  ) ,
                       sigma = ( 0.010 , 0.005 , 0.020 ) ,
                       title = '' 
                       ) :
            
            PEAK  .__init__ ( self , name      , xvar , mean , sigma ) 
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

    ## the function
    def function ( x , m , s ) :
        dx = ( x - m ) / s        
        return math.exp ( -0.5 * dx * dx ) * NORM / s 
    
    with timing ("Using-PyPDF2", logger ) :
        
        pdf.mean  = random.gauss ( 3.100 , 0.010 )
        pdf.sigma = random.gauss ( 0.012 , 0.001 )

        ## construct PDF
        gauss = PyGauss2 ( name  ='G2' , function = function ,
                           xvar = mass , mean = pdf.mean , sigma = pdf.sigma )
        
        ## model 
        model = Fit1D ( signal = gauss , background = background , S = S , B = B  , name = 'M2' )
        
        ##  fit!
        r, _ = model  .fitTo ( dataset , draw = False , silent = True , ncpu=1 )
        with wait ( 1 ) , use_canvas ( "test_PyPDF2" ) : 
            r, f = model  .fitTo ( dataset , draw = True  , silent = True , ncpu=1 )
            logger.info  ("Fit result for `pure python' PDF: PyPDF2 with python function \n%s" % r.table ( prefix = "# " ) )
            
        del model
        del gauss
        del PyGauss2

# =============================================================================
if '__main__' == __name__ :

    
    test_PyPDF2   ()

    del model0
    del pdf
    del background
    
    
# =============================================================================
##                                                                      The END 
# ============================================================================= 
