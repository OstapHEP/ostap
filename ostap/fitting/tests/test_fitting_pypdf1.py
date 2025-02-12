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
from   ostap.utils.utils        import wait, batch_env  
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
## set batch from environment 
batch_env ( logger )
# =============================================================================
#
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
## Test pure python PDF: <code>PyPDF</code>
#  @attention For *OLD* PyROOT only!
#  @see Ostap::Models::PyPdf 
def test_PyPDF() :
    """Test pure python PDF: pyPDF
    - For *OLD* PyROOT only!
    - see Ostap.Models.PyPdf 
    """

    logger = getLogger("test_PyPDF")
    
    if not old_PyROOT :
        logger.warning("test enabled only for *(very)OLD* PyROOT!")
        return
    
    logger.info  ("Test pure python PDF: PyPDF ")

    from   ostap.fitting.pypdf  import PyPDF    
    # =============================================================================
    ## @class PyGauss
    #  local `pure-python' PDF 
    class PyGauss(PEAK,PyPDF) :
        """Local ``pure-python'' PDF
        """    
        def __init__ ( self         ,
                       name         ,
                       xvar         ,
                       mean  = ( 3.080 , 3.05  , 3.15  ) ,
                       sigma = ( 0.010 , 0.005 , 0.020 ) ,
                       pypdf = None ) :
            
            PEAK .__init__ ( self , name      , xvar , mean , sigma ) 
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

    with timing ("Using-PyPDF", logger ) :
        
        pdf.mean  = random.gauss ( 3.100 , 0.010 )
        pdf.sigma = random.gauss ( 0.012 , 0.001 )

        ## create PDF 
        gauss   = PyGauss( 'PyGauss'   , xvar = mass ,  
                           mean = pdf.mean , sigma = pdf.sigma )
        
        ## build fit model 
        model = Fit1D ( signal = gauss , background = background , S = S , B = B  , name = 'M0' )
        
        ##  fit!
        r, _ = model  .fitTo ( dataset , draw = False , silent = True , ncpu=1 )
        with wait ( 1 ) , use_canvas ( "test_PyPDF" ) : 
            r, f = model  .fitTo ( dataset , draw = True  , silent = True , ncpu=1 )            
            logger.info  ("Fit result for `pure python' PDF: PyPDF \n%s" % r.table ( prefix = "# " ) )
        
        del model
        del gauss
            
# =============================================================================
## Test pure python PDF: <code>PyPDF</code> with analytical integral
#  @attention For *OLD* PyROOT only!
#  @see Ostap::Models::PyPdf 
def test_PyPDF_AI() :
    """Test pure python PDF: PyPDF with analytical integral 
    - For *OLD* PyROOT only!
    - see Ostap.Models.PyPdf 
    """

    logger = getLogger("test_PyPDF&AI")
    
    if not old_PyROOT :
        logger.warning("test enabled only for *(very)OLD* PyROOT!")
        return
    
    logger.info  ("Test pure python PDF: PyPDF with analytical integral")


    from   ostap.fitting.pypdf  import PyPDF    
    # =========================================================================
    ## @class PyGaussAI
    #  local ``pure-python'' PDF with analytical integrals 
    class PyGaussAI(PEAK,PyPDF) :
        """Local ``pure-python'' PDF with analytical integrals """    
        def __init__ ( self         ,
                       name         ,
                       xvar         ,
                       mean  = ( 3.080 , 3.05  , 3.15  ) ,
                       sigma = ( 0.010 , 0.005 , 0.020 ) ,
                       pypdf = None ) :
            
            
            PEAK .__init__ ( self , name      , xvar , mean , sigma ) 
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

    with timing ("Using-PyPDF&AI", logger ) :

        pdf.mean  = random.gauss ( 3.100 , 0.010 )
        pdf.sigma = random.gauss ( 0.012 , 0.001 )

        ## create PDF 
        gauss   = PyGaussAI( 'PyGaussAI'   , xvar = mass ,
                             mean = pdf.mean , sigma = pdf.sigma )
        
        ## build fit model 
        model = Fit1D ( signal = gauss , background = background , S = S , B = B  , name = 'M1' )
        
        ##  fit!
        r, _ = model  .fitTo ( dataset , draw = False , silent = True , ncpu=1 )
        with wait ( 1 ) , use_canvas ( "test_PyPDF&AI" ) : 
            r, f = model  .fitTo ( dataset , draw = True  , silent = True , ncpu=1 )        
            logger.info  ("Fit result for `pure python' PDF: PyPDF with analytical integral \n%s" % r.table ( prefix = "# " ) ) 

        del model
        del gauss
            

# =============================================================================
if '__main__' == __name__ :

    
    test_PyPDF    () 
    test_PyPDF_AI ()
    
    del model0
    del pdf
    del background
    
# =============================================================================
##                                                                      The END 
# ============================================================================= 
