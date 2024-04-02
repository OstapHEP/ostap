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
import ostap.fitting.roofit 
from   builtins                 import range
from   ostap.core.meta_info     import root_info 
from   ostap.core.core          import VE, dsID, Ostap
from   ostap.fitting.pdfbasic   import Generic1D_pdf 
from   ostap.fitting.fit1d      import PEAK ,  Fit1D 
from   ostap.utils.utils        import timing
from   ostap.core.meta_info     import old_PyROOT 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait 
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
pdf   = Gauss_pdf( "G"                               ,
                   xvar  = mass                      , 
                   mean  = ( 3.080 , 3.05  , 3.15  ) ,
                   sigma = ( 0.010 , 0.005 , 0.020 ) ) 


my_exp = math.exp 
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
    #  local ``pure-python'' PDF 
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
            return my_exp ( -0.5 * dx * dx ) * NORM / s 

    with timing ("Using-PyPDF", logger ) :
        
        pdf.mean  = random.gauss ( 3.100 , 0.010 )
        pdf.sigma = random.gauss ( 0.012 , 0.001 )

        ## create PDF 
        gauss   = PyGauss( 'PyGauss'   , xvar = mass ,  
                           mean = pdf.mean , sigma = pdf.sigma )
        
        ## build fit model 
        model = Fit1D ( signal = gauss , background = None  , name = 'M0' )
        
        ##  fit!
        r, _ = model  .fitTo ( dataset , draw = False , silent = True , ncpu=1 )
        with wait ( 1 ) , use_canvas ( "test_PyPDF" ) : 
            r, f = model  .fitTo ( dataset , draw = True  , silent = True , ncpu=1 )            
            logger.info  ("Fit result for `pure python' PDF: PyPDF \n%s" % r.table ( prefix = "# " ) )
        

# =============================================================================
## Test pure python PDF: <code>PyPDF</code> with analytical integral
#  @attention For *OLD* PyROOT only!
#  @see Ostap::Models::PyPdf 
def test_PyPDF_AI() :
    """Test pure python PDF: PyPDF with analytical integral 
    - For *OLD* PyROOT only!
    - see Ostap.Models.PyPdf 
    """

    logger = getLogger("test_PyPDF_AI")
    
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

    with timing ("Using-PyPDF+AI", logger ) :

        pdf.mean  = random.gauss ( 3.100 , 0.010 )
        pdf.sigma = random.gauss ( 0.012 , 0.001 )

        ## create PDF 
        gauss   = PyGaussAI( 'PyGaussAI'   , xvar = mass ,
                             mean = pdf.mean , sigma = pdf.sigma )
        
        ## build fit model 
        model = Fit1D ( signal = gauss , background = None  , name = 'M1' )
        
        ##  fit!
        r, _ = model  .fitTo ( dataset , draw = False , silent = True , ncpu=1 )
        with wait ( 1 ) , use_canvas ( "test_PyPDF_AI" ) : 
            r, f = model  .fitTo ( dataset , draw = True  , silent = True , ncpu=1 )        
            logger.info  ("Fit result for `pure python' PDF: PyPDF with analytical integral \n%s" % r.table ( prefix = "# " ) ) 


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
        """Local ``pure-python'' PxDF """    
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
        model = Fit1D ( signal = gauss , background = 'p0' , name = 'M2' )
        
        ##  fit!
        r, _ = model  .fitTo ( dataset , draw = False , silent = True , ncpu=1 )
        with wait ( 1 ) , use_canvas ( "test_PyPDF2" ) : 
            r, f = model  .fitTo ( dataset , draw = True  , silent = True , ncpu=1 )
            logger.info  ("Fit result for `pure python' PDF: PyPDF2 with python function \n%s" % r.table ( prefix = "# " ) )

        del r, model, gauss   

    
# =============================================================================
## Test pure python PDF: <code>PyPdf</code>
#  @attention For *NEW* PyROOT only!
#  @see Ostap::Models::PyPdf 
def test_PyPdf() :
    """Test pure python PDF: pyPDF
    - For *NEW* PyROOT only!
    - see Ostap.Models.PyPdf 
    """

    logger = getLogger("test_PyPdf")
    
    logger.info  ("Test pure python PDF: PyPdf ")

    if old_PyROOT :
        logger.warning("test enabled only for NEW PyROOT! quit!")
        return
    
    if (6,31) <= root_info :
        logger.warning ( 'Test is TEMPORARILY disabled for ROOT>6.31/01 %s' % ROOT.gROOT.GetVersion() ) 
        return

    # =============================================================================
    ## @class MyGauss1
    #  local ``pure-python'' PDF 
    class MyGauss1(Ostap.Models.PyPdf) :
        """Local ``pure-python'' PDF
        """    
        def __init__ ( self , name , xvar = None , mean = None , sigma = None , clone = None  ) :


            if clone and isinstance ( clone, ROOT.RooAbsPdf  ) :

                super(MyGauss1,self).__init__ ( clone , name )

            else  :
                
                vars = ROOT.RooArgList()
                
                vars.add ( xvar  )
                vars.add ( mean  )
                vars.add ( sigma )
                
                super(MyGauss1,self).__init__ ( name , 'title' , vars )
            
        ## the  main method 
        def evaluate ( self ) :
            
            vlist = self.varlist
            
            x = self.variable ( 0 ) 
            m = self.variable ( 1 )  
            s = self.variable ( 2 )
            
            dx = ( x - m ) / s        
            return math.exp ( -0.5 * dx * dx ) * NORM / s 

        def clone ( self , newname ) :

            cl = MyGauss1 ( newname , clone = self  ) 
            ROOT.SetOwnership ( cl  , False )
            
            return cl

        
    with timing ("Using-PyPdf", logger ) :

        pdf.mean  = random.gauss ( 3.100 , 0.010 )
        pdf.sigma = random.gauss ( 0.012 , 0.001 )

        pdf_  = MyGauss1 ( 'MyGauss1' , pdf.xvar , pdf.mean , pdf.sigma )
        gauss = Generic1D_pdf  ( pdf_ , xvar = pdf.xvar )
                
        ## build fit model 
        model = Fit1D ( signal = gauss , background = None  , name = 'M3' )
        
        ##  fit!
        r, _ = model  .fitTo ( dataset , draw = False , silent = True , ncpu=1 )
        with wait ( 1 ) , use_canvas ( "test_PyPdf" ) : 
            r, f = model  .fitTo ( dataset , draw = True  , silent = True , ncpu=1 )        
            logger.info  ("Fit result for `pure python' PDF: PyPdf \n%s" % r.table ( prefix = "# " ) )

        del r,f
        del model, gauss,
        del pdf_    
    
# =============================================================================
## Test pure python PDF: <code>PyPdf</code> + analytical integrals 
#  @attention For *NEW* PyROOT only!
#  @see Ostap::Models::PyPdf 
def test_PyPdf_AI() :
    
    """Test pure python PDF: pyPDF
    - For *NEW* PyROOT only!
    - see Ostap.Models.PyPdf 
    """

    logger = getLogger("test_PyPdf_AI")
    logger.info ( "Test pure python PDF: PyPdf with analytical  integral")
    
    if old_PyROOT :
        logger.warning("test enabled only for NEW PyROOT! quit.")
        return

    if (6,31) <= root_info :
        logger.warning ( 'Test is TEMPORARILY disabled for ROOT>6.31/01 %s' % ROOT.gROOT.GetVersion() ) 
        return
    
    # =============================================================================
    ## @classMyGauss2
    #  local ``pure-python'' PDF 
    class MyGauss2(Ostap.Models.PyPdf) :
        """Local `pure-python' PDF
        """    
        def __init__ ( self , name , xvar = None ,  mean = None  , sigma = None , clone = None ) :

            print ( 'CONSTRUCTOR/1' , name ) 
            if clone and isinstance ( clone , ROOT.RooAbsPdf ) :
                
                print ( 'CONSTRUCTOR/C' , name ) 
                super(MyGauss2,self).__init__ ( clone, name )
                
            else  :
                    
                print ( 'CONSTRUCTOR/R', name  ) 
                vars = ROOT.RooArgList()
                
                vars.add ( xvar  )
                vars.add ( mean  )
                vars.add ( sigma )
                
                super(MyGauss2,self).__init__ ( name , 'title' , vars )
                
        ## the  main method 
        def evaluate ( self ) :
                    
            x = self.variable ( 0 ) 
            m = self.variable ( 1 )  
            s = self.variable ( 2 )
            
            dx = ( x - m ) / s        
            return math.exp ( -0.5 * dx * dx ) * NORM / s 
    
        def clone ( self , newname ) :

            if not newname : newname  = self.name  
            cl = MyGauss2 ( newname , clone = self )
            
            ROOT.SetOwnership ( cl  , False )

            for s in cl.servers():
                print ( 'server' , s )
                
            return cl
        
        ## declare analytical integral 
        def get_analytical_integral ( self ) :
            """Declare the analytical integral"""
            
            x  = self.variables()[0]
            
            if self.matchArgs ( x ) : return 1 ## get the integration code
            
            return 0
    
        ## calculate the analytical integral 
        def analytical_integral ( self ) :
            """Calculate the analytical integral"""
            
            assert 1 == self.intCode () , 'Invalid integration code!'
            
            vlist = self.variables() 
            
            rn    = self.rangeName()        
            xv    = vlist [ 0 ]        
            xmax  = xv.getMax ( rn )
            xmin  = xv.getMin ( rn )
            
            m     = float ( vlist [ 1 ] ) 
            s     = float ( vlist [ 2 ] )        
            
            return CDF ( xmax , m , s  ) - CDF ( xmin , m , s  )
        

    with timing ("Using-PyPdf+AI", logger ) :

        pdf.mean  = random.gauss ( 3.100 , 0.010 )
        pdf.sigma = random.gauss ( 0.012 , 0.001 )

        pdf_  = MyGauss2 ( 'MyGauss2' , pdf.xvar , pdf.mean , pdf.sigma )
        
        ROOT.SetOwnership ( pdf_ , True)

        gauss = Generic1D_pdf  ( pdf_ , xvar = pdf.xvar )
        
        ## build fit model 
        model = Fit1D ( signal = gauss , background = None  , name = 'M3' )
        
        ##  fit!
        r, _ = model  .fitTo ( dataset , draw = False , silent = True , ncpu=1 )
        with wait ( 1 ) , use_canvas ( "test_PyPdf_AI" ) : 
            r, f = model  .fitTo ( dataset , draw = True  , silent = True , ncpu=1 )            
            logger.info  ("Fit result for `pure python' PDF: PyPdf with analytical integral\n%s" % r.table ( prefix = "# " ) )

    del r, f,
    del model, gauss, pdf_
    
# =============================================================================
if '__main__' == __name__ :

    
    test_PyPDF    ()
    test_PyPDF_AI ()
    test_PyPDF2   ()
    test_PyPdf    ()
    test_PyPdf_AI ()


# =============================================================================
##                                                                      The END 
# ============================================================================= 
