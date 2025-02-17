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
from   ostap.utils.timing       import timing
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait, batch_env
from   ostap.utils.basic        import typename 
from   ostap.fitting.models     import Gauss_pdf, Generic1D_pdf,  Fit1D  
from   ostap.fitting.pypdf      import PyPDF, PyPDFLite
from   ostap.io.zipshelve       import tmpdb 
import ostap.fitting.roofit 
import ROOT, math 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_fitting_pypdf1' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
## set batch from environment 
batch_env ( logger )
# =============================================================================

## make simple test mass 
xvar       = ROOT.RooRealVar ( 'test_mass' , 'Some test mass' , 3.0 , 3.2 )

## referece Gaussian
gauss_ref = Gauss_pdf ( "G"                               ,
                        xvar  = xvar                      , 
                        mean  = ( 3.100 , 3.050 , 3.150 ) ,
                        sigma = ( 0.015 , 0.005 , 0.025 ) )
## reference model
model_ref   = Fit1D ( signal = gauss_ref , background = "flat" , suffix =  "R" )
model_ref.S = 900
model_ref.B = 100

## generate dataset 
dataset = model_ref.generate ( 1000 )
NORM    = 1.0/math.sqrt ( 2 * math.pi )

logger.info ('DATASET\n%s' % dataset )

CDF   = Ostap.Math.gauss_cdf

def mg1_factory ( *args ) : return MyGauss1 ( *args )
def mg2_factory ( *args ) : return MyGauss2 ( *args )
    
# ==============================================================================
## @class MyGauss1
#  The simplest pythonic PDF 
#  - see <code>clone</code>
#  - see <code>evaluate</code>
class MyGauss1(PyPDF) :
    """ The simplest pythonic PDF 
    - see `clone`
    - see `evaluate`
    """
    
    def __init__ ( self , name , title , xvar , mean , sigma , clone = None ) :

        assert not clone or isinstance ( clone , MyGauss1 ) , \
            "MyGauss1: nivalid `clone` type %s" % typename ( clone )

        ## order of variables 
        self._variables = ( xvar , mean , sigma )
        
        if clone : super(MyGauss1,self).__init__ ( clone = clone , name = name )
        else     : super(MyGauss1,self).__init__ ( name  = name  , title = title , variables = self._variables )

    ## mandatory clone method 
    def clone ( self , newname = '' ) :
        cloned = MyGauss1 ( newname if newname else self.name , self.title , None , None , None , clone = self )
        ROOT.SetOwnership ( cloned  , False )
        return cloned

    ## mandatory evaluate method 
    def evaluate ( self ) :
        x, mean, sigma = self.variables
        return Ostap.Math.gauss_pdf ( float ( x ) , float ( mean ) , float ( sigma ) )

    def __reduce__  ( self ) :
        return mg1_factory , ( self.name , self.title ) + tuple ( v for v in self.variables ) 
# ==============================================================================
## @class MyGauss3
#  A bit mor ecomplicated methpod using analystical integrals 
#  - see <code>clone</code>
#  - see <code>evaluate</code>
class MyGauss2(PyPDF) :
    """ The simplest pythonic PDF 
    - see `clone`
    - see `evaluate`
    """
    
    def __init__ ( self , name , title , xvar , mean , sigma , clone = None ) :
        
        assert not clone or isinstance ( clone , MyGauss2 ) , \
            "MyGauss2: nivalid `clone` type %s" % typename ( clone )

        ## order of variables 
        self._variables = ( xvar , mean , sigma )
        
        if clone : super(MyGauss2,self).__init__ ( clone = clone , name = name )
        else     : super(MyGauss2,self).__init__ ( name  = name  , title = title , variables = self._variables )

    ## mandatory clone method 
    def clone ( self , newname = '' ) :
        newname = newname if newname else self.name 
        cloned  = MyGauss2 ( newname , self.title , None , None , None , clone = self )
        ROOT.SetOwnership ( cloned  , False )
        return cloned

    ## mandatory evaluate method 
    def evaluate ( self ) :
        x, mean, sigma = self.variables 
        return Ostap.Math.gauss_pdf ( float ( x ) , float ( mean ) , float ( sigma ) )
    
    ## Declare the analystical integrtaio over x 
    def get_analytical_integral ( self ) :
        """ Declare the analystical integrtaio over x"""
        x  = self.variables [ 0 ] 
        return 1 if self.matchArg ( x ) else 0 

    ## Perform analytical integration over x 
    def analytical_integral     ( self ) :
        """ Perform analytical integration over x"""
        
        assert 1 == self.intCode () , 'Invalid integration code!'

        xvar, mean, sigma = self.variables 

        rn    = self.rangeName() 
        xmin, xmax = xvar.getMin ( rn ) , xvar.getMax ( rn )
        
        m = float ( mean  ) 
        s = float ( sigma )        
        
        return CDF ( xmax , m , s  ) - CDF ( xmin , m , s  )

    def __reduce__  ( self ) :
        return mg2_factory , ( self.name , self.title ) + tuple ( v for v in self.variables ) 

# =============================================================================
## Test pure python PDF: <code>PyPDF</code>
#  @see ostap.fitting.pypdf.PyPDF
#  @see Ostap::Models::PyPdf 
def test_PyPDF1() :
    """Test pure python PDF: pyPDF
    - see `ostap.fitting.pypdf.PyPDF`
    - see `Ostap::Models::PyPdf` 
    """

    logger = getLogger("test_PyPDF1")
    logger.info ( 'Use plain PyPDF' ) 
    
    mygauss = MyGauss1 ( 'MyGauss1'                           ,
                          title    = "local pure python PDF" ,
                          xvar     = gauss_ref.xvar          ,
                          mean     = gauss_ref.mean          ,
                          sigma    = gauss_ref.sigma         ) 
    
    signal = Generic1D_pdf ( pdf = mygauss , xvar = xvar )
    model  = Fit1D ( signal = signal  , background = None , suffix = '_P1' )

    ##  fit!
    with use_canvas ( "test_PyPDF1: plain PyPDF" , wait = 2 ) : 
        r, _ = model.fitTo ( dataset , draw = True , nbins = 50 , quiet  = True  )

    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        with tmpdb ( )as db :
            db [ mygauss.name ] = mygauss
            db ['signal'      ] = signal 
            db ['gauss_ref'   ] = gauss_ref
            db ['model_ref'   ] = model_ref        
            db.ls()
        # =====================================================================
    except : # ================================================================
        # =====================================================================
        logger.error ( "Cannot serialize!" , exc_info = True )

# =============================================================================
## Test pure python PDF: <code>PyPDF</code> with analytical integration 
#  @see ostap.fitting.pypdf.PyPdf
#  @see Ostap::Models::PyPdf 
def test_PyPDF2() :
    """Test pure python PDF: pyPDF wth analysticla integration 
    - see `ostap.fitting.pypdf.PyPDF
    - see `Ostap::Models::PyPdf` 
    """

    logger = getLogger("test_PyPDF2")
    logger.info ('Use PyPDF with analytical integrals' )
    
    mygauss = MyGauss2 ( 'MyGauss2'                           ,
                         title    = "local pure python PDF" ,
                         xvar     = gauss_ref.xvar          ,
                         mean     = gauss_ref.mean          ,
                         sigma    = gauss_ref.sigma         ) 
    
    signal = Generic1D_pdf ( pdf = mygauss , xvar = xvar )
    model  = Fit1D ( signal = signal  , background = None , suffix = '_P2' )
    
    ##  fit!
    with use_canvas ( "test_PyPDF2: use analytical integrals" , wait = 2 ) : 
        r, _ = model.fitTo ( dataset , draw = True , nbins = 50 , quiet  = True  )

    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        with tmpdb ( )as db :
            db [ mygauss.name ] = mygauss
            db ['signal'      ] = signal 
            db ['gauss_ref'   ] = gauss_ref
            db ['model_ref'   ] = model_ref        
            db.ls()
        # =====================================================================
    except : # ================================================================
        # =====================================================================
        logger.error ( "Cannot serialize!" , exc_info = True )

# =============================================================================
gauss_cpp = Ostap.Math.gauss_pdf

def gauss_fun ( x , mean , sigma ) :
    return gauss_cpp ( x , mean , sigma ) 

class GaussPdf ( object) :
    def __call__ ( self , x , mean , sigma ) :
        return gauss_fun ( x , mean , sigma )

gauss_obj = GaussPdf()

# =============================================================================
## Test pure python PDF: <code>PyPDFLite</code
#  @see ostap.fitting.pypdf.PyPDFLifgt 
#  @see Ostap::Models::PyPdfLight 
def test_PyPDFLite1() :
    """ Test pure python PDF: pyPDF wth analysticla integration 
    - see `ostap.fitting.pypdf.PyPdf`
    - see `Ostap::Models::PyPdf` 
    """

    logger   = getLogger("test_PyPDFLite1")
    logger.info ( "Use `light' PDF:  global funtion" )
        
    function  = gauss_fun 

    variables = gauss_ref.xvar, gauss_ref.mean, gauss_ref.sigma    
    mygauss   = PyPDFLite ( "MyGauss3"             ,
                             function   = function  ,
                             variables  = variables )
    
    logger.info ( 'Use %s' % mygauss ) 
    
    signal = Generic1D_pdf ( pdf = mygauss , xvar = xvar )
    model  = Fit1D ( signal = signal  , background = None , suffix = '_P2' )
    
    ##  fit!
    with use_canvas ( "test_PyPDFLite1: global function " , wait = 2 ) : 
        r, _ = model.fitTo ( dataset , draw = True , nbins = 50 , quiet  = True  )

    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        with tmpdb ( )as db :
            db [ mygauss.name ] = mygauss
            db ['signal'      ] = signal 
            db ['gauss_ref'   ] = gauss_ref
            db ['model_ref'   ] = model_ref        
            db.ls()
        # =====================================================================
    except : # ================================================================
        # =====================================================================
        logger.error ( "Cannot serialize!" , exc_info = True )


# =============================================================================
## Test pure python PDF: <code>PyPDFLite</code
#  @see ostap.fitting.pypdf.PyPDFLifgt 
#  @see Ostap::Models::PyPdfLight 
def test_PyPDFLite2() :
    """ Test pure python PDF: pyPDF wth analysticla integration 
    - see `ostap.fitting.pypdf.PyPdf`
    - see `Ostap::Models::PyPdf` 
    """

    logger   = getLogger("test_PyPDFLite2")
    logger.info ( "Use `light' PDF:  global c++ function" )
        
    function  = gauss_cpp 

    variables = gauss_ref.xvar, gauss_ref.mean, gauss_ref.sigma    
    mygauss   = PyPDFLite ( "MyGauss4"             ,
                             function   = function  ,
                             variables  = variables )

    logger.info ( 'Use %s' % mygauss ) 
    
    signal = Generic1D_pdf ( pdf = mygauss , xvar = xvar )
    model  = Fit1D ( signal = signal  , background = None , suffix = '_P2' )
    
    ##  fit!
    with use_canvas ( "test_PyPDFLite1: global function " , wait = 2 ) : 
        r, _ = model.fitTo ( dataset , draw = True , nbins = 50 , quiet  = True  )

    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        with tmpdb ( )as db :
            db [ mygauss.name ] = mygauss
            db ['signal'      ] = signal 
            db ['gauss_ref'   ] = gauss_ref
            db ['model_ref'   ] = model_ref        
            db.ls()
        # =====================================================================
    except : # ================================================================
        # =====================================================================
        logger.error ( "Cannot serialize!" , exc_info = True )
            
# =============================================================================
## Test pure python PDF: <code>PyPDFLite</code
#  @see ostap.fitting.pypdf.PyPDFLifgt 
#  @see Ostap::Models::PyPdfLight 
def test_PyPDFLite3() :
    """ Test pure python PDF: pyPDF wth analysticla integration 
    - see `ostap.fitting.pypdf.PyPdf`
    - see `Ostap::Models::PyPdf` 
    """

    logger   = getLogger("test_PyPDFLite2")
    logger.info ( "Use `light' PDF:  global object" )
        
    function  = gauss_obj 

    variables = gauss_ref.xvar, gauss_ref.mean, gauss_ref.sigma    
    mygauss   = PyPDFLite ( "MyGauss5"             ,
                             function   = function  ,
                             variables  = variables )

    logger.info ( 'Use %s' % mygauss ) 
    
    signal = Generic1D_pdf ( pdf = mygauss , xvar = xvar )
    model  = Fit1D ( signal = signal  , background = None , suffix = '_P2' )
    
    ##  fit!
    with use_canvas ( "test_PyPDFLite1: global object" , wait = 2 ) : 
        r, _ = model.fitTo ( dataset , draw = True , nbins = 50 , quiet  = True  )

    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        with tmpdb ( )as db :
            db [ mygauss.name ] = mygauss
            db ['signal'      ] = signal 
            db ['gauss_ref'   ] = gauss_ref
            db ['model_ref'   ] = model_ref        
            db.ls()
        # =====================================================================
    except : # ================================================================
        # =====================================================================
        logger.error ( "Cannot serialize!" , exc_info = True )


# =============================================================================
## Test pure python PDF: <code>PyPDFLite</code
#  @see ostap.fitting.pypdf.PyPDFLifgt 
#  @see Ostap::Models::PyPdfLight 
def test_PyPDFLite4() :
    """ Test pure python PDF: pyPDF wth analysticla integration 
    - see `ostap.fitting.pypdf.PyPdf`
    - see `Ostap::Models::PyPdf` 
    """

    logger   = getLogger("test_PyPDFLite2")
    logger.info ( "Use `light' PDF: local function" )

    def local_fun ( *args ) : return gauss_cpp ( *args )
    
    function  = local_fun 

    variables = gauss_ref.xvar, gauss_ref.mean, gauss_ref.sigma    
    mygauss   = PyPDFLite ( "MyGauss6"             ,
                             function   = function  ,
                             variables  = variables )

    logger.info ( 'Use %s' % mygauss ) 
    
    signal = Generic1D_pdf ( pdf = mygauss , xvar = xvar )
    model  = Fit1D ( signal = signal  , background = None , suffix = '_P2' )
    
    ##  fit!
    with use_canvas ( "test_PyPDFLite1: local function" , wait = 2 ) : 
        r, _ = model.fitTo ( dataset , draw = True , nbins = 50 , quiet  = True  )

    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        with tmpdb ( )as db :
            db [ mygauss.name ] = mygauss
            db ['signal'      ] = signal 
            db ['gauss_ref'   ] = gauss_ref
            db ['model_ref'   ] = model_ref        
            db.ls()
        # =====================================================================
    except : # ================================================================
        # =====================================================================
        logger.error ( "Cannot serialize!" , exc_info = True )
        
# =============================================================================
if '__main__' == __name__ :
    
    ## test_PyPDF1     () 
    ## test_PyPDF2     () 

    test_PyPDFLite1 () 
    test_PyPDFLite2 () 
    test_PyPDFLite3 () 
    test_PyPDFLite4 () 
    
# =============================================================================
##                                                                      The END 
# ============================================================================= 
