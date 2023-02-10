#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/fit3d.py
#  Set of useful basic utilities to build various fit models 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
# =============================================================================
"""Set of useful basic utilities to build various 2D-fit models"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'Flat3D'        , ## the most trivial 3D-pdf - constant
    'Model3D'       , ## trivial class to build 3D model from 1D-components 
    'Sum3D'         , ## non-extended sum two PDFs 
    'H3D_pdf'       , ## convertor of 1D-histo to RooDataPdf
    'Shape3D_pdf'   , ## simple PDF from C++ shape
    ##
    'Fit3D'         , ## the model for                3D-fit
    'Fit3DSym'      , ## the model for      symmetric 3D-fit
    'Fit3DMix'      , ## the model for half-symmetric 3D-fit
    ##
    )
# =============================================================================
from   builtins                 import range
from   ostap.core.core          import Ostap , valid_pointer, roo_silent 
from   ostap.core.ostap_types   import integer_types
from   ostap.fitting.utils      import component_similar , component_clone
from   ostap.fitting.fit2d      import Model2D 
from   ostap.fitting.pdfbasic   import PDF2,PDF3, Generic3D_pdf
from   ostap.fitting.fithelpers import H3D_dset, Fractions  
import ROOT, random
# =============================================================================
from   ostap.logger.logger  import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.fit3d' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
## @class Flat3D
#  The most trivial 3D-model - constant
#  @code 
#  pdf = Flat3D( 'flat' , xvar = ...  , yvar = ... , zvar = ... )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
class Flat3D(PDF3) :
    """The most trival 3D-model - constant
    >>> pdf = Flat3D( 'flat' , xvar = ...  , yvar = ... , zvar = ... )
    """
    def __init__ ( self , xvar , yvar , zvar , name = ''  , title = '' ) :

        name = name if name else self.generate_name ( prefix = 'Flat3D_')                            
        PDF3.__init__ ( self  , name , xvar , yvar , zvar ) 
        
        if not title : title = 'flat3(%s)' % name 
        self.pdf = Ostap.Models.Uniform ( name , title , self.xvar , self.yvar , self.zvar )
        assert 3 == self.pdf.dim() , 'Flat3D: wrong dimensionality!'
        
        ## save configuration
        self.config = {
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'zvar'     : self.zvar ,
            'name'     : self.name ,            
            'title'    : title     ,             
            }

# ===========================================================================
## @class Model3D
#  Trivial class to construct 3D model as a product of split 1D-models
#  actually it is a tiny  wrapper over <code>ROOT.RooProdPdf</code>
#  @code
#  pdfx = ...
#  pdfy = ...
#  pdfz = ...
#  pdf3D = Model3D( 'D3' , xmodel = pdfx , ymodel =  pdfy , zmodel = pdfz )
#  @endcode 
class Model3D(PDF3) :
    """Trivial class to construct 3D model as a product of split 1D-models
    - actually it is a tiny  wrapper over ROOT.RooProdPdf
    >>> pdfx = ...
    >>> pdfy = ...
    >>> pdfz = ...
    >>> pdf3D = Model3D( 'D3' , xmodel = pdfx , ymodel =  pdfy , zmodel = pdfz )
    """
    def __init__ ( self         ,
                   name         ,
                   xmodel       ,
                   ymodel       , 
                   zmodel       ,
                   xvar  = None ,
                   yvar  = None ,
                   zvar  = None ,
                   title = ''   ) :

        if xvar and not xmodel : xmodel = Flat1D ( xvar )
        if yvar and not ymodel : ymodel = Flat1D ( yvar )
        if zvar and not zmodel : zmodel = Flat1D ( zvar )
                
        self.__xmodel , xvar = self.make_PDF1 ( xmodel , xvar = xvar , prefix = 'X' )
        self.__ymodel , yvar = self.make_PDF1 ( ymodel , xvar = yvar , prefix = 'Y' )
        self.__zmodel , zvar = self.make_PDF1 ( zmodel , xvar = zvar , prefix = 'Z' )
        
        name  = name  if name  else self.generate_name ( '(%s)*(%s)*(%s)'  % ( self.xmodel.name ,
                                                                               self.ymodel.name ,
                                                                               self.zmodel.name ) )
        ## initialize the base 
        PDF3.__init__ (  self                    ,
                         name = name             ,
                         xvar = self.xmodel.xvar ,
                         yvar = self.ymodel.xvar ,
                         zvar = self.zmodel.xvar ) 
        
        ## make pdf
        from ostap.fitting.pdf_ops import raw_product 
        self.pdf = raw_product ( self , self.xmodel , self.ymodel , self.zmodel )

        ## save configuration 
        self.config = {
            'name'   : self.name   ,
            'xmodel' : self.xmodel ,
            'ymodel' : self.ymodel ,
            'zmodel' : self.zmodel ,
            'xvar'   : self.xvar   ,
            'yvar'   : self.yvar   ,            
            'zvar'   : self.zvar   ,
            }

    ## redefine the clone 
    def clone ( self , **kwargs ) :
        """ Redefine the clone
        """
        
        name   = kwargs.pop ( 'name' , self.name )        
        xvar   = kwargs.pop ( 'xvar' , self.xvar )
        yvar   = kwargs.pop ( 'yvar' , self.yvar )        
        zvar   = kwargs.pop ( 'zvar' , self.zvar )
        
        xmodel = self.xmodel if self.xvar is xvar else self.xmodel.clone ( xvar = xvar , **kwargs )
        ymodel = self.ymodel if self.yvar is yvar else self.ymodel.clone ( xvar = yvar , **kwargs )
        zmodel = self.zmodel if self.zvar is zvar else self.zmodel.clone ( xvar = zvar , **kwargs )
        
        return PDF3.clone ( self ,
                            name   = name   , 
                            xmodel = xmodel ,
                            ymodel = ymodel ,
                            zmodel = zmodel ,
                            xvar   = xvar   ,
                            yvar   = yvar   , 
                            zvar   = zvar   , **kwargs )
  

    @property
    def xmodel ( self ) :
        """'xmodel' : x-component of M(x)*M(y)*M(z) PDF"""
        return self.__xmodel

    @property
    def ymodel ( self ) :
        """'y-model' : y-component of M(x)*M(y)*M(z) PDF"""
        return self.__ymodel
    
    @property
    def zmodel ( self ) :
        """'zmodel' : z-component of M(x)*M(y)*M(z) PDF"""
        return self.__zmodel




# =============================================================================
## @class Sum3D
#  Non-extended sum of several PDFs
#  It is just a small wrapper for <code>ROOT.RooAddPdf</code>
#  @see RooAddPdf 
class Sum3D (PDF3,Fractions) :
    """Non-extended sum of several PDFs:
    
    It is just a small wrapper for <code>ROOT.RooAddPdf</code>
    - see RooAddPdf 
    
    >>> sum  = Sum3D ( [ pdf1 , pdf2 , pdf3 ]  ) 
    
    """
    def __init__ ( self             ,
                   pdfs             , ## input list of PDFs  
                   xvar      = None , 
                   yvar      = None , 
                   zvar      = None , 
                   name      = ''   ,
                   recursive = True ,
                   prefix    = 'f'  , ## prefix for fraction names 
                   suffix    = ''   , ## suffix for fraction names 
                   fractions = None ) :

        assert 2 <= len ( pdfs ) , 'Sum3D: at least two PDFs are needed!'

        pdf_list = []           
        for i , p in enumerate ( pdfs ) :
            cmp , xvar , yvar , zvar = self.make_PDF3 ( p , xvar = xvar , yvar = yvar , zvar = zvar ) 
            pdf_list.append ( cmp )
            
        ## generic name 
        patname =  '+'.join ( '(%s)' % p.name for p in pdf_list )
        ## check the instance name 
        name    = name if name else self.new_name ( patname ) 

        ## initialize the base class
        PDF3.     __init__ ( self , name , xvar , yvar , zvar ) 
        Fractions.__init__ ( self , pdf_list ,
                             prefix    = prefix    ,
                             suffix    = suffix    ,
                             recursive = recursive ,
                             fractions = fractions ) 

        for p in self.pdfs      : self.alist1.add ( p.pdf )
        for f in self.frac_list : self.alist2.add ( f     )
        
        ## finally build PDF
        self.pdf = ROOT.RooAddPdf ( self.new_roo_name ( patname , suffix ) , 
                                    patname        ,
                                    self.alist1    ,
                                    self.alist2    ,
                                    self.recursive )
        
        self.config = {
            'pdfs'      : self.pdfs      ,
            'xvar'      : self.xvar      ,
            'yvar'      : self.yvar      ,
            'zvar'      : self.zvar      ,
            'name'      : self.name      , 
            'prefix'    : self.prefix    ,
            'suffix'    : self.suffix    ,
            'fractions' : self.fractions ,
            'recursive' : self.recursive        
            }
            
# =============================================================================
## Generic 2D-shape from C++ callable
#  @see Ostap::Models:Shape3D
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2020-07-20
class Shape3D_pdf(PDF3) :
    """ Generic 3D-shape from C++ callable
    - see Ostap::Models:Shape3D
    """
    
    def __init__ ( self , name , shape , xvar , yvar , zvar , tag = 0 ) :

        if isinstance ( shape , ROOT.TH3 ) and not xvar :
            xvar = shape.xminmax()

        if isinstance ( shape , ROOT.TH3 ) and not yvar :
            yvar = shape.yminmax()

        if isinstance ( shape , ROOT.TH3 ) and not zvar :
            zvar = shape.zminmax()
            
        if isinstance ( shape , ROOT.TH3 ) :
            
            self.histo = shape
            shape      = Ostap.Math.Histo3D     ( shape )
            tag        = Ostap.Utils.hash_histo ( shape ) 
            
        elif hasattr ( shape , 'tag' ) and not tag : 
            tag = shape.tag() 

        ##  iniialize the base 
        PDF3.__init__ ( self , name , xvar , yvar , zvar ) 
            
        self.__shape = shape
        self.__tag   = tag

        
        if isinstance ( self.shape , Ostap.Math.Histo2D ) :
            
            ## create the actual pdf
            self.pdf = Ostap.Models.Histo3D ( self.roo_name ( 'histo3_' ) , 
                                              "Histo-3D %s" % self.name   ,
                                              self.xvar                   ,
                                              self.yvar                   ,
                                              self.zvar                   ,
                                              self.shape                  )
        else :
            
            ## create the actual pdf
            self.pdf = Ostap.Models.Shape3D.create  (
                self.roo_name ( 'shape3_' ) , 
                "Shape-3D %s" % self.name   ,
                self.xvar                   ,
                self.yvar                   ,
                self.zvar                   ,
                self.shape                  ,
                self.tag                    ) 

        ## save the configuration
        self.config = {
            'name'    : self.name    , 
            'shape'   : self.shape   , 
            'xvar'    : self.xvar    , 
            'yvar'    : self.yvar    , 
            'zvar'    : self.zvar    , 
            'tag'     : self.tag     , 
            }
        
    @property
    def shape  ( self ) :
        """'shape' : the actual C++ callable shape"""
        return self.__shape  
    @property
    def tag   ( self ) :
        """'tag' : uqnue tag used for cache-integration"""
        return self.__tag 
  
# =============================================================================
## simple convertor of 3D-histogram into PDF
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H3D_pdf(PDF3) :
    """Simple convertor of 3D-histogram into PDF 
    """
    def __init__ ( self            ,
                   name            ,
                   histo           ,
                   xvar    = None  , 
                   yvar    = None  ,
                   zvar    = None  ,
                   density = False ,
                   order   = 0     , 
                   silent  = False ) :
        
        assert isinstance ( order, integer_types ) and 0 <= order ,\
               'Invalid interpolation order: %s/%s' % ( order , type ( order ) )

        self.__ds = H3D_dset ( histo , xvar = xvar , yvar = yvar , zvar = zvar ,
                               density = density ,  silent = silent )
        
        PDF3    .__init__ ( self , name = name  ,
                            xvar = self.ds.xvar ,
                            yvar = self.ds.yvar ,
                            zvar = self.ds.zvar ) 
        
        self.__vset  = ROOT.RooArgSet  ( self.xvar , self.yvar , self.zvar )
        
        #
        ## finally create PDF :
        #
        with roo_silent ( silent ) : 
            self.pdf    = ROOT.RooHistPdf (
                self.roo_name ( 'histo3_' ) , 
                'Histo-3D PDF: %s/%s' % ( histo3.GetName() , histo2.GetTitle() ) , 
                self.__vset  , 
                self.dset    ,
                order        )

        ## and declare it be be a "signal"
        self.signals.add ( self.pdf ) 
            
        ## save the configuration
        self.config = {
            'name'    : self.name    , 
            'histo'   : self.histo   , 
            'xvar'    : self.xvar    , 
            'yvar'    : self.yvar    , 
            'zvar'    : self.zvar    , 
            'density' : self.density , 
            'silent'  : self.silent  ,             
            'order'   : self.order   ,             
            }

    @property
    def ds ( self ) :
        """'ds' : the H3D_dset object"""
        return self.__ds 
    @property     
    def xaxis  ( self ) :
        """The histogram x-axis variable (same as xvar)"""
        return self.xvar
    @property     
    def yaxis  ( self ) :
        """The histogram y-axis variable (same as yvar)"""
        return self.yvar
    @property     
    def zaxis  ( self ) :
        """The histogram z-axis variable (same as zvar)"""
        return self.zvar
    
    @property
    def histo ( self ) :
        """The  histogram itself"""
        return self.ds.histo    
    @property
    def density( self ) :
        """Treat the histo as 'density' histogram?"""
        return self.ds.density    
    @property
    def skip_zero ( self ) :
        """'skip_zero' : skip zero bins for weighted dataset in histo?"""
        return self.ds.skip_zero    
    @property
    def silent( self ) :
        """Use the silent mode?"""
        return self.ds.silent
    @property
    def dset ( self ) :
        """'dset' : ROOT.RooDataHist object"""
        return self.ds.dset
    @property
    def histo_hash ( self ) :
        """Hash value for the histogram"""
        return self.ds.histo_hash
    @property
    def weight ( self ) :
        """'weight' : get weight variable if defined, None otherwise"""
        return self.ds.wvar
        
    @property
    def order  ( self ) :
        """'order' : interpolation order"""
        return self.pdf.getInterpolationOrder () 
    @order.setter
    def order  ( self , value ) :
        assert isinstance ( value , integer_types ) and 0 <= value,\
               'Invalid interpolation order %s/%s' % ( value , type ( value ) )
        self.pdf.setInterpolationOrder ( value )
      
# =============================================================================
# Compound models for 3D-fit
# =============================================================================

# =============================================================================
## @class Fit3D
#  The actual model for 3D-fits
#
#  @code
# 
#  model   = Models.Fit3D (
#      signalx = Models.Gauss_pdf ( 'Gx' , mass = m_x ) ,
#      signaly = Models.Gauss_pdf ( 'Gy' , mass = m_y ) ,
#      signalz = Models.Gauss_pdf ( 'Gz' , mass = m_z ) )
#
#  r = model.fitTo ( dataset ) ## fit dataset 
#
#  print ( r )                   ## get results  
#
#  fx  = model.draw1 ()          ## visualize X-projection
#  fy  = model.draw2 ()          ## visualize Y-projection
#  fz  = model.draw3 ()          ## visualize Z-projection
#
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2017-07-25
class Fit3D (PDF3) :
    """The actual model for 3D-fits
    
    >>>  model   = Models.Fit3D (
    ...      signal_x = Models.Gauss_pdf ( 'Gx' , mass = m_x ) ,
    ...      signal_y = Models.Gauss_pdf ( 'Gy' , mass = m_y ) ,
    ...      signal_x = Models.Gauss_pdf ( 'Gz' , mass = m_z ) ,
    ...      bkg_1x   = 1 , 
    ...      bkg_1y   = 0 ,
    ...      bkg_1z   = 0 )
    >>> r,f = model.fitTo ( dataset ) ## fit dataset 
    >>> print ( r  )                  ## get results  
    >>> fx  = model.draw1 ()          ## visualize X-projection
    >>> fy  = model.draw2 ()          ## visualize Y-projection
    >>> fz  = model.draw3 ()          ## visualize Z-projection

    Parameters
    ----------
    signalx  :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in X-direction
    signaly  :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in Y-direction
    signalz  :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in Z-direction
    suffix   : string
        An optional suffix to be  added to the names of created PDFs and variables
    bkg_1x    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for Bx1(x)* Sy(y)* Sz(z) term
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkg_1y    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for Sx(z)*By1(y)* Sz(z) term
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkg_1x    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for   Sx(x)* Sy(y)*Bz1(z) term 
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkg_2xy    : RooFit/PDF, Ostap/PDF or None
        2D-background for Bxy(x,y)*Sz(z) term  
        Use directly RooFit/PDF or Ostap/PDF otherwise create from bkgX2 and bkgY2
    bkg_2xz    : RooFit/PDF, Ostap/PDF or None
        2D-background for Bxz(x,z)*Sy(y) term 
        Use directly RooFit/PDF or Ostap/PDF otherwise create from bkgX2 and bkgZ2
    bkg_2yz    : RooFit/PDF, Ostap/PDF or None
        2D-background for Bxz(x,z)*Sy(y) term        
        Use directly RooFit/PDF or Ostap/PDF otherwise create from bkgX2 and bkgZ2
    bkg_2x    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for used to create Bxy(x,y) and Bxz(x,z) if they are not specified.
        If None - bkgX1 is used 
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkg_2y    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for used to create Bxy(x,y) and Byz(y,z) if they are not specified.
        If None - bkgY1 is used 
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkg_2z    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for used to create Bxz(x,z) and Byz(y,z) if they are not specified.
        If None - bkgZ1 is used
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkg_3x    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for used to create Bxyz(x,y,z) if it is not specified.
        If None - bkgX2 is used 
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkg_3y   : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for used to create Bxyz(x,y,z) if it is not specified.
        If None - bkgY2 is used 
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkg_3z   : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for used to create Bxyz(x,y,z) if it is not specified.
        If None - bkgZ2 is used         
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkg_3D   : RooFit/PDF or Ostap/PDF 3to descrive 3D-backround component Bxyz(x,y,z)
    

    sss      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of sig(1) * sig(2) * sig(3) component
         Use directly RooRealVar, otherwise create it using self.make_var function
    ssb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of sig(1) * sig(2) * bkg (3) component
         Use directly RooRealVar, otherwise create it using self.make_var function
    sbs      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of sig(1) * bkg(2) * sig (3) component
         Use directly RooRealVar, otherwise create it using self.make_var function
    bss      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of bkg(1) * sig(2) * sig (3) component
         Use directly RooRealVar, otherwise create it using self.make_var function
    sbb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of sig(1) * bkg(2) * bkg (3) component
         Use directly RooRealVar, otherwise create it using self.make_var function
    bsb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of bkg(1) * sig(2) * bkg (3) component
         Use directly RooRealVar, otherwise create it using self.make_var function
    bbs      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of bkg(1) * bkg(2) * sig (3) component
         Use directly RooRealVar, otherwise create it using self.make_var function
    bbb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of bkg(1) * bkg(2) * bkg (3) component
         Use directly RooRealVar, otherwise create it using self.make_var function

    """
    def __init__ ( self               ,
                   #
                   signal_x           , 
                   signal_y           ,
                   signal_z           ,
                   suffix = ''        ,
                   #
                   bkg_1x     = None  , ## 1D-background for Bx(x)*Sy(y)*Sz(z) term
                   bkg_1y     = None  , ## 1D-background for Sx(z)*By(y)*Sz(z) term
                   bkg_1z     = None  , ## 1D-background for Sx(x)*Sy(y)*Bz(z) term 
                   #
                   bkg_2xy    = None  , ## 2D-background for Bxy(x,y)*Sz(z) term 
                   bkg_2xz    = None  , ## 2D-background for Bxz(x,z)*Sy(y) term 
                   bkg_2yz    = None  , ## 2D-background for Byz(y,z)*Sx(z) term  
                   ## *if* no XY,XZ,BC backgrounds are specified, combine them from 
                   bkg_2x      = None  , ## Bxy(x,y) = Bx2(x)*By2(y)
                   bkg_2y      = None  , ## Bxz(x,z) = Bx2(x)*Bz2(z)
                   bkg_2z      = None  , ## Bkg(y,z) = By2(y)*Bz2(z)
                   ##
                   bkg_3D     = None  , ## 3D-backround component  B(x,y,z)
                   ## *if* no 3D-background components is specified, combine it from
                   bkg_3x     = None  , ## Bkg(x,y,z) = Bx3(x)*By3(y)*Bz3(z)
                   bkg_3y     = None  , ## Bkg(x,y,z) = Bx3(x)*By3(y)*Bz3(z)
                   bkg_3z     = None  , ## Bkg(x,y,z) = Bx3(x)*By3(y)*Bz3(z)
                   #
                   ## Yields of the main components :
                   sss        = None  , ## sig(x) * sig(y) * sig(z) 
                   ssb        = None  , ## sig(x) * sig(y) * bkg(z)
                   sbs        = None  , ## sig(x) * bkg(y) * sig(z)
                   bss        = None  , ## bkg(x) * sig(y) * sig(z)
                   sbb        = None  , ## sig(x) * bkg(y,z)
                   bsb        = None  , ## sig(y) * bkg(x,z)
                   bbs        = None  , ## sig(z) * bkg(x,y)
                   bbb        = None  , ## background-3D 
                   ## additional components 
                   components = []    ,
                   xvar       = None  ,
                   yvar       = None  ,
                   zvar       = None  ,                   
                   name       = ''    ) : 

        ## keep all arguments 
        self.__args = {
            ##
            'signal_x'   : signal_x ,
            'signal_y'   : signal_y ,
            'signal_z'   : signal_z ,
            ##
            'bkg_1x'     : bkg_1x   ,
            'bkg_1y'     : bkg_1y   ,
            'bkg_1z'     : bkg_1z   ,
            ## 
            'bkg_2xy'    : bkg_2xy  ,
            'bkg_2xz'    : bkg_2xz  ,
            'bkg_2yz'    : bkg_2yz  ,
            ##
            'bkg_2x'     : bkg_2x   ,
            'bkg_2y'     : bkg_2y   ,
            'bkg_2z'     : bkg_2z   ,
            ## 
            'bkg_3D'     : bkg_3D   ,
            ##
            'bkg_3x'     : bkg_3x   ,
            'bkg_3y'     : bkg_3y   ,
            'bkg_3z'     : bkg_3z   ,
            ## 
            'sss'        : sss      ,            
            'ssb'        : ssb      ,
            'sbs'        : sbs      ,
            'bss'        : bss      ,
            'sbb'        : sbb      ,
            'bsb'        : bsb      ,
            'bbs'        : bbs      ,
            'bbb'        : bbb      ,
            ##
            'components' : components ,
            'xvar'       : xvar     , 
            'yvar'       : yvar     , 
            'zvar'       : zvar     , 
            ##
            'name'     : name     
            }
        
        self.__suffix      = suffix


        self.__signal_x , xvar = self.make_PDF1 ( signal_x , xvar , prefix = 'SX' , suffix = suffix )
        self.__signal_y , yvar = self.make_PDF1 ( signal_y , yvar , prefix = 'SY' , suffix = suffix )
        self.__signal_z , zvar = self.make_PDF1 ( signal_z , zvar , prefix = 'SZ' , suffix = suffix )
            
        #
        ## initialize base class
        #
        if not name : 
            name = "%s&%s&%s" % ( self.__signal_x.name ,
                                  self.__signal_y.name ,
                                  self.__signal_z.name )                                  
            if suffix : name += '_' + suffix 
            
        PDF3.__init__ ( self , name           ,
                        self.__signal_x.xvar  ,
                        self.__signal_y.xvar  ,
                        self.__signal_z.xvar  ) 
     
        # =====================================================================
        ## 1) First component: all   signals
        # =====================================================================
        
        self.__sss_cmp  = Model3D ( name   = self.new_name ( 'SSS' , suffix ) ,
                                    xmodel = self.__signal_x ,
                                    ymodel = self.__signal_y ,
                                    zmodel = self.__signal_z )
        
        # =====================================================================
        ## 2-4) Three terms:  ( 2 signals )  x ( 1 background ) 
        # =====================================================================
        
        self.__bkg_1x    = self.make_bkg ( bkg_1x , self.new_name ( 'Bkg1X_BSS' , suffix ) , self.xvar )
        self.__bkg_1y    = self.make_bkg ( bkg_1y , self.new_name ( 'Bkg1Y_SBS' , suffix ) , self.yvar )
        self.__bkg_1z    = self.make_bkg ( bkg_1z , self.new_name ( 'Bkg1Z_SSB' , suffix ) , self.zvar )
        
        self.__ssb_cmp = Model3D ( self.new_name ( "SSB" , suffix ) ,
                                   self.__signal_x , self.__signal_y , self.__bkg_1z   ,
                                   title = "Signal(x) x Signal(y) x Background1(x)" )
        self.__sbs_cmp = Model3D ( self.new_name ( "SBS" , suffix ) ,
                                   self.__signal_x , self.__bkg_1y   , self.__signal_z , 
                                   title = "Signal(x) x Background1(y) x Signal(z)" ) 
        self.__bss_cmp = Model3D ( self.new_name ( "BSS" , suffix ) ,
                                   self.__bkg_1x   , self.__signal_y , self.__signal_z ,
                                   title = "Background1(x) x Signal(y) x Signal(z)" )

        
        # =====================================================================
        ## (intermezzo-1) Assumptions about SBB-background sub-components 
        # =====================================================================
        
        if   component_clone   ( bkg_2x ) :
            bkg_2x = self.__bkg_1x
            self.debug ( 'bkg_2x set to [CLONE]   %s' % bkg_2x ) 
        elif component_similar ( bkg_2x ) :
            bkg_2x =        bkg_1x
            self.debug ( 'bkg_2x set to [SIMILAR] %s' % bkg_2x ) 

        if   component_clone   ( bkg_2y ) :
            bkg_2y = self.__bkg_1y
            self.debug ( 'bkg_2y set to [CLONE]   %s' % bkg_2y )
        elif component_similar ( bkg_2x ) :
            bkg_2y =        bkg_1y
            self.debug ( 'bkg_2y set to [SIMILAR] %s' % bkg_2y ) 
            
        if   component_clone   ( bkg_2z ) :
            bkg_2z = self.__bkg_1z
            self.debug ( 'bkg_2z set to [CLONE]   %s' % bkg_2z ) 
        elif component_similar ( bkg_2z ) :
            bkg_2z =        bkg_1z
            self.debug ( 'bkg_2z set to [SIMILAR] %s' % bkg_2z ) 
        # =====================================================================
        
        from ostap.fitting.pdf_ops import Prod3D_pdf 
        
        # =====================================================================
        ## 5-7) Three terms: (1 signal) x (2 backgrounds)
        # =====================================================================
                
        self.__bkg_2x = None 
        self.__bkg_2y = None 
        self.__bkg_2z = None 

        bkg_2xy_name = self.new_name ( 'Bkg2XY' , suffix )
        bkg_2xz_name = self.new_name ( 'Bkg2XZ' , suffix )
        bkg_2yz_name = self.new_name ( 'Bkg2YZ' , suffix )
        
        if bkg_2xy and isinstance ( bkg_2xy , ( tuple , list ) ) :            
            from ostap.fitting.models_2d import make_B2D
            self.__bkg_2xy = make_B2D ( bkg2sy_name , self.xvar , self.yvar , *bkg_2xy )
        elif bkg_2xy :
            self.__bkg_2xy = self.make_PDF2 ( bkg_2xy , xvar = self.xvar , yvar = self.yvar , prefix = 'BKg2XY' , suffix = suffix )  [ 0 ]
        else :
            if not self.__bkg_2x : self.__bkg_2x = self.make_bkg ( bkg_2x , self.new_name ( 'Bkg2X_S2B' , suffix ) , self.xvar )    
            if not self.__bkg_2y : self.__bkg_2y = self.make_bkg ( bkg_2y , self.new_name ( 'Bkg2Y_S2B' , suffix ) , self.yvar )                    
            self.__bkg_2xy = Model2D ( bkg_2xy_name    ,
                                       self.__bkg_2x   ,
                                       self.__bkg_2y   ,
                                       title =  'Background2(x) x Background2(y)' )
            
        if bkg_2xz and isinstance ( bkg_2xz , ( tuple , list ) ) :            
            from ostap.fitting.models_2d import make_B2D
            self.__bkg_2xz = make_B2D ( bkg_2xz_name , self.xvar , self.zvar , *bkg_2xz )
        elif bkg_2xz :
            self.__bkg_2xz = self.make_PDF2 ( bkg_2xz , xvar = self.xvar , yvar = self.zvar , prefix = 'BKg2XZ' , suffix = suffix ) [ 0 ]
        else :
            if not self.__bkg_2x : self.__bkg_2x = self.make_bkg ( bkg_2x , self.new_name ( 'Bkg2X_S2B' , suffix ) , self.xvar )    
            if not self.__bkg_2z : self.__bkg_2z = self.make_bkg ( bkg_2z , self.new_name ( 'Bkg2Z_S2B' , suffix ) , self.zvar )                                
            self.__bkg_2xz = Model2D ( bkg_2xz_name   ,                                     
                                       self.__bkg_2x  ,
                                       self.__bkg_2z  , 
                                       title =  'Background2(x) x Background2(z)' )
            
        if bkg_2yz and isinstance ( bkg_2yz , ( tuple , list ) ) :            
            from ostap.fitting.models_2d import make_B2D
            self.__bkg_2yz = make_B2D ( bkg_2yz_name , self.yvar , self.zvar , *bkg_2yz )
        elif bkg_2yz :
            self.__bkg_2yz = self.make_PDF2 ( bkg_2yz , xvar = self.yvar , yvar = self.zvar , prefix = 'Bkg2YZ' , suffix = suffix ) [ 0 ] 
        else :            
            if not self.__bkg_2y : self.__bkg_2y = self.make_bkg ( bkg_2y , self.new_name ( 'Bkg2Y_S2B' , suffix ) , self.yvar )    
            if not self.__bkg_2z : self.__bkg_2z = self.make_bkg ( bkg_2z , self.new_name ( 'Bkg2Z_S2B' , suffix ) , self.zvar )                                
            self.__bkg_2yz = Model2D ( bkg_2yz_name   ,
                                       self.__bkg_2y  ,
                                       self.__bkg_2z  ,
                                       title =  'Background2(y) x Background2(z)' )

            
        ## create components
        self.__sbb_cmp = Prod3D_pdf ( ( self.__signal_x , self.__bkg_2yz ) , 
                                      xvar = self.xvar , yvar = self.yvar , zvar = self.zvar ,
                                      name = self.new_name ( 'SBB' , "S(x)*B(y,z)" ) ) 
        self.__bsb_cmp = Prod3D_pdf ( ( self.__signal_y , self.__bkg_2xz ) ,
                                      xvar = self.xvar , yvar = self.yvar , zvar = self.zvar ,
                                      name = self.new_name ( 'BSB' , "S(y)*B(x,z)" ) ) 
        self.__bbs_cmp = Prod3D_pdf ( ( self.__signal_z , self.__bkg_2xy ) ,
                                      xvar = self.xvar , yvar = self.yvar , zvar = self.zvar ,
                                      name = self.new_name ( 'BBS' , "S(z)*B(x,y)" ) ) 
        
        # =====================================================================
        ## (intermezzo-2) Assumptions about BBB-background sub-components 
        # =====================================================================

        if   component_clone   ( bkg_3x ) :
            bkg_3x = self.__bkg_2x
            self.debug ( 'bkg_3x set to [CLONE]   %s' % bkg_3x ) 
        elif component_similar ( bkg_3x ) :
            bkg_3x =        bkg_2x
            self.debug ( 'bkg_3x set to [SIMILAR] %s' % bkg_3x ) 

        if   component_clone   ( bkg_3y ) :
            bkg_3y = self.__bkg_2y
            self.debug ( 'bkg_3y set to [CLONE]   %s' % bkg_3y  )
        elif component_similar ( bkg_3x ) :
            bkg_3y =        bkg_2y
            self.debug ( 'bkg_3y set to [SIMILAR] %s' % bkg_3y ) 
            
        if   component_clone   ( bkg_3z ) :
            bkg_3z = self.__bkg_2z
            self.debug ( 'bkg_3z set to [CLONE]   %s' % bkg_3z ) 
        elif component_similar ( bkg_3x ) :
            bkg_3z =        bkg_2z
            self.debug ( 'bkg_3z set to [SIMILAR] %s' % bkg_3z ) 

        # =====================================================================
        ## 8) pure background 
        # =====================================================================
        
        self.__bkg_3x = None 
        self.__bkg_3y = None 
        self.__bkg_3z = None 
        
        bbb_name = self.new_name ( 'BBB' , suffix ) 
        if bkg_3D and isinstance ( bkg_3D , (  tuple , list ) ) :
            from ostap.fitting.models_3d import make_B3D 
            self.__bbb_cmp = make_B3D ( bbb_name , self.xvar , self.yvar , self.zvar , *bkg_3D )
        elif bkg_3D :
            self.__bbb_cmp = self.make_PDF3 ( bkg_3D , xvar = self.xvar , yvar = self.yvar , zvar = zvar , prefix = 'BBB' , suffix = suffix ) [ 0 ]
        else :
            
            self.__bkg_3x = self.make_bkg ( bkg_3x , self.new_name ( 'Bkg3X_BBB' , suffix ) , self.xvar )
            self.__bkg_3y = self.make_bkg ( bkg_3y , self.new_name ( 'Bkg3Y_BBB' , suffix ) , self.yvar )
            self.__bkg_3z = self.make_bkg ( bkg_3z , self.new_name ( 'Bkg3Z_BBB' , suffix ) , self.zvar )
            
            self.__bbb_cmp = Model3D ( bbb_name ,
                self.__bkg_3x ,
                self.__bkg_3y ,
                self.__bkg_3z ,
                title = "Background3(x) x Backrgound3(y) x Background3(z)" )
        #
        ## coefficients
        #
        self.__sss = self.make_var ( sss   , "SSS"          + suffix ,
                                     "Signal(x)&Signal(y)&Signal(z)"     + suffix , None , 1000  , 0 ,  1.e+7 )
        self.__ssb = self.make_var ( ssb   , "SSB"          + suffix ,
                                     "Signal(x)&Signal(y)&Background(z)" + suffix , None , 1000  , 0 ,  1.e+7 )
        self.__sbs = self.make_var ( sbs   , "SBS"          + suffix ,
                                     "Signal(x)&Background(y)&Signal(z)" + suffix , None , 1000  , 0 ,  1.e+7 )
        self.__bss = self.make_var ( bss   , "BSS"          + suffix ,
                                     "Background(x)&Signal(y)&Signal(z)" + suffix , None , 1000  , 0 ,  1.e+7 )
        self.__sbb = self.make_var ( sbb  , "SBB"           + suffix ,
                                     "Signal(x)&Background(y,z)"         + suffix , None , 1000  , 0 ,  1.e+7 )
        self.__bsb = self.make_var ( bsb  , "BSB"           + suffix ,
                                     "Signal(y)&Background(x,z)"         + suffix , None , 1000  , 0 ,  1.e+7 )
        self.__bbs = self.make_var ( bbs  ,  "BBS"          + suffix ,
                                     "Signal(z)&Background(x,y)"         + suffix , None , 1000  , 0 ,  1.e+7 )
        self.__bbb = self.make_var ( bbb  , "BBB"           + suffix ,
                                     "Background(x,y,z)"                 + suffix , None , 1000  , 0 ,  1.e+7 )
        
        self.alist1 = ROOT.RooArgList (
            self.__sss_cmp.pdf ,
            self.__ssb_cmp.pdf ,
            self.__sbs_cmp.pdf ,
            self.__bss_cmp.pdf ,
            self.__sbb_cmp.pdf ,
            self.__bsb_cmp.pdf ,
            self.__bbs_cmp.pdf ,
            self.__bbb_cmp.pdf )
        self.alist2 = ROOT.RooArgList (
            self.__sss     ,
            self.__ssb     ,
            self.__sbs     ,
            self.__bss     ,
            self.__sbb     ,
            self.__bsb     ,
            self.__bbs     ,
            self.__bbb     )

        ## treat additional components (if specified)
        self.__nums_components = [] 
        icmp = 0
        self.__more_components = []
        for i , cmp in enumerate ( components ) :
            
            if   isinstance  ( cmp , PDF3           ) : cc = cmp  
            elif isinstance  ( cmp , ROOT.RooAbsPdf ) :
                cc = Generic3D_pdf ( cmp ,  self.xvar , self.yvar, self.zvar , prefix = 'C%d_' % i , suffix = suffix ) 
            else :
                self.error ("unknown 'other' component %s/%s, skip it!" % ( cc , type(cc) ) )
                continue  
            self.__more_components.append ( cc     )
            self.components.add           ( cc.pdf ) 

        nc = len( self.__more_components )
        if 1 == nc :
            cf = self.make_var ( None , "C"+suffix , "Component" + suffix , None , 1 , 0 , 1.e+7 )
            self.alist1.add  ( self.components[0] )
            self.__nums_components.append ( cf ) 
        elif 2 <= nc : 
            fic = self.make_fracs ( nc , 'C_%%d%s' % suffix ,  'C(%%d)%s'  % suffix , fractions  = False )
            for c in self.components : self.alist1.add ( c)
            for f in fic             : self.__nums_components.append ( f )
            
        self.__nums_components  = tuple ( self.__nums_components  ) 
        for c in self.__nums_components  : self.alist2.add ( c )

        ## 
        #
        ## build the final PDF 
        #
        pdfname  = self.new_roo_name ( "fit3d" , suffix )
        pdftitle = "Fit3D %s" % self.name
        pdfargs  = pdfname , pdftitle , self.alist1 , self.alist2
        self.pdf = ROOT.RooAddPdf  ( *pdfargs )

        self.signals     .add ( self.__sss_cmp.pdf )
        self.backgrounds .add ( self.__bbb_cmp.pdf )
        self.crossterms1 .add ( self.__ssb_cmp.pdf ) ## cross-terms
        self.crossterms1 .add ( self.__sbs_cmp.pdf ) ## cross-terms
        self.crossterms1 .add ( self.__bss_cmp.pdf ) ## cross-terms
        self.crossterms2 .add ( self.__sbb_cmp.pdf ) ## cross-terms 
        self.crossterms2 .add ( self.__bsb_cmp.pdf ) ## cross-terms 
        self.crossterms2 .add ( self.__bbs_cmp.pdf ) ## cross-terms 

        ## save the configuration
        self.config = {
            'signal_x'   : self.signal_x ,
            'signal_y'   : self.signal_y ,
            'signal_z'   : self.signal_z ,
            'suffix'     : self.suffix  ,
            ##
            'bkg_1x'     : self.bkg_1x  ,
            'bkg_1y'     : self.bkg_1y  ,
            'bkg_1z'     : self.bkg_1z  ,
            ##
            'bkg_2x'     : self.bkg_2x  ,
            'bkg_2y'     : self.bkg_2y  ,
            'bkg_2z'     : self.bkg_2z  ,
            ##
            'bkg_2xy'    : self.bkg_2xy ,
            'bkg_2xz'    : self.bkg_2xz ,
            'bkg_2yz'    : self.bkg_2yz ,
            ##
            'bkg_3x'     : self.bkg_3x  ,
            'bkg_3y'     : self.bkg_3y  ,
            'bkg_3z'     : self.bkg_3z  ,
            ##
            'bkg_3D'     : self.bkg_3D  ,
            ##
            'sss'        : self.SSS     ,
            'ssb'        : self.SSB     ,
            'sbs'        : self.SBS     ,
            'bss'        : self.BSS     ,
            'sbb'        : self.SBB     ,
            'bsb'        : self.BSB     ,
            'bbs'        : self.BBS     ,
            'bbb'        : self.BBB     ,
            #
            'components' : self.more_components ,            
            'xvar'       : self.xvar    ,
            'yvar'       : self.yvar    ,
            'zvar'       : self.zvar    ,
            'name'       : self.name    ,             
            }

        self.checked_keys.add  ( 'xvar' )
        self.checked_keys.add  ( 'yvar' )
        self.checked_keys.add  ( 'zvar' )
        
    @property
    def SSS ( self ) :
        """The yield of Signal(x)*Signal(y)*Signal(z) component"""
        return self.__sss
    @SSS.setter 
    def SSS ( self , value ) :
        value = float ( value  )
        assert value in self.__sss, "Value %s is out of the allowed range %s " % ( value , self.__sss.minmax() )
        self.__sss.setVal ( value ) 

    @property
    def SSB ( self ) :
        """The yield of Signal(x)*Signal(y)*Background(z) component"""
        return self.__ssb
    @SSB.setter 
    def SSB ( self , value ) :
        value = float ( value  )
        assert value in self.__ssb, "Value %s is out of the allowed range %s " % ( value , self.__ssb.minmax() )
        self.__ssb.setVal ( value ) 

    @property
    def SBS ( self ) :
        """The yield of Signal(x)*Background(y)*Signal(z) component"""
        return self.__sbs
    @SBS.setter 
    def SBS ( self , value ) :
        value = float ( value  )
        assert value in self.__sbs, "Value %s is out of the allowed range %s " % ( value , self.__sbs.minmax() )
        self.__sbs.setVal ( value ) 

    @property
    def BSS ( self ) :
        """The yield of Background(x)*Signal(y)*Signal(z) component"""
        return self.__bss
    @BSS.setter 
    def BSS ( self , value ) :
        value = float ( value  )
        assert value in self.__bss, "Value %s is out of the allowed range %s " % ( value , self.__bss.minmax() )
        self.__bss.setVal ( value ) 

    @property
    def SBB ( self ) :
        """The yield of Signal(x)*Background(y,z) component"""
        return self.__sbb
    @SBB.setter 
    def SBB ( self , value ) :
        value = float ( value  )
        assert value in self.__sbb, "Value %s is out of the allowed range %s " % ( value , self.__sbb.minmax() )
        self.__sbb.setVal ( value ) 

    @property
    def BSB ( self ) :
        """The yield of Background(x,z)*Signal(y) component"""
        return self.__bsb
    @BSB.setter 
    def BSB ( self , value ) :
        value = float ( value  )
        assert value in self.__bsb, "Value %s is out of the allowed range %s " % ( value , self.__bsb.minmax() )
        self.__bsb.setVal ( value ) 

    @property
    def BBS ( self ) :
        """The yield of Background(x,y)*Signal(z) component"""
        return self.__bbs
    @BBS.setter 
    def BBS ( self , value ) :
        value = float ( value  )
        assert value in self.__bbs, "Value %s is out of the allowed range %s " % ( value , self.__bbs.minmax() )
        self.__bbs.setVal ( value ) 

    @property
    def BBB ( self ) :
        """The yield of Background(x,y,z) component"""
        return self.__bbb
    @BBB.setter 
    def BBB ( self , value ) :
        value = float ( value  )
        assert value in self.__bbb, "Value %s is out of the allowed range %s " % ( value , self.__bbb.minmax() )
        self.__bbb.setVal ( value ) 


    @property
    def C ( self ) :
        """Get the  yields of 'other' component(s) 
        For single 'other' component:
        >>> print pdf.C           ## read the single 'other' component 
        >>> pdf.C = 100           ## assign to it 
        For multiple 'other' components:
        >>> print pdf.C[4]        ## read the 4th 'other' component 
        >>> pdf.C = 4,100         ## assign to it 
        ... or, alternatively:
        >>> print pdf.C[4]        ## read the 4th 'other' component 
        >>> pdf.C[4].value 100    ## assign to it         
        """
        return self.component_getter ( self.__nums_components  )     
    @C.setter
    def C (  self , value ) :        
        self.component_setter ( self.__nums_components , value )
   
    @property
    def yields    ( self ) :
        """The list/tuple of the yields of all numeric components"""
        return tuple ( [ i for i in  self.alist2 ] )

    @property 
    def total_yield ( self ) :
        """'total_yield' : get the total yield"""
        if not self.fit_result                   : return None
        if not valid_pointer ( self.fit_result ) : return None
        return self.fit_result.sum ( *self.yields ) 
    
    # =========================================================================
    # components
    # =========================================================================
    
    @property 
    def signal_x ( self  ) :
        """'signal_x' : Signal(x) component/PDF"""
        return self.__signal_x

    @property 
    def signal_y ( self  ) :
        """'signal_y' : Signal(y) component/PDF"""
        return self.__signal_y
    
    @property 
    def signal_z ( self  ) :
        """'signal_z' : Signal(z) component/PDF"""
        return self.__signal_z

    @property
    def bkg_1x ( self ) :
        """'bkg_1x' : B(x) component for B(x)*S(y)*S(z) term"""
        return self.__bkg_1x
    @property
    def bkg_1y ( self ) :
        """'bkg_1y' : B(y) component for S(x)*B(y)*S(z) term"""
        return self.__bkg_1y
    @property
    def bkg_1z ( self ) :
        """'bkg_1z' : B(z) component for S(x)*S(y)*B(z) term"""
        return self.__bkg_1z

    @property
    def bkg_2xy ( self ) :
        """'bkg_2xy' : B(x,y) component for B(x,y)*S(z) term"""
        return self.__bkg_2xy
    @property
    def bkg_2xz ( self ) :
        """'bkg_2xz' : B(x,z) component for B(x,z)*S(y) term"""
        return self.__bkg_2xz
    @property
    def bkg_2yz ( self ) :
        """'bkg_2yz' : B(y,z) component for B(y,z)*S(x) term"""
        return self.__bkg_2yz
    
    @property
    def bkg_2x ( self ) :
        """'bkg_2x' : B(x) component for B(x,y)*S(z) & B(x,z)*S(y) terms"""
        return self.__bkg_2x 
    @property
    def bkg_2y ( self ) :
        """'bkg_2y' : B(y) component for B(y,z)*S(x) & B(x,y)*S(z) terms"""
        return self.__bkg_2y 
    @property
    def bkg_2z ( self ) :
        """'bkg_2z' : B(z) component for B(x,z)*S(y) & B(y,z)*S(x) terms"""
        return self.__bkg_2z 

    @property
    def bkg_3x ( self ) :
        """'bkg_3x': B(x) component for B(x,y,z) term"""
        return self.__bkg_3x 
    @property
    def bkg_3y ( self ) :
        """'bkg_3y' : B(y) component for B(x,y,z) term"""
        return self.__bkg_3y 
    @property
    def bkg_3z ( self ) :
        """'bkg_3z' : B(z) component for B(z,y,z) term"""
        return self.__bkg_3z 

    @property
    def bkg_3D ( self ) :
        """'`bkg_3D' : B(x,y,z) component/PDF for the final PDF"""
        return self.__bbb_cmp 

    @property
    def cmp_SSS ( self ) :
        """'triple-signal':  component/PDF"""
        return self.__sss_cmp

    @property
    def cmp_SSB ( self ) :
        """'signal-signal-background'  component/PDF"""
        return self.__ssb_cmp
    
    @property
    def cmp_SBS ( self ) :
        """'signal-background-signal' component/PDF"""
        return self.__sbs_cmp
    
    @property
    def cmp_BSS ( self ) :
        """'background-signal-signal' component/PDF"""
        return self.__bss_cmp

    @property
    def cpm_SBB ( self ) :
        """'signal-background-background' component/PDF"""
        return self.__sbb_cmp
    
    @property
    def cmp_BSB ( self ) :
        """'background-signal-background' component/PDF"""
        return self.__bsb_cmp
    
    @property
    def cmp_BBS ( self ) :
        """'background-background-signal' component/PDF"""
        return self.__bbs_cmp

    @property
    def cmp_BBB ( self ) :
        """'triple-background' component/PDF"""
        return self.__bbb_cmp

    @property
    def more_components ( self ) :
        """'additional 'other' components"""
        return tuple( self.__more_components  )

    @property
    def suffix ( self ) :
        """'suffix' used to build the name"""
        return self.__suffix



# =============================================================================
## @class Fit3DSym
#  The actual model for fully symmetric 3D-fits
#
#  @param signal_x (RooFit/PDF or Ostap/PDF) PDF to describe (1D)-signal in X-direction
#  @param signal_y (RooFit/PDF, Ostap/PDF or None) PDF to describe (1D)-signal in Y-direction
#  @param signal_z (RooFit/PDF, Ostap/PDF or None) PDF to describe (1D)-signal in Z-direction
#  @param suffix   (string) An optional suffix to be  added to the names of created PDFs and variables
#  @param bkg_1x   (RooFit/PDF, Ostap/PDF,  integer, RooRealVar or None)
#                  1D x-background for SSB-terms
#                  Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
#  @param bkg_2x   (RooFit/PDF, Ostap/PDF,integer, RooRealVar or None)
#                  1D x-background for SBB-terms, if <code>bkg2D</code> is not specified. 
#                  Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
#  @param bkg_2xy  (RooFit/PDF or Ostap/PDF)
#                   2D x,y-background for SBB-terms
#                   Use directly RooFit/PDF or Ostap/PDF
#  @param bkg_3x   (RooFit/PDF, Ostap/PDF,integer, RooRealVar or None)
#                   1D x-background for BBB-term, if <code>bkg3D</code> is not specified. 
#                   Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
#  @param bkg_3D   (RooFit/PDF or Ostap/PDF)
#                   3D x,y,z-background for BBB-term
#                   Use directly RooFit/PDF or Ostap/PDF        
#  @param sss      (None, RooRealVar, non-negative float or tuple)
#                   Variable for the yield of SSS component.
#                   Use directly RooRelaVar, otherwise create it using self.make_var function
#  @param ssb      (None, RooRealVar, non-negative float or tuple)
#                   Variable for the yield of SSB component.
#                   Use directly RooRelaVar, otherwise create it using self.make_var function
#  @param sbb      (None, RooRealVar, non-negative float or tuple)
#                  Variable for the yield of SBB component.
#                  Use directly RooRelaVar, otherwise create it using self.make_var function
#  @param bbb      (None, RooRealVar, non-negative float or tuple)
#                  Variable for the yield of BBB component.
#                  Use directly RooRelaVar, otherwise create it using self.make_var function
#  @param components ([])   the list of additional 3D-PDFs to be used in the fit
#  @param xvar       (None) the x-variable
#  @param yvar       (None) the y-variable
#  @param zvar       (None) the z-variable
#  @param name       ("")   the PDF name 
#  @code
#  model   = Models.Fit3DSym (
#      signal_x = Models.Gauss_pdf ( 'Gx' , mass = m_x ) ,
#      signal_y = Models.Gauss_pdf ( 'Gy' , mass = m_y ) ,
#      signal_z = Models.Gauss_pdf ( 'Gz' , mass = m_z ) ,
#      bkg_1x   = -1 , 
#      bkg_2x   = -1 ,
#      bkg_3x   = -1 )
#
#  r = model.fitTo ( dataset ) ## fit dataset 
#
#  print r                       ## get results  
#
#  fx  = model.draw1 ()          ## visualize X-projection
#  fy  = model.draw2 ()          ## visualize Y-projection
#  fz  = model.draw3 ()          ## visualize Z-projection
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2017-07-25
class Fit3DSym (PDF3) :
    """The actual model for fully symmetric 3D-fits
    
    >>>  model   = Models.Fit3DSym (
    ...      signal_x = Models.Gauss_pdf ( 'Gx' , mass = m_x ) ,
    ...      signal_y = Models.Gauss_pdf ( 'Gy' , mass = m_y ) ,
    ...      signal_z = Models.Gauss_pdf ( 'Gz' , mass = m_z ) ,
    ...      bkg_1x   =  1 , 
    ...      bkg_2x   = -1 ,
    ...      bkg_3x   = -1 )
    >>> r,f = model.fitTo ( dataset ) ## fit dataset 
    >>> print r                       ## get results  
    >>> fx  = model.draw1 ()          ## visualize X-projection
    >>> fy  = model.draw2 ()          ## visualize Y-projection
    >>> fz  = model.draw3 ()          ## visualize Z-projection

    Parameters
    ----------
    signal_x :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in X-direction
    signal_y :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in Y-direction
    signal_z :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in Z-direction
    suffix   : string
        An optional suffix to be  added to the names of created PDFs and variables
    bkg_1x   : RooFit/PDF, Ostap/PDF, integer, RooRealVar or None
        1D x-background for SSB-terms
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
    bkg_2x  : RooFit/PDF, Ostap/PDF, integer, RooRealVar or None
        1D x-background for SBB-terms
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
    bkg_2xy : RooFit/PDF, Ostap/PDF, list/tuple or None 
        2D (x,y)-background for SBB-terms
        Use directly RooFit/PDF or Ostap/PDF
    bkg_3x  : RooFit/PDF, Ostap/PDF, integer integer, RooRealVar or None
        1D x-background for BBB-term
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
    bkg_3D  : RooFit/PDF, Ostap/PDF, list/tuple or None 
        3D x,y,z-background for BBB-term
        Use directly RooFit/PDF or Ostap/PDF
        
    sss      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of SSS component
         Use directly RooRelaVar, otherwise create it using self.make_var function
    ssb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of SSB component
         Use directly RooRelaVar, otherwise create it using self.make_var function
    sbb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of SBB component
         Use directly RooRelaVar, otherwise create it using self.make_var function
    bbb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of BBB component
         Use directly RooRelaVar, otherwise create it using self.make_var function
         
    components : the list of additional 3D-PDFs to be used in the fit
    xvar       : x-variable
    yvar       : y-variable
    zvar       : z-variable
    name       : the PDF name 
         
    """
    def __init__ ( self               ,
                   #
                   signal_x           , 
                   signal_y = None    ,
                   signal_z = None    ,
                   suffix   = ''      ,
                   # background  for SSB-terms:
                   bkg_1x     = None  , ## 1D     x-background for BSS-terms
                   # background for  SBB-terms:
                   bkg_2x     = None  , ## 1D     x-background for BBS-terms, if no bkg2d is specified   
                   bkg_2xy    = None  , ## 2D   x,y-background for BBS-terms     (symmetric)
                   # background for BBB-term
                   bkg_3x     = None  , ## 1D     x-background for BBB-term
                   bkg_3D     = None  , ## 3D x,y,z-background for B(x,y,z) term (symmetric) 
                   #
                   ## Yields of the main components :
                   sss        = None  , ## sig(1) * sig(2) * sig(3) 
                   ssb        = None  , ## SSB components 
                   sbb        = None  , ## SBB components 
                   bbb        = None  , ## background-3D 
                   ## additional components 
                   components = []    ,
                   xvar       = None  ,
                   yvar       = None  ,
                   zvar       = None  ,                   
                   name       = ''    ) : 
        
        ## keep all arguments 
        self.__args = {
            #
            'signal_x'   : signal_x   ,
            'signal_y'   : signal_y   ,
            'signal_z'   : signal_z   ,
            #
            'bkg_1x'     : bkg_1x     ,
            'bkg_2x'     : bkg_2x     ,
            'bkg_2xy'    : bkg_2xy    ,
            'bkg_3x'     : bkg_3x     ,
            'bkg_3D'     : bkg_3D     ,
            ##
            'sss'        : sss        ,            
            'ssb'        : ssb        ,
            'sbb'        : sbb        ,
            'bbb'        : bbb        ,
            ##
            'components' : components ,
            ##
            'xvar'       : xvar       ,
            'yvar'       : yvar       ,
            'zvar'       : zvar       ,
            ##
            'name'       : name                   
            }
        
        self.__suffix      = suffix
 

        self.__signal_x , xvar = self.make_PDF1 ( signal_x , xvar = xvar , prerfix = 'SX' , suffix = suffix )

        if yvar and not signal_y :
            self.__signal_y = self.__signal_x.clone ( xvar = yvar , name_prefix = 'SY_' )
            self.debug('signal y-component is cloned from the signal_x component')
        else :
            self.__signal_y , yvar = self.make_PDF1 ( signal_y , xvar = yvar , prerfix = 'SY' , suffix = suffix )

        if zvar and not signal_z :
            self.__signal_z = self.__signal_x.clone ( xvar = zvar , name_prefix = 'SZ_' )
            self.debug('signal z-component is cloned from the signal_x component')
        else : 
            self.__signal_z , zvar = self.make_PDF1 ( signal_z , xvar = zvar , prerfix = 'SZ' , suffix = suffix )
        
        #
        ## initialize base class
        #
        if not name :
            name = "%s&%s&%s" % ( self.__signal_x.name ,
                                  self.__signal_y.name ,
                                  self.__signal_z.name )                                  
            if suffix : name += '_' + suffix 
            
        PDF3.__init__ ( self , name          ,
                        self.__signal_x.xvar ,
                        self.__signal_y.xvar ,
                        self.__signal_z.xvar ) 
     
        # =====================================================================
        ## 1) First component: all   signals
        # =====================================================================
        
        self.__sss_cmp  = Model3D (
            'SSS_' + self.name , self.__signal_x , self.__signal_y , self.__signal_z )
        
        # =====================================================================
        ## 2) ( 2 signals )  x ( 1 background )
        # =====================================================================
        
        self.__bkg_1x = self.make_bkg (        bkg_1x , self.new_name ( 'Bkg1X_BSS' , suffix ) , self.xvar )
        self.__bkg_1y = self.make_bkg ( self.__bkg_1x , self.new_name ( 'Bkg1Y_SBS' , suffix ) , self.yvar )
        self.__bkg_1z = self.make_bkg ( self.__bkg_1x , self.new_name ( 'Bkg1Z_SSB' , suffix ) , self.zvar )
        
        self.__ssb_cmp_raw = Model3D ( self.new_name ( "SSB_raw" , suffix ) ,
                                       self.__signal_x , self.__signal_y , self.__bkg_1z   ,
                                       title = "Signal(x) x Signal(y) x Background1(z)" )
        self.__sbs_cmp_raw = Model3D ( self.new_name ( "SBS_raw" , suffix ) ,
                                       self.__signal_x , self.__bkg_1y   , self.__signal_z , 
                                       title = "Signal(x) x Background1(y) x Signal(z)" ) 
        self.__bss_cmp_raw = Model3D ( self.new_name ( "BSS_raw" , suffix ) ,
                                       self.__bkg_1x   , self.__signal_y , self.__signal_z ,
                                       title = "Background1(x) x Signal(y) x Signal(z)" )

        
        self.__ssb_cmp     = Generic3D_pdf (
            self.make_sum ( self.new_name ( "SSB" , suffix ) ,
                            "S(x)*S(y)*B(z)+S(x)*B(y)*S(z)+B(x)*S(y)*B(z)" ,
                            self.__ssb_cmp_raw.pdf ,
                            self.__sbs_cmp_raw.pdf ,
                            self.__bss_cmp_raw.pdf ) , 
            self.xvar , self.yvar , self.zvar )
        
        self.__sbs_cmp = self.__ssb_cmp
        self.__bss_cmp = self.__ssb_cmp
        

        # =====================================================================
        ## (intermezzo-1) Assumptions about the SBB-background sub-components 
        # =====================================================================        
        
        if   component_clone   ( bkg_2x ) :
            bkg_2x = self.__bkg_1x
            self.debug ( 'bkg_2x set to [CLONE]   %s' % bkg_2x ) 
        elif component_similar ( bkg_2x ) :
            bkg_2x =        bkg_1x
            self.debug ( 'bkg_2x set to [SIMILAR] %s' % bkg_2x ) 
        # =====================================================================

        
        from ostap.fitting.pdf_ops import Prod3D_pdf 
            
        # =====================================================================
        ## 3 Three terms: (1 signal) x (2 backgrounds)
        # =====================================================================
        
        self.__bkg_2x = None 
        self.__bkg_2y = None 
        self.__bkg_2z = None 

        if   bkg_2xy and isinstance ( bkg_2xy , PDF2 ) :

            self.__bkg_2xy = bkg_2xy
            self.__bkg_2xz = bkg_2xy.clone ( xvar = self.xvar , yvar = self.zvar , name_prefix = 'Bkg2XZ_' , name_suffix = suffix )
            self.__bkg_2yz = bkg_2xy.clone ( xvar = self.yvar , yvar = self.zvar , name_prefix = 'Bkg2YZ_' , name_suffix = suffix )
            
        elif bkg_2xy and isinstance ( bkg_2xy , ( tuple , list ) ) :

            from ostap.fitting.models_2d import make_B2Dsym 
            self.__bkg_2xy = make_B2Dsym ( self.new_name ( 'Bkg2XY' , suffix ) , self.xvar , self.yvar  , *bkg_2xy )
            self.__bkg_2xz = self.__bkg_2xy.clone ( name = self.new_name ( 'Bkg2XZ' , suffix ) ,
                                                    xvar = self.xvar             ,
                                                    yvar = self.zvar             )
            self.__bkg_2yz = self.__bkg_2xy.clone ( name = self.new_name ( 'Bkg2YZ' , suffix ) ,
                                                    xvar = self.yvar             ,
                                                    yvar = self.zvar             )              
        else :

            self.__bkg_2x  = self.make_bkg (        bkg_2x , self.new_name ( 'Bkg2X_S2B' , suffix ) , xvar = self.xvar )        
            self.__bkg_2y  = self.make_bkg ( self.__bkg_2x , self.new_name ( 'Bkg2Y_S2B' , suffix ) , xvar = self.yvar )        
            self.__bkg_2z  = self.make_bkg ( self.__bkg_2x , self.new_name ( 'Bkg2Z_S2B' , suffix ) , xvar = self.zvar )
            
            self.__bkg_2xy = Model2D ( self.new_name ( 'Bkg2XY' , suffix ) ,
                                       self.__bkg_2x     ,
                                       self.__bkg_2y     , title =  'Background2(x,y)' )
            self.__bkg_2xz = Model2D ( self.new_name ( 'Bkg2XZ' , suffix ) ,
                                       self.__bkg_2x     ,
                                       self.__bkg_2z     , title =  'Background2(x,z)' )
            self.__bkg_2yz = Model2D ( self.new_name ( 'Bkg2YZ' , suffix ) ,
                                       self.__bkg_2y     ,
                                       self.__bkg_2z     , title =  'Background2(y,z)' )

        ## make components
        self.__sbb_cmp_raw = Prod3D_pdf ( ( self.__signal_x , self.__bkg_2yz ) , 
                                          xvar = self.xvar , yvar = self.yvar , zvar = self.zvar ,
                                          name = self.new_name ( 'SBB_raw' ,  "S(x)*B(y,z)" ) )                 
        self.__bsb_cmp_raw = Prod3D_pdf ( ( self.__signal_y , self.__bkg_2xz ) ,
                                          xvar = self.xvar , yvar = self.yvar , zvar = self.zvar ,
                                          name = self.new_name ( 'BSB_raw' ,  "S(y)*B(x,z)" ) )
        self.__bbs_cmp_raw = Prod3D_pdf ( ( self.__signal_z , self.__bkg_2xy ) ,
                                          xvar = self.xvar , yvar = self.yvar , zvar = self.zvar ,
                                          name = self.new_name ( 'BBS_raw' ,  "S(z)*B(x,y)" ) )
        
        self.__sbb_cmp     = Generic3D_pdf ( self.make_sum ( self.new_name ( "SBB" , suffix ) ,
                                                             "S(x)*B(y,z)+S(y)*B(x,z)+S(Z)*B(x,y)" ,
                                                             self.__sbb_cmp_raw.pdf ,
                                                             self.__bsb_cmp_raw.pdf ,
                                                             self.__bbs_cmp_raw.pdf ) ,
                                             self.xvar , self.yvar , self.zvar )
        
        self.__bsb_cmp = self.__sbb_cmp
        self.__bbs_cmp = self.__sbb_cmp
        
        # =====================================================================
        ## (intermezzo-2) Assumptions about the BBB-background sub-components 
        # =====================================================================        
        if   component_clone   ( bkg_3x ) :
            bkg_3x = self.__bkg_2x
            self.debug ( 'bkg_3x set to [CLONE]   %s' % bkg_3x ) 
        elif component_similar ( bkg_3x ) :
            bkg_3x =        bkg_2x
            self.debug ( 'bkg_3x set to [SIMILAR] %s' % bkg_3x ) 
        # =====================================================================

        # =====================================================================
        ## 8) pure background 
        # =====================================================================
        
        self.__bkg_3x  = None 
        self.__bkg_3y  = None 
        self.__bkg_3z  = None 

        bbb_name = self.new_name ( 'BBB' , suffix )
        if bkg_3D and isinstance ( bkg_3D , (  tuple , list )  ) :
            from ostap.fitting.models_3d import make_B3Dsym
            self.__bbb_cmp = make_B3Dsym ( bbb_name , self.xvar , self.yvar , self.zvar , *bkg_3D )
        elif bkg_3D :
            self.__bbb_cmp = self.make_PDF3 ( bkg_3D ,
                                              xvar   = self.xvar ,
                                              yvar   = self.yvar ,
                                              zvar   = self.zvar ,
                                              prefix = 'BBB'     ,
                                              suffix = suffix    ) [ 0 ] 
        else :
            
            self.__bkg_3x  = self.make_bkg (        bkg_3x , self.new_name ( 'Bkg3X_BBB' , suffix ) , self.xvar )        
            self.__bkg_3y  = self.make_bkg ( self.__bkg_3x , self.new_name ( 'Bkg3Y_BBB' , suffix ) , self.yvar )        
            self.__bkg_3z  = self.make_bkg ( self.__bkg_3x , self.new_name ( 'Bkg3Z_BBB' , suffix ) , self.zvar )

            self.__bbb_cmp = Model3D ( bbb_name      ,
                                       self.__bkg_3x ,
                                       self.__bkg_3y ,
                                       self.__bkg_3z , title = "Background(x,y,z)" )
        
        #
        ## coefficients
        #
        self.__sss = self.make_var ( sss   , "SSS" + suffix ,
                                     "Signal(x)&Signal(y)&Signal(z)" + suffix , None , 1000 , 0 , 1.e+7 )
        self.__ssb = self.make_var ( ssb   , "SSB" + suffix ,
                                     "Signal*2&Background"           + suffix , None , 1000 , 0 , 1.e+7 )
        self.__sbb = self.make_var ( sbb   , "SBB" + suffix ,
                                     "Signal&Backrgound*2"           + suffix , None , 1000 , 0 , 1.e+7 )
        self.__bbb = self.make_var ( bbb  , "BBB"  + suffix ,
                                     "Background*3"                  + suffix , None , 1000 , 0 , 1.e+7 )
        
        self.__sbs = self.__ssb ## the same 
        self.__bss = self.__ssb ## the same 
        self.__bsb = self.__sbb ## the same 
        self.__bbs = self.__sbb ## the same         

        self.alist1 = ROOT.RooArgList (
            self.__sss_cmp.pdf ,
            self.__ssb_cmp.pdf ,
            self.__sbb_cmp.pdf ,
            self.__bbb_cmp.pdf )
        self.alist2 = ROOT.RooArgList (
            self.__sss     ,
            self.__ssb     ,
            self.__sbb     ,
            self.__bbb     )

        ## treat additional components (if specified)
        self.__nums_components = [] 
        icmp = 0
        self.__more_components = [
            self.make_PDF3 ( cmp ,
                             xvar   = self.xvar ,
                             yvar   = self.yvar ,
                             zvar   = self.zvar ,
                             prefix = 'C%d' % i ,
                             suffix = suffix    ) [ 0 ] for ( i, cmp ) in enumerate ( components ) ]
        
        for cmp in self.__more_components : 
            self.components.add ( cmp.pdf ) 

        nc = len( self.__more_components )
        if 1 == nc :
            cf = self.make_var ( None , "C"+suffix , "Component" + suffix , None , 1 , 0 , 1.e+7 )
            self.alist1.add  ( self.components[0] )
            self.__nums_components.append ( cf ) 
        elif 2 <= nc : 
            fic = self.make_fracs ( nc , 'C_%%d%s' % suffix ,  'C(%%d)%s'  % suffix , fractions  = False )
            for c in self.components : self.alist1.add ( c)
            for f in fic             : self.__nums_components.append ( f )
            
        self.__nums_components  = tuple ( self.__nums_components  ) 
        for c in self.__nums_components  : self.alist2.add ( c )

        ## 
        #
        ## build the final PDF 
        #
        pdfname  = self.new_roo_name ( "fit3ds" , suffix ) 
        pdftitle = "Fit3DSym %s" % self.name
        pdfargs  = pdfname , pdftitle , self.alist1 , self.alist2
        self.pdf = ROOT.RooAddPdf  ( *pdfargs )

        self.signals     .add ( self.__sss_cmp.pdf )
        self.backgrounds .add ( self.__bbb_cmp.pdf )
        self.crossterms1 .add ( self.__ssb_cmp.pdf ) ## cross-terms
        self.crossterms2 .add ( self.__sbb_cmp.pdf ) ## cross-terms 
        
        ## save the configuration
        self.config = {
            ##
            'signal_x'   : self.signal_x ,
            'signal_y'   : self.signal_y ,
            'signal_z'   : self.signal_z ,
            'suffix'     : self.suffix  ,
            ##
            'bkg_1x'     : self.bkg_1x  ,
            'bkg_2x'     : self.bkg_2x  ,
            'bkg_2xy'    : self.bkg_2xy ,
            'bkg_3x'     : self.bkg_3x  ,
            'bkg_3D'     : self.bkg_3D  ,
            ##
            'sss'        : self.SSS     ,
            'ssb'        : self.SSB     ,
            'sbb'        : self.SBB     ,
            'bbb'        : self.BBB     ,
            ##
            'components' : self.more_components ,
            ##
            'xvar'       : self.xvar    ,
            'yvar'       : self.yvar    ,
            'zvar'       : self.zvar    ,
            ##
            'name'       : self.name    ,             
            }
        
        self.checked_keys.add  ( 'xvar' )
        self.checked_keys.add  ( 'yvar' )
        self.checked_keys.add  ( 'zvar' )
        
    @property
    def SSS ( self ) :
        """The yield of Signal(x)*Signal(y)*Signal(z) component"""
        return self.__sss
    @SSS.setter 
    def SSS ( self , value ) :
        value = float ( value  )
        assert value in self.__sss, "Value %s is out of the allowed range %s " % ( value , self.__sss.minmax() )
        self.__sss.setVal ( value ) 

    @property
    def SSB ( self ) :
        """The yield of S(x)*S(y)*B(z)+S(x)*B(y)*S(z)+B(x)*S(y)*S(z) component (same as SBS, BSS)"""
        return self.__ssb
    @SSB.setter 
    def SSB ( self , value ) :
        value = float ( value  )
        assert value in self.__ssb, "Value %s is out of the allowed range %s " % ( value , self.__ssb.minmax() )
        self.__ssb.setVal ( value ) 

    @property
    def SBS ( self ) :
        """The yield of S(x)*S(y)*B(z)+S(x)*B(y)*S(z)+B(x)*S(y)*S(z) component (same as SSB, BSS)"""
        return self.__sbs
    @SBS.setter 
    def SBS ( self , value ) :
        value = float ( value  )
        assert value in self.__sbs, "Value %s is out of the allowed range %s " % ( value , self.__sbs.minmax() )
        self.__sbs.setVal ( value ) 

    @property
    def BSS ( self ) :
        """The yield of S(x)*S(y)*B(z)+S(x)*B(y)*S(z)+B(x)*S(y)*S(z) component (same as SSB, SBS)"""
        return self.__bss
    @BSS.setter 
    def BSS ( self , value ) :
        value = float ( value  )
        assert value in self.__bss, "Value %s is out of the allowed range %s " % ( value , self.__bss.minmax() )
        self.__bss.setVal ( value ) 

    @property
    def SBB ( self ) :
        """The yield of S(x)*B(y,z)+S(y)*B(x,z)+S(z)*B(y,z) component (same as BSB, BBS)"""
        return self.__sbb
    @SBB.setter 
    def SBB ( self , value ) :
        value = float ( value  )
        assert value in self.__sbb, "Value %s is out of the allowed range %s " % ( value , self.__sbb.minmax() )
        self.__sbb.setVal ( value ) 

    @property
    def BSB ( self ) :
        """The yield of S(x)*B(y,z)+S(y)*B(x,z)+S(z)*B(y,z) component (same as SBB, BBS)"""
        return self.__bsb
    @BSB.setter 
    def BSB ( self , value ) :
        value = float ( value  )
        assert value in self.__bsb, "Value %s is out of the allowed range %s " % ( value , self.__bsb.minmax() )
        self.__bsb.setVal ( value ) 

    @property
    def BBS ( self ) :
        """The yield of S(x)*B(y,z)+S(y)*B(x,z)+S(z)*B(y,z) component (same as SBB, BSB )"""
        return self.__bbs
    @BBS.setter 
    def BBS ( self , value ) :
        value = float ( value  )
        assert value in self.__bbs, "Value %s is out of the allowed range %s " % ( value , self.__bbs.minmax() )
        self.__bbs.setVal ( value ) 

    @property
    def BBB ( self ) :
        """The yield of Background(x,y,z) component"""
        return self.__bbb
    @BBB.setter 
    def BBB ( self , value ) :
        value = float ( value  )
        assert value in self.__bbb, "Value %s is out of the allowed range %s " % ( value , self.__bbb.minmax() )
        self.__bbb.setVal ( value ) 

    @property
    def C ( self ) :
        """Get the  yields of 'other' component(s) 
        For single 'other' component:
        >>> print pdf.C           ## read the single 'other' component 
        >>> pdf.C = 100           ## assign to it 
        For multiple 'other' components:
        >>> print pdf.C[4]        ## read the 4th 'other' component 
        >>> pdf.C = 4,100         ## assign to it 
        ... or, alternatively:
        >>> print pdf.C[4]        ## read the 4th 'other' component 
        >>> pdf.C[4].value 100    ## assign to it         
        """
        return self.component_getter ( self.__nums_components  )     
    @C.setter
    def C (  self , value ) :
        self.component_setter ( self.__nums_components , value )
   
    @property
    def yields    ( self ) :
        """The list/tuple of the yields of all numeric components"""
        return tuple ( [ i for i in  self.alist2 ] )

    @property 
    def total_yield ( self ) :
        """'total_yield'' : get the total yield"""
        if not self.fit_result                   : return None
        if not valid_pointer ( self.fit_result ) : return None
        return self.fit_result.sum ( *self.yields ) 
    
    # =========================================================================
    # components
    # =========================================================================
    
    @property 
    def signal_x ( self  ) :
        """'signal_x' : Signal(x) component/PDF"""
        return self.__signal_x

    @property 
    def signal_y ( self  ) :
        """'signal_y' : Signal(y) component/PDF"""
        return self.__signal_y
    
    @property 
    def signal_z ( self  ) :
        """'signal_z' : Signal(z) component/PDF"""
        return self.__signal_z

    @property
    def bkg_1x ( self ) :
        """'bkg_1x' : B(x) component for B(x)*S(y)*S(z) term"""
        return self.__bkg_1x
    @property
    def bkg_1y ( self ) :
        """'bkg_1y' : B(y) component for S(x)*B(y)*S(z) term"""
        return self.__bkg_1y
    @property
    def bkg_1z ( self ) :
        """'bkg_1z' : B(z) component for S(x)*S(y)*B(z) term"""
        return self.__bkg_1z
    
    @property
    def bkg_2xy( self ) :
        """'bkg_2xy' : B(x,y) component for B(x,y)*S(z) term"""
        return self.__bkg_2xy
    @property
    def bkg_2xz( self ) :
        """'bkg_2xz' : B(x,z) component for B(x,z)*S(y) term"""
        return self.__bkg_2xz
    @property
    def bkg_2yz( self ) :
        """'bkg_2yz' : B(y,z) component for B(y,z)*S(x) term"""
        return self.__bkg_2yz 

    @property
    def bkg_2x ( self ) :
        """'bkg_2x' : B(x) component for B(x,y)*S(z) & B(x,z)*S(y) terms"""
        return self.__bkg_2x
    @property
    def bkg_2y ( self ) :
        """'bkg_2y' : B(y) component for B(y,z)*S(x) & B(x,y)*S(z) terms"""
        return self.__bkg_2y
    @property
    def bkg_2z ( self ) :
        """'bkg_2z' : B(z) component for B(x,z)*S(y) & B(y,z)*S(x) terms"""
        return self.__bkg_2z

    @property
    def bkg_3x ( self ) :
        """'bkg_3x' : B(x) component for B(x)*B(y)*B(z) term"""
        return self.__bkg_3x
    @property
    def bkg_3y ( self ) :
        """'bkg_3y' : B(y) component for B(x)*B(y)*B(z) term"""
        return self.__bkg_3y
    @property
    def bkg_3z ( self ) :
        """'bkg_3z' : B(z) component for B(x)*B(y)*B(z) term"""
        return self.__bkg_3z

    @property
    def bkg_3D ( self ) :
        """'bkg_3D' : B(x,y,z) component/PDF for the final PDF"""
        return self.__bbb_cmp 

    @property
    def cmp_SSS ( self ) :
        """'`triple-signal' component/PDF"""
        return self.__sss_cmp

    @property
    def cmp_SSB ( self ) :
        """'`signal-signal-background' symmetrized component/PDF"""
        return self.__ssb_cmp
    
    @property
    def cmp_SBS ( self ) :
        """'`signal-background-signal' symmetrized component/PDF (same as cmp_SSB)"""
        return self.__sbs_cmp
    
    @property
    def cmp_BSS ( self ) :
        """'`background-signal-signal' symmetrized component/PDF (same as cmp_SSB)"""
        return self.__bss_cmp

    @property
    def cpm_SBB ( self ) :
        """'`signal-background-background' symmetrized component/PDF"""
        return self.__sbb_cmp
    
    @property
    def cmp_BSB ( self ) :
        """'`background-signal-background' symmetrized component/PDF (same as cmp_SBB)"""
        return self.__bsb_cmp
    
    @property
    def cmp_BBS ( self ) :
        """'`background-background-signal' symmetrized component/PDF (same as cmp_SBB)"""
        return self.__bbs_cmp

    @property
    def cmp_BBB ( self ) :
        """'`triple-background' component/PDF"""
        return self.__bbb_cmp

    @property
    def more_components ( self ) :
        """additional/'other' components"""
        return tuple( self.__more_components  )

    @property
    def suffix ( self ) :
        """'suffix', used to build the name"""
        return self.__suffix

    # =========================================================================
    ## Raw, non-symmetrized fit components/PDF (for debugging)
    # =========================================================================

    @property
    def cpm_raw_SSB ( self ) :
        """'signal-signal-background' raw,non-symmetrized component/PDF"""
        return self.__ssb_cmp_raw
    
    @property
    def cmp_raw_SBS ( self ) :
        """'signal-background-signal' raw,non-symmetrized component/PDF"""
        return self.__sbs_cmp_raw
    
    @property
    def cmp_raw_BSS ( self ) :
        """'background-signal-signal' raw,non-symmetrized component/PDF"""
        return self.__bss_cmp_raw

    @property
    def cpm_raw_SBB ( self ) :
        """'signal-background-background' raw,non-symmetrized component/PDF"""
        return self.__sbb_cmp_raw
    
    @property
    def cmp_raw_BSB ( self ) :
        """'background-signal-background' raw,non-symmetrized component/PDF"""
        return self.__bsb_cmp_raw
    
    @property
    def cmp_raw_BBS ( self ) :
        """'background-background-signal' raw,non-symmetrized component/PDF"""
        return self.__bbs_cmp_raw


# =============================================================================
## @class Fit3DMix
#  The actual model for fully "mixed-symemtry" 3D-fits  (symmetric for y<-->z)
#
#  @param signal_x (RooFit/PDF or Ostap/PDF) PDF to describe (1D)-signal in X-direction
#  @param signal_y (RooFit/PDF, Ostap/PDF or None) PDF to describe (1D)-signal in Y-direction
#  @param signal_z (RooFit/PDF, Ostap/PDF or None) PDF to describe (1D)-signal in Z-direction
#  @param suffix   (string) An optional suffix to be  added to the names of created PDFs and variables
#  @param bkg_1x   (RooFit/PDF, Ostap/PDF,  integer, RooRealVar or None)
#                  1D x-background for SSB-terms
#                  Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
#  @param bkg_1y   (RooFit/PDF, Ostap/PDF,integer, RooRealVar or None)
#                  1D x-background for SBB-terms, if <code>bkg2D</code> is not specified. 
#                  Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
#  @param bkg_2x   (RooFit/PDF or Ostap/PDF)
#                   2D x,y-background for SBB-terms
#                   Use directly RooFit/PDF or Ostap/PDF
#  @param bkg3     (RooFit/PDF, Ostap/PDF,integer, RooRealVar or None)
#                   1D x-background for BBB-term, if <code>bkg3D</code> is not specified. 
#                   Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
#  @param bkg3D    (RooFit/PDF or Ostap/PDF)
#                   3D x,y,z-background for BBB-term
#                   Use directly RooFit/PDF or Ostap/PDF        
#  @param sss      (None, RooRealVar, non-negative float or tuple)
#                   Variable for the yield of SSS component.
#                   Use directly RooRelaVar, otherwise create it using self.make_var function
#  @param ssb      (None, RooRealVar, non-negative float or tuple)
#                   Variable for the yield of SSB component.
#                   Use directly RooRelaVar, otherwise create it using self.make_var function
#  @param sbb      (None, RooRealVar, non-negative float or tuple)
#                  Variable for the yield of SBB component.
#                  Use directly RooRelaVar, otherwise create it using self.make_var function
#  @param bbb      (None, RooRealVar, non-negative float or tuple)
#                  Variable for the yield of BBB component.
#                  Use directly RooRelaVar, otherwise create it using self.make_var function
#  @param components ([])   the list of additional 3D-PDFs to be used in the fit
#  @param xvar       (None) the x-variable
#  @param yvar       (None) the y-variable
#  @param zvar       (None) the z-variable
#  @param name       ("")   the PDF name 
#  @code
#  model   = Models.Fit3DSym (
#      signal_x = Models.Gauss_pdf ( 'Gx' , mass = m_x ) ,
#      signal_y = Models.Gauss_pdf ( 'Gy' , mass = m_y ) ,
#      signal_z = Models.Gauss_pdf ( 'Gz' , mass = m_z ) ,
#      bkg_1x     = -1 , 
#      bkg_2x     = -1 ,
#      bkg_1y     = -1 )
#
#  r = model.fitTo ( dataset ) ## fit dataset 
#
#  print r                       ## get results  
#
#  fx  = model.draw1 ()          ## visualize X-projection
#  fy  = model.draw2 ()          ## visualize Y-projection
#  fz  = model.draw3 ()          ## visualize Z-projection
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2017-07-25
class Fit3DMix (PDF3) :
    """The actual model for y<->z-symmetric 3D-fits
    
    >>>  model   = Models.Fit3DSym (
    ...      signal_x = Models.Gauss_pdf ( 'Gx' , mass = m_x ) ,
    ...      signal_y = Models.Gauss_pdf ( 'Gy' , mass = m_y ) ,
    ...      signal_z = Models.Gauss_pdf ( 'Gz' , mass = m_z ) ,
    ...      bkg_1x   =  1 , 
    ...      bkg_2x   = -1 ,
    ...      bkg_1y   = -1 )
    >>> r,f = model.fitTo ( dataset ) ## fit dataset 
    >>> print r                       ## get results  
    >>> fx  = model.draw1 ()          ## visualize X-projection
    >>> fy  = model.draw2 ()          ## visualize Y-projection
    >>> fz  = model.draw3 ()          ## visualize Z-projection

    Parameters
    ----------
    signal_x :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in X-direction
    signal_y :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in Y-direction
    signal_z :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in Z-direction
    suffix   : string
        An optional suffix to be  added to the names of created PDFs and variables
    bkg_1x   : RooFit/PDF, Ostap/PDF, integer, RooRealVar or None
        1D x-background for BSS-terms
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
    bkg_1y   : RooFit/PDF, Ostap/PDF, integer, RooRealVar or None
        1D y-background for BSS-terms
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
    bkg_2x   : RooFit/PDF, Ostap/PDF, integer, RooRealVar or None
        1D x-background for BBS-terms
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
    bkg_2y   : RooFit/PDF, Ostap/PDF, integer, RooRealVar or None
        1D y-background for BBS-terms
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
    bkg_2xy  : RooFit/PDF, Ostap/PDF, list/tuple or None
        2D (x,y)-background for BBS-terms
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created 
    bkg_2yz  : RooFit/PDF, Ostap/PDF, list/tuple or None
        2D (y,z)-background for BBS-termsm must be symmetric! 
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created
    bkg_3x   : RooFit/PDF, Ostap/PDF, integer, RooRealVar or None
        1D x-background for BBB-terms
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
    bkg_3y   : RooFit/PDF, Ostap/PDF, integer, RooRealVar or None
        1D y-background for BBB-terms
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
    bkg_3yz  : RooFit/PDF, Ostap/PDF, list/tuple or None
        2D (y,z)-background for BBB-terms (must be symmetric!)
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created 
    bkg_3D   : RooFit/PDF, Ostap/PDF, list/tuple or None
        3D (x,y,z)-background for BBB-terms (must be symmetric for  y<-->z !)
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created 
        
    sss      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of SSS component
         Use directly RooRelaVar, otherwise create it using self.make_var function
    ssb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of SSB component
         Use directly RooRelaVar, otherwise create it using self.make_var function
    sbb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of SBB component
         Use directly RooRelaVar, otherwise create it using self.make_var function
    bbb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of BBB component
         Use directly RooRelaVar, otherwise create it using self.make_var function
         
    components : the list of additional 3D-PDFs to be used in the fit
    xvar       : x-variable
    yvar       : y-variable
    zvar       : z-variable
    name       : the PDF name 
         
    """
    def __init__ ( self               ,
                   #
                   signal_x           , 
                   signal_y           ,
                   signal_z = None    , ## cloned from signal_y if None 
                   suffix   = ''      ,
                   # background  for SSB-terms:
                   bkg_1x  = None  , ## 1D     x-background for BSS-terms
                   bkg_1y  = None  , ## 1D     y-background for BSS-terms,   z-background is cloned  
                   # background for  SBB-terms:
                   bkg_2x  = None  , ## 1D     x-background for SBB-terms
                   bkg_2y  = None  , ## 1D     y-background for SBB-terms,   z-background is cloned
                   # background for  SBB-terms:
                   bkg_2xy = None  , ## 2D   x,y-background for SBB-terms, (x,z)-backgroud is cloned 
                   bkg_2yz = None  , ## 3D   y,z-background for SBB-terms, must be SYMMETRIC!
                   # background for  BBB-terms:
                   bkg_3x  = None  , ## 1D     x-background for SBB-terms
                   bkg_3y  = None  , ## 1D     y-background for SBB-terms,   z-background is cloned
                   # background for  BBB-terms:
                   bkg_3yz = None  , ## 3D   y,z-background for SBB-terms, must be SYMMETRIC!
                   bkg_3D  = None  , ## 3D x,y,z-background for B(x,y,z) term (symmetric) 
                   ## Yields of the main components :
                   sss     = None  , ## S(x)*S(y)*S(z)             component 
                   bss     = None  , ## B(x)*S(y)*S(z)             component 
                   sbb     = None  , ## S(x)*B(y,z)                component 
                   ssb     = None  , ## S(x)*[S(y)*B(z)+B(y)*S(z)] component
                   bbs     = None  , ## B(x)*[S(y)*B(z)+B(y)*S(z)] component 
                   bbb     = None  , ## B(x,y,z)                   component    
                   ## additional components 
                   components = [] , ## the list additional components 
                   xvar    = None  , ## x-variable 
                   yvar    = None  , ## y-variable 
                   zvar    = None  , ## z-variable                   
                   name    = ''    ) : 
        
        ## keep all arguments 
        self.__args = {
            #
            'signal_x'   : signal_x   ,
            'signal_y'   : signal_y   ,
            'signal_z'   : signal_z   ,
            ##
            'bkg_1x'     : bkg_1x     ,
            'bkg_1y'     : bkg_1y     ,
            ##
            'bkg_2x'     : bkg_2x     ,
            'bkg_2y'     : bkg_2y     ,
            'bkg_2xy'    : bkg_2xy    ,
            'bkg_2yz'    : bkg_2yz    ,
            ##
            'bkg_3x'     : bkg_3x     ,
            'bkg_3y'     : bkg_3y     ,
            'bkg_3yz'    : bkg_3yz    ,
            ##
            'bkg_3D'     : bkg_3D     ,
            #
            'sss'        : sss        ,            
            'bss'        : bss        ,
            'ssb'        : ssb        ,
            'sbb'        : sbb        ,
            'bbs'        : bbs        ,
            'bbb'        : bbb        ,
            #
            'components' : components ,
            'xvar'       : xvar       ,
            'yvar'       : yvar       ,
            'zvar'       : zvar       ,
            'name'       : name                   
            }
        
        self.__suffix      = suffix

        self.__signal_x , xvar = self.make_PDF1 ( signal_x , xvar = xvar , prerfix = 'SX' , suffix = suffix )
        self.__signal_y , yvar = self.make_PDF1 ( signal_y , xvar = yvar , prerfix = 'SY' , suffix = suffix )

        if zvar and not signal_z :
            self.__signal_z = self.__signal_y.clone ( xvar = zvar , name_prefix = 'SZ_' )
            self.debug('signal z-component is cloned from the signal_y component')
        else : 
            self.__signal_z , zvar = self.make_PDF1 ( signal_z , xvar = zvar , prerfix = 'SZ' , suffix = suffix )
            
        #
        ## initialize base class
        #
        if not name :
            name = "%s&%s&%s" % ( self.__signal_x.name ,
                                  self.__signal_y.name ,
                                  self.__signal_z.name )                                  
            if suffix : name += '_'+ suffix
            
        PDF3.__init__ ( self , name          ,
                        self.__signal_x.xvar ,
                        self.__signal_y.xvar ,
                        self.__signal_z.xvar ) 

        # =====================================================================
        ## 1) All signals component 
        # =====================================================================
        self.__sss_cmp  = Model3D (
            self.new_name ( 'SSS' , suffix ) ,
            self.__signal_x , self.__signal_y , self.__signal_z )
        
        # =====================================================================
        ## 2) background x signal x  signal 
        # =====================================================================        
        self.__bkg_1x  = self.make_bkg (  bkg_1x , self.new_name ( 'Bkg1X_BSS' , suffix ) , self.xvar )
        self.__bss_cmp = Model3D ( self.new_name ( "BSS" , suffix ) ,
                                    self.__bkg_1x , self.__signal_y , self.__signal_z ,
                                    title = "Background1(x) x Signal(y) x Signal(z)" )
        
        # =====================================================================
        ## 3) signal x (  signal x background + backround x signal ) 
        # =====================================================================

        self.__bkg_1y  = self.make_bkg (        bkg_1y , self.new_name ( 'Bkg1Y_SBS' , suffix ) , self.yvar )
        self.__bkg_1z  = self.make_bkg ( self.__bkg_1y , self.new_name ( 'Bkg1Z_SSB' , suffix ) , self.zvar )
        
        self.__ssb_cmp_raw = Model3D ( self.new_name ( "SSB_raw" , suffix ) ,
                                       self.__signal_x , self.__signal_y , self.__bkg_1z   ,
                                       title = "Signal(x) x Signal(y) x Background1(z)" )
        self.__sbs_cmp_raw = Model3D ( self.new_name ( "SBS_raw" , suffix ) ,
                                       self.__signal_x , self.__bkg_1y   , self.__signal_z , 
                                       title = "Signal(x) x Background1(y) x Signal(z)" ) 

        self.__ssb_sym_cmp  = Generic3D_pdf (
            self.make_sum ( self.new_name ( "SSB" , suffix ) ,
                            "S(x)*S(y)*B(z)+S(x)*B(y)*S(z)" ,
                            self.__ssb_cmp_raw.pdf   ,
                            self.__sbs_cmp_raw.pdf ) ,
            self.xvar , self.yvar , self.zvar )
        
        self.__ssb_cmp = self.__ssb_sym_cmp  ## 
        self.__sbs_cmp = self.__ssb_sym_cmp  ## ditto
        

        # =====================================================================
        ## (intermezzo-1) Assumptions about the SBB-background sub-components 
        # =====================================================================        
        
        if   component_clone   ( bkg_2x ) :
            bkg_2x = self.__bkg_1x
            self.debug ( 'bkg_2x set to [CLONE]   %s' % bkg_2x )
        elif component_similar ( bkg_2x ) :
            bkg_2x =        bkg_1x
            self.debug ( 'bkg_2x set to [SIMILAR] %s' % bkg_2x ) 

        if   component_clone   ( bkg_2y ) :
            bkg_2y = self.__bkg_1y
            self.debug ( 'bkg_2y set to [CLONE]   %s' % bkg_2y ) 
        elif component_similar ( bkg_2y ) :
            bkg_2y =        bkg_1y
            self.debug ( 'bkg_2y set to [SIMILAR] %s' % bkg_2y )
        # =====================================================================


        from ostap.fitting.pdf_ops import Prod3D_pdf 
            
        # =====================================================================
        ## 4) signal x background x background 
        # =====================================================================
        
        self.__bkg_2x = None 
        self.__bkg_2y = None 
        self.__bkg_2z = None 


        bkg_2yz_name = self.new_name ( 'Bkg2YZ' , suffix )
        if bkg_2yz and isinstance ( bkg_2yz , (  tuple , list ) ) :
            from ostap.fitting.models_2d import make_B2Dsym
            self.__bkg_2yz = make_B2Dsym ( bkg_2yz_name , self.yvar , self.zvar , *bkg_2yz  )            
        elif bkg_2yz:
            self.__bkg_2yz = self.make_PDF2 ( bkg_2yz ,
                                              xvar   = self.yvar ,
                                              yvar   = self.zvar ,
                                              prefix = 'Bkg2YZ'  ,
                                              suffix = suffix    ) [ 0 ] 
        else :
            self.__bkg_2y  = self.make_bkg (        bkg_2y , self.new_name ( 'Bkg2Y_S2B' , suffix ) , self.yvar )        
            self.__bkg_2z  = self.make_bkg ( self.__bkg_2y , self.new_name ( 'Bkg2Z_S2B' , suffix ) , self.zvar )            
            self.__bkg_2yz = Model2D ( bkg_2yz_name   ,
                                       self.__bkg_2y  ,
                                       self.__bkg_2z  , title =  'Background2(y,z)' )

        ## make component

            
        self.__sbb_cmp = Prod3D_pdf ( ( self.__signal_x.pdf , self.__bkg_2yz )  ,
                                      xvar = self.xvar , yvar = self.yvar , zvar = self.zvar ,
                                      name = self.new_name ( 'SBB' ,  "S(x)*(y,z)" ) )
            
        # =====================================================================
        ## 5) background x ( signal x background + background x  signal ) 
        # =====================================================================

        bkg_2xz_name = self.new_name ( 'Bkg2XZ' , suffix  )
        
        if   bkg_2xy and isinstance ( bkg_2xy , PDF2 ) :            
            self.__bkg_2xy = bkg_2xy
            self.__bkg_2xz = bkg_2xy.clone ( name = self.new_name ( 'Bkg2XZ' ,  suffix ) ,
                                             xvar = self.xvar             ,
                                             yvar = self.zvar             )            

        elif bkg_2xy and isinstance ( bkg_2xy ,  ( tuple , list )  ) :
            
            from ostap.fitting.models_2d import make_B2D
            self.__bkg_2xy = make_B2D ( self.new_name ( 'Bkg2XY' , suffix  ) , self.xvar , self.yvar , *bkg_2xy  )
            self.__bkg_2xz = self.__bkg_2xy.clone ( name = self.new_name ( 'Bkg2XZ' , suffix ) ,
                                                    xvar = self.xvar         ,
                                                    yvar = self.zvar         )          
        else :

            if not self.__bkg_2x : self.__bkg_2x = self.make_bkg (        bkg_2x , self.new_name ( 'Bkg2X_S2B' , suffix ) , self.xvar )        
            if not self.__bkg_2y : self.__bkg_2y = self.make_bkg (        bkg_2y , self.new_name ( 'Bkg2Y_S2B' , suffix ) , self.yvar )        
            if not self.__bkg_2z : self.__bkg_2z = self.make_bkg ( self.__bkg_2y , self.new_name ( 'Bkg2Z_S2B' , suffix ) , self.zvar )
            
            self.__bkg_2xy = Model2D ( self.new_name ( 'Bkg2XY' , suffix ) , 
                                       self.__bkg_2x         ,
                                       self.__bkg_2y         , title = 'Background2(x,y)' )
            self.__bkg_2xz = Model2D ( self.new_name ( 'Bkg2XZ' , suffix ) ,
                                       self.__bkg_2x         ,
                                       self.__bkg_2z         , title = 'Background2(x,z)' )

        ## make components
            
        self.__bsb_cmp_raw = Prod3D_pdf ( ( self.__signal_y , self.__bkg_2xz ) ,
                                          xvar = self.xvar , yvar = self.yvar , zvar = self.zvar ,
                                          name = self.new_name ( 'BSB_raw' , "S(y)*(x,z)" ) )

        self.__bbs_cmp_raw = Prod3D_pdf ( ( self.__signal_z , self.__bkg_2xy ) ,
                                          xvar = self.xvar , yvar = self.yvar , zvar = self.zvar ,
                                          name = self.new_name ( 'BBS_raw' , "S(z)*(x,y)" ) )

        self.__bbs_sym_cmp = Generic3D_pdf ( self.make_sum ( self.new_name ( "SBB" , suffix ) ,
                                                             "S(x)*B(y,z)+S(y)*B(x,z)+S(Z)*B(x,y)" ,
                                                             self.__bsb_cmp_raw.pdf ,
                                                             self.__bbs_cmp_raw.pdf ) ,
                                                 self.xvar , self.yvar , self.zvar )

        self.__bbs_cmp = self.__bbs_sym_cmp ## ditto
        self.__bsb_cmp = self.__bbs_sym_cmp
        
        # =====================================================================
        ## (intermezzo-2) Assumptions about the BBB-background sub-components 
        # =====================================================================        
        if   component_clone   ( bkg_3x ) :
            bkg_3x = self.__bkg_2x
            self.debug ( 'bkg_3x set to [CLONE]   %s' % bkg_3x ) 
        elif component_similar ( bkg_3x ) :
            bkg_3x =        bkg_2x
            self.debug ( 'bkg_3x set to [SIMILAR] %s' % bkg_3x ) 

        if   component_clone   ( bkg_3y ) :
            bkg_3y = self.__bkg_2y
            self.debug ( 'bkg_3y set to [CLONE]   %s' % bkg_3y ) 
        elif component_similar ( bkg_3y ) :
            bkg_3y =        bkg_2y
            self.debug ( 'bkg_3y set to [SIMILAR] %s' % bkg_3y ) 

        if   component_clone   ( bkg_3yz ) :
            bkg_3yz = self.__bkg_2yz
            self.debug ( 'bkg_3yz set to [CLONE]   %s' % bkg_3yz ) 
        elif component_similar ( bkg_3yz ) :
            bkg_3yz =        bkg_2yz
            self.debug ( 'bkg_3yz set to [SIMILAR] %s' % bkg_3yz ) 
        # =====================================================================

        # =====================================================================
        ## 6) pure background 
        # =====================================================================
        
        self.__bkg_3x  = None 
        self.__bkg_3y  = None 
        self.__bkg_3z  = None
        
        self.__bkg_3yz = None 

        bbb_name = self.new_name ( 'BBB' , suffix ) 
        if bkg_3D and isinstance ( bkg_3D ,  ( tuple , list )  ) :
            from ostap.fitting.models_2d import make_B2DmixYZ 
            self.__bbb_cmp = make_B2DmixYZ ( bbb_name , self.xvar , self.yvar , self.zvar , *bkg_3D )
        elif bkg_3D :
            self.__bbb_cmp = self.make_PDF3 ( bkg_3D ,
                                              xvar   = self.xvar ,
                                              yvar   = self.yvar ,
                                              zvar   = self.zvar ,
                                              prefix = 'BBB'     ,
                                              suffix = suffix    ) [ 0 ] 
        else :

            if   bkg_3yz and isinstance ( bkg_3yz , PDF2 ) :
                self.__bkg_3xy = bkg_3xy
            elif bkg_3yz and isinstance ( bkg_3yz , ROOT.RooAbsPdf ) :
                self.__bkg_3yz = Generic2D_pdf ( bkg_3yz , self.yvar , self.zvar , self.new_name ( 'Bkg3YZ' , suffix ) )
            else :
                self.__bkg_3y  = self.make_bkg (        bkg_3y , self.new_name ( 'Bkg3Y_BBB' , suffix ) , self.yvar )        
                self.__bkg_3z  = self.make_bkg ( self.__bkg_3y , self.new_name ( 'Bkg3Z_BBB' , suffix ), self.zvar )
                self.__bkg_3yz = Model2D ( self.new_name ( 'Bkg3YZ' , suffix ) , self.__bkg_3y , self.__bkg_3z )
                
            self.__bkg_3x  = self.make_bkg ( bkg_3x , self.new_name ( 'Bkg3X_BBB' , suffix ) , self.xvar )

            ## make component
            
            self.__bbb_cmp = Prod3D_pdf ( ( self.__bkg_3x , self.__bkg_3yz ) , 
                                          xvar = self.xvar , yvar = self.yvar , zvar = self.zvar ,
                                          name = self.new_name ( 'BBB' ,  "B(x)*B(y,z)" ) )
            
        # =====================================================================
        ## coefficients
        # =====================================================================
        
        self.__sss = self.make_var ( sss   , "SSS" + suffix ,
                                     "Signal(x)&Signal(y)&Signal(z)"     + suffix , None , 1000 , 0 , 1.e+7 )
        self.__bss = self.make_var ( bss   , "BSS" + suffix ,
                                     "Background(x)&Signal(y)&Signal(z)" + suffix , None , 1000 , 0 , 1.e+7 )
        self.__ssb = self.make_var ( ssb   , "SSB" + suffix ,
                                     "Signal&(SB+BS)"                    + suffix , None , 1000 , 0 , 1.e+7 )
        self.__sbb = self.make_var ( sbb   , "SBB" + suffix ,
                                     "Signal&Background(x,y)"            + suffix , None , 1000 , 0 , 1.e+7 )
        self.__bbs = self.make_var ( bbs   , "BBS" + suffix ,
                                     "Backgroun&(SB+BS)"                 + suffix , None , 1000 , 0 , 1.e+7 )        
        self.__bbb = self.make_var ( bbb  , "BBB"  + suffix ,
                                     "Background(x,y,z)"                 + suffix , None , 1000 , 0 , 1.e+7 )

        self.__sbs = self.__ssb ## the same 
        self.__bsb = self.__bbs ## the same 
        
        self.alist1 = ROOT.RooArgList (
            self.__sss_cmp.pdf ,
            self.__bss_cmp.pdf ,
            self.__ssb_cmp.pdf ,
            self.__sbb_cmp.pdf ,
            self.__bbs_cmp.pdf ,
            self.__bbb_cmp.pdf )
        
        self.alist2 = ROOT.RooArgList (
            self.__sss     ,
            self.__bss     ,
            self.__ssb     ,
            self.__sbb     ,
            self.__bbs     ,
            self.__bbb     )
        
        ## treat additional components (if specified)
        self.__nums_components = [] 
        icmp = 0
        self.__more_components = [
            self.make_PDF3 ( cmp ,
                             xvar   = self.xvar ,
                             yvar   = self.yvar ,
                             zvar   = self.zvar , 
                             prefix = 'C%d' % i ,
                             suffix = suffix    ) [ 0 ] for ( i , cmp ) in enumerate ( components ) ]  
        for cmp in self.__more_components :
            self.components.add ( cmp.pdf ) 

        nc = len( self.__more_components )
        if 1 == nc :
            cf = self.make_var ( None , "C"+suffix , "Component" + suffix , None , 1 , 0 , 1.e+7 )
            self.alist1.add  ( self.components[0] )
            self.__nums_components.append ( cf ) 
        elif 2 <= nc : 
            fic = self.make_fracs ( nc , 'C_%%d%s' % suffix ,  'C(%%d)%s'  % suffix , fractions  = False )
            for c in self.components : self.alist1.add ( c)
            for f in fic             : self.__nums_components.append ( f )
            
        self.__nums_components  = tuple ( self.__nums_components  ) 
        for c in self.__nums_components  : self.alist2.add ( c )

        ## 
        #
        ## build the final PDF 
        # 
        pdfname  = self.new_roo_name ( "fit3dm" , suffix  ) 
        pdftitle = "Fit3DMix %s " % self.name
        pdfargs  = pdfname , pdftitle , self.alist1 , self.alist2
        self.pdf = ROOT.RooAddPdf  ( *pdfargs )

        self.signals     .add ( self.__sss_cmp.pdf )
        self.backgrounds .add ( self.__bbb_cmp.pdf )
        self.crossterms1 .add ( self.__ssb_cmp.pdf ) ## cross-terms
        self.crossterms2 .add ( self.__sbb_cmp.pdf ) ## cross-terms 
        
        ## save the configuration
        self.config = {
            ##
            'signal_x'   : self.signal_x ,
            'signal_y'   : self.signal_y ,
            'signal_z'   : self.signal_z ,
            ##
            'suffix'     : self.suffix  ,
            ## SSB-terms 
            'bkg_1x'     : self.bkg_1x  ,
            'bkg_1y'     : self.bkg_1y  ,
            ## SBB-terms 
            'bkg_2x'     : self.bkg_2x  ,
            'bkg_2y'     : self.bkg_2y  ,
            'bkg_2xy'    : self.bkg_2xy ,
            'bkg_2yz'    : self.bkg_2yz ,
            ## BBB-term
            'bkg_3x'     : self.bkg_3x  ,
            'bkg_3y'     : self.bkg_3y  ,
            ## 'bkg_3xz'    : self.bkg_3xz ,
            'bkg_3yz'    : self.bkg_3yz ,
            'bkg_3D'     : self.bkg_3D  ,
            ## Yields 
            'sss'        : self.SSS     ,
            'bss'        : self.BSS     ,
            'ssb'        : self.SSB     ,
            'sbb'        : self.SBB     ,
            'bbs'        : self.BBS     ,
            'bbb'        : self.BBB     ,            
            #
            'components' : self.more_components ,            
            'xvar'       : self.xvar    ,
            'yvar'       : self.yvar    ,
            'zvar'       : self.zvar    ,
            'name'       : self.name    ,             
            }
        
        self.checked_keys.add  ( 'xvar' )
        self.checked_keys.add  ( 'yvar' )
        self.checked_keys.add  ( 'zvar' )

    @property
    def SSS ( self ) :
        """The yield of Signal(x)*Signal(y)*Signal(z) component"""
        return self.__sss
    @SSS.setter 
    def SSS ( self , value ) :
        value = float ( value  )
        assert value in self.__sss, "Value %s is out of the allowed range %s " % ( value , self.__sss.minmax() )
        self.__sss.setVal ( value ) 

    @property
    def SSB ( self ) :
        """The yield of S(x)*[S(y)*(z)+B(y)*S(z)] symmetrized component (same as SBS)"""
        return self.__bss
    @SSB.setter 
    def SSB ( self , value ) :
        value = float ( value  )
        assert value in self.__bss, "Value %s is out of the allowed range %s " % ( value , self.__bss.minmax() )
        self.__bss.setVal ( value ) 

    @property
    def SBS ( self ) :
        """The yield of S(x)*[S(y)*B(z)+B(y)*S(z)] symmetrized component (same as SSB)"""
        return self.__sbs
    @SBS.setter 
    def SBS ( self , value ) :
        value = float ( value  )
        assert value in self.__sbs, "Value %s is out of the allowed range %s " % ( value , self.__sbs.minmax() )
        self.__sbs.setVal ( value ) 

    @property
    def BSS ( self ) :
        """The yield of B(x)*S(y)*S(z) component"""
        return self.__bss
    @BSS.setter 
    def BSS ( self , value ) :
        value = float ( value  )
        assert value in self.__bss, "Value %s is out of the allowed range %s " % ( value , self.__bss.minmax() )
        self.__bss.setVal ( value ) 

    @property
    def SBB ( self ) :
        """The yield of S(x)*B(y,z) component"""
        return self.__sbb
    @SBB.setter 
    def SBB ( self , value ) :
        value = float ( value  )
        assert value in self.__sbb, "Value %s is out of the allowed range %s " % ( value , self.__sbb.minmax() )
        self.__sbb.setVal ( value ) 

    @property
    def BSB ( self ) :
        """The yield of B(x)*[S(y)*B(z)+B(y)*S(z)] symmetrized component (same as BBS)"""
        return self.__bsb
    @BSB.setter 
    def BSB ( self , value ) :
        value = float ( value  )
        assert value in self.__bsb, "Value %s is out of the allowed range %s " % ( value , self.__bsb.minmax() )
        self.__bsb.setVal ( value ) 

    @property
    def BBS ( self ) :
        """The yield of B(x)*[S(y)*B(z)+B(y)*S(z)] symmetrized component (same as BSB)"""
        return self.__bbs
    @BBS.setter 
    def BBS ( self , value ) :
        value = float ( value  )
        assert value in self.__bbs, "Value %s is out of the allowed range %s " % ( value , self.__bbs.minmax() )
        self.__bbs.setVal ( value ) 

    @property
    def BBB ( self ) :
        """The yield of B(x,y,z) component"""
        return self.__bbb
    @BBB.setter 
    def BBB ( self , value ) :
        value = float ( value  )
        assert value in self.__bbb, "Value %s is out of the allowed range %s " % ( value , self.__bbb.minmax() )
        self.__bbb.setVal ( value ) 

    @property
    def C ( self ) :
        """Get the  yields of 'other' component(s) 
        For single 'other' component:
        >>> print pdf.C           ## read the single 'other' component 
        >>> pdf.C = 100           ## assign to it 
        For multiple 'other' components:
        >>> print pdf.C[4]        ## read the 4th 'other' component 
        >>> pdf.C = 4,100         ## assign to it 
        ... or, alternatively:
        >>> print pdf.C[4]        ## read the 4th 'other' component 
        >>> pdf.C[4].value 100    ## assign to it         
        """
        return self.component_getter ( self.__nums_components  )     
    @C.setter
    def C (  self , value ) :
        self.component_setter ( self.__nums_components , value )
   
    @property
    def yields    ( self ) :
        """The list/tuple of the yields of all numeric components"""
        return tuple ( [ i for i in  self.alist2 ] )

    @property 
    def total_yield ( self ) :
        """'total_yield'' : get the total yield"""
        if not self.fit_result                   : return None
        if not valid_pointer ( self.fit_result ) : return None
        return self.fit_result.sum ( *self.yields ) 
    
    # =========================================================================
    # components
    # =========================================================================
    
    @property 
    def signal_x ( self  ) :
        """'signal_x' : Signal(x) component/PDF"""
        return self.__signal_x
    @property 
    def signal_y ( self  ) :
        """'signal_y' : Signal(y) component/PDF"""
        return self.__signal_y
    @property 
    def signal_z ( self  ) :
        """'signal_z' : Signal(z) component/PDF"""
        return self.__signal_z

    @property
    def bkg_1x ( self ) :
        """'bkg_1x' : B(x) component for B(x)*S(y)*S(z) term"""
        return self.__bkg_1x
    @property
    def bkg_1y ( self ) :
        """'bkg_1y' : B(y) component for S(x)*B(y)*S(z) term"""
        return self.__bkg_1y
    @property
    def bkg_1z ( self ) :
        """'bkg_1z' : B(z) component for S(x)*S(y)*B(z) term"""
        return self.__bkg_1z
    
    @property
    def bkg_2xy( self ) :
        """'bkg_2xy' : B(x,y) component for B(x,y)*S(z) term"""
        return self.__bkg_2xy
    @property
    def bkg_2xz( self ) :
        """'bkg_2xz' : B(x,z) component for B(x,z)*S(y) term"""
        return self.__bkg_2xz
    @property
    def bkg_2yz( self ) :
        """'bkg_2yz' : B(y,z) component for B(y,z)*S(x) term"""
        return self.__bkg_2yz 

    @property
    def bkg_2x ( self ) :
        """'bkg_2x' : B(x) component for B(x,y)*S(z) & B(x,z)*S(y) terms"""
        return self.__bkg_2x
    @property
    def bkg_2y ( self ) :
        """'bkg_2y' : B(y) component for B(y,z)*S(x) & B(x,y)*S(z) terms"""
        return self.__bkg_2y
    @property
    def bkg_2z ( self ) :
        """'bkg_2z' : B(z) component for B(x,z)*S(y) & B(y,z)*S(x) terms"""
        return self.__bkg_2z

    @property
    def bkg_3x ( self ) :
        """'bkg_3x' : B(x) component for B(x)*B(y)*B(z) term"""
        return self.__bkg_3x
    @property
    def bkg_3y ( self ) :
        """'bkg_3y' : B(y) component for B(x)*B(y)*B(z) term"""
        return self.__bkg_3y
    @property
    def bkg_3z ( self ) :
        """'bkg_3z' : B(z) component for B(x)*B(y)*B(z) term"""
        return self.__bkg_3z

    @property
    def bkg_3yz ( self ) :
        """'bkg_3yz' : B(y,z) component for B(x)*B(y,z) term"""
        return self.__bkg_3yz

    @property
    def bkg_3D ( self ) :
        """B(x,y,z) component/PDF for the final PDF"""
        return self.__bbb_cmp 

    @property
    def cmp_SSS ( self ) :
        """'triple-signal' component/PDF"""
        return self.__sss_cmp

    @property
    def cmp_SSB ( self ) :
        """'signal-signal-background' symmetrized component/PDF"""
        return self.__ssb_cmp
    
    @property
    def cmp_SBS ( self ) :
        """'signal-background-signal' symmetrized component/PDF (same as cmp_SSB)"""
        return self.__sbs_cmp
    
    @property
    def cmp_BSS ( self ) :
        """'background-signal-signal' symmetrized component/PDF (same as cmp_SSB)"""
        return self.__bss_cmp

    @property
    def cpm_SBB ( self ) :
        """'signal-background-background'  symmetrized component/PDF"""
        return self.__sbb_cmp
    
    @property
    def cmp_BSB ( self ) :
        """'background-signal-background'  symmetrized component/PDF (same as cmp_SBB)"""
        return self.__bsb_cmp
    
    @property
    def cmp_BBS ( self ) :
        """'background-background-signal' symmetrized component/PDF (same as cmp_SBB)"""
        return self.__bbs_cmp

    @property
    def cmp_BBB ( self ) :
        """'triple-background' component/PDF"""
        return self.__bbb_cmp

    @property
    def more_components ( self ) :
        """additional/'other' components"""
        return tuple( self.__more_components  )

    @property
    def suffix ( self ) :
        """'suffix', used to build the name"""
        return self.__suffix

    # =========================================================================
    ## Raw, non-symmetrized fit components/PDF (for debugging)
    # =========================================================================

    @property
    def cpm_raw_SSB ( self ) :
        """'signal-signal-background' raw,non-symmetrized component/PDF"""
        return self.__ssb_cmp_raw
    
    @property
    def cmp_raw_SBS ( self ) :
        """'signal-background-signal'' raw,non-symmetrized component/PDF"""
        return self.__sbs_cmp_raw
    
    @property
    def cmp_raw_BBS ( self ) :
        """'background-background-signal' raw,non-symmetrized component/PDF"""
        return self.__bss_cmp_raw

    @property
    def cmp_raw_BSB ( self ) :
        """'background-signal-background' raw,non-symmetrized component/PDF"""
        return self.__bsb_cmp_raw


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
