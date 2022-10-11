#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/fit2d.py
#  Set of useful basic utilities to build various fit models 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
# =============================================================================
"""Set of useful basic utilities to build various 2D-fit models"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'Flat2D'        , ## simplest 2D-model: constant function
    'Model2D'       , ## helper class to construct 2D-models. 
    'Sum2D'         , ## non-extended sum of two PDFs
    'H2D_pdf'       , ## convertor of 1D-histo to RooDataPdf
    'Shape2D_pdf'   , ## simple PDF from C++ shape     
    ## 
    'Fit2D'         , ## the model for 2D-fit: signal + background + optional components
    'Fit2DSym'      , ## the model for 2D-fit: signal + background + optional components
    ##
    )
# =============================================================================
from   builtins                 import range
from   ostap.core.core          import Ostap , valid_pointer, roo_silent 
from   ostap.core.ostap_types   import integer_types, num_types, list_types, iterable_types   
from   ostap.fitting.funbasic   import FUN2
from   ostap.fitting.pdfbasic   import PDF2, Generic2D_pdf 
from   ostap.fitting.utils      import component_similar , component_clone 
from   ostap.fitting.fithelpers import H2D_dset, Fractions  
import ROOT, random
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.fit2d' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
## @class Model2D
#  Trivial class to construct 2D model as a product of split 1D-models
#  actually it is a tiny  wrapper over <code>ROOT.RooProdPdf</code>
#  @code
#  pdfx = ...
#  pdfy = ...
#  pdf2D = Model2D( 'D2' , xmodel = pdfx , ymodel =  pdfy )
#  @endcode 
class Model2D(PDF2) :
    """Trivial class to construct 2D model as a product of split 1D-models
    - actually it is a tiny  wrapper over ROOT.RooProdPdf
    >>> pdfx = ...
    >>> pdfy = ...
    >>> pdf2D = Model2D( 'D2' , xmodel = pdfx , ymodel =  pdfy )
    """
    def __init__ ( self         ,
                   name         ,
                   xmodel       ,
                   ymodel       ,
                   xvar  = None ,
                   yvar  = None ,
                   title = ''   ) :

        if xvar and not xmodel : xmodel = Flat1D ( xvar )
        if yvar and not ymodel : ymodel = Flat1D ( yvar )

        self.__xmodel , xvar = self.make_PDF1 ( xmodel , xvar = xvar , prefix = 'X' )
        self.__ymodel , yvar = self.make_PDF1 ( ymodel , xvar = yvar , prefix = 'Y' )

        name  = name  if name  else self.generate_name ( '(%s)*(%s)'  % ( self.xmodel.name , self.ymodel.name ) )
        
        ## initialize the base 
        PDF2.__init__ ( self , name , self.xmodel.xvar , self.ymodel.xvar ) 

        ## make pdf
        from ostap.fitting.pdf_ops import raw_product 
        self.pdf = raw_product ( self , self.xmodel , self.ymodel )
            
        ## save configuration 
        self.config = {
            'name'   :  self.name   ,
            'xmodel' :  self.xmodel ,
            'ymodel' :  self.ymodel ,
            'xvar'   :  self.xvar   ,
            'yvar'   :  self.yvar   ,            
            }
        
    ## redefine the clone 
    def clone ( self , **kwargs ) :
        """ Redefine the clone
        """

        ## xm     = kwargs.pop ( 'xmodel' , None )
        ## assert ( not xm ) or ( xm is self.xmodel ) , "Model2D.clone: invalid usage of 'xmodel'!"            
        ## ym     = kwargs.pop ( 'ymodel' , None )
        ## assert ( not ym ) or ( ym is self.ymodel ) , "Model2D.clone: invalid usage of 'ymodel'!"
        
        name   = kwargs.pop ( 'name' , self.name )        
        xvar   = kwargs.pop ( 'xvar' , self.xvar )
        yvar   = kwargs.pop ( 'yvar' , self.yvar )

        ## transpose 
        if  ( xvar is self.yvar ) and ( yvar is self.xvar ) :
            
            xmodel = self.ymodel 
            ymodel = self.xmodel
            
        else :

            if xvar is self.xvar : xmodel = self.xmodel
            else                 : xmodel = self.xmodel.clone ( xvar = xvar , **kwargs )

            if yvar is self.yvar : ymodel = self.ymodel
            else                 : ymodel = self.ymodel.clone ( xvar = yvar , **kwargs )
            
        return PDF2.clone ( self ,
                            name   = name   , 
                            xmodel = xmodel ,
                            ymodel = ymodel ,
                            xvar   = xvar   ,
                            yvar   = yvar   , **kwargs )
    
    @property
    def xmodel ( self ) :
        """'x-model'' x-component of Model(x)*Model(y) PDF"""
        return self.__xmodel

    @property
    def ymodel ( self ) :
        """'y-model'' y-component of Model(x)*Model(y) PDF"""
        return self.__ymodel


# =============================================================================
## @class Sum2D
#  Non-extended sum of several PDFs
#  It is just a small wrapper for <code>ROOT.RooAddPdf</code>
#  @see RooAddPdf 
class Sum2D (PDF2,Fractions) :
    """Non-extended sum of several PDFs:
    
    It is just a small wrapper for <code>ROOT.RooAddPdf</code>
    - see RooAddPdf 
    
    >>> sum  = Sum2D ( [ pdf1 , pdf2 , pdf3 ]  ) 

    """
    def __init__ ( self             ,
                   pdfs             , ## input list of PDFs  
                   xvar      = None , 
                   yvar      = None , 
                   name      = ''   ,
                   recursive = True ,
                   prefix    = 'f'  , ## prefix for fraction names 
                   suffix    = ''   , ## suffix for fraction names 
                   fractions = None ) :

        assert 2 <= len ( pdfs ) , 'Sum2D: at least two PDFs are needed!'
        
        pdf_list = []           
        for i , p in enumerate ( pdfs ) :
            cmp , xvar , yvar = self.make_PDF2 ( p , xvar = xvar , yvar = yvar ) 
            pdf_list.append ( cmp )

        ## generic name 
        patname =  '+'.join ( '(%s)' % p.name for p in pdf_list )

        ## check the instance name 
        name    = name if name else self.new_name ( patname ) 
        
        ## ininialize the base classes 
        PDF2.     __init__ ( self , name , xvar , yvar )        
        Fractions.__init__ ( self , pdf_list       , 
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
            'name'      : self.name      , 
            'prefix'    : self.prefix    ,
            'suffix'    : self.suffix    ,
            'fractions' : self.fractions ,
            'recursive' : self.recursive        
            }
  
# =============================================================================
## @class Flat2D
#  The most trivial 2D-model - constant
#  @code 
#  pdf = Flat2D( 'flat' , xvar = ...  , yvar = ... )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
class Flat2D(PDF2) :
    """The most trival 2D-model - constant
    >>> pdf = Flat2D( 'flat' , xvar = ...  , yvar = ... )
    """
    def __init__ ( self , xvar , yvar , name = '' ,  title = '' ) :

        name = name if name else self.generate_name ( prefix = 'flat2D_')                            
        PDF2.__init__ ( self  , name , xvar , yvar ) 
                        
        if not title : title = 'flat2(%s)' % name 
        self.pdf = Ostap.Models.Uniform ( name , title , self.xvar , self.yvar )
        assert 2 == self.pdf.dim() , 'Flat2D: wrong dimensionality!'
        
        ## save configuration
        self.config = {
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'name'     : self.name ,            
            'title'    : title     ,             
            }                   

# =============================================================================
## Generic 2D-shape from C++ callable
#  @see Ostap::Models:Shape2D
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2020-07-20
class Shape2D_pdf(PDF2) :
    """ Generic 2D-shape from C++ callable
    - see Ostap::Models:Shape2D
    """
    
    def __init__ ( self , name , shape , xvar , yvar ) :

        if isinstance ( shape , ROOT.TH2 ) and not isinstance ( shape , ROOT.TH3 ) and not xvar :
            xvar = shape.xminmax()
            
        if isinstance ( shape , ROOT.TH2 ) and not isinstance ( shape , ROOT.TH3 ) and not yvar :
            yvar = shape.yminmax()
            
        if isinstance ( shape , ROOT.TH2 ) and not isinstance ( shape , ROOT.TH3 ) :
            self.histo = shape
            shape      = Ostap.Math.Histo2D ( shape )

        ##  iniialize the base 
        PDF2.__init__ ( self , name , xvar , yvar ) 
            
        self.__shape = shape
        
        if isinstance ( self.shape , Ostap.Math.Histo2D ) :
            
            ## create the actual pdf
            self.pdf = Ostap.Models.Histo2D ( self.roo_name ( 'histo2_' ) , 
                                              "Histo-2D %s" % self.name   ,
                                              self.xvar                   ,
                                              self.yvar                   ,
                                              self.shape                  )            
        else : 
            
            ## create the actual pdf
            self.pdf = Ostap.Models.Shape2D.create  (
                self.roo_name  ( 'shape2_' ) , 
                "Shape-2D %s" % self.name    ,
                self.xvar                    ,
                self.yvar                    ,
                self.shape                   )  
            
        ## save the configuration
        self.config = {
            'name'    : self.name    , 
            'shape'   : self.shape   , 
            'xvar'    : self.xvar    , 
            'yvar'    : self.yvar    , 
            }
        
    @property
    def shape  ( self ) :
        """'shape' : the actual C++ callable shape"""
        return self.__shape 
 
# ===================================================it==========================
## simple convertor of 2D-histogram into PDF
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H2D_pdf(PDF2) :
    """Simple convertor of 2D-histogram into PDF 
    """
    def __init__ ( self            ,
                   name            ,
                   histo           ,
                   xvar    = None  , 
                   yvar    = None  ,
                   density = False ,
                   order   = 0     , ## interpolation order 
                   silent  = False ) :

        assert isinstance ( order, integer_types ) and 0 <= order ,\
               'Invalid interpolation order: %s/%s' % ( order , type ( order ) )

        self.__ds = H2D_dset ( histo , xvar = xvar , yvar = yvar , density = density ,  silent = silent )
        
        PDF2    .__init__ ( self , name  , self.ds.xaxis , self.ds.yaxis ) 
        
        self.__vset  = ROOT.RooArgSet  ( self.xvar , self.yvar )

        #
        ## finally create PDF :
        #
        with roo_silent ( silent ) : 
            self.pdf    = ROOT.RooHistPdf (
                self.roo_name  ( 'histo2_' ) , 
                'Histo-2D PDF: %s/%s' % ( self.histo.GetName() , self.histo.GetTitle() ) , 
                self.__vset , 
                self.dset   ,
                order       )
            
        ## and declare it be be a "signal"
        self.signals.add ( self.pdf ) 

        ## save the configuration
        self.config = {
            'name'    : self.name    , 
            'histo'   : self.histo   , 
            'xvar'    : self.xvar    , 
            'yvar'    : self.yvar    , 
            'density' : self.density , 
            'silent'  : self.silent  ,             
            'order'   : self.order   ,             
            }
        
    @property
    def ds ( self ) :
        """'ds' : the H2D_dset object"""
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
# Compound 2D-fit models 
# =============================================================================

# =============================================================================
## @class Fit2D
#  The actual model for 2D-fits. It consists of four main components :
#    - pure signal :        S(x)*S(y)
#    - signal x background: S(x)*B(y)
#    - signal x bakcgronud: B(x)*S(y)
#    - pure backrground:    B(x,y) or B(x)*B(y) 
#  Other 2D-components could be specified in addition
#  @param  signal_x  PDF for the S(x)-signal component
#  @param  signal_y  PDF for the S(y)-signal component
#  @param  suffix    suffix to be used for the PDF and variable names
#  @param  bkg_1x    x-background component for B(x)*S(y) term 
#  @param  bkg_1y    y-background component for S(x)*B(y) term
#  @param  bkg_2x    x-background component for B(x)*B(y) term, if <code>bkg_2D</code> is not specified 
#  @param  bkg_2y    y-background component for B(x)*B(y) term, if <code>bkg_2D</code> is not specified 
#  @param  bkg_2D    PDF for 2D-background component for B(x,y)    term
#  @param  sig_2D    PDF for 2D-signal component for S(x,y)        term
#  @param  ss        the yield of  S(x)*S(y) component
#  @param  sb        the yield of  S(x)*B(y) component
#  @param  bs        the yield of  B(x)*S(y) component
#  @param  bb        the yield of  B(x,y)    component
#  @param  componens list of other 2D-components
#  @param  xvar      the x-variable
#  @param  yvar      the y-variable
#  @param  name      the name of PDF 
#  @code
#  model   = Models.Fit2D (
#      signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
#      signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y ) ,
#      bkg_1x   = 1 , 
#      bkg_1y   = 1 )
#
#  r,f = model.fitTo ( dataset ) ## fit dataset 
#
#  print r                       ## get results  
#
#  fx  = model.draw1 ()          ## visualize X-projection
#  fy  = model.draw2 ()          ## visualize Y-projection#
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class Fit2D (PDF2) :
    """The actual model for 2D-fits
    
    It consists of four main components :
    1. pure signal :        S(x,y)
    2. signal x background: S(x)*B(y)
    3. signal x backgronud: B(x)*S(y)
    4. pure backrground:    B(x,y) or B(x)*B(y) 
    Other 2D-components could be specified in addition
    
    Arguments:
    
    - signal_x        : PDF for the S(x)-signal component
    - signal_y        : PDF for the S(y)-signal component
    - suffix          : the suffix to be used for the PDF and variable names
    - bkg_1x          : x-background component for B(x)*S(y) term 
    - bkg_1y          : y-background component for S(x)*B(y) term
    - bkg_2x          : x-background component for B(x)*B(y) term, if bkg2D is not specified 
    - bkg_2y          : y-background component for B(x)*B(y) term, if bkg2D is not specified 
    - bkg_2D          : PDF for 2D-background component for B(x,y)    term
    - sig_2D          : PDF for 2D-signal component S(x,y) term
    - ss              : the yield of  S(x,y)    component
    - sb              : the yield of  S(x)*B(y) component
    - bs              : the yield of  B(x)*S(y) component
    - bb              : the yield of  B(x,y)    component
    - components      : the list of other 2D-components
    - xvar            : the x-variable
    - yvar            : the y-variable
    - name            : the name of PDF
    
    Example:
    
    >>>  model   = Models.Fit2D (
    ...      signal_x = Models.Gauss_pdf ( 'Gx' , mass = m_x ) ,
    ...      signal_y = Models.Gauss_pdf ( 'Gy' , mass = m_y ) ,
    ...      bkg1x    = 1 , 
    ...      bkg1y    = 1 )
    >>> r,f = model.fitTo ( dataset ) ## fit dataset 
    >>> print r                       ## get results  
    >>> fx  = model.draw1 ()          ## visualize X-projection
    >>> fy  = model.draw2 ()          ## visualize Y-projection

    """
    def __init__ ( self               ,
                   #
                   signal_x           , 
                   signal_y           ,
                   suffix = ''        ,
                   #
                   bkg_1x     = None  ,
                   bkg_1y     = None  ,
                   #
                   bkg_2x     = None  ,
                   bkg_2y     = None  ,
                   #
                   bkg_2D     = None  ,
                   sig_2D     = None  , ## 2D-signal component 
                   #
                   ## main components :
                   ss         = None  , ## signal    (1) * signal     (2)
                   sb         = None  , ## signal    (1) * background (2) 
                   bs         = None  , ## background(1) * signal     (2)
                   bb         = None  , ## background-2D 
                   ## additional components 
                   components = []    ,
                   xvar       = None  ,
                   yvar       = None  ,                   
                   name       = ''    ) :
        
        ## collect all the arguments 
        self.__args = {
            'signal_x'   : signal_x   ,
            'signal_y'   : signal_y   ,
            'bkg_1x'     : bkg_1x     ,
            'bkg_1y'     : bkg_1y     ,
            'bkg_2x'     : bkg_2x     ,
            'bkg_2y'     : bkg_2y     ,
            'bkg_2D'     : bkg_2D     ,
            'sig_2D'     : sig_2D     ,
            'components' : components ,
            ##
            'ss'         : ss ,
            'bb'         : bb ,
            'sb'         : sb ,
            'bs'         : bs ,
            ##
            'suffix'     : suffix   ,
            'name'       : name     ,
            }
        
        
        self.__suffix      = suffix

        self.__signal_x , xvar = self.make_PDF1 ( signal_x , xvar , prexis = 'SX' , suffix = suffix )
        self.__signal_y , yvar = self.make_PDF1 ( signal_y , yvar , prexis = 'SY' , suffix = suffix )
        
        #
        ## initialize base class
        #
        if not name :
            name = self.generate_name ( "fit2:%s&%s" % ( self.__signal_x.name , self.__signal_y.name ) )
            if suffix : name += '_' + suffix 
            
        PDF2.__init__ ( self , name          ,
                        self.__signal_x.xvar ,
                        self.__signal_y.xvar ) 

        
        # =====================================================================
        ## Build components for the  final 2D-PDF
        # =====================================================================
        
        # =====================================================================
        ## First component: Signal(1) and Signal(2)
        # =====================================================================

        ss_name = self.new_name ( 'SS' , suffix ) 
        if sig_2D : self.__ss_cmp = self.make_PDF2 ( sig_2D ,
                                                     xvar = self.xvar ,
                                                     yvar = self.yvar ,
                                                     name = ss_name   ) [ 0 ] 
        else      : self.__ss_cmp  = Model2D ( ss_name         ,
                                               self.__signal_x ,
                                               self.__signal_y , 
                                               title = "Signal(x) x Signal(y)" )
        
        # =====================================================================
        ## Second component: Background(1) and Signal(2)
        # =====================================================================

        self.__bkg_1x = self.make_bkg ( bkg_1x  , self.new_name ( 'Bkg1X_BS', suffix ) , self.xvar )
        self.__bs_cmp = Model2D ( self.new_name ( "BS" ,  suffix ) ,
                                  self.__bkg_1x             ,
                                  self.__signal_y           ,
                                  title = "Backround1(x) x Signal(y)" )
        
        # =====================================================================
        ## Third component:  Signal(1) and Background(2)
        # =====================================================================
        
        self.__bkg_1y = self.make_bkg ( bkg_1y   , self.new_name ( 'Bkg1Y_SB' , suffix ) , self.yvar )
        self.__sb_cmp = Model2D ( self.new_name ( "SB" ,  suffix ) ,
                                  self.__signal_x           ,
                                  self.__bkg_1y             ,
                                  title = "Signal(x) x Background1(y)" )
            
        # =====================================================================
        ## (intermezzo) assumptions about the background sub-components 
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
            
        # =====================================================================
        ## Fourth component: Background(1) and Background(2) 
        # =====================================================================
    
        self.__bkg_2x = None 
        self.__bkg_2y = None 

        bb_name = self.new_name ( "BB" , suffix  )
        if bkg_2D and isinstance ( bkg_2D , ( tuple , list ) ) :
            from ostap.fitting.models_2d import make_B2D
            self.__bb_cmp = make_B2D ( bb_name , self.xvar , self.yvar , *bkg_2D )
        elif bkg_2D  :            
            self.__bb_cmp = self.make_PDF2 ( bkg_2D           ,
                                             xvar = self.xvar ,
                                             yvar = self.yvar ,
                                             name = bb_name   )[0]
        else      : 
            self.__bkg_2x = self.make_bkg ( bkg_2x , self.new_name ( 'Bkg2X_BB' , suffix ) , self.xvar )
            self.__bkg_2y = self.make_bkg ( bkg_2y , self.new_name ( 'Bkg2Y_BB' , suffix ) , self.yvar )            
            self.__bb_cmp = Model2D ( bb_name       ,
                                      self.__bkg_2x ,
                                      self.__bkg_2y ,
                                      title = "Background2(x) x Background2(y)" )
            
        # =====================================================================
        ## coefficients/yields 
        # =====================================================================
    
        self.__ss = self.make_var ( ss   , "SS"               + suffix ,
                                    "Signal(x,y)"             + suffix , None , 1000  , 0 , 1.e+7 )
        self.__sb = self.make_var ( sb   ,  "SB"              + suffix ,
                                    "Signal(x)&Background(y)" + suffix , None ,  100  , 0 , 1.e+7 )
        self.__bs = self.make_var ( bs   , "BS"               + suffix ,
                                    "Background(x)&Signal(y)" + suffix , None ,  100  , 0 , 1.e+7 )        
        self.__bb = self.make_var ( bb   , "BB"               + suffix ,
                                    "Background(x,y)"         + suffix , None ,   10  , 0 , 1.e+7 )
        
        self.alist1 = ROOT.RooArgList (
            self.__ss_cmp.pdf ,
            self.__sb_cmp.pdf ,
            self.__bs_cmp.pdf ,
            self.__bb_cmp.pdf )
        self.alist2 = ROOT.RooArgList (
            self.__ss ,
            self.__sb ,
            self.__bs ,
            self.__bb )

        ## treat additional components (if specified)
        self.__nums_components = [] 
        icmp = 0
        self.__more_components = [
            self.make_PDF2 ( cmp                 ,
                             xvar   = self.xvar  ,
                             yvar   = self.yvar  ,
                             prefix = 'C%d_' % i ,
                             suffix = suffix     ) [0]  for (i,cmp) in enumerate ( components ) ]
        
        for cmp in self.__more_components : self.components.add (  cmp.pdf ) 

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

        #
        ## build the final PDF 
        # 
        pdfname  = self.new_roo_name ( 'fit2d' , suffix ) 
        pdftitle = "Fit2D %s" % self.name
        pdfargs  = pdfname , pdftitle , self.alist1 , self.alist2
        self.pdf = ROOT.RooAddPdf  ( *pdfargs )

        self.signals     .add ( self.__ss_cmp.pdf )
        self.backgrounds .add ( self.__bb_cmp.pdf )
        self.crossterms1 .add ( self.__sb_cmp.pdf ) ## cross-terms 
        self.crossterms2 .add ( self.__bs_cmp.pdf ) ## cross-terms 

        ## save configuration
        self.config = {
            'signal_x'   : self.signal_x        ,
            'signal_y'   : self.signal_y        ,            
            'suffix'     : self.suffix          ,
            'bkg_1x'     : self.bkg_1x          , 
            'bkg_1y'     : self.bkg_1y          , 
            'bkg_2x'     : self.bkg_2x          , 
            'bkg_2y'     : self.bkg_2y          , 
            'bkg_2D'     : self.bkg_2D          ,
            'sig_2D'     : self.sig_2D          ,
            'ss'         : self.SS              ,
            'sb'         : self.SB              ,
            'bs'         : self.BS              ,
            'bb'         : self.BB              ,
            'components' : self.more_components ,
            'xvar'       : self.xvar            , 
            'yvar'       : self.yvar            , 
            'name'       : self.name    
            }
        
        self.checked_keys.add  ( 'xvar' )
        self.checked_keys.add  ( 'yvar' )
        
    @property
    def SS ( self ) :
        """The yield of Signal(x)*Signal(y) component"""
        return self.__ss
    @SS.setter 
    def SS ( self , value ) :
        value = float ( value  )
        assert value in self.__ss, "Value %s is out of the allowed range %s " % ( value , self.__ss.minmax() )
        self.__ss.setVal ( value ) 

    @property
    def SB ( self ) :
        """The yield of Signal(x)*Background(y) component"""
        return self.__sb
    @SB.setter 
    def SB ( self , value ) :
        value = float ( value  )
        assert value in self.__sb, "Value %s is out of the allowed range %s " % ( value , self.__sb.minmax() )
        self.__sb.setVal ( value ) 

    @property
    def BS ( self ) :
        """The yield of Background(x)*Signal(y) component"""
        return self.__bs
    @BS.setter 
    def BS ( self , value ) :
        value = float ( value  )
        assert value in self.__bs, "Value %s is out of the allowed range %s " % ( value , self.__bs.minmax() )
        self.__bs.setVal ( value ) 

    @property
    def BB ( self ) :
        """The yield of Background(x,y) component"""
        return self.__bb
    @BB.setter 
    def BB ( self , value ) :
        value = float ( value  )
        assert value in self.__bb, "Value %s is out of the allowed range %s " % ( value , self.__bb.minmax() )
        self.__bb.setVal ( value ) 
    
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
        if hasattr ( self , 'extended' ) and not self.extended : return None 
        if not self.fit_result                                 : return None
        if not valid_pointer ( self.fit_result )               : return None
        yields = self.yields
        if not yields                                          : return None
        ##  if 1 ==  len ( yields ) : return yield[0]. 
        return self.fit_result.sum ( *yields ) 
 

    # =========================================================================
    # components
    # =========================================================================

    @property 
    def signal_x ( self  ) :
        """Signal(x) component/PDF"""
        return self.__signal_x

    @property 
    def signal_y ( self  ) :
        """Signal(y) component/PDF"""
        return self.__signal_y

    @property
    def bkg_1x( self ) :
        """'bkg_1x' : The background(x) PDF for Backgroud(x)*Signal(y) component/PDF"""
        return self.__bkg_1x
    
    @property
    def bkg_1y( self ) :
        """'bkg_1y' : The background(y) PDF for Signal(x)*Background(y) component/PDF"""
        return self.__bkg_1y

    @property
    def bkg_2x( self ) :
        """'bkg_2x' : The background(x) PDF for Backgroud(x)*Background(y) component/PDF, when bkg2D is not specified"""
        return self.__bkg_2x

    @property
    def bkg_2y( self ) :
        """'bkg_2y' : The background(y) PDF for Backgroud(x)*Background(y) component/PDF, when bkg2D is not specified"""
        return self.__bkg_2y

    @property
    def bkg_2D( self ) :
        """'bkg_2D' : The PDF for Backgroud(x,y) component/PDF (same as cmp_BB)"""
        return self.__bb_cmp

    @property
    def sig_2D( self ) :
        """'sig_2D' : The PDF for Signal(x,y) component/PDF (same as cmp_SS)"""
        return self.__ss_cmp 

    @property
    def more_components ( self ) :
        """additional 'other' components"""
        return tuple( self.__more_components  )

    @property
    def suffix ( self ) :
        """'suffix' , used to build the name"""
        return self.__suffix
    
    @property
    def cmp_SS ( self ) :
        """'cmp_SS'  : Signal(x&y))            component in the fit (PDF)"""
        return self.__ss_cmp
    @property
    def cmp_SB ( self ) :
        """'cmp_SB'  : Signal(x)xBackground(y) component in the fit (PDF)"""
        return self.__sb_cmp
    @property
    def cmp_BS ( self ) :
        """'cmp_BS'  : Background(x)xSignal(y) component in the fit (PDF)"""
        return self.__bs_cmp
    @property
    def cmp_BB ( self ) :
        """'BB_cmp'  : Background(x,y)         component in the fit (PDF)"""
        return self.__bb_cmp
    
# =============================================================================
## @class Fit2DSym
#  The actual model for symmetric 2D-fits. It consists of three main components :
#  1. pure signal :        S(x,y)
#  2. signal x background & background x signal:  S(x)*B(y) + B(x)*S(y)
#  3. pure backrground:    B(x,y) or B(x)*B(y) 
#  Other 2D-components could be specified in addition
#  @param  signal_x  PDF for the S(x)-signal component
#  @param  signal_y  PDF for the S(y)-signal component, cloned from S(x), if None 
#  @param  suffix    suffix to be used for the PDF and variable names
#  @param  bkg_1x    x-background component for B(x)*S(y) term, B(y) is cloned 
#  @param  bkg_2x    x-background component for B(x)*B(y) term, if bkg2D is not specified, B(y) is cloned  
#  @param  bkg_2D     PDF for (symmetric) 2D-background component for B(x,y)    term
#  @param  sig_2D     PDF for (symmetric) 2D-signal     component for S(x,y)    term
#  @param  ss         the yield of  S(x)*S(y) component
#  @param  sb         the yield of  S(x)*B(y)+B(x)*S(y)component
#  @param  bb         the yield of  B(x,y)    component
#  @param  components list of other 2D-components
#  @param  xvar       the x-variable
#  @param  yvar       the y-variable
#  @param  name       the name of PDF 

#  @code
# 
#  model   = Models.Fit2D (
#      signal_x = Models.Gauss_pdf ( 'Gx' , mass = m_x ) ,
#      signal_y = Models.Gauss_pdf ( 'Gy' , mass = m_y ) ,
#      bkg_1x   = 1 , 
#      bkg_2x   = 1 )
#
#  r,f = model.fitTo ( dataset ) ## fit dataset 
#
#  print r                       ## get results  
#
#  fx  = model.draw1 ()          ## visualize X-projection
#  fy  = model.draw2 ()          ## visualize X-projection
#
#  @endcode 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class Fit2DSym (PDF2) :
    """The actual model for *symmetric**2D-fits
    The actual model for symmetric 2D-fits. It consists of three main components :
    - pure signal                               :        S(x)*S(y)
    - signal x background & background x signal :  S(x)*B(y) + B(x)*S(y)
    - pure background:                          :  B(x,y) or B(x)*B(y) 
    Other 2D-components could be specified in addition
    
    Arguments:
    
    - signal_x       : PDF for the S(x)-signal component
    - signal_y       : PDF for the S(y)-signal component; cloned from S(x), if None 
    - suffix         : suffix to be used for the PDF and variable names
    - bkg_1x         : x-background component for B(x)*S(y) term; B(y) is cloned 
    - bkg_2x         : x-background component for B(x)*B(y) term, if bkg2D is not specified; B(y) is   cloned  
    - bkg_2D         : PDF for (symmetric) 2D-background component for B(x,y)    term
    - sig_2D         : PDF for (symmetric) 2D-signal     component for S(x,y)    term
    - ss             : the yield of  S(x)*S(y) component
    - sb             : the yield of  S(x)*B(y)+B(x)*S(y)component
    - bb             : the yield of  B(x,y)    component
    - componens      : the list of other 2D-components
    - xvar           : the x-variable
    - yvar           : the y-variable
    - name           : the name of PDF 
    
    Example:
    
    >>>  model   = Models.Fit2DSym (
    ...      signal_x = Models.Gauss_pdf ( 'Gx' , xvar = m_x ) ,
    ...      signal_y = Models.Gauss_pdf ( 'Gy' , xvar = m_y ) ,
    ...      bkg1x    = 1 , 
    ...      bkg2x    = 1 )
    >>> r,f = model.fitTo ( dataset ) ## fit dataset 
    >>> print r                       ## get results  
    >>> fx  = model.draw1 ()          ## visualize X-projection
    >>> fy  = model.draw2 ()          ## visualize X-projection

    """
    def __init__ ( self               ,
                   #
                   signal_x           , 
                   signal_y   = None  ,
                   suffix     = ''    ,
                   #
                   bkg_1x     = None  ,
                   bkg_2x     = None  ,
                   bkg_2D     = None  ,
                   sig_2D     = None  ,
                   #
                   ## main components :
                   ss         = None  , ## signal (1) * signal     (2)
                   sb         = None  , ## signal     * background 
                   bb         = None  , ## background * background  
                   ## additional components 
                   components = []    ,
                   xvar       = None  ,
                   yvar       = None  ,                   
                   name       = ''    ) :
        
        ## collect all the arguments 
        self.__args = {
            'signal_x'   : signal_x   ,
            'signal_y'   : signal_y   ,
            'bkg_1x'     : bkg_1x     , 
            'bkg_2x'     : bkg_2x     , 
            'bkg_2D'     : bkg_2D     ,
            'sig_2D'     : sig_2D     ,
            'components' : components ,
            ##
            'ss'         : ss         ,
            'sb'         : sb         , 
            'bb'         : bb         ,
            ##
            'suffix'     : suffix     ,
            'name'       : name       ,
            ##
            'xvar'       : xvar       ,
            'yvar'       : yvar       ,            
            }
        
        self.__suffix      = suffix

        # ============================================================================

        self.__signal_x , xvar = self.make_PDF1 ( signal_x , xvar , prefix = 'SX' , suffix = suffix )
        self.__signal_y , yvar = self.make_PDF1 ( signal_y , yvar , prefix = 'SY' , suffix = suffix )
            
        # =====================================================================
        ## initialize base class
        # =====================================================================
        if not name :
            name = "%s&%s" % ( self.__signal_x.name , self.__signal_y.name )
            if suffix : name += '_' + suffix 

        ##  initialize the base class 
        PDF2.__init__ ( self , name          ,
                        self.__signal_x.xvar ,
                        self.__signal_y.xvar ) 
        
        # =====================================================================
        ## First component: Signal(1) and Signal(2)
        # =====================================================================

        ss_name = self.new_name ( 'SS' , suffix  )
        if sig_2D : self.__ss_cmp = self.make_PDF2 ( sig_2D ,
                                                     xvar = self.xvar ,
                                                     yvar = self.yvar ,
                                                     name  = ss_name ) [ 0 ] 
        else      : self.__ss_cmp = Model2D ( ss_name         ,
                                              self.__signal_x ,
                                              self.__signal_y , 
                                              title = "Signal(x) x Signal(y)" )
        
        self.__bkg_1x = self.make_bkg (        bkg_1x , self.new_name ( 'Bkg1X_BS' , suffix ) , self.xvar )
        self.__bkg_1y = self.make_bkg ( self.__bkg_1x , self.new_name ( 'Bkg1Y_SB' , suffix ) , self.yvar )

        # =====================================================================
        ## Second sub-component: Background (1) and Signal     (2)
        ## Third  sub-component: Signal     (1) and Background (2)
        # =====================================================================

        self.__sb_cmp_raw = Model2D ( self.new_name ( "S1B2" , suffix ) ,
                                      self.__signal_x           ,
                                      self.__bkg_1y             ,
                                      title = "Signal(x) x Background(y)" )
        
        self.__bs_cmp_raw = Model2D ( self.new_name ( "B1S2" , suffix ) , 
                                      self.__bkg_1x             ,
                                      self.__signal_y           ,    
                                      title = "Background(x) x Signal(y)" )

        sb_name = self.new_name ( "SB" , suffix ) 
        self.__sb_cmp     = Generic2D_pdf (
            self.make_sum ( sb_name ,
                            "Signal(x) x Background(y) + Background(x) x Signal(y)"   ,
                            self.__sb_cmp_raw.pdf ,
                            self.__bs_cmp_raw.pdf ) , self.xvar , self.yvar , sb_name ) 
        
        ## alias, just for convinience 
        self.__bs_cmp    = self.__sb_cmp
        
        # =====================================================================
        ## (intermezzo) Assumptions about the background sub-components 
        # =====================================================================
        
        if   component_clone   ( bkg_2x ) :
            bkg_2x = self.__bkg_1x
            self.debug ( 'bkg_2x set to [CLONE]   %s' % bkg_2x ) 
        elif component_similar ( bkg_2x ) :
            bkg_2x =        bkg_1x
            self.debug ( 'bkg_2x set to [SIMILAR] %s' % bkg_2x ) 
            
        # =====================================================================
        ## fourth component: Background(1) and Background(2) 
        # =====================================================================
    
        self.__bkg_2x = None
        self.__bkg_2y = None

        bb_name = self.generate_name ( 'BB_' + self.name )        
        if  isinstance  ( bkg_2D , int ) :            
            from ostap.fitting.models_2d import make_B2Dsym
            self.__bb_cmp = make_B2Dsym ( bb_name , self.xvar , self.yvar , bkg_2D )            
        elif bkg_2D  :            
            self.__bb_cmp = self.make_PDF2 ( bkg_2D           ,
                                             xvar = self.xvar ,
                                             yvar = self.yvar ,
                                             name = bb_name   )[0]
        else     :                        
            self.__bkg_2x = self.make_bkg (        bkg_2x , self.generate_name ( 'Bkg2X_BB' + self.name ) , self.xvar )
            self.__bkg_2y = self.make_bkg ( self.__bkg_2x , self.generate_name ( 'Bkg2Y_BB' + self.name ) , self.yvar )
            self.__bb_cmp = Model2D ( bb_name       ,
                                      self.__bkg_2x ,
                                      self.__bkg_2y ,
                                      title = "Background2(x) x Backrgound2(y)" )
            
        # =====================================================================
        ## coefficients
        # =====================================================================
        
        self.__ss = self.make_var ( ss , "SS"             + suffix ,
                                    "Signal(x)&Signal(y)" + suffix , None , 1000  , 0 ,  1.e+7 )
        
        self.__bb = self.make_var ( bb , "BB"             + suffix ,
                                    "Background(x,y)"     + suffix , None ,   10  , 0 ,  1.e+7 )
        
        self.__sb = self.make_var ( sb , "SB"             + suffix ,
                                    "Signal(x)&Background(y)+Background(x)&Signal(y)" + suffix , None ,  100  , 0 ,  1.e+7 )
        
        ## duplicate
        
        self.__bs = self.__sb
        
        self.alist1 = ROOT.RooArgList (
            self.__ss_cmp.pdf ,
            self.__sb_cmp.pdf ,
            self.__bb_cmp.pdf )
        self.alist2 = ROOT.RooArgList (
            self.__ss         ,
            self.__sb         ,
            self.__bb         )

        ## treat additional components (if specified)
        self.__nums_components = [] 
        icmp = 0
        self.__more_components = [
            self.make_PDF2 ( cmp                 ,
                             xvar   = self.xvar  ,
                             yvar   = self.yvar  ,
                             prefix = 'C%d_' % i ,
                             suffix = suffix     ) [ 0 ] for ( i , cmp ) in enumerate ( components ) ]
        for cmp in self.__more_components : self.components.add ( cmp.pdf )

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
            
        #
        ## build the final PDF 
        #
        pdfname  = self.new_roo_name ( 'fit2ds' , suffix ) 
        pdftitle = "Fit2Dsym %s" % self.name
        pdfargs  = pdfname , pdftitle , self.alist1 , self.alist2
        self.pdf = ROOT.RooAddPdf  ( *pdfargs )


        self.signals     .add ( self.__ss_cmp.pdf )
        self.backgrounds .add ( self.__bb_cmp.pdf )
        self.crossterms1 .add ( self.__sb_cmp.pdf ) ## cross-terms 
        self.crossterms2 .add ( self.__bs_cmp.pdf ) ## cross-terms 

        ## save configuration
        self.config = {
            'signal_x'   : self.signal_x        ,
            'signal_y'   : self.signal_y        ,            
            'suffix'     : self.suffix          ,
            'bkg_1x'     : self.bkg_1x          , 
            'bkg_2x'     : self.bkg_2x          , 
            'bkg_2D'     : self.bkg_2D          ,
            'sig_2D'     : self.sig_2D          ,
            'ss'         : self.SS              ,
            'sb'         : self.SB              ,
            'bb'         : self.BB              ,
            'components' : self.more_components ,
            'xvar'       : self.xvar            , 
            'yvar'       : self.yvar            , 
            'name'       : self.name            ,
            }

        self.checked_keys.add ( 'xvar' )
        self.checked_keys.add ( 'yvar' )
        
    @property
    def SS ( self ) :
        """The yield of Signal(x)*Signal(y) component"""
        return self.__ss
    @SS.setter 
    def SS ( self , value ) :
        value = float ( value  )
        assert value in self.__ss, "Value %s is out of the allowed range %s " % ( value , self.__ss.minmax() )
        self.__ss.setVal ( value ) 

    @property
    def SB ( self ) :
        """The yield of Signal(x)*Background(y)+Background(x)*Signal(y) component (same as 'BS')"""
        return self.__sb
    @SB.setter 
    def SB ( self , value ) :
        value = float ( value  )
        assert value in self.__sb, "Value %s is out of the allowed range %s " % ( value , self.__sb.minmax() )
        self.__sb.setVal ( value ) 

    @property
    def BS ( self ) :
        """The yield of Signal(x)*Background(y)+Background(x)*Signal(y) component (same as 'SB')"""
        return self.SB 
    @BS.setter 
    def BS ( self , value ) :
        self.SB = value
        return self.SB.getVal()

    @property
    def BB ( self ) :
        """The yield of Background(x,y) component"""
        return self.__bb
    @BB.setter 
    def BB ( self , value ) :
        value = float ( value  )
        assert value in self.__bb, "Value %s is out of the allowed range %s " % ( value , self.__bb.minmax() )
        self.__bb.setVal ( value ) 
    
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
        return self.component_getter ( self.__nums_components )
    @C.setter
    def C (  self , value ) :
        self.component_setter ( self.__nums_componenents , value )

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
    def bkg_1x( self ) :
        """'bkg_1x' : The background PDF for Backgroud(x)*Signal(y) component/PDF"""
        return self.__bkg_1x
    
    @property
    def bkg_1y( self ) :
        """'bkg_1y' : The background PDF for Signal(x)*Background(y) component/PDF"""
        return self.__bkg_1y 

    @property
    def bkg_2x( self ) :
        """'bkg_2x' : The background(x) PDF for Backgroud(x)*Background(y) component/PDF"""
        return self.__bkg_2x
    
    @property
    def bkg_2y( self ) :
        """'bkg_2y' : The background(y) PDF for Backgroud(x)*Background(y) component/PDF"""
        return self.__bkg_2y 

    @property
    def bkg_2D( self ) :
        """'bkg_2D' : The PDF for Backgroud(x&y) component/PDF (same as cmp_BB)"""
        return self.__bb_cmp

    @property
    def sig_2D( self ) :
        """'sig_2D' : The PDF for Signal(x&y) component/PDF (same as cmp_SS)"""
        return self.__ss_cmp
 
    @property
    def more_components ( self ) :
        """additional 'other' components"""
        return tuple( self.__more_components  )

    @property
    def suffix ( self ) :
        """'suffix' , used to build the name"""
        return self.__suffix
    
    @property
    def cmp_SS ( self ) :
        """'cmp_SS'  : Sig(1&2) component in the fit (PDF)"""
        return self.__ss_cmp
    @property
    def cmp_SB ( self ) :
        """'cmp_SB'  : Sig(1)xBkg(2)+Bkg(1)*Sig(2) component in the fit (PDF)"""
        return self.__sb_cmp
    @property
    def cmp_BS ( self ) :
        """'cmp_BS'  : Sig(1)xBkg(2)+Bkg(1)*Sig(2) component in the fit (PDF) (same as  cmp_SB)"""
        return self.__bs_cmp
    @property
    def cmp_BB ( self ) :
        """'cmp_BB'  : Bkg(1&2)      component in the fit (PDF)"""
        return self.__bb_cmp
 
    # =========================================================================
    ## Raw, non-symmetrized fit components/PDF (for debugging)
    # =========================================================================

    def cmp_raw_SB ( self ) :
        """'cmp_SB'  : Sig(1)xBkg(2) raw, non-symmetrized component in the fit (PDF)"""
        return self.__sb_cmp_raw
    @property
    def cmp_BS ( self ) :
        """'cmp_BS'  :  Bkg(1)*Sig(2) raw, non-symmetrized component in the fit (PDF)"""
        return self.__bs_cmp_raw

 
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
