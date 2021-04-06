#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/modifiers.py
#  Set of useful utilities to modify the certain PDF 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2018-11-29
# =============================================================================
"""Set of useful utilities to modify the certain PDF """
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2018-11-29"
__all__     = (
    'Product1D_pdf'     , ## make a product of 1D PDFs
    'Modify1D_pdf'      , ## modify 1D pdf with a positive polynomials
    'CutOff_pdf'        , ## Cut-off PDF
    'CutOffGauss_pdf'   , ## Gaussian cut-off
    'CutOffStudent_pdf' , ## Student's t-like/power law cut-off    
    )
# =============================================================================
import ROOT, math
from   ostap.fitting.basic import PDF , Generic1D_pdf
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.modifiers' )
else                       : logger = getLogger ( __name__                  )
# =============================================================================
models = []
# =============================================================================
## @class Product1D_pdf
#  Simple product of 1D-pdfs
#  - trivial wrapper for RooProdPdf
#  @code
#  pdf1     = ...
#  pdf2     = ...
#  pdf_prod = Product1D_pdf( pdf1 , pdf2 ) 
#  @endcode
#  @see RooProdPdf
#  @attention it could be rather CPU-inefficient!
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2018-11-29  
class Product1D_pdf(PDF) :
    """Simple product of 1D-PDFs
    - actually it is a trivial wrapper for RooProdPdf
    
    >>> pdf1 = ...
    >>> pdf2 = ...
    >>> pdf  = Product1D_pdf( pdf1 , pdf2 )
    
    - attention: it could be rather CPU-inefficient!
    """
    def __init__ ( self         ,
                   pdf1         ,
                   pdf2         ,
                   xvar  = None ,
                   name  = ''   ,
                   title = ''   ) :
        
        self.__pdf__1 = pdf1
        self.__pdf__2 = pdf2

        self.__pdf1 = None 
        self.__pdf2 = None

        if isinstance ( pdf1 , PDF ) :
            self.__pdf1 = pdf1
            if xvar and not ( xvar is pdf1.xvar ) :
                self.error ("Mismatch in ``xvar''-observable (1) %s/%s" % ( xvar , pdf1.xvar ) ) 
            elif not xvar : xvar = pdf1.xvar 
        elif isinstance ( pdf1 , ROOT.RooAbsPdf ) and xvar :
            self.__pdf1 = Generic1D_pdf ( pdf1 , xvar )
        else :
            raise TypeError("Illegal setting for ``pdf1'': %s/%s" % ( pdf1 , type ( pdf1 ) ) )

        assert isinstance ( self.__pdf1  , PDF ), 'Invalid pdf1 type'

        if isinstance ( pdf2 , PDF ) :
            self.__pdf2 = pdf2
            if xvar and not ( xvar is pdf2.xvar ) :
                self.error ("Mismatch in ``xvar''-observable (2) %s/%s" % ( xvar , pdf2.xvar ) ) 
            elif not xvar : xvar = pdf2.xvar 
        elif isinstance ( pdf2 , ROOT.RooAbsPdf ) and xvar :
            self.__pdf2 = Generic1D_pdf ( pdf2 , xvar )
        else :
            raise TypeError("Illegal setting for ``pdf2'': %s/%s" % ( pdf2 , type ( pdf1 ) ) )
        
        assert isinstance ( self.__pdf2  , PDF ), 'Invalid pdf2 type'
        
        assert isinstance ( xvar , ROOT.RooAbsReal ),\
               "Invalid ``xvar'':%s/%s" % ( xvar , type  ( xvar ) ) 

        name = name if name else self.generate_name ( prefix = "product_%s_%s_"  % ( self.pdf1.name , self.pdf2.name ) )

        ## initialize the base class
        PDF.__init__ ( self , name , xvar =  xvar )

        em1 = self.pdf1.pdf.extendMode()
        em2 = self.pdf2.pdf.extendMode()
        
        if   2 == em1 : self.warning ( "pdf1 ``must-be-extended''" )
        elif 1 == em1 : self.warning ( "pdf1  ``can-be-extended''" )
        
        if   2 == em2 : self.warning ( "pdf2 ``must-be-extended''" )
        elif 1 == em2 : self.warning ( "pdf2  ``can-be-extended''" )

        ## finally build PDF 
        self.pdf = ROOT.RooProdPdf (
            self.roo_name ( 'prod1_' ) ,
            title if title else 'Product of two pdfs %s' % self.name , 
            self.pdf1.pdf  ,
            self.pdf2.pdf )

        ## save configuration for cloning
        self.config = {
            'pdf1'  : self.pdf1  ,
            'pdf2'  : self.pdf2  ,
            'xvar'  : self.xvar  ,
            'name'  : self.name  ,
            'title' : self.title ,            
            }
        
    @property
    def pdf1 ( self ) :
        """``pdf1'' : the first PDF"""
        return self.__pdf1
    
    @property
    def pdf2 ( self ) :
        """``pdf2'' : the second PDF"""
        return self.__pdf2

    ## redefine the clone 
    def clone ( self , **kwargs ) :
        """ Redefine the clone
        """
        
        pdf1 = kwargs.pop ( 'pdf1' , None )
        pdf2 = kwargs.pop ( 'pdf2' , None )
        
        if not pdf1 : pdf1 = self.pdf1.clone ( **kwargs )
        if not pdf2 : pdf2 = self.pdf2.clone ( **kwargs )

        return FUNC.clone ( self , pdf1 = pdf1 , pdf2 = pdf2 , **kwargs ) 
    
        

models.append ( Product1D_pdf ) 
# =============================================================================
## @class Modify1D_pdf
#  Modify the certain PDF with positive polynomial function
#  @code
#  pdf1 = ...
#  pdf  = Modify1D ( pdf1 , power = 1  ) 
#  @endcode
#  @attention it could be rather CPU-inefficient!
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2018-11-29  
class Modify1D_pdf(Product1D_pdf) :
    """Modify the certain PDF with positive polynomial function
    
    >>> pdf1 = ...
    >>> pdf  = Modify1D ( pdf1 , power = 1  )
    
    - attention: it could be rather CPU-inefficient!
    """
    def __init__ ( self         ,
                   pdf          ,
                   power = 1    ,
                   xvar  = None ,
                   name  = ''   ,
                   xmin  = None ,
                   xmax  = None ,
                   title = ''   ) :
        
        assert isinstance ( power , int ) and 0 <= power,\
               "Invalid ``power''   %s" % power

        self.__pdf_1 = pdf
        
        if isinstance ( pdf , PDF ) : 
            if xvar and not ( xvar is pdf.xvar ) :
                self.error ("Mismatch in ``xvar''-observable") 
            elif not xvar : xvar = pdf.xvar            
        elif isinstance ( pdf , ROOT.RooAbsPdf ) and xvar :
            pdf = Generic1D_pdf ( pdf , xvar )            
        else :
            raise TypeError("Illegal setting for ``pdf'': %s/%s" % ( pdf , type ( pdf ) ) )

        assert isinstance ( xvar , ROOT.RooAbsReal ),\
               "Invalid ``xvar'':%s/%s" % ( xvar , type  ( xvar ) ) 

        name = name if name else self.generate_name ( prefix = "modify_%s_%s"  % ( pdf.name , power ) )

        from ostap.fitting.background import PolyPos_pdf
        pdf2 = PolyPos_pdf( self.generate_name ( 'M_%s_%s' % ( pdf.name , power ) ) ,
                            xvar  = xvar  ,
                            power = power ,
                            xmin  = xmin  ,
                            xmax  = xmax  )
        
        self.__pdf_2 = pdf2
        
        ## initialize the base
        Product1D_pdf.__init__ ( self          ,
                                 pdf1  = pdf   ,
                                 pdf2  = pdf2  , 
                                 xvar  = xvar  ,
                                 name  = name  ,
                                 title = title )         

        ## for drawing...
        
        for c in self.pdf1.signals     : self.signals    .add ( c )
        for c in self.pdf1.backgrounds : self.backgrounds.add ( c )
        for c in self.pdf1.components  : self.components .add ( c )
        for c in self.pdf1.crossterms1 : self.crossterms1.add ( c )
        for c in self.pdf1.crossterms2 : self.crossterms2.add ( c )
        
        self.config = {
            'name'  : self.name  ,
            'pdf'   : self.pdf1  ,
            'xvar'  : self.xvar  ,
            'power' : power      ,
            'xmin'  : xmin       ,
            'xmax'  : xmin       ,
            'title' : self.title ,            
            }

    @property
    def old_pdf ( self ):
        """``old_pdf''  : original (non-modifier) PDF"""
        return self.pdf1
    
    ## redirect any other attributes to original PDF
    ## def __getattr__ ( self , attr ) :
    ##    """Get all extra attributes from the original PDF"""
    ##    print 'getting %s attribute' % attr 
    ##    opdf = self.pdf1
    ##    return  getattr ( opdf , attr )


models.append ( Modify1D_pdf )


# =============================================================================
# @class CutOffGauss_pdf
#  Useful function for smooth Gaussian cut-off:
#  \f[ f(x;x_0;\sigma) = \left\{ 
#    \begin{array}{ll}
#    1  & \mathrm{for~} x \le x_0  \\
#    \mathrm{e}^{-\frac{1}{2}\left( \frac{ (x-x_0)^2}{\sigma^2} \right)}
#       & \mathrm{for~} x >   x_0 
#    \end{array}\right. \f] 
# @see Ostap::Math::CutOffGauss
# @see Ostap::Models::CutOffGauss
class CutOffGauss_pdf(PDF) :
    """ Useful function for smooth Gaussian-like  cut-off
    - see Ostap.Math.CutOffGauss
    - see Ostap.Models.CutOffGauss
    """
    def __init__ ( self  ,
                   name  ,
                   xvar  ,
                   right ,
                   x0    ,
                   sigma ) :
        
        PDF.__init__ ( self , name , xvar = xvar ) 
        
        self.__x0   = self.make_var ( x0       ,
                                      'x0_%s'      % name ,
                                      '#x_{0}(%s)' % name , x0 , x0 )
        
        self.__sigma = self.make_var ( sigma   ,
                                       'sigma_%s'        % name ,
                                       '#sigma_{CB}(%s)' % name , sigma , sigma , 0 , 1.e+6 ) 


        self.__right = True if right else False
        
        self.pdf = Ostap.Models.CutOffGauss (
            self.roo_name ( 'cofg_'  )          , 
            'Gaussian cut-off %s' % self.name  ,
            self.xvar                          ,
            self.right                         , 
            self.x0                            ,
            self.sigma                         ) 
        
        ## save the configuration
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'x0'    : self.x0    ,
            'right' : self.right , 
            'sigma' : self.sigma ,
            }
        
    @property
    def x0 ( self ) :
        """``x0'' : threshold/location parameter for Gaussial cut-off"""
        return self.__x0
    @x0.setter 
    def x0 ( self , value ) :
        self.set_value ( self.__x0 , value ) 

    @property
    def sigma ( self ) :
        """``sigma'' : width parameter for Gaussial cut-off"""
        return self.__sigma
    @sigma.setter 
    def sigma ( self , value ) :
        self.set_value ( self.__sigma , value ) 
    
    @property
    def right ( self ) :
        """``right'' : parameter of the Gaussian cut-off"""
        return self.__right 

    @property
    def left ( self ) :
        """``left'' : parameter of the Gaussian cut-off"""
        return not self.right 


models.append ( CutOffGauss_pdf )



# =============================================================================
# @class CutOffStudent_pdf
#  Useful function for smooth Student's t=-like (power-law) cut-off:
#  \f[ f(x;x_0;\sigma) = \left\{ 
#    \begin{array}{ll}
#    1  & \mathrm{for~} x \le x_0  \\
#    \left( \frac{1}{\nu} \left( \frac{(x-x_0)}{\sigma^2} \right)^{ - \frac{\nu+1}{2}} \right) 
#       & \mathrm{for~} x >   x_0 
#    \end{array}\right. \f] 
# @see Ostap::Math::CutOffGauss
# @see Ostap::Models::CutOffGauss
class CutOffStudent_pdf(PDF) :
    """ Useful function for smooth Student's t=-like (power-law) cut-off:
    - see Ostap.Math.CutOffStudent
    - see Ostap.Models.CutOffStudent
    """
    def __init__ ( self  ,
                   name  ,
                   xvar  ,
                   right ,
                   x0    ,
                   nu    ,
                   sigma ) :
        
        PDF.__init__ ( self , name , xvar = xvar ) 
        
        self.__x0   = self.make_var ( x0       ,
                                      'x0_%s'      % name ,
                                      '#x_{0}(%s)' % name , x0 , x0 )
        
        self.__nu   = self.make_var ( nu    ,
                                      'nu_%s'        % name ,
                                      '#nu_{CB}(%s)' % name , nu , nu , 0 , 1000 )
        
        self.__sigma = self.make_var ( sigma   ,
                                       'sigma_%s'        % name ,
                                       '#sigma_{CB}(%s)' % name , sigma , sigma , 0 , 1.e+6 ) 


        self.__right = True if right else False
        
        self.pdf = Ostap.Models.CutOffStudent (
            self.roo_name ( 'cofs_'  )          , 
            "Student's t-like cut-off %s" % self.name  ,
            self.xvar                          ,
            self.right                         , 
            self.x0                            ,
            self.nu                            ,
            self.sigma                         ) 
        
        ## save the configuration
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'right' : self.right , 
            'x0'    : self.x0    ,
            'nu'    : self.nu    ,
            'sigma' : self.sigma ,
            }
        
    @property
    def x0 ( self ) :
        """``x0'' : threshold/location parameter for Student's t-like cut-off"""
        return self.__x0
    @x0.setter 
    def x0 ( self , value ) :
        self.set_value ( self.__x0 , value ) 

    @property
    def sigma ( self ) :
        """``sigma'' : width parameter for Student's t-like cut-off"""
        return self.__sigma
    @sigma.setter 
    def sigma ( self , value ) :
        self.set_value ( self.__sigma , value ) 

    @property
    def nu ( self ) :
        """``nu'' : power parameter for Student's t-like cut-off"""
        return self.__nu
    @nu.setter 
    def nu ( self , value ) :
        self.set_value ( self.__nu , value ) 
    
    @property
    def right ( self ) :
        """``right'' : parameter of the Student's t-like cut-off"""
        return self.__right 

    @property
    def left ( self ) :
        """``left'' : parameter of the Student's t-like cut-off"""
        return not self.right 


models.append ( CutOffStudent_pdf )


# ==============================================================================
# @class CutOff_pdf
# Trivial wrapper for <code>Product1D_pdf</code>
class CutOff_pdf(Product1D_pdf) :
    """Trivial wrapper for `Product1D_pdf`
    """
    def __init__ ( self          ,
                   pdf           ,
                   cutoff        ,
                   xvar   = None , 
                   name   = ''   ) :

        Product1D_pdf.__init__ ( self          ,
                                 pdf1 = pdf    ,
                                 pdf2 = cutoff ,
                                 xvar = xvar   ,
                                 name = name   )

        ## save the configuration
        self.config = {
            'name'   : self.name     ,
            'pdf'    : self.orig_pdf ,
            'cutoff' : self.cutoff   ,
            'xvar'   : self.xvar     ,
            }
        
    @property
    def orig_pdf  ( self ) :
        """``orig_pdf'' : original PDF"""
        return self.pdf1

    @property
    def cutoff   ( self ) :
        """``cutoff'' : PDF ised for cut-off"""
        return self.pdf2
    
# =============================================================================
if '__main__' == __name__ : 

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )

# =============================================================================
##                                                                      The END 
# =============================================================================
