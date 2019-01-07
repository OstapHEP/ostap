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
    'Product1D_pdf' , ## make a product of 1D PDFs
    'Modify1D_pdf'  , ## modify 1D pdf with a positive polynomials 
    )
# =============================================================================
import ROOT, math
from   ostap.fitting.basic import PDF , Generic1D_pdf
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.modifiers' )
else                       : logger = getLogger ( __name__                  )
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
                self.error ("Mismatch in ``xvar''-observable") 
            elif not xvar : xvar = pdf1.xvar 
        elif isinstance ( pdf1 , ROOT.RooAbsPdf ) and xvar :
            self.__pdf1 = Generic1D_pdf ( pdf1 , xvar )
        else :
            raise TypeError("Illegal setting for ``pdf1'': %s/%s" % ( pdf1 , type ( pdf1 ) ) )

        assert isinstance ( self.__pdf1  , PDF ), 'Invalid pdf1 type'

        if isinstance ( pdf2 , PDF ) :
            self.__pdf2 = pdf2
            if xvar and not ( xvar is pdf2.xvar ) :
                self.error ("Mismatch in ``xvar''-observable") 
            elif not xvar : xvar = pdf2.xvar 
        elif isinstance ( pdf2 , ROOT.RooAbsPdf ) and xvar :
            self.__pdf2 = Generic1D_pdf ( pdf2 , xvar )
        else :
            raise TypeError("Illegal setting for ``pdf2'': %s/%s" % ( pdf2 , type ( pdf1 ) ) )
        
        assert isinstance ( self.__pdf2  , PDF ), 'Invalid pdf2 type'
        
        assert isinstance ( xvar , ROOT.RooAbsReal ),\
               "Invalid ``xvar'':%s/%s" % ( xvar , type  ( xvar ) ) 

        
        if not name  : name  = "product_%s_%s"  % ( self.pdf1.name , self.pdf2.name )
        if not title : title = "Product(%s,%s)" % ( self.pdf1.name , self.pdf2.name )

        ## initialize the base class
        PDF.__init__ ( self , name , xvar =  xvar )

        em1 = self.pdf1.pdf.extendMode()
        em2 = self.pdf2.pdf.extendMode()
        
        if   2 == em1 : self.warning ( "pdf1 ``must-be-extended''" )
        elif 1 == em1 : self.warning ( "pdf1  ``can-be-extended''" )
        
        if   2 == em2 : self.warning ( "pdf2 ``must-be-extended''" )
        elif 1 == em2 : self.warning ( "pdf2  ``can-be-extended''" )

        ## finally build PDF 
        self.pdf = ROOT.RooProdPdf ( name , title , self.pdf1.pdf , self.pdf2.pdf )

        ## save configuration for cloning
        self.config = {
            'pdf1'  : self.pdf1  ,
            'pdf2'  : self.pdf2  ,
            'xvar'  : self.xvar  ,
            'name'  : self.name  ,
            'title' : self.title ,            
            }
        print  'I AM PRODUCT'
        print  'PDF1:', self.pdf1
        print  'PDF2:', self.pdf2
        
    @property
    def pdf1 ( self ) :
        """``pdf1'' : the first PDF"""
        print 'I am property pdf1'
        return self.__pdf1
    
    @property
    def pdf2 ( self ) :
        """``pdf2'' : the second PDF"""
        return self.__pdf2
    
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

        if not name  : name  = "modify_%s_%s"  % ( pdf.name , power )
        if not title : title = "Modify(%s,%s)" % ( pdf.name , power )

        from ostap.fitting.background import PolyPos_pdf
        pdf2 = PolyPos_pdf( 'M_%s_%s' % ( pdf.name , power ) ,
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

# =============================================================================
if '__main__' == __name__ : 

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )

# =============================================================================
# The END 
# =============================================================================
