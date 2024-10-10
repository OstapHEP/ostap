#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/adjust.py
#  Simple "adjustment" of PDFs 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
# =============================================================================
"""Simple ``adjustment'' of PDFs"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'Adjust1D'     , ## helper class to adjust 1D-pdf 
    'Adjust2D'     , ## helper class to adjust 2D-pdf 
    'Adjust3D'     , ## helper class to adjust 3D-pdf
    ##
    'Adjust1D_pdf' , ## make adjusted 1D-pdf 
    'Adjust2D_pdf' , ## make adjusted 2D-pdf 
    'Adjust3D_pdf' , ## make adjusted 3D-pdf 
    )
# =============================================================================
from   ostap.fitting.fithelpers import VarMaker
from   ostap.fitting.pdfbasic   import ( PDF1 , Generic1D_pdf ,
                                         PDF2 , Generic2D_pdf ,
                                         PDF3 , Generic3D_pdf )
from   ostap.fitting.fit1d      import Flat1D 
from   ostap.fitting.fit2d      import Flat2D
from   ostap.fitting.fit3d      import Flat3D
import ROOT
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.adjust' )
else                       : logger = getLogger ( __name__               )
# =============================================================================
## @class Adjust1D
#  simple class to adjust certain PDF to avoid zeroes 
class Adjust1D(VarMaker) :
    """ Simple class to ``adjust'' certain PDF to avoid zeroes
    - a small flat component is added and the compound PDF is constructed
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,
                   xvar             , 
                   pdf              ,
                   fraction = 1.e-5 ) : 

        assert isinstance ( pdf  , ROOT.RooAbsPdf  ) , "``pdf''  must be ROOT.RooAbsPdf"
        assert isinstance ( xvar , ROOT.RooAbsReal ) , "``xvar'' must be ROOT.RooAbsReal"
        
        self.name      = name 
        self.__old_pdf = pdf
        
        from ostap.fitting.basic import Flat1D 
        self.__flat    = Flat1D  ( xvar  , name = 'flat_' + name ) 
        self.__frac    = self.make_var ( fraction , 'fracA_%s'                 % name ,
                                         'small fraction of flat component %s' % name ,
                                         fraction , 1.e-5 , 0 , 1 )
        
        self.__alist1 = ROOT.RooArgList ( self.__flat.pdf , self.__old_pdf )        
        self.__alist2 = ROOT.RooArgList ( self.__frac     )
        #
        ## final PDF
        # 
        self.__pdf    = ROOT.RooAddPdf  (
            self.roo_name ( "adjust_" ) ,            
            "Adjust %s" % self.name ,
            self.__alist1 ,
            self.__alist2 )
        
    @property
    def fraction( self ) :
        """``fraction''-parameter: the fraction of flat ``background'' added"""
        return  self.__frac
    @fraction.setter 
    def fraction( self , value ) :
        value = float ( value )
        assert 0 < value < 1 , 'Fraction  must be between 0 and 1'
        self.__frac.setVal ( value )
        
    @property
    def flat ( self ) :
        """new artificial ``flat'' component for the Pdf"""
        return self.__flat

    @property
    def pdf ( self ) :
        """``new'' (adjusted) Pdf"""
        return self.__pdf
    
    @property
    def old_pdf ( self ) :
        """``old'' (non-adjusted) Pdf"""
        return self.__old_pdf

# =============================================================================
## @class Adjust2D
#  simple class to adjust certain PDF to avoid zeroes 
class Adjust2D(VarMaker) :
    """ Simple class to ``adjust'' certain PDF to avoid zeroes
    - a small flat component is added and the compound PDF is constructed
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,
                   xvar             , 
                   yvar             , 
                   pdf              ,
                   fraction = 1.e-5 ) : 
        
        assert isinstance ( pdf  , ROOT.RooAbsPdf  ) , "``pdf''  must be ROOT.RooAbsPdf"
        assert isinstance ( xvar , ROOT.RooAbsReal ) , "``xvar'' must be ROOT.RooAbsReal"
        assert isinstance ( yvar , ROOT.RooAbsReal ) , "``yvar'' must be ROOT.RooAbsReal"
        
        self.name      = name 
        self.__old_pdf = pdf
        
        self.__flat    = Flat2D  ( xvar , yvar , name = 'flat_' + name )
        self.__frac    = self.make_var ( fraction , 'fracA_%s'                  % name ,
                                         'small  fraction of flat component %s' % name ,
                                         fraction , 1.e-5 , 0 , 1 )
        
        self.__alist1  = ROOT.RooArgList ( self.__flat.pdf , self.__old_pdf )        
        self.__alist2  = ROOT.RooArgList ( self.__frac     )
        #
        ## final PDF
        # 
        self.__pdf     = ROOT.RooAddPdf  (
            self.roo_name ( "adjust2_" ) , 
            "Adjust 2D %s" % self.name   ,
            self.__alist1 ,
            self.__alist2 )
        
    @property
    def fraction( self ) :
        """``fraction''-parameter: the fraction of flat ``background'' added"""
        return  self.__frac
    @fraction.setter 
    def fraction( self , value ) :
        value = float ( value )
        assert 0 < value < 1 , 'Fraction  must be between 0 and 1'
        self.__frac.setVal ( value )
        
    @property
    def flat ( self ) :
        """new artificial ``flat'' component for the Pdf"""
        return self.__flat

    @property
    def pdf ( self ) :
        """``new'' (adjusted) Pdf"""
        return self.__pdf
    
    @property
    def old_pdf ( self ) :
        """``old'' (non-adjusted) Pdf"""
        return self.__old_pdf
    
# =============================================================================
## @class Adjust3D
#  simple class to adjust certain PDF to avoid zeroes
class Adjust3D(VarMaker) :
    """ Simple class to ``adjust'' certain PDF to avoid zeroes
    - a small flat component is added and the compound PDF is constructed
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,
                   xvar             , 
                   yvar             , 
                   zvar             , 
                   pdf              ,
                   fraction = 1.e-5 ) : 
        
        assert isinstance ( pdf  , ROOT.RooAbsPdf  ) , "``pdf''  must be ROOT.RooAbsPdf"
        assert isinstance ( xvar , ROOT.RooAbsReal ) , "``xvar'' must be ROOT.RooAbsReal"
        assert isinstance ( yvar , ROOT.RooAbsReal ) , "``yvar'' must be ROOT.RooAbsReal"
        assert isinstance ( zvar , ROOT.RooAbsReal ) , "``zvar'' must be ROOT.RooAbsReal"
        
        self.name      = name 
        self.__old_pdf = pdf
        
        self.__flat    = Flat3D  ( xvar , yvar ,  xvar , name = 'flat_' + name )
        self.__frac    = self.make_var ( fraction , 'fracA_%s'                  % name ,
                                         'small  fraction of flat component %s' % name ,
                                         fraction , 1.e-5 , 0 , 1 )
        
        self.__alist1 = ROOT.RooArgList ( self.__flat.pdf , self.__old_pdf )        
        self.__alist2 = ROOT.RooArgList ( self.__frac     )
        #
        ## final PDF
        # 
        self.__pdf    = ROOT.RooAddPdf  (
            self.roo_name ( "adjust3_" ) , 
            "Adjust 3D %s" % self.name   ,
            self.__alist1 ,
            self.__alist2 )
        
    @property
    def fraction( self ) :
        """``fraction''-parameter: the fraction of flat ``background'' added"""
        return  self.__frac
    @fraction.setter 
    def fraction( self , value ) :
        value = float ( value )
        assert 0 < value < 1 , 'Fraction  must be between 0 and 1'
        self.__frac.setVal ( value )
        
    @property
    def flat ( self ) :
        """new artificial ``flat'' component for the Pdf"""
        return self.__flat

    @property
    def pdf ( self ) :
        """``new'' (adjusted) Pdf"""
        return self.__pdf
    
    @property
    def old_pdf ( self ) :
        """``old'' (non-adjusted) Pdf"""
        return self.__old_pdf

# =============================================================================
## @class Adjust1D_pdf
#  Create adjusted PDF
#  @code
#  pdf  = ....
#  apdf = Adjust1D_pdf ( pdf = pdf , xvar = xvar , fraction = 1.e-4 )
#  @endcode 
#  @see Adjust1D
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
class Adjust1D_pdf(PDF1) :
    """ Simple class to ``adjust'' certain 1D-PDF to avoid zeroes
    - a small flat component is added and the compound 1D-PDF is constructed
    
    >>> opdf = ...
    >>> apdf = Adjust1D_pdf ( pdf = opdf , fraction = 1e-5 )
    """
    ## constructor
    def __init__ ( self             ,
                   pdf              ,
                   xvar     = None  , 
                   fraction = 1.e-5 , 
                   name     = ''    ) :
                
        if isinstance    ( pdf , PDF ) :
            self.__old_pdf = pdf
            if xvar and not xvar is pdf.xvar : self.warning("mismatch in xvar?")
            xvar = pdf.xvar 
        elif insinstance ( pdf , ROOT.RooAbsPdf ) and xvar :
            self.__old_pdf = Generic1D_pdf ( pdf , xvar = xvar )
        else :
            raise TypeError ("Unknown type of pdf %s/%s"  % ( pdf , type ( pdf ) ) )

        assert isinstance ( xvar , ROOT.RooAbsReal ) , "``xvar'' must be ROOT.RooAbsReal"


        name = name if name else self.generate_name ( 'Adjust1D_' + self.old_pdf.name  ) 


        ## initialize the base
        PDF1.__init__ ( self , name = name , xvar = xvar )

        em = self.old_pdf.pdf.extendMode()
        if   1 == em : self.warning("PDF  ``canBeExtended''")
        elif 2 == em : self.warning("PDF ``mustBeExtended''")
        
        ## make the real adjustment 
        self.__adjustment = Adjust1D ( name                        ,
                                       pdf      = self.old_pdf.pdf ,
                                       xvar     = xvar             ,
                                       fraction = fraction         )
        
        self.pdf  = self.adjustment.pdf 
        self.pdf.fixCoefNormalization ( self.vars ) ## VB: added 10/10/2024 to suppress warnings 
        
        self.config = {
            'xvar'     : self.xvar    ,
            'name'     : self.name    ,
            'pdf'      : self.old_pdf ,
            'fraction' : self.fraction 
            }
        
    @property
    def adjustment ( self ) :
        """``adjustment'' : the adjustment object"""
        return self.__adjustment
    
    @property
    def fraction( self ) :
        """``fraction''-parameter: the fraction of flat ``background'' added"""
        return  self.adjustment.fraction
    @fraction.setter 
    def fraction( self , value ) :
        self.adjustment.fraction = value

    @property
    def old_pdf ( self ) :
        """``old'' (non-adjusted) PDF"""
        return self.__old_pdf

    ## redirect any other attributes to original PDF
    def __getattr__ ( self , attr ) :
        """Get all extra attributes from the original PDF"""
        opdf = self.old_pdf 
        return  getattr ( opdf , attr )

# =============================================================================
## @class Adjust2D_pdf
#  Create adjusted 2D-PDF
#  @code
#  pdf  = ....
#  apdf = Adjust2D_pdf ( pdf = pdf , xvar = xvar , yvar =  yvar , fraction = 1.e-4 )
#  @endcode 
#  @see Adjust2D
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
class Adjust2D_pdf(PDF2) :
    """ Simple class to ``adjust'' certain 1D-PDF to avoid zeroes
    - a small flat component is added and the compound 2D-PDF is constructed
    
    >>> opdf = ...
    >>> apdf = Adjust1D_pdf ( pdf = opdf , fraction = 1e-5 )
    """
    ## constructor
    def __init__ ( self             ,
                   pdf              ,
                   xvar     = None  , 
                   yvar     = None  , 
                   fraction = 1.e-5 , 
                   name     = ''    ) :
                
        if isinstance    ( pdf , PDF2 ) :
            self.__old_pdf = pdf 
            if xvar and not xvar is pdf.xvar : self.warning("mismatch in xvar?")
            if yvar and not yvar is pdf.yvar : self.warning("mismatch in yvar?")
            xvar = pdf.xvar 
            yvar = pdf.yvar 
        elif insinstance ( pdf , ROOT.RooAbsPdf ) and xvar and yvar :
            self.__old_pdf = Generic2D_pdf ( pdf , xvar = xvar , yvar = yvar)
        else :
            raise TypeError ("Unknown type of pdf %s/%s"  % ( pdf , type ( pdf ) ) )

        assert isinstance ( xvar , ROOT.RooAbsReal ) , "``xvar'' must be ROOT.RooAbsReal"
        assert isinstance ( yvar , ROOT.RooAbsReal ) , "``yvar'' must be ROOT.RooAbsReal"

        name = name if name else self.generate_name ( 'Adjust2D_' + self.old_pdf.name  ) 

        ## initialize the base
        PDF2.__init__ ( self , name = name , xvar = xvar , yvar = yvar )
        
        em = self.old_pdf.pdf.extendMode()
        if   1 == em : self.warning("PDF  ``canBeExtended''")
        elif 2 == em : self.warning("PDF ``mustBeExtended''")

        ## make the real adjustment 
        self.__adjustment = Adjust2D ( name                        ,
                                       pdf      = self.old_pdf.pdf ,
                                       xvar     = xvar             ,
                                       yvar     = yvar             ,
                                       fraction = fraction         )
        
        self.pdf  = self.adjustment.pdf 
        self.pdf.fixCoefNormalization ( self.vars ) ## VB: added 10/10/2024 to suppress warnings 
        
        self.config = {
            'xvar'     : self.xvar    ,
            'yvar'     : self.yvar    ,
            'name'     : self.name    ,
            'pdf'      : self.old_pdf ,
            'fraction' : self.fraction 
            }
        
    @property
    def adjustment ( self ) :
        """``adjustment'' : the adjustment object"""
        return self.__adjustment
    
    @property
    def fraction( self ) :
        """``fraction''-parameter: the fraction of flat ``background'' added"""
        return  self.adjustment.fraction
    @fraction.setter 
    def fraction( self , value ) :
        self.adjustment.fraction = value
        
    @property
    def old_pdf ( self ) :
        """``old'' (non-adjusted) PDF"""
        return self.__old_pdf

    ## redirect any other attributes to original PDF
    def __getattr__ ( self , attr ) :
        """Get all extra attributes from the original PDF"""
        opdf = self.old_pdf 
        return  getattr ( opdf , attr )

# =============================================================================
## @class Adjust3D_pdf
#  Create adjusted 3D-PDF
#  @code
#  pdf  = ....
#  apdf = Adjust3D_pdf ( pdf  = pdf  ,
#                        xvar = xvar ,
#                        yvar = yvar ,
#                        zvar = zvar ,
#                        fraction = 1.e-4 )
#  @endcode 
#  @see Adjust2D
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
class Adjust3D_pdf(PDF3) :
    """ Simple class to ``adjust'' certain 1D-PDF to avoid zeroes
    - a small flat component is added and the compound 2D-PDF is constructed
    
    >>> opdf = ...
    >>> apdf = Adjust3D_pdf ( pdf = opdf , fraction = 1e-5 )
    """
    ## constructor
    def __init__ ( self             ,
                   pdf              ,
                   xvar     = None  , 
                   yvar     = None  , 
                   zvar     = None  , 
                   fraction = 1.e-5 , 
                   name     = ''    ) :
                
        if isinstance    ( pdf , PDF3 ) :
            self.__old_pdf = pdf 
            if xvar and not xvar is pdf.xvar : self.warning("mismatch in xvar?")
            if yvar and not yvar is pdf.yvar : self.warning("mismatch in yvar?")
            if zvar and not zvar is pdf.zvar : self.warning("mismatch in yvar?")
            xvar = pdf.xvar 
            yvar = pdf.yvar 
            zvar = pdf.zvar 
        elif insinstance ( pdf , ROOT.RooAbsPdf ) and xvar and yvar and zvar :
            self.__old_pdf = Generic3D_pdf ( pdf , xvar = xvar , yvar = yvar , zvar = zvar )
        else :
            raise TypeError ("Unknown type of pdf %s/%s"  % ( pdf , type ( pdf ) ) )

        assert isinstance ( xvar , ROOT.RooAbsReal ) , "``xvar'' must be ROOT.RooAbsReal"
        assert isinstance ( yvar , ROOT.RooAbsReal ) , "``yvar'' must be ROOT.RooAbsReal"
        assert isinstance ( zvar , ROOT.RooAbsReal ) , "``zvar'' must be ROOT.RooAbsReal"

        name = name if name else self.generate_name ( 'Adjust3D_' + self.old_pdf.name  ) 

        ## initialize the base
        PDF3.__init__ ( self , name = name , xvar = xvar , yvar = yvar , zvar = zvar )
        
        em = self.old_pdf.pdf.extendMode()
        if   1 == em : self.warning("PDF  ``canBeExtended''")
        elif 2 == em : self.warning("PDF ``mustBeExtended''")

        ## make the real adjustment 
        self.__adjustment = Adjust3D ( name                        ,
                                       pdf      = self.old_pdf.pdf ,
                                       xvar     = xvar             ,
                                       yvar     = yvar             ,
                                       zvar     = yvar             ,
                                       fraction = fraction         )
        
        self.pdf  = self.adjustment.pdf 
        self.pdf.fixCoefNormalization ( self.vars ) ## VB: added 10/10/2024 to suppress warnings 
        
        self.config = {
            'xvar'     : self.xvar    ,
            'yvar'     : self.yvar    ,
            'zvar'     : self.zvar    ,
            'name'     : self.name    ,
            'pdf'      : self.old_pdf ,
            'fraction' : self.fraction 
            }
        
    @property
    def adjustment ( self ) :
        """``adjustment'' : the adjustment object"""
        return self.__adjustment
    
    @property
    def fraction( self ) :
        """``fraction''-parameter: the fraction of flat ``background'' added"""
        return  self.adjustment.fraction
    @fraction.setter 
    def fraction( self , value ) :
        self.adjustment.fraction = value
        
    @property
    def old_pdf ( self ) :
        """``old'' (non-adjusted) PDF"""
        return self.__old_pdf

    ## redirect any other attributes to original PDF
    def __getattr__ ( self , attr ) :
        """Get all extra attributes from the original PDF"""
        opdf = self.old_pdf 
        return  getattr ( opdf , attr )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
