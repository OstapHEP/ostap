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
    'Adjust1D_pdf' , ## make adjusted 1D-pdf 
    'Adjust2D_pdf' , ## make adjusted 2D-pdf 
    'Adjust3D_pdf' , ## make adjusted 3D-pdf 
    )
# =============================================================================
from   ostap.core.ostap_types   import num_types 
from   ostap.fitting.pdfbasic   import ( Sum1D , Flat1D , 
                                         Sum2D , Flat2D , 
                                         Sum3D , Flat3D ) 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.adjust' )
else                       : logger = getLogger ( __name__               )
# =============================================================================
## @class Adjust1D_pdf
#  Create adjusted PDF
#  @code
#  pdf  = ....
#  apdf = Adjust1D_pdf ( pdf = pdf , xvar = xvar , fraction = 1.e-4 )
#  @endcode 
#  @see Adjust1D
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
class Adjust1D_pdf(Sum1D) :
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
                   name     = ''    ,
                   prefix   = 'f'   , ## prefix for fraction names 
                   suffix   = ''    , ## suffix for fraction names 
                   fix_norm = False ) :

        pdf1 , xvar    = self.make_PDF1 ( pdf , xvar  )
        flat           = Flat1D  ( xvar = xvar )
        Sum1D.__init__ ( self , [ flat , pdf ] ,
                         name     = name      ,
                         fraction = fraction  ,  
                         prefix   = prefix    ,
                         suffix   = suffix    ,
                         fix_norm = fix_norm  )

        self.F.fix()
        
        self.config = {
            'xvar'     : self.xvar     ,
            'name'     : self.name     ,
            'pdf'      : self.old_pdf  ,
            'fraction' : self.fraction ,
            'prefix'   : self.prefix   ,  
            'suffix'   : self.suffix   , 
            'fix_norm' : self.fix_norm 
        }

    @property
    def old_pdf ( self ) :
        """``old'' (non-adjusted) PDF"""
        return self.pdf2
    
    @property
    def flat    ( self ) :
        """`flat` : get flat component"""
        return self.pdf1

    ## redirect any other attributes to original PDF
    def __getattr__ ( self , attr ) :
        """ Get all extra attributes from the original PDF"""
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
class Adjust2D_pdf(Sum2D) :
    """ Simple class to ``adjust'' certain 1D-PDF to avoid zeroes
    - a small flat component is added and the compound 2D-PDF is constructed
    
    >>> opdf = ...
    >>> apdf = Adjust2D_pdf ( pdf = opdf , fraction = 1e-5 )
    """
    ## constructor
    def __init__ ( self             ,
                   pdf              ,
                   xvar     = None  , 
                   yvar     = None  , 
                   fraction = 1.e-5 , 
                   name     = ''    ,
                   prefix   = 'f'   , ## prefix for fraction names 
                   suffix   = ''    , ## suffix for fraction names 
                   fix_norm = False ) :
        
        pdf1 , xvar , yvar = self.make_PDF2 ( pdf , xvar , yvar )
        flat           = Flat2D  ( xvar = xvar , yvar = yvar )
        Sum2D.__init__ ( self , [ flat , pdf ] ,
                         name     = name       ,
                         fraction = fraction   ,  
                         prefix   = prefix     ,
                         suffix   = suffix     ,
                         fix_norm = fix_norm   )

        self.F.fix()
        
        self.config = {
            'xvar'     : self.xvar     ,
            'yvar'     : self.yvar     ,
            'name'     : self.name     ,
            'pdf'      : self.old_pdf  ,
            'fraction' : self.fraction ,
            'prefix'   : self.prefix   ,  
            'suffix'   : self.suffix   , 
            'fix_norm' : self.fix_norm 
        }
        
    @property
    def flat    ( self ) :
        """`flat` : get flat component"""
        return self.pdf2
    
    @property
    def old_pdf ( self ) :
        """``old'' (non-adjusted) PDF"""
        return self.pdf1
    
    ## redirect any other attributes to original PDF
    def __getattr__ ( self , attr ) :
        """ Get all extra attributes from the original PDF"""
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
class Adjust3D_pdf(Sum3D) :
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
                   name     = ''    ,
                   prefix   = 'f'   , ## prefix for fraction names 
                   suffix   = ''    , ## suffix for fraction names 
                   fix_norm = False ) :
        
        pdf1 , xvar , yvar , zvar = self.make_PDF3 ( pdf , xvar , yvar , zvar )
        flat                      = Flat3D  ( xvar = xvar , yvar = yvar , zvar = zvar )
        Sum3D.__init__ ( self , [ flat , pdf ] ,
                         name     = name       ,
                         fraction = fraction   ,  
                         prefix   = prefix     ,
                         suffix   = suffix     ,
                         fix_norm = fix_norm   )
        
        self.F.fix()
        
        self.config = {
            'xvar'     : self.xvar     ,
            'yvar'     : self.yvar     ,
            'zvar'     : self.zvar     ,
            'name'     : self.name     ,
            'pdf'      : self.old_pdf  ,
            'fraction' : self.fraction ,
            'prefix'   : self.prefix   ,  
            'suffix'   : self.suffix   , 
            'fix_norm' : self.fix_norm 
        }
        
    @property
    def flat    ( self ) :
        """`flat` : get flat component"""
        return self.pdf1
    
    @property
    def old_pdf ( self ) :
        """``old'' (non-adjusted) PDF"""
        return self.pdf2
    
    ## redirect any other attributes to original PDF
    def __getattr__ ( self , attr ) :
        """ Get all extra attributes from the original PDF"""
        opdf = self.old_pdf 
        return  getattr ( opdf , attr )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
