#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostap/fitting/pdf_ops.py
#  Introduce helper operators for PDFs
#   - PDF multiplication via operators: easy wrappers for class RooProdPdf
#   - non-extended sum of two PDFs 
#   - non-extended sum of two PDFs 

#  @see RooProdPdf
#  @see ostap.fitting.basic.PDF1
#  @see ostap.fitting.fit2d.PDF2
#  @see ostap.fitting.fit3d.PDF3
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2019-01-28
# =============================================================================
"""Add PDF multiplication  via operators
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2019-01-28"
__all__     = (
    #
    'Prod1D_pdf'  , ## helper class to implement product of two RAW PDFs 
    'Prod2D_pdf'  , ## helper class to implement product of two RAW PDFs 
    'Prod3D_pdf'  , ## helper class to implement product of two RAW PDFs
    #
    'pdf_product' , ## helper fnuction to create a product of PDFs  
    'pdf_sum'     , ## helper fnuction to create a non-extended sum of PDFs  
    )
# =============================================================================
from   ostap.core.ostap_types import sequence_types, sized_types 
from   ostap.fitting.funbasic import constant_types 
from   ostap.fitting.pdfbasic import APDF1, PDF1, PDF2, PDF3
import ROOT
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.pdf_ops' )
else                       : logger = getLogger ( __name__                )
# ========================================================================
## @class Prod1D_pdf
#  Simple product of two pdfs
#  - trivial wrapper for RooProdPdf
#  @code
#  pdf1     = ...
#  pdf2     = ...
#  pdf_prod = Prod1D_pdf( [ pdf1 , pdf2 ] , xvar = ...  ) 
#  @endcode
#  @see RooProdPdf
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2018-11-29  
class Prod1D_pdf(PDF1) :
    """Simple product of 1D-PDFs
    - actually it is a trivial wrapper for RooProdPdf    
    >>> pdf1 = ...
    >>> pdf2 = ...
    >>> pdf  = Prod1D_pdf( [ pdf1 , pdf2 ] , xvar = ... )    
    """
    def __init__ ( self       ,
                   pdfs       ,
                   xvar       ,
                   name  = '' ) :

        assert xvar and isinstance ( xvar , ROOT.RooAbsReal ) , \
               "Invalid type of 'xvar' %s/%s" % ( xvar , type ( xvar ) )
        assert 2 <= len ( pdfs ) , "There must be at least two components!"

        ## keep the argument list locally 
        self.__pdfs = tuple ( pdfs )
        
        ## generic name 
        patname =  '*'.join ( '(%s)' % p.name for p in self.pdfs )
        ## check the instance name 
        name    = name if name else self.new_name ( patname ) 
        
        ## initialize the base class
        PDF1.__init__ ( self , name , xvar = xvar )
        
        self.__pdfs = tuple ( pdfs )
        
        ## the actual product of PDFs
        self.pdf = self.raw_product ( *self.pdfs )
        
        ## save configuration for cloning/pickling 
        self.config = {
            'pdfs'    : self.pdfs ,
            'xvar'    : self.xvar ,
            'name'    : self.name ,
            }
        
    @property
    def pdfs ( self ) :
        """'pdfs' : the list of PDF/RooAbsPdf-like objects"""
        return self.__pdfs

    @property
    def pdf1 ( self ) :
        """'pdf1' : the first object"""
        return self.__pdfs[0]
    
    @property
    def pdf2 ( self ) :
        """'pdf2' : the second object"""
        return self.__pdfs[1]

    @property
    def tail ( self ) :
        """'tail' : get 'other' PDFs"""
        return self.__pdfs[2:]

# =============================================================================
## @class Prod2D_pdf
#  Simple product of two pdfs
#  - trivial wrapper for RooProdPdf
#  @code
#  pdf1     = ...
#  pdf2     = ...
#  pdf_prod = Prod2D_pdf( [ pdf1 , pdf2 ] , xvar = ... ) 
#  @endcode
#  @see RooProdPdf
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2018-11-29  
class Prod2D_pdf(PDF2) :
    """Simple product of 1D-PDFs
    - actually it is a trivial wrapper for RooProdPdf
    
    >>> pdf1 = ...
    >>> pdf2 = ...
    >>> pdf  = Prod2D_pdf( [ pdf1 , pdf2 ) ] , xvar = ... ) 
    """
    def __init__ ( self       ,
                   pdfs       ,
                   xvar       ,
                   yvar       ,
                   name  = '' ) : 
        
        assert xvar and isinstance ( xvar , ROOT.RooAbsReal ) , \
               "Invalid type of 'xvar' %s/%s" % ( xvar , type ( xvar ) )
        assert yvar and isinstance ( yvar , ROOT.RooAbsReal ) , \
               "Invalid type of 'yvar' %s/%s" % ( yvar , type ( xvar ) )
        assert 2 <= len ( pdfs ) , "There must be at least two components!"
        
        self.__pdfs = tuple ( pdfs )

        ## generic name 
        patname =  '*'.join ( '(%s)' % p.name for p in self.pdfs )
        ## check the instance name 
        name    = name if name else self.new_name ( patname ) 
        
        ## initialize the base class
        PDF2.__init__ ( self , name , xvar = xvar , yvar = yvar )
        
        ## the actual product of PDFs
        self.pdf = self.raw_product ( *self.pdfs )
        
        ## save configuration for cloning/pickling 
        self.config = {
            'pdfs'    : self.pdfs    ,
            'xvar'    : self.xvar    ,
            'yvar'    : self.yvar    ,
            'name'    : self.name    ,
            }
        
    @property
    def pdfs ( self ) :
        """'pdfs' : the list of PDF/RooAbsPdf-like objects"""
        return self.__pdfs

    @property
    def pdf1 ( self ) :
        """'pdf1' : the first object"""
        return self.__pdfs[0]
    
    @property
    def pdf2 ( self ) :
        """'pdf2' : the second object"""
        return self.__pdfs[1]

    @property
    def tail ( self ) :
        """'tail' : get 'other' PDFs"""
        return self.__pdfs[2:]

# =============================================================================
## @class Prod3D_pdf
#  Simple product of PDFs
#  - trivial wrapper for RooProdPdf
#  @code
#  pdf1     = ...
#  pdf2     = ...
#  pdf_prod = Prod3D_pdf( [ pdf1 , pdf2 ] , xvar = ... , yvar = ... , zvar = ... ) ) 
#  @endcode
#  @see RooProdPdf
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2018-11-29  
class Prod3D_pdf(PDF3) :
    """Simple product of 1D-PDFs
    - actually it is a trivial wrapper for RooProdPdf
    
    >>> pdf1 = ...
    >>> pdf2 = ...
    >>> pdf  = Prod3D_pdf( [ pdf1 , pdf2 ) ] , xvar = ... , yvar = ... , zvar = ... )  
    """
    def __init__ ( self       ,
                   pdfs       ,
                   xvar       ,
                   yvar       ,
                   zvar       ,
                   name  = '' ) :
    
        assert xvar and isinstance ( xvar , ROOT.RooAbsReal ) , \
               "Invalid type of 'xvar' %s/%s" % ( xvar , type ( xvar ) )
        assert yvar and isinstance ( yvar , ROOT.RooAbsReal ) , \
               "Invalid type of 'yvar' %s/%s" % ( yvar , type ( xvar ) )
        assert zvar and isinstance ( zvar , ROOT.RooAbsReal ) , \
               "Invalid type of 'zvar' %s/%s" % ( zvar , type ( xvar ) )
        
        assert 2 <= len ( pdfs ) , "There must be at least two components!"
        
        self.__pdfs = tuple ( pdfs )

        ## generic name 
        patname =  '*'.join ( '(%s)' % p.name for p in self.pdfs )
        ## check the instance name 
        name    = name if name else self.new_name ( patname ) 
        
        ## initialize the base class
        PDF3.__init__ ( self , name , xvar = xvar , yvar = yvar , zvar = zvar )
        
        ## the actual product of PDFs
        self.pdf = self.raw_product ( *self.pdfs )
        
        ## save configuration for cloning/pickling 
        self.config = {
            'pdfs'    : self.pdfs ,
            'xvar'    : self.xvar ,
            'yvar'    : self.yvar ,
            'zvar'    : self.zvar ,
            'name'    : self.name ,
            }
        
    @property
    def pdfs ( self ) :
        """'pdfs' : the list of PDF/RooAbsPdf-like objects"""
        return self.__pdfs

    @property
    def pdf1 ( self ) :
        """'pdf1' : the first object"""
        return self.__pdfs[0]
    
    @property
    def pdf2 ( self ) :
        """'pdf2' : the second object"""
        return self.__pdfs[1]

    @property
    def tail ( self ) :
        """'tail' : get 'other' PDFs"""
        return self.__pdfs[2:]

# =============================================================================
## helper functon to make a raw product of PDFs or RooAbsPDF objects
def raw_product ( self , *pdfs ) :
    """Make a raw product of PDFs or RooAbsPDF objects
    """
    lpdfs = [] 
    for i , p in enumerate ( pdfs ) :
        if   p and isinstance ( p , APDF1          ) : lpdfs.append ( p.pdf )
        elif p and isinstance ( p , ROOT.RooAbsPdf ) : lpdfs.append ( p     )
        else : raise TypeError ( "Invalid type for %s component %s/%s" % ( i , p , type ( p ) ) ) 
        
    assert 2 <= len ( lpdfs ) , 'raw_product: there should be at leats two elements in the PDF list!'
    
    name  = self.new_roo_name ( '*'.join ( '(%s)' % p.name for p in pdfs ) ) 
    title = 'product: ' +     ( '*'.join ( '(%s)' % p.name for p in pdfs ) )
    
    self.aux_keep.append ( lpdfs ) 
    if 2 == len ( lpdfs ) : return ROOT.RooProdPdf ( name , title , *lpdfs )
    
    plst = ROOT.RooArgList()
    for p in lpdfs : plst.add ( p )
    
    self.aux_keep.append ( plst ) 
    return ROOT.RooProdPdf ( name , title , plst  )
# =============================================================================

Prod1D_pdf.raw_product = raw_product
Prod2D_pdf.raw_product = raw_product
Prod3D_pdf.raw_product = raw_product

# ============================================================================
## Product of two PDFs :
#  @code
#  pdf1 = ...
#  pdf2 = ...
#  pdf = pdf1 * pdf2
#  @endcode
#  Supported argument types and signatures 
#  - PDF1 ( x )         * PDF3 ( ... )       -> process as PDF3 (...) * PDF ( x )
#  - PDF1 ( x )         * PDF2 ( ... )       -> process as PDF2 (...) * PDF ( x )
#  - PDF1 ( x )         * PDF1 ( x )         -> PDF1 ( x )
#  - PDF1 ( x )         * PDF1 ( y )         -> PDF2 ( x , y )
#
#  Other argument types and signatures are not supported
#
def pdf1_product ( pdf1 , pdf2 ) :
    """ Product of two PDFs :
    >>> pdf1 = ...
    >>> pdf2 = ...
    >>> pdf = pdf1 * pdf2    
    Supported argument types and signatures 
    - PDF1 ( x )         * PDF3 ( ... )       -> process as PDF3 (...) * PDF ( x )
    - PDF1 ( x )         * PDF2 ( ... )       -> process as PDF2 (...) * PDF ( x )
    - PDF1 ( x )         * PDF1 ( x )         -> PDF1 ( x )
    - PDF1 ( x )         * PDF1 ( y )         -> PDF2 ( x , y )
    
    Other argument types and signatures are not supported
    """
    if   isinstance ( pdf2 , PDF3 ) : return pdf3_product ( pdf2 , pdf1 )
    elif isinstance ( pdf2 , PDF2 ) : return pdf2_product ( pdf2 , pdf1 )
    elif isinstance ( pdf2 , PDF1 ) :
        
        if pdf2.xvar in pdf1.vars : return Prod1D_pdf ( ( pdf1 , pdf2 ) , xvar = pdf1.xvar )
        else                      : return Prod2D_pdf ( ( pdf1 , pdf2 ) , xvar = pdf1.xvar , yvar = pdf2.xvar )
        
    return NotImplemented 

# =============================================================================
## Product of two PDFs :
#  @code
#  pdf1 = ...
#  pdf2 = ...
#  pdf = pdf1 * pdf2
#  @endcode
#  Supported argument types and signatures 
#  - PDF2 ( x , y )     * PDF3 ( ... )       -> process as PDF3 (...) * PDF2 ( x , y )
#  - PDF2 ( x , y )     * PDF2 ( x , y )     -> PDF2 ( x , y )
#  - PDF2 ( x , y )     * PDF1 ( x )         -> PDF2 ( x , y )
#  - PDF2 ( x , y )     * PDF1 ( y )         -> PDF2 ( x , y )
def pdf2_product ( pdf1 , pdf2 ) :
    """ Product of two PDFs :
    Supported argument  types and signatures:
    - PDF2 ( x , y )     * PDF3 ( ... )       -> process as PDF3 (...) * PDF2 ( x , y )
    - PDF2 ( x , y )     * PDF2 ( x , y )     -> PDF2 ( x , y )
    - PDF2 ( x , y )     * PDF1 ( x )         -> PDF2 ( x , y )
    - PDF2 ( x , y )     * PDF1 ( y )         -> PDF2 ( x , y )
    
    >>> pdf1 = ...
    >>> pdf2 = ...
    >>> pdf = pdf1 * pdf2 
    """
    
    if   isinstance ( pdf2 , PDF3 ) : return pdf3_product ( pdf2 , pdf1 )
    elif isinstance ( pdf2 , PDF2 ) : 

        if   pdf2.xvar in pdf1.vars and pdf2.yvar in pdf1.vars :
            return Prod2D_pdf ( ( pdf1 , pdf2 ) , xvar = pdf1.xvar , yvar = pdf1.yvar )
        elif pdf2.xvar in pdf1.vars :
            return Prod3D_pdf ( ( pdf1 , pdf2 ) , xvar = pdf1.xvar , yvar = pdf1.yvar , zvar = pdf2.yvar )
        elif pdf2.yvar in pdf1.vars :
            return Prod3D_pdf ( ( pdf1 , pdf2 ) , xvar = pdf1.xvar , yvar = pdf1.yvar , zvar = pdf2.xvar )

        return NotImplemented 
            
    elif isinstance ( pdf2 , PDF1 ) :

        if   pdf2.xvar in pdf1.vars : return Prod2D_pdf ( ( pdf1 , pdf2 ) , xvar = pdf1.xvar , yvar = pdf1.yvar )
        else                        : return Prod3D_pdf ( ( pdf1 , pdf2 ) , xvar = pdf1.xvar , yvar = pdf1.yvar, zvar = pdf2.xvar  )
        
    return NotImplemented

# =============================================================================
## Product of two PDFs :
#  @code
#  pdf1 = ...
#  pdf2 = ...
#  pdf = pdf1 * pdf2
#  @endcode
# 
#  Supported argument types and signatures 
#  - PDF3 ( x , y , z ) * PDF3 ( x , y , z ) -> PDF3 ( x , y , z )
#  - PDF3 ( x , y , z ) * PDF2 ( x , y )     -> PDF3 ( x , y , z )
#  - PDF3 ( x , y , z ) * PDF2 ( x , z )     -> PDF3 ( x , y , z )
#  - PDF3 ( x , y , z ) * PDF2 ( y , z )     -> PDF3 ( x , y , z )
#  - PDF3 ( x , y , z ) * PDF1 ( x )         -> PDF3 ( x , y , z )
#  - PDF3 ( x , y , z ) * PDF1 ( y )         -> PDF3 ( x , y , z )
#  - PDF3 ( x , y , z ) * PDF1 ( z )         -> PDF3 ( x , y , z ) 
#  Other argument types and signatures are not supported
def pdf3_product ( pdf1 , pdf2 ) :
    """ Product of two PDFs :
    >>> pdf1 = ...
    >>> pdf2 = ...
    >>> pdf = pdf1 * pdf2 
    Supported argument  types and signatures:
    - PDF3 ( x , y , z ) * PDF3 ( x , y , z ) -> PDF3 ( x , y , z )
    - PDF3 ( x , y , z ) * PDF2 ( x , y )     -> PDF3 ( x , y , z )
    - PDF3 ( x , y , z ) * PDF2 ( x , z )     -> PDF3 ( x , y , z )
    - PDF3 ( x , y , z ) * PDF2 ( y , z )     -> PDF3 ( x , y , z )
    - PDF3 ( x , y , z ) * PDF1 ( x )         -> PDF3 ( x , y , z )
    - PDF3 ( x , y , z ) * PDF1 ( y )         -> PDF3 ( x , y , z )
    - PDF3 ( x , y , z ) * PDF1 ( z )         -> PDF3 ( x , y , z ) 
    
    """

    if   isinstance ( pdf2 , PDF3 )   \
           and pdf2.xvar in self.pdf1 \
           and pdf2.yvar in self.pdf1 \
           and pdf3.zvar in self.pdf1 :        
        return Prod3D_pdf ( ( pdf1 , pdf2 ) , xvar = pdf1.xvar , yvar = pdf1.yvar, zvar = pdf1.zvar  )
    
    elif isinstance ( pdf2 , PDF2 )     \
             and pdf2.xvar in pdf1.vars \
             and pdf2.yars in pdf1.vars :        
        return Prod3D_pdf ( ( pdf1 , pdf2 ) , xvar = pdf1.xvar , yvar = pdf1.yvar, zvar = pdf1.zvar  )
    
    elif isinstance ( pdf2 , PDF1 )       \
             and pdf2.xvar in pdf1.vars  :        
        return Prod3D_pdf ( ( pdf1 , pdf2 ) , xvar = pdf1.xvar , yvar = pdf1.yvar, zvar = pdf1.zvar  )

    return NotImplemented 

# =============================================================================
## Product of two PDFs :
#  @code
#  pdf1 = ...
#  pdf2 = ...
#  pdf = pdf1 * pdf2
#  @endcode
# 
#  Supported argument types and signatures 
#  - PDF3 ( x , y , z ) * PDF3 ( x , y , z ) -> PDF3 ( x , y , z )
#  - PDF3 ( x , y , z ) * PDF2 ( x , y )     -> PDF3 ( x , y , z )
#  - PDF3 ( x , y , z ) * PDF2 ( x , z )     -> PDF3 ( x , y , z )
#  - PDF3 ( x , y , z ) * PDF2 ( y , z )     -> PDF3 ( x , y , z )
#  - PDF3 ( x , y , z ) * PDF1 ( x )         -> PDF3 ( x , y , z )
#  - PDF3 ( x , y , z ) * PDF1 ( y )         -> PDF3 ( x , y , z )
#  - PDF3 ( x , y , z ) * PDF1 ( z )         -> PDF3 ( x , y , z ) 
#  - PDF2 ( x , y )     * PDF3 ( ... )       -> process as PDF3 (...) * PDF2 ( x , y )
#  - PDF2 ( x , y )     * PDF2 ( x , y )     -> PDF2 ( x , y )
#  - PDF2 ( x , y )     * PDF1 ( x )         -> PDF2 ( x , y )
#  - PDF2 ( x , y )     * PDF1 ( y )         -> PDF2 ( x , y )
#  - PDF1 ( x )         * PDF3 ( ... )       -> process as PDF3 (...) * PDF ( x )
#  - PDF1 ( x )         * PDF2 ( ... )       -> process as PDF2 (...) * PDF ( x )
#  - PDF1 ( x )         * PDF1 ( x )         -> PDF1 ( x )
#  - PDF1 ( x )         * PDF1 ( y )         -> PDF2 ( x , y )
#  Other argument types and signatures are not supported
def pdf_product ( pdf1 , pdf2 ) :
    """ Product of two PDFs :
    Supported argument  types and signatures:
    - PDF3 ( x , y , z ) * PDF3 ( x , y , z ) -> PDF3 ( x , y , z )
    - PDF3 ( x , y , z ) * PDF2 ( x , y )     -> PDF3 ( x , y , z )
    - PDF3 ( x , y , z ) * PDF2 ( x , z )     -> PDF3 ( x , y , z )
    - PDF3 ( x , y , z ) * PDF2 ( y , z )     -> PDF3 ( x , y , z )
    - PDF3 ( x , y , z ) * PDF1 ( x )         -> PDF3 ( x , y , z )
    - PDF3 ( x , y , z ) * PDF1 ( y )         -> PDF3 ( x , y , z )
    - PDF3 ( x , y , z ) * PDF1 ( z )         -> PDF3 ( x , y , z ) 
    - PDF2 ( x , y )     * PDF3 ( ... )       -> process as PDF3 (...) * PDF2 ( x , y )
    - PDF2 ( x , y )     * PDF2 ( x , y )     -> PDF2 ( x , y )
    - PDF2 ( x , y )     * PDF1 ( x )         -> PDF2 ( x , y )
    - PDF2 ( x , y )     * PDF1 ( y )         -> PDF2 ( x , y )
    - PDF1 ( x )         * PDF3 ( ... )       -> process as PDF3 (...) * PDF ( x )
    - PDF1 ( x )         * PDF2 ( ... )       -> process as PDF2 (...) * PDF ( x )
    - PDF1 ( x )         * PDF1 ( x )         -> PDF1 ( x )
    - PDF1 ( x )         * PDF1 ( y )         -> PDF2 ( x , y )
    
    >>> pdf1 = ...
    >>> pdf2 = ...
    >>> pdf = pdf1 * pdf2 
    """

    result = NotImplemented 
    if   isinstane  ( pdf1 , PDF3 ) : result = pdf3_product ( pdf1 , pdf2 )
    elif isinstance ( pdf1 , PDF2 ) : result = pdf2_product ( pdf1 , pdf2 )
    elif isinstance ( pdf1 , PDF1 ) : result = pdf1_product ( pdf1 , pdf2 )
    
    if result is NotImplemted and isinstance ( pdf2 , PDF3 ) : result = pdf3_product ( pdf2 , pdf2 )
    if result is NotImplemted and isinstance ( pdf2 , PDF2 ) : result = pdf2_product ( pdf2 , pdf2 )
    if result is NotImplemted and isinstance ( pdf2 , PDF1 ) : result = pdf1_product ( pdf2 , pdf2 )
    
    if result is NotImplemented :
        raise NotImplementedError ( "Product of %s and %s is nudefined!" % ( pdf1 , pdf2 ) )

    return result 

# =============================================================================
## Make an non-extended  sum of the 1D PDFs
#  @see Sum1D
def pdf1_sum ( pdf1 , pdf2 , *other ) :
    """Make an non-extended  sum of the 1D PDFs
    - see Sum1D
    """
    if   isinstance ( pdf2 , PDF1 ) : pass 
    elif isinstance ( pdf2 , sequence_types ) :
        args = tuple ( pdf2 ) + other 
        return pdf1_sum ( pdf1 , *args )

    pall = ( pdf1 , pdf2 , ) + other

    head = pall [0]
    tail = pall [1:]
    if not isinstance ( head , PDF1 )  : return NotImplemented
    
    for p in tail :    
        if not isinstance ( p , PDF1 ) : return NotImplemented
        if not p.xvar in head.vars     : return NotImplemented 

    from ostap.fitting.fit1d import Sum1D 
    return Sum1D ( pall )
    
# =============================================================================
## Make an non-extended  sum of the 2D PDFs
#  @see Sum2D
def pdf2_sum ( pdf1 , pdf2 , *other ) :
    """Make an non-extended  sum of the 2D PDFs
    - see Sum2D
    """
    if   isinstance ( pdf2 , PDF2 ) : pass 
    elif isinstance ( pdf2 , sequence_types ) :
        args = tuple ( pdf2 ) + other 
        return pdf2_sum ( pdf1 , *args )

    pall = ( pdf1 , pdf2 , ) + other
    
    head = pall [0]
    tail = pall [1:]
    if not isinstance ( head , PDF2 )  : return NotImplemented
    
    for p in tail :    
        if not isinstance ( p , PDF2 ) : return NotImplemented
        if not p.xvar in head.vars     : return NotImplemented 
        if not p.yvar in head.vars     : return NotImplemented 
        
    from ostap.fitting.fit2d import Sum2D
    return Sum2D ( pall , xvar = head.xvar , yvar = head.yvar ) 
    
# =============================================================================
## Make an non-extended  sum of the 3D PDFs
#  @see Sum3D
def pdf3_sum ( pdf1 , pdf2 , *other ) :
    """Make an non-extended  sum of the 3D PDFs
    - see Sum3D
    """
    
    if   isinstance ( pdf2 , PDF3 ) : pass 
    elif isinstance ( pdf2 , sequence_types ) :
        args = tuple ( pdf2 ) + other 
        return pdf3_sum ( pdf1 , *args )

    pall = ( pdf1 , pdf2 , ) + other

    head = pall [0]
    tail = pall [1:]
    if not isinstance ( head , PDF3 )  : return NotImplemented
    
    for p in tail :    
        if not isinstance ( p , PDF3 ) : return NotImplemented
        if not p.xvar in head.vars     : return NotImplemented 
        if not p.yvar in head.vars     : return NotImplemented 
        if not p.zvar in head.vars     : return NotImplemented 

    from ostap.fitting.fit3d import Sum3D
    return Sum3D ( pall , xvar = head.xvar , yvar = head.yvar , zvar = head.zvar ) 

     
# =============================================================================
## Non-extended sum of two PDFs
#  @code
#  pdf1 = ...
#  pdf2 = ...
#  pdf = pdf1 + pdf2
#  @endcode
#  @see Sum1D
#  @see Sum2D
#  @see Sum3D
def pdf_sum ( pdf1 , pdf2 , *other ) :
    """ Non-extended sum of two PDFs
    >>> pdf1 = ...
    >>> pdf2 = ...
    >>> pdf = pdf1 + pdf2
    """
    result = NotImplemented 
    
    if   isinstance ( pdf1 , PDF3 ) : result = pdf3_sum ( pdf1 , pdf2 , *other )
    elif isinstance ( pdf1 , PDF2 ) : result = pdf2_sum ( pdf1 , pdf2 , *other )
    elif isinstance ( pdf1 , PDF1 ) : result = pdf1_sum ( pdf1 , pdf2 , *other )
    
    if result is NotImplemented :
        raise NotImplementedError ( "Sum of  %s, %s and %s  is undefined!" % ( pdf1 , pdf2 , list ( others ) ) ) 
    
    return result 


# =============================================================================
## make a convolution (FFT) for a given PDF
#  @code
#  pdf        = ...
#  resolution = ...
#  result1    = pdf % resolution
#  result2    = pdf % Convolution ( ... ) 
#  result3    = pdf % ( 10 * MeV )
#  result4    = pdf % (  5 * MeV , 15 * MeV ) 
#  result5    = pdf % ( 10 * MeV , 5 * MeV , 15 * MeV )
#  ## only for Python 3
#  result1    = pdf @ resolution
#  result2    = pdf @ Convolution ( ... ) 
#  result3    = pdf @ ( 10 * MeV )
#  result4    = pdf @ (  5 * MeV , 15 * MeV ) 
#  result5    = pdf @ ( 10 * MeV , 5 * MeV , 15 * MeV ) 
#  @endcode
#  The resolution can be :
#   - PDF
#   - ostap.fitting.convolution.Convolution 
#   - ROOT.RooAbsPdf
#   - ROOT.RooAbsReal (Gaussian PDF will be constructed with sigma = resoltuion )
#   - float           (Gaussian PDF will be constructed with sigma = resoltuion )
#   - 2/3 tuple       (Gaussian PDF will be constructed with sigma = resoltuion )
#  @see ostap.fitting.convolution.Convolution
#  @see ostap.fitting.convolution.Convolution_pdf 
def pdf_convolution ( pdf , resolution ) :
    """Make a convolution (FFT) for a given PDF
    The resolution can be :
    - PDF
    - ostap.fitting.convolution.Convolution
    - ROOT.RooAbsPdf
    - ROOT.RooAbsReal (Gaussian PDF will be constructed with sigma = resolution )
    - float           (Gaussian PDF will be constructed with sigma = resolution )
    - 2/3 tuple       (Gaussian PDF will be constructed with sigma = resolution )
    
    >>> pdf        = ...
    >>> resolution = ...

    ## python 2 and 3 
    >>> result1    = pdf % resolution
    >>> result2    = pdf % Convolution ( ... ) 
    >>> result3    = pdf % ( 10 * MeV )                      ## Gaussian
    >>> result4    = pdf % (  5 * MeV , 15 * MeV )           ## Gaussian
    >>> result5    = pdf % ( 10 * MeV , 5 * MeV , 15 * MeV ) ## Gaussian

    ## only python 3 
    >>> result1    = pdf @ resolution
    >>> result2    = pdf @ Convolution ( ... ) 
    >>> result3    = pdf @ ( 10 * MeV )                      ## Gaussian
    >>> result4    = pdf @ (  5 * MeV , 15 * MeV )           ## Gaussian
    >>> result5    = pdf @ ( 10 * MeV , 5 * MeV , 15 * MeV ) ## Gaussian

    - see ostap.fitting.convolution.Convolution
    - see ostap.fitting.convolution.Convolution_pdf 
    """
    import ostap.fitting.convolution as CNV

    if not isinstance ( pdf , PDF1 ) : return NotImplemented


    ## Ready to use convolution  PDF 
    if   isinstance ( resolution , PDF1 ) and resolution.xvar in pdf.vars : 
        ## treat it as convolution
        return CNV.Convolution_pdf ( pdf , resolution , **CNV.CnvConfig.config() )
    ## Resolution PDF 
    elif isinstance ( resolution , ROOT.RooAbsPdf  ) : pass
    ## Gaussian sigma 
    elif isinstance ( resolution , ROOT.RooAbsReal ) : pass
    ## Gaussian sigma 
    elif isinstance ( resolution , constant_types  ) and 0 < float ( resolution ) : pass 
    ## Gaussian sigma 
    elif isinstance ( resolution , sized_types     ) and 2 <= len  ( resolution ) <= 4 : pass

    else : return NotImplemented 
            
    
    ## create convolution object 
    cnv = CNV.Convolution ( pdf        = pdf        ,
                            resolution = resolution ,
                            xvar       = pdf.xvar   , **CNV.CnvConfig.config() )
    
    return CNV.Convolution_pdf ( pdf , cnv )
    

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
