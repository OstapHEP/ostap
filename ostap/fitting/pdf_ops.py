#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostap/fitting/pdf_ops.py
#  PDF multiplication  via operators: easy wrappers for class RooProdPdf
#  @see RooProdPdf
#  @see ostap.fitting.basic.PDF 
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
__all__     = ()
# =============================================================================
import ROOT
import ostap.fitting.basic 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.pdf_ops' )
else                       : logger = getLogger ( __name__                )
# =============================================================================        

# =============================================================================
def _prod_ ( pdf1 , pdf2 ) :
    return ROOT.RooProdPdf (
        'Product_%s_%s'    % ( pdf1.name ,  pdf2.name ) ,
        'Product:(%s)x(%s)'% ( pdf1.name ,  pdf2.name ) , pdf1.pdf , pdf2.pdf )

# =============================================================================
def _in_ ( a , *others ) :

    for b in others :
        if a is b : return True
        
    return False 
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
#  - PDF3 ( x , y , z ) * PDF  ( x )         -> PDF3 ( x , y , z )
#  - PDF3 ( x , y , z ) * PDF  ( y )         -> PDF3 ( x , y , z )
#  - PDF3 ( x , y , z ) * PDF  ( z )         -> PDF3 ( x , y , z ) 
#  - PDF2 ( x , y )     * PDF3 ( ... )       -> process as PDF3 (...) * PDF2 ( x , y )
#  - PDF2 ( x , y )     * PDF2 ( x , y )     -> PDF2 ( x , y )
#  - PDF2 ( x , y )     * PDF  ( x )         -> PDF2 ( x , y )
#  - PDF2 ( x , y )     * PDF  ( y )         -> PDF2 ( x , y )
#  - PDF  ( x )         * PDF3 ( ... )       -> process as PDF3 (...) * PDF ( x )
#  - PDF  ( x )         * PDF2 ( ... )       -> process as PDF2 (...) * PDF ( x )
#  - PDF  ( x )         * PDF  ( x )         -> PDF  ( x )
#  - PDF  ( x )         * PDF  ( y )         -> PDF2 ( x , y )
#  Other argument types and signatures are not supported
#  @see ostap.fitting.basic.PDF 
#  @see ostap.fitting.fit2d.PDF2
#  @see ostap.fitting.fit3d.PDF3
#  @see ostap.fitting.modifiers.Product1D            
#  @see ostap.fitting.fit2d.Model2D 
def pdf_product ( pdf1 , pdf2 ) :
    """ Product of two PDFs :
    - see ostap.fitting.basic.PDF 
    - see ostap.fitting.fit2d.PDF2
    - ostap.fitting.fit3d.PDF3
    - ostap.fitting.modifiers.Product1D            
    - ostap.fitting.fit2d.Model2D
    Supported argument  types and signatures:
    - PDF3 ( x , y , z ) * PDF3 ( x , y , z ) -> PDF3 ( x , y , z )
    - PDF3 ( x , y , z ) * PDF2 ( x , y )     -> PDF3 ( x , y , z )
    - PDF3 ( x , y , z ) * PDF2 ( x , z )     -> PDF3 ( x , y , z )
    - PDF3 ( x , y , z ) * PDF2 ( y , z )     -> PDF3 ( x , y , z )
    - PDF3 ( x , y , z ) * PDF  ( x )         -> PDF3 ( x , y , z )
    - PDF3 ( x , y , z ) * PDF  ( y )         -> PDF3 ( x , y , z )
    - PDF3 ( x , y , z ) * PDF  ( z )         -> PDF3 ( x , y , z ) 
    - PDF2 ( x , y )     * PDF3 ( ... )       -> process as PDF3 (...) * PDF2 ( x , y )
    - PDF2 ( x , y )     * PDF2 ( x , y )     -> PDF2 ( x , y )
    - PDF2 ( x , y )     * PDF  ( x )         -> PDF2 ( x , y )
    - PDF2 ( x , y )     * PDF  ( y )         -> PDF2 ( x , y )
    - PDF  ( x )         * PDF3 ( ... )       -> process as PDF3 (...) * PDF ( x )
    - PDF  ( x )         * PDF2 ( ... )       -> process as PDF2 (...) * PDF ( x )
    - PDF  ( x )         * PDF  ( x )         -> PDF  ( x )
    - PDF  ( x )         * PDF  ( y )         -> PDF2 ( x , y )
    
    >>> pdf1 = ...
    >>> pdf2 = ...
    >>> pdf = pdf1 * pdf2 
    """

    import ostap.fitting.basic as _1D 
    import ostap.fitting.fit2d as _2D 
    import ostap.fitting.fit3d as _3D 
    
    ## 1D * ...
    if   isinstance ( pdf1  , _3D.PDF3 ) :
        
        x1 = pdf1.xvar 
        y1 = pdf1.yvar 
        z1 = pdf1.zvar
        v1 = x1 , y1 , z1
        
        if   isinstance  ( pdf2 , _3D.PDF3 )   :

            x2 = pdf2.xvar 
            y2 = pdf2.yvar 
            z2 = pdf2.zvar 

            if _in_ ( x2 , *v1 ) and _in_ ( y2 , *v1 ) and _in_ ( z2 , *v1 ) :
                return _3D.Generic3D_pdf ( _prod_ ( pdf1 , pdf2 ) , *v1 )
            
        elif isinstance  ( pdf2 , _2D.PDF2 )   :
            
            x2 = pdf2.xvar 
            y2 = pdf2.yvar 

            if _in_ ( x2 , *v1 ) and _in_ ( y2 , *v1 ) :
                return _3D.Generic3D_pdf ( _prod_ ( pdf1 , pdf2 ) , *v1 )

        elif isinstance  ( pdf2 , _1D.PDF  )   :
            
            x2 = pdf2.xvar 

            if _in_ ( x2 , *v1 ) : 
                return _3D.Generic3D_pdf ( _prod_ ( pdf1 , pdf2 ) , *v1  )

        return NotImplemented 
        
    elif isinstance ( pdf1  , _2D.PDF2 ) :

        x1 = pdf1.xvar 
        y1 = pdf1.yvar 
        v1 = x1 , y1
        
        if   isinstance ( pdf2 , _3D.PDF3 ) : return pdf_product ( pdf2 , pdf1 )
        elif isinstance ( pdf2 , _2D.PDF2 ) :

            x2 = pdf2.xvar 
            y2 = pdf2.yvar 

            if   _in_ ( x2 , *v1 ) and _in_ ( y2 , *v1 ) : 
                return _2D.Generic2D_pdf ( _prod_ ( pdf1 , pdf2 ) , *v1 )
            elif x1 is x2 : 
                return _3D.Generic3D_pdf ( _prod_ ( pdf1 , pdf2 ) , x1 , y1 , y2 )
            elif x1 is y2 : 
                return _3D.Generic3D_pdf ( _prod_ ( pdf1 , pdf2 ) , x1 , y1 , x2 )
            elif y1 is x2 : 
                return _3D.Generic3D_pdf ( _prod_ ( pdf1 , pdf2 ) , x1 , y1 , y2 )
            elif y1 is y2 : 
                return _3D.Generic3D_pdf ( _prod_ ( pdf1 , pdf2 ) , x1 , y1 , x2 )
                
        elif isinstance ( pdf2 , _1D.PDF  ) :
            
            x2 = pdf2.yvar
            
            if  _in_ ( x2 , *v1 ) :
                return _2D.Generic2D_pdf ( _prod_ ( pdf1 , pdf2 ) , *v1 )
            
            return     _3D.Generic3D_pdf ( _prod_ ( pdf1 , pdf2 ) , x1 , y1 , x2 ) 

        return NotImplemented
    
    elif isinstance ( pdf1  , _1D.PDF  ) :

        if   isinstance ( pdf2  , _3D.PDF3 ) : return pdf_product ( pdf2 , pdf1 ) 
        elif isinstance ( pdf2  , _2D.PDF2 ) : return pdf_product ( pdf2 , pdf1 ) 
        elif isinstance ( pdf2  , _1D.PDF  ) :
            
            x1 = pdf1.xvar
            x2 = pdf2.xvar
            
            if x1 is x2 :
                from   ostap.fitting.modifiers import Product1D_pdf            
                return Product1D_pdf (        pdf1 , pdf2 , x1      )
            else        :
                return _2D.Model2D   ( ''   , pdf1 , pdf2 , x1 , x2 ) 

    return NotImplemented

# =============================================================================
## add multiplication operator for PDF 
ostap.fitting.basic.PDF . __mul__  = pdf_product

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
