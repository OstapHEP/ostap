#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostap/fitting/pdf_ops.py
#  Introduce helper operators for PDFs
#   - PDF multiplication via operators: easy wrappers for class RooProdPdf
#   - non-extended sum of two PDFs 
#   - non-extended sum of two PDFs 

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
        
        
        if   isinstance  ( pdf2 , _3D.PDF3 )   :

            x2 , y2 , z2 = pdf2.xvar , pdf2.yvar , pdf2.zvar 
            if x2 in pdf1.vars and y2 in pdf1.vars and z2 in pdf1.vars :
                return _3D.Generic3D_pdf ( _prod_ ( pdf1 , pdf2 ) , pdf1.xvar , pdf1.yvar , pdf1.zvar )
            
        elif isinstance  ( pdf2 , _2D.PDF2 )   :
            
            x2 , y2 = pdf2.xvar , pdf2.yvar 
            if x2 in pdf1.vars and y2 in pdf1.vars:
                return _3D.Generic3D_pdf ( _prod_ ( pdf1 , pdf2 ) , pdf1.xvar , pdf1.yvar , pdf1.zvar )

        elif isinstance  ( pdf2 , _1D.PDF  )   :
            
            x2 = pdf2.xvar 
            if x2 in pdf1.vars :
                return _3D.Generic3D_pdf ( _prod_ ( pdf1 , pdf2 ) , pdf1.xvar , pdf1.yvar , pdf1.zvar )

        return NotImplemented 
        
    elif isinstance ( pdf1  , _2D.PDF2 ) :

        x1 , y1 = pdf1.xvar , pdf1.yvar 
        
        if   isinstance ( pdf2 , _3D.PDF3 ) : return pdf_product ( pdf2 , pdf1 )
        elif isinstance ( pdf2 , _2D.PDF2 ) :

            x2 , y2  = pdf2.xvar ,  pdf2.yvar 

            if   x2 in pdf1.vars and y2 in pdf1.vars :  
                return _2D.Generic2D_pdf ( _prod_ ( pdf1 , pdf2 ) , pdf1.xvar , pdf1.yvar )
            elif x1 is x2 : 
                return _3D.Generic3D_pdf ( _prod_ ( pdf1 , pdf2 ) , x1 , y1 , y2 )
            elif x1 is y2 : 
                return _3D.Generic3D_pdf ( _prod_ ( pdf1 , pdf2 ) , x1 , y1 , x2 )
            elif y1 is x2 : 
                return _3D.Generic3D_pdf ( _prod_ ( pdf1 , pdf2 ) , x1 , y1 , y2 )
            elif y1 is y2 : 
                return _3D.Generic3D_pdf ( _prod_ ( pdf1 , pdf2 ) , x1 , y1 , x2 )
                
        elif isinstance ( pdf2 , _1D.PDF  ) :
            
            x2 = pdf2.xvar
            
            if  x2 in pdf1.vars :
                return _2D.Generic2D_pdf ( _prod_ ( pdf1 , pdf2 ) , pdf1.xvar , pdf2.yvar )
            
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
## Non-extended sum of two PDFs
#  @code
#  pdf1 = ...
#  pdf2 = ...
#  pdf = pdf1 + pdf2
#  @endcode
#  @see ostap.fitting.basic.Sum1D
#  @see ostap.fitting.fit2d.Sum2D
#  @see ostap.fitting.fit3d.Sum2D
def pdf_sum ( pdf1 , pdf2 ) :
    """ Non-extended sum of two PDFs
    - see ostap.fitting.basic.Sum1D
    - see ostap.fitting.fit2d.Sum2D
    - see ostap.fitting.fit3d.Sum2D
    >>> pdf1 = ...
    >>> pdf2 = ...
    >>> pdf = pdf1 + pdf2
    """
    
    import ostap.fitting.basic as _1D 
    import ostap.fitting.fit2d as _2D 
    import ostap.fitting.fit3d as _3D 
    
    if   isinstance ( pdf1 , _3D.PDF3 ) and isinstance ( pdf2 , _3D.PDF3 ) :
        
        if not pdf1.xvar in pdf2.vars : return NotImplemented
        if not pdf1.yvar in pdf2.vars : return NotImplemented
        if not pdf1.zvar in pdf2.vars : return NotImplemented
        
        return _3D.Sum3D ( pdf1 , pdf2 )
    
    elif isinstance ( pdf1 , _3D.PDF3 ) or  isinstance ( pdf2 , _3D.PDF3 ) : return NotImplemented
    
    elif isinstance ( pdf1 , _2D.PDF2 ) and isinstance ( pdf2 , _2D.PDF2 ) :

        if not pdf1.xvar in pdf2.vars : return NotImplemented
        if not pdf1.yvar in pdf2.vars : return NotImplemented
        
        return _2D.Sum2D ( pdf1 , pdf2 )
            
    elif isinstance ( pdf1 , _2D.PDF2 ) or  isinstance ( pdf2 , _2D.PDF2 ) : return NotImplemented
    
    elif isinstance ( pdf1 , _1D.PDF  ) and isinstance ( pdf2 , _1D.PDF  ) :
        
        if not pdf1.xvar in pdf2.vars : return NotImplemented
        
        return _1D.Sum1D ( pdf1 , pdf2 )

    return NotImplemented 

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

    import ostap.fitting.basic       as      _1D
    import ostap.fitting.convolution as     _CNV
    from   ostap.core.ostap_types          import num_types
    from   ostap.core.core           import VE 
    
    if not isinstance ( pdf , _1D.PDF ) : return NotImplemented

    if   isinstance ( resolution , _1D.PDF          ) and pdf.xvar is resolution.xvar :
        return _CNV.Convolution_pdf ( pdf , resolution )
    elif isinstance ( resolution , _CNV.Convolution ) and pdf.xvar is resolution.xvar :
        return _CNV.Convolution_pdf ( pdf , resolution )
    elif isinstance ( resolution , ROOT.RooAbsPdf   ) :
        return _CNV.Convolution_pdf ( pdf , resolution )
    elif isinstance ( resolution , ROOT.RooAbsReal   ) :
        return _CNV.Convolution_pdf ( pdf , resolution )
    elif isinstance ( resolution , num_types         ) :
        return _CNV.Convolution_pdf ( pdf , resolution )
    elif isinstance ( resolution , VE                ) :
        return _CNV.Convolution_pdf ( pdf , resolution )
    elif isinstance ( resolution , tuple             ) :
        return _CNV.Convolution_pdf ( pdf , resolution )

    return NotImplemented

    
    
        
# =============================================================================
## add multiplication operator for PDFs 
ostap.fitting.basic.PDF . __mul__     = pdf_product
## add addition       operator for PDFs 
ostap.fitting.basic.PDF . __add__     = pdf_sum 
## add convolution operator for PDFs 
ostap.fitting.basic.PDF . __mod__     = pdf_convolution
## Python3: add convolution operator for PDFs 
ostap.fitting.basic.PDF . __matmul__  = pdf_convolution

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
