#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file models_3d.py
#
#  Smooth non-factorizable 3D-models to describe background distribtions
#
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2017-11-20
# =============================================================================
""" Set of useful non-factorizable 3D-models to describe background distribtions
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    'PolyPos3D_pdf'   , ## A positive polynomial in 3D  
    'PolyPos3Dsym_pdf', ## A positive symmetric polynomial in 3D
    'PolyPos3Dmix_pdf', ## A positive partly symmetric (x<-->y) polynomial in 3D
    )
# =============================================================================
import ROOT, math
from   ostap.core.core     import cpp, Ostap
from   ostap.math.base     import iszero
from   ostap.fitting.basic import makeVar
from   ostap.fitting.fit3d import PDF3 
# =============================================================================
from   ostap.logger.logger     import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.models_3d' )
else                       : logger = getLogger ( __name__                  )
# =============================================================================
models = []
# =============================================================================
## @class PolyPos3D_pdf
#  The 3D-polynomial of order Nx*Ny*Nz, that is constrained 
#  to be non-negative over the  defined range      
#  \f[  P(x,y,z) = \sum_{i,j,k} a_{ijk}B^{n_x}_i(x) B^{n_y}_j(y) B^{n_z}_k(z)\f] 
#  where all coefficients \f$a_{ijk}\f$ are non-negative and 
#  \f$ \sum_{i,j,k} a_{ijk}=1 \f$ 
#  @author Vanya BELYAEV Ivan.Belayev@itep.ru
#  @date 2017-11-14
#  @see Ostap::Models::Poly3DPositive
#  @see Ostap::Math::Positive3D
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class PolyPos3D_pdf(PDF3) :
    """Positive (non-factorizable!) polynomial in 3D:
    
    The 3D-polynomial of order Nx*Ny*Nz, that is constrained 
    to be non-negative over the defined range
    
    - P(x,y,z) = \sum_{i,j,k} a_{ijk}B^{n_x}_i(x) B^{n_y}_j(y) B^{n_z}_k(z)
    
    where  B^n_i - are Bernstein polynomials and all coefficients \a_{ijk} are:
    - non-negative: 0<=a_{ijk}
    
    """
    def __init__ ( self   ,
                   name   ,
                   x      ,   ## the first  dimension  
                   y      ,   ## the second dimension
                   z      ,   ## the third  dimension
                   nx = 1 ,   ## polynomial degree in X 
                   ny = 1 ,   ## polynomial degree in Y
                   nz = 1 ) : ## polynomial degree in Z

        ##   inialize the base :
        PDF3.__init__ ( self , name , x , y , z ) 
        #
        if 0 > nx : raise ValueError('PolyPos3D_pdf: Invalid nx=%s' % nx )
        if 0 > ny : raise ValueError('PolyPos3D_pdf: Invalid ny=%s' % ny )
        if 0 > nz : raise ValueError('PolyPos3D_pdf: Invalid nz=%s' % nz )

        #
        ## create parameters
        #
        num = ( nx + 1 ) * ( ny + 1 ) *  ( nz + 1 ) 
        self.makePhis ( num - 1 ) 
            
        #
        ## finally build PDF 
        #
        self.pdf = cpp.Ostap.Models.Poly3DPositive (
            'p3Dp_%s'            % name ,
            'Poly3DPositive(%s)' % name ,
            self.x        ,
            self.y        ,
            self.z        ,
            nx            ,
            ny            , 
            nz            , 
            self.phi_list )
        
models.append ( PolyPos3D_pdf ) 

# =============================================================================
## @class PolyPos3Dsym_pdf
#  The 3D-polynomial of order N*N*N, that is constrained 
#  to be non-negative ans symmetric over the  defined range      
#   \f[  P(x,y,z) = \sum_{i,j,k} a_{ijk}B^{n}_i(x) B^{n}_j(y) B^{n}_k(z)\f] 
#  where all coefficients \f$a_{ijk}\f$ are:
# - non-negative: \f$ a_{ijk}\ge0 \f$
# - symmetric: \f$ a_{ijk}=a_{jik}=a_{ikj}\f$
# - constrainted: \f$ \sum_{i,j,k} a_{ijk}=1 \f$ 
#  @author Vanya BELYAEV Ivan.Belayev@itep.ru
#  @date 2017-11-14
#  @see Analysis::Models::Poly3DSymPositive
#  @see Gaudi::Math::Positive3DSym
class PolyPos3Dsym_pdf(PDF3) :
    """Positive (non-factorizable!) symmetric polynomial in 3D:
    
    The 3D-polynomial of order N*N*N, that is constrained 
    to be non-negative ans symmetric over the  defined range
    
    - P(x,y,z) = \sum_{i,j,k} a_{ijk}B^{n}_i(x) B^{n}_j(y) B^{n}_k(z)
    
    where  B^n_i - are Bernstein polynomials and all coefficients \a_{ijk} are:
    - non-negative:   0<= a_{ijk}
    - symmetric:      a_{ijk}=a_{jik}=a_{ikj}\f$
        
    
    where all coefficients a_{ijk} are non-negative: 0<=a_{ijk}
    
    - 0<=P(x,y,x) for whole 3D-range
    """
    def __init__ ( self   ,
                   name   ,
                   x      ,   ## the first  dimension  
                   y      ,   ## the second dimension
                   z      ,   ## the third  dimension
                   n  = 1 ) : ## polynomial degree in X,Y,Z

        ##   inialize the base :
        PDF3.__init__ ( self , name , x , y , z ) 
        #
        if 0 > n : raise ValueError('PolyPos3Dsym_pdf: Invalid n=%s' % n )

        if self.x.getMin() != self.y.getMin() or self.y.getMin() != self.z.getMin() or \
           self.x.getMax() != self.y.getMax() or self.y.getMax() != self.z.getMax() :
            logger.warning("PolyPos3Dsym_pdf: non-equal ranges for 3D-symmetric model")  
            
        #
        ## create parameters
        #
        num = ( n + 1 ) * ( n + 2 ) *  ( n + 3 ) / 6
        self.makePhis ( num - 1 ) 
        #
            
        #
        ## finally build PDF 
        #
        self.pdf = cpp.Ostap.Models.Poly3DSymPositive (
            'p3Ds_%s'               % name ,
            'Poly3DSymPositive(%s)' % name ,
            self.x        ,
            self.y        ,
            self.z        ,
            n             ,
            self.phi_list )
        
models.append ( PolyPos3Dsym_pdf ) 


# =============================================================================
## @class PolyPos3Dmix_pdf
#  The 3D-polynomial of order N*N*Nz, that is constrained 
#  to be non-negative and x<-->y symmetric over the  defined range      
#   \f[  P(x,y,z) = \sum_{i,j,k} a_{ijk}B^{n}_i(x) B^{n}_j(y) B^{n}_k(z)\f] 
#  where all coefficients \f$a_{ijk}\f$ are:
# - non-negative: \f$ a_{ijk}\ge0 \f$
# - patly symmetric: \f$ a_{ijk}=a_{jik}\f$
# - constrainted: \f$ \sum_{i,j,k} a_{ijk}=1 \f$ 
#  @author Vanya BELYAEV Ivan.Belayev@itep.ru
#  @date 2017-11-14
#  @see Analysis::Models::Poly3DMixPositive
#  @see Gaudi::Math::Positive3DMix
class PolyPos3Dmix_pdf(PDF3) :
    """Positive (non-factorizable!)  x<-->y symmetric polynomial in 3D:
    
    The 3D-polynomial of order N*N*N, that is constrained 
    to be non-negative ans symmetric over the  defined range
    
    - P(x,y,z) = \sum_{i,j,k} a_{ijk}B^{n}_i(x) B^{n}_j(y) B^{n_z}_k(z)
    
    where  B^n_i - are Bernstein polynomials and all coefficients \a_{ijk} are:
    - non-negative:   0<= a_{ijk}
    - partly symmetric:      a_{ijk}=a_{jik}
        
    """
    def __init__ ( self   ,
                   name   ,
                   x      ,   ## the first  dimension  
                   y      ,   ## the second dimension
                   z      ,   ## the third  dimension
                   n  = 1 ,   ## polynomial degree in X,Y
                   nz = 1 ) : ## polynomial degree in Z

        ##   inialize the base :
        PDF3.__init__ ( self , name , x , y , z ) 
        #
        if 0 > n  : raise ValueError('PolyPos3Dmix_pdf: Invalid n=%s'  % n )
        if 0 > nz : raise ValueError('PolyPos3Dmix_pdf: Invalid nz=%s' % nz )
        
        if self.x.getMin() != self.y.getMin() or self.x.getMax() != self.y.getMax() :
            logger.warning("PolyPos3Dmix_pdf: non-equal ranges for 3D-(x<-->y) symmetric model")  
            
        #
        ## create parameters
        #
        num = ( n + 1 ) * ( n + 2 ) *  ( nz + 1 ) / 2 
        self.makePhis ( num - 1 ) 
        #
            
        #
        ## finally build PDF 
        #
        self.pdf = cpp.Ostap.Models.Poly3DMixPositive (
            'p3Dm_%s'               % name ,
            'Poly3DMixPositive(%s)' % name ,
            self.x        ,
            self.y        ,
            self.z        ,
            n             ,
            nz            ,
            self.phi_list )
        
models.append ( PolyPos3Dmix_pdf ) 

           
# =============================================================================
if '__main__' == __name__ : 
         
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )
    
# =============================================================================
# The END 
# =============================================================================
