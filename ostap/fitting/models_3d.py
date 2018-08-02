#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/models_3d.py
#  Smooth non-factorizable 3D-models to describe background distribtions
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
from   ostap.fitting.fit3d import PDF3 
from   ostap.fitting.utils import Phases
# =============================================================================
from   ostap.logger.logger     import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.models_3d' )
else                       : logger = getLogger ( __name__                  )
# =============================================================================
models = []
# =============================================================================
##  @class PolyBase3
#   helper base class to implement various polynomial-like shapes
class PolyBase3(PDF3,Phases) :
    """Helper base class to implement various polynomial-like shapes
    """
    def __init__ ( self , name , xvar , yvar , zvar , power , the_phis = None ) :
        PDF3  .__init__ ( self , name  , xvar , yvar , zvar )
        Phases.__init__ ( self , power , the_phis  )
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
class PolyPos3D_pdf(PolyBase3) :
    """Positive (non-factorizable!) polynomial in 3D:
    
    The 3D-polynomial of order Nx*Ny*Nz, that is constrained 
    to be non-negative over the defined range
    
    - P(x,y,z) = \sum_{i,j,k} a_{ijk}B^{n_x}_i(x) B^{n_y}_j(y) B^{n_z}_k(z)
    
    where  B^n_i - are Bernstein polynomials and all coefficients \a_{ijk} are:
    - non-negative: 0<=a_{ijk}
    
    """
    def __init__ ( self            ,
                   name            ,
                   xvar            ,   ## the first  dimension  
                   yvar            ,   ## the second dimension
                   zvar            ,   ## the third  dimension
                   nx = 1          ,   ## polynomial degree in X 
                   ny = 1          ,   ## polynomial degree in Y
                   nz = 1          ,   ## polynomial degree in Z
                   the_phis = None ) : 
        
        ## check arguments 
        assert isinstance ( nx , int ) and 0 <= nx < 50 , "``nx''-parameter is illegal: %s" % nx 
        assert isinstance ( ny , int ) and 0 <= ny < 50 , "``ny''-parameter is illegal: %s" % ny
        assert isinstance ( nz , int ) and 0 <= nz < 50 , "``nz''-parameter is illegal: %s" % nz
        
        ##   inialize the base :
        PolyBase3.__init__ ( self , name , xvar , yvar , zvar ,
                             ( nx + 1 ) * ( ny + 1 ) *  ( nz + 1 ) - 1 , the_phis )
        #
        self.__nx = nx
        self.__ny = ny
        self.__nz = nz
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.Poly3DPositive (
            'p3Dp_%s'            % name ,
            'Poly3DPositive(%s)' % name ,
            self.x        ,
            self.y        ,
            self.z        ,
            self.nx       ,
            self.ny       , 
            self.nz       , 
            self.phi_list )
        
        ## save  the configurtaion
        self.config = {
            'name'     : self.name ,
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'zvar'     : self.zvar ,
            'nx'       : self.nx   ,
            'ny'       : self.ny   ,
            'nz'       : self.nz   ,
            'the_phis' : self.phis ,
            }

    @property
    def nx ( self ) :
        """``nx''-parameter - order/degree of 3D-polynom in x-direction"""
        return self.__nx
    @property
    def ny ( self ) :
        """``ny''-parameter - order/degree of 3D-polynom in z-direction"""
        return self.__ny
    @property
    def nz ( self ) :
        """``nz''-parameter - order/degree of 3D-polynom in z-direction"""
        return self.__nz
    
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
#  @see Ostap::Models::Poly3DSymPositive
#  @see Ostap::Math::Positive3DSym
class PolyPos3Dsym_pdf(PolyBase3) :
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
    def __init__ ( self            ,
                   name            ,
                   xvar            ,   ## the first  dimension  
                   yvar            ,   ## the second dimension
                   zvar            ,   ## the third  dimension
                   n  = 1          ,  ## polynomial degree in X,Y,Z
                   the_phis = None ) :
        
        ## check arguments 
        assert isinstance ( n , int ) and 0 <= n < 50 , "``n''-parameter is illegal: %s" % n 
        
        ##   inialize the base :
        PolyBase3.__init__ ( self , name , xvar , yvar , zvar ,
                             ( n + 1 ) * ( n + 2 ) *  ( n + 3 ) / 6 - 1 ,   the_phis )
        
        if self.xminmax() != self.yminmax() : 
            logger.warning ( 'PolyPos3Dsym: x&y have different edges %s vs %s' % ( self.xminmax() , self.yminmax() ) )
        if self.xminmax() != self.zminmax() : 
            logger.warning ( 'PolyPos3Dsym: x&z have different edges %s vs %s' % ( self.xminmax() , self.zminmax() ) )
        if self.yminmax() != self.zminmax() : 
            logger.warning ( 'PolyPos3Dsym: y&z have different edges %s vs %s' % ( self.yminmax() , self.zminmax() ) )
            
        self.__n = n 
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.Poly3DSymPositive (
            'p3Ds_%s'               % name ,
            'Poly3DSymPositive(%s)' % name ,
            self.x        ,
            self.y        ,
            self.z        ,
            self.n        ,
            self.phi_list )
        
        ## save  the configurtaion
        self.config = {
            'name'     : self.name ,
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'zvar'     : self.zvar ,
            'n'        : self.n    ,
            'the_phis' : self.phis ,
            }
        
    @property
    def n  ( self ) :
        """``n''-parameter - order/degree of 3D-polynom in x,y&z-directions"""
        return self.__n
    @property
    def nx ( self ) :
        """``nx''-parameter - order/degree of 3D-polynom in x-direction"""
        return self.__n
    @property
    def ny ( self ) :
        """``ny''-parameter - order/degree of 3D-polynom in z-direction"""
        return self.__n
    @property
    def nz ( self ) :
        """``nz''-parameter - order/degree of 3D-polynom in z-direction"""
        return self.__n
 
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
#  @see Ostap::Models::Poly3DMixPositive
#  @see Ostap::Math::Positive3DMix
class PolyPos3Dmix_pdf(PolyBase3) :
    """Positive (non-factorizable!)  x<-->y symmetric polynomial in 3D:
    
    The 3D-polynomial of order N*N*N, that is constrained 
    to be non-negative ans symmetric over the  defined range
    
    - P(x,y,z) = \sum_{i,j,k} a_{ijk}B^{n}_i(x) B^{n}_j(y) B^{n_z}_k(z)
    
    where  B^n_i - are Bernstein polynomials and all coefficients \a_{ijk} are:
    - non-negative:   0<= a_{ijk}
    - partly symmetric:      a_{ijk}=a_{jik}
        
    """
    def __init__ ( self            ,
                   name            ,
                   xvar            ,   ## the first  dimension  
                   yvar            ,   ## the second dimension
                   zvar            ,   ## the third  dimension
                   n  = 1          ,   ## polynomial degree in X,Y
                   nz = 1          , ## polynomial degree in Z
                   the_phis = None ) :
        
        ## check arguments 
        assert isinstance ( n  , int ) and 0 <= n  < 50 , "``n''-parameter is illegal: %s"  % n
        assert isinstance ( nz , int ) and 0 <= nz < 50 , "``nz''-parameter is illegal: %s" % nz
        
        ##   inialize the base :
        PolyBase3.__init__ ( self , name , xvar , yvar , zvar ,
                             ( n + 1 ) * ( n + 2 ) *  ( nz + 1 ) / 2 - 1 , the_phis )
        #
        if self.xminmax() != self.yminmax() : 
            logger.warning( 'PolyPos3Dmix: x&y have different edges %s vs %s' % ( self.xminmax() , self.yminmax() ) )
        
        self.__n = n
        self.__nz = nz
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.Poly3DMixPositive (
            'p3Dm_%s'               % name ,
            'Poly3DMixPositive(%s)' % name ,
            self.x        ,
            self.y        ,
            self.z        ,
            self.n        ,
            self.nz       ,
            self.phi_list )
        
        ## save  the configurtaion
        self.config = {
            'name'     : self.name ,
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'zvar'     : self.zvar ,
            'n'        : self.n    ,
            'nz'       : self.nz   ,
            'the_phis' : self.phis ,
            }
       
    @property
    def n  ( self ) :
        """``n''-parameter - order/degree of 3D-polynom in x&y-directions"""
        return self.__n
    @property
    def nx ( self ) :
        """``nx''-parameter - order/degree of 3D-polynom in x-direction"""
        return self.__n
    @property
    def ny ( self ) :
        """``ny''-parameter - order/degree of 3D-polynom in z-direction"""
        return self.__n
    @property
    def nz ( self ) :
        """``nz''-parameter - order/degree of 3D-polynom in z-direction"""
        return self.__nz

models.append ( PolyPos3Dmix_pdf ) 
           
# =============================================================================
if '__main__' == __name__ : 
         
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )
    
# =============================================================================
# The END 
# =============================================================================
