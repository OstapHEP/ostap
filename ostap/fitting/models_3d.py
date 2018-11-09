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
    'PolyPos3D_pdf'     , ## A positive polynomial in 3D  
    'PolyPos3Dsym_pdf'  , ## A positive symmetric polynomial in 3D
    'PolyPos3DmixXY_pdf', ## A positive partly symmetric (x<-->y) polynomial in 3D 
    'PolyPos3DmixYZ_pdf', ## A positive partly symmetric (y<-->z) polynomial in 3D 
    'PolyPos3DmixXZ_pdf', ## A positive partly symmetric (x<-->z) polynomial in 3D 
    #
    'make_B3D'          , ## create                         3D "background" function 
    'make_B3Dsym'       , ## create symmetric               3D "background" function 
    'make_B3DmixXY'     , ## create mixed symmetry (x<-->y) 3D "background" function 
    'make_B3DmixYZ'     , ## create mixed symmetry (y<-->z) 3D "background" function 
    'make_B3DmixXZ'     , ## create mixed symmetry (x<-->z) 3D "background" function 
    )
# =============================================================================
import ROOT, math
from   ostap.core.core     import cpp, Ostap
from   ostap.math.base     import iszero
from   ostap.fitting.fit3d import PDF3, Flat3D 
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
## @class PolyPos3DmixXY_pdf
#  The 3D-polynomial of order N*N*Nz, that is constrained 
#  to be non-negative and x<-->y symmetric over the  defined range      
#   \f[  P(x,y,z) = \sum_{i,j,k} a_{ijk}B^{n}_i(x) B^{n}_j(y) B^{n}_k(z)\f] 
#  where all coefficients \f$a_{ijk}\f$ are:
# - non-negative: \f$ a_{ijk}\ge0 \f$
# - patly symmetric: \f$ a_{ijk}=a_{jik}\f$
# - constrainted: \f$ \sum_{i,j,k} a_{ijk}=1 \f$ 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2017-11-14
#  @see Ostap::Models::Poly3DMixPositive
#  @see Ostap::Math::Positive3DMix
class PolyPos3DmixXY_pdf(PolyBase3) :
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
                   nz = 1          ,   ## polynomial degree in Z
                   the_phis = None ) :
        
        ## check arguments 
        assert isinstance ( n  , int ) and 0 <= n  < 50 , "``n''-parameter is illegal: %s"  % n
        assert isinstance ( nz , int ) and 0 <= nz < 50 , "``nz''-parameter is illegal: %s" % nz
        
        ##   inialize the base :
        PolyBase3.__init__ ( self , name , xvar , yvar , zvar ,
                             ( n + 1 ) * ( n + 2 ) *  ( nz + 1 ) / 2 - 1 , the_phis )
        #
        if self.xminmax() != self.yminmax() : 
            logger.warning( 'PolyPos3DmixXY: x&y have different edges %s vs %s' % ( self.xminmax() , self.yminmax() ) )
        
        self.__n  = n
        self.__nz = nz
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.Poly3DMixPositive (
            'p3DmXY_%s'               % name ,
            'Poly3DMixPositiveXY(%s)' % name ,
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

models.append ( PolyPos3DmixXY_pdf ) 


# =============================================================================
## @class PolyPos3DmixYZ_pdf
#  Positive (non-factorizable!)  y<-->z symmetric polynomial in 3D.
#
#  The 3D-polynomial of order Nx*N*N, that is constrained 
#  to be non-negative and y<-->z symmetric over the  defined range      
#   \f[  P(x,y,z) = \sum_{i,j,k} \alpha_{ijk}B^{n_x}_i(x) B^{n}_j(y) B^{n}_k(z)\f] 
#  where all coefficients \f$ \alpha_{ijk}\f$ are:
# - non-negative:     \f$ \alpha_{ijk} \ge 0 \f$
# - partly symmetric: \f$ \alpha_{ijk}= \alpha_{ikj}\f$
# - constrainted:     \f$ \sum_{i,j,k} \alpha_{ijk}=1 \f$ 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2017-11-14
#  @see Ostap::Models::Poly3DMixPositive
#  @see Ostap::Math::Positive3DMix
class PolyPos3DmixYZ_pdf(PolyBase3) :
    r"""Positive (non-factorizable!)  y<-->z symmetric polynomial in 3D:
    
    The 3D-polynomial of order Nx*N*N, that is constrained 
    to be non-negative ans symmetric over the  defined range
    
    - P(x,y,z) = \sum_{i,j,k} a_{ijk}B^{n_x}_i(x) B^{n}_j(y) B^{n_z}_k(z)
    
    where  B^n_i - are Bernstein polynomials and all coefficients \a_{ijk} are:
    - non-negative:       0<= a_{ijk}
    - partly symmetric:   a_{ijk}=a_{ikj}
    - constrained:        \sum a_{ijk} = 1 
        
    """
    def __init__ ( self            ,
                   name            ,
                   xvar            ,   ## the first  dimension  
                   yvar            ,   ## the second dimension
                   zvar            ,   ## the third  dimension
                   nx = 1          ,   ## polynomial degree in X
                   n  = 1          ,   ## polynomial degree in Y,Z
                   the_phis = None ) :
        
        ## check arguments 
        assert isinstance ( nx , int ) and 0 <= nx < 50 , "``nx''-parameter is illegal: %s" % nx
        assert isinstance ( n  , int ) and 0 <= n  < 50 , "``n''-parameter  is illegal: %s" % n
        
        ##   inialize the base :
        PolyBase3.__init__ ( self , name , xvar , yvar , zvar ,
                             ( n + 1 ) * ( n + 2 ) *  ( nx + 1 ) / 2 - 1 , the_phis )
        #
        if self.yminmax() != self.zminmax() : 
            logger.warning( 'PolyPos3DmixYZ: y&z have different edges %s vs %s' % ( self.yminmax() , self.zminmax() ) )
        
        self.__nx = nx
        self.__n  = n
        
        # =====================================================================
        ## finally build PDF 
        # =====================================================================
        self.pdf = Ostap.Models.Poly3DMixPositive (
            'p3DmYZ_%s'               % name ,
            'Poly3DMixPositiveYZ(%s)' % name ,
            self.y        , ## note "strange" order!
            self.z        , ## note "strange" order! 
            self.x        , ## note "strange" order! 
            self.n        , ## note "strange" order! 
            self.nx       , ## note "strange" order! 
            self.phi_list )
        
        # =====================================================================
        # ATTENTION! 
        # due to the swap of variables for the underlying function
        # some tricks&shortcuts must be disabled! 
        self.tricks = False                                       ## ATTENTION! 
        # =====================================================================
        
        ## save  the configurtaion
        self.config = {
            'name'     : self.name ,
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'zvar'     : self.zvar ,
            'nx'       : self.nx   ,
            'n'        : self.n    ,
            'the_phis' : self.phis ,
            }
       
    @property
    def n  ( self ) :
        """``n''-parameter - order/degree of 3D-polynom in y&z-directions"""
        return self.__n
    @property
    def nx ( self ) :
        """``nx''-parameter - order/degree of 3D-polynom in x-direction"""
        return self.__nx
    @property
    def ny ( self ) :
        """``ny''-parameter - order/degree of 3D-polynom in y-direction"""
        return self.__n
    @property
    def nz ( self ) :
        """``nz''-parameter - order/degree of 3D-polynom in z-direction"""
        return self.__n

models.append ( PolyPos3DmixYZ_pdf ) 



# =============================================================================
## @class PolyPos3DmixXZ_pdf
#  Positive (non-factorizable!)  x<-->z symmetric polynomial in 3D.
#
#  The 3D-polynomial of order N*Ny*N, that is constrained 
#  to be non-negative and x<-->z symmetric over the  defined range      
#   \f[  P(x,y,z) = \sum_{i,j,k} \alpha_{ijk}B^{n}_i(x) B^{n_y}_j(y) B^{n}_k(z)\f] 
#  where all coefficients \f$ \alpha_{ijk}\f$ are:
# - non-negative:     \f$ \alpha_{ijk} \ge 0 \f$
# - partly symmetric: \f$ \alpha_{ijk}= \alpha_{kji}\f$
# - constrainted:     \f$ \sum_{i,j,k} \alpha_{ijk}=1 \f$ 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2017-11-14
#  @see Ostap::Models::Poly3DMixPositive
#  @see Ostap::Math::Positive3DMix
class PolyPos3DmixXZ_pdf(PolyBase3) :
    r"""Positive (non-factorizable!)  x<-->z symmetric polynomial in 3D:
    
    The 3D-polynomial of order N*Ny*N, that is constrained 
    to be non-negative ans symmetric over the  defined range
    
    - P(x,y,z) = \sum_{i,j,k} a_{ijk}B^{n}_i(x) B^{n_y}_j(y) B^{n_z}_k(z)
    
    where  B^n_i - are Bernstein polynomials and all coefficients \a_{ijk} are:
    - non-negative:       0<= a_{ijk}
    - partly symmetric:   a_{ijk}=a_{kji}
    - constrained:       \sum a_{ijk} = 1 
        
    """
    def __init__ ( self            ,
                   name            ,
                   xvar            ,   ## the first  dimension  
                   yvar            ,   ## the second dimension
                   zvar            ,   ## the third  dimension
                   n  = 1          ,   ## polynomial degree in X,Z
                   ny = 1          ,   ## polynomial degree in Y
                   the_phis = None ) :
        
        ## check arguments 
        assert isinstance ( ny , int ) and 0 <= ny < 50 , "``ny''-parameter is illegal: %s" % ny
        assert isinstance ( n  , int ) and 0 <= n  < 50 , "``n''-parameter  is illegal: %s" % n
        
        ##   inialize the base :
        PolyBase3.__init__ ( self , name , xvar , yvar , zvar ,
                             ( n + 1 ) * ( n + 2 ) *  ( ny + 1 ) / 2 - 1 , the_phis )
        #
        if self.xminmax() != self.zminmax() : 
            logger.warning( 'PolyPos3DmixXZ: x&z have different edges %s vs %s' % ( self.xminmax() , self.zminmax() ) )
        
        self.__n  = n
        self.__ny = ny
        
        # =====================================================================
        ## finally build PDF 
        # =====================================================================
        self.pdf = Ostap.Models.Poly3DMixPositive (
            'p3DmXZ_%s'               % name ,
            'Poly3DMixPositiveXZ(%s)' % name ,
            self.x        , ## note "strange" order!
            self.z        , ## note "strange" order! 
            self.y        , ## note "strange" order! 
            self.n        , ## note "strange" order! 
            self.ny       , ## note "strange" order! 
            self.phi_list )
        
        # =====================================================================
        # ATTENTION! 
        # due to the swap of variables for the underlying function
        # some tricks&shortcuts must be disabled! 
        self.tricks = False                                       ## ATTENTION! 
        # =====================================================================
        
        ## save  the configurtaion
        self.config = {
            'name'     : self.name ,
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'zvar'     : self.zvar ,
            'n'        : self.n    ,
            'ny'       : self.ny   ,
            'the_phis' : self.phis ,
            }
       
    @property
    def n  ( self ) :
        """``n''-parameter - order/degree of 3D-polynom in x&z-directions"""
        return self.__n
    @property
    def nx ( self ) :
        """``nx''-parameter - order/degree of 3D-polynom in x-direction"""
        return self.__n
    @property
    def ny ( self ) :
        """``ny''-parameter - order/degree of 3D-polynom in y-direction"""
        return self.__ny
    @property
    def nz ( self ) :
        """``nz''-parameter - order/degree of 3D-polynom in z-direction"""
        return self.__n

models.append ( PolyPos3DmixXZ_pdf ) 



# ==============================================================================
## Easy creation of  3D function for background
#  @code
#  xvar = ...
#  yvar = ...
#  zvar = ...
#  bkg  = make_B3D ( 'BBB' , xvar , yvar , zvar , 0 , 3 , 2 ) ## create PolyPos3D_pdf  
#  bkg  = make_B3D ( 'BBB' , xvar , yvar , zvar , 0 , 0 , 0 ) ## create Flat3D 
#  endcode
def make_B3D ( name , xvar , yvar , zvar , nx , ny , nz ) :
    """Easy creation of  3D function for background
    >>> xvar = ...
    >>> yvar = ...
    >>> zvar = ...
    >>> bkg  = make_B3D ( 'BBB' , xvar , yvar , zvar , 1 , 1 , 2 ) ## create PolyPos3D_pdf
    >>> bkg  = make_B3D ( 'BBB' , xvar , yvar , zvar , 0 , 0 , 0 ) ## create Flat3D 
    """
    
    if   0 == nx and 0 == ny and 0 == nz :
        return Flat3D    ( name = name , xvar = xvar , yvar = yvar , zvar = zvar )
    
    return PolyPos3D_pdf ( name = name , xvar = xvar , yvar = yvar , zvar = zvar ,
                           nx   = abs ( nx )  ,
                           ny   = abs ( ny )  ,
                           nz   = abs ( nz )  )     

# ==============================================================================
## Easy creation of  symmetric 3D function for background
#  @code
#  xvar = ...
#  yvar = ...
#  zvar = ...
#  bkg  = make_B3Dsym ( 'BBB' , xvar , yvar , zvar , 1 ) ## create PolyPos3Dsym_pdf 
#  bkg  = make_B3Dsym ( 'BBB' , xvar , yvar , zvar , 0 ) ## create Flat3D 
#  endcode
def make_B3Dsym ( name , xvar , yvar , zvar , n ) :
    """Easy creation of  symmetric 3D function for background
    >>> xvar = ...
    >>> yvar = ...
    >>> zvar = ...
    >>> bkg  = make_B3Dsym ( 'BBB' , xvar , yvar , zvar , 1 ) ## create PolyPos3Dsym_pdf
    >>> bkg  = make_B3Dsym ( 'BBB' , xvar , yvar , zvar , 0 ) ## create Flat3D 
    """
    
    if   0 == n :
        return Flat3D       ( name = name , xvar = xvar , yvar = yvar , zvar = zvar )
    
    return PolyPos3Dsym_pdf ( name = name , xvar = xvar , yvar = yvar , zvar = zvar ,
                              n    = abs ( n ) )


# ==============================================================================
## Easy creation of  3D function for background with mixed  x<-->y symmetry 
#  @code
#  xvar = ...
#  yvar = ...
#  zvar = ...
#  bkg  = make_B3DmixXY ( 'BBB' , xvar , yvar , zvar , 0 , 3 ) ## create PolyPos3DmixXY_pdf 
#  bkg  = make_B3DmixXY ( 'BBB' , xvar , yvar , zvar , 0 , 0 ) ## create Flat3D 
#  endcode
def make_B3DmixXY ( name , xvar , yvar , zvar , n , nz ) :
    """Easy creation of mixed symmetry (x<-->y) 3D function for background
    >>> xvar = ...
    >>> yvar = ...
    >>> zvar = ...
    >>> bkg  = make_B3DmixXY ( 'BBB' , xvar , yvar , zvar , 1 , 2 ) ## create PolyPos3DmixXY_pdf  
    >>> bkg  = make_B3DmixXY ( 'BBB' , xvar , yvar , zvar , 0 , 0 ) ## create Flat3D 
    """
    if   0 == n and 0 == nz :
        return Flat3D         ( name = name , xvar = xvar , yvar = yvar , zvar = zvar )
    
    return PolyPos3DmixXY_pdf ( name = name , xvar = xvar , yvar = yvar , zvar = zvar ,
                                n    = abs ( n  )  ,
                                nz   = abs ( nz )  )     


# ==============================================================================
## Easy creation of  3D function for background with mixed  y<-->z symmetry 
#  @code
#  xvar = ...
#  yvar = ...
#  zvar = ...
#  bkg  = make_B3DmixYZ ( 'BBB' , xvar , yvar , zvar , 0 , 3 ) ## create PolyPos3DmixYZ_pdf 
#  bkg  = make_B3DmixYZ ( 'BBB' , xvar , yvar , zvar , 0 , 0 ) ## create Flat3D 
#  endcode
def make_B3DmixYZ ( name , xvar , yvar , zvar , nx , n ) :
    """Easy creation of mixed symmetry (x<-->y) 3D function for background
    >>> xvar = ...
    >>> yvar = ...
    >>> zvar = ...
    >>> bkg  = make_B3DmixYZ ( 'BBB' , xvar , yvar , zvar , 1 , 2 ) ## create PolyPos3DmixYZ_pdf  
    >>> bkg  = make_B3DmixYZ ( 'BBB' , xvar , yvar , zvar , 0 , 0 ) ## create Flat3D 
    """
    if   0 == nx and 0 == n :
        return Flat3D         ( name = name , xvar = xvar , yvar = yvar , zvar = zvar )
    
    return PolyPos3DmixYZ_pdf ( name = name , xvar = xvar , yvar = yvar , zvar = zvar ,
                                nx   = abs ( nx ) ,
                                n    = abs ( n  ) )     

# ==============================================================================
## Easy creation of  3D function for background with mixed  x<-->z symmetry 
#  @code
#  xvar = ...
#  yvar = ...
#  zvar = ...
#  bkg  = make_B3DmixXZ ( 'BBB' , xvar , yvar , zvar , 0 , 3 ) ## create PolyPos3DmixXZ_pdf 
#  bkg  = make_B3DmixZZ ( 'BBB' , xvar , yvar , zvar , 0 , 0 ) ## create Flat3D 
#  endcode
def make_B3DmixXZ ( name , xvar , yvar , zvar , n , ny ) :
    """Easy creation of mixed symmetry (x<-->z) 3D function for background
    >>> xvar = ...
    >>> yvar = ...
    >>> zvar = ...
    >>> bkg  = make_B3DmixXZ ( 'BBB' , xvar , yvar , zvar , 1 , 2 ) ## create PolyPos3DmixXZ_pdf  
    >>> bkg  = make_B3DmixXZ ( 'BBB' , xvar , yvar , zvar , 0 , 0 ) ## create Flat3D 
    """
    if   0 == nx and 0 == n :
        return Flat3D         ( name = name , xvar = xvar , yvar = yvar , zvar = zvar )
    
    return PolyPos3DmixXZ_pdf ( name = name , xvar = xvar , yvar = yvar , zvar = zvar ,
                                n    = abs ( n  ) ,
                                ny   = abs ( ny ) )     


# =============================================================================
if '__main__' == __name__ : 
         
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )
    
# =============================================================================
# The END 
# =============================================================================
