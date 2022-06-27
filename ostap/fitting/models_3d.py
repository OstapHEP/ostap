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
    'Gauss3D_pdf'       , ## 3D Gaussain 
    #
    'RooKeys3D_pdf'     , ## wrapper for native RooNDKeysPdf
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
## @class PolyBase3
#  helper base class to implement various polynomial-like shapes
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
    r"""Positive (non-factorizable!) polynomial in 3D:
    
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
            self.roo_name ( 'p3_' ) ,
            'Positive 3D polynomial %s' % self.name ,
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
    r"""Positive (non-factorizable!) symmetric polynomial in 3D:
    
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
            self.roo_name ( 'p3s_' ) ,
            'Positive symmetric 3D polynomial %s' % self.name ,
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
    r"""Positive (non-factorizable!)  x<-->y symmetric polynomial in 3D:
    
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
            self.roo_name ( 'p3m_' ) ,
            'Positive 3D polynomial with mixed symmetry %s' % self.name ,
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
            self.roo_name ( 'p3m_' ) ,
            'Positive 3D polynomial with mixed symmetry %s' % self.name ,
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
            self.roo_name ( 'p3m_' ) ,
            'Positive 3D polynomial with mixed symmetry %s' % self.name ,
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



# =============================================================================        
## @class Gauss2D_pdf
#  Simple 3D gaussian function
#  @see Ostap::Models::Gauss3D 
#  @see Ostap::Math::Gauss3D 
class Gauss3D_pdf(PDF3) :
    """ Simple 3D gaussian function
    - see Ostap.Models.Gauss3D 
    - see Ostap.Math.Gauss3D 
    """
    
    ## constructor
    def __init__ ( self            ,
                   name            ,   ## the name 
                   xvar            ,   ## the variable
                   yvar            ,   ## the variable
                   zvar            ,   ## the variable
                   muX      = None ,
                   muY      = None ,
                   muZ      = None ,
                   sigmaX   = None ,
                   sigmaY   = None ,
                   sigmaZ   = None ,
                   phi      = None ,
                   theta    = None ,
                   psi      = None ) : 
    
        ## initialize the base class 
        PDF3.__init__ (  self , name , xvar , yvar , zvar )
        
        sx_lims = ()
        mx_lims = ()
        sy_lims = ()
        my_lims = ()
        sz_lims = ()
        mz_lims = ()

        xlims = self.xminmax()        
        if xlims :
            dx = xlims[1] - xlims[0] 
            sx_lims = 0 , 1.0 * dx 
            mx_lims = xlims[0] - 0.1 * dx , xlims[1] + 0.1 * dx 

        ylims = self.yminmax()
        if ylims :
            dy = ylims[1] - ylims[0] 
            sy_lims = 0 , 1.0 * dy 
            my_lims = ylims[0] - 0.1 * dy , ylims[1] + 0.1 * dy 

        zlims = self.zminmax()
        if zlims :
            dz = zlims[1] - zlims[0] 
            sz_lims = 0 , 1.0 * dy 
            mz_lims = zlims[0] - 0.1 * dz , zlims[1] + 0.1 * dz 

            
        self.__muX = self.make_var ( muX  ,
                                     'mu_x_%s'     % self.name ,
                                     '#mu_{x}(%s)' % self.name ,
                                     None , muX , *mx_lims )
        
        self.__muY = self.make_var ( muY  ,
                                     'mu_y_%s'     % self.name ,
                                     '#mu_{y}(%s)' % self.name ,
                                     None , muY , *my_lims )

        self.__muZ = self.make_var ( muZ  ,
                                     'mu_z_%s'     % self.name ,
                                     '#mu_{z}(%s)' % self.name ,
                                     None , muZ , *mz_lims )
        
        
        self.__sigmaX = self.make_var ( sigmaX  ,
                                        'sigma_x_%s'     % self.name ,
                                        '#sigma_{x}(%s)' % self.name ,
                                        None     ,  sigmaX , *sx_lims )
        
        self.__sigmaY = self.make_var ( sigmaY  ,
                                        'sigma_y_%s'     % self.name ,
                                        '#sigma_{y}(%s)' % self.name ,
                                        None     ,  sigmaY , *sy_lims )

        self.__sigmaZ = self.make_var ( sigmaZ  ,
                                        'sigma_z_%s'     % self.name ,
                                        '#sigma_{z}(%s)' % self.name ,
                                        None     ,  sigmaZ , *sy_lims )
        
        self.__phi    = self.make_var ( phi  ,
                                        'phi_%s'   % self.name ,
                                        '#phi(%s)' % self.name ,
                                        None     ,  phi , -10 , +10  )
        
        self.__theta   = self.make_var ( theta  ,
                                         'theta_%s'   % self.name ,
                                         '#theta(%s)' % self.name ,
                                         None     ,  theta , -10 , +10  )
        self.__psi    = self.make_var ( psi  ,
                                        'psi_%s'   % self.name ,
                                        '#psi(%s)' % self.name ,
                                        None     ,  psi , -10 , +10  )
        
        ## make PDF
        self.pdf = Ostap.Models.Gauss3D (
            self.roo_name ( 'gauss3d'  ) , 
            'Gauss3D %s' % self.name     ,
            self.x      ,
            self.y      ,
            self.z      ,
            self.muX    ,
            self.muY    ,
            self.muZ    ,
            self.sigmaX ,
            self.sigmaY ,
            self.sigmaZ ,
            self.phi    ,            
            self.theta  ,            
            self.psi    ,            
            )
        
        ## save configuration
        self.config = {
            'name'     : self.name   ,            
            'xvar'     : self.xvar   ,
            'yvar'     : self.yvar   ,            
            'zvar'     : self.zvar   ,            
            'muX'      : self.muX    ,
            'muY'      : self.muY    ,
            'muZ'      : self.muZ    ,
            'sigmaX'   : self.sigmaX ,
            'sigmaY'   : self.sigmaY ,
            'sigmaZ'   : self.sigmaZ ,
            'phi'      : self.phi    ,
            'theta'    : self.theta  ,
            'psi'      : self.psi  
            }
        
    @property
    def muX ( self ) :
        """``x-location for 3D gaussian"""
        return self.__muX
    @muX.setter
    def muX ( self , value ) :
        self.set_value ( self.__muX , value )

    @property
    def muY ( self ) :
        """``y-location for 3D gaussian"""
        return self.__muY
    @muY.setter
    def muY ( self , value ) :
        self.set_value ( self.__muY , value )

    @property
    def muZ ( self ) :
        """``z-location for 3D gaussian"""
        return self.__muZ
    @muZ.setter
    def muZ ( self , value ) :
        self.set_value ( self.__muZ , value )

    @property
    def sigmaX ( self ) :
        """``sigma-X for 3D gaussian"""
        return self.__sigmaX
    @sigmaX.setter
    def sigmaX ( self , value ) :
        self.set_value ( self.__sigmaX , value )

    @property
    def sigmaY ( self ) :
        """``sigma-Y for 3D gaussian"""
        return self.__sigmaY
    @sigmaY.setter
    def sigmaY ( self , value ) :
        self.set_value ( self.__sigmaY , value )

    @property
    def sigmaZ ( self ) :
        """``sigma-Z for 3D gaussian"""
        return self.__sigmaZ
    @sigmaZ.setter
    def sigmaZ ( self , value ) :
        self.set_value ( self.__sigmaZ , value )
        
    @property
    def phi  ( self ) :
        """``phi'' : Euler angle ``phi'' for 3D gaussian"""
        return self.__phi
    @phi.setter
    def phi  ( self , value ) :
        self.set_value ( self.__phi , value )

    @property
    def theta ( self ) :
        """``theta'' : Euler anlgle ``theta'' for 3D gaussian"""
        return self.__theta
    @theta.setter
    def theta ( self , value ) :
        self.set_value ( self.__theta , value )

    @property
    def psi ( self ) :
        """``psi'' : Euler angle ``psi'' for 3D gaussian"""
        return self.__psi
    @psi.setter
    def psi ( self , value ) :
        self.set_value ( self.__psi , value )


# ==============================================================================
## @class RooKeys3D_pdf
#  Trivial Ostap wrapper for the native <code>RooNDKeysPdf</code> from RooFit
#  @code
#  xvar = ...
#  yvar = ...
#  zvar = ...
#  pdf  = RooKeys3D_pdf ( 'P4' , xvar , yvar ,  zvar , data  ) ;
#  @endcode
#  @see RooNDKeysPdf
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2019-04-27
class RooKeys3D_pdf(PDF3) :
    """Trivial Ostap wrapper for the native RooNDKeysPdf from RooFit
    - see ROOT.RooKeysPdf
    >>> xvar = ...
    >>> yvar = ...
    >>> zvar = ...
    >>> pdf  = RooKeys3D_pdf ( 'Keys' , xvar , yvar ,  zvar , data  )
    """
    ## constructor
    def __init__ ( self           ,
                   name           ,   ## the name 
                   xvar           ,   ## the variable
                   yvar           ,   ## the variable
                   zvar           ,   ## the variable
                   data           ,   ## data set 
                   options = 'ma' ,   ## options 
                   rho     = 1    ,   ## global scale 
                   nsigma  = 3    ,   ##
                   rotate  = True ,   
                   sort    = True ) : 
        
        ## initialize the base class 
        PDF3.__init__ (  self , name , xvar , yvar , zvar )

        self.__data    = data
        self.__options = mirror
        self.__rho     = rho 
        self.__options = options 
        self.__nsigma  = nsigma 
        self.__rotate  = True if rotate else False 
        self.__sort    = True if sort   else False 

        self.__keys_vlst = ROOT.RooArgList()
        self.__keys_vlst.Add ( self.xvar )
        self.__keys_vlst.Add ( self.yvar )
        self.__keys_vlst.Add ( self.zvar )
        
        ## create PDF
        self.pdf = ROOT.RooNDKeysPdf (
            self.roo_name ( 'keys3_' ) ,
            'Keys 3D %s' % self.name ,
            self.__keys_vlst  ,
            self.data    ,
            self.options , 
            self.rho     , 
            self.nsigma  , 
            self.rotate  , 
            self.sort    ) 

        ## save configuration
        self.config = {
            'name'    : self.name    ,
            'xvar'    : self.xvar    ,
            'yvar'    : self.yvar    ,
            'zvar'    : self.yvar    ,
            'data'    : self.data    ,            
            'options' : self.options ,            
            'rho'     : self.rho     ,            
            'nsigma'  : self.nsigma  ,            
            'rotate'  : self.rotate  ,            
            'sort'    : self.sort    ,            
            }

    @property
    def data   ( self ) :
        """``data'' : the actual data set for RooNDKeysPdf"""
        return self.__data
    @property
    def options ( self ) :
        """``options'' : ``ootions'' string for RooNDKeysPdf"""
        return self.__mirror        
    @property
    def rho    ( self )  :
        """``rho'' : ``rho'' parameter for RooNDKeysPdf"""
        return self.__rho 
    @property
    def rotate    ( self )  :
        """``rotate'' : ``rotate'' flag for RooNDKeysPdf"""
        return self.__rotate
    @property
    def sort      ( self )  :
        """``sort'' : ``sort'' flag for RooNDKeysPdf"""
        return self.__sort 
    @property
    def nsigma    ( self )  :
        """``nsigma'' : ``nsigma'' parameter for RooNDKeysPdf"""
        return self.__nsigma 


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
##                                                                      The END 
# =============================================================================
