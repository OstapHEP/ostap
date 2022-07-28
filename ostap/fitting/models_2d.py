#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/models_2d.py
#  Smooth non-factorizable 2D-models to describe background distribtions
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
# =============================================================================
""" Set of useful non-factorizable 2D-models to describe background distribtions
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    'PolyPos2D_pdf'   , ## A positive polynomial in 2D  
    'PolyPos2Dsym_pdf', ## A positive symmetric polynomial in 2D
    'PSPol2D_pdf'     , ## Product of phase spaces, modulated with 2D polynomial
    'PSPol2D2_pdf'    , ## Product of phase spaces, modulated with 2D polynomial
    'PSPol2D3_pdf'    , ## Product of phase spaces, modulated with 2D polynomial
    'PSPol2Dsym_pdf'  , ## Symmetric product of phase spaces, modulated with 2D polynomial
    'PSPol2D2sym_pdf' , ## Symmetric product of phase spaces, modulated with 2D polynomial
    'PSPol2D3sym_pdf' , ## Symmetric product of phase spaces, modulated with 2D polynomial
    'ExpoPSPol2D_pdf' , ## Exponential times  phase space times positive 2D-polynomial
    'ExpoPol2D_pdf'   , ## Product of exponents times positive 2D-polynomial
    'ExpoPol2Dsym_pdf', ## Symmetric version of above
    ##
    'Spline2D_pdf'    , ## 2D generic   positive spline 
    'Spline2Dsym_pdf' , ## 2D symmetric positive spline
    #
    'Gauss2D_pdf'     , ## 2D gaussian 
    # 
    'RooKeys2D_pdf'   , ## wrapper for native RooNDKeysPdf
    #
    'make_B2D'        , ## create           2D "background" function 
    'make_B2Dsym'     , ## create symmetric 2D "background" function 
    )
# =============================================================================
from   ostap.core.core          import cpp, Ostap
from   ostap.math.base          import iszero
from   ostap.fitting.fithelpers import Phases
from   ostap.fitting.fit2d      import PDF2, Flat2D
from   ostap.fitting.signals    import Gauss_pdf, CB2_pdf
from   ostap.core.meta_info     import root_info
import ROOT 
# =============================================================================
from   ostap.logger.logger     import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.models_2d' )
else                       : logger = getLogger ( __name__                  )
# =============================================================================
models = []
# =============================================================================
##  @class PolyBase2
#   helper base class to implement various polynomial-like shapes
class PolyBase2(PDF2,Phases) :
    """Helper base class to implement various polynomial-like shapes
    """
    def __init__ ( self , name , xvar  , yvar , power , the_phis = None ) :
        PDF2  .__init__ ( self , name  = name  , xvar = xvar , yvar = yvar  )
        Phases.__init__ ( self , power = power , the_phis = the_phis  )
# =============================================================================
## @class PolyPos2D_pdf
#  positive polynomial in 2D:
#  \f[ P(x,y) =
#   \sum^{n}_{i=0} \sum{k}_{j=0} a_{ij} B^n_i(x) B^k_j(y)
#  \f] 
#  where
#  - \f$ B^n_i(x)\f$ denotes the basic Bersntein polynomial
#  - \f$ a_{ij} \ge 0 \f$
#  - \f$ \sum_{i,j} a_{ij} = 1 \f$
#  @see Ostap::Models::Poly2DPositive
#  @see Ostap::Math::Poly2DPositive
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class PolyPos2D_pdf(PolyBase2) :
    r"""Positive (non-factorizable!) polynomial in 2D:
    
    P_{n,k}(x,y) = \sum^{n}_{i=0}\sum{k}_{j=0} a_{ij} B^n_i(x) B^k_j(y)
    
    where:
    - B^n_i - are Bernstein polynomials
    - a_{ij} >= 0
    - \sum_{ij} a_{ij} = 1
    
    Note:
    
    - f(x,y)>=0 for whole 2D-range
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,   ##  the first  dimension  
                   yvar             ,   ##  the second dimension
                   nx       = 2     ,   ##  polynomial degree in X 
                   ny       = 2     ,   ##  polynomial degree in Y 
                   the_phis = None  ) : 

        ## check arguments 
        assert isinstance ( nx , int ) and 0 <= nx < 100 , "'nx'-parameter is illegal: %s" % nx 
        assert isinstance ( ny , int ) and 0 <= ny < 100 , "'ny'-parameter is illegal: %s" % ny
        ## 
        PolyBase2.__init__ ( self ,
                             name = name ,
                             xvar = xvar ,
                             yvar = yvar ,
                             power = ( nx + 1 ) * ( ny + 1 ) - 1 , the_phis = the_phis )

        self.__nx = nx 
        self.__ny = ny
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.Poly2DPositive (
            self.roo_name ( 'pol2_'  ) , 
            'Positive 2D polynomial %s' % self.name ,
            self.xvar     ,
            self.yvar     ,
            self.nx       ,
            self.ny       , 
            self.phi_list )
        
        ## save configuration
        self.config = {
            'name'     : self.name ,            
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'nx'       : self.nx   ,            
            'ny'       : self.ny   ,
            'the_phis' : self.phis ,
            }

    @property
    def nx ( self ) :
        """'nx' : order/degree of 2D-polynom in x-direction"""
        return self.__nx
    @property
    def ny ( self ) :
        """'ny' : order/degree of 2D-polynom in y-direction"""
        return self.__ny
    
        
models.append ( PolyPos2D_pdf ) 
# =============================================================================
## @class PolyPos2Dsym_pdf
#  Positive symetric polynomial in 2D. Symmetrized version of PolyPos2D_pdf
#  \f[ P_{n}(x,y) =
#   \sum^{n}_{i=0} \sum{k}_{j=0} a_{ij} B^n_i(x) B^k_j(y)
#  \f] 
#  where
#  - \f$ B^n_i(x)\f$ denotes the basic bersntein polynomial
#  - \f$ a_{ij} \ge 0 \f$
#  - \f$ a_{ji}=a_{i,j} \f$ 
#  - \f$ \sum_{i,j} a_{ij} = 1 \f$
#  @see Ostap::Models::Poly2DSymPositive
#  @see Ostap::Math::Poly2DSymPositive
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class PolyPos2Dsym_pdf(PolyBase2) :
    r"""Positive (non-factorizable!) SYMMETRIC polynomial in 2D:
    
    P_{n}(x,y) = \sum^{n}_{i=0} \sum{k}_{j=0} a_{ij} B^n_i(x) B^k_j(y)
    
    where:
    - B^n_i - are Bernstein polynomials
    - a_{ij} = a_{ji}
    - a_{ij} \ge 0 
    - a_{ji}=a_{i,j} 
    - \sum_{i,j} a_{ij} = 1
    
    Note:
    - f(x,y)>=0 for whole 2D-range
    - f(x,y) = f(y,x)
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,   ##  the first  dimension  
                   yvar             ,   ##  the second dimension
                   n        = 2     ,   ##  polynomial degree
                   the_phis = None  ) : 
        
        ## check arguments 
        assert isinstance ( n , int ) and 0 <= n < 100 , "'n'-parameter is illegal: %s" % n
        ## 
        self.__n = n 
        PolyBase2.__init__ ( self , name , xvar , yvar , ( n + 1 ) * ( n + 2 ) / 2 - 1 , the_phis )

        if self.xminmax() != self.yminmax() :
            logger.warning( 'PolyPos2Dsym: x&y have different edges %s vs %s' % ( self.xminmax() , self.yminmax() ) )

        
        ## finally build PDF 
        self.pdf = Ostap.Models.Poly2DSymPositive (
            self.roo_name ( 'pol2sym_'  ) , 
            'Positive symmetric 2D polynomial %s' % self.name ,
            self.x        ,
            self.y        ,
            self.n        ,
            self.phi_list )
        
        ## save configuration
        self.config = {
            'name'     : self.name ,            
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'n'        : self.n    ,            
            'the_phis' : self.phis ,
            }
    
    @property
    def n  ( self ) :
        """'n'  :  order/degree of 2D-polynom in x&y-directions"""
        return self.__n
    @property
    def nx ( self ) :
        """'nx' : order/degree of 2D-polynom in x-direction"""
        return self.__n
    @property
    def ny ( self ) :
        """'ny' : order/degree of 2D-polynom in y-direction"""
        return self.__n
       
models.append ( PolyPos2Dsym_pdf )

# =============================================================================
## @class PSPol2D_pdf
#  The 2D-function, that represent a cross-product of two phase-space factors,
#
#  The function is:
#    \f[ f(x,y) = 
#        \Phi_{k,l}(x;x_{low}, x_{high})
#        \Phi_{m,n}(y;y_{low}, y_{high})
#        P_{N,M}(x,y) \f]
#   where  
#  - \f$ \Phi_{k,l}(x;x_{low},x_{high}) \f$ is a phase-space function for x-axis
#  - \f$ \Phi_{m,n}(y;y_{low},y_{high}) \f$ is a phase-space function for y-axis
#  - \f$ P_{N,M}(x,y) \f$ is 2D positive Bernstein polynomial
#
#  @see Ostap::Models::PS2DPol
#  @see Ostap::Math::PS2DPol
#  @see Ostap::Math::PhaseSpaceNL 
#  @see Ostap::Math::Positive2D
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class PSPol2D_pdf(PolyBase2) :
    r"""Product of phase space factors, modulated by the positive polynom in 2D
    
    f (x,y) =
    = \Phi_{k,l}(x;x_{low}, x_{high}) 
    * \Phi_{m,n}(y;y_{low}, y_{high})
    *  P_{N,M}(x,y) 
    
    where:
    
    - \Phi_{k,l}(x;x_{low},x_{high}) is a phase-space function for x-axis
    - \Phi_{m,n}(y;y_{low},y_{high}) is a phase-space function for y-axis
    - \P_{N,M}(x,y) is 2D positive Bernstein polynomial
        
    Note:
    - f(x,y)>=0 for whole 2D-range
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,   ##  the first  dimension  
                   yvar             ,   ##  the second dimension
                   psx              ,   ##  phase space in X, Ostap::Math::PhaseSpaceNL 
                   psy              ,   ##  phase space in Y, Ostap::Math::PhaseSpaceNL 
                   nx       = 2     ,   ##  polynomial degree in X 
                   ny       = 2     ,   ##  polynomial degree in Y 
                   the_phis = None  ) :
        
        ## check arguments 
        assert isinstance ( nx , int ) and 0 <= nx < 100 , "'nx'-parameter is illegal: %s" % nx
        assert isinstance ( ny , int ) and 0 <= ny < 100 , "'ny'-parameter is illegal: %s" % ny

        ## the base 
        PolyBase2.__init__ ( self , name , xvar , yvar ,
                             ( nx + 1 ) * ( ny + 1 ) - 1 , the_phis )
        
        self.__phasespacex = psx  
        self.__phasespacey = psy
        
        self.__nx          =  nx
        self.__ny          =  ny

        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.PS2DPol (            
            self.roo_name ( 'ps2_'  ) , 
            'Product of phase space factors (with polynomials) %s' % self.name ,
            self.x           ,
            self.y           ,
            self.phasespacex ,
            self.phasespacey ,
            self.nx          ,
            self.ny          , 
            self.phi_list    )
        
        ## save configuration
        self.config = {
            'name'     : self.name ,            
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'psx'      : self.phasespacex , 
            'psy'      : self.phasespacey ,
            'nx'       : self.nx   , 
            'ny'       : self.ny   , 
            'the_phis' : self.phis ,
            }
    
    @property 
    def mass1 ( self ) :
        """'mass1'-variable for the fit (alias for 'x' and 'xvar')"""
        return self.xvar
    
    @property 
    def mass2 ( self ) :
        """'mass2'-variable for the fit (alias for 'y' or 'yvar')"""
        return self.yvar

    @property
    def phasespacex( self ) :
        """'x-phasespace'  : function for  PSPol2D-function"""
        return self.__phasespacex
    
    @property
    def phasespacey( self ) :
        """'y-phasespace' : function for  PSPol2D-function"""
        return self.__phasespacey 
    
    @property
    def nx ( self ) :
        """'nx' : order/degree of 2D-polynom in x-direction"""
        return self.__nx
    @property
    def ny ( self ) :
        """'ny' : order/degree of 2D-polynom in y-direction"""
        return self.__ny
    
        
models.append ( PSPol2D_pdf ) 


# =============================================================================
## @class PSPol2D2_pdf
#  The 2D-function, that represent non-factorizeable "product" of  
#  phase-space functions modulated by the 2D-positive polynomial.
#  The  function is useful to describe e.g. 2D-distributions of 
#  \f$ m_{23}\f$ vs \f$m_{45}\f$ from 5-body decays. 
#
#  The function is:
#  \f[ f(x,y) = \frac{1}{2}
#   \left( \Phi_{k,n}(x;x_{low},x_{high}) \Phi_{l,m-1}(y,y_{low},m_{max}-x) 
#        + \Phi_{l,m}(y;y_{low},y_{high}) \Phi_{k,n-1}(x,x_{low},m_{max}-y) 
#   \right) P_{N^{x},N^{y}}(x,y) \f]
#  where 
#  - \f$ \Phi_{i,j}(z;z_{low},z_{high}\f$ are normalized phase space functions
#    for mass of \f$i\f$-particles from \f$j\f$-body decays
#  - \f$ P_{N^{x},N^{y}}(x,y) \f$ is 2D positive Bernstein polynomial
#  - \f$m_{max}\f$ is a maximal allowed mass for \f$x+y\f$
#  
#  @see Ostap::Models::PS2DPol2
#  @see Ostap::Math::PS2DPol2
#  @see Ostap::Math::PhaseSpaceNL 
#  @see Ostap::Math::Positive2D
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class PSPol2D2_pdf(PolyBase2) :
    r"""Product of phase space factors, modulated by the positive polynom in 2D
    
    f(x,y) =
    \frac{1}{2} \left(
    \Phi_{k,n}(x;x_{low},x_{high}) \Phi_{l,m-1}(y,y_{low},m_{max}-x) 
    + \Phi_{l,m}(y;y_{low},y_{high}) \Phi_{k,n-1}(x,x_{low},m_{max}-y) \right)
    P_{N^{x},N^{y}}(x,y)
    
    where 
    - \Phi_{i,j}(z;z_{low},z_{high} are normalized phase space functions
    for mass of i-particles from j-body decays
    -  P_{N^{x},N^{y}}(x,y)  is 2D positive Bernstein polynomial
    - m_{max} is a maximal allowed mass for \f$x+y\f$
    
    Note:
    - f(x,y)>=0 for whole 2D-range
    
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,   ##  the first  dimension  
                   yvar             ,   ##  the second dimension
                   psx              ,   ##  phase space in X, Ostap::Math::PhaseSpaceNL 
                   psy              ,   ##  phase space in Y, Ostap::Math::PhaseSpaceNL
                   mmax             ,   ##  max-mass 
                   nx       = 2     ,   ##  polynomial degree in X 
                   ny       = 2     ,   ##  polynomial degree in Y 
                   the_phis = None  ) :
        
        ## check arguments 
        assert isinstance ( nx , int ) and 0 <= nx < 100 , "'nx'-parameter is illegal: %s" % nx
        assert isinstance ( ny , int ) and 0 <= ny < 100 , "'ny'-parameter is illegal: %s" % ny
        assert isinstance ( mmax , float ) , "'mmax'-parameter is illegal: %s" % mmax

        ## the base 
        PolyBase2.__init__ ( self , name , xvar , yvar ,
                             ( nx + 1 ) * ( ny + 1 ) - 1 , the_phis )
        
        self.__phasespacex = psx  
        self.__phasespacey = psy
        
        self.__nx          = nx
        self.__ny          = ny
        self.__mmax        = mmax
        
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.PS2DPol2 (
            self.roo_name ( 'ps22_'  ) , 
            'Product of phase space factors (with polynomials) %s' % self.name ,
            self.x           ,
            self.y           ,
            self.phasespacex ,
            self.phasespacey ,
            self.mmax        ,
            self.nx          ,
            self.ny          , 
            self.phi_list    )
        
        ## save configuration
        self.config = {
            'name'     : self.name ,            
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'psx'      : self.phasespacex , 
            'psy'      : self.phasespacey ,
            'mmax'     : self.mmax , 
            'nx'       : self.nx   , 
            'ny'       : self.ny   , 
            'the_phis' : self.phis ,
            }
    
    @property 
    def mass1 ( self ) :
        """'mass1'-variable for the fit (alias for 'x' or 'xvar')"""
        return self.xvar
    
    @property 
    def mass2 ( self ) :
        """'mass2'-variable for the fit (alias for 'y' or 'yvar')"""
        return self.yvar

    @property
    def phasespacex( self ) :
        """'x-phasespace' : function for  PSPol2D-function"""
        return self.__phasespacex
    
    @property
    def phasespacey( self ) :
        """'y-phasespace' : function for  PSPol2D-function"""
        return self.__phasespacey 

    @property
    def mmax  ( self ) :
        """'mmax' :  the maximal allowed mass"""
        return self.__mmax
    
    @property
    def nx ( self ) :
        """'nx'  : order/degree of 2D-polynom in x-direction"""
        return self.__nx
    @property
    def ny ( self ) :
        """'ny' : order/degree of 2D-polynom in y-direction"""
        return self.__ny
    
models.append ( PSPol2D2_pdf ) 


# =============================================================================
## @class PSPol2D3_pdf
#  The 2D-function, that represent non-factorizeable "product" of  
#  two modulated phase-space functions.
#  It can be considered as a simpler alternative for class PSPol2D2_pdf
#  The  function is useful to describe e.g. 2D-distributions of 
#  \f$ m_{23}\f$ vs \f$m_{45}\f$ from 5-body decays. 
#
#  The function is:
#  \f[ f(x,y) = \frac{1}{2}
#   \left( \Phi^{(N^{x})}_{k,n}(x;x_{low},x_{high}) \Phi_{l,m-1}(y,y_{low},m_{max}-x) 
#        + \Phi^{(N^{y})}_{l,m}(y;y_{low},y_{high}) \Phi_{k,n-1}(x,x_{low},m_{max}-y) 
#    \right) \f]
#  where 
#  - \f$\Phi_{i,j}(z;z_{low},z_{high}\f$ are normalized phase space functions
#    for mass of \f$i\f$-particles from \f$j\f$-body decays
#  - \f$\Phi^{(N)}_{i,j}(z;z_{low},z_{high}\f$ are normalized phase space functions
#    for mass of \f$i\f$-particles from \f$j\f$-body decays, modulated by
#    1D positive benrstein polynomial of degree \f$N\f$
#  - \f$m_{max}\f$ is a maximal allowed mass for \f$x+y\f$
# 
#  @see Ostap::Models::PS2DPol3
#  @see Ostap::Math::PS2DPol3
#  @see Ostap::Math::PhaseSpacePol
#  @see Ostap::Math::PhaseSpaceNL 
#  @see Ostap::Math::Positive2D
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class PSPol2D3_pdf(PolyBase2) :
    r"""Product of phase space factors, modulated by the positive polynom in 2D
    
    The function is:
    f(x,y) = \frac{1}{2}
    \left( \Phi^{(N^{x})}_{k,n}(x;x_{low},x_{high}) \Phi_{l,m-1}(y,y_{low},m_{max}-x) 
    + \Phi^{(N^{y})}_{l,m}(y;y_{low},y_{high}) \Phi_{k,n-1}(x,x_{low},m_{max}-y) 
    \right)
    
    where 
    - \f$\Phi_{i,j}(z;z_{low},z_{high}\f$ are normalized phase space functions
    for mass of \f$i\f$-particles from \f$j\f$-body decays
    - \f$\Phi^{(N)}_{i,j}(z;z_{low},z_{high}\f$ are normalized phase space functions
    for mass of \f$i\f$-particles from \f$j\f$-body decays, modulated by
    1D positive benrstein polynomial of degree \f$N\f$
    - \f$m_{max}\f$ is a maximal allowed mass for \f$x+y\f$
    
    Note:
    - f(x,y)>=0 for whole 2D-range
    
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,   ##  the first  dimension  
                   yvar             ,   ##  the second dimension
                   psx              ,   ##  phase space in X, Ostap::Math::PhaseSpaceNL 
                   psy              ,   ##  phase space in Y, Ostap::Math::PhaseSpaceNL
                   mmax             ,   ##  max-mass 
                   nx       = 2     ,   ##  polynomial degree in X 
                   ny       = 2     ,   ##  polynomial degree in Y 
                   the_phis = None  ) :
        
        ## check arguments 
        assert isinstance ( nx , int ) and 0 <= nx < 100 , "'nx'-parameter is illegal: %s" % nx
        assert isinstance ( ny , int ) and 0 <= ny < 100 , "'ny'-parameter is illegal: %s" % ny
        assert isinstance ( mmax , float ) , "'mmax'-parameter is illegal: %s" % mmax

        ## the base 
        PolyBase2.__init__ ( self , name , xvar , yvar , nx + ny , the_phis )
        
        self.__phasespacex = psx  
        self.__phasespacey = psy
        
        self.__nx          = nx
        self.__ny          = ny
        self.__mmax        = mmax
        
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.PS2DPol3 (
            self.roo_name ( 'ps23_'  ) , 
            'Product of phase space factors (with polynomials) %s' % self.name ,
            self.x           ,
            self.y           ,
            self.phasespacex ,
            self.phasespacey ,
            self.mmax        ,
            self.nx          ,
            self.ny          , 
            self.phi_list    )
        
        ## save configuration
        self.config = {
            'name'     : self.name ,            
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'psx'      : self.phasespacex , 
            'psy'      : self.phasespacey ,
            'mmax'     : self.mmax , 
            'nx'       : self.nx   , 
            'ny'       : self.ny   , 
            'the_phis' : self.phis ,
            }
    
    @property 
    def mass1 ( self ) :
        """'mass1'-variable for the fit (alias for 'x' or 'xvar')"""
        return self.xvar
    
    @property 
    def mass2 ( self ) :
        """'mass2'-variable for the fit (alias for 'y' or 'yvar')"""
        return self.yvar

    @property
    def phasespacex( self ) :
        """'x-phasespace' : function for  PSPol2D-function"""
        return self.__phasespacex
    
    @property
    def phasespacey( self ) :
        """'y-phasespace' : function for  PSPol2D-function"""
        return self.__phasespacey 

    @property
    def mmax  ( self ) :
        """'mmax' : the maximal allowed mass"""
        return self.__mmax
    
    @property
    def nx ( self ) :
        """'nx'  : order/degree of 2D-polynom in x-direction"""
        return self.__nx
    @property
    def ny ( self ) :
        """'ny' : order/degree of 2D-polynom in y-direction"""
        return self.__ny
    
models.append ( PSPol2D3_pdf ) 


# =============================================================================
## @class PSPol2Dsym_pdf
#  The symmetric 2D-function, that represent a cross-product 
#  \f$ \Phi_{k,l}(x)\f$ and \f$ \Phi_{m,n}(y)\f$,  
#  modulated by the 2D-positive symmetric polynomial.
#  It is a "symmetrised" version of class PSPol2D
#
#  The function is:
#  \f[ f(x,y) = 
#      \Phi_{k,l}(x;x_{low}, x_{high})
#      \Phi_{k,l}(y;y_{low}, y_{high})
#       P_{N,N}(x,y) \f]
#  where  
#  - \f$ \Phi_{k,l}(x;x_{low},x_{high}) \f$ is a phase-space function,
#        \f$ y_{low}=x_{low}\f$ and \f$y_{high}=x_{high}\f$
#  - \f$ P_{N,N}(x,y) \f$ is 2D positive symmetric  Bernstein polynomial 
#
#  Clearly the function is symmetric: \f$f(x,y) = f(y,x) \f$ 
#  @see Ostap::Models::PS2DPolSym
#  @see Ostap::Math::PS2DPolSym
#  @see Ostap::Math::PhaseSpaceNL 
#  @see Ostap::Math::Positive2DSym
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class PSPol2Dsym_pdf(PolyBase2) :
    r"""Symmetric product of phase space factors, modulated by the positive polynom in 2D
    
    #  The function is:
    #  \f[ f(x,y) = 
    #      \Phi_{k,l}(x;x_{low}, x_{high})
    #      \Phi_{k,l}(y;y_{low}, y_{high})
    #       P_{N,N}(x,y) \f]
    #  where  
    #  - \f$ \Phi_{k,l}(x;x_{low},x_{high}) \f$ is a phase-space function,
    #        \f$ y_{low}=x_{low}\f$ and \f$y_{high}=x_{high}\f$
    #  - \f$ P_{N,N}(x,y) \f$ is 2D positive symmetric  Bernstein polynomial 
    #
    
    Note:
    - f(x,y)>=0 for whole 2D-range
    - f(x,y)=f(y,x) for whole 2D-range
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,   ##  the first  dimension  
                   yvar             ,   ##  the second dimension
                   ps               ,   ##  phase space in X, Ostap::Math::PhaseSpaceNL 
                   n        = 2     ,   ##  polynomial degree in Y
                   the_phis = None  ) :
        
        ## check arguments 
        assert isinstance ( n , int ) and 0 <= n < 100 , "'n'-parameter is illegal: %s" % n
        ## 

        ## the base 
        PolyBase2.__init__ ( self , name , xvar , yvar , ( n + 1 ) * ( n + 2 ) / 2  - 1 , the_phis )

        if self.xminmax() != self.yminmax() :
            logger.warning( 'PSPol2Dsym_pdf: x&y have different edges %s vs %s' % ( self.xminmax() , self.yminmax() ) )

        
        self.__phasespace = ps  

        self.__n = n 
        
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.PS2DPolSym (
            self.roo_name ( 'ps2s_'  ) , 
            'Symmetric product of phase space factors (with polynomials) %s' % self.name ,
            self.x          ,
            self.y          ,
            self.phasespace ,
            self.n          ,
            self.phi_list   )
        
        ## save configuration
        self.config = {
            'name'     : self.name ,            
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'ps'       : self.phasespace  , 
            'n'        : self.n    , 
            'the_phis' : self.phis , 
            }
        
    @property 
    def mass1 ( self ) :
        """'mass1'-variable for the fit (alias for 'x' and 'xvar')"""
        return self.xvar
    
    @property 
    def mass2 ( self ) :
        """'mass2'-variable for the fit (alias for 'y' or 'yvar')"""
        return self.yvar

    @property
    def phasespace ( self ) :
        """'phasespace' : function for PSPol2DSym-function"""
        return self.__phasespace
    
    @property
    def phasespacex( self ) :
        """'x-phasespace' : function for PSPol2Dsym-function"""
        return self.__phasespace
    
    @property
    def phasespacey( self ) :
        """'y-phasespace' : function for PSPol2Dsym-function"""
        return self.__phasespace 

    @property
    def n  ( self ) :
        """'n'  :  order/degree of 2D-polynom in x&y-directions"""
        return self.__n
    @property
    def nx ( self ) :
        """'nx' : order/degree of 2D-polynom in x-direction"""
        return self.__n
    @property
    def ny ( self ) :
        """'ny' : order/degree of 2D-polynom in y-direction"""
        return self.__n
    
        
models.append ( PSPol2Dsym_pdf ) 


# =============================================================================
## @class PSPol2D2sym_pdf
#
#  The symmetric 2D-function, that represent non-factorizeable "product" of
#  phase-space functions modulated by the 2D-positive polynomial.
#  It is a  symmetrised version of class PSPol2D2.
#  The  function is useful to describe e.g. 2D-distributions of 
#  \f$ m_{23}\f$ vs \f$m_{45}\f$ from 5-body decays. 
#
# The function is:
#  \f[ f(x,y) = \frac{1}{2}
#   \left( \Phi_{k,n}(x;x_{low},x_{high}) \Phi_{k,n-1}(y,y_{low},m_{max}-x) 
#        + \Phi_{k,n}(y;y_{low},y_{high}) \Phi_{k,n-1}(x,x_{low},m_{max}-y) 
#        \right)
#        P_{N,N}(x,y) \f]
#  where 
#  - \f$ \Phi_{i,j}(x;x_{low},x_{high}\f$ are normalized phase space function,
#    for mass of \f$i\f$-particles from \f$j\f$-body decays;
#  - \f$ y_{low}=x_{low}\f$ and \f$y_{high}=x_{high}\f$
#  - \f$ P_{N,N}(x,y) \f$ is 2D positive symmertic Bernstein polynomial
#  - \f$m_{max}\f$ is a maximal allowed mass for \f$x+y\f$
#
#  Clearly the function is symmetric \f$f(x,y) = f(y,x) \f$  
#  @see Ostap::Models::PS2DPolSym
#  @see Ostap::Math::PS2DPolSym
#  @see Ostap::Math::PhaseSpaceNL 
#  @see Ostap::Math::Positive2DSym
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class PSPol2D2sym_pdf(PolyBase2) :
    """Symmetric product of phase space factors, modulated by the positive polynom in 2D
    
    f(x,y) = PS(x) * PS(y) * Pn(x,y) 
    
    where
    - PS(x) is a phase space function  (Ostap::Math::PhaseSpaceNL)
    - Pnk(x,y) is positive symmetric non-factorizable polynom (Ostap::Math::Positive2DSym)
    -- Pn(x,y) = sum^{i=n}_{i=0}sum^{j=n}_{j=0} a^2_{ij} B^n_i(x) B^n_j(y), where 
    --- B^n_i - are Bernstein polynomials
    --- a_{ij}=a_{ji}
    
    Note:
    - f(x,y)>=0 for whole 2D-range
    - f(x,y)=f(y,x) for whole 2D-range
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,   ##  the first  dimension  
                   yvar             ,   ##  the second dimension
                   ps               ,   ##  phase space in X, Ostap::Math::PhaseSpaceNL
                   mmax             ,   ##   max-mass
                   n        = 2     ,   ##  polynomial degree in Y
                   the_phis = None  ) :
        
        ## check arguments 
        assert isinstance ( n , int ) and 0 <= n < 100 , "'n'-parameter is illegal: %s" % n
        assert isinstance ( mmax , float ) , "'mmax'-parameter is illegal: %s" % mmax
        ## 

        ## the base 
        PolyBase2.__init__ ( self , name , xvar , yvar , ( n + 1 ) * ( n + 2 ) / 2  - 1 , the_phis )

        if self.xminmax() != self.yminmax() :
            logger.warning( 'PSPol2Dsym_pdf: x&y have different edges %s vs %s' % ( self.xminmax() , self.yminmax() ) )

        
        self.__phasespace = ps  

        self.__n    = n 
        self.__mmax = mmax 
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.PS2DPol2Sym (
            self.roo_name ( 'ps22s_'  ) , 
            'Symmetric product of phase space factors (with polynomials) %s' % self.name ,
            self.x          ,
            self.y          ,
            self.phasespace ,
            self.mmax       , 
            self.n          ,
            self.phi_list   )
        
        ## save configuration
        self.config = {
            'name'     : self.name ,            
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'ps'       : self.phasespace  , 
            'mmax'     : self.mmax , 
            'n'        : self.n    , 
            'the_phis' : self.phis , 
            }
        
    @property 
    def mass1 ( self ) :
        """'mass1'-variable for the fit (alias for 'x' or 'xvar')"""
        return self.xvar
    
    @property 
    def mass2 ( self ) :
        """'mass2'-variable for the fit (alias for 'y' or 'yvar')"""
        return self.yvar

    @property
    def phasespace ( self ) :
        """'phasespace' : function for PSPol2DSym-function"""
        return self.__phasespace
    
    @property
    def phasespacex( self ) :
        """'x-phasespace' : function for PSPol2Dsym-function"""
        return self.__phasespace
    
    @property
    def phasespacey( self ) :
        """'y-phasespace' : function for PSPol2Dsym-function"""
        return self.__phasespace 

    @property
    def mmax  ( self ) :
        """'mmax' : the maximal allowed mass"""
        return self.__mmax
    
    @property
    def n  ( self ) :
        """'n'  : order/degree of 2D-polynom in x&y-directions"""
        return self.__n
    @property
    def nx ( self ) :
        """'nx' : order/degree of 2D-polynom in x-direction"""
        return self.__n
    @property
    def ny ( self ) :
        """'ny' : order/degree of 2D-polynom in y-direction"""
        return self.__n
    
        
models.append ( PSPol2D2sym_pdf ) 



# =============================================================================
## @class PSPol2D3sym_pdf
#  
#  The symmetric 2D-function, that represent non-factorizeable "product" of  
#  two modulated phase-space functions.
#  It is a  symmetrized version of PSPol2D3_pdf 
#  The  function is useful to describe e.g. 2D-distributions of 
#  \f$ m_{23}\f$ vs \f$m_{45}\f$ from 5-body decays. 
#
#  The function is:
#  \f[ f(x,y) = \frac{1}{2}
#   \left( \Phi^{(N)}_{k,n}(x;x_{low},x_{high}) \Phi_{k,n-1}(y,y_{low},m_{max}-x) 
#        + \Phi^{(N)}_{k,n}(y;y_{low},y_{high}) \Phi_{k,n-1}(x,x_{low},m_{max}-y) 
#   \right) \f]
#  where 
#  - \f$ \Phi_{i,j}(z;z_{low},z_{high}\f$ are normalized phase space functions
#    for mass of \f$i\f$-particles from \f$j\f$-body decays
#  - \f$\Phi^{(N)}_{i,j}(z;z_{low},z_{high}\f$ are normalized phase space functions
#    for mass of \f$i\f$-particles from \f$j\f$-body decays, modulated by
#        1D positive benrstein polynomial of degree \f$N\f$
#  - \f$m_{max}\f$ is a maximal allowed mass for \f$x+y\f$
#
#  Clearly the function is symmetric:  \f$f(y,x)=f(x,y)\f$                             
#  @see Ostap::Models::PS2DPolSym
#  @see Ostap::Math::PS2DPolSym
#  @see Ostap::Math::PhaseSpaceNL 
#  @see Ostap::Math::Positive2DSym
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class PSPol2D3sym_pdf(PolyBase2) :
    """Symmetric product of phase space factors, modulated by the positive polynom in 2D
    
    f(x,y) = PS(x) * PS(y) * Pn(x,y) 
    
    where
    - PS(x) is a phase space function  (Ostap::Math::PhaseSpaceNL)
    - Pnk(x,y) is positive symmetric non-factorizable polynom (Ostap::Math::Positive2DSym)
    -- Pn(x,y) = sum^{i=n}_{i=0}sum^{j=n}_{j=0} a^2_{ij} B^n_i(x) B^n_j(y), where 
    --- B^n_i - are Bernstein polynomials
    --- a_{ij}=a_{ji}
    
    Note:
    - f(x,y)>=0 for whole 2D-range
    - f(x,y)=f(y,x) for whole 2D-range
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,   ##  the first  dimension  
                   yvar             ,   ##  the second dimension
                   ps               ,   ##  phase space in X, Ostap::Math::PhaseSpaceNL
                   mmax             ,   ##   max-mass
                   n        = 2     ,   ##  polynomial degree in Y
                   the_phis = None  ) :
        
        ## check arguments 
        assert isinstance ( n , int ) and 0 <= n < 100 , "'n'-parameter is illegal: %s" % n
        assert isinstance ( mmax , float ) , "'mmax'-parameter is illegal: %s" % mmax
        ## 

        ## the base 
        PolyBase2.__init__ ( self , name , xvar , yvar , n , the_phis )

        if self.xminmax() != self.yminmax() :
            logger.warning( 'PSPol2Dsym_pdf: x&y have different edges %s vs %s' % ( self.xminmax() , self.yminmax() ) )

        
        self.__phasespace = ps  

        self.__n    = n 
        self.__mmax = mmax 
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.PS2DPol3Sym (
            self.roo_name ( 'ps23s_'  ) , 
            'Symmetric product of phase space factors (with polynomials) %s' % self.name ,
            self.x          ,
            self.y          ,
            self.phasespace ,
            self.mmax       , 
            self.n          ,
            self.phi_list   )
        
        ## save configuration
        self.config = {
            'name'     : self.name ,            
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'ps'       : self.phasespace  , 
            'mmax'     : self.mmax , 
            'n'        : self.n    , 
            'the_phis' : self.phis , 
            }
        
    @property 
    def mass1 ( self ) :
        """'mass1'-variable for the fit (alias for 'x' and 'xvar')"""
        return self.xvar
    
    @property 
    def mass2 ( self ) :
        """'mass2'-variable for the fit (alias for 'y' and 'yvar')"""
        return self.yvar

    @property
    def phasespace ( self ) :
        """'phasespace' : function for PSPol2DSym-function"""
        return self.__phasespace
    
    @property
    def phasespacex( self ) :
        """'x-phasespace' : function for PSPol2Dsym-function"""
        return self.__phasespace
    
    @property
    def phasespacey( self ) :
        """'y-phasespace' : function for PSPol2Dsym-function"""
        return self.__phasespace 

    @property
    def mmax  ( self ) :
        """'mmax' : the maximal allowed mass"""
        return self.__mmax
    
    @property
    def n  ( self ) :
        """'n'  : order/degree of 2D-polynom in x&y-directions"""
        return self.__n
    @property
    def nx ( self ) :
        """'nx' : order/degree of 2D-polynom in x-direction"""
        return self.__n
    @property
    def ny ( self ) :
        """'ny' : order/degree of 2D-polynom in y-direction"""
        return self.__n
    
        
models.append ( PSPol2D3sym_pdf ) 


# =============================================================================
## @class ExpoPSPol2D_pdf
#  Product of phase space factors, modulated by positiev polynomial in 2D 
#  \f$  f(x,y) = exp(\tau x) \times \Phi (y) \times P^+(x,y) \f$,
#  where \f$ P^+(x,y)\f$ denotes the positive polynomial,
#  @see Ostap::Models::ExpoPS2DPol
#  @see Ostap::Math::ExpoPS2DPol
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class ExpoPSPol2D_pdf(PolyBase2) :
    """Product of exponential and phase space factor,
    modulated by the positive polynom in 2D

    f(x,y) = exp(tau*x) * PS(y) * Pnk(x,y)
    where
    - PS (y) is a phase space function for y-axis (Ostap::Math::PhaseSpaceNL)
    - Pnk(x,y) is positive non-factorizable polynom
    
    Pnk(x,y) = sum^{i=n}_{i=0}sum{j=k}_{j=0} a^2_{ij} B^n_i(x) B^k_j(y)
    where:
    - B^n_i - are Bernstein polynomials
    
    Note:
    - f(x,y)>=0 for whole 2D-range
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             , ##  the first  dimension  
                   yvar             , ##  the second dimension
                   psy      = None  , ##  phase space in Y, Ostap::Math::PhaseSpaceNL 
                   nx       = 2     , ##  polynomial degree in X 
                   ny       = 2     , ##  polynomial degree in Y 
                   tau      = None  , ##  the exponent 
                   the_phis = None  ) :
        
        ## check arguments 
        assert isinstance ( nx , int ) and 0 <= nx < 100 , "'nx'-parameter is illegal: %s" % nx
        assert isinstance ( ny , int ) and 0 <= ny < 100 , "'ny'-parameter is illegal: %s" % ny

        ## the base 
        PolyBase2.__init__ ( self , name , xvar , yvar , ( nx + 1 ) * ( ny + 1 ) - 1 , the_phis = the_phis )
        
        limits_tau = () 
        if self.xminmax() : 
            mn , mx     = self.xminmax()
            mmax        = max ( abs ( mn ) , abs ( mx ) )
            limits_tau  = -500. / mmax ,  500. / mmax
            
        ## the exponential slope
        self.__tau  = self.make_var ( tau              ,
                                      "tau_%s"  % name ,
                                      "tau(%s)" % name ,
                                      None , 0 , *limits_tau )
        #
        self.__phasespace = psy

        self.__nx = nx
        self.__ny = ny
        
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.ExpoPS2DPol (
            self.roo_name ( 'eps2_'  ) , 
            'Phase space times exponential (with polynomials) %s' % self.name ,
            self.x           ,
            self.y           ,
            self.tau         ,
            self.phasespacey , 
            self.nx          ,
            self.ny          , 
            self.phi_list    )
        
        ## save configuration
        self.config = {
            'name'     : self.name ,            
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'psy'      : self.phasespacey , 
            'nx'       : self.nx   , 
            'ny'       : self.ny   , 
            'tau'      : self.tau  , 
            'the_phis' : self.phis ,
            }

    @property 
    def mass1 ( self ) :
        """'mass1'-variable for the fit (alias for 'x' and  'xvar')"""
        return self.xvar
    
    @property 
    def mass2 ( self ) :
        """'mass2'-variable for the fit (alias for 'y' and 'yvar')"""
        return self.yvar

    @property
    def tau ( self ) :
        """'tau' : the exponential slope for  x-dimension"""
        return   self.__tau 
    @tau.setter
    def tau ( self , value ) :
        self.set_value ( self.__tau , value )
    
    @property
    def phasespace ( self ) :
        """'phasespace' : function for PSPol2DSym-function"""
        return self.__phasespace
        
    @property
    def phasespacey( self ) :
        """'y-phasespace' : function for PSPol2Dsym-function"""
        return self.phasespace
    
    @property
    def nx ( self ) :
        """'nx' : order/degree of 2D-polynom in x-direction"""
        return self.__nx
    @property
    def ny ( self ) :
        """'ny' : order/degree of 2D-polynom in y-direction"""
        return self.__ny
    

models.append ( ExpoPSPol2D_pdf ) 
# =============================================================================
## @class ExpoPol2D_pdf
#  Product of phase space factors, modulated by positive polynomial in 2D 
#  \f$  f(x,y) = exp(\tau_x x) \times exp ( \tau_y y) \times P^+(x,y) \f$,
#  where \f$ P^+(x,y)\f$ denotes the positive polynomial,
#  @see Ostap::Models::Expo2DPol
#  @see Ostap::Math::Expo2DPol
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class ExpoPol2D_pdf(PolyBase2) :
    """Product of exponential factors
    modulated by the positive polynom in 2D

    f(x,y) = exp(tau_x*x) * exp(tau_y*y) * Pnk(x,y)
    where
    - Pnk(x,y) is positive non-factorizable polynom
    
    Pnk(x,y) = sum^{i=n}_{i=0}sum{j=k}_{j=0} a^2_{ij} B^n_i(x) B^k_j(y)
    where:
    - B^n_i - are Bernstein polynomials
    
    Note:
    - f(x,y)>=0 for whole 2D-range

    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,   ##  the first  dimension  
                   yvar             ,   ##  the second dimension
                   nx   = 2         ,   ##  polynomial degree in X 
                   ny   = 2         ,   ##  polynomial degree in Y
                   taux = None      ,   ##  the exponent in X 
                   tauy = None      ,   ##  the exponent in Y
                   the_phis = None  ) : 
        
        ## check arguments 
        assert isinstance ( nx , int ) and 0 <= nx < 100 , "'nx'-parameter is illegal: %s" % nx
        assert isinstance ( ny , int ) and 0 <= ny < 100 , "'ny'-parameter is illegal: %s" % ny

        PolyBase2.__init__ ( self , name , xvar , yvar ,
                             ( nx + 1 ) * ( ny + 1 ) - 1 , the_phis )
        
        limits_taux = () 
        if self.xminmax() : 
            mn , mx     = self.xminmax()
            mmax        = max ( abs ( mn ) , abs ( mx ) )
            limits_taux = -500. / mmax ,  500. / mmax

        limits_tauy = () 
        if self.yminmax() : 
            mn , mx     = self.yminmax()
            mmax        = max ( abs ( mn ) , abs ( mx ) )
            limits_tauy = -500. / mmax ,  500. / mmax

        self.__nx = nx 
        self.__ny = ny
        #
        ## the exponential slopes
        #
        self.__taux  = self.make_var ( taux              ,
                                       "taux_%s"  % name ,
                                       "taux(%s)" % name ,
                                       None , 0 , *limits_taux )
        #
        self.__tauy  = self.make_var ( tauy              ,
                                       "tauy_%s"  % name ,
                                       "tauy(%s)" % name ,
                                       None , 0 , *limits_tauy )
            
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.Expo2DPol (
            self.roo_name ( 'exp2_'  ) , 
            'Exponentials (with polynomials) %s' % self.name ,
            self.x        ,
            self.y        ,
            self.taux     ,
            self.tauy     ,
            nx            ,
            ny            , 
            self.phi_list )
        
        ## save configuration
        self.config = {
            'name'     : self.name ,            
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'nx'       : self.nx   , 
            'ny'       : self.ny   , 
            'taux'     : self.taux , 
            'tauy'     : self.tauy , 
            'the_phis' : self.phis ,
            }

    @property
    def taux ( self ) :
        """'tau-x' : the exponential slope for x-dimension"""
        return   self.__taux 
    @taux.setter
    def taux ( self , value ) :
        self.set_value ( self.__taux , value )

    @property
    def tauy ( self ) :
        """'tau-y'' : the exponential slope for y-dimension"""
        return   self.__tauy
    @tauy.setter
    def tauy ( self , value ) :
        self.set_value ( self.__tauy , value )
    
    @property
    def nx ( self ) :
        """'nx' : order/degree of 2D-polynom in x-direction"""
        return self.__nx
    @property
    def ny ( self ) :
        """'ny' : order/degree of 2D-polynom in y-direction"""
        return self.__ny
    
        
models.append ( ExpoPol2D_pdf ) 
# =============================================================================
## @class ExpoPol2Dsym_pdf
#  Product of phase space factors, modulated by positiev polynomial in 2D 
#  \f$  f(x,y) = exp(\tau x) \times exp ( \tau y) \times P^+_{sym}(x,y) \f$,
#  where \f$ P^+_{sym}(x,y)\f$ denotes the symmetric positive polynomial,
#  @see Ostap::Models::Expo2DPolSym
#  @see Ostap::Math::Expo2DPolSym
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class ExpoPol2Dsym_pdf(PolyBase2) :
    """Symmetric product of exponential factors modulated by the positive polynom in 2D
    
    f(x,y) = exp(tau*x) * exp(tau*y) * Sn(x,y)
    where
    - Sn(x,y) is positive symmetric non-factorizable polynom
    
    Sn(x,y) = sum^{i=n}_{i=0}sum{j=n}_{j=0} a^2_{ij} B^n_i(x) B^n_j(y)
    where:
    - B^n_i - are Bernstein polynomials
    - a_{ij}=a_{ji}
    
    Note:
    - f(x,y)>=0 for whole 2D-range
    - f(x,y)=f(y,x) 
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,   ## the first  dimension  
                   yvar             ,   ## the second dimension
                   n    = 2         ,   ## polynomial degree in X and Y
                   tau  = None      ,   ## the exponent 
                   the_phis = None  ) : 
        
        ## check arguments 
        assert isinstance ( n , int ) and 0 <= n < 100 , "'n'-parameter is illegal: %s" % n
        ## 
        PolyBase2.__init__ ( self , name , xvar , yvar ,
                             ( n + 1 ) * ( n + 2 ) / 2 - 1 , the_phis )
        
        if self.xminmax() != self.yminmax() :
            logger.warning( 'PSPol2Dsym_pdf: x&y have different edges %s vs %s' % ( self.xminmax() , self.yminmax() ) )

             
        limits_tau = () 
        if self.xminmax() : 
            mn , mx    = self.xminmax()
            mmax       = max ( abs ( mn ) , abs ( mx ) )
            limits_tau = -500. / mmax ,  500. / mmax

        self.__n = n 
        #
        ## the exponential slopes
        #
        self.__tau  = self.make_var ( tau              ,
                                      "tau_%s"  % name ,
                                      "tau(%s)" % name ,
                                      None , 0 , *limits_tau )
            
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.Expo2DPolSym (
            self.roo_name ( 'exp2s_'  ) , 
            'Symmetric exponentials (with polynomials) %s' % self.name ,
            self.x        ,
            self.y        ,
            self.tau      ,
            self.n        ,
            self.phi_list )
        
        ## save configuration
        self.config = {
            'name'     : self.name ,            
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'n'        : self.n    , 
            'tau'      : self.tau  , 
            'the_phis' : self.phis 
            }
        
    @property
    def tau ( self ) :
        """'tau' : the exponential slope for x&y-dimensions"""
        return   self.__tau 
    @tau.setter
    def tau ( self , value ) :
        self.set_value ( self.__tau , value )

    @property
    def taux ( self ) :
        """'tau-x' : the exponential slope for x-dimension"""
        return   self.__tau
    @taux.setter
    def taux ( self , value ) :
        self.set_value ( self.__tau , value )

    @property
    def tauy ( self ) :
        """'tau-y' : the exponential slope for y-dimension"""
        return   self.__tau
    @tauy.setter
    def tauy ( self , value ) :
        self.set_value ( self.__tau , value )
    
    @property
    def n  ( self ) :
        """'n'  : order/degree of 2D-polynom in x&y-directions"""
        return self.__n
    @property
    def nx ( self ) :
        """'nx'' : order/degree of 2D-polynom in x-direction"""
        return self.__n
    @property
    def ny ( self ) :
        """'ny' : order/degree of 2D-polynom in y-direction"""
        return self.__n
    

models.append ( ExpoPol2Dsym_pdf ) 
# =============================================================================
## @class Spline2D_pdf
#  positive spline in 2D:
#  \f{displaymath}  f(x,y) = \sum^{i=n}_{i=0}\sum{j=k}_{j=0} a^2_{ij} M^n_i(x) M^k_j(y) \f},
#  where \f$ B^n_i(x)\f$ denotes the M-splines  
#  @see Ostap::Models::Spline2D
#  @see Ostap::Math::Spline2D
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class Spline2D_pdf(PolyBase2) :
    """Positive non-factorizable spline in 2D
    
    f(x,y) = sum_i sum_j a^2_{i,j} Nx_i(x) * Ny_j(y),
    where
    - Nx_i and Ny_j are normailzed B-splines for x and y-axes

    Note:
    - f(x,y)>=0 for whole 2D-range
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,   ##  the first  dimension  
                   yvar             ,   ##  the second dimension
                   spline           ,   ## the spline: Ostap.Math.PositiveSpline2D 
                   the_phis = None  ) : ## 
        
        PolyBase2.__init__ ( self , name , xvar , yvar , spline.npars() , the_phis ) 
        
        self.__spline = spline

        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.Spline2D (
            self.roo_name ( 'spline2_'  ) , 
            'Positive 2D spline %s' % self.name ,
            self.x        ,
            self.y        ,
            self.spline   ,
            self.phi_list )
        
        ## save configuration
        self.config = {
            'name'     : self.name   ,            
            'xvar'     : self.xvar   ,
            'yvar'     : self.yvar   ,
            'spline'   : self.spline , 
            'the_phis' : self.phis   
            }

    @property
    def spline ( self ) :
        """'spline'-function for Spline2D PDF"""
        return self.__spline

models.append ( Spline2D_pdf ) 
# =============================================================================
## @class Spline2Dsym_pdf
#  symmetric positive spline in 2D:
#  \f{displaymath}  f(x,y) = \sum^{i=n}_{i=0}\sum{j=k}_{j=0} a^2_{ij} M^n_i(x) M^k_j(y) \f},
#  where \f$ B^n_i(x)\f$ denotes the M-splines  
#  @see Ostap::Models::Spline2DSym
#  @see Ostap::Math::Spline2DSym
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class Spline2Dsym_pdf(PolyBase2) :
    """SYMMETRIC positive non-factorizable spline in 2D
    
    f(x,y) = sum_i sum_j a^2_{i,j} N_i(x) * N_j(y),
    where
    - N_i are normailzed B-splines 

    Note:
    - f(x,y)>=0     for whole 2D-range
    - f(x,y)=f(y,x) for whole 2D-range 
    """
    def __init__ ( self              ,
                   name              ,
                   xvar              ,   ## the first  dimension  
                   yvar              ,   ## the second dimension
                   spline            ,   ## the spline: Ostap.Math.PositiveSpline2DSym
                   the_phis = None ) :

        PolyBase2.__init__ ( self , name , xvar , yvar , spline.npars() , the_phis )        

        if self.xminmax() != self.yminmax() :
            logger.warning( 'Spline2Dsym_pdf: x&y have different edges %s vs %s' % ( self.xminmax() , self.yminmax() ) )
             
        self.__spline = spline
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.Spline2DSym (
            self.roo_name ( 'spline2s_'  ) , 
            'Positive symmetric 2D spline %s' % self.name ,
            self.x        ,
            self.y        ,
            self.spline   ,
            self.phi_list )
        
        ## save configuration
        self.config = {
            'name'     : self.name   ,            
            'xvar'     : self.xvar   ,
            'yvar'     : self.yvar   ,
            'spline'   : self.spline , 
            'the_phis' : self.phis   
            }
        
    @property
    def spline ( self ) :
        """'spline'-function for Spline2Dsym PDF"""
        return self.__spline

models.append ( Spline2Dsym_pdf )


# =============================================================================        
## @class Gauss2D_pdf
#  Simple 2D gaussian function
#  @see Ostap::Models::Gauss2D 
#  @see Ostap::Math::Gauss2D 
class Gauss2D_pdf(PDF2) :
    """ Simple 2D gaussian function
    - see Ostap.Models.Gauss2D 
    - see Ostap.Math.Gauss2D 
    """
    
    ## constructor
    def __init__ ( self            ,
                   name            ,   ## the name 
                   xvar            ,   ## the variable
                   yvar            ,   ## the variable
                   muX      = None ,
                   muY      = None ,
                   sigmaX   = None ,
                   sigmaY   = None ,
                   theta    = None ) : 
    
        ## initialize the base class 
        PDF2.__init__ (  self , name , xvar , yvar )
        
        sx_lims = ()
        mx_lims = ()
        sy_lims = ()
        my_lims = ()

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

            
        self.__muX    = self.make_var ( muX  ,
                                        'mu_x_%s'     % self.name ,
                                        '#mu_{x}(%s)' % self.name ,
                                        None , *mx_lims )
        
        self.__muY    = self.make_var ( muY  ,
                                        'mu_y_%s'     % self.name ,
                                        '#mu_{y}(%s)' % self.name ,
                                        None , *my_lims )
        
        
        self.__sigmaX = self.make_var ( sigmaX  ,
                                        'sigma_x_%s'     % self.name ,
                                        '#sigma_{x}(%s)' % self.name ,
                                        None    , *sx_lims )
        
        self.__sigmaY = self.make_var ( sigmaY  ,
                                        'sigma_y_%s'     % self.name ,
                                        '#sigma_{y}(%s)' % self.name ,
                                        None    , *sy_lims )
        
        self.__theta   = self.make_var ( theta  ,
                                         'theta_%s'   % self.name ,
                                         '#theta(%s)' % self.name ,
                                         None   , -10 , +10  )
        
        ## make PDF
        self.pdf = Ostap.Models.Gauss2D (
            self.roo_name ( 'gauss2d'  ) , 
            'Gauss2D %s' % self.name     ,
            self.x      ,
            self.y      ,
            self.muX    ,
            self.muY    ,
            self.sigmaX ,
            self.sigmaY ,
            self.theta  ,            
            )
        
        ## save configuration
        self.config = {
            'name'     : self.name   ,            
            'xvar'     : self.xvar   ,
            'yvar'     : self.yvar   ,            
            'muX'      : self.muX    ,
            'muY'      : self.muY    ,
            'sigmaX'   : self.sigmaX ,
            'sigmaY'   : self.sigmaY ,
            'theta'    : self.theta
            }
        
    @property
    def muX ( self ) :
        """x-locaiton for 2D gaussian"""
        return self.__muX
    @muX.setter
    def muX ( self , value ) :
        self.set_value ( self.__muX , value )

    @property
    def muY ( self ) :
        """y-locaiton for 2D gaussian"""
        return self.__muY
    @muY.setter
    def muY ( self , value ) :
        self.set_value ( self.__muY , value )
        
    @property
    def sigmaX ( self ) :
        """'sigma-X' for 2D gaussian"""
        return self.__sigmaX
    @sigmaX.setter
    def sigmaX ( self , value ) :
        self.set_value ( self.__sigmaX , value )

    @property
    def sigmaY ( self ) :
        """'sigma-Y' for 2D gaussian"""
        return self.__sigmaY
    @sigmaY.setter
    def sigmaY ( self , value ) :
        self.set_value ( self.__sigmaY , value )
        
    @property
    def theta ( self ) :
        """'theta' :  rotation for 2D gaussian"""
        return self.__theta
    @theta.setter
    def theta ( self , value ) :
        self.set_value ( self.__theta , value )


# =============================================================================        
## @class RooKeys2D_pdf
#  Trivial Ostap wrapper for the native <code>RooNDKeysPdf</code> from RooFit
#  @code
#  xvar = ...
#  yvar = ...
#  pdf  = RooKeys2D_pdf ( 'P4' , xvar , yvar ,  data  ) ;
#  @endcode
#  @see RooNDKeysPdf
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2019-04-27
class RooKeys2D_pdf(PDF2) :
    """Trivial Ostap wrapper for the native RooNDKeysPdf from RooFit
    - see ROOT.RooKeysPdf
    >>> xvar = ...
    >>> yvar = ...
    >>> pdf  = RooKeys2D_pdf ( 'Keys' , xvar , yvar,  data  )
    """
    ## constructor
    def __init__ ( self           ,
                   name           ,   ## the name 
                   xvar           ,   ## the variable
                   yvar           ,   ## the variable
                   data           ,   ## data set 
                   options = 'ma' ,   ## options 
                   rho     = 1    ,   ## global scale 
                   nsigma  = 3    ,   ##
                   rotate  = True ,   
                   sort    = True ) : 
        
        ## initialize the base class 
        PDF2.__init__ (  self , name , xvar , yvar )

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
        
        ## create PDF
        self.pdf = ROOT.RooNDKeysPdf (
            "rookeys2_%s"        % name ,
            "RooNDKeysPdf(%s,2)" % name ,
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
            'data'    : self.data    ,            
            'options' : self.options ,            
            'rho'     : self.rho     ,            
            'nsigma'  : self.nsigma  ,            
            'rotate'  : self.rotate  ,            
            'sort'    : self.sort    ,            
            }

    @property
    def data   ( self ) :
        """'data' : the actual data set for RooNDKeysPdf"""
        return self.__data
    @property
    def options ( self ) :
        """'options' : 'options'-string for RooNDKeysPdf"""
        return self.__mirror        
    @property
    def rho    ( self )  :
        """'rho' :  'rho'-parameter for RooNDKeysPdf"""
        return self.__rho 
    @property
    def rotate    ( self )  :
        """'rotate' : `'rotate'-flag for RooNDKeysPdf"""
        return self.__rotate
    @property
    def sort      ( self )  :
        """'sort' : 'sort''-flag for RooNDKeysPdf"""
        return self.__sort 
    @property
    def nsigma    ( self )  :
        """'nsigma' : 'nsigma' parameter for RooNDKeysPdf"""
        return self.__nsigma 

# =============================================================================
# some tiny decoration of underlying classes 
# =============================================================================
root_major = root_info.major 
def _2d_get_pars_ ( self ) :
    """
    Get parameters of underlying positive Berstein polynomial
    """
    if   hasattr ( self , 'pdf'      ) :
        return _2d_get_pars_ ( self.pdf         )        
    elif hasattr ( self , 'function' ) :
        return _2d_get_pars_ ( self.function () )
    elif hasattr ( self , 'positive' ) :
        return _2d_get_pars_ ( self.positive () )
    elif hasattr ( self , 'polynom'  ) :
        return _2d_get_pars_ ( self. polynom () )
    elif hasattr ( self , 'bernstein' ) :
        b = self.bernstein()
        m = ROOT.TMatrix ( b.nX() + 1 , b.nY() + 1 )
        for i in range ( 0 , b.nX() + 1 ) :
            for j in range ( 0 , b.nY() + 1 ) :
                
                if  root_major < 6 : m[i][j] = b.par(i,j)
                else               : m[i, j] = b.par(i,j)
                    
        return m 
        
    return ROOT.TMatrix()


for t in ( PolyPos2D_pdf    ,
           PolyPos2Dsym_pdf ,
           PSPol2D_pdf      ,
           PSPol2Dsym_pdf   ,
           ExpoPSPol2D_pdf  ,
           ExpoPol2D_pdf    ,
           ExpoPol2Dsym_pdf ) :

    t.pars = _2d_get_pars_ 
                
# =============================================================================


# ==============================================================================
## Easy creation of  2D function for background
#  @code
#  xvar = ...
#  yvar = ...
#  bkg  = make_B2D ( 'BB' , xvar , yvar , -1 , -1 ) ## create PolyPos2D_pdf 
#  bkg  = make_B2D ( 'BB' , xvar , yvar ,  1 ,  1 ) ## create ExpoPol2D_pdf 
#  bkg  = make_B2D ( 'BB' , xvar , yvar ,  1 , -1 ) ## create ExpoPol2D_pdf, fix tau_y 
#  bkg  = make_B2D ( 'BB' , xvar , yvar , -1 ,  1 ) ## create ExpoPol2D_pdf, fix tau_x  
#  bkg  = make_B2D ( 'BB' , xvar , yvar ,  0 ,  0 ) ## create Flat2D 
#  endcode
def make_B2D ( name , xvar , yvar , nx , ny ) :
    """Easy creation of  2D function for background
    >>> xvar = ...
    >>> yvar = ...
    >>> bkg  = make_B2D ( 'BB' , xvar , yvar , -1 , -1 ) ## create PolyPos2D_pdf
    >>> bkg  = make_B2D ( 'BB' , xvar , yvar ,  1 ,  1 ) ## create ExpoPol2D_pdf 
    >>> bkg  = make_B2D ( 'BB' , xvar , yvar ,  1 , -1 ) ## create ExpoPol2D_pdf, fix tau_y 
    >>> bkg  = make_B2D ( 'BB' , xvar , yvar , -1 ,  1 ) ## create ExpoPol2D_pdf, fix tau_x  
    >>> bkg  = make_B2D ( 'BB' , xvar , yvar ,  0 ,  0 ) ## create Flat2D     
    """
    
    if   0 == nx and 0 == ny :
        return Flat2D        ( name = name , xvar = xvar , yvar = yvar )
    elif 0 >= nx and 0 >= ny : 
        return PolyPos2D_pdf ( name = name , xvar = xvar , yvar = yvar , nx = abs ( nx ) , ny = abs ( ny ) )     

    fun2 = ExpoPol2D_pdf     ( name = name , xvar = xvar , yvar = yvar , nx = abs ( nx ) , ny = abs ( ny ) )
    if 0 > nx : fun2.taux.fix ( 0 )
    if 0 > ny : fun2.tauy.fix ( 0 )
    
    return fun2

# ==============================================================================
## Easy creation of symmetric 2D function for background
#  @code
#  xvar = ...
#  yvar = ...
#  bkg  = make_B2Dsym ( 'BB' , xvar , yvar , -1 ) ## create PolyPol2Dsym_pdf 
#  bkg  = make_B2Dsym ( 'BB' , xvar , yvar ,  1 ) ## create ExpoPol2Dsym_pdf 
#  bkg  = make_B2Dsym ( 'BB' , xvar , yvar ,  0 ) ## create Flat2D 
#  endcode
def make_B2Dsym ( name , xvar , yvar , n ) :
    """Easy creation of symmetric 2D function for background
    >>> xvar = ...
    >>> yvar = ...
    >>> bkg  = make_B2Dsym ( 'BB' , xvar , yvar , -1 ) ## create PolyPol2Dsym_pdf 
    >>> bkg  = make_B2Dsym ( 'BB' , xvar , yvar ,  1 ) ## create ExpoPol2Dsym_pdf 
    >>> bkg  = make_B2Dsym ( 'BB' , xvar , yvar ,  0 ) ## create Flat2D 
    """
    
    if   0 == n :
        return Flat2D           ( name = name , xvar = xvar , yvar = yvar )
    elif 0 >= n : 
        return PolyPos2DSym_pdf ( name = name , xvar = xvar , yvar = yvar , n = abs ( n ) )     

    fun2 = ExpoPol2Dsym_pdf     ( name = name , xvar = xvar , yvar = yvar , n = abs ( n ) )
    
    return fun2



# =============================================================================
if '__main__' == __name__ : 
         
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
