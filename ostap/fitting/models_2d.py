#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
# =============================================================================
## @file models_2d.py
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
    'PSPol2Dsym_pdf'  , ## Symmetric product of phase spaces, modulated with 2D polynomial
    'ExpoPSPol2D_pdf' , ## Exponential times  phase space times positive 2D-polynomial
    'ExpoPol2D_pdf'   , ## Product of exponents times positive 2D-polynomial
    'ExpoPol2Dsym_pdf', ## Symmetric version of above
    ##
    'Spline2D_pdf'    , ## 2D generic   positive spline 
    'Spline2Dsym_pdf' , ## 2D symmetric positive spline 
    )
# =============================================================================
import ROOT, math
from   ostap.core.core     import cpp, Ostap, iszero  
from   ostap.fitting.basic import makeVar, PDF2 
# =============================================================================
from   ostap.logger.logger     import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.models_2d' )
else                       : logger = getLogger ( __name__                  )
# =============================================================================
models = []
# =============================================================================
## @class PolyPos2D_pdf
#  positive polynomial in 2D:
#  \f$  f(x,y) = \sum^{i=n}_{i=0}\sum{j=k}_{j=0} a^2_{\ij} B^n_i(x) B^k_j(y) \f$,
#  where \f$ B^n_i(x)\f$ denotes the basic bersntein polynomial 
#  @see Ostap::Models::Poly2DPositive
#  @see Gaudi::Math::Poly2DPositive
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class PolyPos2D_pdf(PDF2) :
    """Positive (non-factorizable!) polynomial in 2D:

    f(x,y) = sum^{i=n}_{i=0}sum{j=k}_{j=0} a^2_{ij} B^n_i(x) B^k_j(y)
    
    where  B^n_i - are Bernstein polynomials
    
    Note:
    
    - f(x,y)>=0 for whole 2D-range
    """
    def __init__ ( self             ,
                   name             ,
                   x                ,   ##  the first  dimension  
                   y                ,   ##  the second dimension
                   nx = 2           ,   ##  polynomial degree in X 
                   ny = 2           ) : ##  polynomial degree in Y 

        PDF2.__init__ ( self , name , x , y ) 
        
        # 
        self.makePhis ( ( nx + 1 ) * ( ny + 1 ) - 1 ) 
        #
            
        #
        ## finally build PDF 
        #
        self.pdf = cpp.Analysis.Models.Poly2DPositive (
            'p2Dp_%s'            % name ,
            'Poly2DPositive(%s)' % name ,
            self.x        ,
            self.y        ,
            nx            ,
            ny            , 
            self.phi_list )
        
models.append ( PolyPos2D_pdf ) 
# =============================================================================
## @class PolyPos2Dsym_pdf
#  Positive symetric polynomial in 2D:
#  \f$  f(x,y) = \sum^{i=n}_{i=0}\sum{j=n}_{j=0} a^2_{\ij} B^n_i(x) B^n_j(y) \f$,
#  where \f$ B^n_i(x)\f$ denotes the basic bersntein polynomial and
#  \f$a_{ij} = a_{ji}\f$
#  @see Ostap::Models::Poly2DSymPositive
#  @see Gaudi::Math::Poly2DSymPositive
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class PolyPos2Dsym_pdf(PDF2) :
    """Positive (non-factorizable!) SYMMETRIC polynomial in 2D:
    
    f(x,y) = sum^{i=n}_{i=0}sum{j=n}_{j=0} a^2_{ij} B^n_i(x) B^n_j(y)
    
    where:
    - B^n_i - are Bernstein polynomials
    - a_{ij} = a_{ji}
    
    Note:
    - f(x,y)>=0 for whole 2D-range
    - f(x,y) = f(y,x)
    """
    def __init__ ( self             ,
                   name             ,
                   x                ,   ##  the first  dimension  
                   y                ,   ##  the second dimension
                   n  = 2           ) : ##  polynomial degree 

        if x.getMin() != y.getMin() :
            logger.warning( 'PolyPos2Dsym: x&y have different low  edges %s vs %s'
                            % ( x.getMin() , y.getMin() ) )
        if x.getMax() != y.getMax() :
            logger.warning( 'PolyPos2Dsym: x&y have different high edges %s vs %s'
                            % ( x.getMax() , y.getMax() ) )
            
        PDF2.__init__ ( self , name , x , y ) 
        
        # 
        num = ( n + 1 ) * ( n + 2 ) / 2 
        self.makePhis ( num - 1 ) 
        #

        ## finally build PDF 
        #
        self.pdf = Ostap.Models.Poly2DSymPositive (
            'p2Dsp_%s'              % name ,
            'Poly2DSymPositive(%s)' % name ,
            self.x        ,
            self.y        ,
            n             ,
            self.phi_list )


models.append ( PolyPos2Dsym_pdf ) 
# =============================================================================
## @class PSPol2D_pdf
#  Product of phase space factors, modulated by positive polynomial in 2D 
#  \f$  f(x,y) = \Phi_1(x) \times \Phi_2(y) \times P^+(x,y) \f$,
#  where \f$ P^+(x,y)\f$ denotes the positive polynomial,
#  @see Ostap::Models::PS2DPol
#  @see Gaudi::Math::PS2DPol
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class PSPol2D_pdf(PDF2) :
    """Product of phase space factors, modulated by the positive polynom in 2D
    
    f(x,y) = PSX(x) * PSY(y) * Pnk(x,y)
    
    where
    - PSX(x) is a phase space function for x-axis (Ostap::Math::PhaseSpaceNL)
    - PSY(y) is a phase space function for y-axis (Ostap::Math::PhaseSpaceNL)
    - Pnk(x,y) is positive non-factorizable polynom
    
    Pnk(x,y) = sum^{i=n}_{i=0}sum{j=k}_{j=0} a^2_{ij} B^n_i(x) B^k_j(y)
    where:
    - B^n_i - are Bernstein polynomials
    
    Note:
    - f(x,y)>=0 for whole 2D-range

    """
    def __init__ ( self             ,
                   name             ,
                   x                ,   ##  the first  dimension  
                   y                ,   ##  the second dimension
                   psx              ,   ##  phase space in X, Ostap::Math::PhaseSpaceNL 
                   psy              ,   ##  phase space in Y, Ostap::Math::PhaseSpaceNL 
                   nx  = 2          ,   ##  polynomial degree in X 
                   ny  = 2          ) : ##  polynomial degree in Y 
        
        PDF2.__init__ ( self , name , x , y ) 

        self.psx = psx  
        self.psy = psy
                        
        #
        num = ( nx + 1 ) * ( ny + 1 )  - 1 
        self.makePhis ( num ) 
        #
            
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.PS2DPol (
            'ps2D_%s'     % name ,
            'PS2DPol(%s)' % name ,
            self.x        ,
            self.y        ,
            psx           ,
            psy           ,
            nx            ,
            ny            , 
            self.phi_list )
        
models.append ( PSPol2D_pdf ) 
# =============================================================================
## @class PSPol2Dsym_pdf
#  Symmetric product of phase space factors, modulated by positiev polynomial in 2D 
#  \f$  f(x,y) = \Phi(x) \times \Phi(y) \times P^+_{sym}(x,y) \f$,
#  where \f$ P^+_{sym}(x,y)\f$ denotes the symmetric positive polynomial,
#  @see Ostap::Models::PS2DPolSym
#  @see Ostap::Math::PS2DPolSym
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class PSPol2Dsym_pdf(PDF2) :
    """Symmetric product of phase space factors, modulated by the symmetrical
    positive polynom in 2D
    
    f(x,y) = PS(x) * PS(y) * Sn(x,y)
    
    where
    - PS(x) is a phase space function for axis (Ostap::Math::PhaseSpaceNL)
    - Sn(x,y) is positive non-factorizable symemtric polynom
    
    Sn(x,y)= sum^{i=n}_{i=0}sum{j=n}_{j=0} a^2_{ij} B^n_i(x) B^n_j(y)
    
    where:
    
    - B^n_i - are Bernstein polynomials
    - a_{ij} = a_{ji}
    
    Note:
    - f(x,y)>=0 for whole 2D-range
    - f(x,y) = f(y,x)
    """
    def __init__ ( self             ,
                   name             ,
                   x                ,   ##  the first  dimension  
                   y                ,   ##  the second dimension
                   ps               ,   ##  phase space in X, Ostap::Math::PhaseSpaceNL 
                   n   = 2          ) : ##  polynomial degree in Y 
        
        if x.getMin() != y.getMin() :
            logger.warning( 'PSPos2Dsym: x&y have different low  edges %s vs %s'
                            % ( x.getMin() , y.getMin() ) )
            
        if x.getMax() != y.getMax() :
            logger.warning( 'PSPos2Dsym: x&y have different high edges %s vs %s'
                            % ( x.getMax() , y.getMax() ) )
                
        PDF2.__init__ ( self , name , x , y ) 

        self.ps  = ps
        self.psx = ps
        self.psy = ps
            
        #
        num = ( n + 1 ) * ( n + 2 ) / 2  - 1 
        self.makePhis ( num ) 
        #
            
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.PS2DPolSym (
            'ps2Ds_%s'       % name ,
            'PS2DPolSym(%s)' % name ,
            self.x        ,
            self.y        ,
            ps            ,
            n             ,
            self.phi_list )
        
models.append ( PSPol2Dsym_pdf ) 
# =============================================================================
## @class ExpoPSPol2D_pdf
#  Product of phase space factors, modulated by positiev polynomial in 2D 
#  \f$  f(x,y) = exp(\tau x) \times \Phi (y) \times P^+(x,y) \f$,
#  where \f$ P^+(x,y)\f$ denotes the positive polynomial,
#  @see Ostap::Models::ExpoPS2DPol
#  @see Ostap::Math::ExpoPS2DPol
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class ExpoPSPol2D_pdf(PDF2) :
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
                   x                ,   ##  the first  dimension  
                   y                ,   ##  the second dimension
                   psy = None       ,   ##  phase space in Y, Ostap::Math::PhaseSpaceNL 
                   nx  = 2          ,   ##  polynomial degree in X 
                   ny  = 2          ,   ##  polynomial degree in Y 
                   tau = None       ) : ##  the exponent 
        
        PDF2.__init__ ( self , name , x , y ) 

        #
        ## get tau
        taumax  = 100
        mn,mx   = x.minmax()
        mc      = 0.5 * ( mn + mx )
        #
        if not iszero ( mn ) : taumax =                100.0 / abs ( mn ) 
        if not iszero ( mc ) : taumax = min ( taumax , 100.0 / abs ( mc ) )
        if not iszero ( mx ) : taumax = min ( taumax , 100.0 / abs ( mx ) )
        
        #
        ## the exponential slope
        #
        self.tau  = makeVar ( tau              ,
                              "tau_%s"  % name ,
                              "tau(%s)" % name , tau , 0 , -taumax , taumax )
        #
        
        self.psy = psy
                        
        #
        num = ( nx + 1 ) * ( ny + 1 ) - 1 
        self.makePhis ( num ) 
        #
            
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.ExpoPS2DPol (
            'ps2D_%s'     % name ,
            'PS2DPol(%s)' % name ,
            self.x        ,
            self.y        ,
            self.tau      ,
            self.psy      , 
            nx            ,
            ny            , 
            self.phi_list )
        
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
class ExpoPol2D_pdf(PDF2) :
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
                   x                ,   ##  the first  dimension  
                   y                ,   ##  the second dimension
                   nx  = 2          ,   ##  polynomial degree in X 
                   ny  = 2          ,   ##  polynomial degree in Y
                   taux = None      ,   ##  the exponent in X 
                   tauy = None      ) : ##  the exponent in Y
        
        PDF2.__init__ ( self , name , x , y ) 

        #
        ## get tau_x
        #
        xtaumax = 100
        mn,mx   = x.minmax()
        mc      = 0.5 * ( mn + mx )
        #
        if not iszero ( mn ) : xtaumax =                 100.0 / abs ( mn ) 
        if not iszero ( mc ) : xtaumax = min ( xtaumax , 100.0 / abs ( mc ) )
        if not iszero ( mx ) : xtaumax = min ( xtaumax , 100.0 / abs ( mx ) )
        
        ytaumax = 100
        mn,mx   = y.minmax()
        mc      = 0.5 * ( mn + mx )
        #
        if not iszero ( mn ) : ytaumax =                 100.0 / abs ( mn ) 
        if not iszero ( mc ) : ytaumax = min ( ytaumax , 100.0 / abs ( mc ) )
        if not iszero ( mx ) : ytaumax = min ( ytaumax , 100.0 / abs ( mx ) )

        #
        ## the exponential slopes
        #
        self.taux  = makeVar ( taux              ,
                               "taux_%s"  % name ,
                               "taux(%s)" % name , taux , 0 , -xtaumax , xtaumax )
        #
        self.tauy  = makeVar ( tauy              ,
                               "tauy_%s"  % name ,
                               "tauy(%s)" % name , tauy , 0 , -ytaumax , ytaumax )
        #
        
        self.m1  = x ## ditto 
        self.m2  = y ## ditto

        #
        num = ( nx + 1 ) * ( ny + 1 ) - 1 
        self.makePhis ( num ) 
        #
            
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.Expo2DPol (
            'exp2D_%s'      % name ,
            'Expo2DPol(%s)' % name ,
            self.x        ,
            self.y        ,
            self.taux     ,
            self.tauy     ,
            nx            ,
            ny            , 
            self.phi_list )
        
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
class ExpoPol2Dsym_pdf(PDF2) :
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
                   x                ,   ## the first  dimension  
                   y                ,   ## the second dimension
                   n    = 2         ,   ## polynomial degree in X and Y
                   tau  = None      ) : ## the exponent 
        
        if x.getMin() != y.getMin() :
            logger.warning( 'PSPos2Dsym: x&y have different low  edges %s vs %s'
                            % ( x.getMin() , y.getMin() ) )
            
        if x.getMax() != y.getMax() :
            logger.warning( 'PSPos2Dsym: x&y have different high edges %s vs %s'
                            % ( x.getMax() , y.getMax() ) )
                
        PDF2.__init__ ( self , name , x , y ) 

        #
        ## get tau
        # 
        taumax  = 100
        mn,mx   = x.minmax()
        mc      = 0.5 * ( mn + mx )
        #
        if not iszero ( mn ) : taumax =                100.0 / abs ( mn ) 
        if not iszero ( mc ) : taumax = min ( taumax , 100.0 / abs ( mc ) )
        if not iszero ( mx ) : taumax = min ( taumax , 100.0 / abs ( mx ) )

        #
        ## the exponential slopes
        #
        self.tau  = makeVar ( tau              ,
                              "tau_%s"  % name ,
                              "tau(%s)" % name , tau , 0 , -taumax , taumax )
        #
        self.m1  = x ## ditto 
        self.m2  = y ## ditto

        #
        num = ( n + 1 ) * ( n + 2 ) / 2 - 1 
        self.makePhis ( num ) 
        #
            
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.Expo2DPolSym (
            'exp2Ds_%s'        % name ,
            'Expo2DPolSym(%s)' % name ,
            self.x        ,
            self.y        ,
            self.tau      ,
            n             ,
            self.phi_list )
        

models.append ( ExpoPol2Dsym_pdf ) 
# =============================================================================
## @class Spline2D_pdf
#  positive spline in 2D:
#  \f$  f(x,y) = \sum^{i=n}_{i=0}\sum{j=k}_{j=0} a^2_{\ij} M^n_i(x) M^k_j(y) \f$,
#  where \f$ B^n_i(x)\f$ denotes the M-splines  
#  @see Ostap::Models::Spline2D
#  @see Ostap::Math::Spline2D
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class Spline2D_pdf(PDF2) :
    """Positive non-factorizable spline in 2D
    
    f(x,y) = sum_i sum_j a^2_{i,j} Nx_i(x) * Ny_j(y),
    where
    - Nx_i and Ny_j are normailzed B-splines for x and y-axes

    Note:
    - f(x,y)>=0 for whole 2D-range
    """
    def __init__ ( self             ,
                   name             ,
                   x                ,   ##  the first  dimension  
                   y                ,   ##  the second dimension
                   spline           ) : ## the spline: Gaudi.Math.Spline2D 

        PDF2.__init__ ( self , name , x , y )
        
        self.spline = spline

        #
        self.makePhis ( spline.npars() ) 
        #
            
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.Spline2D (
            's2Dp_%s'      % name ,
            'Spline2D(%s)' % name ,
            self.x        ,
            self.y        ,
            self.spline   ,
            self.phi_list )

models.append ( Spline2D_pdf ) 
# =============================================================================
## @class Spline2Dsym_pdf
#  symmetric positive spline in 2D:
#  \f$  f(x,y) = \sum^{i=n}_{i=0}\sum{j=k}_{j=0} a^2_{\ij} M^n_i(x) M^k_j(y) \f$,
#  where \f$ B^n_i(x)\f$ denotes the M-splines  
#  @see Ostap::Models::Spline2DSym
#  @see Ostap::Math::Spline2DSym
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-01-10
class Spline2Dsym_pdf(PDF2) :
    """SYMMETRIC positive non-factorizable spline in 2D
    
    f(x,y) = sum_i sum_j a^2_{i,j} N_i(x) * N_j(y),
    where
    - N_i are normailzed B-splines 

    Note:
    - f(x,y)>=0     for whole 2D-range
    - f(x,y)=f(y,x) for whole 2D-range 
    """
    def __init__ ( self             ,
                   name             ,
                   x                ,   ##  the first  dimension  
                   y                ,   ##  the second dimension
                   spline           ) : ## the spline: Gaudi.Math.Spline2DSym 

        PDF2.__init__ ( self , name , x , y )
        
        self.spline = spline

        #
        self.makePhis ( spline.npars() ) 
        #

        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.Spline2DSym (
            's2Dp_%s'         % name ,
            'Spline2DSym(%s)' % name ,
            self.x        ,
            self.y        ,
            self.spline   ,
            self.phi_list )

models.append ( Spline2Dsym_pdf ) 
# =============================================================================
# some tiny decoration of underlying classes 
# =============================================================================
_rv = ROOT.gROOT.GetVersionInt() // 10000
def _2d_get_pars_ ( self ) :
    """Get parameters of underlying positive Berstein polynomial
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
                
                if _rv < 6 : m[i][j] = b.par(i,j)
                else       : m[i, j] = b.par(i,j)
                    
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
if '__main__' == __name__ : 
         
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )
    
# =============================================================================
# The END 
# =============================================================================
