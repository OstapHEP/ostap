#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/spectra.py
#  A set of predefined ready-to-use
#  functions and shapes for fitting of pt-spectra 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
# =============================================================================
"""A set of predefined ready-to-use shapes and PDFs
for fitting pT-spectra"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2015-07-11"
__all__     = (
    #
    'Tsallis_pdf'   , ## Tsallis   PDF 
    'QGSM_pdf'      , ## QGSM      PDG 
    'Hagedorn_pdf'  , ## Hagedorn's PDF 
    'Tsallis2_pdf'  , ## Tsallis   PDF 
    #
    )
# =============================================================================
from   ostap.core.core             import cpp, Ostap, VE , funID
from   ostap.fitting.distributions import GammaDist_pdf
from   ostap.fitting.pdfbasic      import PDF1, PDF2
import ROOT, math
# =============================================================================
from   ostap.logger.logger     import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.spectra' )
else                       : logger = getLogger ( __name__           )
# =============================================================================
models = []
# =============================================================================
## @class Tsallis_pdf
#  Useful function to describe pT-spectra of particles 
#
#  - C. Tsallis, 
#  "Possible generalization of Boltzmann-Gibbs statistics,
#  J. Statist. Phys. 52 (1988) 479.
#  - C. Tsallis, 
#  Nonextensive statistics: theoretical, experimental and computational 
#  evidences and connections, Braz. J. Phys. 29 (1999) 1.
# 
#  \f[ \frac{d\sigma}{dp_T} \propto  
#    p_T\times \left( 1 + \frac{E_{kin}}{Tn}\right)^{-n}\f],
#  where \f$E_{kin} = \sqrt{p_T^2-M^2}-M\f$ 
#  is transverse kinetic energy 
#
#  @see Ostap::Models::Tsallis
#  @see Ostap::Math::Tsallis
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Tsallis_pdf(PDF1) :
    r""" Useful function to describe pT-spectra of particles 
    
    - C. Tsallis, 
    Possible generalization of Boltzmann-Gibbs statistics,
    J. Statist. Phys. 52 (1988) 479.
    - C. Tsallis, 
    Nonextensive statistics: theoretical, experimental and computational 
    evidences and connections, Braz. J. Phys. 29 (1999) 1.
    
    \f[ \frac{d\sigma}{dp_T} \propto  
    p_T\times \left( 1 + \frac{E_{kin}}{Tn}\right)^{-n}\f],
    
    where \f$E_{kin} = \sqrt{p_T^2-M^2}-M\f$
    
    is transverse kinetic energy 
    """
    def __init__ ( self                   , * , 
                   xvar                   ,   ## pT-variable (for fitting) 
                   name      = ''         , 
                   m0        = 0.135      ,   ## particle mass (may be fixed)
                   n         = None       ,   ## shape parameter
                   T         = None       ) : ## temperature parameter                   

        ## initialize the base 
        PDF1.__init__  ( self , name = name , xvar = xvar )
        
        self.__m0   = self.make_var ( m0                   ,
                                      'm0_%s'  % self.name , 
                                      'm0(%s)' % self.name ,
                                      True , 1.e-10   , 1e+6 )
        
        self.__n    = self.make_var ( n               ,
                                      'n_%s'   % self.name , 
                                      'n(%s) ' % self.name ,
                                      False , 1 , 0.01  , 1000 )  
        
        self.__T    = self.make_var ( T               ,
                                      'T_%s'   % self.name , 
                                      'T(%s) ' % self.name ,
                                      False , 1 , 1.e-4 , 1e+6 )
        
        self.pdf  = Ostap.Models.Tsallis (
            self.roo_name ( 'tsallis_' )   ,
            'Tsallis %s' % self.name  , 
            self.pt               ,
            self.n                ,
            self.T                ,
            self.m0               )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'n'     : self.n     ,            
            'T'     : self.T     ,            
            'm0'    : self.m0    ,            
            }
    
    @property
    def pt ( self ) :
        """`pt'-variable for Tsallis distribution (the same as ``x'')"""
        return self.xvar
    
    @property
    def m0 ( self ) :
        """`m0'-parameter of Tsallis' function"""
        return self.__m0
    @m0.setter
    def m0 ( self , value ) :
        self.set_value ( self.__m0 , value )

    @property
    def n ( self ) :
        """`n'-parameter of Tsallis' function"""
        return self.__n
    @n.setter
    def n ( self , value ) :
        self.set_value ( self.__n , value )

    @property
    def T ( self ) :
        """`T'-parameter of Tsallis' function"""
        return self.__T
    @T.setter
    def T ( self , value ) :
        self.set_value ( self.__T , value )
 
# =============================================================================
## @class QGSM_pdf
#  Useful function to describe pT-spectra of particles 
#
# - A. B. Kaidalov and O. I. Piskunova, Z. Phys. C 30 (1986) 145.
# - O. I. Piskounova, arXiv:1301.6539 [hep-ph]; 
# - O. I. Piskounova, arXiv:1405.4398 [hep-ph].
# - A. A. Bylinkin and O. I. Piskounova, 
#  "Transverse momentum distributions of baryons at LHC energies",
#  arXiv:1501.07706.
#
#  \f[ \frac{d\sigma}{dp_T} \propto 
#  p_T \times \mathrm{e}^{ -b_0 (m_T-m)} \f], 
#  where transverse mass is defined as \f$m_T = \sqrt{p_T^2+m^2}\f$
# 
#  @see Ostap::Models::QGSM
#  @see Ostap::Math::QGSM
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class QGSM_pdf(PDF1) :
    r""" Useful function to describe pT-spectra of particles 
    
    - A. B. Kaidalov and O. I. Piskunova, Z. Phys. C 30 (1986) 145.
    - O. I. Piskounova, arXiv:1301.6539 [hep-ph]; 
    - O. I. Piskounova, arXiv:1405.4398 [hep-ph].
    - A. A. Bylinkin and O. I. Piskounova, 
    'Transverse momentum distributions of baryons at LHC energies',
    arXiv:1501.07706.

    \f[ \frac{d\sigma}{dp_T} \propto p_T \times \mathrm{e}^{ -b_0 (m_T-m)} \f],
    
    where transverse mass is defined as \f$m_T = \sqrt{p_T^2+m^2}\f$
    """
    def __init__ ( self             ,
                   name             , 
                   xvar             ,   ## pT-variable (for fitting) 
                   m0        = 0    ,   ## particle mass (may be fixed)
                   b         = None ) : ## slope parameter
        
        ## initialize the base 
        PDF1.__init__  ( self , name = name , xvar = xvar )

        
        self.__m0   = self.make_var ( m0              ,
                                      'm0_%s'  % self.name , 
                                      'm0(%s)' % self.name ,
                                      True , 0  , 1e+6 )
        
        self.__b    = self.make_var ( b               ,
                                      'b_%s'   % self.name , 
                                      'b(%s) ' % self.name ,
                                      False , 0. , 1e+6 )  
        
        self.pdf  = Ostap.Models.QGSM (
            self.roo_name ( 'qgsm_' ) ,
            'QGSM %s' % self.name , 
            self.pt               ,
            self.b                ,
            self.m0               )

        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'b'     : self.b     ,            
            'm0'    : self.m0    ,            
            }
    
    @property
    def pt ( self ) :
        """`pt'-variable for QGSM distribution (the same as ``x'')"""
        return self.xvar
    
    @property
    def m0 ( self ) :
        """`m0'-parameter of QGSM function"""
        return self.__m0
    @m0.setter
    def m0 ( self , value ) :
        self.set_value ( self.__m0 , value )

    @property
    def b ( self ) :
        """`b'-parameter of QGSM function"""
        return self.__b
    @b.setter
    def b ( self , value ) :
        self.set_value ( self.__b , value )

# =============================================================================
## @class Hagedorn_pdf
#  Useful function to describe pT-spectra of particles 
#  @see R.Hagedorn, "Multiplicities, p_T distributions and the 
#       expected hadron \to Quark - Gluon Phase Transition", 
#       Riv.Nuovo Cim. 6N10 (1983) 1-50
#  @see https://doi.org/10.1007/BF02740917 
#  @see https://inspirehep.net/literature/193590
#  
#  \f[ f(p_T; m, T) \propto 
#   p_T \sqrt{p^2_T + m^2} K_1( \beta \sqrt{ p^2_T+m^2} ) \f] 
#
#  where \f$ \beta \f$ is inverse temporature 
#  \f$ \beta = \frac{1}{T} f$ 
#
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2022-12-06
#  @see Ostap::Models::Hagedorn
#  @see Ostap::Math::Hagedorn
class Hagedorn_pdf(PDF1) :
    r"""Useful function to describe pT-spectra of particles 
    
    - see R.Hagedorn, 'Multiplicities, p_T distributions and the 
    expected hadron \to Quark - Gluon Phase Transition', 
    Riv.Nuovo Cim. 6N10 (1983) 1-50
    - see https://doi.org/10.1007/BF02740917 
    - see https://inspirehep.net/literature/193590
    
    \f[ f(p_T; m, T) \propto 
    p_T \sqrt{p^2_T + m^2} K_1( \beta \sqrt{ p^2_T+m^2} ) \f] 
    
    where \f$ \beta \f$ is inverse temporature 
    \f$ \beta = \frac{1}{T} f$ 
    """
    def __init__ ( self             ,
                   name             , 
                   xvar             ,   ## pT-variable (for fitting) 
                   m0        = 0    ,   ## particle mass (may be fixed)
                   beta      = None ) : ## inverse temperature
        
        ## initialize the base 
        PDF1.__init__  ( self , name  = name  , xvar = xvar )

        
        self.__m0   = self.make_var ( m0              ,
                                      'm0_%s'  % self.name , 
                                      'm0(%s)' % self.name ,
                                      True , 0  , 1e+6 )
        
        self.__beta = self.make_var ( beta                 ,
                                      'beta_%s'    % self.name  , 
                                      '#beta(%s) ' % self.name  ,
                                      False , 1.e-6 , 1e+6 )  
        
        self.pdf  = Ostap.Models.Hagedorn (
            self.roo_name ( 'hage_' ) ,
            'Hagedorn %s' % self.name , 
            self.pt               ,
            self.beta             ,
            self.m0               )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'beta'  : self.beta  ,             
            'm0'    : self.m0    ,            
            }
    
    @property
    def pt ( self ) :
        """'pt'-variable for Hagedorn distribution (the same as ``x'')"""
        return self.xvar
    
    @property
    def m0 ( self ) :
        """'m0'-parameter of Hagedorn function"""
        return self.__m0
    @m0.setter
    def m0 ( self , value ) :
        self.set_value ( self.__m0 , value )

    @property
    def beta ( self ) :
        """'beta'-parameter (inverse temperature) of Hagedorn  function"""
        return self.__beta
    @beta.setter
    def beta ( self , value ) :
        self.set_value ( self.__beta , value )
       
# =============================================================================
## @class Tsallis2_pdf
#  Useful function to describe pT and saidity -spectra of particles 
#
#  2D particle density distribution as function of pt and rapidity 
#  @see L. Marques, J. Cleymans, A. Deppman, 
#       "Description of High-Energy pp Collisions 
#        Using Tsallis Thermodynamics: 
#        Transverse Momentum and Rapidity Distributions", 
#        Phys. Rev. D 91, 054025, 	arXiv:1501.00953 
#  @see https://arxiv.org/abs/1501.00953
#  @see https://doi.org/10.1103/PhysRevD.91.054025
#  @see Ostap::Models::Tsallis2
#  @see Ostap::Models::Tsallis
#  @see Ostap::Math::Tsallis2
#  @see Ostap::Math::Tsallis
class Tsallis2_pdf(PDF2) :
    """ 2D particle density distribution as function of pt and rapidity 
    @see L. Marques, J. Cleymans, A. Deppman, 
    ``Description of High-Energy pp Collisions 
    Using Tsallis Thermodynamics: 
    Transverse Momentum and Rapidity Distributions'', 
    Phys. Rev. D 91, 054025, 	arXiv:1501.00953 
    - see https://arxiv.org/abs/1501.00953
    - see https://doi.org/10.1103/PhysRevD.91.054025
    - see `Ostap.Models.Tsallis2`
    - see `Ostap.Models.Tsallis`
    - see `Ostap.Math.Tsallis2`
    - see `Ostap.Math.Tsallis`
    """
    def __init__ ( self                   ,
                   xvar                   ,   ## pT-observable )
                   yvar                   ,   ## rapidity observable
                   m0        = 1          ,   ## partile mass (presumably constant) 
                   T         = None       ,   ## temperature parameter
                   q         = 1.1        ,   ## q-parameter
                   mu        = 0          ,   ## chemical potential                    
                   name      = ''         ) :
        
        ## initialize the base 
        PDF2.__init__  ( self , name = name , xvar = xvar , yvar =  yvar )
        
        self.__m0   = self.make_var ( m0              ,
                                      'm0_%s'  % self.name , 
                                      'm0(%s)' % self.name ,
                                      True , 0     , 1e+6 )
        
        self.__q    = self.make_var ( q               ,
                                      'q_%s'   % self.name , 
                                      'q(%s) ' % self.name ,
                                      None , 1.1 , 1.e-6, 10 )  
        
        self.__T    = self.make_var ( T               ,
                                      'T_%s'   % self.name , 
                                      'T(%s) ' % self.name ,
                                      None , 1 , 1.e-4 , 1e+6 )
        
        self.__mu   = self.make_var ( mu               ,
                                      'mu_%s'   % self.name , 
                                      '#mu(%s) '% self.name ,
                                      True  , 0 , 0 , 100 )
        
        self.pdf  = Ostap.Models.Tsallis2 (
            self.roo_name ( 'tsallis2_' )   ,
            'Tsallis2 %s' % self.name  , 
            self.pt               ,
            self.y                ,
            self.m0               ,
            self.T                ,
            self.q                ,
            self.mu               )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'yvar'  : self.yvar  ,
            'm0'    : self.m0    ,            
            'q'     : self.q     ,            
            'T'     : self.T     ,            
            'mu'    : self.mu    ,            
            }
        
    @property
    def pt ( self ) :
        """'pt'-observable (transverse momentum) for Tsallis distribution (the same as 'xvar')"""
        return self.xvar
    
    @property
    def y ( self ) :
        """'y'-observable (rapidity) for Tsallis distribution (the same as 'yvar')"""
        return self.yvar

    @property
    def m0 ( self ) :
        """'m0'-parameter (particle mass) of Tsallis' function"""
        return self.__m0
    @m0.setter
    def m0 ( self , value ) :
        self.set_value ( self.__m0 , value )

    @property
    def q ( self ) :
        """'q'-parameter (shape) of Tsallis' function"""
        return self.__q
    @q.setter
    def q ( self , value ) :
        self.set_value ( self.__q , value )

    @property
    def T ( self ) :
        """'T'-parameter (temperature) of Tsallis' function"""
        return self.__T
    @T.setter
    def T ( self , value ) :
        self.set_value ( self.__T , value )

    @property
    def mu ( self ) :
        """'mu'-parameter (chemical potential) of Tsallis' function"""
        return self.__mu
    @mu.setter
    def mu ( self , value ) :
        self.set_value ( self.__mu , value )
     
models  .append ( Tsallis_pdf   )        
models  .append ( QGSM_pdf      )
models  .append ( Hagedorn_pdf  )
models  .append ( Tsallis2_pdf  )
models  .append ( GammaDist_pdf )
# ==============================================================================
## @class PtFitBase
#  helper object  
class PtFitBase(object) :
    """ Helper class to implement pt-spectrum fitter
    """
    def __init__ ( self , ptmax , ptmin = 0 ) :
        self._integral = -1
        ptmax = float ( ptmax ) 
        ptmin = float ( ptmin ) 
        assert ptmin < ptmax , "PtFitBase: wrong ptmin/ptmax: %s/%s" % ( ptmin ,ptmax )
        self._ptmin    = ptmin
        self._ptmax    = ptmax

    ## get the mean-value 
    def mean     ( self ) : return self._fun.mean     ( self._ptmin , self._ptmax )
    ## get variance
    def variance ( self ) : return self._fun.variance ( self._ptmin , self._ptmax )
    ## get rms 
    def rms      ( self ) : return self._fun.rms      ( self._ptmin , self._ptmax )

    ## make final evaluation of the result 
    def result   ( self , pt , norm , changed  ) :
        ## recalculate integral (if needed) 
        if changed or 0 >= self._integral :
            self._integral = self._fun.integral ( self._ptmin , self._ptmax ) 
            
        return norm * self._fun ( pt ) / self._integral 
    
# ==============================================================================
## @class TsallisFun
#  Tsallis'  function for fitting pt-spectra 
class TsallisFun(PtFitBase) :
    """Tsallis'  function for fitting pt-spectra 
    """
    def __init__ ( self , ptmax , ptmin = 0 ) :        
        ## initialize the base 
        PtFitBase.__init__ ( self , ptmax , ptmin ) 
        ## create Tsallis function
        self._fun = Ostap.Math.Tsallis(0,10,1.1)
        
    def __call__ ( self , x , pars ) :
        
        pt   = x    [ 0 ]
        
        norm = pars [ 0 ]  ## normalization (almost arbitrary) 
        m0   = pars [ 1 ]  ## 
        n    = pars [ 2 ]
        T    = pars [ 3 ]
        
        changed = False 
        if self._fun.setMass ( m0 ) : changed = True 
        if self._fun.setN    ( n  ) : changed = True 
        if self._fun.setT    ( T  ) : changed = True 

        return self.result ( pt , norm , changed ) 


# ==============================================================================
## @class QGSMFun
#  QGSM function for fitting pt-spectra 
class QGSMFun(PtFitBase) :
    """QGSM function for fitting pt-spectra
    """
    def __init__ ( self , ptmax , ptmin = 0 ) :
        ## initialize the base 
        PtFitBase.__init__ ( self , ptmax , ptmin ) 
        ## create QGSM function 
        self._fun = Ostap.Math.QGSM(0,1)
        
    ## get the mean-value 
    def mean     ( self ) :
        return self._fun.mean ( self._ptmin , self._ptmax )
    
    def __call__ ( self , x , pars ) :
        
        pt = x    [ 0 ]
        
        norm = pars [ 0 ]  ## normalization (almost arbitrary) 
        m0   = pars [ 1 ]  ## 
        b    = pars [ 2 ]
        
        changed = False 
        if self._fun.setMass ( m0 ) : changed = True 
        if self._fun.setB    ( b  ) : changed = True 

        return self.result ( pt , norm , changed ) 

# ==============================================================================
## @class GammaDistFun
#  Gamma  distribution as  fit-fuction
class GammaDistFun(PtFitBase) :
    """Gamma  distribution as  fit-fuction
    """
    def __init__ ( self , ptmax , ptmin = 0 ) :
        ## initialize the base 
        PtFitBase.__init__ ( self , ptmax , ptmin ) 
        ## create QGSM function 
        self._fun = Ostap.Math.GammaDist(2,1)
        
    def __call__ ( self , x , pars ) :
        
        pt    = x    [ 0 ]
        
        norm  = pars [ 0 ]  ## normalization (almost arbitrary) 
        k     = pars [ 1 ]  ## 
        theta = pars [ 2 ]
        
        changed = False 
        if self._fun.setK     ( k     ) : changed = True 
        if self._fun.setTheta ( theta ) : changed = True 
        
        return self.result ( pt , norm , changed ) 

# =============================================================================
## create Tsallis function (as TF1) for fitting pT-spectra
#  @code
#  >>> histo =
#  >>> fun   = tsallisTF1 ( ptmax = 50 )
#  >>> histo.Fit( fun , 'SI' )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def tsallisTF1( ptmax , ptmin = 0 , mass = None , name = '' ) :
    """Create Tsallis function (as TF1) for fitting pt-spectra:
    
    >>> histo =
    >>> fun   = tsallisFun ( ptmax = 50 )
    >>> histo.Fit( fun , 'SI' )    
    """
    if ptmax <= ptmin :
        raise ArrtibuteError("Tsallis/TF1: wrong ptmin/ptmax: %s/%s" % ( ptmin ,ptmax ) )
    
    obj       = TsallisFun ( ptmax , ptmin )
    if not name : name =  funID() 
    func      = ROOT.TF1( name , obj , ptmin , ptmax , 4 )
    func._obj = obj
    func.SetParNames (
        'norm' , ## 0
        'm'    , ## 1 
        'n'    , ## 2
        'T'    ) ## 3
    ##
    if isinstance ( mass , float ) and 0 <= mass : func.FixParameter ( 1 , mass )
    else                                         : func.SetParameter ( 1 , 0    )
    ## 
    func.SetParameter ( 0 , 1   ) ## norm 
    func.SetParameter ( 2 , 10  )
    func.SetParameter ( 3 , 1.1 )
    ##
    return func

# =============================================================================
## create QGSM function (as TF1) for fitting pT-spectra
#  @code
#  >>> histo =
#  >>> fun   = qgsmTF1 ( ptmax = 50 )
#  >>> histo.Fit( fun , 'SI' )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def qgsmTF1 ( ptmax , ptmin = 0 , mass = None , name = '') :
    """Create QGSM functin (as TF1) for fitting pt-spectra
    
    >>> histo =
    >>> fun   = qgsmTF1 ( ptmax = 50 )
    >>> histo.Fit( fun , 'SI' )
    """
    if ptmax <= ptmin :
        raise ArrtibuteError("QGSM/TF1: wrong ptmin/ptmax: %s/%s" % ( ptmin ,ptmax ) )
    
    obj       = QGSMFun  (                 ptmax , ptmin     )
    if not name : name =  funID() 
    func      = ROOT.TF1 ( name , obj , ptmin , ptmax , 3 )
    func._obj = obj
    func.SetParNames (
        'norm' , ## 0
        'm'    , ## 1 
        'b'    ) ## 2
    ##
    if isinstance ( mass , float ) and 0 <= mass : func.FixParameter ( 1 , mass )
    else                                         : func.SetParameter ( 1 , 0    )
    ## 
    func.SetParameter ( 0 , 1 )  ## norm 
    func.SetParameter ( 2 , 1 )
    ##
    return func

# =============================================================================
## create GammaDist function (as TF1) for fitting pT-spectra
#  @code
#  >>> histo =
#  >>> fun   = gammaDistTF1 ( ptmax = 50 )
#  >>> histo.Fit( fun , 'SI' )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def gammaDistTF1 ( ptmax , ptmin = 0 , mass = None , name = '' ) :
    """ Create GammaDist function (as TF1) for fitting pt-spectra
    >>> histo =
    >>> fun   = gammaDistTF1 ( ptmax = 50 )
    >>> histo.Fit( fun , 'SI' )
    """
    if ptmax <= ptmin :
        raise ArrtibuteError("GammaDist/TF1: wrong ptmin/ptmax: %s/%s" % ( ptmin ,ptmax ) )
    
    obj       = GammaDistFun (                 ptmax , ptmin     )
    if not name : name =  funID() 
    func      = ROOT.TF1     ( name , obj , ptmin , ptmax , 3 )
    func._obj = obj
    func.SetParNames (
        'norm'  , ## 0
        'k'     , ## 1 
        'theta' ) ## 2
    ##
    func.SetParameter ( 0 , 1 )  ## norm 
    func.SetParameter ( 1 , 2 )
    func.SetParameter ( 2 , 2 )
    ##
    return func


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )
 
# =============================================================================
##                                                                     The END 
# =============================================================================
