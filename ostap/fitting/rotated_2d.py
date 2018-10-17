#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/rotated_2d.py
#  Collection of ``rotated'' 2D-models
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2018-10-12
# =============================================================================
"""Collection of ``rotated'' 2D-models
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2018-10-12"
__all__     = (
    'Rotated2Gauss_pdf' , ## Rotated product of two Gaussians 
    'Rotated2CB_pdf'    , ## Rotated product of two double-sided Crystal Ball functions 
    )
# =============================================================================
import ROOT, math
from   ostap.core.core       import Ostap
from   ostap.fitting.fit2d   import PDF2
from   ostap.fitting.signals import Gauss_pdf, CB2_pdf 
# =============================================================================
from   ostap.logger.logger     import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.rotated_2d' )
else                       : logger = getLogger ( __name__                   )
# =============================================================================
models = [] 
# =============================================================================
##  @class Rotated2D
#   helper base class for rotated models
#   - keeps the rotation angle <code>phi</code>
class Rotated2D (PDF2) :
    """Helper base class for rotated 2D-models
    - it keeps the rotation angle ``phi''
    """
    def __init__ ( self , name , xvar  , yvar , phi = None ) :

        ## initialize the base 
        PDF2  .__init__ ( self , name  , xvar , yvar  )

        ## create the rotation phase 
        self.__phi = self.make_var ( phi               ,
                                     "phi_%s"   % name ,
                                     "phi(%s)"  % name , phi , 0 , -3.5 , +6.5 )
    @property
    def phi ( self ) :
        "``phi''-rotation phase"
        return self.__phi
    @phi.setter
    def phi ( self, value ) :
        self.__phi.setVal ( float ( value ) ) 
        
# ==============================================================================
## @class Rotated2Gauss_pdf
#  "Rotated" profduct of two gaussian functions
#  @see Ostap::Models::Rotated2Gaussians
#  @see Ostap::Models::RotatedProduct
#  @see Ostap::Math::RotatedProduct
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2018-10-12
class Rotated2Gauss_pdf (Rotated2D) :
    """ ``Rotated'' product of two gaussian functions
    """
    def __init__ ( self           ,
                   name           ,
                   xvar           ,   ## x-variable 
                   yvar           ,   ## y-variable 
                   phi    = None  ,   ## rotation phase 
                   gauss1 = None  ,   ## 1st gausssian
                   mean1  = None  ,   ## the mean of the first gaussian
                   sigma1 = None  ,   ## the sigma of the first gaussian 
                   gauss2 = None  ,   ## 2nd gaussian 
                   mean2  = None  ,   ## the mean of the second gaussian
                   sigma2 = None  ):  ## the sigma of the second gaussian 
        
        ## inialize the base 
        Rotated2D  .__init__ ( self , name , xvar , yvar , phi )
        
        if   gauss1 and isinstance ( gauss1 , Gauss_pdf ) :
            mean1  = gauss1.mean
            sigma1 = gauss1.sigma
        elif gauss1 and isinstance ( gauss1 , Ostap.Math.Gauss ) :
            mean1  = gauss1.mean  ()
            sigma1 = gauss1.sigma ()
        elif gauss1 and isinstance ( gauss1 , ( tuple , list ) ) :
            mean1  = gauss1[0]
            sigma1 = gauss1[1]

        if   gauss2 and isinstance ( gauss2 , Gauss_pdf ) :
            mean2  = gauss2.mean
            sigma2 = gauss2.sigma
        elif gauss2 and isinstance ( gauss2 , Ostap.Math.Gauss ) :
            mean2  = gauss2.mean  ()
            sigma2 = gauss2.sigma ()
        elif gauss2 and isinstance ( gauss2 , ( tuple , list ) ) :
            mean2  = gauss2[0]
            sigma2 = gauss2[1]
            
        self.__gauss1  = gauss1
        self.__gauss2  = gauss2
            
        self.__gauss_1 = Gauss_pdf ( name + '_X' , xvar , mean1 , sigma1 )
        self.__gauss_2 = Gauss_pdf ( name + '_Y' , yvar , mean2 , sigma2 )

        ## finally create PDF 
        self.pdf = Ostap.Models.Rotated2Gaussians (
            'r2g_%s'            % name ,
            'Rotated2Gauss(%s)' % name ,
            self.x      ,
            self.y      ,
            self.phi    , 
            self.mean1  ,
            self.sigma1 ,
            self.mean2  ,
            self.sigma2 )
            
        ## save configuration
        self.config = {
            'name'   : self.name     ,
            'xvar'   : self.xvar     ,
            'yvar'   : self.yvar     ,
            'phi'    : self.phi      ,
            'gauss1' : self.__gauss1 ,
            'mean1'  : self.mean1    ,
            'sigma1' : self.sigma1   ,
            'gauss2' : self.__gauss2 ,
            'mean2'  : self.mean2    ,
            'sigma2' : self.sigma2   }
        
    @property
    def phi ( self ) :
        "``phi''-rotation phase"
        return self.__phi
    @phi.setter
    def phi ( self, value ) :
        self.__phi.setVal ( float ( value ) ) 

    @property
    def mean1 ( self ) :
        "``mean1''-the mean-value of the first Gaussian"
        return self.__gauss_1.mean 
    @mean1.setter
    def mean1 ( self, value ) :
        self.__gauss_1.mean = value 

    @property
    def sigma1 ( self ) :
        "``sigma1''-the sigma-value of the first Gaussian"
        return self.__gauss_1.sigma
    @sigma1.setter
    def sigma1 ( self, value ) :
        self.__gauss_1.sigma = value 

    @property
    def mean2 ( self ) :
        "``mean2''-the mean-value of the first Gaussian"
        return self.__gauss_2.mean 
    @mean2.setter
    def mean2 ( self, value ) :
        self.__gauss_2.mean = value 

    @property
    def sigma2 ( self ) :
        "``sigma2''-the sigma-value of the second Gaussian"
        return self.__gauss_2.sigma
    @sigma2.setter
    def sigma2 ( self, value ) :
        self.__gauss_2.sigma = value 


# ==============================================================================
## @class Rotated2CB_pdf
#  "Rotated" profduct of two double-sided Crystal Ball functions
#  @see Ostap::Models::Rotated2CrystalBAlls
#  @see Ostap::Models::RotatedProduct
#  @see Ostap::Math::RotatedProduct
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2018-10-12
class Rotated2CB_pdf (PDF2) :
    """ ``Rotated'' product of two double-sided Crystal Ball functions
    """
    def __init__ ( self            ,
                   name            ,
                   xvar            ,   ## x-variable 
                   yvar            ,   ## y-variable 
                   phi     = None  ,   ## rotation phase
                   ## 
                   cb1     = None  ,   ## 1st CB2 
                   mean1   = None  ,   ## the mean  of the first CB2
                   sigma1  = None  ,   ## the sigma of the first CB2
                   alphaL1 = None  ,   ## the alphaL parameter for the first CB2 
                   alphaR1 = None  ,   ## the alphaR parameter for the first CB2
                   nL1     = None  ,   ## the nL parameter     for the first CB2
                   nR1     = None  ,   ## the nR parameter     the the first CB2
                   ## 
                   cb2     = None  ,   ## 2nd CB2 
                   mean2   = None  ,   ## the mean  of the second CB2
                   sigma2  = None  ,   ## the sigma of the second  CB2
                   alphaL2 = None  ,   ## the alphaL parameter for the second CB2 
                   alphaR2 = None  ,   ## the alphaR parameter for the second CB2
                   nL2     = None  ,   ## the nL parameter     for the second CB2
                   nR2     = None  ) : ## the nR parameter     the the second CB2
        
        ## inialize the base 
        PDF2  .__init__ ( self , name  , xvar , yvar  )
        
        self.__phi  = self.make_var ( phi               ,
                                      "phi_%s"   % name ,
                                      "phi(%s)"  % name , phi , 0 , -3.0 , +6.0 )
        
        if cb1 and isinstance ( cb1 , CB2_pdf ) :
            mean1   = cb1.mean
            sigma1  = cb1.sigma
            alphaL1 = cb1.aL
            alphaR1 = cb1.aR
            nL1     = cb1.nL
            nR1     = cb1.nR

        if cb2 and isinstance ( cb2 , CB2_pdf ) :
            mean2   = cb2.mean
            sigma2  = cb2.sigma
            alphaL2 = cb2.aL
            alphaR2 = cb2.aR
            nL2     = cb2.nL
            nR2     = cb2.nR

        ## save areguments 
        self.__cb1  = cb1
        self.__cb2  = cb2

        ## create helper Crystal Balls 
        self.__cb_1 = CB2_pdf ( name + '_X' , xvar , mean1 , sigma1 , alphaL1 , alphaR1 , nL1 , nR1 )
        self.__cb_2 = CB2_pdf ( name + '_Y' , yvar , mean2 , sigma2 , alphaL2 , alphaR2 , nL2 , nR2 )

        ## finally create PDF 
        self.pdf = Ostap.Models.Rotated2CrystalBalls (
            'r2cb_%s'         % name ,
            'Rotated2CBs(%s)' % name ,
            self.x      ,
            self.y      ,
            self.phi    ,
            ##
            self.mean1  ,
            self.sigma1 ,
            self.aL1    ,
            self.nL1    ,
            self.aR1    ,
            self.nR1    ,
            ##
            self.mean2  ,
            self.sigma2 ,
            self.aL2    ,
            self.nL2    ,
            self.aR2    ,
            self.nR2    )
            
        ## save configuration
        self.config = {
            'name'    : self.name   ,
            'xvar'    : self.xvar   ,
            'yvar'    : self.yvar   ,
            'phi'     : self.phi    ,
            ##
            'cb1'     : self.__cb1  ,
            'mean1'   : self.mean1  ,
            'sigma1'  : self.sigma1 ,
            'alphaL1' : self.aL1    ,
            'alphaR1' : self.aR1    ,
            'nL1'     : self.nL1    ,
            'nR1'     : self.nR1    ,
            ##
            'cb2'     : self.__cb2  ,
            'mean2'   : self.mean2  ,
            'sigma2'  : self.sigma2 ,
            'alphaL2' : self.aL2    ,
            'alphaR2' : self.aR2    ,
            'nL2'     : self.nL2    ,
            'nR2'     : self.nR2    }
        
    @property
    def phi ( self ) :
        "``phi''-rotation phase"
        return self.__phi
    @phi.setter
    def phi ( self, value ) :
        self.__phi.setVal ( float ( value ) ) 

    @property
    def mean1 ( self ) :
        "``mean1''-the mean-value of the first double-sided Crystal Ball"
        return self.__cb_1.mean 
    @mean1.setter
    def mean1 ( self, value ) :
        self.__cb_1.mean = value 

    @property
    def sigma1 ( self ) :
        "``sigma1''-the sigma-value of the first double-sided Crystal Ball"
        return self.__cb_1.sigma
    @sigma1.setter
    def sigma1 ( self, value ) :
        self.__cb_1.sigma = value

    @property
    def aL1 ( self ) :
        "``aL1''-the ``alpha-left''-value of the first double-sided Crystal Ball"
        return self.__cb_1.aL
    @aL1.setter
    def aL1 ( self, value ) :
        self.__cb_1.aL = value

    @property
    def aR1 ( self ) :
        "``aR1''-the ``alpha-right''-value of the first double-sided Crystal Ball"
        return self.__cb_1.aR
    @aR1.setter
    def aR1 ( self, value ) :
        self.__cb_1.aR = value

    @property
    def nL1 ( self ) :
        "``nL1''-the ``n-left''-value of the first double-sided Crystal Ball"
        return self.__cb_1.nL
    @nL1.setter
    def nL1 ( self, value ) :
        self.__cb_1.nL = value

    @property
    def nR1 ( self ) :
        "``nR1''-the ``n-right''-value of the first double-sided Crystal Ball"
        return self.__cb_1.nR
    @nR1.setter
    def nR1 ( self, value ) :
        self.__cb_1.nR = value


    @property
    def mean2 ( self ) :
        "``mean2''-the mean-value of the second double-sided Crystal Ball"
        return self.__cb_2.mean 
    @mean2.setter
    def mean2 ( self, value ) :
        self.__cb_2.mean = value 

    @property
    def sigma2 ( self ) :
        "``sigma2''-the sigma-value of the second double-sided Crystal Ball"
        return self.__cb_2.sigma
    @sigma2.setter
    def sigma2 ( self, value ) :
        self.__cb_2.sigma = value

    @property
    def aL2 ( self ) :
        "``aL2''-the ``alpha-left''-value of the second double-sided Crystal Ball"
        return self.__cb_2.aL
    @aL2.setter
    def aL1 ( self, value ) :
        self.__cb_2.aL = value

    @property
    def aR2 ( self ) :
        "``aR2''-the ``alpha-right''-value of the second double-sided Crystal Ball"
        return self.__cb_2.aR
    @aR2.setter
    def aR2 ( self, value ) :
        self.__cb_2.aR = value

    @property
    def nL2 ( self ) :
        "``nL2''-the ``n-left''-value of the second double-sided Crystal Ball"
        return self.__cb_2.nL
    @nL2.setter
    def nL2 ( self, value ) :
        self.__cb_2.nL = value

    @property
    def nR2 ( self ) :
        "``nR2''-the ``n-right''-value of the second double-sided Crystal Ball"
        return self.__cb_2.nR
    @nR2.setter
    def nR2 ( self, value ) :
        self.__cb_2.nR = value


models += [ Rotated2Gauss_pdf ]
models += [ Rotated2CB_pdf    ]
# =============================================================================
if '__main__' == __name__ : 
         
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )
    
# =============================================================================
# The END 
# =============================================================================
