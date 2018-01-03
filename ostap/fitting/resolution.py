#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file resolution.py
## Set of useful resolution models:
#  - single Gaussian
#  - double Gaussian                     (gaussian   tails)
#  - symmetric Apolonious                (exponenial tails)
#  - symmetric double-sided Crystal Ball (power-law  tails)
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2017-07-13
# =============================================================================
"""Set of useful resolution models:
- single Gaussian
- double Gaussian                     (gaussian   tails)
- symmetric Apolonious                (exponenial tails)
- symmetric double-sided Crystal Ball (power-law  tails)
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'ResoGauss'     , ## simple single-Gaussian resolution model,
    'ResoGauss2'    , ## double-Gaussian resolutin model,
    'ResoApo2'      , ## symmetric Apolonios resolution model,
    'ResoCB2'       , ## symmetric double-sided Crystal Ball resolution model,
    )
# =============================================================================
import ROOT
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.resolution' )
else                       : logger = getLogger ( __name__                   )
# =============================================================================
from ostap.fitting.basic import RESOLUTION, makeVar
# =============================================================================    
models = set() 
# =============================================================================
## sigle gaussian model for resolution
# =============================================================================
## @class ResoGauss
#  Trivial single gaussian resolution model
class ResoGauss(RESOLUTION) :
    """Trivial single gaussian resolution model
    """
    def __init__ ( self      ,
                   name      ,   ## the  name 
                   mass      ,   ## the variable 
                   sigma     ,   ## the first sigma 
                   mean  = 0 ) : ## mean-value
        ## initialize the base  
        super(ResoGauss,self).__init__( name  = name  ,
                                        mass  = mass  ,
                                        sigma = sigma ,
                                        mean  = mean  )
        
        ## build gaussian resolution model 
        self.gauss = ROOT.RooGaussModel(
            'ResoGauss_'    + name ,
            'ResoGauss(%s)' % name ,
            self.mass            ,
            self.mean            , 
            self.sigma           )
        
        ## the final PDF 
        self.pdf = self.gauss

models.add ( ResoGauss ) 
# =============================================================================
## @class ResoGauss2
#  Double Gaussian model for  resoltuion
class ResoGauss2(RESOLUTION) :
    """Double-gaussian resolution model
    - sigma of core Gaussian
    - ratio of wide/core widths
    - fraction of core component
    """        
    def __init__ ( self           ,
                   name           ,   ## the name 
                   mass           ,   ## the variable 
                   sigma          ,   ## the core sigma
                   scale    = 1.2 ,   ## sigma2/sigma1 ratio 
                   fraction = 0.5 ,   ## fraction of
                   mean     = 0.0 ) : ## the mean value
        ## initialize the base 
        super(ResoGauss2,self). __init__ ( name  = name  ,
                                           mass  = mass  ,
                                           sigma = sigma ,
                                           mean  = mean  )
        ## fraction of sigma1-component 
        self.__fraction = makeVar (
            fraction                   , 
            'CoreFraction_'     + name ,
            'CoreFraction(%s)'  % name , fraction , 0 ,  1 ) 
        
        ## sigma2/sigma1 width ratio;
        self.__scale = makeVar (
            scale ,
            'SigmaScale_'       + name ,
            'SigmaScale(%s)'    % name , scale    , 1 , 10 ) 
        
        from ostap.core.core import Ostap
        self.pdf = Ostap.Models.DoubleGauss (           
            "Reso2Gauss_"       + name ,
            "Reso2Gauss(%s)"    % name ,
            self.mass     ,
            self.sigma    ,
            self.fraction ,
            self.scale    ,
            self.mean    
            )

    @property
    def fraction ( self  ) :
        """``Fraction'' parameter for double Gaussian resolution function
        """
        return self.__fraction
    @fraction.setter
    def fraction ( self , value ) :
        value = float ( value )
        assert 0<= value <= 1, "``Fraction'' must be in  (0,1) range"
        self.__fraction.setVal ( value )
        return self.__fraction.getVal()

    @property
    def scale ( self  ) :
        """``Scale'' parameter for double Gaussian resolution function
        """
        return self.__scale
    @scale.setter
    def scale ( self , value ) :
        value = float ( value )
        assert 0 < value, "``Value'' must be >1"
        self.__scale.setVal ( value )
        return self.__scale.getVal()
    
models.add ( ResoGauss2 ) 
# =============================================================================
## @class ResoApo2
#  Symmetrical  Apolonios  model for resolution
class ResoApo2(RESOLUTION) :
    """Symmetric  Apolonios model for resolution
    """
    def __init__ ( self      ,
                   name      ,   ## the  name 
                   mass      ,   ## the variable 
                   sigma     ,
                   beta  = 1 ,   ## the first sigma
                   mean  = 0 ) : ## the mean value 

        ##  initlialize the base 
        super(ResoApo2,self).__init__ ( name  = name  ,
                                        mass  = mass  ,
                                        sigma = sigma ,
                                        mean  = mean  )
        self.__beta    = makeVar (
            beta ,
            'ResoBeta_%s'  % name  ,
            'ResoBeta(%s)' % name  , beta , 1.e-6  , 10000 )
        
        ## build (symmetric) resoltuion model
        from ostap.core.core import Ostap
        self.apo2  = Ostap.Models.Apolonios2 (
            "ResoApolonious_"   + name ,
            "ResoApolonios(%s)" % name ,
            self.mass   ,
            self.mean   ,
            self.sigma  ,
            self.sigma  ,
            self.beta   ) 

        self.pdf = self.apo2

    @property
    def beta ( self  ) :
        """``Beta'' parameter for Apolonious resolution function
        """
        return self.__beta
    @beta.setter
    def beta ( self , value ) :
        value = float ( value )
        assert 0<= value , "``Beta'' must be non-negative"
        self.__beta.setVal ( value )
        return self.__beta.getVal()

        
        ## 
models.add ( ResoApo2 ) 
# =============================================================================
## @class ResoCB2
#  Symmetrical double-sided Crystal Ball model for resolution
class ResoCB2(RESOLUTION) :
    """Symmetric double-sided Crystal Ball model for resolution
    """
    def __init__ ( self        ,
                   name        ,   ## the  name 
                   mass        ,   ## the  variable 
                   sigma       ,   ## core r esolution
                   alpha = 1.5 ,   ## alpha  
                   n     = 5   ,   ## power-law exponent
                   mean  = 0   ) : ## the mean value

        ## initialize the base 
        super(ResoCB2,self).__init__ ( name  = name  ,
                                       mass  = mass  ,
                                       sigma = sigma ,
                                       mean  = mean  )
            
        self.__alpha = makeVar (
            alpha                  ,
            'ResoAlpha_'    + name ,
            'ResoAlpha(%s)' % name , alpha , 0.5   ,  5)
        
        self.__n     = makeVar (
            n                  ,
            'ResoN_'        + name ,
            'ResoN(%s)'     % name , n     , 1.e-6 , 50 )
        
        ## 
        from ostap.core.core import Ostap
        self.cb2 = Ostap.Models.CrystalBallDS (  
            'ResoCB2_'   + name ,
            'ResoCB2(%s' % name ,
            self.mass           ,
            self.mean           , 
            self.sigma          ,
            self.alpha          ,
            self.n              ,
            self.alpha          ,
            self.n              )
        
        ## the final PDF 
        self.pdf = self.cb2

    @property
    def alpha ( self  ) :
        """``Alpha'' parameter for Doubble-sided symmetric resolution function
        """
        return self.__alpha
    @alpha.setter
    def alpha ( self , value ) :
        value = float ( value )
        assert 0.1<= value<=10 , "``Alpha'' must be in [0.1,10] interval"
        self.__alpha.setVal ( value )
        return self.__alpha.getVal()

    @property
    def n ( self  ) :
        """``n'' parameter for Doubble-sided symmetric resolution function
        """
        return self.__n
    @n.setter
    def n ( self , value ) :
        value = float ( value )
        assert 1.e-6 <= value <= 100 , "``n'' must be in [1.e-6,100] interval"
        self.__n.setVal ( value )
        return self.__n.getVal()

        
models.add ( ResoCB2 ) 
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
