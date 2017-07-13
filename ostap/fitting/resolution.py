#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file resolution.py
## Set of useful resoltuion models:
#  - single Gaussian
#  - double Gaussian                    (gaussian   tails)
#  - symmetric Apolonious               (exponenial tails)
#  - symmetric double-sided Crytal Ball (power-law  tails)
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2017-07-13
# =============================================================================
"""Set of useful resoltuion models:
- single Gaussian
- double Gaussian                    (gaussian   tails)
- symmetric Apolonious               (exponenial tails)
- symmetric double-sided Crytal Ball (power-law  tails)
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'ResoGauss'     , ## simple single-Gaussian resolutin model,
    'ResoGauss2'    , ## double-Gaussian resolutin model,
    'ResoApo2'      , ## symmetric Apolonios resolutin model,
    'ResoCB2'       , ## symmetric double-sided Crystal Ball resolutin model,
    )
# =============================================================================
import ROOT
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.resolution' )
else                       : logger = getLogger ( __name__                   )
# =============================================================================
models = [] 
# =============================================================================
## sigle gaussian model for resolution
# =============================================================================
import ostap.fitting.signals as OFS
class ResoGauss(OFS.Gauss_pdf) :
    """Trivial single gaussian resolution model
    """
    def __init__ ( self      ,
                   name      ,   ## the  name 
                   mass      ,   ## the variable 
                   sigma     ,   ## the first sigma 
                   mean  = 0 ) : ## mean-value 
        
        from ostap.fitting.basic import RESOLUTION
        with RESOLUTION() :            
            super(ResoGauss,self).__init__( name          ,
                                            mass  = mass  ,
                                            sigma = sigma ,
                                            mean  = mean  ) 
            PDF.__init__ ( self , name )
            
        ## gaussian 
        self.gauss = ROOT.RooGaussModel(
            'CoreGaussian'    + name ,
            'CoreGaussian(%s' % name ,
            self.mass                ,
            self.mean                , 
            self.sigma               )
        
        ## the final PDF 
        self.pdf = self.gauss

models.append ( ResoGauss ) 
# =============================================================================
## double gaussian model for  resoltuion
# =============================================================================
class ResoGauss2(OFS.DoubleGauss_pdf) :
    """Double-gaussian resolution model
    """        
    def __init__ ( self           ,
                   name           ,   ## the name 
                   mass           ,   ## the variable 
                   sigma          ,   ## the core sigma
                   scale    = 1.2 ,   ## sigma2/sigma1 ratio 
                   fraction = 0.5 ,   ## fraction of
                   mean     = 0.0 ) : ## the mean value 
        
        from ostap.fitting.basic import RESOLUTION
        with RESOLUTION() :            
            super(ResoGauss2,self). __init__ ( self                ,
                                               name                ,
                                               mass     = mass     ,
                                               sigma    = sigma    ,
                                               mean     = mean     ,
                                               fraction = fraction ,
                                               scale    = scale    )
            
models.append ( ResoGauss2 ) 
# =============================================================================
## Symmetrical  Apolonios  model for resolution
# =============================================================================
class ResoApo2(OFS.Apolonios2_pdf) :
    """Symmetric  Apolonios model for resolution
    """
    def __init__ ( self      ,
                   name      ,   ## the  name 
                   mass      ,   ## the variable 
                   sigma     ,
                   beta  = 1 ,   ## the first sigma
                   mean  = 0 ) : ## the mean value 

        from ostap.fitting.basic import RESOLUTION
        with RESOLUTION() :            
            super(ResoApo2,self).__init__ ( self               ,
                                            name               ,
                                            mass      = mass   ,
                                            sigma     = sigma  ,
                                            mean      = mean   ,
                                            beta      = beta   ,
                                            asymmetry = 0      )
            
models.append ( ResoApo2 ) 
# =============================================================================
## Symmetrical double-sided Crystal Ball model for resolution
# =============================================================================
class ResoCB2(OFS.Gauss_pdf) :
    """Symmetric double-sided Crystal Ball model for resolution
    """
    def __init__ ( self        ,
                   name        ,   ## the  name 
                   mass        ,   ## the  variable 
                   sigma       ,   ## core r esolution
                   alpha = 1.5 ,   ## alpha  
                   n     = 5   ,   ## power-law exponent
                   mean  = 0   ) : ## the mean value

        from ostap.fitting.basic import RESOLUTION
        with RESOLUTION() :            
            super(ResoCB2,self).__init__ ( self , name ,  mass , sigma , mean )
        
        self.alpha = makeVar (
            alpha                  ,
            'ResoAlpha'     + name ,
            'ResoAlpha(%s)' % name , alpha , 0.5 , 6 )
        
        self.n     = makeVar (
            n                  ,
            'ResoN'     + name ,
            'ResoN(%s)' % name , n , 1.e-4 , 20 )
        
        ## gaussian 
        self.cb2 = Ostap.Models.CrystalBallDS (  
            'ResoCB2'    + name ,
            'ResoCB2(%s' % name ,
            self.mass                ,
            self.mean                , 
            self.sigma               ,
            self.alpha               ,
            self.n                   ,
            self.alpha               ,
            self.n                   )
        
        ## the final PDF 
        self.pdf = self.cb2
        
models.append ( ResoCB2 ) 
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
