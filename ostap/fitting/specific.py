#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/specific.py
#  A set of predefined ready-to-use shapes and PDFs 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
# 
# =============================================================================
"""A set of predefined ready-to-use shapes and PDFs"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    #
    'D0_pdf'     , ## PDF for D0        : Bukin 
    'Dp_pdf'     , ## PDF for D+        : Bukin
    'Ds_pdf'     , ## PDF for Ds+       : Bukin 
    'Lc_pdf'     , ## PDF for Lambda_c+ : Gauss
    #
    'Bd_pdf'     , ## pdf for B0        : double-sided Crystal Ball 
    'B0_pdf'     , ## pdf for B0        : double-sided Crystal Ball 
    'Bu_pdf'     , ## pdf for B+        : double-sided Crystal Ball 
    'Bs_pdf'     , ## pdf for Bs        : double-sided Crystal Ball 
    'Bc_pdf'     , ## pdf for Bc+       : double-sided Crystal Ball 
    #
    'DpDs_pdf'   , ## ready-to-use model for D+ and D_s+ signals 
    'BdBs_pdf'   , ## ready-to-use model for Bd and Bs signals 
    'Manca_pdf'  , ## Manca function to fit Y->mu mu spectrum  [Y(1S),Y(2S),Y(3S)]
    'Manca2_pdf' , ## Manca function to fit Y->mu mu spectrum  [Y(1S),Y(2S),Y(3S)]
    'MancaX_pdf' , ## associative production of Y+X 
    #
    'HORNSdini_pdf' , ## HORNdini PDF (need to ve convolved!)
    'HILLdini_pdf'  , ## HILLdini PDF (need to ve convolved!)
    'HHdini_pdf'    , ## ready to use HORNSdini + HILLdini PDF
    )
# =============================================================================
from   ostap.core.core           import Ostap 
from   ostap.fitting.pdfbasic    import PDF1, PDF2, Sum1D
from   ostap.fitting.fit1d       import PEAK
from   ostap.fitting.signals     import CB2_pdf, Needham_pdf, Bukin_pdf
from   ostap.fitting.convolution import Convolution_pdf 
import ROOT, math
# =============================================================================
from   ostap.logger.logger     import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.specific' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================
models = [] 
# =============================================================================
## @class Bd_pdf
#  simple wrapper over CB2-pdf
#  @see Ostap::Models::CrystalBallDS
#  @see Ostap::Math::CrystalBallDS
#  @attention: mass is mandatory variable! 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Bd_pdf(CB2_pdf) :
    """ B0: double-sided Crystal Ball function
    """
    def __init__ ( self                   ,
                   xvar                   ,   ## mass is mandatory here! 
                   name      = 'Bd'       ,
                   mean      = 5.2791e+00 ,   ## to be released later 
                   sigma     = 7.2938e-03 ,   ## to be released later 
                   alphaL    = 1.4499e+00 ,   ## to be released later 
                   alphaR    = 1.9326e+00 ,   ## to be released later
                   nL        = 8.7234e+00 ,   ## to be released later 
                   nR        = 2.0377e+00 ) : ## to be released later 
        ## 
        CB2_pdf.__init__ ( self             ,
                           name             ,
                           xvar             , 
                           mean             ,
                           sigma            ,
                           alphaL           ,
                           alphaR           ,
                           nL               ,
                           nR               )
        
B0_pdf = Bd_pdf
models.append ( Bd_pdf ) 
# =============================================================================
## @class Bu_pdf
#  simple wrapper over CB2-pdf
#  @see Ostap::Models::CrystalBallDS
#  @see Ostap::Math::CrystalBallDS
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Bu_pdf(CB2_pdf) :
    """ B+: double-sided Crystal Ball function
    """
    def __init__ ( self                   ,
                   xvar                   ,   ## mass is mandatory here! 
                   name      = 'Bu'       ,
                   mean      = 5.2791e+00 ,   ## to be released later 
                   sigma     = 7.2938e-03 ,   ## to be released later 
                   alphaL    = 1.4499e+00 ,   ## to be released later 
                   alphaR    = 1.9326e+00 ,   ## to be released later
                   nL        = 8.7234e+00 ,   ## to be released later 
                   nR        = 2.0377e+00 ) : ## to be released later 
        ## 
        CB2_pdf.__init__ ( self             ,
                           name             ,
                           xvar             , 
                           mean             ,
                           sigma            ,
                           alphaL           ,
                           alphaR           ,
                           nL               ,
                           nR               )

models.append ( Bu_pdf ) 
# =============================================================================
## @class Bs_pdf
#  simple wrapper over CB2-pdf
#  @see Ostap::Models::CrystalBallDS
#  @see Ostap::Math::CrystalBallDS
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Bs_pdf(CB2_pdf) :
    """ Bs: double-sided Crystal Ball function
    """
    def __init__ ( self                   ,
                   xvar                   ,    ## mass is mandatory here! 
                   name      = 'Bs'       ,
                   mean      = 5.3661e+00 ,    ## to be released later 
                   sigma     = 7.2938e-03 ,    ## to be released later 
                   alphaL    = 1.4499e+00 ,    ## to be released later 
                   alphaR    = 1.9326e+00 ,    ## to be released later
                   nL        = 8.7234e+00 ,    ## to be released later 
                   nR        = 2.0377e+00 ) :  ## to be released later 
        ## 
        CB2_pdf.__init__ ( self             ,
                           name             ,
                           xvar             , 
                           mean             ,
                           sigma            ,
                           alphaL           ,
                           alphaR           ,
                           nL               ,
                           nR               )
        ## save configuration
        self.config = {
            'xvar'   : self.xvar  ,
            'name'   : self.name  ,
            'mean'   : self.mean  ,
            'sigma'  : self.sigma , 
            'alphaL' : self.aL    ,
            'alphaR' : self.aR    ,
            'nL'     : self.nL    ,
            'nR'     : self.nR    ,            
            }

models.append ( Bs_pdf )

# =============================================================================
## @class Bc_pdf
#  simple wrapper over CB2-pdf
#  @see Ostap::Models::CrystalBallDS
#  @see Ostap::Math::CrystalBallDS
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Bc_pdf(CB2_pdf) :
    """ Bc: double-sided Crystal Ball function
    """
    def __init__ ( self                   ,
                   xvar                   ,   ## mass is mandatory here! 
                   name      = 'Bc'       ,
                   mean      = 6.267e+00  ,   ## to be released later 
                   sigma     = 7.2938e-03 ,   ## to be released later 
                   alphaL    = 1.4499e+00 ,   ## to be released later 
                   alphaR    = 1.9326e+00 ,   ## to be released later
                   nL        = 8.7234e+00 ,   ## to be released later 
                   nR        = 2.0377e+00 ) : ## to be released later 
        ## 
        CB2_pdf.__init__ ( self             ,
                           name             ,
                           xvar             , 
                           mean             ,
                           sigma            ,
                           alphaL           ,
                           alphaR           ,
                           nL               ,
                           nR               )

models.append ( Bc_pdf ) 
# =============================================================================
# Specializations for Bukin function
# =============================================================================
## @class D0_pdf
#  simple wrapper over Bukin-pdf
#  @see RooBukinPdf
#  @see Ostap::Math::Bukin
#  @see Analusis::Models::Bukin

#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class D0_pdf(Bukin_pdf) :
    """D0: Bukin function 
    """
    def __init__ ( self                 ,
                   xvar                 , ## mass is mandatory here! 
                   name   = 'D0'        ,
                   mean   =  1.8648e+00 , 
                   sigma  =  7.3651e-03 , 
                   xi     = -1.7616e-03 , 
                   rhoL   =  3.7311e-01 ,
                   rhoR   =  4.4033e-01 ) : 
        
        Bukin_pdf.__init__ ( self          ,
                             name          ,
                             xvar          , 
                             mean          ,
                             sigma         ,
                             xi            , 
                             rhoL          ,
                             rhoR          ) 
                             
models.append ( D0_pdf ) 
# =============================================================================
## @class Dp_pdf
#  simple wrapper over Bukin-pdf
#  @see Ostap::Math::Bukin
#  @see Analusis::Models::Bukin
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Dp_pdf(Bukin_pdf) :
    """ D+: Bukin function 
    """
    def __init__ ( self                    ,
                   xvar                    , ## mass is mandatory here 
                   name     = 'Dp'         ,
                   mean     =  1.869       , 
                   sigma    =  7.1183e-03  ,
                   xi       = -7.7344e-03  ,
                   rhoL     =  3.0241e-01  , 
                   rhoR     =  3.7452e-01  ) :

        Bukin_pdf.__init__ ( self          ,
                             name          ,
                             xvar          ,
                             mean          ,
                             sigma         ,
                             xi            ,                            
                             rhoL          ,
                             rhoR          ) 
        
models.append ( Dp_pdf ) 

# =============================================================================
## @class Ds_pdf
#  simple wrapper over Bukin-pdf
#  @see Ostap::Math::Bukin
#  @see Analusis::Models::Bukin
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Ds_pdf(Bukin_pdf) :
    """ Ds: Bukin function 
    """
    def __init__ ( self                    , 
                   xvar                    , ## mass is mandatory 
                   name     = 'Ds'         ,
                   mean     =  1.969       ,
                   sigma    =  0.0068      ,
                   xi       = -6.45755e-04 ,
                   rhoL     =  3.0241e-01  , 
                   rhoR     =  3.7452e-01  ) :
        
        Bukin_pdf.__init__ ( self          ,
                             name          ,
                             xvar          , 
                             mean          ,
                             sigma         ,
                             xi            ,
                             rhoL          ,
                             rhoR          ) 
        
models.append ( Ds_pdf ) 

# =============================================================================
## @class Lc_pdf
#  simple wrapper over Bukin-pdf
#  @see Ostap::Math::Bukin
#  @see Analusis::Models::Bukin
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Lc_pdf(Bukin_pdf) :
    """ Lc: Bukin function 
    """
    def __init__ ( self                     ,
                   xvar                     , 
                   name     = 'Lc'          ,                   
                   mean     =  2.28590e+00  ,
                   sigma    =  5.11874e-03  ,
                   xi       =  1.82493e-03  ,
                   rhoL     =  3.0241e-01   , 
                   rhoR     =  3.7452e-01   ) :
        
        Bukin_pdf.__init__ ( self     ,
                             name     ,
                             xvar     , 
                             mean     ,
                             sigma    ,
                             xi       ,
                             rhoL     ,
                             rhoR     ) 
models.append ( Lc_pdf ) 

# =============================================================================
## Ready-to-use fit model 
# =============================================================================

# =============================================================================
## @class BdBs_pdf
#  Ready-to-use model to fit a distribution with Bd and Bs peaks 
#  - Bs mass is calculated as  <code>m(Bd)+delta_m</code>
#  - Bs resolution is calculated as <code>sigma(Bd)*(m(Bs)/m(Bd)))</code>
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-08-02
class BdBs_pdf (PEAK) :
    """ Ready-to-use model to fit a distribution iwth D+ and D_s+ peaks
      - Bs mass is calculated as  `m(Bd)+delta_m`
      - Bs resolution is calculated as `sigma(Bd)*(m(Bs)/m(Bd))`
    """
    def __init__ ( self                               ,
                   xvar                               ,
                   name       = 'BdBd'                ,
                   mBd        = 5.287                 ,  ## location of Bd peak 
                   sBd        = 7.2938e-03            ,  ## to be released later 
                   dm         = 1.969 - 1.869         ,  ## mass differenc:  m(Bs) - m(Bd)
                   alphaL     = 1.4499e+00            ,  ## to be released later 
                   alphaR     = 1.9326e+00            ,  ## to be released later
                   nL         = 8.7234e+00            ,  ## to be released later 
                   nR         = 2.0377e+00            ,  ## to be released later 
                   background = 'convex decreasing 2' , 
                   NBd        = None                  ,
                   NBs        = None                  ,
                   B          = None                  ,
                   fix_norm   = False                 ) :
        
        ## initialize the base
        PEAK.__init__      ( self                         ,
                             name       = name            ,
                             xvar       = xvar            , 
                             mean       = mBd             , 
                             sigma      = sBd             ,
                             mean_name  = 'mBd_%s' % name ,
                             sigma_name = 'sBd_%s' % name )
        
        self.__Bd = Bd_pdf ( xvar   = self.xvar           , 
                             name   = 'Bd_%s' % self.name ,
                             mean   = self.mean           ,
                             sigma  = self.sigma          ,
                             alphaL = alphaL              ,
                             alphaR = alphaR              ,
                             nL     = nL                  ,
                             nR     = nR                  )
        
        self.__dm  = self.make_var ( dm                      ,
                                     "dm_"     + self.name   ,
                                     "dm(%s)"  % self.name   ,
                                     False                   )
        
        self.__mBs = Ostap.MoreRooFit.Addition ( "mBs"     + self.name ,
                                                 "mBs(%s)" % self.name ,
                                                 self.mBd  , self.dm   ) 
        
        ## resoltuion of Bs is a scaled resolution of Bd by m(Bs)/m(Bd) 
        self.__sBs    = Ostap.MoreRooFit.ABC ( 
            "sBs_"            + self.name  ,
            "#sigma_{Bs}(%s)" % self.name  ,
            self.sBd   ,
            self.__mBs ,
            self.mBd   )        
        
        self.__Bs = Bs_pdf ( xvar  = self.xvar           , 
                             name  = 'Bs_%s' % self.name ,
                             mean  = self.__mBs          ,
                             sigma = self.__sBs          ,
                             alphaL = self.Bd.aL         ,
                             alphaR = self.Bd.aR         ,
                             nL     = self.Bd.nL         ,
                             nR     = self.Bd.nR         )
    
        # ===============================================================================
        ## use helper function to create the background  
        # ===============================================================================
        self.__bkg   = self.make_bkg ( background , 'Bkg%s' % self.name , self.xvar )

        # ===============================================================================
        ## Yields:
        # ===============================================================================        
        self.__nBd = self.make_var ( NBd  ,
                                     "NBd" + name  ,
                                     "Signal(Bd)"  ,  None , 1000 ,  0 ,  5.e+7 )
        self.__nBs = self.make_var ( NBs  ,
                                     "NBs" + name  ,
                                     "Signal(Bs)" ,  None , 1000 ,  0 ,  5.e+7 )
        self.__b   = self.make_var ( B     ,
                                     "B"   + name  ,
                                     "Background"  ,  None ,  100 ,  0 ,  9.e+8 )

        # ===============================================================================
        ## create the final PDF 
        # ===============================================================================
        self.alist1 = ROOT.RooArgList ( self.Bd.pdf , self.Bs.pdf , self.bkg.pdf ) 
        self.alist2 = ROOT.RooArgList ( self.NBd    , self.NBs    , self.B       ) 
        self.pdf    = ROOT.RooAddPdf  (
            "BdBs_%s"  % self.name ,
            "BdBs(%s)" % self.name ,
            self.alist1 ,
            self.alist2 )
        
        if fix_norm : self.pdf.fixCoefNormalization ( self.vars ) ## VB: added 10/10/2024 to suppress warnings 
        
        # ===============================================================================
        ## finally declare the components 
        # ===============================================================================
        self.signals     . add ( self.Bd.pdf  )
        self.signals     . add ( self.Bs.pdf  )
        self.backgrounds . add ( self.bkg.pdf ) 

        ## save configuration
        self.config = {
            'xvar'        : self.xvar     ,
            'name'        : self.name     ,
            'background'  : self.bkg      ,
            'mBd'         : self.mBd      ,
            'sBd'         : self.sBd      ,
            ##
            'alphaL'      : self.Bd.aL    ,
            'alphaR'      : self.Bd.aR    ,
            'nL'          : self.Bd.nL    ,
            'nR'          : self.Bd.nR    ,
            ##
            'dm'          : self.dm       , 
            ## 
            'NBd'         : self.NBd      ,
            'NBs'         : self.NBs      ,
            'B'           : self.B        ,
            'fix_norm'    : self.fix_norm }
    
    @property
    def Bd ( self ) :
        """'Bd' : model for the Bd peak"""
        return self.__Bd

    @property
    def Bs ( self ) :
        """'Bs' : model for the Bs peak"""
        return self.__Bs

    @property
    def mBd ( self ) :
        """''mBd' : mass/location for D+ peak (same as 'mean') """
        return self.mean
    @mBd.setter
    def mBd  ( self , value ) :
        self.mean = value 

    @property
    def sBd ( self ) :
        """''sBd' : resolution for the Bd peak (same as 'sigma') """
        return self.sigma
    @sBd.setter
    def sBd ( self , value ) : 
        self.sigma = value

    @property
    def dm ( self )  :
        """'dm' : mass difference between Bs and Bs mesons"""
        return self.__dm
    @dm.setter
    def dm ( self , value ) :
        self.set_value ( self.__dm , value )
        
    @property
    def NBd ( self ) :
        """'NBd'-parameter: yield of Bd"""
        return self.__nBd
    @NBd.setter
    def NBd ( self , value ) :
        self.set_value ( self.__nBd , value )

    @property
    def NBs ( self ) :
        """'NBs'-parameter: yield of Bs"""
        return self.__nBs
    @NBs.setter
    def NBs ( self , value ) :
        self.set_value ( self.__nBs , value )
        
    @property
    def B ( self ) :
        """'B'-parameter: yield of background component"""
        return self.__b
    @B.setter
    def B ( self , value ) :
        self.set_value ( self.__b , value )

    @property
    def bkg ( self ) :
        """'background' : get the background shape"""
        return self.__bkg 

    @property
    def fix_norm ( self ) :
        """`fix-norm`: 
        - see `ROOT.RooAbsPdf.SetCoefNormalization`
        - see `ROOT.RooAbsPdf.getCoefNormalization`
        """
        pars = self.pdf.getCoefNormalization()
        return True if pars else False 
    
models.append ( BdBs_pdf ) 

# =============================================================================
## @class DpDs_pdf
#  Ready-to-use model to fit a distribution with D+ and Ds+ peaks 
#  - Ds+ mass is calculated as  <code>m(D+)+delta_m</code>
#  - Ds+ resolution is calculates as <code>sigma(D+)*(m(Ds+)/m(D+)))</code>
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-08-02
class DpDs_pdf (PEAK) :
    """ Ready-to-use model to fit a distribution with D+ and Ds+ peaks
      - Ds+ mass is calculated as  `m(D+)+delta_m`
      - Ds+ resolution is calculated as `sigma(D+)*(m(Ds+)/m(D+))`
    """
    def __init__ ( self                               ,
                   xvar                               ,
                   name       =  'DpDs'               ,
                   mDp        =  1.869                ,  ## location of Bd peak
                   dm         =  1.969 - 1.869        ,  ## mass differenc:  m(Ds+) - m(D+)                   
                   sDp        =  7.1183e-03           ,
                   xi         = -7.7344e-03           ,
                   rhoL       =  3.0241e-01           ,
                   rhoR       =  3.7452e-01           ,                    
                   background = 'convex decreasing 2' , 
                   NDp        = None                  ,
                   NDs        = None                  ,
                   B          = None                  ,
                   fix_norm   = False                 ) :
        
        ## initialize the base
        PEAK.__init__      ( self                         ,
                             name       = name            ,
                             xvar       = xvar            , 
                             mean       = mDp             , 
                             sigma      = sDp             ,
                             mean_name  = 'mDp_%s' % name ,
                             sigma_name = 'sDp_%s' % name )
        
        self.__Dp = Dp_pdf ( xvar   = self.xvar           , 
                             name   = 'Dp_%s' % self.name ,
                             mean   = self.mean           ,
                             sigma  = self.sigma          ,
                             xi     = xi                  ,
                             rhoL   = rhoL                ,
                             rhoR   = rhoR                )
        
        self.__dm  = self.make_var ( dm                      ,
                                     "dm_"     + self.name   ,
                                     "dm(%s)"  % self.name   ,
                                     False                   )
        
        self.__mDs = Ostap.MoreRooFit.Addition ( "mDs"     + self.name ,
                                                 "mDs(%s)" % self.name ,
                                                 self.mDp  , self.dm   ) 
        
        ## resolution of Ds+ is a scaled resolution of D+ by m(Ds+)/m(D+) 
        self.__sDs    = Ostap.MoreRooFit.ABC ( 
            "sDs_"            + self.name  ,
            "#sigma_{Ds}(%s)" % self.name  ,
            self.sDp   ,
            self.__mDs ,
            self.mDp   )        
        
        self.__Ds = Ds_pdf ( xvar   = self.xvar           , 
                             name   = 'Ds_%s' % self.name ,
                             mean   = self.__mDs          ,
                             sigma  = self.__sDs          ,
                             xi     = self.Dp.xi          ,
                             rhoL   = self.Dp.rhoL        ,
                             rhoR   = self.Dp.rhoR        )
    
        # ===============================================================================
        ## use helper function to create the background  
        # ===============================================================================
        self.__bkg   = self.make_bkg ( background , 'Bkg%s' % self.name , self.xvar )

        # ===============================================================================
        ## Yields:
        # ===============================================================================        
        self.__nDp = self.make_var ( NDp  ,
                                     "NDp" + name  ,
                                     "Signal(Dp)"  ,  None , 1000 ,  0 ,  5.e+7 )
        self.__nDs = self.make_var ( NDs  ,
                                     "NDs" + name  ,
                                     "Signal(Ds)" ,  None , 1000 ,  0 ,  5.e+7 )
        self.__b   = self.make_var ( B     ,
                                     "B"   + name  ,
                                     "Background"  ,  None ,  100 ,  0 ,  9.e+8 )

        # ===============================================================================
        ## create the final PDF 
        # ===============================================================================
        self.alist1 = ROOT.RooArgList ( self.Dp.pdf , self.Ds.pdf , self.bkg.pdf ) 
        self.alist2 = ROOT.RooArgList ( self.NDp    , self.NDs    , self.B       ) 
        self.pdf    = ROOT.RooAddPdf  (
            "DpDs_%s"  % self.name ,
            "DpDs(%s)" % self.name ,
            self.alist1 ,
            self.alist2 )

        if fix_norm : self.pdf.fixCoefNormalization ( self.vars ) ## VB: added 10/10/2024 to suppress warnings 

        # ===============================================================================
        ## finally declare the components 
        # ===============================================================================
        self.signals     . add ( self.Dp.pdf  )
        self.signals     . add ( self.Ds.pdf  )
        self.backgrounds . add ( self.bkg.pdf ) 

        ## save configuration
        self.config = {
            'xvar'        : self.xvar     ,
            'name'        : self.name     ,
            'background'  : self.bkg      ,
            'mDp'         : self.mDp      ,
            'sDp'         : self.sDp      ,
            ##
            'xi'          : self.Dp.xi    ,
            'rhoL'        : self.Dp.rhoL  ,
            'rhoR'        : self.Dp.rhoR  ,
            ##
            'dm'          : self.dm       , 
            ## 
            'NDp'         : self.NDp      ,
            'NDs'         : self.NDs      ,
            'B'           : self.B        ,
            'fix_norm'    : self.fix_norm }
    
    @property
    def Dp ( self ) :
        """'Dp' : model for the D+ peak"""
        return self.__Dp

    @property
    def Ds ( self ) :
        """'Ds' : model for the Ds+ peak"""
        return self.__Ds

    @property
    def mDp ( self ) :
        """''mDp' : mass/location for D+ peak (same as 'mean') """
        return self.mean
    @mDp.setter
    def mDp  ( self , value ) :
        self.mean = value 

    @property
    def sDp ( self ) :
        """''sDp' : resolution for the D+ peak (same as 'sigma') """
        return self.sigma
    @sDp.setter
    def sDp ( self , value ) : 
        self.sigma = value

    @property
    def dm ( self )  :
        """'dm' : mass difference between Ds+ and D+ mesons"""
        return self.__dm
    @dm.setter
    def dm ( self , value ) :
        self.set_value ( self.__dm , value )
        
    @property
    def NDp ( self ) :
        """'NDp'-parameter: yield of D+"""
        return self.__nDp
    @NDp.setter
    def NDp ( self , value ) :
        self.set_value ( self.__nDp , value )

    @property
    def NDs ( self ) :
        """'NDs'-parameter: yield of Ds+"""
        return self.__nDs
    @NDs.setter
    def NDs ( self , value ) :
        self.set_value ( self.__nDs , value )
        
    @property
    def B ( self ) :
        """'B'-parameter: yield of background component"""
        return self.__b
    @B.setter
    def B ( self , value ) :
        self.set_value ( self.__b , value )

    @property
    def bkg ( self ) :
        """'background' : get the background shape"""
        return self.__bkg
    
    @property
    def fix_norm ( self ) :
        """`fix-norm`: 
        - see `ROOT.RooAbsPdf.SetCoefNormalization`
        - see `ROOT.RooAbsPdf.getCoefNormalization`
        """
        pars = self.pdf.getCoefNormalization()
        return True if pars else False 

    
models.append ( DpDs_pdf ) 

# =============================================================================
## @class MANCA
#  Helper base class for implementaton of PDFs for Y-> mu+mu- peaks
#  It helps to make a connection between the three Upsulon peaks,
#  starting form the Y(1S) parameters 
#  - Y(2S) mass is calculated as  <code>m(Y(1S))+dm21</code>
#  - Y(2S) mass is calculated as  <code>m(Y(2S))+dm32</code>
#  - Y(2S) resolution is calculated as <code>sigma(Y(1S))*(m(Y(2S))/m(Y(1S)))</code>
#  - Y(3S) resolution is calculated as <code>sigma(Y(1S))*(m(Y(3S))/m(Y(1S)))</code>
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-08-02
class MANCA(PEAK) :
    """ Helper base class for implementaton of PDFs for Y-> mu+mu- peaks
    - It helps to make a connection between three Upsilons peaks:
    starting form the Y(1S) parameters 
    - Y(2S) mass is calculated as `m(Y(1S))+dm21`
    - Y(2S) mass is calculated as  `m(Y(2S))+dm32`
    - Y(2S) resolution is calculated as `sigma(Y(1S))*(m(Y(2S))/m(Y(1S)))`
    - Y(3S) resolution is calculated as `sigma(Y(1S))*(m(Y(3S))/m(Y(1S)))`
    """
    ## create MancaBase 
    #  @param xvar observable, mu+mu- mass
    #  @param name       the name of the PDF
    #  @param m1s        mass/location of Y(1S) peak 
    #  @param s1s        resolution parameter for Y(1S) peak  
    #  @param dm21        mass difference between Y(2S) and Y(1S) states 
    #  @param dm32        mass difference between Y(3S) and Y(2S) states
    #  @param  background backrgudn shape: symbols, PDF or ROOT.RooAbsPdf
    #  @param N1S         yield for the Y(1S) signal
    #  @param N2S         yield for the Y(2S) signal
    #  @param N3S         yield for the Y(3S) signal
    #  @param B           yield for the background component    
    def __init__  ( self  ,
                    xvar               ,   ## the observable : mu+mu- mass
                    name       = 'Y'   ,   ## PDF name 
                    m1s        = None  ,   ## mass/location      for the Y(!S) peak 
                    s1s        = None  ,   ## resution parameter for the Y(1S) peak
                    ##
                    dm21       = None  ,   ## mass difference bewween Y(2S) and Y(1S)
                    dm32       = None  ,   ## mass difference bewween Y(3S) and Y(2S)
                    ##
                    background = 'd3'  , ## decreasing 3rd order polyminial 
                    ##
                    N1S        = None  ,   ## Y(1S) signal
                    N2S        = None  ,   ## Y(2S) signal
                    N3S        = None  ,   ## Y(3S) signal
                    B          = None  ) :  ## background signal
        """
        - xvar       : observable, mu+mu- mass
        - name       : the name of the PDF
        - m1s        : mass/location of Y(1S) peak 
        - s1s        : resolution parameter for Y(1S) peak  
        - dm21       : mass difference between Y(2S) and Y(1S) states 
        - dm32       : mass difference between Y(3S) and Y(2S) states 
        - N1S        : yield for the Y(1S) signal
        - N2S        : yield for the Y(2S) signal
        - N3S        : yield for the Y(3S) signal
        - B          : yield for the background component    
        """
        ## initialize the base
        PEAK.__init__ ( self                         ,
                        name       = name            ,
                        xvar       = xvar            ,
                        mean       = m1s             , 
                        sigma      = s1s             ,
                        mean_name  = 'm1s_%s' % name ,
                        sigma_name = 's1s_%s' % name )
        
        if   9460. in self.xvar and 10023. in self.xvar and 10355. in self.xvar : self._gev = 1000        
        elif 9.460 in self.xvar and 10.023 in self.xvar and 10.355 in self.xvar : self._gev = 1
        else :
            raise AttributeError ( "Illegal mass range %s<m<%s" % self.xminmax()  ) 

        m_y1s  =   9.46030     * self._gev
        m_y2s  =  10.02328     * self._gev
        m_y3s  =  10.3552      * self._gev
        
        s_y1s  =  4.03195e-02  * self._gev 
        dm_y21 = m_y2s - m_y1s
        dm_y32 = m_y3s - m_y2s

        ## mass differences for the Y(2S) & Y(3S) states 
    
        self.__dm21  = self.make_var ( dm21                  ,
                                       "dm21"      + name    ,
                                       "dm21(%s)"  % name    ,
                                       False                 , 
                                       dm_y21                ,
                                       dm_y21 - 0.10 * s_y1s , 
                                       dm_y21 + 0.10 * s_y1s )
        
        self.__dm32  = self.make_var ( dm32                  ,
                                       "dm32"      + name    ,
                                       "dm32(%s)"  % name    ,
                                       False                 , 
                                       dm_y32                ,
                                       dm_y32 - 0.10 * s_y1s , 
                                       dm_y32 + 0.10 * s_y1s )
        
        ## masses for the Y(2S) & Y(3S) states 

        self.__m2s    = Ostap.MoreRooFit.Addition ( "m_" + name + '2S'   ,
                                                    "m2s(%s)"  % name    ,
                                                    self.m1s , self.dm21 ) 

        self.__m3s    = Ostap.MoreRooFit.Addition ( "m_" + name + '3S'   ,
                                                    "m3s(%s)"  % name    ,
                                                    self.m2s , self.dm32 ) 

        ## resolution parameters for the Y(2S) and Y(3S) states
       
        ## width of Y(2S) is a scaled width of Y(1S) by m(2S)/m(1S) 
        self.__s2s    = Ostap.MoreRooFit.ABC ( 
            "sigma_"  + name + '2S'    ,
            "#sigma_{Y2S}(%s)" % name  ,
            self.s1s ,
            self.m2s ,
            self.m1s )        

        ## width of Y(3S) is a scaled width of Y(1S) by m(3S)/m(1S) 
        self.__s3s    = Ostap.MoreRooFit.ABC ( 
            "sigma_"  + name + '3S'    ,
            "#sigma_{Y3S}(%s)" % name  ,
            self.s1s ,
            self.m3s ,
            self.m1s )

        ## Yields:
        
        self.__n1s = self.make_var ( N1S  ,
                                     "N1S" + name  ,
                                     "Signal(Y1S)" ,  None , 1000 ,  0 ,  5.e+7 )
        self.__n2s = self.make_var ( N2S  ,
                                     "N2S" + name  ,
                                     "Signal(Y2S)" ,  None ,  300 ,  0 ,  5.e+6 )
        self.__n3s = self.make_var ( N3S  ,
                                     "N3S" + name  ,
                                     "Signal(Y3S)" ,  None ,  100 ,  0 ,  2.e+6 )
        self.__b   = self.make_var ( B     ,
                                     "B"   + name  ,
                                     "Background"  ,  None ,  100 ,  0 ,  9.e+8 )

        # ===============================================================================
        ## use helper function to create the background  
        # ===============================================================================
        self.__bkg   = self.make_bkg ( background , 'Bkg%s' % name , self.mass )

    @property
    def  m1s (  self ) :
        """'m1s'-parameter, location of the  Y(1S) peak, same as 'mean'"""
        return self.mean
    @m1s.setter
    def  m1s (  self , value ) :
        self.mean = value 
    
    @property
    def  s1s (  self ) :
        """'s1s'-parameter (resolution for Y(1S) peak) (same as 'sigma')"""
        return self.sigma
    @s1s.setter
    def  s1s (  self , value ) :
        self.sigma = value 

    @property
    def  dm21 (  self ) :
        """'dm21'-parameter (mass  difference for Y(2S) and Y(1S) peaks)"""
        return self.__dm21
    @dm21.setter
    def  dm21 (  self , value ) :
        self.set_value ( self.__dm21 , value ) 
        
    @property
    def  dm32 (  self ) :
        """'dm32'-parameter (mass  difference for Y(3S) and Y(2S) peaks)"""
        return self.__dm32    
    @dm32.setter
    def  dm32 (  self , value ) :
        self.set_value ( self.__dm32 , value )

    @property
    def  m2s (  self ) :
        """'m2s'-parameter (mass/location for the Y(2S) peak)"""
        return self.__m2s
    @property
    def  m3s (  self ) :
        """'m3s'-parameter (mass/location for the Y(3S) peak)"""
        return self.__m3s

    @property
    def  s2s (  self ) :
        """'s2s'-parameter (resolution for Y(2S) peak)"""
        return self.__s2s
    @property
    def  s3s (  self ) :
        """'s3s'-parameter (resolution for Y(3S) peak)"""
        return self.__s3s
        
    @property
    def N1S ( self ) :
        """'N1S'-parameter: yield of Y(1S)"""
        return self.__n1s
    @N1S.setter
    def N1S ( self , value ) :
        self.set_value ( self.__n1s , value )

    @property
    def N2S ( self ) :
        """'N2S'-parameter: yield of Y(2S)"""
        return self.__n2s
    @N2S.setter
    def N2S ( self , value ) :
        self.set_value ( self.__n2s , value )

    @property
    def N3S ( self ) :
        """'N3S'-parameter: yield of Y(3S)"""
        return self.__n3s
    @N3S.setter
    def N3S ( self , value ) :
        self.set_value ( self.__n3s , value )

    @property
    def B ( self ) :
        """'B'-parameter: yield of background component"""
        return self.__b
    @B.setter
    def B ( self , value ) :
        self.set_value ( self.__b , value )

    @property
    def bkg ( self ) :
        """'background' : get the background shape"""
        return self.__bkg 
    
    @property     
    def fix_norm ( self ) :
        """`fix-norm`: 
        - see `ROOT.RooAbsPdf.SetCoefNormalization`
        - see `ROOT.RooAbsPdf.getCoefNormalization`
        """
        pars = self.pdf.getCoefNormalization()
        return True if pars else False 
       
# =============================================================================
## @class Manca_pdf 
#  the final full PDF for Y->mu+mu- fit
#  This is physically well-motivated function for fits in narrow
#  bins in pt and rapidity
#
#  Parameters of three Upsilon peaks are related:
#  - Y(2S) mass is calculated as  <code>m(Y(1S))+dm21</code>
#  - Y(2S) mass is calculated as  <code>m(Y(2S))+dm32</code>
#  - Y(2S) resolution is calculated as <code>sigma(Y(1S))*(m(Y(2S))/m(Y(1S)))</code>
#  - Y(3S) resolution is calculated as <code>sigma(Y(1S))*(m(Y(3S))/m(Y(1S)))</code>
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-08-02
class Manca_pdf (MANCA) :
    """ Manca: the final fit model for Y->mu+mu- fit
    This is physically well-motivated function for fits in narrow bins in pt and rapidity
    - three Needham functions for Y(1S), Y(2S) and Y(3S) peaks
    - constrants for their resolutions and masses 
    - background: ... 

    Parameters of three Upsilon peaks are related:
    - Y(2S) mass is calculated as `m(Y(1S))+dm21`
    - Y(2S) mass is calculated as  `m(Y(2S))+dm32`
    - Y(2S) resolution is calculated as `sigma(Y(1S))*(m(Y(2S))/m(Y(1S)))`
    - Y(3S) resolution is calculated as `sigma(Y(1S))*(m(Y(3S))/m(Y(1S)))`
    """
    # ===========================================================================
    ## create Manca PDF 
    #  @param  xvar observable, mu+mu- mass
    #  @param  name       the name of the PDF
    #  @param  background backrgudn shape: symbols, PDF or ROOT.RooAbsPdf
    #  @param  m1s        mass/location of Y(1S) peak 
    #  @param  s1s        resolution parameter for Y(1S) peak    
    #  @param a0          a0-parameter for Needham function
    #  @param a1          a1-parameter for Needham function
    #  @param a2          a2-parameter for Needham function
    #  @param dm21        mass difference between Y(2S) and Y(1S) states 
    #  @param dm32        mass difference between Y(3S) and Y(2S) states 
    #  @param N1S         yield for the Y(1S) signal
    #  @param N2S         yield for the Y(2S) signal
    #  @param N3S         yield for the Y(3S) signal
    #  @param B           yield for the background component        
    def __init__ ( self               ,
                   xvar               ,
                   name        = 'Y'  ,
                   background  = 0    ,                   
                   m1s         = None ,
                   s1s         = None ,
                   ## 
                   c0          = 1.91 ,
                   c1          = None ,
                   c2          = None ,
                   ## 
                   dm21       = None  ,   ## m(2S) - m(1S) 
                   dm32       = None  ,   ## m(3S) - m(2S)
                   ##
                   N1S        = None  ,   ## Y(2S) signal
                   N2S        = None  ,   ## Y(2S) signal
                   N3S        = None  ,   ## Y(3S) signal 
                   B          = None  ,   ## background
                   fix_norm   = False ) :
        
        """ Create Manca1 PDF: function to fit Y->mu+mu- signals in narrow kinematic range 
        - xvar       : observable, mu+mu- mass
        - name       : the name of the PDF
        - background : background  shape: symbol, PDF or ROOT.RooAbsPdf        
        - m1s        : mass/location of Y(1S) peak 
        - s1s        : resolution parameter for Y(1S) peak  
        - c0         : c0-parameter for Needham function
        - c1         : c1-parameter for Needham function
        - c2         : cc2-parameter for Needham function
        - dm21       : mass difference between Y(2S) and Y(1S) states 
        - dm32       :  mass difference between Y(3S) and Y(2S) states 
        - N1S        : yield for the Y(1S) signal
        - N2S        : yield for the Y(2S) signal
        - N3S        : yield for the Y(3S) signal
        - B          : yield for the background component    
        
        Parameters of three Upsilon peaks are related:
        - Y(2S) mass is calculated as `m(Y(1S))+dm21`
        - Y(2S) mass is calculated as  `m(Y(2S))+dm32`
        - Y(2S) resolution is calculated as `sigma(Y(1S))*(m(Y(2S))/m(Y(1S)))`
        - Y(3S) resolution is calculated as `sigma(Y(1S))*(m(Y(3S))/m(Y(1S)))`
        """
        ### initialize the base 
        MANCA.__init__ ( self                    ,
                         name       = name       ,
                         xvar       = xvar       ,
                         m1s        = m1s        ,
                         s1s        = s1s        ,
                         dm21       = dm21       ,
                         dm32       = dm32       ,
                         background = background , 
                         N1S        = N1S        , 
                         N2S        = N2S        , 
                         N3S        = N3S        , 
                         B          = B          )
        
        # =====================================================================
        ## Shape parameters
        # =====================================================================

        self.__c0   = self.make_var ( c0                 ,
                                      'c0m_%s' % name    ,
                                      "c0 for Needham's function" ,
                                      True  , 1.91  , 0.1 , 3.0   )
        
        self.__c1   = self.make_var ( c1                 ,
                                      'c1m_%s' % name    ,
                                      "c1 for Needham's function" ,
                                      True , 1.1174 / self._gev ,  -10.0 / self._gev , 10.0 / self._gev )
        
        self.__c2   = self.make_var ( c2                 ,
                                      'c2m_%s' % name    ,
                                      "c2 for Needham's function" , 
                                      True , -5.299 / self._gev**2   , -100.0  / self._gev**2  , 100.0  / self._gev**2 )
            
        # =====================================================================
        ## Y(1S)
        # =====================================================================
        self.__Y1S   = Needham_pdf (
            name + '1S'           ,
            xvar     = self.xvar  ,
            mean     = self.m1s   ,
            sigma    = self.s1s   ,
            c0       = self.c0    ,
            c1       = self.c1    ,
            c2       = self.c2    ) 
        
        # =====================================================================
        ## Y(2S)
        # =====================================================================
        self.__Y2S   = Needham_pdf (
            name + '2S'           ,
            xvar     = self.xvar  ,
            mean     = self.m2s   ,
            sigma    = self.s2s   ,
            c0       = self.c0    ,
            c1       = self.c1    ,
            c2       = self.c2    ) 
        
        # =====================================================================
        ## Y(3S)
        # =====================================================================
        self.__Y3S   = Needham_pdf (
            name + '3S'           ,
            xvar     = self.xvar  ,
            mean     = self.m3s   ,
            sigma    = self.s3s   ,
            c0       = self.c0    ,
            c1       = self.c1    ,
            c2       = self.c2    ) 
        
        # ===============================================================================
        ## create the final PDF 
        # ===============================================================================
        self.alist1 = ROOT.RooArgList ( self.Y1S.pdf , self.Y2S.pdf , self.Y3S.pdf , self.bkg.pdf ) 
        self.alist2 = ROOT.RooArgList ( self.N1S     , self.N2S     , self.N3S     , self.B       ) 
        self.pdf    = ROOT.RooAddPdf  (
            "manca_%s"  % name ,
            "manca(%s)" % name ,
            self.alist1 ,
            self.alist2 )

        if fix_norm : self.pdf.fixCoefNormalization ( self.vars ) ## VB: added 10/10/2024 to suppress warnings 

        # ===============================================================================
        ## finally declare the components 
        # ===============================================================================
        self.signals     . add ( self.Y1S.pdf )
        self.signals     . add ( self.Y2S.pdf )
        self.signals     . add ( self.Y3S.pdf )
        self.backgrounds . add ( self.bkg.pdf ) 

        ## save configurtaion
        self.config = {
            'xvar'        : self.xvar     ,
            'name'        : self.name     ,
            'background'  : self.bkg      ,
            'm1s'         : self.m1s      ,
            's1s'         : self.s1s      ,
            ##
            'c0'          : self.c0       ,
            'c1'          : self.c1       , 
            'c2'          : self.c2       ,
            ##
            'dm21'        : self.dm21     , 
            'dm32'        : self.dm32     ,
            ## 
            'N1S'         : self.N1S      ,
            'N2S'         : self.N2S      ,
            'N3S'         : self.N3S      ,
            'B'           : self.B        ,
            'fix_norm'    : self.fix_norm }
    
    @property
    def c0 ( self ) :
        """'c0' parameter for Needham's function"""
        return self.__c0
    @c0.setter 
    def c0 (  self , value ) :
        self.set_value ( self.__c0 , value )
  
    @property
    def c1 ( self ) :
        """'c1' parameter for Needham's function"""
        return self.__c1
    @c1.setter 
    def c1 (  self , value ) :
        self.set_value ( self.__c1 , value )
          
    @property
    def c2 ( self ) :
        """'c2' parameter for Needham's function"""
        return self.__c2
    @c2.setter 
    def c2 (  self , value ) :
        self.set_value ( self.__c2 , value )
      
    @property
    def Y1S ( self ) :
        """'Y1S' : Y(1S) shape"""
        return self.__Y1S

    @property
    def Y2S ( self ) :
        """'Y2S' : Y(2S) shape"""
        return self.__Y2S

    @property
    def Y3S ( self ) :
        """'Y3S' : Y(3S) shape"""
        return self.__Y3S


models.append ( Manca_pdf ) 

# =============================================================================
## @class Manca2_pdf 
#  the final full PDF for Y->mu+mu- fit
#  This is an effective function for fit in global bin, without pt/y-binning
#  - Y signals are parameterised by the double-sided Cruystal Ball fnuctions
#
#  Parameters of three Upsilon peaks are related:
#  - Y(2S) mass is calculated as  <code>m(Y(1S))+dm21</code>
#  - Y(2S) mass is calculated as  <code>m(Y(2S))+dm32</code>
#  - Y(2S) resolution is calculated as <code>sigma(Y(1S))*(m(Y(2S))/m(Y(1S)))</code>
#  - Y(3S) resolution is calculated as <code>sigma(Y(1S))*(m(Y(3S))/m(Y(1S)))</code>
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-06-24
class Manca2_pdf (MANCA) :
    """ Manca2: the final fit model for Y->mu+mu- fit
    This is an effective function for fit in global bin, without pt/y-binning
    - three double-sided Crystal Ball functions for Y(1S), Y(2S) and Y(3S) peaks
    - constrants for their resolutions and masses 
    - background: ... 
    The Y signals are parameterised by the double-sided Crystal Ball fnuctions

    Parameters of three Upsilon peaks are related:
    - Y(2S) mass is calculated as `m(Y(1S))+dm21`
    - Y(2S) mass is calculated as  `m(Y(2S))+dm32`
    - Y(2S) resolution is calculated as `sigma(Y(1S))*(m(Y(2S))/m(Y(1S)))`
    - Y(3S) resolution is calculated as `sigma(Y(1S))*(m(Y(3S))/m(Y(1S)))`
    """
    # ===========================================================================
    ## create Manca2 PDF 
    #  @param  xvar observable, mu+mu- mass
    #  @param  name       the name of the PDF
    #  @param  background backrgudn shape: symbols, PDF or ROOT.RooAbsPdf
    #  @param  m1s        mass/location of Y(1S) peak 
    #  @param  s1s        resolution parameter for Y(1S) peak  
    #  @param alphaL      left  alpha-parameter for double-sided  Crystal Ball function
    #  @param alphaR      right alpha-parameter for double-sided  Crystal Ball function
    #  @param nL          left  n-parameter for double-sided  Crystal Ball function
    #  @param nR          right n-parameter for double-sided  Crystal Ball function
    #  @param dm21        mass difference between Y(2S) and Y(1S) states 
    #  @param dm32        mass difference between Y(3S) and Y(2S) states 
    #  @param N1S         yield for the Y(1S) signal
    #  @param N2S         yield for the Y(2S) signal
    #  @param N3S         yield for the Y(3S) signal
    #  @param B           yield for the background component    
    def __init__ ( self                    ,
                   xvar                    ,   ## the observable: mu+mu- mass 
                   name       = 'Y'        ,   ## PDF name 
                   background = 'd3'       ,   ## decreasing positive 3rd order polynomial
                   m1s        = None       ,   ## mass ioif Y(1S) state 
                   s1s        = None       ,   ## resoltuion parameter for Y(1S) state 
                   ## 
                   alphaL     = 1.65       ,   ## left  alpha-parameter for double-sided  Crystal Ball function
                   alphaR     = 1.65       ,   ## left  alpha-parameter for double-sided  Crystal Ball function
                   nL         = 1.3110     ,   ## left  n-parameter for double-sided  Crystal Ball function
                   nR         = 1.5751e+01 ,   ## right n-parameter for double-sided  Crystal Ball function
                   ## 
                   dm21       = None       ,   ## m(2S) - m(1S) mass difference 
                   dm32       = None       ,   ## m(3S) - m(2S) mass difference 
                   ## 
                   N1S        = None       ,   ## Y(2S) signal
                   N2S        = None       ,   ## Y(2S) signal
                   N3S        = None       ,   ## Y(3S) signal 
                   B          = None       ,   ## background
                   fix_norm   = False      ) : 
        
        """ Create Manca2 PDF: function to fit Y->mu+mu- signals in wide kinematic ranges 
        - xvar       : observable, mu+mu- mass
        - name       : the name of the PDF
        - background : background  shape: symbol, PDF or ROOT.RooAbsPdf        
        - m1s        : mass/location of Y(1S) peak 
        - s1s        : resolution parameter for Y(1S) peak  
        - alphaL     : left  alpha-parameter for double-sided  Crystal Ball function
        - alphaR     : right alpha-parameter for double-sided  Crystal Ball function
        - nL         : left  n-parameter for double-sided  Crystal Ball function
        - nR         : right n-parameter for double-sided  Crystal Ball function
        - dm21       : mass difference between Y(2S) and Y(1S) states 
        - dm32       :  mass difference between Y(3S) and Y(2S) states 
        - N1S        : yield for the Y(1S) signal
        - N2S        : yield for the Y(2S) signal
        - N3S        : yield for the Y(3S) signal
        - B          : yield for the background component    

        Parameters of three Upsilon peaks are related:
        - Y(2S) mass is calculated as `m(Y(1S))+dm21`
        - Y(2S) mass is calculated as  `m(Y(2S))+dm32`
        - Y(2S) resolution is calculated as `sigma(Y(1S))*(m(Y(2S))/m(Y(1S)))`
        - Y(3S) resolution is calculated as `sigma(Y(1S))*(m(Y(3S))/m(Y(1S)))`
        """
        
        ### initialize the base 
        MANCA.__init__ ( self                    ,
                         name       = name       ,
                         xvar       = xvar       ,
                         m1s        = m1s        ,
                         s1s        = s1s        ,
                         dm21       = dm21       ,
                         dm32       = dm32       ,
                         background = background , 
                         N1S        = N1S        , 
                         N2S        = N2S        , 
                         N3S        = N3S        , 
                         B          = B          )

        # =====================================================================
        ## Double Crystal Ball shape parameters 
        # =====================================================================
        self.__aL    = self.make_var ( alphaL                  ,
                                       "aL_%s"          % name ,
                                       "#alpha_{L}(%s)" % name ,
                                       None   , 1.5462     , 0.1   , 10 )
        self.__nL    = self.make_var ( nL                      ,                     
                                       "nL_%s"          % name ,
                                       "n_{L}(%s)"      % name ,
                                       None   , 1.3119     , 1.e-5 , 25 )
        self.__aR    = self.make_var ( alphaR                  ,
                                       "aR_%s"          % name ,
                                       "#alpha_{R}(%s)" % name ,
                                       None , 1.6952e+00 , 0.1     , 10 )
        self.__nR    = self.make_var ( nR                      ,
                                       "nR_%s"          % name ,
                                       "n_{R}(%s)"      % name ,
                                       None , 1.5751e+01 , 1.e-5   , 100 )

        # =====================================================================
        ## Y(1S) peak 
        # =====================================================================
        self.__Y1S   = CB2_pdf (
            name + '1S'           ,
            xvar     = self.xvar  ,
            mean     = self.m1s   ,
            sigma    = self.s1s   ,
            alphaL   = self.aL    ,
            alphaR   = self.aR    ,
            nL       = self.nL    ,
            nR       = self.nR    )

        # =====================================================================
        ## Y(2S)
        # =====================================================================
        self.__Y2S  = CB2_pdf (
            name + '2S'           ,
            xvar     = self.xvar  ,
            mean     = self.m2s   ,
            sigma    = self.s2s   ,
            alphaL   = self.aL    ,
            alphaR   = self.aR    ,
            nL       = self.nL    ,
            nR       = self.nR    )
                
        # =====================================================================
        ## Y(3S)
        # =====================================================================
        self.__Y3S  = CB2_pdf (
            name + '3S'          ,
            xvar     = self.xvar ,
            mean     = self.m3s  ,
            sigma    = self.s3s  ,
            alphaL   = self.aL   ,
            alphaR   = self.aR   ,
            nL       = self.nL   ,
            nR       = self.nR   )
        
        # ===============================================================================
        ## create the final PDF 
        # ===============================================================================
        self.alist1 = ROOT.RooArgList ( self.Y1S.pdf , self.Y2S.pdf , self.Y3S.pdf , self.bkg.pdf ) 
        self.alist2 = ROOT.RooArgList ( self.N1S     , self.N2S     , self.N3S     , self.B       ) 
        self.pdf    = ROOT.RooAddPdf  (
            "manca_%s"  % name ,
            "manca(%s)" % name ,
            self.alist1 ,
            self.alist2 )

        if fix_norm : self.pdf.fixCoefNormalization ( self.vars ) ## VB: added 10/10/2024 to suppress warnings 

        # ===============================================================================
        ## finally declare the components 
        # ===============================================================================
        self.signals     . add ( self.Y1S.pdf )
        self.signals     . add ( self.Y2S.pdf )
        self.signals     . add ( self.Y3S.pdf )
        self.backgrounds . add ( self.bkg.pdf ) 

        ## save configurtaion
        self.config = {
            'xvar'        : self.xvar     ,
            'name'        : self.name     ,
            'background'  : self.bkg      ,
            'm1s'         : self.m1s      ,
            's1s'         : self.s1s      ,
            ##
            'alphaL'      : self.aL       , 
            'alphaR'      : self.aR       , 
            'nL'          : self.nL       ,
            'nR'          : self.nR       ,
            ##
            'dm21'        : self.dm21     , 
            'dm32'        : self.dm32     ,
            ## 
            'N1S'         : self.N1S      ,
            'N2S'         : self.N2S      ,
            'N3S'         : self.N3S      ,
            'B'           : self.B        ,
            'fix_norm'    : self.fix_norm } 
        
    @property
    def aL (  self ) :
        """(left)'alpha'-parameter for Y-peaks"""
        return self.__aL
    @aL.setter 
    def aL (  self , value ) :
        self.set_value ( self.__aL , value )
        
    @property
    def aR (  self ) :
        """(right)'alpha'-parameter for Y-peaks"""
        return self.__aR
    @aR.setter 
    def aR (  self , value ) :
        self.set_value ( self.__aR , value )
        
    @property
    def nL (  self ) :
        """(left)'n'-parameter for Y-peaks"""
        return self.__nL
    @nL.setter 
    def nL (  self , value ) :
        self.set_value ( self.__nL , value )

    @property
    def nR (  self ) :
        """(right)'n'-parameter for Y-peaks"""
        return self.__nR
    @nR.setter 
    def nR (  self , value ) :
        self.set_value ( self.__nR , value )

    @property
    def Y1S ( self ) :
        """'Y1S' : Y(1S) shape"""
        return self.__Y1S

    @property
    def Y2S ( self ) :
        """'Y2S' : Y(2S) shape"""
        return self.__Y2S

    @property
    def Y3S ( self ) :
        """'Y3S' : Y(3S) shape"""
        return self.__Y3S

models.append ( Manca2_pdf )

# =============================================================================
## Specific model for fitting of Y+X
class MancaX_pdf(PDF2) :
    """ MancaX: 2D-model to study associative production of Upsilon and X 
    """
    def __init__ ( self             ,
                   manca            , ## manca pdf, that defined 3 upsilon peaks  
                   charm            , ## charm pdf
                   bkg1     = 0     ,
                   bkg2     = 0     ,
                   bkgA     = 0     ,
                   bkgB     = 0     ,
                   suffix   = ''    ,
                   fix_norm = False ) :

        PDF2.__init__ ( self , 'MancaX' + suffix , manca.xvar , charm.xvar )  
        self._crossterms1 = ROOT.RooArgSet ()
        self._crossterms2 = ROOT.RooArgSet ()

        self.__suffix = suffix 
        self.__manca  = manca
        self.__charm  = charm
                                 
        self.suffix  = suffix
        self.signal1 = manca  
        self.signal2 = charm

        #
        ## background components
        #
        self.b_Y = self.make_bkg ( bkg1 , 'BkgY' + suffix , self.m1 )   
        self.b_C = self.make_bkg ( bkg2 , 'BkgC' + suffix , self.m2 )   
        self.b_A = self.make_bkg ( bkgA , 'BkgA' + suffix , self.m1 )   
        self.b_B = self.make_bkg ( bkgB , 'BkgB' + suffix , self.m2 )
        
        #
        ## pure signal components: 3 
        #
        RPP  = ROOT.RooProdPdf 
        
        self.y1s_c   = RPP     ( 'Y1SC' + suffix , 'Y(1S)(x)Charm' , manca.y1s , charm.pdf )
        self.y2s_c   = RPP     ( 'Y2SC' + suffix , 'Y(2S)(x)Charm' , manca.y2s , charm.pdf )
        self.y3s_c   = RPP     ( 'Y3SC' + suffix , 'Y(3S)(x)Charm' , manca.y3s , charm.pdf )
        
        ## charm + background
        self.bs_pdf  = RPP     ( 'BC'   + suffix , 'B(2mu)(x)Charm' , self.b_Y.pdf , charm.pdf ) 
        
        ## Y     + background
        self.y1s_b   = RPP     ( 'Y1SB' + suffix , 'Y(1S)(x)Bkg'    , manca.y1s , self.b_C.pdf ) 
        self.y2s_b   = RPP     ( 'Y2SB' + suffix , 'Y(2S)(x)Bkg'    , manca.y2s , self.b_C.pdf ) 
        self.y3s_b   = RPP     ( 'Y3SB' + suffix , 'Y(3S)(x)Bkg'    , manca.y3s , self.b_C.pdf ) 

        ## background + background
        self.bb_pdf  = RPP ( 'BB'   + suffix , 'Bkg(x)Bkg'     , self.b_A.pdf , self.b_B.pdf ) 
        
        ## coefficients
        self.s1s = self.make_var ( None , 'NY1C' + suffix , 'N(Y1S+C)' + suffix , None , 200 , 0 , 1.e+5 )
        self.s2s = self.make_var ( None , 'NY2C' + suffix , 'N(Y2S+C)' + suffix , None , 100 , 0 , 1.e+4 )
        self.s3s = self.make_var ( None , 'NY3C' + suffix , 'N(Y3S+C)' + suffix , None ,  20 , 0 , 1.e+4 )
        self.bs  = self.make_var ( None , 'NBC'  + suffix , 'N(Bkg+C)' + suffix , None , 500 , 0 , 1.e+5 )
        self.s1b = self.make_var ( None , 'NY1B' + suffix , 'N(Y1S+B)' + suffix , None , 500 , 0 , 1.e+5 )
        self.s2b = self.make_var ( None , 'NY2B' + suffix , 'N(Y2S+B)' + suffix , None , 200 , 0 , 1.e+4 )
        self.s3b = self.make_var ( None , 'NY3B' + suffix , 'N(Y3S+B)' + suffix , None , 200 , 0 , 1.e+4 )
        self.bb  = self.make_var ( None , 'NBB'  + suffix , 'N(Bkg+B)' + suffix , None , 500 , 0 , 1.e+5 )

        self.S1S_name = self.s1s.GetName()
        self.S2S_name = self.s2s.GetName()
        self.S3S_name = self.s3s.GetName()
        self.S1B_name = self.s1b.GetName()
        self.S2B_name = self.s2b.GetName()
        self.S3B_name = self.s3b.GetName()
        self.BS_name  = self.bs .GetName()
        self.BB_name  = self.bb .GetName()
        
        self.alist1 = ROOT.RooArgList ( self.y1s_c  , self.y2s_c , self.y3s_c ,
                                        self.bs_pdf ,
                                        self.y1s_b  , self.y2s_b , self.y3s_b ,
                                        self.bb_pdf )
        self.alist2 = ROOT.RooArgList ( self.s1s    , self.s2s   , self.s3s   ,
                                        self.bs     ,
                                        self.s1b    , self.s2b   , self.s3b   ,
                                        self.bb     )
        #
        ## build final PDF 
        # 
        self.pdf  = ROOT.RooAddPdf  ( "model2D"      + suffix ,
                                      "Model2D(%s)"  % suffix ,
                                      self.alist1 ,
                                      self.alist2 )
        
        if fix_norm : self.pdf.fixCoefNormalization ( self.vars ) ## VB: added 10/10/2024 to suppress warnings 

        self.name    = self.pdf.GetName()
        ##
        self.signals     () . add ( self.y1s_c  )
        self.signals     () . add ( self.y2s_c  )
        self.signals     () . add ( self.y3s_c  )
        ##
        self.backgrounds () . add ( self.bb_pdf )
        ##
        self.crossterms1 () . add ( self.y1s_b  )
        self.crossterms1 () . add ( self.y2s_b  )
        self.crossterms1 () . add ( self.y3s_b  )
        ##
        self.crossterms2 () . add ( self.bs_pdf )                              
    
        ## save configurtaion
        self.config = {
            'mass'     : self.mass     ,
            'manca'    : self.name     ,
            'charm'    : self.charma   ,
            'bkg1'     : self.b_Y      , 
            'bkg2'     : self.b_C      , 
            'bkgA'     : self.b_A      , 
            'bkgB'     : self.b_B      , 
            'suffix'   : self.__siffix ,
            'fix_norm' : self.fix_norm  
            }

    ## get all declared components 
    def crossterms1 ( self ) : return self._crossterms1
    ## get all declared components 
    def crossterms2 ( self ) : return self._crossterms2

    @property
    def manca (  self ) :
        """'manca'-function to decribe Upsilon peaks"""
        return self.__manca
    @property
    def charm (  self ) :
        """'charm'-function to decribe Charm peak"""
        return self.__charm
    
    @property
    def fix_norm ( self ) :
        """`fix-norm`: 
        - see `ROOT.RooAbsPdf.SetCoefNormalization`
        - see `ROOT.RooAbsPdf.getCoefNormalization`
        """
        pars = self.pdf.getCoefNormalization()
        return True if pars else False 

models.append ( MancaX_pdf ) 

# =============================================================================
## @class HORNSdini_pdf
#  \f[ f(x;a,\delta, \phi) = 
#  \frac{3}{2\delta}\left( z \right)^2
#  \left( \cos^2( \phi + \frac{\pi}{4}) ( 1 + z ) +
#         \sin^2( \phi + \frac{\pi}{4}) ( 1 - z ) \right) \f]
#  where  \f$ z = \frac{ x - ( a - \delta ) } { \delta } \f$ 
#  for \f$ a \le x \le a + 2\delta\$ and zero otherwise 
#  
# The first factor accound for two-horn parabolic shape, 
# and the second factor accouns for the linear correction factor 
# ("efficiency")
#
#  @attention: For the actual use it needs to be convoluted with resolution function 
#  @see Ostap::Math::HORNSdini
#  @see Ostap::Modls::HORNSdini
class HORNSdini_pdf(PEAK) :
    """ HORNSdini PDF
    """
    def __init__ ( self              ,
                   name              ,
                   xvar              ,   ## mass is mandatory here! 
                   a                 , 
                   delta             ,
                   phi        = None ,
                   resolution = None ,
                   cnvpars    = {}   ) :


        PEAK.__init__ ( self           ,
                        name        = name   ,
                        xvar        = xvar   ,
                        mean        = a      ,
                        sigma       = delta  ,
                        mean_name   = 'a_%s'               % name ,
                        mean_title  = 'a_{HORNS}(%s)'      % name ,
                        sigma_name  = 'delta_%s'           % name ,
                        sigma_title = '#delta_{HORNS}(%s)' % name )
                        
        self.__delta = self.sigma
        self.__a     = self.mean
        self.__phi   = self.make_var ( phi                           ,
                                       'phi_%s'          % self.name ,
                                       '#phi_{HORN}(%s)' % self.name ,
                                       None , 0 , -12 , 12 )
        
        self.__horns = Ostap.Models.HORNSdini (
            self.new_roo_name ( 'horns_' ) ,
            "HORNSdini %s" % self.name     ,
            self.xvar                      ,
            self.a                         ,
            self.delta                     ,
            self.phi                       )


        self.__resolution = None
        self.__cnvpars    = {}
        self.__cnvpars.update ( cnvpars )
        
        if resolution :

            self.__resolution = resolution 
            cname = self.generate_name ( prefix = self.name , suffix = 'cnv' )  
            from ostap.fitting.convolution import Convolution_pdf 
            self.__convolved = Convolution_pdf ( name       = cname       , 
                                                 pdf        = self.horns  ,
                                                 xvar       = self.xvar   ,
                                                 resolution = resolution  ,
                                                 **self.cnvpars           )
            
            ## save the resolution 
            self.__resolution = self.__convolved.resolution
            
            ## finally get the convolved PDF 
            self.pdf = self.__convolved.pdf
            
        else :
            
            self.pdf = self.horns
            
        
        ## save configuration
        self.config = {
            'xvar'       : self.xvar       ,
            'name'       : self.name       ,
            'a'          : self.a          ,
            'delta'      : self.delta      ,
            'phi'        : self.phi        ,
            'resolution' : self.resolution ,
            'cnvpars'    : self.cnvpars    }

    @property
    def a ( self ) :
        """'a' : position of the left horn  (same as 'mean')"""
        return self.__a
    @a.setter
    def a ( self ) :
        self.set_value ( self.__a , value )

    @property
    def delta ( self ) :
        """'delta' : half-distance from the left to right horns (same as 'sigma')"""
        return self.__delta
    @delta.setter
    def delta ( self ) :
        self.set_value ( self.__delta , value )

    @property
    def phi ( self ) :
        """'phi' : correction/modificaiton parameter"""
        return self.__phi
    @phi.setter
    def phi ( self , value ) :
        self.set_value ( self.__phi , value )

    @property
    def horns ( self )  :
        """'horns' : get the HORNSdini PDF"""
        return self.__horns

    @property
    def resolution ( self )  :
        """'resolution' : the resolution function"""
        return self.__resolution 

    @property
    def cnvpars ( self )  :
        """'cnvpars' : parameters used for convolution"""
        return self.__cnvpars 

models.append ( HORNSdini_pdf ) 



# =============================================================================
## @class HILLdini_pdf
#  \f[ f(x;a,\delta, \phi) = 
#  \frac{3}{4\delta}\left( 1 - z^2 \right)
#  \left( \cos^2( \phi + \frac{\pi}{4}) ( 1 + z ) +
#         \sin^2( \phi + \frac{\pi}{4}) ( 1 - z ) \right) \f]
#  where  \f$ z = \frac{ x - ( a - \delta ) } { \delta } \f$ 
#  for \f$ a \le x \le a + 2\delta\$ and zero otherwise 
#  
# The first factor accound for two-horn parabolic shape, 
# and the second factor accouns for the linear correction factor 
# ("efficiency")
#
#  @attention: For the actual use it needs to be convoluted with resolution function 
#  @see Ostap::Math::HILLdini
#  @see Ostap::Modls::HILLdini
class HILLdini_pdf(PEAK) :
    """HORNSdini PDF
    """
    def __init__ ( self              ,
                   name              ,
                   xvar              ,   ## mass is mandatory here! 
                   a                 , 
                   delta             ,
                   phi        = None ,
                   resolution = None ,
                   cnvpars    = {}   ) :


        PEAK.__init__ ( self           ,
                        name        = name   ,
                        xvar        = xvar   ,
                        mean        = a      ,
                        sigma       = delta  ,
                        mean_name   = 'a_%s'               % name ,
                        mean_title  = 'a_{HILL}(%s)'       % name ,
                        sigma_name  = 'delta_%s'           % name ,
                        sigma_title = '#delta_{HILL}(%s)'  % name )
                        
        self.__delta = self.sigma
        self.__a     = self.mean
        self.__phi   = self.make_var ( phi                           ,
                                       'phi_%s'          % self.name ,
                                       '#phi_{HILL}(%s)' % self.name ,
                                       None , 0 , -12 , 12 )
        
        self.__hill = Ostap.Models.HILLdini (
            self.new_roo_name ( 'hill_' ) ,
            "HILLdini %s" % self.name     ,
            self.xvar                     ,
            self.a                        ,
            self.delta                    ,
            self.phi                      )


        self.__resolution = None
        self.__cnvpars    = {}
        self.__cnvpars.update ( cnvpars )
        
        if resolution :

            self.__resolution = resolution 
            cname = self.generate_name ( prefix = self.name , suffix = 'cnv' )  
            from ostap.fitting.convolution import Convolution_pdf 
            self.__convolved = Convolution_pdf ( name       = cname       , 
                                                 pdf        = self.hill   ,
                                                 xvar       = self.xvar   ,
                                                 resolution = resolution  ,
                                                 **self.cnvpars           )
            
            ## save the resolution 
            self.__resolution = self.__convolved.resolution
            
            ## finally get the convolved PDF 
            self.pdf = self.__convolved.pdf
            
        else :
            
            self.pdf = self.hill
            
        
        ## save configuration
        self.config = {
            'xvar'       : self.xvar       ,
            'name'       : self.name       ,
            'a'          : self.a          ,
            'delta'      : self.delta      ,
            'phi'        : self.phi        ,
            'resolution' : self.resolution ,
            'cnvpars'    : self.cnvpars    }

    @property
    def a ( self ) :
        """'a' : position of the left horn  (same as 'mean')"""
        return self.__a
    @a.setter
    def a ( self ) :
        self.set_value ( self.__a , value )

    @property
    def delta ( self ) :
        """'delta' : half-distance from the left to right horns (same as 'sigma')"""
        return self.__delta
    @delta.setter
    def delta ( self ) :
        self.set_value ( self.__delta , value )

    @property
    def phi ( self ) :
        """'phi' : correction/modificaiton parameter"""
        return self.__phi
    @phi.setter
    def phi ( self , value ) :
        self.set_value ( self.__phi , value )

    @property
    def hill ( self )  :
        """'hill' : get the HILLdini PDF"""
        return self.__hill

    @property
    def resolution ( self )  :
        """'resolution' : the resolution function"""
        return self.__resolution 

    @property
    def cnvpars ( self )  :
        """'cnvpars' : parameters used for convolution"""
        return self.__cnvpars 

models.append ( HORNSdini_pdf ) 




# =============================================================================
## @class HHdini_pdf
#  fL*HORNS + ( 1-fl) * HILL
#  @see Ostap::Math::HORNSdini
#  @see Ostap::Math::HILLdini
#  @see Ostap::Modls::HORNSdini
#  @see Ostap::Modls::HILLdini
class HHdini_pdf(PEAK) :
    """HHdini PDF: fl*HORNS+(1-f)*HILL
    """
    def __init__ ( self              ,
                   name              ,
                   xvar              ,   ## mass is mandatory here! 
                   a                 , 
                   delta             ,
                   phi        = None ,
                   fL         = None ,
                   resolution = None ,
                   cnvpars    = {}   ) :
        
        PEAK.__init__ ( self           ,
                        name        = name   ,
                        xvar        = xvar   ,
                        mean        = a      ,
                        sigma       = delta  ,
                        mean_name   = 'a_%s'             % name ,
                        mean_title  = 'a_{HH}(%s)'       % name ,
                        sigma_name  = 'delta_%s'         % name ,
                        sigma_title = '#delta_{HH}(%s)'  % name )
                        
        self.__delta = self.sigma
        self.__a     = self.mean
        self.__phi   = self.make_var ( phi                           ,
                                       'phi_%s'          % self.name ,
                                       '#phi_{HH}(%s)'   % self.name ,
                                       None , 0 , -12 , 12 )

        self.__fL = self.make_var    ( fL ,
                                       'fL_%s'           % self.name ,
                                       'f_{L}(%s)'       % self.name ,
                                       None , 0.5 , 0 , 1 )
        
        
        self.__horns = HORNSdini_pdf ( self.new_name ( 'HORNS' )          ,
                                       xvar       = self.xvar             ,
                                       a          = self.a                ,
                                       delta      = self.delta            ,
                                       phi        = self.phi              ,
                                       resolution = resolution            ,
                                       cnvpars    = cnvpars               )
        
        self.__hill  = HILLdini_pdf  ( self.new_name ( 'HILL'  )          ,
                                       xvar       = self.horns.xvar       ,
                                       a          = self.horns.a          ,
                                       delta      = self.horns.delta      ,
                                       phi        = self.horns.phi        ,
                                       resolution = self.horns.resolution ,
                                       cnvpars    = self.horns.cnvpars    )
        
        self.__sum = Sum1D ( ( self.horns , self.hill ) ,
                             recursive = True           ,
                             xvar      =   self.xvar    ,
                             fractions = ( self.fL , )  )
        
        ## finally get the PDF
        self.pdf = self.__sum.pdf

        
        ## save configuration
        self.config = {
            'xvar'       : self.xvar       ,
            'name'       : self.name       ,
            'a'          : self.a          ,
            'delta'      : self.delta      ,
            'fL'         : self.fL         ,
            'phi'        : self.phi        ,
            'resolution' : self.resolution ,
            'cnvpars'    : self.cnvpars    }

    @property
    def a ( self ) :
        """'a' : position of the left horn (same as 'mean')"""
        return self.__a
    @a.setter
    def a ( self ) :
        self.set_value ( self.__a , value )

    @property
    def delta ( self ) :
        """'delta' : half-distance from the left to right horns (same as 'sigma')"""
        return self.__delta
    @delta.setter
    def delta ( self ) :
        self.set_value ( self.__delta , value )

    @property
    def fL ( self ) :
        """'fL' : fraction of HORNS component"""
        return self.__fL 
    @fL.setter
    def fL ( self , value ) :
        self.set_value ( self.__fL , value )

    @property
    def phi ( self ) :
        """'phi' : correction/modificaiton parameter"""
        return self.__phi
    @phi.setter
    def phi ( self , value ) :
        self.set_value ( self.__phi , value )

    @property
    def horns ( self )  :
        """'horns' : get the HORNSdini PDF"""
        return self.__horns

    @property
    def hill ( self )  :
        """'hill' : get the HILLdini PDF"""
        return self.__hill

    @property
    def resolution ( self )  :
        """'resolution' : the resolution function"""
        return self.horns.resolution 

    @property
    def cnvpars ( self )  :
        """'cnvpars' : parameters used for convolution"""
        return self.horns.cnvpars 

models.append ( HHdini_pdf ) 


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )

 
# =============================================================================
##                                                                      The END 
# =============================================================================
