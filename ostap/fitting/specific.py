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
    'D0_pdf'  , ## PDF for D0        : Bukin 
    'Dp_pdf'  , ## PDF for D+        : Bukin
    'Ds_pdf'  , ## PDF for Ds+       : Bukin 
    'Lc_pdf'  , ## PDF for Lambda_c+ : Gauss
    #
    'Bd_pdf'  , ## pdf for B0        : double-sided Crystal Ball 
    'B0_pdf'  , ## pdf for B0        : double-sided Crystal Ball 
    'Bu_pdf'  , ## pdf for B+        : double-sided Crystal Ball 
    'Bs_pdf'  , ## pdf for Bs        : double-sided Crystal Ball 
    'Bc_pdf'  , ## pdf for Bc+       : double-sided Crystal Ball 
    #
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
from   ostap.fitting.pdfbasic    import PDF1, PDF2 
from   ostap.fitting.fit1d       import PEAK, Sum1D 
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
    """B0: double-sided Crystal Ball function
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
    """B+: double-sided Crystal Ball function
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

models.append ( Bu_pdf ) 
# =============================================================================
## @class Bs_pdf
#  simple wrapper over CB2-pdf
#  @see Ostap::Models::CrystalBallDS
#  @see Ostap::Math::CrystalBallDS
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Bs_pdf(CB2_pdf) :
    """Bs: double-sided Crystal Ball function
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
    """Bc: double-sided Crystal Ball function
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


models.append ( Bc_pdf ) 
# =============================================================================
# Specializations for Bukin function
# =============================================================================
## @class D0_pdf
#  simple wrapper over Bukin-pdf
#  @see RooBukinPdf
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
        ## save configuration
        self.config = {
            'xvar'   : self.xvar  ,
            'name'   : self.name  ,
            'mean'   : self.mean  ,
            'sigma'  : self.sigma ,
            'xi'     : self.xi    ,
            'rhoL'   : self.rhoL  ,
            'rhoR'   : self.rhoR  ,
            }

                             
models.append ( D0_pdf ) 
# =============================================================================
## @class Dp_pdf
#  simple wrapper over Bukin-pdf
#  @see RooBukinPdf
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Dp_pdf(Bukin_pdf) :
    """D+: Bukin function 
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
        ## save configuration
        self.config = {
            'xvar'   : self.xvar  ,
            'name'   : self.name  ,
            'mean'   : self.mean  ,
            'sigma'  : self.sigma , 
            'xi'     : self.xi    ,
            'rhoL'   : self.rhoL  ,
            'rhoR'   : self.rhoR  ,
            }
        
models.append ( Dp_pdf ) 
# =============================================================================
## @class Ds_pdf
#  simple wrapper over Bukin-pdf
#  @see RooBukinPdf
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Ds_pdf(Bukin_pdf) :
    """Ds: Bukin function 
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
        ## save configuration
        self.config = {
            'xvar'   : self.xvar  ,
            'name'   : self.name  ,
            'mean'   : self.mean  ,
            'sigma'  : self.sigma ,
            'xi'     : self.xi    ,
            'rhoL'   : self.rhoL  ,
            'rhoR'   : self.rhoR  ,
            }

        
models.append ( Ds_pdf ) 
# =============================================================================
## @class Lc_pdf
#  simple wrapper over Bukin-pdf
#  @see RooBukinPdf
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Lc_pdf(Bukin_pdf) :
    """Lc: Bukin function 
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
        ## save configuration
        self.config = {
            'xvar'   : self.xvar  ,
            'name'   : self.name  ,
            'mean'   : self.mean  ,
            'sigma'  : self.sigma ,
            'xi'     : self.xi    ,
            'rhoL'   : self.rhoL  ,
            'rhoR'   : self.rhoR  ,
            }

        
models.append ( Lc_pdf ) 
# =============================================================================
## @class Manca_pdf 
#  the final full PDF for Y->mu+mu- fit
#  This is physically well-motivated function for fits in narrow
#  bins in pt and rapidity  
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-08-02
class Manca_pdf (PDF1) :
    """Manca: the final fit model for Y->mu+mu- fit
    This is physically well-motivated function for fits in narrow bins in pt and rapidity
    - three Needham functions for Y(1S), Y(2S) and Y(3S) peaks
    - constrants for their resoltuions and masses 
    - background: exponent modulated by positive polynomial 
    
    """
    def __init__ ( self          ,
                   mass          ,
                   name   = 'Y'  ,
                   power  = 0    ,
                   m1s    = None ,
                   sigma  = None ,
                   a0     = 1.91 ,
                   a1     = None ,
                   a2     = None ) : 

        #
        PDF1.__init__ ( self , name , mass )
        #
        if   9460.  in self.mass and 10023. in self.mass and 10355. in self.mass :  gev_ = 1000        
        elif 9.460  in self.mass and 10.023 in self.mass and 10.355 in self.mass : gev_ = 1
        else :
            raise AttributeError ( "Illegal mass range %s<m<%s" % self.xminmax()  ) 
        
        m_y1s  =  9.46030     * gev_ 
        s_y1s  =  4.3679e-02  * gev_ 
        dm_y2s = 10.02326     * gev_ - m_y1s
        dm_y3s = 10.3552      * gev_ - m_y1s

        # =====================================================================
        ## Y(1S)
        # =====================================================================
        
        self.__a0   = self.make_var ( a0                 ,
                                      'a0m_%s' % name    ,
                                      "a0 for Needham's function" ,
                                      True  , 1.91  , 0.1 , 3.0   )
        
        self.__a1   = self.make_var ( a1                 ,
                                      'a1m_%s' % name    ,
                                      "a1 for Needham's function" ,
                                      True , 1.1174 / gev_ ,  -10.0 / gev_ , 10.0 /gev_ )
        
        self.__a2   = self.make_var ( a2                 ,
                                      'a2m_%s' % name    ,
                                      "a2 for Needham's function" , 
                                      True , -5.299 / gev_**2   , -100.0  / gev_**2  , 100.0  /gev_**2 )
        
        ## default logic does not work nicely here, therefore we need to be explicit:
        if a0 is None and not self.a0.isConstant () : self.a0.fix (  1.91             ) 
        if a1 is None and not self.a1.isConstant () : self.a1.fix (  1.1174 / gev_    )
        if a2 is None and not self.a2.isConstant () : self.a2.fix ( -5.299  / gev_**2 ) 
        
        self.__m1s   = self.make_var ( m1s                   ,
                                       "m1S_%s"       % name ,
                                       "mass Y1S(%s)" % name , 
                                       None , m_y1s , m_y1s - 0.15 * s_y1s , m_y1s + 0.15 * s_y1s ) 
        
        self.__s1s   = self.make_var ( sigma                  ,
                                       "s1S_%s"        % name ,
                                       "sigma Y1S(%s)" % name ,
                                       None , s_y1s , 0.3 * s_y1s , 4 * s_y1s )
        
        self.__sigma = self.s1s 
        self.__Y1S   = Needham_pdf (
            name + '1S'           ,
            xvar     = self.mass  ,
            mean     = self.m1s   ,
            sigma    = self.s1s   ,
            a0       = self.a0    ,
            a1       = self.a1    ,
            a2       = self.a2    ) 
        
        # =====================================================================
        ## Y(2S)
        # =====================================================================
        self.__dm2s  = self.make_var ( None ,
                                       "dm2s"      + name    ,
                                       "dm2s(%s)"  % name    ,
                                       True ,
                                       dm_y2s ,
                                       dm_y2s - 0.20 * s_y1s , 
                                       dm_y2s + 0.20 * s_y1s )
        
        self.__aset11 = ROOT.RooArgList ( self.__m1s , self.__dm2s )

        ## self.__m2s    = ROOT.RooFormulaVar (
        ## "m_" + name + '2S'   ,
        ## "m2s(%s)"  % name    ,
        ## "%s+%s" % ( self.__m1s.GetName() , self.__dm2s.GetName()  ) , 
        ## self.__aset11       )
        self.__m2s    = Ostap.MoreRooFit.Addition ( "m_" + name + '2S' ,
                                                    "m2s(%s)"  % name  ,
                                                    self.__m1s , self.__dm2s ) 
        
        self.__aset12 = ROOT.RooArgList ( self.__sigma , self.__m1s , self.__m2s ) 
        self.__s2s    = Ostap.FormulaVar (
            "sigma_"  + name + '2S'    ,
            "#sigma_{Y2S}(%s)" % name  ,
            "%s*(%s/%s)"  % ( self.__sigma.GetName() ,
                              self.__m2s  .GetName() ,
                              self.__m1s  .GetName() ) ,
            self.__aset12  )
        
        self.__Y2S   = Needham_pdf (
            name + '2S'           ,
            xvar     = self.mass  ,
            mean     = self.m2s   ,
            sigma    = self.s2s   ,
            a0       = self.a0    ,
            a1       = self.a1    ,
            a2       = self.a2    ) 
        
        # =====================================================================
        ## Y(3S)
        # =====================================================================
        self.__dm3s  = self.make_var ( None                  ,
                                       "dm3s"      + name    ,
                                       "dm3s(%s)"  % name    ,
                                       True                  , 
                                       dm_y3s                ,
                                       dm_y3s - 0.20 * s_y1s , 
                                       dm_y3s + 0.20 * s_y1s )
        
        self.__aset21 = ROOT.RooArgList ( self.__m1s , self.__dm3s )
        ## self.__m3s    = ROOT.RooFormulaVar (
        ## "m_"       + name + '(3S)' ,
        ## "m3s(%s)"  % name          ,
        ## "%s+%s" % ( self.__m1s.GetName() , self.__dm3s.GetName() ) ,
        ## self.__aset21       )
        self.__m3s    = Ostap.MoreRooFit.Addition ( "m_" + name + '3S' ,
                                                    "m3s(%s)"  % name  ,
                                                    self.__m1s , self.__dm3s ) 
        
        
        self.__aset22 = ROOT.RooArgList ( self.__sigma , self.__m1s , self.__m3s ) 
        self.__s3s    = Ostap.FormulaVar (
            "sigma_"  + name + '3S'    ,
            "#sigma_{Y3S}(%s)" % name  ,
            "%s*(%s/%s)"  % ( self.__sigma.GetName() ,
                              self.__m3s  .GetName() ,
                              self.__m1s  .GetName() ) , 
            self.__aset22       )
        
        self.Y3S   = Needham_pdf (
            name + '3S'           ,
            xvar     = self.mass  ,
            mean     = self.m3s   ,
            sigma    = self.s3s   ,
            a0       = self.a0    ,
            a1       = self.a1    ,
            a2       = self.a2    ) 
        
        #
        ## the actual signal PDFs
        # 
        self.__y1s   = self.Y1S.pdf
        self.__y2s   = self.Y2S.pdf
        self.__y3s   = self.Y3S.pdf

        self.__power = power 
        ## use helper function to create background  
        self.background = self.make_bkg ( power , 'Bkg%s' % name , self.mass )

        self.__n1s = self.make_var ( None ,
                                     "N1S" + name  ,
                                     "Signal(Y1S)" ,  None , 1000 ,  0 ,  1.e+8 )
        self.__n2s = self.make_var ( None ,
                                     "N2S" + name  ,
                                     "Signal(Y2S)" ,  None ,  300 ,  0 ,  1.e+7 )
        self.__n3s = self.make_var ( None ,
                                     "N3S" + name  ,
                                     "Signal(Y3S)" ,  None ,  100 ,  0 ,  1.e+6 )
        self.__b   = self.make_var ( None ,
                                     "B"   + name  ,
                                     "Background"  ,  None ,  100 ,  0 ,  1.e+8 )
        
        self.alist1 = ROOT.RooArgList ( self.__y1s , self.__y2s , self.__y3s ) 
        self.alist2 = ROOT.RooArgList ( self.__n1s , self.__n2s , self.__n3s ) 
        
        self.alist1 . add ( self.background.pdf )
        self.alist2 . add ( self.__b            )
        
        self.pdf  = ROOT.RooAddPdf  (
            "manca_%s"  % name ,
            "manca(%s)" % name ,
            self.alist1 ,
            self.alist2 )
        
        self.__dm2s.setConstant ( True )
        self.__dm3s.setConstant ( True )
        
        self._splots = []
        
        self.s1_name = self.__n1s.GetName ()
        self.s2_name = self.__n2s.GetName ()
        self.s3_name = self.__n3s.GetName ()

        # 
        ## finally declare components 
        #
        self.signals    () . add ( self.__y1s )
        self.signals    () . add ( self.__y2s )
        self.signals    () . add ( self.__y3s )
        self.backgrounds() . add ( self.background.pdf ) 

        ## save configurtaion
        self.config = {
            'mass'  : self.mass  ,
            'name'  : self.name  ,
            'power' : self.power ,
            'm1s'   : self.m1s   ,
            'sigma' : self.sigma ,
            'a0'    : self.a0    ,
            'a1'    : self.a1    ,
            'a2'    : self.a2    ,
            }
    
    def alpha_1S ( self ) : return self.Y1S.pdf.alpha ()
    def alpha_2S ( self ) : return self.Y2S.pdf.alpha ()
    def alpha_3S ( self ) : return self.Y3S.pdf.alpha ()

    @property
    def mass ( self ) :
        """'mass'-variable, ailas for 'x' and 'xvar'"""
        return self.xvar
    
    @property
    def power ( self ) :
        """'power'-variable for background"""
        return self.__power

    @property
    def  s1s (  self ) :
        """'s1s'-parameter (resoltuion for Y(1S) peak) (same as 'sigma')"""
        return self.__sigma
    @s1s.setter
    def  s1s (  self , value ) :
        self.set_value ( self.__s1s , value )

    @property
    def  sigma (  self ) :
        """'sigma'-parameter (resolution for Y(1S) peak) (same as 's1s')"""
        return self.s1s
    @sigma.setter
    def  sigma (  self , value ) :
        self.s1s = value

    @property
    def a0 ( self ) :
        """'a0' parameter for Needham function"""
        return self._a0
    @property
    def a1 ( self ) :
        """'a1' parameter for Needham function"""
        return self._a1
    @property
    def a2 ( self ) :
        """'a2' parameter for Needham function"""
        return self._a0
    
    @property
    def  m1s (  self ) :
        """'m1s'-parameter (mass for Y(1S) peak)"""
        return self.__m1s
    @m1s.setter
    def  m1s (  self , value ) :
        self.set_value ( self.__m1s , value )

    @property
    def  m2s (  self ) :
        """'m2s'-parameter (mass for Y(2S) peak)"""
        return self.__m2s
    @property
    def  m3s (  self ) :
        """'m3s'-parameter (mass for Y(3S) peak)"""
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
    def  dm2s (  self ) :
        """'dm2s'-parameter (mass  difference for Y(2S) and Y(1S) peaks)"""
        return self.__dm2s    
    @dm2s.setter
    def  dm2s (  self , value ) :
        self.set_value ( self.__dm2s , value )

    @property
    def  dm3s (  self ) :
        """'dm3s'-parameter (mass  difference for Y(3S) and Y(1S) peaks)"""
        return self.__dm3s    
    @dm3s.setter
    def  dm3s (  self , value ) :
        self.set_avlue ( self.__dm3s , value )
        
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
        self.set_avlue ( self.__n2s , value )

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
        
    
models.append ( Manca_pdf ) 
# =============================================================================
## @class Manca2_pdf 
#  the final full PDF for Y->mu+mu- fit
#  This is an effective function for fit in global bin, without pt/y-binning 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-06-24
class Manca2_pdf (PDF1) :
    """Manca2: the final fit model for Y->mu+mu- fit
    This is an effective function for fit in global bin, without pt/y-binning
    - three double-sided Crystal Ball functions for Y(1S), Y(2S) and Y(3S) peaks
    - constrants for their resoltuions and masses 
    - background: exponent modulated by positive polynomial 
    """
    def __init__ ( self            ,
                   mass            ,
                   name   = 'Y'    ,
                   power  = 0      ,
                   m1s    = None   ,
                   sigma  = None   ,
                   alphaL = 1.5462 ,
                   alphaR = 1.6952 ,
                   nL     = 1.3110 ,
                   nR     = 1.5751e+01 ) :

        #
        PDF1.__init__ ( self , name , mass )
        #
        if   9460.  in self.mass and 10023. in self.mass and 10355. in self.mass :  gev_ = 1000        
        elif 9.460  in self.mass and 10.023 in self.mass and 10.355 in self.mass : gev_ = 1
        else :
            raise AttributeError ( "Illegal mass range %s<m<%s" % self.xminmax()  ) 

        m_y1s  =  9.46030     * gev_
        s_y1s  =  4.03195e-02 * gev_ 
        dm_y2s = 10.02326     * gev_ - m_y1s
        dm_y3s = 10.3552      * gev_ - m_y1s

        # =====================================================================
                
        # =====================================================================
        ## Y(1S)
        # =====================================================================
        self.__aL    = self.make_var ( alphaL                  ,
                                       "aL_%s"          % name ,
                                       "#alpha_{L}(%s)" % name ,
                                       True , 1.5462     , 0.1   , 10 )
        self.__nL    = self.make_var ( nL                      ,                     
                                       "nL_%s"          % name ,
                                       "n_{L}(%s)"      % name ,
                                       True , 1.3119     , 1.e-5 , 25 )
        self.__aR    = self.make_var ( alphaR                  ,
                                       "aR_%s"          % name ,
                                       "#alpha_{R}(%s)" % name ,
                                       True , 1.6952e+00 , 0.1   , 10 )
        self.__nR    = self.make_var ( nR                      ,
                                       "nR_%s"          % name ,
                                       "n_{R}(%s)"      % name ,
                                       True , 1.5751e+01 , 1.e-5 , 25 )
        
        self.__m1s   = self.make_var ( m1s                   ,
                                       "m1S_%s"       % name ,
                                       "mass Y1S(%s)" % name , m1s ,
                                       None  , 
                                       m_y1s , m_y1s - 0.15 * s_y1s , m_y1s + 0.15 * s_y1s ) 
        
        self.__s1s   = self.make_var ( sigma                  ,
                                       "s1S_%s"        % name ,
                                       "sigma Y1S(%s)" % name ,
                                       None  , s_y1s , 0.3 * s_y1s , 4 * s_y1s )
        self.__sigma = self.__s1s
        
        self.__Y1S   = CB2_pdf (
            name + '1S'             ,
            xvar     = self.mass    ,
            mean     = self.__m1s   ,
            sigma    = self.__s1s   ,
            alphaL   = self.__aL    ,
            alphaR   = self.__aR    ,
            nL       = self.__nL    ,
            nR       = self.__nR    )
        
        # =====================================================================
        ## Y(2S)
        # =====================================================================
        self.__dm2s  = self.make_var ( None                  ,
                                       "dm2s"      + name    ,
                                       "dm2s(%s)"  % name    ,
                                       True                  , 
                                       dm_y2s                ,
                                       dm_y2s - 0.20 * s_y1s , 
                                       dm_y2s + 0.20 * s_y1s )
        
        self.__aset11 = ROOT.RooArgList ( self.__m1s , self.__dm2s )
        
        ## self.__m2s    = ROOT.RooFormulaVar (
        ##     "m_" + name + '2S'   ,
        ##     "m2s(%s)"  % name    ,
        ##     "%s+%s" % ( self.__m1s.GetName() , self.__dm2s.GetName()  ) , 
        ##    self.__aset11       )
        self.__m2s    = Ostap.MoreRooFit.Addition ( "m_" + name + '2S' ,
                                                    "m2s(%s)"  % name  ,
                                                    self.__m1s , self.__dm2s ) 
        
        
        self.__aset12 = ROOT.RooArgList ( self.__sigma , self.__m1s , self.__m2s ) 
        self.s2s    = Ostap.FormulaVar (
            "sigma_"  + name + '2S'    ,
            "#sigma_{Y2S}(%s)" % name  ,
            "%s*(%s/%s)"  % ( self.__sigma.GetName() ,
                              self.__m2s  .GetName() ,
                              self.__m1s  .GetName() ) ,
            self.__aset12  )
        
        self.__Y2S  = CB2_pdf (
            name + '2S'             ,
            xvar     = self.mass    ,
            mean     = self.__m2s   ,
            sigma    = self.__s2s   ,
            alphaL   = self.__aL    ,
            alphaR   = self.__aR    ,
            nL       = self.__nL    ,
            nR       = self.__nR    )
                
        # =====================================================================
        ## Y(3S)
        # =====================================================================
        self.__dm3s  = self.make_var ( None                  ,
                                       "dm3s"      + name    ,
                                       "dm3s(%s)"  % name    ,
                                       True                  , 
                                       dm_y3s                ,
                                       dm_y3s - 0.20 * s_y1s , 
                                       dm_y3s + 0.20 * s_y1s )
        
        self.__aset21 = ROOT.RooArgList ( self.__m1s , self.__dm3s )
        ## self.__m3s    = ROOT.RooFormulaVar (
        ##    "m_"       + name + '(3S)' ,
        ##    "m3s(%s)"  % name          ,
        ##    "%s+%s" % ( self.__m1s.GetName() , self.__dm3s.GetName() ) ,
        ##    self.__aset21       )
        self.__m3s    = Ostap.MoreRooFit.Addition ( "m_" + name + '3S' ,
                                                    "m3s(%s)"  % name  ,
                                                    self.__m1s , self.__dm3s ) 
        
        
        self.__aset22 = ROOT.RooArgList ( self.__sigma , self.__m1s , self.__m3s ) 
        self.__s3s    = Ostap.FormulaVar (
            "sigma_"  + name + '3S'    ,
            "#sigma_{Y3S}(%s)" % name  ,
            "%s*(%s/%s)"  % ( self.__sigma.GetName() ,
                              self.__m3s  .GetName() ,
                              self.__m1s  .GetName() ) , 
            self.__aset22       )
        
        self.__Y3S  = CB2_pdf (
            name + '3S'           ,
            xvar     = self.mass  ,
            mean     = self.__m3s ,
            sigma    = self.__s3s ,
            alphaL   = self.__aL  ,
            alphaR   = self.__aR  ,
            nL       = self.__nL  ,
            nR       = self.__nR  )
        
        #
        ## the actual signal PDFs
        # 
        self.__y1s   = self.__Y1S.pdf
        self.__y2s   = self.__Y2S.pdf
        self.__y3s   = self.__Y3S.pdf
        
        
        ## use helper function to create background  
        self.background = self.make_bkg ( power , 'Bkg%s' % name , self.mass )

        
        self.__n1s = self.make_var ( None ,
                                     "N1S" + name  ,
                                     "Signal(Y1S)" ,  None , 1000 ,  0 ,  1.e+7 )
        self.__n2s = self.make_var ( None ,
                                     "N2S" + name  ,
                                     "Signal(Y2S)" ,  None ,  300 ,  0 ,  1.e+6 )
        self.__n3s = self.make_var ( None ,
                                     "N3S" + name  ,
                                     "Signal(Y3S)" ,  None ,  100 ,  0 ,  1.e+6 )
        self.__b   = self.make_var ( None ,
                                     "B"   + name  ,
                                     "Background"  ,  None ,  100 ,  0 ,  1.e+8 )
        
        self.alist1 = ROOT.RooArgList ( self.__y1s , self.__y2s , self.__y3s ) 
        self.alist2 = ROOT.RooArgList ( self.__n1s , self.__n2s , self.__n3s ) 
        
        self.alist1 . add ( self.background.pdf )
        self.alist2 . add ( self.__b            )
        
        self.pdf  = ROOT.RooAddPdf  (
            "manca_%s"  % name ,
            "manca(%s)" % name ,
            self.alist1 ,
            self.alist2 )
        
        self.__dm2s.setConstant ( True )
        self.__dm3s.setConstant ( True )
        
        self._splots = []
        
        self.s1_name = self.N1S.GetName ()
        self.s2_name = self.N2S.GetName ()
        self.s3_name = self.N3S.GetName ()
        self.b_name  = self.B  .GetName ()

        # 
        ## finally declare components 
        #
        self.signals    () . add ( self.__y1s )
        self.signals    () . add ( self.__y2s )
        self.signals    () . add ( self.__y3s )
        self.backgrounds() . add ( self.background.pdf ) 

        ## save configurtaion
        self.config = {
            'mass'   : self.mass  ,
            'name'   : self.name  ,
            'power'  : self.power ,
            'm1s'    : self.m1s   ,
            'sigma'  : self.sigma ,
            'alphaL' : self.aL    ,
            'alphaR' : self.aR    ,
            'nL'     : self.nL    ,
            'nR'     : self.nR    ,
            }

    @property
    def mass ( self ) :
        """'mass'-variable, ailas for 'x' and 'xvar'"""
        return self.xvar
    
    @property
    def  s1s (  self ) :
        """'s1s'-parameter (resoltuion for Y(1S) peak) (same as 'sigma')"""
        return self.__sigma
    @s1s.setter
    def  s1s (  self , value ) :
        self.set_value ( self.__s1s , value )
    
    @property
    def  sigma (  self ) :
        """'sigma'-parameter (resoltuion for Y(1S) peak) (same as 's1s')"""
        return self.s1s
    @sigma.setter
    def  sigma (  self , value ) :
        self.s1s = value

    @property
    def  m1s (  self ) :
        """'m1s'-parameter (mass for Y(1S) peak)"""
        return self.__m1s
    @m1s.setter
    def  m1s (  self , value ) :
        self.set_value ( self.__m1s , value )

    @property
    def  m2s (  self ) :
        """'m2s'-parameter (mass for Y(2S) peak)"""
        return self.__m2s
    @property
    def  m3s (  self ) :
        """'m3s'-parameter (mass for Y(3S) peak)"""
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
    def  dm2s (  self ) :
        """'dm2s'-parameter (mass  difference for Y(2S) and Y(1S) peaks)"""
        return self.__dm2s    
    @dm2s.setter
    def  dm2s (  self , value ) :
        self.set_value ( self.__dm2s , value ) 
        
    @property
    def  dm3s (  self ) :
        """'dm3s'-parameter (mass  difference for Y(3S) and Y(1S) peaks)"""
        return self.__dm3s    
    @dm3s.setter
    def  dm3s (  self , value ) :
        self.Set_value ( self.__dm3s , value )
        
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
        self.set_avlue ( self.__aR , value )
        
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
        
models.append ( Manca2_pdf ) 
# =============================================================================
## Specific model for fitting of Y+X
class MancaX_pdf(PDF2) :
    """MancaX: 2D-model to study associative production of Upsilon and X 
    """
    def __init__ ( self         ,
                   manca        , ## manca pdf, that defined 3 upsilon peaks  
                   charm        , ## charm pdf
                   bkg1   = 0   ,
                   bkg2   = 0   ,
                   bkgA   = 0   ,
                   bkgB   = 0   ,
                   suffix = ''  ) :

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
            'mass'   : self.mass     ,
            'manca'  : self.name     ,
            'charm'  : self.charma   ,
            'bkg1'   : self.b_Y      , 
            'bkg2'   : self.b_C      , 
            'bkgA'   : self.b_A      , 
            'bkgB'   : self.b_B      , 
            'suffix' : self.__siffix
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
