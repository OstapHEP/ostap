#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
# =============================================================================
## @file specific.py
#  A set of predefined ready-to-use shapes and PDFs 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
# 
#                    $Revision$
#  Last modification $Date$
#                 by $Author$
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
    )
# =============================================================================
import ROOT, math
# =============================================================================
from   ostap.logger.logger     import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.specific' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================
from   ostap.core.core           import cpp, Ostap, VE 
# =============================================================================
# Specializations of double-sided Crystal Ball function 
# =============================================================================
from   ostap.fitting.basic   import PDF , PDF2 
from   ostap.fitting.signals import CB2_pdf
models = [] 
# =============================================================================
## @class Bd_pdf
#  simple wrapper over CB2-pdf
#  @see Analysis::Models::CrystalBallDS
#  @see Gaudi::Math::CrystalBallDS
#  @attention: mass is mandatory variable! 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Bd_pdf(CB2_pdf) :
    """B0: double-sided Crystal Ball function
    """
    def __init__ ( self                   ,
                   mass                   ,   ## mass is mandatory here! 
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
                           mass.getMin()    ,
                           mass.getMax()    ,
                           mass             , 
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
#  @see Analysis::Models::CrystalBallDS
#  @see Gaudi::Math::CrystalBallDS
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Bu_pdf(CB2_pdf) :
    """B+: double-sided Crystal Ball function
    """
    def __init__ ( self                   ,
                   mass                   ,   ## mass is mandatory here! 
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
                           mass.getMin()    ,
                           mass.getMax()    ,
                           mass             , 
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
#  @see Analysis::Models::CrystalBallDS
#  @see Gaudi::Math::CrystalBallDS
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Bs_pdf(CB2_pdf) :
    """Bs: double-sided Crystal Ball function
    """
    def __init__ ( self                   ,
                   mass                   ,    ## mass is mandatory here! 
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
                           mass.getMin()    ,
                           mass.getMax()    ,
                           mass             , 
                           mean             ,
                           sigma            ,
                           alphaL           ,
                           alphaR           ,
                           nL               ,
                           nR               )

models.append ( Bs_pdf ) 
# =============================================================================
## @class Bc_pdf
#  simple wrapper over CB2-pdf
#  @see Analysis::Models::CrystalBallDS
#  @see Gaudi::Math::CrystalBallDS
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Bc_pdf(CB2_pdf) :
    """Bc: double-sided Crystal Ball function
    """
    def __init__ ( self                   ,
                   mass                   ,   ## mass is mandatory here! 
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
                           mass.getMin()    ,
                           mass.getMax()    ,
                           mass             , 
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
from   ostap.fitting.signals   import Bukin_pdf 
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
                   mass                 , ## mass is mandatory here! 
                   name   = 'D0'        ,
                   mean   =  1.8648e+00 , 
                   sigma  =  7.3651e-03 , 
                   xi     = -1.7616e-03 , 
                   rhoL   =  3.7311e-01 ,
                   rhoR   =  4.4033e-01 ) : 
        
        Bukin_pdf.__init__ ( self          ,
                             name          ,
                             mass.getMin() ,
                             mass.getMax() ,
                             mass          , 
                             mean          ,
                             sigma         ,
                             xi            , 
                             rhoL          ,
                             rhoR          ) 
                             
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
                   mass                    , ## mass is mandatory here 
                   name     = 'Dp'         ,
                   mean     =  1.869       , 
                   sigma    =  7.1183e-03  ,
                   xi       = -7.7344e-03  ,
                   rhoL     =  3.0241e-01  , 
                   rhoR     =  3.7452e-01  ) :

        Bukin_pdf.__init__ ( self          ,
                             name          ,
                             mass.getMin() , 
                             mass.getMax() , 
                             mass          ,
                             mean          ,
                             sigma         ,
                             xi            ,                            
                             rhoL          ,
                             rhoR          ) 
        
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
                   mass     = None         , ## mass is mandatory 
                   name     = 'Ds'         ,
                   mean     =  1.969       ,
                   sigma    =  0.0068      ,
                   xi       = -6.45755e-04 ,
                   rhoL     =  3.0241e-01  , 
                   rhoR     =  3.7452e-01  ) :
        
        Bukin_pdf.__init__ ( self          ,
                             name          ,
                             mass.getMin() , 
                             mass.getMax() , 
                             mass          , 
                             mean          ,
                             sigma         ,
                             xi            ,
                             rhoL          ,
                             rhoR          ) 
        
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
                   mass                     , 
                   name     = 'Lc'          ,                   
                   mean     =  2.28590e+00  ,
                   sigma    =  5.11874e-03  ,
                   xi       =  1.82493e-03  ,
                   rhoL     =  3.0241e-01   , 
                   rhoR     =  3.7452e-01   ) :
        
        Bukin_pdf.__init__ ( self     ,
                             name     ,
                             mass.getMin() , 
                             mass.getMax() , 
                             mass     , 
                             mean     ,
                             sigma    ,
                             xi       ,
                             rhoL     ,
                             rhoR     ) 
        
models.append ( Lc_pdf ) 
# =============================================================================
## @class Manca_pdf 
#  the final full PDF for Y->mu+mu- fit
#  This is physically well-motivated function for fits in narrow
#  bins in pt and rapidity  
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-08-02
class Manca_pdf (PDF) :
    """ Manca function
    Fhe final fit model for Y->mu+mu- fit
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
        PDF.__init__ ( self , name )
        #
        if     mass.getMin() <  9.460 and   9.60  <= mass.getMax()  : gev_ =    1
        elif   mass.getMin() < 10.    and  10.500 <= mass.getMax()  : gev_ =    1
        elif   mass.getMin() < 10.0   and  10.200 <= mass.getMax()  : gev_ =    1
        elif   mass.getMin() <  9460  and  10355  <= mass.getMax()  : gev_ = 1000
        elif   mass.getMin() < 10000  and  10500  <= mass.getMax()  : gev_ = 1000
        elif   mass.getMin() < 10000  and  10200  <= mass.getMax()  : gev_ = 1000
        else : raise TypeError ( "Illegal mass range %s<m<%s"
                                 % ( mass.getMin() , mass.getMax() ) ) 
        
        m_y1s  =  9.46030     * gev_ 
        s_y1s  =  4.3679e-02  * gev_ 
        dm_y2s = 10.02326     * gev_ - m_y1s
        dm_y3s = 10.3552      * gev_ - m_y1s

        # 
        self.mass = mass

        # =====================================================================
        from   ostap.fitting.basic   import makeVar
        from   ostap.fitting.signals import Needham_pdf
        # =====================================================================
        ## Y(1S)
        # =====================================================================
        
        self.a0   = makeVar ( a0                 ,
                              'a0m_%s' % name    ,
                              "a0 for Needham's function" , a0 , 
                              1.91               ,    0.1             ,   3.0            )

        self.a1   = makeVar ( a1                 ,
                              'a1m_%s' % name    ,
                              "a1 for Needham's function" , a1 ,  
                              1.1174 / gev_      ,  -10.0  / gev_     ,  10.0  /gev_     )
        
        self.a2   = makeVar ( a2                 ,
                              'a2m_%s' % name    ,
                              "a2 for Needham's function" , a2 , 
                              -5.299 / gev_**2   , -100.0  / gev_**2  , 100.0  /gev_**2  )

        ## default logic does not work nicely here, therefore we need to be explicit:
        if a0 is None and not self.a0.isConstant () : self.a0.fix (  1.91             ) 
        if a1 is None and not self.a1.isConstant () : self.a1.fix (  1.1174 / gev_    )
        if a2 is None and not self.a2.isConstant () : self.a2.fix ( -5.299  / gev_**2 ) 
        
        self.m1s   = makeVar ( m1s                   ,
                               "m1S_%s"       % name ,
                               "mass Y1S(%s)" % name , m1s ,
                               m_y1s , m_y1s - 0.15 * s_y1s , m_y1s + 0.15 * s_y1s ) 
        
        self.s1s   = makeVar ( sigma                  ,
                               "s1S_%s"        % name ,
                               "sigma Y1S(%s)" % name , sigma ,
                               s_y1s , 0.3 * s_y1s , 4 * s_y1s )
        
        self.sigma = self.s1s 
        self.Y1S   = Needham_pdf (
            name + '1S'           ,
            mass.getMin()         ,
            mass.getMax()         ,
            mass     = self.mass  ,
            mean     = self.m1s   ,
            sigma    = self.s1s   ,
            a0       = self.a0    ,
            a1       = self.a1    ,
            a2       = self.a2    ) 
        
        # =====================================================================
        ## Y(2S)
        # =====================================================================
        self.dm2s  = makeVar ( None ,
                               "dm2s"      + name    ,
                               "dm2s(%s)"  % name    ,
                               dm_y2s                ,
                               dm_y2s - 0.20 * s_y1s , 
                               dm_y2s + 0.20 * s_y1s )
        
        self.aset11 = ROOT.RooArgList ( self.m1s , self.dm2s )
        self.m2s    = ROOT.RooFormulaVar (
            "m_" + name + '2S'   ,
            "m2s(%s)"  % name    ,
            "%s+%s" % ( self.m1s.GetName() , self.dm2s.GetName()  ) , 
            self.aset11       )
        
        self.aset12 = ROOT.RooArgList ( self.sigma , self.m1s , self.m2s ) 
        self.s2s    = ROOT.RooFormulaVar (
            "sigma_"  + name + '2S'    ,
            "#sigma_{Y2S}(%s)" % name  ,
            "%s*(%s/%s)"  % ( self.sigma.GetName() ,
                              self.m2s  .GetName() ,
                              self.m1s  .GetName() ) ,
            self.aset12  )
        
        self.Y2S   = Needham_pdf (
            name + '2S'           ,
            mass.getMin()         ,
            mass.getMax()         ,
            mass     = self.mass  ,
            mean     = self.m2s   ,
            sigma    = self.s2s   ,
            a0       = self.a0    ,
            a1       = self.a1    ,
            a2       = self.a2    ) 
        
        # =====================================================================
        ## Y(3S)
        # =====================================================================
        self.dm3s  = makeVar ( None ,
                               "dm3s"      + name ,
                               "dm3s(%s)"  % name    ,
                               dm_y3s                ,
                               dm_y3s - 0.20 * s_y1s , 
                               dm_y3s + 0.20 * s_y1s )
        
        self.aset21 = ROOT.RooArgList ( self.m1s , self.dm3s )
        self.m3s    = ROOT.RooFormulaVar (
            "m_"       + name + '(3S)' ,
            "m3s(%s)"  % name          ,
            "%s+%s" % ( self.m1s.GetName() , self.dm3s.GetName() ) ,
            self.aset21       )
        
        
        self.aset22 = ROOT.RooArgList ( self.sigma , self.m1s , self.m3s ) 
        self.s3s    = ROOT.RooFormulaVar (
            "sigma_"  + name + '3S'    ,
            "#sigma_{Y3S}(%s)" % name  ,
            "%s*(%s/%s)"  % ( self.sigma.GetName() ,
                              self.m3s  .GetName() ,
                              self.m1s  .GetName() ) , 
            self.aset22       )
        
        self.Y3S   = Needham_pdf (
            name + '3S'           ,
            mass.getMin()         ,
            mass.getMax()         ,
            mass     = self.mass  ,
            mean     = self.m3s   ,
            sigma    = self.s3s   ,
            a0       = self.a0    ,
            a1       = self.a1    ,
            a2       = self.a2    ) 
        
        #
        ## the actual signal PDFs
        # 
        self.y1s   = self.Y1S.pdf
        self.y2s   = self.Y2S.pdf
        self.y3s   = self.Y3S.pdf
        
        ## use helper function to create background  
        from ostap.fitting.models_bkg import makeBkg
        self.background = makeBkg ( power , 'Bkg%s' % name , self.mass )

        self.n1s = makeVar ( None ,
                             "N1S" + name  ,
                             "Signal(Y1S)" ,  None , 1000 ,  0 ,  1.e+7 )
        self.n2s = makeVar ( None ,
                             "N2S" + name  ,
                             "Signal(Y2S)" ,  None ,  300 ,  0 ,  1.e+6 )
        self.n3s = makeVar ( None ,
                             "N3S" + name  ,
                             "Signal(Y3S)" ,  None ,  100 ,  0 ,  1.e+5 )
        self.b   = makeVar ( None ,
                             "B"   + name  ,
                             "Background"  ,  None ,  100 ,  0 ,  1.e+8 )
        
        self.alist1 = ROOT.RooArgList ( self.y1s , self.y2s , self.y3s ) 
        self.alist2 = ROOT.RooArgList ( self.n1s , self.n2s , self.n3s ) 
        
        self.alist1 . add ( self.background.pdf )
        self.alist2 . add ( self.b              )
        
        self.pdf  = ROOT.RooAddPdf  (
            "manca_%s"  % name ,
            "manca(%s)" % name ,
            self.alist1 ,
            self.alist2 )
        
        self.dm2s.setConstant ( True )
        self.dm3s.setConstant ( True )
        
        self._splots = []
        
        self.s1_name = self.n1s.GetName ()
        self.s2_name = self.n2s.GetName ()
        self.s3_name = self.n3s.GetName ()

        # 
        ## finally declare components 
        #
        self.signals    () . add ( self.y1s )
        self.signals    () . add ( self.y2s )
        self.signals    () . add ( self.y3s )
        self.backgrounds() . add ( self.background.pdf ) 
        
    def alpha_1S ( self ) : return self.Y1S.pdf.alpha ()
    def alpha_2S ( self ) : return self.Y2S.pdf.alpha ()
    def alpha_3S ( self ) : return self.Y3S.pdf.alpha ()

    
models.append ( Manca_pdf ) 
# =============================================================================
## @class Manca2_pdf 
#  the final full PDF for Y->mu+mu- fit
#  This is an effective function for fit in global bin, without pt/y-binning 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-06-24
class Manca2_pdf (PDF) :
    """
    Manca2: the final fit model for Y->mu+mu- fit
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
        PDF.__init__ ( self , name )
        # 
        if     mass.getMin() <  9.460 and   9.60  <= mass.getMax()  : gev_ =    1
        elif   mass.getMin() < 10.    and  10.500 <= mass.getMax()  : gev_ =    1
        elif   mass.getMin() < 10.0   and  10.200 <= mass.getMax()  : gev_ =    1
        elif   mass.getMin() <  9460  and  10355  <= mass.getMax()  : gev_ = 1000
        elif   mass.getMin() < 10000  and  10500  <= mass.getMax()  : gev_ = 1000
        elif   mass.getMin() < 10000  and  10200  <= mass.getMax()  : gev_ = 1000
        else : raise TypeError ( "Illegal mass range %s<m<%s"
                                 % ( mass.getMin() , mass.getMax() ) ) 

        m_y1s  =  9.46030     * gev_
        s_y1s  =  4.03195e-02 * gev_ 
        dm_y2s = 10.02326     * gev_ - m_y1s
        dm_y3s = 10.3552      * gev_ - m_y1s

        # 
        self.mass = mass

        # =====================================================================
        from   ostap.fitting.basic   import makeVar
        from   ostap.fitting.signals import CB2_pdf
                
        # =====================================================================
        ## Y(1S)
        # =====================================================================
        self.aL    = makeVar ( alphaL                  ,
                               "aL_%s"          % name ,
                               "#alpha_{L}(%s)" % name , alphaL , 1.5462     , 0 , 10 )
        self.nL    = makeVar ( nL                      ,                     
                               "nL_%s"          % name ,
                               "n_{L}(%s)"      % name , nL     , 1.3119     , 0 , 10 )
        self.aR    = makeVar ( alphaR                  ,
                               "aR_%s"          % name ,
                               "#alpha_{R}(%s)" % name , alphaR , 1.6952e+00 , 0 , 10 )
        self.nR    = makeVar ( nR                      ,
                               "nR_%s"          % name ,
                               "n_{R}(%s)"      % name , nR     , 1.5751e+01 , 0 , 25 )
        
        self.m1s   = makeVar ( m1s                   ,
                               "m1S_%s"       % name ,
                               "mass Y1S(%s)" % name , m1s ,
                               m_y1s , m_y1s - 0.15 * s_y1s , m_y1s + 0.15 * s_y1s ) 
        
        self.s1s   = makeVar ( sigma                  ,
                               "s1S_%s"        % name ,
                               "sigma Y1S(%s)" % name , sigma ,
                               s_y1s , 0.3 * s_y1s , 4 * s_y1s )
        self.sigma = self.s1s
        
        self.Y1S   = CB2_pdf (
            name + '1S'           ,
            mass.getMin()         ,
            mass.getMax()         ,
            mass     = self.mass  ,
            mean     = self.m1s   ,
            sigma    = self.s1s   ,
            alphaL   = self.aL    ,
            alphaR   = self.aR    ,
            nL       = self.nL    ,
            nR       = self.nR    )
        
        # =====================================================================
        ## Y(2S)
        # =====================================================================
        self.dm2s  = makeVar ( None                  ,
                               "dm2s"      + name    ,
                               "dm2s(%s)"  % name    ,
                               dm_y2s                ,
                               dm_y2s - 0.20 * s_y1s , 
                               dm_y2s + 0.20 * s_y1s )
        
        self.aset11 = ROOT.RooArgList ( self.m1s , self.dm2s )
        self.m2s    = ROOT.RooFormulaVar (
            "m_" + name + '2S'   ,
            "m2s(%s)"  % name    ,
            "%s+%s" % ( self.m1s.GetName() , self.dm2s.GetName()  ) , 
            self.aset11       )
        
        self.aset12 = ROOT.RooArgList ( self.sigma , self.m1s , self.m2s ) 
        self.s2s    = ROOT.RooFormulaVar (
            "sigma_"  + name + '2S'    ,
            "#sigma_{Y2S}(%s)" % name  ,
            "%s*(%s/%s)"  % ( self.sigma.GetName() ,
                              self.m2s  .GetName() ,
                              self.m1s  .GetName() ) ,
            self.aset12  )
        
        self.Y2S  = CB2_pdf (
            name + '2S'           ,
            mass.getMin()         ,
            mass.getMax()         ,
            mass     = self.mass  ,
            mean     = self.m2s   ,
            sigma    = self.s2s   ,
            alphaL   = self.aL    ,
            alphaR   = self.aR    ,
            nL       = self.nL    ,
            nR       = self.nR    )
                
        # =====================================================================
        ## Y(3S)
        # =====================================================================
        self.dm3s  = makeVar ( None                  ,
                               "dm3s"      + name    ,
                               "dm3s(%s)"  % name    ,
                               dm_y3s                ,
                               dm_y3s - 0.20 * s_y1s , 
                               dm_y3s + 0.20 * s_y1s )
        
        self.aset21 = ROOT.RooArgList ( self.m1s , self.dm3s )
        self.m3s    = ROOT.RooFormulaVar (
            "m_"       + name + '(3S)' ,
            "m3s(%s)"  % name          ,
            "%s+%s" % ( self.m1s.GetName() , self.dm3s.GetName() ) ,
            self.aset21       )
        
        
        self.aset22 = ROOT.RooArgList ( self.sigma , self.m1s , self.m3s ) 
        self.s3s    = ROOT.RooFormulaVar (
            "sigma_"  + name + '3S'    ,
            "#sigma_{Y3S}(%s)" % name  ,
            "%s*(%s/%s)"  % ( self.sigma.GetName() ,
                              self.m3s  .GetName() ,
                              self.m1s  .GetName() ) , 
            self.aset22       )
        
        self.Y3S  = CB2_pdf (
            name + '3S'           ,
            mass.getMin()         ,
            mass.getMax()         ,
            mass     = self.mass  ,
            mean     = self.m3s   ,
            sigma    = self.s3s   ,
            alphaL   = self.aL    ,
            alphaR   = self.aR    ,
            nL       = self.nL    ,
            nR       = self.nR    )
        
        #
        ## the actual signal PDFs
        # 
        self.y1s   = self.Y1S.pdf
        self.y2s   = self.Y2S.pdf
        self.y3s   = self.Y3S.pdf
        
        
        ## use helper function to create background  
        from ostap.fitting.models_bkg import makeBkg
        self.background = makeBkg ( power , 'Bkg%s' % name , self.mass )

        
        self.n1s = makeVar ( None ,
                             "N1S" + name  ,
                             "Signal(Y1S)" ,  None , 1000 ,  0 ,  1.e+7 )
        self.n2s = makeVar ( None ,
                             "N2S" + name  ,
                             "Signal(Y2S)" ,  None ,  300 ,  0 ,  1.e+6 )
        self.n3s = makeVar ( None ,
                             "N3S" + name  ,
                             "Signal(Y3S)" ,  None ,  100 ,  0 ,  1.e+6 )
        self.b   = makeVar ( None ,
                             "B"   + name  ,
                             "Background"  ,  None ,  100 ,  0 ,  1.e+8 )
        
        self.alist1 = ROOT.RooArgList ( self.y1s , self.y2s , self.y3s ) 
        self.alist2 = ROOT.RooArgList ( self.n1s , self.n2s , self.n3s ) 
        
        self.alist1 . add ( self.background.pdf )
        self.alist2 . add ( self.b              )
        
        self.pdf  = ROOT.RooAddPdf  (
            "manca_%s"  % name ,
            "manca(%s)" % name ,
            self.alist1 ,
            self.alist2 )
        
        self.dm2s.setConstant ( True )
        self.dm3s.setConstant ( True )
        
        self._splots = []
        
        self.s1_name = self.n1s.GetName ()
        self.s2_name = self.n2s.GetName ()
        self.s3_name = self.n3s.GetName ()

        # 
        ## finally declare components 
        #
        self.signals    () . add ( self.y1s )
        self.signals    () . add ( self.y2s )
        self.signals    () . add ( self.y3s )
        self.backgrounds() . add ( self.background.pdf ) 


models.append ( Manca2_pdf ) 
# =============================================================================
## Specific model for fitting of Y+X
class MancaX_pdf(PDF2) :
    """
    MancaX: 2D-model to study associative production of Upsilon and X 
    """
    def __init__ ( self         ,
                   manca        , ## manca pdf, that defined 3 upsolon peaks  
                   charm        , ## charm pdf
                   bkg1   = 0   ,
                   bkg2   = 0   ,
                   bkgA   = 0   ,
                   bkgB   = 0   ,
                   suffix = ''  ) :

        PDF2.__init__ ( self , 'MancaX' + suffix , manca.mass , charm.mass )  
        self._crossterms1 = ROOT.RooArgSet ()
        self._crossterms2 = ROOT.RooArgSet ()
                         
        self.suffix  = suffix
        self.signal1 = manca  
        self.signal2 = charm

        ## use helper function to create background  
        from ostap.fitting.models_bkg import makeBkg
        #
        ## background components
        #
        self.b_Y = makeBkg ( bkg1 , 'BkgY' + suffix , self.m1 )   
        self.b_C = makeBkg ( bkg2 , 'BkgC' + suffix , self.m2 )   
        self.b_A = makeBkg ( bkgA , 'BkgA' + suffix , self.m1 )   
        self.b_B = makeBkg ( bkgB , 'BkgB' + suffix , self.m2 )
        
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
        from Ostap.FitBasic import makeVar 
        self.s1s = makeVar ( None , 'NY1C' + suffix , 'N(Y1S+C)' + suffix , None , 200 , 0 , 1.e+5 )
        self.s2s = makeVar ( None , 'NY2C' + suffix , 'N(Y2S+C)' + suffix , None , 100 , 0 , 1.e+4 )
        self.s3s = makeVar ( None , 'NY3C' + suffix , 'N(Y3S+C)' + suffix , None ,  20 , 0 , 1.e+4 )
        self.bs  = makeVar ( None , 'NBC'  + suffix , 'N(Bkg+C)' + suffix , None , 500 , 0 , 1.e+5 )
        self.s1b = makeVar ( None , 'NY1B' + suffix , 'N(Y1S+B)' + suffix , None , 500 , 0 , 1.e+5 )
        self.s2b = makeVar ( None , 'NY2B' + suffix , 'N(Y2S+B)' + suffix , None , 200 , 0 , 1.e+4 )
        self.s3b = makeVar ( None , 'NY3B' + suffix , 'N(Y3S+B)' + suffix , None , 200 , 0 , 1.e+4 )
        self.bb  = makeVar ( None , 'NBB'  + suffix , 'N(Bkg+B)' + suffix , None , 500 , 0 , 1.e+5 )

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
    
    ## get all declared components 
    def crossterms1 ( self ) : return self._crossterms1
    ## get all declared components 
    def crossterms2 ( self ) : return self._crossterms2

models.append ( MancaX_pdf ) 
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )

 
# =============================================================================
# The END 
# =============================================================================
