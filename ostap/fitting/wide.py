#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/wide .py
#
#  Set of useful PDFs for various `wide' signals 
#  It includes
#  - soeme empricial PDFs to describe narrow peaks: Gauss, CrystalBall, ....
#  - some PDF to describe "wide" peaks: BreitWigner,LASS, Bugg, Flatte, ...
#  - some useful PDFs to describe smooth background: phase space ;
#    expo times polynomial; phase space times polynomial, ...
#  - set of smooth non-facrorizable model for 2D fits 
#
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
# 
# =============================================================================
""" Set of useful PDFs for various `wide` signals

PDF to describe 'wide' peaks

  - BreitWigner
  - BreitWigner with interference 
  - LASS
  - Bugg
  - Flatte
  - Swanson's S=wave cusp 
  - ...

Special stuff:

  - Voigt & PseudoVoigt
  - BW3L

"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
# =============================================================================
__all__ = (
    # =========================================================================
    'Voigt_pdf'              , ## Voigt-profile
    'PseudoVoigt_pdf'        , ## PseudoVoigt-profile
    'BreitWigner_pdf'        , ## (relativistic) 2-body Breit-Wigner
    'BWMC_pdf'               , ## Multi-channel version of Breit-Wigner function
    'BWI_pdf'                , ## (relativistic) Breit-Wigner with interference
    'BWPS_pdf'               , ## Breit-Wigner function modulated with extra phase-space and polynomial factors
    'BW3L_pdf'               , ## Breit-Wigner function modulated with p**(2L+1) factor 
    'Flatte_pdf'             , ## Flatte-function  (pipi/KK)
    'FlattePS_pdf'           , ## Flatte-PS function     
    'FlatteBugg_pdf'         , ## Flatte-Bugg function  (pipi)
    'LASS_pdf'               , ## kappa-pole
    'Bugg_pdf'               , ## sigma-pole
    ##
)
# =============================================================================
from   ostap.core.core          import Ostap
from   ostap.fitting.funbasic   import FUN1   , Fun1D 
from   ostap.fitting.pdfbasic   import PDF1   , all_args
from   ostap.fitting.fit1d      import PEAK   , PEAKMEAN , CheckMean
from   ostap.fitting.fithelpers import Phases , ZERO     
from   ostap.fitting.variables  import SETVAR
import ostap.math.dalitz
import ostap.math.models
import ROOT, math
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.wide' )
else                       : logger = getLogger ( __name__             )
# =============================================================================
models = [] 
# =============================================================================
## @class Voigt_pdf
#  Voigt-pdf distribution
#  @see Ostap::Models::Voigt
#  @see Ostap::Math::Voigt
#  The implementation relied on Faddeeva function 
#  @see http://en.wikipedia.org/wiki/Faddeeva_function
#  @see http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class Voigt_pdf(PEAK) :
    """ Voigt function:
    Convolution of non-relativistic Breit-Wigner with Gaussian resolution
    
    The implementation relied on Faddeeva function 
    http://en.wikipedia.org/wiki/Faddeeva_function
    http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package
    
    Parameters
    - mean  : location 
    - gamma : gamma for breight-wigner pole
    - sigma : resolution parameter for gaussian 
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   m0        = None ,
                   sigma     = None ,
                   gamma     = None ) :
        #
        ## initialize the base
        # 
        PEAK.__init__  ( self , name , xvar ,
                              mean        = m0                  ,
                              sigma       = sigma               ,
                              mean_name   = 'm0_%s'      % name ,
                              mean_title  = '#m_{0}(%s)' % name ) 
                         
        limits_gamma = ()
        if  self.xminmax() :
            mn , mx = self.xminmax() 
            dm = mx - mn
            limits_gamma = 1.e-5 * dm , dm
            
        self.__gamma  = self.make_var ( gamma                ,
                                        'gamma_%s'   % name  ,   
                                        '#gamma(%s)' % name  ,
                                        None , *limits_gamma )
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.Voigt (
            self.roo_name ( 'voigt_' ) , 
            "Voigt %s" % self.name ,
            self.xvar   ,
            self.m0     ,
            self.gamma  ,
            self.sigma  )

        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'm0'        : self.m0    ,
            'sigma'     : self.sigma ,
            'gamma'     : self.gamma ,
            }


    @property
    def gamma ( self ) :
        """'gamma'-parameter for Voigt function"""
        return self.__gamma
    @gamma.setter
    def gamma ( self, value ) :
        value = float ( value )
        assert 0 < value , "'gamma'-parameter must be positive"
        self.__gamma.setVal ( value ) 
        return self.__gamma.getVal ()

    ## ALIASES 
    m0    = PEAK.mean 
    Gamma = gamma 

models.append ( Voigt_pdf )                          

# =============================================================================
## @class PseudoVoigt_pdf
#  PseudoVoigt-pdf distribution
#  @see Ostap::Models::PseudoVoigt
#  @see Ostap::Math::PseudoVoigt
#  CPU-efficient Approximation of Voight profile
#  @see T. Ida, M. Ando and H. Toraya, 
#       "Extended pseudo-Voigt function for approximating the Voigt profile"
#       J. Appl. Cryst. (2000). 33, 1311-1316
#  @see doi:10.1107/S0021889800010219
#  @see https://doi.org/10.1107/S0021889800010219
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2016-06-15
class PseudoVoigt_pdf(Voigt_pdf) :
    """ Voigt function:
    Convolution of non-relativistic Breit-Wigner with gaussian resolution
    
    CPU-efficient Approximation of Voight profile
    -@see T. Ida, M. Ando and H. Toraya, 
    'Extended pseudo-Voigt function for approximating the Voigt profile'
    J. Appl. Cryst. (2000). 33, 1311-1316
    - see doi:10.1107/S0021889800010219
    - see https://doi.org/10.1107/S0021889800010219
    
    Parameters
    - mean  : location 
    - gamma : gamma for Breight-Wigner pole
    - sigma : resolution parameter for gaussian 
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   m0        = None ,
                   sigma     = None ,
                   gamma     = None ) :
        #
        ## initialize the base
        # 
        Voigt_pdf.__init__  ( self , name , xvar , m0 , sigma , gamma ) 

        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.PseudoVoigt (
            self.roo_name ( 'pvoigt_' ) , 
            "Pseudo Voigt %s" % self.name ,
            self.xvar   ,
            self.mean   ,
            self.gamma  ,
            self.sigma  )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'm0'        : self.mean  ,
            'sigma'     : self.sigma ,
            'gamma'     : self.gamma ,
            }
    
models.append ( PseudoVoigt_pdf )                          




# =============================================================================
## @class BreitWigner_pdf 
#  Relativistic Breit-Wigner function using Jackson's parameterization
#  J.D.Jackson,
#  "Remarks on the Phenomenological Analysis of Resonances",
#  In Nuovo Cimento, Vol. XXXIV, N.6
#  http://www.springerlink.com/content/q773737260425652/
#  - Blatt-Weisskopf forfactors  are also possible
#  @see Ostap::Models::BreitWigner 
#  @see Ostap::Math::BreitWigner
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class BreitWigner_pdf(PEAK) :
    """ Relativistic Breit-Wigner function using Jackson's parameterization
    J.D.Jackson, 'Remarks on the Phenomenological Analysis of Resonances',
    In Nuovo Cimento, Vol. XXXIV, N.6
    http://www.springerlink.com/content/q773737260425652/

    >>> bw    = Ostap.Math.BreitWigner( m_X , g_X , m_Jpsi , m_pipi , 0 )
    >>> breit = Models.BreitWigner_pdf ( 'BW'          ,
    ...                                  bw            ,
    ...                                  xvar  = mass  ,
    ...                                  m0    = m_X   ,
    ...                                  gamma = g_X   )
    
    Parameters:
    - m0          : location Breigt-Wigner function
    - gamma       : width of Breigt-Wigner function
    
    """
    def __init__ ( self               ,
                   name               ,
                   breitwigner        , ## Ostap::Math::BreitWeigner object
                   xvar               ,
                   m0          = None , 
                   gamma       = None ) :        
        #
        ## initialize the base
        #
        PEAK.__init__  ( self  , name  , xvar ,
                         mean        = m0                  ,
                         sigma       = gamma               ,
                         mean_name   = 'm0_%s'      % name ,
                         mean_title  = '#m_{0}(%s)' % name ,                         
                         sigma_name  = 'gamma_%s'   % name ,
                         sigma_title = '#Gamma(%s)' % name )
        
        bw = breitwigner
        assert isinstance ( bw , Ostap.Math.BW ), \
               'Invalid  type of the Breit-Wigner object: %s/%s' % ( bw   , type ( bw ) )
        #
        ## define the actual BW-shape using
        #      Ostap::Math::BreitWeigner object
        self.__breitwigner = breitwigner  ## Ostap::Math::BreitWeigner object

        ## create PDF 
        self.pdf = Ostap.Models.BreitWigner ( 
            self.roo_name ( 'rbw_' ) , 
            "Relativistic Breit-Wigner %s" % self.name ,
            self.xvar        ,
            self.m0          ,
            self.gamma       ,
            self.breitwigner )

        ## save the configuration
        self.config = {
            'name'        : self.name          ,
            'breitwigner' : self.breitwigner   ,
            'xvar'        : self.xvar          ,
            'm0'          : self.m0            ,
            'gamma'       : self.gamma         ,
            }

    ## ALIASES 
    m0    = PEAK.mean
    gamma = PEAK.sigma
    Gamma = PEAK.sigma
        
    @property
    def breitwigner ( self ) :
        """`breitwigner` : The Breit-Wigner' function  itself"""
        return self.__breitwigner

    # =========================================================================
    ## prepare Argand plot as <code>TGraph</code>
    #  @code
    #  bw = ...
    #  argand = bw.argand ( npx = 1000 )
    #  argand.draw ( 'al')  
    #  @endcode
    #  @see  TGraph 
    def argand ( self , x_min =  None , x_max = None , npx = 1000 ) :
        """ prepare Argand plot as `TGraph`
        >>> bw = ...
        >>> argand = bw.argand ( npx = 1000 )
        >>> argand.draw ( 'al')  
        """
        bw    = self.pdf.function()
        xmnmx = self.xminmax() 
        if x_min is None and xmnmx : x_min = xmnmx [ 0 ] 
        if x_max is None and xmnmx : x_max = xmnmx [ 1 ] 
        ## make Argand plot 
        return bw.argand ( xmin = x_min , xmax = x_max , npx = npx ) 
            
models.append ( BreitWigner_pdf )

# =============================================================================
## @class BWMC_pdf
#  Multi-channel version of Breit-Wigner function
#  @see Ostap::Models::BreitWignerMC
#  @see Ostap::Math::BreitWignerMC
#  @see Ostap::Math::Channel
#  @code
#  m_Kp  =  493.677 * MeV
#  m_Kz  =  497.614 * MeV
#  m_phi = 1019.462 * MeV
#  g_phi =    4.249 * MeV 
#  br_pm = 0.492
#  br_00 = 0.340 
#  br_xx = 1 - br_pm - br_00
#  ## define three  channels 
#  ch_pm = Ostap.Math.Channel ( g_phi * br_pm , m_Kp , m_Kp , 1 , 1 )
#  ch_00 = Ostap.Math.Channel ( g_phi * br_00 , m_Kz , m_Kz , 1 , 1 )
#  ch_xx = Ostap.Math.Channel ( g_phi * br_xx , 0    , 0    )
#  ## define the Breit-wigner function
#  bw    = Ostap.Math.BreitWignerMC ( m_phi , ch_pm , ch_00 , ch_xx )
#
#  mKK   = ROOT.RooRealVar ( ... ) 
#  pdf = BWMC_pdf ( 'BW' , breitwigner = bw ,
#                    xvar = mKK , mean = m_phi , gamma = g_phi ,
#                    fractions = [ br_pm , br_00 , br_xx ] ) 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2018-11-25
class BWMC_pdf(PEAK) :
    """ Multi-channel version of Relativistic Breit-Wigner function

    >>> m_Kp  =  493.677 * MeV    ## mass of K+ 
    >>> m_Kz  =  497.614 * MeV    ## mass of K0
    >>> m_phi = 1019.462 * MeV    ## mass of phi(1020)
    >>> g_phi =    4.249 * MeV    ## width of phi(1020)
    >>> br_pm = 0.492             ## Br ( phi -> K+K-) 
    >>> br_00 = 0.340             ## Br ( phi -> K0s K0L  )
    >>> br_xx = 1 - br_pm - br_00 ## Br ( phi ->  anything else )
    
    Define three decay channels:
    
    >>> ch_pm = Ostap.Math.Channel ( g_phi * br_pm , m_Kp , m_Kp , 1 , 1 )  ## K+ K- 
    >>> ch_00 = Ostap.Math.Channel ( g_phi * br_00 , m_Kz , m_Kz , 1 , 1 )  ## K0L K0L
    >>> ch_xx = Ostap.Math.Channel ( g_phi * br_xx , 0    , 0    )          ## light stuff 
    
    Define the Breit-Wigner function:
    
    >>> bw    = Ostap.Math.BreitWignerMC ( m_phi , ch_pm , ch_00 , ch_xx )

    >>> mKK   = ROOT.RooRealVar ( ... ) 
    >>> pdf = BWMC_pdf ( 'BW' , breitwigner = bw ,
    ...                  xvar = mKK , mean = m_phi , gamma = g_phi ,
    ...                  fractions = [ br_pm , br_00 , br_xx ] ) 
        
    Parameters:
    - xvar        : fitting variable/observable  
    - mean        : location of Breit-Wigner pole
    - gamma       : total width of Breit-Wigner pole
    - widths      : partial width for the channels  (mutually exclusive with gamma and fractions)
    - fractions   : branching fractions 
    
    """
    def __init__ ( self               ,
                   name               ,
                   breitwigner        , ## Ostap::Math::BreitWignerMC object
                   xvar               ,
                   m0          = None , 
                   gamma       = None ,
                   widths      = []   ,
                   fractions   = []   ) : 

        ## correct type of Breit-Wigner function?
        bw = breitwigner 
        assert isinstance ( bw , Ostap.Math.BreitWignerMC ), \
               'Invalid  type of the Breit-Wigner object: %s/%s' % ( bw   , type ( bw ) )

        ## number of channels 
        nc      = bw.nChannels    ()
        
        self.__brfrs  = ROOT.RooArgList()
        self.__widths = ROOT.RooArgList()
        
        case         = None
        self.__trash = []  ## keep the trash

        ## Valid  cases: 
        if   widths    and nc == len ( width     ) and gamma is None and not fractions : case = 1
        elif fractions and nc == len ( fractions )                   and not widths    : case = 2
        else : raise TypeError ('Gamma/widths/fraction mismatch!')

        ## partial  widths are specified:
        if 1 == case :
            
            for i in range ( nc ) : 
                gi = widths[i] 
                gg = self.make_var ( gi ,
                                     'gamma_%d_%s'     % ( i+1 , name ) ,
                                     '#Gamma_{%d}(%s)' % ( i+1 , name ) ,
                                     None )
                
                self.__trash.append ( gg )
                self.widths.add     ( gg )

            self.__trash.append ( self.widths ) 
            ## construct the total gamma
            formula = '%s ' % self.widths[0].GetName()
            for i in range ( 1 , nc ) : formula += ' + %s' % self.widths[i].GetName()  
            ## for g in self.widths[1:] :

            ## create gamma 
            gamma = Ostap.FormulaVar ( 'gamma_%s'    % name ,
                                       '#Gamma_(%s)' % name , formula , self.widths )
            self.__trash.append ( gamma )
            
        # =====================================================================
        ## initialize the base 
        # =====================================================================
        PEAK.__init__ ( self , name , xvar                ,
                        mean        =  m0                 ,
                        siga        =  gamma              ,
                        mean_name   = 'm0_%s'      % name ,
                        mean_title  = '#m_{0}(%s)' % name ,
                        sigma_name  = 'gamma_%s'   % name ,
                        sigma_title = '#Gamma(%s)' % name )
                
        ## create branching fractions 
        if 1 == case :
            
            for i in range ( nc ) :
                
                gi  = self.widths[i]
                lst = ROOT.RooArgList ( self.gamma ,  gi )
                self.__formulas_lists.append ( lst ) 
                ## br  = ROOT.RooFormulaVar ( 'brfr_%d_%s'  % ( i + 1 , name ) ,
                ##                            'Br_{%d}(%s)' % ( i + 1 , name ) ,
                ##                            '%s / %s'     % ( gi.GetName() , self.gamma.GetName() ) , lst )
                br  =  Ostap.MoreRooFit.Division ( 'brfr_%d_%s'  % ( i + 1 , name ) ,
                                                   'Br_{%d}(%s)' % ( i + 1 , name ) , self.gamma , gi )
                self.brfrs.add      ( br )
                self.__trash.append ( br )
                
        ##  branching fractions are specified 
        elif 2 == case : 
            
            for i in range ( nc ) :
                
                bi = fractions [i] 
                br = self.make_var       ( bi ,
                                           'brfr_%d_%s'  % ( i+1 , name ) ,
                                           'Br_{%d}(%s)' % ( i+1 , name ) ,None ) 
                self.brfrs.add       ( br )    
                self.__trash.append  ( br ) 
                
                ls = ROOT.RooArgList ( self.gamma ,  br )
                self.__trash.append  ( ls ) 

                ## gg  = ROOT.RooFormulaVar ( 'gamma_%d_%s'     % ( i + 1 , name ) ,
                ##                           '#Gamma_{%d}(%s)' % ( i + 1 , name ) ,
                ##                           '%s * %s'         % ( br.GetName() , self.gamma.GetName() ) , ls )
                
                gg  = Ostap.MoreRooFit.Product ( 'gamma_%d_%s'     % ( i + 1 , name ) ,
                                                 '#Gamma_{%d}(%s)' % ( i + 1 , name ) , self.gamma , br ) 
                self.widths.add      ( gg )
                self.__trash.append  ( gg ) 
            

        self.__gammas    = tuple ( [ i for i in self.widths ] )
        self.__fractions = tuple ( [ i for i in self.brfrs  ] )
                       
        #
        ## define the actual BW-shape using
        #      Ostap::Math::BreitWeignerMC object
        #
        self.__breitwigner = breitwigner  ## Ostap::Math::BreitWeignerMC object
        
        ## create PDF 
        self.pdf = Ostap.Models.BreitWignerMC ( 
            self.roo_name ( 'rbwmc_' ) , 
            "Multi-channel relativistic Breit-Wigner %s" % self.name ,
            self.xvar        ,
            self.mean        ,
            self.widths      , 
            self.breitwigner )
            
        ## save the configuration
        self.config = {
            'name'        : self.name          ,
            'breitwigner' : self.breitwigner   ,
            'xvar'        : self.xvar          ,
            'mean'        : self.mean          }
        
        if   1 == case : self.config.update (  { 'widths'    : self.gammas    ,
                                                 'fractions' : ()             } )
        elif 2 == case : self.config.update (  { 'fractions' : self.fractions ,
                                                 'gamma'     : self.gamma     ,
                                                 'widths'    : ()             } )

    m0    = PEAK.mean
    gamma = PEAK.sigma
    Gamma = PEAK.sigma

    @property
    def breitwigner ( self ) :
        """The Breit-Wigner function  itself"""
        return self.__breitwigner

    @property
    def widths  ( self ) :
        """'widths'  : partial widths for different decay channels"""
        return self.__widths

    @property
    def brfrs   ( self ) :
        """'brfrs'  : branching fractions for different decay channels"""
        return self.__brfrs 

    @property
    def gammas ( self ) :
        """'gammas'  : partial widths for different decay channels"""
        return self.__gammas  

    @property
    def fractions ( self ) :
        """'fractions'  : branching fractions for different decay channels"""
        return self.__fractions

    # =========================================================================
    ## prepare Argand plot as <code>TGraph</code>
    #  @code
    #  bw = ...
    #  argand = bw.argand ( npx = 1000 )
    #  argand.draw ( 'al')  
    #  @endcode
    #  @see  TGraph 
    def argand ( self , x_min =  None , x_max = None , npx = 1000 ) :
        """ prepare Argand plot as `TGraph`
        >>> bw = ...
        >>> argand = bw.argand ( npx = 1000 )
        >>> argand.draw ( 'al')  
        """
        bw = self.pdf.function()
        xmnmx = self.xminmax() 
        if x_min is None and xmnmx : x_min = xmnmx [ 0 ] 
        if x_max is None and xmnmx : x_max = xmnmx [ 1 ] 
        ## make Argand plot 
        return bw.argand ( xmin = x_min , xmax = x_max , npx = npx ) 
            

models.append ( BWMC_pdf )

# =============================================================================
## @class BWI_pdf
#  (relativistic) Breit-Wigner function + some interference
#   Breit-Wigner with some embedded interference: 
#   \f[ f(x) = \left| \upalpha b(x) + A(x)_{\mathrm{BW}} \right|^2 \f], 
#   where \f$b(x)\f$ - any smooth function and 
#   \f$ A(x)_{\mathrm{BW}} \f$ is Breit-Wigner amplitude 
#  @see Ostap.Models.BWI
class BWI_pdf (BreitWigner_pdf) :
    """ (Relativistic) Breit-Wigner function + some interference
    Breit-Wigner with some embedded interference
    - see Ostap.Models.BWI
    """
    def __init__ ( self             ,
                   name             ,
                   breitwigner      ,   ## Breit-Wigner function 
                   xvar             ,   ## observable 
                   m0        = None ,   ## m0-parameter
                   gamma     = None ,   ## gamma-parameter(s)  
                   magnitude = -1   ,   ## background magnitude 
                   phase     = -1   ,   ## background phase 
                   scale1    = None ,   ## background magnitide scale (if/when needed)
                   scale2    = None ) : ## background phase scale (if/when needed)
        
        ## initialize the base 
        BreitWigner_pdf.__init__ ( self , name , breitwigner , xvar , m0 , gamma )

        ## save the constructed PDF  (we'll need it)
        self.__bw = self.pdf
        
        self.__scale1 = ROOT.RooFit.RooConst ( 1 )
        self.__scale2 = ROOT.RooFit.RooConst ( 1 )
        # =====================================================================
        ## take care on background and scale 
        # =====================================================================

        if isinstance ( magnitude , PDF1 ) and self.xvar is magnitude.xvar :
            
            ## almost ideal case , but the scale factor is needed 
            self.__magnitude     = magnitude
            self.__scale1        = self.make_var ( scale1 ,
                                                   "scale1_%s"  % name ,
                                                   "scale1(%s)" % name ,
                                                   False , 1 , 1e-6 , 1e+6 )
            self.__magnitude_tot = self.magnitude.as_FUN() * self.scale1 
            
        elif   isinstance ( magnitude , FUN1 ) and self.xvar is magnitude.xvar :
            
            ## ideal case
            self.__magnitude     =      magnitude
            self.__magnitude_tot = self.magnitude
            
        elif isinstance ( magnitude  , ROOT.RooRealVar ) :
            
            ## simple case: background is a kind of a simple variable 
            self.__magnitude     = Fun1D ( magnitude , xvar = self.xvar , name = self.new_name ( 'magnitude' ) ) 
            self.__magnitude_tot = self.magnitude

        elif isinstance ( magnitude , ROOT.RooAbsPdf ) :
            
            ## simple case: background is RooAbsPdf, scale factor is needed
            self.__magnitude    = Fun1D ( magnitude , xvar = self.xvar , name = self.new_name ( 'bkg' ) ) 
            self.__scale1       = self.make_var ( scale1 ,
                                                  "scale1_%s"  % name ,
                                                  "scale1(%s)" % name ,
                                                  False , 1 , 1e-6 , 1e+6 )
            self.__magnitude_tot = self.magnitude * self.scale1 

        elif isinstance ( magnitude , ROOT.RooAbsReal ) :
            
            ## simple case: background is RooAbsReal, scale factor is NOT needed 
            self.__magnitude     = Fun1D ( magnitude , xvar = self.xvar , name = self.new_name ( 'bkg' ) ) 
            self.__magnitude_tot = self.magnitude

        else :
            
            ## use the predefined shapes (PDFs), scale factor is needed  
            from ostap.fitting.background import make_bkg as MKB
            self.__magnitude = MKB ( magnitude ,  name = 'B4'+ self.name , xvar  = self.xvar )
            self.__scale1    = self.make_var ( scale1 ,
                                               "scale1_%s"  % name ,
                                               "scale1(%s)" % name ,
                                               False , 1 , 1e-6 , 1e+6 )            
            self.__magnitude_tot = self.magnitude.as_FUN() * self.scale1 
            
        # =====================================================================
        ## now take care on the phase
        # =====================================================================

        if   isinstance ( phase , PDF1 ) and self.xvar is phase.xvar :
            
            ## almost ideal case , but the scale factor is needed 
            self.__phase     = phase 
            self.__scale2    = self.make_var ( scale2 ,
                                               "scale2_%s"  % name ,
                                              "scale3(%s)" % name ,
                                              False , 1 , 1e-6 , 1e+6 )
            self.__phase_tot = self.phase.as_FUN () * self.scale2
            
        elif isinstance ( phase , FUN1 ) and self.xvar is phase.xvar :
            
            ## the ideal case 
            self.__phase     =      phase 
            self.__phase_tot = self.phase 
            
        elif isinstance ( phase , ROOT.RooRealVar ) :
            
            ## simple case: phase is a kind of a simple constant/variable
            self.__phase     = Fun1D ( phase , name = self.new_name ( 'phase' ) , xvar = self.xvar ) 
            self.__phase_tot = self.phase 
            
        elif isinstance ( phase , ROOT.RooAbsPdf ) :
            
            ## simple case: phase is RooAbsPdf, scale factor is needed
            self.__phase     = Fun1D ( phase , name = self.new_name ( 'phase' ) , xvar = self.xvar ) 
            self.__scale2    = self.make_var ( scale2 ,
                                               "scale2_%s"  % name ,
                                               "scale2(%s)" % name ,
                                               False , 1 , 1e-6 , 1e+6 )
            self.__phase_tot = self.phase * self.scale2 
            
        elif isinstance ( phase , ROOT.RooAbsReal ) :
            
            ## simple case: phase is a kind of a simple constant/variable
            self.__phase     = Fun1D ( phase , name = self.new_name ( 'phase' ) ) 
            self.__phase_tot = self.phase 
            
        else :  

            ## use the predefined shapes (PDFs), scale factor is needed  
            from ostap.fitting.background import make_bkg as MKB
            self.__phase     = MKB ( phase ,  name = 'P4'+ self.name , xvar  = self.xvar )
            self.__scale2    = self.make_var ( scale2 ,
                                               "scale2_%s"  % name ,
                                               "scale2(%s)" % name ,
                                               False , 1 , 1e-6 , 1e+6 )
            self.__phase_tot = self.phase.as_FUN () * self.scale2

        ## finally create PDF
        self.pdf = Ostap.Models.BWI (
            self.roo_name ( 'bwi_' ) ,
            "Breit-Wigner with interference %s" % self.name  ,
            self.bw            ,
            self.magnitude.fun ,
            self.phase    .fun ,
            self.scale1        ,
            self.scale2        )

        ## save configuration
        self.config = {
            'name'        : self.name        ,
            'xvar'        : self.xvar        , 
            'breitwigner' : self.breitwigner , 
            'm0'          : self.m0          ,
            'gamma'       : self.gamma       ,
            'magnitude'   : self.magnitude   ,
            'phase'       : self.phase       ,
            'scale1'      : self.scale1      , 
            'scale2'      : self.scale2      } 

    @property
    def bw     ( self ) :
        """'bw' : original BreitWigner PDF (no interference)"""
        return self.__bw
    
    @property
    def magnitude ( self ) :
        """The magnitude of the coherent background"""
        return self.__magnitude

    @property
    def phase ( self ) :
        """The phase of the coherent background"""
        return self.__phase
    
    @property
    def magnitude_total ( self ) :
        """The total magnitude of the coherent background : magnitide * scale1 """
        return self.__magnitude_tot 

    @property
    def phase_total  ( self ) :
        """The total phase of the coherent background: phase * scale2 """
        return self.__phase_tot 
    
    @property
    def scale1 ( self ) :
        """The scale-factor for the coherent background"""
        return self.__scale1
    @scale1.setter
    def scale1 ( self , value ) :
        self.set_value ( self.__scale1 , value ) 

    @property
    def scale2 ( self ) :
        """The scale-factor for the coherent background"""
        return self.__scale2
    @scale2.setter
    def scale2 ( self , value ) :
        self.set_value ( self.__scale2 , value ) 

    # =========================================================================
    ## Get the "signal" component 
    def cmp_signal ( self , x )  :
        """get the 'signal' component"""
        
        with SETVAR ( self.xvar ) :
            self.xvar.setVal ( x )
            amp = self.pdf.bw_amplitude ()
            return self.pdf.breit_wigner ( x , amp )

    # =========================================================================
    ## Get the "background" component 
    def cmp_bkg    ( self , x )  :
        """Get the 'signal' component"""
        
        with SETVAR ( self.xvar ) :
            self.xvar.setVal ( x )
            amp = self.pdf.amplitude () - self.pdf.bw_amplitude () 
            return self.pdf.breit_wigner ( x , amp )
        
    # =========================================================================
    ## Get the "total" component 
    def cmp_total ( self , x )  :
        """Get the 'total' component"""
        
        with SETVAR ( self.xvar ) :
            self.xvar.setVal ( x )
            amp = self.pdf.amplitude ()
            return self.pdf.breit_wigner ( x , amp )

    # =========================================================================
    ## Get the "interference" component 
    def cmp_interference ( self , x )  :
        """Get the 'interference' component"""

        s = self.cmp_signal ( x )
        b = self.cmp_bkg    ( x )
        t = self.cmp_total  ( x )

        return t - ( s + b ) 
    
    # =========================================================================
    ## Get the signal, background & interference fit fractions
    #  @code
    #  pdf = ...
    #  pdf.fitTo  ( ... )
    #  sF , bF , iF = pdf.fit_fractions () 
    #  @endcode 
    def ffs  ( self ) :
        """ Get the signal, background & interference fit fractions 
        >>> pdf = ...
        >>> pdf.fitTo  ( ... )
        >>> sF , bF, iF  = pdf.fit_fractions () 
        """

        self.pdf.setPars()
        
        xmn, xmx = self.xminmax() 
        
        import ostap.math.integral as I
        
        iS = I.integral ( self.cmp_signal       , xmin = xmn , xmax = xmx )
        iB = I.integral ( self.cmp_bkg          , xmin = xmn , xmax = xmx )
        iI = I.integral ( self.cmp_interference , xmin = xmn , xmax = xmx )
        iT = I.integral ( self.cmp_total        , xmin = xmn , xmax = xmx )
                
        return iS/iT, iB/iT, iI/iT  
            
# =============================================================================
## @class BWPS
#  Breit-Wigner function modulated with extra phase-space and polynomial factors
#  @see Ostap::Models::BWPS
#  @see Ostap::Math::BWPS
class BWPS_pdf(BreitWigner_pdf,Phases) :
    """Breit-Wigner function modulated with extra phase-space and polynomial factors
    - see Ostap.Models.BWPS
    - see Ostap.Math.BWPS
    """    
    def __init__ ( self             ,
                   name             ,
                   breitwigner      , ## Ostap::Math::BWPS object
                   xvar             ,
                   m0        = None ,
                   gamma     = None ,
                   the_phis  = None ) :

        if   isinstance ( breitwigner , Ostap.Math.BWPS ) : pass
        elif isinstance ( breitwigner , tuple           ) :
            breitwigner = Ostap.Math.BWPS  ( *breitwigner )
        else :
            raise ArgumentError("BWPS_pdf: Invalid type of `breitwigner`") 
        
        ## initialize the base classes 
        BreitWigner_pdf.__init__  ( self ,
                                    name ,
                                    breitwigner = breitwigner.breit_wigner () , 
                                    xvar        = xvar   ,
                                    m0          = m0     ,
                                    gamma       = gamma  )
        
        Phases.__init__ ( self , breitwigner.npars () , the_phis ) 

        ## make "original" BW-pdf 
        self.__bw_pdf = BreitWigner_pdf ( name        = self.name + '_orig' ,
                                          breitwigner = self.breitwigner    ,
                                          xvar        = self.xvar           ,
                                          m0          = self.m0             ,
                                          gamma       = self.gamma          )
        self.__bwps = breitwigner
        
        ## finally create PDF        
        self.pdf = Ostap.Models.BWPS (
            self.roo_name ( 'bwps_' ) ,
            "Breit-Wigner with phase space %s" % self.name  ,
            self.xvar         ,
            self.m0           ,
            self.gamma        ,
            self.phi_list     ,
            self.bwps         ) 
            
        ## save configuration
        self.config = {
            'name'        : self.name  ,
            'xvar'        : self.xvar  , 
            'breitwigner' : ( self.bwps.breit_wigner () ,
                              self.bwps.phase_space  () ,
                              self.bwps.use_rho      () ,
                              self.bwps.use_N2       () ) ,                               
            'm0'          : self.mean  ,
            'gamma'       : self.gamma ,
            'the_phis'    : self.phis  }
        
    @property
    def bw_pdf  ( self ) :
        """'bw_pdf' : 'original' Breit-Wigner pdf (no additional phase space  factors)"""
        return self.__bw_pdf
    
    @property
    def bwps ( self ) :
        """The Breit-Wigner function (BWPS) itself"""
        return self.__bwps

    
models.append ( BWPS_pdf )

# =============================================================================
## @class BW3L_pdf
#  Breit-Wigner function modulated with \f$ p^{2L+1}\f$ factor
#   - it can approximate the mass distrbition from 3-body decays
# e.g.  \f$ \eta^{\prime)  \rigtharrow \left(\rho^0 
#               \rigtharrow \pi^+ \pi^-\right)\gamma \f$~decays
#    or similar  configurations  
# 
# \f[ f(x) \equiv F_{\mathrm{BW}}(x) p(x|M_0,m_3)^{2L+1} \f]
#  - \f$ p(x|M,m_3) \f$ is a momentumm of the 3rd particle, \f$P_3\f$ 
#       in the \f$ P \rightarrow \left( P_{\mathrm{BW}} \rightharrow 
#      P_1 P_2 \right) P_3 \f$ decay chain
#  - \f$ M \f$ is a (fixed) mass of "mother" particle \f$P\f$
#  - \f$ m_1\f$ is a (fixed) mass of 1st particle \f$P_1\f$
#  - \f$ m_2\f$ is a (fixed) mass of 2nd particle \f$P_2\f$
#  - \f$ m_3\f$ is a (fixed) mass of 3rd particle \f$P_3\f$
#  - \f$ x \equiv m_{23} \f$ is a mass intermediate Breit-Wigner particle \f$P_{\mathrm{BW}}\f$
#  - \f$ L \f$  is an orbital momentum between \f$ P_{\mathrm{BW}}\f$ and \f$ P_3\f$
# 
#  It is assumed that  \f$ m_1\f$  and \f$ m_2\f$ parameters 
#  are in agreement with the Breit-Wigner definition 
#  @see Ostap::Models::BW3L
#  @see Ostap::Math::BW3L
class BW3L_pdf(BreitWigner_pdf) :
    """ Breit-Wigner function modulated with  p^{2L+1} factor
    - it can approximate the mass distrbition from 3-body decays
    e.g. ( eta'  ->  ( rho0 -> pi+ pi- ) gamma ) decay or similar  configurations
    
    - see Ostap.Models.BW3L
    - see Ostap.Math.BW3L
    """
    
    def __init__ ( self             ,
                   name             ,
                   breitwigner      , ## Ostap::Math::BW3L object
                   xvar             ,
                   m0        = None ,
                   gamma     = None ) : 

        if   isinstance ( breitwigner , Ostap.Math.BW3L ) : pass
        elif isinstance ( breitwigner , tuple ) :
            breitwigner = Ostap.Math.BW3L  ( *breitwigner )
        else :
            raise ArgumentError("BW3L_pdf: Invalidd type of breitwigner") 
        
        ## initialize the base classes 
        BreitWigner_pdf.__init__  ( self ,
                                    name ,
                                    breitwigner = breitwigner.breit_wigner () , 
                                    xvar        = xvar   ,
                                    m0          = m0     ,
                                    gamma       = gamma  )
        
        ## make "original" BW-pdf 
        self.__bw_pdf = BreitWigner_pdf ( name        = self.name + '_orig' ,
                                          breitwigner = self.breitwigner    ,
                                          xvar        = self.xvar           ,
                                          m0          = self.m0             ,
                                          gamma       = self.gamma          )
        self.__bw3l = breitwigner
        
        ## finally create PDF        
        self.pdf = Ostap.Models.BW3L (
            self.roo_name ( 'bw3l_' ) ,
            "Breit-Wigner form 3-body decay  %s" % self.name  ,
            self.xvar         ,
            self.m0           ,
            self.gamma        ,
            self.bw3l         ) 
            
        ## save configuration
        self.config = {
            'name'        : self.name  ,
            'xvar'        : self.xvar  , 
            'breitwigner' : ( self.bw3l.breit_wigner () ,
                              self.bw3l.M  () ,
                              self.bw3l.m1 () ,
                              self.bw3l.m2 () ,
                              self.bw3l.m3 () ,
                              self.bw3l.L  () ) , 
            'm0'          : self.mean  ,
            'gamma'       : self.gamma }
        
    @property
    def bw_pdf  ( self ) :
        """'bw_pdf' : 'original' Breit-Wigner pdf (no additional factors)"""
        return self.__bw_pdf
    
    @property
    def bw3l ( self ) :
        """The Breit-Wigner function (BW3L) itself"""
        return self.__bw3l

    
models.append ( BW3L_pdf )

# =============================================================================
## @class Flatte_pdf
#  Flatte function
#  S.M.Flatte, "Coupled-channel analysis of the \f$\pi\eta\f$ 
#    and \f$K\bar{K}\f$ systems near \f$K\bar{K}\f$ threshold  
#    Phys. Lett. B63, 224 (1976)
#  Well suitable for \f$f_0(980)\rightarrow \pi^+ \pi^-\f$
#  @see http://www.sciencedirect.com/science/article/pii/0370269376906547
#  @see Ostap::Models::Flatte
#  @see Ostap::Math::Flatte
#  @see Ostap::Math::Flatte2
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-01-18
class Flatte_pdf(PEAKMEAN) :
    """Flatte function:
    S.M.Flatte, 'Coupled-channel analysis of the (pi eta)
    and (KbarK) systems near (KbarK) threshold' 
    Phys. Lett. B63, 224 (1976
    http://www.sciencedirect.com/science/article/pii/0370269376906547

    Typical case:    f0(980) -> pi+ pi- & K+ K- shapes 
    """
    def __init__ ( self              ,
                   name              ,
                   flatte            ,    ## Ostap::Math::Flatte/Flatte2
                   xvar              ,
                   m0       = None   ,    ## the pole 
                   m0g1     = None   ,    ## m0*g1 
                   g2og1    = None   ,    ## g2/g1 
                   g1       = None   ,    ## g1 
                   g2       = None   ,    ## g2 
                   gamma0   = None   ) :  ## gamma0 

        assert isinstance ( flatte , Ostap.Math.Flatte ), \
               'Invalid type for flatte %s' %  type ( flatte )

        ## initialize the base
        with CheckMean ( False ) :
            # for Flatte-function m0 can be outside the interesting interval 
            PEAKMEAN.__init__  ( self , name , xvar ,
                                 mean       = m0  ,
                                 mean_name  = 'm0_%s'      % name ,
                                 mean_title = '#m_{0}(%s)' % name )

        self.__my_case = 0 
        if   all_args ( self.mean , m0g1 , g2og1 ) : self.__my_case = 1 
        elif all_args ( self.mean ,   g1 , g2    ) : self.__my_case = 2
            
        assert self.case in  ( 1 , 2 ), 'Invalid combination of (m0g1,g2og1:g1,g2)arguments!'
        
        self.__flatte = flatte
            
        self.__gamma0 = self.make_var  ( gamma0                  ,
                                         'gamma0_%s'      % name ,
                                         '#Gamma_{0}(%s)' % name ,
                                         True    , 0 , 10  )
        if  1 == self.case : 
            
            self.__m0g1 = self.make_var  ( m0g1                     ,
                                           'm0g1_%s'          % name ,
                                           '#m_{0}#g_{1}(%s)' % name ,
                                           None  )
            
            self.__g2og1 = self.make_var ( g2og1 ,
                                           'g2og1_%s'          % name ,
                                           '#g_{2}/#g_{1}(%s)' % name ,
                                           None , 0.001 , 200  ) 

            self.__g1 = self.vars_divide   ( self.m0g1  , self.m0 , name = 'g1_%s' % name , title = "g_1(%s)" % name )
            self.__g2 = self.vars_multiply ( self.g2og1 , self.g1 , name = 'g2_%s' % name , title = "g_2(%s)" % name )
            
        elif 2 == self.case :
            
            self.__g1 =  self.make_var  ( g1                 ,
                                          'g1_%s'     % name ,
                                          'g_{1}(%s)' % name ,
                                          None )
            self.__g2 =  self.make_var  ( g2                 ,
                                          'g2_%s'     % name ,
                                          'g_{2}(%s)' % name ,
                                          None )
            
            self.__m0g1  = self.vars_multiply ( self.m0 , self.g1 , name = 'm0g1_%s'  % name , title = "m_0g_1(%s)"  % name )
            self.__g2og1 = self.vars_divide   ( self.g2 , self.g1 , name = 'g2og1_%s' % name , title = "g_2/g_1(%s)" % name )
                
        ## create PDF 
        self.pdf = Ostap.Models.Flatte (
            self.roo_name ( 'flatte_' ) ,
            "Flatte %s" % self.name  ,
            self.xvar    ,
            self.m0      ,
            self.g1      ,
            self.g2      ,
            self.gamma0  ,
            self.flatte  )

        ## save the configuration
        cnf = {
            'name'        : self.name    ,
            'flatte'      : self.flatte  ,
            'xvar'        : self.xvar    ,
            'm0'          : self.m0      ,
            'gamma0'      : self.gamma0  ,
            }
        
        if   1 == self.case : cnf.update ( { 'm0g1' : self.m0g1 , 'g2og1' : self.g2og1 } )
        elif 2 == self.case : cnf.update ( { 'g1'   : self.g1   , 'g2'    : self.g2    } )

        self.config = cnf
        
    @property
    def case  (self ) :
        """'case' : How the input argument are  specified: 1: (m0g1,g2og1) vs 2: (g1,g2) """
        return self.__my_case

    ## ALIAS 
    m0 = PEAK.mean

    @property
    def m0g1 ( self ) :
        """'m0*g1'-parameter for Flatte-function"""
        return self.__m0g1
    @m0g1.setter
    def m0g1 ( self, value ) :
        assert gasattr ( self.__m0g1 , 'setVal' ),"'m0g1'-parameter can't be set!"
        value = float ( value )
        self.__m0g1.setVal ( value ) 

    @property
    def g2og1 ( self ) :
        """'g2/g1'-parameter for Flatte-function"""
        return self.__g2og1
    @g2og1.setter
    def g2og1 ( self, value ) :
        assert hasattr ( self.__g2og1 , 'setVal'),"'g2og1'-parameter can't be set!"        
        value = float ( value )
        assert 0 < value, "'g2/g1'-parameter for Flatte-function must be positive"
        self.__g2og1.setVal ( value )

    @property
    def g1 ( self ) :
        "'g1'-parameter for Flatte-function"
        return self.__g1
    @g1.setter
    def g1 ( self , value ) :
        assert hasattr ( self.__g1 , 'setVal' ),"'g1'-parameter can't be set!"
        value = float ( value )
        self.__g1.setVal ( value ) 

    @property
    def g2 ( self ) :
        "'g2'-parameter for Flatte-function"
        return self.__g2
    @g2.setter
    def g2 ( self , value ) :
        assert hasattr ( self.__g2 , 'setVal' ),"'g2'-parameter can't be set!"
        value = float ( value )
        self.__g2.setVal ( value ) 

    @property
    def gamma0 ( self ) :
        "'gamma0'-parameter for Flatte-function"
        return self.__gamma0
    @gamma0.setter
    def gamma0 ( self , value ) :
        value = float ( value )
        self.__gamma0.setVal ( value ) 

    @property
    def flatte ( self ) :
        """The Flatte function itself"""
        return self.__flatte

    # =========================================================================
    ## prepare Argand plot as <code>TGraph</code>
    #  @code
    #  bw = ...
    #  argand = bw.argand ( npx = 1000 )
    #  argand.draw ( 'al')  
    #  @endcode
    #  @see  TGraph 
    def argand ( self , x_min =  None , x_max = None , npx = 1000 ) :
        """ prepare Argand plot as `TGraph`
        >>> bw = ...
        >>> argand = bw.argand ( npx = 1000 )
        >>> argand.draw ( 'al')  
        """
        bw    = self.pdf.function()
        xmnmx = self.xminmax() 
        if x_min is None and xmnmx : x_min = xmnmx [ 0 ] 
        if x_max is None and xmnmx : x_max = xmnmx [ 1 ] 
        ## make Argand plot 
        return bw.argand ( xmin = x_min , xmax = x_max , npx = npx ) 
            

models.append ( Flatte_pdf )                          

# ============================================================================
class FlattePS_pdf(Flatte_pdf,Phases) :
    """Flatte function:
    S.M.Flatte, 'Coupled-channel analysis of the (pi eta)
    and (KbarK) systems near (KbarK) threshold' 
    Phys. Lett. B63, 224 (1976
    http://www.sciencedirect.com/science/article/pii/0370269376906547

    Typical case:    f0(980) -> pi+ pi- & K+ K- shapes 
    """
    def __init__ ( self              ,
                   name              ,
                   flatte            ,   ## Ostap::Math::BWPS 
                   xvar              ,
                   m0       = None   ,   ## the pole 
                   m0g1     = None   ,   ## m0*g1 
                   g2og1    = None   ,   ## g2/g1 
                   g1       = None   ,   ## g1 
                   g2       = None   ,   ## g2 
                   gamma0   = None   ,   ## gamma0 
                   the_phis = None   ) : ##
        
        if   isinstance ( flatte , Ostap.Math.BWPS ) : pass
        elif isinstance ( flatte , tuple ) :
            flatte = Ostap.Math.BWPS  ( *flatte )
        else :
            raise ArgumentError("FlattePS_pdf: Invalidd type of flatte") 
        

        assert isinstance ( flatte, Ostap.Math.BWPS ),\
               'Invalid type for breitwigner %s' %  type ( flatte )
        
        ## initialize the base classes 
        Flatte_pdf.__init__  ( self ,
                               name ,
                               flatte = flatte.breit_wigner () ,
                               xvar   = xvar     , 
                               m0     = m0       ,
                               m0g1   = m0g1     ,
                               g2og1  = g2og1    ,
                               g1     = g1       ,
                               g2     = g2       ,
                               gamma0 = gamma0   )
        
        Phases.__init__ ( self , flatte.npars () , the_phis )

        ## store the "original"" Flatte PDF 
        if  1 ==  self.case : 
            self.__flatte_pdf = Flatte_pdf ( name   = self.name + "_orig" ,
                                             flatte = self.flatte ,  
                                             xvar   = self.xvar   ,
                                             m0     = self.m0     ,
                                             m0g1   = self.m0g1   ,
                                             g2og1  = self.g2og1  ,
                                             gamma0 = self.gamma0 )
        elif 2 ==  self.case : 
            self.__flatte_pdf = Flatte_pdf ( name   = self.name + "_orig" ,
                                             flatte = self.flatte ,  
                                             xvar   = self.xvar   ,
                                             m0     = self.m0     ,
                                             g1     = self.g1     ,
                                             g2     = self.g2     ,
                                             gamma0 = self.gamma0 )
            
        self.__bwps   = flatte
  
        self.__g_list = ROOT.RooArgList ( self.g1 , self.g2 , self.gamma0 )
        
        ## finally create PDF
        self.pdf = Ostap.Models.BWPS (
            self.roo_name ( 'flatteps_' ) ,
            "Flatte with phase space %s" % self.name  ,
            self.xvar         ,
            self.m0           ,
            self.g_list       ,
            self.phi_list     ,
            self.bwps         ) 
        
        ## save the configuration
        cnf = {
            'name'        : self.name    ,
            'flatte'      : ( self.bwps.breit_wigner () ,
                              self.bwps.phase_space  () ,
                              self.bwps.use_rho      () ,
                              self.bwps.use_N2       () ) ,                               
            'xvar'        : self.xvar    ,
            'm0'          : self.m0      ,
            'gamma0'      : self.gamma0  ,
            'the_phis'    : self.phis    , 
            }
        
        if   1 == self.case : cnf.update ( { 'm0g1' : self.m0g1 , 'g2og1' : self.g2og1 } )
        elif 2 == self.case : cnf.update ( { 'g1'   : self.g1   , 'g2'    : self.g2    } )

        self.config = cnf

    @property
    def flatte_pdf  ( self ) :
        """'flatte_pdf' : 'original' Flatte pdf (no additional phase space  factors)"""
        return self.__flatte_pdf
    
    @property
    def bwps ( self ) :
        """The Breit-Wigner function (BWPS) itself"""
        return self.__bwps
    
    @property
    def g_list ( self ) :
        """'g_list' list of gammas for Breit-Wigner"""
        return self.__g_list

models.append ( FlattePS_pdf )

# =============================================================================
## @class FlatteBugg_pdf
#  Bugg's modification of Flatte channel 
#  @see D.V. Bugg, "Re-analysis of data on a(0)(1450) and a(0)(980)"
#           Phys.Rev.D 78 (2008) 074023
#  @see https://doi.org/10.1103/PhysRevD.78.074023
#  @see https://arxiv.org/abs/0808.2706
#  Well suitable for \f$f_0(980)\rightarrow \pi^+ \pi^-\f$
#  @see Ostap::Models::FlatteBugg
#  @see Ostap::Math::FlatteBugg
#  @see Ostap::Math::Flatte
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-01-18
class FlatteBugg_pdf(PEAKMEAN) :
    """Bugg's modification of Flatte channel 
    - see D.V. Bugg, "Re-analysis of data on a(0)(1450) and a(0)(980)"
    Phys.Rev.D 78 (2008) 074023
    - see https://doi.org/10.1103/PhysRevD.78.074023
    - see https://arxiv.org/abs/0808.2706
    
    Typical case:    f0(980) -> pi+ pi- & K+ K- shapes 
    """            
    def __init__ ( self              ,
                   name              ,
                   flatte            ,    ## Ostap::Math::FlatteBugg
                   xvar              ,
                   m0       = None   ,    ## the pole 
                   g1       = None   ,    ## g1 
                   g2og1    = None   ,    ## g2/g1 
                   gamma0   = None   ) :  ## gamma0 

        assert isinstance ( flatte , Ostap.Math.FlatteBugg ), \
               'Invalid type for flatte %s' %  type ( flatte )

        ## initialize the base
        with CheckMean ( False ) :
            # for Flatte-function m0 can be outside the interesting interval 
            PEAKMEAN.__init__  ( self , name , xvar ,
                                 mean       = m0  ,
                                 mean_name  = 'm0_%s'      % name ,
                                 mean_title = '#m_{0}(%s)' % name )
            
        self.__flatte = flatte
            
        self.__gamma0 = self.make_var  ( gamma0                  ,
                                         'gamma0_%s'      % name ,
                                         '#Gamma_{0}(%s)' % name ,
                                         True  )
        
        self.__g1     = self.make_var  ( g1                     ,
                                         'g1_%s'          % name ,
                                         '#g_{1}(%s)'     % name ,
                                         None  )
        self.__g2og1  = self.make_var  ( g2og1 ,
                                       'g2og1_%s'          % name ,
                                         '#g_{2}/#g_{1}(%s)' % name ,
                                         None , 0.001 , 200  ) 
        self.__g2     = self.vars_multiply ( self.g2og1 , self.g1 , name = 'g2_%s' % name , title = "g_2(%s)" % name )
        

        ## create PDF 
        self.pdf = Ostap.Models.FlatteBugg (
            self.roo_name ( 'fb_' ) ,
            "FlatteBugg %s" % self.name  ,
            self.xvar    ,
            self.m0      ,
            self.g1      ,
            self.g2      ,
            self.gamma0  ,
            self.flatte  )

        ## save the configuration
        self.config = {
            'name'        : self.name    ,
            'flatte'      : self.flatte  ,
            'xvar'        : self.xvar    ,
            'm0'          : self.m0      ,
            'g1'          : self.g1      ,
            'g2og1'       : self.g2og1   , 
            'gamma0'      : self.gamma0  ,
            }
        
    @property
    def case  (self ) :
        """'case' : How the input argument are  specified: 1: (m0g1,g2og1) vs 2: (g1,g2) """
        return self.__my_case
    
    @property
    def m0 ( self ) :
        """'m0'-parameter for Flatte-function (same as 'mean')"""
        return self.mean 
    @m0.setter
    def m0  ( self, value ) :
        self.set_value ( self.__mean , value ) 

    @property
    def g1 ( self ) :
        """'g1'-parameter for Flatte-function"""
        return self.__g1
    @g1.setter
    def g1 ( self, value ) :
        self.set_value ( self.__g1, value ) 

    @property
    def g2 ( self ) :
        """'g2'-parameter for Flatte-function"""
        return self.__g2

    @property
    def g2og1 ( self ) :
        """'g2/g1'-parameter for Flatte-function"""
        return self.__g2og1
    @g2og1.setter
    def g2og1 ( self, value ) :
        self.set_value ( self.__g2og1, value ) 

    @property
    def gamma0 ( self ) :
        "'gamma0'-parameter for Flatte-function"
        return self.__gamma0
    @gamma0.setter
    def gamma0 ( self , value ) :
        self.set_value ( self.__gamma0, value ) 

    @property
    def flatte ( self ) :
        """The FlatteBugg function itself"""
        return self.__flatte

    # =========================================================================
    ## prepare Argand plot as <code>TGraph</code>
    #  @code
    #  bw = ...
    #  argand = bw.argand ( npx = 1000 )
    #  argand.draw ( 'al')  
    #  @endcode
    #  @see  TGraph 
    def argand ( self , x_min =  None , x_max = None , npx = 1000 ) :
        """ prepare Argand plot as `TGraph`
        >>> bw = ...
        >>> argand = bw.argand ( npx = 1000 )
        >>> argand.draw ( 'al')  
        """
        bw    = self.pdf.function()
        xmnmx = self.xminmax() 
        if x_min is None and xmnmx : x_min = xmnmx [ 0 ] 
        if x_max is None and xmnmx : x_max = xmnmx [ 1 ] 
        ## make Argand plot 
        return bw.argand ( xmin = x_min , xmax = x_max , npx = npx ) 


models.append ( FlatteBugg_pdf )                          

# =============================================================================
## @class LASS_pdf
#  The LASS parameterization (Nucl. Phys. B296, 493 (1988))
#  describes the 0+ component of the Kpi spectrum.
#  It consists of the K*(1430) resonance together with an
#  effective range non-resonant component
#  @see Ostap::Models::LASS
#  @see Ostap::Math::LASS
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class LASS_pdf(PEAK) :
    """Kappa pole:
    The LASS parameterization (Nucl. Phys. B296, 493 (1988))
    describes the 0+ component of the Kpi spectrum.
    It consists of the K*(1430) resonance together with an
    effective range non-resonant component
    """
    def __init__ ( self               ,
                   name               ,
                   xvar               ,
                   m0       = None    ,    ## mass  of K*(1430)
                   g0       = None    ,    ## width of K*(1430)
                   a        = 1.94e-3 , 
                   r        = 1.76e-3 ,
                   e        = 1.0     ,    ## elasticity                    
                   mKaon    = 493.7   ,    ## kaon mass 
                   mPion    = 139.6   ) :  ## pion mass 

        #
        ## initialize the base
        # 
        PEAK.__init__  ( self , name , xvar , 
                         mean        = m0 ,
                         sigma       = g0 ,
                         mean_name   = 'm0_%s'      % name ,
                         mean_title  = '#m_{0}(%s)' % name ,
                         sigma_name  = 'g_%s'       % name ,
                         sigma_title = 'g_0(%s)'    % name )
        
        
        self.__g0 = self.sigma
        self.__m0 = self.mean
            
        self.__a = self.make_var ( a                  ,
                                   'aLASS_%s'  % name ,
                                   "aLASS(%s)" % name ,
                                   None               , 
                                   1.94e-3            ,
                                   1.94e-3            ,
                                   1.94e-3            ) 
        self.__r = self.make_var ( r             ,
                                   'rLASS_%s'  % name ,
                                   "rLASS(%s)" % name ,
                                   None               ,
                                   1.76e-3            ,
                                   1.76e-3            ,
                                   1.76e-3            ) 
        self.__e = self.make_var ( e            ,
                                   'eLASS_%s'  % name ,
                                   "eLASS(%s)" % name ,
                                   None               , 
                                   1.0                , 
                                   0.5                ,
                                   2.0                )
        
        self.__mKaon = mKaon
        self.__mPion = mPion
        
        ## create PDF 
        self.pdf = Ostap.Models.LASS (
            self.roo_name ( 'lass_' ) ,
            "LASS/kappa %s" % self.name  ,
            self.xvar    ,
            self.m0      ,
            self.g0      ,
            self.a       ,
            self.r       ,
            self.e       ,
            self.__mKaon ,
            self.__mPion )

        ## save the configuration
        self.config = {
            'name'        : self.name    ,
            'xvar'        : self.xvar    ,
            'm0'          : self.m0      ,
            'g0'          : self.g0      ,
            'a'           : self.a       ,
            'r'           : self.r       ,
            'e'           : self.e       ,
            'mKaon'       : self.__mKaon ,
            'mPion'       : self.__mPion ,
            }
    
    @property
    def g0 ( self ) :
        """'g0'-parameter for LASS-function (same as 'sigma')"""
        return self.sigma
    @g0.setter
    def g0 ( self, value ) :
        self.sigma = value 

    @property
    def m0 ( self ) :
        """'m0'-parameter for LASS-function (same as 'mean')"""
        return self.__gamma
    @m0.setter
    def m0 ( self, value ) :
        self.mean = value 

    @property
    def a ( self ) :
        """'a'-parameter for LASS-function"""
        return self.__a_
    @a.setter
    def a ( self, value ) :
        value = float ( value )
        self.__a.setVal ( value ) 

    @property
    def r ( self ) :
        """'r'-parameter for LASS-function"""
        return self.__r
    @r.setter
    def r ( self, value ) :
        value = float ( value )
        self.__r.setVal ( value ) 

    @property
    def e ( self ) :
        """'e'-parameter for LASS-function"""
        return self.__e
    @e.setter
    def e ( self, value ) :
        value = float ( value )
        self.__e.setVal ( value ) 

models.append ( LASS_pdf )                          


# =============================================================================
## @class Bugg_pdf
#  The parameterization of sigma pole by B.S.Zou and D.V.Bugg, Phys.Rev. D48 (1993) R3948.
#  @see https://doi.org/10.1103/PhysRevD.48.R3948
#  @see Ostap::Models::Bugg
#  @see Ostap::Math::Bugg
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class Bugg_pdf(PEAK) :
    """ The parameterization of sigma pole by
    B.S.Zou and D.V.Bugg, Phys.Rev. D48 (1993) R3948.
    https://doi.org/10.1103/PhysRevD.48.R3948
    """
    def __init__ ( self           ,
                   name           ,
                   xvar           ,
                   m     = 0.9264 , 
                   g2    = 0.0024 , ## g2-parameter
                   b1    = 0.5848 , ## b1-parameter [GeV]
                   b2    = 1.6663 , ## b2-parameter [GeV^-1]
                   a     = 1.082  , ##  a-parameter [GeV^2]
                   s1    = 2.8    , ## s1-parameter [GeV^2]
                   s2    = 3.5    , ## s2-parameter
                   mPion = 0.1396 ) :  ## pion mass 
        #
        ## initialize the base
        # 
        PEAK.__init__  ( self , name , xvar ,
                         mean        = m    ,
                         sigma       = g2   ,
                         mean_name   = 'mBugg_%s'        % name ,
                         mean_title  = '#m_{Bugg}_2(%s)' % name ,
                         sigma_name  = 'gamma2_%s'       % name ,
                         sigma_title = '#gamma_2(%s)'    % name )

        
        self.__bugg_g2 = self.sigma
        self.__gamma   = self.sigma
        ##
        self.__bugg_m = self.mean

        self.__bugg_b1 = self.make_var ( b1                  ,
                                         'b1Bugg_%s'  % name ,
                                         "b1Bugg(%s)" % name ,
                                         None , 0.5848 , 0 , 2  )
        
        self.__bugg_b2 = self.make_var ( b2             ,
                                         'b2Bugg_%s'  % name ,
                                         "b2Bugg(%s)" % name ,
                                         None ,  1.6663 , 1 , 2  ) 
        
        self.__bugg_a  = self.make_var ( a             ,
                                         'aBugg_%s'  % name ,
                                         "aBugg(%s)" % name ,
                                         None , 1.082 , 0.5 , 5  ) 
        
        self.__bugg_s1  = self.make_var ( s1           ,
                                          's1Bugg_%s'  % name ,
                                          "s1Bugg(%s)" % name ,
                                          None , 2.8 , 1 , 5  ) 
        
        self.__bugg_s2  = self.make_var ( s2           ,
                                          's2Bugg_%s'  % name ,
                                          "s2Bugg(%s)" % name ,
                                          None , 3.5 , 1 , 5  ) 
        
        self.__mPion = mPion 
        ## create PDF 
        self.pdf = Ostap.Models.Bugg ( 
            self.roo_name ( 'bugg_' ) ,
            "Bugg/sigma %s" % self.name  ,
            self.xvar      ,
            self.__bugg_m  ,
            self.__bugg_g2 ,
            self.__bugg_b1 ,
            self.__bugg_b2 ,
            self.__bugg_a  ,
            self.__bugg_s1 ,
            self.__bugg_s2 ,
            self.__mPion )
        
        ## save the configuration
        self.config = {
            'name'        : self.name      ,
            'xvar'        : self.xvar      ,
            'm'           : self.__bugg_m  ,
            'g2'          : self.__bugg_g2 ,
            'b1'          : self.__bugg_b1 ,
            'b2'          : self.__bugg_b2 ,
            'a'           : self.__bugg_a  ,
            's1'          : self.__bugg_s1 ,
            's2'          : self.__bugg_s2 ,
            'mPion'       : self.__mPion   ,
            }

    @property
    def gamma ( self ) :
        """'gamma'-parameter ('g2') for Bugg function"""
        return self.__gamma
    @gamma.setter
    def gamma ( self, value ) :
        value = float ( value )
        self.__gamma.setVal ( value ) 

    @property
    def g2 ( self ) :
        """'g2'-parameter for Bugg function"""
        return self.__bugg_g2
    @g2.setter
    def g2 ( self, value ) :
        value = float ( value )
        self.__bugg_g2.setVal ( value ) 

    @property
    def b1 ( self ) :
        """'b1'-parameter for Bugg function"""
        return self.__bugg_b1
    @b1.setter
    def b1 ( self, value ) :
        value = float ( value )
        self.__bugg_b1.setVal ( value ) 

    @property
    def b2 ( self ) :
        """'b2'-parameter for Bugg function"""
        return self.__bugg_b2
    @b2.setter
    def b1 ( self, value ) :
        value = float ( value )
        self.__bugg_b2.setVal ( value ) 

    @property
    def a ( self ) :
        """'a'-parameter for Bugg function"""
        return self.__bugg_a
    @a.setter
    def a ( self, value ) :
        value = float ( value )
        self.__bugg_a.setVal ( value ) 

    @property
    def s1 ( self ) :
        """'s1'-parameter for Bugg function"""
        return self.__bugg_s1
    @s1.setter
    def s1 ( self, value ) :
        value = float ( value )
        self.__bugg_s1.setVal ( value ) 

    @property
    def s2 ( self ) :
        """'s2'-parameter for Bugg function"""
        return self.__bugg_s2
    @s2.setter
    def s2 ( self, value ) :
        value = float ( value )
        self.__bugg_s2.setVal ( value ) 

        
models.append ( Bugg_pdf )



## # =============================================================================
## ## @class Swanson_pdf
## #  S-wave cusp
## #  @see LHCb-PAPER-2016-019, Appendix
## #  @see E. S. Swanson, Cusps and exotic charmonia, arXiv:1504.07952
## #  @see http://arxiv.org/abs/1504.07952
## #  @see Ostap::Models::Swanson
## #  @see Ostap::Math::Swanson
## #  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
## #  @date 2016-06-12
## class Swanson_pdf(PDF) :
##     """ S-wave cusp
##     - LHCb-PAPER-2016-019, Appendix
##     - E. S. Swanson, Cusps and exotic charmonia, arXiv:1504.07952
##     - http://arxiv.org/abs/1504.07952
##     """
##     def __init__ ( self              ,
##                    name              ,
##                    swanson           , ## Ostap::Math::Swanson objects 
##                    xvar              ,
##                    beta0    = None   ) : 
##         #
##         ## initialize the base
##         #
##         if   isinstance ( xvar , ROOT.TH1   ) :
##             m_title = xvar.GetTitle ()            
##             xvar    = xvar.xminmax  ()
##         elif isinstance ( xvar , ROOT.TAxis ) :
##             xvar    = xvar.GetXmin() , mass.GetXmax()
            
##         ## create the variable 
##         if isinstance ( xvar , tuple ) and 2 == len(mass) :  
##             xvar = self.make_var ( xvar       , ## var 
##                              m_name     , ## name 
##                              m_title    , ## title/comment
##                              *mass      , ## min/max 
##                              fix = None ) ## fix ? 
##         elif isinstance ( xvar , ROOT.RooAbsReal ) :
##             xvar = self.make_var ( xvar       , ## var 
##                              m_name     , ## name 
##                              m_title    , ## title/comment
##                              fix = None ) ## fix ? 
##         else :
##             raise AttributeError("Swanson: Unknown type of 'xvar' parameter %s/%s" % ( type ( xvar ) , xvar ) )

            
##         PDF.__init__  ( self , name , xvar , None , None    ) 
        
##         self.__swanson = swanson 
##         beta_max       = max ( swanson.mmin() , swanson.cusp() )
##         self.__beta0   = self.make_var ( beta0 , 
##                                    'b0_swanson_%s'   % name ,
##                                    'b0_swanson(%s)'  % name ,
##                                    beta0 , 
##                                    0 , beta_max )
##         ## create PDF 
##         self.pdf = Ostap.Models.Swanson ( 
##             "Swanson_"    + name ,
##             "Swanson(%s)" % name ,
##             self.xvar    ,
##             self.beta0   ,
##             self.swanson )
        
##         ## save the configuration
##         self.config = {
##             'name'        : self.name      ,
##             'swanson'     : self.swanson   ,
##             'xvar'        : self.xvar      ,
##             'beta0'       : self.beta0     ,
##             }
    
##     @property
##     def beta0 ( self ) :
##         """'beta0'-parameter for Swanson function"""
##         return self.__beta0
##     @beta0.setter
##     def beta0 ( self, value ) :
##         value = float ( value )
##         self.__beta0.setVal ( value ) 

##     @property
##     def swanson ( self ) :
##         """'swanson'-function itself for Swanson PDF"""
##         return self.__swanson

## models.append ( Swanson_pdf )

# =============================================================================
if '__main__' == __name__ : 
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
