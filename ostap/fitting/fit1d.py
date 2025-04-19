#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostap/fitting/fit1d.py
#  Set of useful basic utilities to build various 1D-fit models 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
# =============================================================================
"""Set of useful basic utilities to build various 1D-fit models"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'PEAKMEAN'      , ## useful base class to create "signal" PDFs for peak-like fits
    'PEAK'          , ## useful base class to create "signal" PDFs for peak-like fits
    'RESOLUTION'    , ## useful base class to create "resolution" PDFs
    ##
    'Fit1D'         , ## the basic compound 1D-fit model 
    ##
    )
# =============================================================================
from   ostap.core.core          import Ostap , VE , valid_pointer, roo_silent 
from   ostap.core.ostap_types   import ( is_integer     , string_types   , 
                                         integer_types  , num_types      ,
                                         list_types     , all_numerics   ) 
from   ostap.fitting.funbasic   import FUN1
from   ostap.fitting.pdfbasic   import PDF1, APDF1, Sum1D
from   ostap.fitting.utils      import make_name
import ROOT, math,  random
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.basic' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
##  helper utilities to imlement resolution models.
# =============================================================================
class _CHECKMEAN(object) : check = True
def checkMean() :
    return True if  _CHECKMEAN.check else False
# =============================================================================
## @class CheckMean 
#  Helper contex manager to enable/disable check for the mean/location-values
class CheckMean(object) :
    """ Helper contex manager to enable/disable check for the mean/location-values
    """
    def __init__  ( self , check ) :
        self.__check = True if check else False 
    def __enter__ ( self ) :
        self.__old       = _CHECKMEAN.check 
        _CHECKMEAN.check =  self.__check
    def __exit__  ( self , *_ ) :
        _CHECKMEAN.check =  self.__old
    @property
    def check ( self ) :
        """'check'  : check the mean/location?"""
        return self.__check
    
# =============================================================================
## helper base class for implementation  of various helper pdfs
#  - it defines alias <code>mass</code> for <code>xvar</code>
#  - it defiens a variable <code>mean</code> alias <code>location</code>
#  - optionally it checks that this variable is withing the specified range  
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class PEAKMEAN(PDF1) :
    """ Helper base class for implementation of various pdfs
    It is useful for 'peak-like' distributions, where one can talk about
    - 'mean/location'
    - it defines alias `mass` for `xvar`
    - it defiens a variable `mean` (alias `location`)
    - optionally it checks that this variable is withing the specified range  
    """
    def __init__ ( self              ,
                   name              ,
                   xvar              ,
                   mean      = None  ,
                   mean_name  = ''   , 
                   mean_title = ''   ) : 

        
        if   isinstance ( xvar , ROOT.TH1   ) : xvar = xvar.xminmax ()
        elif isinstance ( xvar , ROOT.TAxis ) : xvar = xvar.GetXmin () , xvar.GetXmax ()

        ## intialize the base 
        PDF1.__init__ ( self , name , xvar = xvar )

        ## check mean/location values ? 
        self.__check_mean = self.xminmax () and checkMean () 
        
        self.__limits_mean  = ()
        if  self.check_mean and self.xminmax () and not isinstance ( mean , ROOT.RooAbsReal ) :      
            mn , mx = self.xminmax()
            dm      =  mx - mn
            self.__limits_mean  = mn - 0.35 * dm , mx + 0.35 * dm

        ## mean-value
        m_name  = mean_name  if mean_name  else "mean_%s"  % name
        m_title = mean_title if mean_title else "mean(%s)" % name
        self.__mean = self.make_var ( mean , m_name , m_title , False , *self.limits_mean )
        
        ##
        if self.limits_mean :  
            mn , mx = self.limits_mean  
            dm      =  mx - mn
            if   self.mean.isConstant() :
                if not mn <= self.mean.getVal() <= mx : 
                    self.error ( 'PEAKMEAN(%s): Fixed mass %s is not in mass-range (%s,%s)' % ( name , self.mean.getVal() , mn , mx  ) )
            elif self.mean.minmax() :
                mmn , mmx = self.mean.minmax()
                self.mean.setMin ( max ( mmn , mn ) )
                self.mean.setMax ( min ( mmx , mx ) )
                self.debug ( 'mean range is adjusted  to be %s' % list ( self.mean.minmax() ) )

        ## save the configuration
        self.config = {
            'name'        : self.name  ,
            'xvar'        : self.xvar  ,
            'mean'        : self.mean  ,
            'mean_name'   : mean_name  ,
            'mean_title'  : mean_title ,
            }

    @property 
    def mass ( self ) :
        """'mass'-variable (the same as 'x' and 'xvar')"""
        return self.xvar
    
    @property
    def mean ( self ):
        """'mean/location''-variable (the same as 'location')"""
        return self.__mean
    @mean.setter
    def mean ( self , value ) :
        value =  float ( value )
        mn , mx = self.mean.minmax()
        if not mn <= value <= mx :
            self.warning( "'%s'': %s is outside the interval (%s,%s)/1" % ( self.mean.name , value , mn , mx ) )
        if self.check_mean and self.limits_mean  :  
            mn , mx = self.limits_mean 
            if not mn <= value <= mx :
                self.error ("'%s'': %s is outside the interval (%s,%s)/2"  % ( self.mean.name , value , mn , mx ) )                
        self.mean.setVal ( value )
        
    @property
    def location ( self ):
        """'location/mean'-variable (the same as 'mean')"""
        return self.mean
    @location.setter
    def location ( self , value ) :
        self.mean =  value

    @property
    def check_mean ( self ) :
        """'check_mean' : Is mean-value to be checked?"""
        return self.__check_mean
    @check_mean.setter
    def check_mean ( self, value ) :
        self.__check_mean = True if  value else False
        
    @property
    def limits_mean ( self ) :
        """'limits_mean' : reasonable limits for mean/location"""
        return self.__limits_mean
    
# =============================================================================
## @class PEAK
#  helper base class for implementation  of various helper peak-like pdfs 
#  - mean/location
#  - sigma/width/scale
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class PEAK(PEAKMEAN) :
    """ Helper base class for implementation of various signal-like pdfs
    It is useful for 'peak-like' distributions, where one can talk about
    - 'mean/location'
    - 'sigma/width/scale' 
    """
    def __init__ ( self               ,
                   name               ,
                   xvar               ,
                   mean        = None ,
                   sigma       = None , 
                   mean_name   = ''   , 
                   mean_title  = ''   ,
                   sigma_name  = ''   , 
                   sigma_title = ''   ) : 

        ## base class 
        PEAKMEAN.__init__ ( self                    ,
                            name       = name       ,
                            xvar       = xvar       , 
                            mean       = mean       ,
                            mean_name  = mean_name  ,
                            mean_title = mean_title )
        
        self.__limits_sigma = ()        
        if  self.xminmax() and not isinstance ( sigma , ROOT.RooAbsReal ) :            
            mn , mx   = self.xminmax()
            dm        = mx - mn
            sigma_max = 3 * dm / math.sqrt(12)  
            self.__limits_sigma = 1.e-4 * sigma_max , sigma_max 

        ## sigma
        s_name  = sigma_name  if sigma_name  else "sigma_%s"   % name
        s_title = sigma_title if sigma_title else "#sigma(%s)" % name
        #
        self.__check_sigma = True 
        self.__sigma = self.make_var ( sigma  , s_name , s_title , False , *self.limits_sigma )
        
        ## save the configuration
        self.config = {
            'name'        : self.name   ,
            'xvar'        : self.xvar   ,
            'mean'        : self.mean   ,
            'sigma'       : self.sigma  ,
            'mean_name'   : mean_name   ,
            'mean_title'  : mean_title  ,
            'sigma_name'  : sigma_name  ,
            'sigma_title' : sigma_title ,
            }
            
    @property
    def sigma ( self ):
        """'sigma/width/scale/spread'-variable"""
        return self.__sigma
    @sigma.setter
    def sigma ( self , value ) :
        value =   float ( value )
        mn , mx = self.sigma.minmax()
        if not mn <= value <= mx :
            self.warning ("'%s': %s is outside the interval (%s,%s)/1" % ( self.sigma.name , value , mn , mx ) )
        if self.limits_sigma and self.check_sigma  : 
            mn , mx = self.limits_sigma 
            if not mn <= value <= mx :
                self.error ("'%s': %s is outside the interval (%s,%s)/2" % ( self.sigma.name , value , mn , mx ) )
        self.sigma.setVal ( value )

    @property
    def check_sigma ( self ) :
        """'check_mean' : Is mean-value to be checked?"""
        return self.__check_sigma 
    @check_sigma.setter
    def check_sigma ( self, value ) :
        self.__check_sigma = True if  value else False 
    
    @property
    def limits_sigma ( self ) :
        """'limits_sigma' : reasonable limits for sigma/width"""
        return self.__limits_sigma

# =============================================================================
## @class RESOLUTION
#  helper base class  to parameterize the resolution
#  - It allows setting of the <code>mean</code> to zero,
#  - It containg "fudge-factor" for the resolution parameter <code>sigma</code>
#  - It simplify creation of the soft/gaussian constraint for the "fudge-factor"
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2017-07-13
class RESOLUTION(PEAK) :
    """ Helper base class  to parameterize the resolution
    - It allows setting of the 'mean' to zero,
    - It contains 'fudge-factor' for the resolution parameter 'sigma'
    - It simplify creation of the soft/gaussian constraint for the 'fudge-factor'    
    """
    ## constructor
    #  @param name   the name of PDF
    #  @param xvar   the variable/observable
    #  @param sigma  sigma/resoltuion parameter 
    #  @param mean   "mean"-variable
    #  @param fudge  "fudge-factor" to be aplied to sigma
    def __init__ ( self               ,
                   name               ,
                   xvar        = None ,
                   sigma       = None , 
                   mean        = None ,
                   fudge       = 1.0  ,
                   mean_name   = ''   ,
                   mean_title  = ''   ,
                   sigma_name  = ''   ,
                   sigma_title = ''   ) :
        
        ## mean-value
        if mean is None : mean = ROOT.RooFit.RooConst ( 0 ) 
            
        with CheckMean ( False ) :
            super(RESOLUTION,self).__init__ ( name        = name        ,
                                              xvar        = xvar        ,
                                              sigma       = sigma       ,
                                              mean        = mean        ,
                                              mean_name   = mean_name   ,
                                              mean_title  = mean_title  ,
                                              sigma_name  = sigma_name  ,
                                              sigma_title = sigma_title )
            
        self.__fudge            = fudge        
        self.__fudge_constraint = None
        self.__sigma_corr       = None
        
        if isinstance ( fudge , VE ) :
            
            assert 0 < fudge.value() and 0 < fudge.cov2(),\
                   "Invalid value for 'fudge-factor': %s" % s 
            
            value  = fudge.value()
            error  = fudge.error()
            vmin   = max ( 1.e-3 , value - 10 * error )
            vmax   =               value + 10 * error
            
            ## make fudge-factor to be a variable 
            self.__fudge = self.make_var ( value ,
                                           'fudge_factor_%s'  % self.name ,
                                           'fudge_factor(%s)' % self.name ,
                                           True , vmin , vmax            )
            
            ## create soft/gaussian constraint for fudge-factor
            self.__fudge_constraint = self.soft_constraint (
                self.fudge ,
                fudge      ,
                name  = 'Fudge_constraint_%s'  % self.name ,
                title = 'Fudge_constraint(%s)' % self.name )
            
        elif isinstance ( fudge , ROOT.RooAbsReal ) :
            
            ## make fudge-factor to be a variable 
            self.__fudge = self.make_var ( fudge  ,
                                           'fudge_factor_%s'  % self.name ,
                                           'fudge_factor(%s)' % self.name )
            
        elif isinstance ( fudge , num_types ) and 1 == fudge :

            ## fudge is trivial 
            self.__fudge = ROOT.RooFit.RooConst ( fudge )

            ## corrected sigma is trivial 
            self.__sigma_corr = self.sigma
            
        elif isinstance ( fudge , num_types ) :
            
            ## fudge is trivial 
            self.__fudge = ROOT.RooFit.RooConst ( fudge )

        else :
                
            ## make fudge-factor to be a variable 
            self.__fudge = self.make_var ( fudge  ,
                                           'fudge_factor_%s'  % self.name ,
                                           'fudge_factor(%s)' % self.name ,
                                           False  , 0.2 , 3.0 ) 
            
        ## create corrected sigma 
        if self.__sigma_corr is None :            
            ## corrected sigma 
            self.__sigma_corr = self.vars_multiply ( self.sigma ,
                                                     self.fudge ,
                                                     'Corrected_%s' % self.sigma.name )
            
        ## save the configuration
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'mean'  : self.mean  ,
            'sigma' : self.sigma ,
            'fudge' : self.fudge ,
            }

    @property 
    def fudge ( self ) :
        """'fudge' : fudge factor for resolution"""
        return self.__fudge
    @property 
    def fudge_constraint ( self ) :
        """'fudge_constraint' : constraint for fudge factor for the resolution"""
        return self.__fudge_constraint
    @property 
    def sigma_corr ( self ) :
        """'sigma_corr' : the corrected sigma parameter: sigma*fudge """
        return self.__sigma_corr

    # =========================================================================
    ## right modulus: convoluttion
    #  @code
    #  pdf        = ...
    #  resolution = ...
    #  result  = pdf % resolution ## python 2&3
    #  result  = pdf @ resolution ## python 3 only
    #  @endcode 
    #  The configuration can be specified via `ConvolutionConfig`
    #    context manager:
    #  @code
    #  pdf        = ...
    #  resolution = ... 
    #  with ConvolutionConfig ( buffer = 0.25 , nbins = 1000 ):  
    #      result = pdf % other ## python 2 and 3 
    #      result = pdf @ other ## python 3 only
    #  @endcode 
    def __rmod__ ( self , other ) :
        """ Right modulus: convolution 
        >>> pdf = ...
        >>> resolution = ...
        >>> result  = pdf % resolution ## python 2&3
        >>> result  = pdf @ resolution ## python 3 only        
        The configuration can be specified via `ConvolutionConfig` context manager:
        >>> pdf        = ...
        >>> resolution = ... 
        >>> with ConvolutionConfig ( buffer = 0.25 , nbins = 1000 ):  
        ...     result = pdf % other ## python 2 an d3 
        ...     result = pdf @ other ## python 3 only
        """        
        from ostap.fitting.pdf_ops import pdf_convolution
        return pdf_convolution ( other , self )
    
    __rmatmult__ = __rmod__     
    
    

# =============================================================================
## @class Fit1D
#  The actual model for 1D-mass fits
#  @param signal                PDF for 'signal'     component                 (Ostap/PDF or RooAbsPdf)
#  @param background            PDF for 'background' component                 (Ostap/PDF or RooAbsPdf)
#  @param othersignals          list of PDFs for other 'signal' components     (Ostap/PDF or RooAbsPdf)
#  @param otherbackgrouds       list of PDFs for other 'background' components (Ostap/PDF or RooAbsPdf)
#  @param others                list of 'other' components                     (Ostap/PDF or RooAbsPdf)
#  @param signals               list 'signal'     components
#  @param backgrounds           list of 'background' component               
#  @param suffix                ... add this  suffix for the PDF name
#  @param name                  the name of compound PDF 
#  @param xvar                  the fitting variable, must be specified if components are given as RooAbsPdf
#  @param extended              build 'extended' PDF
#  @param combine_signals       combine all signal components into single SIGNAL?
#  @param combine_backgrounds   combine all background components into single BACKGROUND?
#  @param combine_others        combine all other components into single COMPONENT?
#  @param recursive             use recursive fractions for compound PDF
#  @param recursive_signals     use recursive fractions for compound signal
#  @param recursive_backgriunds use recursive fractions for compound background
#  @param recursive_others      use recursive fractions for compound others
#  @param S                     yields of signal components 
#  @param B                     yields of background components 
#  @param C                     yields of others components
#  @param F                     component fractions for non-extended fit
#  @param fS                    fractions for compound signal
#  @param fB                    fractions for compound background
#  @param fC                    fractions for compound others 
#  @code 
#  gauss = Gauss_pdf( ... ) 
#  pdf   = Fit1D ( signal = gauss , background = 0 ) ## Gauss as signal and exponent as background
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-08-02
class Fit1D (PDF1) :
    """The actual fit-model for generic 1D-fits
    Parameters
    - signal              : PDF for 'signal'     component                 (Ostap/PDF or RooAbsPdf)
    - background          : PDF for 'background' component                 (Ostap/PDF or RooAbsPdf)
    - othersignals        : list of PDFs for other 'signal' components     (Ostap/PDF or RooAbsPdf)
    - otherbackgrouds     : list of PDFs for other 'background' components (Ostap/PDF or RooAbsPdf)
    - others              : list of 'other' components                     (Ostap/PDF or RooAbsPdf)
    - name                : The name of compound PDF 
    - suffix              : ... add this  suffix for the PDF name
    - extended            : build 'extended' PDF
    - combine_signals     : combine all signal components into single SIGNAL?
    - combine_backgrounds : combine all background components into single BACKGROUND?
    - combine_others      : combine all other components into single COMPONENT?
    - recursive           : use recursive fractions for compound PDF
    - xvar                : the fitting variable, must be specified if components are given as RooAbsPdf

    >>> gauss = Gauss_pdf( ... ) 
    >>> pdf   = Fit1D ( signal = gauss , background = 0 ) ## Gauss as signal and exponent as background 
    """
    def __init__ ( self                            , 
                   signal                = None    ,    ## the main signal 
                   background            = None    ,    ## the main background 
                   othersignals          = ()      ,    ## additional signal         components
                   otherbackgrounds      = ()      ,    ## additional background     components
                   others                = ()      ,    ## additional non-classified components
                   signals               = ()      ,    ## alternative : all signals 
                   backgrounds           = ()      ,    ## alternative : all backgrounds  
                   suffix                = ''      ,    ## the suffix 
                   name                  = ''      ,    ## the name
                   xvar                  = None    ,    ## x-variable 
                   extended              = True    ,    ## extended fits ?
                   combine_signals       = False   ,    ## combine signal     components into single "SIGNAL"     ? 
                   combine_backgrounds   = False   ,    ## combine background components into single "BACKGROUND" ?            
                   combine_others        = False   ,    ## combine other      components into single "COMPONENT"  ?             
                   recursive             = True    ,    ## recursive fractions for NON-extended models?
                   recursive_signals     = True    ,    ## recursive fractions for combined signal?
                   recursive_backgrounds = True    ,    ## recursive fractions for combined background?
                   recursive_others      = True    ,    ## recursive fractions for combined other components?
                   S                     = ()      ,    ## yields for 'signals'
                   B                     = ()      ,    ## yields for 'background'
                   C                     = ()      ,    ## yields for 'components'
                   F                     = ()      ,    ## fractions for noin-extended fit 
                   fS                    = ()      ,    ## fraction for combined signal
                   fB                    = ()      ,    ## fraction for combined background
                   fC                    = ()      ,    ## fraction for combined components
                   fix_norm              = False   ,    ## SetCoefNoralization.getCoefNormalization
                   **kwargs                        ) :  ## other arguments (e.g. drawing)

        
        self.__suffix                = suffix
        self.__extended              = True if extended              else False
        
        self.__combine_signals       = True if combine_signals       else False
        self.__combine_backgrounds   = True if combine_backgrounds   else False
        self.__combine_others        = True if combine_others        else False

        self.__recursive             = True if recursive             else False
        self.__recursive_signals     = True if recursive_signals     else False
        self.__recursive_backgrounds = True if recursive_backgrounds else False
        self.__recursive_others      = True if recursive_others      else False

        self.__signal_components     = ()
        self.__background_components = () 
        self.__other_components      = ()
        
        # =====================================================================
        ## Signals
        # =====================================================================
        
        assert       signals   or       signal         , "Fit1D:'signals' or 'signal' must be specified!"        
        assert ( not signals ) or ( not signal       ) , "Fit1D:'signals' and 'signal' are mutually exclusive!"        
        assert ( not signals ) or ( not othersignals ) , "Fit1D:'signals' and 'othersignals' are mutually exclusive!"  
        
        sig_lst = list ( othersignals ) + list ( signals ) 
        if signal : sig_lst.insert ( 0 , signal )

        pdfs = [] 
        for i , signal in enumerate ( sig_lst ) :
            s , xvar = self.make_PDF1 ( signal , xvar , prefix = "Sig%d_" % i , suffix = self.suffix )
            pdfs.append ( s )


        ## sinal components 
        self.__signal_components  = tuple ( pdfs )
        
        # =====================================================================
        ## initialize the base class
        # =====================================================================
        name = name if name else self.generate_name ( prefix = 'Fit%s' % self.signal_components[0].name , suffix = self.suffix ) 
        PDF1.__init__ ( self , name , xvar , **kwargs ) 
                            
        # =====================================================================
        ## Backgrounds
        # =====================================================================

        assert ( not backgrounds ) or ( not otherbackgrounds ) ,\
               "Fit1D:'backgrounds' and 'otherbackgrounds' are mutually exclusive!"        
        assert ( not backgrounds ) or ( not background ) ,\
               "Fit1D:'backgrounds' and 'background' are mutually exclusive!"        
        
        if not backgrounds :
            ## create background
            bkg_name = 'Background_%s' % self.name if not self.suffix else 'Background_%s' % self.suffix 
            background = self.make_bkg ( background , name = bkg_name , xvar = self.xvar ) 
            
        bkg_lst  = list ( otherbackgrounds ) + list ( backgrounds )
        if background : bkg_lst.insert ( 0 , background ) 

        pdfs = [ self.make_PDF1 (
            cmp  ,
            xvar ,
            prefix = "Bkg%d_" % i ,
            suffix = suffix       ) [ 0 ]  for ( i , cmp ) in enumerate ( bkg_lst ) ]

        ## background components 
        self.__background_components = tuple ( pdfs ) 
        
        # =====================================================================
        ## Other fit components
        # =====================================================================
        
        pdfs = [ self.make_PDF1 (
            cmp  ,
            xvar ,
            prefix = "Cmp%d_" % i ,
            suffix = self.suffix ) [ 0 ]  for ( i, cmp ) in enumerate ( others ) ]
                     
        ## all other componnets 
        self.__other_components = tuple ( pdfs ) 
        
        # =====================================================================
        ## Merge them if requested 
        # =====================================================================
        self.__combined_signal     = None
        self.__combined_background = None
        self.__combined_others     = None

        sigs = list ( self.signal_components     )
        bkgs = list ( self.background_components ) 
        cmps = list ( self.other_components      ) 
        
        if combine_signals     and 2 <= len ( self.signal_components     ) :            
            combined = Sum1D ( self.signal_components ,
                               xvar      = self.xvar     , 
                               prefix    = 'fS'          ,
                               suffix    = suffix        ,
                               recursive = recursive     ,                                              
                               fractions = fS            ) ## read  fS from arguments
            
            self.__combined_signal = combined 
            sigs = [ combined ] 
            
        if combine_backgrounds and 2 <= len ( self.background_components ) :                
            combined = Sum1D ( self.background_components ,
                               xvar      = self.xvar      , 
                               prefix    = 'fB'           ,
                               suffix    = suffix         ,
                               recursive = recursive      ,                                             
                               fractions = fB             ) ## read fB from arguments 
            
            self.__combined_background = combined
            bkgs = [ combined ]  
            
        if combine_others      and 2 <= len ( self.other_components           ) :                        
            combined = Sum1D ( self.other_components     ,
                               xvar      = self.xvar     , 
                               prefix    = 'fC'          ,
                               suffix    = suffix        ,
                               recursive = recursive     ,                                              
                               fractions = fC            ) ## read fC from arguments 
            
            self.__combined_others = combined
            cmps = [ combined ] 

        self.__fit_signals     = tuple ( sigs )
        self.__fit_backgrounds = tuple ( bkgs )
        self.__fit_others      = tuple ( cmps )
        
        ## final list of fit components 
        self.__fit_components  = self.fit_signals + self.fit_backgrounds + self.fit_others 
        for p in self.fit_components  : self.alist1.add ( p.pdf )  

        ## Yields/fractions 
        
        self.__S = () 
        self.__B = () 
        self.__C = ()
        self.__F = ()
        
        if self.extended :
            
            ns        = len  ( self.fit_signals ) 
            fname     = make_name ( 'S' , '%d' if 1 != ns else '' , suffix )
            title     = "Yield(s) for 'signal' component(s)/%s" % self.name
            title     = title if 1 != ns else title.replace ( '(s)', '' )

            self.__S  = self.make_yields ( ns , fname , title , yields = S ) ## read S from arguments
            
            nb        = len ( self.fit_backgrounds ) 
            fname     = make_name ( 'B' , '%d' if 1 != nb else '' , suffix )
            title     = "Yield(s) for 'background' component(s)/%s" % self.name
            title     = title if 1 != nb else title.replace ( '(s)', '' )             
            self.__B  = self.make_yields ( nb , fname , title , yields = B ) ## read S from arguments

            nc        = len ( self.fit_others  ) 
            fname     = make_name ( 'C' , '%d' if 1 != nc else '' , suffix )
            title     = "Yield(s) for 'other' component(s)/%s" % self.name
            title     = title if 1 != nc else title.replace ( '(s)', '' )             
            self.__C  = self.make_yields ( nc , fname , title , yields = C ) ## read S from arguments

            for y in self.yields : self.alist2.add ( y )   

            assert len ( self.alist2 ) == len ( self.alist1 ) ,\
                   'Fit1D: inconsistent parameters for ROOT.RooAddPdf' 

        else :
            
            assert 1 < len ( self.fit_components ) ,\
                   'Fot1D: At least two components are required to build proper non-extended PDF!'

            nt       = len ( self.fit_components )
            nf       = nt - 1 
            fname    = make_name ( 'F' , '%d' if 1 != nf else '' , suffix )
            title    = "Fraction(s) for various component(s)/%s" % self.name
            title    = title if 1 != nf else title.replace ( '(s)', '' )                         
            self.__F = self.make_fractions ( nt , fname ,  title , fractions = F ) ## read F from arguments
                        
            for f in self.__F : self.alist2.add ( f )

            assert len ( self.alist2 ) + 1 == len ( self.alist1 ) ,\
                   'Fit1D: inconsistent parameters for ROOT.RooAddPdf' 

        ## now we finally can create PDF    
        pdf_name  = self.new_roo_name ( 'fit1d' , suffix ) 
        pdf_title = "Fit1D %s" % self.name
        pdf_args  = pdf_name , pdf_title , self.alist1 , self.alist2

        if not self.extended :
            pdf_args = pdf_args + ( True if recursive else False , ) ## RECURSIVE ?
            
        self.pdf = ROOT.RooAddPdf     ( *pdf_args )
        
        if fix_norm : self.pdf.fixCoefNormalization ( self.vars ) ## VB: added 10/10/2024 to suppress warnings 

        ## sanity checks

        ## drawing stuff
        if self.combined_background         : self.combined_backgrounds.add ( self.combined_background.pdf ) ## for drawing 
        if self.combined_signal             : self.combined_signals    .add ( self.combined_signal    .pdf ) ## for drawing 
        if self.combined_others             : self.combined_components .add ( self.combined_others    .pdf ) ## for drawing 

        for p in self.signal_components     : self.signals    .add ( p.pdf ) 
        for p in self.background_components : self.backgrounds.add ( p.pdf ) 
        for p in self.other_components      : self.components .add ( p.pdf ) 

        ## save the configuration
        self.config = {
            ## 
            'signals'               : self.signal_components     ,
            'backgrounds'           : self.background_components ,
            'others'                : self.other_components      ,
            ##
            'suffix'                : self.suffix                ,
            'name'                  : self.name                  ,            
            'extended'              : self.extended              ,
            'xvar'                  : self.xvar                  ,
            ## 
            'combine_signals'       : self.combine_signals       ,
            'combine_backgrounds'   : self.combine_backgrounds   ,
            'combine_others'        : self.combine_others        ,
            ## 
            'recursive'             : self.recursive             ,
            'recursive_signals'     : self.recursive_signals     ,
            'recursive_backgrounds' : self.recursive_backgrounds ,
            'recursive_others'      : self.recursive_others      ,
            ##
            'fS'                    : self.fS                    ,
            'fB'                    : self.fB                    ,
            'fC'                    : self.fC                    ,
            ##                      
            'S'                     : self.S                     ,
            'B'                     : self.B                     ,
            'C'                     : self.C                     ,
            'F'                     : self.F                     ,
            ##
            'fix_norm'              : self.fix_norm
            ##
            }

        self.checked_keys.add  ( 'xvar' )

    @property
    def extended ( self ) :
        """'extended' : build extended PDF?"""
        return  self.__extended 
    @property
    def suffix   ( self ) :
        """'suffix'  : append the names  with the specified suffix"""
        return self.__suffix

    @property
    def combine_signals ( self ) :
        """Combine all 'signal'-components into single 'signal' componet?"""
        return self.__combine_signals
    @property
    def combine_backgrounds ( self ) :
        """Combine all 'background'-components into single 'background'  componet?"""
        return self.__combine_backgrounds
    @property
    def combine_others ( self ) :
        """Combine all 'others'-components into single 'other' componet?"""
        return self.__combine_others 
    @property
    def recursive ( self ) :
        """'recursive':  use recursive fit fractions for non-extended fit?"""
        return  self.__recursive
    @property
    def recursive_signals ( self ) :
        """'recursive_signals' :  use recursive fractions for combined signal?"""
        return  self.__recursive_signals
    @property
    def recursive_backgrounds ( self ) :
        """'recursive_backgrounds' :  use recursive fractions for combined background?"""
        return  self.__recursive_backgrounds
    @property
    def recursive_others      ( self ) :
        """'recursive_backgrounds' :  use recursive fractions for combined other components?"""
        return  self.__recursive_others

    @property
    def signal_components ( self )  :
        """'signal_components' : all 'signal'' components"""
        return self.__signal_components         
    @property
    def background_components ( self )  :
        """'background_components'' : all 'background'' components"""
        return self.__background_components 
    @property
    def other_components ( self )  :
        """'other_components'  : all 'other'' components"""
        return self.__other_components 

    @property
    def combined_signal     ( self ) :
        """'combined_signal' : PDF for combined 'signal'' component"""
        return self.__combined_signal
    @property
    def combined_background ( self ) :
        """'combined_background' : PDF for combined 'background'' component"""
        return self.__combined_background
    @property
    def combined_others    ( self ) :
        """'combined_background' : PDF for combined 'others'' component"""
        return self.__combined_others

    @property
    def fS ( self ) :
        """'fS'' : fractions (possible recursive) of components in combined signal"""
        return () if not self.combined_signal else self.combined_signal.F
    @fS.setter
    def fS ( self , value ) :
        assert  self.combined_signal, "'fS'': no combined signal is defined!"
        self.combined_signal.F = value

    @property
    def fB ( self ) :
        """'fB'' : fractions (possible recursive) of components in combined background"""
        return () if not self.combined_background else self.combined_background.F
    @fB.setter
    def fB ( self , value ) :
        assert  self.combined_background, "'fB'': no combined background is defined!"
        self.combined_background.F = value

    @property
    def fC ( self ) :
        """'fC'' : fractions (possible recursive) of components in combined 'others''"""
        return () if not self.combined_others else self.combined_others.F
    @fC.setter
    def fC ( self , value ) :
        assert  self.combined_others, "'fC'': no combined 'others'' is defined!"
        self.combined_others.F = value

    @property
    def fit_components  ( self ) :
        """'fit_components'' : list of fit components"""
        return self.__fit_components
    @property
    def fit_signals      ( self ) :
        """'fit_signals'': list of (the 1st order) 'signal'' components in the model"""
        return self.__fit_signals 
    @property
    def fit_backgrounds  ( self ) :
        """'fit_backgrounds'': list of (the 1st order) 'background'' components in the model"""
        return self.__fit_backgrounds
    @property
    def fit_others       ( self ) :
        """'fit_others'': list of (the 1st order) 'others'' components in the model"""
        return self.__fit_others 

    @property
    def signal ( self ) :
        """'signal'' : get 'signal'' PDF ('combined_signal'' or the 1st from 'signal_components'')"""
        if self.__combined_signal : return self.__combined_signal
        assert self.signal_components, "signal: empty lst of 'signal'' components!"    
        if 1 != len ( self.signal_components ) :
            self.warning ("signal: get the 1st 'signal'' component")
        return self.signal_components[0]
    @property
    def background ( self ) :
        """'background'' : get 'background'' PDF ('combined_background'' or the 1st from 'background_components'')"""
        if self.combined_background : return self.combined_background
        assert self.background_components, "background: empty lst of 'background'' components!"    
        if 1 != len ( self.background_components ) :
            self.warning ("background: get the 1st 'background'' component")
        return self.background_components[0]
        
    @property
    def S ( self ) :
        """Get the  yields of signal component(s) (empty for non-extended fits)
        For single signal component:
        >>> print pdf.S          ## read the single single component 
        >>> pdf.S = 100          ## assign to it
        For multiple signal components:
        >>> print pdf.S[4]       ## read the 4th signal component 
        >>> pdf.S = (1,2,3,4,5,6)## assign to it 
        ... or, alternatively:
        >>> print pdf.S[4]       ## read the 4th signal component 
        >>> pdf.S[4].value = 100 ## assign to it         
        """
        return () if not self.extended else self.component_getter ( self.__S )    
    @S.setter
    def S (  self , value ) :
        assert self.extended, "'S'' cannot be set for non-exteded model!"
        self.component_setter ( self.__S , value )

    @property
    def B ( self ) :
        """Get the  yields of background  component(s) (empty for non-extended fits)
        For single background component:
        >>> print pdf.B          ## read the single background component 
        >>> pdf.B = 100          ## assign to it 
        For multiple background components:
        >>> print pdf.B[4]            ## read the 4th background component 
        >>> pdf.B = ( 1, 2, 3, 4, 5 ) ## assign to it 
        ... or, alternatively:
        >>> print pdf.B[4]       ## read the 4th background component 
        >>> pdf.B[4].value = 100 ## assign to it 
        """
        return () if not self.extended else self.component_getter ( self.__B ) 
    @B.setter
    def B (  self , value ) :
        assert self.extended, "'B'' cannot be set for non-extended model!"
        self.component_setter ( self.__B , value )

    @property
    def C ( self ) :
        """Get the  yields of 'other'' component(s) (empty for non-extended fits)
        For single 'other'' component:
        >>> print pdf.C           ## read the single 'other'' component 
        >>> pdf.C = 100           ## assign to it 
        For multiple 'other'' components:
        >>> print pdf.C[4]            ## read the 4th 'other'' component 
        >>> pdf.C = ( 1, 2, 3, 4, 5 ) ## assign to it 
        ... or, alternatively:
        >>> print pdf.C[4]        ## read the 4th 'other'' component 
        >>> pdf.C[4].value 100    ## assign to it         
        """
        return () if not self.extended else self.component_getter ( self.__C )     
    @C.setter
    def C (  self , value ) : 
        assert self.extended, "'C'' cannot be set for non-extended model!"
        self.component_setter ( self.__C , value )

    @property
    def yields ( self ) :
        """'yields'' : yields for 'all'' components, same as 'S+B+C'', emtpy for non-extended fits"""
        return () if not self.extended else self.__S + self.__B + self.__C 

    @property 
    def F ( self ) :
        """Get fit fractions for non-expended fits (empty for extended fits)
        For single fraction (2 fit components):
        >>> print pdf.F           ## read the single fraction 
        >>> pdf.F = 0.1           ## assign to it 
        For multiple fractions (>2 fit components):
        >>> print pdf.F[4]        ## read the 4th fraction
        >>> pdf.F = (0.1,0.2,0.3,0.4,0.6) ## assign to it 
        ... or, alternatively:
        >>> print pdf.F[4]        ## read the 4th fraction
        >>> pdf.F[4].value = 0.1  ## assign to it         
        """
        return () if self.extended else self.component_getter ( self.__F )     
    @F.setter
    def F (  self , value ) :
        assert not self.extended, "'F'' cannot be set for extended model!"        
        self.component_setter ( self.__F , value )

    @property
    def natural ( self ) :
        """Are all yields natural? """
        if not self.yeilds : return False
        for y in self.yields :
            if not isinstance  ( y , ROOT.RooRealVar ) : return False 
        return True 
        
    def total_yield ( self ) :
        """'total_yield''' : get the total yield if/when possible"""
        if not self.extended                                   : return None 
        if not self.fit_result                                 : return None
        if not valid_pointer ( self.fit_result )               : return None
        yields = self.yields
        if not yields                                          : return None
        if not self.natural                                    : return None 
        if 1 ==  len ( yields )                                : return yields[0].value  
        return self.fit_result.sum ( *yields )
    
    @property
    def fix_norm ( self ) :
        """`fix-norm`: 
        - see `ROOT.RooAbsPdf.SetCoefNormalization`
        - see `ROOT.RooAbsPdf.getCoefNormalization`
        """
        pars = self.pdf.getCoefNormalization()
        return True if pars else False 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END 
# =============================================================================
