#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/convolution.py
#  Set of useful basic utilities to use convolution 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
# =============================================================================
""" Set of useful basic utilities to use convoluion"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'ConvolutionConfig', ## context manager  to provide convolution parmaeters 
    'Convolution'      , ## helper utility to build convolution
    'Convolution_pdf'  , ## ``ready-to-use'' PDF for convolution 
    )
# =============================================================================
from   ostap.fitting.pdfbasic  import PDF1, Generic1D_pdf
from   ostap.core.ostap_types  import num_types, integer_types, string_types 
from   ostap.fitting.rooreduce import root_store_factory
from   ostap.utils.basic       import typename 
import ostap.logger.table      as     T 
import ROOT, math
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.convolution' )
else                       : logger = getLogger ( __name__         )
# =============================================================================
## helper class to keep convolution configuration
class CnvConfig(object) :
    """ Helper class to keep convolution configuration"""
    CONFIG = {}
    ## retrive (copy) of current configuration
    @classmethod 
    def config ( klass ) :
        """ Return currrent configuration
        """
        return klass.CONFIG
    # ========================================================================
    @classmethod
    def clean  ( klass ) :
        ## clean the global configuraiton 
        while klass.CONFIG : klass.CONFIG.popitem ()
        
# ============================================================================
## Context manager to treat the convoluton nconfiguration
class ConvolutionConfig(object):
    """ Context manager to treat the convolution configuration
    """
    def __init__ ( self , **config ) :
        self.__config = config
        self.__backup = {} 
    # ========================================================================
    ## context manager 
    def __enter__ ( self ) :
        """ ENTER: update global convolution configuration 
        """
        self.__backup.update ( CnvConfig.config() )
        CnvConfig.clean () 
        CnvConfig.CONFIG.update ( self.config )                                   
        return self
    # ========================================================================
    ## context manager 
    def __exit__ ( self , *_ ) : 
        """ EXIT: restore global convolution configuration 
        """
        ## restore previous configration
        CnvConfig.clean () 
        CnvConfig.CONFIG.update ( self.backup )                                   
    @property
    def config ( self ) :
        """`config`: new configuration parameters for convolution"""
        return self.__config
    @property
    def backup ( self ) :
        """`backuo`: old configuration parameters for convolution"""
        return self.__backup
    
# =============================================================================
## @class Convolution
#  Helper class to perform the convolution
#  @code
#  >>> pdf  = .... ##  original PDF (the one from Ostap or bare ROOT.RooAbsPdf) 
#  >>> resolution = 3 * MeV                                     ## fixed number ``sigma''
#  >>> # resolution = ROOT.RooRealVar( 'sigma','',2*MeV,4*MeV)  ## ``sigma'' as ROOT.RooAbsReal 
#  >>> # resolution = ResoGauss_pdf ( ... )                     ## pdf from Ostap
#  >>> # resolution = ...                                       ## bare ROOT.RooAbsPdf
#  >>> cnv = Convolution ('CNV' , pdf , xvar  =  xvar , resolution = resolution )
#  >>> cnv_pdf = cnv.pdf 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-07-13
class Convolution(object):
    """ Helper class to make a convolution PDF:
    
    >>> pdf  = .... ##  original PDF (the one from Ostap or bare ROOT.RooAbsPdf) 
    >>> resolution = 3 * MeV                                     ## fixed number ``sigma''
    >>> # resolution = ROOT.RooRealVar( 'sigma','',2*MeV,4*MeV)  ## ``sigma'' as ROOT.RooAbsReal 
    >>> # resolution = ResoGauss_pdf ( ... )                     ## pdf from Ostap
    >>> # resolution = ...                                       ## bare ROOT.RooAbsPdf
    >>> cnv = Convolution ('CNV' , pdf , xvar  =  xvar , resolution = resolution )
    >>> cnv_pdf = cnv.pdf 
    """
    def __init__ ( self                ,
                   pdf                 ,  ## the PDF to be convoluted 
                   xvar                ,  ## the axis variable
                   resolution          ,  ## the resolution
                   name     = ''       ,  ## name for Gaussian resolution pdf 
                   useFFT   = True     ,  ## use FFT ? 
                   nbins    = 10000    ,  ## number of bins for FFT
                   buffer   = 0.25     ,  ## buffer fraction use for setBufferFraction
                   bufstrat = 'Extend' ,  ## "Buffer strategy" : (0,1,2)
                   shift1   = None     ,  ## shift1 parameter
                   shift2   = None     ,  ## shift2 parameter
                   nsigmas  = 6        ,  ## number of sigmas for setConvolutionWindow
                   silent   = True     ) : 
        
        ## the axis 
        assert isinstance ( xvar , ROOT.RooAbsReal ) or not xvar , "`xvar' must be ROOT.RooAbsReal"
        
        self.__xvar   = xvar
        self.__useFFT = True if useFFT else False

        self.__arg_pdf        = pdf
        self.__arg_resolution = resolution
        self.__arg_xvar       = xvar
        
        ## pdf itself 
        if   isinstance ( pdf , PDF1           ) :
            self.__old_pdf = pdf
            assert not self.xvar or xvar == pdf.xvar , "Invalid 'xvar/pdf' setting!"
            self.__xvar = pdf.xvar             
        elif isinstance ( pdf , ROOT.RooAbsPdf ) and self.xvar :
            self.__old_pdf = Generic1D_pdf ( pdf , xvar = self.__xvar )
        else :
            raise TypeError("Convolution: invalid ``pdf'' %s/%s"  % ( pdf , type ( pdf ) ) )

        ## resolution  function 
        if   isinstance ( resolution , PDF1           ) :
            assert resolution.xvar == self.xvar, "Invalid 'xvar/resolution' setting!"
            self.__resolution = resolution
        elif isinstance ( resolution , ROOT.RooAbsPdf ) :
            self.__resolution = Generic1D_pdf ( resolution , xvar = self.__xvar ) 
        else :
            ## use   Gaussian resolution
            import ostap.fitting.resolution as OFR
            rname  = name if name else 'ResoGauss'
            rname  = self.old_pdf.generate_name ( prefix = rname )
            self.__resolution = OFR.ResoGauss ( rname              ,
                                                self.__xvar        ,
                                                sigma = resolution )

        self.__nbins    = nbins
        self.__buffer   = buffer
        self.__bufstrat = bufstrat 
        self.__nsigmas  = nsigmas
        self.__shift1   = shift1 
        self.__shift2   = shift2

        name = name if name else PDF1.generate_name ( prefix = 'cnv_%s@%s' % ( pdf.name , self.resolution.name ) )
        self.__name = name

        # =====================================================================
        if self.useFFT : ## Use Fast Fourier transform  (fast) # ==============
            # =================================================================

            bs = self.__bufstrat
            if isinstance   ( bs , string_types ) :
                bs = bs.lower ()
                if   'extend' == bs : self.__bufstrat = ROOT.RooFFTConvPdf.Extend
                elif 'flat'   == bs : self.__bufstrat = ROOT.RooFFTConvPdf.Flat
                elif 'mirror' == bs : self.__bufstrat = ROOT.RooFFTConvPdf.Mirror
                else :
                    logger.error ( "Unknown bufstrat %s" % bs )
                    self.__bufstrat = ROOT.RooFFTConvPdf.Extend
            elif isinstance   ( bs , integer_types ) :
                if not bs in ( ROOT.RooFFTConvPdf.Extend ,
                               ROOT.RooFFTConvPdf.Flat   ,
                               ROOT.RooFFTConvPdf.Mirror ) :
                    logger.error ( "Unknown bufstrat %s" % bs )
                    self.__bufstrat = ROOT.RooFFTConvPdf.Extend
            else :
                logger.error ( "Unknown bufstrat %s" % bs )
                self.__bufstrat = ROOT.RooFFTConvPdf.Extend

            assert isinstance ( nbins  , integer_types ) and 100  <= abs ( nbins  )  , \
                   "Invalid `nbins'  parameter %s/%s for fast Fourier transform"  % ( nbins  , type ( nbins  ) )
            assert isinstance ( buffer , float         ) and 0.05 <= buffer <= 0.95  , \
                   "Invalid `buffer' parameter %s/%s for `setBufferFraction'"   % ( buffer , type ( buffer ) )

            ## adjust #bins if positive. keep it as it is if negative 
            if hasattr ( self.__resolution , 'sigma' ) and self.__xvar.minmax() and  0 < self.__nbins :
                mn , mx = self.xvar.minmax()
                dm  = mx - mn
                sv  = self.__resolution.sigma.getVal() 
                dm /= sv
                nb  = min ( 50 * ( int ( dm ) + 1  ) , 2**14 )
                nb  = 2 ** math.frexp(nb)[1]
                if nb > self.nbinsFFT : 
                    self.__nbins  = nb       
                    logger.info('Convolution: choose #bins %d' % self.__nbins )

            self.__xvar.setBins ( self.nbinsFFT , 'cache' )

            self.__pdf = ROOT.RooFFTConvPdf (
                self.old_pdf.new_roo_name ( 'fft' ) ,
                'FFT convolution: %s (*) %s' %  ( pdf.name , self.resolution.name ) ,
                self.__xvar              ,
                self.__old_pdf    .pdf   ,
                self.__resolution .pdf   )            
            self.__pdf.setBufferFraction ( self.buffer )

            ## buffer strategy 
            self.__pdf.setBufferStrategy ( self.__bufstrat )

            ## set shift-parameters
            if ( not self.shift1 is None ) and ( not self.shift2 is None ) :             
                self.__pdf.setShift ( self.shift1 , self.shift2 )

            # =================================================================
        else : ##  Use plain numerical integration (could be VERY slow) # =====
            # =================================================================
            
            assert isinstance ( nsigmas  , num_types ) and 2.5 <= nsigmas , \
                   "Invalid `nsigmas'  parameter  %s/%s for `setConvolutionWindow'"  % ( nsigmas , type ( nsigmas ) )
            
            self.__pdf = ROOT.RooNumConvPdf (
                self.old_pdf.new_roo_name ( 'numcnv' ) ,
                'NUM convolution: %s (*) %s' %  ( pdf.name , self.resolution.name ) ,
                self.__xvar            ,
                self.__old_pdf    .pdf ,
                self.__resolution .pdf )
            
            if hasattr ( self.__resolution , 'sigma' ) :                
                if hasattr ( self.__resolution , 'mean' ) :
                    self.__pdf.setConvolutionWindow ( self.__resolution.mean  ,
                                                      self.__resolution.sigma , self.__nsigmas  )
                    logger.debug('Convolution: choose window of %s' % self.__nsigmas )

        if not silent :
            title = 'Convolution'
            logger.info ( 'Convolution configuration:\n%s' % self.table ( title = title , prefix = '# ' ) ) 
                        
    @property
    def name ( self ) :
        """'name' : the bname of convolution object"""
        return self.__name    
    @property
    def xvar (self ) :
        """The axis variable for  convolution"""
        return self.__xvar
    @property
    def useFFT  ( self ) :
        """``useFFT'' :    use Fast Fourier Transform?"""
        return self.__useFFT
    @property
    def resolution ( self  ) :
        """`resolution': pdf for resolution function"""
        return self.__resolution 
    @property
    def old_pdf ( self ) :
        """`old' - the original pdf before convolution"""
        return self.__old_pdf
    @property
    def pdf ( self ) :
        """`new' (convoluted) PDF"""
        return self.__pdf
    @property
    def nbinsFFT ( self ) :
        """`nbinsFFT`: number of cache bins for Fast Fourier Transform"""
        return abs ( self.__nbins )
    @property
    def nbins    ( self ) :
        """`nbins`: number of cache bins for Fast Fourier Transform"""
        return abs ( self.__nbins ) 
    @property
    def buffer ( self ) :
        """`buffer' : buffer fraction for Fast Fourier Transform"""
        return self.__buffer
    @property
    def bufstrat ( self ) :
        """`bufstrat' : buffer strategy:
        - 'Extend/0' means is that the input p.d.f convolution observable range is widened to include the buffer range
        - 'Flat/1'   means that the buffer is filled with the p.d.f. value at the boundary of the observable range
        - 'Mirror/2' means that the buffer is filled with a mirror image of the p.d.f. around the convolution observable boundary
        """
        return ( 'Extend'     if ROOT.RooFFTConvPdf.Extend == self.__bufstrat else
                 ( 'flat'     if ROOT.RooFFTConvPdf.Flat   == self.__bufstrat else
                   ( 'Mirror' if ROOT.RooFFTConvPdf.Mirror == self.__bufstrat else '????' ) ) ) 
    @property
    def shift1 ( self ) :
        """`shift1` : parameter for RooFFTConvPdf"""
        return self.__shift1
    @property
    def shift2 ( self ) :
        """`shift2` : parameter for RooFFTConvPdf"""
        return self.__shift2    
    @property
    def nsigmas ( self ) :
        """`nsigmas' : convolution window for RooNumConvPdf"""
        return self.__nsigmas
    @property
    def name    ( self ) :
        """`name' : name of this convoltuoon object/name of pdf"""
        return self.__pdf.name
    def __str__ ( self ) :
        part1 = "pdf=%s,resolution=%s"                      % ( self.old_pdf    ,
                                                                self.resolution ) 
        part2  = "useFFT=%s,nbins=%s,buffer=%s,bufstrat=%s" % ( self.useFFt     , 
                                                                self.nbins      ,
                                                                self.buffer     , 
                                                                self.bufstrat   )  
        return "Convolution(%s,%s,nsigmas=%s)" %  ( part1 , part2 , self.nsigmas)
    __repr__ = __str__

    ## serialize the convolution object object 
    def __reduce__ ( self ) :
        """ Serialize the convolution object object"""
        return root_store_factory  , ( type ( self )    ,
                                       self.old_pdf     , 
                                       self.xvar        , 
                                       self.resolution  ,
                                       self.name        , 
                                       self.useFFT      ,
                                       self.nbinsFFT    ,
                                       self.buffer      ,
                                       self.bufstrat    ,
                                       self.shift1      ,
                                       self.shift2      ,
                                       self.nsigmas     )
    # =========================================================================
    ## Get the convolution result as table 
    def table ( self , title = '' , prefix = '' ) :
        """ Get the convolution configuration as table """
        
        rows = [ ( 'Parameter' , 'value' ) ]        
        row  = 'name'       , '%s' % self.name
        rows.append ( row )        
        row  = 'pdf'        , '%s' % typename ( self.old_pdf )
        rows.append ( row )        
        row  = 'xvar'       , '%s' % self.xvar.name             
        rows.append ( row )        
        row  = 'resolution' , '%s' % typename ( self.resolution )
        rows.append ( row )        
        row  = 'use FFT?'   , '%s' % self.useFFT 
        rows.append ( row )
        if self.useFFT : 
            row  = 'nbins/cache' , '%s' % self.nbinsFFT 
            rows.append ( row )        
            row  = 'buffer'   , '%s' % self.buffer
            rows.append ( row )
            if   self.__bufstrat == ROOT.RooFFTConvPdf.Extend : bs = 'Extend'
            elif self.__bufstrat == ROOT.RooFFTConvPdf.Flat   : bs = 'Flat'
            elif self.__bufstrat == ROOT.RooFFTConvPdf.Mirror : bs = 'Mirror'
            else : bs = self.bufstrat 
            ## row  = 'bufstrat' , '%s' % self.bufstrat
            row  = 'bufstrat' , '%s' % bs 
            rows.append ( row )        
            row  = 'shift1'   , '%s' % self.shift1
            rows.append ( row )        
            row  = 'shift2'   , '%s' % self.shift2
            rows.append ( row )
        else :
            row  = 'nsigmas' , '%s' % self.sigmas
            rows.append ( row )
            
        title = title if title else 'Convolution config'
        return T.table ( rows , title = title , prefix = prefix , alignment = 'll' )

    __repr__ = table 

        
# =============================================================================
## @class Convolution_pdf
#  Helper class to simplify the convolutions
#  @code
#  pdf = ...
#  pdfc = Convolution_pdf( pdf  , xvar = ... , resolution = ... , useFFT = True )
#  @endcode
class Convolution_pdf(PDF1) :
    """ Helper class/PDF to simplify the convolution:
    >>> pdf = ...
    >>> pdfc = Convolution_pdf( pdf  , xvar = ... , resolution = ... , useFFT = True )
    """
    def __init__ ( self             ,
                   pdf              ,   ## the PDF to be convoluted 
                   resolution       ,   ## the convolution/resolution
                   xvar     = None  ,   ## the axis varable
                   name     = ''    ,   ## name
                   **kwargs         ) : ## convolution configuration 
        
        self.__arg_pdf        = pdf
        self.__arg_resolution = resolution 
        self.__arg_xvar       = xvar
        
        if   isinstance ( pdf , PDF1 ) :
            self.__old_pdf = pdf 
            xvar           = pdf.xvar
        elif isinstance ( pdf , ROOT.RooAbsPdf ) and xvar :
            self.__old_pdf = Generic1D_pdf  ( pdf , xvar )
        else :
            raise TypeError ("Convolution_pdf: invalid pdf/xvar %s/%s"  % ( pdf , xvar ) ) 

        em = self.old_pdf.pdf.extendMode ()
        if   1 == em : self.warning ( "PDF  ``canBeExtended''" )
        elif 2 == em : self.error   ( "PDF ``mustBeExtended''" )

        ## make the actual convolution
        if isinstance ( resolution , Convolution ) :
            assert resolution.xvar is xvar, "Mismatch in ``xvar'': %s vs %s" % ( xvar , resolution.xvar )
            self.__cnv = resolution
        else :
            self.__cnv = Convolution ( pdf        = self.old_pdf ,
                                       xvar       = xvar         ,
                                       resolution = resolution   , **kwargs ) 

        name = name if name else self.generate_name ( prefix = '(%s)@(%s)' % ( pdf.name , self.resolution.name ) ) 
                            
        ## initialize the base 
        PDF1.__init__ ( self , name , xvar = xvar )


        ## the  actual convoluted PDF 
        self.pdf = self.__cnv.pdf 

        ## save configuration 
        self.config = {
            'name'       : self.name              ,
            'xvar'       : self.xvar              ,
            'pdf'        : self.old_pdf           ,
            'resolution' : self.__arg_resolution  , ## attention! 
            'useFFT'     : self.cnv.useFFT        ,
            'nbins'      : self.cnv.nbins         ,
            'buffer'     : self.cnv.buffer        ,
            'bufstrat'   : self.cnv.bufstrat      ,
            'nsigmas'    : self.cnv.nsigmas       }

    @property
    def convolution ( self ) :
        """`convolution' : the actual convolution object (same as ``cnv'')"""
        return self.__cnv
    @property
    def cnv         ( self ) :
        """`cnv' : the actual convolution object (same as ``convolution'')"""
        return self.__cnv
    @property
    def old_pdf ( self ):
        """`old_pdf'  : original (non-convolved) PDF"""
        return self.__old_pdf
    @property
    def original_pdf ( self ):
        """`original_pdf'  : original (non-convolved) PDF"""
        return self.__old_pdf
    
    @property
    def resolution ( self ) :
        """`resolution' :  the actual resolution function/PDF"""
        return self.cnv.resolution
        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================

