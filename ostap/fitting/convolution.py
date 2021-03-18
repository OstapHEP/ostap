#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/convolution.py
#  Set of useful basic utilities to use convolution 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
# =============================================================================
"""Set of useful basic utilities to use convoltuion"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'Convolution'      , ## helper utility to build convolution
    'Convolution_pdf'  , ## ``ready-to-use'' PDF for convolution 
    )
# =============================================================================
import ROOT, math
from   ostap.fitting.basic import PDF, Generic1D_pdf
from   ostap.core.ostap_types    import num_types ,  integer_types 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.convolution' )
else                       : logger = getLogger ( __name__         )
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
    """Helper class to make a convolution PDF:
    
    >>> pdf  = .... ##  original PDF (the one from Ostap or bare ROOT.RooAbsPdf) 
    >>> resolution = 3 * MeV                                     ## fixed number ``sigma''
    >>> # resolution = ROOT.RooRealVar( 'sigma','',2*MeV,4*MeV)  ## ``sigma'' as ROOT.RooAbsReal 
    >>> # resolution = ResoGauss_pdf ( ... )                     ## pdf from Ostap
    >>> # resolution = ...                                       ## bare ROOT.RooAbsPdf
    >>> cnv = Convolution ('CNV' , pdf , xvar  =  xvar , resolution = resolution )
    >>> cnv_pdf = cnv.pdf 
    """
    def __init__ ( self              ,
                   name              ,
                   pdf               ,   ## the PDF to be convoluted 
                   xvar              ,   ## the axis variable
                   resolution        ,   ## the resolution
                   useFFT   = True   ,   ## use FFT ? 
                   nbins    = 10000  ,   ## number of bins for FFT
                   buffer   = 0.25   ,   ## buffer fraction use for setBufferFraction
                   bufstrat = None   ,   ## "Buffer strategy" : (0,1,2)
                   nsigmas  = 6      ) : ## number of sigmas for setConvolutionWindow

        ## the axis 
        assert isinstance ( xvar , ROOT.RooAbsReal ) , "``xvar'' must be ROOT.RooAbsReal"
        
        self.__xvar   = xvar
        self.__useFFT = True if useFFT else False

        self.__arg_pdf        = pdf
        self.__arg_resolution = resolution
        self.__arg_xvar       = xvar
        
        ## pdf itself 
        if   isinstance ( pdf , PDF            ) : self.__old_pdf = pdf
        elif isinstance ( pdf , ROOT.RooAbsPdf ) :
            self.__old_pdf = Generic1D_pdf ( pdf , xvar = self.__xvar )
        else :
            raise AttributeError("Convolution: invalid ``pdf'' %s/%s"  % ( pdf , type ( pdf ) ) )
        
        ## resolution  function 
        if   isinstance ( resolution , PDF            ) :
            self.__resolution = resolution
        elif isinstance ( resolution , ROOT.RooAbsPdf ) :
            self.__resolution = Generic1D_pdf ( resolution , xvar = self.__xvar ) 
        else :
            ## use   Gaussial resolution
            import ostap.fitting.resolution as OFR
            rname  = ( 'Reso%s_' % name ) if name else 'ResoGauss_'
            rname  = PDF.generate_name ( prefix = rname ) 
            self.__resolution = OFR.ResoGauss ( rname              ,
                                                self.__xvar        ,
                                                sigma = resolution ,
                                                mean  = None       )
        self.__nbins    = nbins
        self.__buffer   = buffer
        self.__bufstrat = bufstrat 
        self.__nsigmas  = nsigmas

        name = name if name else PDF.generate_name ( prefix = 'cnv_%s@%s' % ( pdf.name , self.resolution.name ) )
        
        if self.useFFT : ## Use Fast Fourier transform  (fast)
            
            assert isinstance ( nbins  , integer_types ) and 500   <  abs ( nbins  )  , \
                   "Invalid ``nbins''  parameter %s/%s for fast Fourier transform"  % ( nbins  , type ( nbins  ) )
            assert isinstance ( buffer ,  float        ) and 0.03  <= buffer <=0.9    , \
                   "Invalid ``buffer'' parameter %s/%s for ``setBufferFraction''"   % ( buffer , type ( buffer ) )

            ## adjust #bins if positive. keep it as it is if negavtive 
            if hasattr ( self.__resolution , 'sigma' ) and self.__xvar.minmax() and  self.__nbins > 0 :
                mn , mx = self.xvar.minmax()
                dm  = mx - mn
                sv  = self.__resolution.sigma.getVal() 
                dm /= sv
                nb  = min ( 50 * ( int ( dm ) + 1  ) , 2**14 )
                nb  = 2**math.frexp(nb)[1]
                if nb > self.nbinsFFT : 
                    self.__nbins  = nb       
                    logger.info('Convolution: choose #bins %d' % self.__nbins )

            self.__xvar.setBins ( self.nbinsFFT , 'cache' )

            pdf_name  = PDF.generate_name ( prefix = 'FFT_' + name ) 
            pdf_title = 'FFT convolution: %s (*) %s' %  ( pdf.name , self.resolution.name )
            self.__pdf = ROOT.RooFFTConvPdf (
                pdf_name               ,
                pdf_title              ,
                self.__xvar            ,
                self.__old_pdf    .pdf ,
                self.__resolution .pdf )            
            self.__pdf.setBufferFraction ( self.buffer )
            
            if isinstance ( self.bufstrat , int ) and 0 <= self.bufstrat <= 2 : 
                self.__pdf.setBufferStrategy ( self.bufstrat )
                
        else :           ##  Use plain numerical integration (could be slow)
            
            assert isinstance ( nsigmas  , num_types ) and 2.5 <= nsigmas , \
                   "Invalid ``nsigmas''  parameter  %s/%s for ``setConvolutionWindow''"  % ( nsigmas , type ( nsigmas ) )
            
            pdf_name  = PDF.generate_name ( prefix = 'NUMCNV_' + name ) 
            pdf_title = 'NUM convolution: %s (*) %s' %  ( pdf.name , self.resolution.name )
            self.__pdf = ROOT.RooNumConvPdf (
                pdf_name               ,
                pdf_title              ,
                self.__xvar            ,
                self.__old_pdf    .pdf ,
                self.__resolution .pdf )
            
            if hasattr ( self.__resolution , 'sigma' ) :                
                if hasattr ( self.__resolution , 'mean' ) :
                    self.__pdf.setConvolutionWindow ( self.__resolution.mean  ,
                                                      self.__resolution.sigma , self.__nsigmas  )
                    logger.debug('Convolution: choose window of %s' % self.__nsigmas ) 
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
        """``resolution'': pdf for resolution function"""
        return self.__resolution 
    @property
    def old_pdf ( self ) :
        """``old'' - the original pdf before convolution"""
        return self.__old_pdf
    @property
    def pdf ( self ) :
        """``new'' (convoluted) PDF"""
        return self.__pdf
    @property
    def nbinsFFT ( self ) :
        """number of cache bins for Fast Fourier Transform"""
        return abs ( self.__nbins )
    @property
    def nbins    ( self ) :
        """number of cache bins for Fast Fourier Transform"""
        return abs ( self.__nbins ) 
    @property
    def buffer ( self ) :
        """``buffer'' : buffer fraction for Fast Fourier Transform"""
        return self.__buffer
    @property
    def bufstrat ( self ) :
        """``bufstrat'' : buffer strategy:
        - 'Extend/0' means is that the input p.d.f convolution observable range is widened to include the buffer range
        - 'Flat/1'   means that the buffer is filled with the p.d.f. value at the boundary of the observable range
        - 'Mirror/2' means that the buffer is filled with a mirror image of the p.d.f. around the convolution observable boundary
        """
        return self.__bufstrat
    @property
    def nsigmas ( self ) :
        """``nsigmas'' : convolution window for RooNumConvPdf"""
        return self.__nsigmas
    @property
    def name    ( self ) :
        """``name'' : name of this convoltuoon object/name of pdf"""
        return self.__pdf.name
    
# =============================================================================
## @class Convolution_pdf
#  Helper class to simplify the convolutions
#  @code
#  pdf = ...
#  pdfc = Convolution_pdf( pdf  , xvar = ... , resolution = ... , useFFT = True )
#  @endcode
class Convolution_pdf(PDF) :
    """Helper class/PDF to simplify the convolution:
    >>> pdf = ...
    >>> pdfc = Convolution_pdf( pdf  , xvar = ... , resolution = ... , useFFT = True )
    """
    def __init__ ( self             ,
                   pdf              ,   ## the PDF to be convoluted 
                   resolution       ,   ## the convolution/resolution
                   xvar     = None  ,   ## the axis varable
                   useFFT   = True  ,   ## use  FastFourierTransform?
                   nbins    = 2**14 ,   ## #bins for FFT
                   buffer   = 0.25  ,   ## buffer fraction ## setBufferFraction
                   bufstrat = None  ,   ## "Buffer strategy" : (0,1,2)
                   nsigmas  = 6     ,   ## number of sigmas for setConvolutionWindow
                   name     = ''    ) : ## the name 

        self.__arg_pdf        = pdf
        self.__arg_resolution = resolution 
        self.__arg_xvar       = xvar
        
        if   isinstance ( pdf , PDF ) :
            self.__old_pdf = pdf 
            xvar = pdf.xvar
        elif isinstance ( pdf , ROOT.RooAbsPdf ) and xvar :
            self.__old_pdf = Generic1D_pdf  ( pdf , xvar )
        else :
            raise AttributeError ("Convolution_pdf: invalid pdf/xvar %s/%s"  % ( pdf , xvar ) ) 

        em = pdf.pdf.extendMode()
        if   1 == em : self.warning ( "PDF  ``canBeExtended''" )
        elif 2 == em : self.error   ( "PDF ``mustBeExtended''" )

        ## make the actual convolution
        if isinstance ( resolution , Convolution ) :
            assert resolution.xvar is xvar, "Mismatch in ``xvar'': %s vs %s" % ( xvar , resolution.xvar )
            self.__cnv = resolution
        else :
            self.__cnv = Convolution ( name       = ''               ,
                                       pdf        = self.old_pdf.pdf ,
                                       xvar       = xvar             ,
                                       resolution = resolution       ,
                                       useFFT     = useFFT           ,
                                       nbins      = nbins            ,
                                       buffer     = buffer           ,
                                       bufstrat   = bufstrat         ,
                                       nsigmas    = nsigmas          )

        name = name if name else self.generate_name ( prefix = 'Cnv_%s@%s_' %  ( pdf.name , self.resolution.name ) ) 
                            
        ## initialize the base 
        PDF.__init__ ( self , name , xvar )


        ## the  actual convoluted PDF 
        self.pdf = self.__cnv.pdf 

        ## save configuration 
        self.config = {
            'name'       : self.name           ,
            'xvar'       : self.xvar           ,
            'pdf'        : self.old_pdf        ,
            'resolution' : self.cnv.resolution ,
            'useFFT'     : self.cnv.useFFT     ,
            'nbins'      : self.cnv.nbinsFFT   ,
            'buffer'     : self.cnv.buffer     ,
            'bufstrat'   : self.cnv.bufstrat   ,
            'nsigmas'    : self.cnv.nsigmas    ,
            }

    @property
    def convolution ( self ) :
        """``convolution'' : the actual convolution object (same as ``cnv'')"""
        return self.__cnv
    @property
    def cnv         ( self ) :
        """``cnv'' : the actual convolution object (same as ``convolution'')"""
        return self.__cnv
    @property
    def old_pdf ( self ):
        """``old_pdf''  : original (non-convolved) PDF"""
        return self.__old_pdf
    @property
    def original_pdf ( self ):
        """``original_pdf''  : original (non-convolved) PDF"""
        return self.__old_pdf
    
    @property
    def resolution ( self ) :
        """``resolution'' :  the actual resolution function/PDF"""
        return self.cnv.resolution
    

    ## ## redirect any other attributes to original PDF
    ## def __getattr__ ( self , attr ) :
    ##     """Get all extra attributes from the original PDF"""

    ##     try :
    ##         return getattr ( self.old_pdf     , attr )
    ##     except AttributeError :
    ##         pass

    ##     try :
    ##         return getattr ( self.convolution , attr )
    ##     except AttributeError :
    ##         pass
    ##     raise AttributeError("Unknown attribute %s:" % attr )
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================

