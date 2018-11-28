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
import ROOT
from   ostap.fitting.basic import PDF, Generic1D_pdf  
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
                   resolution        ,   ## the    resolution
                   useFFT  = True    ,   ## use FFT ? 
                   nbins   = 50000   ,   ##  #bins for FFT
                   buffer  = 0.25    ,   ## buffer fraction ## setBufferFraction
                   nsigmas = 6       ) : ## number of sigmas for setConvolutionWindow

        ## the axis 
        assert isinstance ( xvar , ROOT.RooAbsReal ) , "``xvar'' must be ROOT.RooAbsReal"
        self.__xvar   = xvar
        self.__useFFT = True if  useFFT else False
        
        ## pdf itself 
        if   isinstance ( pdf , PDF            ) : self.__old_pdf = pdf
        elif isinstance ( pdf , ROOT.RooAbsPdf ) :
            self.__old_pdf = Generic1D_pdf ( pdf , xvar = self.__xvar )
        else :
            raise AttributeError("Convolution: invalid ``pdf'' %s/%s"  % ( pdf , type ( pdf ) ) )
        
        ## resolution  function 
        if   isinstance ( resolution , PDF            ) : self.__resolution = resolution
        elif isinstance ( resolution , ROOT.RooAbsPdf ) :
            self.__resolution = Generic1D_pdf ( resolution , xvar = self.__xvar ) 
        else :
            ## use   Gaussial resolution
            import ostap.fitting.resolution as OFR 
            self.__resolution = OFR.ResoGauss ( 'ResoGauss' + name ,
                                                self.__xvar        ,
                                                sigma = resolution ,
                                                mean  = None       )

        self.__nbins   = nbins
        self.__buffer  = buffer
        self.__nsigmas = nsigmas
        
        if self.useFFT : ## Use Fast Fourier transform  (fast)
            
            assert isinstance ( nbins  , (long,int) ) and 100  < nbins          , \
                   "Invalid ``nbins''  parameter  %s/%s for fast Fourier transform"  % ( nbins  , type ( nbins  ) )
            assert isinstance ( buffer ,  float     ) and 0.1  < buffer < 0.9   , \
                   "Invalid ``buffer'' parameter  %s/%s for ``setBufferFraction''"     % ( buffer , type ( buffer ) )
                   
            if hasattr ( self.__resolution , 'sigma' ) and self.__xvar.minmax() : 
                mn , mx = xvar.minmax()
                dm  = mx - mn
                sv  = self.__resolution.sigma.getVal() 
                dm /= sv 
                self.__nbins  = max ( self.__nbins , 100 * int ( dm ) )               
                logger.debug('Convolution: choose #bins %d' % self.__nbins ) 
                
            self.__xvar.setBins ( self.__nbins , 'cache' )
            
            self.__pdf = ROOT.RooFFTConvPdf (
                'FFT'     + name       , 
                'FFT(%s)' % name       ,
                self.__xvar            ,
                self.__old_pdf    .pdf ,
                self.__resolution .pdf )            
            self.__pdf.setBufferFraction ( 0.25 )
            
        else :           ##  Use plain numerical integration (could be slow)

            
            assert isinstance ( nsigmas  , ( long , int , float ) ) and 2 < nsigmas , \
                   "Invalid ``nsigmas''  parameter  %s/%s for ``setConvolutionWindow''"  % ( nsigmas , type ( nsigmas ) )
            
            self.__pdf = ROOT.RooNumConvPdf (
                'CNV'     + name       ,
                'CNV(%s)' % name       ,
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
        """``resoltuion'': pdf for resolution function"""
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
        return self.__nbins 

    @property
    def buffer ( self ) :
        """``buffer'' : buffer fraction for Fast Fourier Transform"""
        return self.__buffer

    @property
    def nsigmas ( self ) :
        """``nsigmas'' : convolution window for RooNumConvPdf"""
        return self.__nsigmas
        
# =============================================================================
## @class Convolution_pdf
#  Helper class to simplify the convolutions
#  @code
#  pdf = ...
#  pdfc = Convolution_pdf(  'C' , pdf  , xvar = ... , resolution = ... , useFFT = True )
#  @endcode
class Convolution_pdf(PDF) :
    """Helper class/PDF to simplify the convolution:
    >>> pdf = ...
    >>> pdfc = Convolution_pdf(  'C' , pdf  , xvar = ... , resolution = ... , useFFT = True )
    """
    def __init__ ( self              ,
                   name              ,   ## the name 
                   pdf               ,   ## the PDF to be convoluted 
                   resolution        ,   ## the resolution
                   xvar    = None    ,   ## the axis varable
                   useFFT  = True    ,   ## use  FastFourierTransform?
                   nbins   = 1000000 ,   ## #bins for FFT
                   buffer  = 0.25    ,   ## buffer fraction ## setBufferFraction
                   nsigmas = 6       ) : ## number of sigmas for setConvolutionWindow

        if   isinstance ( pdf , PDF ) :
            xvar = pdf.xvar
        elif isinstance ( pdf , ROOT.RooAbsPdf ) and isinstance ( xvar , ROOT.RooAbsReal ) :
            pdf  = Generic1D_pdf  ( pdf , xvar )
        else :
            raise AttributeError ("Convolution_pdf: invalid pdf/xvar    %s/%s"  % ( pdf , xvar ) ) 

        ## initialize the base 
        PDF.__init__ ( self , name , xvar )

        ## make the actual convolution
        self.__cnv = Convolution ( name       = name       ,
                                   pdf        = pdf        ,
                                   xvar       = xvar       ,
                                   resolution = resolution ,
                                   useFFT     = useFFT     ,
                                   nbins      = nbins      ,
                                   buffer     = buffer     ,
                                   nsigmas    = nsigmas    )

        self.pdf = self.__cnv.pdf 

        ## save configuration 
        self.config = {
            'name'       : self.name           ,
            'xvar'       : self.xvar           ,
            'pdf'        : self.cnv.old_pdf    ,
            'resolution' : self.cnv.resolution ,
            'useFFT'     : self.cnv.useFFT     ,
            'nbins'      : self.cnv.nbinsFFT   ,
            'buffer'     : self.cnv.buffer     ,
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
        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
