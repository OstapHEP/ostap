#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostap/fitting/morphing_pdf.py
#  Morphing PDF 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2020-07-22
# =============================================================================
"""Morphing PDF 
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2020-07-22"
__all__     = (
    ##
    'Morphing1D_pdf' , ## 1D-morphing/1D PDF
    'Morphing2D_pdf' , ## 1D-morphing/2D PDF
    ##
    )
# =============================================================================
import ROOT
# =============================================================================
from   ostap.fitting.pdfbasic import PDF1, Generic1D_pdf
from   ostap.core.ostap_types import integer_types 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.morphing_pdf' )
else                       : logger = getLogger ( __name__                     )
# =============================================================================
## @class Morphing1D_pdf
#  Wrapper for <code>RooMomentMorph</code> PDF
#  - 1D morphing/1D PDF
#  @see RooMomentMorph
#  @see Baak, M., Gadatsch, S., Harrington, R., & Verkerke, W. (2015).
#       "Interpolation between multi-dimensional histograms using
#       a new non-linear moment morphing method".
#       Nuclear Instruments & Methods in Physics Research.
#       Section A - Accelerators Spectrometers Detectors and Associated Equipment, 771, 39-48.
#  @see https://doi.org/10.1016/j.nima.2014.10.033
class Morphing1D_pdf (PDF1) :
    """ Wrapper for ROOT.RooMomentMorph PDF
    - 1D morphing/1D PDF
    - see ROOT.RooMomentMorph
    - see Baak, M., Gadatsch, S., Harrington, R., & Verkerke, W. (2015).
    'Interpolation between multi-dimensional histograms using
    a new non-linear moment morphing method'.
    Nuclear Instruments & Methods in Physics Research.
    Section A - Accelerators Spectrometers Detectors and Associated Equipment, 771, 39-48.
    - see https://doi.org/10.1016/j.nima.2014.10.033
    """
    def __init__ ( self                ,
                   name                , ## PDF name 
                   pdfs                , ## dictionary {mu1,mu2 -> pdf }
                   setting    = None   , ## morphing setting 
                   morph_var  = None   , ## morphing variable mu 
                   xvar       = None ) : ## observable (1D) 

        assert pdfs and 2 <= len ( pdfs ) , \
               "Invalid dictionary of morphing PDFs!"

        if setting is None : setting = ROOT.RooMomentMorph.Linear 
            
        assert isinstance ( setting , integer_types ) and 0 <= setting < 5,\
               'Invalid value for the setting %s/%s' %  ( setting , type ( setting ) )

        for k in pdfs :
            if xvar : break
            p = pdfs [ k ]
            if isinstance ( p , PDF ) :
                xvar = p.xvar
                break
        else :
            raise TypeError("Morphing_pdf: cannot identify xvar!")
                
        ## initialize the base class 
        PDF1.__init__ ( self , name , xvar )

        ## convert the dictionary of PDFs into  ordered list/tuple of pairs (mu,pdf)        
        self.__pdflist = []
        for k in sorted ( pdfs.keys()  ) :
            
            pdfk = pdfs [ k ]

            if   isinstance ( pdfk , PDF            ) and pdfk.xvar is self.xvar : pass 
            elif isinstance ( pdfk , ROOT.RooAbsPdf ) : 
                pdfk = Generic1D_pdf ( pdfk , xvar = self.xvar ) 
            else :
                pass 

            pair = k , pdfk 
            self.pdflist.append ( pair )
            
        ## convert to tuple 
        self.__pdflist = tuple ( self.__pdflist ) 

        ## save setting 
        self.__setting = setting
        
        ## min and maximal value of morhing parameter
        mu_min = self.pdflist[ 0][0]
        mu_max = self.pdflist[-1][0]

        ## create morphing variable 
        self.__mu = self.make_var (
            morph_var                ,
            "mu_%s"           % name ,
            "morphing mu(%s)" % name , morph_var , mu_min , mu_max )

        ## vector of morphing values 
        muvct  = ROOT.TVectorD ( len ( self.pdflist ) )
        pdflst = ROOT.RooArgList()
        for i , p in enumerate ( self.pdflist ) :
            muvct [ i ] = p [ 0 ]
            pdflst.add (  p [ 1 ].pdf ) 

        observables = ROOT.RooArgList (  self.xvar )
        self.aux_keep.append  ( observables )
        
        ## create the PDF  
        self.pdf = ROOT.RooMomentMorph (
            self.roo_name ( 'morph_' ) ,
            "Morphing %s" % self.name  , 
            self.mu               , ## morphing variable
            observables           , ## observables 
            pdflst                , ## ordered list of PDFs  
            muvct                 , ## values of morhing parameter 
            self.setting          ) ## morphing setting 

        #
        self.config = {
            'name'      : self.name             ,
            'setting'   : self.setting          ,
            'xvar'      : self.xvar             ,
            'morph_var' : self.mu               ,
            'pdfs'      : dict ( self.pdflist ) ,
            }
        
    @property
    def mu ( self ) :
        """'mu' : morphing variable"""
        return self.__mu
    @mu.setter
    def mu ( self , value ) :
        self.set_value ( self.__mu , value ) 

    @property
    def pdflist ( self ) :
        """'pdflist' : (sorted) tuple of (morphing parameter, pdf) pairs"""
        return self.__pdflist

    @property 
    def setting ( self ) :
        """'setting': morphing setting"""
        return self.__setting 


# =============================================================================
## @class Morphing2D_pdf
#  Wrapper for <code>RooMomentMorph</code> PDF
#  - 1D morphing/2D  PDF   
#  @see RooMomentMorphND
#  @see Baak, M., Gadatsch, S., Harrington, R., & Verkerke, W. (2015).
#       "Interpolation between multi-dimensional histograms using
#       a new non-linear moment morphing method".
#       Nuclear Instruments & Methods in Physics Research.
#       Section A - Accelerators Spectrometers Detectors and Associated Equipment, 771, 39-48.
#  @see https://doi.org/10.1016/j.nima.2014.10.033
class Morphing2D_pdf (PDF1) :
    """ Wrapper for ROOT.RooMomentMorphND PDF for N = 2 
    - 1D morphing/2D PDF
    - see ROOT.RooMomentMorphND
    - see Baak, M., Gadatsch, S., Harrington, R., & Verkerke, W. (2015).
    'Interpolation between multi-dimensional histograms using
    a new non-linear moment morphing method'.
    Nuclear Instruments & Methods in Physics Research.
    Section A - Accelerators Spectrometers Detectors and Associated Equipment, 771, 39-48.
    - see https://doi.org/10.1016/j.nima.2014.10.033
    """
    def __init__ ( self                ,
                   name                , ## PDF name 
                   pdfs                , ## dictionary {mu1,mu2 -> pdf }
                   setting    = None   , ## morphing setting 
                   morph_var1 = None   , ## morphing variable mu1 
                   morph_var2 = None   , ## morphing variable mu2 
                   xvar       = None ) : ## observable (1D) 

        assert pdfs and 2 <= len ( pdfs ) , \
               "Invalid dictionary of morphing PDFs!"
        
        if setting is None : setting = ROOT.RooMomentMorphND.Linear 

        assert isinstance ( setting , integer_types ) and 0 <= setting < 5,\
               'Invalid value for the setting %s/%s' %  ( setting , type ( setting ) )

        v1ps = set ()
        v2ps = set () 
        
        for k in pdfs :
            v1 , v2 = k
            v1ps.add ( v1 )
            v2ps.add ( v2 )
            p = pdfs [ k ]
            if not xvar and isinstance ( p , PDF1 ) :
                xvar = p.xvar

        assert xvar and isinstance ( xvar , ROOT.RooAbsReal ) , 'Cannot deduce xvar!'

        assert 2 <= len ( v1ps ) and 2 <= len ( v2ps ) , 'Invalid number of bins!'

        v1ps = list ( v1ps ) 
        v2ps = list ( v2ps )

        v1ps.sort ()
        v2ps.sort ()

        assert  len ( pdfs ) ==  len ( v1ps ) * len ( v2ps ) ,\
               'Invalid table/dict structure!'
        
        ## initialize the base class 
        PDF1.__init__ ( self , name , xvar )

        ## create morphing variables 
        self.__mu1 = self.make_var (
            morph_var1                ,
            "mu1_%s"           % name ,
            "morphing mu1(%s)" % name , morph_var1 , v1ps[0] , v1ps[-1] )
        
        ## create morphing variables 
        self.__mu2 = self.make_var (
            morph_var2                ,
            "mu2_%s"           % name ,
            "morphing mu2(%s)" % name , morph_var2 , v2ps[0] , v2ps[-1] )
        
        self.__pdfdict = {}
        for k in sorted ( pdfs.keys()  ) :

            v1 , v2 = k
            
            pdfk = pdfs [ k ] 
            if   isinstance ( pdfk , PDF            ) and pdfk.xvar is self.xvar : pass 
            elif isinstance ( pdfk , ROOT.RooAbsPdf ) : 
                pdfk = Generic1D_pdf ( pdfk , xvar = self.xvar ) 
            else :
                raise TypeError( "Invalid component type: %s/%s" % ( pdfk , type ( pdfk ) ) ) 

            pair = k , pdfk 
            self.pdfdict[ ( v1, v2 ) ] = pdfk
            
        ## save setting 
        self.__setting = setting

        
        ## fill morphing grid
        from ostap.fitting.variables import binning
        
        bins_v1     = binning ( v1ps , name = 'morph1' ) 
        bins_v2     = binning ( v2ps , name = 'morph2' ) 
        self.__grid = ROOT.RooMomentMorphND.Grid ( bins_v1 , bins_v2 ) 
        
        for k in self.pdfdict :

            p       = self.pdfdict [ k ]
            
            v1 , v2 = k

            assert v1 in v1ps , 'Morphing2D_pdf: Invalid v1 value %s' % v1 
            assert v2 in v2ps , 'Morphing2D_pdf: Invalid v2 value %s' % v2 
            
            ib1 = v1ps.index ( v1 ) 
            ib2 = v2ps.index ( v2 )
            
            ## ib1 = bins_v1.binNumber ( v1 )
            ## ib2 = bins_v2.binNumber ( v2 )

            self.__grid.addPdf ( p.pdf , ib1 , ib2 )

        morph_vars  = ROOT.RooArgList ( self.mu1 , self.mu2 )
        observables = ROOT.RooArgList ( self.xvar )

        self.aux_keep.append  ( morph_vars  )
        self.aux_keep.append  ( observables )
        
        ## create the PDF  
        self.pdf = ROOT.RooMomentMorphND (
            self.roo_name ( 'morph2_' )   ,
            "Morphing 2D %s" % self.name  , 
            morph_vars              , ## morphing variables 
            observables             , ## observables 
            self.grid               , ## morphing grid 
            self.setting            ) ## morphing setting 

        #
        self.config = {
            'name'       : self.name    ,
            'setting'    : self.setting ,
            'xvar'       : self.xvar    ,
            'morph_var1' : self.mu1     ,
            'morph_var2' : self.mu2     ,
            'pdfs'       : self.pdfdict ,
            }
        
    @property
    def mu1 ( self ) :
        """'mu1' : the first morphing variable"""
        return self.__mu1
    @mu1.setter
    def mu1 ( self , value ) :
        v = float ( value ) 
        mm = self.__mu1.minmax()
        if mm and not mm[0] <= v <= mm[1] :
            self.error ( "Morphing parameter %s is outside [%s,%s]" % ( v , mm[0] , mm[1] ) )
        self.__mu1.setVal ( v )

    @property
    def mu2 ( self ) :
        """'mu2' : the second morphing variable"""
        return self.__mu2
    @mu2.setter
    def mu2 ( self , value ) :
        v = float ( value ) 
        mm = self.__mu2.minmax()
        if mm and not mm[0] <= v <= mm[1] :
            self.error ( "Morphing parameter %s is outside [%s,%s]" % ( v , mm[0] , mm[1] ) )
        self.__mu2.setVal ( v )

    @property
    def grid    ( self ) :
        """'grid' : morphing grid"""
        return self.__grid
    
    @property
    def pdfdict ( self ) :
        """'pdfdict' : Dictionary { morphing parameters : pdf } """
        return self.__pdfdict

    @property 
    def setting ( self ) :
        """'setting': morphing setting"""
        return self.__setting 
                       
    
        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
