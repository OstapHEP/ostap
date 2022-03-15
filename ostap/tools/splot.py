#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==========================================================================================
## @file ostap/tools/splot.py
#  Helper utilities to get sWeights in a form of a function or histogram
#  (often needed in practice, e.g to add these values to TTree, 
#  avoiding the direct usage of ROOT.RooStat.SPlot)
# 
#  @see RooStat::SPlot
#  @see M.Pivk, F.R. Le Deberder,
#      "SPlot: A Statistical tool to unfold data distributions" 
#       Published in Nucl.Instrum.Meth. A555 (2005) 356
#  @see http://arxiv.org/abs/physics/0402083
#  @see https://doi.org/10.1016/j.nima.2005.08.106
#  @date   2019-05-14
#  @author Vanya  BELYAEV Ivan.Belyaev@itep.ru

# =============================================================================
""" Helper utilities to get sWeights in a form of a function or histogram
(often needed in practice, e.g to add these values to TTree,
avoiding the direct usage of ROOT.RooStat.SPlot)

- see RooStat::SPlot
- see M.Pivk, F.R. Le Deberder,
... ``SPlot: A Statistical tool to unfold data distributions''
...  Published in Nucl.Instrum.Meth. A555 (2005) 356
- see http://arxiv.org/abs/physics/0402083
- see https://doi.org/10.1016/j.nima.2005.08.106
"""
# =============================================================================
__author__  = 'Vanya BELYAEV  Ivan.Belyaev@itep.ru'
__date__    = "2019-05-14"
__version__ = '$Revision$'
__all__     = (
    'sPlot1D'  , ## 1D-splot
    )
# =============================================================================
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
# =============================================================================
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.tools.splot' )
else                       : logger = getLogger ( __name__            )
# =============================================================================
import ostap.fitting.roofitresult
from   ostap.fitting.basic        import PDF,  Generic1D_pdf 
from   ostap.fitting.variables    import FIXVAR
from   ostap.histos.histos        import Histo1DFun 
# =============================================================================
#  @class sPlot1D
#  Helper class to get <code>sWeigts</code> in a form of a historgams  or function objects.
#  It is often useful to avoid the direct usage of ROOT.RooStat.SPlot 
#  @see RooStat::SPlot
#  @see M.Pivk, F.R. Le Deberder,
#      "SPlot: A Statistical tool to unfold data distributions" 
#       Published in Nucl.Instrum.Meth. A555 (2005) 356
#  @see http://arxiv.org/abs/physics/0402083
#  @see https://doi.org/10.1016/j.nima.2005.08.106
#  @date   2019-05-14
#  @code
#  dataset = ...
#  pdf     = ...
#  r , f   = pdf.fitTo ( dataset , ... )
#  s = sPlot1D ( pdf , dataset )
#  weigths  = s. weights
#  hweigths = s.hweights
#  @endcode
#  @author Vanya  BELYAEV Ivan.Belyaev@itep.ru
class sPlot1D(object) :
    """ Helper class to get sWeigtts in form of historgams/function objects.
    It is often useful to avoid the direct usage of ROOT.RooStat.SPlot 
    - see ROOT.RooStat.SPlot
    - see M.Pivk, F.R. Le Deberder,
    ...``SPlot: A Statistical tool to unfold data distributions'' 
    ...    Published in Nucl.Instrum.Meth. A555 (2005) 356
    - see http://arxiv.org/abs/physics/0402083
    - see https://doi.org/10.1016/j.nima.2005.08.106
    
    >>> pdf     = ...
    >>> r , f   = pdf.fitTo ( dataset , ... )
    >>> s = sPlot1D ( pdf , dataset )
    >>> weigths  = s. weights
    >>> hweigths = s.hweights
    """
    def __init__ ( self             ,
                   pdf              ,
                   dataset   = None , 
                   fitresult = None ,  
                   fast      = True ,   ## fast histogram filling ? (bin centers)
                   nbins     = 100  ,   ## histogram bining 
                   access    = {}   ,   ## histogram access options 
                   fitopts   = {}   ) : ## PDF.fitTo options 
                
        assert dataset or fitresult, 'Either dataset or fitresult must be specified!'
        
        assert isinstance ( pdf     , PDF              ) and \
               isinstance ( pdf.pdf , ROOT.RooAddPdf   ) and \
               len ( pdf.alist1 ) ==  len ( pdf.alist2 )     , 'Invalid type of PDF!'

        cmps   = pdf.alist2 
        names  = [ c.name for c in cmps ] 

        ## if  datset is specified - perform the fit
        if  dataset :
            
            vars   = pdf.pdf.getParameters ( dataset )

            ## make a proper (re)fit fixing everything  but yields
            with FIXVAR ( [ v  for v in vars if not v in cmps ] ) :
                logger.info ('Refit with the fixed parameters') 
                fitresult , f = pdf.fitTo ( dataset , silent = True , draw = False , **fitopts )
                
        elif fitresult :
            
            pars   = fitresult.floatParsFinal()
            pnames = set( [ p.name for p in pars ] )
            
            if set ( names )  != pnames :
                raise RuntimeError("Rerun fit with with only %s floating " % names )

        ## temlate historgam 
        template = pdf.make_histo ( nbins )

        ## dictionary of components 
        hcomponents  = {}

        ## the list of PDFs 
        cpdfs        = [ Generic1D_pdf ( p , xvar = pdf.xvar ) for p in pdf.alist1 ]
        
        for p , n in zip  ( cpdfs , names ) :
            
            if fast  : hc = p.roo_histo ( histo = template , events = False )
            else     : hc = p.    histo ( histo = template , errors = False )
            
            hcomponents [ n ] = hc

        ## sum of all histograms 
        hsum = template.clone()
        hsum.Reset() ;
        if not hsum.GetSumw2() : hsum . Sumw2() 

        for k in hcomponents : hsum += hcomponents[k] * fitresult( k )[0].value() 

        hweights = {}
        l = len ( names )
        for i in range ( l ) :

            cmp = template.clone() ;
            cmp.Reset() ;
            if not cmp.GetSumw2() : cmp.Sumw2() 
            for j in range ( l ) :
                cmp += fitresult.cov ( names[i] , names[j] ) * hcomponents [ names[j] ] 
                
            cmp /= hsum
            hweights [ names [ i ] ] = cmp 
            
        del hsum
        del template 

        components = {}
        for k in hcomponents : components [k] = Histo1DFun ( hcomponents [k] , **access )
        
        weights    = {} 
        for k in hweights    : weights    [k] = Histo1DFun ( hweights    [k] , **access )  
        
        self.__hcomponents = hcomponents
        self.__components  =  components 
        self.__hweights    = hweights 
        self.__weights     =  weights         
                
    @property
    def components  ( self ) :
        """``components''  :  get fit components (as functions)"""
        return self.__components

    @property
    def hcomponents ( self ) :
        """``hcomponents'' : get fit components  (as histograms)"""
        return self.__hcomponents
        
    @property
    def weights     ( self ) :
        """``weights''     : get sWeights (as functions)"""
        return self.__weights
    
    @property
    def hweights    ( self ) :
        """``hweights''    : get sWeights (as histograms)"""
        return self.__hweights
    

## =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


    
# =============================================================================
# The END 
# =============================================================================
