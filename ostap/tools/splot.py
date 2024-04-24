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
... `SPlot: A Statistical tool to unfold data distributions'
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
from   ostap.core.core            import Ostap 
from   ostap.core.ostap_types     import string_types 
from   ostap.fitting.pdfbasic     import PDF1, Generic1D_pdf 
from   ostap.fitting.variables    import FIXVAR
from   ostap.histos.histos        import Histo1DFun 
import ostap.fitting.roofitresult
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
# =============================================================================
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.tools.splot' )
else                       : logger = getLogger ( __name__            )
# =============================================================================
#  @class sPlot1D
#  Helper class to get <code>sWeigts</code> in a form of a histograms  or function objects.
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
    """ Helper class to get sWeigtts in form of histograms/function objects.
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
        
        assert isinstance ( pdf     , PDF1             ) and \
               isinstance ( pdf.pdf , ROOT.RooAddPdf   ) and \
               len ( pdf.alist1 ) ==  len ( pdf.alist2 )     , 'Invalid type of PDF!'

        cmps   = pdf.alist2 
        names  = sorted ( set ( c.name for c in cmps ) ) 

        ## if  datset is specified - perform the fit
        if  dataset :
            
            vars   = pdf.pdf.getParameters ( dataset )

            ## make a proper (re)fit fixing everything  but yields
            to_fix =  [ v  for v in vars if not v in cmps ]
            with FIXVAR ( to_fix ) :
                ns = ','.join ( sorted ( v.name for v in to_fix ) ) 
                logger.info ( 'Refit with the fixed parameters: %s' % ns )
                
                fitresult , f = pdf.fitTo ( dataset , silent = True , draw = False , **fitopts )
                ## fitresult , f = pdf.fitTo ( dataset , silent = True , draw = True , **fitopts )
                                
        elif fitresult :
            
            pars   = fitresult.floatParsFinal()
            pnames = sorted ( set ( p.name for p in pars ) ) 
            
            if names != pnames :                
                ns = ','.join ( n for n in names  ) 
                raise RuntimeError ( "Rerun fit with with only floating %s" % ns )

        ## template histogram 
        template = pdf.make_histo ( nbins )

        ## dictionary of components 
        hcomponents  = {}

        ## the list of PDFs 
        cpdfs        = [ Generic1D_pdf ( p , xvar = pdf.xvar ) for p in pdf.alist1 ]
        
        for p , v in zip  ( cpdfs , pdf.alist2  ) :
            
            if fast  : hc = p.roo_histo ( histo = template , events = False )
            else     : hc = p.    histo ( histo = template , events = False , errors = False )
            
            hcomponents [ v.name ] = hc

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
        self.__components  = components 
        self.__hweights    = hweights 
        self.__weights     =  weights         
        self.__access      =  access


    # =========================================================================
    ## Add sPlot results to TTree
    #  @code
    #  splot = ...
    #  tree  = ...
    #  splot.add_to_tree ( tree , parallel = True , suffix = '_sw'
    #  @endcode 
    def add_to_tree ( self , tree , var , suffix = '_sw' , prefix = '' , parallel = False ) :
        """ Add sPlot results to TTree
        >>> splot = ...
        >>> tree  = ...
        >>> splot.add_to_tree ( tree , parallel = True , suffix = '_sw'
        """
        import ostap.trees.trees 
        assert tree and isinstance ( tree , ROOT.TTree   ) , "sPlot: Ivalid tree!"
        assert var  and isinstance ( var  , string_types ) and var  in tree, "sPlot: Variable '%s' not in the tree" % var 
        
        if parallel :
            import ostap.parallel.parallel_add_branch
            
        suffix = suffix.strip   ()
        suffix = suffix.replace ( ' ' , '_' )
        prefix = prefix.strip   ()
        prefix = prefix.replace ( ' ' , '_' )
        
        fmap   = {}
        for ww in self.hweights  :            
            hw     = self.hweights [ ww ]
            fw     = Ostap.Functions.FuncTH1 ( hw , var )
            newvar = '%s%s%s' % ( prefix , ww , suffix )
            assert not newvar in tree , "sPlot: Variable '%s' is already in the Tree!" % newvar 
            fmap [ newvar ] = fw
            
        vvs = sorted ( fmap.keys() )
        bbs = ','.join ( vvs ) 
        logger.info ( "Adding sPlot results %s to TTree" % vvs ) 
        
        if parallel : result = tree.padd_new_branch ( None, fmap ) 
        else        : result = tree. add_new_branch ( None, fmap ) 

        return result 
            
    @property
    def components  ( self ) :
        """'components'  :  get fit components (as 1D-functions)"""
        return self.__components 

    @property
    def hcomponents ( self ) :
        """'hcomponents' : get fit components  (as 1D-histograms)"""
        return self.__hcomponents
        
    @property
    def weights     ( self ) :
        """'weights'     : get sWeights (as 1D-functions)"""
        return self.__weights
    
    @property
    def hweights    ( self ) :
        """'hweights'    : get sWeights (as 1D-histograms)"""
        return self.__hweights

    @property
    def access ( self ) :
        """'access' : parameters for histogram->function transition"""
        return self.__access

    ## serialisation/pickling 
    def __getstate__  ( self ) :
        """Serializarion/pickling
        """
        return { 'hcomponents' : self.hcomponents ,
                 'hweights'    : self.hweights    ,
                 'access'      : self.access      } 
                 
    ## deserialisation/unpickling
    def __setstate__  ( self , state ) :
        """Deserialisation/unpickling
        """

        self.__hcomponents = {}
        self.__hweights    = {}
        self.__access      = {}

        self.__hcomponents.update ( state.get ( 'hcomponents' , {} ) )
        self.__hweights   .update ( state.get ( 'hweights'    , {} ) )
        self.__access     .update ( state.get ( 'access'      , {} ) )
        
        self.__components = {}
        for k in self.__hcomponents :
            self.__components [k] = Histo1DFun ( self.__hcomponents [k] , **self.__access )
        
        self.__weights    = {} 
        for k in self.__hweights    :
            self.__weights    [k] = Histo1DFun ( self.__hweights    [k] , **self.__access )  

## =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


    
# =============================================================================
##                                                                      The END 
# =============================================================================
