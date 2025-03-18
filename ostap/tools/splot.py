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
    'COWs'     , ## the basic (unbinned)  COWs object
    'sPLOT'    , ## (unbinned) sPlot machinery 
    'hPlot1D'  , ## binned version of 1D-splot
    'hPlot2D'  , ## binned version of 2D-splot
    'hPlot3D'  , ## binned version of 3D-splot
    )
# =============================================================================
from   ostap.core.ostap_types     import string_types, integer_types  
from   ostap.core.core            import Ostap
from   ostap.fitting.pdfbasic     import ( APDF1 ,
                                           PDF1  , Generic1D_pdf ,
                                           PDF2  , Generic2D_pdf ,
                                           PDF3  , Generic3D_pdf )
from   ostap.fitting.variables    import FIXVAR
from   ostap.histos.histos        import Histo1DFun, Histo2DFun 
from   ostap.utils.basic          import typename
from   ostap.utils.progress_conf  import progress_conf 
from   ostap.utils.progress_bar   import progress_bar
import ostap.trees.trees 
import ostap.fitting.roofitresult
import ROOT, abc
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
# =============================================================================
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.tools.splot' )
else                       : logger = getLogger ( __name__            )
# =============================================================================
pdf_types = ROOT.RooAddPdf , ROOT.RooRealSumPdf
# =============================================================================
# @class COWs
# (over) simplified version of COWS
# @see H.Dembinski, M.Kenzie, C.Langenbruch and M.Schmelling, 
#      "Custom Orthogonal Weight functions (COWs) for Event Classification",
#      NIMA 1040 (2022) 167270, arXiv:2112.04574
# @see https://doi.org/10.48550/arXiv.2112.04574
# @see https://doi.org/10.1016/j.nima.2022.167270   
class COWs ( object ) :
    """ (over) simplified version of COWS
    - see H.Dembinski, M.Kenzie, C.Langenbruch and M.Schmelling, 
           "Custom Orthogonal Weight functions (COWs) for Event Classification",
           NIMA 1040 (2022) 167270, arXiv:2112.04574
    - see https://doi.org/10.48550/arXiv.2112.04574
    - see https://doi.org/10.1016/j.nima.2022.167270
    """
    def __init__ ( self            ,
                   pdf             , 
                   dataset         ,
                   progress = True ) :

        assert isinstance ( pdf     , APDF1     ) , 'Invalid type of PDF!'
        assert isinstance ( pdf.pdf , pdf_types ) , 'Invalid type of pdf.pdf: %s' % typename ( pdf.pdf )
        
        self.__pdf      = pdf
        self.__cows     = None
        self.__progress = True if progress else False

        ## if dataset provides, prepare COWs object
        if dataset :
            ## 
            obsset   = pdf.pdf.getObservables( dataset )
            assert obsset , "Empy set of observables!"
            ## 
            progress = progress_conf    ( self.progress ) 
            self.__cows = Ostap.Utils.COWs ( self.pdf.pdf  , ## model 
                                             dataset       , ## data 
                                             self.pdf.vars , ## normalization 
                                             progress      ) 
            
    # =========================================================================
    ## Add COWs weights to TTree
    #  @code
    #  cows = ...
    #  tree = ...
    #  cows.cows2tree ( tree , parallel = True , .. )'
    #  @endcode 
    def cows2tree ( self , tree , *vars , names = [] , parallel = False ) :
        """ Add COWs results to TTree
        >>> cows = ...
        >>> tree = ...
        >>> cows.cows2tree ( tree , parallel = True , suffix = '_sw'
        """
        assert tree and isinstance ( tree , ROOT.TTree   ) , "COWs: Ivalid tree!"
        for var in vars : 
            assert var and isinstance ( var , string_types ) and var in tree, \
                "COWS: Variable '%s' not in the tree" % var 

        assert self.cows , "COWs: `Ostap.Utils.COWs` object is None!"
        
        if not names : names = tuple (  'CMP%02d_cw' % i for i in range ( self.cows.size ) ) 
        assert len ( names ) == self.cows.size() and all ( isinstance ( n , string_types ) for n in names ) , \
            "Invalid `names` %s " % str ( names )
        
        if parallel :
            import ostap.parallel.parallel_add_branch

        kwargs  = { 'names' : names , 'progress' : self.progress , 'report' : True }
        
        mapping = {} 
        for a , b in zip ( vars , self.pdf.vars ) : mapping [ b.name ] = a
        kwargs [ 'mapping' ]  = mapping

        if parallel : result = tree.padd_new_branch ( self.cows , **kwargs ) 
        else        : result = tree. add_new_branch ( self.cows , **kwargs )
        
        return result 

    @property
    def progress ( self  ) :
        """`progress` : show progress?"""
        return self.__progress

    @property
    def pdf      ( self ) :
        """`pdf` : actual pdf for COWs/SPLOT"""
        return self.__pdf

    @property
    def cows     ( self ) :
        """`cows` : get he actual `Ostap.Utils.COWs` object"""
        return self.__cows

    def __getstate__ ( self ) :
        return { 'pdf'      : self.pdf      ,
                 'cows'     : self.cows     ,
                 'progress' : self.progress }
    
    def __setstate__ ( self, state ) :
        self.__pdf      = state.pop ( 'pdf'      )
        self.__cows     = state.pop ( 'cows'     )
        self.__progress = state.pop ( 'progress' )
        
# =============================================================================
# @class sPLOT
# Helper class to get <code>sWeigts</code> in a form of a histograms  or function objects.
# It is often useful to avoid the direct usage of ROOT.RooStat.SPlot 
# @see RooStat::SPlot
# @see M.Pivk, F.R. Le Deberder,
#      "SPlot: A Statistical tool to unfold data distributions" 
#       Published in Nucl.Instrum.Meth. A555 (2005) 356
# @see http://arxiv.org/abs/physics/0402083
# @see https://doi.org/10.1016/j.nima.2005.08.106
# @date   2019-05-14
# @code
# dataset = ...
# pdf     = ...
# r , f   = pdf.fitTo ( dataset , ... )
# s = sPlot1D ( pdf , dataset )
# weigths  = s. weights
# hweigths = s.hweights
# @endcode
# @author Vanya  BELYAEV Ivan.Belyaev@itep.ru
class sPLOT(COWs) :
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
    def __init__ ( self              ,
                   pdf               ,                   
                   dataset    = None , 
                   fitresult  = None ,  
                   fitopts    = {}   ,   ## PDF.fitTo options 
                   progress   = True ) : 
        
        assert dataset or fitresult, 'Either dataset or fitresult must be specified!'        
        assert isinstance ( pdf     , APDF1     ) , 'Invalid type of PDF!'
        assert isinstance ( pdf.pdf , pdf_types ) , 'Invalid type of pdf.pdf: %s' % typename ( pdf.pdf )
        ## there are some limits for sPlot:
        ## (1) pdf must be extended 
        assert pdf.pdf.canBeExtended()  , 'PDF cannot be extedned!'        
        assert len ( pdf.alist1 ) ==  len ( pdf.alist2 ) , \
            'PDF is not extended: %s != %ss' % ( len ( pdf.alist1 ) , len ( pdf.alist2 ) )

        if dataset and isinstance ( dataset , ROOT.TH1 ) :
            assert dataset.GetDimension() == len ( pdf.vars ) , 'Invalid histogram/dataset dimension!'
            params = pdf.pdf.getParameters ( pdf.vars )
            to_fix = tuple ( p for p in params if not p in pdf.alist2 )
            with FIXVAR ( to_fix ) :
                ns = ', '.join ( sorted ( v.name for v in to_fix ) ) 
                logger.info ( '(Re)fit the input histogram with fixed parameters: %s' % ns )                
                fitresult , _ = pdf.fitTo ( dataset , silent = True , draw = False , **fitopts )
                dataset = pdf.histo_data.dset 
 
        ## initialize the base
        super(sPLOT,self).__init__( pdf      = pdf                          ,
                                    dataset  = dataset if dataset else None , 
                                    progress = progress                     ) 

        # =====================================================================
        cmps   = self.pdf.alist2 
        names  = tuple ( sorted ( set ( c.name for c in cmps ) ) ) 
        # =====================================================================
        ## if dataset is specified - perform the fit
        if dataset : # ========================================================
            # =================================================================
            ## get all parameters 
            if isinstance ( dataset , ROOT.RooAbsData ) : vars = self.pdf.pdf.getParameters ( dataset )
            else : vars = tuple ( v for v in self.pdf.pdf.getParameters ( ROOT.nullptr ) if not v in self.pdf.vars )
            ## make a proper (re)fit fixing everything but the yields
            to_fix =  [ v for v in vars if not v in cmps ]
            with FIXVAR ( to_fix ) :
                ns = ', '.join ( sorted ( v.name for v in to_fix ) ) 
                logger.info ( '(Re)fit with fixed parameters: %s' % ns )
                fitresult , _ = self.pdf.fitTo ( dataset , silent = True , draw = False , **fitopts )
                title = '(Re)fit with fixed parameters'
                logger.info ( '%s:\n%s' % ( title , fitresult.table ( title = title , prefix = '# ' ) ) ) 
            # =================================================================
        elif fitresult : # ====================================================
            # =================================================================            
            pars   = fitresult.floatParsFinal()
            if pars and all ( ( p in self.pdf.alist2 ) for p in pars ) : pass
            else :
                ns = ', '.join ( n for n in names  )                 
                raise TypeError ( "Please re-run fit with with (at most) floating: %s" % ns )
            
        self.__fitresult = fitresult     
        self.__splot     = Ostap.Utils.SPLOT ( self.pdf.pdf   ,
                                               self.pdf.vars  ,
                                               self.fitresult ,
                                               self.pdf.vars  )     
    @property
    def fitresult ( self ) :
        """`fitresult` : fit results"""
        return self.__fitresult 

    @property
    def splot     ( self ) :
        """`splot` : get the actual `Ostap.Utils.SPLOT` object 
        """
        return self.__splot

    def __getstate__ ( self ) :
        state = super(sPLOT,self).__getstate__ ()
        state [ 'fitresult' ] = self.fitresult
        state [ 'splot' ]     = self.splot 
        return state

    def __setstate__ ( self , state ) :
        super(sPLOT,self).__setstate__ ( state )
        self.__fitresult = state.pop ( 'fitresult' )
        self.__splot     = state.pop ( 'splot'     )
        
    # =========================================================================
    ## Add sPLOT weights to TTree
    #  @code
    #  splot = ...
    #  tree  = ...
    #  splot.splot2tree ( tree , parallel = True , suffix = '_sw'
    #  @endcode 
    def splot2tree ( self , tree , *vars , prefix = '' , suffix = '_sw' , parallel = False ) :
        """ Add SPLOT results to TTree
        >>> splot = ...
        >>> tree  = ...
        >>> cows.splot2tree ( tree , parallel = True , suffix = '_sw'
        """
        assert tree and isinstance ( tree , ROOT.TTree   ) , "COWs: Ivalid tree!"
        for var in vars : 
            assert var and isinstance ( var , string_types ) and var in tree, "COWs: Variable '%s' not in the tree" % var 

        assert self.splot , "sPlot: `Ostap.Utils.SPLOT` object is None!"

        if parallel :
            import ostap.parallel.parallel_add_branch
            
        kwargs  = { 'prefix' : prefix , 'suffix' : suffix , 'progress' : self.progress , 'report' : True }
        
        mapping = {} 
        for a , b in zip ( vars , self.pdf.vars ) : mapping [ b.name ] =  a
        kwargs [ 'mapping' ]  = mapping
    
        if parallel : result = tree.padd_new_branch ( self.splot , **kwargs ) 
        else        : result = tree. add_new_branch ( self.splot , **kwargs )

        return result

# ==========================================================================
## @class hPlot
#  Binned versio of sPLOT
# Helper class to get <code>sWeigts</code> in a form of a histograms  or function objects.
# It is often useful to avoid the direct usage of ROOT.RooStat.SPlot 
# @see RooStat::SPlot
# @see M.Pivk, F.R. Le Deberder,
#      "SPlot: A Statistical tool to unfold data distributions" 
#       Published in Nucl.Instrum.Meth. A555 (2005) 356
# @see http://arxiv.org/abs/physics/0402083
# @see https://doi.org/10.1016/j.nima.2005.08.106
# @date   2019-05-14
# @code
# dataset = ...
# pdf     = ...
# r , f   = pdf.fitTo ( dataset , ... )
# s = hPlot1D ( pdf , dataset )
# weigths  = s. weights
# hweigths = s.hweights
# @endcode
# @author Vanya  BELYAEV Ivan.Belyaev@itep.ru
class hPlot(sPLOT) :
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
    >>> s = hPlot1D ( pdf , dataset )
    >>> weigths  = s. weights
    >>> hweigths = s.hweights
    """
    def __init__ ( self              ,
                   pdf               ,                   
                   template          , 
                   dataset    = None , 
                   fitresult  = None ,  
                   fast       = True ,   ## fast histogram filling ? (bin centers)
                   access     = {}   ,   ## histogram access options
                   fitopts    = {}   ,   ## PDF.fitTo options 
                   progress   = True ) : 
        
        ## initialize the base
        super(hPlot,self).__init__( pdf       = pdf       ,
                                    dataset   = dataset   , 
                                    fitresult = fitresult ,
                                    fitopts   = fitopts   , 
                                    progress  = progress  ) 

        assert isinstance ( template , ROOT.TH1 ) and template.GetDimension() == len ( pdf.vars ) , \
            "Invalid type of template histogram: %s" % typename ( template )

        self.__template    = template 
        self.__fast        = True if fast else False 
        self.__access      = access

        self.__hcomponents = {} 
        self.__hweights    = {} 
        self.__components  = {} 
        self.__weights     = {} 

    @property
    def access ( self ) :
        """'access' : parameters for histogram->function transition"""
        return self.__access
    
    @property
    def fast    ( self  ) :
        """`fast` : use bin centers of intergal over bins? """
        return self.__fast

    @property
    def hcomponents ( self ) :
        """'hcomponents' : get fit components  (as 1D-histograms)"""
        if not self.__hcomponents : self.calculate () 
        return self.__hcomponents
        
    @property
    def hweights    ( self ) :
        """'hweights'    : get sWeights (as 1D-histograms)"""
        if not self.__hweights : self.calculate () 
        return self.__hweights

    @property
    def components  ( self ) :
        """'components'  :  get fit components (as nD-functions)"""
        if not self.__components : self.calculate() 
        return self.__components 

    @property
    def weights ( self ) :
        """'weights'     : get sWeights (as nD-functions)"""
        if not self.__weights : self.calculate() 
        return self.__weights

    def __getstate__ ( self ) :
        state = super(hPlot,self).__getstate__ ()
        state [ 'hcomponents' ] = self.hcomponents 
        state [ 'hweights'    ] = self.hweights 
        return state

    def __setstate__ ( self , state ) :
        super(hPlot,self).__setstate__ ( state )
        self.__hcomponents = state.pop ( 'hcomponents' )
        self.__hweights    = state.pop ( 'hweights'    )
        self.__components  = {} 
        self.__weights     = {} 
        
    ## calculate sPlot histograms 
    def calculate ( self , progress = True )  :
        """ Calculate sPlot histograms 
        """
        
        ## dictionary of components 
        hcomponents  = {}
        
        ## get the list of PDFs 
        cpdfs        = self.pdf.alist_cmp 

        if not self.fast : logger.warning  ( "sPlot(fast=False) is specified: it could be biased!" )
        
        ## total (extended) sum of all components 
        total = None

        ## 1) get the fit components in a form of histograms & calcualte the total (extended) sum
        nmax = min ( len ( cpdfs ) , len ( self.pdf.alist2 ) ) 
        for p , v in progress_bar ( zip ( cpdfs , self.pdf.alist2 ) ,
                                    max_value   = nmax                ,
                                    description = '#Components:' , silent = not self.progress ) : 
            
            if self.fast : hc = p.roo_histo ( histo = self.__template , events = False )
            else         : hc = p.    histo ( histo = self.__template , events = False , errors = False , integral = False )

            key = v.name
            
            hcomponents [ key ] = hc

            ## calculate the total (extended) sum 
            factor = float ( self.fitresult ( key ) [ 0 ] )
            if not total : total  = hc * factor 
            else         : total += hc * factor 

        # =======================================================================
        ## calculate the weights
        # ========================================================================
        cmps   = self.pdf.alist2 
        names  = tuple ( sorted ( set ( c.name for c in cmps ) ) ) 

        hweights = {}
        for i , iname  in enumerate ( names ) :
            
            cmp = None 
            for j, jname in enumerate ( names ) : 
                cov_ij = self.fitresult.cov ( iname , jname ) 
                if not cmp : cmp  = cov_ij * hcomponents [ jname ] 
                else       : cmp += cov_ij * hcomponents [ jname ] 
                
            cmp /= total 
            hweights [ names [ i ] ] = cmp
            
        del total 

        self.__hcomponents = hcomponents
        self.__hweights    = hweights 

        self.__components = {}
        self.__weights    = {}        
        for k in self.__hcomponents : self.__components [ k ] = self.HFUN ( self.__hcomponents [ k ] , **self.access )
        for k in self.__hweights    : self.__weights    [ k ] = self.HFUN ( self.__hweights    [ k ] , **self.access )
        
    # =========================================================================
    ## Add sPlot results to TTree
    #  @code
    #  splot = ...
    #  tree  = ...
    #  splot.add_to_tree ( tree , parallel = True , suffix = '_sw'
    #  @endcode 
    def _add_to_tree ( self , tree , *vars , prefix = '' , suffix = '_hw' , parallel = False ) :
        """ Add sPlot results to TTree
        >>> splot = ...
        >>> tree  = ...
        >>> splot.add_to_tree ( tree , parallel = True , suffix = '_hw'
        """
        assert tree and isinstance ( tree , ROOT.TTree   ) , "sPlot: Ivalid tree!"
        for var in vars : 
            assert var and isinstance ( var , string_types ) and var in tree, \
                "sPlot: Variable '%s' not in the tree" % var
            
        if parallel :
            import ostap.parallel.parallel_add_branch
            
        suffix = suffix.replace ( ' ' , '' ).replace ( '-' , '_' ).strip() 
        prefix = prefix.replace ( ' ' , '' ).replace ( '-' , '_' ).strip() 

        fmap   = {}
        for ww in self.hweights  :            
            hw     = self.hweights [ ww ]
            
            newvar = '%s%s%s' % ( prefix , ww , suffix )
            fw     = self.make_thfun ( hw , *vars )        
            
            assert not newvar in tree , "sPlot: Variable '%s' is already in the Tree!" % newvar 
            fmap [ newvar ] = fw
            
            vvs = sorted ( fmap.keys() )
            bbs = ','.join ( vvs )
            logger.info ( "Adding sPlot results %s to TTree" % vvs )
            
        kwargs = { 'progress' : self.progress , 'report' : True }
        if parallel : result = tree.padd_new_branch ( fmap , **kwargs ) 
        else        : result = tree. add_new_branch ( fmap , **kwargs )
        
        return result 

    # =========================================================================
    ## Create proper TH1/TH2-function 
    @abc.abstractmethod
    def make_thfun ( self , histo , *vars ) :
        """ Create propeer TH1/TH2-function
        """
        pass 

# =============================================================================
## @class hPlot1D
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
class hPlot1D(hPlot) :
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
                   nbins     = 200  ,   ## histogram binning 
                   access    = {}   ,   ## histogram access options 
                   fitopts   = {}   ,   ## PDF.fitTo options 
                   progress  = True ) :

        ## 
        assert isinstance ( pdf   , PDF1          ) ,               'Invalid type of PDF!'
        assert isinstance ( nbins , integer_types ) and 1 < nbins , "Invalid nhins!"
        
        if dataset and isinstance ( dataset , ROOT.TH1 ) :
            assert dataset.GetDimension() == 1 , 'Invalid dataset/histogram dimension!'
        
        ## template histogram 
        template = pdf.make_histo ( nbins )
        
        ## initialize the base class 
        hPlot.__init__ ( self ,
                         pdf       = pdf       ,
                         template  = template  ,
                         dataset   = dataset   ,
                         fitresult = fitresult ,
                         fast      = fast      ,
                         access    = access    ,
                         fitopts   = fitopts   ,
                         progress  = progress  )
        
        ## we do not need it anymore 
        del template
        
    ## H-function type 
    def HFUN  ( self , *args , **kwargs ) : return Histo1DFun ( *args , **kwargs ) 
        
    # =========================================================================
    ## Add sPlot results to TTree
    #  @code
    #  splot = ...
    #  tree  = ...
    #  splot.add_to_tree ( tree , xvar , parallel = True , suffix = '_sw'
    #  @endcode 
    def add_to_tree ( self , tree , xvar , * , suffix = '_hw' , prefix = '' , parallel = False ) :
        """ Add sPlot results to TTree
        >>> splot = ...
        >>> tree  = ...
        >>> splot.add_to_tree ( tree , xvar , parallel = True , suffix = '_sw'
        """
        return self._add_to_tree ( tree              ,
                                   xvar              ,
                                   suffix   = suffix ,
                                   prefix   = prefix ,
                                   parallel = parallel )
    
    # =========================================================================
    ## make proper TH1 function 
    def make_thfun ( self , histo , xvar ) :
        """ Create proper TH1-funcion
        """
        return Ostap.Functions.FuncTH1 ( histo , xvar )

# =============================================================================
## @class hPlot2D
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
#  s = sPlot2D ( pdf , dataset )
#  weigths  = s. weights
#  hweigths = s.hweights
#  @endcode
#  @author Vanya  BELYAEV Ivan.Belyaev@itep.ru
class hPlot2D(hPlot) :
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
    >>> s = sPlot2D ( pdf , dataset )
    >>> weigths  = s. weights
    >>> hweigths = s.hweights
    """
    def __init__ ( self                      ,
                   pdf                       ,
                   dataset   = None          , 
                   fitresult = None          ,  
                   fast      = True          ,   ## fast histogram filling ? (bin centers)                   
                   xbins     = 100           ,   ## histogram binning 
                   ybins     = 100           ,   ## histogram binning 
                   access    = {}            ,   ## histogram access options 
                   fitopts   = {}            ,   ## PDF.fitTo options 
                   progress  = True          ) :

        ## 
        assert isinstance ( pdf   , PDF2          )               , 'Invalid type of PDF!'
        assert isinstance ( xbins , integer_types ) and 1 < xbins , "Invalid xhins!"
        assert isinstance ( ybins , integer_types ) and 1 < ybins , "Invalid xhins!"
        
        if dataset and isinstance ( dataset , ROOT.TH1 ) :
            assert isinstance ( dataset , ROOT.TH2 ) and dataset.GetDimension() == 2 , \
                'Invalid dataset/histogram dimension!'
            
        ## template histogram 
        template = pdf.make_histo ( xbins , ybins )

        ## initialize the base class 
        hPlot.__init__ ( self ,
                         pdf       = pdf       ,
                         template  = template  ,
                         dataset   = dataset   ,
                         fitresult = fitresult ,
                         fast      = fast      ,
                         access    = access    ,
                         fitopts   = fitopts   ,
                         progress  = progress  )

        ## we do not need it anymore 
        del template
        
    ## H-function type 
    def HFUN  ( self , *args , **kwargs ) : return Histo2DFun ( *args , **kwargs ) 

    # =========================================================================
    ## Add sPlot results to TTree
    #  @code
    #  splot = ...
    #  tree  = ...
    #  splot.add_to_tree ( tree , xvar , yvar , parallel = True , suffix = '_sw'
    #  @endcode 
    def add_to_tree ( self , tree , xvar , yvar , * , suffix = '_hw' , prefix = '' , parallel = False ) :
        """ Add sPlot results to TTree
        >>> splot = ...
        >>> tree  = ...
        >>> splot.add_to_tree ( tree , xvar , yvar , parallel = True , suffix = '_sw'
        """
        return self._add_to_tree ( tree                ,
                                   xvar , yvar         ,
                                   suffix   = suffix   ,
                                   prefix   = prefix   ,
                                   parallel = parallel )
    
    # =========================================================================
    ## make proper TH2-function 
    def make_thfun ( self , histo , xvar , yvar ) :
        """ Create propeer TH2-funcion
        """
        return Ostap.Functions.FuncTH2 ( histo , xvar , yvar )
    
# =============================================================================
## @class hPlot3D
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
#  s = sPlot3D ( pdf , dataset )
#  weigths  = s. weights
#  hweigths = s.hweights
#  @endcode
#  @author Vanya  BELYAEV Ivan.Belyaev@itep.ru
class hPlot3D(hPlot) :
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
    >>> s = sPlot3D ( pdf , dataset )
    >>> weigths  = s. weights
    >>> hweigths = s.hweights
    """
    def __init__ ( self                      ,
                   pdf                       ,
                   dataset   = None          , 
                   fitresult = None          ,  
                   fast      = True          , ## fast histogram filling ? (bin centers)                   
                   xbins     = 100           , ## histogram binning 
                   ybins     = 100           , ## histogram binning 
                   zbins     = 100           , ## histogram binning 
                   access    = {}            , ## histogram access options 
                   fitopts   = {}            , ## PDF.fitTo options 
                   progress  = True          ) :

        ## 
        assert isinstance ( pdf   , PDF3          )               , 'Invalid type of PDF!'
        assert isinstance ( xbins , integer_types ) and 1 < xbins , "Invalid xhins!"
        assert isinstance ( ybins , integer_types ) and 1 < ybins , "Invalid yhins!"
        assert isinstance ( zbins , integer_types ) and 1 < zbins , "Invalid zhins!"
        
        if dataset and isinstance ( dataset , ROOT.TH1 ) :
            assert isinstance ( dataset , ROOT.TH3 ) and dataset.GetDimension() == 3 , \
                'Invalid dataset/histogram dimension!'
            
        ## template histogram 
        template = pdf.make_histo ( xbins , ybins , zbins )

        ## H-function type 
        self.HFUN = Histo3DFun
        
        ## initialize the base class 
        hPlot.__init__ ( self ,
                         pdf       = pdf       ,
                         template  = template  ,
                         dataset   = dataset   ,
                         fitresult = fitresult ,
                         fast      = fast      ,
                         access    = access    ,
                         fitopts   = fitopts   ,
                         progress  = progress  )

        ## we do not need it anymore 
        del template
        
    ## H-function type 
    def HFUN  ( self , *args , **kwargs ) : return Histo3DFun ( *args , **kwargs ) 

    # =========================================================================
    ## Add sPlot results to TTree
    #  @code
    #  splot = ...
    #  tree  = ...
    #  splot.add_to_tree ( tree , xvar , yvar , xvar , parallel = True , suffix = '_sw'
    #  @endcode 
    def add_to_tree ( self , tree , xvar , yvar , zvar , * , suffix = '_hw' , prefix = '' , parallel = False ) :
        """ Add sPlot results to TTree
        >>> splot = ...
        >>> tree  = ...
        >>> splot.add_to_tree ( tree , xvar , yvar , zvar , parallel = True , suffix = '_sw'
        """
        return self._add_to_tree ( tree                ,
                                   xvar , yvar , zvar  ,
                                   suffix   = suffix   ,
                                   prefix   = prefix   ,
                                   parallel = parallel )
    
    # =========================================================================
    ## make proper TH3 function 
    def make_thfun ( self , histo , xvar , yvar , zvar ) :
        """ Create proper TH3-funcion
        """
        return Ostap.Functions.FuncTH3 ( histo , xvar , yvar , zvar )
    
# =============================================================================
## reconstruct COWs object
def _cows_factory_ ( *args ) :
    """ Reconstruct SPLOT object
    """
    result = Ostap.Utils.COWs ( *args )
    result._args = args
    return result
# =============================================================================
## reduce COWS object 
def _cows_reduce_  ( cows   ) :
    """ Reduce COWS object
    """
    return _cows_factory_ , ( cows.pdf           () ,
                              cows.observables   () ,                            
                              cows.normalization () , 
                              cows. A            () )

Ostap.Utils.COWs.__reduce__ = _cows_reduce_

# =============================================================================
## reconstruct SPLOT object
def _splot_factory_ ( *args ) :
    """ Reconstruct SPLOT object
    """
    result = Ostap.Utils.SPLOT ( *args )
    result._args = args
    return result
# =============================================================================
## reduce SPlit3Tree object 
def _splot_reduce_  ( sp4t   ) :
    """ Reduce SPLOT object
    """
    return _splot_factory_ , ( sp4t.pdf           () ,
                              sp4t.observables   () ,
                              sp4t.fitresult     () ,                            
                              sp4t.normalization () )

Ostap.Utils.SPLOT.__reduce__ = _splot_reduce_

## =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


    
# =============================================================================
##                                                                      The END 
# =============================================================================
