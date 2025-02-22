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
    'sPlot2D'  , ## 2D-splot
    'sPlot3D'  , ## 3D-splot
    )
# =============================================================================
from   ostap.core.meta_info       import root_info
from   ostap.core.ostap_types     import string_types, integer_types  
from   ostap.core.core            import Ostap
from   ostap.fitting.pdfbasic     import ( APDF1 ,
                                           PDF1  , Generic1D_pdf ,
                                           PDF2  , Generic2D_pdf )
from   ostap.fitting.variables    import FIXVAR
from   ostap.histos.histos        import Histo1DFun, Histo2DFun 
from   ostap.utils.basic          import typename
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
#  @class sPlot
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
class sPlot(object) :
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
                   template          , 
                   dataset    = None , 
                   fitresult  = None ,  
                   fast       = True ,   ## fast histogram filling ? (bin centers)
                   access     = {}   ,   ## histogram access options
                   fitopts    = {}   ,   ## PDF.fitTo options 
                   progress   = True ) : 
        
        assert dataset or fitresult, 'Either dataset or fitresult must be specified!'

        assert isinstance ( pdf     , APDF1            ) , 'Invalid type of PDF!'
        assert isinstance ( pdf.pdf , pdf_types        ) , 'Invalid type of pdf.pdf: %s' % typename ( pdf.pdf )
        assert pdf.pdf.canBeExtended()                   , 'PDF cannot be extedned!'        
        assert len ( pdf.alist1 ) ==  len ( pdf.alist2 ) , 'PDF is not extedned: %s != %ss' %  ( len ( pdf.alist1 ) ,
                                                                                                 len ( pdf.alist2 ) )
        if dataset and isinstance ( dataset , ROOT.TH1 ) :
            assert dataset.GetDimension() == len ( pdf.vars ) , 'Invalid histogram/dataset dimension!'

        assert isinstance ( template , ROOT.TH1 ) and \
            template.GetDimension()       == len ( pdf.vars ) , 'Invalid template histogram typee!'
        
        cmps   = pdf.alist2 
        names  = tuple ( sorted ( set ( c.name for c in cmps ) ) ) 

        # =====================================================================
        ## if  dataset is specified - perform the fit
        if  dataset : # =======================================================
            # =================================================================
            ## get all parameters 
            if isinstance ( dataset , ROOT.RooAbsData ) : vars = pdf.pdf.getParameters ( dataset )
            else : vars = tuple ( v for v in pdf.pdf.getParameters ( ROOT.nullptr ) if not v in pdf.vars )
            ## make a proper (re)fit fixing everything  but yields
            to_fix =  [ v for v in vars if not v in cmps ]
            with FIXVAR ( to_fix ) :
                ns = ', '.join ( sorted ( v.name for v in to_fix ) ) 
                logger.info ( '(Re)fit with fixed parameters: %s' % ns )
                fitresult , _ = pdf.fitTo ( dataset , silent = True , draw = False , **fitopts )
                title = '(Re)fit with fixed parameters'
                logger.info ( '%s:\n%s' % ( title , fitresult.table ( title = title , prefix = '# ' ) ) ) 
            # =================================================================
        elif fitresult : # ====================================================
            # =================================================================            
            pars   = fitresult.floatParsFinal()
            if pars and all ( ( p in pdf.alist2 ) for p in pars ) : pass
            else :
                ns = ', '.join ( n for n in names  )                 
                raise TypeError ( "Please re-run fit with with (at most) floating: %s" % ns )

        self.__pdf         = pdf
        self.__fitresult   = fitresult 
        self.__template    = template 
        self.__fast        = fast
        self.__access      = access
        self.__progress    = True if progress else False

        self.__hcomponents = {} 
        self.__hweights    = {} 
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
                                    max_value = nmax                ,
                                    description = '#Components:' , silent = not self.progress ) : 
            
            if self.fast : hc = p.roo_histo ( histo = self.__template , events = False )
            else         : hc = p.    histo ( histo = self.__template , events = False , errors = False , integral = False )

            key = v.name
            
            hcomponents [ key ] = hc

            ## calcuate the total (extendd) sum 
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
                cov_ij = self.__fitresult.cov ( iname , jname ) 
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

    @property
    def pdf         ( self ) :
        """`pdf` : actual fit pdf"""
        return self.__pdf

    @property
    def fitresult  ( self ) :
        """`fitresult` : fit results"""
        return self.__fitresult 

    @property
    def access ( self ) :
        """'access' : parameters for histogram->function transition"""
        return self.__access
    
    @property
    def fast    ( self  ) :
        """`fast` : use bin centers of intergal over bins? """
        return self.__fast

    @property
    def progress ( self  ) :
        """`progress` : show progress?"""
        return self.__progress

    @property
    def hcomponents ( self ) :
        """'hcomponents' : get fit components  (as 1D-histograms)"""
        if not self.__hcomponent : self.calculate () 
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

    ## serialisation/pickling 
    def __getstate__  ( self ) :
        """ Serializarion/pickling
        """
        state = { 'pdf'       : self.pdf        ,
                  'fitresult' : self.fitresult  ,
                  'template'  : self.__template ,
                  'fast'      : self.fast       , 
                  'access'    : self.access     , 
                  'progress'  : self.progress   }
        if self.__hcomponents :  state ['hcomponents'] = self.__hcomponents
        if self.__hweights    :  state ['hweights']    = self.__hweights         
        return state 
    # =========================================================================
    ## deserialisation/unpickling
    def __setstate__  ( self , state ) :
        """ Deserialisation/unpickling
        """
        self.__pdf         = state.get ( 'pdf'         , None )
        self.__fitresult   = state.get ( 'fitresult'   , None )
        self.__template    = state.get ( 'template'    , None )
        self.__access      = state.get ( 'access'      , {}   )
        self.__fast        = state.get ( 'fast'        , True )
        self.__progress    = state.get ( 'progress'    , True )
        self.__hcomponents = state.get ( 'hcomponents' , {}   )
        self.__hweights    = state.get ( 'hweights'    , {}   )
        
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
    def _add_to_tree ( self , tree , *vars , prefix = '' , suffix = '_sw' , unbinned = True , parallel = False ) :
        """ Add sPlot results to TTree
        >>> splot = ...
        >>> tree  = ...
        >>> splot.add_to_tree ( tree , parallel = True , suffix = '_sw'
        """
        assert tree and isinstance ( tree , ROOT.TTree   ) , "sPlot: Ivalid tree!"
        for var in vars : 
            assert var  and isinstance ( var  , string_types ) and var in tree, "sPlot: Variable '%s' not in the tree" % var 

        if parallel :
            import ostap.parallel.parallel_add_branch
            
        suffix = suffix.replace ( ' ' , '_' ).replace ( '-' , '_' ).strip() 
        prefix = prefix.replace ( ' ' , '_' ).replace ( '-' , '_' ).strip() 

        if root_info < ( 6 , 24 , 6 ) and not unbinned :
            unbinned = True
            logger.warning ( 'Switch to *UNBINNED* processing for ROOT version %s.%s/%s' % ( root_info.major , 
                                                                                             root_info.minor , 
                                                                                             root_info.patch ) ) 

        ## process!! 
        if unbinned :
            
            ## add to tree using Splti4tree machinery :
            splot = Ostap.MoreRooFit.SPlot4Tree (
                self.pdf.pdf   ,
                self.pdf.vars  ,
                self.fitresult , 
                self.pdf.vars  ,
            )

            kwargs = { 'prefix' : prefix , 'suffix' : suffix , 'progress' : self.progress , 'report' : True }
            if parallel : result = tree.padd_new_branch ( splot , **kwargs ) 
            else        : result = tree. add_new_branch ( splot , **kwargs )

        else : 
        
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
## @class sPlot1D
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
class sPlot1D(sPlot) :
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
        sPlot.__init__ ( self ,
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
    def add_to_tree ( self , tree , xvar , * , suffix = '_sw' , prefix = '' , unbinned = True , parallel = False ) :
        """ Add sPlot results to TTree
        >>> splot = ...
        >>> tree  = ...
        >>> splot.add_to_tree ( tree , xvar , parallel = True , suffix = '_sw'
        """
        return self._add_to_tree ( tree , xvar , suffix = suffix , prefix = prefix , unbinned = unbinned , parallel = parallel )
    
    # =========================================================================
    ## make proper TH1 function 
    def make_thfun ( self , histo , xvar ) :
        """ Create proper TH1-funcion
        """
        return Ostap.Functions.FuncTH1 ( histo , xvar )

# =============================================================================
## @class sPlot2D
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
class sPlot2D(sPlot) :
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
            assert isinstance ( dataset , ROOT.TH2 ) and dataset.GetDimension() == 2 , 'Invalid dataset/histogram dimension!'
            
        ## template histogram 
        template = pdf.make_histo ( xbins , ybins )

        ## initialize the base class 
        sPlot.__init__ ( self ,
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
    def add_to_tree ( self , tree , xvar , yvar , * , suffix = '_sw' , prefix = ''  , unbinned = True , parallel = False ) :
        """ Add sPlot results to TTree
        >>> splot = ...
        >>> tree  = ...
        >>> splot.add_to_tree ( tree , xvar , yvar , parallel = True , suffix = '_sw'
        """
        return self._add_to_tree ( tree , xvar, yvar , suffix = suffix , prefix = prefix , unbinned = unbinned , parallel = parallel )
    
    # =========================================================================
    ## make proper TH2-function 
    def make_thfun ( self , histo , xvar , yvar ) :
        """ Create propeer TH2-funcion
        """
        return Ostap.Functions.FuncTH2 ( histo , xvar , yvar )
    
# =============================================================================
## @class sPlot2D
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
class sPlot3D(sPlot) :
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
            assert isinstance ( dataset , ROOT.TH3 ) and dataset.GetDimension() == 3 , 'Invalid dataset/histogram dimension!'
            
        ## template histogram 
        template = pdf.make_histo ( xbins , ybins , zbins )

        ## H-function type 
        self.HFUN = Histo3DFun
        
        ## initialize the base class 
        sPlot.__init__ ( self ,
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
    def add_to_tree ( self , tree , xvar , yvar , zvar , * , suffix = '_sw' , prefix = ''  , unbinned = True , parallel = False ) :
        """ Add sPlot results to TTree
        >>> splot = ...
        >>> tree  = ...
        >>> splot.add_to_tree ( tree , xvar , yvar , zvar , parallel = True , suffix = '_sw'
        """
        return self._add_to_tree ( tree , xvar, yvar , zvar , suffix = suffix , prefix = prefix  , unbinned = unbinned , parallel = parallel )
    
    # =========================================================================
    ## make proper TH3 function 
    def make_thfun ( self , histo , xvar , yvar , zvar ) :
        """ Create proper TH3-funcion
        """
        return Ostap.Functions.FuncTH3 ( histo , xvar , yvar , zvar )
    




# =============================================================================
## reconstrucct SPlot4Tree object
def _sp4t_factory_ ( *args ) :
    """ Reconstrucct SPlot4Tree object
    """
    result = Ostap.MoreRooFit.SPlot4Tree ( *args )
    result._args = args
    return result
# =============================================================================
## reduce SPlit3Tree object 
def _sp4t_reduce_  ( sp4t   ) :
    """ Reduce SPlit3Tree object
    """
    return _sp4t_factory_ , ( sp4t.pdf           () ,
                              sp4t.observables   () ,
                              sp4t.fitresult     () ,                            
                              sp4t.normalization () )

Ostap.MoreRooFit.SPlot4Tree.__reduce__ = _sp4t_reduce_

## =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


    
# =============================================================================
##                                                                      The END 
# =============================================================================
