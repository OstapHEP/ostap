#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostap/fitting/simfit.py
#  Collection of utilities that simplify the Simultaneous fit
#  @see RooSimultaneous
#  @see https://github.com/OstapHEP/ostap/blob/master/ostap/fitting/SIMFIT.md
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2018-11-23
# =============================================================================
""" Collection of utilities that simplify the Simultaneous fit
- see https://github.com/OstapHEP/ostap/blob/master/ostap/fitting/SIMFIT.md
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    'SimFit'         , ## fit model for simultaneon     fit
    'Sim1D'          , ## fit model for simultaneous 1D-fit (obsolete) 
    'combined_data'  , ## prepare combined dataset for the simultaneous fit
    'combined_hdata' , ## prepare combined binned dataset for the simultaneous fit
    )
# =============================================================================
from   ostap.core.core          import std , Ostap , dsID , items_loop
from   ostap.utils.utils        import chunked 
from   ostap.fitting.fithelpers import VarMaker, ConfigReducer
from   ostap.fitting.pdfbasic   import ( PDF1 , Generic1D_pdf , 
                                         PDF2 , Generic2D_pdf , 
                                         PDF3 , Generic3D_pdf )
import ostap.fitting.variables 
import ostap.fitting.roocmdarg
import ROOT, math,  random , warnings 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.simfit' )
else                       : logger = getLogger ( __name__               )

# =============================================================================
## Create combined dataset for simultaneous fit
#  @code
#  sample = ROOT.RooCategory ( 'sample' , 'sample' , 'cc' , 'zz' )
#  vars   = ROOT.RooArgSet   ( m2c )
#  ds_cmb = combined_data    ( sample  ,
#                vars    , { 'cc' : ds_cc ,  'zz' : ds_00 } )
#  @endcode
#  - weighted variant: combine unweighted datasets and then apply weight 
#  @code
#  wvars = ROOT.RooArgSet ( m2c , SS_sw ) 
#  dsw_cmb   = combined_data ( sample ,
#                 wvars  , { 'cc' : dsn_cc ,  'zz' : dsn_00 } ,
#                 args = ( ROOT.RooFit.WeightVar( 'SS_sw' ) , ) )
#  @endcode
def combined_data ( sample          ,
                    varset          , 
                    datasets        ,
                    name     = ''   ,
                    title    = ''   ,
                    args     = ()   ) :
    """
     Create combined  dataset for simultaneous fit

     >>> sample = ROOT.RooCategory ( 'sample' , 'sample' , 'cc' , 'zz' )
     >>> vars   = ROOT.RooArgSet   ( m2c )
     >>> ds_cmb = combined_data ( sample  ,
     ...          vars    , { 'cc' : ds_cc ,  'zz' : ds_00 } )
     
     Weighted variant:
     
     >>> wvars = ROOT.RooArgSet ( m2c , SS_sw ) 
     >>> dsw_cmb   = combined_data ( sample ,
     ...             wvars  , { 'cc' : dsn_cc ,  'zz' : dsn_00 } ,
     ...             args = ( ROOT.RooFit.WeightVar( 'SS_sw' ) , ) )
     
     """
    
    labels  = sample.labels()
    
    largs   = [ ROOT.RooFit.Index ( sample ) ] 

    weights = set() 
    ds_keep = []
    
    store_error      = False 
    store_asym_error = False 

    ## collect all vars 
    all_vars = ROOT.RooArgSet()
    for label in labels :

        dset = None 
        if isinstance ( datasets , dict ) : dset = datasets[label]
        else :
            for ds in datasets :
                if label == ds [ 0 ] :
                    dset =  ds [ 1 ]
                    break
                
        assert dset and isinstance ( dset , ROOT.RooAbsData ),\
               'Invalid data set for label %s' % label

        for v in dset.vars :
            if not v in all_vars : all_vars.add ( v ) 

        ##assert not dset.isNonPoissonWeighted () ,\
        ##       'Weighted data cannot be combined!'

        if not dset.isWeighted () :
            largs.append (  ROOT.RooFit.Import ( label , dset ) )
        else :
            store_error      = dset.store_error      () 
            store_asym_error = dset.store_asym_error ()            
            uwdset , wnam    = dset.unWeighted ()
            assert uwdset and wnam, "Cannot 'unweight' dataset!"
            largs.append   ( ROOT.RooFit.Import ( label , uwdset ) )
            ds_keep.append ( uwdset ) 
            weights.add    ( wnam   )
            
    ## we can do with only one weight variable!
    assert len ( weights ) < 2 , 'Invalid number of weights %s' % list ( weights )

    weight = weights.pop() if weights else None 
    
    name  = name  if name  else dsID()
    title = title if title else 'Data for simultaneous fit/%s' % sample.GetName()

    args = args + tuple ( largs )
        
    vars = ROOT.RooArgSet()
    if    isinstance ( varset , ROOT.RooAbsReal     ) : vars.add ( varset )
    elif  isinstance ( varset , ROOT.RooAbsCategory ) : vars.add ( varset )
    else :
        for i , v in enumerate ( varset ) :
            if   isinstance ( v , ROOT.RooAbsReal     ) : vars.add ( v )
            elif isinstance ( v , ROOT.RooAbsCategory ) : vars.add ( v )
            elif v in all_vars :
                vv = all_vars [ v ]
                vars.add ( vv )
            else : 
                assert TypeError ( 'Invalid variable [%d] %s/%s' % ( i , v , type ( v ) ) )

    wset = ROOT.RooArgSet() 
    if weight :
        args = args + ( ROOT.RooFit.WeightVar ( weight ) , )
        if not weight in vars : 
            wvar = ROOT.RooRealVar ( weight , 'weight variable' , 1 , -1.e+100 , 1.e+100 ) 
            vars.add ( wvar )

        wvar = vars [ weight ]
        wset.add ( wvar )

        if   store_asym_error : args = args + ( ROOT.RooFit.StoreAsymError ( wset ) , )
        elif store_error      : args = args + ( ROOT.RooFit.StoreError     ( wset ) , )

    NARGS = 8
    if len ( args ) < NARGS :  
        ds = ROOT.RooDataSet ( name , title , vars , *args )
    else :
        keep  = set()
        for a in args : keep.add ( a )
        cargs = tuple ( a for a in args )
        while NARGS <= len ( cargs )  :
            keep.add ( cargs )
            nargs  = []
            for chunk in chunked ( cargs , min ( 5 , NARGS ) ) :
                newarg = ROOT.RooFit.MultiArg ( *chunk )
                nargs.append ( newarg )
            cargs = tuple ( nargs )
        ds = ROOT.RooDataSet ( name , title , vars , *cargs )

    ## delete helper intermediate datasets 
    while ds_keep :
        d = ds_keep.pop()
        d.reset()
        del d
            
    return ds


# =============================================================================
## create combined binned dataset for simultaneous fit as RooDataHist 
#  - combine 2D histograms:
#  @code
#  sample = ROOT.RooCategory ( 'sample' , 'fitting sample' , 'A' , 'B' )
#  hA  = ...
#  hB  = ...
#  var = ROOT.RooRealVar ( ... )
#  ds  = combined_hdata ( sample , var , { 'A' : hA , 'B' : hB } )  
#  @endcode
#  - combine 2D histograms:
#  @code
#  sample = ROOT.RooCategory ( 'sample' , 'fitting sample' , 'A' , 'B' )
#  hA   = ...
#  hB   = ...
#  xvar = ROOT.RooRealVar ( ... )
#  yvar = ROOT.RooRealVar ( ... )
#  ds   = combined_hdata ( sample , (xvar, yvar) , { 'A' : hA , 'B' : hB } )  
#  @endcode
#  - combine 3D histograms:
#  @code
#  sample = ROOT.RooCategory ( 'sample' , 'fitting sample' , 'A' , 'B' )
#  hA   = ...
#  hB   = ...
#  xvar = ROOT.RooRealVar ( ... )
#  yvar = ROOT.RooRealVar ( ... )
#  zvar = ROOT.RooRealVar ( ... )
#  ds   = combined_hdata ( sample , (xvar, yvar, zvar) , { 'A' : hA , 'B' : hB } )  
#  @endcode
#  Since the output is <code>RooDataHist</code> the specifications for
#  the individual histograms must be the same  
def combined_hdata ( sample        ,
                     varset        ,
                     histograms    ,
                     name     = '' ,
                     title    = '' ) :
    """Create combined binned dataset for simultaneous fit
    - combine 2D histograms:

    >>> sample = ROOT.RooCategory ( 'sample' , 'fitting sample' , 'A' , 'B' )
    >>> hA  = ...
    >>> hB  = ...
    >>> var = ROOT.RooRealVar ( ... )
    >>> ds  = combined_hdata ( sample , var , { 'A' : hA , 'B' : hB } )  

    - combine 2D histograms:

    >>> sample = ROOT.RooCategory ( 'sample' , 'fitting sample' , 'A' , 'B' )
    >>> hA   = ...
    >>> hB   = ...
    >>> xvar = ROOT.RooRealVar ( ... )
    >>> yvar = ROOT.RooRealVar ( ... )
    >>> ds   = combined_hdata ( sample , (xvar, yvar) , { 'A' : hA , 'B' : hB } )  

    - combine 3D histograms:

    >>> sample = ROOT.RooCategory ( 'sample' , 'fitting sample' , 'A' , 'B' )
    >>> hA   = ...
    >>> hB   = ...
    >>> xvar = ROOT.RooRealVar ( ... )
    >>> yvar = ROOT.RooRealVar ( ... )
    >>> zvar = ROOT.RooRealVar ( ... )
    >>> ds   = combined_hdata ( sample , (xvar, yvar, zvar) , { 'A' : hA , 'B' : hB } )  

    - Since the output is `RooDataHist` the specifications for the individual histograms must be the same  
    """

    MAP  = std.map  ( 'std::string'       , 'TH1*' )
    PAIR = std.pair ( 'const std::string' , 'TH1*' )
    mm   = MAP()
    
    nd1   = 0 
    nd2   = 0 
    nd3   = 0
    
    labels = sample.labels ()
    
    for label in labels :
        
        histo = histograms.pop  ( label ) 
        
        if   isinstance ( histo , ROOT.TH3 ) : nd3 += 1
        elif isinstance ( histo , ROOT.TH2 ) : nd2 += 1
        elif isinstance ( histo , ROOT.TH1 ) : nd1 += 1
        
        mm.insert ( PAIR ( label , histo ) )

    assert not histograms, 'Unknown histograms: %s' % histograms.keys() 

    d1 = 0 < nd1
    d2 = 0 < nd2
    d3 = 0 < nd3
    
    ## assert ( not d1 and not d2 ) or \
    ##       ( not d2 and not d3 ) or \
    ##       ( not d3 and not d1 ) , 'Mismatch in histogram dimensions!'
            
    name  = name  if name  else dsID()
    title = title if title else 'Data for simultaneous fit/%s' % sample.GetName()

    varlst = ROOT.RooArgList()
    if   isinstance ( varset , ROOT.RooAbsReal     ) : varlst.add ( varset )
    elif isinstance ( varset , ROOT.RooAbsCategory ) : varlst.add ( varset )
    else :  
        for i , v in enumerate ( varset ) :
            varlst.add ( v )
            
    ## assert ( d3 and 3 == len ( varlst ) ) or \
    ##       ( d2 and 2 == len ( varlst ) ) or \
    ##       ( d1 and 1 == len ( varlst ) )  , \
    ##       'Invalid dimension of dataset!'
    
    return ROOT.RooDataHist ( name , title , varlst , sample  , mm ) 
    
# =============================================================================        
## @class Sim1D
#  Helper class to simplify creation and usage of simultaneous PDF
#  @code
#  @endcode 
#  @see RooSimultaneous
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2018-11-23
class Sim1D(PDF1) :
    """Helper class to simplify the creation and usage of simultaneous PDF
    
    - see RooSimultaneous 
    """
    
    def __init__ ( self              ,
                   sample            , 
                   categories        ,
                   xvar       = None ,
                   name       = None , 
                   title      = ''   ) :
        
        warnings.warn ( "Usage of obsolete Sim1D. Use SimFit instead!" )
        
        if isinstance ( sample , ( tuple , list ) ) :
            _cat = ROOT.RooCategory ( 'sample' , 'sample' )
            for i in sample : _cat.defineType ( i ) 
            sample =  _cat
            
        assert isinstance ( sample , ROOT.RooCategory ),\
               'Invalid type for "sample":' % ( sample ,  type ( sample ) )
        
        name = name if name else self.generate_name ( prefix = 'sim1D_%s_' % sample.GetName() )
        
        self.__sample     = sample 
        self.__categories = {}

        # =====================================================================
        ## components
        # =====================================================================
        labels = sample.labels()

        _xvar = xvar 
        for label in labels :
            
            cmp   = None 
            if isinstance ( categories , dict ) : cmp = categories [ label ]
            else :
                for ii in categories :
                    if ii[0] == label :
                        cmp = ii[1]
                        break
                    
            if   isinstance ( cmp , PDF1 ) :
                
                _xv = cmp.xvar
                if _xvar and not ( _xvar is _xv ) : self.error('Mismatch in XVAR!')
                elif not _xvar                    : _xvar = _xv
                
                self.__categories [ label ] = cmp
                
            elif isinstance ( cmp , ROOT.RooAbsPdf ) and _xvar :
                
                self.__categories [ label ] = Generic1D_pdf ( pdf = cmp , cmp = _xvar )
                
            else :
                
                raise TypeError ( 'Can not find the category "%s"' % label ) 

        # =====================================================================
        assert _xvar, 'Unable to define "xvar"'

        ## initialize the base 
        PDF1.__init__ ( self , name , xvar = _xvar ) 
        
        self.warning  ( "Usage of obsolete Sim1D. Use SimFit instead!" )


        self.pdf = ROOT.RooSimultaneous (
            self.roo_name ( 'sim1d_' ) ,
            title if title else "Simultaneous %s" % self.name , 
            self.sample )

        keys = self.categories.keys()
        for key in sorted ( keys ) :
            self.pdf.addPdf ( self.categories[key].pdf , key )

        for k , pdf in items_loop ( self.categories ) :

            self.copy_structures ( pdf )
            
            ## for c in pdf.signals              : self.signals             .add ( c ) 
            ## for c in pdf.backgrounds          : self.backgrounds         .add ( c ) 
            ## for c in pdf.crossterms1          : self.crossterms1         .add ( c ) 
            ## for c in pdf.crossterms2          : self.crossterms2         .add ( c ) 
            ## for c in pdf.combined_signals     : self.combined_signals    .add ( c ) 
            ## for c in pdf.combined_backgrounds : self.combined_background .add ( c ) 
            ## for c in pdf.combined_components  : self.combined_components .add ( c ) 

            ## ## copy "components" (needed for e.g. splot) 
            ## for c in pdf.alist1      : self.alist1.add ( c ) 
            ## for c in pdf.alist2      : self.alist2.add ( c ) 

            ## copy draw options 
            for k in pdf.draw_options :
                self.draw_options [ k ] = pdf.draw_options [ k ]


        self.config = {
            'name'       : self.name       ,
            'title'      : title           ,            
            'sample'     : self.sample     ,
            'categories' : self.categories ,
            'xvar'       : self.xvar       ,
            }

    
    @property
    def sample  ( self ) :
        "'sample': RooCategory object for simultaneous PDF"
        return self.__sample

    @property
    def samples ( self ) :
        "'samples' : list/tuple of the known categories"
        return tuple ( self.__categories.keys() ) 

    @property
    def categories ( self ) :
        "'categories' : map { category : pdf }"
        return self.__categories

    # =========================================================================
    ## delegate  attribute search to the components 
    def __getattr__ ( self , attr ) :
        """Delegate attribte search to the category components
        """
        if attr in self.samples : return self.components[attr]
        
        raise  AttributeError('Unknown attibute %s' % attr )

    # =========================================================================
    ## make the actual fit (and optionally draw it!)
    #  @code
    #  r,f = model.fitTo ( dataset )
    #  r,f = model.fitTo ( dataset , weighted = True )    
    #  r,f = model.fitTo ( dataset , ncpu     = 10   )    
    #  r,f = model.fitTo ( dataset , draw = 'signal' , nbins = 300 )    
    #  @endcode 
    def fitTo ( self           ,
                dataset        ,
                draw   = False ,
                nbins  = 100   ,
                silent = False ,
                refit  = False ,
                timer  = False ,
                args   = ()    , **kwargs ) :
        """
        Perform the actual fit (and draw it)
        >>> r,f = model.fitTo ( dataset )
        >>> r,f = model.fitTo ( dataset , weighted = True )    
        >>> r,f = model.fitTo ( dataset , ncpu     = 10   )    
        >>> r,f = model.fitTo ( dataset , draw = 'signal' , nbins = 300 )    
        """
        assert self.sample in dataset      ,\
               'Category %s is not in dataset' % self.sample.GetName()

        res , frame = PDF1.fitTo ( self ,
                                   dataset = dataset ,
                                   draw    = False   , ## ATTENTION! False is here! 
                                   nbins   = nbins   ,
                                   silent  = silent  ,
                                   refit   = refit   ,
                                   timer   = timer   , 
                                   args    = args    , **kwargs )
        
        if   not draw                 : return res , None 
        elif not draw in self.samples :
            self.error ('Unknown category for drawing %s' % draw )
            return res , None
        
        from ostap.plotting.fit_draw import draw_options
        draw_opts = draw_options ( **kwargs )

        frame = self.draw ( category = draw    ,
                            dataset  = dataset ,
                            nbins    = nbins   ,
                            silent   = silent  , **draw_opts )
        
        return res , frame
    
    # ========================================================================
    ## Draw the PDF&data for the given category
    #  @code
    #  pdf.fitTo ( dataset )
    #  pdf.draw ( 'signal' , dataset , nbins = 100 ) 
    #  @endcode 
    def draw ( self           ,
               category       ,  ## must be specified! 
               dataset        ,  ## must be specified!
               nbins   =  100 ,
               silent  = True ,
               args    = ()   , 
               **kwargs       ) :
        """
        Draw the PDF&data for the given   category
        >>> pdf.fitTo ( dataset )
        >>> pf.draw ( 'signal' , dataset , nbins = 100 ) 
        """
        
        assert category    in self.samples ,\
               'Category %s is not in %s' % ( category , self.samples )
        assert self.sample in dataset      ,\
               'Category %s is not in dataset' % self.sample.GetName()
                
        ## 
        sname = self.sample.GetName() 
        dcut  = ROOT.RooFit.Cut ( "%s==%s::%s"  % ( sname , sname , category ) )
        
        data_options = self.draw_option ( 'data_options' , **kwargs ) +  ( dcut , ) 
        
        self._tmp_vset = ROOT.RooArgSet ( self.sample ) 
        _proj  = ROOT.RooFit.ProjWData  ( self._tmp_vset , dataset  ) 
        _slice = ROOT.RooFit.Slice      ( self.sample    , category )

        bkgoptions   = self.draw_option ( 'backrground_options' , **kwargs ) + ( _slice , _proj )
        ct1options   = self.draw_option ( 'crossterm1_options'  , **kwargs ) + ( _slice , _proj )
        ct2options   = self.draw_option ( 'crossterm2_options'  , **kwargs ) + ( _slice , _proj )        
        cmpoptions   = self.draw_option (  'component_options'  , **kwargs ) + ( _slice , _proj )
        sigoptions   = self.draw_option (     'signal_options'  , **kwargs ) + ( _slice , _proj )
        totoptions   = self.draw_option (  'total_fit_options'  , **kwargs ) + ( _slice , _proj )
        
        kwargs [ 'data_options'       ] = data_options
        kwargs [ 'signal_options'     ] = sigoptions 
        kwargs [ 'background_options' ] = bkgoptions 
        kwargs [ 'crossterm1_options' ] = ct1options 
        kwargs [ 'crossterm2_options' ] = ct2options
        kwargs [ 'component_options'  ] = cmpoptions 
        kwargs [ 'total_fit_options'  ] = totoptions 

        from ostap.fitting.roocollections import KeepArgs

        cat_pdf = self.categories[ category ]
        with KeepArgs     ( self . signals     , cat_pdf . signals     ) as _k1 ,\
                 KeepArgs ( self . backgrounds , cat_pdf . backgrounds ) as _k2 ,\
                 KeepArgs ( self . components  , cat_pdf . components  ) as _k3 ,\
                 KeepArgs ( self . crossterms1 , cat_pdf . crossterms1 ) as _k4 ,\
                 KeepArgs ( self . crossterms2 , cat_pdf . crossterms2 ) as _k5 :
            
            return PDF1.draw ( self ,
                               dataset = dataset ,
                               nbins   = nbins   ,
                               silent  = silent  ,
                               args    = args    , **kwargs )
        
        
# =============================================================================
## suppress methods specific for 1D-PDFs only
for _a in (
    ##'_get_stat_'   ,
    'integral'       , 
    'histo'          , 
    'roo_histo'      , 
    'residual'       , 
    'pull'           ,
    #
    'rms'            , 
    'fwhm'           , 
    'skewness'       , 
    'kurtosis'       , 
    'mode'           , 
    'mode'           , 
    'median'         , 
    'get_mean'       , 
    'moment'         , 
    'central_moment' , 
    'quantile'       , 
    'cl_symm'        , 
    'cl_asymm'       ,
    'derivative'     ) :

    if hasattr ( Sim1D , _a ) :
        ## method is suppressed
        def _suppress_ ( self , *args , **kwargs ) :
            """Method is suppressed"""
            raise AttributeError ( "'%s' object has no attribute '%s'" % ( type ( self ) , _a ) )
        setattr ( Sim1D , _a , _suppress_ ) 
        logger.verbose ( 'Remove attribute %s from Sim1D' ) 


# =============================================================================        
## @class SimFit
#  Helper pdf-like class to simplify creation and usage of simultaneous PDF
#  @code
#  sample = ROOT.RooCategory( 'sample', 'fitting sample' , 'A' , 'B' )
#  pdfA   = ... ## pdf for the sample 'A'
#  pdfB   = ... ## pdf for the sample 'B'
#  simfit = SimFit (  sample , { 'A' : pdfA , 'B' : pdfB } )
#  @endcode 
#  @see RooSimultaneous
#  Note that this class is *not* PDF, but it behaves rather similar to PDF,
#  and, in partcualr has such methods like
#  - <code>fitTo</code>
#  - <code>draw</code>
#  - <code>nll</code>
#  - <code>draw_nll</code>
#  - <code>wilks</code>
#  - <code>minuit</code> 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2018-11-23
class SimFit (VarMaker,ConfigReducer) :
    """Helper PDF-like class to simplify the creation and usage of simultaneous PDF
    
    >>> sample = ROOT.RooCategory( 'sample', 'fitting sample' , 'A' , 'B' )
    >>> pdfA   = ... ## pdf for the sample 'A'
    >>> pdfB   = ... ## pdf for the sample 'B'
    >>> simfit = SimFit (  sample , { 'A' : pdfA , 'B' : pdfB } )

    - see RooSimultaneous

    Note that this class is *not* PDF, but it behaves rather similar to PDF,
    and, in partcular has such methods as 
    - fitTo
    - draw
    - nll
    - draw_nll
    - wilks
    - minuit
    """    
    def __init__ ( self              ,
                   sample            , 
                   categories        ,
                   name       = None , 
                   title      = ''   ) :
        """Helper PDF-like class to simplify the creation and usage of simultaneous PDF
        
        >>> sample = ROOT.RooCategory( 'sample', 'fitting sample' , 'A' , 'B' )
        >>> pdfA   = ... ## pdf for the sample 'A'
        >>> pdfB   = ... ## pdf for the sample 'B'
        >>> simfit = SimFit (  sample , { 'A' : pdfA , 'B' : pdfB } )
        """

        ## initialize the base class 
        ConfigReducer.__init__ ( self )
        
        self.__categories   = {}
        
        if isinstance ( sample , ( tuple , list ) ) :
            _cat = ROOT.RooCategory ( 'sample' , 'sample' )
            for i in sample : _cat.defineType ( i ) 
            sample =  _cat
            
        assert isinstance ( sample , ROOT.RooCategory ),\
               'Invalid type for "sample":' % ( sample ,  type ( sample ) )
        
        name = name if name else self.generate_name ( 'simfit' , '' , sample.GetName() )
        
        ## propagate the name 
        self.name = name
        
        self.__sample       = sample 

        ## observables 
        self.__observables  = ROOT.RooArgSet( self.__sample )
        
        # =====================================================================
        ## components
        # =====================================================================
        labels = sample.labels()

        _xv = None 
        for label in labels :
            
            cmp   = None 
            if isinstance ( categories , dict ) : cmp = categories [ label ]
            else :
                for ii in categories :
                    if ii[0] == label :
                        cmp = ii[1]
                        break

            if not isinstance  ( cmp , ( PDF1 , PDF2 , PDF3 ) ) :
                raise TypeError ( 'Can not find the proper category component: "%s"' % label ) 
            
            self.__categories [ label ] = cmp
            _xv = cmp.xvar

        sim_pdf     = PDF1 ( self.name + '_Sim' , xvar = _xv )            
        sim_pdf.pdf = ROOT.RooSimultaneous (
            self.roo_name ( 'simfit_' ) ,
            title if title else "Simultaneous %s" % self.name , 
            self.sample )
        
        ## get all observables  
        for label in self.categories :
            pdf = self.categories[ label ]
            for v in pdf.observables :
                if not v in self.__observables :
                    self.__observables.add ( v )
              
        keys = self.categories.keys()
        for key in sorted ( keys ) :
            sim_pdf.pdf.addPdf ( self.categories[key].pdf , key )

        self.__pdf = sim_pdf 
        
        for k , cmp in items_loop ( self.categories ) :

            self.pdf.copy_structures ( cmp )
            
            ## for c in cmp.signals              : self.pdf.signals             .add ( c ) 
            ## for c in cmp.backgrounds          : self.pdf.backgrounds         .add ( c ) 
            ## for c in cmp.crossterms1          : self.pdf.crossterms1         .add ( c ) 
            ## for c in cmp.crossterms2          : self.pdf.crossterms2         .add ( c )
            ## for c in cmp.combined_signals     : self.pdf.combined_signals    .add ( c )
            ## for c in cmp.combined_backgrounds : self.pdf.combined_backgrounds.add ( c )
            ## for c in cmp.combined_components  : self.pdf.combined_components .add ( c )
            
            ## ## copy "components" (needed for e.g. splot) 
            ## for c in cmp.alist1               : self.pdf.alist1 . add ( c ) 
            ## for c in cmp.alist2               : self.pdf.alist2 . add ( c )
            
            self.pdf.draw_options.update ( cmp.draw_options )

        # =====================================================================
        ##  drawing helpers
        # =====================================================================
        
        self.__drawpdfs   = {}
        for key in sorted ( keys ) :

            cmp = self.categories [ key ] 
            if isinstance  ( cmp , PDF3 ) :
                dpdf = Generic3D_pdf ( sim_pdf.pdf , 
                                       cmp.xvar    ,
                                       cmp.yvar    ,
                                       cmp.zvar    ,
                                       name           = sim_pdf.name + '_draw_' + key , 
                                       add_to_signals = False )
            elif isinstance  ( cmp , PDF2 ) :
                dpdf = Generic2D_pdf ( sim_pdf.pdf , 
                                       cmp.xvar    ,
                                       cmp.yvar    ,
                                       name           = sim_pdf.name + '_draw_' + key , 
                                       add_to_signals = False )
            elif isinstance  ( cmp , PDF1 ) :
                dpdf = Generic1D_pdf ( sim_pdf.pdf , 
                                       cmp.xvar    ,
                                       name           = sim_pdf.name + '_draw_' + key , 
                                       add_to_signals = False  )
                
            ## copy structural elements 
            dpdf.copy_structures ( cmp )
            ## update copy draw options 
            dpdf.draw_options.update ( cmp.draw_options )
                        
            self.__drawpdfs [ key ]  = dpdf

        self.config = {
            'name'       : self.name       ,
            'title'      : title           ,            
            'sample'     : self.sample     ,
            'categories' : self.categories ,
            }

    
    @property
    def sample  ( self ) :
        "'sample' : RooCategory object for simultaneous PDF"
        return self.__sample

    @property
    def pdf     ( self ) :
        "'pdf'  : the actual PDF with RooSimultaneous "
        return self.__pdf

    @property
    def simpdf  ( self ) :
        """'simpdf' : the actual RooSimultaneous instance """
        return self.__pdf.pdf 
        
    @property
    def roo_pdf ( self ) :
        """'roo_pdf' : get the underlying RooFit RooAbsPdf object (same as `simpdf` here)"""
        return self.pdf.pdf 

    @property
    def samples ( self ) :
        "'samples' : list/tuple of known categories"
        return tuple ( self.__categories.keys() ) 

    @property
    def categories ( self ) :
        "'categories' : map { category : pdf } "
        return self.__categories

    @property
    def observables ( self ) :
        """'observables' : get the observables (as `ROOT.RooArgSet`, same as `vars`)"""
        return self.__observables
    
    @property
    def vars        ( self ) :
        """'vars' : get the observables (as `ROOT.RooArgSet`, same as `obsevables`)"""
        return self.observables
            
    @property
    def drawpdfs ( self ) :
        "'drawpdfs' dictionary with PDFs (mainly for drawing)"
        return self.__drawpdfs

    @property
    def draw_options ( self ) :
        """`draw_options' : dictionary with predefined draw-options for this PDF
        """
        return self.pdf.draw_options

    # =========================================================================
    ## get the certain predefined drawing option
    #  @code
    #  options = ROOT.RooFit.LineColor(2), ROOT.RooFit.LineWidth(4)
    #  pdf = ...
    #  pdf.draw_options['signal_style'] = [ options ]
    #  ## and later:
    #  options = pdf.draw_option ( 'signal_style' )
    #  @endcode 
    def draw_option ( self , key , default = () , **kwargs ) :
        """Get the certain predefined drawing option
        >>> options = ROOT.RooFit.LineColor(2), ROOT.RooFit.LineWidth(4)
        >>> pdf = ...
        >>> pdf.draw_options['signal_style'] = options 
        - and later:
        >>> options = pdf.draw_option ( 'signal_style' )
        """
        return self.pdf.draw_option ( key , default , **kwargs ) 

    ## # =========================================================================
    ## ## delegate  attribute search to the components 
    ## def __getattr__ ( self , attr ) :
    ##     """Delegate attribute search to the category components
    ##     """
    ##     print ( 'CHECKIING', attr ) 
    ##     if attr in self.__categories :
    ##         return self.__categories[attr]
        
    ##     raise  AttributeError('Unknown attibute %s' % attr )

    # =========================================================================
    ## make the actual fit (and optionally draw it!)
    #  @code
    #  r,f = model.fitTo ( dataset )
    #  r,f = model.fitTo ( dataset , weighted = True )    
    #  r,f = model.fitTo ( dataset , ncpu     = 10   )    
    #  r,f = model.fitTo ( dataset , draw = 'signal' , nbins = 300 )    
    #  @endcode 
    def fitTo ( self           ,
                dataset        ,
                draw   = False ,
                nbins  = 100   ,
                silent = False ,
                refit  = False ,
                timer  = False ,
                args   = ()    , **kwargs ) :
        """
        Perform the actual fit (and draw it optionally)
        >>> r,f = model.fitTo ( dataset )
        >>> r,f = model.fitTo ( dataset , weighted = True )    
        >>> r,f = model.fitTo ( dataset , ncpu     = 10   )    
        >>> r,f = model.fitTo ( dataset , draw = 'signal' , nbins = 300 )    
        """
        assert self.sample in dataset      ,\
               'Category %s is not in dataset' % self.sample.GetName()

        res , frame = self.pdf.fitTo ( 
            dataset = dataset ,
            draw    = False   , ## ATTENTION! False is here! 
            nbins   = nbins   ,
            silent  = silent  ,
            refit   = refit   ,
            timer   = timer   , 
            args    = args    , **kwargs )
        
        if   not draw                  : return res , None
        elif draw in self.samples      : pass
        elif isinstance ( draw , str ) :
            draw , s , dvar = draw.partition('/')
            if not draw in self.samples : 
                self.error ('Unknown category for drawing %s' % draw )
                return res , None
        
        from ostap.plotting.fit_draw import draw_options
        draw_opts = draw_options ( **kwargs )

        frame = self.draw ( category = draw    ,
                            dataset  = dataset ,
                            nbins    = nbins   ,
                            silent   = silent  , **draw_opts )
        
        return res , frame

    
    # ========================================================================
    ## Draw the PDF&data for the given category
    #  @code
    #  pdf.fitTo ( dataset )
    #  pdf.draw ( 'signal' , dataset , nbins = 100 ) 
    #  @endcode 
    def draw ( self           ,
               category       ,  ## must be specified! 
               dataset        ,  ## must be specified!
               nbins   =  100 ,
               silent  = True ,
               args    = ()   , 
               **kwargs       ) :
        """
        Draw the PDF&data for the given category
        >>> pdf.fitTo ( dataset )
        >>> pf.draw ( 'signal' , dataset , nbins = 100 ) 
        """

        dvar = None
        if   isinstance ( category , ( tuple , list ) ) and 2 == len ( category ) :
            category , dvar = category 
        elif isinstance ( category , str ) :
            category , _ , dvar = category.partition('/')
            
        if dvar :
            try :
                ivar = int ( dvar )
                dvar = ivar 
            except ValueError :
                pass

        assert category    in self.samples ,\
               'Category %s is not in %s' % ( category , self.samples )
        assert self.sample in dataset      ,\
               'Category %s is not in dataset' % self.sample.GetName()
                
        ## 
        sname = self.sample.GetName() 
        dcut  = ROOT.RooFit.Cut ( "%s==%s::%s"  % ( sname , sname , category ) )

        kwargs [ 'data_options' ] = self.draw_option ( 'data_options' , **kwargs ) +  ( dcut , )
        
        self._tmp_vset = ROOT.RooArgSet ( self.sample ) 

        _proj  = ROOT.RooFit.ProjWData  ( self._tmp_vset , dataset  ) 
        _slice = ROOT.RooFit.Slice      ( self.sample    , category )
        for key in  ( 'total_fit_options'           ,
                      #
                      'signal_options'              ,
                      'background_options'          ,
                      'component_options'           ,
                      'crossterm1_options'          ,
                      'crossterm2_options'          ,
                      #
                      'combined_signal_options'     ,
                      'combined_background_options' ,
                      'combined_component_options'  ) :

            kwargs [ key ] = self.draw_option ( key , **kwargs ) + ( _slice , _proj ) 
        
        from ostap.fitting.roocollections import KeepArgs

        cat_pdf  = self.categories [ category ]
        draw_pdf = self.drawpdfs   [ category ]

        if 1 < 2 : 
        ## with KeepArgs     ( draw_pdf . signals     , cat_pdf . signals     ) as _k1 ,\
        ##          KeepArgs ( draw_pdf . backgrounds , cat_pdf . backgrounds ) as _k2 ,\
        ##          KeepArgs ( draw_pdf . components  , cat_pdf . components  ) as _k3 ,\
        ##          KeepArgs ( draw_pdf . crossterms1 , cat_pdf . crossterms1 ) as _k4 ,\
        ##          KeepArgs ( draw_pdf . crossterms2 , cat_pdf . crossterms2 ) as _k5 :

            if   isinstance ( draw_pdf , PDF3 ) :

                if   3 == dvar or dvar in  ( 'z' , 'Z' , '3' , draw_pdf.zvar.name ) :    
                    return draw_pdf.draw3 ( dataset = dataset ,
                                            nbins   = nbins   ,
                                            silent  = silent  ,
                                            args    = args    , **kwargs )
                elif 2 == dvar or dvar in  ( 'y' , 'Y' , '2' , draw_pdf.yvar.name ) : 
                    return draw_pdf.draw2 ( dataset = dataset ,
                                            nbins   = nbins   ,
                                            silent  = silent  ,
                                            args    = args    , **kwargs )
                elif 1 == dvar or dvar in  ( 'x' , 'X' , '1' , draw_pdf.xvar.name ) : 
                    return draw_pdf.draw1 ( dataset = dataset ,
                                            nbins   = nbins   ,
                                            silent  = silent  ,
                                            args    = args    , **kwargs )
                else :
                    self.error("Unknown 'dvar' for 3D-draw pdf!")
                    return None
                
            elif isinstance ( draw_pdf , PDF2 ) :
                
                if   2 == dvar or dvar in  ( 'y' , 'Y' , '2' , draw_pdf.yvar.name ) :
                    return draw_pdf.draw2 ( dataset = dataset ,
                                            nbins   = nbins   ,
                                            silent  = silent  ,
                                            args    = args    , **kwargs )
                elif 1 == dvar or dvar in  ( 'x' , 'X' , '1' , draw_pdf.xvar.name ) : 
                    return draw_pdf.draw1 ( dataset = dataset ,
                                            nbins   = nbins   ,
                                            silent  = silent  ,
                                            args    = args    , **kwargs )
                else :
                    self.error("Unknown 'dvar' for 2D-draw pdf! %s" %  dvar )
                    return None 

            elif isinstance ( draw_pdf , PDF1 ) :                
                return draw_pdf.draw ( dataset = dataset ,
                                       nbins   = nbins   ,
                                       silent  = silent  ,
                                       args    = args    , **kwargs )

            self.error ("draw: inconsistent combination of draw_pdf '%s' and dvar '%s'" % (
                type ( draw_pdf ) , dvar ) )
            
    # =========================================================================
    ## create NLL
    #  @code
    #  model.fitTo ( dataset , ... )
    #  nll, sfactor  = model.nll ( 'dataset )
    #  @endcode
    #  @see RooAbsPdf::createNLL
    #  @attention Due to the bug/typo in<c> RooAbsPdf.clreateNLL</c>, line 817 
    #  <c>CloneData</c> depends on <c>Optimize</c>
    #  @todo report problem to RooFit and fix it! 
    def nll ( self            ,
              dataset         ,
              silent  = True  ,
              args    = ()    , **kwargs ) :
        """Get NLL object from the pdf
        >>> model.fitTo ( dataset , ... )
        >>> nll, sf = model.nll ( dataset )
        - see RooAbsPdf::createNLL 
        """
        return self.pdf.nll ( dataset , silent = silent , args = args , **kwargs )

    # =========================================================================
    ## draw/prepare NLL or LL-profiles for selected variable
    #  @code
    #  model.fitTo ( dataset , ... )
    #  nll  , f1 = model.draw_nll ( 'B' ,  dataset )
    #  prof , f2 = model.draw_nll ( 'B' ,  dataset , profile = True )
    #  @endcode    
    def draw_nll ( self            ,
                   var             ,
                   dataset         ,
                   profile = False ,
                   draw    = True  ,
                   silent  = True  , 
                   args    = ()    , **kwargs ) :        
        """Draw/prepare NLL or LL-profile for seleted variable:        
        >>> model.fitTo ( dataset , ... )
        >>> nll  , f1 = model.draw_nll ( 'B' ,  dataset )
        >>> prof , f2 = model.draw_nll ( 'B' ,  dataset , profile = True )
        """
        return self.pdf.draw_nll ( var , dataset ,
                                   profile = profile ,
                                   draw    = draw    ,
                                   silent  = silent  ,
                                   args    = args    , **kwargs )  
        
    # ========================================================================
    ## evaluate "significance" using Wilks' theorem via NLL
    #  @code
    #  data = ...
    #  pdf  = ...
    #  pdf.fitTo ( data , ... )
    #  sigmas = pdf.wilks ( 'S' , data )
    #  @endcode
    def wilks ( self                     ,
                var                      ,
                dataset                  ,
                range    = ( 0 , None )  ,
                silent   = True          ,
                args     = () , **kwargs ) :
        """Evaluate 'significance' using Wilks' theorem via NLL
        >>> data = ...
        >>> pdf  = ...
        >>> pdf.fitTo ( data , ... )
        >>> sigmas = pdf.wilks ( 'S' , data )
        """
        return self.pdf.wilks ( var , dataset   ,
                                range  = range  ,
                                silent = silent ,
                                args   = args   , **kwargs )


    # ========================================================================
    ## evaluate "significance" using Wilks' theorem via NLL
    #  @code
    #  data = ...
    #  pdf  = ...
    #  pdf.fitTo ( data , ... )
    #  sigmas = pdf.wilks2 ( 'S' , data , fix = [ 'mean' , 'gamma' ] )
    #  @endcode
    def wilks2 ( self                           ,
                 var                            ,
                 dataset                        ,
                 fix                            , ## variables to fix 
                 range          = ( 0 , None )  ,
                 silent         = True          ,
                 ## opt_const      = True          ,
                 ## max_calls      = 10000         ,
                 ## max_iterations = -1            ,
                 ## strategy       = None          ,
                 args           = () , **kwargs ) :
        """Evaluate 'significance' using Wilks' theorem via NLL
        >>> data = ...
        >>> pdf  = ...
        >>> pdf.fitTo ( data , ... )
        >>> sigmas = pdf.wilks2 ( 'S' , data , fix = [ 'mean' , 'gamma'] )
        """
        return self.pdf.wilks2 ( var                              ,
                                 dataset        = dataset         ,
                                 fix            = fix             ,
                                 range          = range           , 
                                 silent         = silent          , 
                                 ## opt_const      = opt_const       ,
                                 ## max_calls      = max_calls       ,
                                 ## max_iterations = max_iterations  ,
                                 ## strategy       = strategy        ,
                                 args           = args , **kwargs )
    
    # ========================================================================
    ## get the actual minimizer for the explicit manipulations
    #  @code
    #  data = ...
    #  pdf  = ...
    #  m    = pdf.minuit  ( data )
    #  m.migrad()
    #  m.hesse ()
    #  m.minos ( param )
    #  @endcode
    #  @see RooMinimizer
    def minuit ( self , dataset        ,
                 max_calls      = -1   ,
                 max_iterations = -1   , 
                 opt_const      = True , ## optimize const 
                 strategy       = None ,
                 nLL            = None , 
                 args           =   () , **kwargs  ):
        """Get the actual minimizer for the explicit manipulations
        >>> data = ...
        >>> pdf  = ...
        >>> m    = pdf.minuit ( data )
        >>> m.migrad()
        >>> m.hesse ()
        >>> m.minos ( param )
        - see ROOT.RooMinimizer
        """
        return self.pdf.minuit ( dataset                         ,
                                 max_calls      = max_calls      ,
                                 max_iterations = max_iterations ,
                                 opt_const      = opt_const      ,
                                 strategy       = strategy       ,
                                 nLL            = nLL            ,
                                 args           = args           , **kwargs )

    
    # =========================================================================
    ## generate toy-sample according to PDF
    #  @code
    #  model  = ....
    #  data   = model.generate ( 10000 ) ## generate dataset with 10000 events
    #  varset = ....
    #  data   = model.generate ( 100000 , varset , sample = False )
    #  data   = model.generate ( 100000 , varset , sample = True  )     
    #  @endcode
    def generate ( self                  , 
                   nEvents               , 
                   varset        = None  ,
                   args          = ()    ,
                   binning       = {}    ,
                   sample        = True  ,
                   storage       = None  ,
                   category_args = {}    ) :
        """Generate toy-sample according to PDF
        >>> model  = ....
        >>> data   = model.generate ( 10000 ) ## generate dataset with 10000 events
        
        >>> varset = ....
        >>> data   = model.generate ( 100000 , varset , sample = False )
        >>> data   = model.generate ( 100000 , varset , sample = True  )
        """

        labels = self.sample.labels()
        
        assert len ( labels ) == len ( nEvents ), 'Invalid length of nEvents array'
        
        vars   = ROOT.RooArgSet()
        data   = {}

        weight = None
        wvar   = None
        ## generate all categories separately:        
        for l , n in zip ( labels , nEvents ) :
            
            cargs = []
            for a in args                         : cargs.append ( a )
            for a in category_args.get ( l , () ) : cargs.append ( a )
            cargs = tuple ( cargs )
                
            pdf   = self.categories [ l ]
            ds    = pdf.generate ( n                   ,
                                   varset   = varset   ,
                                   binning  = binning  ,
                                   sample   = sample   ,
                                   storage  = storage  , 
                                   args     = cargs    )

            if ds.isWeighted() :
                ds , weight  = ds.unWeighted ()
                if weight : wvar = getattr ( ds , weight )
                
            data [ l ]  = ds
            for v in ds.varset() :
                if not v in vars : vars.add ( v ) 
                ## vars       |= ds.varset()                 

        ## combine generated datasets
        args = () if not weight else ( ROOT.RooFit.WeightVar ( weight ) , )
        return combined_data ( self.sample ,
                               vars        , 
                               data        ,
                               args = args )
    
    # =========================================================================
    ## perform sPlot-analysis 
    #  @code
    #  r,f = model.fitTo ( dataset )
    #  model.sPlot ( dataset ) 
    #  @endcode 
    def sPlot ( self , dataset , silent = False , **kwargs ) : 
        """ Make sPlot analysis
        >>> r,f = model.fitTo ( dataset )
        >>> model.sPlot ( dataset ) 
        """
        return self.pdf.sPlot ( dataset , silent = silent , **kwargs )

    # =========================================================================
    ## Load parameters from external dictionary <code>{ name : value }</code>
    #  or sequence of <code>RooAbsReal</code> objects
    #  @code
    #  pdf     = ...
    #  dataset = ...
    #  params  = { 'A' : 10 , 'B' : ... }
    #  pdf.load_params ( dataset , params ) 
    #  params  = ( A , B , C , ... )
    #  pdf.load_params ( dataset , params )  
    #  @endcode 
    def load_params ( self , params = {}  , dataset = None , silent = False , **kwargs ) :
        """Load parameters from external dictionary <code>{ name : value }</code>
        #  or sequence of <code>RooAbsReal</code> objects
        >>> pdf      = ...
        >>> dataset = ... 
        >>> params = { 'A' : 10 , 'B' : ... }
        >>> pdf.load_params ( params , dataset ) 
        >>> params = ( A , B , C , ... )
        >>> pdf.load_params ( params , dataset )  
        """
        return self.pdf.load_params ( params = params , dataset = dataset , silent = silent , **kwargs )

    # =========================================================================
    ##  Does this function depend on this variable,
    #   @code
    #   fun = ...
    #   var = ...
    #   if var in fun :
    #      ... 
    #   @endcode  
    def __contains__ ( self , var ) : 
        """Does this function depend on this variable?
        >>> fun = ...
        >>> var = ...
        >>> if var in fun :
        ...      ... 
        """
        return var in self.pdf
    
    # =========================================================================
    ## Get the parameters
    #  @code
    #  fun = ...
    #  parameters = fun.params ( )
    #  @endcode
    #  Or
    #  @code  
    #  pdf       = ...
    #  dataset   = ...
    #  parameters = pdf.params ( dataset)
    #  @endcode
    #  @see RooAbsReal::getParameters
    def params ( self , dataset = None  ) :
        """Get the parameters
        >>>  fun = ...
        >>> parameters = fun.params ( )
        Or
        >>>  pdf       = ...
        >>> dataset   = ...
        >>> parameters = pdf.params ( dataset)
        - see RooAbsReal::getParameters
        """
        return self.pdf.params  ( dataset )
    
    # =========================================================================
    ## get all parameters/variables in form of dictionary
    #  @code
    #  pdf    = ...
    #  params = pdf.parameters ( dataset ) 
    #  @endcode
    def parameters ( self , dataset = None ) :
        """ Get all parameters/varibales in form of dictionary
        >>> pdf    = ...
        >>> params = pdf.parameters ( dataset ) 
        """
        return self.pdf.parameters ( dataset )
    
    # ========================================================================
    ## get the parameter value by name
    #  @code
    #  pdf = ...
    #  p   = pdf.parameter  ( 'A' )
    #  @endcode
    def parameter ( self , param , dataset = None ) :
        """Get the parameter value by name
        >>> pdf = ...
        >>> p   = pdf.parameter  ( 'A' )
        """
        return self.pdf.parameter ( param , dataset )

    # ==========================================================================
    ## get parameter by name 
    #  @code
    #  pdf = ...
    #  a   = pdf['A']
    #  @endcode
    def __getitem__ ( self , param ) :
        """Get parameter by name 
        >>> pdf = ...
        >>> a   = pdf['A']
        """
        return self.pdf.__getitem__ ( param )

    # =========================================================================
    ## get NLL/profile-graph for the variable, using the specified bscissas
    #  @code
    #  pdf   = ...
    #  graph = pdf.graph_nll ( 'S'               ,
    #                          range ( 0 , 100 ) ,
    #                          dataset           )
    #  @endcode
    def graph_nll ( self , *args , **kwargs ) :
        """Get NLL/profile-graph for the variable, using the specified abscissas
        >>> pdf   = ...
        >>> graph = pdf.graph_nll ( 'S'               ,
        ...                          range ( 0 , 100 ) ,
        ...                          dataset           )
        """
        return self.pdf.graph_nll ( *args , **kwargs )
    
    # =========================================================================
    ## get NLL-profile-graph for the variable, using the specified abscissas
    #  @code
    #  pdf   = ...
    #  graph = pdf.graph_profile ( 'S'                       ,
    #                              vrange ( 0 , 12.5 , 10  ) ,
    #                              dataset                   )
    #  @endcode
    def graph_profile ( self , *args , **kwargs ) :
        """Get profile-graph for the variable, using the specified abscissas
        >>> pdf   = ...
        >>> graph = pdf.graph_profile ( 'S'                     ,
        ...                             range ( 0 , 12.5 , 20 ) ,
        ...                             dataset                 )
        """
        return self.pdf.graph_profile ( *args , **kwargs )
         
    # =========================================================================
    ## access to all signal components
    @property
    def signals ( self ) :
        """'signals' : get all signal components"""
        return self.pdf.signals
    # =========================================================================
    ## access to all background  components
    @property
    def backgrounds ( self ) :
        """'backgrounds' : get all background components"""
        return self.pdf.backgrounds
    # =========================================================================
    ## access to all "other" componets
    @property
    def components  ( self ) :
        """'components' : get all 'other' components"""
        return self.pdf.components
    # =========================================================================
    ## access to all 'crossterms1' components
    @property
    def crossterms1 ( self ) :
        """'crossterms1' : get all 'crossterms1' components"""
        return self.pdf.crossterms1
    # =========================================================================
    ## access to all 'crossterms2' componets
    @property
    def crossterms2 ( self ) :
        """'crossterms2' : get all 'crossterms2' components"""
        return self.pdf.crossterms2
    # =========================================================================
    ## access to all 'combined-signal' components
    @property
    def combined_signals ( self ) :
        """'combined_signals' : get all combined-signal components"""
        return self.pdf.combined_signals
    # =========================================================================
    ## access to all 'combined-backgrounds' components
    @property
    def combined_backgrunds ( self ) :
        """'combined_backgrounds' : get all combined-background components"""
        return self.pdf.combined_background
    # =========================================================================
    ## access to all 'combined-backgrounds' components
    @property
    def combined_components ( self ) :
        """'combined_components' : get all combined  components"""
        return self.pdf.combined_components
    # =========================================================================
    ## access to 'alist1' structure (crucial for sPlot)
    @property
    def alist1 ( self ) :
        """'alist1' : get 'alist1' structure (crucial for sPlot)"""
        return self.pdf.alist1 
    # =========================================================================
    ## access to 'alist2' structure (crucial for sPlot)
    @property
    def alist2 ( self ) :
        """'alist2' : get 'alist2' structure (crucial for sPlot)"""
        return self.pdf.alist2
    # =========================================================================
    @property
    def fit_options ( self ) :
        """'fit_options' : the predefined 'fitTo'-options for this PDF
        - tuple of ROOT.RooArgCmd
        pdf = ...
        pdf.fit_options = ROOT.RooFit.Optimize ( 1 )
        pdf.fit_options = ROOT.RooFit.Optimize ( 1 ) , ROOT.RooFit.PrintEvalError ( 2 ) 
        """
        return self.pdf.fit_options 
    @fit_options.setter
    def fit_options ( self , value )  :
        self.pdf.fit_options = value 
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END 
# =============================================================================
