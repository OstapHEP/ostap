#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostap/fitting/simfit.py
#  Collection of utilities that simplify Simultaneous fit
#  @see RooSimultaneous
#  @see https://github.com/OstapHEP/ostap/blob/master/ostap/fitting/SIMFIT.md
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2018-11-23
# =============================================================================
""" Collection of utilities that simplify Simultaneous fit
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
import ROOT, math,  random , warnings 
from   ostap.core.core     import std , Ostap , dsID , items_loop 
from   ostap.fitting.utils import MakeVar
from   ostap.fitting.basic import PDF , Generic1D_pdf
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.simfit' )
else                       : logger = getLogger ( __name__               )
# =============================================================================
## redefine ROOT.RooCategory constructor to define the categories
#  @code
#  sample = ROOT.RooCategory('sample','fitting sample', 'signal' , 'control' )
#  @endcode
def _rc_init_ ( self , name , title , *categories ) :
    """Modified ROOT.RooCategory constructor to define the categories
    >>> sample = ROOT.RooCategory('sample','fitting sample', 'signal' , 'control' )
    """
    ROOT.RooCategory._old_init_ ( self , name ,  title )
    for c in categories : self.defineType ( c )

if not hasattr ( ROOT.RooCategory , '_old_init_' ) :
    ROOT.RooCategory._old_init_ = ROOT.RooCategory.__init__
    ROOT.RooCategory.__init__   = _rc_init_     
    
# =============================================================================
## Get the list/tuple of categories 
#  @code
#  cat = ....
#  labels = cat.labels()
#  @endcode
#  @see RooCategory
def _rc_labels_ ( self ) :
    """Get the list/tuple of categories
    >>> cat = ....
    >>> labels = cat.labels()
    """
    _iter = Ostap.Utils.Iterator ( self.typeIterator() )
    _icat = _iter.Next()

    labs = [] 
    while _icat  :

        labs.append ( _icat.GetName() )
        _icat = _iter.Next() 
        
    del _iter

    return tuple ( labs ) 

ROOT.RooCategory.labels = _rc_labels_

# =============================================================================
## Create combined dataset for simultaneous fit
#  @code
#  sample = ROOT.RooCategory ( 'sample' , 'sample' , 'cc' , 'zz' )
#  vars   = ROOT.RooArgSet   ( m2c )
#  ds_cmb = combined_data ( sample  ,
#                vars    , { 'cc' : ds_cc ,  'zz' : ds_00 } )
#  @endcode
#  - weighted variant:
#  @code
#  wvars = ROOT.RooArgSet ( m2c , SS_sw ) 
#  dsw_cmb   = combined_data ( sample ,
#                 wvars  , { 'cc' : dsn_cc ,  'zz' : dsn_00 } ,
#                 args = ( ROOT.RooFit.WeightVar( 'SS_sw' ) , ) )
#  @endcode
def combined_data ( sample        ,
                    varset        , 
                    datasets      ,
                    name     = '' ,
                    title    = '' ,
                    args     = () ) :
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
    
    labels = sample.labels()
    
    largs   = [ ROOT.RooFit.Index ( sample ) ] 

    for label in labels :

        dset = None 
        if isinstance ( datasets , dict ) : dset = datasets[label]
        else :
            for ds in dataset :
                if label == ds[0] :
                    dset =  ds[1]
                    break
                
        assert isinstance ( dset , ROOT.RooAbsData ),\
               'Invalid data set for label %s' % label
        
        largs.append (  ROOT.RooFit.Import ( label , dset ) )

    name  = name  if name  else dsID()
    title = title if title else 'Data for simultaneous fit/%s' % sample.GetName()

    args = args + tuple ( largs )

    vars = ROOT.RooArgSet()
    if   isinstance ( varset , ROOT.RooArgSet  ) : vars = varset
    elif isinstance ( varset , ROOT.RooAbsReal ) : vars.add ( varset )
    else :
        for v in varset : vars.add ( v )
        
    return ROOT.RooDataSet ( name , title , vars , *args )

# =============================================================================
## create combined binned dataset for simultaneous fit
def combined_hdata ( sample        ,
                     varset        ,
                     histograms    ,
                     name     = '' ,
                     title    = '' ) :

    MAP  = std.map  ( 'std::string'       , 'TH1*' )
    PAIR = std.pair ( 'const std::string' , 'TH1*' )
    mm   = MAP()
    for key in histograms : 
        mm.insert ( PAIR ( key , histograms [ key ] ) ) 
        
    name  = name  if name  else dsID()
    title = title if title else 'Data for simultaneous fit/%s' % sample.GetName()

    varlst = ROOT.RooArgList()
    if isinstance ( varset , ROOT.RooAbsReal ) : varlst.add ( varset )
    else :  
        for v in varset : varlst.add ( v )
    
    return ROOT.RooDataHist ( name , title , varlst , sample  , mm ) 
    
# =============================================================================        
## @class Sim1D
#  Helper class to simplify creation and usage of simultaneous PDF
#  @code
#  @endcode 
#  @see RooSimultaneous
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2018-11-23
class Sim1D(PDF) :
    """Helper class to simplify the creation and usage of simultaneous PDF
    
    - see RooSimultaneous 
    """
    
    def __init__ ( self              ,
                   sample            , 
                   categories        ,
                   xvar       = None ,
                   name       = None , 
                   title      = ''   ) :

        warnings.warn("Usage of obsolete Sim1D. Use SimFit instead")
        
        if isinstance ( sample , ( tuple , list ) ) :
            _cat = ROOT.RooCategory ( 'sample' , 'sample' )
            for i in sample : _cat.defineType ( i ) 
            sample =  _cat
            
        assert isinstance ( sample , ROOT.RooCategory ),\
               'Invalid type for "sample":' % ( sample ,  type ( sample ) )
        
        if not name  : name  = 'SimFit_'                 +          sample.GetName()
        if not title : title = 'Simultaneous PDF(%s,%s)' % ( name , sample.GetName() )
        
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
                    
            if   isinstance ( cmp , PDF ) :
                
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
        PDF.__init__ ( self , name , xvar = _xvar ) 
        
        self.pdf = ROOT.RooSimultaneous ( 'sim_' + self.name , 
                                          title              , 
                                          self.sample        )

        keys = self.categories.keys()
        for key in sorted ( keys ) :
            self.pdf.addPdf ( self.categories[key].pdf , key )

        for k , pdf in items_loop ( self.categories ) :
            
            for c in pdf.signals     : self.signals    .add ( c ) 
            for c in pdf.backgrounds : self.backgrounds.add ( c ) 
            for c in pdf.crossterms1 : self.crossterms1.add ( c ) 
            for c in pdf.crossterms2 : self.crossterms2.add ( c ) 

            ## for c in pdf.alist1      : self.alist1.add ( c ) 
            ## for c in pdf.alist2      : self.alist2.add ( c ) 

        self.config = {
            'name'       : self.name       ,
            'title'      : title           ,            
            'sample'     : self.sample     ,
            'categories' : self.categories ,
            'xvar'       : self.xvar       ,
            }

    
    @property
    def sample  ( self ) :
        "``sample'' : RooCategory object for simultaneous PDF"
        return self.__sample

    @property
    def samples ( self ) :
        "``samples'' : list/tuple of known categories"
        return tuple ( self.__categories.keys() ) 

    @property
    def categories ( self ) :
        "``categories'' : map { category : pdf } "
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

        res , frame = PDF.fitTo ( self ,
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
            
            return PDF.draw ( self ,
                              dataset = dataset ,
                              nbins   = nbins   ,
                              silent  = silent  , **kwargs )
        
        
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
#  Helper class to simplify creation and usage of simultaneous PDF
#  @code
#  @endcode 
#  @see RooSimultaneous
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2018-11-23
class SimFit ( MakeVar ) :
    """Helper class to simplify the creation and usage of simultaneous PDF
    
    - see RooSimultaneous 
    """
    
    def __init__ ( self              ,
                   sample            , 
                   categories        ,
                   name       = None , 
                   title      = ''   ) :
        
        if isinstance ( sample , ( tuple , list ) ) :
            _cat = ROOT.RooCategory ( 'sample' , 'sample' )
            for i in sample : _cat.defineType ( i ) 
            sample =  _cat
            
        assert isinstance ( sample , ROOT.RooCategory ),\
               'Invalid type for "sample":' % ( sample ,  type ( sample ) )
        
        if not name  : name  = 'SimFit_'                 +          sample.GetName()
        if not title : title = 'Simultaneous PDF(%s,%s)' % ( name , sample.GetName() )

        ## propagatethe name 
        self.name = name
        
        self.__sample       = sample 
        self.__categories   = {}

        # =====================================================================
        ## components
        # =====================================================================
        labels = sample.labels()

        from ostap.fitting.basic import PDF 
        from ostap.fitting.fit2d import PDF2
        from ostap.fitting.fit3d import PDF3

        _xv = None 
        for label in labels :
            
            cmp   = None 
            if isinstance ( categories , dict ) : cmp = categories [ label ]
            else :
                for ii in categories :
                    if ii[0] == label :
                        cmp = ii[1]
                        break

            if not isinstance  ( cmp , PDF  ) :
                raise TypeError ( 'Can not find the proper category component: "%s"' % label ) 
            
            self.__categories [ label ] = cmp
            _xv = cmp.xvar
            
        sim_pdf     = PDF ( self.name , xvar = _xv )
        sim_pdf.pdf = ROOT.RooSimultaneous ( 'Sim_' + self.name , title , self.sample )
        
        keys = self.categories.keys()
        for key in sorted ( keys ) :
            sim_pdf.pdf.addPdf ( self.categories[key].pdf , key )

        self.__pdf = sim_pdf 
        
        for k , cmp in items_loop ( self.categories ) :
            
            for c in cmp.signals      : self.pdf.signals    .add ( c ) 
            for c in cmp.backgrounds  : self.pdf.backgrounds.add ( c ) 
            for c in cmp.crossterms1  : self.pdf.crossterms1.add ( c ) 
            for c in cmp.crossterms2  : self.pdf.crossterms2.add ( c )
            
            self.pdf.draw_options.update ( cmp.draw_options )
            
        # =====================================================================
        ##  drawing helpers
        # =====================================================================
        
        self.__drawpdfs   = {}
        for key in sorted ( keys ) :

            cmp = self.categories [ key ] 
            if isinstance  ( cmp , PDF3 ) :
                from ostap.fitting.fit3d import Generic3D_pdf                
                dpdf = Generic3D_pdf ( sim_pdf.pdf ,
                                       cmp.xvar    ,
                                       cmp.yvar    ,
                                       cmp.zvar    ,
                                       add_to_signals = False )
            elif isinstance  ( cmp , PDF2 ) :
                from ostap.fitting.fit2d import Generic2D_pdf                                
                dpdf = Generic2D_pdf ( sim_pdf.pdf ,
                                       cmp.xvar    ,
                                       cmp.yvar    ,
                                       add_to_signals = False )
            elif isinstance  ( cmp , PDF  ) :
                from ostap.fitting.basic import Generic1D_pdf   
                dpdf = Generic1D_pdf ( sim_pdf.pdf ,
                                       cmp.xvar    ,
                                       add_to_signals = False )
                
            for c in cmp.signals     : dpdf.signals    .add ( c ) 
            for c in cmp.backgrounds : dpdf.backgrounds.add ( c ) 
            for c in cmp.crossterms1 : dpdf.crossterms1.add ( c ) 
            for c in cmp.crossterms2 : dpdf.crossterms2.add ( c )


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
        "``sample'' : RooCategory object for simultaneous PDF"
        return self.__sample

    @property
    def pdf     ( self ) :
        "``pdf''  : the actual PDF with RooSimultaneous "
        return self.__pdf
    
    @property
    def samples ( self ) :
        "``samples'' : list/tuple of known categories"
        return tuple ( self.__categories.keys() ) 

    @property
    def categories ( self ) :
        "``categories'' : map { category : pdf } "
        return self.__categories

    @property
    def drawpdfs ( self ) :
        "``drawpdfs'' dictionary with PDFs for drawing"
        return self.__drawpdfs

    @property
    def draw_options ( self ) :
        """``draw_options'' : disctionary with predefined draw-options for this PDF
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
        >>> pdf.draw_options['signal_style'] = [ options ]
        - and later:
        >>> options = pdf.draw_option ( 'signal_style' )
        """
        return self.pdf.draw_option ( key , default , **kwargs ) 

    # =========================================================================
    ## delegate  attribute search to the components 
    def __getattr__ ( self , attr ) :
        """Delegate attribute search to the category components
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
               **kwargs       ) :
        """
        Draw the PDF&data for the given   category
        >>> pdf.fitTo ( dataset )
        >>> pf.draw ( 'signal' , dataset , nbins = 100 ) 
        """
        dvar = None
        if   isinstance ( category , ( tuple , list ) ) and 2 == len ( category ) :
            category, dvar = category 
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

        cat_pdf  = self.categories [ category ]
        draw_pdf = self.drawpdfs   [ category ]

        if 1 < 2 : 
        ## with KeepArgs     ( draw_pdf . signals     , cat_pdf . signals     ) as _k1 ,\
        ##          KeepArgs ( draw_pdf . backgrounds , cat_pdf . backgrounds ) as _k2 ,\
        ##          KeepArgs ( draw_pdf . components  , cat_pdf . components  ) as _k3 ,\
        ##          KeepArgs ( draw_pdf . crossterms1 , cat_pdf . crossterms1 ) as _k4 ,\
        ##          KeepArgs ( draw_pdf . crossterms2 , cat_pdf . crossterms2 ) as _k5 :

            from ostap.fitting.basic import PDF 
            from ostap.fitting.fit2d import PDF2
            from ostap.fitting.fit3d import PDF3

            if   isinstance ( draw_pdf , PDF3 ) :

                if   3 == dvar or dvar in  ( 'z' , 'Z' , '3' , draw_pdf.zvar.name ) : 
                    return draw_pdf.draw3 ( dataset = dataset ,
                                            nbins   = nbins   ,
                                            silent  = silent  , **kwargs )
                elif 2 == dvar or dvar in  ( 'y' , 'Y' , '2' , draw_pdf.yvar.name ) : 
                    return draw_pdf.draw2 ( dataset = dataset ,
                                            nbins   = nbins   ,
                                            silent  = silent  , **kwargs )
                elif 1 == dvar or dvar in  ( 'x' , 'X' , '1' , draw_pdf.xvar.name ) : 
                    return draw_pdf.draw1 ( dataset = dataset ,
                                            nbins   = nbins   ,
                                            silent  = silent  , **kwargs )
                else :
                    self.error('Unknown ``dvar'' for 3D-draw pdf!')
                    return None
                
            elif isinstance ( draw_pdf , PDF2 ) :

                if   2 == dvar or dvar in  ( 'y' , 'Y' , '2' , draw_pdf.yvar.name ) :
                    print 'here2'
                    return draw_pdf.draw2 ( dataset = dataset ,
                                            nbins   = nbins   ,
                                            silent  = silent  , **kwargs )
                elif 1 == dvar or dvar in  ( 'x' , 'X' , '1' , draw_pdf.xvar.name ) : 
                    print 'here1'
                    return draw_pdf.draw1 ( dataset = dataset ,
                                            nbins   = nbins   ,
                                            silent  = silent  , **kwargs )
                else :
                    self.error('Unknown ``dvar'' for 2D-draw pdf! %s' %  dvar )
                    return None 

            elif isinstance ( draw_pdf , PDF  ) :
                
                return draw_pdf.draw ( dataset = dataset ,
                                       nbins   = nbins   ,
                                       silent  = silent  , **kwargs )
            
        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
# The END 
# =============================================================================
