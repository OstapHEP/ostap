#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file basic.py
#  Set of useful basic utilities to build various fit models 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
# =============================================================================
"""Set of useful basic utilities to build various 2D-fit models"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'PDF2'          , ## useful base class for 2D-models
    'Fit2D'         , ## the model for 2D-fit: signal + background + optional components
    'Fit2DSym'      , ## the model for 2D-fit: signal + background + optional components
    ##
    'H2D_dset'      , ## convertor of 2D-histo to RooDataHist 
    'H2D_pdf'       , ## convertor of 1D-histo to RooDataPdf
    ##
    'Generic2D_pdf' , ## wrapper over imported RooFit (2D)-pdf  
    )
# =============================================================================
import ROOT
from   ostap.fitting.basic import PDF, makeVar, makeBkg, H2D_dset
from   ostap.logger.utils  import roo_silent 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.fit2d' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
# @class PDF2
#  The helper base class for implementation of 2D-pdfs 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-08-21
class PDF2 (PDF) :
    """ Useful helper base class for implementation of PDFs for 2D-fit
    """
    def __init__ ( self , name , xvar = None , yvar = None ) : 
        
        PDF.__init__ ( self , name , xvar )
        
        self.__yvar = None
        
        ## create the variable 
        if isinstance ( yvar , tuple ) and 2 == len(yvar) :  
            self.__yvar = makeVar ( yvar         , ## var 
                                    'y'          , ## name 
                                    'y-variable' , ## title/comment
                                    *yvar        , ## min/max 
                                    fix = None   ) ## fix ? 
        elif isinstance ( yvar , ROOT.RooAbsReal ) :
            self.__yvar = makeVar ( yvar         , ## var 
                                    'y'          , ## name 
                                    'y-variable' , ## title/comment
                                    fix = None   ) ## fix ? 
        else :
            logger.warning('PDF2: ``y-variable''is not specified properly %s/%s' % ( yvar , type ( yvar ) ) )
            self.__yvar = makeVar( yvar , 'y' , 'y-variable' )

        ## save the configuration
        self.config = {
            'name' : self.name ,
            'xvar' : self.xvar ,
            'yvar' : self.yvar ,            
            }
        
    def yminmax ( self ) :
        """Min/max values for y-varibale"""
        return self.__yvar.minmax()
    
    @property 
    def yvar ( self ) :
        """``y''-variable for the fit (same as ``y'')"""
        return self.__yvar

    @property 
    def y    ( self ) :
        """``y''-variable for the fit (same as ``yvar'')"""
        return self.__yvar

    # =========================================================================
    ## make the actual fit (and optionally draw it!)
    #  @code
    #  r,f = model.fitTo ( dataset )
    #  r,f = model.fitTo ( dataset , weighted = True )    
    #  r,f = model.fitTo ( dataset , ncpu     = 10   )    
    #  r,f = model.fitTo ( dataset , draw = True , nbins = 300 )    
    #  @endcode 
    def fitTo ( self           , 
                dataset        ,
                draw   = False ,
                nbins  =    50 ,
                ybins  =  None , 
                silent = False ,
                refit  = False , *args , **kwargs ) :
        """
        Perform the actual fit (and draw it)
        >>> r,f = model.fitTo ( dataset )
        >>> r,f = model.fitTo ( dataset , weighted = True )    
        >>> r,f = model.fitTo ( dataset , ncpu     = 10   )    
        >>> r,f = model.fitTo ( dataset , draw = True , nbins = 300 )    
        """
        if isinstance ( dataset , ROOT.TH2 ) :
            density = kwargs.pop ( 'density' , True  ) 
            chi2    = kwargs.pop ( 'chi2'    , False ) 
            return self.fitHisto ( dataset   , draw , silent , density , chi2 , *args , **kwargs )

        
        result,f = PDF.fitTo ( self    ,
                               dataset ,
                               False   , ## false here!
                               nbins   ,
                               silent  ,
                               refit   , *args , **kwargs ) 
        if not draw :
            return result , None
        
        ## 2D 
        if 1< nbins and isinstance ( ybins , ( int , long ) ) and 1<ybins :
            return result, self.draw ( None , dataset , nbins , ybins , silent = silent )
        
        if     1<= nbins : return result, self.draw1 ( dataset ,  nbins , silent = silent )
        elif  -1>= nbins : return result, self.draw2 ( dataset , -nbins , silent = silent )

        ## return 2D 
        return result, self.draw ( None , dataset , silent = silent )
    
    # =========================================================================
    ## draw the projection over 1st variable
    #
    #  @code
    #  r,f = model.fitTo ( dataset ) ## fit dataset
    #  fx  = model.draw1 ( dataset , nbins = 100 ) ## draw results
    #
    #  f1  = model.draw1 ( dataset , nbins = 100 , in_range = (2,3) ) ## draw results
    #
    #  model.yvar.setRange ( 'QUQU2' , 2 , 3 ) 
    #  f1  = model.draw1 ( dataset , nbins = 100 , in_range = 'QUQU2') ## draw results
    #
    #  @endcode 
    def draw1 ( self            ,
                dataset  = None ,
                nbins    = 100  ,
                silent   = True ,
                in_range = None , *args , **kwargs ) :
        """ Draw the projection over 1st variable
        
        >>> r,f = model.fitTo ( dataset ) ## fit dataset
        >>> fx  = model.draw1 ( dataset , nbins = 100 ) ## draw results
        
        >>> f1  = model.draw1 ( dataset , nbins = 100 , in_range = (2,3) ) ## draw results

        >>> model.yvar.setRange ( 'QUQU2' , 2 , 3 ) 
        >>> f1  = model.draw1 ( dataset , nbins = 100 , in_range = 'QUQU2') ## draw results
        
        """
        if in_range and isinstance ( in_range , tuple ) and 2 == len ( in_range ) :
            self.yvar.setRange ( 'aux_rng2' , in_range[0] , in_range[1] )
            in_range = 'aux_rng2'
            
        return self.draw ( self.xvar , 
                           dataset   ,
                           nbins     ,
                           20        , ## fake 
                           silent    ,
                           in_range  , *args , **kwargs )
    
    # =========================================================================
    ## draw the projection over 2nd variable
    #
    #  @code
    #  r,f = model.fitTo ( dataset ) ## fit dataset
    #  fy  = model.draw2 ( dataset , nbins = 100 ) ## draw results
    #
    #  f2  = model.draw2 ( dataset , nbins = 100 , in_range = (2,3) ) ## draw results
    #
    #  model.xvar.setRange ( 'QUQU1' , 2 , 3 ) 
    #  f2  = model.draw2 ( dataset , nbins = 100 , in_range = 'QUQU1') ## draw results
    #
    #  @endcode 
    def draw2 ( self            ,
                dataset  = None ,
                nbins    = 100  ,
                silent   = True ,
                in_range = None , *args , **kwargs ) :
        """
        Draw the projection over 2nd variable
        
        >>> r,f = model.fitTo ( dataset ) ## fit dataset
        >>> fy  = model.draw2 ( dataset , nbins = 100 ) ## draw results
        
        >>> f2  = model.draw2 ( dataset , nbins = 100 , in_range = (2,3) ) ## draw results

        >>> model.xvar.setRange ( 'QUQU1' , 2 , 3 ) 
        >>> f2  = model.draw2 ( dataset , nbins = 100 , in_range = 'QUQU1') ## draw results

        """
        if in_range and isinstance ( in_range , tuple ) and 2 == len ( in_range ) :
            self.xvar.setRange ( 'aux_rng1' , in_range[0] , in_range[1] )
            in_range = 'aux_rng1'

        return self.draw ( self.yvar ,
                           dataset   ,
                           nbins     ,
                           20        , ## fake
                           silent    , in_range , *args , **kwargs )

    # =========================================================================
    ## draw as 2D-histograms 
    def draw_H2D ( self           ,
                   dataset = None ,  
                   xbins   = 20   ,
                   ybins   = 20   ) :
        """
        Make/draw 2D-histograms 
        """
        
        _xbins = ROOT.RooFit.Binning ( xbins ) 
        _ybins = ROOT.RooFit.Binning ( ybins ) 
        _yvar  = ROOT.RooFit.YVar    ( self.yvar , _ybins )
        _clst  = ROOT.RooLinkedList  ()
        hdata  = self.pdf.createHistogram ( hID() , self.xvar , _xbins , _yvar )
        hpdf   = self.pdf.createHistogram ( hID() , self.xvar , _xbins , _yvar )
        hdata.SetTitle(';;;')
        hpdf .SetTitle(';;;')
        _lst   = ROOT.RooArgList ( self.xvar , self.yvar )  
        if dataset : dataset.fillHistogram( hdata , _lst ) 
        self.pdf.fillHistogram  ( hpdf , _lst )
        
        if not ROOT.gROOT.IsBatch() :
            from Ostap.Utils import  rootWarning
            with rootWarning ():
                hdata.lego ()
                hpdf .Draw ( 'same surf')
        
        return hpdf , hdata 
    
    # =========================================================================
    ## make 1D-plot
    def draw ( self                         ,
               drawvar               = None ,
               dataset               = None ,
               nbins                 =  100 ,
               ybins                 =   20 ,
               silent                = True ,
               in_range              = None ,
               **kwargs                     ) : 
        """
        Make 1D-plot:
        """
        
        #
        ## special case:  do we need it? 
        # 
        if not drawvar : return self.draw_H2D( dataset , nbins , ybins )

        ## copy arguments:
        args = kwargs.copy ()
        
        import Ostap.FitDraw as FD
        if not isinstance ( in_range , (list,tuple) ) : in_range = in_range ,  
        if in_range :
            data_options        = args.pop (       'data_options' , FD.         data_options )
            background_options  = args.pop ( 'background_options' , FD. background2D_options )
            signal_options      = args.pop (     'signal_options' , FD.       signal_options )
            component_options   = args.pop (  'component_options' , FD.    component_options )
            crossterm1_options  = args.pop ( 'crossterm1_options' , FD.   crossterm1_options )
            crossterm2_options  = args.pop ( 'crossterm2_options' , FD.   crossterm2_options )
            total_fit_options   = args.pop (  'total_fit_options' , FD.    total_fit_options )

            for i in in_range :  
                data_options       += ROOT.RooFit.CutRange        ( i ) , 
                signal_options     += ROOT.RooFit.ProjectionRange ( i ) , 
                background_options += ROOT.RooFit.ProjectionRange ( i ) , 
                component_options  += ROOT.RooFit.ProjectionRange ( i ) , 
                crossterm1_options += ROOT.RooFit.ProjectionRange ( i ) , 
                crossterm2_options += ROOT.RooFit.ProjectionRange ( i ) , 
                total_fit_options  += ROOT.RooFit.ProjectionRange ( i ) , 
            
            args [       'data_options' ] =       data_options
            args [     'signal_options' ] =     signal_options
            args [ 'background_options' ] = background_options
            args [  'component_options' ] =  component_options
            args [ 'crossterm1_options' ] = crossterm1_options
            args [ 'crossterm2_options' ] = crossterm2_options
            args [  'total_fit_options' ] =  total_fit_options
            
        background_options    = args.pop ( 'background_options'    , FD.background2D_options    )
        base_background_color = args.pop ( 'base_background_color' , FD.base_background2D_color )
        args [ 'background_options'    ] = background_options
        args [ 'base_background_color' ] = base_background_color
        
        
        #
        ## redefine the drawing variable:
        # 
        self.draw_var = drawvar
        
        #
        ## delegate the actual drawing to the base class
        # 
        return PDF.draw ( self    ,
                          dataset ,
                          nbins   ,
                          silent  ,  **args ) 
    
    # =========================================================================
    ## fit the 2D-histogram (and draw it)
    #
    #  @code
    #
    #  histo = ...
    #  r,f = model.fitHisto ( histo )
    #
    #  @endcode
    def fitHisto ( self            ,
                   histo           ,
                   draw    = False ,
                   silent  = False ,
                   density = True  ,
                   chi2    = False , *args , **kwargs ) :
        """Fit the histogram (and draw it)
        
        >>> histo = ...
        >>> r,f = model.fitHisto ( histo , draw = True )
        
        """

        xminmax = histo.xminmax()
        yminmax = histo.yminmax()        
        with RangeVar( self.xvar , *xminmax ) , RangeVar ( self.yvar , *yminmax ): 
            ## convert it!
            self.histo_data = H2D_dset ( histo , self.xvar , self.yvar  , density , silent )
            data = self.histo_data.dset 
            
            ## fit it!!
            if chi2 : return self.chi2fitTo ( data              ,
                                              draw    = draw    ,
                                              silent  = False   ,
                                              density = density , *args , **kwargs )
            else     : return self.fitTo    ( data              ,
                                              draw    = draw    ,
                                              nbins   = histo.nbinsx() ,
                                              ybins   = histo.nbinsy() ,
                                              silent  = silent  , *args , **kwargs )            
    ## adjust PDF a little bit to avoid zeroes
    #  A tiny  ``flat'' component is added and the orginal PDF is replaced by a new compound PDF.
    #  The fraction of added  component is fixed and defined by ``value''
    #  @code
    #  >>> pdf = ...
    #  >>> pdf.adjust ( 1.e-6 )
    #  @endcode
    #  The  fraction can be changed and/or relesed:
    #  @code
    #  >>> pdf.adjustment.fraction = 1.e-4    ## release it
    #  >>> pdf.adjustment.fraction.release()  ## allow to  vary in the fit 
    #  @endcode 
    #  The original PDF is stored as:
    #  @code
    #  >>> orig_pdf = pdf.adjustment.old_pdf 
    #  @endcode
    def adjust ( self , value =  1.e-5 ) :
        """``adjust'' PDF a little bit to avoid zeroes
        A tiny  ``flat'' component is added and the orginal PDF is replaced by a new compound PDF.
        The fraction of added  component is fixed and defined by ``value''

        >>> pdf = ...
        >>> pdf.adjust ( 1.e-6 )

        The  fraction can be changed/relesed

        >>> pdf.adjustment.fraction = 1.e-4    ## change the value 
        >>> pdf.adjustment.fraction.release()  ## release it, allow to vary in the fit 
        
        The original PDF is stored as:
        
        >>> orig_pdf = pdf.adjustment.old_pdf 
        
        """
        if self.adjustment :
            logger.warning ( "PDF is already adjusted, skip it!")
            return

        ## create adjustment object and  use it to adjust PDF:
        self.__adjustment = Adjust2D ( self.name , self.xvar , self.yvar , self.pdf , value )
        ## replace original PDF  with  adjusted one:
        self.pdf          = self.__adjustment.pdf

    # =========================================================================
    ## simple 'function-like' interface 
    def __call__ ( self , x , y ) :
        """ Simple  function-like interface
        >>>  pdf = ...
        >>>  print pdf(0.1,0.5) 
        """
        if     isinstance ( self.xvar , ROOT.RooRealVar ) and \
               isinstance ( self.yvar , ROOT.RooRealVar ) :
            
            from Ostap.RooFitDeco import SETVAR
            if x in self.xvar and y in self.yvar : 
                with SETVAR ( self.xvar ) , SETVAR( self.yvar ) :
                    self.xvar.setVal ( x )
                    self.yvar.setVal ( y )
                    return self.pdf.getVal()
            else :
                return 0.0
            
        raise AttributeError, 'something wrong goes here'

    # =========================================================================
    ## get integral over (xmin,xmax,ymin,ymax) region
    #  @code
    #  pdf = ...
    #  print pdf.integral( 0,1,0,2)
    #  @endcode
    def integral ( self, xmin , xmax , ymin , ymax ) :
        """Get integral over (xmin,xmax,ymin,ymax) region
        >>> pdf = ...
        >>> print pdf.integral( 0,1,0,2)
        """
        xmn , xmx = self.xminmax()
        ymn , ymx = self.yminmax()

        xmin = max ( xmin , xmn )
        xmax = min ( xmax , xmx )
        ymin = max ( ymin , ymn )
        ymax = min ( ymax , ymx )

        ## make a try to use analytical integral (could be fast)
        if hasattr ( self , 'pdf' ) :
            _pdf = self.pdf 
            if hasattr ( _pdf , 'setPars'  ) : _pdf.setPars() 
            try: 
                if hasattr ( _pdf , 'function' ) :
                    _func = _pdf.function() 
                    if hasattr ( _func , 'integral' ) :
                        return _func.integral ( xmin , xmax , ymin , ymax )
            except:
                pass
            
        ## use numerical integration 
        from scipy import integrate 
        result = integrate.dblquad ( self ,
                                     ymin ,
                                     ymax ,
                                     lambda x : xmin ,
                                     lambda x : xmax , 
                                     *args , **kwargs )
        return result[0]

# =============================================================================
## suppress methods specific for 1D-PDFs only
for _a in (
    ##'_get_stat_'     ,
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

    if hasattr ( PDF2 , _a ) :
        def _suppress_ ( self , *args , **kwargs ) :
            raise AttributeError ( "'%s' object has no attribute '%s'" % ( type(self) , _a ) )
        setattr ( PDF2 , _a , _suppress_ ) 
        logger.verbose ( 'Remove attribute %s from PDF2' ) 


# =============================================================================
## @class Flat2D
#  The most trivial 2D-model - constant
#  @code 
#  pdf = Flat2D( 'flat' , xvar = ...  , yvar = ... )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
class Flat2D(PDF2) :
    """The most trival 2D-model - constant
    >>> pdf = Flat2D( 'flat' , xvar = ...  , yvar = ... )
    """
    def __init__ ( self , xvar , yvar , name = 'Flat2D') :

        PDF2.__init__ ( self  , name , xvar , yvar ) 
                        
        self.__xp0 = ROOT.RooPolynomial( 'xp0_%s'   % name , 'xpoly0(%s)'   % name , xvar )        
        self.__yp0 = ROOT.RooPolynomial( 'yp0_%s'   % name , 'ypoly0(%s)'   % name , yvar )
        
        self.pdf   = ROOT.RooProdPDF ( name , 'poly0_2D(%s)' % name , self.__xp0 , self.__yp0  )
        ## save configuration
        self.config = {
            'name'     : self.name ,            
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            }                   

# =============================================================================
## simple class to adjust certaint PDF to avoid zeroes 
class Adjust2D(object) :
    """Simple class to ``adjust'' certain PDF to avoid zeroes
    - a small flat component is added and the compound PDF is constructed
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,
                   xvar             , 
                   yvar             , 
                   pdf              ,
                   value    = 1.e-5 ) : 
        
        assert isinstance ( pdf  , ROOT.RooAbsPdf  ) , "``pdf''  must be ROOT.RooAbsPdf"
        assert isinstance ( xvar , ROOT.RooAbsReal ) , "``xvar'' must be ROOT.RooAbsReal"
        assert isinstance ( yvar , ROOT.RooAbsReal ) , "``yvar'' must be ROOT.RooAbsReal"
        
        self.name      = name 
        self.__old_pdf = pdf
        
        self.__flat    = Flat2D  ( xvar , yvar , name = 'flat_' + name )
        self.__frac    = makeVar ( value , 'fracA_%s'                     % name ,
                                   'small  fraction of flat component %s' % name ,
                                   value , 1.e-4 , 0 , 1 )
        
        self.__alist1  = ROOT.RooArgList ( self.__flat.pdf , self.__old_pdf )        
        self.__alist2  = ROOT.RooArgList ( self.__frac     )
        #
        ## final PDF
        # 
        self.__pdf     = ROOT.RooAddPdf  ( "adjust_"    + name ,
                                           "Adjust(%s)" % name ,
                                           self.__alist1 ,
                                           self.__alist2 )        
    @property
    def fraction( self ) :
        """``fraction''-parameter: the fraction of flat background added"""
        return  self.__frac
    @fraction.setter 
    def fraction( self , value ) :
        value = float ( value )
        assert 0 < value < 1 , 'Fraction  must be between 0 and 1'
        self.__frac.setVal ( value )
        
    @property
    def flat ( self ) :
        """new artificial ``flat'' component for the PDF"""
        return self.__flat

    @property
    def pdf ( self ) :
        """``new'' (adjusted) PDF"""
        return self.__pdf
    
    @property
    def old_pdf ( self ) :
        """``old'' (non-adjusted) PDF"""
        return self.__old_pdf
               
# =============================================================================
## @class Model2D
#  Trivial class to construct 2D model as a product of split 1D-models
#  actually it is a tiny  wrapper over <code>ROOT.RooProdPdf</code>
#  @code
#  pdfx = ...
#  pdfy = ...
#  pdf2D = Model2D( 'D2' , xmodel = pdfx , ymodel =  pdfy )
#  @endcode 
class Model2D(PDF2) :
    """Trivial class to construct 2D model as a product of split 1D-models
    - actually it is a tiny  wrapper over ROOT.RooProdPdf
    >>> pdfx = ...
    >>> pdfy = ...
    >>> pdf2D = Model2D( 'D2' , xmodel = pdfx , ymodel =  pdfy )
    """
    def __init__ ( self         ,
                   name         ,
                   xmodel       ,
                   ymodel       ,
                   xvar  = None ,
                   yvar  = None ,
                   title = ''   ) :

        if   isinstance ( xmodel , PDF            ) : self.__xmodel = xmodel
        elif isinstance ( xmodel , ROOT.RooAbsPdf ) and xvar :
            self.__xmodel = Generic1D_pdf  ( xmodel , xvar )
        else : raise AttributeError ( "Invalid ``x-model'' attribute" )

        if   isinstance ( ymodel , PDF            ) : self.__ymodel = ymodel
        elif isinstance ( ymodel , ROOT.RooAbsPdf ) and xvar :
            self.__ymodel = Generic1D_pdf  ( ymodel , yvar )
        else : raise AttributeError ( "Invalid ``y-model'' attribute" )

        ## initialize the base 
        PDF2.__init__ (  self , name , self.__xmodel.xvar , self.__ymodel.xvar ) 

        ## check the title 
        if not title : title = '%s x %s' % ( self.__xmodel.name , self.__ymodel.name )
        
        ## build the final PDF 
        self.pdf = ROOT.RooProdPdf (
            name  ,
            title ,
            self.__xmodel.pdf ,
            self.__ymodel.pdf )

        ## save configuration 
        self.config = {
            'name'   :  self.name   ,
            'xmodel' :  self.xmodel ,
            'xmodel' :  self.ymodel ,
            'xvar'   :  self.xvar   ,
            'yvar'   :  self.yvar   ,            
            'title'  :  self.pdf.GetTitle() 
            }

    @property
    def xmodel ( self ) :
        """``x-model'' x-component of Model(x)*Model(y) PDF"""
        return self.__xmodel

    @property
    def ymodel ( self ) :
        """``y-model'' y-component of Model(x)*Model(y) PDF"""
        return self.__ymodel

# =============================================================================
## @class Fit2D
#  The actual model for 2D-fits
#
#  @code
# 
#  model   = Models.Fit2D (
#      signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
#      signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y ) ,
#      bkg1     = 1 , 
#      bkg2     = 1 )
#
#  r,f = model.fitTo ( dataset ) ## fit dataset 
#
#  print r                       ## get results  
#
#  fx  = model.draw1 ()          ## visualize X-projection
#  fy  = model.draw2 ()          ## visualize Y-projection
#
#  @endcode 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class Fit2D (PDF2) :
    """The actual model for 2D-fits
    
    >>>  model   = Models.Fit2D (
    ...      signal_1 = Models.Gauss_pdf ( 'Gx' , mass = m_x ) ,
    ...      signal_2 = Models.Gauss_pdf ( 'Gy' , mass = m_y ) ,
    ...      bkg1     = 1 , 
    ...      bkg2     = 1 )
    >>> r,f = model.fitTo ( dataset ) ## fit dataset 
    >>> print r                       ## get results  
    >>> fx  = model.draw1 ()          ## visualize X-projection
    >>> fy  = model.draw2 ()          ## visualize Y-projection

    """
    def __init__ ( self               ,
                   #
                   signal_1           , 
                   signal_2           ,
                   suffix = ''        ,
                   #
                   bkg1       = None  ,
                   bkg2       = None  ,
                   #
                   bkgA       = None  ,
                   bkgB       = None  ,
                   #
                   bkg2D      = None  ,
                   #
                   ## main components :
                   ss         = None  , ## signal    (1) * signal     (2)
                   sb         = None  , ## signal    (1) * background (2) 
                   bs         = None  , ## background(1) * signal     (2)
                   bb         = None  , ## background-2D 
                   ## additional components 
                   components = []    ,
                   xvar       = None  ,
                   yvar       = None  ,                   
                   name       = ''    ) :
        
        ## collect all the arguments 
        self.__args = {
            'signal_1'   : signal_1 , 'signal_2' : signal_2 ,
            'bkg1'       : bkg1     , 'bkg2'     : bkg2     ,
            'bkgA'       : bkgA     , 'bkgB'     : bkgB     ,
            'bkg2D'      : bkg2D ,
            'components' : bkg2D ,
            ##
            'ss'         : ss , 'bb'         : bb ,
            'sb'         : sb ,'bs'         : bs ,
            ##
            'suffix'     : suffix   ,
            'name'       : name     ,
            }
        
        self.__crossterms1 = ROOT.RooArgSet()
        self.__crossterms2 = ROOT.RooArgSet()
        
        self.__suffix      = suffix

        if   isinstance ( signal_1 , PDF            )          : self.__signal1 = signal_1
        elif isinstance ( signal_1 , ROOT.RooAbsPdf ) and xvar :
            self.__signal1 = Generic1D_pdf ( signal_1 , xvar , 'SIGNAL-X' )
        else : raise AttributeError("Invalid ``signal1'' attribute" )
            
        if   isinstance ( signal_2 , PDF            )          : self.__signal2 = signal_2
        elif isinstance ( signal_2 , ROOT.RooAbsPdf ) and yvar :
            self.__signal2 = Generic1D_pdf ( signal_2 , yvar , 'SIGNAL-Y' )
        else : raise AttributeError("Invalid ``signal2'' attribute" )
            
        #
        ## initialize base class
        #
        if not name and self.__signal1.name and self.__signal2.name :
            name = '%s_and_%s_%s' % ( self.__signal1.name ,
                                      self.__signal2.name , suffix )
            
        PDF2.__init__ ( self , name , self.__signal1.xvar , self.__signal2.xvar ) 

        # =====================================================================
        ## Build components for the  final 2D-PDF
        # =====================================================================
        
        #
        ## First component: Signal(1) and Signal(2)
        # 
        self.__ss_cmp = Model2D ( "SS_pdf" + suffix ,
                                  self.__signal1            ,
                                  self.__signal2            , 
                                  title = "Sig(1) x Sig(2)" )
        #
        ## Second component: Background(1) and Signal(2)
        # 
        self.__bkg1   = makeBkg ( bkg1   , 'Bkg(1)' + suffix  , self.xvar )
        self.__bs_cmp = Model2D ( "BS_pdf" + suffix         ,
                                  self.__bkg1               ,
                                  self.__signal2            ,
                                  title = "Bkg(1) x Sig(2)" )
        
        #
        ## Third component:  Signal(1) and Background(2)
        # 
        self.__bkg2   = makeBkg ( bkg2   , 'Bkg(2)' + suffix  , self.yvar )
        self.__sb_cmp = Model2D ( "SB_pdf" + suffix         ,
                                  self.__signal1            ,
                                  self.__bkg2               ,
                                  title = "Sig(1) x Bkg(2)" )
            
        #
        ## fourth component: Background(1) and Background(2) 
        #
        if bkg2D : self.__bb2D  = bkg2D
        #
        self.__bkgA = None 
        self.__bkgB = None 

        if   isinstance ( bkg2D , PDF2           ) : self.__bb_cmp = bkg2D  
        elif isinstance ( bkg2D , ROOT.RooAbsPdf ) :
            self.__bb_cmp  = Generic2D_pdf  ( bkg2D , self.xvar , self.yvar )
        else     :            
            
            if bkgA is None : bkgA = bkg1
            if bkgB is None : bkgB = bkg2
            
            self.__bkgA = makeBkg ( bkgA   , 'Bkg(A)' + suffix , self.xvar )
            self.__bkgB = makeBkg ( bkgB   , 'Bkg(B)' + suffix , self.yvar )
            
            self.__bb_cmp = Model2D ( "BB_pdf" + suffix         ,
                                      self.__bkgA               ,
                                      self.__bkgB               ,
                                      title = "Bkg(A) x Bkg(B)" )
        #
        ## coefficients
        #
        self.__ss = makeVar ( ss   ,
                              "SS"            + suffix ,
                              "Sig(1)&Sig(2)" + suffix , None , 1000  , 0 , 1.e+7 )
        self.__sb = makeVar ( sb   ,
                              "SB"            + suffix ,
                              "Sig(1)&Bkg(2)" + suffix , None ,  100  , 0 , 1.e+7 )
        self.__bs = makeVar ( bs   ,
                              "BS"            + suffix ,
                              "Bkg(1)&Sig(2)" + suffix , None ,  100  , 0 , 1.e+7 )        
        self.__bb = makeVar ( bb   ,
                              "BB"            + suffix ,
                              "Bkg(1)&Bkg(2)" + suffix , None ,   10  , 0 , 1.e+7 )
        
        self.alist1 = ROOT.RooArgList (
            self.__ss_cmp.pdf ,
            self.__sb_cmp.pdf ,
            self.__bs_cmp.pdf ,
            self.__bb_cmp.pdf )
        self.alist2 = ROOT.RooArgList (
            self.__ss ,
            self.__sb ,
            self.__bs ,
            self.__bb )

        ## treat additional components (if specified)
        self.__nums_components = [] 
        icmp = 0
        self.__more_components = []
        for cmp in components :
            
            if   isinstance  ( c , PDF2           ) : cc = c  
            elif isinstance  ( c , ROOT.RooAbsPdf ) : cc = Generic2D_pdf ( cs ,  self.xvar , self.yvar ) 
            else :
                logger.error ("Fit2D(%s): Unknown ``other''component %s/%s, skip it!" % ( self.name , cc , type(cc) ) )
                continue  
            self.__more_components.append ( cc     )
            self.components.add           ( cc.pdf ) 

        nc = len( self.__more_components )
        if 1 == nc :
            cf = makeVar ( None , "C"+suffix , "Component" + suffix , None , 1 , 0 , 1.e+7 )
            self.alist1.add  ( self.components[0] )
            self.__num_components.append ( cf ) 
        elif 2 <= nc : 
            fic = makeFracs ( nc , 'C_%%d%s' % suffix ,  'C(%%d)%s'  % suffix , fractions  = False , model = self )
            for c in self.components : self.alist1.add ( c)
            for f in fic             : self.__num_components.append ( f )
            
        self.__nums_components  = tuple ( self.__nums_components  ) 
        for c in self.__nums_components  : self.alist2.add ( c )

        #
        ## build the final PDF 
        # 
        self.pdf  = ROOT.RooAddPdf  ( "model2D"      + suffix ,
                                      "Model2D(%s)"  % suffix ,
                                      self.alist1 ,
                                      self.alist2 )

        self.signals     .add ( self.__ss_cmp.pdf )
        self.backgrounds .add ( self.__bb_cmp.pdf )
        self.crossterms1 .add ( self.__sb_cmp.pdf ) ## cross-terms 
        self.crossterms2 .add ( self.__bs_cmp.pdf ) ## cross-terms 

        ## save configuration
        self.config = {
            'signal_1'   : self.signal1         ,
            'signal_2'   : self.signal2         ,            
            'suffix'     : self.suffix          ,
            'bkg1'       : self.bkg1            , 
            'bkg2'       : self.bkg2            , 
            'bkgA'       : self.bkgA            , 
            'bkgB'       : self.bkgB            , 
            'bkg2D'      : self.bkg2D           ,
            'ss'         : self.SS              ,
            'sb'         : self.SB              ,
            'bs'         : self.BS              ,
            'BB'         : self.BB              ,
            'components' : self.more_components ,
            'xvar'       : self.xvar            , 
            'yvar'       : self.yvar            , 
            'name'       : self.name    
            }
    
    @property
    def SS ( self ) :
        """The yield of Signal(x)*Signal(y) component"""
        return self.__ss
    @SS.setter 
    def SS ( self , value ) :
        value = float ( value  )
        assert value in self.__ss, "Value %s is out of the allowed range %s " % ( value , self.__ss.minmax() )
        self.__ss.setVal ( value ) 

    @property
    def SB ( self ) :
        """The yield of Signal(x)*Background(y) component"""
        return self.__sb
    @SB.setter 
    def SB ( self , value ) :
        value = float ( value  )
        assert value in self.__sb, "Value %s is out of the allowed range %s " % ( value , self.__sb.minmax() )
        self.__sb.setVal ( value ) 

    @property
    def BS ( self ) :
        """The yield of Background(x)*Signal(y) component"""
        return self.__bs
    @BS.setter 
    def BS ( self , value ) :
        value = float ( value  )
        assert value in self.__bs, "Value %s is out of the allowed range %s " % ( value , self.__bs.minmax() )
        self.__bs.setVal ( value ) 

    @property
    def BB ( self ) :
        """The yield of Background(x,y) component"""
        return self.__bb
    @BB.setter 
    def BB ( self , value ) :
        value = float ( value  )
        assert value in self.__bb, "Value %s is out of the allowed range %s " % ( value , self.__bb.minmax() )
        self.__bb.setVal ( value ) 
    
    @property
    def C ( self ) :
        """Get the  yields of ``other'' component(s) 
        For single ``other'' component:
        >>> print pdf.C           ## read the single ``other'' component 
        >>> pdf.C = 100           ## assign to it 
        For multiple ``other'' components:
        >>> print pdf.C[4]        ## read the 4th ``other'' component 
        >>> pdf.C = 4,100         ## assign to it 
        ... or, alternatively:
        >>> print pdf.C[4]        ## read the 4th ``other'' component 
        >>> pdf.C[4].value 100    ## assign to it         
        """
        lst = [ i for i in self.__nums_components ]
        if not lst          : return ()     ## extended fit? no other components?
        elif  1 == len(lst) : return lst[0] ## single component?
        return tuple ( lst )
    @C.setter
    def C (  self , value ) :
        _n = len ( self.__nums_components )
        assert 1 <= _n , "No ``other'' components are defined, assignement is impossible"
        if 1 ==  _n :
            _c    = self.C 
            value = float ( value )
        else : 
            index = value [0]
            assert isinstance ( index , int ) and 0 <= index < _n, "Invalid ``other'' index %s/%d" % ( index , _n ) 
            value = float ( value[1] )
            _c    = self.C[index]
        ## assign 
        assert value in _c , "Value %s is outside the allowed region %s"  % ( value , _c.minmax() )
        _c.setVal ( value )

    @property
    def  yields    ( self ) :
        """The list/tuple of the yields of all numeric components"""
        return tuple ( [ i for i in  self.alist2 ] )

    # =========================================================================
    # components
    # =========================================================================

    @property 
    def signal1 ( self  ) :
        """Signal(x) component/PDF"""
        return self.__signal1

    @property 
    def signal2 ( self  ) :
        """Signal(y) component/PDF"""
        return self.__signal2

    @property
    def bkg1( self ) :
        """ The background PDF for Backgroud(x)*Signal(y) component/PDF"""
        return self.__bkg1
    
    @property
    def bkg2( self ) :
        """ The background PDF for Signal(x)*Background(y) component/PDF"""
        return self.__bkg2

    @property
    def bkgA( self ) :
        """ The background(x) PDF for Backgroud(x)*Background(y) component/PDF"""
        return self.__bkgA

    @property
    def bkgB( self ) :
        """ The background(y) PDF for Backgroud(x)*Background(y) component/PDF"""
        return self.__bkgB

    @property
    def bkg2D( self ) :
        """ The PDF for Backgroud(x,y) component/PDF"""
        return self.__bb_cmp
 
    @property 
    def crossterms1 ( self ) :
        """``cross-terms'': pdfs for signal(x) and backgrond(y)"""        
        return self.__crossterms1

    @property
    def crossterms2 ( self ) :
        """``cross-terms'':  pdfs for background(x) and signal(y) """        
        return self.__crossterms2

    @property
    def more_components ( self ) :
        """additional ``other'' components"""
        return tuple( self.__more_components  )

    @property
    def suffix ( self ) :
        """``suffix'', used to build the name"""
        return self.__suffix
    
# =============================================================================
## @class Fit2DSym
#  The actual model for 2D-fits
#
#  @code
# 
#  model   = Models.Fit2D (
#      signal_1 = Models.Gauss_pdf ( 'Gx' , mass = m_x ) ,
#      signal_2 = Models.Gauss_pdf ( 'Gy' , mass = m_y ) ,
#      bkg1     = 1 , 
#      bkg2     = 1 )
#
#  r,f = model.fitTo ( dataset ) ## fit dataset 
#
#  print r                       ## get results  
#
#  fx  = model.draw1 ()          ## visualize X-projection
#  fy  = model.draw2 ()          ## visualize X-projection
#
#  @endcode 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class Fit2DSym (PDF2) :
    """The actual model for *symmetric**2D-fits
    
    >>>  model   = Models.Fit2D (
    ...      signal_1 = Models.Gauss_pdf ( 'Gx' , xvar = m_x ) ,
    ...      signal_2 = Models.Gauss_pdf ( 'Gy' , xvar = m_y ) ,
    ...      bkg1     = 1 , 
    ...      bkg2     = 1 )
    >>> r,f = model.fitTo ( dataset ) ## fit dataset 
    >>> print r                       ## get results  
    >>> fx  = model.draw1 ()          ## visualize X-projection
    >>> fy  = model.draw2 ()          ## visualize X-projection

    """
    def __init__ ( self               ,
                   #
                   signal_1           , 
                   signal_2           ,
                   suffix = ''        ,
                   #
                   bkg1       = None  ,
                   bkgA       = None  ,
                   bkg2D      = None  ,
                   #
                   ## main components :
                   ss         = None  , ## signal (1) * signal     (2)
                   sb         = None  , ## signal     * background 
                   bb         = None  , ## background * background  
                   ## additional components 
                   components = []    ,
                   xvar       = None  ,
                   yvar       = None  ,                   
                   name       = ''    ) :
        
        ## collect all the arguments 
        self.__args = {
            'signal_1'   : signal_1   , 'signal_2' : signal_2 ,
            'bkg1'       : bkg1       , 
            'bkgA'       : bkgA       , 
            'bkg2D'      : bkg2D      ,
            'components' : components ,
            ##
            'ss'         : ss         , 'bb'       : bb      , 'sb'        : sb , 
            ##
            'suffix'     : suffix     ,
            'name'       : name       ,
            ##
            'xvar'       : xvar , 'yvar'       : yvar ,            
            }
        
        self.__crossterms1 = ROOT.RooArgSet()
        self.__crossterms2 = ROOT.RooArgSet()
        
        self.__suffix      = suffix

        if   isinstance ( signal_1 , PDF            )          : self.__signal1 = signal_1
        elif isinstance ( signal_1 , ROOT.RooAbsPdf ) and xvar :
            self.__signal1 = Generic1D_pdf ( signal_1 , xvar , 'SIGNAL-X' )
        else : raise AttributeError ( "Invalid ``signal1'' attribute" )
            
        if   isinstance ( signal_2 , PDF            )          : self.__signal2 = signal_2
        elif isinstance ( signal_2 , ROOT.RooAbsPdf ) and yvar :
            self.__signal2 = Generic1D_pdf ( signal_2 , yvar , 'SIGNAL-Y' )
        else : raise AttributeError ( "Invalid ``signal2'' attribute" )
            
        #
        ## initialize base class
        #
        if not name and self.__signal1.name and self.__signal2.name :
            name = self.__signal1.name +'_AND_' + self.__signal2.name + '_'+ suffix

        ##  initialize the base class 
        PDF2.__init__ ( self , name , self.__signal1.xvar , self.__signal2.xvar ) 
        
        #
        ## First component: Signal(1) and Signal(2)
        # 
        self.__ss_cmp = Model2D ( "SS_pdf" + suffix         ,
                                  self.__signal1            ,
                                  self.__signal2            , 
                                  title = "Sig(1) x Sig(2)" )
        
        self.__bkg1  = makeBkg (      bkg1 , 'Bkg(1)' + suffix , self.xvar )
        self.__bkg2  = makeBkg ( self.bkg1 , 'Bkg(2)' + suffix , self.yvar )

        #
        ## Second sub-component: Background (1) and Signal     (2)
        ## Third  sub-component: Signal     (1) and Background (2)
        # 
        self.__sb_cmp_raw = Model2D ( "S1B2_pdf" + suffix       ,
                                      self.__signal1            ,
                                      self.__bkg2               ,
                                      title = "Sig(1) x Bkg(2)" )
        
        self.__bs_cmp_raw = Model2D ( "B1S2_pdf" + suffix       ,                                      
                                      self.__bkg1               ,
                                      self.__signal2            ,    
                                      title = "Bkg(1) x Sig(2)" )
        
        self.__f_cross    = ROOT.RooConstVar ( 'SxB_fraction'   + suffix  ,
                                               '(S1B2-vs-B1S2) fraction' , 0.5 )
        
        self.__sb_cmp     = Generic2D_pdf (
            ROOT.RooAddPdf ( "SB_pdf" + suffix ,
                             "Sig(1) x Bkg(2) + Bkg(1) x Sig(2)"   ,
                             self.__sb_cmp_raw.pdf ,
                             self.__bs_cmp_raw.pdf ,
                             self.__f_cross        ) , self.xvar , self.yvar )
        
        ## alias, just for convinience 
        self.__bs_cmp    = self.__sb_cmp
        
        #
        ## fourth component: Background(1) and Background(2) 
        #
        self.__bkgA = None
        self.__bkgB = None
        
        if   isinstance ( bkg2D , PDF2           ) : self.__bb_cmp = bkg2D  
        elif isinstance ( bkg2D , ROOT.RooAbsPdf ) :
            self.__bb_cmp  = Generic2D_pdf  ( bkg2D , self.xvar , self.yvar )
        else     :            
            
            if bkgA is None : bkgA = bkg1
            
            self.__bkgA = makeBkg (      bkgA , 'Bkg(A)' + suffix , self.xvar )
            self.__bkgB = makeBkg ( self.bkgA , 'Bkg(B)' + suffix , self.yvar )
            
            self.__bb_cmp = Model2D ( "BB_pdf" + suffix         ,
                                      self.__bkgA               ,
                                      self.__bkgB               ,
                                      title = "Bkg(1) x Bkg(2)" )
        #
        ## coefficients
        #
        self.__ss = makeVar ( ss                       ,
                              "SS"            + suffix ,
                              "Sig(1)&Sig(2)" + suffix , None , 1000  , 0 ,  1.e+7 )
        
        self.__bb = makeVar ( bb                       ,
                              "BB"            + suffix ,
                              "Bkg(1)&Bkg(2)" + suffix , None ,   10  , 0 ,  1.e+7 )
        
        self.__sb = makeVar ( sb                       ,
                              "SB"            + suffix ,
                              "Sig(1)&Bkg(2)+Bkg(1)&Sig(2)" + suffix , None ,  100  , 0 ,  1.e+7 )
        
        ## duplicate 
        self.__bs = self.__sb
        
        self.alist1 = ROOT.RooArgList (
            self.__ss_cmp.pdf ,
            self.__sb_cmp.pdf ,
            self.__bb_cmp.pdf )
        self.alist2 = ROOT.RooArgList (
            self.__ss         ,
            self.__sb         ,
            self.__bb         )

        ## treat additional components (if specified)
        self.__nums_components = [] 
        icmp = 0
        self.__more_components = []
        for cmp in components :
            
            if   isinstance  ( c , PDF2           ) : cc = c  
            elif isinstance  ( c , ROOT.RooAbsPdf ) : cc = Generic2D_pdf ( cs ,  self.xvar , self.yvar ) 
            else :
                logger.error ("Fit2D(%s): Unknown ``other''component %s/%s, skip it!" % ( self.name , cc , type(cc) ) )
                continue  
            self.__more_components.append ( cc     )
            self.components.add           ( cc.pdf ) 

        nc = len( self.__more_components )
        if 1 == nc :
            cf = makeVar ( None , "C"+suffix , "Component" + suffix , None , 1 , 0 , 1.e+7 )
            self.alist1.add  ( self.components[0] )
            self.__num_components.append ( cf ) 
        elif 2 <= nc : 
            fic = makeFracs ( nc , 'C_%%d%s' % suffix ,  'C(%%d)%s'  % suffix , fractions  = False , model = self )
            for c in self.components : self.alist1.add ( c)
            for f in fic             : self.__num_components.append ( f )
            
        self.__nums_components  = tuple ( self.__nums_components  ) 
        for c in self.__nums_components  : self.alist2.add ( c )
            
        #
        ## build the final PDF 
        # 
        self.pdf  = ROOT.RooAddPdf  ( "model2D"      + suffix ,
                                      "Model2D(%s)"  % suffix ,
                                      self.alist1 ,
                                      self.alist2 )


        self.signals     .add ( self.__ss_cmp.pdf )
        self.backgrounds .add ( self.__bb_cmp.pdf )
        self.crossterms1 .add ( self.__sb_cmp.pdf ) ## cross-terms 
        self.crossterms2 .add ( self.__bs_cmp.pdf ) ## cross-terms 

        ## save configuration
        self.config = {
            'signal_1'   : self.signal1         ,
            'signal_2'   : self.signal2         ,            
            'suffix'     : self.suffix          ,
            'bkg1'       : self.bkg1            , 
            'bkgA'       : self.bkgA            , 
            'bkg2D'      : self.bkg2D           ,
            'ss'         : self.SS              ,
            'sb'         : self.SB              ,
            'bb'         : self.BB              ,
            'components' : self.more_components ,
            'xvar'       : self.xvar            , 
            'yvar'       : self.yvar            , 
            'name'       : self.name            ,
            }
    
    @property
    def SS ( self ) :
        """The yield of Signal(x)*Signal(y) component"""
        return self.__ss
    @SS.setter 
    def SS ( self , value ) :
        value = float ( value  )
        assert value in self.__ss, "Value %s is out of the allowed range %s " % ( value , self.__ss.minmax() )
        self.__ss.setVal ( value ) 

    @property
    def SB ( self ) :
        """The yield of Signal(x)*Background(y)+Background(x)*Signal(y) component (same as ``BS'')"""
        return self.__sb
    @SB.setter 
    def SB ( self , value ) :
        value = float ( value  )
        assert value in self.__sb, "Value %s is out of the allowed range %s " % ( value , self.__sb.minmax() )
        self.__sb.setVal ( value ) 

    @property
    def BS ( self ) :
        """The yield of Signal(x)*Background(y)+Background(x)*Signal(y) component (same as ``SB'')"""
        return self.SB 
    @BS.setter 
    def BS ( self , value ) :
        self.SB = value
        return self.SB.getVal()

    @property
    def BB ( self ) :
        """The yield of Background(x,y) component"""
        return self.__bb
    @BB.setter 
    def BB ( self , value ) :
        value = float ( value  )
        assert value in self.__bb, "Value %s is out of the allowed range %s " % ( value , self.__bb.minmax() )
        self.__bb.setVal ( value ) 
    
    @property
    def C ( self ) :
        """Get the  yields of ``other'' component(s) 
        For single ``other'' component:
        >>> print pdf.C           ## read the single ``other'' component 
        >>> pdf.C = 100           ## assign to it 
        For multiple ``other'' components:
        >>> print pdf.C[4]        ## read the 4th ``other'' component 
        >>> pdf.C = 4,100         ## assign to it 
        ... or, alternatively:
        >>> print pdf.C[4]        ## read the 4th ``other'' component 
        >>> pdf.C[4].value 100    ## assign to it         
        """
        lst = [ i for i in self.__nums_components ]
        if not lst          : return ()     ## extended fit? no other components?
        elif  1 == len(lst) : return lst[0] ## single component?
        return tuple ( lst )
    @C.setter
    def C (  self , value ) :
        _n = len ( self.__nums_components )
        assert 1 <= _n , "No ``other'' components are defined, assignement is impossible"
        if 1 ==  _n :
            _c    = self.C 
            value = float ( value )
        else : 
            index = value [0]
            assert isinstance ( index , int ) and 0 <= index < _n, "Invalid ``other'' index %s/%d" % ( index , _n ) 
            value = float ( value[1] )
            _c    = self.C[index]
        ## assign 
        assert value in _c , "Value %s is outside the allowed region %s"  % ( value , _c.minmax() )
        _c.setVal ( value )

    @property
    def  yields    ( self ) :
        """The list/tuple of the yields of all numeric components"""
        return tuple ( [ i for i in  self.alist2 ] )

    # =========================================================================
    # components
    # =========================================================================
    
    @property 
    def signal1 ( self  ) :
        """Signal(x) component/PDF"""
        return self.__signal1

    @property 
    def signal2 ( self  ) :
        """Signal(y) component/PDF"""
        return self.__signal2

    @property
    def bkg1( self ) :
        """ The background PDF for Backgroud(x)*Signal(y) component/PDF"""
        return self.__bkg1
    
    @property
    def bkg2( self ) :
        """ The background PDF for Signal(x)*Background(y) component/PDF"""
        return self.__bkg2

    @property
    def bkgA( self ) :
        """ The background(x) PDF for Backgroud(x)*Background(y) component/PDF"""
        return self.__bkgA

    @property
    def bkgB( self ) :
        """ The background(y) PDF for Backgroud(x)*Background(y) component/PDF"""
        return self.__bkgB

    @property
    def bkg2D( self ) :
        """ The PDF for Backgroud(xmy) component/PDF"""
        return self.__bb_cmp
 
    @property 
    def crossterms1 ( self ) :
        """``cross-terms'': pdfs for signal(x) and backgrond(y)"""        
        return self.__crossterms1

    @property
    def crossterms2 ( self ) :
        """``cross-terms'':  pdfs for background(x) and signal(y) """        
        return self.__crossterms2

    @property
    def more_components ( self ) :
        """additional ``other'' components"""
        return tuple( self.__more_components  )

    @property
    def suffix ( self ) :
        """``suffix'', used to build the name"""
        return self.__suffix
 
# =============================================================================
## @class Generic2D_pdf
#  "Wrapper" over generic RooFit (2D)-pdf
#  @code
#  raw_pdf = 
#  pdf     = Generic2D_pdf ( raw_pdf , xvar = ... , yvar = ... )  
#  @endcode 
#  If more functionality is required , more actions are possible:
#  @code
#  ## for sPlot 
#  pdf.alist2 = ROOT.RooArgList ( n1 , n2 , n3 ) ## for sPlotting 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-03-29
class Generic2D_pdf(PDF2) :
    """ Wrapper for generic (2D) RooFit pdf    
    >>> raw_pdf = 
    >>> pdf     = Generic2D_pdf ( raw_pdf , xvar = ... , yvar = ... )
    """
    ## constructor 
    def __init__ ( self , pdf , xvar , yvar , name = None ) :

        assert isinstance ( pdf  , ROOT.RooAbsPdf  ) , "``pdf''  must be ROOT.RooAbsPdf"
        assert isinstance ( xvar , ROOT.RooAbsReal ) , "``xvar'' must be ROOT.RooAbsReal"
        assert isinstance ( yvar , ROOT.RooAbsReal ) , "``yvar'' must be ROOT.RooAbsReal"
        
        if not name : name = pdf.GetName()
        PDF2  . __init__ ( self , name , xvar , yvar )
        self.pdf = pdf
        
        self.signals.add ( self.pdf )
       
        ## save the configuration
        self.config = {
            'pdf'  : self.pdf  ,
            'xvar' : self.xvar ,
            'yvar' : self.yvar ,
            'name' : self.name 
            }

    ## redefine the clone method, allowing only the name to be changed
    #  @attention redefinition of parameters and variables is disabled,
    #             since it can't be done in a safe way                  
    def clone ( self , name = '' , xvar  = None , yvar = None ) :
        """Redefine the clone method, allowing only the name to be changed
         - redefinition of parameters and variables is disabled,
         since it can't be done in a safe way          
        """
        if xvar and not xvar is self.xvar :
            raise AttributeError("Generic2D_pdf can not be cloned with different `xvar''")
        if yvar and not yvar is self.yvar :
            raise AttributeError("Generic2D_pdf can not be cloned with different `yvar''")
        return PDF.clone ( self , name = name ) if name else PDF.clone( self )
        
# =============================================================================
## simple convertor of 2D-histogram into PDF
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H2D_pdf(H2D_dset,PDF2) :
    """Simple convertor of 2D-histogram into PDF 
    """
    def __init__ ( self            ,
                   name            ,
                   histo           ,
                   xvar  = None    , 
                   yvar  = None    ,
                   density = True  ,
                   silent  = False ) :
        
        H2D_dset.__init__ ( self , histo , xvar , yvar , density , silent )
        PDF2    .__init__ ( self , name  , self.xaxis , self.yaxis ) 
        
        self.__vset  = ROOT.RooArgSet  ( self.xvar , self.yvar )
        
        #
        ## finally create PDF :
        #
        from   Ostap.Utils      import roo_silent 
        with roo_silent ( silent ) : 
            self.pdf    = ROOT.RooHistPdf (
                'hpdf_%s'            % name ,
                'Histo2PDF(%s/%s/%s)' % ( name , self.histo.GetName() , self.histo.GetTitle() ) , 
                self.__vset , 
                self.dset   )
            
        ## and declare it be be a "signal"
        self.signals.add ( self.pdf ) 

        ## save the configuration
        self.config = {
            'name'    : self.name    , 
            'histo'   : self.histo   , 
            'xvar'    : self.xvar    , 
            'yvar'    : self.yvar    , 
            'density' : self.density , 
            'silent'  : self.silent  ,             
            }

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
