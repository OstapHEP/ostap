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
from   ostap.fitting.basic import PDF, makeVar, makeBkg
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
    """Useful helper base class for implementation of PDFs for 2D-fit
    """
    def __init__ ( self , name , xvar = None , yvar = None ) : 


        PDF.__init__ ( self , name , xvar )
        
        self.__yvar = None 
        ## create the variable 
        if isinstance ( yvar , tuple ) and 2 == len(yvar) :  
            self.__yvar = makeVar ( yvar               , ## var 
                                    'y'                , ## name 
                                    'y-varibale(mass)' , ## title/comment
                                    *yvar              , ## min/max 
                                    fix = None         ) ## fix ? 
        elif isinstance ( yvar , ROOT.RooAbsReal ) :
            self.__yvar = makeVar ( yvar               , ## var 
                                    'y'                , ## name 
                                    'y-variable/mass'  , ## title/comment
                                    fix = None         ) ## fix ? 
        else :
            ##logger.warning('x-varibale is not specified (yet)')
            self.__yvar = makeVar( yvar , 'y' , 'y-variable' )


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
        """Perform the actual fit (and draw it)
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
        """Draw the projection over 2nd variable
        
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
        """Make/draw 2D-histograms 
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
            hdata.lego ()
            hpdf .Draw ( 'same surf')
        
        return hpdf , hdata 
    
    # =========================================================================
    ## make 1D-plot
    #  @code
    #  model.fitTo ( dataset , ... )
    #  frame = model.draw( var1 , dataset , ... )  
    #  @endcode 
    def draw ( self                         ,
               drawvar               = None ,
               dataset               = None ,
               nbins                 =  100 ,
               ybins                 =   20 ,
               silent                = True ,
               in_range              = None ,
               **kwargs                     ) : 
        """ Make 1D-plot:
        >>> model.fitTo ( dataset , ... )
        >>> frame = model.draw( var1 , dataset , ... )  
        """
                #
        ## special case:  do we need it? 
        # 
        if not drawvar : return self.draw_H2D( dataset , nbins , ybins )

        ## copy arguments:
        args = kwargs.copy ()
        
        import ostap.plotting.fit_draw as FD
        if in_range :
            data_options        = args.pop (       'data_options' , FD.         data_options )
            background_options  = args.pop ( 'background_options' , FD. background2D_options )
            signal_options      = args.pop (     'signal_options' , FD.       signal_options )
            component_options   = args.pop (  'component_options' , FD.    component_options )
            crossterm1_options  = args.pop ( 'crossterm1_options' , FD.   crossterm1_options )
            crossterm2_options  = args.pop ( 'crossterm2_options' , FD.   crossterm2_options )
            total_fit_options   = args.pop (  'total_fit_options' , FD.    total_fit_options )
            
            data_options       += ROOT.RooFit.CutRange        ( in_range ) , 
            signal_options     += ROOT.RooFit.ProjectionRange ( in_range ) , 
            background_options += ROOT.RooFit.ProjectionRange ( in_range ) , 
            component_options  += ROOT.RooFit.ProjectionRange ( in_range ) , 
            crossterm1_options += ROOT.RooFit.ProjectionRange ( in_range ) , 
            crossterm2_options += ROOT.RooFit.ProjectionRange ( in_range ) , 
            total_fit_options  += ROOT.RooFit.ProjectionRange ( in_range ) , 
            
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
    #  @code
    #  histo = ...
    #  r,f = model.fitHisto ( histo )
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
            self.hdset = H2D_dset ( histo , self.xvar , self.yvar  , density , silent )
            self.hset  = self.hdset.dset
                
            ## fit it!!
            return self.fitTo ( self.hset      ,
                                draw           ,
                                histo.nbinsx() ,
                                histo.nbinsy() ,
                                silent         , *args , **kwargs )

    # =================================================================================
    ## simple 'function-like' interface
    #  @code
    #  pdf = ...
    #  x, y = 0.45, 0.88 
    #  print 'Value of PDF at x=%f,y=%s is %f' % ( x , y , pdf ( x , y ) ) 
    #  @endcode
    def __call__ ( self , x , y ) :
        """Simple 'function-like' interface
        >>> pdf = ...
        >>> x, y = 0.45, 0.88 
        >>> print 'Value of PDF at x=%f,y=%s is %f' % ( x , y , pdf ( x , y ) ) 
        """        
        if isinstance ( self.xvar , ROOT.RooRealVar ) and isinstance ( self.yvar , ROOT.RooRealVar ) :
            from ostap.fitting.roofit import SETVAR
            if x in self.xvar and y in  self.xvar : 
                with SETVAR( self.xvar ) , SETVAR( self.yvar ) :
                    self.xvar.setVal ( x )
                    self.yvar.setVal ( y )
                    return self.pdf.getVal()
            else : return 0.0
            
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

        ## 1) make a try to use analytical integral (could be fast)
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
            
        ## 2) use numerical integration 
        from ostap.math.integral import integral2 as _integral2 
        return _integral2 ( self ,
                            xmin , xmax ,
                            ymin , ymax )
      
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
#  fy  = model.draw2 ()          ## visualize X-projection
#
#  @endcode 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class Fit2D (PDF2) :
    """The actual model for 2D-fits
    
    >>>  model   = Models.Fit2D (
    ...      signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
    ...      signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y ) ,
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
                   bkg2       = None  ,
                   bkgA       = None  ,
                   bkgB       = None  ,
                   bkg2D      = None  ,
                   #
                   ## main components :
                   ss         = None  , ## signal    (1) * signal     (2)
                   sb         = None  , ## signal    (1) * background (2) 
                   bs         = None  , ## background(1) * signal     (2)
                   bb         = None  , ## background-2D 
                   ## additional components 
                   components = []    ,
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
            'suffix'   : suffix   ,
            'name'     : name     ,
            }
        
        self.__crossterms1 = ROOT.RooArgSet()
        self.__crossterms2 = ROOT.RooArgSet()
        
        self.suffix    = suffix 
        self.__signal1 = signal_1
        self.__signal2 = signal_2

        #
        ## initialize base class
        #
        if not name and signal_1.name and signal_2.name :
            name = '%s_and_%s_%s' % ( signal_1.name , signal_2.name , suffix )
            
        PDF2.__init__ ( self , name , signal_1.mass , signal_2.mass ) 
        
        #
        ## First component: Signal(1) and Signal(2)
        # 
        self.__ss_pdf = ROOT.RooProdPdf ( "S1S2pdf" + suffix ,
                                          "Sig(1) x Sig(2)"  ,
                                          self.signal1.pdf   ,
                                          self.signal2.pdf   )
        
        self.__bkg1 = makeBkg ( bkg1   , 'Bkg(1)' + suffix , self.xvar )
        
        #
        ## Second component: Background(1) and Signal(2)
        # 
        self.__bs_pdf = ROOT.RooProdPdf ( "B1S2pdf" + suffix  ,
                                          "Bkg(1) x Sig(2)"   ,
                                          self.__bkg1.pdf     ,
                                          self.signal2.pdf    )
        
        self.__bkg2 = makeBkg ( bkg2   , 'Bkg(2)' + suffix , self.yvar )
        
        #
        ## Third component:  Signal(1) and Background(2)
        # 
        self.__sb_pdf = ROOT.RooProdPdf ( "S1B2pdf" + suffix  ,
                                          "Sig(1) x Bkg(2)"   ,
                                          self.signal1.pdf    ,
                                          self.__bkg2.pdf     )
        
        ## 
        self.__bkgs = ( bkg1 , bkg2 , bkgA , bkgB ) 
        #
        ## fourth component: Background(1) and Background(2) 
        #
        if bkg2D : self.__bb2D  = bkg2D
        #
        self.__bkgA = None 
        self.__bkgB = None 

        if   bkg2D and isinstance ( bkg2D , ROOT.RooAbsPdf ) : self.__bb_pdf = bkg2D 
        elif bkg2D and hasattr    ( bkg2D , 'pdf'          ) : self.__bb_pdf = bkg2D.pdf
        else     :            

            if bkgA is None : bkgA = bkg1
            if bkgB is None : bkgB = bkg2
            
            self.__bkgA = makeBkg ( bkgA   , 'Bkg(A)' + suffix , self.xvar )
            self.__bkgB = makeBkg ( bkgB   , 'Bkg(B)' + suffix , self.yvar )
            
            self.__bb_pdf = ROOT.RooProdPdf ( "B1B2pdf" + suffix ,
                                              "Bkg(A) x Bkg(B)"  ,
                                              self.__bkgA.pdf    ,
                                              self.__bkgB.pdf    )
        #
        ## coefficients
        #
        self.__ss = makeVar ( ss   ,
                              "S1S2"          + suffix ,
                              "Sig(1)&Sig(2)" + suffix , None , 1000  , 0 , 1.e+7 )
        self.__sb = makeVar ( sb   ,
                              "S1B2"          + suffix ,
                              "Sig(1)&Bkg(2)" + suffix , None ,  100  , 0 , 1.e+7 )
        self.__bs = makeVar ( bs   ,
                              "B1S2"          + suffix ,
                              "Bkg(1)&Sig(2)" + suffix , None ,  100  , 0 , 1.e+7 )        
        self.__bb = makeVar ( bb   ,
                              "B1B2"          + suffix ,
                              "Bkg(1)&Bkg(2)" + suffix , None ,   10  , 0 , 1.e+7 )
        
        self.SS_name = self.S1S2.GetName()
        self.SB_name = self.S1B2.GetName()
        self.BS_name = self.B1S2.GetName()
        self.BB_name = self.B1B2.GetName()
        
        self.alist1 = ROOT.RooArgList (
            self.__ss_pdf ,
            self.__sb_pdf ,
            self.__bs_pdf ,
            self.__bb_pdf )
        self.alist2 = ROOT.RooArgList (
            self.__ss ,
            self.__sb ,
            self.__bs ,
            self.__bb )

        #
        ## treat additional components (if specified)
        # 
        self.other       = components
        self._cmps       = []
        icmp = 0 
        for cmp in self.other :

            icmp += 1
            
            if   isinstance ( cmp , ROOT.RooAbsPdf         ) : pass 
            elif hasattr    ( cmp , 'pdf'                  ) :
                self._cmps += [ cmp ] 
                cmp = cmp.pdf 
            elif isinstance ( cmp , ( float , int , long ) ) and not isinstance ( cmp , bool ) :
                px  = ROOT.RooPolynomial ( 'Px%d'    % icmp + suffix ,
                                           'Px(%d)'  % icmp + suffix , self.xvar ) 
                py  = ROOT.RooPolynomial ( 'Py%d'    % icmp + suffix ,
                                           'Py(%d)'  % icmp + suffix , self.yvar) 
                cmp = ROOT.RooProdPdf    ( "Pxy%d"   % icmp + suffix ,
                                           "Pxy(%d)" % icmp + suffix , px , py )  
                self._cmps += [ px,py,cmp]
            else :
                logger.error( 'Unknown type of component %d %d ' % ( imcp , type(cmp) ) )

                
            nn = makeVar ( None ,
                           'Cmp%d'   % icmp + suffix ,
                           'Cmp(%d)' % icmp + suffix ,
                           None ,  100  ,  0 , 1.e+6 )  
            self._cmps += [ nn ]

            self.alist1.add ( cmp )
            self.alist2.add ( nn  )
            
            self.components ().add ( cmp ) 
            
        #
        ## build the final PDF 
        # 
        self.pdf  = ROOT.RooAddPdf  ( "model2D"      + suffix ,
                                      "Model2D(%s)"  % suffix ,
                                      self.alist1 ,
                                      self.alist2 )


        self.signals     ().add ( self.__ss_pdf )
        self.backgrounds ().add ( self.__bb_pdf )
        self.crossterms1 ().add ( self.__sb_pdf      ) ## cross-terms 
        self.crossterms2 ().add ( self.__bs_pdf      ) ## cross-terms 


    ## get all declared components 
    def crossterms1 ( self ) : return self.__crossterms1
    ## get all declared components 
    def crossterms2 ( self ) : return self.__crossterms2

    @property
    def S1S2 ( self ) :
        """The yield of Signal(x)*Signal(y) component"""
        return self.__ss
    @S1S2.setter 
    def S1S2 ( self , value ) :
        value = float ( value  )
        assert value in self.__ss, "Value %s is out of the allowed range %s " % ( value , self.__ss.minmax() )
        self.__ss.setVal ( value ) 
        return self.__ss.getValue() 

    @property
    def S1B2 ( self ) :
        """The yield of Signal(x)*Background(y) component"""
        return self.__sb
    @S1B2.setter 
    def S1B2 ( self , value ) :
        value = float ( value  )
        assert value in self.__sb, "Value %s is out of the allowed range %s " % ( value , self.__sb.minmax() )
        self.__sb.setVal ( value ) 
        return self.__sb.getValue() 

    @property
    def B1S2 ( self ) :
        """The yield of Background(x)*Signal(y) component"""
        return self.__bs
    @B1S2.setter 
    def B1S2 ( self , value ) :
        value = float ( value  )
        assert value in self.__bs, "Value %s is out of the allowed range %s " % ( value , self.__bs.minmax() )
        self.__bs.setVal ( value ) 
        return self.__bs.getValue() 

    @property
    def B1B2 ( self ) :
        """The yield of Background(x,y) component"""
        return self.__bb
    @B1B2.setter 
    def B1B2 ( self , value ) :
        value = float ( value  )
        assert value in self.__bb, "Value %s is out of the allowed range %s " % ( value , self.__bb.minmax() )
        self.__bb.setVal ( value ) 
        return self.__bb.getValue() 

    @property 
    def signal1 ( self  ) :
        """Signal(x) component/pdf"""
        return self.__signal1

    @property 
    def signal2 ( self  ) :
        """Signal(y) component/pdf"""
        return self.__signal2

    @property
    def bkg1( self ) :
        """ The background PDF for Backgroud(x)*Signal(y) component"""
        return self.__bkg1
    
    @property
    def bkg2( self ) :
        """ The background PDF for Signal(x)*Background(y) component"""
        return self.__bkg2

    @property
    def bkgA( self ) :
        """ The background(x) PDF for Backgroud(x)*Background(y) component"""
        return self.__bkgA

    @property
    def bkgB( self ) :
        """ The background(y) PDF for Backgroud(x)*Background(y) component"""
        return self.__bkgB
    


    
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
## @class Fit2DSym
#  The actual model for 2D-fits
#
#  @code
# 
#  model   = Models.Fit2D (
#      signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
#      signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y ) ,
#      bkg1     = 1 ) 
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
    """The actual model for SYMMETRIC 2D-fits
    
    >>>  model   = Models.Fit2D (
    ...      signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
    ...      signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y ) ,
    ...      bkg1     = 1 )
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
                   bkgA       = None  ,
                   bkg2D      = None  ,
                   #
                   ## main components :
                   ss         = None  , ## signal (1) * signal     (2)
                   sb         = None  , ## signal (*) background 
                   bb         = None  , ## background-2D 
                   ## additional components 
                   components = []    ,
                   name       = ''    ) :
        
        self._crossterms1 = ROOT.RooArgSet()
        self._crossterms2 = ROOT.RooArgSet()
        
        self.suffix    = suffix 
        self.__signal1   = signal_1
        self.__signal2   = signal_2

        #
        ## initialize base class
        #
        if not name and signal_1.name and signal_2.name :
            name = '%s_and_%s_%s' % ( signal_1.name , signal_2.name , suffix )
            
        PDF2.__init__ ( self , name , signal_1.mass , signal_2.mass ) 
        

        #
        ## First component: Signal(1) and Signal(2)
        # 
        self.__ss_pdf = ROOT.RooProdPdf ( "S1S2pdf" + suffix ,
                                          "Sig(1) x Sig(2)"  ,
                                          self.signal1.pdf   ,
                                          self.signal2.pdf   )
        self.__arg_bkg1 = bkg1 
        self.__bkg1  = makeBkg ( bkg1   , 'Bkg(1)' + suffix , self.xvar )
        
        if bkg1 :
            if hasattr ( self.bkg1, 'tau' )  :
                self.__bkg2 = makeBkg ( bkg1   , 'Bkg(2)' + suffix , self.yvar , the_phis = self.bkg1 , tau = self.bkg1.tau )
            else :
                self.__bkg2 = makeBkg ( bkg1   , 'Bkg(2)' + suffix , self.yvar , the_phis = self.bkg1 )
        else    :
            if hasattr ( self.bkg1, 'tau' )  :
                self.__bkg2 = makeBkg ( bkg1   , 'Bkg(2)' + suffix , self.yvar , the_phis = self.bkg1 , tau = self.bkg1.tau )
            else :
                self.__bkg2 = makeBkg ( bkg1   , 'Bkg(2)' + suffix , self.yvar , the_phis = self.bkg1 )
        
        
        #
        ## Second sub-component: Background (1) and Signal     (2)
        ## Third  sub-component: Signal     (1) and Background (2)
        # 
        self.__sub_bs_pdf = ROOT.RooProdPdf ( "B1S2pdf" + suffix  ,
                                         "Bkg(1) x Sig(2)"   ,
                                         self.bkg1.pdf       ,
                                         self.signal2.pdf    )
        self.__sub_sb_pdf = ROOT.RooProdPdf ( "S1B2pdf" + suffix  ,
                                              "Sig(1) x Bkg(2)"   ,
                                              self.signal1.pdf    ,
                                              self.bkg2.pdf       )
        
        self.__f_cross = ROOT.RooConstVar ( 'SxBfraction'   + suffix  , '(S1B2-vs-B1S2) fraction' , 0.5 )
        # 
        self.__sb_pdf   = ROOT.RooAddPdf ( "SxB_pdf" + suffix  ,
                                           "Sig(1) x Bkg(2) + Bkg(1) x Sig(2)"   ,
                                           self.__sub_sb_pdf  ,
                                           self.__sub_bs_pdf  ,
                                           self.__f_cross ) 
        
        ## just for convinience 
        self.__bs_pdf   = self.__sb_pdf 
        ## 
        self._bkgs = ( bkg1 , bkg1 , bkgA , bkgA ) 

        #
        ## fourth component: Background(1) and Background(2) 
        #
        if bkg2D : self.__bb2D  = bkg2D
        # 
        self.__bkgA = None  
        self.__bkgB = None 
        
        if   bkg2D and isinstance ( bkg2D , ROOT.RooAbsPdf ) : self.__bb_pdf = bkg2D 
        elif bkg2D and hasattr    ( bkg2D , 'pdf'          ) : self.__bb_pdf = bkg2D.pdf
        else     :            

            self.__arg_bkgA = bkgA
            if bkgA is None : bkgA = bkg1
            
            self.__bkgA = makeBkg ( bkgA   , 'Bkg(A)' + suffix , self.xvar )
            self.__bkgB = makeBkg ( bkgA   , 'Bkg(B)' + suffix , self.yvar , the_phis = self.bkgA )
            
            self.__bb_pdf = ROOT.RooProdPdf ( "B1B2pdf" + suffix ,
                                              "Bkg(A) x Bkg(B)"  ,
                                              self.bkgA.pdf      ,
                                              self.bkgB.pdf      )
        #
        ## coefficients
        #
        self.__ss = makeVar ( ss   ,
                              "S1S2"          + suffix ,
                              "Sig(1)&Sig(2)" + suffix , None , 1000  , 0 ,  1.e+7 )
        
        self.__bb = makeVar ( bb   ,
                              "B1B2"          + suffix ,
                              "Bkg(1)&Bkg(2)" + suffix , None ,   10  , 0 ,  1.e+7 )
        
        self.__sb = makeVar ( sb   ,
                              "SxB"           + suffix ,
                              "Sig(1)&Bkg(2)+Bkg(1)&Sig(2)" + suffix , None ,  100  , 0 ,  1.e+7 )
        
        self.__bs = self.__sb
        
        self.__sb_  = ROOT.RooFormulaVar (
            'S1B2' + suffix ,
            'Sig(1)&Bkg(2)' ,
            '0.5*%s' % self.__sb.GetName() , ROOT.RooArgList ( self.__sb ) )
        
        self.__bs_  = ROOT.RooFormulaVar (
            'B1S2' + suffix ,
            'Bkg(1)&Sig(2)' ,
            '0.5*%s' % self.__sb.GetName() , ROOT.RooArgList ( self.__sb ) )        
        
        self.SS_name = self.S1S2.GetName()
        self.BB_name = self.B1B2.GetName()
        self.SB_name = self.S1B2.GetName()
        self.BS_name = self.B1S2.GetName()
        
        self.alist1 = ROOT.RooArgList ( self.__ss_pdf , self.__sb_pdf , self.__bb_pdf )
        self.alist2 = ROOT.RooArgList ( self.__ss     , self.__sb     , self.__bb     )

        #
        ## treat additional components (if specified)
        # 
        self.other       = components
        self._cmps       = []
        icmp = 0 
        for cmp in self.other :

            icmp += 1

            if   isinstance ( cmp , ROOT.RooAbsPdf         ) : pass 
            elif hasattr    ( cmp , 'pdf'                  ) :
                self._cmps += [ cmp ] 
                cmp = cmp.pdf 
            elif isinstance ( cmp , ( float , int , long ) ) and not isinstance ( cmp , bool ) :
                px  = ROOT.RooPolynomial ( 'Px%d'    % icmp + suffix ,
                                           'Px(%d)'  % icmp + suffix , self.xvar ) 
                py  = ROOT.RooPolynomial ( 'Py%d'    % icmp + suffix ,
                                           'Py(%d)'  % icmp + suffix , self.yvar ) 
                cmp = ROOT.RooProdPdf    ( "Pxy%d"   % icmp + suffix ,
                                           "Pxy(%d)" % icmp + suffix , px , py )  
                self._cmps += [ px,py,cmp]
            else :
                logger.error( 'Unknown type of component %d %d ' % ( imcp , type(cmp) ) )

                
            nn = makeVar ( None ,
                           'Cmp%d'   % icmp + suffix ,
                           'Cmp(%d)' % icmp + suffix ,
                           None ,  100  ,  0 , 1.e+7 )  
            self._cmps += [ nn ]

            self.alist1.add ( cmp )
            self.alist2.add ( nn  )
            
            self.components ().add ( cmp ) 
            
        #
        ## build the final PDF 
        # 
        self.pdf  = ROOT.RooAddPdf  ( "model2D"      + suffix ,
                                      "Model2D(%s)"  % suffix ,
                                      self.alist1 ,
                                      self.alist2 )


        self.signals     ().add ( self.__ss_pdf )
        self.backgrounds ().add ( self.__bb_pdf )
        self.crossterms1 ().add ( self.__sb_pdf ) ## cross-terms 
        self.crossterms2 ().add ( self.__bs_pdf ) ## cross-terms 


        

    ## get all declared components 
    def crossterms1 ( self ) : return self._crossterms1
    ## get all declared components 
    def crossterms2 ( self ) : return self._crossterms2

    @property
    def S1S2 ( self ) :
        """The yield of Signal(x)*Signal(y) component"""
        return self.__ss
    @S1S2.setter 
    def S1S2 ( self , value ) :
        value = float ( value  )
        assert value in self.__ss, "Value %s is out of the allowed range %s " % ( value , self.__ss.minmax() )
        self.__ss.setVal ( value ) 
        return self.__ss.getValue() 

    @property
    def S1B2 ( self ) :
        """The yield of Signal(x)*Background(y)+BAckgroun(1)*Signal(2) component (the same as ``S2B1'')"""
        return self.__sb
    @S1B2.setter 
    def S1B2 ( self , value ) :
        value = float ( value  )
        assert value in self.__sb, "Value %s is out of the allowed range %s " % ( value , self.__sb.minmax() )
        self.__sb.setVal ( value ) 
        return self.__sb.getValue() 
    
    @property
    def B1S2 ( self ) :
        """The yield of Signal(x)*Background(y)+BAckgroun(1)*Signal(2) component (the same as ``S1B2'')"""
        return self.__bs
    @B1S2.setter 
    def B1S2 ( self , value ) :
        value = float ( value  )
        assert value in self.__bs, "Value %s is out of the allowed range %s " % ( value , self.__bs.minmax() )
        self.__bs.setVal ( value ) 
        return self.__bs.getValue() 

    @property
    def B1B2 ( self ) :
        """The yield of Background(x,y) component"""
        return self.__bb
    @B1B2.setter 
    def B1B2 ( self , value ) :
        value = float ( value  )
        assert value in self.__bb, "Value %s is out of the allowed range %s " % ( value , self.__bb.minmax() )
        self.__bb.setVal ( value ) 
        return self.__bb.getValue() 

    @property 
    def signal1 ( self  ) :
        """Signal(x) component/pdf"""
        return self.__signal1

    @property 
    def signal2 ( self  ) :
        """Signal(y) component/pdf"""
        return self.__signal2

    @property
    def bkg1( self ) :
        """ The background PDF for Backgroud(x)*Signal(y) component"""
        return self.__bkg1
    
    @property
    def bkg2( self ) :
        """ The background PDF for Signal(x)*Background(y) component"""
        return self.__bkg2

    @property
    def bkgA( self ) :
        """ The background(x) PDF for Backgroud(x)*Background(y) component"""
        return self.__bkgA

    @property
    def bkgB( self ) :
        """ The background(y) PDF for Backgroud(x)*Background(y) component"""
        return self.__bkgB
    

# =============================================================================
## @class Generic2D_pdf
#  "Wrapper" over generic RooFit (2D)-pdf
#  @code
#     
#  raw_pdf = 
#  pdf     = Generic1D_pdf ( raw_pdf )  
# 
#  @endcode 
#  If more functionality is required , more actions are possible:
#  @code
#  ## for sPlot 
#  pdf.alist2 = ROOT.RooArgList ( n1 , n2 , n3 ) ## for sPlotting 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-03-29
class Generic2D_pdf(PDF2) :
    """Wrapper for generic RooFit pdf
    
    >>> raw_pdf = 
    >>> pdf     = Generic2D_pdf ( raw_pdf )
    """
    ## constructor 
    def __init__ ( self , pdf , varx , vary , name = None ) :
        if not name : name = pdf.GetName()
        PDF2  . __init__ ( self , name , varx , vary )
        self.pdf = pdf

            
# =============================================================================
## simple convertor of 2D-histo to data set
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H2D_dset(object) :
    """Simple convertor of 2D-histogram into data set
    """
    def __init__ ( self            ,
                   histo2          ,
                   mass    = None  ,
                   mass2   = None  ,
                   density = True  ,
                   silent  = False ) :
        #
        ## use mass-variable
        #
        name         = histo3.GetName() 
        self.mass    = makeVar ( mass  , 'm_%s'  % name , 'mass (%s)' % name , None , *(histo3.xminmax()) )
        self.mass1   = self.mass 
        self.mass2   = makeVar ( mass2 , 'm2_%s' % name , 'mass2(%s)' % name , None , *(histo3.yminmax()) )

        self.impDens = density 
        self.var1    = self.mass1
        self.var2    = self.mass2
        self.x       = self.var1 
        self.y       = self.var2
        
        with roo_silent ( silent ) : 

            self.vlst  = ROOT.RooArgList    ( self.mass1 , self.mass2 )
            self.vimp  = ROOT.RooFit.Import ( histo2 , density )
            self.dset  = ROOT.RooDataHist   (
                rootID ( 'hds_' ) ,
                "Data set for histogram '%s'" % histo2.GetTitle() ,
                self.vlst  ,
                self.vimp  )


# =============================================================================
## simple convertor of 2D-histogram into PDF
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H2D_pdf(H2D_dset,PDF2) :
    """Simple convertor of 2D-histogram into PDF 
    """
    def __init__ ( self            ,
                   name            ,
                   histo2          ,
                   mass    = None  , 
                   mass2   = None  ,
                   density = True  ,
                   silent  = False ) :
        
        H2D_dset.__init__ ( self , histo2 ,      mass  ,      mass2 , density , silent )
        PDF2    .__init__ ( self , name   , self.mass1 , self.mass2 ) 

        self.vset  = ROOT.RooArgSet  ( self.mass , self.mass2 )
        
        #
        ## finally create PDF :
        #
        with roo_silent ( silent ) : 
            self.pdf    = ROOT.RooHistPdf (
                'hpdf_%s'            % name ,
                'HistoPDF(%s/%s/%s)' % ( name , histo2.GetName() , histo2.GetTitle() ) , 
                self.vset  , 
                self.dset  )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
