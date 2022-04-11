#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/fit2d.py
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
    'Generic2D_pdf' , ## wrapper over imported RooFit (2D)-pdf
    'Flat2D'        , ## simplest 2D-model: constant function
    'Model2D'       , ## helper class to construct 2D-models. 
    'Combine2D'     , ## non-extended sum of several PDFs
    'Sum2D'         , ## non-extended sum of two PDFs
    'H2D_pdf'       , ## convertor of 1D-histo to RooDataPdf
    'Shape2D_pdf'   , ## simple PDF from C++ shape     
    ## 
    )
# =============================================================================
import ROOT, random 
from   builtins               import range
from   ostap.core.meta_info   import root_info 
from   ostap.core.core        import dsID , VE , Ostap, hID , iszero, valid_pointer
from   ostap.core.ostap_types import integer_types, num_types, list_types, iterable_types   
from   ostap.fitting.roofit   import SETVAR
from   ostap.logger.utils     import roo_silent, rooSilent, rootWarning 
from   ostap.fitting.basic    import PDF , Flat1D
from   ostap.fitting.funbasic import FUNC2 
from   ostap.fitting.utils    import H2D_dset , component_similar , component_clone 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.fit2d' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
# @class PDF2
# The helper base class for implementation of 2D-pdfs 
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date 2014-08-21
class PDF2 (PDF,FUNC2) :
    """ Useful helper base class for implementation of PDFs for 2D-fit
    """
    def __init__ ( self , name , xvar , yvar , special = False , **kwargs ) : 
        
        PDF  .__init__ ( self ,      name ,      xvar , special = special , **kwargs )
        FUNC2.__init__ ( self , self.name , self.xvar , yvar  )
                    
        ## save the configuration
        self.config = {
            'name' : self.name ,
            'xvar' : self.xvar ,
            'yvar' : self.yvar ,            
            }
        
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
                refit  = False ,
                timer  = False ,
                args   = ()    , **kwargs ) :
        """
        Perform the actual fit (and draw it)
        >>> r,f = model.fitTo ( dataset )
        >>> r,f = model.fitTo ( dataset , weighted = True )    
        >>> r,f = model.fitTo ( dataset , ncpu     = 10   )    
        >>> r,f = model.fitTo ( dataset , draw = True , nbins = 300 )    
        """
        if   isinstance ( dataset , H2D_dset ) : dataset = dataset.dset        
        elif isinstance ( dataset , ROOT.TH2 ) :
            density = kwargs.pop ( 'density' , False ) 
            chi2    = kwargs.pop ( 'chi2'    , False ) 
            return self.fitHisto ( dataset   ,
                                   draw    = draw    ,
                                   silent  = silent  ,
                                   density = density ,
                                   chi2    = chi2    , args = args , **kwargs )
        
        ## play a bit with binning cache for convolutions 
        if self.yvar.hasBinning ( 'cache' ) :
            nb1 = self.yvar.getBins( 'cache' ) 
            yv  = getattr ( dataset , self.yvar.name , None )
            if   yv and yv.hasBinning ( 'cache' ) :
                nb2 = yv.getBins('cache')
                if  nb1 != nb2 :
                    yv.setBins ( max (  nb1 , nb2 ) , 'cache' )
                    self.info ('Adjust binning cache %s->%s for variable %s in dataset' % ( nb2 , nb1 , yv.name ) )
            elif yv :
                yv.setBins (        nb1         , 'cache' )
                self    .info ('Set binning cache %s for variable %s in dataset' %  ( nb1 , yv.name )  )
                                
        result , f = PDF.fitTo ( self            ,
                                 dataset         ,
                                 draw   = False  , ## false here!
                                 nbins  = nbins  ,
                                 silent = silent ,
                                 refit  = refit  ,
                                 timer  = timer  , 
                                 args   = args   , **kwargs ) 
        if not draw :
            return result , None

        
        ## 2D 
        if 1 < nbins and isinstance ( ybins , integer_types ) and 1 < ybins :
            return result, self.draw ( None , dataset , nbins , ybins , silent = silent )
        
        if isinstance ( draw , str ) :
            if   draw.upper() in ( '1' , 'X' ) :
                return result, self.draw1 ( dataset , nbins = nbins , silent = silent )
            elif draw.upper() in ( '2' , 'Y' ) :
                return result, self.draw2 ( dataset , nbins = nbins , silent = silent )

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
                in_range = None ,
                args     = ()   , **kwargs ) :
        """ Draw the projection over 1st variable
        
        >>> r,f = model.fitTo ( dataset ) ## fit dataset
        >>> fx  = model.draw1 ( dataset , nbins = 100 ) ## draw results
        
        >>> f1  = model.draw1 ( dataset , nbins = 100 , in_range = (2,3) ) ## draw results

        >>> model.yvar.setRange ( 'QUQU2' , 2 , 3 ) 
        >>> f1  = model.draw1 ( dataset , nbins = 100 , in_range = 'QUQU2') ## draw results
        
        """
        if in_range and isinstance ( in_range , tuple ) and 2 == len ( in_range ) :
            range_name = 'aux2_rng2_%s' % self.name 
            with rooSilent ( 3 ) : 
              self.yvar.setRange ( range_name , in_range[0] , in_range[1] )
              if dataset:
                dataset.get_var(self.yvar.GetName()).setRange ( range_name , in_range[0] , in_range[1] )

            in_range = range_name 

        return self.draw ( drawvar  = self.xvar , 
                           dataset  = dataset   ,
                           nbins    = nbins     ,
                           ybins    = 20        , ## fake 
                           silent   = silent    ,
                           in_range = in_range  ,
                           args     = args      , **kwargs )
    
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
                in_range = None ,
                args     = ()   , **kwargs ) :
        """
        Draw the projection over 2nd variable
        
        >>> r,f = model.fitTo ( dataset ) ## fit dataset
        >>> fy  = model.draw2 ( dataset , nbins = 100 ) ## draw results
        
        >>> f2  = model.draw2 ( dataset , nbins = 100 , in_range = (2,3) ) ## draw results

        >>> model.xvar.setRange ( 'QUQU1' , 2 , 3 ) 
        >>> f2  = model.draw2 ( dataset , nbins = 100 , in_range = 'QUQU1') ## draw results

        """
        if in_range and isinstance ( in_range , tuple ) and 2 == len ( in_range ) :
            range_name = 'aux2_rng1_%s' % self.name 
            with rooSilent ( 3 ) : 
              self.xvar.setRange ( range_name , in_range[0] , in_range[1] )
              if dataset:
                dataset.get_var(self.xvar.GetName()).setRange ( range_name , in_range[0] , in_range[1] )

            in_range = range_name

        return self.draw ( drawvar  = self.yvar ,
                           dataset  = dataset   ,
                           nbins    = nbins     ,
                           ybins    = 20        , ## fake
                           silent   = silent    ,
                           in_range = in_range  ,
                           args     = args      , **kwargs )

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

        groot = ROOT.ROOT.GetROOT()
        if not groot.IsBatch() :
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
               args                  = ()   , 
               **kwargs                     ) : 
        """
        Make 1D-plot:
        """
        if   drawvar in ( 'x'  , 'X' , '1' , 1 , self.xvar.name ) : drawvar = self.xvar
        elif drawvar in ( 'y'  , 'Y' , '2' , 2 , self.yvar.name ) : drawvar = self.yvar

        # 
        ## special case:  do we need it? 
        #

        if drawvar is None : return self.draw_H2D( dataset , nbins , ybins )
                
        newargs = kwargs.copy ()
        
        if in_range and isinstance ( in_range , list_types ) and 2 == len ( in_range ) :
            low  = in_range [ 0 ]
            high = in_range [ 1 ]
            if isinstance ( low , num_types ) and isinstance ( high , num_types ) and low < high :
                range_name = 'aux2_range_%s' % self.name 
                with rooSilent ( 3 ) : 
                  drawvar.setRange ( range_name , low , high )
                  dataset.get_var(drawvar.GetName()).setRange ( range_name , low , high )
                in_range = range_name
    
  #      if in_range and not isinstance ( in_range , list_types ) :
   #         in_range = in_range ,
        
        if in_range :
            options_cut =tuple ( [ ROOT.RooFit.CutRange ( in_range ) , ] )
            newargs [ 'data_options' ] = self.draw_option ( 'data_options' , **newargs ) + options_cut
            
        if in_range : 
            options_project =  tuple (  [ROOT.RooFit.ProjectionRange ( in_range ) ,] )
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
                newargs [ key ] =  self.draw_option ( key , **newargs ) + options_project
                
        #
        ## redefine the drawing variable:
        # 
        self.draw_var = drawvar
        ##
        if in_range and root_info < ( 6 , 24 ) :
            from itertools import chain 
            for p in chain ( self.signals               ,
                             self.backgrounds           ,
                             self.components            ,
                             self.crossterms1           ,
                             self.crossterms2           ,
                             self.combined_signals      ,
                             self.combined_backgrounds  ,
                             self.combined_components   ) :
                
                if   isinstance ( p , ( ROOT.RooHistPdf , ROOT.RooParamHistFunc ) ) :
                    self.warning ("``in_range'' is specified, it does not work properly for ROOT<6.24")
                elif isinstance ( p , ( ROOT.RooAddPdf , ROOT.RooProdPdf ) ) : 
                    for pp in p.pdfList() :
                        if isinstance  ( pp , ( ROOT.RooHistPdf , ROOT.RooParamHistFunc ) ) :
                            self.warning ("``in_range'' is specified, it does not work properly for ROOT<6.24")
                            
        #
        ## delegate the actual drawing to the base class
        #
        result = PDF.draw ( self            ,
                            dataset         ,
                            nbins  = nbins  ,
                            silent = silent ,
                            args   = args   , **newargs )

        self.draw_var = None
        return result 
    
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
                   density = False ,
                   chi2    = False ,
                   args    = ()    , **kwargs ) :
        """Fit the histogram (and draw it)
        
        >>> histo = ...
        >>> r,f = model.fitHisto ( histo , draw = True )
        
        """

        xminmax = histo.xminmax()
        yminmax = histo.yminmax()        
        with RangeVar( self.xvar , *xminmax ) , RangeVar ( self.yvar , *yminmax ):
            
            hdata = getattr ( self , 'histo_data' , None )
            if hdata and isinstance ( hdata  , H2D_dset ) and \
                   hdata.histo      is histo              and \
                   hdata.density    == density            and \
                   hdata.histo_hash == hash ( histo ) :
                ## reuse the existing dataset
                self.debug ('Reuse the existing H2D_dset') 
                data = hdata.dset
            else :                
                ## convert it!
                self.debug ('Create new H2D_dset'        ) 
                self.histo_data = H2D_dset ( histo , self.xvar , self.yvar  , density , silent )
                data = self.histo_data.dset 
            
            ## fit it!!
            if chi2 : return self.chi2fitTo ( data                     ,
                                              draw    = draw           ,
                                              silent  = False          ,
                                              density = density        ,
                                              args    = args           , **kwargs )
            else     : return self.fitTo    ( data                     ,
                                              draw    = draw           ,
                                              nbins   = histo.nbinsx() ,
                                              ybins   = histo.nbinsy() ,
                                              silent  = silent         ,
                                              args    = args           , **kwargs )
            
    # =========================================================================
    ## generate toy-sample according to PDF
    #  @code
    #  model  = ....
    #  data   = model.generate ( 10000 ) ## generate dataset
    #  varset = ....
    #  data   = model.generate ( 100000 , varset , sample = False )
    #  data   = model.generate ( 100000 , varset , sample = True  )     
    #  @endcode
    def generate ( self             ,  
                   nEvents          ,
                   varset   = None  ,
                   binning  = {}    ,
                   sample   = True  , 
                   args     = ()    ) :
        """Generate toy-sample according to PDF
        >>> model  = ....
        >>> data   = model.generate ( 10000 ) ## generate dataset
        
        >>> varset = ....
        >>> data   = model.generate ( 100000 , varset , sample = False )
        >>> data   = model.generate ( 100000 , varset , sample = True  )
        """
        nEvents = self.gen_sample ( nEvents ) if sample else nEvents 
        assert 0 <= nEvents , 'Invalid number of Events %s' % nEvents  
        
        args = args + ( ROOT.RooFit.Name ( dsID() ) , ROOT.RooFit.NumEvents ( nEvents ) )
        
        if binning is True :
            args    = args + ( ROOT.AllBinned() , ) 
            binning = {}

        if   not varset :
            varset = ROOT.RooArgSet( self.xvar , self.yvar )
        elif isinstance ( varset , ROOT.RooAbsReal ) :
            varset = ROOT.RooArgSet( varset )

        if not self.xvar in varset :
            vs = ROOT.RooArgSet()
            vs . add ( self.xvar )
            for  v in varset : vs.add ( v )
            varset = vs

        if not self.yvar in varset :
            vs = ROOT.RooArgSet()
            vs . add ( self.yvar )
            for  v in varset : vs.add ( v )
            varset = vs

        from ostap.fitting.variables import KeepBinning        
        with KeepBinning ( self.xvar ) , KeepBinning ( self.yvar ) : 

            if binning :
                
                xbins = binning.get ( self.xvar.name , None )
                ybins = binning.get ( self.yvar.name , None )

                if xbins : self.xvar.bins = xbins
                if ybins : self.yvar.bins = ybins
                                        
            return self.pdf.generate ( varset , *args )

    # =========================================================================
    ## simple 'function-like' interface 
    def __call__ ( self , x , y , error = False , normalized = True ) :
        """ Simple  function-like interface
        >>>  pdf = ...
        >>>  print pdf(0.1,0.5) 
        """
        if     isinstance ( self.xvar , ROOT.RooRealVar ) and \
               isinstance ( self.yvar , ROOT.RooRealVar ) :
            
            if x in self.xvar and y in self.yvar : 
                with SETVAR ( self.xvar ) , SETVAR( self.yvar ) :
                    self.xvar.setVal ( x )
                    self.yvar.setVal ( y )
                    
                    v = self.pdf.getVal ( self.vars ) if normalized else self.pdf.getValV ()
                    
                    if error and self.fit_result :
                        e = self.pdf.getPropagatedError ( self.fit_result )
                        if 0<= e : return  VE ( v ,  e * e )
                    return v 
            else :
                return 0.0
            
        raise AttributeError('something wrong goes here')


    # ========================================================================
    ## check minmax of the PDF using the random shoots
    #  @code
    #  pdf     = ....
    #  mn , mx = pdf.minmax()            
    #  @endcode 
    def minmax ( self , nshoots =  100000 ) :
        """Check min/max for the PDF using  random shoots 
        >>> pdf     = ....
        >>> mn , mx = pdf.minmax()        
        """
        ## try to get minmax directly from pdf/function 
        if self.tricks and hasattr ( self.pdf , 'function' ) :
            if hasattr ( self.pdf , 'setPars' ) : self.pdf.setPars() 
            f = self.pdf.function()
            if hasattr ( f , 'minmax' ) :
                try :
                    mn , mx = f.minmax()
                    if  0<= mn and mn <= mx and 0 < mx :   
                        return mn , mx
                except :
                    pass
            if hasattr ( f , 'max' ) :
                try :
                    mx = f.max()
                    if 0 < mx : return 0 , mx
                except :
                    pass

        ## check RooAbsReal functionality
        code = self.pdf.getMaxVal( ROOT.RooArgSet ( self.xvar , self.yvar ) )
        if 0 < code :
            mx = self.pdf.maxVal ( code )
            if 0 < mx : return 0 , mx
            
        ## not try  to use random
                
        mn , mx = -1 , -10
        if hasattr ( self.pdf , 'min' ) : mn = self.pdf.min()
        if hasattr ( self.pdf , 'max' ) : mx = self.pdf.max()
        if 0 <= mn and mn <= mx and 0 < mx : return mn , mx
        
        if not self.xminmax() : return ()
        if not self.yminmax() : return ()
        
        mn  , mx = -1 , -10
        xmn , xmx = self.xminmax()
        ymn , ymx = self.yminmax()
        for i in range ( nshoots ) : 
            xx = random.uniform ( xmn , xmx )
            yy = random.uniform ( ymn , ymx )
            with SETVAR ( self.xvar ) :
                with SETVAR ( self.yvar ) :
                    self.xvar.setVal ( xx )
                    self.yvar.setVal ( yy )
                    vv = self.pdf.getVal()
                    if mn < 0 or vv < mn : mn = vv
                    if mx < 0 or vv > mx : mx = vv
                    
        return mn , mx 
        

    # =========================================================================
    ## get integral over (xmin,xmax,ymin,ymax) region
    #  @code
    #  pdf = ...
    #  print pdf.integral( 0,1,0,2)
    #  @endcode
    def integral ( self, xmin , xmax , ymin , ymax , nevents = True ) :
        """Get integral over (xmin,xmax,ymin,ymax) region
        >>> pdf = ...
        >>> print pdf.integral( 0,1,0,2)
        """
        if self.xminmax() :
            xmn , xmx = self.xminmax()
            xmin = max ( xmin , xmn )
            xmax = min ( xmax , xmx )

        if self.yminmax() : 
            ymn , ymx = self.yminmax()            
            ymin = max ( ymin , ymn )
            ymax = min ( ymax , ymx )

        value , todo  = 0 , True 
        
        ## 1) make a try to use analytical integral (could be fast)
        if self.tricks :
            try:
                if hasattr ( self.pdf , 'setPars'  ) : self.pdf.setPars() 
                fun          = self.pdf.function()
                value , todo = fun.integral ( xmin , xmax , ymin , ymax ) , False 
            except:
                pass

        ## use numerical integration 
        from ostap.math.integral import integral2 as _integral2

        extended =  self.pdf.canBeExtended() or isinstance ( self.pdf , ROOT.RooAddPdf )

        if   todo and extended : value   = _integral2 ( self , xmin , xmax , ymin , ymax )
        elif todo  :
            
            ## use unormalized PDF here to speed up the integration 
            ifun   = lambda x, y  :  self ( x , y , error = False , normalized = False )
            value  = _integral2 ( ifun , xmin , xmax , ymin , ymax )
            norm   = self.pdf.getNorm ( self.vars )
            value /= norm

        if nevents and self.pdf.mustBeExtended () :
            evts = self.pdf.expectedEvents( self.vars )
            if evts  <= 0 or iszero ( evts ) :
                self.warning ( "integral: expectedEvents is %s" % evts )
            value *= evts 

        return value


    # ==========================================================================
    ## get a minimum of PDF for certain interval
    #  @code
    #  pdf2 = ...
    #  x ,y = pdf2.minimum() 
    #  @endcode 
    def minimum ( self ,
                  xmin = None , xmax = None ,
                  ymin = None , ymax = None , x0 = () ) :
        """Get a minimum of PDF for certain interval
        >>> pdf2 = ...
        >>> x, y = pdf2.minimum()
        """
        
        if xmin is None : xmin = self.xminmax()[0]
        if xmax is None : xmax = self.xminmax()[1]
        if self.xminmax() :
            xmin =  max ( xmin , self.xminmax()[0] )
            xmax =  min ( xmax , self.xminmax()[1] )

        if ymin is None : ymin = self.yminmax()[0]
        if ymax is None : ymax = self.yminmax()[1]
        if self.yminmax() :
            ymin =  max ( ymin , self.yminmax()[0] )
            ymax =  min ( ymax , self.yminmax()[1] )
            
        if not x0 : x0 = 0.5 * ( xmin + xmax ) , 0.5 * ( ymin + ymax )
        
        if not xmin <= x0[0] <= xmax :
            self.error("Wrong xmin/x0[0]/xmax: %s/%s/%s"   % ( xmin , x0[0] , xmax ) )

        if not ymin <= x0[1] <= ymax : 
            self.error("Wrong ymin/x0[1]/ymax: %s/%s/%s"   % ( ymin , x0[1] , ymax ) )
        
        from ostap.math.minimize import sp_minimum_2D
        return sp_minimum_2D (  self ,
                                xmin , xmax ,
                                ymin , ymax , x0 )

    # ==========================================================================
    ## get a maximum of PDF for certain interval
    #  @code
    #  pdf2 = ...
    #  x,y  = pdf2.maximum() 
    #  @endcode 
    def maximum ( self , xmin = None , xmax = None , x0 = None ) :
        """Get a maximum of PDF for certain interval
        >>> pdf2  = ...
        >>> x , y  = pdf2.maximum()
        """
        if xmin is None : xmin = self.xminmax()[0]
        if xmax is None : xmax = self.xminmax()[1]
        if self.xminmax() :
            xmin =  max ( xmin , self.xminmax()[0] )
            xmax =  min ( xmax , self.xminmax()[1] )

        if ymin is None : ymin = self.yminmax()[0]
        if ymax is None : ymax = self.yminmax()[1]
        if self.yminmax() :
            ymin =  max ( ymin , self.yminmax()[0] )
            ymax =  min ( ymax , self.yminmax()[1] )
            
        if not x0 : x0 = 0.5 * ( xmin + xmax ) , 0.5 * ( ymin + ymax )

        if not xmin <= x0[0] <= xmax :
            self.error("Wrong xmin/x0[0]/xmax: %s/%s/%s"   % ( xmin , x0[0] , xmax ) )

        if not ymin <= x0[1] <= ymax : 
            self.error("Wrong ymin/x0[1]/ymax: %s/%s/%s"   % ( ymin , x0[1] , ymax ) )

        from ostap.math.minimize import sp_maximum_2D
        return sp_maximum_2D (  self ,
                                xmin , xmax ,
                                ymin , ymax , x0 )
    
    # ==========================================================================
    ## convert PDF into TF2 object, e.g. to profit from TF2::Draw options
    #  @code
    #  pdf = ...
    #  tf2 = pdf.tf()
    #  tf2.Draw('colz')
    #  @endcode
    def tf ( self , xmin = None , xmax = None , ymin = None , ymax = None ) :
        """Convert PDF to TF2 object, e.g. to profit from TF2::Draw options
        >>> pdf = ...
        >>> tf2 = pdf.tf()
        >>> tf1.Draw('colz')
        """
        def _aux_fun_ ( x , pars = [] ) :
            return self ( x[0] , x[1] , error = False )

        if xmin == None and self.xminmax() : xmin = self.xminmax()[0]
        if xmax == None and self.xminmax() : xmax = self.xminmax()[1]
        if ymin == None and self.yminmax() : ymin = self.yminmax()[0]
        if ymax == None and self.yminmax() : ymax = self.yminmax()[1]
        
        if xmin == None : xmin = 0.0
        if xmax == None : xmin = 1.0
        if ymin == None : ymin = 0.0
        if ymax == None : ymin = 1.0
        
        from ostap.core.core import fID
        return ROOT.TF2 ( fID() , _aux_fun_ , xmin , xmax , ymin , ymax ) 

    
    # ==========================================================================
    ## Create the histo according to specifications 
    def make_histo ( self , 
                     xbins    = 20    , xmin = None , xmax = None ,
                     ybins    = 20    , ymin = None , ymax = None ,
                     hpars    = ()    , 
                     histo    = None  ) :
        """Create the histogram accordig to specifications
        """
        
        import ostap.histos.histos

        # histogram is provided 
        if histo :
            
            assert isinstance ( histo , ROOT.TH2 ) and not isinstance ( histo , ROOT.TH3 ) , \
                   "Illegal type of ``histo''-argument %s" % type( histo )
            
            histo = histo.clone()
            histo.Reset()

        # arguments for the histogram constructor 
        elif hpars :
            
            histo = ROOT.TH2F ( hID () , 'PDF%s' % self.name , *hpars  )
            if not histo.GetSumw2() : histo.Sumw2()

        # explicit construction from (#bins,min,max)-triplet  
        else :
            
            assert isinstance ( xbins , integer_types ) and 0 < xbins, \
                   "Wrong ``xbins''-argument %s" % xbins 
            assert isinstance ( ybins , integer_types ) and 0 < ybins, \
                   "Wrong ``ybins''-argument %s" % ybins 
            if xmin == None and self.xminmax() : xmin = self.xminmax()[0]
            if xmax == None and self.xminmax() : xmax = self.xminmax()[1]
            if ymin == None and self.yminmax() : ymin = self.yminmax()[0]
            if ymax == None and self.yminmax() : ymax = self.yminmax()[1]
            
            histo = ROOT.TH2F ( hID() , 'PDF%s' % self.name ,
                                xbins , xmin , xmax ,
                                ybins , ymin , ymax )
            if not histo.GetSumw2() : histo.Sumw2()

        return histo 
                     
    # ==========================================================================
    ## Convert PDF to the 2D-histogram
    #  @code
    #  pdf = ...
    #  h1  = pdf.histo ( 100 , 0. , 10. , 20 , 0. , 10 ) ## specify histogram parameters
    #  histo_template = ...
    #  h2  = pdf.histo ( histo = histo_template ) ## use historgam template
    #  h3  = pdf.histo ( ... , integral = True  ) ## use PDF integral within the bin  
    #  h4  = pdf.histo ( ... , density  = True  ) ## convert to "density" histogram 
    #  @endcode
    def histo ( self             ,
                xbins    = 20    , xmin = None , xmax = None ,
                ybins    = 20    , ymin = None , ymax = None ,
                hpars    = ()    , 
                histo    = None  ,
                integral = False ,
                errors   = False , 
                density  = False ) :
        """Convert PDF to the 2D-histogram
        >>> pdf = ...
        >>> h1  = pdf.histo ( 100 , 0. , 10. , 20 , 0. , 10 ) ## specify histogram parameters
        >>> histo_template = ...
        >>> h2  = pdf.histo ( histo = histo_template ) ## use historgam template
        >>> h3  = pdf.histo ( ... , integral = True  ) ## use PDF integral within the bin  
        >>> h4  = pdf.histo ( ... , density  = True  ) ## convert to 'density' histogram 
        """
        
        
        histos = self.make_histo ( xbins = xbins , xmin = xmin , xmax = xmax ,
                                   ybins = ybins , ymin = ymin , ymax = ymax ,
                                   hpars = hpars ,
                                   histo = histo )

        # loop over the historgam bins 
        for ix,iy,x,y,z in histo.items() :

            xv , xe = x.value() , x.error()
            yv , ye = y.value() , y.error()
            
            # value at the bin center 
            c = self ( xv , yv , error = errors ) 

            if not integral : 
                histo[ix,iy] = c
                continue

            # integral over the bin 
            v  = self.integral( xv - xe , xv + xe , yv - ye , yv + ye )
            
            if errors :
                if    0 == c.cov2 () : pass
                elif  0 != c.value() and 0 != v : 
                    v = c * ( v / c.value() )
                    
            histo[ix,iy] = v 

        ## coovert to density historgam, if requested 
        if density : histo =  histo.density()
        
        return histo


    # ==========================================================================
    ## Convert PDF to the 2D-histogram, taking taking PDF-values at bin-centres
    #  @code
    #  pdf = ...
    #  h1  = pdf.roo_histo ( 100 , 0. , 10. , 20 , 0. , 10 ) 
    #  histo_template = ...
    #  h2  = pdf.roo_histo ( histo = histo_template ) ## use historgam template
    #  h3  = pdf.roo_histo ( ... , density  = True  ) ## convert to "density" histogram 
    #  @endcode
    def roo_histo ( self             ,
                   xbins   = 20    , xmin = None , xmax = None ,
                   ybins   = 20    , ymin = None , ymax = None ,
                   hpars   = ()    , 
                   histo   = None  , 
                   events  = True  ) : 
        """Convert PDF to the 2D-histogram, taking PDF-values at bin-centres
        >>> pdf = ...
        >>> h1  = pdf.as_histo ( 100 , 0. , 10. , 20 , 0. , 10 ) 
        >>> histo_template = ...
        >>> h2  = pdf.as_histo ( histo = histo_template ) ## use historgam template
        >>> h3  = pdf.as_histo ( ... , density  = True  ) ## convert to 'density' histogram 
        """
        
        histo = self.make_histo ( xbins = xbins , xmin = xmin , xmax = xmax ,
                                  ybins = ybins , ymin = ymin , ymax = ymax ,
                                  hpars = hpars ,
                                  histo = histo )
        
        hh = self.pdf.createHistogram (
            hID()     ,
            self.xvar ,                    self.binning ( histo.GetXaxis() , 'histo2x' )   ,
            ROOT.RooFit.YVar ( self.yvar , self.binning ( histo.GetYaxis() , 'histo2y' ) ) , 
            ROOT.RooFit.Scaling  ( False ) , 
            ROOT.RooFit.Extended ( False ) ) 

        for i in hh : hh.SetBinError ( i , 0 ) 
        
        if events and self.pdf.mustBeExtended() :
            
            for ix , iy , x , y , z in hh.items() :
                volume          = 4 * x.error() * y.error() 
                hh [ ix , iy ] *= volume
                
            hh *= self.pdf.expectedEvents ( self.vars ) / hh.sum() 
                
        histo += hh
        
        return histo 

    
    # ==========================================================================
    ## get the residual histogram : (data-fit) 
    #  @see PDF.as_histo
    #  @see PDF.residual_histo
    #  @see PDF.make_histo
    #  @code
    #  data = ...
    #  pdf  = ...
    #  pdf.fitTo ( data )
    #  residual = pdf.residual ( data , nbins = 100 ) 
    #  @endcode 
    def residual ( self  , dataset , **kwargs ) :
        """Get the residual histogram
        - see PDF.as_histo
        - see PDF.residual_histo
        - see PDF.make_histo

        >>> data = ...
        >>> pdf  = ...
        >>> pdf.fitTo ( data )
        >>> residual = pdf.residual ( data , nbins = 100 ) 
        """
        hdata = self.make_histo ( **kwargs )
        dataset.project ( hdata , ( self.yvar.name , self.xvar.name )  )
        return self.residual_histo ( hdata ) 
        
    # ==========================================================================
    ## get the pull histogram : (data-fit)/data_error 
    #  @see PDF.as_histo
    #  @see PDF.residual_histo
    #  @see PDF.make_histo
    #  @code
    #  data = ...
    #  pdf  = ...
    #  pdf.fitTo ( data )
    #  residual = pdf.pull ( data , nbins = 100 ) 
    #  @endcode 
    def pull ( self  , dataset , **kwargs ) :
        """Get the pull  histogram: (data-fit)/data_error
        - see PDF.as_histo
        - see PDF.residual_histo
        - see PDF.make_histo

        >>> data = ...
        >>> pdf  = ...
        >>> pdf.fitTo ( data )
        >>> residual = pdf.residual ( data , nbins = 100 ) 
        """
        hdata = self.make_histo ( **kwargs )
        dataset.project ( hdata , ( self.yvar.name , self.xvar.name ) ) 
        return self.pull_histo ( hdata ) 
        
    ## conversion to string 
    def __str__ (  self ) :
        return '%s(%s,xvar=%s,yvar=%s)' % (
            self.__class__.__name__ , self.name , self.xvar.name , self.yvar.name )
    __repr__ = __str__ 
    

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
    'derivative'     ,
    ##
    'roo_moment'     , 
    'roo_mean'       ,
    'roo_variance'   ,
    'roo_dispersion' ,
    'roo_rms'        ,
    'roo_skewness'   ,
    'roo_kurtosis'   ,    
    ) :

    if hasattr ( PDF2 , _a ) :
        def _suppress_ ( self , *args , **kwargs ) :
            raise AttributeError ( "'%s' object has no attribute '%s'" % ( type(self) , _a ) )
        setattr ( PDF2 , _a , _suppress_ ) 
        logger.verbose ( 'Remove attribute %s from PDF2' ) 


 
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
    def __init__ ( self , pdf , xvar , yvar ,
                   name           = None    ,
                   special        = False   ,                  
                   add_to_signals = True    ,
                   prefix         = ''      ,
                   suffix         = ''      ) :

        assert isinstance ( xvar , ROOT.RooAbsReal ) , "``xvar'' must be ROOT.RooAbsReal"
        assert isinstance ( yvar , ROOT.RooAbsReal ) , "``yvar'' must be ROOT.RooAbsReal"        
        assert isinstance ( pdf  , ROOT.RooAbsReal  ) , "``pdf'' must be ROOT.RooAbsReal"
        
        name = name if name else self.generate_name ( prefix = prefix + '%s_' % pdf.GetName() , suffix = suffix ) 
        PDF2  . __init__ ( self , name , xvar , yvar , special = special )

        if not self.special : 
            assert isinstance ( pdf  , ROOT.RooAbsPdf  ) , "``pdf''  must be ROOT.RooAbsPdf"

        ## PDF! 
        self.pdf = pdf

        if not self.xvar in self.params () : 
            self.warning ( "Function/PDF does not depend on xvar=%s" % self.xvar.name )
        if not self.yvar in self.params () : 
            self.warning ( "Function/PDF does not depend on yvar=%s" % self.yvar.name )

        ## add it to the list of signal components ?
        self.__add_to_signals = True if add_to_signals else False
        
        if self.add_to_signals :
            self.signals.add ( self.pdf )
       
        ## save the configuration
        self.config = {
            'pdf'            : self.pdf            ,
            'xvar'           : self.xvar           ,
            'yvar'           : self.yvar           ,
            'name'           : self.name           ,
            'special'        : self.special        ,  
            'add_to_signals' : self.add_to_signals , 
            'prefix'         : prefix              ,
            'suffix'         : suffix              ,            
            }

        self.checked_keys.add ( 'pdf'     )
        self.checked_keys.add ( 'xvar'    )
        self.checked_keys.add ( 'yvar'    )
        self.checked_keys.add ( 'special' )
        
    @property
    def add_to_signals ( self ) :
        """``add_to_signals'' : should PDF be added into list of signal components?"""
        return self.__add_to_signals 


# =============================================================================
## @class Combine2D
#  Non-extended sum of several PDFs
#  It is just a small wrapper for <code>ROOT.RooAddPdf</code>
#  @see RooAddPdf 
class Combine2D (PDF2) :
    """Non-extended sum of several PDFs:
    
    It is just a small wrapper for <code>ROOT.RooAddPdf</code>
    - see RooAddPdf 
    
    >>> sum  = Combine2D ( [ pdf1 , pdf2 , pdf3 ]  ) 

    """
    def __init__ ( self             ,
                   pdfs             , ## input list of PDFs  
                   xvar      = None , 
                   yvar      = None , 
                   name      = ''   ,
                   recursive = True ,
                   prefix    = 'f'  , ## prefix for fraction names 
                   suffix    = 'f'  , ## suffix for fraction names 
                   fractions = None ) :

        assert 2 <= len ( pdfs ) , 'Combine1D: at least two PDFs are needed!'

        pdf_list = []        
        for i , p in enumerate ( pdfs ) :

            if isinstance ( p , PDF2 ) :
                
                assert ( not xvar ) or xvar is p.xvar, "Invalid xvar/pdf%d.xvar: %s/%s" % ( i , xvar , p.xvar ) 
                assert ( not yvar ) or yvar is p.yvar, "Invalid yvar/pdf%d.yvar: %s/%s" % ( i , yvar , p.yvar ) 
                xvar = p.xvar
                yvar = p.yvar
                pdf_list.append ( p )
                
            elif isinstance ( pdf1 , ROOT.RooAbsPdf )                   \
                     and xvar and isinstance ( xvar , ROOT.RooAbsReal ) \
                     and yvar and isinstance ( xvar , ROOT.RooAbsReal ) :
                
                pdf_list.append ( Generic2D_pdf ( p , xvar , yvar ) )
                
            else :
                raise TypeError ( "Invalid type: pdf%d xvar %s/%s xvar=%s, yvar=%s" % ( i , p , type(p) , xvar , yvar ) )
        
        ## check the name 
        name = name if name else self.generate_name ( prefix = 'sum2' )
        
        ## ininialize the base class
        PDF2.__init__ ( self , name , xvar , yvar ) 

        for i , p in enumerate ( pdf_list )  :
            if p.pdf.canBeExtended() : self.warning ("``pdf%f'' can be extended!" % i ) 
                

        while prefix.endswith  ('_') : prefix = prefix[:-1]
        while suffix.startswith('_') : suffix = suffix[1:]
        
        self.__prefix    = prefix if prefix else 'f'
        self.__suffix    = suffix
        self.__recursive = True if recursive else False
        
        if 2 < len ( pdf_list ) : 
            if self.prefix and self.suffix : fr_name = '%s_%%d_%s' % ( self.prefix , self.suffix )
            else                           : fr_name = '%s_%%d'    %   self.prefix
        else :
            if self.prefix and self.suffix : fr_name = '%s_%s'     % ( self.prefix , self.suffix )
            else                           : fr_name =                 self.prefix

        ## make list of fractions 
        fraction_list = self.make_fractions  ( len ( pdf_list )           ,
                                               name      = fr_name        , 
                                               recursive = self.recursive ,
                                               fractions = fractions      )

        ## keep them
        self.__pdfs      = tuple ( pdf_list )
        self.__fractions = tuple ( fraction_list ) 
        
        for p in self.__pdfs      : self.alist1.add ( p.pdf )
        for f in self.__fractions : self.alist2.add ( f     )
        
        ## finally build PDF 
        self.pdf = ROOT.RooAddPdf ( self.roo_name ( 'combine2' ) ,
                                    ' + '.join ( '(%s)' % p.name for p in self.pdfs  ) ,
                                    self.alist1    ,
                                    self.alist2    ,
                                    self.recursive )
        
        self.config = {
            'pdfs'      : self.pdfs      ,
            'xvar'      : self.xvar      ,
            'yvar'      : self.yvar      ,
            'name'      : self.name      , 
            'prefix'    : self.prefix    ,
            'suffix'    : self.suffix    ,
            'fractions' : self.fractions ,
            'recursive' : self.recursive        
            }
  
        
    @property
    def prefix ( self ) :
        """``prefix'' : prefix for fraction names"""
        return self.__prefix

    @property
    def suffix ( self ) :
        """``suffix'' : suffix for fraction names"""
        return self.__suffix

    @property
    def recursive ( self ) :
        """``recursive'' : recursive fractions?"""
        return self.__recursive
    
    @property
    def pdfs ( self ) :
        """``pdfs'' : get list/tuple of involved PDFs (same as ``components'')"""
        return self.__pdfs
    @property
    def components ( self ) :
        """``components'' : get list/tuple of involved PDFs (same as ``pdfs'')"""
        return self.pdfs
        
    @property
    def fractions ( self ) :
        """``fractions'' : get involved fractions (same as ``F'')"""
        return self.component_getter ( self.__fractions )
    @fractions.setter
    def fractions ( self , values ) :
        self.component_setter ( self.__fractions , values )

    @property
    def F         ( self ) :
        """``F'' : get involved fractions (same as ``fractions'')"""
        return self.component_getter ( self.__fractions )
    @F.setter
    def F         ( self , values ) :
        self.fractions = values 




# =============================================================================
## @class Sum2D
#  Non-extended sum of two PDFs
#  @code
#  pdf1 = ...
#  pdf2 = ...
#  sum  = Sum2D ( pdf1 , pdf2 ) 
#  @endcode
#  It is just a small wrapper for <code>ROOT.RooAddPdf</code>
#  @see RooAddPdf 
class Sum2D (Combine2D) :
    """Non-extended sum of two PDFs:    
    It is just a small wrapper for `ROOT.RooAddPdf`    
    >>> pdf1 = ...
    >>> pdf2 = ...
    >>> sum  = Sum2D ( pdf1 , pdf2 ) 
    - see `ROOT.RooAddPdf` 
    """
    def __init__ ( self             ,
                   pdf1             ,
                   pdf2             ,  
                   xvar      = None , 
                   yvar      = None ,
                   name      = ''   , 
                   prefix    = 'f'  ,
                   suffix    = ''   ,
                   fraction  = None ,
                   others    = []   ,
                   recursive = True ) :                    
        
        ## check the name 
        name = name if name else self.generate_name ( prefix = 'sum2' )

        ## initialize the base class 
        Combine2D.__init__ ( self                  ,
                             name      = name      , 
                             pdfs      =  [ pdf1 , pdf2] + others ,
                             xvar      = xvar      ,
                             yvar      = yvar      ,
                             recursive = recursive ,         
                             prefix    = prefix    ,                             
                             suffix    = suffix    ,
                             fractions = fraction  )


        self.config = {
            'pdf1'      : self.pdf1      ,
            'pdf2'      : self.pdf2      ,
            'xvar'      : self.xvar      ,
            'yvar'      : self.yvar      ,
            'name'      : self.name      ,
            'prefix'    : self.prefix    ,
            'suffix'    : self.suffix    , 
            'fraction'  : self.fraction  ,
            'others'    : self.others    ,            
            'recursive' : self.recursive ,
            }
        
    @property
    def pdf1 ( self ) :
        """``pdf1'' : the first PDF"""
        return self.pdfs[0]
    
    @property
    def pdf2 ( self ) :
        """``pdf2'' : the second PDF"""
        return self.pdfs[1]

    @property 
    def others ( self ) :
        """``others'' : other PDFs (if any)"""
        return self.pdfs[2:]
    
    @property
    def fraction ( self ) :
        """``fraction'' : the fraction of the first PDF in the sum (same as ``fractions'')"""
        return self.fractions 
    @fraction.setter
    def fraction ( self , value ) :
        self.fractions = value 

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
    def __init__ ( self , xvar , yvar , name = '' ,  title = '' ) :

        name = name if name else self.generate_name ( prefix = 'flat2D_')                            
        PDF2.__init__ ( self  , name , xvar , yvar ) 
                        
        if not title : title = 'flat2(%s)' % name 
        self.pdf = Ostap.Models.Uniform ( name , title , self.xvar , self.yvar )
        assert 2 == self.pdf.dim() , 'Flat2D: wrong dimensionality!'
        
        ## save configuration
        self.config = {
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'name'     : self.name ,            
            'title'    : title     ,             
            }                   

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
        else : raise AttributeError ( "Invalid ``x-model'' argument: %s" % xmodel )

        if   isinstance ( ymodel , PDF            ) : self.__ymodel = ymodel
        elif isinstance ( ymodel , ROOT.RooAbsPdf ) and xvar :
            self.__ymodel = Generic1D_pdf  ( ymodel , yvar )
        else : raise AttributeError ( "Invalid ``y-model'' argument: %s" % ymodel )

        name  = name  if name  else self.generate_name ( 'Model2D_%s_%s'  % ( self.xmodel.name , self.ymodel.name ) )
        
        ## initialize the base 
        PDF2.__init__ (  self , name , self.__xmodel.xvar , self.__ymodel.xvar ) 

        ## check the title 
        if not title : title = '%s x %s' % ( self.__xmodel.name , self.__ymodel.name )
        
        def _triv_ ( m ) :
            _U = Ostap.Models.Uniform 
            if     isinstance ( m , Flat1D ) : return True 
            return isinstance ( m.pdf , _U ) and 1 == m.pdf.dim()

        ## trivial case: 
        if _triv_ ( self.xmodel ) and _triv_ ( self.ymodel ) :
            
            self.debug ('use Flat2D-model for the trivial product')
            self.__flat = Flat2D ( self.xvar , self.yvar , name = self.generate_name ( name ) , title = title )
            self.pdf    = self.__flat.pdf  
            
        else :
            
            ## build the final PDF 
            self.pdf = ROOT.RooProdPdf (
                self.roo_name ( 'model2_' ) , 
                title             ,
                self.__xmodel.pdf ,
                self.__ymodel.pdf )

        ## save configuration 
        self.config = {
            'name'   :  self.name   ,
            'xmodel' :  self.xmodel ,
            'ymodel' :  self.ymodel ,
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
## Generic 2D-shape from C++ callable
#  @see Ostap::Models:Shape2D
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2020-07-20
class Shape2D_pdf(PDF2) :
    """ Generic 2D-shape from C++ callable
    - see Ostap::Models:Shape2D
    """
    
    def __init__ ( self , name , shape , xvar , yvar ) :

        if isinstance ( shape , ROOT.TH2 ) and not isinstance ( shape , ROOT.TH3 ) and not xvar :
            xvar = shape.xminmax()
            
        if isinstance ( shape , ROOT.TH2 ) and not isinstance ( shape , ROOT.TH3 ) and not yvar :
            yvar = shape.yminmax()
            
        if isinstance ( shape , ROOT.TH2 ) and not isinstance ( shape , ROOT.TH3 ) :
            self.histo = shape
            shape      = Ostap.Math.Histo2D ( shape )

        ##  iniialize the base 
        PDF2.__init__ ( self , name , xvar , yvar ) 
            
        self.__shape = shape
        
        if isinstance ( self.shape , Ostap.Math.Histo2D ) :
            
            ## create the actual pdf
            self.pdf = Ostap.Models.Histo2D ( self.roo_name ( 'histo2_' ) , 
                                              "Histo-2D %s" % self.name   ,
                                              self.xvar                   ,
                                              self.yvar                   ,
                                              self.shape                  )            
        else : 
            
            ## create the actual pdf
            self.pdf = Ostap.Models.Shape2D.create  (
                self.roo_name  ( 'shape2_' ) , 
                "Shape-2D %s" % self.name    ,
                self.xvar                    ,
                self.yvar       `            ,
                self.shape                   )  
            
        ## save the configuration
        self.config = {
            'name'    : self.name    , 
            'shape'   : self.shape   , 
            'xvar'    : self.xvar    , 
            'yvar'    : self.yvar    , 
            }
        
    @property
    def shape  ( self ) :
        """``shape'': the actual C++ callable shape"""
        return self.__shape 
 
# ===================================================it==========================
## simple convertor of 2D-histogram into PDF
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H2D_pdf(H2D_dset,PDF2) :
    """Simple convertor of 2D-histogram into PDF 
    """
    def __init__ ( self            ,
                   name            ,
                   histo           ,
                   xvar    = None  , 
                   yvar    = None  ,
                   density = False ,
                   order   = 0     , ## interpolation order 
                   silent  = False ) :
        
        H2D_dset.__init__ ( self , histo , xvar , yvar , density , silent )
        PDF2    .__init__ ( self , name  , self.xaxis , self.yaxis ) 
        
        self.__vset  = ROOT.RooArgSet  ( self.xvar , self.yvar )

        assert isinstance ( order, integer_types ) and 0 <= order ,\
               'Invalid interpolation order: %s/%s' % ( order , type ( order ) )
        
        #
        ## finally create PDF :
        #
        with roo_silent ( silent ) : 
            self.pdf    = ROOT.RooHistPdf (
                self.roo_name  ( 'histo2_' ) , 
                'Histo-2D PDF: %s/%s' % ( self.histo.GetName() , self.histo.GetTitle() ) , 
                self.__vset , 
                self.dset   ,
                order       )
            
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
            'order'   : self.order   ,             
            }
        
    @property
    def order  ( self ) :
        """``order'': interpolation order"""
        return self.pdf.getInterpolationOrder () 
    @order.setter
    def order  ( self , value ) :
        assert isinstance ( value , integer_types ) and 0 <= value,\
               'Invalid interpolation order %s/%s' % ( value , type ( value ) )
        self.pdf.setInterpolationOrder ( value )


# =============================================================================
# Compound 2D-fit models 
# =============================================================================

# =============================================================================
## @class Fit2D
#  The actual model for 2D-fits. It consists of four main components :
#    - pure signal :        S(x)*S(y)
#    - signal x background: S(x)*B(y)
#    - signal x bakcgronud: B(x)*S(y)
#    - pure backrground:    B(x,y) or B(x)*B(y) 
#  Other 2D-components could be specified in addition
#  @param  signal_x  PDF for the S(x)-signal component
#  @param  signal_y  PDF for the S(y)-signal component
#  @param  suffix    suffix to be used for the PDF and variable names
#  @param  bkg_1x    x-background component for B(x)*S(y) term 
#  @param  bkg_1y    y-background component for S(x)*B(y) term
#  @param  bkg_2x    x-background component for B(x)*B(y) term, if <code>bkg_2D</code> is not specified 
#  @param  bkg_2y    y-background component for B(x)*B(y) term, if <code>bkg_2D</code> is not specified 
#  @param  bkg_2D    PDF for 2D-background component for B(x,y)    term
#  @param  sig_2D    PDF for 2D-signal component for S(x,y)        term
#  @param  ss        the yield of  S(x)*S(y) component
#  @param  sb        the yield of  S(x)*B(y) component
#  @param  bs        the yield of  B(x)*S(y) component
#  @param  bb        the yield of  B(x,y)    component
#  @param  componens list of other 2D-components
#  @param  xvar      the x-variable
#  @param  yvar      the y-variable
#  @param  name      the name of PDF 
#  @code
#  model   = Models.Fit2D (
#      signal_1 = Models.Gauss_pdf ( 'Gx' , m_x.getMin () , m_x.getMax () , mass = m_x ) ,
#      signal_2 = Models.Gauss_pdf ( 'Gy' , m_y.getMin () , m_y.getMax () , mass = m_y ) ,
#      bkg_1x   = 1 , 
#      bkg_1y   = 1 )
#
#  r,f = model.fitTo ( dataset ) ## fit dataset 
#
#  print r                       ## get results  
#
#  fx  = model.draw1 ()          ## visualize X-projection
#  fy  = model.draw2 ()          ## visualize Y-projection#
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class Fit2D (PDF2) :
    """The actual model for 2D-fits
    
    It consists of four main components :
    1. pure signal :        S(x,y)
    2. signal x background: S(x)*B(y)
    3. signal x backgronud: B(x)*S(y)
    4. pure backrground:    B(x,y) or B(x)*B(y) 
    Other 2D-components could be specified in addition
    
    Arguments:
    
    - signal_x        : PDF for the S(x)-signal component
    - signal_y        : PDF for the S(y)-signal component
    - suffix          : the suffix to be used for the PDF and variable names
    - bkg_1x          : x-background component for B(x)*S(y) term 
    - bkg_1y          : y-background component for S(x)*B(y) term
    - bkg_2x          : x-background component for B(x)*B(y) term, if bkg2D is not specified 
    - bkg_2y          : y-background component for B(x)*B(y) term, if bkg2D is not specified 
    - bkg_2D          : PDF for 2D-background component for B(x,y)    term
    - sig_2D          : PDF for 2D-signal component S(x,y) term
    - ss              : the yield of  S(x,y)    component
    - sb              : the yield of  S(x)*B(y) component
    - bs              : the yield of  B(x)*S(y) component
    - bb              : the yield of  B(x,y)    component
    - components      : the list of other 2D-components
    - xvar            : the x-variable
    - yvar            : the y-variable
    - name            : the name of PDF
    
    Example:
    
    >>>  model   = Models.Fit2D (
    ...      signal_x = Models.Gauss_pdf ( 'Gx' , mass = m_x ) ,
    ...      signal_y = Models.Gauss_pdf ( 'Gy' , mass = m_y ) ,
    ...      bkg1x    = 1 , 
    ...      bkg1y    = 1 )
    >>> r,f = model.fitTo ( dataset ) ## fit dataset 
    >>> print r                       ## get results  
    >>> fx  = model.draw1 ()          ## visualize X-projection
    >>> fy  = model.draw2 ()          ## visualize Y-projection

    """
    def __init__ ( self               ,
                   #
                   signal_x           , 
                   signal_y           ,
                   suffix = ''        ,
                   #
                   bkg_1x     = None  ,
                   bkg_1y     = None  ,
                   #
                   bkg_2x     = None  ,
                   bkg_2y     = None  ,
                   #
                   bkg_2D     = None  ,
                   sig_2D     = None  , ## 2D-signal component 
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
            'signal_x'   : signal_x   ,
            'signal_y'   : signal_y   ,
            'bkg_1x'     : bkg_1x     ,
            'bkg_1y'     : bkg_1y     ,
            'bkg_2x'     : bkg_2x     ,
            'bkg_2y'     : bkg_2y     ,
            'bkg_2D'     : bkg_2D     ,
            'sig_2D'     : sig_2D     ,
            'components' : components ,
            ##
            'ss'         : ss ,
            'bb'         : bb ,
            'sb'         : sb ,
            'bs'         : bs ,
            ##
            'suffix'     : suffix   ,
            'name'       : name     ,
            }
        
        
        self.__suffix      = suffix

        if   isinstance ( signal_x , PDF            )          : self.__signal_x = signal_x
        elif isinstance ( signal_x , ROOT.RooAbsPdf ) and xvar :
            self.__signal_x = Generic1D_pdf ( signal_x , xvar , prefix = 'SX' , suffix = suffix )
        else : raise AttributeError("Invalid ``signal_x'' argument: %s" % signal_x )
            
        if   isinstance ( signal_y , PDF            )          : self.__signal_y = signal_y
        elif isinstance ( signal_y , ROOT.RooAbsPdf ) and yvar :
            self.__signal_y = Generic1D_pdf ( signal_y , yvar , prefix = 'SY' , suffix = suffix )
        else : raise AttributeError("Invalid ``signal_y'' argument: %s" % signal_y )
            
        #
        ## initialize base class
        #
        if not name :
            name = self.generate_name ( "fit2:%s&%s" % ( self.__signal_x.name , self.__signal_y.name ) )
            if suffix : name += '_' + suffix 
            
        PDF2.__init__ ( self , name          ,
                        self.__signal_x.xvar ,
                        self.__signal_y.xvar ) 

        
        # =====================================================================
        ## Build components for the  final 2D-PDF
        # =====================================================================
        
        # =====================================================================
        ## First component: Signal(1) and Signal(2)
        # =====================================================================

        if   sig_2D and isinstance ( sig_2D , PDF2 ) :
            self.__ss_cmp = sig_2D
        elif sig_2D and isinstance ( sig_2D , ROOT.RooAbsPdf ) :
            self.__ss_cmp = Generic2D_pdf ( sig_2D , self.xvar , self.yvar , self.generate_name ( 'SS_' + self.name ) )
        elif not sig_2D : 
            self.__ss_cmp = Model2D ( self.generate_name ( "SS_" + self.name ) ,
                                      self.__signal_x   ,
                                      self.__signal_y   , 
                                      title = "Signal(x) x Signal(y)" )
        else :
            raise TypeError("Fit2D: can't create Signal(x,y)-component!")
        
        # =====================================================================
        ## Second component: Background(1) and Signal(2)
        # =====================================================================

        self.__bkg_1x = self.make_bkg ( bkg_1x  , self.generate_name ( 'Bkg1X_BS' + self.name ) , self.xvar )
        self.__bs_cmp = Model2D ( self.generate_name ( "BS_" + self.name ) ,
                                  self.__bkg_1x             ,
                                  self.__signal_y           ,
                                  title = "Backround1(x) x Signal(y)" )
        
        # =====================================================================
        ## Third component:  Signal(1) and Background(2)
        # =====================================================================
        
        self.__bkg_1y = self.make_bkg ( bkg_1y   , self.generate_name ( 'Bkg1Y_SB' + self.name ) , self.yvar )
        self.__sb_cmp = Model2D ( self.generate_name ( "SB_" + self.name ) ,
                                  self.__signal_x           ,
                                  self.__bkg_1y             ,
                                  title = "Signal(x) x Background1(y)" )
            
        # =====================================================================
        ## (intermezzo) assumptions about the background sub-components 
        # =====================================================================
        
        if   component_clone   ( bkg_2x ) :
            bkg_2x = self.__bkg_1x
            self.debug ( 'bkg_2x set to [CLONE]   %s' % bkg_2x ) 
        elif component_similar ( bkg_2x ) :
            bkg_2x =        bkg_1x
            self.debug ( 'bkg_2x set to [SIMILAR] %s' % bkg_2x )

        if   component_clone   ( bkg_2y ) :
            bkg_2y = self.__bkg_1y
            self.debug ( 'bkg_2y set to [CLONE]   %s' % bkg_2y )
        elif component_similar ( bkg_2x ) :
            bkg_2y =        bkg_1y
            self.debug ( 'bkg_2y set to [SIMILAR] %s' % bkg_2y ) 
            
        # =====================================================================
        ## Fourth component: Background(1) and Background(2) 
        # =====================================================================
    
        self.__bkg_2x = None 
        self.__bkg_2y = None 

        if   isinstance ( bkg_2D , PDF2             ) : self.__bb_cmp = bkg_2D        
        elif isinstance ( bkg_2D , ROOT.RooAbsPdf   ) : ## generic PDF 
            self.__bb_cmp  = Generic2D_pdf  ( bkg_2D , self.xvar , self.yvar , self.generate_name ( 'BB_' + self.name ) )
        elif isinstance ( bkg_2D , ( tuple , list ) ) : ## polynomials 
            from ostap.fitting.models_2d import make_B2D
            self.__bb_cmp = make_B2D ( self.generate_name ( "BB_" + self.name  ) , self.xvar , self.yvar , *bkg_2D )
        else     :                       
            self.__bkg_2x = self.make_bkg ( bkg_2x , self.generate_name ( 'Bkg2X_BB' + self.name ) , self.xvar )
            self.__bkg_2y = self.make_bkg ( bkg_2y , self.generate_name ( 'Bkg2Y_BB' + self.name ) , self.yvar )            
            self.__bb_cmp = Model2D ( self.generate_name ( "BB_" + self.name ) ,
                                      self.__bkg_2x             ,
                                      self.__bkg_2y             ,
                                      title = "Background2(x) x Background2(y)" )
            
        # =====================================================================
        ## coefficients/yields 
        # =====================================================================
    
        self.__ss = self.make_var ( ss   , "SS"               + suffix ,
                                    "Signal(x,y)"             + suffix , None , 1000  , 0 , 1.e+7 )
        self.__sb = self.make_var ( sb   ,  "SB"              + suffix ,
                                    "Signal(x)&Background(y)" + suffix , None ,  100  , 0 , 1.e+7 )
        self.__bs = self.make_var ( bs   , "BS"               + suffix ,
                                    "Background(x)&Signal(y)" + suffix , None ,  100  , 0 , 1.e+7 )        
        self.__bb = self.make_var ( bb   , "BB"               + suffix ,
                                    "Background(x,y)"         + suffix , None ,   10  , 0 , 1.e+7 )
        
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
        for i , cmp in enumerate ( components ) :
            
            if   isinstance  ( cmp , PDF2           ) : cc = cmp  
            elif isinstance  ( cmp , ROOT.RooAbsPdf ) : cc = Generic2D_pdf ( cmp ,  self.xvar , self.yvar ,  prefix = 'C%d_' % i , suffix = suffix ) 
            else :
                self.error ("unknown ``other''component %s/%s, skip it!" % ( cc , type(cc) ) )
                continue  
            self.__more_components.append ( cc     )
            self.components.add           ( cc.pdf ) 

        nc = len( self.__more_components )
        if 1 == nc :
            cf = self.make_var ( None , "C"+suffix , "Component" + suffix , None , 1 , 0 , 1.e+7 )
            self.alist1.add  ( self.components[0] )
            self.__nums_components.append ( cf ) 
        elif 2 <= nc : 
            fic = self.make_fracs ( nc , 'C_%%d%s' % suffix ,  'C(%%d)%s'  % suffix , fractions  = False )
            for c in self.components : self.alist1.add ( c)
            for f in fic             : self.__nums_components.append ( f )
            
        self.__nums_components  = tuple ( self.__nums_components  ) 
        for c in self.__nums_components  : self.alist2.add ( c )

        #
        ## build the final PDF 
        # 
        pdfname  = self.roo_name ( 'fit2d_' ) 
        pdftitle = "Fit2D %s" % self.name
        pdfargs  = pdfname , pdftitle , self.alist1 , self.alist2
        self.pdf = ROOT.RooAddPdf  ( *pdfargs )

        self.signals     .add ( self.__ss_cmp.pdf )
        self.backgrounds .add ( self.__bb_cmp.pdf )
        self.crossterms1 .add ( self.__sb_cmp.pdf ) ## cross-terms 
        self.crossterms2 .add ( self.__bs_cmp.pdf ) ## cross-terms 

        ## save configuration
        self.config = {
            'signal_x'   : self.signal_x        ,
            'signal_y'   : self.signal_y        ,            
            'suffix'     : self.suffix          ,
            'bkg_1x'     : self.bkg_1x          , 
            'bkg_1y'     : self.bkg_1y          , 
            'bkg_2x'     : self.bkg_2x          , 
            'bkg_2y'     : self.bkg_2y          , 
            'bkg_2D'     : self.bkg_2D          ,
            'sig_2D'     : self.sig_2D          ,
            'ss'         : self.SS              ,
            'sb'         : self.SB              ,
            'bs'         : self.BS              ,
            'bb'         : self.BB              ,
            'components' : self.more_components ,
            'xvar'       : self.xvar            , 
            'yvar'       : self.yvar            , 
            'name'       : self.name    
            }
        
        self.checked_keys.add  ( 'xvar' )
        self.checked_keys.add  ( 'yvar' )
        
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
        return self.component_getter ( self.__nums_components  )     
    @C.setter
    def C (  self , value ) :
        self.component_setter ( self.__nums_components , value )

    @property
    def yields    ( self ) :
        """The list/tuple of the yields of all numeric components"""
        return tuple ( [ i for i in  self.alist2 ] )
    
    @property 
    def total_yield ( self ) :
        """``total_yield''' : get the total yield"""
        if hasattr ( self , 'extended' ) and not self.extended : return None 
        if not self.fit_result                                 : return None
        if not valid_pointer ( self.fit_result )               : return None
        yields = self.yields
        if not yields                                          : return None
        ##  if 1 ==  len ( yields ) : return yield[0]. 
        return self.fit_result.sum ( *yields ) 
 

    # =========================================================================
    # components
    # =========================================================================

    @property 
    def signal_x ( self  ) :
        """Signal(x) component/PDF"""
        return self.__signal_x

    @property 
    def signal_y ( self  ) :
        """Signal(y) component/PDF"""
        return self.__signal_y

    @property
    def bkg_1x( self ) :
        """``bkg_1x'': The background(x) PDF for Backgroud(x)*Signal(y) component/PDF"""
        return self.__bkg_1x
    
    @property
    def bkg_1y( self ) :
        """``bkg_1y'': The background(y) PDF for Signal(x)*Background(y) component/PDF"""
        return self.__bkg_1y

    @property
    def bkg_2x( self ) :
        """``bkg_2x'': The background(x) PDF for Backgroud(x)*Background(y) component/PDF, when bkg2D is not specified"""
        return self.__bkg_2x

    @property
    def bkg_2y( self ) :
        """``bkg_2y'': The background(y) PDF for Backgroud(x)*Background(y) component/PDF, when bkg2D is not specified"""
        return self.__bkg_2y

    @property
    def bkg_2D( self ) :
        """``bkg_2D'': The PDF for Backgroud(x,y) component/PDF (same as cmp_BB)"""
        return self.__bb_cmp

    @property
    def sig_2D( self ) :
        """``sig_2D'': The PDF for Signal(x,y) component/PDF (same as cmp_SS)"""
        return self.__ss_cmp 

    @property
    def more_components ( self ) :
        """additional ``other'' components"""
        return tuple( self.__more_components  )

    @property
    def suffix ( self ) :
        """``suffix'', used to build the name"""
        return self.__suffix
    
    @property
    def cmp_SS ( self ) :
        """``cmp_SS'' : Signal(x&y))            component in the fit (PDF)"""
        return self.__ss_cmp
    @property
    def cmp_SB ( self ) :
        """``cmp_SB'' : Signal(x)xBackground(y) component in the fit (PDF)"""
        return self.__sb_cmp
    @property
    def cmp_BS ( self ) :
        """``cmp_BS'' : Background(x)xSignal(y) component in the fit (PDF)"""
        return self.__bs_cmp
    @property
    def cmp_BB ( self ) :
        """``BB_cmp'' : Background(x,y)         component in the fit (PDF)"""
        return self.__bb_cmp
    
# =============================================================================
## @class Fit2DSym
#  The actual model for symmetric 2D-fits. It consists of three main components :
#  1. pure signal :        S(x,y)
#  2. signal x background & background x signal:  S(x)*B(y) + B(x)*S(y)
#  3. pure backrground:    B(x,y) or B(x)*B(y) 
#  Other 2D-components could be specified in addition
#  @param  signal_x  PDF for the S(x)-signal component
#  @param  signal_y  PDF for the S(y)-signal component, cloned from S(x), if None 
#  @param  suffix    suffix to be used for the PDF and variable names
#  @param  bkg_1x    x-background component for B(x)*S(y) term, B(y) is cloned 
#  @param  bkg_2x    x-background component for B(x)*B(y) term, if bkg2D is not specified, B(y) is cloned  
#  @param  bkg_2D     PDF for (symmetric) 2D-background component for B(x,y)    term
#  @param  sig_2D     PDF for (symmetric) 2D-signal     component for S(x,y)    term
#  @param  ss         the yield of  S(x)*S(y) component
#  @param  sb         the yield of  S(x)*B(y)+B(x)*S(y)component
#  @param  bb         the yield of  B(x,y)    component
#  @param  components list of other 2D-components
#  @param  xvar       the x-variable
#  @param  yvar       the y-variable
#  @param  name       the name of PDF 

#  @code
# 
#  model   = Models.Fit2D (
#      signal_x = Models.Gauss_pdf ( 'Gx' , mass = m_x ) ,
#      signal_y = Models.Gauss_pdf ( 'Gy' , mass = m_y ) ,
#      bkg_1x   = 1 , 
#      bkg_2x   = 1 )
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
    The actual model for symmetric 2D-fits. It consists of three main components :
    - pure signal                               :        S(x)*S(y)
    - signal x background & background x signal :  S(x)*B(y) + B(x)*S(y)
    - pure background:                          :  B(x,y) or B(x)*B(y) 
    Other 2D-components could be specified in addition
    
    Arguments:
    
    - signal_x       : PDF for the S(x)-signal component
    - signal_y       : PDF for the S(y)-signal component; cloned from S(x), if None 
    - suffix         : suffix to be used for the PDF and variable names
    - bkg_1x         : x-background component for B(x)*S(y) term; B(y) is cloned 
    - bkg_2x         : x-background component for B(x)*B(y) term, if bkg2D is not specified; B(y) is   cloned  
    - bkg_2D         : PDF for (symmetric) 2D-background component for B(x,y)    term
    - sig_2D         : PDF for (symmetric) 2D-signal     component for S(x,y)    term
    - ss             : the yield of  S(x)*S(y) component
    - sb             : the yield of  S(x)*B(y)+B(x)*S(y)component
    - bb             : the yield of  B(x,y)    component
    - componens      : the list of other 2D-components
    - xvar           : the x-variable
    - yvar           : the y-variable
    - name           : the name of PDF 
    
    Example:
    
    >>>  model   = Models.Fit2DSym (
    ...      signal_x = Models.Gauss_pdf ( 'Gx' , xvar = m_x ) ,
    ...      signal_y = Models.Gauss_pdf ( 'Gy' , xvar = m_y ) ,
    ...      bkg1x    = 1 , 
    ...      bkg2x    = 1 )
    >>> r,f = model.fitTo ( dataset ) ## fit dataset 
    >>> print r                       ## get results  
    >>> fx  = model.draw1 ()          ## visualize X-projection
    >>> fy  = model.draw2 ()          ## visualize X-projection

    """
    def __init__ ( self               ,
                   #
                   signal_x           , 
                   signal_y   = None  ,
                   suffix     = ''    ,
                   #
                   bkg_1x     = None  ,
                   bkg_2x     = None  ,
                   bkg_2D     = None  ,
                   sig_2D     = None  ,
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
            'signal_x'   : signal_x   ,
            'signal_y'   : signal_y   ,
            'bkg_1x'     : bkg_1x     , 
            'bkg_2x'     : bkg_2x     , 
            'bkg_2D'     : bkg_2D     ,
            'sig_2D'     : sig_2D     ,
            'components' : components ,
            ##
            'ss'         : ss         ,
            'sb'         : sb         , 
            'bb'         : bb         ,
            ##
            'suffix'     : suffix     ,
            'name'       : name       ,
            ##
            'xvar'       : xvar       ,
            'yvar'       : yvar       ,            
            }
        
        self.__suffix      = suffix

        # ============================================================================
        if   isinstance ( signal_x , PDF            )          : self.__signal_x = signal_x
        elif isinstance ( signal_x , ROOT.RooAbsPdf ) and xvar :
            self.__signal_x = Generic1D_pdf ( signal_x , xvar , prefix = 'SX' , suffix = suffix )
        else : raise AttributeError ( "Invalid ``signal_x'' argument: %s" %  signal_x )
            
        if   isinstance ( signal_y , PDF            )          : self.__signal_y = signal_y
        elif isinstance ( signal_y , ROOT.RooAbsPdf ) and yvar :
            self.__signal_y = Generic1D_pdf ( signal_y , yvar , prefix = 'SX' , suffix = suffix )
        elif yvar and not signal_y :
            self.__signal_y = self.__signal_x.clone ( xvar = yvar , name_prefix = 'SY_' )
            self.debug('signal y-component is cloned from the signal_x component')
        else : raise AttributeError ( "Invalid ``signal_y'' argument: %s" % signal_y )
            
        # =====================================================================
        ## initialize base class
        # =====================================================================
        if not name :
            name = "%s&%s" % ( self.__signal_x.name , self.__signal_y.name )
            if suffix : name += '_' + suffix 

        ##  initialize the base class 
        PDF2.__init__ ( self , name          ,
                        self.__signal_x.xvar ,
                        self.__signal_y.xvar ) 
        
        # =====================================================================
        ## First component: Signal(1) and Signal(2)
        # =====================================================================

        if   sig_2D and isinstance ( sig_2D , PDF2 ) :
            self.__ss_cmp = sig_2D
        elif sig_2D and isinstance ( sig_2D , ROOT.RooAbsPdf ) :
            self.__ss_cmp = Generic2D_pdf ( sig_2D , self.xvar , self.yvar , self.generate_name ( 'SS_' + self.name  ) ) 
        elif not sig_2D : 
            self.__ss_cmp = Model2D ( self.generate_name ( 'SS_' + self.name ) ,
                                      self.__signal_x    ,
                                      self.__signal_y    , 
                                      title = "Signal(x) x Signal(y)" )
        else :
            raise TypeError("Fit2D: can't create Signal(x,y)-component!")
        
        self.__bkg_1x = self.make_bkg (        bkg_1x , self.generate_name ( 'Bkg1X_BS' + self.name ) , self.xvar )
        self.__bkg_1y = self.make_bkg ( self.__bkg_1x , self.generate_name ( 'Bkg1Y_SB' + self.name ) , self.yvar )

        # =====================================================================
        ## Second sub-component: Background (1) and Signal     (2)
        ## Third  sub-component: Signal     (1) and Background (2)
        # =====================================================================

        self.__sb_cmp_raw = Model2D ( self.generate_name ( "S1B2_" + self.name ) ,
                                      self.__signal_x           ,
                                      self.__bkg_1y             ,
                                      title = "Signal(x) x Background(y)" )
        
        self.__bs_cmp_raw = Model2D ( self.generate_name ( "B1S2_" + self.name ) , 
                                      self.__bkg_1x             ,
                                      self.__signal_y           ,    
                                      title = "Background(x) x Signal(y)" )
        
        self.__sb_cmp     = Generic2D_pdf (
            self.make_sum ( self.generate_name ( "SB_" + self.name ) ,
                            "Signal(x) x Background(y) + Background(x) x Signal(y)"   ,
                            self.__sb_cmp_raw.pdf ,
                            self.__bs_cmp_raw.pdf ) , self.xvar , self.yvar , 'SB_' + self.name )
        
        ## alias, just for convinience 
        self.__bs_cmp    = self.__sb_cmp
        
        # =====================================================================
        ## (intermezzo) Assumptions about the background sub-components 
        # =====================================================================
        
        if   component_clone   ( bkg_2x ) :
            bkg_2x = self.__bkg_1x
            self.debug ( 'bkg_2x set to [CLONE]   %s' % bkg_2x ) 
        elif component_similar ( bkg_2x ) :
            bkg_2x =        bkg_1x
            self.debug ( 'bkg_2x set to [SIMILAR] %s' % bkg_2x ) 
            
        # =====================================================================
        ## fourth component: Background(1) and Background(2) 
        # =====================================================================
    
        self.__bkg_2x = None
        self.__bkg_2y = None

        bb_name = self.generate_name ( 'BB_' + self.name )
        if   isinstance ( bkg_2D , PDF2           ) : self.__bb_cmp = bkg_2D  
        elif isinstance ( bkg_2D , ROOT.RooAbsPdf ) :
            self.__bb_cmp  = Generic2D_pdf  ( bkg_2D , self.xvar , self.yvar , bb_name )
        elif isinstance ( bkg_2D , int ) :

            from ostap.fitting.models_2d import make_B2Dsym
            self.__bb_cmp = make_B2Dsym ( bb_name , self.xvar , self.yvar , bkg_2D )

        else     :                        
            self.__bkg_2x = self.make_bkg (        bkg_2x , self.generate_name ( 'Bkg2X_BB' + self.name ) , self.xvar )
            self.__bkg_2y = self.make_bkg ( self.__bkg_2x , self.generate_name ( 'Bkg2Y_BB' + self.name ) , self.yvar )
            self.__bb_cmp = Model2D ( bb_name       ,
                                      self.__bkg_2x ,
                                      self.__bkg_2y ,
                                      title = "Background2(x) x Backrgound2(y)" )
            
        # =====================================================================
        ## coefficients
        # =====================================================================
        
        self.__ss = self.make_var ( ss , "SS"             + suffix ,
                                    "Signal(x)&Signal(y)" + suffix , None , 1000  , 0 ,  1.e+7 )
        
        self.__bb = self.make_var ( bb , "BB"             + suffix ,
                                    "Background(x,y)"     + suffix , None ,   10  , 0 ,  1.e+7 )
        
        self.__sb = self.make_var ( sb , "SB"             + suffix ,
                                    "Signal(x)&Background(y)+Background(x)&Signal(y)" + suffix , None ,  100  , 0 ,  1.e+7 )
        
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
        for i , cmp in enumerate ( components ) :
            
            if   isinstance  ( cmp , PDF2           ) : cc = cmp  
            elif isinstance  ( cmp , ROOT.RooAbsPdf ) : cc = Generic2D_pdf ( cmp ,  self.xvar , self.yvar , prefix = 'C%d_' % i , suffix = suffix )
            else :
                self.error ("unknown ``other''component %s/%s, skip it!" % ( cc , type(cc) ) )
                continue  
            self.__more_components.append ( cc     )
            self.components.add           ( cc.pdf ) 

        nc = len( self.__more_components )
        if 1 == nc :
            cf = self.make_var ( None , "C"+suffix , "Component" + suffix , None , 1 , 0 , 1.e+7 )
            self.alist1.add  ( self.components[0] )
            self.__nums_components.append ( cf ) 
        elif 2 <= nc : 
            fic = self.make_fracs ( nc , 'C_%%d%s' % suffix ,  'C(%%d)%s'  % suffix , fractions  = False )
            for c in self.components : self.alist1.add ( c)
            for f in fic             : self.__nums_components.append ( f )
            
        self.__nums_components  = tuple ( self.__nums_components  ) 
        for c in self.__nums_components  : self.alist2.add ( c )
            
        #
        ## build the final PDF 
        #
        pdfname  = self.roo_name ( 'fit2ds_' ) 
        pdftitle = "Fit2Dsym %s" % self.name
        pdfargs  = pdfname , pdftitle , self.alist1 , self.alist2
        self.pdf = ROOT.RooAddPdf  ( *pdfargs )


        self.signals     .add ( self.__ss_cmp.pdf )
        self.backgrounds .add ( self.__bb_cmp.pdf )
        self.crossterms1 .add ( self.__sb_cmp.pdf ) ## cross-terms 
        self.crossterms2 .add ( self.__bs_cmp.pdf ) ## cross-terms 

        ## save configuration
        self.config = {
            'signal_x'   : self.signal_x        ,
            'signal_y'   : self.signal_y        ,            
            'suffix'     : self.suffix          ,
            'bkg_1x'     : self.bkg_1x          , 
            'bkg_2x'     : self.bkg_2x          , 
            'bkg_2D'     : self.bkg_2D          ,
            'sig_2D'     : self.sig_2D          ,
            'ss'         : self.SS              ,
            'sb'         : self.SB              ,
            'bb'         : self.BB              ,
            'components' : self.more_components ,
            'xvar'       : self.xvar            , 
            'yvar'       : self.yvar            , 
            'name'       : self.name            ,
            }

        self.checked_keys.add ( 'xvar' )
        self.checked_keys.add ( 'yvar' )
        
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
        return self.component_getter ( self.__nums_components )
    @C.setter
    def C (  self , value ) :
        self.component_setter ( self.__nums_componenents , value )

    @property
    def yields    ( self ) :
        """The list/tuple of the yields of all numeric components"""
        return tuple ( [ i for i in  self.alist2 ] )

    @property 
    def total_yield ( self ) :
        """``total_yield''' : get the total yield"""
        if not self.fit_result                   : return None
        if not valid_pointer ( self.fit_result ) : return None
        return self.fit_result.sum ( *self.yields ) 
 
    # =========================================================================
    # components
    # =========================================================================
    
    @property 
    def signal_x ( self  ) :
        """``signal_x'': Signal(x) component/PDF"""
        return self.__signal_x

    @property 
    def signal_y ( self  ) :
        """``signal_y'': Signal(y) component/PDF"""
        return self.__signal_y

    @property
    def bkg_1x( self ) :
        """``bkg_1x'': The background PDF for Backgroud(x)*Signal(y) component/PDF"""
        return self.__bkg_1x
    
    @property
    def bkg_1y( self ) :
        """``bkg_1y'': The background PDF for Signal(x)*Background(y) component/PDF"""
        return self.__bkg_1y 

    @property
    def bkg_2x( self ) :
        """``bkg_2x'': The background(x) PDF for Backgroud(x)*Background(y) component/PDF"""
        return self.__bkg_2x
    
    @property
    def bkg_2y( self ) :
        """``bkg_2y'': The background(y) PDF for Backgroud(x)*Background(y) component/PDF"""
        return self.__bkg_2y 

    @property
    def bkg_2D( self ) :
        """``bkg_2D'': The PDF for Backgroud(x&y) component/PDF (same as cmp_BB)"""
        return self.__bb_cmp

    @property
    def sig_2D( self ) :
        """``sig_2D'': The PDF for Signal(x&y) component/PDF (same as cmp_SS)"""
        return self.__ss_cmp
 
    @property
    def more_components ( self ) :
        """additional ``other'' components"""
        return tuple( self.__more_components  )

    @property
    def suffix ( self ) :
        """``suffix'', used to build the name"""
        return self.__suffix
    
    @property
    def cmp_SS ( self ) :
        """``cmp_SS'' : Sig(1&2) component in the fit (PDF)"""
        return self.__ss_cmp
    @property
    def cmp_SB ( self ) :
        """``cmp_SB'' : Sig(1)xBkg(2)+Bkg(1)*Sig(2) component in the fit (PDF)"""
        return self.__sb_cmp
    @property
    def cmp_BS ( self ) :
        """``cmp_BS'' : Sig(1)xBkg(2)+Bkg(1)*Sig(2) component in the fit (PDF) (same as  cmp_SB)"""
        return self.__bs_cmp
    @property
    def cmp_BB ( self ) :
        """``cmp_BB'' : Bkg(1&2)      component in the fit (PDF)"""
        return self.__bb_cmp
 
    # =========================================================================
    ## Raw, non-symmetrized fit components/PDF (for debugging)
    # =========================================================================

    def cmp_raw_SB ( self ) :
        """``cmp_SB'' : Sig(1)xBkg(2) raw, non-symmetrized component in the fit (PDF)"""
        return self.__sb_cmp_raw
    @property
    def cmp_BS ( self ) :
        """``cmp_BS'' :  Bkg(1)*Sig(2) raw, non-symmetrized component in the fit (PDF)"""
        return self.__bs_cmp_raw

 
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
