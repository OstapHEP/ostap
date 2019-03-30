#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/fit3d.py
#  Set of useful basic utilities to build various fit models 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
# =============================================================================
"""Set of useful basic utilities to build various 2D-fit models"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'PDF3'          , ## useful base class for 3D-models
    'Fit3D'         , ## the model for                3D-fit
    'Fit3DSym'      , ## the model for      symmetric 3D-fit
    'Fit3DMix'      , ## the model for half-symmetric 3D-fit
    ##
    'Generic3D_pdf' , ## wrapper over imported RooFit (2D)-pdf  
    'Flat3D'        , ## the most trivial 3D-pdf - constant
    'Model3D'       , ## trivial class to build 3D model from 1D-components 
    'H3D_pdf'       , ## convertor of 1D-histo to RooDataPdf
    )
# =============================================================================
import ROOT, random
from   ostap.core.core      import dsID , hID ,  VE , Ostap 
from   ostap.core.types     import integer_types
from   ostap.logger.utils   import roo_silent , rooSilent
from   ostap.fitting.utils  import H3D_dset , component_similar , component_clone
from   ostap.fitting.basic  import PDF  , Flat1D 
from   ostap.fitting.fit2d  import PDF2 , Model2D 
from   ostap.fitting.roofit import SETVAR
from   builtins             import range
# =============================================================================
from   ostap.logger.logger  import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.fit3d' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
# @class PDF3
# The helper base class for implementation of 3D-pdfs 
# @author Vanya BELYAEV Ivan.Belyaev@itep.ru
# @date 2017-11-11
class PDF3 (PDF2) :
    """ Useful helper base class for implementation of PDFs for 3D-fit
    """
    def __init__ ( self , name , xvar = None , yvar = None , zvar  = None , special = False ) : 
        
        PDF2.__init__ ( self , name , xvar , yvar , special = special ) 
        
        ## create the variable 
        if isinstance ( zvar , tuple ) and 2 == len(zvar) :  
            self.__zvar = self.make_var ( zvar               , ## var 
                                    'z'                , ## name 
                                    'z-variable(mass)' , ## title/comment
                                    None               , ## fix ?
                                    *zvar              ) ## min/max 
        elif isinstance ( zvar , ROOT.RooAbsReal ) :
            self.__zvar = self.make_var ( zvar               , ## var 
                                    'z'                , ## name 
                                    'z-variable/mass'  , ## title/comment
                                    fix = None         ) ## fix ? 
        else :
            self.warning('``z-variable''is not specified properly %s/%s' % ( zvar , type ( zvar ) ) )
            self.__zvar = self.make_var( zvar , 'z' , 'z-variable' )


        self.vars.add ( self.__zvar )
        
        ## save the configuration
        self.config = {
            'name' : self.name ,
            'xvar' : self.xvar ,
            'yvar' : self.yvar ,            
            'zvar' : self.zvar ,            
            }
        
    def zminmax ( self ) :
        """Min/max values for z-variable"""
        return self.__zvar.minmax()
    
    @property 
    def zvar ( self ) :
        """``z''-variable for the fit (same as ``z'')"""
        return self.__zvar

    @property 
    def z    ( self ) :
        """``z''-variable for the fit (same as ``zvar'')"""
        return self.__zvar

    # =========================================================================
    ## make the actual fit 
    #  @code
    #  r,f = model.fitTo ( dataset )
    #  r,f = model.fitTo ( dataset , weighted = True )    
    #  r,f = model.fitTo ( dataset , ncpu     = 10   )    
    #  r,f = model.fitTo ( dataset )    
    #  @endcode 
    def fitTo ( self           , 
                dataset        ,
                silent = False ,
                refit  = False ,
                timer  = False ,
                args   = ()    , **kwargs ) :
        """
        Perform the actual fit (and draw it)
        >>> r,f = model.fitTo ( dataset )
        >>> r,f = model.fitTo ( dataset , weighted = True )    
        >>> r,f = model.fitTo ( dataset , ncpu     = 10   )    
        >>> r,f = model.fitTo ( dataset )    
        """
        if   isinstance ( dataset , H3D_dset ) : dataset = dataset.dset        
        elif isinstance ( dataset , ROOT.TH3 ) :
            density = kwargs.pop ( 'density' , True  ) 
            chi2    = kwargs.pop ( 'chi2'    , False ) 
            return self.fitHisto ( dataset   ,
                                   draw    = draw    ,
                                   silent  = silent  ,
                                   density = dentity ,
                                   chi2    = chi2    , args = args , **kwargs )
        
        result,f = PDF2.fitTo ( self    ,
                                dataset = dataset ,
                                draw    = False   , ## False here!
                                nbins   = 50      , ## fake  here!
                                ybins   = 20      , ## fake  here!
                                silent  = silent  ,
                                refit   = refit   ,
                                timer   = timer   ,
                                args    = args    , **kwargs )
        
        return result
    
    # =========================================================================
    ## draw the projection over 1st variable
    #
    #  @code
    #  r,f = model.fitTo ( dataset ) ## fit dataset
    #  fx  = model.draw1 ( dataset , nbins = 100 ) ## draw results
    #
    #  fx  = model.draw1 ( dataset , nbins = 100 , in_range2 = (2,3) ) ## draw results
    #
    #  model.yvar.setRange ( 'QUQU2' , 2 , 3 ) 
    #  fx  = model.draw1 ( dataset , nbins = 100 , in_range2 = 'QUQU2') ## draw results
    #
    #  @endcode 
    def draw1 ( self            ,
                dataset   = None ,
                nbins     = 100  ,
                silent    = True ,
                in_range2 = None ,
                in_range3 = None , **kwargs ) :
        """ Draw the projection over 3rd variable
        
        >>> r,f = model.fitTo ( dataset ) ## fit dataset
        >>> fx  = model.draw1 ( dataset , nbins = 100 ) ## draw results
        
        >>> fx  = model.draw1 ( dataset , nbins = 100 , in_range2 = (2,3) ) ## draw results

        >>> model.yvar.setRange ( 'QUQU2' , 2 , 3 ) 
        >>> fx  = model.draw1 ( dataset , nbins = 100 , in_range2 = 'QUQU2') ## draw results
        
        """
        if in_range2 and isinstance ( in_range2 , tuple ) and 2 == len ( in_range2 ) :
            with rooSilent ( 3 ) : self.yvar.setRange ( 'aux_rng2' , in_range2[0] , in_range2[1] )
            in_range2 = 'aux_rng2'

        if in_range3 and isinstance ( in_range3 , tuple ) and 2 == len ( in_range3 ) :
            with rooSilent ( 3 ) : self.zvar.setRange ( 'aux_rng3' , in_range3[0] , in_range3[1] )
            in_range3 = 'aux_rng3'

        in_range = []
        if in_range2 : in_range.append( in_range2 )
        if in_range3 : in_range.append( in_range3 )
        in_ranage = tuple( in_range ) 
        return self.draw ( drawvar  = self.xvar , 
                           dataset  = dataset   ,
                           nbins    = nbins     ,
                           ybins    = 20        , ## fake 
                           silent   = silent    ,
                           in_range = in_range  , **kwargs )


    # =========================================================================
    ## draw the projection over 2nd variable
    #
    #  @code
    #  r,f = model.fitTo ( dataset ) ## fit dataset
    #  fy  = model.draw1 ( dataset , nbins = 100 ) ## draw results
    #
    #  fy  = model.draw1 ( dataset , nbins = 100 , in_range1 = (2,3) ) ## draw results
    #
    #  model.xvar.setRange ( 'QUQU1' , 2 , 3 ) 
    #  fy  = model.draw1 ( dataset , nbins = 100 , in_range1 = 'QUQU1') ## draw results
    #
    #  @endcode 
    def draw2 ( self            ,
                dataset   = None ,
                nbins     = 100  ,
                silent    = True ,
                in_range1 = None ,
                in_range3 = None , **kwargs ) :
        """ Draw the projection over 2nd variable
        
        >>> r,f = model.fitTo ( dataset ) ## fit dataset
        >>> fy  = model.draw2 ( dataset , nbins = 100 ) ## draw results
        
        >>> fx  = model.draw2 ( dataset , nbins = 100 , in_range1 = (2,3) ) ## draw results

        >>> model.xvar.setRange ( 'QUQU1' , 2 , 3 ) 
        >>> fx  = model.draw2 ( dataset , nbins = 100 , in_range1 = 'QUQU1') ## draw results
        
        """
        if in_range1 and isinstance ( in_range1 , tuple ) and 2 == len ( in_range1 ) :
            with rooSilent ( 3 ) : self.xvar.setRange ( 'aux_rng1' , in_range1[0] , in_range1[1] )
            in_range1 = 'aux_rng1'

        if in_range3 and isinstance ( in_range3 , tuple ) and 2 == len ( in_range3 ) :
            with rooSilent ( 3 ) : self.zvar.setRange ( 'aux_rng3' , in_range3[0] , in_range3[1] )
            in_range3 = 'aux_rng3'

        in_range = []
        if in_range1 : in_range.append( in_range1 )
        if in_range3 : in_range.append( in_range3 )
        in_ranage = tuple( in_range ) 
        return self.draw ( drawvar  = self.yvar , 
                           dataset  = dataset   ,
                           nbins    = nbins     ,
                           ybins    = 20        , ## fake 
                           silent   = silent    ,
                           in_range = in_range  , **kwargs )


    # =========================================================================
    ## draw the projection over 3rd variable
    #
    #  @code
    #  r,f = model.fitTo ( dataset ) ## fit dataset
    #  fz  = model.draw3 ( dataset , nbins = 100 ) ## draw results
    #
    #  fz  = model.draw3 ( dataset , nbins = 100 , in_range2 = (2,3) ) ## draw results
    #
    #  model.yvar.setRange ( 'QUQU2' , 2 , 3 ) 
    #  f  = model.draw3 ( dataset , nbins = 100 , in_range2 = 'QUQU2') ## draw results
    #  @endcode 
    def draw3 ( self             ,
                dataset   = None ,
                nbins     = 100  ,
                silent    = True ,
                in_range1 = None ,
                in_range2 = None , **kwargs ) :
        """ Draw the projection over 3rd variable
        
        >>> r,f = model.fitTo ( dataset ) ## fit dataset
        >>> fx  = model.draw3 ( dataset , nbins = 100 ) ## draw results
        
        >>> fx  = model.draw3 ( dataset , nbins = 100 , in_range2 = (2,3) ) ## draw results

        >>> model.yvar.setRange ( 'QUQU2' , 2 , 3 ) 
        >>> fx  = model.draw3 ( dataset , nbins = 100 , in_range2 = 'QUQU2') ## draw results
        
        """
        if in_range1 and isinstance ( in_range1 , tuple ) and 2 == len ( in_range1 ) :
            with rooSilent ( 3 ) : self.xvar.setRange ( 'aux_rng1' , in_range1[0] , in_range1[1] )
            in_range1 = 'aux_rng1'

        if in_range2 and isinstance ( in_range2 , tuple ) and 2 == len ( in_range2 ) :
            with rooSilent ( 3 ) : self.yvar.setRange ( 'aux_rng2' , in_range2[0] , in_range2[1] )
            in_range2 = 'aux_rng2'

        in_range = []
        if in_range1 : in_range.append( in_range1 )
        if in_range2 : in_range.append( in_range2 )
        in_ranage = tuple( in_range ) 
        return self.draw ( drawvar  = self.zvar , 
                           dataste  = dataset   ,
                           nbins    = nbins     ,
                           ybins    = 20        , ## fake 
                           silent   = silent    ,
                           in_range = in_range  , **kwargs )

    # =========================================================================
    ## fit the 3D-histogram (and draw it)
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
                   chi2    = False , 
                   args    = ()    , **kwargs ) :
        """Fit the histogram (and draw it)
        
        >>> histo = ...
        >>> r,f = model.fitHisto ( histo , draw = True )
        
        """
        
        xminmax = histo.xminmax()
        yminmax = histo.yminmax()
        zminmax = histo.zminmax()
        
        with     RangeVar ( self.xvar , *xminmax ) , \
                 RangeVar ( self.yvar , *yminmax ) , \
                 RangeVar ( self.xvar , *zminmax ): 
            
            ## convert it! 
            self.histo_data = H3D_dset ( histo , self.xvar , self.yvar  , self.zvar ,
                                         density , silent )
            data = self.histo_data
            if chi2 : return self.chi2fitTo ( data              ,
                                              draw    = draw    ,
                                              silent  = False   ,
                                              density = density ,
                                              args    = args    , **kwargs )
            else    : return self.fitTo     ( data              ,
                                              silent  = silent  ,
                                              args    =  args   , **kwargs ) 

    # =========================================================================
    ## generate toy-sample according to PDF
    #  @code
    #  model  = ....
    #  data   = model.generate ( 10000 ) ## generate dataset with 10000 events
    #  varset = ....
    #  data   = model.generate ( 100000 , varset )
    #  data   = model.generate ( 100000 , varset , extended =  =   True )     
    #  @endcode
    def generate ( self ,  nEvents , varset = None , extended = False ,  *args ) :
        """Generate toy-sample according to PDF
        >>> model  = ....
        >>> data   = model.generate ( 10000 ) ## generate dataset with 10000 events
        
        >>> varset = ....
        >>> data   = model.generate ( 100000 , varset )
        >>> data   = model.generate ( 100000 , varset , extended = True )
        """
        args = args + ( ROOT.RooFit.Name ( dsID() ) , ROOT.RooFit.NumEvents ( nEvents ) )
        if  extended :
            args = args + ( ROOT.RooFit.Extended () , )
        if   not varset :
            varset = ROOT.RooArgSet( self.xvar , self.yvar , self.zvar )
        elif isinstance ( varset , ROOT.RooAbsReal ) :
            varset = ROOT.RooArgSet( varser )

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

        if not self.zvar in varset :
            vs = ROOT.RooArgSet()
            vs . add ( self.zvar )
            for  v in varset : vs.add ( v )
            varset = vs
            
        return self.pdf.generate ( varset , *args )


    # ====================================================================================
    ## simple 'function-like' interface 
    def __call__ ( self , x , y , z , error = False , normalized = True ) :
        
        if     isinstance ( self.xvar , ROOT.RooRealVar ) and \
               isinstance ( self.yvar , ROOT.RooRealVar ) and \
               isinstance ( self.zvar , ROOT.RooRealVar ) :
           
            if x in self.xvar and y in self.yvar and z in self.zvar : 
                with SETVAR( self.xvar ) , SETVAR( self.yvar ) ,  SETVAR( self.zvar ) :
                    self.xvar.setVal ( x )
                    self.yvar.setVal ( y )
                    self.zvar.setVal ( z )
                    
                    v = self.pdf.getVal ( self.vars ) if normalized else self.pdf.getValV ()
 
                    if error and self.fit_result :
                        e = self.pdf.getPropagatedError ( self.fit_result )
                        if 0<= e : return  VE ( v ,  e * e )
                    return v 
            else : return 0.0
            
        raise AttributeError ( 'something wrong goes here' )


    # ========================================================================
    ## check minmax of the PDF using the random shoots
    #  @code
    #  pdf     = ....
    #  mn , mx = pdf.minmax()            
    #  @endcode 
    def minmax ( self , nshoots = 200000 ) :
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
        code = self.pdf.getMaxVal( ROOT.RooArgSet ( self.xvar , self.yvar , self.zvar ) )
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
        if not self.zminmax() : return ()
        
        mn  , mx = -1 , -10
        xmn , xmx = self.xminmax()
        ymn , ymx = self.yminmax()
        zmn , zmx = self.zminmax()
        for i in range ( nshoots ) : 
            xx = random.uniform ( xmn , xmx )
            yy = random.uniform ( ymn , ymx )
            zz = random.uniform ( zmn , zmx )
            with SETVAR ( self.xvar ) :
                with SETVAR ( self.yvar ) :
                    with SETVAR ( self.zvar ) :
                        self.xvar.setVal ( xx )
                        self.yvar.setVal ( yy )
                        self.zvar.setVal ( zz )
                        vv = self.pdf.getVal()
                        if mn < 0 or vv < mn : mn = vv
                        if mx < 0 or vv > mx : mx = vv
                        
        return mn , mx 

    # =========================================================================
    ## get integral over (xmin,xmax,ymin,ymax,zmin,zmax) region
    #  @code
    #  pdf = ...
    #  print pdf.integral( 0,1,0,2,0,5)
    #  @endcode
    def integral ( self, xmin , xmax , ymin , ymax , zmin , zmax , nevents = True ) :
        """Get integral over (xmin,xmax,ymin,ymax,zmin,zmax) region
        >>> pdf = ...
        >>> print pdf.integral( 0,1,0,2,0,5)
        """
        if self.xminmax() :            
            xmn , xmx = self.xminmax()
            xmin = max ( xmin , xmn )
            xmax = min ( xmax , xmx )
            
        if self.yminmax() : 
            ymn , ymx = self.yminmax()
            ymin = max ( ymin , ymn )
            ymax = min ( ymax , ymx )
            
        if self.zminmax() :
            zmn , zmx = self.zminmax()
            zmin = max ( zmin , zmn )
            zmax = min ( zmax , zmx )

        
        value , todo = 0 , True 
        
         ## 1) make a try to use analytical integral (could be fast)
        if self.tricks :
            try:
                if hasattr ( self.pdf , 'setPars'  ) : self.pdf.setPars() 
                fun          = self.pdf.function()
                value , todo = fun.integral ( xmin , xmax ,
                                              ymin , ymax ,
                                              zmin , zmax ) , False 
            except:
                pass


        ## for numerical integration 
        from ostap.math.integral import integral3 as _integral3
        
        extended = self.pdf.canBeExtended() or isinstance ( self.pdf , ROOT.RooAddPdf )
        if todo  and extended  :
            value = _integral3 ( self , xmin , xmax , ymin , ymax , zmin , zmax )
        elif todo : 
                        
            ## use unormalized PDF here to speed up the integration 
            ifun   = lambda x , y , z : self ( x , y , z , error = False , normalized = False )
            value  = _integral3 ( ifun , xmin , xmax , ymin , ymax , zmin , zmax )
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
    #  x , y , z = pdf3.minimum() 
    #  @endcode 
    def minimum ( self ,
                  xmin = None , xmax = None ,
                  ymin = None , ymax = None ,
                  zmin = None , zmax = None , x0 = () ) :
        """Get a minimum of PDF for certain interval
        >>> pdf3 = ...
        >>> x, y , z = pdf3.minimum()
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

        if zmin is None : zmin = self.zminmax()[0]
        if zmax is None : zmax = self.zminmax()[1]
        if self.zminmax() :
            zmin =  max ( zmin , self.zminmax()[0] )
            zmax =  min ( zmax , self.zminmax()[1] )
            
        if not x0 : x0 = 0.5 * ( xmin + xmax ) , 0.5 * ( ymin + ymax ) , 0.5 * ( zmin + zmax )

        if not xmin <= x0[0] <= xmax :
            logger.error("Wrong xmin/x0[0]/xmax: %s/%s/%s"   % ( xmin , x0[0] , xmax ) )

        if not ymin <= x0[1] <= ymax : 
            logger.error("Wrong ymin/x0[1]/ymax: %s/%s/%s"   % ( ymin , x0[1] , ymax ) )
        
        if not zmin <= x0[2] <= zmax : 
            logger.error("Wrong zmin/x0[2]/zmax: %s/%s/%s"   % ( zmin , x0[2] , zmax ) )

        from ostap.math.minimize import sp_minimum_3D
        return sp_minimum_3D (  self ,
                                xmin , xmax ,
                                ymin , ymax ,
                                zmin , zmax , x0 )
    

    # ==========================================================================
    ## get a maximum of PDF for certain interval
    #  @code
    #  pdf2 = ...
    #  x , y , z = pdf3.maximum() 
    #  @endcode 
    def minimum ( self ,
                  xmin = None , xmax = None ,
                  ymin = None , ymax = None ,
                  zmin = None , zmax = None , x0 = () ) :
        """Get a maximum of PDF for certain interval
        >>> pdf3 = ...
        >>> x, y , z = pdf3.maximum()
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

        if zmin is None : zmin = self.zminmax()[0]
        if zmax is None : zmax = self.zminmax()[1]
        if self.zminmax() :
            zmin =  max ( zmin , self.zminmax()[0] )
            zmax =  min ( zmax , self.zminmax()[1] )
            
        if not x0 : x0 = 0.5 * ( xmin + xmax ) , 0.5 * ( ymin + ymax ) , 0.5 * ( zmin + zmax )

        if not xmin <= x0[0] <= xmax :
            logger.error("Wrong xmin/x0[0]/xmax: %s/%s/%s"   % ( xmin , x0[0] , xmax ) )

        if not ymin <= x0[1] <= ymax : 
            logger.error("Wrong ymin/x0[1]/ymax: %s/%s/%s"   % ( ymin , x0[1] , ymax ) )
        
        if not zmin <= x0[2] <= zmax : 
            logger.error("Wrong zmin/x0[2]/zmax: %s/%s/%s"   % ( zmin , x0[2] , zmax ) )

        from ostap.math.minimize import sp_maximum_3D
        return sp_maximum_3D ( self ,
                               xmin , xmax ,
                               ymin , ymax ,
                               zmin , zmax , x0 )
    

    # ==========================================================================
    ## convert PDF into TF2 object, e.g. to profit from TF3::Draw options
    #  @code
    #  pdf = ...
    #  tf3 = pdf.tf()
    #  tf3.Draw( options )
    #  @endcode
    def tf ( self                      ,
             xmin = None , xmax = None ,
             ymin = None , ymax = None ,
             zmin = None , zmax = None ) :
        """Convert PDF  to TF3 object, e.g. to profit from TF3::Draw options
        >>> pdf = ...
        >>> tf3 = pdf.tf()
        >>> tf3.Draw('colz')
        """
        def _aux_fun_ ( x , pars = [] ) :
            return self ( x[0] , x[1] , x[2] , error = False )

        if xmin == None and self.xminmax() : xmin = self.xminmax()[0]
        if xmax == None and self.xminmax() : xmax = self.xminmax()[1]
        if ymin == None and self.yminmax() : ymin = self.yminmax()[0]
        if ymax == None and self.yminmax() : ymax = self.yminmax()[1]
        if zmin == None and self.zminmax() : zmin = self.zminmax()[0]
        if zmax == None and self.zminmax() : zmax = self.zminmax()[1]
        
        if xmin == None : xmin = 0.0
        if xmax == None : xmin = 1.0
        if ymin == None : ymin = 0.0
        if ymax == None : ymin = 1.0
        if zmin == None : zmin = 0.0
        if zmax == None : zmin = 1.0
        
        from ostap.core.core import fID
        return ROOT.TF3 ( fID() , _aux_fun_ , xmin , xmax , ymin , ymax , zmin , zmax ) 


    # ==========================================================================
    ## create the histogram accoring to specifications 
    def make_histo ( self ,
                     xbins    = 10    , xmin = None , xmax = None ,
                     ybins    = 10    , ymin = None , ymax = None ,
                     zbins    = 10    , zmin = None , zmax = None ,
                     hpars    = ()    , 
                     histo    = None  ) :
        """Create the histogram accoring to specifications"""
        
        import ostap.histos.histos

        # histogram is provided 
        if histo :
            
            assert isinstance ( histo  , ROOT.TH3 ), \
                   "Illegal type of ``histo''-argument %s" % type( histo )
            
            histo = histo.clone()
            histo.Reset()

        # arguments for the histogram constructor 
        elif hpars :
            
            from ostap.core.core import hID
            histo = ROOT.TH3F ( hID () , 'PDF%s' % self.name , *hpars  )
            if not histo.GetSumw2() : histo.Sumw2()

        # explicit contruction from (#bins,min,max)-triplet  
        else :
            
            assert isinstance ( xbins , integer_types ) and 0 < xbins, \
                   "Wrong ``xbins''-argument %s" % xbins 
            assert isinstance ( ybins , integer_types ) and 0 < ybins, \
                   "Wrong ``ybins''-argument %s" % ybins 
            assert isinstance ( zbins , integer_types ) and 0 < zbins, \
                   "Wrong ``zbins''-argument %s" % zbins 
            if xmin == None and self.xminmax() : xmin = self.xminmax()[0]
            if xmax == None and self.xminmax() : xmax = self.xminmax()[1]
            if ymin == None and self.yminmax() : ymin = self.yminmax()[0]
            if ymax == None and self.yminmax() : ymax = self.yminmax()[1]
            if zmin == None and self.zminmax() : zmin = self.zminmax()[0]
            if zmax == None and self.zminmax() : zmax = self.zminmax()[1]
            
            from ostap.core.core import hID
            histo = ROOT.TH3F ( hID() , 'PDF%s' % self.name ,
                                xbins , xmin , xmax ,
                                ybins , ymin , ymax ,
                                zbins , zmin , zmax )
            if not histo.GetSumw2() : histo.Sumw2()

        return histo 


    # ==========================================================================
    ## Convert PDF to the 3D-histogram in correct  way 
    #  @code
    #  pdf = ...
    #  h1  = pdf.histo ( 10 , 0. , 10. , 10 , 0. , 4. , 10 , 0. , 3 ) ## specify histogram parameters
    #  histo_template = ...
    #  h2  = pdf.histo ( histo = histo_template ) ## use historgam template
    #  h3  = pdf.histo ( ... , integral = True  ) ## use PDF integral within the bin  
    #  @endcode
    def histo ( self             ,
                xbins    = 10    , xmin = None , xmax = None ,
                ybins    = 10    , ymin = None , ymax = None ,
                zbins    = 10    , zmin = None , zmax = None ,
                hpars    = ()    , 
                histo    = None  ,
                intergal = True  ,
                errors   = False ) :
        """Convert PDF to the 3D-histogram in correct way
        >>> pdf = ...
        >>> h1  = pdf.histo ( 10 , 0. , 10. , 10 , 0. , 4. , 10 , 0. , 3 ) ## specify histogram parameters
        >>> histo_template = ...
        >>> h2  = pdf.histo ( histo = histo_template ) ## use historgam template
        >>> h3  = pdf.histo ( ... , integral = True  ) ## use PDF integral within the bin  
        """


        histo = self.make_histo ( xbins = xbins , xmin = xmin , xmax = xmax ,
                                  ybins = ybins , ymin = ymin , ymax = ymax ,
                                  zbins = zbins , zmin = zmin , zmax = zmax ,
                                  hpars = hpars ,
                                  histo = histo )
        
        
        # loop over the historgam bins 
        for ix , iy , iz , x , y , z , w in histo.items() :

            xv , xe = x.value() , x.error()
            yv , ye = y.value() , y.error()
            zv , ze = z.value() , z.error()
            
            # value at the bin center 
            c = self ( xv , yv , zv , error = errors ) 

            if not integral : 
                histo[ix,iy,iz] = c
                continue

            # integral over the bin 
            v  = self.integral( xv - xe , xv + xe ,
                                yv - ye , yv + ye ,
                                zv - ze , zv + ze )
            
            if errors :
                if    0 == c.cov2 () : pass
                elif  0 != c.value() and 0 != v : 
                    v = c * ( v / c.value() )
                    
            histo[ix,iy,iz] = v 

        return histo
    
    # ==========================================================================
    ## Convert PDF to the 3D-histogram, taking PDF-values at bin-centres
    #  @code
    #  pdf = ...
    #  h1  = pdf.roo_histo ( 10 , 0. , 10. , 10 , 0. , 4. , 10 , 0. , 3 ) 
    #  histo_template = ...
    #  h2  = pdf.roo_histo ( histo = histo_template ) ## use historgam template
    #  @endcode
    def roo_histo ( self            ,
                   xbins    = 10    , xmin = None , xmax = None ,
                   ybins    = 10    , ymin = None , ymax = None ,
                   zbins    = 10    , zmin = None , zmax = None ,
                   hpars    = ()    , 
                   histo    = None  ,
                   events   = True) :
        """Convert PDF to the 3D-histogram, taking PDF-values at bin-centres
        >>> pdf = ...
        >>> h1  = pdf.roo_histo ( 10 , 0. , 10. , 10 , 0. , 4. , 10 , 0. , 3 )
        >>> histo_template = ...
        >>> h2  = pdf.roo_histo ( histo = histo_template ) ## use historgam template
        """
        
        histo = self.make_histo ( xbins = xbins , xmin = xmin , xmax = xmax ,
                                  ybins = ybins , ymin = ymin , ymax = ymax ,
                                  zbins = zbins , zmin = zmin , zmax = zmax ,
                                  hpars = hpars ,
                                  histo = histo )

        hh = self.pdf.createHistogram (
            hID()     ,
            self.xvar ,                    self.binning ( histo.GetXaxis() , 'histo3x' )   ,
            ROOT.RooFit.YVar ( self.yvar , self.binning ( histo.GetYaxis() , 'histo3y' ) ) , 
            ROOT.RooFit.ZVar ( self.zvar , self.binning ( histo.GetZaxis() , 'histo3z' ) ) , 
            ROOT.RooFit.Scaling  ( False ) , 
            ROOT.RooFit.Extended ( False ) )
        
        for i in hh : hh.SetBinError ( i , 0 ) 
        
        if events and self.pdf.mustBeExtended() :            
            for ix , iy , iz , x , y , z , v  in hh.items() :
                volume               = 8 * x.error()  * y.error() * z.error() 
                hh [ iz , iy , iz ] *= volume
                
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
        dataset.project ( hdata , ( self.xvar.name , self.yvar.name , self.xvar.name )  )
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
        dataset.project ( hdata , ( self.zvar.name , self.yvar.name , self.xvar.name ) ) 
        return self.pull_histo ( hdata ) 

    ## conversion to string 
    def __str__ (  self ) :
        return '%s(%s,xvar=%s,yvar=%s,zvar=%s)' % (
            self.__class__.__name__ , self.name ,
            self.xvar.name , self.yvar.name , self.zvar.name )
    __repr__ = __str__ 

# =============================================================================
## @class Generic3D_pdf
#  "Wrapper" over generic RooFit (3D)-pdf
#  @code
#     
#  raw_pdf = 
#  pdf     = Generic3D_pdf ( raw_pdf )  
# 
#  @endcode 
#  If more functionality is required , more actions are possible:
#  @code
#  ## for sPlot 
#  pdf.alist2 = ROOT.RooArgList ( n1 , n2 , n3 ) ## for sPlotting 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-03-29
class Generic3D_pdf(PDF3) :
    """ Wrapper for generic (3D) RooFit pdf:    
    >>> raw_pdf =
    >>> x,y,z   = 
    >>> pdf     = Generic3D_pdf ( raw_pdf ,  xvar = x ,   yvar = y , zvar =  z) 
    """
    ## constructor 
    def __init__ ( self , pdf , xvar , yvar , zvar ,
                   name           = None  ,
                   special        = False ,
                   add_to_signals = True  ) :
        
        assert isinstance ( xvar , ROOT.RooAbsReal ) , "``xvar'' must be ROOT.RooAbsReal"
        assert isinstance ( yvar , ROOT.RooAbsReal ) , "``yvar'' must be ROOT.RooAbsReal"
        assert isinstance ( zvar , ROOT.RooAbsReal ) , "``zvar'' must be ROOT.RooAbsReal"
        assert isinstance ( pdf  , ROOT.RooAbsReal ) , "``pdf''  must be ROOT.RooAbsReal"

        name = name if name else pdf.GetName ()
        
        PDF3  . __init__ ( self , name , xvar , yvar , zvar , special = special )

        if not self.special : 
            assert isinstance ( pdf  , ROOT.RooAbsPdf  ) , "``pdf'' must be ROOT.RooAbsPdf"

        ## PDF! 
        self.pdf = pdf

        ## add it to the list of signal components ?
        self.__add_to_signals = True if add_to_signals else False
        
        if self.add_to_signals :
            self.signals.add ( self.pdf )

        ## save the configuration
        self.config = {
            'pdf'            : self.pdf            ,
            'xvar'           : self.xvar           ,
            'yvar'           : self.yvar           ,
            'zvar'           : self.zvar           ,
            'name'           : self.name           ,
            'special'        : self.special        ,  
            'add_to_signals' : self.add_to_signals , 
            }

    
    @property
    def add_to_signals ( self ) :
        """``add_to_signals'' : shodul PDF be added into list of signal components?"""
        return self.__add_to_signals 
    
    ## redefine the clone method, allowing only the name to be changed
    #  @attention redefinition of parameters and variables is disabled,
    #             since it can't be done in a safe way                  
    def clone ( self , name = '' , xvar = None , yvar = None , zvar = None  ) :
        """Redefine the clone method, allowing only the name to be changed
         - redefinition of parameters and variables is disabled,
         since it can't be done in a safe way          
        """
        if xvar and not xvar is self.xvar :
            raise AttributeError("Generic2D_pdf can not be cloned with different `xvar''")
        if yvar and not yvar is self.yvar :
            raise AttributeError("Generic2D_pdf can not be cloned with different `yvar''")
        if zvar and not zvar is self.zvar :
            raise AttributeError("Generic2D_pdf can not be cloned with different `zvar''")
        return PDF.clone ( self , name = name ) if name else PDF.clone( self )

# =============================================================================
## @class Flat3D
#  The most trivial 3D-model - constant
#  @code 
#  pdf = Flat3D( 'flat' , xvar = ...  , yvar = ... , zvar = ... )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
class Flat3D(PDF3) :
    """The most trival 3D-model - constant
    >>> pdf = Flat3D( 'flat' , xvar = ...  , yvar = ... , zvar = ... )
    """
    def __init__ ( self , xvar , yvar , zvar , name = 'Flat3D'  , title = '' ) :
        
        PDF3.__init__ ( self  , name , xvar , yvar , zvar ) 
        
        if not title : title = 'flat3(%s)' % name 
        self.pdf = Ostap.Models.Uniform ( name , title , self.xvar , self.yvar , self.zvar )
        assert 3 == self.pdf.dim() , 'Flat3D: wrong dimensionality!'
        
        ## save configuration
        self.config = {
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'zvar'     : self.zvar ,
            'name'     : self.name ,            
            'title'    : title     ,             
            }
        
# =============================================================================
## @class Model3D
#  Trivial class to construct 3D model as a product of split 1D-models
#  actually it is a tiny  wrapper over <code>ROOT.RooProdPdf</code>
#  @code
#  pdfx = ...
#  pdfy = ...
#  pdfz = ...
#  pdf3D = Model3D( 'D3' , xmodel = pdfx , ymodel =  pdfy , zmodel = pdfz )
#  @endcode 
class Model3D(PDF3) :
    """Trivial class to construct 3D model as a product of split 1D-models
    - actually it is a tiny  wrapper over ROOT.RooProdPdf
    >>> pdfx = ...
    >>> pdfy = ...
    >>> pdfz = ...
    >>> pdf3D = Model3D( 'D3' , xmodel = pdfx , ymodel =  pdfy , zmodel = pdfz )
    """
    def __init__ ( self         ,
                   name         ,
                   xmodel       ,
                   ymodel       , 
                   zmodel       ,
                   xvar  = None ,
                   yvar  = None ,
                   zvar  = None ,
                   title = ''   ) :

        if   isinstance ( xmodel , PDF            ) : self.__xmodel = xmodel
        elif isinstance ( xmodel , ROOT.RooAbsPdf ) and xvar :
            self.__xmodel = Generic1D_pdf  ( xmodel , xvar )
        else : raise AttributeError ( "Invalid ``x-model'' attribute" )

        if   isinstance ( ymodel , PDF            ) : self.__ymodel = ymodel
        elif isinstance ( ymodel , ROOT.RooAbsPdf ) and yvar :
            self.__ymodel = Generic1D_pdf  ( ymodel , yvar )
        else : raise AttributeError ( "Invalid ``y-model'' attribute" )

        if   isinstance ( zmodel , PDF            ) : self.__zmodel = zmodel
        elif isinstance ( zmodel , ROOT.RooAbsPdf ) and zvar :
            self.__zmodel = Generic1D_pdf  ( zmodel , zvar )
        else : raise AttributeError ( "Invalid ``z-model'' attribute" )
        
        ## initialize the base 
        PDF3.__init__ (  self , name ,
                         self.__xmodel.xvar ,
                         self.__ymodel.xvar ,
                         self.__zmodel.xvar ) 

        self.__plst = ROOT.RooArgList (
            self.__xmodel.pdf ,
            self.__ymodel.pdf ,
            self.__zmodel.pdf ,
            )

        if not title : title = '%s x %s x %s' % ( self.__xmodel.name ,
                                                  self.__ymodel.name ,
                                                  self.__zmodel.name )
        
       
        def _triv_ ( m ) :
            
            _U = Ostap.Models.Uniform 
            if     isinstance ( m , Flat1D ) : return True 
            return isinstance ( m.pdf , _U ) and 1 == m.pdf.dim()
        
        ## trivial case: 
        if _triv_ ( self.xmodel ) and _triv_ ( self.ymodel ) and _triv_ ( self.zmodel ) :
            
            self.debug ('use Flat3D-model for the trivial product')
            self.__flat = Flat3D ( self.xvar , self.yvar , self.zvar , name = name , title = title )
            self.pdf    = self.__flat.pdf
            
        else :
            
            ## build pdf         
            self.__plst = ROOT.RooArgList (
                self.__xmodel.pdf ,
                self.__ymodel.pdf ,
                self.__zmodel.pdf ,
                )
            self.pdf = ROOT.RooProdPdf ( name , title , self.__plst )

        ## save configuration 
        self.config = {
            'name'   : self.name   ,
            'xmodel' : self.xmodel ,
            'ymodel' : self.ymodel ,
            'zmodel' : self.zmodel ,
            'xvar'   : self.xvar   ,
            'yvar'   : self.yvar   ,            
            'zvar'   : self.xvar   ,
            'title'  : self.pdf.GetTitle() 
            }

    @property
    def xmodel ( self ) :
        """``x-model'' x-component of M(x)*M(y)*M(z) PDF"""
        return self.__xmodel

    @property
    def ymodel ( self ) :
        """``y-model'' y-component of M(x)*M(y)*M(z) PDF"""
        return self.__ymodel
    
    @property
    def zmodel ( self ) :
        """``z-model'' z-component of M(x)*M(y)*M(z) PDF"""
        return self.__zmodel
    
# =============================================================================
## simple convertor of 3D-histogram into PDF
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H3D_pdf(H3D_dset,PDF3) :
    """Simple convertor of 3D-histogram into PDF 
    """
    def __init__ ( self            ,
                   name            ,
                   histo           ,
                   xvar    = None  , 
                   yvar    = None  ,
                   zvar    = None  ,
                   density = True  ,
                   silent  = False ) :
        
        H3D_dset.__init__ ( self , histo3 ,      xvar  ,      yvar  ,      zvar  ,  density , silent )
        PDF3    .__init__ ( self , name   , self.xaxis , self.yaxis , self.zaxis ) 
        
        self.__vset  = ROOT.RooArgSet  ( self.xvar , self.yvar , self.zvar )
        
        #
        ## finally create PDF :
        #
        with roo_silent ( silent ) : 
            self.pdf    = ROOT.RooHistPdf (
                'hpdf_%s'            % name ,
                'Histo3PDF(%s/%s/%s)' % ( name , histo3.GetName() , histo2.GetTitle() ) , 
                self.__vset  , 
                self.dset    )

        ## and declare it be be a "signal"
        self.signals.add ( self.pdf ) 
            
        ## save the configuration
        self.config = {
            'name'    : self.name    , 
            'histo'   : self.histo   , 
            'xvar'    : self.xvar    , 
            'yvar'    : self.yvar    , 
            'zvar'    : self.zvar    , 
            'density' : self.density , 
            'silent'  : self.silent  ,             
            }
             
# =============================================================================
# Compound models for 3D-fit
# =============================================================================

# =============================================================================
## @class Fit3D
#  The actual model for 3D-fits
#
#  @code
# 
#  model   = Models.Fit3D (
#      signal_x = Models.Gauss_pdf ( 'Gx' , mass = m_x ) ,
#      signal_y = Models.Gauss_pdf ( 'Gy' , mass = m_y ) ,
#      signal_z = Models.Gauss_pdf ( 'Gz' , mass = m_z ) )
#
#  r = model.fitTo ( dataset ) ## fit dataset 
#
#  print r                       ## get results  
#
#  fx  = model.draw1 ()          ## visualize X-projection
#  fy  = model.draw2 ()          ## visualize Y-projection
#  fz  = model.draw3 ()          ## visualize Z-projection
#
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2017-07-25
class Fit3D (PDF3) :
    """The actual model for 3D-fits
    
    >>>  model   = Models.Fit3D (
    ...      signal_x = Models.Gauss_pdf ( 'Gx' , mass = m_x ) ,
    ...      signal_y = Models.Gauss_pdf ( 'Gy' , mass = m_y ) ,
    ...      signal_x = Models.Gauss_pdf ( 'Gz' , mass = m_z ) ,
    ...      bkg_1x   = 1 , 
    ...      bkg_1y   = 0 ,
    ...      bkg_1z   = 0 )
    >>> r,f = model.fitTo ( dataset ) ## fit dataset 
    >>> print r                       ## get results  
    >>> fx  = model.draw1 ()          ## visualize X-projection
    >>> fy  = model.draw2 ()          ## visualize Y-projection
    >>> fz  = model.draw3 ()          ## visualize Z-projection

    Parameters
    ----------
    signal_x :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in X-direction
    signal_y :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in Y-direction
    signal_z :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in Z-direction
    suffix   : string
        An optional suffix to be  added to the names of created PDFs and variables
    bkg_1x    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for Bx1(x)* Sy(y)* Sz(z) term
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkg_1y    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for Sx(z)*By1(y)* Sz(z) term
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkg_1x    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for   Sx(x)* Sy(y)*Bz1(z) term 
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkg_2xy    : RooFit/PDF, Ostap/PDF or None
        2D-background for Bxy(x,y)*Sz(z) term  
        Use directly RooFit/PDF or Ostap/PDF otherwise create from bkgX2 and bkgY2
    bkg_2xz    : RooFit/PDF, Ostap/PDF or None
        2D-background for Bxz(x,z)*Sy(y) term 
        Use directly RooFit/PDF or Ostap/PDF otherwise create from bkgX2 and bkgZ2
    bkg_2yz    : RooFit/PDF, Ostap/PDF or None
        2D-background for Bxz(x,z)*Sy(y) term        
        Use directly RooFit/PDF or Ostap/PDF otherwise create from bkgX2 and bkgZ2
    bkg_2x    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for used to create Bxy(x,y) and Bxz(x,z) if they are not specified.
        If None - bkgX1 is used 
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkg_2y    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for used to create Bxy(x,y) and Byz(y,z) if they are not specified.
        If None - bkgY1 is used 
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkg_2z    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for used to create Bxz(x,z) and Byz(y,z) if they are not specified.
        If None - bkgZ1 is used
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkg_3x    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for used to create Bxyz(x,y,z) if it is not specified.
        If None - bkgX2 is used 
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkg_3y   : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for used to create Bxyz(x,y,z) if it is not specified.
        If None - bkgY2 is used 
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkg_3z   : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for used to create Bxyz(x,y,z) if it is not specified.
        If None - bkgZ2 is used         
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkg_3D   : RooFit/PDF or Ostap/PDF 3to descrive 3D-backround component Bxyz(x,y,z)
    

    sss      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of sig(1) * sig(2) * sig(3) component
         Use directly RooRealVar, otherwise create it using self.make_var function
    ssb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of sig(1) * sig(2) * bkg (3) component
         Use directly RooRealVar, otherwise create it using self.make_var function
    sbs      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of sig(1) * bkg(2) * sig (3) component
         Use directly RooRealVar, otherwise create it using self.make_var function
    bss      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of bkg(1) * sig(2) * sig (3) component
         Use directly RooRealVar, otherwise create it using self.make_var function
    sbb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of sig(1) * bkg(2) * bkg (3) component
         Use directly RooRealVar, otherwise create it using self.make_var function
    bsb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of bkg(1) * sig(2) * bkg (3) component
         Use directly RooRealVar, otherwise create it using self.make_var function
    bbs      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of bkg(1) * bkg(2) * sig (3) component
         Use directly RooRealVar, otherwise create it using self.make_var function
    bbb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of bkg(1) * bkg(2) * bkg (3) component
         Use directly RooRealVar, otherwise create it using self.make_var function

    """
    def __init__ ( self               ,
                   #
                   signal_x           , 
                   signal_y           ,
                   signal_z           ,
                   suffix = ''        ,
                   #
                   bkg_1x     = None  , ## 1D-background for Bx(x)*Sy(y)*Sz(z) term
                   bkg_1y     = None  , ## 1D-background for Sx(z)*By(y)*Sz(z) term
                   bkg_1z     = None  , ## 1D-background for Sx(x)*Sy(y)*Bz(z) term 
                   #
                   bkg_2xy    = None  , ## 2D-background for Bxy(x,y)*Sz(z) term 
                   bkg_2xz    = None  , ## 2D-background for Bxz(x,z)*Sy(y) term 
                   bkg_2yz    = None  , ## 2D-background for Byz(y,z)*Sx(z) term  
                   ## *if* no XY,XZ,BC backgrounds are specified, combine them from 
                   bkg_2x      = None  , ## Bxy(x,y) = Bx2(x)*By2(y)
                   bkg_2y      = None  , ## Bxz(x,z) = Bx2(x)*Bz2(z)
                   bkg_2z      = None  , ## Bkg(y,z) = By2(y)*Bz2(z)
                   ##
                   bkg_3D     = None  , ## 3D-backround component  B(x,y,z)
                   ## *if* no 3D-background components is specified, combine it from
                   bkg_3x     = None  , ## Bkg(x,y,z) = Bx3(x)*By3(y)*Bz3(z)
                   bkg_3y     = None  , ## Bkg(x,y,z) = Bx3(x)*By3(y)*Bz3(z)
                   bkg_3z     = None  , ## Bkg(x,y,z) = Bx3(x)*By3(y)*Bz3(z)
                   #
                   ## Yields of the main components :
                   sss        = None  , ## sig(x) * sig(y) * sig(z) 
                   ssb        = None  , ## sig(x) * sig(y) * bkg(z)
                   sbs        = None  , ## sig(x) * bkg(y) * sig(z)
                   bss        = None  , ## bkg(x) * sig(y) * sig(z)
                   sbb        = None  , ## sig(x) * bkg(y,z)
                   bsb        = None  , ## sig(y) * bkg(x,z)
                   bbs        = None  , ## sig(z) * bkg(x,y)
                   bbb        = None  , ## background-3D 
                   ## additional components 
                   components = []    ,
                   xvar       = None  ,
                   yvar       = None  ,
                   zvar       = None  ,                   
                   name       = ''    ) : 

        ## keep all arguments 
        self.__args = {
            ##
            'signal_x'   : signal_x ,
            'signal_y'   : signal_y ,
            'signal_z'   : signal_z ,
            ##
            'bkg_1x'     : bkg_1x   ,
            'bkg_1y'     : bkg_1y   ,
            'bkg_1z'     : bkg_1z   ,
            ## 
            'bkg_2xy'    : bkg_2xy  ,
            'bkg_2xz'    : bkg_2xz  ,
            'bkg_2yz'    : bkg_2yz  ,
            ##
            'bkg_2x'     : bkg_2x   ,
            'bkg_2y'     : bkg_2y   ,
            'bkg_2z'     : bkg_2z   ,
            ## 
            'bkg_3D'     : bkg_3D   ,
            ##
            'bkg_3x'     : bkg_3x   ,
            'bkg_3y'     : bkg_3y   ,
            'bkg_3z'     : bkg_3z   ,
            ## 
            'sss'        : sss      ,            
            'ssb'        : ssb      ,
            'sbs'        : sbs      ,
            'bss'        : bss      ,
            'sbb'        : sbb      ,
            'bsb'        : bsb      ,
            'bbs'        : bbs      ,
            'bbb'        : bbb      ,
            ##
            'components' : components ,
            'xvar'       : xvar     , 
            'yvar'       : yvar     , 
            'zvar'       : zvar     , 
            ##
            'name'     : name     
            }
        
        self.__suffix      = suffix

        if   isinstance ( signal_x , PDF            )          : self.__signal_x = signal_x
        elif isinstance ( signal_x , ROOT.RooAbsPdf ) and xvar :
            self.__signal_x = Generic1D_pdf ( signal_1 , xvar , 'SX' )
        else : raise AttributeError ( "Invalid ``signal_x'' argument: %s" % signal_x )
            
        if   isinstance ( signal_y , PDF            )          : self.__signal_y = signal_y
        elif isinstance ( signal_y , ROOT.RooAbsPdf ) and yvar :
            self.__signal_y = Generic1D_pdf ( signal_y , yvar , 'SY' )
        else : raise AttributeError ( "Invalid ``signal_y'' argument: %s" % signal_y )

        if   isinstance ( signal_z , PDF            )          : self.__signal_z = signal_z
        elif isinstance ( signal_z , ROOT.RooAbsPdf ) and zvar :
            self.__signal_z = Generic1D_pdf ( signal_z , zvar , 'SZ' )
        else : raise AttributeError ( "Invalid ``signal_z'' argument: %s" % signal_z )
            
        #
        ## initialize base class
        #
        if not name : 
            name = "%s&%s&%s" % ( self.__signal_x.name ,
                                  self.__signal_y.name ,
                                  self.__signal_z.name )                                  
            if suffix : name += '_'+ suffix
            
        PDF3.__init__ ( self , name ,
                        self.__signal_x.xvar ,
                        self.__signal_y.xvar ,
                        self.__signal_z.xvar ) 
     
        # =====================================================================
        ## 1) First component: all   signals
        # =====================================================================
        
        self.__sss_cmp  = Model3D (
            'SSS_pdf' + suffix , self.__signal_x , self.__signal_y , self.__signal_z )
        
        # =====================================================================
        ## 2-4) Three terms:  ( 2 signals )  x ( 1 background ) 
        # =====================================================================
        
        self.__bkg_1x    = self.make_bkg ( bkg_1x , 'Bkg1X_BSS' + suffix , self.xvar )
        self.__bkg_1y    = self.make_bkg ( bkg_1y , 'Bkg1Y_SBS' + suffix , self.yvar )
        self.__bkg_1z    = self.make_bkg ( bkg_1z , 'Bkg1Z_SSB' + suffix , self.zvar )
        
        self.__ssb_cmp = Model3D ( "SSB_pdf" + suffix ,
                                   self.__signal_x , self.__signal_y , self.__bkg_1z   ,
                                   title = "Signal(x) x Signal(y) x Background1(x)" )
        self.__sbs_cmp = Model3D ( "SBS_pdf" + suffix ,
                                   self.__signal_x , self.__bkg_1y   , self.__signal_z , 
                                   title = "Signal(x) x Background1(y) x Signal(z)" ) 
        self.__bss_cmp = Model3D ( "BSS_pdf" + suffix ,
                                   self.__bkg_1x   , self.__signal_y , self.__signal_z ,
                                   title = "Background1(x) x Signal(y) x Signal(z)" )
        
        # =====================================================================
        ## (intermezzo-1) Assumptions about SBB-background sub-components 
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
            
        if   component_clone   ( bkg_2z ) :
            bkg_2z = self.__bkg_1z
            self.debug ( 'bkg_2z set to [CLONE]   %s' % bkg_2z ) 
        elif component_similar ( bkg_2z ) :
            bkg_2z =        bkg_1z
            self.debug ( 'bkg_2z set to [SIMILAR] %s' % bkg_2z ) 
        # =====================================================================
        
        # =====================================================================
        ## 5-7) Three terms: (1 signal) x (2 backgrounds)
        # =====================================================================
                
        self.__bkg_2x = None 
        self.__bkg_2y = None 
        self.__bkg_2z = None 

        if   bkg_2xy and isinstance ( bkg_2xy , PDF2 ) :
            self.__bkg_2xy = bkg_2xy
        elif bkg_2xy and isinstance ( bkg_2xy , ROOT.RooAbsPdf ) : 
            self.__bkg_2xy = Generic2D_pdf ( bkg_2xy  , self.xvar , self.yvar )
        elif bkg_2xy and isinstance ( bkg_2xy , ( tuple , list ) ) :            
            from ostap.fitting.models_2d import make_B2D
            self.__bkg_2xy = make_B2D ( 'Bkg2XY' + suffix , self.xvar , self.yvar , *bkg_2xy )
        else :

            if not self.__bkg_2x : self.__bkg_2x = self.make_bkg ( bkg_2x , 'Bkg2X_S2B' + suffix , self.xvar )    
            if not self.__bkg_2y : self.__bkg_2y = self.make_bkg ( bkg_2y , 'Bkg2Y_S2B' + suffix , self.yvar )                    
            self.__bkg_2xy = Model2D ( 'Bkg2XY' + suffix ,
                                       self.__bkg_2x     ,
                                       self.__bkg_2y     ,
                                       title =  'Backrgound2(x) x Background2(y)' )
            
        if   bkg_2xz and isinstance ( bkg_2xz , PDF2 ) :
            self.__bkg_2xz = bkg_2xz
        elif bkg_2xz and isinstance ( bkg_2xz , ROOT.RooAbsPdf ) : 
            self.__bkg_2xz = Generic2D_pdf ( bkg_2xz  , self.xvar , self.zvar )
        elif bkg_2xz and isinstance ( bkg_2xz , ( tuple , list ) ) :            
            from ostap.fitting.models_2d import make_B2D
            self.__bkg_2xz = make_B2D ( 'Bkg2XZ' + suffix , self.xvar , self.zvar , *bkg_2xz )
        else :
            if not self.__bkg_2x : self.__bkg_2x = self.make_bkg ( bkg_2x , 'Bkg2X_S2B' + suffix , self.xvar )    
            if not self.__bkg_2z : self.__bkg_2z = self.make_bkg ( bkg_2z , 'Bkg2Z_S2B' + suffix , self.zvar )                                
            self.__bkg_2xz = Model2D ( 'Bkg2XZ' + suffix ,                                     
                                       self.__bkg_2x     ,
                                       self.__bkg_2z       , 
                                       title =  'Backrgound2(x) x Background2(z)' )
            
        if   bkg_2yz and isinstance ( bkg_2yz , PDF2 ) :
            self.__bkg_2yz = bkg_2yz
        elif bkg_2yz and isinstance ( bkg_2yz , ROOT.RooAbsPdf ) : 
            self.__bkg_2yz = Generic2D_pdf ( bkg_2yz  , self.yvar , self.zvar)
        elif bkg_2yz and isinstance ( bkg_2yz , ( tuple , list ) ) :            
            from ostap.fitting.models_2d import make_B2D
            self.__bkg_2yz = make_B2D ( 'Bkg2YZ' + suffix , self.yvar , self.zvar , *bkg_2yz )
        else :            
            if not self.__bkg_2y : self.__bkg_2y = self.make_bkg ( bkg_2y , 'Bkg2Y_S2B' + suffix , self.yvar )    
            if not self.__bkg_2z : self.__bkg_2z = self.make_bkg ( bkg_2z , 'Bkg2Z_S2B' + suffix , self.zvar )                                
            self.__bkg_2yz = Model2D ( 'Bkg2YZ' + suffix ,
                                       self.__bkg_2y     ,
                                       self.__bkg_2z     ,
                                       title =  'Backrgound2(y) x Background2(z)' )

        self.__sbb_cmp = Generic3D_pdf ( 
            ROOT.RooProdPdf ( "SBB_pdf" + suffix , "Signal(x) x Background(y,z)" , self.__signal_x.pdf , self.__bkg_2yz.pdf ) ,
            self.xvar , self.yvar , self.zvar )
        self.__bsb_cmp = Generic3D_pdf (
            ROOT.RooProdPdf ( "BSB_pdf" + suffix , "Signal(y) x Background(x,z)" , self.__signal_y.pdf , self.__bkg_2xz.pdf ) ,
            self.xvar , self.yvar , self.zvar )
        self.__bbs_cmp = Generic3D_pdf ( 
            ROOT.RooProdPdf ( "BBS_pdf" + suffix , "Signal(z) x Background(x,y)" , self.__signal_z.pdf , self.__bkg_2xy.pdf ) ,
            self.xvar , self.yvar , self.zvar )
        

        # =====================================================================
        ## (intermezzo-2) Assumptions about BBB-background sub-components 
        # =====================================================================

        if   component_clone   ( bkg_3x ) :
            bkg_3x = self.__bkg_2x
            self.debug ( 'bkg_3x set to [CLONE]   %s' % bkg_3x ) 
        elif component_similar ( bkg_3x ) :
            bkg_3x =        bkg_2x
            self.debug ( 'bkg_3x set to [SIMILAR] %s' % bkg_3x ) 

        if   component_clone   ( bkg_3y ) :
            bkg_3y = self.__bkg_2y
            self.debug ( 'bkg_3y set to [CLONE]   %s' % bkg_3y  )
        elif component_similar ( bkg_3x ) :
            bkg_3y =        bkg_2y
            self.debug ( 'bkg_3y set to [SIMILAR] %s' % bkg_3y ) 
            
        if   component_clone   ( bkg_3z ) :
            bkg_3z = self.__bkg_2z
            self.debug ( 'bkg_3z set to [CLONE]   %s' % bkg_3z ) 
        elif component_similar ( bkg_3x ) :
            bkg_3z =        bkg_2z
            self.debug ( 'bkg_3z set to [SIMILAR] %s' % bkg_3z ) 


        # =====================================================================
        ## 8) pure background 
        # =====================================================================
        
        self.__bkg_3x = None 
        self.__bkg_3y = None 
        self.__bkg_3z = None 
        
        if   bkg_3D and isinstance ( bkg_3D , PDF3 ) :
            self.__bbb_cmp = bkg_3D
        elif bkg_3D and isinstance ( bkg_3D , ROOT.RooAbsPdf ) : 
            self.__bbb_cmp = Generic3D_pdf ( bkg_3D , self.xvar , self.yvar , self.zvar )
        elif bkg_3D and isinstance ( bkg_3D , (  tuple , list ) ) :
            from ostap.fitting.models_3d import make_B3D 
            self.__bbb_cmp = make_B3D ( 'BBB_pdf' + suffix ,
                                        self.xvar , self.yvar , self.zvar , *bkg_3D )
        else :
            
            self.__bkg_3x = self.make_bkg ( bkg_3x , 'Bkg3X_BBB' + suffix , self.xvar )
            self.__bkg_3y = self.make_bkg ( bkg_3y , 'Bkg3Y_BBB' + suffix , self.yvar )
            self.__bkg_3z = self.make_bkg ( bkg_3z , 'Bkg3Z_BBB' + suffix , self.zvar )
            
            self.__bbb_cmp = Model3D (
                "BBB_pdf" + suffix ,
                self.__bkg_3x ,
                self.__bkg_3y ,
                self.__bkg_3z ,
                title = "Background3(x) x Backrgound3(y) x Background3(z)" )
        #
        ## coefficients
        #
        self.__sss = self.make_var ( sss   , "SSS"          + suffix ,
                               "Signal(x)&Signal(y)&Signal(z)"     + suffix , sss , 1000  , 0 ,  1.e+7 )
        self.__ssb = self.make_var ( ssb   , "SSB"          + suffix ,
                               "Signal(x)&Signal(y)&Background(z)" + suffix , ssb , 1000  , 0 ,  1.e+7 )
        self.__sbs = self.make_var ( sbs   , "SBS"          + suffix ,
                               "Signal(x)&Background(y)&Signal(z)" + suffix , sbs , 1000  , 0 ,  1.e+7 )
        self.__bss = self.make_var ( bss   , "BSS"          + suffix ,
                               "Background(x)&Signal(y)&Signal(z)" + suffix , bss , 1000  , 0 ,  1.e+7 )
        self.__sbb = self.make_var ( sbb  , "SBB"           + suffix ,
                               "Signal(x)&Background(y,z)"         + suffix , sbb , 1000  , 0 ,  1.e+7 )
        self.__bsb = self.make_var ( bsb  , "BSB"           + suffix ,
                               "Signal(y)&Background(x,z)"         + suffix , bsb , 1000  , 0 ,  1.e+7 )
        self.__bbs = self.make_var ( bbs  ,  "BBS"          + suffix ,
                               "Signal(z)&Background(x,y)"         + suffix , bbs , 1000  , 0 ,  1.e+7 )
        self.__bbb = self.make_var ( bbb  , "BBB"           + suffix ,
                               "Background(x,y,z)"                 + suffix , bbb , 1000  , 0 ,  1.e+7 )
        
        self.alist1 = ROOT.RooArgList (
            self.__sss_cmp.pdf ,
            self.__ssb_cmp.pdf ,
            self.__sbs_cmp.pdf ,
            self.__bss_cmp.pdf ,
            self.__sbb_cmp.pdf ,
            self.__bsb_cmp.pdf ,
            self.__bbs_cmp.pdf ,
            self.__bbb_cmp.pdf )
        self.alist2 = ROOT.RooArgList (
            self.__sss     ,
            self.__ssb     ,
            self.__sbs     ,
            self.__bss     ,
            self.__sbb     ,
            self.__bsb     ,
            self.__bbs     ,
            self.__bbb     )

        ## treat additional components (if specified)
        self.__nums_components = [] 
        icmp = 0
        self.__more_components = []
        for cmp in components :
            
            if   isinstance  ( c , PDF2           ) : cc = c  
            elif isinstance  ( c , ROOT.RooAbsPdf ) : cc = Generic2D_pdf ( cs ,  self.xvar , self.yvar ) 
            else :
                self.error ("unknown ``other''component %s/%s, skip it!" % ( cc , type(cc) ) )
                continue  
            self.__more_components.append ( cc     )
            self.components.add           ( cc.pdf ) 

        nc = len( self.__more_components )
        if 1 == nc :
            cf = self.make_var ( None , "C"+suffix , "Component" + suffix , None , 1 , 0 , 1.e+7 )
            self.alist1.add  ( self.components[0] )
            self.__num_components.append ( cf ) 
        elif 2 <= nc : 
            fic = self.make_fracs ( nc , 'C_%%d%s' % suffix ,  'C(%%d)%s'  % suffix , fractions  = False )
            for c in self.components : self.alist1.add ( c)
            for f in fic             : self.__num_components.append ( f )
            
        self.__nums_components  = tuple ( self.__nums_components  ) 
        for c in self.__nums_components  : self.alist2.add ( c )

        ## 
        #
        ## build the final PDF 
        # 
        self.pdf  = ROOT.RooAddPdf  ( "model3D"      + suffix ,
                                      "Model3D(%s)"  % suffix ,
                                      self.alist1 ,
                                      self.alist2 )
        
        self.signals     .add ( self.__sss_cmp.pdf )
        self.backgrounds .add ( self.__bbb_cmp.pdf )
        self.crossterms1 .add ( self.__ssb_cmp.pdf ) ## cross-terms
        self.crossterms1 .add ( self.__sbs_cmp.pdf ) ## cross-terms
        self.crossterms1 .add ( self.__bss_cmp.pdf ) ## cross-terms
        self.crossterms2 .add ( self.__sbb_cmp.pdf ) ## cross-terms 
        self.crossterms2 .add ( self.__bsb_cmp.pdf ) ## cross-terms 
        self.crossterms2 .add ( self.__bbs_cmp.pdf ) ## cross-terms 

        ## save the configuration
        self.config = {
            'signal_x'   : self.signal_x ,
            'signal_y'   : self.signal_y ,
            'signal_z'   : self.signal_z ,
            'suffix'     : self.suffix  ,
            ##
            'bkg_1x'     : self.bkg_1x  ,
            'bkg_1y'     : self.bkg_1y  ,
            'bkg_1z'     : self.bkg_1z  ,
            ##
            'bkg_2x'     : self.bkg_2x  ,
            'bkg_2y'     : self.bkg_2y  ,
            'bkg_2z'     : self.bkg_2z  ,
            ##
            'bkg_2xy'    : self.bkg_2xy ,
            'bkg_2xz'    : self.bkg_2xz ,
            'bkg_2yz'    : self.bkg_2yz ,
            ##
            'bkg_3x'     : self.bkg_3x  ,
            'bkg_3y'     : self.bkg_3y  ,
            'bkg_3z'     : self.bkg_3z  ,
            ##
            'bkg_3D'     : self.bkg_3D  ,
            ##
            'sss'        : self.SSS     ,
            'ssb'        : self.SSB     ,
            'sbs'        : self.SBS     ,
            'bss'        : self.BSS     ,
            'sbb'        : self.SBB     ,
            'bsb'        : self.BSB     ,
            'bbs'        : self.BBS     ,
            'bbb'        : self.BBB     ,
            #
            'components' : self.more_components ,            
            'xvar'       : self.xvar    ,
            'yvar'       : self.yvar    ,
            'zvar'       : self.zvar    ,
            'name'       : self.name    ,             
            }
        
    ## redefine the clone method, allowing only the name to be changed
    #  @attention redefinition of parameters and variables is disabled,
    #             since it can't be done in a safe way                  
    def clone ( self , name = '' , xvar = None , yvar = None , zvar = None  ) :
        """Redefine the clone method, allowing only the name to be changed
         - redefinition of parameters and variables is disabled,
         since it can't be done in a safe way          
        """
        if xvar and not xvar is self.xvar :
            raise AttributeError("Fit3D can not be cloned with different `xvar''")
        if yvar and not yvar is self.yvar :
            raise AttributeError("Fit3D can not be cloned with different `yvar''")
        if zvar and not zvar is self.zvar :
            raise AttributeError("Fit3D can not be cloned with different `zvar''")
        return PDF.clone ( self , name = name ) if name else PDF.clone( self )

    @property
    def SSS ( self ) :
        """The yield of Signal(x)*Signal(y)*Signal(z) component"""
        return self.__sss
    @SSS.setter 
    def SSS ( self , value ) :
        value = float ( value  )
        assert value in self.__sss, "Value %s is out of the allowed range %s " % ( value , self.__sss.minmax() )
        self.__sss.setVal ( value ) 

    @property
    def SSB ( self ) :
        """The yield of Signal(x)*Signal(y)*Background(z) component"""
        return self.__ssb
    @SSB.setter 
    def SSB ( self , value ) :
        value = float ( value  )
        assert value in self.__ssb, "Value %s is out of the allowed range %s " % ( value , self.__ssb.minmax() )
        self.__ssb.setVal ( value ) 

    @property
    def SBS ( self ) :
        """The yield of Signal(x)*Background(y)*Signal(z) component"""
        return self.__sbs
    @SBS.setter 
    def SBS ( self , value ) :
        value = float ( value  )
        assert value in self.__sbs, "Value %s is out of the allowed range %s " % ( value , self.__sbs.minmax() )
        self.__sbs.setVal ( value ) 

    @property
    def BSS ( self ) :
        """The yield of Background(x)*Signal(y)*Signal(z) component"""
        return self.__bss
    @BSS.setter 
    def BSS ( self , value ) :
        value = float ( value  )
        assert value in self.__bss, "Value %s is out of the allowed range %s " % ( value , self.__bss.minmax() )
        self.__bss.setVal ( value ) 

    @property
    def SBB ( self ) :
        """The yield of Signal(x)*Background(y,z) component"""
        return self.__sbb
    @SBB.setter 
    def SBB ( self , value ) :
        value = float ( value  )
        assert value in self.__sbb, "Value %s is out of the allowed range %s " % ( value , self.__sbb.minmax() )
        self.__sbb.setVal ( value ) 

    @property
    def BSB ( self ) :
        """The yield of Background(x,z)*Signal(y) component"""
        return self.__bsb
    @BSB.setter 
    def BSB ( self , value ) :
        value = float ( value  )
        assert value in self.__bsb, "Value %s is out of the allowed range %s " % ( value , self.__bsb.minmax() )
        self.__bsb.setVal ( value ) 

    @property
    def BBS ( self ) :
        """The yield of Background(x,y)*Signal(z) component"""
        return self.__bbs
    @BBS.setter 
    def BBS ( self , value ) :
        value = float ( value  )
        assert value in self.__bbs, "Value %s is out of the allowed range %s " % ( value , self.__bbs.minmax() )
        self.__bbs.setVal ( value ) 

    @property
    def BBB ( self ) :
        """The yield of Background(x,y,z) component"""
        return self.__bbb
    @BBB.setter 
    def BBB ( self , value ) :
        value = float ( value  )
        assert value in self.__bbb, "Value %s is out of the allowed range %s " % ( value , self.__bbb.minmax() )
        self.__bbb.setVal ( value ) 


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
    def signal_z ( self  ) :
        """``signal_z'': Signal(z) component/PDF"""
        return self.__signal_z

    @property
    def bkg_1x ( self ) :
        """``bkg_1x'': B(x) component for B(x)*S(y)*S(z) term"""
        return self.__bkg_1x
    @property
    def bkg_1y ( self ) :
        """``bkg_1y'': B(y) component for S(x)*B(y)*S(z) term"""
        return self.__bkg_1y
    @property
    def bkg_1z ( self ) :
        """``bkg_1z'': B(z) component for S(x)*S(y)*B(z) term"""
        return self.__bkg_1z

    @property
    def bkg_2xy ( self ) :
        """``bkg_2xy'': B(x,y) component for B(x,y)*S(z) term"""
        return self.__bkg_2xy
    @property
    def bkg_2xz ( self ) :
        """``bkg_2xz'': B(x,z) component for B(x,z)*S(y) term"""
        return self.__bkg_2xz
    @property
    def bkg_2yz ( self ) :
        """``bkg_2yz'': B(y,z) component for B(y,z)*S(x) term"""
        return self.__bkg_2yz
    
    @property
    def bkg_2x ( self ) :
        """``bkg_2x'': B(x) component for B(x,y)*S(z) & B(x,z)*S(y) terms"""
        return self.__bkg_2x 
    @property
    def bkg_2y ( self ) :
        """``bkg_2y'': B(y) component for B(y,z)*S(x) & B(x,y)*S(z) terms"""
        return self.__bkg_2y 
    @property
    def bkg_2z ( self ) :
        """``bkg_2z'': B(z) component for B(x,z)*S(y) & B(y,z)*S(x) terms"""
        return self.__bkg_2z 

    @property
    def bkg_3x ( self ) :
        """``bkg_3x'': B(x) component for B(x,y,z) term"""
        return self.__bkg_3x 
    @property
    def bkg_3y ( self ) :
        """``bkg_3y'': B(y) component for B(x,y,z) term"""
        return self.__bkg_3y 
    @property
    def bkg_3z ( self ) :
        """``bkg_3z'': B(z) component for B(z,y,z) term"""
        return self.__bkg_3z 

    @property
    def bkg_3D ( self ) :
        """```bkg_3D'': B(x,y,z) component/PDF for the final PDF"""
        return self.__bbb_cmp 

    @property
    def cmp_SSS ( self ) :
        """```triple-signal'' component/PDF"""
        return self.__sss_cmp

    @property
    def cpm_SSB ( self ) :
        """```signal-signal-background'' component/PDF"""
        return self.__ssb_cmp
    
    @property
    def cmp_SBS ( self ) :
        """```signal-background-signal'' component/PDF"""
        return self.__sbs_cmp
    
    @property
    def cmp_BSS ( self ) :
        """```background-signal-signal'' component/PDF"""
        return self.__bss_cmp

    @property
    def cpm_SBB ( self ) :
        """```signal-background-background'' component/PDF"""
        return self.__sbb_cmp
    
    @property
    def cmp_BSB ( self ) :
        """```background-signal-background'' component/PDF"""
        return self.__bsb_cmp
    
    @property
    def cmp_BBS ( self ) :
        """```background-background-signal'' component/PDF"""
        return self.__bbs_cmp

    @property
    def cmp_BBB ( self ) :
        """```triple-background'' component/PDF"""
        return self.__bbb_cmp

    @property
    def more_components ( self ) :
        """additional/``other'' components"""
        return tuple( self.__more_components  )

    @property
    def suffix ( self ) :
        """``suffix'', used to build the name"""
        return self.__suffix



# =============================================================================
## @class Fit3DSym
#  The actual model for fully symmetric 3D-fits
#
#  @param signal_x (RooFit/PDF or Ostap/PDF) PDF to describe (1D)-signal in X-direction
#  @param signal_y (RooFit/PDF, Ostap/PDF or None) PDF to describe (1D)-signal in Y-direction
#  @param signal_z (RooFit/PDF, Ostap/PDF or None) PDF to describe (1D)-signal in Z-direction
#  @param suffix   (string) An optional suffix to be  added to the names of created PDFs and variables
#  @param bkg_1x   (RooFit/PDF, Ostap/PDF,  integer, RooRealVar or None)
#                  1D x-background for SSB-terms
#                  Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
#  @param bkg_2x   (RooFit/PDF, Ostap/PDF,integer, RooRealVar or None)
#                  1D x-background for SBB-terms, if <code>bkg2D</code> is not specified. 
#                  Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
#  @param bkg_2xy  (RooFit/PDF or Ostap/PDF)
#                   2D x,y-background for SBB-terms
#                   Use directly RooFit/PDF or Ostap/PDF
#  @param bkg_3x   (RooFit/PDF, Ostap/PDF,integer, RooRealVar or None)
#                   1D x-background for BBB-term, if <code>bkg3D</code> is not specified. 
#                   Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
#  @param bkg_3D   (RooFit/PDF or Ostap/PDF)
#                   3D x,y,z-background for BBB-term
#                   Use directly RooFit/PDF or Ostap/PDF        
#  @param sss      (None, RooRealVar, non-negative float or tuple)
#                   Variable for the yield of SSS component.
#                   Use directly RooRelaVar, otherwise create it using self.make_var function
#  @param ssb      (None, RooRealVar, non-negative float or tuple)
#                   Variable for the yield of SSB component.
#                   Use directly RooRelaVar, otherwise create it using self.make_var function
#  @param sbb      (None, RooRealVar, non-negative float or tuple)
#                  Variable for the yield of SBB component.
#                  Use directly RooRelaVar, otherwise create it using self.make_var function
#  @param bbb      (None, RooRealVar, non-negative float or tuple)
#                  Variable for the yield of BBB component.
#                  Use directly RooRelaVar, otherwise create it using self.make_var function
#  @param components ([])   the list of additional 3D-PDFs to be used in the fit
#  @param xvar       (None) the x-variable
#  @param yvar       (None) the y-variable
#  @param zvar       (None) the z-variable
#  @param name       ("")   the PDF name 
#  @code
#  model   = Models.Fit3DSym (
#      signal_x = Models.Gauss_pdf ( 'Gx' , mass = m_x ) ,
#      signal_y = Models.Gauss_pdf ( 'Gy' , mass = m_y ) ,
#      signal_z = Models.Gauss_pdf ( 'Gz' , mass = m_z ) ,
#      bkg_1x   = -1 , 
#      bkg_2x   = -1 ,
#      bkg_3x   = -1 )
#
#  r = model.fitTo ( dataset ) ## fit dataset 
#
#  print r                       ## get results  
#
#  fx  = model.draw1 ()          ## visualize X-projection
#  fy  = model.draw2 ()          ## visualize Y-projection
#  fz  = model.draw3 ()          ## visualize Z-projection
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2017-07-25
class Fit3DSym (PDF3) :
    """The actual model for fully symmetric 3D-fits
    
    >>>  model   = Models.Fit3DSym (
    ...      signal_x = Models.Gauss_pdf ( 'Gx' , mass = m_x ) ,
    ...      signal_y = Models.Gauss_pdf ( 'Gy' , mass = m_y ) ,
    ...      signal_z = Models.Gauss_pdf ( 'Gz' , mass = m_z ) ,
    ...      bkg_1x   =  1 , 
    ...      bkg_2x   = -1 ,
    ...      bkg_3x   = -1 )
    >>> r,f = model.fitTo ( dataset ) ## fit dataset 
    >>> print r                       ## get results  
    >>> fx  = model.draw1 ()          ## visualize X-projection
    >>> fy  = model.draw2 ()          ## visualize Y-projection
    >>> fz  = model.draw3 ()          ## visualize Z-projection

    Parameters
    ----------
    signal_x :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in X-direction
    signal_y :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in Y-direction
    signal_z :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in Z-direction
    suffix   : string
        An optional suffix to be  added to the names of created PDFs and variables
    bkg_1x   : RooFit/PDF, Ostap/PDF, integer, RooRealVar or None
        1D x-background for SSB-terms
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
    bkg_2x  : RooFit/PDF, Ostap/PDF, integer, RooRealVar or None
        1D x-background for SBB-terms
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
    bkg_2xy : RooFit/PDF, Ostap/PDF, list/tuple or None 
        2D (x,y)-background for SBB-terms
        Use directly RooFit/PDF or Ostap/PDF
    bkg_3x  : RooFit/PDF, Ostap/PDF, integer integer, RooRealVar or None
        1D x-background for BBB-term
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
    bkg_3D  : RooFit/PDF, Ostap/PDF, list/tuple or None 
        3D x,y,z-background for BBB-term
        Use directly RooFit/PDF or Ostap/PDF
        
    sss      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of SSS component
         Use directly RooRelaVar, otherwise create it using self.make_var function
    ssb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of SSB component
         Use directly RooRelaVar, otherwise create it using self.make_var function
    sbb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of SBB component
         Use directly RooRelaVar, otherwise create it using self.make_var function
    bbb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of BBB component
         Use directly RooRelaVar, otherwise create it using self.make_var function
         
    components : the list of additional 3D-PDFs to be used in the fit
    xvar       : x-variable
    yvar       : y-variable
    zvar       : z-variable
    name       : the PDF name 
         
    """
    def __init__ ( self               ,
                   #
                   signal_x           , 
                   signal_y = None    ,
                   signal_z = None    ,
                   suffix   = ''      ,
                   # background  for SSB-terms:
                   bkg_1x     = None  , ## 1D     x-background for BSS-terms
                   # background for  SBB-terms:
                   bkg_2x     = None  , ## 1D     x-background for BBS-terms, if no bkg2d is specified   
                   bkg_2xy    = None  , ## 2D   x,y-background for BBS-terms     (symmetric)
                   # background for BBB-term
                   bkg_3x     = None  , ## 1D     x-background for BBB-term
                   bkg_3D     = None  , ## 3D x,y,z-background for B(x,y,z) term (symmetric) 
                   #
                   ## Yields of the main components :
                   sss        = None  , ## sig(1) * sig(2) * sig(3) 
                   ssb        = None  , ## SSB components 
                   sbb        = None  , ## SBB components 
                   bbb        = None  , ## background-3D 
                   ## additional components 
                   components = []    ,
                   xvar       = None  ,
                   yvar       = None  ,
                   zvar       = None  ,                   
                   name       = ''    ) : 
        
        ## keep all arguments 
        self.__args = {
            #
            'signal_x'   : signal_x   ,
            'signal_y'   : signal_y   ,
            'signal_z'   : signal_z   ,
            #
            'bkg_1x'     : bkg_1x     ,
            'bkg_2x'     : bkg_2x     ,
            'bkg_2xy'    : bkg_2xy    ,
            'bkg_3x'     : bkg_3x     ,
            'bkg_3D'     : bkg_3D     ,
            ##
            'sss'        : sss        ,            
            'ssb'        : ssb        ,
            'sbb'        : sbb        ,
            'bbb'        : bbb        ,
            ##
            'components' : components ,
            ##
            'xvar'       : xvar       ,
            'yvar'       : yvar       ,
            'zvar'       : zvar       ,
            ##
            'name'       : name                   
            }
        
        self.__suffix      = suffix
 
           
        if   isinstance ( signal_x , PDF            )          : self.__signal_x= signal_x
        elif isinstance ( signal_x , ROOT.RooAbsPdf ) and xvar :
            self.__signal_x = Generic1D_pdf ( signal_x , xvar , 'SX' )
        else : raise AttributeError ( "Invalid ``signal_x'' argument: %s" % signal_x )
            
        if   isinstance ( signal_y , PDF            )          : self.__signal_y = signal_y
        elif isinstance ( signal_y , ROOT.RooAbsPdf ) and yvar :
            self.__signal_y = Generic1D_pdf ( signal_y , yvar , 'SY' )
        elif yvar and not signal_y :
            self.__signal_y = self.__signal_x.clone ( xvar = yvar , name = 'SY' )
            self.debug('signal y-component is cloned from the signal_x component')
        else : raise AttributeError ( "Invalid ``signal_y'' argument: %s" % signal_y )

        if   isinstance ( signal_z , PDF            )          : self.__signal_z = signal_z
        elif isinstance ( signal_z , ROOT.RooAbsPdf ) and zvar :
            self.__signal_z = Generic1D_pdf ( signal_z , zvar , 'SZ' )
        elif zvar and not signal_z :
            self.__signal_z = self.__signal_x.clone ( xvar = zvar , name = 'SZ' )
            self.debug('signal z-component is cloned from the signal_x component')
        else : raise AttributeError ( "Invalid ``signal_z'' argument: %s" % signal_z )
            
        #
        ## initialize base class
        #
        if not name :
            name = "%s&%s&%s" % ( self.__signal_x.name ,
                                  self.__signal_y.name ,
                                  self.__signal_z.name )                                  
            if suffix : name += '_'+ suffix
            
        PDF3.__init__ ( self , name ,
                        self.__signal_x.xvar ,
                        self.__signal_y.xvar ,
                        self.__signal_z.xvar ) 
     

        # =====================================================================
        ## 1) First component: all   signals
        # =====================================================================
        
        self.__sss_cmp  = Model3D (
            'SSS_pdf' + suffix , self.__signal_x , self.__signal_y , self.__signal_z )
        
        # =====================================================================
        ## 2) ( 2 signals )  x ( 1 background )
        # =====================================================================
        
        self.__bkg_1x = self.make_bkg (        bkg_1x , 'Bkg1X_BSS' + suffix , self.xvar )
        self.__bkg_1y = self.make_bkg ( self.__bkg_1x , 'Bkg1Y_SBS' + suffix , self.yvar )
        self.__bkg_1z = self.make_bkg ( self.__bkg_1x , 'Bkg1Z_SSB' + suffix , self.zvar )
        
        self.__ssb_cmp_raw = Model3D ( "SSB_raw" + suffix ,
                                       self.__signal_x , self.__signal_y , self.__bkg_1z   ,
                                       title = "Signal(x) x Signal(y) x Background1(z)" )
        self.__sbs_cmp_raw = Model3D ( "SBS_raw" + suffix ,
                                       self.__signal_x , self.__bkg_1y   , self.__signal_z , 
                                       title = "Signal(x) x Background1(y) x Signal(z)" ) 
        self.__bss_cmp_raw = Model3D ( "BSS_raw" + suffix ,
                                       self.__bkg_1x   , self.__signal_y , self.__signal_z ,
                                       title = "Background1(x) x Signal(y) x Signal(z)" )
        
        self.__ssb_cmp     = Generic3D_pdf (
            self.make_sum ( "SSB_pdf" + suffix ,
                            "S(x)*S(y)*B(z)+S(x)*B(y)*S(z)+B(x)*S(y)*B(z)" ,
                            self.__ssb_cmp_raw.pdf ,
                            self.__sbs_cmp_raw.pdf ,
                            self.__bss_cmp_raw.pdf ) , 
            self.xvar , self.yvar , self.zvar )
        
        self.__sbs_cmp = self.__ssb_cmp
        self.__bss_cmp = self.__ssb_cmp
        

        # =====================================================================
        ## (intermezzo-1) Assumptions about the SBB-background sub-components 
        # =====================================================================        
        
        if   component_clone   ( bkg_2x ) :
            bkg_2x = self.__bkg_1x
            self.debug ( 'bkg_2x set to [CLONE]   %s' % bkg_2x ) 
        elif component_similar ( bkg_2x ) :
            bkg_2x =        bkg_1x
            self.debug ( 'bkg_2x set to [SIMILAR] %s' % bkg_2x ) 
        # =====================================================================

        
        # =====================================================================
        ## 3 Three terms: (1 signal) x (2 backgrounds)
        # =====================================================================
        
        self.__bkg_2x = None 
        self.__bkg_2y = None 
        self.__bkg_2z = None 

        if   bkg_2xy and isinstance ( bkg_2xy , PDF2 ) :
            
            self.__bkg_2xy = bkg_2xy
            self.__bkg_2xz = bkg_2xy.clone ( xvar = self.xvar , yvar = self.zvar , name = 'Bkg2XZ_pdf' + suffix ) 
            self.__bkg_3yz = bkg_2xy.clone ( xvar = self.yvar , yvar = self.yvar , name = 'Bkg2YZ_pdf' + suffix )
            
        elif bkg_2xy and isinstance ( bkg_2xy , ROOT.RooAbsPdf ) :
            
            self.__bkg_2xy = Generic2D_pdf ( bkg_2xy  , self.xvar , self.yvar )
            self.__bkg_2xz = self.__bkg_2xy.clone ( name = 'Bkg2XZ_pdf' + suffix ,
                                                    xvar = self.xvar             ,
                                                    yvar = self.zvar             )
            self.__bkg_2xz = self.__bkg_2xy.clone ( name = 'Bkg2YZ_pdf' + suffix ,
                                                    xvar = self.yvar             ,
                                                    yvar = self.zvar             )
        elif bkg_2xy and isinstance ( bkg_2xy , ( tuple , list ) ) :

            from ostap.fitting.models_2d import make_B2Dsym 
            self.__bkg_2xy = make_B2Dsym ( 'Bkg2XY_pdf' , self.xvar , self.yvar  , *bkg_2xy )
            self.__bkg_2xz = self.__bkg_2xy.clone ( name = 'Bkg2XZ_pdf' + suffix ,
                                                    xvar = self.xvar             ,
                                                    yvar = self.zvar             )
            self.__bkg_2xz = self.__bkg_2xy.clone ( name = 'Bkg2YZ_pdf' + suffix ,
                                                    xvar = self.yvar             ,
                                                    yvar = self.zvar             )              
        else :

            self.__bkg_2x  = self.make_bkg (        bkg_2x , 'Bkg2X_S2B' + suffix , self.xvar )        
            self.__bkg_2y  = self.make_bkg ( self.__bkg_2x , 'Bkg2Y_S2B' + suffix , self.yvar )        
            self.__bkg_2z  = self.make_bkg ( self.__bkg_2x , 'Bkg2Z_S2B' + suffix , self.zvar )
            
            self.__bkg_2xy = Model2D ( 'Bkg2XY' + suffix ,
                                       self.__bkg_2x     ,
                                       self.__bkg_2y     , title =  'Background2(x,y)' )
            self.__bkg_2xz = Model2D ( 'Bkg2XZ' + suffix ,
                                       self.__bkg_2x     ,
                                       self.__bkg_2z     , title =  'Background2(x,z)' )
            self.__bkg_2yz = Model2D ( 'Bkg2YZ' + suffix ,
                                       self.__bkg_2y     ,
                                       self.__bkg_2z     , title =  'Background2(y,z)' )
            
        self.__sbb_cmp_raw = Generic3D_pdf ( 
            ROOT.RooProdPdf ( "SBB_raw" + suffix , "Signal(x) x Background2(y,z)" , self.__signal_x.pdf , self.__bkg_2yz.pdf ) ,
            self.xvar , self.yvar , self.zvar )
        self.__bsb_cmp_raw = Generic3D_pdf (
            ROOT.RooProdPdf ( "BSB_raw" + suffix , "Signal(y) x Background2(x,z)" , self.__signal_y.pdf , self.__bkg_2xz.pdf ) ,
            self.xvar , self.yvar , self.zvar )
        self.__bbs_cmp_raw = Generic3D_pdf ( 
            ROOT.RooProdPdf ( "BBS_raw" + suffix , "Signal(z) x Background2(x,y)" , self.__signal_z.pdf , self.__bkg_2xy.pdf ) ,
            self.xvar , self.yvar , self.zvar )

        self.__sbb_cmp     = Generic3D_pdf (
            self.make_sum( "SBB_pdf" + suffix ,
                           "S(x)*B(y,z)+S(y)*B(x,z)+S(Z)*B(x,y)" ,
                           self.__sbb_cmp_raw.pdf ,
                           self.__bsb_cmp_raw.pdf ,
                           self.__bbs_cmp_raw.pdf ) ,
            self.xvar , self.yvar , self.zvar )

        self.__bsb_cmp = self.__sbb_cmp
        self.__bbs_cmp = self.__sbb_cmp
        
        # =====================================================================
        ## (intermezzo-2) Assumptions about the BBB-background sub-components 
        # =====================================================================        
        if   component_clone   ( bkg_3x ) :
            bkg_3x = self.__bkg_2x
            self.debug ( 'bkg_3x set to [CLONE]   %s' % bkg_3x ) 
        elif component_similar ( bkg_3x ) :
            bkg_3x =        bkg_2x
            self.debug ( 'bkg_3x set to [SIMILAR] %s' % bkg_3x ) 
        # =====================================================================

        # =====================================================================
        ## 8) pure background 
        # =====================================================================
        
        self.__bkg_3x  = None 
        self.__bkg_3y  = None 
        self.__bkg_3z  = None 

        if   bkg_3D and isinstance ( bkg_3D , PDF3 ) :
            self.__bbb_cmp = bkg_3D
        elif bkg_3D and isinstance ( bkg_3D , ROOT.RooAbsPdf ) : 
            self.__bbb_cmp = Generic3D_pdf ( bkg_3D , self.xvar , self.yvar , self.zvar )
        elif bkg_3D and isinstance ( bkg_3D , (  tuple , list )  ) :
            from ostap.fitting.models_3d import make_B3Dsym
            self.__bbb_cmp = make_B3Dsym ( 'BBB_pdf' , self.xvar , self.yvar , self.zvar , *bkg_3D )
        else :
            
            self.__bkg_3x  = self.make_bkg (        bkg_3x , 'Bkg3X_BBB' + suffix , self.xvar )        
            self.__bkg_3y  = self.make_bkg ( self.__bkg_3x , 'Bkg3Y_BBB' + suffix , self.yvar )        
            self.__bkg_3z  = self.make_bkg ( self.__bkg_3x , 'Bkg3Z_BBB' + suffix , self.zvar )
            
            self.__bbb_cmp = Model3D ( "BBB_pdf" + suffix ,
                                       self.__bkg_3x      ,
                                       self.__bkg_3y      ,
                                       self.__bkg_3z      , title = "Background(x,y,z)" )
        #
        ## coefficients
        #
        self.__sss = self.make_var ( sss   , "SSS" + suffix ,
                               "Signal(x)&Signal(y)&Signal(z)" + suffix , sss , 1000 , 0 , 1.e+7 )
        self.__ssb = self.make_var ( ssb   , "SSB" + suffix ,
                               "Signal*2&Background"           + suffix , ssb , 1000 , 0 , 1.e+7 )
        self.__sbb = self.make_var ( sbb   , "SBB" + suffix ,
                               "Signal&Backrgound*2"           + suffix , sbb , 1000 , 0 , 1.e+7 )
        self.__bbb = self.make_var ( bbb  , "BBB"  + suffix ,
                               "Background*3"                  + suffix , bbb , 1000 , 0 , 1.e+7 )
        
        self.__sbs = self.__ssb ## the same 
        self.__bss = self.__ssb ## the same 
        self.__bsb = self.__sbb ## the same 
        self.__bbs = self.__sbb ## the same         

        self.alist1 = ROOT.RooArgList (
            self.__sss_cmp.pdf ,
            self.__ssb_cmp.pdf ,
            self.__sbb_cmp.pdf ,
            self.__bbb_cmp.pdf )
        self.alist2 = ROOT.RooArgList (
            self.__sss     ,
            self.__ssb     ,
            self.__sbb     ,
            self.__bbb     )

        ## treat additional components (if specified)
        self.__nums_components = [] 
        icmp = 0
        self.__more_components = []
        for cmp in components :
            
            if   isinstance  ( c , PDF2           ) : cc = c  
            elif isinstance  ( c , ROOT.RooAbsPdf ) : cc = Generic2D_pdf ( cs ,  self.xvar , self.yvar ) 
            else :
                self.error ("unknown ``other''component %s/%s, skip it!" % ( cc , type(cc) ) )
                continue  
            self.__more_components.append ( cc     )
            self.components.add           ( cc.pdf ) 

        nc = len( self.__more_components )
        if 1 == nc :
            cf = self.make_var ( None , "C"+suffix , "Component" + suffix , None , 1 , 0 , 1.e+7 )
            self.alist1.add  ( self.components[0] )
            self.__num_components.append ( cf ) 
        elif 2 <= nc : 
            fic = self.make_fracs ( nc , 'C_%%d%s' % suffix ,  'C(%%d)%s'  % suffix , fractions  = False )
            for c in self.components : self.alist1.add ( c)
            for f in fic             : self.__num_components.append ( f )
            
        self.__nums_components  = tuple ( self.__nums_components  ) 
        for c in self.__nums_components  : self.alist2.add ( c )

        ## 
        #
        ## build the final PDF 
        # 
        self.pdf  = ROOT.RooAddPdf  ( "model3D"      + suffix ,
                                      "Model3D(%s)"  % suffix ,
                                      self.alist1 ,
                                      self.alist2 )

        self.signals     .add ( self.__sss_cmp.pdf )
        self.backgrounds .add ( self.__bbb_cmp.pdf )
        self.crossterms1 .add ( self.__ssb_cmp.pdf ) ## cross-terms
        self.crossterms2 .add ( self.__sbb_cmp.pdf ) ## cross-terms 
        
        ## save the configuration
        self.config = {
            ##
            'signal_x'   : self.signal_x ,
            'signal_y'   : self.signal_y ,
            'signal_z'   : self.signal_z ,
            'suffix'     : self.suffix  ,
            ##
            'bkg_1x'     : self.bkg_1x  ,
            'bkg_2x'     : self.bkg_2x  ,
            'bkg_2xy'    : self.bkg_2xy ,
            'bkg_3x'     : self.bkg_3x  ,
            'bkg_3D'     : self.bkg_3D  ,
            ##
            'sss'        : self.SSS     ,
            'ssb'        : self.SSB     ,
            'sbb'        : self.SBB     ,
            'bbb'        : self.BBB     ,
            ##
            'components' : self.more_components ,
            ##
            'xvar'       : self.xvar    ,
            'yvar'       : self.yvar    ,
            'zvar'       : self.zvar    ,
            ##
            'name'       : self.name    ,             
            }
        
    ## redefine the clone method, allowing only the name to be changed
    #  @attention redefinition of parameters and variables is disabled,
    #             since it can't be done in a safe way                  
    def clone ( self , name = '' , xvar = None , yvar = None , zvar = None  ) :
        """Redefine the clone method, allowing only the name to be changed
         - redefinition of parameters and variables is disabled,
         since it can't be done in a safe way          
        """
        if xvar and not xvar is self.xvar :
            raise AttributeError("Fit3DSym can not be cloned with different `xvar''")
        if yvar and not yvar is self.yvar :
            raise AttributeError("Fit3DSym can not be cloned with different `yvar''")
        if zvar and not zvar is self.zvar :
            raise AttributeError("Fit3DSym can not be cloned with different `zvar''")
        return PDF.clone ( self , name = name ) if name else PDF.clone( self )

    @property
    def SSS ( self ) :
        """The yield of Signal(x)*Signal(y)*Signal(z) component"""
        return self.__sss
    @SSS.setter 
    def SSS ( self , value ) :
        value = float ( value  )
        assert value in self.__sss, "Value %s is out of the allowed range %s " % ( value , self.__sss.minmax() )
        self.__sss.setVal ( value ) 

    @property
    def SSB ( self ) :
        """The yield of S(x)*S(y)*B(z)+S(x)*B(y)*S(z)+B(x)*S(y)*S(z) component (same as SBS, BSS)"""
        return self.__ssb
    @SSB.setter 
    def SSB ( self , value ) :
        value = float ( value  )
        assert value in self.__ssb, "Value %s is out of the allowed range %s " % ( value , self.__ssb.minmax() )
        self.__ssb.setVal ( value ) 

    @property
    def SBS ( self ) :
        """The yield of S(x)*S(y)*B(z)+S(x)*B(y)*S(z)+B(x)*S(y)*S(z) component (same as SSB, BSS)"""
        return self.__sbs
    @SBS.setter 
    def SBS ( self , value ) :
        value = float ( value  )
        assert value in self.__sbs, "Value %s is out of the allowed range %s " % ( value , self.__sbs.minmax() )
        self.__sbs.setVal ( value ) 

    @property
    def BSS ( self ) :
        """The yield of S(x)*S(y)*B(z)+S(x)*B(y)*S(z)+B(x)*S(y)*S(z) component (same as SSB, SBS)"""
        return self.__bss
    @BSS.setter 
    def BSS ( self , value ) :
        value = float ( value  )
        assert value in self.__bss, "Value %s is out of the allowed range %s " % ( value , self.__bss.minmax() )
        self.__bss.setVal ( value ) 

    @property
    def SBB ( self ) :
        """The yield of S(x)*B(y,z)+S(y)*B(x,z)+S(z)*B(y,z) component (same as BSB, BBS)"""
        return self.__sbb
    @SBB.setter 
    def SBB ( self , value ) :
        value = float ( value  )
        assert value in self.__sbb, "Value %s is out of the allowed range %s " % ( value , self.__sbb.minmax() )
        self.__sbb.setVal ( value ) 

    @property
    def BSB ( self ) :
        """The yield of S(x)*B(y,z)+S(y)*B(x,z)+S(z)*B(y,z) component (same as SBB, BBS)"""
        return self.__bsb
    @BSB.setter 
    def BSB ( self , value ) :
        value = float ( value  )
        assert value in self.__bsb, "Value %s is out of the allowed range %s " % ( value , self.__bsb.minmax() )
        self.__bsb.setVal ( value ) 

    @property
    def BBS ( self ) :
        """The yield of S(x)*B(y,z)+S(y)*B(x,z)+S(z)*B(y,z) component (same as SBB, BSB )"""
        return self.__bbs
    @BBS.setter 
    def BBS ( self , value ) :
        value = float ( value  )
        assert value in self.__bbs, "Value %s is out of the allowed range %s " % ( value , self.__bbs.minmax() )
        self.__bbs.setVal ( value ) 

    @property
    def BBB ( self ) :
        """The yield of Background(x,y,z) component"""
        return self.__bbb
    @BBB.setter 
    def BBB ( self , value ) :
        value = float ( value  )
        assert value in self.__bbb, "Value %s is out of the allowed range %s " % ( value , self.__bbb.minmax() )
        self.__bbb.setVal ( value ) 

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
    def signal_z ( self  ) :
        """``signal_z'': Signal(z) component/PDF"""
        return self.__signal_z

    @property
    def bkg_1x ( self ) :
        """``bkg_1x'': B(x) component for B(x)*S(y)*S(z) term"""
        return self.__bkg_1x
    @property
    def bkg_1y ( self ) :
        """``bkg_1y'': B(y) component for S(x)*B(y)*S(z) term"""
        return self.__bkg_1y
    @property
    def bkg_1z ( self ) :
        """``bkg_1z'': B(z) component for S(x)*S(y)*B(z) term"""
        return self.__bkg_1z
    
    @property
    def bkg_2xy( self ) :
        """``bkg_2xy'': B(x,y) component for B(x,y)*S(z) term"""
        return self.__bkg_2xy
    @property
    def bkg_2xz( self ) :
        """``bkg_2xz'': B(x,z) component for B(x,z)*S(y) term"""
        return self.__bkg_2xz
    @property
    def bkg_2yz( self ) :
        """``bkg_2yz'': B(y,z) component for B(y,z)*S(x) term"""
        return self.__bkg_2yz 

    @property
    def bkg_2x ( self ) :
        """``bkg_2x'': B(x) component for B(x,y)*S(z) & B(x,z)*S(y) terms"""
        return self.__bkg_2x
    @property
    def bkg_2y ( self ) :
        """``bkg_2y'': B(y) component for B(y,z)*S(x) & B(x,y)*S(z) terms"""
        return self.__bkg_2y
    @property
    def bkg_2z ( self ) :
        """``bkg_2z'': B(z) component for B(x,z)*S(y) & B(y,z)*S(x) terms"""
        return self.__bkg_2z

    @property
    def bkg_3x ( self ) :
        """``bkg_3x'': B(x) component for B(x)*B(y)*B(z) term"""
        return self.__bkg_3x
    @property
    def bkg_3y ( self ) :
        """``bkg_3y'': B(y) component for B(x)*B(y)*B(z) term"""
        return self.__bkg_3y
    @property
    def bkg_3z ( self ) :
        """``bkg_3z'': B(z) component for B(x)*B(y)*B(z) term"""
        return self.__bkg_3z

    @property
    def bkg_3D ( self ) :
        """``bkg_3D'': B(x,y,z) component/PDF for the final PDF"""
        return self.__bbb_cmp 

    @property
    def cmp_SSS ( self ) :
        """```triple-signal'' component/PDF"""
        return self.__sss_cmp

    @property
    def cpm_SSB ( self ) :
        """```signal-signal-background'' symmetrized component/PDF"""
        return self.__ssb_cmp
    
    @property
    def cmp_SBS ( self ) :
        """```signal-background-signal'' symmetrized component/PDF (same as cmp_SSB)"""
        return self.__sbs_cmp
    
    @property
    def cmp_BSS ( self ) :
        """```background-signal-signal'' symmetrized component/PDF (same as cmp_SSB)"""
        return self.__bss_cmp

    @property
    def cpm_SBB ( self ) :
        """```signal-background-background'' symmetrized component/PDF"""
        return self.__sbb_cmp
    
    @property
    def cmp_BSB ( self ) :
        """```background-signal-background'' symmetrized component/PDF (same as cmp_SBB)"""
        return self.__bsb_cmp
    
    @property
    def cmp_BBS ( self ) :
        """```background-background-signal'' symmetrized component/PDF (same as cmp_SBB)"""
        return self.__bbs_cmp

    @property
    def cmp_BBB ( self ) :
        """```triple-background'' component/PDF"""
        return self.__bbb_cmp

    @property
    def more_components ( self ) :
        """additional/``other'' components"""
        return tuple( self.__more_components  )

    @property
    def suffix ( self ) :
        """``suffix'', used to build the name"""
        return self.__suffix

    # =========================================================================
    ## Raw, non-symmetrized fit components/PDF (for debugging)
    # =========================================================================

    @property
    def cpm_raw_SSB ( self ) :
        """```signal-signal-background'' raw,non-symmetrized component/PDF"""
        return self.__ssb_cmp_raw
    
    @property
    def cmp_raw_SBS ( self ) :
        """```signal-background-signal'' raw,non-symmetrized component/PDF"""
        return self.__sbs_cmp_raw
    
    @property
    def cmp_raw_BSS ( self ) :
        """```background-signal-signal'' raw,non-symmetrized component/PDF"""
        return self.__bss_cmp_raw

    @property
    def cpm_raw_SBB ( self ) :
        """```signal-background-background'' raw,non-symmetrized component/PDF"""
        return self.__sbb_cmp_raw
    
    @property
    def cmp_raw_BSB ( self ) :
        """```background-signal-background'' raw,non-symmetrized component/PDF"""
        return self.__bsb_cmp_raw
    
    @property
    def cmp_raw_BBS ( self ) :
        """```background-background-signal'' raw,non-symmetrized component/PDF"""
        return self.__bbs_cmp_raw


# =============================================================================
## @class Fit3DMix
#  The actual model for fully "mixed-symemtry" 3D-fits  (symmetric for y<-->z)
#
#  @param signal_x (RooFit/PDF or Ostap/PDF) PDF to describe (1D)-signal in X-direction
#  @param signal_y (RooFit/PDF, Ostap/PDF or None) PDF to describe (1D)-signal in Y-direction
#  @param signal_z (RooFit/PDF, Ostap/PDF or None) PDF to describe (1D)-signal in Z-direction
#  @param suffix   (string) An optional suffix to be  added to the names of created PDFs and variables
#  @param bkg_1x   (RooFit/PDF, Ostap/PDF,  integer, RooRealVar or None)
#                  1D x-background for SSB-terms
#                  Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
#  @param bkg_1y   (RooFit/PDF, Ostap/PDF,integer, RooRealVar or None)
#                  1D x-background for SBB-terms, if <code>bkg2D</code> is not specified. 
#                  Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
#  @param bkg_2x   (RooFit/PDF or Ostap/PDF)
#                   2D x,y-background for SBB-terms
#                   Use directly RooFit/PDF or Ostap/PDF
#  @param bkg3     (RooFit/PDF, Ostap/PDF,integer, RooRealVar or None)
#                   1D x-background for BBB-term, if <code>bkg3D</code> is not specified. 
#                   Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
#  @param bkg3D    (RooFit/PDF or Ostap/PDF)
#                   3D x,y,z-background for BBB-term
#                   Use directly RooFit/PDF or Ostap/PDF        
#  @param sss      (None, RooRealVar, non-negative float or tuple)
#                   Variable for the yield of SSS component.
#                   Use directly RooRelaVar, otherwise create it using self.make_var function
#  @param ssb      (None, RooRealVar, non-negative float or tuple)
#                   Variable for the yield of SSB component.
#                   Use directly RooRelaVar, otherwise create it using self.make_var function
#  @param sbb      (None, RooRealVar, non-negative float or tuple)
#                  Variable for the yield of SBB component.
#                  Use directly RooRelaVar, otherwise create it using self.make_var function
#  @param bbb      (None, RooRealVar, non-negative float or tuple)
#                  Variable for the yield of BBB component.
#                  Use directly RooRelaVar, otherwise create it using self.make_var function
#  @param components ([])   the list of additional 3D-PDFs to be used in the fit
#  @param xvar       (None) the x-variable
#  @param yvar       (None) the y-variable
#  @param zvar       (None) the z-variable
#  @param name       ("")   the PDF name 
#  @code
#  model   = Models.Fit3DSym (
#      signal_x = Models.Gauss_pdf ( 'Gx' , mass = m_x ) ,
#      signal_y = Models.Gauss_pdf ( 'Gy' , mass = m_y ) ,
#      signal_z = Models.Gauss_pdf ( 'Gz' , mass = m_z ) ,
#      bkg_1x     = -1 , 
#      bkg_2x     = -1 ,
#      bkg_1y     = -1 )
#
#  r = model.fitTo ( dataset ) ## fit dataset 
#
#  print r                       ## get results  
#
#  fx  = model.draw1 ()          ## visualize X-projection
#  fy  = model.draw2 ()          ## visualize Y-projection
#  fz  = model.draw3 ()          ## visualize Z-projection
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2017-07-25
class Fit3DMix (PDF3) :
    """The actual model for y<->z-symmetric 3D-fits
    
    >>>  model   = Models.Fit3DSym (
    ...      signal_x = Models.Gauss_pdf ( 'Gx' , mass = m_x ) ,
    ...      signal_y = Models.Gauss_pdf ( 'Gy' , mass = m_y ) ,
    ...      signal_z = Models.Gauss_pdf ( 'Gz' , mass = m_z ) ,
    ...      bkg_1x   =  1 , 
    ...      bkg_2x   = -1 ,
    ...      bkg_1y   = -1 )
    >>> r,f = model.fitTo ( dataset ) ## fit dataset 
    >>> print r                       ## get results  
    >>> fx  = model.draw1 ()          ## visualize X-projection
    >>> fy  = model.draw2 ()          ## visualize Y-projection
    >>> fz  = model.draw3 ()          ## visualize Z-projection

    Parameters
    ----------
    signal_x :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in X-direction
    signal_y :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in Y-direction
    signal_z :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in Z-direction
    suffix   : string
        An optional suffix to be  added to the names of created PDFs and variables
    bkg_1x   : RooFit/PDF, Ostap/PDF, integer, RooRealVar or None
        1D x-background for BSS-terms
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
    bkg_1y   : RooFit/PDF, Ostap/PDF, integer, RooRealVar or None
        1D y-background for BSS-terms
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
    bkg_2x   : RooFit/PDF, Ostap/PDF, integer, RooRealVar or None
        1D x-background for BBS-terms
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
    bkg_2y   : RooFit/PDF, Ostap/PDF, integer, RooRealVar or None
        1D y-background for BBS-terms
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
    bkg_2xy  : RooFit/PDF, Ostap/PDF, list/tuple or None
        2D (x,y)-background for BBS-terms
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created 
    bkg_2yz  : RooFit/PDF, Ostap/PDF, list/tuple or None
        2D (y,z)-background for BBS-termsm must be symmetric! 
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created
    bkg_3x   : RooFit/PDF, Ostap/PDF, integer, RooRealVar or None
        1D x-background for BBB-terms
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
    bkg_3y   : RooFit/PDF, Ostap/PDF, integer, RooRealVar or None
        1D y-background for BBB-terms
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created via self.make_bkg
    bkg_3yz  : RooFit/PDF, Ostap/PDF, list/tuple or None
        2D (y,z)-background for BBB-terms (must be symmetric!)
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created 
    bkg_3D   : RooFit/PDF, Ostap/PDF, list/tuple or None
        3D (x,y,z)-background for BBB-terms (must be symmetric for  y<-->z !)
        Use directly RooFit/PDF or Ostap/PDF, otherwise PDF is created 
        
    sss      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of SSS component
         Use directly RooRelaVar, otherwise create it using self.make_var function
    ssb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of SSB component
         Use directly RooRelaVar, otherwise create it using self.make_var function
    sbb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of SBB component
         Use directly RooRelaVar, otherwise create it using self.make_var function
    bbb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of BBB component
         Use directly RooRelaVar, otherwise create it using self.make_var function
         
    components : the list of additional 3D-PDFs to be used in the fit
    xvar       : x-variable
    yvar       : y-variable
    zvar       : z-variable
    name       : the PDF name 
         
    """
    def __init__ ( self               ,
                   #
                   signal_x           , 
                   signal_y           ,
                   signal_z = None    , ## cloned from signal_y if None 
                   suffix   = ''      ,
                   # background  for SSB-terms:
                   bkg_1x  = None  , ## 1D     x-background for BSS-terms
                   bkg_1y  = None  , ## 1D     y-background for BSS-terms,   z-background is cloned  
                   # background for  SBB-terms:
                   bkg_2x  = None  , ## 1D     x-background for SBB-terms
                   bkg_2y  = None  , ## 1D     y-background for SBB-terms,   z-background is cloned
                   # background for  SBB-terms:
                   bkg_2xy = None  , ## 2D   x,y-background for SBB-terms, (x,z)-backgroud is cloned 
                   bkg_2yz = None  , ## 3D   y,z-background for SBB-terms, must be SYMMETRIC!
                   # background for  BBB-terms:
                   bkg_3x  = None  , ## 1D     x-background for SBB-terms
                   bkg_3y  = None  , ## 1D     y-background for SBB-terms,   z-background is cloned
                   # background for  BBB-terms:
                   bkg_3yz = None  , ## 3D   y,z-background for SBB-terms, must be SYMMETRIC!
                   bkg_3D  = None  , ## 3D x,y,z-background for B(x,y,z) term (symmetric) 
                   ## Yields of the main components :
                   sss     = None  , ## S(x)*S(y)*S(z)             component 
                   bss     = None  , ## B(x)*S(y)*S(z)             component 
                   sbb     = None  , ## S(x)*B(y,z)                component 
                   ssb     = None  , ## S(x)*[S(y)*B(z)+B(y)*S(z)] component
                   bbs     = None  , ## B(x)*[S(y)*B(z)+B(y)*S(z)] component 
                   bbb     = None  , ## B(x,y,z)                   component    
                   ## additional components 
                   components = [] , ## the list additional components 
                   xvar    = None  , ## x-variable 
                   yvar    = None  , ## y-variable 
                   zvar    = None  , ## z-variable                   
                   name    = ''    ) : 
        
        ## keep all arguments 
        self.__args = {
            #
            'signal_x'   : signal_x   ,
            'signal_y'   : signal_y   ,
            'signal_z'   : signal_z   ,
            ##
            'bkg_1x'     : bkg_1x     ,
            'bkg_1y'     : bkg_1y     ,
            ##
            'bkg_2x'     : bkg_2x     ,
            'bkg_2y'     : bkg_2y     ,
            'bkg_2xy'    : bkg_2xy    ,
            'bkg_2yz'    : bkg_2yz    ,
            ##
            'bkg_3x'     : bkg_3x     ,
            'bkg_3y'     : bkg_3y     ,
            'bkg_3yz'    : bkg_3yz    ,
            ##
            'bkg_3D'     : bkg_3D     ,
            #
            'sss'        : sss        ,            
            'bss'        : bss        ,
            'ssb'        : ssb        ,
            'sbb'        : sbb        ,
            'bbs'        : bbs        ,
            'bbb'        : bbb        ,
            #
            'components' : components ,
            'xvar'       : xvar       ,
            'yvar'       : yvar       ,
            'zvar'       : zvar       ,
            'name'       : name                   
            }
        
        self.__suffix      = suffix


        if   isinstance ( signal_x , PDF            )          : self.__signal_x = signal_x
        elif isinstance ( signal_x , ROOT.RooAbsPdf ) and xvar :
            self.__signal_x = Generic1D_pdf ( signal_x , xvar , 'SX' )
        else : raise AttributeError ( "Invalid ``signal_x'' argumnent: %s" % signal_x )
            
        if   isinstance ( signal_y , PDF            )          : self.__signal_y = signal_y
        elif isinstance ( signal_y , ROOT.RooAbsPdf ) and yvar :
            self.__signal_y = Generic1D_pdf ( signal_y , yvar , 'SY' )
        else : raise AttributeError ( "Invalid ``signal_y'' argument: %s" %   signal_y )

        if   isinstance ( signal_z , PDF            )          : self.__signal_z = signal_z
        elif isinstance ( signal_z , ROOT.RooAbsPdf ) and zvar :
            self.__signal_z = Generic1D_pdf ( signal_z , zvar , 'SZ' )
        elif zvar and not signal_z :
            self.__signal_z = self.__signal_y.clone ( xvar = zvar , name = 'SZ' )
            self.debug('signal z-component is cloned from the signal_y component')
        else : raise AttributeError ( "Invalid ``signal_z'' argument: %s" % signal_z )
            
        #
        ## initialize base class
        #
        if not name :
            name = "%s&%s&%s" % ( self.__signal_x.name ,
                                  self.__signal_y.name ,
                                  self.__signal_z.name )                                  
            if suffix : name += '_'+ suffix
            
        PDF3.__init__ ( self , name ,
                        self.__signal_x.xvar ,
                        self.__signal_y.xvar ,
                        self.__signal_z.xvar ) 

        # =====================================================================
        ## 1) All signals component 
        # =====================================================================
        self.__sss_cmp  = Model3D (
            'SSS_pdf' + suffix , self.__signal_x , self.__signal_y , self.__signal_z )
        
        # =====================================================================
        ## 2) background x signal x  signal 
        # =====================================================================        
        self.__bkg_1x  = self.make_bkg (  bkg_1x , 'Bkg1X_BSS' + suffix , self.xvar )
        self.__bss_cmp = Model3D ( "BSS_pdf" + suffix ,
                                    self.__bkg_1x , self.__signal_y , self.__signal_z ,
                                    title = "Background1(x) x Signal(y) x Signal(z)" )
        
        # =====================================================================
        ## 3) signal x (  signal x background + backround x signal ) 
        # =====================================================================

        self.__bkg_1y  = self.make_bkg (        bkg_1y , 'Bkg1Y_SBS' + suffix , self.yvar )
        self.__bkg_1z  = self.make_bkg ( self.__bkg_1y , 'Bkg1Z_SSB' + suffix , self.zvar )
        
        self.__ssb_cmp_raw = Model3D ( "SSB_raw" + suffix ,
                                       self.__signal_x , self.__signal_y , self.__bkg_1z   ,
                                       title = "Signal(x) x Signal(y) x Background1(z)" )
        self.__sbs_cmp_raw = Model3D ( "SBS_raw" + suffix ,
                                       self.__signal_x , self.__bkg_1y   , self.__signal_z , 
                                       title = "Signal(x) x Background1(y) x Signal(z)" ) 

        self.__ssb_sym_cmp  = Generic3D_pdf (
            self.make_sum ( "SSB_pdf" + suffix ,
                            "S(x)*S(y)*B(z)+S(x)*B(y)*S(z)" ,
                            self.__ssb_cmp_raw.pdf   ,
                            self.__sbs_cmp_raw.pdf ) ,
            self.xvar , self.yvar , self.zvar )
        
        self.__ssb_cmp = self.__ssb_sym_cmp  ## 
        self.__sbs_cmp = self.__ssb_sym_cmp  ## ditto
        

        # =====================================================================
        ## (intermezzo-1) Assumptions about the SBB-background sub-components 
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
        elif component_similar ( bkg_2y ) :
            bkg_2y =        bkg_1y
            self.debug ( 'bkg_2y set to [SIMILAR] %s' % bkg_2y )
        # =====================================================================


        # =====================================================================
        ## 4) signal x background x background 
        # =====================================================================
        
        self.__bkg_2x = None 
        self.__bkg_2y = None 
        self.__bkg_2z = None 

        if   bkg_2yz and isinstance ( bkg_2yz , PDF2 ) :            
            self.__bkg_2yz = bkg_2yz            
        elif bkg_2yz and isinstance ( bkg_2yz , ROOT.RooAbsPdf ) :            
            self.__bkg_2yz = Generic2D_pdf ( bkg_2yz  , self.yvar , self.zvar , name = bkg2D.name )
        elif bkg_2yz and isinstance ( bkh_2yz , (  tuple , list ) ) :
            from ostap.fitting.models_2d import make_B2Dsym
            self.__bkg_2yz = make_B2Dsym ( 'Bkg2YZ' + suffix , self.yvar , self.zvar , *bkg_2yz  )            
        else :
            self.__bkg_2y = self.make_bkg (        bkg_2y , 'Bkg2Y_S2B' + suffix , self.yvar )        
            self.__bkg_2z = self.make_bkg ( self.__bkg_2y , 'Bkg2Z_S2B' + suffix , self.zvar )            
            self.__bkg_2yz = Model2D ( 'Bkg2YZ_pdf' + suffix ,
                                       self.__bkg_2y     ,
                                       self.__bkg_2z     , title =  'Background2(y,z)' )
            
        self.__sbb_cmp = Generic3D_pdf ( 
            ROOT.RooProdPdf ( "SBB_raw_pdf" + suffix , "Signal(x) x Background2(y,z)" , self.__signal_x.pdf , self.__bkg_2yz.pdf ) ,
            self.xvar , self.yvar , self.zvar )

        # =====================================================================
        ## 5) background x ( signal x background + background x  signal ) 
        # =====================================================================
        
        if   bkg_2xy and isinstance ( bkg_2xy , PDF2 ) :            
            self.__bkg_2xy = bkg_2xy
            self.__bkg_2xz = bkg_2xy.clone ( name = 'Bkg2XZ' + suffix ,
                                             xvar = self.xvar         ,
                                             yvar = self.zvar         )            
        elif bkg_2xy and isinstance ( bkg_2xy , ROOT.RooAbsPdf ) :            
            self.__bkg_2xy = Generic2D_pdf ( bkg_2xy  , self.xvar , self.yvar , name = bkg2D.name )
            self.__bkg_2xz = self.__bkg_2xy.clone ( name = 'Bkg2XZ' + suffix ,
                                                    xvar = self.xvar         ,
                                                    yvar = self.zvar         )  
        elif bkg_2xy and isinstance ( bkg_2xy ,  ( tuple , list )  ) :
            from ostap.fitting.models_2d import make_B2D
            self.__bkg_2xy = make_B2D ( 'Bkg2XY' + suffix , self.xvar , self.yvar , *bkg_2xy  )
            self.__bkg_2xz = self.__bkg_2xy.clone ( name = 'Bkg2XZ' + suffix ,
                                                    xvar = self.xvar         ,
                                                    yvar = self.zvar         )          
        else :
            
            if not self.__bkg_2x : self.__bkg_2x = self.make_bkg (        bkg_2x , 'Bkg2X_S2B' + suffix , self.xvar )        
            if not self.__bkg_2y : self.__bkg_2y = self.make_bkg (        bkg_2y , 'Bkg2Y_S2B' + suffix , self.yvar )        
            if not self.__bkg_2z : self.__bkg_2z = self.make_bkg ( self.__bkg_2y , 'Bkg2Z_S2B' + suffix , self.zvar )
            
            self.__bkg_2xy = Model2D ( 'Bkg2XY' + suffix , 
                                       self.__bkg_2x     ,
                                       self.__bkg_2y     , title = 'Background2(x,y)' )
            self.__bkg_2xz = Model2D ( 'Bkg2XZ' + suffix ,
                                       self.__bkg_2x     ,
                                       self.__bkg_2z     , title = 'Background2(x,z)' )
            
        self.__bsb_cmp_raw = Generic3D_pdf (
            ROOT.RooProdPdf ( "BSB_raw" + suffix , "Signal(y) x Background2(x,z)" , self.__signal_y.pdf , self.__bkg_2xz.pdf ) ,
            self.xvar , self.yvar , self.zvar )
        self.__bbs_cmp_raw = Generic3D_pdf ( 
            ROOT.RooProdPdf ( "BBS_raw" + suffix , "Signal(z) x Background2(x,y)" , self.__signal_z.pdf , self.__bkg_2xy.pdf ) ,
            self.xvar , self.yvar , self.zvar )

        self.__bbs_sym_cmp     = Generic3D_pdf (
            self.make_sum ( "SBB_pdf" + suffix ,
                            "S(x)*B(y,z)+S(y)*B(x,z)+S(Z)*B(x,y)" ,
                            self.__bsb_cmp_raw.pdf ,
                            self.__bbs_cmp_raw.pdf ) ,
            self.xvar , self.yvar , self.zvar )

        self.__bbs_cmp = self.__bbs_sym_cmp ## ditto
        self.__bsb_cmp = self.__bbs_sym_cmp
        
        # =====================================================================
        ## (intermezzo-2) Assumptions about the BBB-background sub-components 
        # =====================================================================        
        if   component_clone   ( bkg_3x ) :
            bkg_3x = self.__bkg_2x
            self.debug ( 'bkg_3x set to [CLONE]   %s' % bkg_3x ) 
        elif component_similar ( bkg_3x ) :
            bkg_3x =        bkg_2x
            self.debug ( 'bkg_3x set to [SIMILAR] %s' % bkg_3x ) 

        if   component_clone   ( bkg_3y ) :
            bkg_3y = self.__bkg_2y
            self.debug ( 'bkg_3y set to [CLONE]   %s' % bkg_3y ) 
        elif component_similar ( bkg_3y ) :
            bkg_3y =        bkg_2y
            self.debug ( 'bkg_3y set to [SIMILAR] %s' % bkg_3y ) 

        if   component_clone   ( bkg_3yz ) :
            bkg_3yz = self.__bkg_2yz
            self.debug ( 'bkg_3yz set to [CLONE]   %s' % bkg_3yz ) 
        elif component_similar ( bkg_3yz ) :
            bkg_3yz =        bkg_2yz
            self.debug ( 'bkg_3yz set to [SIMILAR] %s' % bkg_3yz ) 
        # =====================================================================

        # =====================================================================
        ## 6) pure background 
        # =====================================================================
        
        self.__bkg_3x  = None 
        self.__bkg_3y  = None 
        self.__bkg_3z  = None
        
        self.__bkg_3yz = None 

        if   bkg_3D and isinstance ( bkg_3D , PDF3 ) :
            self.__bbb_cmp = bkg_3D
        elif bkg_3D and isinstance ( bkg_3D , ROOT.RooAbsPdf ) : 
            self.__bbb_cmp = Generic3D_pdf ( bkg_3D , self.xvar , self.yvar , self.zvar )
        elif bkg_3D and isinstance ( bkg_2D ,  ( tuple , list )  ) :
            from ostap.fitting.models_2d import make_B2DmixYZ 
            self.__bkg_3D = make_B2DmixYZ ( 'BBB_pdf' + suffix , self.xvar , self.yvar , self.zvar , *bkg_3D )
        else :

            if   bkg_3yz and isinstance ( bkg_3yz , PDF2 ) :
                self.__bkg_3xy = bkg_3xy
            elif bkg_3yz and isinstance ( bkg_3yz , ROOT.RooAbsPdf ) :
                self.__bkg_3yz = Generic2D_pdf ( bkg_3yz , self.yvar , self.zvar )
            else :
                self.__bkg_3y  = self.make_bkg (        bkg_3y , 'Bkg3Y_BBB' + suffix , self.yvar )        
                self.__bkg_3z  = self.make_bkg ( self.__bkg_3y , 'Bkg3Z_BBB' + suffix , self.zvar )
                self.__bkg_3yz = Model2D ( 'Bkg3YZ' + suffix , self.__bkg_3y , self.__bkg_3z )
                
            self.__bkg_3x  = self.make_bkg ( bkg_3x , 'Bkg3X_BBB' + suffix , self.xvar )
            self.__bbb_cmp = Generic3D_pdf (
                ROOT.RooProdPdf ( "BBB_pdf" + suffix ,
                                  "Background3(x) x Background3(y,z)" ,
                                  self.__bkg_3x.pdf , self.__bkg_3yz.pdf ) ,
                self.xvar , self.yvar , self.zvar ) 
            
        # =====================================================================
        ## coefficients
        # =====================================================================
        
        self.__sss = self.make_var ( sss   , "SSS" + suffix ,
                                     "Signal(x)&Signal(y)&Signal(z)"     + suffix , sss , 1000 , 0 , 1.e+7 )
        self.__bss = self.make_var ( bss   , "BSS" + suffix ,
                                     "Background(x)&Signal(y)&Signal(z)" + suffix , bss , 1000 , 0 , 1.e+7 )
        self.__ssb = self.make_var ( ssb   , "SSB" + suffix ,
                                     "Signal&(SB+BS)"                    + suffix , sbb , 1000 , 0 , 1.e+7 )
        self.__sbb = self.make_var ( sbb   , "SBB" + suffix ,
                                     "Signal&Background(x,y)"            + suffix , sbb , 1000 , 0 , 1.e+7 )
        self.__bbs = self.make_var ( bbs   , "BBS" + suffix ,
                                     "Backgroun&(SB+BS)"                 + suffix , bbs , 1000 , 0 , 1.e+7 )        
        self.__bbb = self.make_var ( bbb  , "BBB"  + suffix ,
                                     "Background(x,y,z)"                 + suffix , bbb , 1000 , 0 , 1.e+7 )

        self.__sbs = self.__ssb ## the same 
        self.__bsb = self.__bbs ## the same 
        
        self.alist1 = ROOT.RooArgList (
            self.__sss_cmp.pdf ,
            self.__bss_cmp.pdf ,
            self.__ssb_cmp.pdf ,
            self.__sbb_cmp.pdf ,
            self.__bbs_cmp.pdf ,
            self.__bbb_cmp.pdf )
        self.alist2 = ROOT.RooArgList (
            self.__sss     ,
            self.__bss     ,
            self.__ssb     ,
            self.__sbb     ,
            self.__bbs     ,
            self.__bbb     )
        
        
        ## treat additional components (if specified)
        self.__nums_components = [] 
        icmp = 0
        self.__more_components = []
        for cmp in components :
            
            if   isinstance  ( c , PDF2           ) : cc = c  
            elif isinstance  ( c , ROOT.RooAbsPdf ) : cc = Generic2D_pdf ( cs ,  self.xvar , self.yvar ) 
            else :
                self.error ("unknown ``other''component %s/%s, skip it!" % ( cc , type(cc) ) )
                continue  
            self.__more_components.append ( cc     )
            self.components.add           ( cc.pdf ) 

        nc = len( self.__more_components )
        if 1 == nc :
            cf = self.make_var ( None , "C"+suffix , "Component" + suffix , None , 1 , 0 , 1.e+7 )
            self.alist1.add  ( self.components[0] )
            self.__num_components.append ( cf ) 
        elif 2 <= nc : 
            fic = self.make_fracs ( nc , 'C_%%d%s' % suffix ,  'C(%%d)%s'  % suffix , fractions  = False )
            for c in self.components : self.alist1.add ( c)
            for f in fic             : self.__num_components.append ( f )
            
        self.__nums_components  = tuple ( self.__nums_components  ) 
        for c in self.__nums_components  : self.alist2.add ( c )

        ## 
        #
        ## build the final PDF 
        # 
        self.pdf  = ROOT.RooAddPdf  ( "model3D"      + suffix ,
                                      "Model3D(%s)"  % suffix ,
                                      self.alist1 ,
                                      self.alist2 )

        self.signals     .add ( self.__sss_cmp.pdf )
        self.backgrounds .add ( self.__bbb_cmp.pdf )
        self.crossterms1 .add ( self.__ssb_cmp.pdf ) ## cross-terms
        self.crossterms2 .add ( self.__sbb_cmp.pdf ) ## cross-terms 
        
        ## save the configuration
        self.config = {
            ##
            'signal_x'   : self.signal_x ,
            'signal_y'   : self.signal_y ,
            'signal_z'   : self.signal_z ,
            ##
            'suffix'     : self.suffix  ,
            ## SSB-terms 
            'bkg_1x'     : self.bkg_1x  ,
            'bkg_1y'     : self.bkg_1y  ,
            ## SBB-terms 
            'bkg_2x'     : self.bkg_2x  ,
            'bkg_2y'     : self.bkg_2y  ,
            'bkg_2xz'    : self.bkg_2xz ,
            'bkg_2yz'    : self.bkg_2yz ,
            ## BBB-term
            'bkg_3x'     : self.bkg_3x  ,
            'bkg_3y'     : self.bkg_3y  ,
            ## 'bkg_3xz'    : self.bkg_3xz ,
            'bkg_3yz'    : self.bkg_3yz ,
            'bkg_3D'     : self.bkg_3D  ,
            ## Yields 
            'sss'        : self.SSS     ,
            'bss'        : self.BSS     ,
            'ssb'        : self.SSB     ,
            'sbb'        : self.SBB     ,
            'bbs'        : self.BBS     ,
            'bbb'        : self.BBB     ,            
            #
            'components' : self.more_components ,            
            'xvar'       : self.xvar    ,
            'yvar'       : self.yvar    ,
            'zvar'       : self.zvar    ,
            'name'       : self.name    ,             
            }
        
    ## redefine the clone method, allowing only the name to be changed
    #  @attention redefinition of parameters and variables is disabled,
    #             since it can't be done in a safe way                  
    def clone ( self , name = '' , xvar = None , yvar = None , zvar = None  ) :
        """Redefine the clone method, allowing only the name to be changed
         - redefinition of parameters and variables is disabled,
         since it can't be done in a safe way          
        """
        if xvar and not xvar is self.xvar :
            raise AttributeError("Fit3DMix can not be cloned with different `xvar''")
        if yvar and not yvar is self.yvar :
            raise AttributeError("Fit3DMix can not be cloned with different `yvar''")
        if zvar and not zvar is self.zvar :
            raise AttributeError("Fit3DMix can not be cloned with different `zvar''")
        return PDF.clone ( self , name = name ) if name else PDF.clone( self )

    @property
    def SSS ( self ) :
        """The yield of Signal(x)*Signal(y)*Signal(z) component"""
        return self.__sss
    @SSS.setter 
    def SSS ( self , value ) :
        value = float ( value  )
        assert value in self.__sss, "Value %s is out of the allowed range %s " % ( value , self.__sss.minmax() )
        self.__sss.setVal ( value ) 

    @property
    def SSB ( self ) :
        """The yield of S(x)*[S(y)*(z)+B(y)*S(z)] symmetrized component (same as SBS)"""
        return self.__bss
    @SSB.setter 
    def SSB ( self , value ) :
        value = float ( value  )
        assert value in self.__bss, "Value %s is out of the allowed range %s " % ( value , self.__bss.minmax() )
        self.__bss.setVal ( value ) 

    @property
    def SBS ( self ) :
        """The yield of S(x)*[S(y)*B(z)+B(y)*S(z)] symmetrized component (same as SSB)"""
        return self.__sbs
    @SBS.setter 
    def SBS ( self , value ) :
        value = float ( value  )
        assert value in self.__sbs, "Value %s is out of the allowed range %s " % ( value , self.__sbs.minmax() )
        self.__sbs.setVal ( value ) 

    @property
    def BSS ( self ) :
        """The yield of B(x)*S(y)*S(z) component"""
        return self.__bss
    @BSS.setter 
    def BSS ( self , value ) :
        value = float ( value  )
        assert value in self.__bss, "Value %s is out of the allowed range %s " % ( value , self.__bss.minmax() )
        self.__bss.setVal ( value ) 

    @property
    def SBB ( self ) :
        """The yield of S(x)*B(y,z) component"""
        return self.__sbb
    @SBB.setter 
    def SBB ( self , value ) :
        value = float ( value  )
        assert value in self.__sbb, "Value %s is out of the allowed range %s " % ( value , self.__sbb.minmax() )
        self.__sbb.setVal ( value ) 

    @property
    def BSB ( self ) :
        """The yield of B(x)*[S(y)*B(z)+B(y)*S(z)] symmetrized component (same as BBS)"""
        return self.__bsb
    @BSB.setter 
    def BSB ( self , value ) :
        value = float ( value  )
        assert value in self.__bsb, "Value %s is out of the allowed range %s " % ( value , self.__bsb.minmax() )
        self.__bsb.setVal ( value ) 

    @property
    def BBS ( self ) :
        """The yield of B(x)*[S(y)*B(z)+B(y)*S(z)] symmetrized component (same as BSB)"""
        return self.__bbs
    @BBS.setter 
    def BBS ( self , value ) :
        value = float ( value  )
        assert value in self.__bbs, "Value %s is out of the allowed range %s " % ( value , self.__bbs.minmax() )
        self.__bbs.setVal ( value ) 

    @property
    def BBB ( self ) :
        """The yield of B(x,y,z) component"""
        return self.__bbb
    @BBB.setter 
    def BBB ( self , value ) :
        value = float ( value  )
        assert value in self.__bbb, "Value %s is out of the allowed range %s " % ( value , self.__bbb.minmax() )
        self.__bbb.setVal ( value ) 

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
    def signal_z ( self  ) :
        """``signal_z'': Signal(z) component/PDF"""
        return self.__signal_z

    @property
    def bkg_1x ( self ) :
        """``bkg_1x'': B(x) component for B(x)*S(y)*S(z) term"""
        return self.__bkg_1x
    @property
    def bkg_1y ( self ) :
        """``bkg_1y'': B(y) component for S(x)*B(y)*S(z) term"""
        return self.__bkg_1y
    @property
    def bkg_1z ( self ) :
        """``bkg_1z'': B(z) component for S(x)*S(y)*B(z) term"""
        return self.__bkg_1z
    
    @property
    def bkg_2xy( self ) :
        """``bkg_2xy'': B(x,y) component for B(x,y)*S(z) term"""
        return self.__bkg_2xy
    @property
    def bkg_2xz( self ) :
        """``bkg_2xz'': B(x,z) component for B(x,z)*S(y) term"""
        return self.__bkg_2xz
    @property
    def bkg_2yz( self ) :
        """``bkg_2yz'': B(y,z) component for B(y,z)*S(x) term"""
        return self.__bkg_2yz 

    @property
    def bkg_2x ( self ) :
        """``bkg_2x'': B(x) component for B(x,y)*S(z) & B(x,z)*S(y) terms"""
        return self.__bkg_2x
    @property
    def bkg_2y ( self ) :
        """``bkg_2y'': B(y) component for B(y,z)*S(x) & B(x,y)*S(z) terms"""
        return self.__bkg_2y
    @property
    def bkg_2z ( self ) :
        """``bkg_2z'': B(z) component for B(x,z)*S(y) & B(y,z)*S(x) terms"""
        return self.__bkg_2z

    @property
    def bkg_3x ( self ) :
        """``bkg_3x'': B(x) component for B(x)*B(y)*B(z) term"""
        return self.__bkg_3x
    @property
    def bkg_3y ( self ) :
        """``bkg_3y'': B(y) component for B(x)*B(y)*B(z) term"""
        return self.__bkg_3y
    @property
    def bkg_3z ( self ) :
        """``bkg_3z'': B(z) component for B(x)*B(y)*B(z) term"""
        return self.__bkg_3z

    @property
    def bkg_3yz ( self ) :
        """``bkg_3yz'': B(y,z) component for B(x)*B(y,z) term"""
        return self.__bkg_3yz

    @property
    def bkg_3D ( self ) :
        """B(x,y,z) component/PDF for the final PDF"""
        return self.__bbb_cmp 

    @property
    def cmp_SSS ( self ) :
        """```triple-signal'' component/PDF"""
        return self.__sss_cmp

    @property
    def cpm_SSB ( self ) :
        """```signal-signal-background'' symmetrized component/PDF"""
        return self.__ssb_cmp
    
    @property
    def cmp_SBS ( self ) :
        """```signal-background-signal'' symmetrized component/PDF (same as cmp_SSB)"""
        return self.__sbs_cmp
    
    @property
    def cmp_BSS ( self ) :
        """```background-signal-signal'' symmetrized component/PDF (same as cmp_SSB)"""
        return self.__bss_cmp

    @property
    def cpm_SBB ( self ) :
        """```signal-background-background'' symmetrized component/PDF"""
        return self.__sbb_cmp
    
    @property
    def cmp_BSB ( self ) :
        """```background-signal-background'' symmetrized component/PDF (same as cmp_SBB)"""
        return self.__bsb_cmp
    
    @property
    def cmp_BBS ( self ) :
        """```background-background-signal'' symmetrized component/PDF (same as cmp_SBB)"""
        return self.__bbs_cmp

    @property
    def cmp_BBB ( self ) :
        """```triple-background'' component/PDF"""
        return self.__bbb_cmp

    @property
    def more_components ( self ) :
        """additional/``other'' components"""
        return tuple( self.__more_components  )

    @property
    def suffix ( self ) :
        """``suffix'', used to build the name"""
        return self.__suffix

    # =========================================================================
    ## Raw, non-symmetrized fit components/PDF (for debugging)
    # =========================================================================

    @property
    def cpm_raw_SSB ( self ) :
        """```signal-signal-background'' raw,non-symmetrized component/PDF"""
        return self.__ssb_cmp_raw
    
    @property
    def cmp_raw_SBS ( self ) :
        """```signal-background-signal'' raw,non-symmetrized component/PDF"""
        return self.__sbs_cmp_raw
    
    @property
    def cmp_raw_BBS ( self ) :
        """```background-background-signal'' raw,non-symmetrized component/PDF"""
        return self.__bss_cmp_raw

    @property
    def cmp_raw_BSB ( self ) :
        """```background-signal-background'' raw,non-symmetrized component/PDF"""
        return self.__bsb_cmp_raw


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
