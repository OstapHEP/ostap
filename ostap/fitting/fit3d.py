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
    'PDF3'          , ## useful base class for 2D-models
    'Fit3D'         , ## the model for 2D-fit: signal + background + optional components
    ##
    'H3D_dset'      , ## convertor of 2D-histo to RooDataHist 
    'H3D_pdf'       , ## convertor of 1D-histo to RooDataPdf
    ##
    'Generic3D_pdf' , ## wrapper over imported RooFit (2D)-pdf  
    'Flat3D'        , ## the most trivial 3D-pdf - constant
    'Model3D'       , ## trivial class to build 3D model from 1D-components 
    )
# =============================================================================
import ROOT
from   ostap.fitting.basic import makeVar, makeBkg, H3D_dset
from   ostap.fitting.fit2d import PDF2
from   ostap.logger.utils  import roo_silent, rooSilent 
# =============================================================================
from   ostap.logger.logger import getLogger
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
            self.__zvar = makeVar ( zvar               , ## var 
                                    'z'                , ## name 
                                    'z-variable(mass)' , ## title/comment
                                    None               , ## fix ?
                                    *zvar              ) ## min/max 
        elif isinstance ( zvar , ROOT.RooAbsReal ) :
            self.__zvar = makeVar ( zvar               , ## var 
                                    'z'                , ## name 
                                    'z-variable/mass'  , ## title/comment
                                    fix = None         ) ## fix ? 
        else :
            logger.warning('PDF3: ``z-variable''is not specified properly %s/%s' % ( zvar , type ( zvar ) ) )
            self.__zvar = makeVar( zvar , 'z' , 'z-variable' )

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
                refit  = False , *args , **kwargs ) :
        """
        Perform the actual fit (and draw it)
        >>> r,f = model.fitTo ( dataset )
        >>> r,f = model.fitTo ( dataset , weighted = True )    
        >>> r,f = model.fitTo ( dataset , ncpu     = 10   )    
        >>> r,f = model.fitTo ( dataset )    
        """
        if isinstance ( dataset , ROOT.TH3 ) :
            density = kwargs.pop ( 'density' , True  ) 
            chi2    = kwargs.pop ( 'chi2'    , False ) 
            return self.fitHisto ( dataset   , draw , silent , density , chi2 , *args , **kwargs )

        
        result,f = PDF2.fitTo ( self    ,
                                dataset ,
                                False   , ## false here!
                                20      , ## fake here..
                                silent  ,
                                refit   , *args , **kwargs )
        
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
                in_range3 = None ,
                *args , **kwargs ) :
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
        return self.draw ( self.xvar , 
                           dataset  ,
                           nbins    ,
                           20       , ## fake 
                           silent   ,
                           in_range = in_range , *args , **kwargs )


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
                in_range3 = None ,
                *args , **kwargs ) :
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
        return self.draw ( self.yvar , 
                           dataset   ,
                           nbins     ,
                           20        , ## fake 
                           silent    ,
                           in_range = in_range , *args , **kwargs )


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
    def draw3 ( self            ,
                dataset   = None ,
                nbins     = 100  ,
                silent    = True ,
                in_range1 = None ,
                in_range2 = None ,
                *args , **kwargs ) :
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
        return self.draw ( self.zvar , 
                           dataset   ,
                           nbins     ,
                           20        , ## fake 
                           silent    ,
                           in_range = in_range , *args , **kwargs )

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
                   chi2    = False , *args , **kwargs ) :
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
                                              density = density , *args , **kwargs )
            else    : return self.fitTo     ( data              ,
                                              silent = silent   , *args , **kwargs ) 

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
        self.__adjustment = Adjust3D ( self.name ,
                                       self.xvar , self.yvar , self.zvar ,
                                       self.pdf  , value )
        ## replace the original PDF  with  adjusted one:
        self.pdf          = self.__adjustment.pdf

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
        from ostap.core.core import dsID
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
    def __call__ ( self , x , y , z ) :
        
        if     isinstance ( self.xvar , ROOT.RooRealVar ) and \
               isinstance ( self.yvar , ROOT.RooRealVar ) and \
               isinstance ( self.zvar , ROOT.RooRealVar ) :
           
            from ostap.fitting.roofit import SETVAR
            from ostap.math.ve        import VE
            if x in self.xvar and y in self.yvar and z in self.zvar : 
                with SETVAR( self.xvar ) , SETVAR( self.yvar ) ,  SETVAR( self.zvar ) :
                    self.xvar.setVal ( x )
                    self.yvar.setVal ( y )
                    self.zvar.setVal ( z )
                    if error and self.fit_result :
                        e = self.eff_fun.getPropagatedError ( self.fit_result )
                        if 0<= e : return  VE ( v ,  e * e )
                    return v 
            else : return 0.0
            
        raise AttributeError, 'something wrong goes here'


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
        if hasattr ( self.pdf , 'function' ) :
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
        from ostap.fitting.roofit import SETVAR
        import random
        for i in xrange ( nshoots ) : 
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
    def integral ( self, xmin , xmax , ymin , ymax , zmin , zmax ) :
        """Get integral over (xmin,xmax,ymin,ymax,zmin,zmax) region
        >>> pdf = ...
        >>> print pdf.integral( 0,1,0,2,0,5)
        """
        xmn , xmx = self.xminmax()
        ymn , ymx = self.yminmax()
        zmn , zmx = self.zminmax()

        xmin = max ( xmin , xmn )
        xmax = min ( xmax , xmx )
        ymin = max ( ymin , ymn )
        ymax = min ( ymax , ymx )
        zmin = max ( zmin , zmn )
        zmax = min ( zmax , zmx )
        
        ## make a try to use analytical integral (could be fast)
        if hasattr ( self , 'pdf' ) :
            _pdf = self.pdf 
            if hasattr ( _pdf , 'setPars'  ) : _pdf.setPars() 
            try: 
                if hasattr ( _pdf , 'function' ) :
                    _func = _pdf.function() 
                    if hasattr ( _func , 'integral' ) :
                        return _func.integral ( xmin , xmax , ymin , ymax , zmin , zmax )
            except:
                pass
            
        ## use numerical integration 
        from ostap.math.integral import integral3 as _integral3 
        return _integral3 ( self ,
                            xmin , xmax ,
                            ymin , ymax ,
                            zmin , zmax )

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
        
        if not name : name = pdf.GetName()
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
    def __init__ ( self , xvar , yvar , zvar , name = 'Flat3D') :
        
        PDF3.__init__ ( self  , name , xvar , yvar , zvar ) 
        
        self.__xp0  = ROOT.RooPolynomial( 'xp0_%s'   % name , 'xpoly0(%s)'   % name , xvar )        
        self.__yp0  = ROOT.RooPolynomial( 'yp0_%s'   % name , 'ypoly0(%s)'   % name , yvar )
        self.__zp0  = ROOT.RooPolynomial( 'zp0_%s'   % name , 'zpoly0(%s)'   % name , zvar )
        
        self.__vlst = ROOT.RooArgList ( self.__xp0 ,  self.__yp0 ,  self.__zp0 )
        
        self.pdf    = ROOT.RooProdPdf ( name , 'poly0_3D(%s)' % name , self.__vlst )
        ## save configuration
        self.config = {
            'name'     : self.name ,            
            'xvar'     : self.xvar ,
            'yvar'     : self.yvar ,
            'zvar'     : self.zvar ,
            }                   

# =============================================================================
## simple class to adjust certaint PDF to avoid zeroes 
class Adjust3D(object) :
    """Simple class to ``adjust'' certain PDF to avoid zeroes
    - a small flat component is added and the compound PDF is constructed
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,
                   xvar             , 
                   yvar             , 
                   zvar             , 
                   pdf              ,
                   value    = 1.e-5 ) : 
        
        assert isinstance ( pdf  , ROOT.RooAbsPdf  ) , "``pdf''  must be ROOT.RooAbsPdf"
        assert isinstance ( xvar , ROOT.RooAbsReal ) , "``xvar'' must be ROOT.RooAbsReal"
        assert isinstance ( yvar , ROOT.RooAbsReal ) , "``yvar'' must be ROOT.RooAbsReal"
        assert isinstance ( zvar , ROOT.RooAbsReal ) , "``zvar'' must be ROOT.RooAbsReal"
        
        self.name      = name 
        self.__old_pdf = pdf

        ##
        self.__flat    = Flat3D  ( xvar , yvar ,  xvar , name = 'flat_' + name )
        self.__frac    = makeVar ( value , 'fracA_%s'                     % name ,
                                   'small  fraction of flat component %s' % name ,
                                   value , 1.e-4 , 0 , 1 )
        
        self.__alist1 = ROOT.RooArgList ( self.__flat.pdf , self.__old_pdf )        
        self.__alist2 = ROOT.RooArgList ( self.__frac     )
        #
        ## final PDF
        # 
        self.__pdf    = ROOT.RooAddPdf  ( "adjust_"    + name ,
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
        
        ## build pdf 
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
## @class Fit3D
#  The actual model for 3D-fits
#
#  @code
# 
#  model   = Models.Fit3D (
#      signal_1 = Models.Gauss_pdf ( 'Gx' , mass = m_x ) ,
#      signal_2 = Models.Gauss_pdf ( 'Gy' , mass = m_y ) ,
#      signal_3 = Models.Gauss_pdf ( 'Gz' , mass = m_z ) ,
#      power1   = 1 , 
#      power2   = 1 ,
#      power3   = 1 )
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
    ...      signal_1 = Models.Gauss_pdf ( 'Gx' , mass = m_x ) ,
    ...      signal_2 = Models.Gauss_pdf ( 'Gy' , mass = m_y ) ,
    ...      signal_3 = Models.Gauss_pdf ( 'Gz' , mass = m_z ) ,
    ...      bkgX1    = 1 , 
    ...      bkgY1    = 0 ,
    ...      bkgZ1    = 0 )
    >>> r,f = model.fitTo ( dataset ) ## fit dataset 
    >>> print r                       ## get results  
    >>> fx  = model.draw1 ()          ## visualize X-projection
    >>> fy  = model.draw2 ()          ## visualize Y-projection
    >>> fz  = model.draw3 ()          ## visualize Z-projection

    Parameters
    ----------
    signal_1 :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in X-direction
    signal_2 :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in Y-direction
    signal_3 :  RooFit/PDF or Ostap/PDF 
        PDF to describe (1D)-signal in Z-direction
    suffix   : string
        An optional suffix to be  added to the names of created PDFs and variables
    bkgX1    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for Bx1(x)* Sy(y)* Sz(z) term
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkgY1    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for Sx(z)*By1(y)* Sz(z) term
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkgZ1    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for   Sx(x)* Sy(y)*Bz1(z) term 
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkgXY    : RooFit/PDF, Ostap/PDF or None
        2D-background for Bxy(x,y)*Sz(z) term  
        Use directly RooFit/PDF or Ostap/PDF otherwise create from bkgX2 and bkgY2
    bkgXZ    : RooFit/PDF, Ostap/PDF or None
        2D-background for Bxz(x,z)*Sy(y) term 
        Use directly RooFit/PDF or Ostap/PDF otherwise create from bkgX2 and bkgZ2
    bkgYZ    : RooFit/PDF, Ostap/PDF or None
        2D-background for Bxz(x,z)*Sy(y) term        
        Use directly RooFit/PDF or Ostap/PDF otherwise create from bkgX2 and bkgZ2
    bkgX2    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for used to create Bxy(x,y) and Bxz(x,z) if they are not specified.
        If None - bkgX1 is used 
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkgY2    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for used to create Bxy(x,y) and Byz(y,z) if they are not specified.
        If None - bkgY1 is used 
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkgZ2    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for used to create Bxz(x,z) and Byz(y,z) if they are not specified.
        If None - bkgZ1 is used
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkgX3    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for used to create Bxyz(x,y,z) if it is not specified.
        If None - bkgX2 is used 
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkgY3    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for used to create Bxyz(x,y,z) if it is not specified.
        If None - bkgY2 is used 
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
    bkgZ3    : RooFit/PDF, Ostap/PDF, non-negative integer, RooRealVar or None
        1D-background for used to create Bxyz(x,y,z) if it is not specified.
        If None - bkgZ2 is used         
        Use directly RooFit/PDF or Ostap/PDF, otherwise create and use Bkg_pdf
        
    sss      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of sig(1) * sig(2) * sig(3) component
         Use directly RooRelaVar, otherwise create it using makeVar function
    ssb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of sig(1) * sig(2) * bkg (3) component
         Use directly RooRelaVar, otherwise create it using makeVar function
    sbs      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of sig(1) * bkg(2) * sig (3) component
         Use directly RooRelaVar, otherwise create it using makeVar function
    bss      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of bkg(1) * sig(2) * sig (3) component
         Use directly RooRelaVar, otherwise create it using makeVar function
    sbb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of sig(1) * bkg(2) * bkg (3) component
         Use directly RooRelaVar, otherwise create it using makeVar function
    bsb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of bkg(1) * sig(2) * bkg (3) component
         Use directly RooRelaVar, otherwise create it using makeVar function
    bbs      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of bkg(1) * bkg(2) * sig (3) component
         Use directly RooRelaVar, otherwise create it using makeVar function
    bbb      : None, RooRealVar, non-negative float or tuple
         Variable for the yield of bkg(1) * bkg(2) * bkg (3) component
         Use directly RooRelaVar, otherwise create it using makeVar function

    """
    def __init__ ( self               ,
                   #
                   signal_1           , 
                   signal_2           ,
                   signal_3           ,
                   suffix = ''        ,
                   #
                   bkgX1      = None  , ## 1D-background for Bx1(x)* Sy(y)* Sz(z) term
                   bkgY1      = None  , ## 1D-background for  Sx(z)*By1(y)* Sz(z) term
                   bkgZ1      = None  , ## 1D-background for  Sx(x)* Sy(y)*Bz1(z) term 
                   #
                   bkgXY      = None  , ## 2D-background for Bxy(x,y)*Sz(z) term 
                   bkgXZ      = None  , ## 2D-background for Bxz(x,z)*Sy(y) term 
                   bkgYZ      = None  , ## 2D-background for Byz(y,z)*Sx(z) term  
                   #
                   ## *if* no XY,XZ,BC backgrounds are specified, combine them from 
                   bkgX2      = None  , ## Bxy(x,y) = Bx2(x)*By2(y)
                   bkgY2      = None  , ## Bxz(x,z) = Bx2(x)*Bz2(z)
                   bkgZ2      = None  , ## Bkg(y,z) = By2(y)*Bz2(z)
                   ##
                   bkg3D      = None  , ## 3D-backround component  B(x,y,z)
                   #
                   ## *if* no 3D-background components is specified, combine it from
                   bkgX3      = None  , ## Bkg(x,y,z) = Bx3(x)*By3(y)*Bz3(z)
                   bkgY3      = None  , ## Bkg(x,y,z) = Bx3(x)*By3(y)*Bz3(z)
                   bkgZ3      = None  , ## Bkg(x,y,z) = Bx3(x)*By3(y)*Bz3(z)
                   #
                   ## Yields of the main components :
                   sss        = None  , ## sig(1) * sig(2) * sig(3) 
                   ssb        = None  , ## sig(1) * sig(2) * bkg(3)
                   sbs        = None  , ## sig(1) * bkg(2) * sig(3)
                   bss        = None  , ## bkg(1) * sig(2) * sig(3)
                   sbb        = None  , ## sig(1) * bkg(2) * bks(3)
                   bsb        = None  , ## bkg(1) * sig(2) * bks(3)
                   bbs        = None  , ## bkg(1) * bkg(2) * sig(3)
                   bbb        = None  , ## background-3D 
                   ## additional components 
                   components = []    ,
                   xvar       = None  ,
                   yvar       = None  ,
                   zvar       = None  ,                   
                   name       = ''    ) : 

        ## keep all arguments 
        self.__args = {
            'signal_1' : signal_1 , 'signal_2' : signal_2 , 'signal_3' : signal_3 ,
            'bkgX1'    : bkgX1    , 'bkgY1'    : bkgY1    , 'bkgZ1'    : bkgZ1    ,
            'bkgXY'    : bkgXY    , 'bkgXZ'    : bkgXZ    , 'bkgYZ'    : bkgYZ    ,
            'bkgX2'    : bkgX2    , 'bkgY2'    : bkgY2    , 'bkgZ2'    : bkgZ2    ,
            'bkg3D'    : bkg3D    ,
            'bkgX3'    : bkgX3    , 'bkgY3'    : bkgY3    , 'bkgZ3'    : bkgZ3    ,
            'sss'      : sss      ,            
            'ssb'      : ssb      ,
            'sbs'      : sbs      ,
            'bss'      : bss      ,
            'sbb'      : sbb      ,
            'bsb'      : bsb      ,
            'bbs'      : bbs      ,
            'bbb'      : bbb      ,
            'name'     : name     
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

        if   isinstance ( signal_3 , PDF            )          : self.__signal3 = signal_3
        elif isinstance ( signal_3 , ROOT.RooAbsPdf ) and zvar :
            self.__signal3 = Generic1D_pdf ( signal_3 , zvar , 'SIGNAL-Z' )
        else : raise AttributeError ( "Invalid ``signal3'' attribute" )
            
        #
        ## initialize base class
        #
        if not name and self.__signal1.name and self.__signal2.name and self.__signal3.name :
            name = '%s_and_%s_and_%s_%s' % ( self.__signal1.name ,
                                             self.__signal2.name ,
                                             self.__signal3.name , suffix )
            
        PDF3.__init__ ( self , name ,
                        self.__signal1.xvar ,
                        self.__signal2.xvar ,
                        self.__signal3.xvar ) 
     
        # =====================================================================
        ## 1) First component: all   signals
        # =====================================================================
        
        self.__sss_cmp = Model3D ( 'SSS_pdf' + suffix ,
                                   self.__signal1 , self.__signal2 , self.__signal3 ,
                                   title = "Sig(1) x Sig(2) x Sig(3)" )
                            
        # =====================================================================
        ## 2-4) Three terms:  ( 2 signals )  x ( 1 background ) 
        # =====================================================================
        
        self.__bkgX1   = makeBkg ( bkgX1 , 'Bkg_BSS' + suffix , self.xvar )
        self.__bkgY1   = makeBkg ( bkgY1 , 'Bkg_SBS' + suffix , self.yvar )
        self.__bkgZ1   = makeBkg ( bkgZ1 , 'Bkg_SSB' + suffix , self.zvar )
        
        self.__ssb_cmp = Model3D ( "SSB_pdf" + suffix ,
                                   self.__signal1 , self.__signal2 , self.__bkgZ1   ,
                                   title = "Sig(1) x Sig(2) x Bkg(3)" )
        self.__sbs_cmp = Model3D ( "SBS_pdf" + suffix ,
                                   self.__signal1 , self.__bkgY1   , self.__signal3 , 
                                   title = "Sig(1) x Bkg(2) x Sig(3)" ) 
        self.__bss_cmp = Model3D ( "BSS_pdf" + suffix ,
                                   self.__bkgX1   , self.__signal2 , self.__signal3 ,
                                   title = "Bkg(1) x Sig(2) x Sig(3)" )
        
        # =====================================================================
        ## 5-7) Three terms: (1 signal) x (2 backgrounds)
        # =====================================================================
        
        if bkgX2 is None : bkgX2 = self.__bkgX1
        if bkgY2 is None : bkgY2 = self.__bkgY1
        if bkgZ2 is None : bkgZ2 = self.__bkgZ1
        
        self.__bkgX2  = makeBkg ( bkgX2 , 'BkgX_S2B' + suffix , self.xvar )        
        self.__bkgY2  = makeBkg ( bkgY2 , 'BkgY_S2B' + suffix , self.yvar )        
        self.__bkgZ2  = makeBkg ( bkgZ2 , 'BkgZ_S2B' + suffix , self.zvar )

        if   bkgXY and isinstance ( bkgXY , PDF2 ) :
            self.__bkgXY = bkgXY
        elif bkgXY and isinstance ( bkgXY , ROOT.RooAbsPdf ) : 
            self.__bkgXY = Generic2D_pdf ( bkgXY  , self.xvar , self.yvar , name = bkgXY.name )
        else :
            self.__bkgXY = Model2D ( 'BkgXY_pdf' + suffix       ,
                                     self.__bkgX2               ,
                                     self.__bkgY2               ,
                                     title =  'Bkg(1) x Bkg(2)' )
        if   bkgXZ and isinstance ( bkgXZ , PDF2 ) :
            self.__bkgXZ = bkgXZ
        elif bkgXZ and isinstance ( bkgXZ , ROOT.RooAbsPdf ) : 
            self.__bkgXZ = Generic2D_pdf ( bkgXZ  , self.xvar , self.zvar , name = bkgXZ.name )
        else :
            self.__bkgXZ = Model2D ( 'BkgXZ_pdf' + suffix      ,                                     
                                     self.__bkgX2              ,
                                     self.__bkgZ2              , 
                                     title = 'Bkg(1) x Bkg(3)' )

        if   bkgYZ and isinstance ( bkgYZ , PDF2 ) :
            self.__bkgYZ = bkgYZ
        elif bkgYZ and isinstance ( bkgYZ , ROOT.RooAbsPdf ) : 
            self.__bkgYZ = Generic2D_pdf ( bkgYZ  , self.yvar , self.zvar , name = bkgYZ.name )
        else :            
            self.__bkgYZ = Model2D ( 'BkgYZ_pdf' + suffix      ,
                                     self.__bkgY2              ,
                                     self.__bkgZ2              ,
                                     title = 'Bkg(2) x Bkg(3)' )

        self.__sbb_cmp = Generic3D_pdf ( 
            ROOT.RooProdPdf ( "SBB_pdf" + suffix , "Sig(1) x Bkg(2,3)" , self.__signal1.pdf , self.__bkgYZ.pdf ) ,
            self.xvar , self.yvar , self.zvar )
        self.__bsb_cmp = Generic3D_pdf (
            ROOT.RooProdPdf ( "BSB_pdf" + suffix , "Sig(2) x Bkg(1,3)" , self.__signal2.pdf , self.__bkgXZ.pdf ) ,
            self.xvar , self.yvar , self.zvar )
        self.__bbs_cmp = Generic3D_pdf ( 
            ROOT.RooProdPdf ( "BBS_pdf" + suffix , "Sig(3) x Bkg(1,2)" , self.__signal3.pdf , self.__bkgXY.pdf ) ,
            self.xvar , self.yvar , self.zvar )
        
        # =====================================================================
        ## 8) pure background 
        # =====================================================================
        
        self.__bkgX3 = None 
        self.__bkgY3 = None 
        self.__bkgZ3 = None 
        
        if   bkg3D and isinstance ( bkg3D , PDF3 ) :
            self.__bbb_cmp = bkg3D
        elif bkg3D and isinstance ( bkgYZ , ROOT.RooAbsPdf ) : 
            self.__bbb_cmp = Generic3D_pdf ( bkg3D , self.xvar , self.yvar , self.zvar )
        else :

            if bkgX3 is None : bkgX3 = self.__bkgX2
            if bkgY3 is None : bkgY3 = self.__bkgY2
            if bkgZ3 is None : bkgZ3 = self.__bkgZ2
            
            self.__bkgX3 = makeBkg ( bkgX3 , 'BkgX_BBB' + suffix , self.xvar )
            self.__bkgY3 = makeBkg ( bkgY3 , 'BkgY_BBB' + suffix , self.yvar )
            self.__bkgZ3 = makeBkg ( bkgZ3 , 'BkgZ_BBB' + suffix , self.zvar )
            
            self.__bbb_cmp = Model3D (
                "BBB_pdf" + suffix ,
                self.__bkgX3 ,
                self.__bkgY3 ,
                self.__bkgZ3 ,
                title = "Bkg(1) x Bkg(2) x Bkg(3)" )
        #
        ## coefficients
        #
        self.__sss = makeVar ( sss   , "SSS"          + suffix ,
                               "Sig(1)&Sig(2)&Sig(3)" + suffix , sss , 1000  , 0 ,  1.e+8 )
        self.__ssb = makeVar ( ssb   , "SSB"          + suffix ,
                               "Sig(1)&Sig(2)&Bkg(3)" + suffix , ssb , 1000  , 0 ,  1.e+8 )
        self.__sbs = makeVar ( sbs   , "SBS"          + suffix ,
                               "Sig(1)&Bkg(2)&Sig(3)" + suffix , sbs , 1000  , 0 ,  1.e+8 )
        self.__bss = makeVar ( bss   , "BSS"          + suffix ,
                               "Bkg(1)&Sig(2)&Sig(3)" + suffix , bss , 1000  , 0 ,  1.e+8 )
        self.__sbb = makeVar ( sbb  , "SBB"           + suffix ,
                               "Sig(1)&Bkg(2)&Bkg(3)" + suffix , sbb , 1000  , 0 ,  1.e+8 )
        self.__bsb = makeVar ( bsb  , "BSB"           + suffix ,
                               "Bkg(1)&Sig(2)&Bkg(3)" + suffix , bsb , 1000  , 0 ,  1.e+8 )
        self.__bbs = makeVar ( bbs  ,  "BBS"          + suffix ,
                               "Bkg(1)&Bkg(2)&Sig(3)" + suffix , bbs , 1000  , 0 ,  1.e+8 )
        self.__bbb = makeVar ( bbb  , "BBB"           + suffix ,
                               "Bkg(1)&Bkg(2)&Bkg(3)" + suffix , bbb , 1000  , 0 ,  1.e+8 )
        
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
            'signal_1'   : self.signal1 ,
            'signal_2'   : self.signal2 ,
            'signal_3'   : self.signal3 ,
            'suffix'     : self.suffix  ,
            ##
            'bkgX1'      : self.bkgX1   ,
            'bkgY1'      : self.bkgY1   ,
            'bkgZ1'      : self.bkgZ1   ,
            ##
            'bkgXY'      : self.bkgXY   ,
            'bkgXZ'      : self.bkgXZ   ,
            'bkgYZ'      : self.bkgYZ   ,
            ##
            'bkgX2'      : self.bkgX2   ,
            'bkgY2'      : self.bkgY2   ,
            'bkgZ2'      : self.bkgZ2   ,
            ##
            'bkg3D'      : self.bkg3D   ,
            ##
            'bkgX3'      : self.bkgX3   ,
            'bkgY3'      : self.bkgY3   ,
            'bkgZ3'      : self.bkgZ3   ,
            #
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
        
    @property 
    def crossterms1 ( self ) :
        """``cross-terms'': pdfs for two signals andone background"""        
        return self.__crossterms1

    @property
    def crossterms2 ( self ) :
        """``cross-terms'':  pdfs for two backgrounds and one signal"""        
        return self.__crossterms2

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
        """The yield of Signal(x)*Background(y)*Background(z) component"""
        return self.__sbb
    @SBB.setter 
    def SBB ( self , value ) :
        value = float ( value  )
        assert value in self.__sbb, "Value %s is out of the allowed range %s " % ( value , self.__sbb.minmax() )
        self.__sbb.setVal ( value ) 

    @property
    def BSB ( self ) :
        """The yield of Background(x)*Signal(y)*Background(z) component"""
        return self.__bsb
    @BSB.setter 
    def BSB ( self , value ) :
        value = float ( value  )
        assert value in self.__bsb, "Value %s is out of the allowed range %s " % ( value , self.__bsb.minmax() )
        self.__bsb.setVal ( value ) 

    @property
    def BBS ( self ) :
        """The yield of Background(x)*Background(y)*Signal(z) component"""
        return self.__bbs
    @BBS.setter 
    def BBS ( self , value ) :
        value = float ( value  )
        assert value in self.__bbs, "Value %s is out of the allowed range %s " % ( value , self.__bbs.minmax() )
        self.__bbs.setVal ( value ) 

    @property
    def BBB ( self ) :
        """The yield of Background(x)*Background(y)*Background(z) component"""
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
    def signal3 ( self  ) :
        """Signal(z) component/PDF"""
        return self.__signal3

    @property
    def bkgX1 ( self ) :
        """B(x) component for B(x)*S(y)*S(z) term"""
        return self.__bkgX1
    @property
    def bkgY1 ( self ) :
        """B(y) component for S(x)*B(y)*S(z) term"""
        return self.__bkgY1
    @property
    def bkgZ1 ( self ) :
        """B(z) component for S(x)*S(y)*B(z) term"""
        return self.__bkgZ1
    
    @property
    def bkgXY ( self ) :
        """B(x,y) component for B(x,y)*S(z) term"""
        return self.__bkgXY

    @property
    def bkgXZ ( self ) :
        """B(x,z) component for B(x,z)*S(y) term"""
        return self.__bkgXZ
    
    @property
    def bkgYZ ( self ) :
        """B(y,z) component for S(x)*B(y,z) term"""
        return self.__bkgYZ

    @property
    def bkgX2 ( self ) :
        """B(x) component for B(x)*B(y)*S(z) and B(x)*S(y)*B(z) terms,
        when B(x,z) and B(x,y) are not specified
        """
        return self.__bkgX2
    
    @property
    def bkgY2 ( self ) :
        """B(y) component for B(x)*B(y)*S(z) and S(x)*B(y)*B(z) terms,
        when B(x,y) and B(y,z) are not specified
        """
        return self.__bkgY2
    
    @property
    def bkgZ2 ( self ) :
        """B(z) component for B(x)*S(y)*B(z) and S(x)*B(y)*B(z) terms,
        when B(x,z) and B(y,z) are not specified
        """
        return self.__bkgZ2

    @property
    def bkgX3 ( self ) :
        """B(x) component for B(x)*B(y)*B(z) term,
        when B(x,y,z) is not specified)
        """
        return self.__bkgX3
    @property
    def bkgY3 ( self ) :
        """B(y) component for S(x)*B(y)*B(z) term,
        when B(x,y,z) is not specified)
        """
        return self.__bkgY3
    
    @property
    def bkgZ3 ( self ) :
        """B(z) component for B(x)*B(y)*B(z) term,
        when B(x,y,z) is not specified"""
        return self.__bkgZ3

    @property
    def bkg3D ( self ) :
        """B(x,y,z) component/PDF for the final PDF"""
        return self.__bbb_cmp 

    @property
    def cmp_sss ( self ) :
        """```triple-signal'' component/PDF"""
        return self.__sss_cmp

    @property
    def cpm_ssb ( self ) :
        """```signal-signal-background'' component/PDF"""
        return self.__ssb_cmp
    
    @property
    def cmp_sbs ( self ) :
        """```signal-background-signal'' component/PDF"""
        return self.__sbs_cmp
    
    @property
    def cmp_bss ( self ) :
        """```background-signal-signal'' component/PDF"""
        return self.__bss_cmp

    @property
    def cpm_sbb ( self ) :
        """```signal-background-background'' component/PDF"""
        return self.__sbb_cmp
    
    @property
    def cmp_bsb ( self ) :
        """```background-signal-background'' component/PDF"""
        return self.__bsb_cmp
    
    @property
    def cmp_bbs ( self ) :
        """```background-background-signal'' component/PDF"""
        return self.__bbs_cmp

    @property
    def cmp_bbb ( self ) :
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
        from ostap.logger.utils import roo_silent 
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
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
