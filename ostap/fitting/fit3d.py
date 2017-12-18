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
    )
# =============================================================================
import ROOT
from   ostap.fitting.basic import makeVar, makeBkg
from   ostap.fitting.fit2d import PDF2
from   ostap.logger.utils  import roo_silent 
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
    def __init__ ( self , name , xvar = None , yvar = None , zvar  = None ) : 
        
        
        PDF2.__init__ ( self , name , xvar , yvar ) 
        
        var3  = makeVar ( zvar , 'var3' , '3rd-variable' )
        
        self.varz         = var3
        self.z            = var3 ## ditto
        self.m3           = var3 ## ditto 

    # =========================================================================
    ## make the actual fit (and optionally draw it!)
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
    #  model.m2.setRange ( 'QUQU2' , 2 , 3 ) 
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

        >>> model.m2.setRange ( 'QUQU2' , 2 , 3 ) 
        >>> fx  = model.draw1 ( dataset , nbins = 100 , in_range2 = 'QUQU2') ## draw results
        
        """
        if in_range2 and isinstance ( in_range2 , tuple ) and 2 == len ( in_range2 ) :
            self.m3.setRange ( 'aux_rng2' , in_range2[0] , in_range2[1] )
            in_range2 = 'aux_rng2'

        if in_range3 and isinstance ( in_range3 , tuple ) and 2 == len ( in_ran3e2 ) :
            self.m3.setRange ( 'aux_rng3' , in_range3[0] , in_range3[1] )
            in_range3 = 'aux_rng3'

        in_range = []
        if in_range2 : in_range.append( in_range2 )
        if in_range3 : in_range.append( in_range3 )
        in_ranage = tuple( in_range ) 
        return self.draw ( self.m1  , 
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
    #  model.m1.setRange ( 'QUQU1' , 2 , 3 ) 
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

        >>> model.m1.setRange ( 'QUQU1' , 2 , 3 ) 
        >>> fx  = model.draw2 ( dataset , nbins = 100 , in_range1 = 'QUQU1') ## draw results
        
        """
        if in_range1 and isinstance ( in_range1 , tuple ) and 2 == len ( in_range1 ) :
            self.m1.setRange ( 'aux_rng1' , in_range1[0] , in_range1[1] )
            in_range1 = 'aux_rng1'

        if in_range3 and isinstance ( in_range3 , tuple ) and 2 == len ( in_ran3e2 ) :
            self.m3.setRange ( 'aux_rng3' , in_range3[0] , in_range3[1] )
            in_range3 = 'aux_rng3'

        in_range = []
        if in_range1 : in_range.append( in_range1 )
        if in_range3 : in_range.append( in_range3 )
        in_ranage = tuple( in_range ) 
        return self.draw ( self.m2  , 
                           dataset  ,
                           nbins    ,
                           20       , ## fake 
                           silent   ,
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
    #  model.m2.setRange ( 'QUQU2' , 2 , 3 ) 
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

        >>> model.m2.setRange ( 'QUQU2' , 2 , 3 ) 
        >>> fx  = model.draw3 ( dataset , nbins = 100 , in_range2 = 'QUQU2') ## draw results
        
        """
        if in_range1 and isinstance ( in_range1 , tuple ) and 2 == len ( in_range1 ) :
            self.m1.setRange ( 'aux_rng1' , in_range1[0] , in_range1[1] )
            in_range1 = 'aux_rng1'

        if in_range2 and isinstance ( in_range2 , tuple ) and 2 == len ( in_range2 ) :
            self.m3.setRange ( 'aux_rng2' , in_range2[0] , in_range2[1] )
            in_range2 = 'aux_rng2'

        in_range = []
        if in_range1 : in_range.append( in_range1 )
        if in_range2 : in_range.append( in_range2 )
        in_ranage = tuple( in_range ) 
        return self.draw ( self.m3  , 
                           dataset  ,
                           nbins    ,
                           20       , ## fake 
                           silent   ,
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
        """
        Fit the histogram (and draw it)
        
        >>> histo = ...
        >>> r,f = model.fitHisto ( histo , draw = True )
        
        """
        
        xminmax = histo.xminmax()
        yminmax = histo.yminmax()
        zminmax = histo.zminmax()
        
        with     RangeVar ( self.m1 , *xminmax ) , \
                 RangeVar ( self.m2 , *yminmax ) , \
                 RangeVar ( self.m3 , *zminmax ): 
            
            ## convert it! 
            self.hdset = H3D_dset ( histo , self.m1 , self.m2  , self.m3 , density , silent )
            self.hset  = self.hdset.dset
                
            ## fit it!!
            return self.fitTo ( self.hset      ,
                                draw           ,
                                histo.nbinsx() ,
                                histo.nbinsy() ,
                                histo.nbinsz() ,
                                silent         , *args , **kwargs ) 

    # ====================================================================================
    ## simple 'function-like' interface 
    def __call__ ( self , x , y , z ) :
        
        if     isinstance ( self.m1 , ROOT.RooRealVar ) and \
               isinstance ( self.m2 , ROOT.RooRealVar ) and \
               isinstance ( self.m3 , ROOT.RooRealVar ) :
           
            from ostap.fitting.roofit import SETVAR
            if x in  self.m1 and y in self.m2  and  z in self.m3 : 
                with SETVAR( self.m1 ) , SETVAR( self.m2 ) ,  SETVAR( self.m3 ) :
                    self.m1.setVal ( x )
                    self.m2.setVal ( y )
                    self.m3.setVal ( z )
                    return self.pdf.getVal()
            else : return 0.0
            
        raise AttributeError, 'something wrong goes here'

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
        xmn , xmx = self.m1.minmax()
        ymn , ymx = self.m2.minmax()
        zmn , zmx = self.m3.minmax()

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
                   ## ## additional components 
                   ## components = []    ,
                   name       = ''    ) : 
        
        self._crossterms1 = ROOT.RooArgSet()
        self._crossterms2 = ROOT.RooArgSet()
        
        self.suffix    = suffix 
        self.signal1   = signal_1
        self.signal2   = signal_2
        self.signal3   = signal_3

        
        #
        ## initialize base class
        #
        if not name and signal_1.name and signal_2.name and signal_3.name : 
            name  = signal_1.name
            name += '_x_'
            name += signal_2.name
            name += '_x_'
            name += signal_3.name + '_'+ suffix
            name  = name.replace (' ','_')
            
        PDF3.__init__ ( self          ,
                        name          ,
                        signal_1.mass ,
                        signal_2.mass ,
                        signal_3.mass )

        ## 1) First component: all   signals
        self._ss_list = ROOT.RooArgList ( self.signal1.pdf ,
                                          self.signal2.pdf ,
                                          self.signal3.pdf )
        
        self.sss_pdf = ROOT.RooProdPdf ( "S1S2S3pdf" + suffix ,
                                         "Sig(1) x Sig(2) x Sig(3) "  ,
                                         self._ss_list )
        
        ## 2-4) Three terms:  ( 2 signals )  x ( 1 background ) 

        self._bkgX1 = bkgX1
        self._bkgY1 = bkgY1
        self._bkgZ1 = bkgZ1
        
        self.bkgX1 = makeBkg ( bkgX1 , 'BkgX_S2S3' + suffix , self.m1 )
        self.bkgY1 = makeBkg ( bkgY1 , 'BkgY_S1S3' + suffix , self.m2 )
        self.bkgZ1 = makeBkg ( bkgZ1 , 'BkgZ_S1S2' + suffix , self.m3 )
        
        self._ssb_list = ROOT.RooArgList (
            self.signal1.pdf , self.signal2.pdf , self.bkgZ1  .pdf )
        self._sbs_list = ROOT.RooArgList (
            self.signal1.pdf , self.bkgY1  .pdf , self.signal3.pdf )
        self._bss_list = ROOT.RooArgList (
            self.bkgX1  .pdf , self.signal2.pdf , self.signal3.pdf )

        self.ssb_pdf = ROOT.RooProdPdf ( "S1S2B3pdf" + suffix ,
                                         "Sig(1) x Sig(2) x Bkg(3) "  ,
                                         self._ssb_list )
        self.sbs_pdf = ROOT.RooProdPdf ( "S1B2S3pdf" + suffix ,
                                         "Sig(1) x Bkg(2) x Sig(3) "  ,
                                         self._sbs_list )
        self.bss_pdf = ROOT.RooProdPdf ( "B1S2S3pdf" + suffix ,
                                         "Bkg(1) x Sig(2) x Sig(3) "  ,
                                         self._bss_list )

        ## 5-7) Three terms: (1 signal) x ( 2 backgrounds )
        
        
        self._bkgX2 = bkgX2
        self._bkgY2 = bkgY2
        self._bkgZ2 = bkgZ2

        if bkgX2 is None : bkgX2 = bkgX1
        if bkgY2 is None : bkgY2 = bkgY1
        if bkgZ2 is None : bkgZ2 = bkgZ1
        
        self.bkgX2  = makeBkg ( bkgX2 , 'BkgX_S2' + suffix , self.m1 )        
        self.bkgY2  = makeBkg ( bkgY2 , 'BkgY_S2' + suffix , self.m2 )        
        self.bkgZ2  = makeBkg ( bkgZ2 , 'BkgZ_S2' + suffix , self.m3 )
        
        self._bkgXY = bkgXY
        self._bkgXZ = bkgXZ
        self._bkgYZ = bkgYZ

        if   bkgXY and isinstance ( bkgXY , ROOT.RooAbsPdf ) :
            self.bkgXY = bkgXY            
        elif bkgXY and hasattr ( bkgXY , 'pdf' ) and isintance ( bkgXY.pdf  , ROOT.RooAbsPdf ) :
            self.bkgXY = bkgXY.pdf
        else :
            self.bkgXY = ROOT.RooProdPdf ( 'BkgXY_pdf' + suffix , 'Bkg(1)xBkg(2)' , 
                                           self.bkgX2.pdf , self.bkgY2.pdf )

        if   bkgXZ and isinstance ( bkgXZ , ROOT.RooAbsPdf ) :
            self.bkgXZ = bkgXZ            
        elif bkgXZ and hasattr ( bkgXZ , 'pdf' ) and isintance ( bkgXZ.pdf  , ROOT.RooAbsPdf ) :
            self.bkgXZ = bkgXZ.pdf
        else :
            self.bkgXZ = ROOT.RooProdPdf ( 'BkgXZ_pdf' + suffix , 'Bkg(1)xBkg(3)' , 
                                           self.bkgX2.pdf , self.bkgZ2.pdf )

        if   bkgYZ and isinstance ( bkgYZ , ROOT.RooAbsPdf ) :
            self.bkgYZ = bkgYZ            
        elif bkgYZ and hasattr ( bkgYZ , 'pdf' ) and isintance ( bkgYZ.pdf  , ROOT.RooAbsPdf ) :
            self.bkgYZ = bkgYZ.pdf
        else :
            self.bkgYZ = ROOT.RooProdPdf ( 'BkgYZ_pdf' + suffix , 'Bkg(2)xBkg(3)' , 
                                           self.bkgY2.pdf , self.bkgZ2.pdf )


        self.sbb_pdf = ROOT.RooProdPdf ( "S1B2B3pdf" + suffix ,
                                         "Sig(1) x Bkg(2) x Bkg(3) "   ,
                                         self.signal1.pdf , self.bkgYZ )
        self.bsb_pdf = ROOT.RooProdPdf ( "B1S2B3pdf" + suffix ,
                                         "Bkg(1) x Sig(2) x Bkg(3) "   ,
                                         self.signal2.pdf , self.bkgXZ )
        self.bbs_pdf = ROOT.RooProdPdf ( "B1B2S3pdf" + suffix ,
                                         "Bkg(1) x Bkg(2) x Sig(3) "   ,
                                         self.signal3.pdf , self.bkgXY )
        
        ##  8) the last term: all background 
        self._bkg3D_ = bkg3D
        
        if   bkg3D and isinstance ( bkg3D , ROOT.RooAbsPdf ) : self.bbb_pdf = bkg3D 
        elif bkg3D and hasattr    ( bkg3D , 'pdf'          ) : self.bbb_pdf = bkg3D.pdf
        else     :            

            
            if bkgX3 is None : bkgX3 = bkgX2
            if bkgY3 is None : bkgY3 = bkgY2
            if bkgZ3 is None : bkgZ3 = bkgZ2
            
            self.bkgX3 = makeBkg ( bkgX3 , 'BkgX_S0' + suffix , self.m1 )
            self.bkgY3 = makeBkg ( bkgY3 , 'BkgY_S0' + suffix , self.m2 )
            self.bkgZ3 = makeBkg ( bkgZ3 , 'BkgZ_S0' + suffix , self.m3 )
            
            self._bbb_list = ROOT.RooArgList ( self.bkgX3.pdf ,
                                               self.bkgY3.pdf ,
                                               self.bkgZ3.pdf )
            
            self.bbb_pdf = ROOT.RooProdPdf ( "B1B2B3pdf" + suffix ,
                                             "Bkg(1) x Bkg(2) x Bkg(3) "  ,
                                             self._bbb_list )

        #
        ## coefficients
        #
        self.sss = makeVar ( sss   ,
                             "S1S2S3"               + suffix ,
                             "Sig(1)&Sig(2)&Sig(3)" + suffix , None , 1000  , 0 ,  1.e+8 )
        self.ssb = makeVar ( ssb   ,
                             "S1S2B3"               + suffix ,
                             "Sig(1)&Sig(2)&Bkg(3)" + suffix , None , 1000  , 0 ,  1.e+8 )
        self.sbs = makeVar ( sbs   ,
                             "S1B2S3"               + suffix ,
                             "Sig(1)&Bkg(2)&Sig(3)" + suffix , None , 1000  , 0 ,  1.e+8 )
        self.bss = makeVar ( bss   ,
                             "B1S2S3"               + suffix ,
                             "Bkg(1)&Sig(2)&Sig(3)" + suffix , None , 1000  , 0 ,  1.e+8 )
        self.sbb = makeVar ( sbb  ,
                             "S1B2B3"               + suffix ,
                             "Sig(1)&Bkg(2)&Bkg(3)" + suffix , None , 1000  , 0 ,  1.e+8 )
        self.bsb = makeVar ( bsb  ,
                             "B1S2B3"               + suffix ,
                             "Bkg(1)&Sig(2)&Bkg(3)" + suffix , None , 1000  , 0 ,  1.e+8 )
        self.bbs = makeVar ( bbs  ,
                             "B1B2S3"               + suffix ,
                             "Bkg(1)&Bkg(2)&Sig(3)" + suffix , None , 1000  , 0 ,  1.e+8 )
        self.bbb = makeVar ( bbb  ,
                             "B1B2B3"               + suffix ,
                             "Bkg(1)&Bkg(2)&Bkg(3)" + suffix , None , 1000  , 0 ,  1.e+8 )
        
        self.SSS_name = self.sss.GetName()
        self.SSB_name = self.ssb.GetName()
        self.SBS_name = self.sbs.GetName()
        self.BSS_name = self.bss.GetName()
        self.SBB_name = self.sbb.GetName()
        self.BSB_name = self.bsb.GetName()
        self.BBS_name = self.bbs.GetName()
        self.BBB_name = self.bbb.GetName()
        
        self.alist1 = ROOT.RooArgList (
            self.sss_pdf ,
            self.ssb_pdf ,
            self.sbs_pdf ,
            self.bss_pdf ,
            self.sbb_pdf ,
            self.bsb_pdf ,
            self.bbs_pdf ,
            self.bbb_pdf )
        self.alist2 = ROOT.RooArgList (
            self.sss     ,
            self.ssb     ,
            self.sbs     ,
            self.bss     ,
            self.sbb     ,
            self.bsb     ,
            self.bbs     ,
            self.bbb     )
        
        
        #
        ## build the final PDF 
        # 
        self.pdf  = ROOT.RooAddPdf  ( "model3D"      + suffix ,
                                      "Model3D(%s)"  % suffix ,
                                      self.alist1 ,
                                      self.alist2 )
        
        self.signals     ().add ( self.sss_pdf )
        self.backgrounds ().add ( self.bbb_pdf )
        self.crossterms1 ().add ( self.ssb_pdf ) ## cross-terms
        self.crossterms1 ().add ( self.sbs_pdf ) ## cross-terms
        self.crossterms1 ().add ( self.bss_pdf ) ## cross-terms
        self.crossterms2 ().add ( self.ssb_pdf ) ## cross-terms 
        self.crossterms2 ().add ( self.sbs_pdf ) ## cross-terms 
        self.crossterms2 ().add ( self.bss_pdf ) ## cross-terms 
        
    ## get all declared components 
    def crossterms1 ( self ) : return self._crossterms1
    ## get all declared components 
    def crossterms2 ( self ) : return self._crossterms2


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
    """ Wrapper for generic (3D) RooFit pdf
    # 
    # raw_pdf = 
    # pdf     = Generic3D_pdf ( raw_pdf )
    # 
    """
    ## constructor 
    def __init__ ( self , pdf , varx = None , vary = None , varz = None , name = None ) :
        if not name : name = pdf.GetName()
        PDF3  . __init__ ( self , name , varx , vary , varz )
        self.pdf = pdf

# =============================================================================
## simple convertor of 3D-histo to data set
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H3D_dset(object) :
    """Simple convertor of 3D-histogram into data set
    """
    def __init__ ( self            ,
                   histo3          ,
                   mass    = None  ,
                   mass2   = None  ,
                   mass3   = None  ,
                   density = True  ,
                   silent  = False ) :
        #
        ## use mass-variable
        #
        name         = histo3.GetName() 
        self.mass    = makeVar ( mass  , 'm_%s'  % name , 'mass (%s)' % name , None , *(histo3.xminmax()) )
        self.mass1   = self.mass 
        self.mass2   = makeVar ( mass2 , 'm2_%s' % name , 'mass2(%s)' % name , None , *(histo3.yminmax()) )
        self.mass3   = makeVar ( mass3 , 'm3_%s' % name , 'mass2(%s)' % name , None , *(histo3.zminmax()) )

        self.impDens = density 
        self.var1    = self.mass1
        self.var2    = self.mass2
        self.var3    = self.mass3
        self.x       = self.var1 
        self.y       = self.var2
        self.z       = self.varz
        
        with roo_silent ( silent ) : 

            self.vlst  = ROOT.RooArgList    ( self.mass1 , self.mass2 , self.mass3 )
            self.vimp  = ROOT.RooFit.Import ( histo3 , density )
            self.dset  = ROOT.RooDataHist   (
                rootID ( 'hds_' ) ,
                "Data set for histogram '%s'" % histo3.GetTitle() ,
                self.vlst  ,
                self.vimp  )


# =============================================================================
## simple convertor of 3D-histogram into PDF
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H3D_pdf(H3D_dset,PDF3) :
    """Simple convertor of 3D-histogram into PDF 
    """
    def __init__ ( self            ,
                   name            ,
                   histo3          ,
                   mass    = None  , 
                   mass2   = None  ,
                   mass3   = None  ,
                   density = True  ,
                   silent  = False ) :
        
        H3D_dset.__init__ ( self , histo3 ,      mass  ,      mass2 , mass3 ,  density , silent )
        PDF3    .__init__ ( self , name   , self.mass1 , self.mass2 , self.mass3 ) 

        self.vset  = ROOT.RooArgSet  ( self.mass , self.mass2 , self.mass3 )
        
        #
        ## finally create PDF :
        #
        from   Ostap.Utils      import roo_silent 
        with roo_silent ( silent ) : 
            self.pdf    = ROOT.RooHistPdf (
                'hpdf_%s'            % name ,
                'Histo3PDF(%s/%s/%s)' % ( name , histo3.GetName() , histo2.GetTitle() ) , 
                self.vset  , 
                self.dset  )
        
# =============================================================================
if '__main__' == __name__ :
    
    from Ostap.Line import line 
    logger.info ( __file__ + '\n' + line  )
    logger.info ( 80*'*' )
    logger.info ( __doc__  )
    logger.info ( 80*'*' )
    logger.info ( ' Author  : %s' %         __author__    ) 
    logger.info ( ' Version : %s' %         __version__   ) 
    logger.info ( ' Date    : %s' %         __date__      )
    logger.info ( ' Symbols : %s' %  list ( __all__     ) )
    logger.info ( 80*'*' ) 

# =============================================================================
# The END 
# =============================================================================
