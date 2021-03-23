#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/efficiency.py
#  Set of useful basic utilities to fit "efficiency" 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
# =============================================================================
"""Set of useful basic utilities to fit ``efficiency''"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'Efficiency1D', ## helper utility to get the efficiency (1D-case)
    'Efficiency2D', ## helper utility to get the efficiency (2D-case)
    'Efficiency3D', ## helper utility to get the efficiency (3D-case)
    )
# =============================================================================
import ROOT
from   ostap.fitting.funbasic import FUNC, FUNC2 , FUNC3 , Fun1D , Fun2D , Fun3D 
from   ostap.fitting.basic    import PDF , Generic1D_pdf
from   ostap.fitting.fit2d    import PDF2, Generic2D_pdf
from   ostap.fitting.fit3d    import PDF3, Generic3D_pdf
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.efficiency' )
else                       : logger = getLogger( __name__                   )
# ==============================================================================
logger.debug ( 'Utilities to get efficiency using RooFit machinery')
# ==============================================================================
## @class Efficiency1D
#  Helper class  to get get the efficiency using unbinned fit and ROOT.RooEfficiency class
#  @attention this an internal class, not to be used directly!
class Efficiency ( object ) :
    """ Helper class to get the efficiency using unbinned fit and ROOT.RooEfficiency class
    - attention : this an internal class, not to be used directly!
    """
    def __init__  ( self              ,
                    name              ,
                    eff_pdf           ,                    
                    eff_fun           ,
                    vars              ,
                    cut               ,
                    accept = 'accept' ,
                    scale  = None     ) :

        self.__name   = str(name)
        self.__cut    = cut
        self.__accept = str(accept) 

        assert eff_pdf or eff_fun ,'Function or PDF must be specified!'
        
        assert isinstance ( cut , ROOT.RooCategory ) , "``Cut'' is not RooCategory!"

        self.__eff_pdf = eff_pdf
        self.__eff_fun = eff_fun
        self.__vars    = ROOT.RooArgSet( *vars )
        
        self.__scale   = 1
        
        if not self.eff_fun :

            if scale is None :
                scale = 0.001 , 1.e-10 , 1.e+5 
                mnmx  = self.eff_pdf.minmax ()
                if mnmx :
                    mn , mx = mnmx
                    scale   = 0.25 / mx , 1.e-9 / mx , 10 / mx
                    
            if   isinstance ( scale , ROOT.RooAbsReal ) :
                self.__scale = scale
            else : 
                self.__scale = ROOT.RooRealVar ( 'effscale_%s' % self.name , 'scale factor for efficiency (%s)' % self.name , *scale )

            self.__lst     = ROOT.RooArgList ( self.__scale , self.__eff_pdf.pdf )
            _s = self.scale.GetName()
            _p = self.eff_pdf.pdf.GetName() 
            self.__eff_fun = ROOT.RooFormulaVar (
                'Eff_%s' % self.name , '%s*%s'  % ( _s , _p ) , self.__lst )

        ## create the main PDF: RooEfficiency 
        self.__pdf =  ROOT.RooEfficiency (
            PDF.roo_name ( 'eff_' )       ,
            "Efficiency  %s"  % self.name ,
            self.eff_fun                  ,
            self.cut                      ,
            self.accept                   )

        if 3 == len ( vars ) : 
            ## pdf-object for fit 
            self.__pdf_fit = Generic3D_pdf ( pdf   = self.pdf ,
                                             xvar  = vars[0]  ,
                                             yvar  = vars[1]  ,
                                             zvar  = vars[2]  ,
                                             name  = PDF.generate_name ( 'eff_fit_%s'   % self.name ) ,
                                             special        = True  ,
                                             add_to_signals = False )
            ## pdf-object for drawing
            self.__pdf_draw = Generic3D_pdf ( pdf   = self.eff_fun   ,
                                              xvar  = vars[0]        ,
                                              yvar  = vars[1]        ,
                                              zvar  = vars[2]        ,
                                              name  = PDF.generate_name ( 'eff_draw_%s'  % self.name ) ,
                                              special        = True  ,
                                              add_to_signals = False )
        elif 2 == len (  vars ) :
            ## pdf-object for fit 
            self.__pdf_fit  = Generic2D_pdf ( pdf   = self.pdf ,
                                              xvar  = vars[0]  ,
                                              yvar  = vars[1]  ,
                                              name  = PDF.generate_name ( 'eff_fit_%s'   % self.name ) ,
                                              special        = True  ,
                                              add_to_signals = False )
            ## pdf-object for drawing
            self.__pdf_draw = Generic2D_pdf ( pdf   = self.eff_fun   ,
                                              xvar  = vars[0]        ,
                                              yvar  = vars[1]        ,
                                              name  = PDF.generate_name ( 'eff_draw_%s'  % self.name ) ,
                                              special        = True  ,
                                              add_to_signals = False )
        elif 1 == len (  vars ) :
            ## pdf-object for fit 
            self.__pdf_fit  = Generic1D_pdf ( pdf   = self.pdf ,
                                              xvar  = vars[0]  ,
                                              name  = PDF.generate_name ( 'eff_fit_%s'   % self.name ) ,
                                              special        = True  ,
                                              add_to_signals = False )
            ## pdf-object for drawing
            self.__pdf_draw = Generic1D_pdf ( pdf   = self.eff_fun   ,
                                              xvar  = vars[0]        ,
                                              name  = PDF.generate_name ( 'eff_draw_%s'  % self.name ) ,
                                              special        = True  ,
                                              add_to_signals = False )
        else :
            raise AttributeError("Invalid length of vars: %s" % str( vars ) )
        
        ## th fit results from the last fit
        self.__fit_result = None
        
    @property 
    def name ( self ) :
        """``name'': the name of Efficiency object"""
        return self.__name

    @property
    def vars  ( self ) :
        """``vars'' - list of variables (RooArgSet)"""
        return self.__vars
    
    @property 
    def scale ( self ) :
        """``scale'': the scale factor to be applied to eff_pdf"""
        return self.__scale
    @scale.setter
    def scale ( self , value ) :
        val = float ( value )
        self.__scale.setVal( val )
    
    @property 
    def cut ( self ) :
        """``cut'' : RooCategory to describe ``accepted'' and ``rejected'' samples"""
        return self.__cut
    
    @property 
    def accept ( self ) :
        """``accept'' : the name of category for ``accepted'' sample"""
        return str(self.__accept)
        
    @property
    def eff_pdf ( self ) :
        """``eff_pdf'' : efficiency PDF (could be fictive)"""
        return  self.__eff_pdf
    
    @property
    def eff_fun ( self ) :
        """``eff_fun'' : efficiency function"""
        return  self.__eff_fun

    @property
    def pdf    ( self ) :
        """``pdf''  : the actual   RooEfficiency object"""
        return self.__pdf

    @property
    def pdf_fit  ( self ) :
        """```pdf_fit'' : the actual PDF-object used for fit"""
        return self.__pdf_fit
    
    @property
    def pdf_draw ( self ) :
        """```pdf_draw'' : the actual PDF-object used for drawing"""
        return self.__pdf_draw 

    @property
    def fit_result ( self ) :
        """``fit_result'' :  fit  results  from the last fit"""
        return self.__fit_result
    
    # =========================================================================
    ## fit it!
    #  @code
    #  dataset = ...
    #  eff     = Efficiency1D ( ... ) 
    #  result  = eff.fitTo ( dataset ) 
    #  @endcode 
    def  fitTo ( self , dataset , silent = True , refit = 2 , args = () , **kwargs )  :
        """Fit it!

        >>> dataset = ...
        >>> eff     = Efficiency1D( ... ) 
        >>> result  = eff.fitTo ( dataset ) 
        """
        
        vargs = [] ##  i for i in args ]
        vargs.append  ( ROOT.RooFit.ConditionalObservables ( self.vars ) ) 
        vargs = tuple ( vargs ) 
        ##
        draw  = kwargs.pop ( 'draw' , False ) 
        from ostap.logger.utils import roo_silent
        with roo_silent ( silent ) : 
            result , frame = self.pdf_fit.fitTo ( dataset ,                            
                                                  draw   = False  ,
                                                  nbins  = 100    ,
                                                  silent = silent ,
                                                  refit  = refit  ,
                                                  args   = vargs  , **kwargs )
            if  draw : self.draw ( dataset ) 
        self.__fit_result = result 
        return result

    # =========================================================================
    ## draw the efficiency (and the dataset)
    #  @code
    #  dataset = ... 
    #  eff = Efficiency1D( ... )
    #  eff.fitTo ( dataset ) 
    #  eff.draw  ( dataset ) 
    #  @endcode 
    def draw ( self ,
               dataset = None  ,
               nbins   = 100   ,
               silent  = False ,
               args    = ()    , **kwargs ) :
        """Draw the efficiency (and the dataset)
        >>> dataset = ... 
        >>> eff     = Efficiency1D( ... )
        >>> eff.fitTo ( dataset ) 
        >>> eff.draw  ( dataset ) 
        """
        
        ddopts  = kwargs.get('data_options',() )
        ddopts  = [  o for o in ddopts ] 
        ddopts.append ( ROOT.RooFit.Efficiency ( self.cut ) )
        ddopts  =  tuple ( ddopts )
        kwargs['data_options'] = ddopts
        
        return self.pdf_draw.draw ( dataset         ,
                                    nbins  = nbins  ,
                                    silent = silent ,
                                    args   = args   , **kwargs )
    

# =============================================================================
## @class Efficiency1D
#  Get the efficiency using unbinned fit and ROOT.RooEfficiency class
#  @code
#  cut = ROOT.RooCategory( .... )
#  cut.defineType ( 'accept' , 1 )
#  cut.defineType ( 'reject' , 0 )
#  ...
#  efficiency = ...  ## RooFit function or pdf (RooAbsReal or RooAbdPdf of PDF )
#  eff = Efficiency1D ( 'eff' , efficiency , cut , xvar = x )
#  dataset = .... ##  dataset
#  r = eff.fitTo ( dataset ) ## fit it! 
#  eff.draw( dataset )       ## draw it!
#  v = eff ( 0.115 )         ## use it! 
#  @endcode
class Efficiency1D (Efficiency) :
    """ Get the efficiency using unbinned fit and ROOT.RooEfficiency class
    
    >>> cut = ROOT.RooCategory( .... )
    >>> cut.defineType ( 'accept' , 1 )
    >>> cut.defineType ( 'reject' , 0 )

    >>> efficiency = ...  ## RooFit function or pdf (RooAbsReal or RooAbdPdf of PDF )
    
    >>> eff = Efficiency1D ( 'eff' , efficiency , cut , xvar = x )

    >>> dataset = .... ##  dataset

    Fit it:
    >>> r = eff.fitTo ( dataset ) ## fit it! 

    Draw it: 
    >>> eff.draw( dataset )       ## draw it!

    Use it:
    >>> value = eff ( 0.115 )     ## use it! 
    """
    
    def __init__  ( self              ,
                    name              ,
                    efficiency        , ## the function, FUNC or PDF 
                    cut               ,
                    xvar   = None     ,
                    accept = 'accept' ,
                    scale  = None     ) :

        self.__eff = efficiency
        
        if   isinstance   ( efficiency , PDF  ) :
            eff_pdf = efficiency
            xvar    = efficiency.xvar
            eff_fun = None
        elif isinstance   ( efficiency , FUNC ) :
            eff_fun = efficiency.fun 
            xvar    = efficiency.xvar
            eff_pdf = None  
        elif isinstance ( efficiency , ROOT.RooAbsPdf   ) and xvar and isinstance ( xvar , ROOT.RooAbsReal ) :
            eff_pdf = Generic1D_pdf ( efficiency , xvar )                
            eff_fun = None 
        elif isinstance ( efficiency , ROOT.RooAbsReal  ) and xvar and isinstance ( xvar , ROOT.RooAbsReal ) :            
            eff_pdf = Fun1D        ( efficiency , xvar =  xvar )                
            eff_fun = efficiency            
        else :
            raise AttributeError('Invalid efficiency/xvar combination  %s/%s;%s/%s'  %
                                 ( efficiency ,  type(efficiency) , xvar , type(xvar) ) )

        self.__xvar = xvar 
        Efficiency.__init__ ( self , name , eff_pdf , eff_fun ,  ( xvar,) , cut , accept , scale )
    
    @property 
    def xvar ( self ) :
        """``xvar'': the x-variable """
        return self.__xvar

    # =========================================================================
    ## draw the efficiency (and the dataset)
    #  @code
    #  dataset = ... 
    #  eff = Efficiency1D( ... )
    #  eff.fitTo ( dataset ) 
    #  eff.draw  ( dataset ) 
    #  @endcode 
    def draw ( self , dataset = None , *args ,  **kwargs ) :
        """Draw the efficiency (and the dataset)
        >>> dataset = ... 
        >>> eff     = Efficiency1D( ... )
        >>> eff.fitTo ( dataset ) 
        >>> eff.draw  ( dataset ) 
        """
        
        ddopts  = kwargs.get('data_options',() )
        ddopts  = [  o for o in ddopts ] 
        ddopts.append ( ROOT.RooFit.Efficiency ( self.cut ) )
        ddopts  =  tuple ( ddopts )
        kwargs['data_options'] = ddopts
        
        return self.pdf_draw.draw ( dataset , *args , **kwargs )
    
    # =========================================================================
    ## get the efficiency
    #  @code
    #  dataset = ... 
    #  eff = Efficiency1D( ... )
    #  eff.fitTo ( dataset ) 
    #  x = 0.15
    #  value =  eff(x) 
    #  @endcode     
    def __call__ (  self , x , error = False ) :
        """Get the efficiency
        >>> dataset = ... 
        >>> eff = Efficiency1D( ... )
        >>> eff.fitTo ( dataset ) 
        >>> x = 0.15
        >>> value = eff(x) 
        """
        from ostap.fitting.roofit import SETVAR
        from ostap.math.ve        import VE 
        xx = float ( x ) 
        if xx in self.xvar : 
            with  SETVAR ( self.xvar ) :
                self.xvar.setVal ( xx )                
                v = self.eff_fun.getVal ()
                if error and self.fit_result :
                    e = self.eff_fun.getPropagatedError ( self.fit_result )
                    if 0<= e : return  VE ( v ,  e * e )
                return v 
        logger.error ('Invalid efficiency, return -1 ') 
        return -1 

# =============================================================================
## @class Efficiency2D
#  Get the efficiency using unbinned fit and ROOT.RooEfficiency class
#  @code
#  cut = ROOT.RooCategory( .... )
#  cut.defineType ( 'accept' , 1 )
#  cut.defineType ( 'reject' , 0 )
#  ...
#  efficiency = ...  ## RooFit function or pdf (RooAbsReal or RooAbdPdf of PDF )
#  eff = Efficiency2D ( 'eff' , efficiency , cut , xvar = x , yvar = ... )
#  dataset = .... ##  dataset
#  r = eff.fitTo ( dataset ) ## fit it! 
#  eff.draw( dataset )       ## draw it!
#  v = eff ( 0.115 , 90.0 )         ## use it! 
#  @endcode
class Efficiency2D (Efficiency) :
    """ Get the efficiency using unbinned fit and ROOT.RooEfficiency class
    
    >>> cut = ROOT.RooCategory( .... )
    >>> cut.defineType ( 'accept' , 1 )
    >>> cut.defineType ( 'reject' , 0 )

    >>> efficiency = ...  ## RooFit function or pdf (RooAbsReal or RooAbdPdf of PDF )
    
    >>> eff = Efficiency1D ( 'eff' , efficiency , cut , xvar = x )

    >>> dataset = .... ##  dataset

    Fit it:
    >>> r = eff.fitTo ( dataset ) ## fit it! 

    Draw it: 
    >>> eff.draw( dataset )       ## draw it!

    Use it:
    >>> value = eff ( 0.115 , 90.0 )     ## use it! 
    """
    
    def __init__  ( self              ,
                    name              ,
                    efficiency        , ## the function or PDF 
                    cut               ,
                    xvar   = None     ,
                    yvar   = None     ,
                    accept = 'accept' ) :

        if isinstance   ( efficiency , PDF2  ) :            
            eff_pdf = efficiency
            xvar    = efficiency.xvar
            eff_fun = None            
        elif isinstance   ( efficiency , FUNC2 ) :
            eff_fun = efficiency.fun 
            xvar    = efficiency.xvar
            yvar    = efficiency.yvar
            eff_pdf = None  
        elif isinstance ( efficiency , ROOT.RooAbsReal  ) :            
            okx = xvar and isinstance ( xvar , ROOT.RooAbsReal )
            oky = yvar and isinstance ( yvar , ROOT.RooAbsReal )
            assert oix and oky, 'Invalid efficiency/xvar/yvar setting!'            
            eff_pdf = Generic2D_pdf ( efficiency , xvar , yvar , special = True )                
            eff_fun = None  if isinstance ( efficiency , ROOT.RooAbsPdf ) else efficiency            
        elif isinstance ( efficiency , ROOT.RooAbsPdf   ) :            
            okx = xvar and isinstance ( xvar , ROOT.RooAbsReal )
            oky = yvar and isinstance ( yvar , ROOT.RooAbsReal )
            assert oix and oky, 'Invalid efficiency/xvar/yvar setting!'            
            eff_pdf = Generic2D_pdf ( efficiency , xvar , yvar , special = True )                
            eff_fun = None       
        elif isinstance ( efficiency , ROOT.RooAbsReal  ) :            
            okx = xvar and isinstance ( xvar , ROOT.RooAbsReal )
            oky = yvar and isinstance ( yvar , ROOT.RooAbsReal )
            assert oix and oky, 'Invalid efficiency/xvar/yvar setting!'            
            eff_pdf = Fun2D ( efficiency , xvar = xvar , yvar = yvar )  
            eff_fun = efficiency 
        else :
            raise AttributeError('Invalid efficiency/xvar/yvat combination  %s/%s/%s/%s'  %
                                 ( efficiency ,  type(efficiency) , xvar , yvar ) )

        
        self.__xvar = xvar 
        self.__yvar = yvar 
        Efficiency.__init__ ( self , name , eff_pdf , eff_fun ,  ( xvar , yvar) , cut , accept )

    @property 
    def xvar ( self ) :
        """``xvar'': the x-variable """
        return self.__xvar
    @property 
    def yvar ( self ) :
        """``yvar'': the y-variable """
        return self.__yvar

    # =========================================================================
    ## draw the efficiency (and the dataset)
    #  @code
    #  dataset = ... 
    #  eff     = ... 
    #  eff.fitTo ( dataset ) 
    #  eff.draw1 ( dataset ) 
    #  @endcode 
    def draw1 ( self           ,
                dataset = None ,
                nbins   = 100  ,
                silent  = True ,  *args , **kwargs ) :
        """Draw the efficiency (and the dataset)
        >>> dataset = ... 
        >>> eff     = ... 
        >>> eff.fitTo ( dataset ) 
        >>> eff.draw  ( dataset ) 
        """
        
        ddopts  = kwargs.get('data_options',() )
        ddopts  = [  o for o in ddopts ] 
        ddopts.append ( ROOT.RooFit.Efficiency ( self.cut ) )
        ddopts  =  tuple ( ddopts )
        kwargs['data_options'] = ddopts
        
        return self.pdf_draw.draw1 ( dataset , nbins , silent , *args , **kwargs )
    
    # =========================================================================
    ## draw the efficiency (and the dataset)
    #  @code
    #  dataset = ... 
    #  eff     = ... 
    #  eff.fitTo ( dataset ) 
    #  eff.draw2 ( dataset ) 
    #  @endcode 
    def draw2 ( self           ,
                dataset = None ,
                nbins   = 100  ,
                silent  = True ,  *args , **kwargs ) :
        """Draw the efficiency (and the dataset)
        >>> dataset = ... 
        >>> eff     = ... 
        >>> eff.fitTo ( dataset ) 
        >>> eff.draw  ( dataset ) 
        """
        
        ddopts  = kwargs.get('data_options',() )
        ddopts  = [  o for o in ddopts ] 
        ddopts.append ( ROOT.RooFit.Efficiency ( self.cut ) )
        ddopts  =  tuple ( ddopts )
        kwargs['data_options'] = ddopts
        
        return self.pdf_draw.draw2 ( dataset , nbins , silent , *args , **kwargs )
        
    # =========================================================================
    ## get the efficiency
    #  @code
    #  dataset = ... 
    #  eff     = ... 
    #  eff.fitTo ( dataset ) 
    #  x = 0.15, y = 90. 
    #  value =  eff(x,y) 
    #  @endcode     
    def __call__ (  self , x , y , error = False ) :
        """Get the efficiency
        >>> dataset = ... 
        >>> eff     = ... 
        >>> eff.fitTo ( dataset ) 
        >>> x = 0.15 , y = 90 
        >>> value = eff(x,y) 
        """
        from ostap.fitting.roofit import SETVAR
        from ostap.math.ve        import VE 
        xx = float ( x ) 
        yy = float ( y ) 
        if xx in self.xvar and yy in self.yvar : 
            with  SETVAR ( self.xvar ) :
                with SETVAR ( self.yvar ) :
                    self.xvar.setVal ( xx )
                    self.yvar.setVal ( yy )
                    v = self.eff_fun.getVal ()
                    if error and self.fit_result :
                        e = self.eff_fun.getPropagatedError ( self.fit_result )
                        if 0<= e : return  VE ( v ,  e * e )
                    return v 
        logger.error ('Invalid efficiency, return -1 ') 
        return -1 


# =============================================================================
## @class Efficiency3D
#  Get the efficiency using unbinned fit and ROOT.RooEfficiency class
#  @code
#  cut = ROOT.RooCategory( .... )
#  cut.defineType ( 'accept' , 1 )
#  cut.defineType ( 'reject' , 0 )
#  ...
#  efficiency = ...  ## RooFit function or pdf (RooAbsReal or RooAbdPdf of PDF )
#  eff = Efficiency2D ( 'eff' , efficiency , cut , xvar = x , yvar = ... ,  zvar  =  ...  )
#  dataset = .... ##  dataset
#  r = eff.fitTo ( dataset ) ## fit it! 
#  eff.draw( dataset )       ## draw it!
#  v = eff ( 0.115 , 90.0 , 15 )         ## use it! 
#  @endcode
class Efficiency3D (Efficiency) :
    """ Get the efficiency using unbinned fit and ROOT.RooEfficiency class
    
    >>> cut = ROOT.RooCategory( .... )
    >>> cut.defineType ( 'accept' , 1 )
    >>> cut.defineType ( 'reject' , 0 )

    >>> efficiency = ...  ## RooFit function or pdf (RooAbsReal or RooAbdPdf of PDF )
    
    >>> eff = Efficiency3D ( 'eff' , efficiency , cut , xvar = x , yvar = y , zvar = z)

    >>> dataset = .... ##  dataset

    Fit it:
    >>> r = eff.fitTo ( dataset ) ## fit it! 

    Draw it: 
    >>> eff.draw( dataset )       ## draw it!

    Use it:
    >>> value = eff ( 0.115 , 90.0 , 15 )     ## use it! 
    """
    
    def __init__  ( self              ,
                    name              ,
                    efficiency        , ## the function or PDF 
                    cut               ,
                    xvar   = None     ,
                    yvar   = None     ,
                    zvar   = None     ,
                    accept = 'accept' ) :

        if isinstance   ( efficiency , PDF3  ) :            
            eff_pdf = efficiency
            xvar    = efficiency.xvar
            yvar    = efficiency.yvar
            zvar    = efficiency.zvar
            eff_fun = None            
        elif isinstance   ( efficiency , FUNC3 ) :
            eff_fun = efficiency.fun 
            xvar    = efficiency.xvar
            yvar    = efficiency.yvar
            zvar    = efficiency.zvar
            eff_pdf = None  
        elif isinstance ( efficiency , ROOT.RooAbsPdf  ) :            
            okx = xvar and isinstance ( xvar , ROOT.RooAbsReal )
            oky = yvar and isinstance ( yvar , ROOT.RooAbsReal )
            okz = zvar and isinstance ( zvar , ROOT.RooAbsReal )
            assert oix and oky and okz, 'Invalid efficiency/xvar/yvar/zvar setting!'            
            eff_pdf = Generic3D_pdf ( efficiency , xvar , yvar , zvar , special = True )
            eff_fun = None  
        elif isinstance ( efficiency , ROOT.RooAbsReal  ) :            
            okx = xvar and isinstance ( xvar , ROOT.RooAbsReal )
            oky = yvar and isinstance ( yvar , ROOT.RooAbsReal )
            okz = zvar and isinstance ( zvar , ROOT.RooAbsReal )
            assert oix and oky and okz, 'Invalid efficiency/xvar/yvar/zvar setting!'            
            eff_pdf = Fun3D ( efficiency , xvar = xvar , yvar = yvar , zvar = zvar )                
            eff_fun = efficiency            
        else :
            raise AttributeError('Invalid efficiency/xvar/yvar/zvar combination  %s/%s/%s/%s/%s'  %
                                 ( efficiency ,  type(efficiency) , xvar , yvar ,  zvar ) )
        
        
        self.__xvar = xvar 
        self.__yvar = yvar 
        self.__zvar = zvar 
        Efficiency.__init__ ( self , name , eff_pdf , eff_fun ,  ( xvar , yvar , zvar ) , cut , accept )

    @property 
    def xvar ( self ) :
        """``xvar'': the x-variable """
        return self.__xvar
    @property 
    def yvar ( self ) :
        """``yvar'': the y-variable """
        return self.__yvar
    @property 
    def zvar ( self ) :
        """``zvar'': the z-variable """
        return self.__zvar

    # =========================================================================
    ## draw the efficiency (and the dataset)
    #  @code
    #  dataset = ... 
    #  eff     = ... 
    #  eff.fitTo ( dataset ) 
    #  eff.draw1 ( dataset ) 
    #  @endcode 
    def draw1 ( self            ,
                dataset   = None ,
                nbins     = 100  ,
                silent    = True ,
                in_range2 = None ,
                in_range3 = None ,
                *args , **kwargs ) :
        """Draw the efficiency (and the dataset)
        >>> dataset = ... 
        >>> eff     = ... 
        >>> eff.fitTo ( dataset ) 
        >>> eff.draw  ( dataset ) 
        """
        
        ddopts  = kwargs.get('data_options',() )
        ddopts  = [  o for o in ddopts ] 
        ddopts.append ( ROOT.RooFit.Efficiency ( self.cut ) )
        ddopts  =  tuple ( ddopts )
        kwargs['data_options'] = ddopts
        
        return self.pdf_draw.draw1 ( dataset , nbins , silent ,
                                     in_range2 , in_range3,
                                     *args , **kwargs )

    # =========================================================================
    ## draw the efficiency (and the dataset)
    #  @code
    #  dataset = ... 
    #  eff     = ... 
    #  eff.fitTo ( dataset ) 
    #  eff.draw2 ( dataset ) 
    #  @endcode 
    def draw2 ( self            ,
                dataset   = None ,
                nbins     = 100  ,
                silent    = True ,
                in_range1 = None ,
                in_range3 = None ,
                *args , **kwargs ) :
        """Draw the efficiency (and the dataset)
        >>> dataset = ... 
        >>> eff     = ... 
        >>> eff.fitTo ( dataset ) 
        >>> eff.draw  ( dataset ) 
        """
        
        ddopts  = kwargs.get('data_options',() )
        ddopts  = [  o for o in ddopts ] 
        ddopts.append ( ROOT.RooFit.Efficiency ( self.cut ) )
        ddopts  =  tuple ( ddopts )
        kwargs['data_options'] = ddopts
        
        return self.pdf_draw.draw2 ( dataset , nbins , silent ,
                                     in_range1 , in_range3,
                                     *args , **kwargs )

    # =========================================================================
    ## draw the efficiency (and the dataset)
    #  @code
    #  dataset = ... 
    #  eff     = ... 
    #  eff.fitTo ( dataset ) 
    #  eff.draw1 ( dataset ) 
    #  @endcode 
    def draw3 ( self            ,
                dataset   = None ,
                nbins     = 100  ,
                silent    = True ,
                in_range1 = None ,
                in_range2 = None ,
                *args , **kwargs ) :
        """Draw the efficiency (and the dataset)
        >>> dataset = ... 
        >>> eff     = ... 
        >>> eff.fitTo ( dataset ) 
        >>> eff.draw  ( dataset ) 
        """
        
        ddopts  = kwargs.get('data_options',() )
        ddopts  = [  o for o in ddopts ] 
        ddopts.append ( ROOT.RooFit.Efficiency ( self.cut ) )
        ddopts  =  tuple ( ddopts )
        kwargs['data_options'] = ddopts
        
        return self.pdf_draw.draw3 ( dataset , nbins , silent ,
                                     in_range1 , in_range2,
                                     *args , **kwargs )

    # =========================================================================
    ## get the efficiency
    #  @code
    #  dataset = ... 
    #  eff     = ... 
    #  eff.fitTo ( dataset ) 
    #  x = 0.15
    #  y = 90
    #  z = 12 
    #  value =  eff(x,y,z) 
    #  @endcode      
    def __call__ (  self , x , y , z , error = False ) :
        """Get the efficiency
        >>> dataset = ... 
        >>> eff     = ... 
        >>> eff.fitTo ( dataset ) 
        >>> x = 0.15
        >>> y = 90 
        >>> z = 15
        >>> value = eff(x,y,z) 
        """
        from ostap.fitting.roofit import SETVAR
        from ostap.math.ve        import VE 
        xx = float ( x ) 
        yy = float ( y ) 
        zz = float ( z ) 
        if xx in self.xvar and yy in self.yvar and zz in self.zvar :
            with SETVAR ( self.xvar ) :
                with SETVAR ( self.yvar ) :
                    with SETVAR ( self.zvar ) :
                        self.xvar.setVal ( xx )
                        self.yvar.setVal ( yy )
                        self.zvar.setVal ( zz )
                        if error and self.fit_result :
                            e = self.eff_fun.getPropagatedError ( self.fit_result )
                            if 0<= e : return  VE ( v ,  e * e )
                        return v 
        logger.error ('Invalid efficiency, return -1 ') 
        return -1 
                        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
