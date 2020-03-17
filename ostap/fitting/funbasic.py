#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/utils.py
#  Set of useful technical utilities to build various fit models 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
# =============================================================================
"""Set of useful technical utilities to build various fit models"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2018-08-14"
__all__     = (
    'FUNC'              , ## base class for    fuuction/PDF-like objects 
    'FUNC2'             , ## base class for 2D-function-like objects
    'FUNC3'             , ## base class for 3D-function-like objects
    'Fun1D'             , ## wrapper for 1D-function
    'Fun2D'             , ## wrapper for 2D-function
    'Fun3D'             , ## wrapper for 3D-function
    )
# =============================================================================
import ROOT, math
from   ostap.core.ostap_types        import list_types , num_types, is_good_number      
from   ostap.core.core               import Ostap , valid_pointer
from   ostap.fitting.variables       import SETVAR
from   ostap.logger.utils            import roo_silent , rootWarning
from   ostap.fitting.roofit          import PDF_fun 
from   ostap.fitting.utils           import MakeVar
import ostap.fitting.variables
import ostap.fitting.roocollections 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.funbasic' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================
## @class FUNC
#  Helper base class for impolementation of various (Roo)Function-wrappers
class FUNC(MakeVar) :
    """Helper base class for implementation of various (Roo)Function-wrappers
    """    
    def __init__ ( self , name , xvar = None ) :
        
        ## name is defined via base class MakeVar 
        self.name  = name ## name is defined via the base class MakeVar 
     
        if   isinstance ( xvar , ROOT.TH1   ) : xvar = xvar.xminmax()
        elif isinstance ( xvar , ROOT.TAxis ) : xvar = xvar.GetXmin() , xvar.GetXmax()

        self.__xvar = None
        
        ## create the variable 
        if isinstance ( xvar , tuple ) and 2 == len ( xvar ) :  
            self.__xvar = self.make_var ( xvar         , ## var 
                                          'x'          , ## name 
                                          'x-variable' , ## title/comment
                                          None         , ## fix ? 
                                          *xvar        ) ## min/max 
        elif isinstance ( xvar , ROOT.RooAbsReal ) :
            self.__xvar = self.make_var ( xvar         , ## var 
                                          'x'          , ## name 
                                          'x-variable' , ## title/comment
                                          fix = None   ) ## fix ? 
        else :
            self.warning('FUNC: ``x-variable''is not specified properly %s/%s' % ( xvar , type ( xvar ) ) )
            self.__xvar = self.make_var( xvar , 'x' , 'x-variable' )
            

        self.__vars = ROOT.RooArgSet  ()
        self.vars.add ( self.__xvar )
        
        self.__config       = {}
        self.__fun          = None
        
        self.__tricks       = True
        self.__fit_result   = None
        
        self.__draw_var     = None
        self.__draw_options = {} ## predefined drawing options for this FUNC/PDF

        self.config = { 'name' : self.name , 'xvar' : self.xvar  }
                
    ## conversion to string 
    def __str__ (  self ) :
        return '%s(%s,xvar=%s)' % ( self.__class__.__name__ , self.name , self.xvar.name )
    __repr__ = __str__ 

    @property 
    def xvar ( self ) :
        """``x''-variable for the fit (same as ``x'')"""
        return self.__xvar
    @property 
    def x    ( self ) :
        """``x''-variable for the fit (same as ``xvar'')"""
        return self.__xvar

    @property 
    def vars ( self ) :
        """``vars'' : variables/observables (as ROOT.RooArgSet)"""
        return self.__vars    

    ## Min/max values for x-variable (when applicable)
    def xminmax ( self ) :
        """Min/max values for x-variable (when applicable)"""
        return self.__xvar.minmax() if self.__xvar else () 

    ## get the proper xmin/xmax range 
    def xmnmx    ( self , xmin , xmax ) :
        """Get the proper xmin/xmax range
        """
        if self.xminmax() :
            
            xmn , xmx = self.xminmax ()
            
            if   is_good_number ( xmin ) : xmin = max ( xmin , xmn )
            else                         : xmin = xmn
            
            if   is_good_number ( xmax ) : xmax = min ( xmax , xmx )
            else                         : xmax = xmx
            
        assert is_good_number ( xmin ),\
               'Invalid type of ``xmin'' %s/%s'  %  ( xmin , type ( xmin ) )
        assert is_good_number ( xmax ),\
               'Invalid type of ``xmin'' %s/%s'  %  ( xmin , type ( xmin ) )

        assert xmin < xmax, 'Invalid xmin/xmax range: %s/%s' % ( xmin , xmax )

        return xmin , xmax
    
    @property
    def fun  ( self ) :
        """The actual function (ROOT.RooAbsReal)"""
        return self.__fun
    @fun.setter
    def fun  ( self , value ) :
        if value is None :
            self.__fun = value
            return        
        assert isinstance ( value , ROOT.RooAbsReal ) , "``pdf'' is not ROOT.RooAbsReal"
        self.__fun = value
        
    @property
    def fun_name ( self ) :
        """``fun_name'' : get the name of the underlying RooAbsReal"""
        return  self.fun.GetName() if self.fun else ''

    @property
    def fit_result ( self ) :
        """``fit_result'' : (the latest) fit resut (TFitResult)"""
        return self.__fit_result
    @fit_result.setter
    def fit_result ( self , value ) :
        assert value is None or isinstance ( value , ROOT.RooFitResult ) , \
               "Invalid value: %s/%s" % ( value , type ( value ) )
        self.__fit_result = None
        if isinstance ( value , ROOT.RooFitResult ) and valid_pointer ( value ) : 
            self.__fit_result = value    
        
    @property
    def config ( self ) :
        """The full configuration info for the PDF"""
        conf = {}
        conf.update ( self.__config )
        return conf
    @config.setter
    def config ( self , value ) :
        conf = {}
        conf.update ( value )
        self.__config = conf

    @property
    def draw_var ( self ) :
        """``draw_var''  :  variable to be drawn if not specified explicitely"""
        return self.__draw_var
    @draw_var.setter
    def draw_var ( self , value ) :
        assert value is None or isinstance ( value , ROOT.RooAbsReal ) , \
               "``draw_var'' has invalid type %s/%s" % ( value , type(value) )
        self.__draw_var = value 

    @property    
    def tricks ( self ) :
        """``tricks'' : flag to allow some  tricks&shortcuts """
        return self.__tricks
    @tricks.setter
    def tricks ( self , value ) :
        val = True if value else False 
        if val and not self.__tricks :
            raise ValueError("Can't allow tricks&shortcuts!")
        self.__tricks = val

    @property
    def draw_options ( self ) :
        """``draw_options'' : dictionary with predefined draw-options for this PDF
        """
        return self.__draw_options

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
        import ostap.plotting.fit_draw as FD

        key_transform = lambda s : s.lower().replace('_','')

        the_key = key_transform ( key ) 

        ## 1. check the explicitely provided arguments
        for k in kwargs :
            if key_transform ( k ) == the_key :
                return kwargs[ k ]
            
        ## check the predefined drawing options for this PDF 
        for k in self.draw_options :
            if key_transform ( k ) == the_key :
                return self.draw_options.get ( k)
            
        ## check the default options
        for k in dir ( FD ) :
            if k.startswith ( '__' ) or k.endswith ( '__' ) : continue 
            if key_transform ( k ) == the_key :
                return getattr ( FD , k ) 

        ##  use the default value 
        return default 

    # ==========================================================================
    ## Add/define new default draw option
    #  @code
    #  pdf = ...
    #  pdf.add_draw_option( 'background_style' ) = Line ( 4 , 2 , 1 ) 
    #  pdf.add_draw_option( 'components_style' ) = Styles ( Line (... ), Area ( ...) , ... ] ) 
    #  pdf.add_draw_option( 'signal_style'     ) = ROOT.RooFit.LineColor ( 2 )
    #  @endcode
    #  @see ostap.plotting.fit_draw
    #  @see ostap.plotting.fit_draw.Style
    #  @see ostap.plotting.fit_draw.Styles
    #  @see ostap.plotting.fit_draw.Styles
    #  @see ostap.plotting.fit_draw.Line
    #  @see ostap.plotting.fit_draw.Area
    def add_draw_option ( self , key , options = () ) :
        """Add/define new default draw option
        - see ostap.plotting.fit_draw
        - see ostap.plotting.fit_draw.Style
        - see ostap.plotting.fit_draw.Styles
        - see ostap.plotting.fit_draw.Styles
        - see ostap.plotting.fit_draw.Line
        - see ostap.plotting.fit_draw.Area
        >>> pdf = ...
        >>> pdf.add_draw_option ( 'data_options'     ) = ROOT.RooFit.MarkerStyle ( 20 ) , ROOT.RooFit.DrawOption  ( 'zp' )
        >>> pdf.add_draw_option ( 'background_style' ) = Line ( 4 , 2 , 1 ) 
        >>> pdf.add_draw_option ( 'components_style' ) = Styles ( Line (... ), Area ( ...) , ... ] 
        >>> pdf.add_draw_option ( 'signal_style'     ) = ROOT.RooFit.LineColor ( 2 )
        """

        key = key.lower() 

        import ostap.plotting.fit_draw as FD
        if not key in FD.keys :
            self.warning("Unknown draw_option '%s'" % key )
            
        option = key.endswith ( '_options' ) 
        style  = key.endswith ( '_style'   )
        if   options :
            if   isinstance ( options , list_types     ) : options = tuple ( options )
            elif isinstance ( options , ROOT.RooCmdArg ) : options = options , 
            else                                         : options = options ,  
        elif style   :
            if   isinstance ( options , FD.Styles      ) : pass
            elif isinstance ( options , FD.Style       ) : options = options , 
            elif isinstance ( options , ROOT.RooCmdArg ) :
                args    = tuple ( 5*[None] + [ options ] )
                options = FD.Styles ( [ FD.Style ( *args ) ] )
            elif isinstance ( options , list_types     ) : options = tuple ( options )
        else :
            self.warning("Neither ``options'' nor ``style''...")

        self.draw_options [ key ] = options 
                    
    # =========================================================================
    ## make a clone for the given function/PDF with the optional
    #  replacement of certain parameters
    #  @code
    #  >>> xpdf = ...
    #  >>> ypdf = xpdf.clone ( xvar = yvar ,  name = 'PDFy' ) 
    #  @endcode 
    def clone ( self , **kwargs ) :
        """Make a clone for the given fuuction/PDF with
        the optional replacement of the certain parameters
        >>> xpdf = ...
        >>> ypdf = xpdf.clone ( xvar = yvar ,  name = 'PDFy' ) 
        """

        ## get config 
        conf = {}
        conf.update ( self.config ) 

        ## modify the name if the name is in config  
        if 'name' in conf :
            name_prefix = kwargs.pop ( 'name_prefix' , '' )
            name_suffix = kwargs.pop ( 'name_suffix' , '' )
            if name_prefix or name_suffix :
                conf['name']  = name_prefix + conf [ 'name' ] + name_suffix
            else :
                conf['name'] += '_copy'
                
        ## update (if needed)
        conf.update ( kwargs )

        KLASS = self.__class__
        cloned = KLASS ( **conf )
        
        return cloned 

    # =========================================================================
    ## make a copy/clone for the given function/PDF 
    #  @code
    #  >>> import copy 
    #  >>> xpdf = ...
    #  >>> ypdf = copy.copy ( xpdf ) 
    #  @endcode 
    def __copy__ ( self ) :
        """Make a copy/clone for the given function/PDF 
        >>> import copy 
        >>> xpdf = ...
        >>> ypdf = copy.copy ( xpdf ) 
        """
        return self.clone()

    # ========================================================================
    # some generic stuff 
    # ========================================================================
    ## helper  function to implement some math stuff 
    def _get_stat_ ( self , funcall , *args , **kwargs ) :
        """Helper  function to implement some math stuff 
        """
        fun         = self.fun
        xmin , xmax = self.xminmax()

        xmin = kwargs.pop ( 'xmin' , xmin )
        xmax = kwargs.pop ( 'xmax' , xmax )
        
        if self.tricks and hasattr ( fun , 'function' ) :    
            ff = fun.function()
            if   hasattr  ( fun , 'setPars'   ) : fun.setPars()             
        else :
            ff = PDF_fun  ( fun , self.xvar , xmin , xmax )
            
        return funcall ( ff , xmin , xmax , *args , **kwargs )

    # ========================================================================
    ## get the effective Full Width at Half Maximum
    def fwhm ( self , **kwargs ) :
        """Get the effective Full Width at  Half Maximum
        >>>  fun = ...
        >>>  print 'FWHM: %s ' % fun.fwhm()
        """
        ## use generic machinery 
        from ostap.stats.moments import width as _width
        w = self._get_stat_ ( _width , **kwargs )
        return  w [ 1 ] - w [ 0 ]

    # ========================================================================
    ## Get the effective midpoint/location:
    #  a mid point of an interval \f$ [x_{low}, x_{high}]\f$,
    #  \f$ x_{mid} = \frac{ x_{low} + x_{high}}{2}\f$, where  
    #  where \f$ f(x_{low}) = f(x_{high}) = \frac{f_{max}(x)}{2} \f$ 
    # (the same points are used for FWHM)
    def mid_point ( self , **kwargs ) :
        """Get the effective midpoint/location:
        - a mid point of an interval  x_low, x_high
        x_mid = (x_low + x_high)/2 {2},
        where  f(x_low) = f(x_high) = f_{max}(x)
        - the same points are used for FWHM
        """
        ## use generic machinery 
        from ostap.stats.moments import width as _width
        w = self._get_stat_ ( _width , **kwargs )
        return 0.5 * ( w [ 1 ] + w [ 0 ] )

    # =========================================================================
    ## get the effective mode 
    def mode ( self , **kwargs ) :
        """Get the effective mode
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'MODE: %s ' % pdf.mode()
        """
        from ostap.stats.moments import mode as _mode
        return self._get_stat_ ( _mode , **kwargs )

    # =========================================================================
    ## simple 'function-like' interface 
    def __call__ ( self , x , error = False , normalized = False ) :
        """ Function as a ``function''
        >>> fun  = ...
        >>> x = 1
        >>> y = fun ( x ) 
        """
        if error and not normalized :
            self.error("Can't get error for non-normalized call" )
            error = False

        if isinstance ( self.xvar , ROOT.RooRealVar ) :
            
            mn , mx = self.xminmax()            
            if not mn <= x <= mx : return 0             ## RETUR N
                
            with SETVAR( self.xvar ) :
                
                self.xvar.setVal ( x )
                
                v = self.fun.getVal ( self.vars ) if normalized else self.fun.getVal ()  

                e = -1 
                if   error and isinstance ( error , ROOT.TFitResult ) : 
                    e = self.fun.getPropagatedError ( error           )
                elif error and self.fit_result : 
                    e = self.fun.getPropagatedError ( self.fit_result )
                    
                if 0 <= e : return  VE ( v ,  e * e )
                
                return v
                
        raise AttributeError('Something wrong goes here')

    # ========================================================================
    ## convert to float 
    def __float__ ( self ) :
        """Convert to float
        >>> pdf = ...
        >>> v = float ( pdf )
        """
        return self.fun.getVal ( self.vars ) 

    # =========================================================================
    ## get the derivative at  point x 
    def derivative ( self , x ) :
        """Get derivative at point x 
        >>> pdf = ...
        >>> print pdf.derivative ( 0 ) 
        """
        ## check limits 
        if hasattr ( self.xvar , 'getMin' ) and x < self.xvar.getMin() : return 0.
        if hasattr ( self.xvar , 'getMax' ) and x > self.xvar.getMax() : return 0.

        ## make a try to use analytical derivatives 
        if self.tricks  and hasattr ( self , 'pdf' ) :
            _fun = self.fun 
            if hasattr ( _fun , 'setPars'  ) : _fun.setPars() 
            try: 
                if hasattr ( _fun , 'function' ) :
                    _ff = _fun.function() 
                    if hasattr ( _ff , 'derivative' ) :
                        return _ff.derivative ( x )
            except:
                pass
            
        ## use numerical derivatives 
        from ostap.math.derivative import derivative as _derivatve
        return _derivative ( self , x )

    # ==========================================================================
    ## get a minimum of PDF for certain interval
    #  @code
    #  pdf = ...
    #  x   = pdf.minimum() 
    #  @endcode 
    def minimum ( self , xmin = None , xmax = None , x0 = None ) :
        """Get a minimum of PDF for certain interval
        >>> pdf = ...
        >>> x = pdf.minimum()
        """
        if xmin is None : xmin = self.xminmax()[0]
        if xmax is None : xmax = self.xminmax()[1]
        if self.xminmax() :
            xmin =  max ( xmin , self.xminmax()[0] )
            xmax =  min ( xmax , self.xminmax()[1] )
            
        if x0 is None           : x0 = 0.5 * ( xmin + xmax )
        
        if not xmin <= x0 <= xmax :
            logger.error("Wrong xmin/x0/xmax: %s/%s/%s"   % ( xmin , x0 , xmax ) )
        
        from ostap.math.minimize import sp_minimum_1D
        return sp_minimum_1D (  self , xmin , xmax , x0 )

    # ==========================================================================
    ## get a maximum of PDF for certain interval
    #  @code
    #  pdf = ...
    #  x   = pdf.maximum() 
    #  @endcode 
    def maximum ( self , xmin = None , xmax = None , x0 = None ) :
        """Get a maximum of PDF for certain interval
        >>> pdf = ...
        >>> x = pdf.maximum()
        """
        if xmin is None : xmin = self.xminmax()[0]
        if xmax is None : xmax = self.xminmax()[1]
        if self.xminmax() :
            xmin =  max ( xmin , self.xminmax()[0] )
            xmax  = min ( xmax , self.xminmax()[1] )
            
        if x0 is None           : x0 = 0.5 * ( xmin + xmax )

        if not xmin <= x0 <= xmax :
            logger.error("Wrong xmin/x0/xmax: %s/%s/%s"   % ( xmin , x0 , xmax ) )
        
        from ostap.math.minimize import sp_maximum_1D
        return sp_maximum_1D (  self , xmin , xmax , x0 )


    # ================================================================================
    ## visualise the function 
    #  @code
    #  fun.draw ) 
    #  @endcode
    #  @see ostap.plotting.fit_draw
    #
    #  Drawing options can be specified as keyword arguments:
    #  @code
    #  curve = ROOT.RooFit.LineColor ( ROOT.kRed ) , ROOT.RooFit.LineWidth ( 3 )
    #  f = fun.draw ( ... , curve_options = curve  , )
    #  @endcode
    #  When options are not provided explicitly, the options defined in the PDF are looked for:
    #  @code
    #  curve = ROOT.RooFit.LineColor ( ROOT.kRed ) , ROOT.RooFit.LineWidth ( 3 )
    #  fun.draw_opptions['curve_options'] = curve 
    #  f = pdf.draw ( ...)
    #  @endcode
    #  Otherwise the default options,  defined in ostap.plotting.fit_draw module, are used 
    #  @see ostap.plotting.fit_draw
    def draw ( self                         ,
               drawvar               = None ,
               silent                = True ,   ## silent mode ?
               style                 = None ,   ## use another style ?
               args                  = ()   , 
               **kwargs                     ) :
        """  Visualize the function
        >>> fun.draw ()

        - Drawing options can be specified as keyword arguments:
        
        >>> curve = ROOT.RooFit.LineColor ( ROOT.kRed ) , ROOT.RooFit.LineWidth ( 3 )
        >>> f = pdf.draw ( ... , curve_options = curve  , ... )
        
        - when options are not provided explicitly, the options defined in the PDF are looked for:
        
        >>> curve = ROOT.RooFit.LineColor ( ROOT.kRed ) , ROOT.RooFit.LineWidth ( 3 )
        >>> pdf.draw_options['curve_options'] = curve 
        >>> f = fun.draw ( ...)
        
        - otherwise the default options, defined in ostap.plotting.fit_draw module, are used 

        """
        #        
        from ostap.plotting.style import useStyle 

        #
        ## again the context
        # 
        with roo_silent ( silent ) , useStyle ( style ) :

            drawvar = drawvar if drawvar else ( self.draw_var if self.draw_var else self.xvar )  

            frame = drawvar.frame ()
            
            #
            ## the total fit curve
            #
            coptions   = self.draw_option ( 'curve_options' , **kwargs )
            self.fun .plotOn ( frame ) ##  , *coptions )
            kwargs.pop ( 'curve_options' , () )            
            #
            #
            ## suppress ugly axis labels
            #
            if not kwargs.get ( 'draw_axis_title' , False ) :  
                frame.SetXTitle ( '' )
                frame.SetYTitle ( '' )
                frame.SetZTitle ( '' )
                
            #
            ## Draw the frame!
            #
            if not ROOT.gROOT.IsBatch() :
                with rootWarning (): frame.draw ( kwargs.pop ( 'draw_options','' ) )
            
            if kwargs :
                self.warning("draw: ignored unknown options: %s" % list( kwargs.keys() ) ) 

            return frame

# =============================================================================
## @class Fun1D
#  Simple wrapper for 1D-function
#  @code
#  func = ...
#  xvar = ...
#  f1d  = Fun1D ( func , xvar = xvar ) 
#  @endcode 
class Fun1D ( FUNC ) :
    """Simple wrapper for 1D-function
    >>> func = ...
    >>> xvar = ...
    >>> f1d  = Fun1D ( func , xvar = xvar ) 
    """
    def __init__ ( self ,  fun , xvar , name = '' ) :

        if isinstance ( fun , FUNC ) :
            self.__argfun = fun 
            fun           = fun.fun
        
        assert xvar and isinstance ( xvar , ROOT.RooAbsReal ) , "``xvar'' must be ROOT.RooAbsReal"
        assert fun  and isinstance ( fun  , ROOT.RooAbsReal ) , "``fun''  must be ROOT.RooAbsReal"

        if fun is xvar : fun = Ostap.MoreRooFit.Id ( "" , "" , xvar )            
            
        if not name : name = 'Fun1D_%s' % fun.GetName() 

        FUNC.__init__ ( self , name , xvar = xvar )

        if not self.xvar in fun.getParameters ( 0 ) and not self.xvar is fun : 
            self.warning ("Function does not depends on xvar=%s" % self.xvar.name )
            
        self.fun = fun
        
        self.config = {
            'fun'  : self.fun  ,
            'xvar' : self.xvar ,
            'name' : self.name ,            
            }

    # =========================================================================
    ## redefine the clone method, allowing only the name to be changed
    #  @attention redefinition of parameters and variables is disabled,
    #             since it can't be done in a safe way                  
    def clone ( self , fun = None , xvar = None , **kwargs ) :
        """Redefine the clone method, allowing only the name to be changed
         - redefinition of parameters and variables is disabled,
         since it can't be done in a safe way          
        """
        if fun  and not  fun is self.fun  :
            raise AttributeError("Fun1D cannot be cloned with different `fun''" )
        if xvar and not xvar is self.xvar :
            raise AttributeError("Fun1D cannot be cloned with different ``xvar''")

        return FUNC.clone ( self , **kwargs ) 

    
# =============================================================================
## @class FUNC2
#  The base class for 2D-function
class FUNC2(FUNC) :
    """Base class for 2D-function
    """
    def __init__ ( self , name , xvar = None , yvar = None ) :

        FUNC.__init__ ( self , name , xvar = xvar )
        
        if   isinstance ( yvar , ROOT.TH1   ) : yvar = yvar.xminmax()
        elif isinstance ( yvar , ROOT.TAxis ) : yvar = yvar.GetXmin() , yvar.GetXmax()

        self.__yvar = None
        
        ## create the variable 
        if isinstance ( yvar , tuple ) and 2 == len ( yvar ) :  
            self.__yvar = self.make_var ( yvar         , ## var 
                                          'y'          , ## name 
                                          'y-variable' , ## title/comment
                                          None         , ## fix ? 
                                          *yvar        ) ## min/max 
        elif isinstance ( yvar , ROOT.RooAbsReal ) :
            self.__yvar = self.make_var ( yvar         , ## var 
                                          'y'          , ## name 
                                          'y-variable' , ## title/comment
                                          fix = None   ) ## fix ? 
        else :
            self.warning ( 'FUNC: ``y-variable''is not specified properly %s/%s' % ( yvar , type ( yvar ) ) )
            self.__yvar = self.make_var( yvar , 'y' , 'y-variable' )
            

        self.vars.add ( self.__yvar )
        
        ## save the configuration
        self.config = {
            'name' : self.name ,
            'xvar' : self.xvar ,
            'yvar' : self.yvar ,            
            }

    ## conversion to string 
    def __str__ (  self ) :
        return '%s(%s,xvar=%s,yvar=%s)' % ( self.__class__.__name__ ,
                                            self.name      ,
                                            self.xvar.name ,
                                            self.yvar.name )
    __repr__ = __str__ 


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
    ## simple 'function-like' interface 
    def __call__ ( self , x , y , error = False , normalized = False ) :
        """ Function as a ``function''
        >>> fun  = ...
        >>> x = 1
        >>> y = 2 
        >>> v = fun ( x , y ) 
        """
        if error and not normalized :
            self.error("Can't get error for non-normalized call" )
            error = False
            
        if isinstance ( self.xvar , ROOT.RooRealVar ) and \
           isinstance ( self.yvar , ROOT.RooRealVar )     :
            
            mn , mx = self.xminmax()            
            if not mn <= x <= mx : return 0             ## RETUR N

            mn , mx = self.yminmax()            
            if not mn <= y <= mx : return 0             ## RETUR N

                
            with SETVAR( self.xvar ), SETVAR ( self.yvar )  :
                
                self.xvar.setVal ( x )
                self.yvar.setVal ( y )
                
                v = self.fun.getVal ( self.vars ) if normalized else self.fun.getVal ()  

                e = -1 
                if   error and isinstance ( error , ROOT.TFitResult ) : 
                    e = self.fun.getPropagatedError ( error           )
                elif error and self.fit_result : 
                    e = self.fun.getPropagatedError ( self.fit_result )
                    
                if 0 <= e : return  VE ( v ,  e * e )
                
                return v
                
        raise AttributeError('Something wrong goes here')


    # =========================================================================
    ## make 1D-plot
    def draw ( self                         ,
               drawvar               = None ,
               silent                = True ,
               in_range              = None ,
               args                  = ()   , 
               **kwargs                     ) : 
        """
        Make 1D-plot:
        """
        if   drawvar in ( 'x'  , 'X' , '1' , 1 , self.xvar.name ) : drawvar = self.xvar
        elif drawvar in ( 'y'  , 'Y' , '2' , 2 , self.yvar.name ) : drawvar = self.yvar

        newargs = kwargs.copy ()
        
        if in_range and isinstance ( in_range , list_types ) and 2 == len ( in_range ) :
            low  = in_range [ 0 ]
            high = in_range [ 1 ]
            if isinstance ( low , num_types ) and isinstance ( high , num_types ) and low < high :
                range_name = 'aux2_range_%s' % self.name 
                with rooSilent ( 3 ) : drawvar.setRange ( range_name , low , high )
                in_range = range_name
    
        if in_range and not isinstance ( in_range , list_types ) :
            in_range = in_range ,
            
        if in_range : 
            options_project = tuple ( [  ROOT.RooFit.ProjectionRange ( i ) for i in in_range ] )
            for key in  ( 'curve_options'               , ) : 
                newargs [ key ] =  self.draw_option ( key , **newargs ) + options_project
                
        ## redefine the drawing variable:
        self.draw_var = drawvar
        
        ## delegate the actual drawing to the base class
        result = FUNC.draw ( self            ,
                             silent = silent ,
                             args   = args   , **newargs )
        
        self.draw_var = None
        return result 

    # =========================================================================
    ## draw function 
    #  @code
    #  fx  = fun.draw1 () 
    #  f1  = fun.draw1 (in_range = (2,3) ) 
    #  model.yvar.setRange ( 'QUQU2' , 2 , 3 ) 
    #  f2  = fun.draw1 ( in_range = 'QUQU2')
    #  @endcode 
    def draw1 ( self            ,
                silent   = True ,
                in_range = None ,
                args     = ()   , **kwargs ) :
        """ Draw the function 
        
        >>> fx  = fun.draw1 ()
        >>> f1  = fun.draw1 ( in_range = (2,3) ) 
        >>> f1  = fun.draw1 ( in_range = 'QUQU2')
        
        """
        if in_range and isinstance ( in_range , tuple ) and 2 == len ( in_range ) :
            range_name = 'aux2_rng2_%s' % self.name 
            with rooSilent ( 3 ) : self.yvar.setRange ( range_name , in_range[0] , in_range[1] )
            in_range = range_name 

        return self.draw ( drawvar  = self.xvar , 
                           silent   = silent    ,
                           in_range = in_range  ,
                           args     = args      , **kwargs )

    # =========================================================================
    ## draw function 
    #  @code
    #  fx  = fun.draw2 () 
    #  f1  = fun.draw2 (in_range = (2,3) ) 
    #  fun.xvar.setRange ( 'QUQU2' , 2 , 3 ) 
    #  f2  = fun.draw2 ( in_range = 'QUQU2')
    #  @endcode 
    def draw2 ( self            ,
                silent   = True ,
                in_range = None ,
                args     = ()   , **kwargs ) :
        """ Draw the function 
        
        >>> fx  = fun.draw2 () ## draw results
        >>> f1  = fun.draw2 (in_range = (2,3) ) 
        >>> fun.xvar.setRange ( 'QUQU2' , 2 , 3 ) 
        >>> f2  = fun.draw2 ( in_range = 'QUQU2')
        """
        if in_range and isinstance ( in_range , tuple ) and 2 == len ( in_range ) :
            range_name = 'aux2_rng2_%s' % self.name 
            with rooSilent ( 3 ) : self.xvar.setRange ( range_name , in_range[0] , in_range[1] )
            in_range = range_name 

        return self.draw ( drawvar  = self.yvar , 
                           silent   = silent    ,
                           in_range = in_range  ,
                           args     = args      , **kwargs )
    
# =============================================================================
## @class Fun2D
#  Simple wrapper for 2D-function
#  @code
#  func = ...
#  xvar = ...
#  yvar = ...
#  f2d  = Fun2D ( func , xvar = xvar , yvar = yvar ) 
#  @endcode 
class Fun2D ( FUNC2 ) :
    """Simple wrapper for 2D-function
    >>> func = ...
    >>> xvar = ...
    >>> yvar = ...
    >>> f2d  = Fun2D ( func , xvar = xvar , yvar = yvar ) 
    """
    def __init__ ( self ,  fun , xvar , yvar , name = '' ) :

        if isinstance ( fun , FUNC ) :
            self.__argfun = fun 
            fun           = fun.fun
        
        assert xvar and isinstance ( xvar , ROOT.RooAbsReal ) , "``xvar'' must be ROOT.RooAbsReal"
        assert yvar and isinstance ( yvar , ROOT.RooAbsReal ) , "``yvar'' must be ROOT.RooAbsReal"
        assert fun  and isinstance ( fun  , ROOT.RooAbsReal ) , "``fun''  must be ROOT.RooAbsReal"

        assert not xvar is yvar, "xvar and yvar must be different!"

        if fun is xvar : fun = Ostap.MoreRooFit.Addition  ( "" , "" , xvar )
        if fun is yvar : fun = Ostap.MoreRooFit.Addition  ( "" , "" , yvar )
            
        if not name : name = 'Fun2D_%s' % fun.GetName() 

        FUNC2.__init__ ( self , name , xvar = xvar , yvar = yvar )

        if not self.xvar in fun.getParameters ( 0 ) and not self.xvar is fun :
            self.warning ("Function does not depends on xvar=%s" % self.xvar.name ) 
        if not self.yvar in fun.getParameters ( 0 ) and not self.yvar is fun :
            self.warning ("Function does not depends on yvar=%s" % self.yvar.name ) 

        self.fun = fun
        
        self.config = {
            'fun'  : self.fun  ,
            'xvar' : self.xvar ,
            'yvar' : self.yvar ,
            'name' : self.name ,            
            }

    # =========================================================================
    ## redefine the clone method, allowing only the name to be changed
    #  @attention redefinition of parameters and variables is disabled,
    #             since it can't be done in a safe way                  
    def clone ( self , fun = None , xvar = None , yvar = None , **kwargs ) :
        """Redefine the clone method, allowing only the name to be changed
         - redefinition of parameters and variables is disabled,
         since it can't be done in a safe way          
        """
        if fun  and not  fun is self.fun  :
            raise AttributeError("Fun2D cannot be cloned with different `fun''" )
        if xvar and not xvar is self.xvar :
            raise AttributeError("Fun2D cannot be cloned with different ``xvar''")
        if yvar and not yvar is self.yvar :
            raise AttributeError("Fun2D cannot be cloned with different ``yvar''")

        return FUNC2.clone ( self , **kwargs ) 

    
# =============================================================================
## @class FUNC3
#  The base class for 3D-function
class FUNC3(FUNC2) :
    """Base class for 3D-function
    """
    def __init__ ( self , name , xvar = None , yvar = None , zvar = None ) :

        FUNC2.__init__ ( self , name , xvar = xvar , yvar = yvar )
        
        if   isinstance ( zvar , ROOT.TH1   ) : zvar = zvar.xminmax()
        elif isinstance ( zvar , ROOT.TAxis ) : zvar = zvar.GetXmin() , zvar.GetXmax()

        self.__zvar = None
        
        ## create the variable 
        if isinstance ( zvar , tuple ) and 2 == len ( zvar ) :  
            self.__zvar = self.make_var ( zvar         , ## var 
                                          'z'          , ## name 
                                          'z-variable' , ## title/comment
                                          None         , ## fix ? 
                                          *zvar        ) ## min/max 
        elif isinstance ( zvar , ROOT.RooAbsReal ) :
            self.__zvar = self.make_var ( zvar         , ## var 
                                          'z'          , ## name 
                                          'z-variable' , ## title/comment
                                          fix = None   ) ## fix ? 
        else :
            self.warning ( 'FUNC: ``z-variable''is not specified properly %s/%s' % ( zvar , type ( zvar ) ) )
            self.__zvar = self.make_var( zvar , 'z' , 'z-variable' )
            

        self.vars.add ( self.__zvar )
        
        ## save the configuration
        self.config = {
            'name' : self.name ,
            'xvar' : self.xvar ,
            'yvar' : self.yvar ,            
            'zvar' : self.zvar ,            
            }

    ## conversion to string 
    def __str__ (  self ) :
        return '%s(%s,xvar=%s,yvar=%s,zvar=%s)' % ( self.__class__.__name__ ,
                                                    self.name      ,
                                                    self.xvar.name ,
                                                    self.yvar.name ,
                                                    self.zvar.name )
    __repr__ = __str__ 

    def zminmax ( self ) :
        """Min/max values for z-varibale"""
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
    ## simple 'function-like' interface 
    def __call__ ( self , x , y , z , error = False , normalized = False ) :
        """ Function as a ``function''
        >>> fun  = ...
        >>> x = 1
        >>> y = 2 
        >>> z = 0 
        >>> v = fun ( x , y , z ) 
        """
        if error and not normalized :
            self.error("Can't get error for non-normalized call" )
            error = False
            
        if isinstance ( self.xvar , ROOT.RooRealVar ) and \
           isinstance ( self.yvar , ROOT.RooRealVar ) and \
           isinstance ( self.zvar , ROOT.RooRealVar )     : 
           
            mn , mx = self.xminmax()            
            if not mn <= x <= mx : return 0             ## RETUR N

            mn , mx = self.yminmax()            
            if not mn <= y <= mx : return 0             ## RETUR N

            mn , mx = self.zminmax()            
            if not mn <= z <= mx : return 0             ## RETUR N

                
            with SETVAR( self.xvar ), SETVAR ( self.yvar ) , SETVAR ( self.zvar ) :
                
                self.xvar.setVal ( x )
                self.yvar.setVal ( y )
                self.zvar.setVal ( z )
                
                v = self.fun.getVal ( self.vars ) if normalized else self.fun.getVal ()  

                e = -1 
                if   error and isinstance ( error , ROOT.TFitResult ) : 
                    e = self.fun.getPropagatedError ( error           )
                elif error and self.fit_result : 
                    e = self.fun.getPropagatedError ( self.fit_result )
                    
                if 0 <= e : return  VE ( v ,  e * e )
                
                return v
                
        raise AttributeError('Something wrong goes here')


    # =========================================================================
    ## draw the 1st variable
    #  @code
    #  fx  = fun.draw1 ( dataset ) 
    #  fx  = fun.draw1 ( in_range2 = (2,3) ) 
    #  fun.yvar.setRange ( 'QUQU2' , 2 , 3 ) 
    #  fx  = fun.draw1 ( in_range2 = 'QUQU2')
    #  @endcode 
    def draw1 ( self             ,
                silent    = True ,
                in_range2 = None ,
                in_range3 = None ,
                args      = ()   , **kwargs ) :
        """ Draw the 1st variable        
        >>> fx  = fun.draw1 ( dataset ) 
        >>> fx  = fun.draw1 ( in_range2 = (2,3) ) 
        >>> fun.yvar.setRange ( 'QUQU2' , 2 , 3 ) 
        >>> fx  = fun.draw1 ( in_range2 = 'QUQU2') 
        """
        if in_range2 and isinstance ( in_range2 , tuple ) and 2 == len ( in_range2 ) :
            range_name = 'aux3_rng12_%s' % self.name
            with rooSilent ( 3 ) : self.yvar.setRange ( range_name , in_range2[0] , in_range2[1] )
            in_range2  = range_name 

        if in_range3 and isinstance ( in_range3 , tuple ) and 2 == len ( in_range3 ) :
            range_name = 'aux3_rng13_%s' % self.name
            with rooSilent ( 3 ) : self.zvar.setRange ( range_name , in_range3[0] , in_range3[1] )
            in_range3  = range_name 

        in_range = []
        if in_range2 : in_range.append ( in_range2 )
        if in_range3 : in_range.append ( in_range3 )
        in_range = tuple ( in_range )
        
        return self.draw ( drawvar  = self.xvar , 
                           silent   = silent    ,
                           in_range = in_range  ,
                           args     = args      , **kwargs )


    # =========================================================================
    ## draw the 2nd variable
    #  @code
    #  fy  = fun.draw2 ( dataset ) 
    #  fy  = fun.draw2 ( in_range3 = (2,3) ) 
    #  fun.zvar.setRange ( 'QUQU2' , 2 , 3 ) 
    #  fy  = fun.draw1 ( in_range3 = 'QUQU2')  
    #  @endcode 
    def draw2 ( self             ,
                silent    = True ,
                in_range1 = None ,
                in_range3 = None ,
                args      = ()   , **kwargs ) :
        """Draw the 2nd variable
        >>> fy  = fun.draw2 ( dataset ) 
        >>> fy  = fun.draw2 ( in_range3 = (2,3) ) 
        >>> fun.zvar.setRange ( 'QUQU2' , 2 , 3 ) 
        >>> fy  = fun.draw1 ( in_range3 = 'QUQU2')  
        """
        
        if in_range1 and isinstance ( in_range1 , tuple ) and 2 == len ( in_range1 ) :
            range_name = 'aux3_rng21_%s' % self.name
            with rooSilent ( 3 ) : self.xvar.setRange ( range_name , in_range1[0] , in_range1[1] )
            in_range1  = range_name

        if in_range3 and isinstance ( in_range3 , tuple ) and 2 == len ( in_range3 ) :
            range_name = 'aux3_rng23_%s' % self.name
            with rooSilent ( 3 ) : self.zvar.setRange ( range_name , in_range3[0] , in_range3[1] )
            in_range3  = range_name 

        in_range = []
        if in_range1 : in_range.append ( in_range1 )
        if in_range3 : in_range.append ( in_range3 )
        in_range = tuple ( in_range )
        
        return self.draw ( drawvar  = self.yvar , 
                           silent   = silent    ,
                           in_range = in_range  ,
                           args     = args      , **kwargs )

    # =========================================================================
    ## draw the 3rd variable
    #  @code
    #  fz  = fun.draw3 ( dataset ) 
    #  fz  = fun.draw3 ( in_range1 = (2,3) ) 
    #  fun.xvar.setRange ( 'QUQU2' , 2 , 3 ) 
    #  fz  = fun.draw1 ( in_range1 = 'QUQU2') 
    #  @endcode 
    def draw3 ( self             ,
                silent    = True ,
                in_range1 = None ,
                in_range2 = None ,
                args      = ()   , **kwargs ) :
        """Draw the 3rd variable
        >>> fz  = fun.draw3 ( dataset ) 
        >>> fz  = fun.draw3 ( in_range1 = (2,3) ) 
        >>> fun.xvar.setRange ( 'QUQU2' , 2 , 3 ) 
        >>> fz  = fun.draw1 ( in_range1 = 'QUQU2') 
        """
        
        if in_range1 and isinstance ( in_range1 , tuple ) and 2 == len ( in_range1 ) :
            range_name = 'aux3_rng31_%s' % self.name
            with rooSilent ( 3 ) : self.xvar.setRange ( range_name ,  in_range1[0] , in_range1[1] )
            in_range1  = range_name 
            
        if in_range2 and isinstance ( in_range2 , tuple ) and 2 == len ( in_range2 ) :
            range_name = 'aux3_rng32_%s' % self.name
            with rooSilent ( 3 ) : self.yvar.setRange ( range_name , in_range2[0] , in_range2[1] )
            in_range2  = range_name 

        in_range = []
        if in_range1 : in_range.append ( in_range1 )
        if in_range2 : in_range.append ( in_range2 )
        in_range = tuple ( in_range )
        
        return self.draw ( drawvar  = self.zvar , 
                           silent   = silent    ,
                           in_range = in_range  ,
                           args     = args      , **kwargs )
    
    # =========================================================================
    ## make 1D-plot
    def draw ( self            ,
               drawvar  = None ,
               silent   = True ,
               in_range = None ,
               args     = ()   , 
               **kwargs        ) : 
        """
        Make 1D-plot:
        """
        if drawvar in ( 'z'  , 'Z' , '3' , 3 , self.zvar.name ) :
            drawvar = self.zvar
            
        return FUNC2.draw  ( self ,
                             drawvar  = drawvar  ,
                             silent   =  silent  ,
                             in_range = in_range ,
                             args     = args     , **kwargs )


    
# =============================================================================
## @class Fun3D
#  Simple wrapper for 3D-function
#  @code
#  func = ...
#  xvar = ...
#  yvar = ...
#  zvar = ...
#  f3d  = Fun3D ( func , xvar = xvar , yvar = yvar , zvar = zvar ) 
#  @endcode 
class Fun3D ( FUNC3 ) :
    """Simple wrapper for 3D-function
    >>> func = ...
    >>> xvar = ...
    >>> yvar = ...
    >>> zvar = ...
    >>> f3d  = Fun3D ( func , xvar = xvar , yvar = yvar , zvar = zvar ) 
    """
    def __init__ ( self ,  fun , xvar , yvar , zvar , name = '' ) :

        if isinstance ( fun , FUNC ) :
            self.__argfun = fun 
            fun           = fun.fun

        assert xvar and isinstance ( xvar , ROOT.RooAbsReal ) , "``xvar'' must be ROOT.RooAbsReal"
        assert yvar and isinstance ( yvar , ROOT.RooAbsReal ) , "``yvar'' must be ROOT.RooAbsReal"
        assert zvar and isinstance ( zvar , ROOT.RooAbsReal ) , "``zvar'' must be ROOT.RooAbsReal"
        assert fun  and isinstance ( fun  , ROOT.RooAbsReal ) , "``fun''  must be ROOT.RooAbsReal"

        assert not xvar is yvar, "xvar and yvar must be different!"
        assert not xvar is zvar, "xvar and zvar must be different!"
        assert not yvar is zvar, "yvar and zvar must be different!"
        
        if fun is xvar : fun = Ostap.MoreRooFit.Id ( xvar )
        if fun is yvar : fun = Ostap.MoreRooFit.Id ( yvar )
        if fun is zvar : fun = Ostap.MoreRooFit.Id ( zvar )
            
        if not name : name = 'Fun3D_%s' % fun.GetName() 

        FUNC3.__init__ ( self , name , xvar = xvar , yvar = yvar , zvar = zvar )

        if not self.xvar in fun.getParameters ( 0 ) and not self.xvar is fun :
            self.warning ("Function does not depends on xvar=%s" % self.xvar.name ) 
        if not self.yvar in fun.getParameters ( 0 ) and not self.yvar is fun :
            self.warning ("Function does not depends on yvar=%s" % self.yvar.name ) 
        if not self.zvar in fun.getParameters ( 0 ) and not self.zvar is fun :
            self.warning ("Function does not depends on zvar=%s" % self.zvar.name ) 


        self.fun = fun
        
        self.config = {
            'fun'  : self.fun  ,
            'xvar' : self.xvar ,
            'yvar' : self.yvar ,
            'zvar' : self.zvar ,
            'name' : self.name ,            
            }

    # =========================================================================
    ## redefine the clone method, allowing only the name to be changed
    #  @attention redefinition of parameters and variables is disabled,
    #             since it can't be done in a safe way                  
    def clone ( self , fun = None , xvar = None , yvar = None , zvar = None , **kwargs ) :
        """Redefine the clone method, allowing only the name to be changed
         - redefinition of parameters and variables is disabled,
         since it can't be done in a safe way          
        """
        if fun  and not  fun is self.fun  :
            raise AttributeError("Fun3D cannot be cloned with different `fun''" )
        if xvar and not xvar is self.xvar :
            raise AttributeError("Fun3D cannot be cloned with different ``xvar''")
        if yvar and not yvar is self.yvar :
            raise AttributeError("Fun3D cannot be cloned with different ``yvar''")
        if zvar and not zvar is self.zvar :
            raise AttributeError("Fun3D cannot be cloned with different ``zvar''")

        return FUNC3.clone ( self , **kwargs ) 


# =============================================================================
## Operator for `3D-function (op) other`:
def _f3_op_ ( fun1 , fun2 , ftype , fname ) :
    """ Operator for `3D-function (op) other`:
    """
    
    xvar, yvar, zvar = fun1.xvar , fun1.yvar , fun1.zvar
    
    if   isinstance ( fun2 , FUNC3 ) :
        
        if  fun2.xvar in fun1.vars and fun2.yvar in fun1.vars and fun2.zvar in fun1.vars :
            ## sum 
            s      = ftype ( fun1.fun , fun2.fun )
            result = Fun3D ( s ,
                             xvar = xvar ,
                             yvar = yvar ,
                             zvar = zvar ,
                             name = fname % ( fun1.name , fun2.name ) )
            
            result.aux_keep.append ( fun1 )
            result.aux_keep.append ( fun2 )
            return result 
            
    elif isinstance ( fun2 , FUNC2 ) :
        
        if fun2.xvar in fun1.vars and fun2.yvar in fun1.vars :
            ## sum 
            s      = ftype ( self.fun , other.fun )
            result = Fun3D ( s ,
                             xvar = xvar ,
                             yvar = yvar ,
                             zvar = zvar ,
                             name = fname % ( fun1.name , fun2.name ) )
            
            result.aux_keep.append ( fun1 )
            result.aux_keep.append ( fun2 )
            return result 

    elif isinstance ( fun2 , FUNC  ) :
        
        if fun2.xvar in fun1.vars :
            ## sum 
            s      = ftype ( fun1.fun , fun2.fun )
            result = Fun3D ( s ,
                             xvar = xvar ,
                             yvar = yvar ,
                             zvar = zvar ,
                             name = fname % ( fun1.name , fun2.name ) )
            
            result.aux_keep.append ( fun1 )
            result.aux_keep.append ( fun2 )
            return result
        
    elif isinstance ( fun2 , num_types ) :
        
        fun2   = ROOT.RooRealConstant.value ( float ( fun2 )  ) 
        s      = ftype ( fun1.fun , fun2  )
        result = Fun3D ( s ,
                         xvar = xvar ,
                         yvar = yvar ,
                         zvar = zvar ,
                         name = fname % ( fun1.name , fun2.name  ) )
        
        result.aux_keep.append ( fun1 )
        result.aux_keep.append ( fun2 )
        return result
    
    elif isinstance ( fun2 , ROOT.RooAbsReal ) :
        
        s      = ftype ( fun1.fun , fun2 )
        result = Fun3D ( s ,
                         xvar = xvar  ,
                         yvar = yvar  ,
                         zvar = zvar  ,
                         name = fname % ( fun1.name , fun2.name ) )
        
        result.aux_keep.append ( fun1 )
        result.aux_keep.append ( fun2 )
        return result
    
    return NotImplemented 

# =============================================================================
## operator for "3D-function + other"
def _f3_add_ ( self , other ) :
    """Operator for ``3D-function + other''"""
    return _f3_op_ ( self , other , Ostap.MoreRooFit.Addition    , "Add_%s_%s"     )

# =============================================================================
## operator for "3D-function - other"
def _f3_sub_ ( self , other ) :
    """Operator for ``3D-function - other''"""
    return _f3_op_ ( self , other , Ostap.MoreRooFit.Subtraction , "Subtract_%s_%s" )

# =============================================================================
## operator for "3D-function * other"
def _f3_mul_ ( self , other ) :
    """Operator for ``3D-function * other''"""
    return _f3_op_ ( self , other , Ostap.MoreRooFit.Product     , "Product_%s_%s"  )

# =============================================================================
## operator for "3D-function / other"
def _f3_div_ ( self , other ) :
    """Operator for ``3D-function / other''"""
    return _f3_op_ ( self , other , Ostap.MoreRooFit.Division    , "Divide_%s_%s"  )

# =============================================================================
## operator for "3D-function ** other"
def _f3_pow_ ( self , other ) :
    """Operator for ``3D-function **  other''"""
    return _f3_op_ ( self , other , Ostap.MoreRooFit.Power       , "Pow_%s_%s"  )
        

FUNC3.__add__     = _f3_add_
FUNC3.__mul__     = _f3_mul_
FUNC3.__sub__     = _f3_sub_
FUNC3.__div__     = _f3_div_
FUNC3.__truediv__ = _f3_div_
FUNC3.__pow__     = _f3_pow_


# =============================================================================
## Operator for `3D-function (op) other`:
def _f3_rop_ ( fun1 , fun2 , ftype , fname ) :
    """ Operator for `3D-function (op) other`:
    """
    
    xvar, yvar, zvar = fun1.xvar , fun1.yvar , fun1.zvar

    if isinstance ( fun2 , num_types ) :
        
        fun2   = ROOT.RooRealConstant.value ( float ( fun2 )  ) 
        s      = ftype ( fun2 , fun1.fun  )
        result = Fun3D ( s ,
                         xvar = xvar ,
                         yvar = yvar ,
                         zvar = zvar ,
                         name = fname % ( fun2.name , fun1.name  ) )
        
        result.aux_keep.append ( fun1 )
        result.aux_keep.append ( fun2 )
        return result
    
    elif isinstance ( fun2 , ROOT.RooAbsReal ) :
        
        s      = ftype ( fun2 , fun1.fun )
        result = Fun3D ( s ,
                         xvar = xvar  ,
                         yvar = yvar  ,
                         zvar = zvar  ,
                         name = fname % ( fun2.name , fun1.name ) )
        
        result.aux_keep.append ( fun1 )
        result.aux_keep.append ( fun2 )
        return result
    
    return NotImplemented 


# =============================================================================
## operator for "3D-function + other"
def _f3_radd_ ( self , other ) :
    """Operator for ``3D-function + other''"""
    return _f3_rop_ ( self , other , Ostap.MoreRooFit.Addition    , "Add_%s_%s"     )

# =============================================================================
## operator for "3D-function - other"
def _f3_rsub_ ( self , other ) :
    """Operator for ``3D-function - other''"""
    return _f3_rop_ ( self , other , Ostap.MoreRooFit.Subtraction , "Subtract_%s_%s" )

# =============================================================================
## operator for "3D-function * other"
def _f3_rmul_ ( self , other ) :
    """Operator for ``3D-function * other''"""
    return _f3_rop_ ( self , other , Ostap.MoreRooFit.Product     , "Product_%s_%s"  )

# =============================================================================
## operator for "3D-function / other"
def _f3_rdiv_ ( self , other ) :
    """Operator for ``3D-function / other''"""
    return _f3_rop_ ( self , other , Ostap.MoreRooFit.Division    , "Divide_%s_%s"  )

# =============================================================================
## operator for "3D-function ** other"
def _f3_rpow_ ( self , other ) :
    """Operator for ``3D-function **  other''"""
    return _f3_rop_ ( self , other , Ostap.MoreRooFit.Power       , "Pow_%s_%s"  )
        

FUNC3.__radd__     = _f3_radd_
FUNC3.__rmul__     = _f3_rmul_
FUNC3.__rsub__     = _f3_rsub_
FUNC3.__rdiv__     = _f3_rdiv_
FUNC3.__rtruediv__ = _f3_rdiv_
FUNC3.__rpow__     = _f3_rpow_


# =============================================================================
## Operator for `2D-function (op) other`:
def _f2_op_ ( fun1 , fun2 , ftype , fname ) :
    """ Operator for `2D-function (op) other`:
    """
    
    xvar, yvar, zvar  = fun1.xvar , fun1.yvar, None 
    
    if   isinstance ( fun2 , FUNC3 ) :

        if xvar in fun2.vars and yvar in fun2.vars :

            ## operation 
            s      = ftype ( fun1.fun , fun2.fun )
            result = Fun3D ( s ,
                             xvar = fun2.xvar ,
                             yvar = fun2.yvar ,
                             zvar = fun2.zvar ,
                             name = fname     % ( fun1.name , fun2.name ) )
            
            result.aux_keep.append ( fun1 )
            result.aux_keep.append ( fun2 )
            return result
        
    elif isinstance ( fun2 , FUNC2 ) :

        if xvar in fun2.vars and yvar in fun2.vars :
            
            ## operation  
            s       = ftype ( self.fun , other.fun )
            retsult = Fun2D ( s ,
                              xvar = xvar  ,
                              yvar = yvar  ,
                              name = fname % ( fun1.name , fun2.name ) )
            
            result.aux_keep.append ( fun1 )
            result.aux_keep.append ( fun2 )
            return result

        elif ( xvar in fun2.vars ) or ( yvar in fun2.vars ) : 
            
            if zvar is None : zvar = fun2.xvar if not fun2.xvar in fun1.vars else None 
            if zvar is None : zvar = fun2.yvar if not fun2.yvar in fun1.vars else None
            
            ## operation  
            s      = ftype ( self.fun , other.fun )
            result = Fun3D ( s ,
                             xvar = xvar  ,
                             yvar = yvar  ,
                             zvar = zvar  ,
                             name = fname % ( fun1.name , fun2.name ) )
        
            result.aux_keep.append ( fun1 )
            result.aux_keep.append ( fun2 )
            return result

    elif isinstance ( fun2 , FUNC  ) :

        if fun2.xvar in fun1.vars :
            
            ## operation 
            s      = ftype ( fun1.fun , fun2.fun )
            result = Fun2D ( s ,
                             xvar = xvar ,
                             yvar = yvar ,
                             name = fname % ( fun1.name , fun2.name ) )
            
            result.aux_keep.append ( fun1 )
            result.aux_keep.append ( fun2 )
            return result
            
        else :
            
            ## operation
            s      = ftype ( fun1.fun , fun2.fun )
            result = Fun3D ( s ,
                             xvar =      xvar ,
                             yvar =      yvar ,
                             zvar = fun2.xvar , 
                             name = fname     % ( fun1.name , fun2.name ) )

            result.aux_keep.append ( fun1 )
            result.aux_keep.append ( fun2 )
            return result

    elif isinstance ( fun2 , num_types ) :
        
        fun2   = ROOT.RooRealConstant.value ( float ( fun2 )  ) 
        s      = ftype ( fun1.fun , fun2 )
        result = Fun2D ( s ,
                         xvar = xvar  ,
                         yvar = yvar  ,
                         name = fname % ( fun1.name , fun2.name ) )
        
        result.aux_keep.append ( fun1 )
        result.aux_keep.append ( fun2 )
        return result

    elif isinstance ( fun2 , ROOT.RooAbsReal ) :
        
        s      = ftype  ( fun1.fun , fun2 )
        result = Fun2D ( s ,
                       xvar = xvar  ,
                       yvar = yvar  ,
                       name = fname % ( fun1.name , fun2.name ) )
        
        result.aux_keep.append ( fun1 )
        result.aux_keep.append ( fun2 )
        return result

    return NotImplemented 


# =============================================================================
## operator for "2D-function + other"
def _f2_add_ ( self , other ) :
    """Operator for ``2D-function + other''"""
    return _f2_op_ ( self , other , Ostap.MoreRooFit.Addition    , "Add_%s_%s"      )

# =============================================================================
## operator for "2D-function - other"
def _f2_sub_ ( self , other ) :
    """Operator for ``2D-function - other''"""
    return _f2_op_ ( self , other , Ostap.MoreRooFit.Subtraction , "Subtract_%s_%s" )

# =============================================================================
## operator for "3D-function * other"
def _f2_mul_ ( self , other ) :
    """Operator for ``2D-function * other''"""
    return _f2_op_ ( self , other , Ostap.MoreRooFit.Product     , "Product_%s_%s"  )

# =============================================================================
## operator for "2D-function / other"
def _f2_div_ ( self , other ) :
    """Operator for ``2D-function / other''"""
    return _f2_op_ ( self , other , Ostap.MoreRooFit.Division    , "Divide_%s_%s"  )
        
# =============================================================================
## operator for "2D-function ** other"
def _f2_pow_ ( self , other ) :
    """Operator for ``2D-function **  other''"""
    return _f2_op_ ( self , other , Ostap.MoreRooFit.Power       , "Pow_%s_%s"     )
        

FUNC2.__add__     = _f2_add_
FUNC2.__mul__     = _f2_mul_
FUNC2.__sub__     = _f2_sub_
FUNC2.__div__     = _f2_div_
FUNC2.__truediv__ = _f2_div_
FUNC2.__pow__     = _f2_pow_


# =============================================================================
## Operator for `2D-function (op) other`:
def _f2_rop_ ( fun1 , fun2 , ftype , fname ) :
    """ Operator for `2D-function (op) other`:
    """
    
    xvar, yvar, zvar  = fun1.xvar , fun1.yvar, None 
    
    if isinstance ( fun2 , num_types ) :
        
        fun2   = ROOT.RooRealConstant.value ( float ( fun2 )  ) 
        s      = ftype ( fun2 ,  fun1.fun )
        result = Fun2D ( s ,
                         xvar = xvar  ,
                         yvar = yvar  ,
                         name = fname % ( fun2.name , fun1.name ) )
        
        result.aux_keep.append ( fun1 )
        result.aux_keep.append ( fun2 )
        return result

    elif isinstance ( fun2 , ROOT.RooAbsReal ) :
        
        s      = ftype ( fun2 , fun1.fun )
        result = Fun2D ( s ,
                       xvar = xvar  ,
                       yvar = yvar  ,
                       name = fname % ( fun2.name , fun1.name ) )
        
        result.aux_keep.append ( fun1 )
        result.aux_keep.append ( fun2 )
        return result

    return NotImplemented 

# =============================================================================
## operator for "2D-function + other"
def _f2_radd_ ( self , other ) :
    """Operator for ``2D-function + other''"""
    return _f2_rop_ ( self , other , Ostap.MoreRooFit.Addition    , "Add_%s_%s"      )

# =============================================================================
## operator for "2D-function - other"
def _f2_rsub_ ( self , other ) :
    """Operator for ``2D-function - other''"""
    return _f2_rop_ ( self , other , Ostap.MoreRooFit.Subtraction , "Subtract_%s_%s" )

# =============================================================================
## operator for "3D-function * other"
def _f2_rmul_ ( self , other ) :
    """Operator for ``2D-function * other''"""
    return _f2_rop_ ( self , other , Ostap.MoreRooFit.Product     , "Product_%s_%s"  )

# =============================================================================
## operator for "2D-function / other"
def _f2_rdiv_ ( self , other ) :
    """Operator for ``2D-function / other''"""
    return _f2_rop_ ( self , other , Ostap.MoreRooFit.Division    , "Divide_%s_%s"  )
        
# =============================================================================
## operator for "2D-function ** other"
def _f2_rpow_ ( self , other ) :
    """Operator for ``2D-function **  other''"""
    return _f2_rop_ ( self , other , Ostap.MoreRooFit.Power       , "Pow_%s_%s"     )
        

FUNC2.__radd__     = _f2_radd_
FUNC2.__rmul__     = _f2_rmul_
FUNC2.__rsub__     = _f2_rsub_
FUNC2.__rdiv__     = _f2_rdiv_
FUNC2.__rtruediv__ = _f2_rdiv_
FUNC2.__rpow__     = _f2_rpow_

# =============================================================================
## Operator for `1D-function (op) other`:
def _f1_op_ ( fun1 , fun2 , ftype , fname ) :
    """ Operator for `1D-function (op) other`:"""
    
    xvar, yvar, zvar  = fun1.xvar , None , None 
    
    if   isinstance ( fun2 , FUNC3 ) :

        if xvar in fun2.vars :
            
            ## operation 
            s      = ftype ( fun1.fun , fun2.fun )
            result = Fun3D ( s ,
                             xvar = fun2.xvar ,
                             yvar = fun2.yvar ,
                             zvar = fun2.zvar ,
                             name = fname     % ( fun1.name , fun2.name ) )
            
            result.aux_keep.append ( fun1 )
            result.aux_keep.append ( fun2 )
            return result
    
    elif isinstance ( fun2 , FUNC2 ) :

        if xvar in fun2.vars : 
            
            ## operation  
            s      = ftype ( self.fun , other.fun )
            result = Fun2D ( s ,
                           xvar = fun2.xvar ,
                           yvar = fun2.yvar ,
                           name = fname     % ( fun1.name , fun2.name ) )
            
            result.aux_keep.append ( fun1 )
            result.aux_keep.append ( fun2 )
            return result
        
        else :
            
            ## operation  
            s      = ftype ( self.fun , other.fun )
            result = Fun3D ( s ,
                             xvar = fun2.xvar ,
                             yvar = fun2.yvar ,
                             zvar =      xvar ,
                             name = fname     % ( fun1.name , fun2.name ) )
            
            result.aux_keep.append ( fun1 )
            result.aux_keep.append ( fun2 )
            return result

    elif isinstance ( fun2 , FUNC  ) :

        if fun2.xvar in fun1.vars :
            
            ## operation
            s      = ftype ( fun1.fun , fun2.fun )
            result = Fun1D ( s ,
                             xvar = xvar      ,
                             name = fname     % ( fun1.name , fun2.name ) )
            
            result.aux_keep.append ( fun1 )
            result.aux_keep.append ( fun2 )
            return result

        else :
            
            ## operation 
            s      = ftype ( fun1.fun , fun2.fun )
            result = Fun2D ( s ,
                             xvar =      xvar ,
                             yvar = fun2.xvar ,
                             name = fname     % ( fun1.name , fun2.name ) )
            
            result.aux_keep.append ( fun1 )
            result.aux_keep.append ( fun2 )
            return result

    elif isinstance ( fun2 , num_types ) :
        
        fun2   = ROOT.RooRealConstant.value ( float ( fun2 )  ) 
        s      = ftype ( fun1.fun , fun2 )
        result = Fun1D ( s ,
                         xvar = xvar  ,
                         name = fname % ( fun1.name , fun2.name ) )
        
        result.aux_keep.append ( fun1 )
        result.aux_keep.append ( fun2 )
        return result


    elif isinstance ( fun2 , ROOT.RooAbsReal ) :
        
        s      = ftype ( fun1.fun , fun2 )
        result = Fun1D ( s ,
                         xvar = xvar  ,
                         name = fname % ( fun1.name , fun2.name ) )
        
        result.aux_keep.append ( fun1 )
        result.aux_keep.append ( fun2 )
        return result

    return NotImplemented 

# =============================================================================
## operator for "1D-function + other"
def _f1_add_ ( self , other ) :
    """Operator for ``1D-function + other''"""
    return _f1_op_ ( self , other , Ostap.MoreRooFit.Addition    , "Add_%s_%s"     )

# =============================================================================
## operator for "1D-function - other"
def _f1_sub_ ( self , other ) :
    """Operator for ``2D-function - other''"""
    return _f1_op_ ( self , other , Ostap.MoreRooFit.Subtraction , "Subtract_%s_%s" )

# =============================================================================
## operator for "1D-function * other"
def _f1_mul_ ( self , other ) :
    """Operator for ``2D-function * other''"""
    return _f1_op_ ( self , other , Ostap.MoreRooFit.Product     , "Product_%s_%s"  )

# =============================================================================
## operator for "1D-function / other"
def _f1_div_ ( self , other ) :
    """Operator for ``1D-function / other''"""
    return _f1_op_ ( self , other , Ostap.MoreRooFit.Division    , "Divide_%s_%s"  )
        
# =============================================================================
## operator for "1D-function ** other"
def _f1_pow_ ( self , other ) :
    """Operator for ``1D-function **  other''"""
    return _f1_op_ ( self , other , Ostap.MoreRooFit.Power       , "Pow_%s_%s"  )
        

FUNC.__add__     = _f1_add_
FUNC.__mul__     = _f1_mul_
FUNC.__sub__     = _f1_sub_
FUNC.__div__     = _f1_div_
FUNC.__truediv__ = _f1_div_
FUNC.__pow__     = _f1_pow_

# =============================================================================
## Operator for `1D-function (op) other`:
def _f1_rop_ ( fun1 , fun2 , ftype , fname ) :
    """ Operator for `1D-function (op) other`:"""
    
    xvar, yvar, zvar  = fun1.xvar , None , None 
    
    if isinstance ( fun2 , num_types ) :
        
        fun2   = ROOT.RooRealConstant.value ( float ( fun2 )  ) 
        s      = ftype ( fun2 ,  fun1.fun )
        result = Fun1D ( s ,
                         xvar = xvar  ,
                         name = fname % ( fun2.name , fun1.name ) )
        
        result.aux_keep.append ( fun1 )
        result.aux_keep.append ( fun2 )
        return result


    elif isinstance ( fun2 , ROOT.RooAbsReal ) :
        
        s      = ftype ( fun2 , fun1.fun )
        result = Fun1D ( s ,
                         xvar = xvar  ,
                         name = fname % ( fun2.name , fun1.name ) )
        
        result.aux_keep.append ( fun1 )
        result.aux_keep.append ( fun2 )
        return result

    return NotImplemented 


# =============================================================================
## operator for "1D-function + other"
def _f1_radd_ ( self , other ) :
    """Operator for ``1D-function + other''"""
    return _f1_rop_ ( self , other , Ostap.MoreRooFit.Addition    , "Add_%s_%s"     )

# =============================================================================
## operator for "1D-function - other"
def _f1_rsub_ ( self , other ) :
    """Operator for ``2D-function - other''"""
    return _f1_rop_ ( self , other , Ostap.MoreRooFit.Subtraction , "Subtract_%s_%s" )

# =============================================================================
## operator for "1D-function * other"
def _f1_rmul_ ( self , other ) :
    """Operator for ``2D-function * other''"""
    return _f1_rop_ ( self , other , Ostap.MoreRooFit.Product     , "Product_%s_%s"  )

# =============================================================================
## operator for "1D-function / other"
def _f1_rdiv_ ( self , other ) :
    """Operator for ``1D-function / other''"""
    return _f1_rop_ ( self , other , Ostap.MoreRooFit.Division    , "Divide_%s_%s"  )
        
# =============================================================================
## operator for "1D-function ** other"
def _f1_rpow_ ( self , other ) :
    """Operator for ``1D-function **  other''"""
    return _f1_rop_ ( self , other , Ostap.MoreRooFit.Power       , "Pow_%s_%s"  )
        

FUNC.__radd__     = _f1_radd_
FUNC.__rmul__     = _f1_rmul_
FUNC.__rsub__     = _f1_rsub_
FUNC.__rdiv__     = _f1_rdiv_
FUNC.__rtruediv__ = _f1_rdiv_
FUNC.__rpow__     = _f1_rpow_

        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
#                                                                       The END 
# =============================================================================
