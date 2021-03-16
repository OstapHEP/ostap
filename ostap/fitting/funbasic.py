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
    'make_fun'          , ## helper functno to reate FUNC-objects
    'SETPARS'           , ## context manager to keep/preserve parameters 
    )
# =============================================================================
import ROOT, math, sys 
from   sys                           import version_info as python_version 
from   ostap.core.ostap_types        import ( integer_types  , num_types   ,
                                              dictlike_types , list_types  ,
                                              is_good_number )     
from   ostap.core.core               import Ostap , valid_pointer
from   ostap.fitting.variables       import SETVAR
from   ostap.logger.utils            import roo_silent , rootWarning
from   ostap.fitting.roofit          import PDF_fun 
from   ostap.fitting.utils           import MakeVar, XVar, YVar, ZVar, NameDuplicates  
import ostap.fitting.variables
import ostap.fitting.roocollections
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.funbasic' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================
py2 = 2 >= sys.version_info.major  
# =============================================================================
## helper factory function
def func_factory ( klass , config ) :
    """Helper factory function, used for unpickling"""
    with NameDuplicates ( True ) : 
        return klass ( **config ) 
# =============================================================================
## @class FUNC
#  Helper base class for impolementation of various (Roo)Function-wrappers
class FUNC(XVar) :
    """Helper base class for implementation of various (Roo)Function-wrappers
    """
    def __new__( cls, *args, **kwargs):
        if  python_version.major > 2 : obj = super(FUNC, cls).__new__( cls )
        else                         : obj = super(FUNC, cls).__new__( cls , *args , **kwargs )
        ##
        obj.__func_init = False  
        return obj
        
    def __init__ ( self , name , xvar ) :

        if self.__func_init : return 
        else                : self.__func_init = True  
        
        ## name is defined via the base class MakeVar 
        self.name  = name ## name is defined via the base class MakeVar 
        
        ##  super(FUNC,self).__init__ ( xvar )
        XVar .__init__ ( self , xvar )

        self.__vars      = ROOT.RooArgSet  ()
        self.__variables = [] 

        self.__config       = {}
        self.__fun          = None
        
        self.__tricks       = True
        self.__fit_result   = None
        
        self.__draw_var     = None
        self.__draw_options = {} ## predefined drawing options for this FUNC/PDF

        self.__checked_keys = set()
        
        self.__dfdx = None
        self.__intx = None
        
        self.vars.add         ( self.xvar )
        self.variables.append ( self.xvar )

        self.config = { 'name' : self.name , 'xvar' : self.xvar  }
        self.__func_init = True  

        ## derived functions/objects
        self.__derived = {}
        
    ## pickling via reduce 
    def __reduce__ ( self ) :
        if py2 : return func_factory , ( type ( self ) , self.config, )
        else   : return type ( self ).factory , ( self.config, )
    ## factory method 
    @classmethod
    def factory ( klass , config ) :
        """Factory method, used for unpickling"""
        with NameDuplicates ( True ) : 
            return klass ( **config ) 
        
    ## conversion to string 
    def __str__ (  self ) :
        return '%s(%s,xvar=%s)' % ( self.__class__.__name__ , self.name , self.xvar.name )
    __repr__ = __str__ 

    @property 
    def vars ( self ) :
        """``vars'' : variables/observables (as ROOT.RooArgSet)"""
        return self.__vars

    @property 
    def variables ( self ) :
        """``variables'' : variables/observables at list"""
        return self.__variables

    @property
    def fun  ( self ) :
        """The actual function (ROOT.RooAbsReal)"""
        return self.__fun
    @fun.setter
    def fun  ( self , value ) :
        if value is None :
            self.__fun = value
            return        
        assert isinstance ( value , ROOT.RooAbsReal ) , "``pdf/fun'' is not ROOT.RooAbsReal"
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
        if not self.fun : return False
        ##
        if var and isinstance ( var , ROOT.RooAbsReal ) :
            return self.fun.depends_on ( var ) 

        ## 
        if isinstance ( var , FUNC ) :
            if self.fun.depends_on ( var.fun  ) : return True
            if self.fun.depends_on ( var.vars ) : return True
                
        return False 
            
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
    def checked_keys ( self ) :
        """``checked keys'' : special keys for clone-method
        """
        return self.__checked_keys

    # =========================================================================
    ##  Convert function to PDF
    #   @code
    #   fun = ...
    #   pdf = fun.as_PDF () 
    #   @endcode
    #   - for PDC, <code>self</code> is returned
    #   - otherwise <code>RooWrappedPdf</code> is used
    #   @see RooWrapperPdf 
    def as_PDF ( self  , name ) :
        """Convert function to PDF
        >>> fun = ...
        >>> pdf = fun.as_PDF () 
        - for PDC, `self` is returned
        - otherwise `RooWrappedPdf` is used
        - see ROOT.RooWrapperPdf 
        """
        from ostap.fit.basic import PDF, make_pdf 
        if    isinstance ( self     , PDF           ) : return self

        name = name if name else "PDF_from_%s" % self.name
        ##
        return make_pdf ( self.fun , name = name , *self.variables )

    # =========================================================================
    @property 
    def dfdx ( self ) :
        """``dfdx'' : derivative dF/dX"""
        return self.__dfdx
    @dfdx.setter
    def dfdx ( self , value ) :
        assert value is None or isinstance ( value , FUNC ) , \
               "``dfdx'' has invalid type %s/%s" % ( value , type ( value ) )
        self.__dfdx = value 

    # =========================================================================
    @property 
    def intx ( self ) :
        """``intx'' : running integral over X"""
        return self.__intx
    @intx.setter
    def intx ( self , value ) :
        assert value is None or isinstance ( value , FUNC ) , \
               "``intx'' has invalid type %s/%s" % ( value , type ( value ) )
        self.__intx = value 

    @property
    def draw_options ( self ) :
        """``draw_options'' : dictionary with predefined draw-options for this PDF/FUN
        """
        return self.__draw_options

    # =========================================================================

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
        assert self.fun, "FUNC::params: Function is invalid!"        
        return self.fun.getParameters ( 0 ) if dataset is None else self.fun.getParameters ( dataset )
    
    # =========================================================================
    ## get the parameter value by name
    #  @code
    #  fun = ...
    #  p   = fun.parameter  ( 'A' )
    #  @endcode
    def parameter ( self , param , dataset = None ) :
        """Get the parameter value by name
        >>> pdf = ...
        >>> p   = pdf.parameter  ( 'A' )
        """
        ## get the list of the actual parameters 
        pars = self.params ( dataset )

        for p in pars :
            if p.name == param : return p
            
        self.error ( "No parameter %s defined" % param )
        raise KeyError ( "No parameter %s defined" % param )

    # =========================================================================
    ## get all parameters/variables in form of dictionary
    #  @code
    #  fun    = ...
    #  params = fun.parameters ( dataset ) 
    #  @endcode
    def parameters ( self , dataset = None ) :
        """ Get all parameters/variables in form of dictionary
        >>> fun    = ...
        >>> params = fun.parameters ( dataset ) 
        """
        
        ## get the list of the actual parameters 
        pars = self.params ( dataset ) 

        tmp    = {}
        for p in pars :
            if not isinstance ( p, ROOT.RooAbsCategory ) :
                tmp [ p.name ] = p.value
                
        keys   = tmp.keys()
        result = {} 
        for key in sorted ( keys ) : result [ key ] = tmp [ key ] 
            
        return result 

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
        ## get the list of the actual parameters 
        pars = self.params ( )
        for p in pars :
            if p.name == param : return p
        raise KeyError ( "No parameter %s defined" % param )


    # =========================================================================
    ## Load parameters from:
    #    - external dictionary <code>{ name : value }</code>
    #    - sequence of <code>RooAbsReal</code> object
    #    - <code>ROOT.RooFitResult</code> object 
    #  @code
    #  pdf     = ...
    #  dataset = ...
    #  params  = { 'A' : 10 , 'B' : ... }
    #  pdf.load_params ( dataset , params ) 
    #  params  = ( A , B , C , ... )
    #  pdf.load_params ( dataset , params )  
    #  @endcode 
    def load_params ( self , params = {} , dataset = None , silent = False  ) :
        """Load parameters from
        - external dictionary `{ name : value }`
        - sequence of `RooAbsReal` objects
        - `RooFitResult` object
        
        >>> pdf      = ...
        >>> dataset = ... 
        >>> params = { 'A' : 10 , 'B' : ... }
        >>> pdf.load_params ( dataset , params ) 
        >>> params = ( A , B , C , ... )
        >>> pdf.load_params ( dataset , params )  
        """

        if dataset :
            assert     isinstance ( dataset , ROOT.RooAbsData ) , "load_params: invalid type of ``dataset':%s'" % type ( dataset ) 
        else :
            dataset = ROOT.nullptr
            
        if params :
            assert not isinstance ( params  , ROOT.RooAbsData ) , "load_params: invalid type of ``params'':%s" % type ( params ) 
            
        if isinstance ( params , ROOT.RooFitResult ) :
            params = params.dct_params () 
        
        ## get the list of the actual parameters 
        pars = self.params ( dataset ) 

        table = [] 
        if isinstance ( params , dictlike_types ) :
            keys   = set () 
            for key in params :
                for p in pars :
                    if not hasattr ( p  , 'setVal' ) : continue
                    if p.name != key                 : continue
                    
                    v  = params[key]
                    vv = float ( v  )
                    pv = p.getVal ()   
                    if vv != pv : 
                        p.setVal   ( vv )
                        item = p.name , "%-14.6g" % pv , "%-+14.6g" % vv 
                        table.append ( item ) 
                    keys.add ( key )

            not_used = set ( params.keys() ) - keys 

        ## list of objects 
        else :
            
            keys = set()        
            for i , pp in enumerate ( params ) :  
                if not isinstance ( pp , ROOT.RooAbsReal ) : continue
                for p in pars :
                    if not hasattr ( p  , 'setVal' )       : continue
                    if p.name != pp.name                   : continue
                    
                    vv = float ( pp )
                    pv = p.getVal () 
                    if vv != pv :
                        p.setVal   ( vv )
                        item = p.name , "%-14.6g" % pv , "%-+14.6g" % vv 
                        table.append ( item ) 
                    keys.add  ( i )

            not_used = []
            for i , pp in enumerate ( params ) :  
                if i in keys : continue
                not_used.append ( pp )

        if not silent :
            
            table.sort()
            npars = len ( table )
            
            if npars :            
                title = 'Parameters loaded: %s' % npars 
                table = [ ('Parameter' ,'old value' , 'new value' ) ] + table
                import ostap.logger.table
                table = ostap.logger.table.table ( table , title , prefix = "# " )
                
            self.info ( "%s parameters loaded:\n%s" % ( npars , table ) ) 
            
            not_used = list ( not_used )
            not_used.sort() 
            if not_used :
                self.warning ("Following keys are unused %s" % not_used ) 
        
        return 
    
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
        """Make a clone for the given fucction/PDF with
        the optional replacement of the certain parameters
        >>> xpdf = ...
        >>> ypdf = xpdf.clone ( xvar = yvar ,  name = 'PDFy' ) 
        """

        ## get config 
        conf = {}
        
        ## check for the special keys 
        for k in kwargs :
            if k in self.checked_keys :
                if not hasattr ( self , k ) or getattr ( self, k ) != kwargs [ k ] :
                    raise AttributeError ("Class %s cannot be cloned with ``%s'' key" % ( self.__class__.__name__ , k ) )

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

        ## KLASS = self.__class__
        ## cloned = KLASS ( **conf )

        cloned = self.factory ( conf )
            
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
        
        if self.xminmax() : 
            xmin , xmax = self.xminmax()
            xmin = kwargs.pop ( 'xmin' , xmin )
            xmax = kwargs.pop ( 'xmax' , xmax )
        else : 
            xmin = kwargs.pop ( 'xmin' , None )
            xmax = kwargs.pop ( 'xmax' , None )

        assert not xmin is None  , 'xmin is not defined!' 
        assert not xmax is None  , 'xmax is not defined!' 

        
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
                if   error and isinstance ( error , ROOT.RooFitResult ) : 
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
        if self.xminmax() :
            xmin , xmax = self.xminmax ()
            if x < xmin or x > xmax : return 0

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
        if self.xminmax () :
            xmn , xmx = self.xminmax() 
            xmin = xmn if xmin is None else max ( xmin , xmn )
            xmax = xmx if xmax is None else min ( xmax , xmx )
            
        assert is_good_number ( xmin ) , 'xmin is not specified!'
        assert is_good_number ( xmax ) , 'xmax is not specified!'
        assert xmin < xmax             , 'Invalid xmin/xmax : %s/%s' % ( xmin , xmax )

        if x0 is None  :  x0 = 0.5 * ( xmin + xmax )
        
        if not xmin <= x0 <= xmax :
            self.error("Wrong xmin/x0/xmax: %s/%s/%s"   % ( xmin , x0 , xmax ) )
        
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
        
        if self.xminmax () :
            xmn , xmx = self.xminmax() 
            xmin = xmn if xmin is None else max ( xmin , xmn )
            xmax = xmx if xmax is None else min ( xmax , xmx )
            
        assert is_good_number ( xmin ) , 'xmin is not specified!'
        assert is_good_number ( xmax ) , 'xmax is not specified!'
        assert xmin < xmax             , 'Invalid xmin/xmax : %s/%s' % ( xmin , xmax )

        if x0 is None  :  x0 = 0.5 * ( xmin + xmax )

        if not xmin <= x0 <= xmax :
            self.error("Wrong xmin/x0/xmax: %s/%s/%s"   % ( xmin , x0 , xmax ) )
        
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
            groot = ROOT.ROOT.GetROOT()
            if not groot.IsBatch() :
                with rootWarning (): frame.draw ( kwargs.pop ( 'draw_options','' ) )
            
            if kwargs :
                self.warning("draw: ignored unknown options: %s" % list( kwargs.keys() ) ) 

            return frame

# =============================================================================
## Context manager to keep/preserve the parameters for function/pdf
#  @code
#  pdf = ...
#  with SETPARS ( pdf ) :
#  ...   <do something here with pdf>
#
#  @endcode 
class SETPARS(object) :
    """Context manager to keep/preserve the parameters for function/pdf
    >>> pdf = ...
    >>> with SETPARS ( pdf ) :
    ...   <do something here with pdf>
    """

    ## constructor with the function and dataset 
    def __init__ ( self , fun , dataset = None ) :

        self.__params  = {}
        self.__fun     = fun        
        if dataset is None : dataset = ROOT.nullptr         
        self.__dataset = dataset 

    ## context manager: ENTER 
    def __enter__ ( self ) :

        params = self.__fun.parameters ( self.__dataset )
        for par in params :
            self.__params [ par ] = float ( params [ par ] )
            
        return self
    
    ## context manager: EXIT
    def __exit__ ( self , *_ ) :
        
        if self.__params : 
            self.__fun.load_params ( params = self.__params , dataset = self.__dataset , silent = True )
            
        self.__fun     = None
        self.__params  = {}
        self.__dataset = ROOT.nullptr 
        
    @property
    def fun ( self ) :
        """``fun'': the actual function/pdf"""
        return self.__fun
    
    @property
    def dataset ( self ) :
        """``dataset'': the dataset"""
        return self.__dataset

    @property
    def params( self ) :
        """``params'': dictionary of parameters"""
        return self.__params

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
        
        self.checked_keys.add  ( 'fun'  ) 
        self.checked_keys.add  ( 'xvar' ) 

    
# =============================================================================
## @class FUNC2
#  The base class for 2D-function
class FUNC2(FUNC,YVar) :
    """Base class for 2D-function
    """
    def __new__( cls, *args, **kwargs):
        if  python_version.major > 2 : obj = super(FUNC2, cls).__new__( cls )
        else                         : obj = super(FUNC2, cls).__new__( cls , *args , **kwargs )
        ##
        obj.__func2_init = False  
        return obj
    
    def __init__ ( self , name , xvar , yvar ) :

        if self.__func2_init : return 
        else                 : self.__func2_init = True  
        
        FUNC .__init__ ( self , name , xvar )
        YVar .__init__ ( self , yvar )
        
        self.__dfdy = None 
        self.__inty = None 
        
        self.vars     .add    ( self.yvar )
        self.variables.append ( self.yvar )
        
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


    # =========================================================================
    @property
    def dfdy ( self ) :
        """``dfdy'': derivative dF/dY"""
        return self.__dfdy
    @dfdy.setter
    def dfdy ( self , value ) :
        assert value is None or isinstance ( value , FUNC2 ) , \
               "``dfdy'' has invalid type %s/%s" % ( value , type ( value ) )
        self.__dfdy = value 
    
    # =========================================================================
    @property 
    def inty ( self ) :
        """``inty'' : running integral over Y"""
        return self.__inty
    @inty.setter
    def inty ( self , value ) :
        assert value is None or isinstance ( value , FUNC2 ) , \
               "``inty'' has invalid type %s/%s" % ( value , type ( value ) )
        self.__inty = value 


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


    # =========================================================================
    ## disabled 
    def __disabled ( self , name ) :
        """Disabled method"""
        self.error ('Method "%s" is disabled for FUNC2/DUNC3' % name )
        raise NotImplementedError('Method "%s"is disabled for FUNC2/DUNC3' % name  )
    
    # =========================================================================
    def minimum    ( self , *args , **kwargs ) : return self.__disabled ( "minimum"    ) 
    def maximum    ( self , *args , **kwargs ) : return self.__disabled ( "maximum"    ) 
    def fwhm       ( self , *args , **kwargs ) : return self.__disabled ( "fwhm"       ) 
    def mid_point  ( self , *args , **kwargs ) : return self.__disabled ( "mid_point"  ) 
    def mode       ( self , *args , **kwargs ) : return self.__disabled ( "mode"       ) 
    def derivative ( self , *args , **kwargs ) : return self.__disabled ( "derivative" ) 

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

        if fun is xvar : fun = Ostap.MoreRooFit.Id ( "" , "" , xvar )
        if fun is yvar : fun = Ostap.MoreRooFit.Id ( "" , "" , yvar )
            
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
        
        self.checked_keys.add  ( 'fun'  ) 
        self.checked_keys.add  ( 'xvar' ) 
        self.checked_keys.add  ( 'yvar' ) 

    
# =============================================================================
## @class FUNC3
#  The base class for 3D-function
class FUNC3(FUNC2,ZVar) :
    """Base class for 3D-function
    """
    def __new__( cls, *args, **kwargs):
        if  python_version.major > 2 : obj = super(FUNC3, cls).__new__( cls )
        else                         : obj = super(FUNC3, cls).__new__( cls , *args , **kwargs )
        ##
        obj.__func3_init = False  
        return obj

    def __init__ ( self , name , xvar , yvar , zvar ) :

        
        if self.__func3_init : return 
        else                 : self.__func3_init = True  

        FUNC2.__init__ ( self , name , xvar , yvar )
        ZVar .__init__ ( self , zvar )
        
        self.__dfdz = None 
        self.__intz = None 
        
        self.vars     .add    ( self.zvar )
        self.variables.append ( self.zvar )
        
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

    # =========================================================================
    @property
    def dfdz ( self ) :
        """``dfdz'': derivative dF/dZ"""
        return self.__dfdz
    @dfdz.setter
    def dfdz ( self , vallue ) :
        assert value is None or isinstance ( value , FUNC3 ) , \
               "``dfdz'' has invalid type %s/%s" % ( value , type ( value ) )
        self.__dfdz = value 
    
    # =========================================================================
    @property 
    def intz ( self ) :
        """``intz'' : (running) integral over Z"""
        return self.__intz
    @intz.setter
    def intz ( self , value ) :
        assert value is None or isinstance ( value , FUNC3 ) , \
               "``intz'' has invalid type %s/%s" % ( value , type ( value ) )
        self.__intz = value 
    
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
        
        if fun is xvar : fun = Ostap.MoreRooFit.Id ( "" , "" , xvar )
        if fun is yvar : fun = Ostap.MoreRooFit.Id ( "" , "" , yvar )
        if fun is zvar : fun = Ostap.MoreRooFit.Id ( "" , "" , zvar )
            
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

        self.checked_keys.add  ( 'fun'  ) 
        self.checked_keys.add  ( 'xvar' ) 
        self.checked_keys.add  ( 'yvar' ) 
        self.checked_keys.add  ( 'zvar' ) 

# =============================================================================
## Helper function to create the function
def make_fun ( fun , args , name ) :
    """Helper function to create the function
    """
    
    assert fun and isinstance ( fun , ROOT.RooAbsReal ), 'make_fun: Invalid type %s' % type ( fun ) 
    
    num = len ( args )
    if   1 == num : return Fun1D ( fun , name = name , *args )
    elif 2 == num : return Fun2D ( fun , name = name , *args )
    elif 3 == num : return Fun3D ( fun , name = name , *args )
    
    raise TypeError ( "Invalid length of arguments %s " % num ) 


# =============================================================================
## Get the derivative  dF/dx for the 1D-function
#  \f[ f(x) = \frac{dF(x)}{dx}\f]
#  @code
#  F    = ...
#  dFdx = F.dFdX ( ) 
#  @endcode 
#  @see RooAsbReal::derivative 
def _f1_deriv_x_ ( self , *args ) :
    """Get the derivative dF/dx for 1D-fuction
    >>> F    = ...
    >>> dFdx = F.dFdX ( ) 
    - see ROOT.RooAbsReal.derivative    
    """
    if not self.dfdx :
        d = self.fun.derivative ( self.xvar , 1 , *args )
        self.dfdx = Fun1D ( d , self.xvar )
    ##
    return self.dfdx 

# =============================================================================
## Get the partial derivative  dF/dx for 2D-function
#  \f[ f(x,y) = \frac{\partial F(x,y)}{\partial x} \f]
#  @code
#  F    = ...
#  dFdx = F.dFdX ( ) 
#  @endcode 
#  @see RooAbsReal::derivative 
def _f2_deriv_x_ ( self , *args ) :
    """Get the partial derivative dF/dx for 2D-function
    >>> F    = ...
    >>> dFdx = F.dFdX ( ) 
    = see ROOT.RooAbsReal.derivative
    """
    if not self.dfdx :
        d = self.fun.derivative ( self.xvar , 1 , *args )
        self.dfdx = Fun2D ( d , self.xvar , self.yvar )
    ##
    return func.dfdx 


# =============================================================================
## Get the partial derivative  dF/dy for 2D-function
#  \f[ f(f,y) = \frac{\partial F(x,y)}{\partial y} \f]
#  @code
#  F    = ...
#  dFdy = F.dFdY ( ) 
#  @endcode 
#  @see RooAbsReal::derivative 
def _f2_deriv_y_ ( self , *args ) :
    """Get the partial derivative dF/dx for 2D-function
    >>> F    = ...
    >>> dFdy = F.dFdY ( ) 
    = see ROOT.RooAbsReal.derivative
    """
    if not self.dfdy :
        d = self.fun.derivative ( self.yvar , 1 , *args )
        self.dfdy = Fun2D ( d , self.xvar , self.yvar )
    ##
    return self.dfdy 

# =============================================================================
## Get the partial derivative  dF/dx for 3D-function
#  \f[ f(x,y,z) = \frac{\partial F(x,y,z)}{\partial x}\f]
#  @code
#  F    = ...
#  dFdx = F.dFdX ( ) 
#  @endcode 
#  @see RooAbsReal::derivative 
def _f3_deriv_x_ ( self , *args ) :
    """Get the partial derivative dF/dx for 3D-function
    >>> F    = ...
    >>> dFdx = F.dFdX ( ) 
    = see ROOT.RooAbsReal.derivative
    """
    if not self.dfdx :
        d = self.fun.derivative ( self.xvar , 1 , *args )
        self.dfdx = Fun3D ( d , self.xvar , self.yvar , self.yvar )
    ##
    return self.dfdx 

# =============================================================================
## Get the partial derivative  dF/dy for 3D-function
#  \f[ f(x,y,z) = \frac{\partial F(x,y,z)}{\partial y} \f]
#  @code
#  F    = ...
#  dFdy = F.dFdY ( ) 
#  @endcode 
#  @see RooAbsReal::derivative 
def _f3_deriv_y_ ( self , *args ) :
    """Get the partial derivative dF/dy for 3D-function
    >>> F    = ...
    >>> dFdy = F.dFdY ( ) 
    = see ROOT.RooAbsReal.derivative
    """
    if not self.dfdy :
        d = self.fun.derivative ( self.yvar , 1 , *args )
        self.dfdy = Fun3D ( d , self.xvar , self.yvar , self.yvar )
    ##
    return self.dfdy 

# =============================================================================
## Get the partial derivative  df/dz for 3D-function
#  \f[ f(x,y,z) = \frac{\partial F(x,y,z)}{\partial z}\f]
#  @code
#  F    = ...
#  dFdz = F.dFdZ ( ) 
#  @endcode 
#  @see RooAbsReal::derivative 
def _f3_deriv_z_ ( self , *args ) :
    """Get the partial derivative dF/dy for 3D-function
    >>> F    = ...
    >>> dFdz = F.dFdZ ( ) 
    = see ROOT.RooAbsReal.derivative
    """
    if not self.dfdz :
        d = self.fun.derivative ( self.zvar , 1 , *args )
        self.dfdz = Fun3D ( d , self.xvar , self.yvar , self.yvar )
    ##
    return self.dfdz 


FUNC .dFdX = _f1_deriv_x_
FUNC2.dFdX = _f2_deriv_x_
FUNC2.dFdY = _f2_deriv_y_
FUNC3.dFdX = _f3_deriv_x_
FUNC3.dFdY = _f3_deriv_y_
FUNC3.dFdZ = _f3_deriv_z_

# =============================================================================
## @var num_bins
#  default number of cache-bins for RooNumRuningInt
#  @see RooNumRuningInt
num_bins = 5000

# =============================================================================
## check cache binnig scheme
def check_bins ( var , *args ) :
    """Check cache bining scheme
    """
    
    if not args :
        if not var.hasBinning('cache') : var.setBins ( 5000 , 'cache' ) 
        return args

    if ininstance ( bins , string_types  ) and var.hasBinning ( bins ) : return bins ,
    if isinstance ( bins , integer_types ) and 100 < bins :
        if var.hasBininig ( 'cache' ) :
            cbins = var.getBinning ( 'cache' ).numBins()
            if cbins < bins :  var.setBins ( bins , 'cache' )
        else : var.setBins ( bins , 'cache' )
        return ()

    raise TypeError ("Invalid binning %s/%s" % ( bins , type ( bins ) ) ) 

        
    
# =============================================================================
## Get the running integral for 1D-function
#  \f[ f(x) = \int_{x_{low}}^{x} F(t) dt \f]
#  @code
#  f = ...
#  g = f.integral_x ( ) 
#  @endcode 
#  @see RooNumRunnigInt 
def _f1_rint_x_ ( self , *args  ) :
    """ Get the running integral for 1D-function
    >>> f = ...
    >>> g = f.integral_x ( ) 
    - see ROOT.RooNumRunnigInt 
    """
    if not self.intx :
        
        kargs = check_bins ( self.xvar , *args )
        
        i = ROOT,RooNumRunnigInt ( "IntX_%s"    % self.name ,
                                   "IntegralX " + self.name ,
                                   self.fun     , self.xvar , *kargs )
        
        self.intx = Fun1D ( i , self.xvar )
    ##
    return func.intx 

# =============================================================================
## Get the running integral for 2D-function
#  \f[ f(x,y) = \int_{x_{low}}^{x} F(t,y) dt \f]
#  @code
#  f = ...
#  g = f.integral_x ( ) 
#  @endcode 
#  @see RooNumRunnigInt 
def _f2_rint_x_ ( self , *args ) :
    """ Get the running integral for 2D-function
    >>> f = ...
    >>> g = f.integral_x ( ) 
    - see ROOT.RooNumRunnigInt 
    """
    if not self.intx :
        
        kargs = check_bins ( self.xvar , *args )
        
        i = ROOT,RooNumRunnigInt ( "IntX_%s"    % self.name ,
                                   "IntegralX " + self.name ,
                                   self.fun     , self.xvar , *kargs )
        
        self.intx = Fun2D ( i , self.xvar , self.yvar )
    ##
    return func.intx 

# =============================================================================
## Get the running integral for 2D-function
#  \f[ f(x,y) = \int_{y_{low}}^{y} F(x,t) dt \f]
#  @code
#  f = ...
#  g = f.integral_y ( ) 
#  @endcode 
#  @see RooNumRunnigInt 
def _f2_rint_y_ ( self , *args  ) :
    """ Get the running integral for 2D-function
    >>> f = ...
    >>> g = f.integral_y ( ) 
    - see ROOT.RooNumRunnigInt 
    """
    if not self.inty :
        
        kargs = check_bins ( self.yvar , *args )

        i = ROOT,RooNumRunnigInt ( "IntY_%s"    % self.name ,
                                   "IntegralY " + self.name ,
                                   self.fun     , self.yvar , *kargs )
        
        self.inty = Fun2D ( i , self.xvar , self.yvar )
    ##
    return func.inty 

# =============================================================================
## Get the running integral for 3D-function
#  \f[ f(x,y,z) = \int_{x_{low}}^{x} F(t,y,z) dt \f]
#  @code
#  f = ...
#  g = f.integral_x ( ) 
#  @endcode 
#  @see RooNumRunnigInt 
def _f3_rint_x_ ( self , *args ) :
    """ Get the running integral for 3D-function
    >>> f = ...
    >>> g = f.integral_x ( ) 
    - see ROOT.RooNumRunnigInt 
    """
    if not self.intx :
        
        kargs = check_bins ( self.xvar , *args )

        i = ROOT,RooNumRunnigInt ( "IntX_%s"    % self.name ,
                                   "IntegralX " + self.name ,
                                   self.fun     , self.xvar , *kargs )
        
        self.intx = Fun3D ( i , self.xvar , self.yvar , self.zvar )
    ##
    return func.intx 

# =============================================================================
## Get the running integral for 3D-function
#  \f[ f(x,y,z) = \int_{y_{low}}^{y} F(x,t,z) dt \f]
#  @code
#  f = ...
#  g = f.integral_y ( ) 
#  @endcode 
#  @see RooNumRunnigInt 
def _f3_rint_y_ ( self , *args  ) :
    """ Get the running integral for 2D-function
    >>> f = ...
    >>> g = f.integral_y ( ) 
    - see ROOT.RooNumRunnigInt 
    """
    if not self.inty :

        kargs = check_bins ( self.yvar , *args )

        i = ROOT,RooNumRunnigInt ( "IntY_%s"    % self.name ,
                                   "IntegralY " + self.name ,
                                   self.fun     , self.yvar , kargs )
        self.inty = Fun3D ( i , self.xvar , self.yvar , self.zvar )
    ##
    return func.inty 


# =============================================================================
## Get the running integral for 3D-function
#  \f[ f(x,y,z) = \int_{z_{low}}^{z} F(x,y,t) dt \f]
#  @code
#  f = ...
#  g = f.integral_z ( ) 
#  @endcode 
#  @see RooNumRunnigInt 
def _f3_rint_z_ ( self , *args ) :
    """ Get the running integral for 3D-function
    >>> f = ...
    >>> g = f.integral_z ( ) 
    - see ROOT.RooNumRunnigInt 
    """
    if not self.intz :
        
        kargs = check_bins ( self.zvar , *args )

        i = ROOT,RooNumRunnigInt ( "IntZ_%s"    % self.name ,
                                   "IntegralZ " + self.name ,
                                   self.fun     , self.zvar , *kargs )
        
        self.intz = Fun3D ( i , self.xvar , self.yvar , self.zvar )
    ##
    return func.intz 


FUNC .integral_x = _f1_rint_x_
FUNC2.integral_x = _f2_rint_x_
FUNC2.integral_y = _f2_rint_y_
FUNC3.integral_x = _f3_rint_x_
FUNC3.integral_y = _f3_rint_y_
FUNC3.integral_y = _f3_rint_z_


# ==============================================================================
## Integrate 2D function over x
#  \f[ f(y) = \int F(x,y) dx \f]
#  @code
#  f = ...
#  g = f.integrate_x ( 'x-range' ) 
#  @endcode 
#  @see RooAbsReal::createIntegral
def _f2_int_x_ ( self , *args ) :
    """ Integrate 2D-function over x  
    >>> f = ...
    >>> g = f.integrate_x ( 'x-range' ) 
    - see ROOT.RooAbsReal.createIntegral
    """
    ##
    vset = ROOT.RooArgSet ( self.xvar ) 
    i = self.fun.createIntegral ( vset , *args )
    ## 
    return Fun1D ( i , self.yvar )

# ==============================================================================
## Integrate 2D function over y
#  \f[ f(x) = \int F(x,y) dy \f]
#  @code
#  f = ...
#  g = f.integrate_y ( 'y-range' ) 
#  @endcode 
#  @see RooAbsReal::createIntegral
def _f2_int_y_ ( self , *args ) :
    """ Integrate 2D-function over y  
    >>> f = ...
    >>> g = f.integrate_x ( 'x-range' ) 
    - see ROOT.RooAbsReal.createIntegral
    """
    ##
    vset = ROOT.RooArgSet ( self.yvar ) 
    i = self.fun.createIntegral ( vset , *args )
    ## 
    return Fun1D ( i , self.xvar )


# ==============================================================================
## Integrate 3D function over x
#  \f[ f(y,z) = \int F(x,y,z) dx \f]
#  @code
#  f = ...
#  g = f.integrate_x ( 'x-range' ) 
#  @endcode 
#  @see RooAbsReal::createIntegral
def _f3_int_x_ ( self , *args ) :
    """ Integrate 3D-function over x  
    >>> f = ...
    >>> g = f.integrate_x ( 'x-range' ) 
    - see ROOT.RooAbsReal.createIntegral
    """
    ##
    vset = ROOT.RooArgSet ( self.xvar ) 
    i = self.fun.createIntegral ( vset , *args )
    ## 
    return Fun2D ( i , self.yvar , self.zvar )

# ==============================================================================
## Integrate 3D function over y
#  \f[ f(x,z) = \int F(x,y,z) dy \f]
#  @code
#  f = ...
#  g = f.integrate_y ( 'y-range' ) 
#  @endcode 
#  @see RooAbsReal::createIntegral
def _f3_int_y_ ( self , *args ) :
    """ Integrate 3D-function over y  
    >>> f = ...
    >>> g = f.integrate_y ( 'y-range' ) 
    - see ROOT.RooAbsReal.createIntegral
    """
    ##
    vset = ROOT.RooArgSet ( self.yvar ) 
    i = self.fun.createIntegral ( vset , *args )
    ## 
    return Fun2D ( i , self.xvar , self.zvar )

# ==============================================================================
## Integrate 3D function over z
#  \f[ f(x,y) = \int F(x,y,z) dz \f]
#  @code
#  f = ...
#  g = f.integrate_z ( 'z-range' ) 
#  @endcode 
#  @see RooAbsReal::createIntegral
def _f3_int_z_ ( self , *args ) :
    """ Integrate 3D-function over z  
    >>> f = ...
    >>> g = f.integrate_z ( 'z-range' ) 
    - see ROOT.RooAbsReal.createIntegral
    """
    ##
    vset = ROOT.RooArgSet ( self.zvar ) 
    i = self.fun.createIntegral ( vset , *args )
    ## 
    return Fun2D ( i , self.xvar , self.yvar )


FUNC2.integrate_x = _f2_int_x_
FUNC2.integrate_y = _f2_int_y_
FUNC3.integrate_x = _f3_int_x_
FUNC3.integrate_y = _f3_int_y_
FUNC3.integrate_z = _f3_int_z_

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
        result = Fun1D ( s            ,
                         xvar = xvar  ,
                         name = fname % ( fun2.name , fun1.name ) )
        
        result.aux_keep.append ( fun1 )
        result.aux_keep.append ( fun2 )
        return result

    elif isinstance ( fun2 , ROOT.RooAbsReal ) :
        
        s      = ftype ( fun2 , fun1.fun )
        result = Fun1D ( s            ,
                         xvar = xvar  ,
                         name = fname % ( fun2.name , fun1.name ) )
        
        result.aux_keep.append ( fun1 )
        result.aux_keep.append ( fun2 )
        return result
    
    return NotImplemented 


# =============================================================================
## operator for "other + 1D-function"
def _f1_radd_ ( self , other ) :
    """Operator for ``other + 1D-function''"""
    return _f1_rop_ ( self , other , Ostap.MoreRooFit.Addition    , "Add_%s_%s"     )

# =============================================================================
## operator for "other - 1D-function"
def _f1_rsub_ ( self , other ) :
    """Operator for ``other - 1D-function''"""
    return _f1_rop_ ( self , other , Ostap.MoreRooFit.Subtraction , "Subtract_%s_%s" )

# =============================================================================
## operator for "other * 1D-function"
def _f1_rmul_ ( self , other ) :
    """Operator for ``other * 1D-function''"""
    return _f1_rop_ ( self , other , Ostap.MoreRooFit.Product     , "Product_%s_%s"  )

# =============================================================================
## operator for "other/1D-function"
def _f1_rdiv_ ( self , other ) :
    """Operator for ``other / 1D-function ''"""
    return _f1_rop_ ( self , other , Ostap.MoreRooFit.Division    , "Divide_%s_%s"  )
        
# =============================================================================
## operator for "other ** 1D-function"
def _f1_rpow_ ( self , other ) :
    """Operator for ``other ** 1D-function''"""
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
