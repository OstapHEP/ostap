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
    'AFUN1'             , ## base   class for 1D-function/PDF-like objects 
    'AFUN2'             , ## base   class for 2D-function/PDF-like objects 
    'AFUN3'             , ## base   class for 3D-function/PDF-like objects 
    'F1AUX'             , ## MIXIN  class for 1D-function-like objects
    'FUN1'              , ## base   class for 1D-function-like objects
    'FUN2'              , ## base   class for 2D-function-like objects
    'FUN3'              , ## base   class for 3D-function-like objects
    'Fun1D'             , ## wrapper for 1D-function
    'Fun2D'             , ## wrapper for 2D-function
    'Fun3D'             , ## wrapper for 3D-function
    'Id'                , ## the most trivial function 
    'make_fun1'         , ## reate popular FUN1 objects from the short description 
    )
# =============================================================================
from   sys                           import version_info as python_version 
from   ostap.core.ostap_types        import ( integer_types  , num_types    ,
                                              dictlike_types , list_types   ,
                                              is_good_number , string_types )
from   ostap.core.core               import ( Ostap         ,
                                              valid_pointer ,
                                              roo_silent    ,
                                              rootWarning   )   
from   ostap.math.base               import iszero , isequal  
from   ostap.fitting.variables       import SETVAR
from   ostap.fitting.roofit          import PDF_fun 
from   ostap.fitting.fithelpers      import ( VarMaker , FitHelper ,
                                              XVar, YVar, ZVar, NameDuplicates ) 
from   ostap.utils.cidict            import cidict
from   ostap.plotting.fit_draw       import key_transform, draw_options  
import ostap.fitting.variables
import ostap.fitting.roocollections
import ROOT, math, sys, abc  
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.funbasic' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================
py2 = 2 >= sys.version_info.major
constant_types = num_types + ( ROOT.RooConstVar , ) 
# =============================================================================
## is valuer equal to 1?
isone = lambda x : isequal ( float ( x ) , 1 )
# =============================================================================
## helper factory function
def func_factory ( klass , config ) :
    """Helper factory function, used for unpickling"""
    with NameDuplicates ( True ) : 
        return klass ( **config ) 
# =============================================================================
## @class AFUN1
#  Helper base class for implementation of various (Roo)Function-wrappers
class AFUN1(XVar,FitHelper) : ## VarMaker) :
    """Helper base class for implementation of various (Roo)Function-wrappers
    """
        
    def __init__ ( self , name , xvar , tricks = True , **kwargs ) :

        ## name is defined via the base class VarMaker 
        self.name  = name ## name is defined via the base class VarMaker 

        XVar .__init__ ( self , xvar )

        self.__vars      = ROOT.RooArgSet  ()
        self.__variables = [] 

        self.__config       = {}
        self.__fun          = None
        
        self.__tricks       = True if tricks else False 
        
        self.__draw_var     = None
        ## predefined drawing options for this FUN/PDF
        self.__draw_options = cidict ( transform = key_transform )
        
        self.__checked_keys = set()
                
        self.vars.add         ( self.xvar )
        self.variables.append ( self.xvar )

        self.config = { 'name' : self.name , 'xvar' : self.xvar  }
        self.__func_init = True  

        ## derived functions/objects
        self.__derived = {}

        ## decode the keyword arguments 
        dropts = draw_options ( **kwargs )
        self.__draw_options.update ( dropts )

        ## check for extra arguments 
        extra  = {}
        for k in kwargs :
            if not k in dropts : extra [ k ] = kwargs [ k ] 
        if extra : self.error ("Unknown arguments %s" % extra )

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
    def title ( self ) :
        """'title' : the title for RooAbsReal/RooAbsPdf"""
        return self.fun.title if self.fun else self.name

    @property 
    def vars ( self ) :
        """'vars' : variables/observables (as ROOT.RooArgSet)"""
        return self.__vars

    @property 
    def variables ( self ) :
        """'variables' : variables/observables at list"""
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
        assert isinstance ( value , ROOT.RooAbsReal ) , "'pdf/fun' is not ROOT.RooAbsReal"
        self.__fun = value

    @property
    def fun_name ( self ) :
        """'fun_name' : the name of the underlying RooAbsReal"""
        return  self.fun.GetName() if self.fun else ''

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
    def tricks ( self ) :
        """'tricks' : flag to allow some  tricks&shortcuts """
        return self.__tricks
    @tricks.setter
    def tricks ( self , value ) :
        val = True if value else False 
        if val and not self.__tricks :
            raise ValueError("Can't allow tricks&shortcuts!")
        self.__tricks = val

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
        if isinstance ( var , AFUN1 ) :
            if self.fun.depends_on ( var.fun  ) : return True
            if self.fun.depends_on ( var.vars ) : return True
                
        return False 
            
    @property
    def draw_var ( self ) :
        """'draw_var'  :  variable to be drawn if not specified explicitely"""
        return self.__draw_var
    @draw_var.setter
    def draw_var ( self , value ) :
        assert value is None or isinstance ( value , ROOT.RooAbsReal ) , \
               "'draw_var' has invalid type %s/%s" % ( value , type(value) )
        self.__draw_var = value 

    @property
    def checked_keys ( self ) :
        """'checked keys' : special keys for clone-method
        """
        return self.__checked_keys

    # =========================================================================
    ## simple 'function-like' interface
    @abc.abstractmethod 
    def __call__ ( self , *args , **kwargs ) :
        """ Function as a 'function'
        """
        return None 
    
    # =========================================================================
    @property
    def draw_options ( self ) :
        """'draw_options' : dictionary with predefined draw-options for this PDF/FUN
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
        assert self.fun, "AFUN1::params: Function is invalid!"        
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
    #  pdf.load_params ( params, dataset ) 
    #  params  = ( A , B , C , ... )
    #  pdf.load_params ( params , dataset )  
    #  @endcode 
    def load_params ( self , params = {} , dataset = None , silent = False , **kwargs ) :
        """Load parameters from
        - external dictionary `{ name : value }`
        - sequence of `RooAbsReal` objects
        - `RooFitResult` object
        
        >>> pdf      = ...
        >>> dataset = ... 
        >>> params = { 'A' : 10 , 'B' : ... }
        >>> pdf.load_params ( params , dataset  ) 
        >>> params = ( A , B , C , ... )
        >>> pdf.load_params ( params , dataset )  
        """

        if dataset :
            assert     isinstance ( dataset , ROOT.RooAbsData ) , "load_params: invalid type of 'dataset:%s'" % type ( dataset ) 
        else :
            dataset = ROOT.nullptr
            
        if params :
            assert not isinstance ( params  , ROOT.RooAbsData ) , "load_params: invalid type of 'params':%s" % type ( params ) 
            
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
                        item = p.name , "%-15.7g" % pv , "%-+15.7g" % vv 
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
                        item = p.name , "%-15.7g" % pv , "%-+15.7g" % vv 
                        table.append ( item ) 
                    keys.add  ( i )

            not_used = []
            for i , pp in enumerate ( params ) :  
                if i in keys : continue
                not_used.append ( pp )

        ## explicit parameters 
        keys = set()        
        for key in kwargs :
            for p in pars :
                if not hasattr ( p  , 'setVal' ) : continue
                if p.name != key                 : continue
                
                v  = kwargs [key]
                vv = float ( v  )
                pv = p.getVal ()   
                if vv != pv : 
                    p.setVal   ( vv )
                    item = p.name , "%-15.7g" % pv , "%-+15.7g" % vv 
                    table.append ( item ) 
                keys.add ( key )
                
            not_used |= set ( kwargs.keys() ) - keys 

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

        key_transform = self.draw_options.transform

        the_key = key_transform ( key ) 

        ## 1. check the explicitely provided arguments
        for k in kwargs :
            if key_transform ( k ) == the_key :
                return kwargs [ k ]
            
        ## check the predefined drawing options for this PDF
        if key in self.draw_options :
            return self.draw_options.get ( key )
            
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
            self.warning("Neither 'options' nor 'style'...")

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
                    raise AttributeError ("Class %s cannot be cloned with '%s' key" % ( self.__class__.__name__ , k ) )

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
            self.plot_on ( self.fun , frame , *coptions )
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

    # =========================================================================
    ## invoke <cdoe>what.plotOn (frame , *options)</code> command
    #  - merge arguments using <code>RooFit::MultiArg</code> to shorted list
    def plot_on ( self , what , frame , *options ) :
        """Invoke `what.plotOn (frame , *options)` command
        - merge arguments using `ROOT.RooFit::MultiArg` to shorted list
        """
        
        NARGS = 8

        assert all ( isinstance ( o , ROOT.RooCmdArg ) for o in options  ), \
               "plot_on: invalid argument types: %s" % list ( options  ) 

        ## for "small' number of arguments use the standard function 
        if len ( options ) <= NARGS :
            return what.plotOn ( frame  , *options )
        
        from ostap.fitting.roocmdarg import command 
        cmd = command ( *options )

        return what.plotOn ( frame , cmd  )

        ## ## merge arguments to get shorter list        

        ## head = options [            : NARGS - 1 ]
        ## tail = options [ NARGS - 1  :           ]
        
        ## from   ostap.utils.utils  import chunked
        ## if 1 == len ( tail ) % NARGS  : chunks = chunked ( tail , NARGS - 1  )
        ## else                          : chunks = chunked ( tail , NARGS      )
        
        ## new_options = head + tuple ( ROOT.RooFit.MultiArg ( *chunk ) for chunk in chunks )

        ## self.debug ( 'plotOn: merged options: %s' % str ( new_options ) ) 
        ## return self.plotOn ( what , frame , *new_options ) 

# =============================================================================
## @class F1AUX
#  Helper MIXIN class for impleementation of various 1D (Roo)Function-wrappers
# 
#  It relies on following atttribute2:
#  - <code>fun</cpode>
#  - <code>tricks</cpode>
#
#  It relies on follownig methods
#  - <code>__call__</code>
#  - <code>xvar</code>
#  - <code>error</code>
#  - <code>xminmax</code>
class F1AUX(object) :
    """Helper MIXIN class for implementation of various 1D (Roo)Function-wrappers
    
    It relies on following atttribute2:
    - `fun`
    - `tricks`
    
    It relies on following methods
    - `__call__`
    - `xvar`
    - `error`
    - `xminmax`
    """
    
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
            ff = PDF_fun  ( fun , self.xvar , xmin = xmin , xmax = xmax )

        return funcall ( ff , xmin = xmin , xmax = xmax , *args , **kwargs )

    # =========================================================================
    ## get the effective mean
    def get_mean ( self , **kwargs ) :
        """Get the effective Mean
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print('MEAN: %s ' % pdf.get_mean() )
        """
        from ostap.stats.moments import mean as _mean
        return self._get_stat_ ( _mean , **kwargs )

    # ========================================================================
    ## get the effective RMS 
    def rms ( self , **kwargs ) :
        """Get the effective RMS
        >>>  fun = ...
        >>>  print('RMS: %s ' % fun.rms())
        """
        
        from ostap.stats.moments import rms      as sp_rms
        from ostap.stats.moments import variance as sp_variance 
        
        fun   = self.fun
        ftype = type ( fun ) 
        if   hasattr ( ftype , 'rms' ) and not ftype.rms is sp_rms :
            return fun.rms ()        
        elif hasattr ( ftype , 'Rms' )                             :
            return fun.Rms ()        
        elif hasattr ( ftype , 'RMS' )                             :
            return fun.RMS ()
        
        elif self.tricks and hasattr ( fun , 'function' ) :

            if   hasattr ( fun , 'setPars'    ) : fun.setPars()
            
            ff = fun.function()
            ftype = type ( ff ) 
            if   hasattr ( ftype , 'rms'        ) and not ftype.rms        is sp_rms      :
                return ff.rms()
            elif hasattr ( ftype , 'variance'   ) and not ftype.variance   is sp_variance :
                return ff.variance   () ** 0.5  
            elif hasattr ( ftype , 'dispersion' ) and not ftype.dispersion is sp_variance :
                return ff.dispersion ()**0.5 

        return  self._get_stat_ ( sp_rms , **kwargs )

    # ========================================================================
    ## get the effective Full Width at Half Maximum
    def fwhm ( self , **kwargs ) :
        """Get the effective Full Width at  Half Maximum
        >>>  fun = ...
        >>>  print('FWHM: %s ' % fun.fwhm())
        """
        ## use generic machinery 
        from ostap.stats.moments import width as _width
        w = self._get_stat_ ( _width , **kwargs )
        return  w [ 1 ] - w [ 0 ]

    # ========================================================================
    ## get the skewness 
    def skewness ( self , **kwargs ) :
        """Get the skewness 
        >>>  fun = ...
        >>>  print('Skewness: %s ' % fun.skewness())
        """
        ## use generic machinery 
        from ostap.stats.moments import skewness as _calc
        return self._get_stat_ ( _calc , **kwargs )

    # ========================================================================
    ## get the kurtosis
    def kurtosis ( self , **kwargs ) :
        """Get the kurtosis
        >>>  fun = ...
        >>>  print('Kurtosis: %s ' % fun.kurtosis())
        """
        ## use generic machinery 
        from ostap.stats.moments import kurtosis as _calc
        return self._get_stat_ ( _calc , **kwargs )

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
        >>>  print('MODE: %s ' % pdf.mode())
        """
        from ostap.stats.moments import mode as _calc
        return self._get_stat_ ( _calc , **kwargs )

    # =========================================================================
    ## get the mediane 
    def median ( self , **kwargs ) :
        """Get the effective mode
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print ( 'MODE: %s ' % pdf.mode() )
        """
        from ostap.stats.moments import median as _calc
        return self._get_stat_ ( _calc , **kwargs )

    # =========================================================================
    ## get the effective moment for the distribution
    def moment ( self , N , **kwargs ) :
        """Get the effective moment
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print ( 'MOMENT: %s ' % pdf.moment( 10 ) )
        """
        ## use generic machinery 
        from ostap.stats.moments import moment as _moment
        return self._get_stat_ ( _moment , N , **kwargs ) 
    
    # =========================================================================
    ## get the effective central moment for the distribution
    def central_moment ( self , N , **kwargs ) :
        """Get the effective central moment
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print ( 'CENTRAL MOMENT: %s ' % pdf.central_moment( 10 ) )
        """
        from ostap.stats.moments import central_moment as _moment
        return self._get_stat_ ( _moment , N , **kwargs ) 

    # =========================================================================
    ## get the standartized moment for the distribution
    def std_moment ( self , N , **kwargs ) :
        """Get the standartized moment
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print('STD-MOMENT: %s ' % pdf.std_moment( 10 ))
        """
        from ostap.stats.moments import std_moment as _moment
        return self._get_stat_ ( _moment , N , **kwargs ) 

    # =========================================================================
    ## get the effective quantile 
    def quantile ( self , prob  , **kwargs ) :
        """Get the effective quantile
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print ( 'QUANTILE: %s ' % pdf.quantile ( 0.10 ) )
        """
        from ostap.stats.moments import quantile as _quantile
        return self._get_stat_ ( _quantile , prob , **kwargs ) 

    # =========================================================================
    ## get the symmetric confidence interval 
    def cl_symm ( self , prob , x0 =  None , **kwargs ) :
        """Get the symmetric confidence interval 
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print ( 'CL :  ',  pdf.cl_symm ( 0.10 ) ) 
        """
        from ostap.stats.moments import cl_symm as _cl
        return self._get_stat_ ( _cl , prob , x0 , **kwargs ) 

    # =========================================================================
    ## get the asymmetric confidence interval 
    def cl_asymm ( self , prob , **kwargs ) :
        """Get the asymmetric confidence interval 
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print ( 'CL :  ',  pdf.cl_asymm ( 0.10 ) )
        """
        from ostap.stats.moments import cl_asymm as _cl
        return self._get_stat_ ( _cl , prob , **kwargs )

    # =========================================================================
    ## get moment using <code>RooAbsPdf::moment</code> method
    #  @see RooAbsPdf::moment
    #  @code
    #  pdf = ...
    #  v4  = pdf.roo_moment ( 5 , central = True )
    #  @endcode
    def roo_moment ( self , order , central ) :
        """Get moment using <code>RooAbsPdf::moment</code> method
        >>> pdf = ...
        >>> v5  = pdf.roo_moment ( 5 , central = True )
        - see `ROOT.RooAbsPdf.moment`
        """
        assert isinstance ( order , integer_types ) and 0<= order , \
               'roo_moment: invalid moment order %s !' % order

        if   central and 0 == order : return 1
        elif central and 1 == order : return 0
        
        with rootWarning() , roo_silent ( True ) : 
            mom    = self.fun.moment ( self.xvar , order , central , False  )
            result = mom.getVal ()
            mom.Delete() 
            del mom
            return result 

    # =========================================================================
    ## get mean using RooAbsPdf method
    #  @see RooAbsReal::mean 
    def roo_mean ( self ) :
        """get mean using RooAbsPdf method
         -see `ROOT.RooAbsReal.mean`
         """
        with rootWarning() , roo_silent ( True ) : 
            mom    = self.fun.mean ( self.xvar )
            result = mom.getVal ()
            mom.Delete() 
            del mom
            return result 

    # =========================================================================
    ## get variance using RooAbdPdf method
    #  @see RooAbdPdf::moment 
    def roo_variance ( self ) :
        """get variance  using the RooAbdReal method
         -see `ROOT.RooAbdPdf.moment`
         """
        return self.roo_moment ( 2 , central = True )

    # =========================================================================
    ## get dispersion  using RooAbdPdf method
    #  @see RooAbdPdf::moment 
    def roo_dispersion  ( self ) :
        """get dispersion using RooAbdPdf method
         -see `ROOT.RooAbdPdf.moment`
         """
        return self.roo_variance () 

    # =========================================================================
    ## get RMS using RooAbdPdf method
    #  @see RooAbdPdf::sigma 
    def roo_rms ( self ) :
        """get RMS using RooAbsPdf method
         -see `ROOT.RooAbdPdf.sigma`
         """
        with rootWarning() , roo_silent ( True ) : 
            mom    = self.fun.sigma ( self.xvar )
            result = mom.getVal ()
            mom.Delete() 
            del mom
            return result 


    # =========================================================================
    ## get skewness using RooAbdPdf method
    #   \f$  k = \frac{\mu_3}{\mu_2^{3/2}}\f$ 
    #  @see RooAbdPdf::moment 
    def roo_skewness ( self ) :
        """get skewness using RooAbdPdf method
        -see `ROOT.RooAbdPdf.moment`
        """
        m2 = self.roo_moment ( 2 , central = True )
        m3 = self.roo_moment ( 3 , central = True )        
        return m3 / ( m2 ** ( 3.0 / 2 ) )

    # =========================================================================
    ## get (excessive) kurtosis  using RooAbdPdf method
    #   \f$  k = \frac{\mu_4}{\mu_2^2} -3 \f$ 
    #  @see RooAbdPdf::moment 
    def roo_kurtosis ( self ) :
        """get (excessive) kurtosis using RooAbdPdf method
        -see `ROOT.RooAbdPdf.moment`
        """
        m2 = self.roo_moment ( 2 , central = True )
        m4 = self.roo_moment ( 4 , central = True )    
        return m4 / ( m2 * m2 ) - 3.0 
    
    # =========================================================================
    ## get the derivative at  point x 
    def derivative ( self , x ) :
        """Get derivative at point x 
        >>> pdf = ...
        >>> print ( pdf.derivative ( 0 ) )
        """
        ## check limits
        if self.xminmax() :
            xmin , xmax = self.xminmax ()
            if x < xmin or x > xmax : return 0

        ## make a try to use analytical derivatives 
        if self.tricks  and hasattr ( self.fun , 'function' ) :
            if hasattr ( _fun , 'setPars'  ) : _fun.setPars() 
            try: 
                funfun = sefl.fun.function() 
                if hasattr ( funfun , 'derivative' ) :
                    return funfun.derivative ( x )
            except:
                pass
            
        ## use numerical derivatives 
        from ostap.math.derivative import derivative as _derivative
        return _derivative ( self , x )

    # ==========================================================================
    ## get a minimum of FUN for certain interval
    #  @code
    #  fun = ...
    #  x   = fun.minimum() 
    #  @endcode 
    def minimum ( self , xmin = None , xmax = None , x0 = None ) :
        """Get a minimum of FUN for certain interval
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
    ## get a maximum of FUN for certain interval
    #  @code
    #  fun = ...
    #  x   = fun.maximum() 
    #  @endcode 
    def maximum ( self , xmin = None , xmax = None , x0 = None ) :
        """Get a maximum of FUN for certain interval
        >>> fun = ...
        >>> x   = fun.maximum()
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

    # ========================================================================
    ## check minmax of the PDF using the random shoots
    #  @code
    #  pdf     = ....
    #  mn , mx = pdf.minmax()            
    #  @endcode 
    def minmax ( self , nshoots =  50000 ) :
        """Check min/max for the PDF using  random shoots 
        >>> pdf     = ....
        >>> mn , mx = pdf.minmax()        
        """
        ## try to get minmax directly from pdf/function

        fun = self.fun
        
        if self.tricks and hasattr ( fun , 'function' ) :
            if hasattr ( fun , 'setPars' ) : fun.setPars() 
            f = self.fun.function()
            if   hasattr ( f , 'minmax' ) :
                try :
                    mn , mx = f.minmax()
                    if mn <= mx : return mn , mx
                except :
                    pass
            elif hasattr ( f , 'max' ) and  hasattr ( f , 'min' ) : 
                try :
                    mx = f.max()
                    mn = f.min()
                    if mn <= mx : return mn , mx 
                except :
                    pass
            elif hasattr ( f , 'max' ) and not hasattr ( f , 'min' ) : 
                try :
                    mx = f.max()
                    if 0 <= mx : return 0 , mx 
                except :
                    pass                
            elif hasattr ( f , 'mode' ) :
                try :
                    mode = f.mode()
                    if mode in self.xvar :
                        mx = f ( mode ) 
                        if 0 < mx : return 0 , mx
                except :
                    pass 
                    
        ## check RooAbsReal functionality
        if hasattr ( fun , 'getMaxVal' ) : 
            code = fun.getMaxVal( ROOT.RooArgSet ( self.xvar ) )
            if 0 < code :
                mx = self.pdf.maxVal ( code )
                if 0 < mx : return 0 , mx
                
        mn , mx = -1 , -10
        if hasattr ( fun , 'min' ) : mn = fun.min()
        if hasattr ( fun , 'max' ) : mx = fun.max()
        if 0 <= mn and mn <= mx and 0 < mx : return mn , mx
        
        ## now try to use brute force and random shoots 
        if not self.xminmax() : return ()
        
        mn  , mx = -1 , -10
        xmn , xmx = self.xminmax()
        for i in range ( nshoots ) : 
            xx = random.uniform ( xmn , xmx )
            with SETVAR ( self.xvar ) :
                self.xvar.setVal ( xx )
                vv = fun.getVal()
                if mn < 0 or vv < mn : mn = vv
                if mx < 0 or vv > mx : mx = vv 
        return mn , mx 


    # =========================================================================
    ## get the integral between xmin and xmax 
    def integral ( self , xmin , xmax , nevents = True ) :
        """Get integral between xmin and xmax
        >>> pdf = ...
        >>> print ( pdf.integral ( 0 , 10 ) ) 
        """
        ## check limits
        if self.xminmax() :
            mn , mx = self.xminmax() 
            xmin = max ( xmin , mn )  
            xmax = max ( xmax , mn )

        ## initialize the value and the flag 
        value , todo = 0 , True

        fun = self.fun
        
        ## 1) make a try to use 'analytical' integral 
        if self.tricks and hasattr ( fun , 'function' ) : 
            try:
                if hasattr ( fun , 'setPars'  ) : fun.setPars() 
                funfun       = fun.function()
                value , todo = funfun.integral ( xmin , xmax ) , False 
            except:
                pass

        extended = hasattr ( fun , 'canBeExtended' ) and fun.canBeExtended()
        if not extended : extended = isinstance ( fun , ROOT.RooAbsPdf )

        ## 2) use numerical integration
        from ostap.math.integral import integral as _integral
                
        if   todo and extended : value = _integral ( self , xmin , xmax )

        elif todo :
                        
            ## use unormalized PDF here to speed up the integration 
            ifun   = lambda x :  self ( x , normalized = False )
            value  = _integral ( ifun , xmin , xmax )
            
            if hasattr ( fun , 'getNorm' ) : 
                norm   = fun.getNorm ( self.vars )
                value /= norm
                
        if nevents and hasattr ( fun , 'mustBeExtended' ) and fun.mustBeExtended () :
            evts = fun.expectedEvents( self.vars )
            if evts <= 0 or iszero ( evts ) :
                self.warning ( "integral: expectedEvents is %s" % evts )
            value *= evts 

        return value


# =============================================================================
## @class FUN1
#  Helper base class for impleementation of various 1D (Roo)Function-wrappers
class FUN1(AFUN1,F1AUX) :
    """Helper base class for implementation of various 1D (Roo)Function-wrappers
    """
    def __init__ ( self , name , xvar , tricks = True , **kwargs ) :
        
        AFUN1 .__init__ ( self , name = name , xvar = xvar , tricks = tricks , **kwargs )

        ## save the configuration
        self.config = {
            'name'   : self.name   ,
            'xvar'   : self.xvar   ,
            'tricks' : self.tricks , 
            }
        self.config.update ( kwargs )
        
        self.__call_OK = isinstance ( self.xvar , ROOT.RooAbsRealLValue ) 
        
    # =========================================================================
    ## simple 'function-like' interface 
    def __call__ ( self , x , normalized = False ) :
        """ Function as a 'function'
        >>> fun  = ...
        >>> x = 1
        >>> y = fun ( x ) 
        """
        
        assert self.__call_OK , "Invalid type for xvar!"
        
        xmnmx = self.xminmax()
        if xmnmx :
            xmn , xmx = xmnmx 
            if not xmn <= x <= xmx : return 0
            
        with SETVAR( self.xvar ) :
            
            self.xvar.setVal ( x )
            
            return self.fun.getVal ( self.vars ) if normalized else self.fun.getVal ()  
        
    # ========================================================================
    ## convert to float 
    def __float__ ( self ) :
        """Convert to float
        >>> fun = ...
        >>> v   = float ( fun )
        """
        return self.fun.getVal () ## self.vars ) 

    # ==========================================================================
    ## convert PDF into TF1 object, e.g. to profit from TF1::Draw options
    #  @code
    #  fun = ...
    #  tf1 = pdf.tf1()
    #  tf1.Draw('colz')
    #  @endcode
    def tf1 ( self , xmin = None , xmax = None ) :
        """Convert FUN  to TF1 object, e.g. to profit from TF1::Draw options
        >>> fun = ...
        >>> tf1 = fun.tf1()
        >>> tf1.Draw('colz')
        """

        if xmin is None or xmax is None :
            xmnmn = self.xminmax ()
            if xmin is None and xmnmx : xmin = xmnmx[0]
            if xmax is None and xmnmx : xmax = xmnmx[1]

        from ostap.core.core import fID
        def _aux_fun_ ( x , pars = [] ) : return self ( x[0]  )
        self.aux_keep.append ( _aux_fun_ )
        
        return ROOT.TF1 ( fID() , _aux_fun_ , xmin , xmax ) 

    # =========================================================================
    ## operations with FUN1 object: self (op) fun2 
    def __f1_op__ ( self , fun2 , ftype , pattern ) :
        """Operations with FUN1 object: self (op) fun2
        """

        fun_name = lambda  name1 , name2 : \
                   self.generate_name ( pattern % ( name1 , name2 ) )
        
        if isinstance ( fun2 , AFUN3 ) and self.xvar in fun2.vars :
            return Fun3op ( self , fun2 , operation = ftype, 
                            xvar    = fun2.xvar ,
                            yvar    = fun2.yvar ,
                            zvar    = fun2.zvar ,                            
                            name    = fun_name  ( self.name , fun2.name ) ,
                            strname = pattern % ( self.name , fun2.name ) )

        elif isinstance ( fun2 , AFUN3 ) :
            self.warning ( '1D (op) 3D : keep only 3D arguments %s ' % [ a.name for a in fun2.vars ] )        
            return Fun3op ( self , fun2 , operation = ftype, 
                            xvar    = fun2.xvar ,
                            yvar    = fun2.yvar ,
                            zvar    = fun2.zvar ,                            
                            name    = fun_name  ( self.name , fun2.name ) ,
                            strname = pattern % ( self.name , fun2.name ) )
       
        elif isinstance ( fun2 , AFUN2 ) and self.xvar in fun2.vars :        
            return Fun2op ( self , fun2 , operation = ftype, 
                            xvar    = fun2.xvar ,
                            yvar    = fun2.yvar ,
                            name    = fun_name  ( self.name , fun2.name ) ,
                            strname = pattern % ( self.name , fun2.name ) )

        elif isinstance ( fun2 , AFUN2 ) :
            return Fun3op ( self , fun2 , operation = ftype, 
                            xvar    = fun2.xvar ,
                            yvar    = fun2.yvar ,
                            zvar    = self.xvar ,                            
                            name    = fun_name  ( self.name , fun2.name ) ,
                            strname = pattern % ( self.name , fun2.name ) )

        elif isinstance ( fun2 , AFUN1 ) and self.xvar in fun2.vars :            
            return Fun1op ( self , fun2 , operation = ftype , 
                            xvar    = self.xvar ,
                            name    = fun_name  ( self.name , fun2.name ) ,
                            strname = pattern % ( self.name , fun2.name ) )
        
        elif isinstance ( fun2 , AFUN1 ) : 
            return Fun2op ( self , fun2 , operation = ftype , 
                            xvar    = self.xvar ,
                            yvar    = fun2.xvar ,
                            name    = fun_name  ( self.name , fun2.name ) ,
                            strname = pattern % ( self.name , fun2.name ) )
        
        elif  isinstance ( fun2 , constant_types ) :            
            fun2   = ROOT.RooFit.RooConst ( float ( fun2 )  )            
            return Fun1op ( self , fun2 , operation = ftype ,
                            xvar    = self.xvar ,
                            name    = fun_name  ( self.name , fun2.name ) ,
                            strname = pattern % ( self.name , fun2.name ) )
        
        elif  isinstance ( fun2 , ROOT.RooAbsReal ) :
            return Fun1op ( self , fun2 , operation = ftype ,
                            xvar    = self.xvar ,
                            name    = fun_name  ( self.name , fun2.name ) ,
                            strname = pattern % ( self.name , fun2.name ) )
        
        return NotImplemented 

    # =========================================================================
    ## operations with FUN1 object: fun2 (op) self 
    def __f1_rop__ ( self , fun2 , ftype , pattern  ) :
        """operations with FUN1 object: fun2 (op) self 
        """
        
        fun_name = lambda  name1 , name2 : self.generate_name ( pattern % ( name1 , name2 ) )
        
        if  isinstance ( fun2 , constant_types ) :
            fun2 = ROOT.RooFit.RooConst ( float ( fun2 )  )
                
        if  isinstance ( fun2 , ROOT.RooAbsReal ) :            
            return Fun1op ( fun2 , self , operation = ftype ,
                            xvar    = self.xvar  ,                           
                            name    = fun_name  ( fun2.name , self.name ) , 
                            strname = pattern % ( fun2.name , self.name ) )
        
        return NotImplemented 

    # =========================================================================
    ## make a constant function 
    def __constant ( self , value ) :
        """ake a constant function"""
        return Fun1D ( value , xvar = self.xvar , name = self.new_name ( 'const1(%s)' % value ) ) 
    
    # =========================================================================
    ## Add two functions
    #  @code
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 + f2
    #  @endcode 
    def __add__ ( self , other ) :
        """Add two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 + f2 
        """
        if isinstance ( other , constant_types ) and iszero ( other ) : return self
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Addition , pattern = '(%s+%s)' )
    
    # =========================================================================
    ## Subtract  two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 - f2 
    #  @endcode 
    def __sub__ ( self , other ) :
        """Subtract two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 + f2 
        """
        if   self is other : return self.__constant ( 0 )
        elif isinstance ( other , constant_types ) and iszero ( other ) : return self
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Subtraction , pattern = "(%s-%s)" )  
        
    # =========================================================================
    ## Multiply two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 * f2 
    #  @endcode 
    def __mul__ ( self , other ) :
        """Multiply two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 * f2 
        """
        if   isinstance ( other , constant_types ) and isone  ( other ) : return self
        elif isinstance ( other , constant_types ) and iszero ( other ) : return self.__constant ( 0 )
        elif self is other                                                   : return self ** 2 
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Product    , pattern = "(%s*%s)" )  
        
    # =========================================================================
    ## Divide two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 / f2 
    #  @endcode 
    def __div__ ( self , other ) :
        """Divide two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 / f2 
        """
        if self is other                                                     : return self.__constant ( 1 )
        elif isinstance ( other , constant_types ) and isequal ( other , 1 ) : return self
        elif isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 )
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Division    , pattern = "(%s/%s)" )  
    __truediv__ = __div__

    # =========================================================================
    ## Power function for two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 ** f2 
    #  @endcode 
    def __pow__ ( self , other ) :
        """Power function  for two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 ** f2 
        """
        if   isinstance ( other , constant_types ) and isequal ( other , 1 ) : return self
        elif isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 )
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Power    , pattern = 'pow(%s,%s)')  

    # =========================================================================
    ## Add two functions
    #  @code
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 + f2
    #  @endcode 
    def __radd__ ( self , other ) :
        """Add two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 + f2 
        """
        if isinstance ( other , constant_types ) and iszero ( other ) : return self
        return self.__f1_rop__ ( other ,  Ostap.MoreRooFit.Addition    , pattern = "(%s+%s)" )

    # =========================================================================
    ## Subtraction for two functions
    #  @code
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 - f2
    #  @endcode 
    def __rsub__ ( self , other ) :
        """Subtract two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 - f2 
        """
        if   self is other : return self.__constant ( 0 )
        return self.__f1_rop__ ( other ,  Ostap.MoreRooFit.Subtraction , pattern = "(%s-%s)" )

    # =========================================================================
    ## Multiply two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 * f2 
    #  @endcode 
    def __rmul__ ( self , other ) :
        """Multiply two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 * f2 
        """
        if   isinstance ( other , constant_types ) and isequal ( other , 1 ) : return self
        elif isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 )
        return self.__f1_rop__ ( other ,  Ostap.MoreRooFit.Product    , pattern = "(%s*%s)" )  
        
    # =========================================================================
    ## Divide two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 / f2 
    #  @endcode 
    def __rdiv__ ( self , other ) :
        """Divide two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 / f2 
        """
        if   self is other                                                   : return self.__constant ( 1 )
        elif isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 )
        return self.__f1_rop__ ( other ,  Ostap.MoreRooFit.Division    , pattern = "(%s/%s)" )  
    __rtruediv__ = __rdiv__

    # =========================================================================
    ## Power function for two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 ** f2 
    #  @endcode 
    def __rpow__ ( self , other ) :
        """Power function  for two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 ** f2 
        """
        if   isinstance ( other , constant_types ) and isequal ( other , 1 ) : return self.__constant ( 1 )
        elif isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 ) 
        return self.__f1_rop__ ( other ,  Ostap.MoreRooFit.Power    , pattern = "pow(%s,%s)"  )  

    # =========================================================================
    ## absolute value : \f$ \left| bf(x) \right| \f$
    #  @code
    #  f1 = ...
    #  ff = abs ( f ) 
    #  @endcode
    def __abs__ ( self , other = 1 ) :
        """Absolute value :  |f(x)|
        >>> f1 = ...
        >>> ff = abs ( f ) 
        """
        if   isinstance ( other    , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 )
        elif isinstance ( self.fun , ROOT.RooAbsPdf ) and \
             isinstance ( other    , constant_types ) and isequal ( other , 1 ) : return self
        
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Abs , pattern = "abs(%s*%s)" )  
        
    # ==============================================================================
    ## Exponent \f$ f = {\mathrm{e}}^{ab} \f$
    #  @code
    #  from ostap.math.math_ve import exp 
    #  f1 = ...
    #  ff = exp ( f1 ) 
    #  @endcode 
    def  __exp__ ( self , other = 1 ) :
        """ Exponent f = exp(ab)
        >>> from ostap.math.math_ve import exp 
        >>> f1 = ...
        >>> ff = exp ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 )
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Exp    , pattern = "exp(%s*%s)" )  

    # ==============================================================================
    ## Logarithm \f$ f = \log ab \f$
    #  @code
    #  from ostap.math.math_ve import log
    #  f1 = ...
    #  ff = log ( f1 ) 
    #  @endcode 
    def  __log__ ( self , other = 1 ) :
        """ Exponent f = exp(ab)
        >>> from ostap.math.math_ve import log
        >>> f1 = ...
        >>> ff = log ( f1 ) 
        """
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Log    , pattern = "log(%s*%s)" )  

    # ==============================================================================
    ## Decimal logarithm \f$ f = \log_{10} ab \f$
    #  @code
    #  from ostap.math.math_ve import log10
    #  f1 = ...
    #  ff = log10 ( f1 ) 
    #  @endcode 
    def  __log10__ ( self , other = 1 ) :
        """ Decimal logarithm f = exp(ab)
        >>> from ostap.math.math_ve import log10
        >>> f1 = ...
        >>> ff = log10 ( f1 ) 
        """
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Log10 , pattern = "log10(%s*%s)" )  
        
    # ==============================================================================
    ## Error function \f$ f = erf ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import erf
    #  f1 = ...
    #  ff = erf( f1 ) 
    #  @endcode 
    def  __erf__ ( self , other = 1 ) :
        """ Error function f = erf(ab)
        >>> from ostap.math.math_ve import erf
        >>> f1 = ...
        >>> ff = erf ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 )
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Exp    , pattern = "erf(%s*%s)" )  

    # ==============================================================================
    ## Complementary error function \f$ f = erfc ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import erfc
    #  f1 = ...
    #  ff = erfc ( f1 ) 
    #  @endcode 
    def  __erfc__ ( self , other = 1 ) :
        """ Complementary error function f = erfc(ab)
        >>> from ostap.math.math_ve import erfc
        >>> f1 = ...
        >>> ff = erfc ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 )
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Erfc   , pattern = "erfc(%s*%s)" )  

    # ==============================================================================
    ## Sine function \f$ f =  \sin  ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import sin
    #  f1 = ...
    #  ff = sin ( f1 ) 
    #  @endcode 
    def  __sin__ ( self , other = 1 ) :
        """ Sine function f = sin(ab)
        >>> from ostap.math.math_ve import sin
        >>> f1 = ...
        >>> ff = sin ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 )
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Sin   , pattern = "sin(%s*%s)" )  

    # ==============================================================================
    ## Cosine function \f$ f =  \cos  ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import cos
    #  f1 = ...
    #  ff = cos ( f1 ) 
    #  @endcode 
    def  __cos__ ( self , other = 1 ) :
        """ Cosine function f = cos(ab)
        >>> from ostap.math.math_ve import cos 
        >>> f1 = ...
        >>> ff = cos ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 )
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Cos    , pattern = "cos(%s*%s)" )  

    # ==============================================================================
    ## Tangent function \f$ f =  \tan  ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import tan
    #  f1 = ...
    #  ff = tan ( f1 ) 
    #  @endcode 
    def  __tan__ ( self , other = 1 ) :
        """ Tangent function f = tan(ab)
        >>> from ostap.math.math_ve import tan
        >>> f1 = ...
        >>> ff = tan ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 )
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Tan    , pattern = "tan(%s*%s)" )  

    # ==============================================================================
    ## Hyperbolic sine function \f$ f =  \sinh ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import sinh
    #  f1 = ...
    #  ff = sinh ( f1 ) 
    #  @endcode 
    def  __sinh__ ( self , other = 1 ) :
        """ Hyperbolic sinh function f = sinh(ab)
        >>> from ostap.math.math_ve import sinh
        >>> f1 = ...
        >>> ff = sinh ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 )
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Sinh  , pattern = "sinh(%s*%s)" )  

    # ==============================================================================
    ## Hyperboilic sosine function \f$ f =  \cosh ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import cosh
    #  f1 = ...
    #  ff = cosh ( f1 ) 
    #  @endcode 
    def  __cosh__ ( self , other = 1 ) :
        """ Hyperbolic cosine function f = cosh(ab)
        >>> from ostap.math.math_ve import cosh 
        >>> f1 = ...
        >>> ff = cosh ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 )
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Cosh   , pattern = "cosh(%s*%s)" )  

    # ==============================================================================
    ## Hyperbolic tangent function \f$ f =  \tanh ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import tanh
    #  f1 = ...
    #  ff = tanh ( f1 ) 
    #  @endcode 
    def  __tanh__ ( self , other = 1 ) :
        """ Hyperbolic tangent function f = tanh(ab)
        >>> from ostap.math.math_ve import tanh
        >>> f1 = ...
        >>> ff = tanh ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 )
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Tanh   , pattern = "tanh(%s*%s)" )  

    # ==============================================================================
    ## Hyperbolic secant function \f$ f =  sech ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import sech
    #  f1 = ...
    #  ff = sech ( f1 ) 
    #  @endcode 
    def  __sech__ ( self , other = 1 ) :
        """ Hyperbolic secant function f = sech(ab)
        >>> from ostap.math.math_ve import sech
        >>> f1 = ...
        >>> ff = sech ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 )
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Sech   , pattern = "sech(%s*%s)" )  

    
    # ==============================================================================
    ## Inverse tangent function \f$ f =  \atan2 ( a , b )  \f$
    #  @code
    #  from ostap.math.math_ve import atan2
    #  f1 = ...
    #  ff = atan2 ( f1 ) 
    #  @endcode 
    def  __atan2__ ( self , other = 1 ) :
        """ Inverse tangent function f = atan2(a,b)
        >>> from ostap.math.math_ve import atan2
        >>> f1 = ...
        >>> ff = atan2 ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0.5 * math.pi )
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Atan2   , pattern = "atan2(%s,%s)" )  
    __atan__ = __atan2__


    # ==============================================================================
    ## maximal for two functions \f$ f =  \max ( a , b )  \f$
    #  @code
    #  from ostap.math.math_ve import maxv
    #  f1 = ...
    #  ff = maxv ( f1 ) 
    #  @endcode 
    def  __maxv__ ( self , other ) :
        """ Maximal of two functions f  = max(a,b)
        >>> from ostap.math.math_ve import maxv
        >>> f1 = ...
        >>> ff = maxv ( f1 ) 
        """
        if self is other : return self 
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.MaxV   , pattern = "maxv(%s,%s)" )  
    
    # ==============================================================================
    ## minimal for two functions \f$ f =  \min ( a , b )  \f$
    #  @code
    #  from ostap.math.math_ve import minv
    #  f1 = ...
    #  ff = minv ( f1 ) 
    #  @endcode 
    def  __minv__ ( self , other ) :
        """ Maximal of two functions f  = min(a,b)
        >>> from ostap.math.math_ve import minv
        >>> f1 = ...
        >>> ff = minv ( f1 ) 
        """
        if self is other : return self 
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.MinV   , pattern = "minv(%s,%s)" )  

    # ==============================================================================
    ## Gamma function \f$ f =  \Gamma( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import  gamma
    #  f1 = ...
    #  ff = gamma( f1 ) 
    #  @endcode 
    def  __tgamma__ ( self , other = 1 ) :
        """ Gamma function f = gamma(ab)
        >>> from ostap.math.math_ve import gamma
        >>> f1 = ...
        >>> ff = gamma( f1 ) 
        """
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Gamma   , pattern = "gamma(%s*%s)" )  
    __gamma__ = __tgamma__

    # ==============================================================================
    ## logarith of Gamma function \f$ f =  \ln \Gamma( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import  lgamma
    #  f1 = ...
    #  ff = lgamma( f1 ) 
    #  @endcode 
    def  __lgamma__ ( self , other = 1 ) :
        """ logarithm of Gamma function f = log(gamma(ab))
        >>> from ostap.math.math_ve import lgamma
        >>> f1 = ...
        >>> ff = lgamma( f1 ) 
        """
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.LGamma   , pattern = "lgamma(%s*%s)" )  
    
    # ==============================================================================
    ## 1/Gamma function \f$ f =  1/\Gamma( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import  igamma
    #  f1 = ...
    #  ff = igamma( f1 ) 
    #  @endcode 
    def  __igamma__ ( self , other = 1 ) :
        """ 1/Gamma function f = 1/gamma(ab))
        >>> from ostap.math.math_ve import igamma
        >>> f1 = ...
        >>> ff = igamma( f1 ) 
        """
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.IGamma   , pattern = "igamma(%s*%s)" )  

    # =========================================================================
    ## Fraction: \f$ f =  a / ( a + b ) \f$
    #  @code
    #  a =
    #  b = 
    #  ff = a.fraction ( b )
    #  @endcode 
    def fraction ( self , other ) :
        """ Fraction:  f =  a / ( a + b ) 
        >>> a =
        >>> b = 
        >>> ff = a.fraction ( b )
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 )
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Fraction , pattern = "fraction(%s,%s)" )
    
    # ==============================================================================
    ## Asymmetry  \f$ f =  ( a - b ) / ( a + b ) \f$
    #  @code
    #  a =
    #  b = 
    #  ff = a.asymmetry ( b )
    #  @endcode 
    def asymmetry ( self , b ) :
        """ Asymmetry  : f = ( a - b ) / ( a+ b )   )
        >>> f =
        >>> a = f.asymmetry ( b )
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 )
        return self.__f1_op__ ( other ,  Ostap.MoreRooFit.Asymmetry , pattern = "asymmetry(%s,%s)" )
    

    # =========================================================================
    ## helper function to create the histogram
    def make_histo ( self  ,
                     nbins    = 100   , xmin = None , xmax = None ,
                     hpars    = ()    , 
                     histo    = None  ) :
        """Create the histogram according to specifications
        """
        
        import ostap.histos.histos
        # histogram is provided 
        if histo :
            
            assert isinstance ( histo , ROOT.TH1 ) and not isinstance ( histo , ROOT.TH2 ) , \
                   "Illegal type of 'histo'-argument %s" % type( histo )
            
            histo = histo.clone()
            histo.Reset()

        # arguments for the histogram constructor 
        elif hpars :
            
            histo = ROOT.TH1F ( hID() , 'PDF%s' % self.name , *hpars  )
            if not histo.GetSumw2() : histo.Sumw2()

        # explicit construction from (#bins,min,max)-triplet  
        else :
            
            assert is_integer ( nbins ) and 0 < nbins, \
                   "Wrong 'nbins'-argument %s" % nbins 
            if xmin == None and self.xminmax() : xmin = self.xminmax()[0]
            if xmax == None and self.xminmax() : xmax = self.xminmax()[1]
            
            histo = ROOT.TH1F ( hID() , 'PDF%s' % self.name , nbins , xmin , xmax )
            if not histo.GetSumw2() : histo.Sumw2()

        return histo 

    # =============================================================================
    ## Get the derivative  dF/dx for the 1D-function
    #  \f[ f(x) = \frac{dF(x)}{dx}\f]
    #  @code
    #  F    = ...
    #  dFdx = F.dFdX ( ) 
    #  @endcode 
    #  @see RooAsbReal::derivative 
    def dFdX ( self , *args ) :
        """Get the derivative dF/dx for 1D-fuction
        >>> F    = ...
        >>> dFdx = F.dFdX ( ) 
        - see ROOT.RooAbsReal.derivative    
        """
        d = self.fun.derivative ( self.xvar , 1 , *args )
        return Fun1D ( d , self.xvar , name = self.new_name ( 'dFdX' ) )
    
    
# =============================================================================
## @class Fun1D
#  Simple wrapper for 1D-function
#  @code
#  func = ...
#  xvar = ...
#  f1d  = Fun1D ( func , xvar = xvar ) 
#  @endcode 
class Fun1D ( FUN1 ) :
    """Simple wrapper for 1D-function
    >>> func = ...
    >>> xvar = ...
    >>> f1d  = Fun1D ( func , xvar = xvar ) 
    """
    def __init__ ( self ,  fun , xvar , name = '' ) :

        self.__argfun = fun 
        if   isinstance ( fun , num_types ) :
            value = float ( fun ) 
            vname = 'C(%s,%s)' %  ( value , xvar.name ) 
            fun   = Ostap.MoreRooFit.Constant ( self.generate_name ( vname , name ) , 'constant %s' % value , value , xvar ) 
            
        elif isinstance ( fun , AFUN1     ) : fun = fun.fun
        
        assert xvar and isinstance ( xvar , ROOT.RooAbsReal ) , "'xvar' must be ROOT.RooAbsReal"
        assert fun  and isinstance ( fun  , ROOT.RooAbsReal ) , "'fun'  must be ROOT.RooAbsReal"

        if fun is xvar : fun = Ostap.MoreRooFit.Id ( "Id_%s" % xvar.name , "Id:%s" % xvar.title , xvar )            
            
        if not name : name = 'Fun1D_%s' % fun.GetName() 

        FUN1.__init__ ( self , name , xvar = xvar )

        if isinstance ( fun , ( ROOT.RooConstVar , Ostap.MoreRooFit.Id ) ) : pass 
        elif not self.xvar in fun.getParameters ( 0 ) and not self.xvar is fun : 
            self.warning ("Function does not depend on xvar=%s" % self.xvar.name )
            
        self.fun = fun
        
        self.config = {
            'fun'  : self.fun  ,
            'xvar' : self.xvar ,
            'name' : self.name ,            
            }
        
        self.checked_keys.add  ( 'fun'  ) 
        self.checked_keys.add  ( 'xvar' ) 


# =============================================================================
## @class Id
#  The most trivial function \f$ f(x) = x \f$
class Id ( FUN1 ) :
    """The most trivial function
    f(x) = x
    """
    def __init__ ( self , xvar ,  name = '' ) :
        
        FUN1.__init__ ( self , name = name if name else self.generate_name ( 'Id(%s)' % xvar.name ) , xvar = xvar )
        
        ## create the function itself
        self.fun = Ostap.MoreRooFit.Id ( self.roo_name ( 'Id' ) ,  '' , self.xvar )
        
        self.config = {
            'xvar' : self.xvar ,
            'name' : self.name ,            
            }
        
        self.checked_keys.add  ( 'xvar' ) 

# =============================================================================
## @class AFUN2
#  The base class for 2D-function
class AFUN2(AFUN1,YVar) :
    """Base class for 2D-function
    """
    def __init__ ( self , name , xvar , yvar , tricks = True , **kwargs ) :

        AFUN1 .__init__ ( self , name , xvar , tricks = tricks , **kwargs )
        YVar  .__init__ ( self , yvar )
        
        
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
        result = AFUN1.draw ( self            ,
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
## @class FUN2
#  The base class for 2D-function
class FUN2(AFUN2) :
    """Base class for 2D-function
    """
    def __init__ ( self , name , xvar , yvar , tricks = True , **kwargs ) :

        AFUN2 .__init__ ( self , name , xvar , yvar , tricks = tricks , **kwargs )
        ## save the configuration
        self.config = {
            'name'   : self.name   ,
            'xvar'   : self.xvar   ,
            'yvar'   : self.yvar   ,
            'tricks' : self.tricks , 
            }
        self.config.update ( kwargs )
        
        self.__call_OK = isinstance ( self.xvar , ROOT.RooAbsRealLValue ) and \
                         isinstance ( self.yvar , ROOT.RooAbsRealLValue ) 
        

    # =========================================================================
    ## simple 'function-like' interface 
    def __call__ ( self , x , y , normalized = False ) :
        """ Function as a 'function'
        >>> fun  = ...
        >>> x = 1
        >>> y = 2 
        >>> v = fun ( x , y ) 
        """
        assert self.__call_OK , "Invalid types for xvar/yvar!"


        xmnmx = self.xminmax()
        if xmnmx :
            xmn , xmx = xmnmx 
            if not xmn <= x <= xmx : return 0

        ymnmx = self.yminmax()
        if ymnmx :
            ymn , ymx = ymnmx 
            if not ymn <= y <= ymx : return 0
            
        with SETVAR( self.xvar ), SETVAR ( self.yvar )  :
            
            self.xvar.setVal ( x )
            self.yvar.setVal ( y )
            
            return self.fun.getVal ( self.vars ) if normalized else self.fun.getVal ()  

    # ========================================================================
    ## convert to float 
    def __float__ ( self ) :
        """Convert to float
        >>> pdf = ...
        >>> v = float ( pdf )
        """
        return self.fun.getVal () ## self.vars ) 


    # ==========================================================================
    ## convert FUN into TF3 object, e.g. to profit from TF2::Draw options
    #  @code
    #  fun = ...
    #  tf2 = fun.tf2 ()
    #  tf2.Draw('colz')
    #  @endcode
    def tf2 ( self , xmin = None , xmax = None , ymin = None , ymax = None  ) :
        """Convert PDF  to TF2 object, e.g. to profit from TF2::Draw options
        >>> fun = ...
        >>> tf2 = fun.tf2()
        >>> tf2.Draw('colz')
        """

        if xmin is None or xmax is None :
            xmnmn = self.xminmax ()
            if xmin is None and xmnmx : xmin = xmnmx[0]
            if xmax is None and xmnmx : xmax = xmnmx[1]

        if ymin is None or ymax is None :
            ymnmn = self.yminmax ()
            if ymin is None and ymnmx : ymin = ymnmx[0]
            if ymax is None and ymnmx : ymax = ymnmx[1]
        
        from ostap.core.core import fID
        def _aux_fun_ ( x , pars = [] ) : return self ( x[0] , x[1] )
        self.aux_keep.append ( _aux_fun_ )
        
        return ROOT.TF2 ( fID() , _aux_fun_ , xmin , xmax , ymin , ymax ) 


    # =========================================================================
    ## operations with FUN2 object: self (op) fun2 
    def __f2_op__ ( self , fun2 , ftype , pattern  ) :
        """Operations with FUN2 object: self (op) fun2
        """
        
        fun_name = lambda  name1 , name2 : self.generate_name ( pattern % ( name1 , name2 ) )
        
        if isinstance ( fun2 , AFUN3 ) and self.xvar in fun2.vars and self.yvar in fun2.vars :
            
            return Fun3op ( self , fun2 , operation = ftype , 
                            xvar    = fun2.xvar ,
                            yvar    = fun2.yvar ,
                            zvar    = fun2.zvar ,
                            name    = fun_name  ( self.name , fun2.name ) ,
                            strname = pattern % ( self.name , fun2.name ) ) 
        
        elif isinstance ( fun2 , AFUN3 ) :
            self.warning ( '1D (op) 3D : keep only 3D arguments %s ' % [ a.name for a in fun2.vars ] )                          
            return Fun3op ( self , fun2 , operation = ftype , 
                            xvar    = fun2.xvar ,
                            yvar    = fun2.yvar ,
                            zvar    = fun2.zvar ,
                            name    = fun_name  ( self.name , fun2.name ) ,
                            strname = pattern % ( self.name , fun2.name ) )
        
        elif isinstance ( fun2 , AFUN2 ) and self.xvar in fun2.vars and self.xvar in fun2.vars :        
            return Fun2op ( self , fun2 , operation = ftype , 
                            xvar    = self.xvar ,
                            yvar    = self.yvar ,
                            name    = fun_name   ( self.name , fun2.name ) ,
                            strname = pattern %  ( self.name , fun2.name ) )
        
        elif isinstance ( fun2 , AFUN2 ) and fun2.xvar in self.vars :            
            return Fun3op ( self , fun2 , operation = ftype , 
                            xvar    = self.xvar ,
                            yvar    = self.yvar ,
                            zvar    = self.yvar ,                            
                            name    = fun_name  ( self.name , fun2.name ) ,
                            strname = pattern % ( self.name , fun2.name ) )

        elif isinstance ( fun2 , AFUN2 ) and fun2.yvar in self.vars :
            return Fun3op ( self , fu2 , operation = ftype , 
                            xvar    = self.xvar ,
                            yvar    = self.yvar ,
                            zvar    = self.xvar ,
                            name    = fun_name  ( self.name , fun2.name ) ,
                            strname = pattern % ( self.name , fun2.name ) )

        elif isinstance ( fun2 , AFUN2 ) :
            self.warning ( '2D (op) 2D : keep only 2D arguments %s ' % [ a.name for a in self.vars ] )            
            return Fun2op ( self , fun2 , operation = ftype , 
                            xvar    = self.xvar ,
                            yvar    = self.yvar ,
                            name    = fun_name  ( self.name , fun2.name ) ,
                            strname = pattern % ( self.name , fun2.name ) )
        
        elif isinstance ( fun2 , AFUN1 ) and self.xvar in fun2.vars :        
            return Fun2op ( self , fun2 , operation = ftype , 
                            xvar    = self.xvar ,
                            yvar    = self.yvar ,
                            name    = fun_name  ( self.name , fun2.name ) ,
                            strname = pattern % ( self.name , fun2.name ) )
        
        elif isinstance ( fun2 , AFUN1 ) :             
            return Fun3op ( self , fun2 , operation = ftype , 
                            xvar    = self.xvar ,
                            yvar    = self.yvar ,
                            zvar    = fun2.xvar ,                            
                            name    = fun_name  ( self.name , fun2.name ) ,
                            strname = pattern % ( self.name , fun2.name ) )
        
        elif  isinstance ( fun2 , constant_types ) :
            fun2   = ROOT.RooFit.RooConst ( float ( fun2 )  )            
            return Fun2op ( self , fun2 , operaiton = ftype ,
                            xvar    = self.xvar ,
                            yvar    = self.yvar ,
                            name    = fun_name  ( self.name , fun2.name ) ,
                            strname = pattern % ( self.name , fun2.name ) )
        
        elif  isinstance ( fun2 , ROOT.RooAbsReal ) :            
            return Fun2op ( self . fun2 , operation = ftype ,
                            xvar    = self.xvar ,
                            yvar    = self.yvar ,
                            name    = fun_name  ( self.name , fun2.name ) ,
                            strname = pattern % ( self.name , fun2.name ) )
        return NotImplemented 

    # =========================================================================
    ## operations with FUN2 object: fun2 (op) self 
    def __f2_rop__ ( self , fun2 , ftype , pattern = '' ) :
        """operations with FUN2 object: fun2 (op) self 
        """
        
        fun_name = lambda  name1 , name2 : self.generate_name ( pattern % ( name1 , name2 ) )

        if  isinstance ( fun2 , constant_types ) :            
            fun2 = ROOT.RooFit.RooConst ( float ( fun2 )  )
            
        if  isinstance ( fun2 , ROOT.RooAbsReal ) :
            return Fun2op ( fun2 , self , operation = ftype ,
                            xvar     = self.xvar  ,
                            yvar     = self.yvar  ,
                            name     = fun_name ( fun2.name , self.name ) ,
                            str_name = pattern % ( fun2.name , self.name ) )
        
        return NotImplemented 
 
    # =========================================================================
    ## make a constant function 
    def __constant ( self , value ) :
        """ake a constant function"""
        return Fun2D ( value , xvar = self.xvar , yvar = self.yvar ,
                       name = self.new_name ( 'const2(%s)' % value ) ) 

    # =========================================================================
    ## Add two functions
    #  @code
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 + f2
    #  @endcode 
    def __add__ ( self , other ) :
        """Add two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 + f2 
        """
        if isinstance ( other , constant_types ) and iszero ( other ) : return self
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Addition    , pattern = "(%s+%s)" )
    
    # =========================================================================
    ## Subtract  two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 - f2 
    #  @endcode 
    def __sub__ ( self , other ) :
        """Subtract two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 + f2 
        """
        if   self is other                                              : return self.__constant ( 0 )
        elif isinstance ( other , constant_types ) and iszero ( other ) : return self
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Subtraction    , pattern = "(%s-%s)" )  
        
    # =========================================================================
    ## Multiply two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 * f2 
    #  @endcode 
    def __mul__ ( self , other ) :
        """Multiply two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 * f2 
        """
        if   isinstance ( other , constant_types ) and isequal ( other , 1 ) : return self
        elif isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 )
        elif self is other                                                   : return self ** 2 
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Product    , pattern = "(%s*%s)" )  
        
    # =========================================================================
    ## Divide two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 / f2 
    #  @endcode 
    def __div__ ( self , other ) :
        """Divide two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 / f2 
        """
        if self is other                                                     : return self.__constant ( 1 )
        elif isinstance ( other , constant_types ) and isequal ( other , 1 ) : return self
        elif isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 )
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Division    , pattern = "(%s/%s)" )  
    __truediv__ = __div__

    # =========================================================================
    ## Power function for two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 ** f2 
    #  @endcode 
    def __pow__ ( self , other ) :
        """Power function  for two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 ** f2 
        """
        if   isinstance ( other , constant_types ) and isequal ( other , 1 ) : return self
        elif isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 )
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Power    , pattern = 'pow(%s,%s)')    

    # =========================================================================
    ## Add two functions
    #  @code
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 + f2
    #  @endcode 
    def __radd__ ( self , other ) :
        """Add two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 + f2 
        """
        if isinstance ( other , constant_types ) and iszero ( other ) : return self
        return self.__f2_rop__ ( other ,  Ostap.MoreRooFit.Addition    , pattern = "(%s+%s)" )

    # =========================================================================
    ## Subtraction for two functions
    #  @code
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 - f2
    #  @endcode 
    def __rsub__ ( self , other ) :
        """Subtract two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 - f2 
        """
        if   self is other : return self.__constant ( 0 ) 
        return self.__f2_rop__ ( other ,  Ostap.MoreRooFit.Subtraction    , pattern = "(%s-%s)" )

    # =========================================================================
    ## Multiply two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 * f2 
    #  @endcode 
    def __rmul__ ( self , other ) :
        """Multiply two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 * f2 
        """
        if   isinstance ( other , constant_types ) and isequal ( other , 1 ) : return self
        elif isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 )
        return self.__f2_rop__ ( other ,  Ostap.MoreRooFit.Product    , pattern = "(%s*%s)" )  
        
    # =========================================================================
    ## Divide two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 / f2 
    #  @endcode 
    def __rdiv__ ( self , other ) :
        """Divide two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 / f2 
        """
        if   self is other                                                   : return self.__constant ( 1 )
        elif isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 )
        return self.__f2_rop__ ( other ,  Ostap.MoreRooFit.Division    , pattern = "(%s/%s)" )  
    __rtruediv__ = __rdiv__

    # =========================================================================
    ## Power function for two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 ** f2 
    #  @endcode 
    def __rpow__ ( self , other ) :
        """Power function  for two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 ** f2 
        """
        if   isinstance ( other , constant_types ) and isequal ( other , 1 ) : return self.__constant ( 1 )
        elif isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 )
        return self.__f2_rop__ ( other ,  Ostap.MoreRooFit.Power    , pattern = "pow(%s,%s)"  )  

    # =========================================================================
    ## absolute value : \f$ \left| bf(x) \right| \f$
    #  @code
    #  f1 = ...
    #  ff = abs ( f ) 
    #  @endcode
    def __abs__ ( self , other = 1 ) :
        """Absolute value :  |f(x)|
        >>>  f1 = ...
        >>> ff = abs ( f ) 
        """
        if   self is other                                                      : return self ** 2 
        elif isinstance ( other    , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 )
        elif isinstance ( self.fun , ROOT.RooAbsPdf ) and \
             isinstance ( other    , constant_types ) and isequal ( other , 1 ) : return self        
        ##
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Abs , pattern = "abs(%s*%s)" )  
        
    # ==============================================================================
    ## Exponent \f$ f = {\mathrm{e}}^{ab} \f$
    #  @code
    #  from ostap.math.math_ve import exp 
    #  f1 = ...
    #  ff = exp ( f1 ) 
    #  @endcode 
    def  __exp__ ( self , other = 1 ) :
        """ Exponent f = exp(ab)
        >>> from ostap.math.math_ve import exp 
        >>> f1 = ...
        >>> ff = exp ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 )
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Exp    , pattern = "exp(%s*%s)" )  

    # ==============================================================================
    ## Logarithm \f$ f = \log ab \f$
    #  @code
    #  from ostap.math.math_ve import log
    #  f1 = ...
    #  ff = log ( f1 ) 
    #  @endcode 
    def  __log__ ( self , other = 1 ) :
        """ Exponent f = exp(ab)
        >>> from ostap.math.math_ve import log
        >>> f1 = ...
        >>> ff = log ( f1 ) 
        """
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Log    , pattern = "log(%s*%s)" )  

    # ==============================================================================
    ## Decimal logarithm \f$ f = \log_{10} ab \f$
    #  @code
    #  from ostap.math.math_ve import log10
    #  f1 = ...
    #  ff = log10 ( f1 ) 
    #  @endcode 
    def  __log10__ ( self , other = 1 ) :
        """ Decimal logarithm f = exp(ab)
        >>> from ostap.math.math_ve import log10
        >>> f1 = ...
        >>> ff = log10 ( f1 ) 
        """
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Log10 , pattern = "log10(%s*%s)" )  
        
    # ==============================================================================
    ## Error function \f$ f = erf ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import erf
    #  f1 = ...
    #  ff = erf( f1 ) 
    #  @endcode 
    def  __erf__ ( self , other = 1 ) :
        """ Error function f = erf(ab)
        >>> from ostap.math.math_ve import erf
        >>> f1 = ...
        >>> ff = erf ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 )
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Exp    , pattern = "erf(%s*%s)" )  

    # ==============================================================================
    ## Complementary error function \f$ f = erfc ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import erfc
    #  f1 = ...
    #  ff = erfc ( f1 ) 
    #  @endcode 
    def  __erfc__ ( self , other = 1 ) :
        """ Complementary error function f = erfc(ab)
        >>> from ostap.math.math_ve import erfc
        >>> f1 = ...
        >>> ff = erfc ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 )
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Erfc   , pattern = "erfc(%s*%s)" )  

    # ==============================================================================
    ## Sine function \f$ f =  \sin  ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import sin
    #  f1 = ...
    #  ff = sin ( f1 ) 
    #  @endcode 
    def  __sin__ ( self , other = 1 ) :
        """ Sine function f = sin(ab)
        >>> from ostap.math.math_ve import sin
        >>> f1 = ...
        >>> ff = sin ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 )
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Sin   , pattern = "sin(%s*%s)" )  

    # ==============================================================================
    ## Cosine function \f$ f =  \cos  ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import cos
    #  f1 = ...
    #  ff = cos ( f1 ) 
    #  @endcode 
    def  __cos__ ( self , other = 1 ) :
        """ Cosine function f = cos(ab)
        >>> from ostap.math.math_ve import cos 
        >>> f1 = ...
        >>> ff = cos ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 )
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Cos    , pattern = "cos(%s*%s)" )  

    # ==============================================================================
    ## Tangent function \f$ f =  \tan  ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import tan
    #  f1 = ...
    #  ff = tan ( f1 ) 
    #  @endcode 
    def  __tan__ ( self , other = 1 ) :
        """ Tangent function f = tan(ab)
        >>> from ostap.math.math_ve import tan
        >>> f1 = ...
        >>> ff = tan ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 )
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Tan    , pattern = "tan(%s*%s)" )  

    # ==============================================================================
    ## Hyperbolic sine function \f$ f =  \sinh ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import sinh
    #  f1 = ...
    #  ff = sinh ( f1 ) 
    #  @endcode 
    def  __sinh__ ( self , other = 1 ) :
        """ Hyperbolic sinh function f = sinh(ab)
        >>> from ostap.math.math_ve import sinh
        >>> f1 = ...
        >>> ff = sinh ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 )
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Sinh  , pattern = "sinh(%s*%s)" )  

    # ==============================================================================
    ## Hyperboilic sosine function \f$ f =  \cosh ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import cosh
    #  f1 = ...
    #  ff = cosh ( f1 ) 
    #  @endcode 
    def  __cosh__ ( self , other = 1 ) :
        """ Hyperbolic cosine function f = cosh(ab)
        >>> from ostap.math.math_ve import cosh 
        >>> f1 = ...
        >>> ff = cosh ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 )
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Cosh   , pattern = "cosh(%s*%s)" )  

    # ==============================================================================
    ## Hyperbolic tangent function \f$ f =  \tanh ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import tanh
    #  f1 = ...
    #  ff = tanh ( f1 ) 
    #  @endcode 
    def  __tanh__ ( self , other = 1 ) :
        """ Hyperbolic tangent function f = tanh(ab)
        >>> from ostap.math.math_ve import tanh
        >>> f1 = ...
        >>> ff = tanh ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 )
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Tanh   , pattern = "tanh(%s*%s)" )  

    # ==============================================================================
    ## Hyperbolic secant function \f$ f =  sech ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import sech
    #  f1 = ...
    #  ff = sech ( f1 ) 
    #  @endcode 
    def  __sech__ ( self , other = 1 ) :
        """ Hyperbolic secant function f = sech(ab)
        >>> from ostap.math.math_ve import sech
        >>> f1 = ...
        >>> ff = sech ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 )
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Sech   , pattern = "sech(%s*%s)" )  

    # ==============================================================================
    ## Inverse tangent function \f$ f =  \atan2 ( a , b )  \f$
    #  @code
    #  from ostap.math.math_ve import atan2
    #  f1 = ...
    #  ff = atan2 ( f1 ) 
    #  @endcode 
    def  __atan2__ ( self , other = 1 ) :
        """ Inverse tangent function f = atan2(a,b)
        >>> from ostap.math.math_ve import atan2
        >>> f1 = ...
        >>> ff = atan2 ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0.5 * math.pi )
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Atan2   , pattern = "atan2(%s,%s)" )  
    __atan__ = __atan2__
    
    # ==============================================================================
    ## maximal for two functions \f$ f =  \max ( a , b )  \f$
    #  @code
    #  from ostap.math.math_ve import maxv
    #  f1 = ...
    #  ff = maxv ( f1 ) 
    #  @endcode 
    def  __maxv__ ( self , other ) :
        """ Maximal of two functions f  = max(a,b)
        >>> from ostap.math.math_ve import maxv
        >>> f1 = ...
        >>> ff = maxv ( f1 ) 
        """
        if self is other : return self 
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.MaxV   , pattern = "maxv(%s,%s)" )  
    
    # ==============================================================================
    ## minimal for two functions \f$ f =  \min ( a , b )  \f$
    #  @code
    #  from ostap.math.math_ve import minv
    #  f1 = ...
    #  ff = minv ( f1 ) 
    #  @endcode 
    def  __minv__ ( self , other ) :
        """ Maximal of two functions f  = min(a,b)
        >>> from ostap.math.math_ve import minv
        >>> f1 = ...
        >>> ff = minv ( f1 ) 
        """
        if self is other : return self 
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.MinV   , pattern = "minv(%s,%s)" )  

    # ==============================================================================
    ## Gamma function \f$ f =  \Gamma( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import  gamma
    #  f1 = ...
    #  ff = gamma( f1 ) 
    #  @endcode 
    def  __tgamma__ ( self , other = 1 ) :
        """ Gamma function f = gamma(ab)
        >>> from ostap.math.math_ve import gamma
        >>> f1 = ...
        >>> ff = gamma( f1 ) 
        """
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Gamma   , pattern = "gamma(%s*%s)" )  
    __gamma__ = __tgamma__

    # ==============================================================================
    ## logarith of Gamma function \f$ f =  \ln \Gamma( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import  lgamma
    #  f1 = ...
    #  ff = lgamma( f1 ) 
    #  @endcode 
    def  __lgamma__ ( self , other = 1 ) :
        """ logarithm of Gamma function f = log(gamma(ab))
        >>> from ostap.math.math_ve import lgamma
        >>> f1 = ...
        >>> ff = lgamma( f1 ) 
        """
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.LGamma   , pattern = "lgamma(%s*%s)" )  
    
    # ==============================================================================
    ## 1/Gamma function \f$ f =  1/\Gamma( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import  igamma
    #  f1 = ...
    #  ff = igamma( f1 ) 
    #  @endcode 
    def  __igamma__ ( self , other = 1 ) :
        """ 1/Gamma function f = 1/gamma(ab))
        >>> from ostap.math.math_ve import igamma
        >>> f1 = ...
        >>> ff = igamma( f1 ) 
        """
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.IGamma   , pattern = "igamma(%s*%s)" )  

    # =========================================================================
    ## Fraction: \f$ f =  a / ( a + b ) \f$
    #  @code
    #  a =
    #  b = 
    #  ff = a.fraction ( b )
    #  @endcode 
    def fraction ( self , other ) :
        """ Fraction:  f =  a / ( a + b ) 
        >>> a =
        >>> b = 
        >>> ff = a.fraction ( b )
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 )
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Fraction , pattern = "fraction(%s,%s)" )
    
    # ==============================================================================
    ## Asymmetry  \f$ f =  ( a - b ) / ( a + b ) \f$
    #  @code
    #  a =
    #  b = 
    #  ff = a.asymmetry ( b )
    #  @endcode 
    def asymmetry ( self , b ) :
        """ Asymmetry  : f = ( a - b ) / ( a+ b )   )
        >>> f =
        >>> a = f.asymmetry ( b )
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 )
        return self.__f2_op__ ( other ,  Ostap.MoreRooFit.Asymmetry , pattern = "asymmetry(%s,%s)" )
    
    # ==============================================================================
    ## Integrate 2D function over x
    #  \f[ f(y) = \int F(x,y) dx \f]
    #  @code
    #  f = ...
    #  g = f.integrate_x ( 'x-range' ) 
    #  @endcode 
    #  @see RooAbsReal::createIntegral
    def integrate_x ( self , *args ) :
        """ Integrate 2D-function over x  
        >>> f = ...
        >>> g = f.integrate_x ( 'x-range' ) 
        - see ROOT.RooAbsReal.createIntegral
        """
        ##
        vset   = ROOT.RooArgSet ( self.xvar ) 
        i      = self.fun.createIntegral ( vset , *args )
        return Fun1D ( i , self.yvar , name = self.new_name( 'integrateX' ) )
    
    # ==============================================================================
    ## Integrate 2D function over y
    #  \f[ f(x) = \int F(x,y) dy \f]
    #  @code
    #  f = ...
    #  g = f.integrate_y ( 'y-range' ) 
    #  @endcode 
    #  @see RooAbsReal::createIntegral
    def integrate_y ( self , *args ) :
        """ Integrate 2D-function over y  
        >>> f = ...
        >>> g = f.integrate_x ( 'x-range' ) 
        - see ROOT.RooAbsReal.createIntegral
        """
        vset = ROOT.RooArgSet ( self.yvar ) 
        i    = self.fun.createIntegral ( vset , *args )
        return Fun1D ( i , self.xvar , name = self.new_name( 'integrateY' ) )


    # =============================================================================
    ## Get the derivative  dF/dx for the 2D-function
    #  \f[ f(x,y) = \frac{dF(x,y)}{dx}\f]
    #  @code
    #  F    = ...
    #  dFdx = F.dFdX ( ) 
    #  @endcode 
    #  @see RooAsbReal::derivative 
    def dFdX ( self , *args ) :
        """Get the derivative dF/dx for 2D-fuction
        >>> F    = ...
        >>> dFdx = F.dFdX ( ) 
        - see ROOT.RooAbsReal.derivative    
        """
        d = self.fun.derivative ( self.xvar , 1 , *args )
        return Fun2D ( d , self.xvar , self.yvar , name = self.new_name ( 'dFdX' ) )

    # =============================================================================
    ## Get the derivative  dF/dy for the 2D-function
    #  \f[ f(x,y) = \frac{dF(x,y)}{dy}\f]
    #  @code
    #  F    = ...
    #  dFdx = F.dFdY ( ) 
    #  @endcode 
    #  @see RooAsbReal::derivative 
    def dFdY ( self , *args ) :
        """Get the derivative dF/dy for 2D-fuction
        >>> F    = ...
        >>> dFdy = F.dFdY ( ) 
        - see ROOT.RooAbsReal.derivative    
        """
        d = self.fun.derivative ( self.yvar , 1 , *args )
        return Fun2D ( d , self.xvar , self.yvar , name = self.new_name ( 'dFdY' ) )

# =============================================================================
## @class Fun2D
#  Simple wrapper for 2D-function
#  @code
#  func = ...
#  xvar = ...
#  yvar = ...
#  f2d  = Fun2D ( func , xvar = xvar , yvar = yvar ) 
#  @endcode 
class Fun2D ( FUN2 ) :
    """Simple wrapper for 2D-function
    >>> func = ...
    >>> xvar = ...
    >>> yvar = ...
    >>> f2d  = Fun2D ( func , xvar = xvar , yvar = yvar ) 
    """
    def __init__ ( self ,  fun , xvar , yvar , name = '' ) :
        
        self.__argfun = fun 
        if   isinstance ( fun , num_types ) :
            value = float ( fun ) 
            vname = 'C(%s,%s,%s)' %  ( value , xvar.name , yvar.name ) 
            fun   = Ostap.MoreRooFit.Constant ( self.generate_name ( vname , name ) , 'constant' , value , xvar , yvar ) 
            
        elif isinstance ( fun , AFUN1     ) : fun = fun.fun
            
        assert xvar and isinstance ( xvar , ROOT.RooAbsReal ) , "'xvar' must be ROOT.RooAbsReal"
        assert yvar and isinstance ( yvar , ROOT.RooAbsReal ) , "'yvar' must be ROOT.RooAbsReal"
        assert fun  and isinstance ( fun  , ROOT.RooAbsReal ) , "'fun'  must be ROOT.RooAbsReal"

        assert not xvar is yvar, "xvar and yvar must be different!"

        if fun is xvar : fun = Ostap.MoreRooFit.Id ( "Id_%s" % xvar.name , "Id:%s" % xvar.title , xvar )            
        if fun is yvar : fun = Ostap.MoreRooFit.Id ( "Id_%s" % yvar.name , "Id:%s" % yvar.title , yvar )            

        if not name : name = 'Fun2D_%s' % fun.GetName() 

        FUN2.__init__ ( self , name , xvar = xvar , yvar = yvar )
        
        if not isinstance ( fun , ( ROOT.RooConstVar , Ostap.MoreRooFit.Id ) ) : 
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
        self.checked_keys.add  ( 'fun'  ) 
        self.checked_keys.add  ( 'xvar' ) 
        self.checked_keys.add  ( 'yvar' ) 

# =============================================================================
## @class AFUN3
#  The base class for 3D-function
class AFUN3(AFUN2,ZVar) :
    """Base class for 3D-function
    """

    def __init__ ( self , name , xvar , yvar , zvar , tricks = True , **kwargs ) :

        AFUN2.__init__ ( self , name , xvar , yvar , tricks = tricks , **kwargs )
        ZVar .__init__ ( self , zvar )

                
        self.vars     .add    ( self.zvar )
        self.variables.append ( self.zvar )
        
        ## save the configuration
        self.config = {
            'name' : self.name ,
            'xvar' : self.xvar ,
            'yvar' : self.yvar ,            
            'zvar' : self.zvar ,            
            }
        self.config.update ( kwargs )
        
    ## conversion to string 
    def __str__ (  self ) :
        return '%s(%s,xvar=%s,yvar=%s,zvar=%s)' % ( self.__class__.__name__ ,
                                                    self.name      ,
                                                    self.xvar.name ,
                                                    self.yvar.name ,
                                                    self.zvar.name )
    __repr__ = __str__ 

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
            
        return FUN2.draw  ( self ,
                            drawvar  = drawvar  ,
                            silent   =  silent  ,
                            in_range = in_range ,
                            args     = args     , **kwargs )




# =============================================================================
## @class FUN3
#  The base class for 3D-function
class FUN3(AFUN3) :
    """Base class for 3D-function
    """
    def __init__ ( self , name , xvar , yvar , zvar , tricks = True , **kwargs ) :

        AFUN3.__init__ ( self , name , xvar , yvar , zvar , tricks = tricks , **kwargs )
        
        self.vars     .add    ( self.zvar )
        self.variables.append ( self.zvar )
        
        ## save the configuration
        self.config = {
            'name'   : self.name   ,
            'xvar'   : self.xvar   ,
            'yvar'   : self.yvar   ,            
            'zvar'   : self.zvar   ,            
            'tricks' : self.tricks ,            
            }
        
        self.__call_OK = isinstance ( self.xvar , ROOT.RooAbsRealLValue ) and \
                         isinstance ( self.yvar , ROOT.RooAbsRealLValue ) and \
                         isinstance ( self.zvar , ROOT.RooAbsRealLValue )

        
    # =========================================================================
    ## simple 'function-like' interface 
    def __call__ ( self , x , y , z , normalized = False ) :
        """ Function as a 'function'
        >>> fun  = ...
        >>> x = 1
        >>> y = 2 
        >>> z = 0 
        >>> v = fun ( x , y , z ) 
        """
        assert self.__call_OK , "Invalid types for xvar/yvar/zvar!"
        
        xmnmx = self.xminmax()
        if xmnmx :
            xmn , xmx = xmnmx 
            if not xmn <= x <= xmx : return 0

        ymnmx = self.yminmax()
        if ymnmx :
            ymn , yxmx = ymnmx 
            if not ymn <= y <= ymx : return 0

        zmnmx = self.zminmax()
        if zmnmx :
            zmn , zxmx = zmnmx 
            if not zmn <= z <= zmx : return 0

        with SETVAR( self.xvar ), SETVAR ( self.yvar ) , SETVAR ( self.zvar ) :
            
            self.xvar.setVal ( x )
            self.yvar.setVal ( y )
            self.zvar.setVal ( z )
            
            return self.fun.getVal ( self.vars ) if normalized else self.fun.getVal ()  

    # ==========================================================================
    ## convert FUN into TF3 object, e.g. to profit from TF3::Draw options
    #  @code
    #  fun = ...
    #  tf3 = fun.tf3 ()
    #  tf3.Draw('colz')
    #  @endcode
    def tf3 ( self ,
              xmin = None , xmax = None ,
              ymin = None , ymax = None ,
              zmin = None , zmax = None ) :
        """Convert PDF  to TF3 object, e.g. to profit from TF3::Draw options
        >>> fun = ...
        >>> tf3 = fun.tf3()
        >>> tf3.Draw('colz')
        """

        if xmin is None or xmax is None :
            xmnmn = self.xminmax ()
            if xmin is None and xmnmx : xmin = xmnmx[0]
            if xmax is None and xmnmx : xmax = xmnmx[1]

        if ymin is None or ymax is None :
            ymnmn = self.yminmax ()
            if ymin is None and ymnmx : ymin = ymnmx[0]
            if ymax is None and ymnmx : ymax = ymnmx[1]

        if zmin is None or zmax is None :
            zmnmn = self.zminmax ()
            if zmin is None and zmnmx : zmin = zmnmx[0]
            if zmax is None and zmnmx : zmax = zmnmx[1]
        
        from ostap.core.core import fID
        def _aux_fun_ ( x , pars = [] ) : return self ( x[0] , x[1] , x[2]  )
        self.aux_keep.append ( _aux_fun_ )
        
        return ROOT.TF3 ( fID() , _aux_fun_ ,
                          xmin , xmax ,
                          ymin , ymax ,
                          zmin , zmax ) 
    
    # =========================================================================
    ## operations with FUN3 object: self (op) fun2 
    def __f3_op__ ( self , fun2 , ftype ,pattern ) :
        """Operations with FUN3 object: s elf (op) fun2
        """
        
        fun_name = lambda  name1 , name2 : self.generate_name ( pattern % ( name1 , name2 ) )
        
        if isinstance ( fun2 ,  ( AFUN3 , AFUN2 , AFUN1 ) )  : 
            return Fun3op ( self , fun2 , operation = ftype , 
                            xvar    = self.xvar ,
                            yvar    = self.yvar ,
                            zvar    = self.zvar ,
                            name    = fun_name  ( self.name , fun2.name ) ,
                            strname = pattern % ( self.name , fun2.name ) )
       
        elif  isinstance ( fun2 , constant_types ) :            
            fun2   = ROOT.RooFit.RooConst ( float ( fun2 )  )
            return Fun3op ( self , fun2 , operation = ftype , 
                            xvar    = self.xvar ,
                            yvar    = self.yvar ,
                            zvar    = self.zvar ,
                            name    = fun_name  ( self.name , fun2.name ) ,
                            strname = pattern % ( self.name , fun2.name ) )
        
        elif  isinstance ( fun2 , ROOT.RooAbsReal ) :
            return Fun3op ( self , fun2 , operation = ftype , 
                            xvar    = self.xvar ,
                            yvar    = self.yvar ,
                            zvar    = self.zvar ,
                            name    = fun_name  ( self.name , fun2.name ) ,
                            strname = pattern % ( self.name , fun2.name ) )
        
        return NotImplemented 

    # =========================================================================
    ## operations with FUN3 object: fun2 (op) self 
    def __f3_rop__ ( self , fun2 , ftype , opname = '' , pattern = '' ) :
        """operations with FUN3 object: fun2 (op) self 
        """
        
        fun_name = lambda  name1 , name2 : self.generate_name ( pattern % ( name1 , name2 ) )
        
        if  isinstance ( fun2 , constant_types ) :            
            fun2 = ROOT.RooFit.RooConst ( float ( fun2 )  )
            
        if  isinstance ( fun2 , ROOT.RooAbsReal ) :
            return Fun3op ( fun2 , self , operation = ftype ,
                            xvar    = self.xvar  ,
                            yvar    = self.yvar  ,
                            zvar    = self.zvar  ,
                            name    = fun_name   ( fun2.name , self.name ) ,
                            strname = pattern  % ( fun2.name , self.name ) )
        
        return NotImplemented 
    
    # =========================================================================
    ## make a constant function 
    def __constant ( self , value ) :
        """ake a constant function"""
        return Fun3D ( value , xvar = self.xvar , yvar = self.yvar , zvar = self.zvar , 
                       name = self.new_name ( 'const3(%s)' % value ) ) 

    # =========================================================================
    ## Add two functions
    #  @code
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 + f2
    #  @endcode 
    def __add__ ( self , other ) :
        """Add two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 + f2 
        """
        if isinstance ( other , constant_types ) and iszero ( other ) : return self
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Addition    , pattern = "(%s+%s)" )
    
    # =========================================================================
    ## Subtract  two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 - f2 
    #  @endcode 
    def __sub__ ( self , other ) :
        """Subtract two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 + f2 
        """
        if   self is other                                              : return self.__constant ( 0 ) 
        elif isinstance ( other , constant_types ) and iszero ( other ) : return self
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Subtraction    , pattern = "(%s-%s)" )  
        
    # =========================================================================
    ## Multiply two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 * f2 
    #  @endcode 
    def __mul__ ( self , other ) :
        """Multiply two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 * f2 
        """
        if   isinstance ( other , constant_types ) and isequal ( other , 1 ) : return self
        elif isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 ) 
        elif self is other                                                   : return self ** 2 
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Product    , pattern = "(%s*%s)" )  
        
    # =========================================================================
    ## Divide two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 / f2 
    #  @endcode 
    def __div__ ( self , other ) :
        """Divide two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 / f2 
        """
        if self is other                                                     : return self.__constant ( 1 ) 
        elif isinstance ( other , constant_types ) and isequal ( other , 1 ) : return self
        elif isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 ) 
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Division    , pattern = "(%s/%s)" )  
    __truediv__ = __div__

    # =========================================================================
    ## Power function for two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 ** f2 
    #  @endcode 
    def __pow__ ( self , other ) :
        """Power function  for two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 ** f2 
        """
        if   isinstance ( other , constant_types ) and isequal ( other , 1 ) : return self
        elif isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 ) 
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Power    , pattern = 'pow(%s,%s)' )  

    # =========================================================================
    ## Add two functions
    #  @code
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 + f2
    #  @endcode 
    def __radd__ ( self , other ) :
        """Add two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 + f2 
        """
        if isinstance ( other , constant_types ) and iszero ( other ) : return self
        return self.__f3_rop__ ( other ,  Ostap.MoreRooFit.Addition    , pattern = "(%s+%s)" )

    # =========================================================================
    ## Subtraction for two functions
    #  @code
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 - f2
    #  @endcode 
    def __rsub__ ( self , other ) :
        """Subtract two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 - f2 
        """
        if   self is other : return self.__constant ( 0 ) 
        return self.__f3_rop__ ( other ,  Ostap.MoreRooFit.Subtraction    , pattern = "(%s-%s)" )

    # =========================================================================
    ## Multiply two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 * f2 
    #  @endcode 
    def __rmul__ ( self , other ) :
        """Multiply two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 * f2 
        """
        if   isinstance ( other , constant_types ) and isequal ( other , 1 ) : return self
        elif isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 ) 
        return self.__f3_rop__ ( other ,  Ostap.MoreRooFit.Product    , pattern = "(%s*%s)" )  
        
    # =========================================================================
    ## Divide two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 / f2 
    #  @endcode 
    def __rdiv__ ( self , other ) :
        """Divide two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 / f2 
        """
        if   self is other                                                   : return self.__constant ( 1 ) 
        elif isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 ) 
        return self.__f3_rop__ ( other ,  Ostap.MoreRooFit.Division    , pattern = "(%s/%s)" )  
    __rtruediv__ = __rdiv__

    # =========================================================================
    ## Power function for two functions
    #  @code 
    #  f1 = ... 
    #  f2 = ...
    #  ff = f1 ** f2 
    #  @endcode 
    def __rpow__ ( self , other ) :
        """Power function  for two functions
        >>> f1 = ... 
        >>> f2 = ...
        >>> ff = f1 ** f2 
        """
        if   isinstance ( other , constant_types ) and isequal ( other , 1 ) : return self.__constant ( 1 ) 
        return self.__f3_rop__ ( other ,  Ostap.MoreRooFit.Power    , pattern = 'pow(%s,%s)' )  


    # =========================================================================
    ## absolute value : \f$ \left| bf(x) \right| \f$
    #  @code
    #  f1 = ...
    #  ff = abs ( f ) 
    #  @endcode
    def __abs__ ( self , other = 1 ) :
        """Absolute value : |f(x)|
        >>>  f1 = ...
        >>> ff = abs ( f ) 
        """
        if   isinstance ( other    , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 ) 
        elif isinstance ( self.fun , ROOT.RooAbsPdf ) and \
             isinstance ( other    , constant_types ) and isequal ( other , 1 ) : return self        
        ##
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Abs , pattern = "abs(%s*%s)" )  
        
    # ==============================================================================
    ## Exponent \f$ f = {\mathrm{e}}^{ab} \f$
    #  @code
    #  from ostap.math.math_ve import exp 
    #  f1 = ...
    #  ff = exp ( f1 ) 
    #  @endcode 
    def  __exp__ ( self , other = 1 ) :
        """ Exponent f = exp(ab)
        >>> from ostap.math.math_ve import exp 
        >>> f1 = ...
        >>> ff = exp ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 )  
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Exp    , pattern = "exp(%s*%s)" )  

    # ==============================================================================
    ## Logarithm \f$ f = \log ab \f$
    #  @code
    #  from ostap.math.math_ve import log
    #  f1 = ...
    #  ff = log ( f1 ) 
    #  @endcode 
    def  __log__ ( self , other = 1 ) :
        """ Exponent f = exp(ab)
        >>> from ostap.math.math_ve import log
        >>> f1 = ...
        >>> ff = log ( f1 ) 
        """
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Log    , pattern = "log(%s*%s)" )  

    # ==============================================================================
    ## Decimal logarithm \f$ f = \log_{10} ab \f$
    #  @code
    #  from ostap.math.math_ve import log10
    #  f1 = ...
    #  ff = log10 ( f1 ) 
    #  @endcode 
    def  __log10__ ( self , other = 1 ) :
        """ Decimal logarithm f = exp(ab)
        >>> from ostap.math.math_ve import log10
        >>> f1 = ...
        >>> ff = log10 ( f1 ) 
        """
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Log10 , pattern = "log10(%s*%s)" )  
        
    # ==============================================================================
    ## Error function \f$ f = erf ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import erf
    #  f1 = ...
    #  ff = erf( f1 ) 
    #  @endcode 
    def  __erf__ ( self , other = 1 ) :
        """ Error function f = erf(ab)
        >>> from ostap.math.math_ve import erf
        >>> f1 = ...
        >>> ff = erf ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 ) 
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Exp    , pattern = "erf(%s*%s)" )  

    # ==============================================================================
    ## Complementary error function \f$ f = erfc ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import erfc
    #  f1 = ...
    #  ff = erfc ( f1 ) 
    #  @endcode 
    def  __erfc__ ( self , other = 1 ) :
        """ Complementary error function f = erfc(ab)
        >>> from ostap.math.math_ve import erfc
        >>> f1 = ...
        >>> ff = erfc ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 ) 
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Erfc   , pattern = "erfc(%s*%s)" )  

    # ==============================================================================
    ## Sine function \f$ f =  \sin  ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import sin
    #  f1 = ...
    #  ff = sin ( f1 ) 
    #  @endcode 
    def  __sin__ ( self , other = 1 ) :
        """ Sine function f = sin(ab)
        >>> from ostap.math.math_ve import sin
        >>> f1 = ...
        >>> ff = sin ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 ) 
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Sin   , pattern = "sin(%s*%s)" )  

    # ==============================================================================
    ## Cosine function \f$ f =  \cos  ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import cos
    #  f1 = ...
    #  ff = cos ( f1 ) 
    #  @endcode 
    def  __cos__ ( self , other = 1 ) :
        """ Cosine function f = cos(ab)
        >>> from ostap.math.math_ve import cos 
        >>> f1 = ...
        >>> ff = cos ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 ) 
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Cos    , pattern = "cos(%s*%s)" )  

    # ==============================================================================
    ## Tangent function \f$ f =  \tan  ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import tan
    #  f1 = ...
    #  ff = tan ( f1 ) 
    #  @endcode 
    def  __tan__ ( self , other = 1 ) :
        """ Tangent function f = tan(ab)
        >>> from ostap.math.math_ve import tan
        >>> f1 = ...
        >>> ff = tan ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 ) 
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Tan    , pattern = "tan(%s*%s)" )  

    # ==============================================================================
    ## Hyperbolic sine function \f$ f =  \sinh ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import sinh
    #  f1 = ...
    #  ff = sinh ( f1 ) 
    #  @endcode 
    def  __sinh__ ( self , other = 1 ) :
        """ Hyperbolic sinh function f = sinh(ab)
        >>> from ostap.math.math_ve import sinh
        >>> f1 = ...
        >>> ff = sinh ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 ) 
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Sinh  , pattern = "sinh(%s*%s)" )  

    # ==============================================================================
    ## Hyperboilic sosine function \f$ f =  \cosh ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import cosh
    #  f1 = ...
    #  ff = cosh ( f1 ) 
    #  @endcode 
    def  __cosh__ ( self , other = 1 ) :
        """ Hyperbolic cosine function f = cosh(ab)
        >>> from ostap.math.math_ve import cosh 
        >>> f1 = ...
        >>> ff = cosh ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 ) 
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Cosh   , pattern = "cosh(%s*%s)" )  

    # ==============================================================================
    ## Hyperbolic tangent function \f$ f =  \tanh ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import tanh
    #  f1 = ...
    #  ff = tanh ( f1 ) 
    #  @endcode 
    def  __tanh__ ( self , other = 1 ) :
        """ Hyperbolic tangent function f = tanh(ab)
        >>> from ostap.math.math_ve import tanh
        >>> f1 = ...
        >>> ff = tanh ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 0 ) 
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Tanh   , pattern = "tanh(%s*%s)" )  

    # ==============================================================================
    ## Hyperbolic secant function \f$ f =  sech ( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import sech
    #  f1 = ...
    #  ff = sech ( f1 ) 
    #  @endcode 
    def  __sech__ ( self , other = 1 ) :
        """ Hyperbolic secant function f = sech(ab)
        >>> from ostap.math.math_ve import sech
        >>> f1 = ...
        >>> ff = sech ( f1 ) 
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 )
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Sech   , pattern = "sech(%s*%s)" )  

    # ==============================================================================
    ## Inverse tangent function \f$ f =  \atan2 ( a , b )  \f$
    #  @code
    #  from ostap.math.math_ve import atan2
    #  f1 = ...
    #  ff = atan2 ( f1 ) 
    #  @endcode 
    def  __atan2__ ( self , other = 1 ) :
        """ Inverse tangent function f = atan2(a,b)
        >>> from ostap.math.math_ve import atan2
        >>> f1 = ...
        >>> ff = atan2 ( f1 ) 
        """
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Atan2   , pattern = "atan2(%s,%s)" )  
    __atan__ = __atan2__

    # ==============================================================================
    ## maximal for two functions \f$ f =  \max ( a , b )  \f$
    #  @code
    #  from ostap.math.math_ve import maxv
    #  f1 = ...
    #  ff = maxv ( f1 ) 
    #  @endcode 
    def  __maxv__ ( self , other ) :
        """ Maximal of two functions f  = max(a,b)
        >>> from ostap.math.math_ve import maxv
        >>> f1 = ...
        >>> ff = maxv ( f1 ) 
        """
        if self is other : return self 
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.MaxV   , pattern = "maxv(%s,%s)" )  
    
    # ==============================================================================
    ## minimal for two functions \f$ f =  \min ( a , b )  \f$
    #  @code
    #  from ostap.math.math_ve import minv
    #  f1 = ...
    #  ff = minv ( f1 ) 
    #  @endcode 
    def  __minv__ ( self , other ) :
        """ Maximal of two functions f  = min(a,b)
        >>> from ostap.math.math_ve import minv
        >>> f1 = ...
        >>> ff = minv ( f1 ) 
        """
        if self is other : return self 
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.MinV   , pattern = "minv(%s,%s)" )  
    
    # ==============================================================================
    ## Gamma function \f$ f =  \Gamma( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import  gamma
    #  f1 = ...
    #  ff = gamma( f1 ) 
    #  @endcode 
    def  __tgamma__ ( self , other = 1 ) :
        """ Gamma function f = gamma(ab)
        >>> from ostap.math.math_ve import gamma
        >>> f1 = ...
        >>> ff = gamma( f1 ) 
        """
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Gamma   , pattern = "gamma(%s*%s)" )  
    __gamma__ = __tgamma__

    # ==============================================================================
    ## logarith of Gamma function \f$ f =  \ln \Gamma( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import  lgamma
    #  f1 = ...
    #  ff = lgamma( f1 ) 
    #  @endcode 
    def  __lgamma__ ( self , other = 1 ) :
        """ logarithm of Gamma function f = log(gamma(ab))
        >>> from ostap.math.math_ve import lgamma
        >>> f1 = ...
        >>> ff = lgamma( f1 ) 
        """
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.LGamma   , pattern = "lgamma(%s*%s)" )  
    
    # ==============================================================================
    ## 1/Gamma function \f$ f =  1/\Gamma( ab )  \f$
    #  @code
    #  from ostap.math.math_ve import  igamma
    #  f1 = ...
    #  ff = igamma( f1 ) 
    #  @endcode 
    def  __igamma__ ( self , other = 1 ) :
        """ 1/Gamma function f = 1/gamma(ab))
        >>> from ostap.math.math_ve import igamma
        >>> f1 = ...
        >>> ff = igamma( f1 ) 
        """
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.IGamma   , pattern = "igamma(%s*%s)" )  

    # =========================================================================
    ## Fraction: \f$ f =  a / ( a + b ) \f$
    #  @code
    #  a =
    #  b = 
    #  ff = a.fraction ( b )
    #  @endcode 
    def fraction ( self , other ) :
        """ Fraction:  f =  a / ( a + b ) 
        >>> a =
        >>> b = 
        >>> ff = a.fraction ( b )
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 ) 
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Fraction , pattern = "fraction(%s,%s)" )
    
    # ==============================================================================
    ## Asymmetry  \f$ f =  ( a - b ) / ( a + b ) \f$
    #  @code
    #  a =
    #  b = 
    #  ff = a.asymmetry ( b )
    #  @endcode 
    def asymmetry ( self , b ) :
        """ Asymmetry  : f = ( a - b ) / ( a+ b )   )
        >>> f =
        >>> a = f.asymmetry ( b )
        """
        if isinstance ( other , constant_types ) and iszero  ( other     ) : return self.__constant ( 1 ) 
        return self.__f3_op__ ( other ,  Ostap.MoreRooFit.Asymmetry , pattern = "asymmetry(%s,%s)" )

    # ==============================================================================
    ## Integrate 3D function over x
    #  \f[ f(y,z) = \int F(x,y,z) dx \f]
    #  @code
    #  f = ...
    #  g = f.integrate_x ( 'x-range' ) 
    #  @endcode 
    #  @see RooAbsReal::createIntegral
    def integrate_x ( self , *args ) :
        """ Integrate 3D-function over x  
        >>> f = ...
        >>> g = f.integrate_x ( 'x-range' ) 
        - see ROOT.RooAbsReal.createIntegral
        """
        ##
        vset   = ROOT.RooArgSet ( self.xvar ) 
        i      = self.fun.createIntegral ( vset , *args )
        return Fun2D ( i , self.yvar , self.zvar , name = self.new_name( 'integrateX' ) )

    # ==============================================================================
    ## Integrate 3D function over y
    #  \f[ f(x,z) = \int F(x,y,z) dy \f]
    #  @code
    #  f = ...
    #  g = f.integrate_y ( 'y-range' ) 
    #  @endcode 
    #  @see RooAbsReal::createIntegral
    def integrate_y ( self , *args ) :
        """ Integrate 3D-function over y  
        >>> f = ...
        >>> g = f.integrate_y ( 'y-range' ) 
        - see ROOT.RooAbsReal.createIntegral
        """
        ##
        vset   = ROOT.RooArgSet ( self.yvar ) 
        i      = self.fun.createIntegral ( vset , *args )
        return Fun2D ( i , self.xvar , self.zvar , name = self.new_name( 'integrateY' ) )

    # ==============================================================================
    ## Integrate 3D function over z
    #  \f[ f(x,y) = \int F(x,y,z) dz \f]
    #  @code
    #  f = ...
    #  g = f.integrate_z ( 'z-range' ) 
    #  @endcode 
    #  @see RooAbsReal::createIntegral
    def integrate_z ( self , *args ) :
        """ Integrate 3D-function over z  
        >>> f = ...
        >>> g = f.integrate_z ( 'z-range' ) 
        - see ROOT.RooAbsReal.createIntegral
        """
        ##
        vset   = ROOT.RooArgSet ( self.zvar ) 
        i      = self.fun.createIntegral ( vset , *args )
        return Fun2D ( i , self.xvar , self.yvar , name = self.new_name( 'integrateXY' ) )
    
    # ==============================================================================
    ## Integrate 3D function over x,y
    #  \f[ f(z) = \int\int F(x,y,z) dx dy\f]
    #  @code
    #  f = ...
    #  g = f.integrate_xy ( 'xy-range' ) 
    #  @endcode 
    #  @see RooAbsReal::createIntegral
    def integrate_xy ( self , *args ) :
        """ Integrate 3D-function over x,y  
        >>> f = ...
        >>> g = f.integrate_xy ( 'xy-range' ) 
        - see ROOT.RooAbsReal.createIntegral
        """
        ##
        vset   = ROOT.RooArgSet ( self.xvar , self.yvar ) 
        i      = self.fun.createIntegral ( vset , *args )
        return Fun1D ( i , self.zvar , name = self.new_name( 'integrateXY' ) )

    # ==============================================================================
    ## Integrate 3D function over x,z
    #  \f[ f(y) = \int\int F(x,y,z) dx dz\f]
    #  @code
    #  f = ...
    #  g = f.integrate_xz ( 'xz-range' ) 
    #  @endcode 
    #  @see RooAbsReal::createIntegral
    def integrate_xz ( self , *args ) :
        """ Integrate 3D-function over x,z  
        >>> f = ...
        >>> g = f.integrate_xz ( 'xz-range' ) 
        - see ROOT.RooAbsReal.createIntegral
        """
        ##
        vset   = ROOT.RooArgSet ( self.xvar , self.zvar ) 
        i      = self.fun.createIntegral ( vset , *args )
        return Fun1D ( i , self.yvar , name = self.new_name( 'integrateXZ' ) )

    # ==============================================================================
    ## Integrate 3D function over y,z
    #  \f[ f(x) = \int\int F(x,y,z) dy dz\f]
    #  @code
    #  f = ...
    #  g = f.integrate_yz ( 'yz-range' ) 
    #  @endcode 
    #  @see RooAbsReal::createIntegral
    def integrate_yz ( self , *args ) :
        """ Integrate 3D-function over x,z  
        >>> f = ...
        >>> g = f.integrate_yz ( 'yz-range' ) 
        - see ROOT.RooAbsReal.createIntegral
        """
        ##
        vset   = ROOT.RooArgSet ( self.yvar , self.zvar ) 
        i      = self.fun.createIntegral ( vset , *args )
        return Fun1D ( i , self.zvar , name = self.new_name( 'integrateYZ' ) )


    # =============================================================================
    ## Get the derivative  dF/dx for the 3D-function
    #  \f[ f(x,y,z) = \frac{dF(x,y,z)}{dx}\f]
    #  @code
    #  F    = ...
    #  dFdx = F.dFdX ( ) 
    #  @endcode 
    #  @see RooAsbReal::derivative 
    def dFdX ( self , *args ) :
        """Get the derivative dF/dx for 3D-fuction
        >>> F    = ...
        >>> dFdx = F.dFdX ( ) 
        - see ROOT.RooAbsReal.derivative    
        """
        d = self.fun.derivative ( self.xvar , 1 , *args )
        return Fun3D ( d , self.xvar , self.yvar , self.zvar , name = self.new_name ( 'dFdX' ) )

    # =============================================================================
    ## Get the derivative  dF/dy for the 3D-function
    #  \f[ f(x,y,z) = \frac{dF(x,y,z)}{dy}\f]
    #  @code
    #  F    = ...
    #  dFdx = F.dFdY ( ) 
    #  @endcode 
    #  @see RooAsbReal::derivative 
    def dFdY ( self , *args ) :
        """Get the derivative dF/dy for 3D-fuction
        >>> F    = ...
        >>> dFdy = F.dFdY ( ) 
        - see ROOT.RooAbsReal.derivative    
        """
        d = self.fun.derivative ( self.yvar , 1 , *args )
        return Fun3D ( d , self.xvar , self.yvar , self.zvar , name = self.new_name ( 'dFdY' ) )

    # =============================================================================
    ## Get the derivative  dF/dz for the 3D-function
    #  \f[ f(x,y,z) = \frac{dF(x,y,z)}{dz}\f]
    #  @code
    #  F    = ...
    #  dFdx = F.dFdZ ( ) 
    #  @endcode 
    #  @see RooAsbReal::derivative 
    def dFdZ ( self , *args ) :
        """Get the derivative dF/dz for 3D-fuction
        >>> F    = ...
        >>> dFdy = F.dFdD ( ) 
        - see ROOT.RooAbsReal.derivative    
        """
        d = self.fun.derivative ( self.zvar , 1 , *args )
        return Fun3D ( d , self.xvar , self.yvar , self.zvar , name = self.new_name ( 'dFdZ' ) )

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
class Fun3D ( FUN3 ) :
    """Simple wrapper for 3D-function
    >>> func = ...
    >>> xvar = ...
    >>> yvar = ...
    >>> zvar = ...
    >>> f3d  = Fun3D ( func , xvar = xvar , yvar = yvar , zvar = zvar ) 
    """
    def __init__ ( self , fun , xvar , yvar , zvar , name = '' ) :

        self.__argfun == fun        
        if   isinstance ( fun , num_types ) :
            value = float ( fun ) 
            vname = 'C(%s,%s,%s)' %  ( value , xvar.name , yvar.name , zvar.name ) 
            fun   = Ostap.MoreRooFit.Constant ( self.generate_name ( vname , name ) , 'constant' , value , xvar , yvar , zvar ) 
            
        elif isinstance ( fun , AFUN1     ) : fun = fun.fun
            
        assert xvar and isinstance ( xvar , ROOT.RooAbsReal ) , "'xvar' must be ROOT.RooAbsReal"
        assert yvar and isinstance ( yvar , ROOT.RooAbsReal ) , "'yvar' must be ROOT.RooAbsReal"
        assert zvar and isinstance ( zvar , ROOT.RooAbsReal ) , "'zvar' must be ROOT.RooAbsReal"
        assert fun  and isinstance ( fun  , ROOT.RooAbsReal ) , "'fun'  must be ROOT.RooAbsReal"

        assert not xvar is yvar, "xvar and yvar must be different!"
        assert not xvar is zvar, "xvar and zvar must be different!"
        assert not yvar is zvar, "yvar and zvar must be different!"
        
        if fun is xvar : fun = Ostap.MoreRooFit.Id ( "Id_%s" % xvar.name , "Id:%s" % xvar.title , xvar )            
        if fun is yvar : fun = Ostap.MoreRooFit.Id ( "Id_%s" % yvar.name , "Id:%s" % yvar.title , yvar )            
        if fun is zvar : fun = Ostap.MoreRooFit.Id ( "Id_%s" % zvar.name , "Id:%s" % zvar.title , zvar )            
        
        if not name : name = 'Fun3D_%s' % fun.GetName() 

        FUN3.__init__ ( self , name , xvar = xvar , yvar = yvar , zvar = zvar )

        if not isinstance ( fun , ( ROOT.RooConstVar , Ostap.MoreRooFit.Id ) ) :            
            if not self.xvar in fun.getParameters ( 0 ) and not self.xvar is fun :
                self.warning ("Function does not depends on xvar=%s" % self.xvar.name )                
            if not self.yvar in fun.getParameters ( 0 ) and not self.yvar is fun :
                self.warning ("Function does not depends on yvar=%s" % self.yvar.name )                
            if not self.zvar in fun.getParameters ( 0 ) and not self.zvar is fun :
                self.warning ("Function does not depends on zvar=%s" % self.zvar.name ) 

        self.fun = fun

        self.config = {
            'fun'  : self.fun   ,
            'xvar' : self.xvar  ,
            'yvar' : self.yvar  ,
            'zvar' : self.zvar  ,
            'name' : self.name  ,            
            }
        
        self.checked_keys.add  ( 'fun'  )        
        self.checked_keys.add  ( 'xvar' ) 
        self.checked_keys.add  ( 'yvar' ) 
        self.checked_keys.add  ( 'zvar' ) 

# =============================================================================
## @class Fun1op
#  Helper function to define operations 
class Fun1op(FUN1) :
    """Helper function to define operations
    """    
    def __init__  ( self         ,
                    fun1         ,
                    fun2         ,
                    operation    , 
                    xvar         ,
                    name         ,
                    strname = '' ) :
        
        ## initialize the base class 
        FUN1.__init__ ( self , name , xvar = xvar )
        
        self.__fun1 = fun1
        self.__fun2 = fun2

        ## Fix for new problem/feature appearing in dev3 at 20221/08/11 
        if operation is Ostap.MoreRooFit.Addition :
            if isinstance ( fun1 , ROOT.RooAbsRealLValue ) :
                fun1 = Ostap.MoreRooFit.Id ( self.new_roo_name ( 'Id_' + fun1.name ) , fun1 )
            if isinstance ( fun2 , ROOT.RooAbsRealLValue ) :
                fun2 = Ostap.MoreRooFit.Id ( self.new_roo_name ( 'Id_' + fun2.name ) , fun2 )
                
        if   isinstance ( fun1 , AFUN1           ) : self.__raw1  = fun1.fun
        elif isinstance ( fun1 , ROOT.RooAbsReal ) : self.__raw1  = fun1
        elif isinstance ( fun1 , constant_types  ) : self.__raw1  = ROOT.RooConst ( float ( fun1 ) )
        else : raise TypeError ( "Invalid 'fun1': %s/%s" % ( self.fun1 , type ( self.fun1 ) ) )
        
        if   isinstance ( fun2 , AFUN1           ) : self.__raw2  = fun2.fun
        elif isinstance ( fun2 , ROOT.RooAbsReal ) : self.__raw2  = fun2
        elif isinstance ( fun2 , constant_types  ) : self.__raw2  = ROOT.RooConst ( float ( fun2 ) )
        else : raise TypeError ( "Invalid 'fun2': %s/%s" % ( self.fun2 , type ( self.fun2 ) ) )
        
        self.__strname   = strname 
        self.__operation = operation

        ## create the actual function 
        self.fun         = operation ( self.raw1 , self.raw2 ) 
            
        self.config = {
            'fun1'      : self.fun1      ,
            'fun2'      : self.fun2      ,
            'operation' : self.operation ,
            'xvar'      : self.xvar      ,
            'name'      : self.name      ,            
            'strname'   : self.strname   ,            
            }
        
    def __str__ ( self ) :
        return self.__strname if self.__strname else AFUN1.__str__ ( self ) 
    __repr__ = __str__
    
    @property
    def fun1 ( self ) :
        """'fun1' : the first argument"""
        return self.__fun1
    @property
    def fun2 ( self ) :
        """'fun2' : the second argument"""
        return self.__fun2

    @property
    def raw1 ( self ) :
        """'raw1' : true 'raw' object for the first argument"""
        return self.__raw1
    @property
    def raw2 ( self ) :
        """'raw2' : true 'raw' object for the second argument"""
        return self.__raw2
    
    @property
    def operation ( self ) :
        """'operation' : the actual operation"""
        return self.__operation
    @property
    def strname   ( self ) :
        """'strname' : name for __str__ function"""
        return self.__strname

# =============================================================================
## @class Fun2op
#  Helper function to define operations 
class Fun2op(FUN2) :
    """Helper function to define operations
    """    
    def __init__  ( self         ,
                    fun1         ,
                    fun2         ,
                    operation    , 
                    xvar         ,
                    yvar         ,
                    name         ,
                    strname = '' ) :
        
        ## initialize the base class 
        FUN2.__init__ ( self , name , xvar = xvar , yvar = yvar )
        
        self.__fun1 = fun1
        self.__fun2 = fun2
        
        ## Fix for new problem/feature appearing in dev3 at 20221/08/11 
        if operation is Ostap.MoreRooFit.Addition :
            if isinstance ( fun1 , ROOT.RooAbsRealLValue ) :
                fun1 = Ostap.MoreRooFit.Id ( self.new_roo_name ( 'Id_' + fun1.name ) , fun1 )
            if isinstance ( fun2 , ROOT.RooAbsRealLValue ) :
                fun2 = Ostap.MoreRooFit.Id ( self.new_roo_name ( 'Id_' + fun2.name ) , fun2 )
                
        if   isinstance ( fun1 , AFUN1           ) : self.__raw1  = fun1.fun
        elif isinstance ( fun1 , ROOT.RooAbsReal ) : self.__raw1  = fun1
        elif isinstance ( fun1 , constant_types  ) : self.__raw1  = ROOT.RooConst ( float ( fun1 ) )
        else : raise TypeError ( "Invalid 'fun1': %s/%s" % ( self.fun1 , type ( self.fun1 ) ) )
        
        if   isinstance ( fun2 , AFUN1           ) : self.__raw2  = fun2.fun
        elif isinstance ( fun2 , ROOT.RooAbsReal ) : self.__raw2  = fun2
        elif isinstance ( fun2 , constant_types  ) : self.__raw2  = ROOT.RooConst ( float ( fun2 ) )
        else : raise TypeError ( "Invalid 'fun2': %s/%s" % ( self.fun2 , type ( self.fun2 ) ) )
        
        self.__strname   = strname 
        self.__operation = operation

        ## create the actual function 
        self.fun         = operation ( self.raw1 , self.raw2 ) 
            
        self.config = {
            'fun1'      : self.fun1      ,
            'fun2'      : self.fun2      ,
            'operation' : self.operation ,
            'xvar'      : self.xvar      ,
            'yvar'      : self.yvar      ,
            'name'      : self.name      ,            
            'strname'   : self.strname   ,            
            }
        
    def __str__ ( self ) : return self.__strname if self.__strname else AFUN2.__str__ ( self ) 
    __repr__ = __str__
    
    @property
    def fun1 ( self ) :
        """'fun1' : the first argument"""
        return self.__fun1
    @property
    def fun2 ( self ) :
        """'fun2' : the second argument"""
        return self.__fun2

    @property
    def raw1 ( self ) :
        """'raw1' : true 'raw' object for the first argument"""
        return self.__raw1
    @property
    def raw2 ( self ) :
        """'raw2' : true 'raw' object for the second argument"""
        return self.__raw2
        
    @property
    def operation ( self ) :
        """'operation' : the actual operation"""
        return self.__operation
    @property
    def strname   ( self ) :
        """'strname' : name for __str__ function"""
        return self.__strname

# =============================================================================
## @class Fun3op
#  Helper function to define operations 
class Fun3op(FUN3) :
    """Helper function to define operations
    """    
    def __init__  ( self         ,
                    fun1         ,
                    fun2         ,
                    operation    , 
                    xvar         ,
                    yvar         ,
                    zvar         ,
                    name         ,
                    strname = '' ) :
        
        ## initialize the base class 
        FUN3.__init__ ( self , name , xvar = xvar , yvar = yvar , zvar = zvar )
        
        self.__fun1 = fun1
        self.__fun2 = fun2
        
        ## Fix for new problem/feature appearing in dev3 at 20221/08/11 
        if operation is Ostap.MoreRooFit.Addition :
            if isinstance ( fun1 , ROOT.RooAbsRealLValue ) :
                fun1 = Ostap.MoreRooFit.Id ( self.new_roo_name ( 'Id_' + fun1.name ) , fun1 )
            if isinstance ( fun2 , ROOT.RooAbsRealLValue ) :
                fun2 = Ostap.MoreRooFit.Id ( self.new_roo_name ( 'Id_' + fun2.name ) , fun2 )

        if   isinstance ( fun1 , AFUN1           ) : self.__raw1  = fun1.fun
        elif isinstance ( fun1 , ROOT.RooAbsReal ) : self.__raw1  = fun1
        elif isinstance ( fun1 , constant_types  ) : self.__raw1  = ROOT.RooConst ( float ( fun1 ) )
        else : raise TypeError ( "Invalid 'fun1': %s/%s" % ( self.fun1 , type ( self.fun1 ) ) )
        
        if   isinstance ( fun2 , AFUN1           ) : self.__raw2  = fun2.fun
        elif isinstance ( fun2 , ROOT.RooAbsReal ) : self.__raw2  = fun2
        elif isinstance ( fun2 , constant_types  ) : self.__raw2  = ROOT.RooConst ( float ( fun2 ) )
        else : raise TypeError ( "Invalid 'fun2': %s/%s" % ( self.fun2 , type ( self.fun2 ) ) )
        
        self.__strname   = strname 
        self.__operation = operation

        ## create the actual function 
        self.fun         = operation ( self.raw1 , self.raw2 ) 
            
        self.config = {
            'fun1'      : self.fun1      ,
            'fun2'      : self.fun2      ,
            'operation' : self.operation , 
            'xvar'      : self.xvar      ,
            'yvar'      : self.yvar      ,
            'zvar'      : self.zvar      ,
            'name'      : self.name      ,            
            'strname'   : self.strname   ,            
            }
        
    def __str__ ( self ) : return self.__strname if self.__strname else AFUN3.__str__ ( self ) 
    __repr__ = __str__
    
    @property
    def fun1 ( self ) :
        """'fun1' : the first argument"""
        return self.__fun1
    @property
    def fun2 ( self ) :
        """'fun2' : the second argument"""
        return self.__fun2
    @property
    
    def raw1 ( self ) :
        """'raw1' : true 'raw' object for the first argument"""
        return self.__raw1
    @property
    def raw2 ( self ) :
        """'raw2' : true 'raw' object for the second argument"""
        return self.__raw2

    @property
    def operation ( self ) :
        """'operation' : the actual operation"""
        return self.__operation
    @property
    def strname   ( self ) :
        """'strname' : name for __str__ function"""
        return self.__strname

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
## make some popular special function from (presumably) 'short' description
#  @param fun  the type of background function/PDF
#  @param name the name of background function/PDF
#  @param xvar the observable
#  Possible values for <code>fun</code>:
#  - any Ostap/PDF      : PDF will be copied or cloned  
#  - any RooAbsReal     : <code>Fun1D</code> 
#  - 'pN', 'polN' , 'polyN'               : <code>BernsteinPoly(power=N)</code>
#  - 'iN', 'incN' , 'incrN','increasingN' : <code>MonotonicPoly(power=N,increasing=True)</code>
#  - 'dN', 'decN' , 'decrN','decreasingN' : <code>MonotonicPoly(power=N,increasing=False)</code>
#  - 'cxN' , 'convexN'                    : <code>ConvexOnlyPoly(power=N,convex=True)</code>
#  - 'cvN' , 'concaveN'                   : <code> ConvexOnlyPoly(power=N,convex=False)</code>
#  - 'roopolyN','rpN'                     : <code>RooPoly(power=N)</code>> 
def make_fun1 ( fun , name , xvar , logger = None , **kwargs ) :
    """Make some popular special function from presumably 'short' description
    Possible values for <code>fun</code>:
    - any Ostap/PDF      : PDF will be copied or cloned  
    - any RooAbsReal     : `Fun1D` 
    - 'bN' , 'pN', 'polN' , 'polyN'        : `BernsteinPoly(power=N)`
    - 'iN', 'incN' , 'incrN','increasingN' : `MonotonicPoly(power=N,increasing=True)`
    - 'dN', 'decN' , 'decrN','decreasingN' : `MonotonicPoly(power=N,increasing=False)`
    - 'cxN' , 'convexN'                    : `ConvexOnlyPoly(power=N,convex=True)`
    - 'cvN' , 'concaveN'                   : `ConvexOnlyPoly(power=N,convex=False)`
    - 'roopolyN','rpN', 'rN'               : `RooPoly(power=N)`
    """

    if not logger : logger = globals()['logger'] 

    prefix = kwargs.pop ( 'prefix' , '' )
    suffix = kwargs.pop ( 'suffix' , '' )

    from ostap.fitting.utils import make_name 
    name   = make_name  ( prefix , name , suffix )

    result = None

    ## some Ostap-based functions ?
    if   isinstance ( fun , AFUN1 ) : 
        
        ## return the same model
        if   xvar is bkg.xvar and not  kwargs :
            result = fun 
            logger.debug ( 'make_fun1: %s model is copied to %s' % ( bkg , model ) )
        else :           ## make a clone : 
            result = fun.clone ( name = name , xvar = xvar , **kwargs )
            logger.debug ( 'make_fun1: %s model is cloned to %s' % ( bkg , model ) )

    elif isinstance ( fun , ROOT.RooAbsReal ) and xvar and isinstance ( xvar , ROOT.RooAbsReal ) :
        result = Fun1D ( fun , xvar = xvar , name = name )

    elif isinstance ( fun , string_types ) :

        sfun = fun.strip().lower() 

        import re        
        ok = re.search ( r'(poly|pol|p|b)(( *)|(_*))(?P<degree>\d)' , sfun , re.IGNORECASE )
        if ok :
            degree = int ( ok.group ( 'degree' ) ) 
            from ostap.fitting.roofuncs import BernsteinPoly as BP 
            result =  BP ( name , xvar = xvar , power = degree )

        ok = re.search ( r'(increasing|increase|incr|inc|i)(( *)|(_*))(?P<degree>\d)' , sfun , re.IGNORECASE )
        if ok : 
            degree = int ( ok.group ( 'degree' ) )
            from ostap.fitting.roofuncs import MonotonicPoly as MP 
            result = MP ( name , xvar = xvar , power = degree , increasing = True )

        ok = re.search ( r'(decreasing|decrease|decr|dec|d)(( *)|(_*))(?P<degree>\d)' , sfun, re.IGNORECASE )
        if ok : 
            degree = int ( ok.group ( 'degree' ) )
            from ostap.fitting.roofuncs import MonotonicPoly as MP 
            result = MP ( name , xvar = xvar , power = degree , increasing = False )

        ok = re.search ( r'(convex|cx)(( *)|(_*))(?P<degree>\d)' , sfun  , re.IGNORECASE )
        if ok :
            degree = int ( ok.group ( 'degree' ) )
            from ostap.fitting.roofuncs import ConvexOnlyPoly as CO 
            result = CO ( name , xvar = xvar , power = degree , convex = True  )

        ok = re.search ( r'(concave|cv)(( *)|(_*))(?P<degree>\d)' , sfun, re.IGNORECASE )
        if ok :
            degree = int ( ok.group ( 'degree' ) )
            from ostap.fitting.roofuncs import ConvexOnlyPoly as CO 
            result = CO ( name , xvar = xvar , power = degree , convex = False  )
        
        ok = re.search ( r'(roopoly|rp|r)(( *)|(_*))(?P<degree>\d)' , sfun, re.IGNORECASE )
        if ok :
            degree = int ( ok.group ( 'degree' ) )
            from ostap.fitting.roofuncs import RooPoly as RP
            result = RP ( name , xvar = xvar , power = degree )

    if result :
        logger.debug ( 'make_fun1: created functiob is %s' % result  ) 
        return result 

    raise  TypeError("make_fun1: Wrong type of 'fun' object: %s/%s " % ( fun , type ( fun ) ) ) 


    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    
# =============================================================================
#                                                                       The END 
# =============================================================================
