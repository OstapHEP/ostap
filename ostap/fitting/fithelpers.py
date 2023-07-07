#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/fithelpers.py
#  Set of useful helpers to build various fit models 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
# =============================================================================
"""Set of useful technical utilities to build various fit models"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2018-08-14"
__all__     = (
    ##
    'VarMaker'          , ## Helper bade class that allow storage of newly created RooFit objects
    'FitHelper'         , ## Helper base class for fit-relate objects
    'ConfigReducer'     , ## Helper base class for pickliing/unpickling via  reduce
    ## 
    'XVar'              , ## helper MIXIN to deal with x-variable
    'YVar'              , ## helper MIXIN to deal with y-variable
    'ZVar'              , ## helper MIXIN to deal with z-variable
    #
    "H1D_dset"          , ## 1D-histogram to RooDataHist converter 
    "H2D_dset"          , ## 2D-histogram to RooDataHist converter 
    "H3D_dset"          , ## 3D-histogram to RooDataHist converter
    #
    'Phases'            , ##  helper class for Ostap polynomial/PDFs
    'ParamsPoly'        , ##  helper class for RooFit polynomials
    'ShiftScalePoly'    , ##  helper class for RooFit polynomials
    #
    "NameDuplicates"    , ## allow/disallow name duplicates
    'SETPARS'           , ## context manager to keep/preserve parameters 
    )
# =============================================================================
from   ostap.core.meta_info    import root_info 
from   ostap.core.core         import ( Ostap, rootID, VE,
                                        items_loop, isequal , roo_silent ) 
from   ostap.core.ostap_types  import ( string_types   , num_types   ,
                                        integer_types  , list_types  , 
                                        is_good_number , is_integer  ,
                                        sequence_types , sized_types )
from   ostap.utils.utils       import make_iterable
from   ostap.fitting.variables import SETVAR, make_formula 
from   ostap.fitting.utils     import make_name, numcpu, ncpu, get_i  
from   ostap.fitting.roocmdarg import check_arg 
from   ostap.math.random_ext   import ve_gauss, poisson
import ROOT, sys, random
# =============================================================================
from   ostap.logger.logger   import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.fithelpers')
else                       : logger = getLogger ( __name__                  )
# =============================================================================
try                : from string import            ascii_letters
except ImportError : from string import letters as ascii_letters
# =============================================================================
_py2 = sys.version_info.major < 3 
# =============================================================================
## @class NameDuplicates
#  Are name duplicates allowed?
#  @code
#  if NameDuplicates.allowed()
#  ...
#  NameDuplicates.allow ( True ) 
#  @endcode
#  @attention allowing name duplication is a dangerous action!
class NameDuplicates(object) :
    """ Are name duplicates allowed?
    >>> if NameDuplicates.allowed() : .... 
    ...
    >>> NameDuplicates.allow ( True )
    - Attention:  allowing name duplication is a dangerous action!
    """
    __allowed = False
    @classmethod
    def allowed ( cls ) : return cls.__allowed
    @classmethod
    def allow   ( cls , allowed ) :
        cls.__allowed = True if allowed else False
    # =========================================================================
    def __init__   ( self , allow  ) :
        self.__allow  = True if allow else False
        self.__backup = self.allowed()
    def __enter__  ( self ) :        
        self.__backup = self.allowed()
        self.allow ( self.__allow )
        return self
    def __exit__   ( self , *_     ) :
        self.allow ( self.__backup )
# =============================================================================
## helper factory function
def config_factory ( klass , config ) :
    """Helper factory function, used for unpickling"""
    with NameDuplicates ( True ) : 
        return klass ( **config ) 
# =============================================================================
## @class ConfigReducer
#  Helper base class for pickling/unpickling FUNs/PDFs
class ConfigReducer(object) :
    """Helper base class for pickling/unpickling FUNs/PDFs
    """
    def __init__ ( self , **kwargs ) :
        self.__config = {}
        self.__config.update ( kwargs )
        
    ## pickling via reduce 
    def __reduce__ ( self ) :
        if _py2 : return config_factory , ( type ( self ) , self.config, )
        else    : return type ( self ).factory , ( self.config, )

    ## factory method 
    @classmethod
    def factory ( klass , config ) :
        """Factory method, used for unpickling"""
        with NameDuplicates ( True ) :
            return klass ( **config )
        
    @property
    def config ( self ) :
        """The full configuration info for the FUN/PDF"""
        conf = {}
        conf.update ( self.__config )
        return conf
    @config.setter
    def config ( self , value ) :
        conf = {}
        conf.update ( value )
        self.__config = conf

    ## oversimplified placeholder for the clone method 
    def clone ( self , **kwargs ) :
        """Oversimplified placeolder for clone method
        - actual cloning for PDFs are a bit more sophistocated 
        """        
        ## get config 
        conf = {}
        conf.update ( self.config ) 
        conf.update ( kwargs      )
        ## 
        return self.factory ( conf ) 

# =============================================================================
## @class VarMaker
#  Helper class that allows implement several purely  technical methods:
#   - creation of <code>ROOT.RooRealVar objects</code>
#   - store newly created RooFit objects
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2018-07-14
class VarMaker (object) :
    """Helper class that allows implement several purely  technical methods:
    - creation of <code>ROOT.RooRealVar objects</code>
    - store newly created RooFit objects
    """
    __pdf_names = set () ## all known/generated PDF/FUN names 
    __var_names = set () ## all known/generated VAR names
    __loggers   =     {} ## the list of local loggers  
    __numnames  = 0
    
    ## @attention ensure that important attributes are available even before __init__
    def __new__( cls, *args, **kwargs):
        if _py2 : obj = super(VarMaker, cls).__new__( cls , *args , **kwargs )
        else    : obj = super(VarMaker, cls).__new__( cls )
        ##
        obj.__aux_keep     = []     ## ATTENTION!        
        obj.__name        = None    ## ATTENTION!
        obj.__local_names = set()
        return obj
    
    ##  produce ERROR    message using the local logger 
    def fatal   ( self , message , *args , **kwargs ) :
        """Produce FATGAL   message using the local logger"""
        return self.logger.fatal   ( message , *args , **kwargs )
    ##  produce ERROR    message using the local logger 
    def error   ( self , message , *args , **kwargs ) :
        """Produce ERROR    message using the local logger"""
        return self.logger.error   ( message , *args , **kwargs )
    ##  produce WARNIING message using the local logger 
    def warning ( self , message , *args , **kwargs ) :
        """Produce WARNIING message using the local logger"""
        return self.logger.warning ( message , *args , **kwargs )
    ##  produce INFO     message using the local logger 
    def info    ( self , message , *args , **kwargs ) :
        """Produce INFO     message using the local logger"""
        return self.logger.info    ( message , *args , **kwargs )
    ##  produce DEBUG    message using the local logger 
    def debug   ( self , message , *args , **kwargs ) :
        """Produce DEBUG    message using the local logger"""
        return self.logger.debug   ( message , *args , **kwargs )
    ##  produce VERBOSE  message using the local logger 
    def verbose ( self , message , *args , **kwargs ) :
        """Produce VERBOSE  message using the local logger"""
        return self.logger.verbose ( message , *args , **kwargs )

    # =========================================================================
    ## The list of objects that are kept 
    @property
    def aux_keep ( self ) :
        """'aux_keep' -  the list of objects to be kept by this PDF/FUN"""
        return self.__aux_keep
    # =========================================================================
    ## The actual local logger object
    @property
    def logger   ( self ) :
        """'logger': get the local logger object"""
        name    = self.name
        logname = str ( self.__class__.__name__ ) 
        if name : logname += '(' + name + ')'
        _logger = self.__loggers.get ( logname , None )
        if not _logger :
            _logger  = getLogger ( logname )
            self.__loggers [ logname ] = _logger
        return _logger
    @property
    def name ( self ) :
        """The name of the object"""
        return self.__name if self.__name else '' 
    @name.setter
    def name ( self , value ) :
        assert isinstance ( value , str ) , "'name' must  be a string, %s/%s is given" % ( value , type ( value ) )
        if self.__name == value : return 
        if value in self.__pdf_names and not NameDuplicates.allowed() :
            self.warning ( 'The name "%s" for PDF/FUN already defined!' % value )
        self.__pdf_names.add ( value )     
        self.__name = value

    # =============================================================================
    ## generate some unique name for PDF/FUN and objects
    @staticmethod 
    def generate_name ( prefix = '', suffix = '' , name = '' ) :
        """Generate some unique name for PDF/FUN and objects
        """

        def good_name ( nam ) :
            if  not nam                      : return False 
            elif nam in VarMaker.__pdf_names : return False
            elif nam in VarMaker.__var_names : return False
            return True 

        ## is name already in prefix/suffix? 
        name_twice = name and ( prefix or suffix ) and ( not name in prefix ) and ( not name in suffix ) 

        if name_twice : newname = make_name ( prefix , name , suffix )
        else          : newname = make_name ( prefix , ''   , suffix )        
        VarMaker.__numnames += 1

        while not good_name ( newname ) :            
            part = ''.join ( ( random.choice ( ascii_letters ) for i in range ( 6 ) )  )
            if name_twice : newname = make_name ( prefix , name + '_' + part , suffix )
            else          : newname = make_name ( prefix ,              part , suffix )            
            VarMaker.__numnames += 1

        return newname
    
    # =============================================================================
    ## generate some unique name for <code>RooFit</code>
    #  @see TNamed 
    #  @see RooNameReg 
    #  @see RooAbsArg 
    @staticmethod
    def roo_name ( prefix = 'roo' , suffix = '' , name = '' ) :
        """Generate some unique name for <code>RooFit</code>
        - see `ROOT.RooNameReg` 
        - see `ROOT.TNamed`
        - see `ROOT.RooAbsArg`
        """
        
        regname = ROOT.RooNameReg.instance()
        
        ## RooFit does no tline coma in the names
        ## RooAbsCollection::selectByName treat commas special case 
        
        to_replace = ( ( ','  ,'_'      ) ,
                       ( ';'  ,'_'      ) ,
                       ( ':'  ,'_'      ) ,
                       ( '?'  ,'_'      ) ,
                       ( '/'  ,'_'      ) ,
                       ( '\\' ,'_'      ) ,
                       ( ' '  ,'_'      ) ,
                       ( '+'  ,'_plus_' ) ,
                       ( '______'  ,'_' ) ,
                       ( '_____'   ,'_' ) ,
                       ( '____'    ,'_' ) ,
                       ( '___'     ,'_' ) ,
                       ( '__'      ,'_' ) ,
                       ( '__'      ,'_' ) ,
                       ( '__'      ,'_' ) )
                       
        def rootify  ( nam ) :
            for a , b in to_replace : nam = nam.replace ( a , b )
            return nam 

        prefix = rootify ( prefix.strip ( '_ ' ) ) 
        suffix = rootify ( suffix.strip ( '_ ' ) ) 
        name   = rootify ( name  .strip ( '_ ' ) ) 

        def good_name ( nam ) :
            if  not nam                      : return False 
            elif nam in VarMaker.__pdf_names : return False
            elif nam in VarMaker.__var_names : return False
            elif regname.known ( nam )       : return False
            return True 

        ## is name already in prefix/suffix? 
        name_twice = name and ( prefix or suffix ) and ( not name in prefix ) and ( not name in suffix ) 

        if name_twice : newname = make_name ( prefix , name , suffix )
        else          : newname = make_name ( prefix , ''   , suffix )
        
        VarMaker.__numnames += 1
        while not good_name ( newname ) :            
            part = ''.join ( ( random.choice ( ascii_letters ) for i in range ( 6 ) )  )
            if name_twice : newname = make_name ( prefix , name + '_' + part , suffix )
            else          : newname = make_name ( prefix ,              part , suffix )            
            VarMaker.__numnames += 1

        return newname

    # =============================================================================
    ## generate unique name 
    def new_name ( self , prefix = '' , suffix = '' ) :
        """generate unique name"""
        return self.generate_name ( prefix = prefix , suffix = suffix , name = self.name )  

    # =============================================================================
    ## generate unique name 
    def new_roo_name ( self , prefix = '' , suffix = '' ) :
        """generate unique name"""
        return self.roo_name ( prefix = prefix , suffix = suffix , name = self.name )  
    
    # =============================================================================
    ## create/modify  the variable
    #  Helper function for creation/modification/adjustment of variable
    #  @code
    #  v = self.make_var ( 10   , 'myvar' , 'mycomment' )
    #  v = self.make_var ( 10   , 'myvar' , 'mycomment' , ' ,     -1 , 1 )
    #  v = self.make_var ( 10   , 'myvar' , 'mycomment' , ' , 0 , -1 , 1 )
    #  v = self.make_var ( None , 'myvar' , 'mycomment' , ' , 0 , -1 , 1 )
    #  v = self.make_var ( None , 'myvar' , 'mycomment' , 10 , 0 , -1 , 1 )
    #  v = self.make_var ( v    , 'myvar' , 'mycomment' , 10 )
    #  @endcode
    #  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
    #  @date 2013-12-01
    def make_var ( self           ,
                   var            ,
                   name           , ## name 
                   title   = ''   , ## title 
                   fix     = None , *args ) :
        """Make/modify  the variable:
        
        v = self.make_var ( 10   , 'myvar' , 'mycomment' )
        v = self.make_var ( 10   , 'myvar' , 'mycomment' , ' ,     -1 , 1 )
        v = self.make_var ( 10   , 'myvar' , 'mycomment' , ' , 0 , -1 , 1 )
        v = self.make_var ( None , 'myvar' , 'mycomment' , ' , 0 , -1 , 1 )
        v = self.make_var ( None , 'myvar' , 'mycomment' , 10 , 0 , -1 , 1 )
        v = self.make_var ( v    , 'myvar' , 'mycomment' , 10 )        
        """
        # var = ( value )
        # var = ( min , max )
        # var = ( value , min , max )
        
        if   isinstance ( var , ROOT.RooAbsReal ) : pass 
        elif var is None : 
            assert name and args , "make_var: 'name' and 'args' must be specified when 'var' is None!"
            var = name , title if title else name            
        elif isinstance ( var , num_types ) :
            assert name , "make_var: 'name' must be specified when 'var' is of numerical type!"
            var = name , title if title else name , float ( var ) 
            
        ## convert sequence of values 
        if   isinstance ( var , ROOT.RooAbsReal ) : pass 
        elif isinstance ( var , sequence_types  ) : var = tuple ( var )

        if   var and isinstance ( var , ROOT.RooAbsReal ) : pass 
        elif var and isinstance ( var , tuple ) :

            ## content of args 

            ## at most four argumentss:  (value, min, max, unit)
            assert len ( args ) <= 4 , "make_var: invalid length of 'args' %s" % len ( args ) 
            
            ## strip out the uniuts from args 
            vunit = ''
            vargs = args 
            if args and isinstance ( args[-1] , string_types ) :
                vunit = args[ -1]
                vargs = args[:-1]
                
            ## convert to numbers 
            vargs = tuple ( float ( v ) for v in vargs )

            ## content of var 
        
            vvars = var   [0:]

            ## pickup the name (if specified...) 
            vname = vvars [0]
            if isinstance ( vname , string_types ) :
                assert vname and ( ( not name ) or vname == name ) , \
                       "make_var: invalid specification of 'name' : %s/%s" % ( vname , name )
                name  = vname 
                vvars = vvars [1:]

            assert name, "make_var: 'name' must be specified for creation of new variable!"

            if vvars :
                ## pickup the title (if specified...) 
                vtitle = vvars [ 0 ] 
                if isinstance ( vtitle , string_types ) :
                    if title and vtitle and title != vtitle : 
                        self.warning ("make_var('%s'): title is redefined: %s -> %s" % ( name , title , vtitle ) ) 
                    if vtitle : title = vtitle
                    vvars = vvars [1:]
                    
            ## at most four  argumentss:  (value, min, max, unit )
            assert len ( vvars ) <= 4 , "make_var('%s'): invalid length of 'vvars' %s" % ( name , len ( vvars ) ) 

            ## strip out the units from vvars 
            unit = '' 
            if vvars and isinstance ( vvars [ -1 ] , string_types ) :
                unit = vvars [ -1]
                vvar = vvars [:-1]

            ## convert to numbers 
            vvars = tuple ( float ( v ) for v in vvars )

            if unit and vunit and unit != unit :
                self.warning ("make_var(%s): unit is redefined: %s -> %s" % ( name , vunit , unit ) ) 
            elif vunit : unit = vunit
            
            len_vvars = len ( vvars )
            len_vargs = len ( vargs )

            vvars_ = vvars
            vargs_ = vargs
            vvars  = tuple ( sorted ( vvars ) ) 
            vargs  = tuple ( sorted ( vargs ) ) 
            ## attention! for 3-tupel the order is (value, min, max)
            if 3 == len ( vvars ) : vvars = vvars [ 1 ] , vvars [ 0 ] , vvars [ 2 ] 
            if 3 == len ( vargs ) : vargs = vargs [ 1 ] , vargs [ 0 ] , vargs [ 2 ] 
            if vvars != vvars_ : self.warning ("make_var('%s'):  %s -> %s " % ( name , str ( vvars_ ) , str ( vvars ) ) )  
            if vargs != vargs_ : self.warning ("make_var('%s'):  %s -> %s " % ( name , str ( vargs_ ) , str ( vargs ) ) )  

            ## there should be at least one useful number! 
            assert 1<= len_vvars + len_vargs , "make_var: empty 'vvars' and 'vargs'!"

            fixed = False
            
            if   0 == len_vvars :
                ## OK 
                params = vargs + ( unit , )
                
            elif 1 == len_vvars and 0 == len_vargs :
                ## OK: use only vvar 
                params = vvars + ( unit , )
                
            elif 1 == len_vvars and 2 == len_vargs and vargs [ 0 ] <= vvars [ 0 ] <= vargs [ 1 ] :
                ## OK, get min/max from vargs 
                params = vvars + vargs + ( unit , )

            elif 1 == len_vvars and 2 == len_vargs :
                ## OK, get min/max from vargs 
                params = vvars + ( min ( vvars [ 0 ] , vargs [ 0 ] ) , max ( vvars[0] , vargs [ 1 ] ) , unit )
                fixed = True 

            elif 1 == len_vvars and 3 == len_vargs and vargs [ 1 ] <= vvars [ 0 ] <= vargs [ 2 ] :
                ## OK 
                params = vvars + ( min ( vvars [ 0 ] , vargs [ 1 ] ) , max ( vvars [ 0 ] , vargs [ 2 ]  ) , unit ) 
                
            elif 1 == len_vvars and 3 == len_vargs :
                ## FIX IT
                params = vvars + ( min ( vvars[0] , vargs [ 1 ] ) , max ( vvars[0] , vargs [ 2 ] ) , unit ) 
                fixed  = True 
                
            elif 1 == len_vvars :
                ## 'vargs" are ignored 
                params = vvars + ( unit , )
                fixed  = True
                
            elif 2 == len_vvars and 0 == len_vargs :
                ## OK: min/max specified via vvar
                params = vvars  + (  unit , ) 

            elif 2 == len_vvars and 1 == len_vargs and vvars [ 0 ] <= vargs [ 0 ] <= vvars [ 1 ] :
                ## get the value from vargs 
                params = vargs + vvars + ( unit , ) 

            elif 2 == len_vvars and 1 == len_vargs :
                ## get the value from vargs 
                params = vvars + ( unit , )
                fixed  = True 

            elif 2 == len_vvars and 2 == len_vargs :
                ## igore vargs 
                params = vvars + ( unit , ) 
                
            elif 2 == len_vvars and 3 == len_vargs and vvars [ 0 ] <= vargs [ 0 ] <= vvars [ 1 ] :
                ## get the value from vargs 
                params = vargs[:1]  + vvars + ( unit , )
            
            elif 2 == len_vvars :
                ## 'vargs" are ignored
                params = vvars + ( unit , )
                fixed  = True 

            elif 3 == len_vvars :
                ## ignore vargs 
                params = vvars + ( unit , )

            else : 
                
                args   = vvar + vargs
                vmin   = vmin ( *args )
                vmax   = vmax ( *args )
                params = vmin , vmax , unit
                fixed  = True
                
            if fixed : self.warning ( "make_var('%s'):  %s + %s -> %s" % ( name , str ( vvars_ ) , str ( vargs_ ) , str ( params[:-1] ) ) )
            else     : self.debug   ( "make_var('%s'):  %s + %s -> %s" % ( name , str ( vvars_ ) , str ( vargs_ ) , str ( params[:-1] ) ) )
            
            ## create the variable!
            var = ROOT.RooRealVar ( name , title , *params ) 
            self.debug ( "New variable is created: name:'%s',title:'%s', args=%s" % ( var.name , var.title , str ( params ) ) )
            
            if   fix is None                     : pass ## no actions 
            elif isinstance ( fix , bool       ) : var.setConstant ( fix )
            elif isinstance ( fix , num_types  ) :                
                vfix = float ( fix )
                if var.hasMin() and vfix < var.getMin () : self.warning ("make_var('%s'): fix value %s < min %s" % ( var.name , vfix , vmin ) ) 
                if var.hasMax() and vfix > var.getMax () : self.warning ("make_var('%s'): fix value %s > max %s" % ( var.name , vfix , vmax ) )
                var.setVal      ( vfix )
                value = float ( var ) 
                if ( value != vfix ) : self.warning ("make_var('%s'): the value for %s variable %s is not fixed at the value of %s" % ( var.name , value , vfix ) )                    
                var.setConstant ( True )
            else :
                self.warning ("make_var('%s'): ignore 'fix' of %s" % ( var.name , str ( fix ) ) )
                
            ## raise TypeError("make_var('%s'): Invalid 'var/fix' setting: %s/%s" % ( var.name , var , fix ) )                

        ## var is either newly created or provided as RooAbsArg 
        assert isinstance ( var , ROOT.RooAbsArg ) , \
               "make_var: invalid 'var' type! %s/%s" % ( var , type ( var ) ) 

        ## store the variable
        self.aux_keep.append ( var ) 
        
        return var

    # ==========================================================================
    ## check the possible name  duplication
    def var_name  ( self , name ) :
        """Check the possible name duplication
        """
        if name in self.__var_names and not NameDuplicates.allowed() :
            self.warning ( 'The variable name "%s" is already defined!' % name )
            
        self.__var_names.add   ( name )
        self.__local_names.add ( name )
        return name

    # =========================================================================
    ## delete the obejct
    #  - remove the registered names in storages
    #  - clear the local storage of names 
    def __del__ ( self ) :
        """Delete the obejct
        - remove the registered names in storages
        - clear the local storage of names 
        """
        
        if self.name and self.name in self.__pdf_names :
            self.__pdf_names.remove ( self.name ) 
            while self.__local_names :
                a = self.__local_names.pop ()
                if a in self.__var_names :
                    self.__var_names.remove  ( a ) 
        
    # =========================================================================
    ## create ROOT.RooFit.Binning from TAxis
    #  @see RooFit::Binning
    #  @see RooBinning
    #  @see TAxis
    #  @code
    #  axis = ...
    #  bins = bining ( axis ) 
    #  @endocode
    def binning ( self , axis , name = '' ) :
        """Create ROOT.RooFit.Binning from TAxis
        >>> axis = ...
        >>> bins = binning ( axis )
        """
        assert isinstance ( axis , ROOT.TAxis ),\
               'Invalid axis %s/%s' % ( axis , type ( axis ) )

        ## uniform binning?
        if not  axis.IsVariableBinSize() : 
            return ROOT.RooFit.Binning ( axis.GetNbins() , axis.GetXmin() , axis.GetXmax() )
        ##
        xbins = axis.GetXbins().GetArray()
        rb = ROOT.RooBinning ( axis.GetNbins() , xbins , name )
        ##
        self.aux_keep.append ( rb )
        ##
        return ROOT.RooFit.Binning ( rb )

# =============================================================================
## @class FitHelper
#  Helper base class for fit-rrleated objects
#  - storage of newly created RooFit objects
#  - loggin
#  - a lot of helper methdos for fitting 
class FitHelper(VarMaker) :
    """ Helper base class for fit-rrleated objects
    - storage of newly created RooFit objects
    - logging
    - a lot of helper methdos for fitting 
    """
            
    # =========================================================================
    ## technical function to parse arguments for <code>fitTo</code> and 
    #  <code>nll</code>  methods
    def parse_args ( self ,  dataset = None , *args , **kwargs ) :
        """Technical function to parse arguments for fitTo/nll/.. methods
        """
        _args = []
        for a in args :
            if not isinstance ( a , ROOT.RooCmdArg ) :
                self.error ( 'parse_args: unknown argument type %s/%s, skip' % ( a , type ( a ) ) )
            else : _args.append ( a ) 

        from ostap.plotting.fit_draw import keys  as drawing_options

        silent  = None
        verbose = None

        from ostap.core.core import cidict_fun as key_transform 
        transformed_draw_options = tuple ( key_transform ( k ) for k in drawing_options )
        
        for k , a in items_loop ( kwargs ) :
            
            key  = key_transform ( k )             
            
            ## skip "drawing" options 
            if   k   in drawing_options          : continue
            if   key in drawing_options          : continue
            elif k   in transformed_draw_options : continue
            elif key in transformed_draw_options : continue
            
            elif key in ( 'draw'                 ,
                          'drawoption'           ,
                          'drawoptions'          ) : continue 
            
            if   isinstance ( a , ROOT.RooCmdArg ) : _args.append ( a )            
            elif key in ( 'verbose' ,       ) and isinstance ( a , bool ) :
                
                if not verbose is None :
                    if a != verbose : 
                        self.warning ( 'parse_args: Redefine VERBOSE to %s' %  a ) 
                        verbose = a                        
                if not silent is None :
                    if a == silent :
                        self.warning ( 'parse_args: confusing VERBOSE/SILENT %s/%s' % ( a , silent ) )
                        silent = not a 
                _args.append ( ROOT.RooFit.Verbose (     a ) )

            elif key in ( 'silent'  , 'silence') and isinstance ( a , bool ) :
                
                if not silent is None :
                    if a != silent : 
                        self.warning ( 'parse_args: Redefine SILENT to %s' %  a ) 
                        verbose = a                        
                if not verbose is None :
                    if a == verbose :
                        self.warning ( 'parse_args: confusing SILENT/VERBOSE %s/%s' % ( a , verbose ) )
                        verbose = not a
                _args.append ( ROOT.RooFit.Verbose ( not a ) )
                
            elif key in ( 'strategy'         , 
                          'miniutstrategy'   ,
                          'strategyminuit'   ) and isinstance ( a , integer_types ) and 0 <= a <= 2 :
                
                _args.append ( ROOT.RooFit.Strategy (    a ) )
                
            elif key in ( 'printlevel'       ,
                          'minuitprint'      ,
                          'minuitlevel'      ) and isinstance ( a , integer_types ) and -1 <= a <= 3 :
                
                _args.append ( ROOT.RooFit.PrintLevel ( a ) )
                
            elif key in ( 'printevalerrors'  ,
                          'printerrors'      ,
                          'errorspront'      ) and isinstance ( a , integer_types ) and -1 <= a :
                
                _args.append ( ROOT.RooFit.PrintEvalErrors ( a ) )
                
            elif key in ( 'timer' , 'timing' ) and isinstance ( a , bool ) :
                
                _args.append ( ROOT.RooFit.Timer    ( a ) )
                
            elif key in ( 'warning' , 'warnings' ) and isinstance ( a , bool ) :
                
                _args.append ( ROOT.RooFit.Warnings ( a ) ) 
            
            elif key in ( 'sumw2'            ,
                          'sumw2err'         ,
                          'sumw2errs'        ,
                          'sumw2error'       ,
                          'sumw2errors'      ) and isinstance ( a , bool ) :
                
                if   dataset and isinstance ( dataset , ROOT.RooDataHist ) :
                    self.warning ('parse_args: no need in SumW2-flag for RooDataHist')
                elif a and dataset and     dataset.isWeighted()           : pass 
                elif a and dataset and not dataset.isWeighted()           :
                    self.warning ('parse_args: SumW2-flag is True  for non-weighted dataset')
                elif       dataset and not dataset.isWeighted() and not a : pass 
                elif       dataset and     dataset.isWeighted() and not a :
                    self.warning ('parse_args: SumW2-flag is False for     weighted dataset')                    

                _args.append (  ROOT.RooFit.SumW2Error( a ) )
                                    
            elif key in ( 'asymptotic'       ,
                          'asymptoticerr'    ,
                          'asymptoticerrs'   ,
                          'asymptoticerror'  ,
                          'asymptoticerrors' ) and isinstance ( a , bool ) and (6,19) <= root_info :
                
                if   dataset and isinstance ( dataset , ROOT.RooDataHist ) : pass
                    self.warning ('parse_args: no need in AsymptoticErrors-flag for RooDataHist')
                elif a and dataset and     dataset.isWeighted()           : pass 
                elif a and dataset and not dataset.isWeighted()           :
                    self.warning ('parse_args: AsymptoticError-flag is True  for non-weighted dataset')
                elif       dataset and not dataset.isWeighted() and not a : pass 
                elif       dataset and     dataset.isWeighted() and not a :
                    self.warning ('parse_args: AsymptoticError-flag is False for     weighted dataset')                    

                if a and root_info <  ( 6 , 27 ) :
                    self.warning ("'AsymptoticError=True' is buggy for this version of ROOT, (ROOT-PR-11282)")
                if a and root_info <= ( 6 , 26 ) :
                    self.warning ("'AsymptoticError=True' will crash if Title!=Name (ROOT-10668)")
                    
                _args.append (  ROOT.RooFit.AsymptoticError ( a ) )
                    
            elif key in ( 'batch'            ,
                          'batchmode'        ,
                          'modebatch'        ) and isinstance ( a , bool ) and  ( 6, 20) <= root_info :
                
                _args.append (  ROOT.RooFit.BatchMode ( a ) )
                
            elif key in ( 'extended' ,       ) and isinstance ( a , bool ) :
                
                _args.append   (  ROOT.RooFit.Extended ( a ) )
                
            elif key in ( 'cpu'              ,
                          'cpus'             ,
                          'ncpu'             ,
                          'ncpus'            ,
                          'numcpu'           ,
                          'numcpus'          ) and isinstance ( a , int ) and 1<= a :
                
                _args.append   (  ROOT.RooFit.NumCPU( a  ) )
                
            elif key in ( 'cpu'              ,
                          'cpus'             ,
                          'ncpu'             ,
                          'ncpus'            ,
                          'numcpu'           ,
                          'numcpus'          ) and \
                 isinstance ( a , list_types ) and 2 == len ( a )  and \
                 isinstance ( a[0] , integer_types ) and 1 <= a[1] and \
                 isinstance ( a[1] , integer_types ) and 0 <= a[1] <=3 :
                
                _args.append   (  ROOT.RooFit.NumCPU( a[0] ,  a[1] ) ) 
                
            elif key in ( 'range'            ,
                          'fitrange'         ,
                          'rangefit'         ,
                          'ranges'           ,
                          'fitranges'        ) and isinstance ( a , string_types ) :
                
                _args.append   (  ROOT.RooFit.Range ( a ) )
                
            elif key in ( 'range'            ,
                          'fitrange'         ,
                          'rangefit'         ) and isinstance ( a , list_types   ) and 2 == len ( a ) \
                 and isinstance ( a[0] ,  num_types ) \
                 and isinstance ( a[1] ,  num_types ) \
                 and a[0] < a[1]  :
                
                _args.append   (  ROOT.RooFit.Range ( a[0] , a[1] ) )
                
            elif key in ( 'minimizer'       ,
                          'minimiser'       ) and isinstance ( a , list_types   ) and 2 == len ( 2 ) \
                 and isinstance ( a[0] ,  string_types ) \
                 and isinstance ( a[1] ,  string_types ) :
                
                _args.append   (  ROOT.RooFit.Minimizer ( a[0] , a[1] ) )
                
            elif key in  ( 'hesse'    ,     ) and isinstance ( a , bool ) :
                
                _args.append   (  ROOT.RooFit.Hesse ( a )  )
                
            elif key in  ( 'initialhesse'    ,
                           'inithesse'       ,
                           'hesseinit'       ,
                           'hesseinitial'    ) and isinstance ( a , bool ) :
                
                _args.append   (  ROOT.RooFit.InitialHesse ( a )  )
                
            elif key in ( 'optimize'         ,
                          'optimise'         ) and isinstance ( a , integer_types  ) :
                
                _args.append   (  ROOT.RooFit.Optimize     ( a )  )
                
            elif key in ( 'minos'    ,       ) and isinstance ( a , bool           ) :
                
                _args.append   (  ROOT.RooFit.Minos        ( a )  )
                
            elif key in ( 'minos'    ,       ) and isinstance ( a , ROOT.RooArgSet ) :
                
                _args.append   (  ROOT.RooFit.Minos        ( a )  )

            elif key in ( 'minos'    ,       ) and isinstance ( a , ROOT.RooAbsReal ) and a in self : 

                vset = ROOT.RooArgSet( a    )
                self.aux_keep.append ( vset ) 
                _args.append   (  ROOT.RooFit.Minos        ( vset )  )
                
            elif key in ( 'minos'    ,       ) and isinstance ( a , string_types   ) \
                     and hasattr  ( self , 'params' ) and a in self.params ( dataset ) :
                
                _v = self.params()[ a ]
                _s = ROOT.RooArgSet ( _v )
                self.aux_keep.append ( _s ) 
                _args.append   (  ROOT.RooFit.Minos        ( _s )  )
                
            elif key in ( 'minos'  ,       ) and not isinstance ( a , string_types ) :

                _s     = ROOT.RooArgSet()
                _pars = self.params ( dataset ) if hasattr  ( self , 'params' ) else ROOT.RooArgSet() 
                for v in a :
                    if   v in _pars and isinstance ( v , string_types ):
                        _v = _pars [ v ] 
                        _s.add ( _v )
                    elif v in _pars and isinstance ( v , ROOT.RooAbsArg  ) :
                        _s.add (  v )
                    else :
                        self.error ( "Can not find %s in parameetrs" % v )

                self.aux_keep.append ( _s ) 
                _args.append   (  ROOT.RooFit.Minos ( _s )  )
                                
            elif key in ( 'save'     ,       ) and isinstance ( a , bool           ) :
                
                _args.append   (  ROOT.RooFit.Save         ( a )  )
                
            elif key in ( 'clone'            ,
                          'clonedata'        ,
                          'dataclone'        ) and isinstance ( a , bool           ) :
                
                _args.append   (  ROOT.RooFit.CloneData    ( a )  )
                
            elif key in ( 'offset' ,         ) and isinstance ( a , bool           ) :
                
                _args.append   (  ROOT.RooFit.Offset       ( a )  )
                
            elif key in ( 'fitoption'  ,
                          'fitoptions' ,
                          'optionfit'  ,
                          'potionsfit' ) and \
                          isinstance ( a , string_types   ) and root_info < ( 6 , 28 ) :
                
                _args.append   (  ROOT.RooFit.FitOptions   ( a )  )
                
            elif key in ( 'constraint'       ,
                          'constraints'      ,
                          'pars'             ,
                          'params'           ,
                          'parameters'       ) :
                c = self.parse_constraints ( a )
                if c is None : self.error ('parse_args: Invalid constraint specification: %s/%s' % ( a , type ( a ) ) )
                else         : _args.append ( c ) 

            elif key in ( 'integratebin'  ,
                          'integratebins' ,                          
                          'binintegrate'  ,
                          'binsintegrate' ,
                          'integrate'     )  and \
                          isinstance ( a , num_types )  and (6,24) <= root_info :
                
                _args.append   (  ROOT.RooFit.IntegrateBins ( a ) )
                
            elif key in ( 'newstyle' ,
                          'stylenew' ,
                          'new'      ) and \
                          isinstance ( a , bool   ) and (6,27) <= root_info :
                
                _args.append   (  ROOT.RooFit.NewStyle ( a ) )

            elif key in ( 'parallelize' ,
                          'parallelise' ,
                          'parallel' ) and \
                          isinstance ( a , sized_types )       and \
                          1<= len ( a ) <= 3                   and (6,27) <= root_info :
                
                _args.append   (  ROOT.RooFit.Parallelize ( *a ) )
                
            elif key in ( 'recover'                     ,
                          'recoverfromundefined'        ,
                          'recoverfromundefinedregions' ) and \
                          isinstance ( a , num_types )    and ( 6, 24 ) <= root_info :
                
                _args.append   ( ROOT.RooFit.RecoverFromUndefinedRegions ( 1.0 * float ( a ) ) ) 
                          
            elif key in ( 'maxcall'  ,
                          'maxcalls' ,
                          'callmax'  ,
                          'callsmax' ) and \
                     isinstance ( a , integer_types )    and (6,27) <= root_info :
                
                _args.append   ( ROOT.RooFit.MaxCalls ( a ) ) 

            else :
               
                self.error ( 'parse_args: Unknown/illegal keyword argument: %s/%s, skip it ' % ( k , type ( a ) ) )
            
        
        if not check_arg ( 'numcpu' , *_args ) :
            if  dataset and not isinstance ( dataset , ROOT.RooDataHist ) :
                _args.append ( ncpu ( len ( dataset ) ) )
            else :
                nc = numcpu()
                if  1 < nc : _args.append ( ROOT.RooFit.NumCPU ( nc ) )

                
        # =============================================================
        ## check options for the weighted datasets 
        if dataset and not isinstance ( dataset , ROOT.RooDataHist ) :
            
            weighted = dataset.isWeighted ()            
            sw2      = check_arg  ( 'SumW2Error'      , *_args  )
            aer      = check_arg  ( 'AsymptoticError' , *_args  )
            binned   = isinstance ( dataset  , ROOT.RooDataHist )
            
            if ( 6 , 25 )<= root_info and binned and ( not sw2 ) and ( not aer ) :
                _args.append ( ROOT.RooFit.SumW2Error ( True ) )                
            elif sw2 and aer :
                self.warning ( "parse_args: Both 'SumW2Error' and 'AsymptoticError' are specified" )                
            elif weighted   and sw2 :
                value = bool ( sw2.getInt( 0 ) )
                if not value : self.warning ("parse_args: 'SumW2=False' is specified for the weighted  dataset!")
            elif weighted and aer : 
                value = bool ( aer.getInt( 0 ) )
                if not value : self.warning ("parse_args: 'AsymptoticError=False' is specified for the weighted  dataset!")
            elif weighted :                
                if binned : self.debug   ( "parse_args: Neither 'SumW2Error' nor 'AsymptoticError' are specified for weighted/binned dataset! 'SumW2=True' is added" )
                else      : self.warning ( "parse_args: Neither 'SumW2Error' nor 'AsymptoticError' are specified for weighted dataset! 'SumW2=True' is added" )
                _args.append ( ROOT.RooFit.SumW2Error ( True ) )                
            elif not weighted and sw2 :
                self.warning ( "parse_args:'SumW2Error' is specified for non-weighted dataset" )
            elif not weighted and aer :
                self.warning ( "parse_args:'AsymptoticError' is specified for non-weighted dataset" )

        keys = [ str ( a ) for a in _args ]
        keys.sort ()
        
        ## check presence of "non-trivial"  keys
        kset = set( keys ) 
        kset.discard  ( 'Save'       ) ## trivial
        kset.discard  ( 'NumCPU'     ) ## trivial
        kset.discard  ( 'Verbose'    ) ## trivial 
        kset.discard  ( 'Timer'      ) ## trivial 
        kset.discard  ( 'PrintLevel' ) ## trivial

        ## duplicates? 
        if len ( kset ) != len ( keys ) :
            self.warning ("duplicated options!")            
        #
        if kset : self.debug ( 'parse_args: Parsed arguments %s' % keys )
        else    : self.debug ( 'parse_args: Parsed arguments %s' % keys )


        ## store them 
        self.aux_keep.append ( _args ) 
        
        return tuple ( _args ) 

    # =========================================================================
    ## technical method to parse the constraint argument
    #  @code
    #  pdf.fiTo ( ..  , constraints = ... , ... )
    #  @endcode
    def parse_constraints ( self , arg ) :
        """Technical method to parse the constraints argument
        >>>  pdf.fiTo ( ..  , constraints = ... , ... )
        """
        
        if   isinstance ( arg , ROOT.RooCmdArg   ) :
            return arg        
        elif isinstance ( arg , ROOT.RooArgSet   ) :
            return ROOT.RooFit.ExternalConstraints  ( arg )            
        elif isinstance ( arg , ROOT.RooAbsReal  ) :
            return ROOT.RooFit.ExternalConstraints  ( ROOT.RooArgSet ( arg ) )
        elif isinstance ( arg , ( tuple , list ) ) and 2 == len ( arg ) \
                 and isinstance ( arg[0] , ROOT.RooAbsReal ) and isinstance ( arg[1] , VE ) :
            return self.make_constraint ( arg[0] , arg[1] ) 
        
        elif isinstance ( arg , ( tuple , list ) ) :
            c = ROOT.RooArgSet ( )
            for ia in arg :
                if   isinstance ( ia , ROOT.RooAbsReal ) : c.add ( ia )
                elif isinstance ( ia , ( tuple , list  ) ) and 2 == len ( ia ) \
                         and isinstance ( ia[0] , ROOT.RooAbsReal ) and isinstance ( ia[1] , VE ) :
                    c.add( self.soft_constraint ( ia[0] , ia[1] ) )
                else  : self.error ('parse_constraint: invalid constraint: %s/%s from %s'  % ( ia , type(ia) , arg ) )  
            self.aux_keep.append ( c )
            return ROOT.RooFit.ExternalConstraints  ( c )

        self.error ('parse_constraint: Invalid specification: %s/%s' % ( arg , type ( arg) ) )
        return None 

    # =========================================================================
    ## set value to a given value with the optional check
    #  @code
    #  pdf = ...
    #  pdf.set_value ( my_var1 , 10 )
    #  pdf.set_value ( my_var2 , 10 , lambda a,b : b>0 ) 
    #  @endcode 
    @staticmethod 
    def set_value ( var , value , ok = lambda a , b : True ) :
        """set value to a given value with the optional check
        pdf = ...
        pdf.set_value ( my_var1 , 10 )
        pdf.set_value ( my_var2 , 10 , lambda a,b : b>0 ) 
        """

        ## must be roofit variable! 
        assert isinstance ( var , ROOT.RooAbsRealLValue ) , \
               "set_value: invalid type of 'var' %s" % type ( var )
        
        if not hasattr ( var ,  'setVal' ) :
            raise ValueError ( "No value can be set for %s/%s" % ( var , type ( var ) ) )  

        ## convert to float 
        value = float ( value )

        ## check for the range, if range is defined 
        minmax = var.minmax ()
        if minmax :
            mn , mx = minmax
            if not ( mn <= value <= mx or isequal ( mn , value ) or isequal ( mx , value ) ) :
                raise ValueError ( "set_value: value %s is outside of the [%s,%s] region" % ( value , mn , mx ) ) 
            
        ## check for external conditions, if specified  
        if not ok ( var , value ) :
            raise ValueError ( "Value %s is not OK" % value ) 

        ## finally set the value 
        var.setVal ( value )

        return isequal ( value , var.getVal () ) 

    # =========================================================================
    ## gettter for certain fit components from the provided list 
    def component_getter ( self , components ) :
        """Gettter for certain fit components from the provided list
        """
        nc = len ( components )
        if   0 == nc : return ()
        elif 1 == nc : return components [ 0 ] 
        return tuple ( c for c in components )

    # ======================================================
    ## setter for certian fit components form provided list 
    def component_setter ( self , components , value ) :
        """Setter for certian fit components form provided list
        """
        assert 0 < len ( components ) , 'Empty list of components, setting is not possible!'
        
        if   isinstance ( value , num_types          ) : values = float ( value ) , 
        elif isinstance ( value , VE                 ) : values = value.value ()  , 
        elif isinstance ( value , ROOT.RooAbsReal    ) : values = float ( value ) , 
        elif isinstance ( value , ROOT.RooArgList    ) : values = tuple ( float ( v ) for v in value ) 
        elif isinstance ( value , sequence_types     ) : values = tuple ( float ( v ) for v in value ) 
        else :
            self.warning ("component setter: unknown type for 'value':%s/%s" % ( str( value) , type ( value ) ) )
            values = value 

        for c , v in zip  ( components , values )  :
            self.set_value (  c , v ) 

    # =========================================================================
    ## Create some expressions with variables
    # =========================================================================
    
    # =========================================================================
    ## construct (on-flight) variable for the product of
    #  <code>var1</code> and <code>var2</code> \f$ v \equiv  v_1 v_2\f$ 
    #  @code
    #  var1 = ...
    #  var2 = ...
    #  var3 = xxx.vars_multiply ( var1 , var2 )
    #  var4 = xxx.vars_multiply ( var1 , 2.0  )    
    #  var3 = xxx.vars_product  ( var2 )
    #  var4 = xxx.vars_product  ( var1 , 2.0  )    
    #  @endcode 
    def vars_multiply ( self , var1 , var2 , name = '' , title = '' ) :
        """Construct (on-flight) variable for var1*var2 
        >>> var1 = ...
        >>> var2 = ...
        >>> var3 = xxx.vars_multiply ( var1 , var2   )
        >>> var4 = xxx.vars_multiply ( var1 , 2 , 'sigma2' , title = 'Scaled sigma' )
        >>> var3 = xxx.vars_product  ( var1 , var2   )
        >>> var4 = xxx.vars_product  ( var1 , 2 , 'sigma2' , title = 'Scaled sigma' )
        """
        
        f1 = isinstance ( var1 , num_types )
        f2 = isinstance ( var2 , num_types )

        if f1 and f2 :
            res  = float ( var1 ) * float ( var2 )
            return ROOT.RooFit.RooConst ( res ) 
        elif f1 :
            # shortcut 
            if   1 == var1 : return var2                       ## SHORTCUT
            elif 0 == var1 : return ROOT.RooFit.RooConst ( 0 ) ## SHORTCUT
            # 
            var1 = ROOT.RooFit.RooConst ( var1 ) 
            return self.vars_multiply ( var1 , var2 , name , title )
        elif f2 : 
            # shortcut 
            if   1 == var2 : return var1                       ## SHORTCUT
            elif 0 == var2 : return ROOT.RooFit.RooConst ( 0 ) ## SHORTCUT
            # 
            var2 = ROOT.RooFit.RooConst ( var2 ) 
            return self.vars_multiply ( var1 , var2 , name , title )
        
        self.aux_keep.append ( var1 )
        self.aux_keep.append ( var2 )

        name  = name  if name   else self.roo_name  ( "mult_%s_%s" % ( var1.name , var2.name ) )
        title = title if title  else                  "(%s)*(%s)"  % ( var1.name , var2.name )
        
        result = Ostap.MoreRooFit. Product  ( var1 , var2 , name , title )
        self.aux_keep.append ( result  )
        
        return result

    # =============================================================================
    ## construct (on-flight) variable for the sum  of
    #  <code>var1</code> and <code>var2</code> \f$ v\equiv  c_1v_1 + c_2v_2\f$ 
    #  @code
    #  var1 = ...
    #  var2 = ...
    #  var3 = xxx.vars_add ( var1 , var2 )
    #  var4 = xxx.vars_add ( var1 , 2.0  )    
    #  var5 = xxx.vars_sum ( var1 , var2 )
    #  var6 = xxx.vars_sum ( var1 , 2.0  )    
    #  @endcode 
    def vars_add ( self , var1 , var2 , c1 = 1 , c2 = 1 , name = '' , title = '' ) :
        """Construct (on-flight) variable for var1*c1+var2*c2 
        >>> var1 = ...
        >>> var2 = ...
        >>> var3 = xxx.vars_add ( var1 , var2   )
        >>> var4 = xxx.vars_add ( var1 , 2 , 'sigma2' , title = 'Scaled sigma' )
        >>> var5 = xxx.vars_sum ( var1 , var2   )
        >>> var6 = xxx.vars_sum ( var1 , 2 , 'sigma2' , title = 'Scaled sigma' )
        """
        
        assert isinstance ( c1 , num_types ) and isinstance ( c2 , num_types ),\
               "vars_add: c1 and c2 must be numeric types!"

        c1 = float ( c1 )
        c2 = float ( c2 )
        
        if   0 == c1 and 0 == c2 :
            return ROOT.RooFit.RooConst ( 0 )   ## shortcut 
        elif 0 == c1 : 
            return self.vars_multiply ( var2 , c2 , name = name , title = title )
        elif 0 == c2 :
            return self.vars_multiply ( var1 , c1 , name = name , title = title )
        
        f1 = isinstance ( var1 , num_types )
        f2 = isinstance ( var2 , num_types )

        if f1 and f2 :
            res  = float ( var1 ) * float ( c1 ) + float ( var2 ) * float ( c2 )
            return ROOT.RooFit.RooConst ( res ) 
        elif f1 :
            ## shortcut 
            if 0 == var1 :
                return self.var_multiply ( var2 , c2 , name = name , title = title )  ## SHORTCUT 
            #
            var1 = ROOT.RooFit.RooConst ( float ( var1 ) * float ( c1 ) )                         
            return self.vars_add ( var1 , var2 , name = name , title = title )
        elif f2 :
            ## shortcut 
            if 0 == var2 :
                return self.var_multiply ( var1 , c1 , name = name , title = title )  ## SHORTCUT 
            #
            var2 = ROOT.RooFit.RooConst ( float ( var2 ) * float ( c2 ) ) 
            return self.vars_add ( var1 , var2 , name =name , title = title  )
        
        self.aux_keep.append ( var1 )
        self.aux_keep.append ( var2 )

        name  = name  if name   else self.roo_name  ( "add_%s_%s" % ( var1.name , var2.name ) )
        title = title if title  else                  "(%s)+(%s)" % ( var1.name , var2.name )

        if c1 == 1 and c2 == 1 : 
            result = Ostap.MoreRooFit.Addition  ( var1 , var2 , name , title )
        else :
            result = Ostap.MoreRooFit.Addition2 ( name , title , var1 , var2 , c1 , c2 )
                
        self.aux_keep.append ( result  )
                
        return result

    # =============================================================================
    ## construct (on-flight) variable for the subtraction  of
    #  <code>var1</code> and <code>var1</code>: \f$ v\equiv v_1 - v_2\f$ 
    #  @code
    #  var1 = ...
    #  var2 = ...
    #  var3 = xxx.vars_subtract   ( var1 , var2 )
    #  var4 = xxx.vars_subtract   ( var1 , 2.0  )    
    #  var5 = xxx.vars_difference ( var1 , var2 )
    #  var6 = xxx.vars_difference ( var1 , 2.0  )    
    #  @endcode 
    def vars_subtract ( self , var1 , var2 , name = '' , title = '' ) :
        """Construct (on-flight) variable  for var1-var2 
        >>> var1 = ...
        >>> var2 = ...
        >>> var3 = xxx.vars_subtract   ( var1 , var2   )
        >>> var4 = xxx.vars_subtract   ( var1 , 2 , 'sigma2' , title = 'Scaled sigma' )
        >>> var5 = xxx.vars_difference ( var1 , var2   )
        >>> var6 = xxx.vars_difference ( var1 , 2 , 'sigma2' , title = 'Scaled sigma' )
        """

        f1 = isinstance ( var1 , num_types )
        f2 = isinstance ( var2 , num_types )

        if f1 and f2 :
            ##
            res  = float ( var1 ) - float ( var2 )
            return ROOT.RooFit.RooConst ( res ) 
        elif f1 :
            ## 
            var1 = ROOT.RooFit.RooConst ( var1 )                         
            return self.vars_subtract ( var1 , var2 , name , title )
        elif f2 :
            ## shortcut 
            if 0 == var2 : return var1                      ## SHORTCUT
            #
            var2 = ROOT.RooFit.RooConst ( var2 ) 
            return self.vars_subtract ( var1 , var2 , name , title )

        self.aux_keep.append ( var1 )
        self.aux_keep.append ( var2 )

        name  = name  if name   else self.roo_name  ( "subtract_%s_%s" % ( var1.name , var2.name ) )
        title = title if title  else                  "(%s)-(%s)"      % ( var1.name , var2.name )

        result = Ostap.MoreRooFit.Subtraction  ( var1 , var2 )
        self.aux_keep.append ( result  )
                
        return result
            
    # =============================================================================
    ## construct (on-flight) variable for the ratio of
    #  <code>var1</code> and <code>var2</code>: \f$   v\equiv \frac{v_1}{v_2} \f$ 
    #  @code
    #  var1 = ...
    #  var2 = ...
    #  var3 = xxx.vars_divide ( var1 , var2 )
    #  var4 = xxx.vars_divide ( var1 , 2.0  )    
    #  var5 = xxx.vars_ratio  ( var1 , var2 )
    #  var6 = xxx.vars_ratio  ( var1 , 2.0  )    
    #  @endcode 
    def vars_divide ( self , var1 , var2 , name = '' , title = '' ) :
        """Construct (on-flight) variable for var1/var2 
        >>> var1 = ...
        >>> var2 = ...
        >>> var3 = xxx.vars_divide ( var1 , var2   )
        >>> var4 = xxx.vars_divide ( var1 , 2 , 'sigma2' , title = 'Scaled sigma' )
        >>> var5 = xxx.vars_ratio  ( var1 , var2   )
        >>> var6 = xxx.vars_ratio  ( var1 , 2 , 'sigma2' , title = 'Scaled sigma' )
        """
        
        f1 = isinstance ( var1 , num_types )
        f2 = isinstance ( var2 , num_types )

        if f1 and f2 :
            res   = float ( var1 ) / float ( var2 )
            return ROOT.RooFit.RooConst ( res ) 
        elif f1 :
            var1 = ROOT.RooFit.RooConst ( var1 ) 
            return self.vars_divide   ( var1 , var2     , name , title )
        elif f2 :
            return self.vars_multiply ( var1 , 1.0/var2 , name , title )
        
        self.aux_keep.append ( var1 )
        self.aux_keep.append ( var2 )
        
        name  = name  if name   else self.roo_name  ( "divide_%s_%s" % ( var1.name , var2.name ) )
        title = title if title  else                  "(%s)/(%s)"    % ( var1.name , var2.name )
        
        result = Ostap.MoreRooFit.Division  ( var1 , var2 , name , title )
        self.aux_keep.append ( result  )
                
        return result

    # =============================================================================
    ## construct (on-flight) variable  for the fraction  of
    #  <code>var1</code> and <code>var2</code>: \f$ v \equiv \frac{v_1}{v_1+v_2}\f$
    #  @code
    #  var1 = ...
    #  var2 = ...
    #  var3 = xxx.vars_fraction ( var1 , var2 )
    #  var4 = xxx.vars_fraction ( var1 , 2.0  )    
    #  @endcode 
    def vars_fraction ( self , var1 , var2 , name = '' , title = '' ) :
        """Construct (on-flight) variable  for var1/(var2+var1)
        >>> var1 = ...
        >>> var2 = ...
        >>> var3 = xxx.vars_fraction ( var1 , var2   )
        >>> var4 = xxx.vars_fraction ( var1 , 2 , 'sigma2' , title = 'exression' )
        """
        
        f1 = isinstance ( var1 , num_types )
        f2 = isinstance ( var2 , num_types )

        if f1 and f2 :
            res  = float ( var1 ) / ( float ( var2 ) + float ( var1 ) )
            return ROOT.RooFit.RooConst ( res ) 
        elif f1 :
            ## shortcut 
            if 0 == var1  : return ROOT.RooFit.RooConst ( 0 ) ## SHORTCUT
            var1 = ROOT.RooFit.RooConst ( var1 ) 
            return self.vars_fraction ( var1 , var2 , name , title )
        elif f2 :
            ## shortcut
            if 0 == var2  : return ROOT.RooFit.RooConst ( 1 ) ## SHORTCUT
            var2 = ROOT.RooFit.RooConst ( var2 ) 
            return self.vars_fraction ( var1 , var2 , name , title )
        
        self.aux_keep.append ( var1 )
        self.aux_keep.append ( var2 )

        name  = name  if name   else self.roo_name  ( "fraction_%s_%s"  % ( var1.name , var2.name ) )
        title = title if title  else                  "fraction(%s,%s)" % ( var1.name , var2.name )

        result = Ostap.MoreRooFit.Fraction  ( var1 , var2 , name , title )
        self.aux_keep.append ( result  )
        
        return result

    # =============================================================================
    ## construct (on-flight) variable  for the asymmetry of
    #  <code>var1</code> and <code>var2</code>: \f$ v \equiv \frac{v_1-v_2}{v_1+v_2}\f$
    #  @code
    #  var1 = ...
    #  var2 = ...
    #  var3 = xxx.vars_asymmetry     ( var1 , var2 )
    #  var4 = xxx.vars_asymmetry     ( var1 , 2.0  )    
    #  var3 = xxx.vars_reldifference ( var1 , var2 )
    #  var4 = xxx.vars_reldifference ( var1 , 2.0  )    
    #  @endcode 
    def vars_asymmetry ( self , var1 , var2 , scale = 1 , name = '' , title = '' ) :
        """Construct (on-flight) variable for (var1-var2)/(var2+var1)
        >>> var1 = ...
        >>> var2 = ...
        >>> var3 = xxx.vars_asymmetry     ( var1 , var2   )
        >>> var4 = xxx.vars_asymmetry     ( var1 , 2 , 'sigma2' , title = 'exression' )
        >>> var3 = xxx.vars_reldifference ( var1 , var2   )
        >>> var4 = xxx.vars_reldifference ( var1 , 2 , 'sigma2' , title = 'exression' )
        """
        
        f1 = isinstance ( var1 , num_types )
        f2 = isinstance ( var2 , num_types )

        if f1 and f2 :
            res  = scale * ( float ( var1 ) - float ( var2 ) ) / ( float ( var2 ) + float ( var1 ) )
            return ROOT.RooFit.RooConst ( res ) 
        elif f1 :
            ## shortcut 
            if 0 == var1 : return ROOT.RooFit.RooConst ( -1.0 * scale ) ## shortcut
            var1 = ROOT.RooFit.RooConst ( var1 ) 
            return self.vars_asymmetry ( var1 , var2 , scale = scale , name = name , title = title )
        elif f2 :
            ## shortcut 
            if 0 == var2 : return ROOT.RooFit.RooConst (  1.0 * scale ) ## shortcut
            var2 = ROOT.RooFit.RooConst ( var2 ) 
            return self.vars_asymmetry ( var1 , var2 , scale = scale , name = name , title = title )
        
        self.aux_keep.append ( var1 )
        self.aux_keep.append ( var2 )

        name  = name  if name   else self.roo_name  ( "asymmetry_%s_%s"  % ( var1.name , var2.name ) )
        title = title if title  else                  "asymmetry(%s,%s)" % ( var1.name , var2.name )

        result = Ostap.MoreRooFit.Asymmetry ( var1 , var2 , name , title , scale )
        self.aux_keep.append ( result  )
        
        return result

    # =============================================================================
    ## construct (on-flight) variable  for \f$ a^b \f$ 
    #  <code>var1</code> and <code>var2</code>: 
    #  @code
    #  var1 = ...
    #  var2 = ...
    #  var3 = xxx.vars_power        ( var1 , var2 )
    #  var4 = xxx.vars_power        ( var1 , 2.0  )    
    #  var4 = xxx.vars_power        ( 2.0  , var2 )    
    #  @endcode 
    def vars_power ( self , var1 , var2 , name = '' , title = '' ) :
        """ construct (on-flight) variable  for \f$ a^b \f$ 
        >>> var1 = ...
        >>> var2 = ...
        >>> var3 = xxx.vars_power        ( var1 , var2 )
        >>> var4 = xxx.vars_power        ( var1 , 2.0  )    
        >>> var4 = xxx.vars_power        ( 2.0  , var2 )    
        """
        f1 = isinstance ( var1 , num_types )
        f2 = isinstance ( var2 , num_types )

        if f1 and f2 :
            res  = float ( var1 ) ** float ( var2 )
            return ROOT.RooFit.RooConst ( res ) 
        elif f1 :
            ## shortcut 
            if 1 == var1 : return ROOT.RooFit.RooConst ( 1 ) ## shortcut
            var1 = ROOT.RooFit.RooConst ( var1 ) 
            return self.vars_power ( var1 , var2 , name , title )
        elif f2 :
            ## shortcut 
            if 0 == var2 : return ROOT.RooFit.RooConst (  1 ) ## shortcut
            var2 = ROOT.RooFit.RooConst ( var2 ) 
            return self.vars_power ( var1 , var2 , name , title )
        
        self.aux_keep.append ( var1 )
        self.aux_keep.append ( var2 )

        name  = name  if name   else self.roo_name  ( "pow_%s_%s"  % ( var1.name , var2.name ) )
        title = title if title  else                  "pow(%s,%s)" % ( var1.name , var2.name )

        result = Ostap.MoreRooFit.Power ( var1 , var2 , name . title )
        self.aux_keep.append ( result  )
        
        return result

    # =============================================================================
    ## construct (on-flight) variable  for \f$ exp(ab) \f$ 
    #  <code>var1</code> and <code>var2</code>: 
    #  @code
    #  var1 = ...
    #  var2 = ...
    #  var3 = xxx.vars_exp ( var1 , var2 )
    #  var4 = xxx.vars_exp ( var1 , -1   )    
    #  var4 = xxx.vars_exp ( -1   , var2 )    
    #  @endcode 
    def vars_exp ( self , var1 , var2 = 1 , name = '' , title = '' ) :
        """ construct (on-flight) variable  for \f$ exp(a*b) \f$ 
        >>> var1 = ...
        >>> var2 = ...
        >>> var3 = xxx.vars_exp ( var1 , var2 )
        >>> var4 = xxx.vars_exp ( var1 , 2.0  )    
        >>> var4 = xxx.vars_exp ( 2.0  , var2 )    
        """
        f1 = isinstance ( var1 , num_types )
        f2 = isinstance ( var2 , num_types )

        if f1 and f2 :
            res  = math.exp ( float ( var1 ) * float ( var2 ) )
            return ROOT.RooFit.RooConst ( res ) 
        elif f1 :
            ## shortcut 
            if 0 == var1 : return ROOT.RooFit.RooConst ( 1 ) ## shortcut
            var1 = ROOT.RooFit.RooConst ( var1 ) 
            return self.vars_exp ( var1 , var2 , name , title )
        elif f2 :
            ## shortcut 
            if 0 == var2 : return ROOT.RooFit.RooConst ( 1 ) ## shortcut
            var2 = ROOT.RooFit.RooConst ( var2 ) 
            return self.vars_exp ( var1 , var2 , name , title )
        
        self.aux_keep.append ( var1 )
        self.aux_keep.append ( var2 )

        name  = name  if name   else self.roo_name  ( "exp_%s_%s"  % ( var1.name , var2.name ) )
        title = title if title  else                  "exp(%s*%s)" % ( var1.name , var2.name )

        result = Ostap.MoreRooFit.Exp ( var1 , var2 , name , title )
        self.aux_keep.append ( result  )
        
        return result

    # =============================================================================
    ## construct (on-flight) variable  for \f$ |ab| \f$ 
    #  <code>var1</code> and <code>var2</code>: 
    #  @code
    #  var1 = ...
    #  var2 = ...
    #  var3 = xxx.vars_abs ( var1 , var2 )
    #  var4 = xxx.vars_abs ( var1 , -1   )    
    #  var4 = xxx.vars_abs ( -1   , var2 )    
    #  @endcode 
    def vars_abs ( self , var1 , var2 = 1 , name = '' , title = '' ) :
        """ construct (on-flight) variable  for \f$ |a*b| \f$ 
        >>> var1 = ...
        >>> var2 = ...
        >>> var3 = xxx.vars_abs ( var1 , var2 )
        >>> var4 = xxx.vars_abs ( var1 , 2.0  )    
        >>> var4 = xxx.vars_abs ( 2.0  , var2 )    
        """
        f1 = isinstance ( var1 , num_types )
        f2 = isinstance ( var2 , num_types )

        if f1 and f2 :
            res  = abs ( float ( var1 ) * float ( var2 ) )
            return ROOT.RooFit.RooConst ( res ) 
        elif f1 :
            ## shortcut 
            if 0 == var1 : return ROOT.RooFit.RooConst ( 0 ) ## shortcut
            var1 = ROOT.RooFit.RooConst ( var1 ) 
            return self.vars_abs ( var1 , var2 , name , title )
        elif f2 :
            ## shortcut 
            if 0 == var2 : return ROOT.RooFit.RooConst ( 0 ) ## shortcut
            var2 = ROOT.RooFit.RooConst ( var2 ) 
            return self.vars_abs ( var1 , var2 , name , title )
        
        self.aux_keep.append ( var1 )
        self.aux_keep.append ( var2 )

        name  = name  if name   else self.roo_name  ( "abs_%s_%s"  % ( var1.name , var2.name ) )
        title = title if title  else                  "abs(%s*%s)" % ( var1.name , var2.name )

        result = Ostap.MoreRooFit.Abs ( var1 , var2 , name , title )
        self.aux_keep.append ( result  )
        
        return result

    vars_sum           = vars_add 
    vars_product       = vars_multiply
    vars_ratio         = vars_divide
    vars_difference    = vars_subtract
    vars_reldifference = vars_asymmetry
    vars_pow           = vars_power
    vars_expo          = vars_exp

    # =========================================================================
    ## helper function to create <code>RooFormulaVar</code>
    #  @code
    #  OBJ = ...
    #  f   = OBJ.vars_formula ( '%s*%s/%s' , [ a , b, c ] , name = 'myvar' )
    #  @endcode 
    def vars_formula ( self , formula , vars , name = '' , title = '' ) :
        """helper function to create <code>RooFormulaVar</code>
        >>> OBJ = ...
        >>> f   = OBJ.vars_formula ( '%s*%s/%s' , [ a , b, c ] , name = 'myvar' )
        """

        assert vars and len ( vars ) , 'Variables must be specified!'

        vvars = []
        for v in  vars :
            if isinstance    ( v , ROOT.RooAbsArg ) :
                vvars.append ( v )
            elif isinstance  ( v , string_types   ) :
                try :
                    vv = self.parameter ( v )
                    vvars.append ( vv ) 
                except :
                    raise TypeError ( "Unknown parameter name %s" % v)
            else :
                raise TypeError( "Unknown parameter type %s/%s" % ( v , type ( v ) ) ) 

        vlst = ROOT.RooArgList()
        for v in vvars : vlst.add ( v )

        has_at      = '@' in formula
        has_percent = '%' in formula
        import re
        has_index   = re.search ( r'\[( *)(?P<degree>\d*)( *)\]' , formula )
        has_format1 = re.search ( r'\{( *)(?P<degree>\d*)( *)\}' , formula )
        has_format2 = re.search ( r'\{( *)(?P<degree>\w*)( *)\}' , formula )

        formula_ = formula 
        if   has_at      : pass 
        elif has_index   : pass 
        elif has_percent :  
            vnames   = tuple ( [ p.name for p in vlst ] )
            formula_ = formula % vnames
        elif has_format1 :            
            vnames   = tuple ( [ p.name for p in vlst ] )
            formula_ = formula.format ( *vnames ) 
        elif has_format2 :
            kw  = {}
            for p in vlist : kw [ p.name ] = p.name
            formula_ = formula.format ( *kw )
            
        name  = name  if name  else 'Formula_%s '    % self.name 
        title = title if title else 'Formula:%s/%s'  % ( formula , self.name )
        
        rfv = make_formula ( self.var_name ( name ) , formula_ , vlst )
        
        self.aux_keep.append ( vlst )
        self.aux_keep.append ( rvf  )
        
        return rfv 

    # =========================================================================
    ## make very specific combination of variables:  alpha*var1*(bets+gamma*var2)    
    #  \f$ r = \alpha v_1 ( \beta + \gamma * v_2 ) \    
    def vars_combination ( self ,
                           var1 ,
                           var2 ,
                           alpha  = 1   ,
                           beta   = 1   ,
                           gamma  = 1   ,
                           name   = ''  , 
                           title  = ''  ) :
        """Make very specific combination of variables:  alpha*var1*(bets+gamma*var2)    
        r = alpha * v_1 ( beta + gamma * v_2 ) 
        """
        
        f1 = isinstance ( var1 , num_types )
        f2 = isinstance ( var2 , num_types )

        assert isinstance ( alpha , num_types ) , "vars_combination: 'alpha' must be numeric types!"
        assert isinstance ( beta  , num_types ) , "vars_combination: 'alpha' must be numeric types!"
        assert isinstance ( gamma , num_types ) , "vars_combination: 'alpha' must be numeric types!"

        if   0 == alpha               : return ROOT.RooFit.RooConst ( 0 ) 
        elif 0 == beta and 0 == gamma : return ROOT.RooFit.RooConst ( 0 ) 

        alpha = float ( alpha )
        beta  = float ( beta  )
        gamma = float ( gamma )        
        
        if f1 and f2 :

            res  = beta   + gamma * float ( var2 )
            res *=          alpha * float ( var1 )
            
            return ROOT.RooFit.RooConst ( 0 )
        
        elif f1 :
            
            return self.sum_add      ( float ( var1 ) * alpha * beta        ,
                                       var2                                 ,
                                       c2    = alpha * gamma * float ( v1 ) , 
                                       name  = name                         ,
                                       title = title                        ) 
        elif f2 :
            
            return self.sum_multiply ( var1 ,
                                       alpha * ( beta + gamma * float ( v2 ) ) ,
                                       name  = name  ,
                                       title = title ) 
        
        
        self.aux_keep.append ( var1 )
        self.aux_keep.append ( var2 )

        if not name :
            if   1 == alpha and 1 == beta and +1 == gamma :
                name = '%s*(1+%s)' % ( var1.name , var2.name )
            elif 1 == alpha and 1 == beta and -1 == gamma :
                name = '%s*(1-%s)' % ( var1.name , var2.name )
            else :
                name = 'Var_%g*%s*(%s%+g*%s)' % ( alpha , var1.name , beta , gamma , var2.name )
                
        if not title : title = name
                    
        result = Ostap.MoreRooFit.Combination ( var1  , var2  ,
                                                name  , title ,
                                                alpha , beta  , gamma )
        
        self.aux_keep.append ( result  )
        
        return result
        
                           
    # =========================================================================
    ## (var1,var2) <--> (var, asymmetry)
    # =========================================================================

    # =========================================================================
    ## Convert pair of variables into 'sum' & 'asymmetry' pair:
    # - 'sum'       :  sum_scale  * ( var1 + var2  )
    # - 'asymmetry' :  asym_scale * ( var1 - var2 ) / (var1+ + var2 ) 
    # @code
    # var1 = ...
    # var2 = ...
    # hsum , asum = self.vars_to_asymmetry ( var1 , var2 ) 
    # @endcode 
    def vars_to_asymmetry ( self             ,
                            var1             ,   ## the first variable 
                            var2             ,   ## the second variable
                            asym_scale = 1   ,   ## scale factor asymmetry
                            sum_scale  = 0.5 ,   ## scale factor for 'sum'
                            asym_name  = ''  ,   ## name for asymmetry variable
                            sum_name   = ''  ,   ## name for 'sum' variable                              
                            asym_title = ''  ,   ## title for asymmetry variable 
                            sum_title  = ''  ) : ## title for 'sum' variable 
        """Convert pair of variables into 'sum' & 'asymmetry' pair
        - 'sum'       :  sum_scale  * ( var1 + var2  )
        - 'asymmetry' :  asym_scale * ( var1 - var2 ) / (var1+ + var2 ) 
        >>> var1 = ...
        >>> var2 = ...
        >>> hsum , asum = self.vars_to_asymmetry ( var1 , var2 ) 
        """
        
        ## asym_scale * (var1-var2)/(var1+var2)
        asym_var = self.vars_asymmetry (
            var1               , ## first vatiable 
            var2               , ## second variable 
            scale = scale      , ## scale
            name  = asym_name  ,
            title = asym_title ) 

        ## sum_scale * ( var1 + var2 ) 
        sum_var  = self.vars_add (
            var1               , ## first variable 
            var2               , ## second variable
            c1    = sum_scale  , ## factor c1 
            c2    = sum_scale  , ## factor c2 
            name  = sum_name   , ## name 
            title = sum_title  ) ## title  
        
        return sum_var , asym_var
    
    # ========================================================================= 
    ## convert a pair of variables 'half-sum'&'asymmetry' into 'var1', 'var2'
    #  @code
    #  halfsum   = ...
    #  asymmetry = ...
    #  var1 , var2 = self.vars_from_asymmetry ( halfsum , asymmetry )
    #  @endcode 
    def vars_from_asymmetry ( self         ,
                              hsumvar      , ## half-sum 
                              asymvar      , ## asymmetry 
                              v1name  = '' , 
                              v2name  = '' ,
                              v1title = '' ,
                              v2title = '' ) :
        """Convert a pair of variables 'half-sum'&'asymmetry' into 'var1', 'var2'
        >>> halfsum   = ...
        >>> asymmetry = ...
        >>> var1 , var2 = self.vars_from_asymmetry ( halfsum , asymmetry )
        """

        if hsumvar is None :
            return hsumvar , hsumvar
        
        if isinstance ( hsumvar , ROOT.RooAbsArg ) and isinstance ( asymvar , ROOT.RooAbsArg )  :
            s = hsumvar.name
            a = asymvar.name
            if not v1name  : v1name  = '%sL' % s
            if not v2name  : v2name  = '%sR' % s            
            if not v1title : v1title = '%s_{L} : %s #times (1 + %s_{%s}) ' %  ( s , s , a , s ) 
            if not v2title : v2title = '%s_{R} : %s #times (1 - %s_{%s}) ' %  ( s , s , a , s ) 
                 
        var1 = self.vars_combination ( hsumvar           ,
                                       asymvar           ,
                                       alpha   =  1      ,
                                       beta    =  1      ,
                                       gamma   = +1      , 
                                       name    = v1name  ,
                                       title   = v1title )
        
        var2 = self.vars_combination ( hsumvar           ,
                                       asymvar           ,
                                       alpha   =  1      ,
                                       beta    =  1      ,
                                       gamma   = -1      , 
                                       name    = v2name  ,
                                       title   = v2title )
        return var1 , var2 
        
        
    # =========================================================================
    ## Soft/Gaussian constraints 
    # =========================================================================
    
    # =========================================================================
    ## Helper function to prepare 'soft' Gaussian constraint for the variable
    #  @attention the constraint is prepared, but not applied!
    #  @code
    #  sigma      = ...
    #  constraint = pdf.soft_constraint( sigma , VE ( 0.15 , 0.01**2 ) )
    #  @endcode 
    def soft_constraint ( self , var , value , name = '' , title = '' , error = None ) :
        """Prepare 'soft' Gaussian constraint for the variable
        -  consraint is prepared but not applied!
        >>> sigma      = ...
        >>> constraint = pdf.make_constraint( sigma , VE ( 0.15 , 0.01**2 ) )
        """
        
        assert isinstance ( var   , ROOT.RooAbsReal ) ,\
               "Invalid 'v': %s/%s"  % ( var , type ( var ) )

        if isinstance ( value , VE ) and 0 < value.cov2() and error is None : 
            error = value.error ()
            value = value.value ()

        assert isinstance ( value , ROOT.RooAbsReal ) or   isinstance ( value , num_types ) , \
               "Invalid 'value': %s/%s"  % ( value , type ( value ) )
               
        assert isinstance ( error , ROOT.RooAbsReal ) or ( isinstance ( error , num_types ) and 0 < error ) , \
               "Invalid 'error': %s/%s"  % ( error , type ( error ) )
        
        name  = name  if name  else self.generate_name ( name = 'Gauss_%s_%s' % ( var.GetName() , self.name ) )
        title = title if title else 'Gaussian Constraint(%s,%s) at %s' % ( var.GetName() , self.name , value )
        
        # value & error as RooFit objects:
        val = value if isinstance ( value , ROOT.RooAbsReal ) else ROOT.RooFit.RooConst ( value )
        err = error if isinstance ( error , ROOT.RooAbsReal ) else ROOT.RooFit.RooConst ( error )
        
        # Gaussian constrains 
        gauss = ROOT.RooGaussian ( self.roo_name ( name ) , title , var , val , err )
        
        # keep all the created technical stuff  
        self.aux_keep.append ( val   )
        self.aux_keep.append ( err   )
        self.aux_keep.append ( gauss )

        v = float ( value )
        e = float ( error ) 
        self.debug  ("Constraint is created '%s' : %s" % ( var.name , VE ( v , e * e ) ) ) 
        return  gauss 

    # =========================================================================
    ## Helper function to prepare 'soft' asymmetric Gaussian constraint for the variable
    #  @attention the constraint is prepared, but not applied!
    #  @code
    #  sigma      = ...
    #  constraint = pdf.soft_constraint2 ( sigma , 0.15 , -0.01 , +0.05 )
    #  @endcode 
    def soft_constraint2 ( self       ,
                           var        ,
                           value      ,
                           neg_error  ,
                           pos_error  , 
                           name  = '' ,
                           title = '' ) :
        """Prepare 'soft' asymetric Gaussian constraint for the variable
        -  consraint is prepared but not applied!
        >>> sigma      = ...
        >>> constraint = pdf.make_constraint2 ( sigma , 0.15 , -0.01 , +0.05 ) 
        """
        
        assert isinstance ( var   , ROOT.RooAbsReal ) ,\
               "Invalid 'v': %s/%s"  % ( var , type ( var ) )        
        assert isinstance ( value     , num_types ) ,\
               "Invalid 'value': %s/%s"  % ( value , type ( value ) )
        assert isinstance ( neg_error , num_types ) ,\
               "Invalid 'neg_error': %s/%s"  % ( neg_error , type ( neg_error ) )        
        assert isinstance ( pos_error , num_types ) and 0 < pos_error ,\
               "Invalid 'pos_error': %s/%s"  % ( pos_error , type ( pos_error ) )
        
        if abs ( neg_error ) == pos_error :
            return self.soft_constraint ( var , VE ( value , pos_error * pos_error ) ,
                                          name = name , title = title )
        
        name  = name  if name  else 'Gauss_%s_%s'                      % ( var.GetName() , self.name ) 
        title = title if title else 'Gaussian Constraint(%s,%s) at %s' % ( var.GetName() , self.name , value )
        
        # value & error as RooFit objects: 
        val    = ROOT.RooFit.RooConst (       float ( value     )   )
        poserr = ROOT.RooFit.RooConst (       float ( pos_error )   )
        negerr = ROOT.RooFit.RooConst ( abs ( float ( neg_error ) ) ) 
        
        # asymmetric/bifurcated  Gaussian constrains 
        agauss = Ostap.Models.BifurcatedGauss (
            self.var_name ( name ) , title , var , val , negerr , poserr )
        
        # keep all the created technical stuff  
        self.aux_keep.append ( val    )
        self.aux_keep.append ( poserr )
        self.aux_keep.append ( negerr )
        self.aux_keep.append ( agauss )
        
        self.debug ("Constraint is created '%s': %+.g+%-.g-%-.g" % ( var.name , value , pos_error , neg_error ) )        
        return agauss 

    # ==========================================================================
    ## Create multivariate Gaussian constraint
    #  @attention the constraint is prepared, but not applied!
    #  - use fit result 
    #  @code
    #  vars       = ...
    #  fit_result = ...
    #  constraint = pdf.soft_multivar_constraint ( vars , fitresult )
    #  @endcode
    #  - use values and covariance matrix 
    #  @code
    #  vars       = ...
    #  values     = ...
    #  cov2       = ... 
    #  constraint = pdf.soft_multivar_constraint ( vars , ( values , cov2 ) )    
    #  @see RooMultiVarGaussian 
    def soft_multivar_constraint ( self , vars , config , name = '' , title = '' ) :
        """Create multivariate Gaussian constraint
        - attention the constraint is prepared, but not applied!
        
        - use fit result:
        
        >>> vars       = ...
        >>> fit_result = ...
        >>> constraint = pdf.soft_multivar_constraint ( vars , fitresult )
        
        - use values and covariance matrix 
        >>> vars       = ...
        >>> values     = ...
        >>> cov2       = ... 
        >>> constraint = pdf.soft_multivar_constraint ( vars , ( values , cov2 ) )    
        
        - see `ROOT.RooMultiVarGaussian`
        """

        if isinstance ( config , ROOT.RooFitResult ) :
            return self.soft_multivar_constraint ( vars , ( config , True ) , name , title )

        variables = []
        for i, v in enumerate ( vars ) :
            assert v in self, 'Invalid/unknown parameter #%d %s/%s' % ( i , v , type(v) ) 
            variables.append ( self.parameter ( v ) )
        variables = tuple ( variables )


        vlist = ROOT.RooArgList()
        for v in variables : vlist.add ( v )
        
        name  = name  if name  else self.roo_name ( 'MVGauss' ) 
        title = title if title else 'MVGaussian Constraint(%s,%s) ' % ( self.name , self.name  )

        N = len ( vlist )

        import ostap.math.linalg

        if isinstance ( config , ROOT.RooFitResult ) :
            first , second = config , True
        elif isinstance ( config , Ostap.Math.VectorE ( N ) ) :
            first , second = config.value() , config.cov2() 
        elif isinstance ( config , sequence_types ) and \
             isinstance ( config , sized_types    ) and 2 == len ( config ) :            
            first , second = config
        else :
            "Invalid type of 'config' %s" % str ( config )

        if   isinstance ( first  , ROOT.RooFitResult          ) : pass 
        elif isinstance ( first  , ROOT.TVectorD              ) : pass 
        elif isinstance ( first  , Ostap.Math.Vector    ( N ) ) : first  = first .tvector()
        elif isinstance ( first  , sequence_types             ) and \
             isinstance ( first  , sized_types                ) and N == len ( first ) :
            vct = ROOT.TVectorD ( N ) 
            for i , v in enumerate ( first ) : vct [ i ] = float ( v )
            first = vct 

        if   isinstance ( second , Ostap.Math.SymMatrix ( N ) ) : second = second.tmatrix() 

        ## 1) the most useful case : values and covariances are provdeid throught RooFitResult obkect
        if isinstance ( first , ROOT.RooFitResult ) :
            
            fitres     = first 
            float_pars = fitres.floatParsFinal()
            for v in variables :
                assert v in float_pars , \
                       'FitResult object does not depend on parameter %s/%s' % ( v , type ( v ) )

        ## 2) values and covariances are provided explicitely 
        elif isinstance ( first , ROOT.TVectorD ) and isinstance ( second , ROOT.TMatrixDSym ) :

            vals = first
            cov2 = second 
            assert vals.IsValid() and cov2.IsValid()             , 'Invalid values/covariances'
            assert vals.GetNrows() == N and cov2.GetNrows() == N , 'Invalid seze for values/covariances'

        else :
            
            raise TypeError ( "Invalid 'config' (%s,%s)/(%s,%s)" % (
                first , second , type ( first ) , type ( second ) ) )

        ## create multivariate gaussian constraint  
        mgauss = ROOT.RooMultiVarGaussian ( name , title , vlist , first , second )
        
        self.aux_keep.append ( vlist  )
        self.aux_keep.append ( mgauss )
        
        self.debug ("Multivariate (%s) constraint is created %s" % ( [ v.name for v in variables ] , mgauss ) ) 
        return mgauss 
                        
    # ==========================================================================
    ## Helper function to  create soft Gaussian constraint
    #  to the ratio of <code>a</code> and <code>b</code>
    #  @code
    #  N1 = ...
    #  N2 = ...
    #  rC = pdf.soft_ratio_contraint ( N1 , N2 , VE (10,1**2) )
    #  @endcode
    def soft_ratio_constraint ( self , a , b , value , name = '' , title = '' ) :
        """Helper function to  create soft Gaussian constraint
        to the ratio of te variables 
        >>> N1 = ...
        >>> N2 = ...
        >>> rC = pdf.soft_ratio_constraint ( N1 , N2 , VE (10,1**2) )
        """
        assert isinstance ( value , VE ) and 0 < value.cov2() ,\
               "Invalid 'value': %s/%s"  % ( value , type ( value ) )
        
        if 1 < abs ( value.value() )  :
            return self.soft_ratio_constraint ( b , a , 1.0 / value )

        fa = isinstance ( a , ROOT.RooAbsReal )
        fb = isinstance ( b , ROOT.RooAbsReal )

        if   fa and fb :
            
            var = self.vars_ratio ( a , b )        
            return self.soft_constraint ( var , value , name , title )
        
        elif fa and isinstance ( b , num_types + (VE,) ) :
            
            return self.soft_constraint ( a , value * b , name , title )
        
        elif fb and isinstance ( a , num_types + (VE,) ) :
            
            return self.soft_constraint ( b , a / value  , name , title )
        
        raise TypeError('Unknown types a&b: %s/%s' % ( type ( a ) , type ( b ) ) )
    
    # ==========================================================================
    ## Helper function to  create soft Gaussian constraint
    #  to the product of <code>a</code> and <code>b</code>
    #  @code
    #  N1 = ...
    #  N2 = ...
    #  rC = pdf.soft_product_contraint ( N1 , N2 , VE (10,1**2) )
    #  @endcode
    def soft_product_constraint ( self , a , b , value , name = '' , title = '' ) :
        """Helper function to  create soft Gaussian constraint
        to the product of the variables 
        >>> N1 = ...
        >>> N2 = ...
        >>> rC = pdf.soft_product_constraint ( N1 , N2 , VE (10,1**2) )
        """
        assert isinstance ( value , VE ) and 0 < value.cov2() ,\
               "Invalid 'value': %s/%s"  % ( value , type ( value ) )

        fa = isinstance ( a , ROOT.RooAbsReal )
        fb = isinstance ( b , ROOT.RooAbsReal )

        if   fa and fb :
            
            var = self.vars_multiply ( a , b )        
            return self.soft_constraint ( var , value , name , title )
        
        elif fa and isinstance ( b , num_types + (VE,) ) :
            
            return self.soft_constraint ( a , value / b , name , title )
        
        elif fb and isinstance ( a , num_types + (VE,) ) :
            
            return self.soft_constraint ( b , value / a  , name , title )
        
        raise TypeError('Unknown types a&b: %s/%s' % ( type ( a ) , type ( b ) ) )

                        
    # ==========================================================================
    ## Helper function to  create soft Gaussian constraint
    #  to the fraction of <code>a</code> and <code>b</code> : a/(a+b)
    #  @code
    #  N1 = ...
    #  N2 = ...
    #  rC = pdf.soft_fraction_contraint ( N1 , N2 , VE (10,1**2) )
    #  @endcode
    def soft_fraction_constraint ( self , a , b , value , name = '' , title = '' ) :
        """Helper function to  create soft Gaussian constraint
        to the fraction of the variables:  a/(a+b)
        >>> N1 = ...
        >>> N2 = ...
        >>> rC = pdf.soft_product_constraint ( N1 , N2 , VE (10,1**2) )
        """
        assert isinstance ( value , VE ) and 0 < value.cov2() ,\
               "Invalid 'value': %s/%s"  % ( value , type ( value ) )

        return self.soft_ratio_constraint ( b , a , 1.0/value - 1.0 , name , value )
        
    # ==========================================================================
    ## Helper function to  create soft Gaussian constraint
    #  to the sum of <code>a</code> and <code>b</code>
    #  @code
    #  N1 = ...
    #  N2 = ...
    #  rC = pdf.soft_sum_contraint ( N1 , N2 , VE (10,1**2) )
    #  @endcode
    def soft_sum_constraint ( self , a , b , value , name = '' , title = '' ) :
        """Helper function to  create soft Gaussian constraint
        to the sum of the variables 
        >>> N1 = ...
        >>> N2 = ...
        >>> rC = pdf.soft_sum_constraint ( N1 , N2 , VE (10,1**2) )
        """

        fa = isinstance ( a , ROOT.RooAbsReal )
        fb = isinstance ( b , ROOT.RooAbsReal )

        if   fa and fb :
            
            var = self.vars_add ( a , b )        
            return self.soft_constraint ( var , value , name , title )
        
        elif fa and isinstance ( b , num_types + (VE,) ) :
            
            return self.soft_constraint ( a , value - b , name , title )
        
        elif fb and isinstance ( a , num_types + (VE,) ) :
            
            return self.soft_constraint ( b , value - a  , name , title )
        
        raise TypeError('Unknown types a&b: %s/%s' % ( type ( a ) , type ( b ) ) )
    
    # ==========================================================================
    ## Helper function to  create soft Gaussian constraint
    #  to the difference of <code>a</code> and <code>b</code>
    #  @code
    #  N1 = ...
    #  N2 = ...
    #  rC = pdf.soft_difference_contraint ( N1 , N2 , VE (10,1**2) )
    #  @endcode
    def soft_difference_constraint ( self , a , b , value , name = '' , title = '' ) :
        """Helper function to  create soft Gaussian constraint
        to the difference  of the variables 
        >>> N1 = ...
        >>> N2 = ...
        >>> rC = pdf.soft_sum_constraint ( N1 , N2 , VE (10,1**2) )
        """
        assert isinstance ( value , VE ) and 0 < value.cov2() ,\
               "Invalid 'value': %s/%s"  % ( value , type ( value ) )
        
        fa = isinstance ( a , ROOT.RooAbsReal )
        fb = isinstance ( b , ROOT.RooAbsReal )

        if   fa and fb :
            
            var = self.vars_subtract ( a , b )        
            return self.soft_constraint ( var , value , name , title )
        
        elif fa and isinstance ( b , num_types + (VE,) ) :
            
            return self.soft_constraint ( a , value + b , name , title )
        
        elif fb and isinstance ( a , num_types + (VE,) ) :
            
            return self.soft_constraint ( b , a - value  , name , title )
        
        raise TypeError('Unknown types a&b: %s/%s' % ( type ( a ) , type ( b ) ) )

    
    # =========================================================================
    ## create ready-to-use soft Gaussian constraint
    #       and wrap it to ROOT.RooFit.ExternalConstraint
    #  @see RooFit::ExternalConstraint
    #  @code
    #  sigma      = ...
    #  constraint = pdf.make_constraint( sigma , VE ( 0.15 , 0.01**2 ) )
    #  pdf.fitTo ( ... ,  constraints =  constraint ) 
    #  @endcode 
    def make_constraint ( self , var , value , name = '' ,  title = '' ) :
        """Create ready-to-use 'soft' gaussian constraint for the variable
        
        >>> var     = ...                              ## the variable 
        >>> extcntr = xxx.constraint ( VE(1,0.1**2 ) ) ## create constrains 
        >>> model.fitTo ( ... , constraint = extcntr ) ## use it in the fit 
        """
        
        ## create the gaussian constraint
        gauss  = self.soft_constraint ( var , value , name ,  title ) 
        
        cnts   = ROOT.RooArgSet ( gauss )
        
        result = ROOT.RooFit.ExternalConstraints ( cnts )
        
        self.aux_keep.append ( cnts   )
        
        return result 

    # =========================================================================
    ## sample 'random' positive number of events
    #  @code
    #  n =  pdf.gen_sample ( 10            ) ## get poissonian 
    #  n =  pdf.gen_sample ( VE ( 10 , 3 ) ) ## get gaussian stuff
    #  @endcode
    def gen_sample ( self , nevents ) :
        """Sample 'random' positive number of events
        >>> n =  pdf.gen_sample ( 10            ) ## get poissonian 
        >>> n =  pdf.gen_sample ( VE ( 10 , 3 ) ) ## get gaussian stuff
        """
        if   isinstance ( nevents , num_types ) and 0 < nevents :
            return poisson ( nevents )
        elif isinstance ( nevents , VE ) and \
                 ( ( 0 <= nevents.cov2 () and 0 < nevents                       ) or 
                   ( 0 <  nevents.cov2 () and 0 < nevents + 3 * nevents.error() ) ) :
            for i in range ( 20000 ) :
                n = int ( ve_gauss ( nEvents ) )
                if 0 < n : return n 
            else :
                self.error ( "Can't generate positive number from %s" % events )
                return
            
        self.error ( "Can't generate positive number from %s/%s" % ( events , type ( events ) ) )
        return 

    # =========================================================================
    ## get the proper min/max range for the variable  
    def vmnmx    ( self , var , vmin , vmax ) :
        """Get the proper min/max range for the variable 
        """
        if var.xminmax() :
            vmn , vmx = var.xminmax ()
            if   is_good_number ( vmin ) : vmin = max ( vmin , vmn )
            else                         : vmin = vmn
            if   is_good_number ( vmax ) : vmax = min ( vmax , vmx )
            else                         : vmax = vmx

        assert is_good_number ( vmin ), "Invalid type of 'min' %s/%s" % ( vmin , type ( vmin ) )
        assert is_good_number ( vmax ), "Invalid type of 'max' %s/%s" % ( vmax , type ( vmax ) )
        assert vmin < vmax, 'Invalid min/max range for %s: %s/%s' % ( var.name  , vmin , vmax )
        
        return vmin , vmax



    # ==========================================================================
    # Several purely technical methods 
    # ==========================================================================
        
    # ==========================================================================
    ## Pure technical methods to make a sum of PDFs with *constant* equal fractions
    #  @code 
    #  sum = self.make_sum ( 'A' , 'B' , pdf1 , pdf2 , pdf3 )
    #  @endcode
    #  It is very useful for e.g. creation of ""symmetrized" PDFs
    #  f(x,y) = 0.5 * [f1(x)*f2(y)] + 0.5* [ f2(x)*f1(y) ]    
    def make_sum ( self , name , title , pdf1 , pdf2 , *pdfs ) :
        """ Pure technical methods to make a sum of PDFs with *constant* equal fractions
        >>> sum = self.make_sum ( 'A' , 'B' , pdf1 , pdf2 , pdf3 )
        It is very useful for e.g. creation of 'symmetrized'PDF:
        f(x,y) = 0.5 * [f1(x)*f2(y)] + 0.5* [ f2(x)*f1(y) ]    
        """
        if self.name : name = name + '_' + self.name 
        _pdfs  = ROOT.RooArgList()
        _pdfs.add ( pdf1 )
        _pdfs.add ( pdf2 )
        for p in pdfs : _pdfs.add ( p )
        n = len(_pdfs)
        _fracs = [] 
        for i in  range(n) :
            f = ROOT.RooConstVar ( "Fraction%d_%s"     % ( i+1 , name         ) ,
                                   "fraction%d(%s,%s)" % ( i+1 , name , title ) , 1.0 / n )
            _fracs.append ( f )
            
        _rlst = ROOT.RooArgList()
        for f in _fracs : _rlst.add ( f ) 
        ## create PDF 
        result = ROOT.RooAddPdf ( name , title , _pdfs , _rlst , False )
        ##
        self.aux_keep.append ( _pdfs  )
        self.aux_keep.append ( _fracs )
        self.aux_keep.append ( _rlst  )
        #
        return result
    
    # =============================================================================
    ## make list of variables/fractions for compound PDFs 
    def make_fracs ( self             ,
                     N                ,
                     pname            ,
                     ptitle           ,
                     fractions = True ,
                     recursive = True ,
                     fracs     = []   )  :
        """Make list of variables/fractions for compound PDF
        """
        assert is_integer ( N ) and 2 <= N , \
               "PDF.make_fracs: Invalid N=%s/%s" % ( N, type ( N ) )
        ##
        ufracs   = [] 
        n        =  ( N - 1 ) if fractions else N
        NN       = n + 1
        vminmax  =  ( 0 , 1 ) if fractions else ( 0 , 1.e+7 )
        value    = 1 
        prod     = 1.0
        for i in range ( 0 , n ) :            
            if not fractions :                
                fv = 1.0 / NN
                if recursive :    
                    fv   /=  prod
                    prod *= ( 1.0 - fv )
                value = fv 
            ## finally create the fraction
            fi = get_i ( fracs , i , None )
            
            if isinstance ( fi , num_types ) and not vminmax[0] <= fi <= vminmax[1] :
                self.warning ("make_fracs: fraction %s is outside the interval %s, ignore" % ( fi , vminmax ) )
                fi = None
                
            f  = self.make_var ( fi , pname % i , ptitle % i , False , value , *vminmax ) 
            ufracs.append ( f )
            
        return ufracs
    
    # =============================================================================
    ## Create list of variables that can be used as 'fractions' for
    #  non-extedned N-component fit
    #  @code
    #  MV = ...
    #  fractions = MV.make_fractions ( 5 )
    #  fractions = MV.make_fractions ( 5 , name = 'F%d_A' , suffix = 'D' )
    #  fractions = MV.make_fractions ( 5 , name = 'F%d_B' , recursive = False )    
    #  fractions = MV.make_fractions ( 5 , name = 'F%d_C' , fractions = (0.4, 0.1 , 0.3 , 0.1 ) )   
    #  @endcode 
    def make_fractions ( self              ,
                         N                 ,
                         name      = 'f%d' , ## pattern to contruct the fraction name from index 
                         title     = ''    , 
                         recursive = True  ,
                         fractions = ()    ) :
        
        assert is_integer ( N ) and 2 <= N ,\
               "make_fractions: there must be at least two components!"

        my_fractions = make_iterable ( fractions , None )
            
        fracs = []
        value = 1.0
        prod  = 1.0


        for i , ff in zip ( range ( N - 1  ) , my_fractions ) :

            value = 1.0 / N
            
            if recursive :
                
                value /= prod
                prod  *= ( 1.0 - value ) 
                
            if   isinstance  ( ff , num_types       ) and not 0 <= ff <= 1 :
                self.error   ("make_fractions: fraction %s is outside [0,1] interval, ignore it!" %  ff )
                ff = value
                
            elif isinstance  ( ff , ROOT.RooAbsReal ) and not 0 <= float ( ff ) <= 1 :
                self.warning ("make_fractions: fraction %s is outside [0,1] interval" %  ff )
                
            if ff is None : ff = value

            fname = name if ( 2 == N and  '%s' not in name and '%d' not in name ) else name % i  
            
            tit   =  title if title else 'Fraction #%d: %s %s' % ( i , fname , self.name )                 
            fvar  = self.make_var ( ff   , fname  , tit , False , value , 0 , 1 )
            
            fracs.append ( fvar  )

        return tuple ( fracs )
          
    # =============================================================================
    ## create a list of variables that can be used as 'yields' for N-component fit
    #  @code
    #  M = ...
    #  nums = M.make_yields ( 5 , name = 'S_%d_A' , minmax = (0, 1000) )
    #  nums = M.make_yields ( 5 , name = 'S_%d_B' , minmax = (0, 1000) ,  yields =  ( 1 , 100 , 50 , 10 , 10 ) )
    #  @endcode     
    def make_yields ( self                   ,
                      N                      ,
                      name  = 'N%d'          , ## pattern to construct yield names from index
                      title = ''             ,
                      minmax = ( 0 , 1.e+6 ) , 
                      yields = ()            ) :
        """create a list of variables that can be used as 'yields' for N-component fit
        >>> M = ...
        >>> nums = M.make_yields ( 5 , name = 'S_%d_A' , minmax = (0, 1000) )
        >>> nums = M.make_yields ( 5 , name = 'S_%d_B' , minmax = (0, 1000) ,  yields =  ( 1 , 100 , 50 , 10 , 10 ) )
        """
        my_yields = make_iterable ( yields , None )

        nums = []        
        for i , n in zip ( range ( N ) , my_yields ) : 

            fname = ( name  % i )  if 1 != N else name           ## ATTENTION!  
            tit   = title if title else 'Yield #%d: %s %s' % ( i , fname , self.name )                 
            if not n is None : nvar  = self.make_var ( n , fname , tit , False , n , *minmax )
            else             : nvar  = self.make_var ( n , fname , tit , False ,     *minmax )
            
            nums.append ( nvar )

        return tuple ( nums ) 

    # =============================================================================
    ## helper function to build composite (non-extended) PDF from components 
    def add_pdf ( self             ,
                  pdflist          ,
                  name             ,
                  title            ,
                  fname            ,
                  ftitle           ,
                  recursive = True ,
                  fractions = []   ) :
        """Helper function to build composite (non-extended) PDF from components 
        """
        ##
        pdfs   = ROOT.RooArgList() 
        for pdf in pdflist : pdfs.add  ( pdf )
        fs     = self.make_fracs ( len ( pdfs )          ,
                                   fname                 ,
                                   ftitle                ,
                                   fractions = True      ,
                                   recursive = recursive ,
                                   fracs     = fractions )
        fracs  = ROOT.RooArgList()
        for f in fs : fracs.add ( f ) 
        pdf    = ROOT.RooAddPdf ( self.roo_name ( name ) , title , pdfs , fracs , recursive )
        ##
        self.aux_keep.append ( pdf   )
        self.aux_keep.append ( pdfs  )
        self.aux_keep.append ( fracs )
        ##
        return pdf , fracs , pdfs

# =============================================================================


# =============================================================================
## @class XVar
#  Helper MIXIN  class to keep all properties of the x-variable
class XVar(object) :
    """Helper MIXIN class to keep all properteis the x-variable
    """
    def  __init__ ( self , xvar , name = 'x' , title = 'x-observable' ) :

        self.__xvar = None

        if   isinstance ( xvar , ROOT.TH1   ) : xvar = xvar.xminmax()
        elif isinstance ( xvar , ROOT.TAxis ) : xvar = xvar.GetXmin() , xvar.GetXmax()
        
        ## create the variable 
        if isinstance ( xvar , tuple ) and ( 2 <= len ( xvar ) <= 6 ) :  
            self.__xvar = self.make_var ( xvar         , ## var 
                                          name         , ## name 
                                          title        , ## title/comment
                                          False        , ## fix ? 
                                          *xvar        ) ## min/max 
        elif isinstance ( xvar , ROOT.RooAbsReal ) :
            self.__xvar = self.make_var ( xvar         , ## var 
                                          name         , ## name 
                                          title        ) ## title/comment
        else :
            raise TypeError( "'x-variable'is not specified properly %s/%s" % ( xvar , type ( xvar ) ) )
                
    @property 
    def xvar ( self ) :
        """'x'-variable  (same as 'x')"""
        return self.__xvar
    @property 
    def x    ( self ) :
        """'x'-variable (same as 'xvar')"""
        return self.xvar
    def xminmax ( self ) :
        """Min/max values for x-variable (if/when specified)"""
        return self.xvar.minmax()

    ## get the proper xmin/xmax range
    def xmnmx    ( self , xmin , xmax ) :
        """Get the proper xmin/xmax range
        """
        return self.vmnmx ( self.xvar , xmin , xmax )



# =============================================================================
## @class YVar
#  Helper MIXIN class to keep all properties of the y-variable
#  - it requires the method <code>make_var</code>
class YVar(object) :
    """Helper class to keep all properteis the y-variable
    - it required the method `make_var`
    """
    def  __init__ ( self , yvar , name = 'y' , title = 'y-observable' ) :
         
        self.__yvar = None

        if   isinstance ( yvar , ROOT.TH1   ) : yvar = yvar.xminmax()
        elif isinstance ( yvar , ROOT.TAxis ) : yvar = yvar.GetXmin() , yvar.GetXmax()
        
        ## create the variable 
        if isinstance ( yvar , tuple ) and ( 2 <= len ( yvar ) <= 6 ) :  
            self.__yvar = self.make_var ( yvar         , ## var 
                                          name         , ## name 
                                          title        , ## title/comment
                                          False        , ## fix ? 
                                          *yvar        ) ## min/max 
        elif isinstance ( yvar , ROOT.RooAbsReal ) :
            self.__yvar = self.make_var ( yvar         , ## var 
                                          name         , ## name 
                                          title        ) ## title/comment
        else :
            raise TypeError( "'y-variable'is not specified properly %s/%s" % ( yvar , type ( yvar ) ) )
                
    @property 
    def yvar ( self ) :
        """'y'-variable (same as 'y')"""
        return self.__yvar
    @property 
    def y    ( self ) :
        """'y'-variable (same as 'yvar')"""
        return self.yvar
    def yminmax ( self ) :
        """Min/max values for y-variable (if/when specified)"""
        return self.yvar.minmax()

    ## get the proper ymin/ymax range
    def ymnmx    ( self , ymin , ymax ) :
        """Get the proper ymin/ymax range
        """
        return self.vmnmx ( self.yvar , ymin , ymax )

# =============================================================================
## @class ZVar
#  Helper MIXIN class to keep all properties of the z-variable
#  - it requires the method <code>make_var</code>
class ZVar(object) :
    """Helper class to keep all properteis the z-variable
    - it requires the method <code>make_var</code>
    """
    def  __init__ ( self , zvar , name = 'z' , title = 'z-observable' ):
         
        self.__zvar = None

        if   isinstance ( zvar , ROOT.TH1   ) : zvar = zvar.xminmax()
        elif isinstance ( zvar , ROOT.TAxis ) : zvar = zvar.GetXmin() , yvar.GetXmax()
        
        ## create the variable 
        if isinstance ( zvar , tuple ) and ( 2 <= len ( zvar ) <= 6 ) :   
            self.__zvar = self.make_var ( zvar         , ## var                                          
                                          name         , ## name 
                                          title        , ## title/comment
                                          False        , ## fix ? 
                                          *zvar        ) ## min/max 
        elif isinstance ( zvar , ROOT.RooAbsReal ) :
            self.__zvar = self.make_var ( zvar         , ## var 
                                          name         , ## name 
                                          title        ) ## title/comment
        else :
            raise TypeError( "'z-variable' is not specified properly %s/%s" % ( zvar , type ( zvar ) ) )
                
    @property 
    def zvar ( self ) :
        """'y'-variable (same as 'z')"""
        return self.__zvar
    @property 
    def z    ( self ) :
        """'z'-variable (same as 'zvar')"""
        return self.zvar
    def zminmax ( self ) :
        """Min/max values for y-variable (if/when specified)"""
        return self.zvar.minmax()

    ## get the proper xmin/xmax range
    def zmnmx    ( self , xmin , xmax ) :
        """Get the proper zmin/zmax range
        """
        return self.vmnmx ( self.zvar , zmin , zmax )


# =============================================================================
## simple convertor of 1D-histo to weighted or binned data set
#  @code
#  h   = ...
#  dset = H1D_dset ( h )
#  @endcode
#  One can create binned (default) or weighted datasets depending
#  on the value of <code>weighted</code> parameter
#  @code
#  h   = ...
#  wset = H1D_dset ( h , weighted = True  )
#  bset = H1D_dset ( h , weighted = False ) ## default 
#  @endcode
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H1D_dset(XVar,VarMaker) :
    """Simple convertor of 1D-histogram into weighted or binned data set
    >>> h   = ...
    >>> dset = H1D_dset ( h )
    One can create `binned` (default) or `weighted` data set
    >>> wset = H1D_dset ( h , weighted = True  )
    >>> bset = H1D_dset ( h , weighted = False ) ## default 
    """
    
    w_min = -1.e+100 
    w_max =  1.e+100
    
    def __init__ ( self              , 
                   histo             ,
                   xaxis     = None  , ## predefined axis/variable ? 
                   density   = False ,
                   weighted  = False , ## weighted or binned?
                   skip_zero = False , ## skip zero bins for weighted dataset? 
                   silent    = False ) :

        import ostap.histos.histos
        
        #
        ## use mass-variable
        #
        assert histo and isinstance ( histo , ROOT.TH1 ) and 1 == histo.dim () , "'histo' is not ROOT.TH1"

        ## set the base from VarMaker 
        self.name = self.new_name ( 'H1D_dset(%s)' % histo.GetName() )
        
        if not xaxis : xaxis = histo.xminmax() 

        ## initialize the XVar base class 
        XVar.__init__ ( self  , xaxis  ,
                        name  = 'x_%s'           % self.name ,
                        title = 'x-observable%s' % self.name )

        
        self.__histo      = histo 
        self.__histo_hash = hash ( histo )        
        self.__skip_zero  = True if skip_zero else False
        self.__density    = True if density   else False 
        self.__silent     = silent 
        
        self.__wvar       = None
        
        with roo_silent ( self.silent ) :  

            ## create weighted dataset ?
            if weighted :
                
                wname = weighted if isinstance ( weighted , string_types ) else 'h1weight'
                
                self.__wvar  = ROOT.RooRealVar  ( wname , "weight-variable" , 1 , self.w_min , self.w_max )
                self.__wname = self.__wvar.GetName() 
                self.__vset  = ROOT.RooArgSet   ( self.xvar   ,  self.__wvar )
                self.__wset  = ROOT.RooArgSet   ( self.__wvar )
                self.__warg  = ROOT.RooFit.WeightVar ( self.__wvar ) , ROOT.RooFit.StoreError ( self.__wset )
                self.__dset  = ROOT.RooDataSet  (
                    rootID ( 'whds_' )  , "Weighted data set for the histogram '%s'" % histo.GetTitle() ,
                    self.__vset , *self.__warg )

                xvar = self.xvar   
                wvar = self.__wvar 
                with SETVAR ( xvar ) :
                    for i, x , v in histo.items () :
                        if skip_zero and 0 == v.value() and 0 == v.error () : continue 
                        xvar.setVal     ( x.value () )
                        self.__dset.add ( self.__vset , v.value() , v.error() ) 
                        
            ## create binned dataset 
            else :
                
                self.__vlst = ROOT.RooArgList    ( self.xvar  )
                self.__vimp = ROOT.RooFit.Import ( self.histo , self.density )
                self.__dset = ROOT.RooDataHist   (
                    rootID ( 'bhds_' ) , "Binned data set for histogram '%s'" % histo.GetTitle() ,
                    self.__vlst  ,
                    self.__vimp  )
                
    @property     
    def xaxis ( self ) :
        """The histogram x-axis variable (same as 'xvar')"""
        return self.xvar 
    @property
    def histo ( self ) :
        """The  histogram itself"""
        return self.__histo
    @property
    def density( self ) :
        """Treat the histo as 'density' histogram?"""
        return self.__density
    @property
    def skip_zero ( self ) :
        """'skip_zero' : skip zero bins for weighted dataset in histo?"""
        return self.__skip_zero    
    @property
    def silent( self ) :
        """Use the silent mode?"""
        return self.__silent
    @property
    def dset ( self ) :
        """'dset' : ROOT.RooDataHist object"""
        return self.__dset
    @property
    def histo_hash ( self ) :
        """Hash value for the histogram"""
        return self.__histo_hash
    @property
    def weight ( self ) :
        """'weight' : get weight variable if defined, None otherwise"""
        return self.__wvar
    @property
    def wname  ( self ) :
        """'wname' : het weight name (if defined)"""
        return None if ( self.__wvar is None ) else self.__wvar.GetName() 
    
# =============================================================================
## simple convertor of 2D-histo to weighted or binned data set
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H2D_dset(XVar,YVar,VarMaker) :
    """Simple convertor of 2D-histogram into weighted or binned data set
    """
    
    w_min = -1.e+100 
    w_max =  1.e+100
    
    def __init__ ( self              ,
                   histo             ,
                   xaxis     = None  ,
                   yaxis     = None  ,
                   density   = False ,
                   weighted  = False ,
                   skip_zero = False , ## skip zero bins for weighted dataset? 
                   silent    = False ) :
        #
        import ostap.histos.histos
        
        assert histo and isinstance ( histo , ROOT.TH2 ) and 2 == histo.dim() , "'histo' is not ROOT.TH2"

        ## set the base from VarMaker 
        self.name = self.new_name ( 'H2D_dset(%s)' % histo.GetName() )

        if not xaxis : xaxis = histo.xminmax() 
        if not yaxis : yaxis = histo.yminmax() 

        ## initialize the XVar base class 
        XVar.__init__ ( self   , xaxis  ,
                        name  = 'x_%s'           % self.name ,
                        title = 'x-observable%s' % self.name )
        
        ## initialize the YVar base class 
        YVar.__init__ ( self  , yaxis  ,
                        name  = 'y_%s'           % self.name ,
                        title = 'y-observable%s' % self.name )

        
        self.__histo      =        histo
        self.__histo_hash = hash ( histo )

        self.__skip_zero = True if skip_zero else False
        self.__density   = True if density else False 
        self.__silent    = silent
        
        self.__wvar    = None

        with roo_silent ( silent ) : 

            ## create weighted dataset 
            if weighted :
                
                wname = weighted if isinstance ( weighted , string_types ) else 'h2weight'
                
                self.__wvar = ROOT.RooRealVar  ( wname , "weight-variable" , 1 , self.w_min , self.w_max )
                self.__vset = ROOT.RooArgSet   ( self.xaxis ,  self.yaxis , self.__wvar )
                self.__wset = ROOT.RooArgSet   ( self.__wvar )
                self.__warg = ROOT.RooFit.WeightVar ( self.__wvar ) , ROOT.RooFit.StoreError ( self.__wset )
                self.__dset = ROOT.RooDataSet (
                    rootID ( 'whds_' )  , "Weighted data set for the histogram '%s'" % histo.GetTitle() ,
                    self.__vset , *self.__warg )

                xvar = self.xvar
                yvar = self.yvar
                wvar = self.__wvar 
                with SETVAR ( xvar ) :
                    for ix , iy , x , y , v in histo.items () :
                        if skip_zero and 0 == v.value() and 0 == v.error () : continue 
                        xvar.setVal     ( x.value () )
                        yvar.setVal     ( y.value () )
                        self.__dset.add ( self.__vset , v.value() , v.error() ) 

            ## create binned dataset 
            else : 
                self.__vlst  = ROOT.RooArgList    ( self.xaxis , self.yaxis )
                self.__vimp  = ROOT.RooFit.Import ( histo , density )
                self.__dset  = ROOT.RooDataHist   (
                    rootID ( 'bhds_' ) , "Binned sata set for histogram '%s'" % histo.GetTitle() ,
                    self.__vlst  ,
                    self.__vimp  )
            
    @property     
    def xaxis  ( self ) :
        """The histogram x-axis variable (same as 'xvar')"""
        return self.xvar
    @property     
    def yaxis  ( self ) :
        """The histogram y-axis variable (same as 'yvar')"""
        return self.yvar
    @property
    def histo ( self ) :
        """The  histogram itself"""
        return self.__histo
    @property
    def density( self ) :
        """Treat the histo as 'density' histogram?"""
        return self.__density    
    @property
    def skip_zero ( self ) :
        """'skip_zero' : skip zero bins for weighted dataset in histo?"""
        return self.__skip_zero    
    @property
    def silent( self ) :
        """Use the silent mode?"""
        return self.__silent
    @property
    def dset ( self ) :
        """'dset' : ROOT.RooDataHist object"""
        return self.__dset
    @property
    def histo_hash ( self ) :
        """Hash value for the histogram"""
        return self.__histo_hash
    @property
    def weight ( self ) :
        """'weight' : get weight variable if defined, None otherwise"""
        return self.__wvar
    
# =============================================================================
## simple convertor of 3D-histo to weighted or binned data set
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H3D_dset(XVar,YVar,ZVar,VarMaker) :
    """Simple convertor of 3D-histogram into weighted or binned data set
    """
    
    w_min = -1.e+100 
    w_max =  1.e+100
    
    def __init__ ( self              ,
                   histo             ,
                   xaxis     = None  ,
                   yaxis     = None  ,
                   zaxis     = None  ,
                   density   = False ,
                   weighted  = False ,
                   skip_zero = False , ## skip zero bins for weighted dataset? 
                   silent    = False ) :
        
        import ostap.histos.histos

        assert histo and isinstance ( histo , ROOT.TH3 ) and 3 == histo.dim () , "'histo' is not ROOT.TH3"

        ## set the base from VarMaker 
        self.name = self.new_name ( 'H3D_dset(%s)' % histo.GetName() )

        if not xaxis : xaxis = histo.xminmax()
        if not yaxis : yaxis = histo.yminmax()
        if not zaxis : zaxis = histo.zminmax()
        
        ## initialize the XVar base class 
        XVar.__init__ ( self  , xaxis  ,
                        name  = 'x_%s'           % self.name ,
                        title = 'x-observable%s' % self.name )
        
        ## initialize the YVar base class 
        YVar.__init__ ( self  , yaxis  ,
                        name  = 'y_%s'           % self.name ,
                        title = 'y-observable%s' % self.name )
        
        ## initialize the YVar base class 
        ZVar.__init__ ( self  , zaxis  ,
                        name  = 'z_%s'           % self.name ,
                        title = 'z-observable%s' % self.name )


        self.__histo      =        histo
        self.__histo_hash = hash ( histo )
        #
        
        self.__skip_zero = True if skip_zero else False
        self.__density   = True if density else False 
        self.__silent    = silent
        
        self.__wvar    = None
        
        with roo_silent ( silent ) : 

            ## create weighted dataset 
            if weighted :
                
                wname = weighted if isinstance ( weighted , string_types ) else 'h2weight'
                
                self.__wvar = ROOT.RooRealVar  ( wname , "weight-variable" , 1 , self.w_min , self.w_max )
                self.__vset = ROOT.RooArgSet   ( self.xaxis ,  self.yaxis , self.zaxis, self.__wvar )
                self.__wset = ROOT.RooArgSet   ( self.__wvar )
                self.__warg = ROOT.RooFit.WeightVar ( self.__wvar ) , ROOT.RooFit.StoreError ( self.__wset )
                self.__dset = ROOT.RooDataSet  (
                    rootID ( 'whds_' )  , "Weighted data set for the histogram '%s'" % histo.GetTitle() ,
                    self.__vset , *self.__warg )

                xvar = self.xvar
                yvar = self.yvar
                zvar = self.zvar 
                wvar = self.__wvar 
                with SETVAR ( xvar ) :
                    for ix , iy , iz , x , y , z , v in histo.items () :
                        if skip_zero and 0 == v.value() and 0 == v.error () : continue 
                        xvar.setVal     ( x.value () )
                        yvar.setVal     ( y.value () )
                        zvar.setVal     ( z.value () )
                        self.__dset.add ( self.__vset , v.value() , v.error() ) 

            ## create binned dataset 
            else : 
                                
                self.__vlst  = ROOT.RooArgList    ( self.xaxis , self.yaxis , self.zaxis )
                self.__vimp  = ROOT.RooFit.Import ( histo , density )
                self.__dset  = ROOT.RooDataHist   (
                    rootID ( 'bhds_' ) , "Binned data set for histogram '%s'" % histo.GetTitle() ,
                    self.__vlst  ,
                    self.__vimp  )
            
    @property     
    def xaxis  ( self ) :
        """The histogram x-axis variable (same as 'xvar')"""
        return self.xvar
    @property     
    def yaxis  ( self ) :
        """The histogram y-axis variable (same as 'yvar')"""
        return self.yvar    
    @property     
    def zaxis  ( self ) :
        """The histogram z-axis variable (same as 'zvar')"""
        return self.zvar
    @property
    def histo ( self ) :
        """The  histogram itself"""
        return self.__histo
    @property
    def density( self ) :
        """Treat the histo as 'density' histogram?"""
        return self.__density    
    @property
    def skip_zero ( self ) :
        """'skip_zero' : skip zero bins for weighted dataset in histo?"""
        return self.__skip_zero    
    @property
    def silent( self ) :
        """Use the silent mode?"""
        return self.__silent
    @property
    def dset ( self ) :
        """'dset' : ROOT.RooDataHist object"""
        return self.__dset
    @property
    def histo_hash ( self ) :
        """Hash value for the histogram"""
        return self.__histo_hash
    @property
    def weight ( self ) :
        """'weight' : get weight variable if defined, None otherwise"""
        return self.__wvar
    
# =============================================================================
## @class ParamsPoly
#  Helper MIXIN class to implement polynomials
#  - it requres the method <code>make_var</code>
#  - it requres the method <code>component_setter</code>
class ParamsPoly(object) :
    """Helper MIXIN class to implement polynomials 
    - it requres the method `make_var`
    - it requres the method `component_setter`
    """
    def __init__ ( self , npars = 1 , pars = None , limits = ( 0 , -1.e+6 , 1.e+6 ) ) :
        
        from ostap.math.base import isint as _isint 
        assert pars or ( ( isinstance ( npars , integer_types ) or _isint ( npars ) ) and 0 <= npars ) ,\
               "ParamsPoly: Inconsistent 'npars' setting"
        
        npars = int ( npars )

        self.__pars     = [] 

        if not limits : limits = 0 , -1.e+6 , 1.e+6
            
        if isinstance ( pars , ParamsPoly      ) :
            
            self.__pars = [ p for p in pars.pars ]
            npars       = len ( self.__pars ) 
            
        else : 

            newpars     = make_iterable ( pars , size = npars )
            self.__pars = [ self.make_var ( p ,
                                            'par%d_%s'            % ( i , self.name ) ,
                                            'parameter %d for %s' % ( i , self.name ) ,
                                            False , *limits ) for ( i , p ) in enumerate ( newpars ) ]

        self.__pars     = tuple ( self.__pars ) 
        self.__pars_lst = ROOT.RooArgList()
        for p in self.pars : self.__pars_lst.add ( p )
        
        self.config = {
            'name'   : self.name ,
            'xvar'   : self.xvar ,
            'pars'   : self.pars ,
            'limits' : limits 
            }
        
    ## release parameter
    #  - 0 <= index < N : release  parameter 
    #  - index      < 0 : release all pamametrs
    #  - N < index      : no actions 
    def release_par ( self , index = -1 ) :
        """Release parameter
        - 0 <= index < N : release  parameter 
        - index      < 0 : release all pamametrs
        - N < index      : no actions 
        """
        n = len ( self.__pars )
        if index < 0        :
            for p in self.__pars : p.release () 
        elif 0 <= index < n : self.__pars[index].release()

    @property
    def pars  ( self ) :
        """'pars' : the polynomial coefficients/parameters"""
        return self.__pars    
    @pars.setter
    def pars ( self , values ) :
        self.component_setter ( self.__pars , values )

    def reset_pars ( self , value = 0 ) :
        """Set all pars to be value 
        >>> pdf = ...
        >>> pdf.reset_pars() 
        """
        for f in self.__pars : f.setVal ( value )

    @property
    def pars_lst ( self ) :
        """'pars_lst' : the polynomial coefficients/parameters as RooArgList"""
        return self.__pars_lst
    
    @property
    def npars ( self ) :
        """'npars'  : number of parameters """
        return len  ( self.pars ) 

# =============================================================================
## @class Phases
#  helper MIXIN class to build/keep the list of 'phi'-arguments
#   - it is needed for polynomial functions
#   - it requires mehtod <code>make_var</code>
#   - it requires mehtod <code>error</code>
#   - it requires mehtod <code>component_setter</code>
class Phases(object) :
    """Helper MIXIN class to build/keep the list of 'phi'-arguments,
    - it is needed for polynomial functions
    - it requires method `make_var`
    - it requires method `error`
    - it requires method `component_setter`
    """
    ## Create vector of phases (needed for various polynomial forms)
    def __init__( self  , power , the_phis = None ) :
        """Create vector of phases (needed for various polynomial forms)
        """

        from ostap.math.base import isint as _isint 
        ## check  the arguments 
        assert ( isinstance ( power , integer_types ) or _isint ( power ) ) and 0 <= power, \
               "Phases: invalid type/value for 'power'-parameter: %s/%s"  % (  power , type ( power ) ) 
        
        power = int ( power ) 

        from math import pi
        limits = 0 , -5 * pi , 5 * pi

        if isinstance ( the_phis , Phases          ) :
            ## copy phases 
            self.__phis = [ i for i in the_phis.phis ]  
            power       = len ( self.__phis )            
        else :
            ## mer new phases or copy them 
            newphis     = make_iterable ( the_phis , size = power )
            self.__phis = [ self.make_var ( phi ,
                                            'phi%d%s'       % ( i , self.name ) ,
                                            '#phi_{%d}(%s)' % ( i , self.name ) ,
                                            False , *limits ) for ( i , phi ) in enumerate ( newphis )  ]
            
        self.__phis     = tuple ( self.__phis )
        self.__phi_list = ROOT.RooArgList ()
        for p in self.phis : self.__phi_list.add ( p ) 
        
    # =========================================================================
    ## release parameter
    #  - 0 <= index < N : release  parameter 
    #  - index      < 0 : release all pamametrs
    #  - N < index      : no actions 
    def release_par ( self , index = -1 ) :
        """Release parameter
        - 0 <= index < N : release  parameter 
        - index      < 0 : release all pamametrs
        - N < index      : no actions 
        """
        n = len ( self.__pars )
        if index < 0        :
            for p in self.__pars : p.release () 
        elif 0 <= index < n : self.__pars[index].release()


    ## set all phis to be 0
    def reset_phis ( self ) :
        """Set all phases to be zero
        >>> pdf = ...
        >>> pdf.reset_phis() 
        """
        for f in self.__phis : f.setVal(0)
        
    @property
    def phis ( self ) :
        """The list/tuple of 'phases', used to parameterize various polynomial-like shapes
        >>> pdf = ...
        >>> for phi in pdf.phis :
        ...    print phi       ## get phase  
        ...    print phi.value ## get phase value 
        ...    phi.value = 0.1 ## set phase value 
        ...    phi.fix(0)      ## fix phase 
        
        """
        return tuple ( self.__phis )

    @phis.setter
    def phis ( self , values ) :
        self.component_setter ( self.__phis , values )

    @property
    def phi_list ( self ) :
        """The list/ROOT.RooArgList of 'phases', used to parameterize polynomial-like shapes
        """
        return self.__phi_list
    
    @property
    def phis_lst ( self ) :
        """The list/ROOT.RooArgList of 'phases', used to parameterize polynomial-like shapes
        """
        return self.__phi_list

    @property
    def power ( self ) :
        """'power'  : polynomial degree """
        return len ( self.__phis ) 
    

# =============================================================================
## @class ShiftScalePoly
#  Helper MIXIN class to implement polynomials
#  \f$ f(x) = a + b P(x) \f$,
#  where \f$P(x)\f$ some special polynomial
#  @see Phases
#  - it requires the method <code>make_var</code>
#  - it requires the method <code>error</code>
class ShiftScalePoly ( Phases ) :
    """Helper MIXIN class to implemnet polynomials
    - see Phases 
    - it requires the method `make_var`
    - it requires the method `error`
    """
    def __init__ ( self         ,
                   a     = 0.0  , ## shift/bias
                   b     = 1.0  , ## scale 
                   power = 1    ,
                   pars  = None ) :
        
        ## initialize the base class 
        Phases.__init__ ( self , power , pars ) 
                          
        ## parameter a 
        self.__a  = self.make_var ( a ,
                                    'a_%s'              % self.name ,
                                    'bias/shift for %s' % self.name , False )
        ## parameter a 
        self.__b  = self.make_var ( b ,
                                    'b_%s'              % self.name ,
                                    'scale for %s'      % self.name , False )
        
        self.config = {
            'name'  : self.name         ,
            'xvar'  : self.xvar         ,
            'pars'  : self.pars         ,
            'a'     : self.a            ,
            'power' : len ( self.phis ) , 
            'b'     : self.b            }
        
    @property
    def a ( self ) :
        """'a' : bias parameter for polynomial:  f(x) = a + b*M(x)"""
        return self.__a
    @a.setter
    def a ( self , value ) :
        vv = float ( value )
        if self.__a.minmax () and not vv in self.__a  :
            self.error ("Value %s is outside the allowed region %s for %s"  % ( vv , self.__a.minmax() , self.__a.name ) )
        self.__a.setVal ( vv )

    @property
    def b ( self ) :
        """'scale' : bias parameter for polynomial:  f(x) = a + b*M(x)"""
        return self.__b
    @b.setter
    def b ( self , value ) :
        vv = float ( value )
        if self.__b.minmax () and not vv in self.__b  :
            self.error ("Value %s is outside the allowed region %s for %s"  % ( vv , self.__b.minmax() . self.__b.name ) )
        self.__b.setVal ( vv )
    
    @property
    def shift  ( self ) :
        """'shift' : bias/shift parameter  (same as 'a')"""
        return self.__a
    @shift.setter
    def shift  ( self , value ) : self.a = value 
    @property
    def bias  ( self ) :
        """'bias' : bias/shift parameter  (same as 'a')"""
        return self.__a
    @bias.setter
    def bias   ( self , value ) : self.a = value
    
    @property
    def scale ( self ) :
        """'scale' : bias/shift parameter  (same as 'b')"""
        return self.__b
    @scale.setter
    def scale ( self , value ) : self.b = value

    @property
    def pars  ( self ) :
        """'pars' :  polynomial parameters (same as 'phis')"""
        return self.phis
    @pars.setter 
    def pars  ( self , values ) :
        self.phis = values 
    @property
    def pars_lst ( self ) :
        """'pars_lst' :  polynomial parameters as RooArgList"""
        return self.phis_lst
    

# =============================================================================
## @class Fractions
#  Helper MIXIN class for implementatiorn of SumXD objects
class Fractions(object) :
    """Helper MIXIN class for implementatiorn of SumXD objects
    """
    def __init__  ( self             ,
                    pdfs             , ## list of PDF objects 
                    prefix    = 'f'  ,                    
                    suffix    = ''   , 
                    recursive = True ,
                    fractions = None )  :

        assert 2 <= len ( pdfs ) , 'Fractions: at least two PDFs are needed!'
        
        self.__pdfs = tuple ( pdfs )

        ## check
        for i , p in enumerate ( self.pdfs )  :
            if   p.pdf.mustBeExtended() : self.warning ("'pdf%f' is 'must be extended'!" % i ) 
            elif p.pdf. canBeExtended() : self.warning ("'pdf%f' is 'can  be extended'!" % i ) 
                    
        while prefix.endswith  ('_') : prefix = prefix[:-1]
        while suffix.startswith('_') : suffix = suffix[1:]
        
        self.__prefix    = prefix if prefix    else 'f' 
        self.__suffix    = suffix
        self.__recursive = True   if recursive else False 
        
        fr_name = make_name ( self.prefix , '%d' if 2 < len ( self.pdfs ) else '' , self.suffix )

        ## make list of fractions 
        self.__fractions = tuple ( self.make_fractions  ( len ( self.pdfs )           ,
                                                          name       = fr_name        , 
                                                          recursive  = self.recursive ,
                                                          fractions  = fractions      ) ) 
    @property
    def pdfs ( self ) :
        """'pdfs' : get list/tuple of involved PDFs (same as 'components')"""
        return self.__pdfs

    @property
    def pdf1 ( self ) :
        """'pdf1' : the first PDF"""
        return self.__pdfs[0] if 1<= len ( self.__pdfs ) else None 
    
    @property
    def pdf2 ( self ) :
        """'pdf2' : the second PDF"""
        return self.__pdfs[1] if 2<= len ( self.__pdfs ) else None 
    
    @property 
    def tail( self ) :
        """'tail' : other PDFs (if any)"""
        return self.pdfs[2:] if 2<= len ( self.__pdfs ) else () 

    @property
    def components ( self ) :
        """'components' : get list/tuple of involved PDFs (same as 'pdfs')"""
        return self.pdfs
        
    @property
    def prefix ( self ) :
        """'fr_prefix' : prefix for fraction names"""
        return self.__prefix

    @property
    def suffix ( self ) :
        """'suffix' : suffix for fraction names"""
        return self.__suffix

    @property
    def recursive ( self ) :
        """'recursive' : recursive fractions?"""
        return self.__recursive

    @property
    def frac_list ( self ) :
        """'frac_list': get fractions as list"""
        return self.__fractions
    
    @property
    def fractions ( self ) :
        """'fractions' : get involved fractions (same as 'F' or 'fraction')"""
        return self.component_getter ( self.__fractions )
    @fractions.setter
    def fractions ( self , values ) :
        self.component_setter ( self.__fractions , values )

    @property
    def fraction ( self ) :
        """'fraction' : get involved fractions (same as 'F' of 'fractions')"""
        return self.fractions 
    @fraction.setter
    def fraction ( self , values ) :
        self.fractions = values 

    @property
    def F         ( self ) :
        """'F' : get involved fractions (same as 'fractions' or 'fraction')"""
        return self.fractions
    @F.setter
    def F         ( self , values ) :
        self.fractions = values 
       
        
        
                    
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
    def __init__ ( self , fun , dataset = None , silent = True ) :

        self.__params  = {}
        self.__fun     = fun        
        self.__dataset = ROOT.nullptr if ( dataset is None ) else dataset 
        self.__silent  = True if silent else False 
        
    ## context manager: ENTER 
    def __enter__ ( self ) :

        params = self.__fun.parameters ( self.__dataset )
        for par in params :
            self.__params [ par ] = float ( params [ par ] )
            
        return self
    
    ## context manager: EXIT
    def __exit__ ( self , *_ ) :
        
        if self.__params :
            self.__fun.load_params ( params = self.__params , dataset = self.__dataset , silent = self.__silent )
            
        self.__fun     = None
        self.__params  = {}
        self.__dataset = ROOT.nullptr 
        
    @property
    def fun ( self ) :
        """'fun': the actual function/pdf"""
        return self.__fun
    
    @property
    def dataset ( self ) :
        """'dataset': the dataset"""
        return self.__dataset

    @property
    def params( self ) :
        """'params': dictionary of parameters"""
        return self.__params


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
#                                                                       The END 
# =============================================================================
