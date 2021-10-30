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
    'cov_qual'          , ## Quality of covariance matrix from Minuit
    'fit_status'        , ## Fit status from Minuit
    #
    'RangeVar'          , ## Helper class to temporary change a range for the variable 
    'MakeVar'           , ## Helper bade class that allow storage of newly created RooFit objects
    'XVar'              , ## helper to deal with x-variable
    'YVar'              , ## helper to deal with y-variable
    'ZVar'              , ## helper to deal with z-variable
    #
    "H1D_dset"          , ## 1D-histogram to RooDataHist converter 
    "H2D_dset"          , ## 2D-histogram to RooDataHist converter 
    "H3D_dset"          , ## 3D-histogram to RooDataHist converter
    #
    'component_similar' , ## Should one use ``similar'' component?
    'component_clone'   , ## Should one use ``cloned'' component?
    # 
    'numcpu'            , ## number of CPUs
    'ncpu'              , ## fuction to build ROOT.RooFit.NumCPU
    ## add/remove RooFit topic
    'remove_topic'      , ## remove topic from RooMsgService
    'add_topic'         , ## add    topic from RooMsgService
    'suppress_topics'   , ## suppress topics from RooMsgService 
    #
    'Phases'            , ##  helper class for Ostap polynomial/PDFs
    'ParamsPoly'        , ##  helper class for RooFit polynomials
    'ShiftScalePoly'    , ##  helper class for RooFit polynomials
    #
    "NameDuplicates"    , ## allow/disallow name duplicates
    )
# =============================================================================
import ROOT, math, random 
import ostap.fitting.variables 
import ostap.fitting.roocollections
from   builtins                import range 
from   ostap.core.core         import Ostap, rootID, VE, items_loop, isequal 
from   ostap.core.ostap_types  import ( num_types      , list_types     ,
                                        integer_types  , string_types   ,
                                        is_good_number , sequence_types )
from   ostap.logger.utils      import roo_silent
from   sys                     import version_info as python_version 
from   ostap.math.random_ext   import ve_gauss, poisson
from   ostap.core.meta_info    import root_version_int
from   ostap.fitting.variables import SETVAR 
from   ostap.fitting.roocmdarg import flat_args, check_arg 
# =============================================================================
from   ostap.logger.logger     import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.utils' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
try :
    from string import ascii_letters
except ImportError :
    from string import letters as ascii_letters
# =============================================================================
## make a name from prefix, name and suffix 
def make_name ( prefix , name , suffix ) :
    """Make a name from prefix, name and suffix
    """
    
    prefix = prefix.replace ( ' ' , '' ) 
    suffix = suffix.replace ( ' ' , '' )
    name   = name  .replace ( ' ' , '' )

    while prefix.endswith   ( '_' ) : prefix = prefix[:-1]
    while suffix.startswith ( '_' ) : suffix = suffix[1:] 

    if   prefix and name and suffix : return "%s_%s_%s" % ( prefix , name , suffix ) 
    elif prefix and suffix          : return "%s_%s"    % ( prefix ,        suffix ) 
    elif prefix and name            : return "%s_%s"    % ( prefix , name          ) 
    elif suffix and name            : return "%s_%s"    % (          name , suffix )

    return "%s" % ( name or prefix or suffix ) 

    
# =============================================================================
## MINUIT covariance matrix status:
# - status = -1 :  not available (inversion failed or Hesse failed)
# - status =  0 : available but not positive defined
# - status =  1 : covariance only approximate
# - status =  2 : full matrix but forced pos def
# - status =  3 : full accurate matrix
_cov_qual_ = {
    -1 :  '-1/not available (inversion failed or Hesse failed or externally provided)' ,
    0  :  ' 0/available but not positive defined',
    1  :  ' 1/covariance only approximate',
    2  :  ' 2/full matrix but forced pos def',
    3  :  ' 3/full accurate matrix',
    }
# =============================================================================
## MINUIT covariance matrix status:
# - status = -1 : not available (inversion failed or Hesse failed)
# - status =  0 : available but not positive defined
# - status =  1 : covariance only approximate
# - status =  2 : full matrix but forced pos def
# - status =  3 : full accurate matrix
def cov_qual ( status ) : return _cov_qual_.get( status , "%s" % status )
# =============================================================================
## Miniut::minimize status code
# - status = 1    : Covariance was made pos defined
# - status = 2    : Hesse is invalid
# - status = 3    : Edm is above max
# - status = 4    : Reached call limit
# - status = 5    : Any other failure
_fit_status_ = {
    0    : ' 0/success' ,
    1    : ' 1/Covariance was made pos defined',
    2    : ' 2/Hesse is invalid',
    3    : ' 3/Edm is above max',
    4    : ' 4/Reached call limit',
    5    : ' 5/Any other failure',
    }
# =============================================================================
## Miniut::minimize status code
# - status = 1    : Covariance was made pos defined
# - status = 2    : Hesse is invalid
# - status = 3    : Edm is above max
# - status = 4    : Reached call limit
# - status = 5    : Any other failure
def fit_status ( status ) : return _fit_status_.get( status ,"%s" % status )
# =============================================================================
_nemax = 1000 ## number of events per CPU-core 
_ncmax =   16 ## maximal number of CPUs: there are some problems with >= 7
              ## @see https://sft.its.cern.ch/jira/browse/ROOT-4897
# ==============================================================================
_ncpus = []
## Get number of cores/CPUs
def  numcpu () :
    """Get number of cores/CPUs
    """
    import multiprocessing
    return multiprocessing.cpu_count()

# =============================================================================
## prepare "NumCPU" argument with reasonable choice of #cpu, depending on
#  number of events in dataset 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-03-31
def ncpu ( events ) :
    """Prepare 'NumCPU' argument with reasonable choice of #cpu, depending on
    the number of events in dataset 
    """
    #
    n_cores = numcpu() 
    if n_cores <= 1 : return ROOT.RooFit.NumCPU ( 1 ) ## fake!!! 
    #
    n  = events // _nemax
    if n       <= 1 : return ROOT.RooFit.NumCPU ( 1 ) ## fake!!! 
    #
    num = min ( n , n_cores , _ncmax )
    if not _ncpus : _ncpus.append ( num )   
    #
    return ROOT.RooFit.NumCPU ( num )

# =============================================================================
## helper class to temporary change a range for the variable 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class RangeVar(object) :
    """Helper class to temporary change a range for the variable 
    """
    def __init__( self , var , vmin , vmax ) :
        self.var  = var
        self.vmin = min ( vmin , vmax ) 
        self.vmax = max ( vmin , vmax )
        self.omin = self.var.getMin ()
        self.omax = self.var.getMax ()
        
    def __enter__ ( self ) :
        self.omin = self.var.getMin ()
        self.omax = self.var.getMax ()
        self.var.setMin ( self.vmin ) 
        self.var.setMax ( self.vmax )
        return self
    
    def __exit__  ( self , *_ ) :        
        self.var.setMin ( self.omin ) 
        self.var.setMax ( self.omax )

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


# =============================================================================
## keep the list of local loggers  
_loggers  = {}           
# =============================================================================
## @class MakeVar
#  Helper class that allows implement several purely  technical methods:
#   - creation of <code>ROOT.RooRealVar objects</code>
#   - store newly created RooFit objects
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2018-07-14
class MakeVar ( object ) :
    """Helper class that allows implement several purely  technical methods:
    - creation of <code>ROOT.RooRealVar objects</code>
    - store newly created RooFit objects
    """
    __pdf_names = set()
    __var_names = set()
    __numnames  = 0
    
    ## @attention ensure that important attributes are available even before __init__
    def __new__( cls, *args, **kwargs):
        if  python_version.major > 2 : obj = super(MakeVar, cls).__new__( cls )
        else                         : obj = super(MakeVar, cls).__new__( cls , *args , **kwargs )
        ##
        obj.__aux_keep     = []                     ## ATTENTION!        
        obj.__name        = None                    ## ATTENTION!
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

    @property
    def aux_keep ( self ) :
        """``aux_keep'' -  the list of objects to be kept by this PDF/FUN"""
        return self.__aux_keep
    @property
    def logger   ( self ) :
        """``logger'': get the local logger object"""
        name    = self.name
        logname = str ( self.__class__.__name__ ) 
        if name : logname += '(' + name + ')'
        _logger = _loggers.get ( logname , None )
        if not _logger :
            _logger              = getLogger ( logname )
            _loggers [ logname ] = _logger
        return _logger
    @property
    def name ( self ) :
        """The name of the object"""
        return self.__name if self.__name else '' 
    @name.setter
    def name ( self , value ) :
        assert isinstance ( value , str ) , "``name'' must  be a string, %s/%s is given" % ( value , type ( value ) )
        if self.__name == value : return 
        if value in self.__pdf_names and not NameDuplicates.allowed() :
            self.warning ( 'The name "%s" for PDF/FUN already defined!' % value )
        self.__pdf_names.add ( value )     
        self.__name = value

    ## # =============================================================================
    ## ## generate some unique name
    ## @classmethod 
    ## def generate_name ( cls , prefix = '' , suffix = '' ) :
    ##     name = prefix + suffix 
    ##     while name in cls.__pdf_names or name in cls.__var_names or not name :
    ##         name = prefix + ''.join ( ( random.choice ( ascii_letters ) for i in range ( 6 ) )  ) + suffix 
    ##     return name

    # =============================================================================
    ## generate some unique name for PDF/FUN and objects
    @staticmethod 
    def generate_name ( prefix = '' , suffix = '' ) :
        """Generate some unique name for PDF/FUN and objects
        """
                
        prefix = prefix.replace ( ' ' , '' )
        if prefix.endswith('_') : prefix = prefix[:-1]
        
        suffix = suffix.replace ( ' ' , '' )
        if suffix.startswith('_') : suffix = suffix[1:] 

        name = make_name ( prefix , '' , suffix ) 
                                                   
        MakeVar.__numnames += 1            
        while ( name in MakeVar.__pdf_names ) or ( name in MakeVar.__var_names ) or ( not name ) :
            
            part = ''.join ( ( random.choice ( ascii_letters ) for i in range ( 6 ) )  ) + suffix

            name = make_name ( prefix , part , suffix )
            
            MakeVar.__numnames += 1

        return name
    
    # =============================================================================
    ## generate some unique name for <code>RooFit</code>
    #  @see TNamed 
    #  @see RooNameReg 
    #  @see RooAbsArg 
    @staticmethod
    def roo_name ( prefix = 'roo' , suffix = '' ) :
        """Generate some unique name for <code>RooFit</code>
        - see `ROOT.RooNameReg` 
        - see `ROOT.TNamed`
        - see `ROOT.RooAbsArg`
        """

        regname = ROOT.RooNameReg.instance()

        prefix = prefix.replace ( ' ' , '' )
        if prefix.endswith('_') : prefix = prefix[:-1]
        
        suffix = suffix.replace ( ' ' , '' )
        if suffix.startswith('_') : suffix = suffix[1:] 

        name = make_name ( prefix , '' , suffix ) 
                                                   
        MakeVar.__numnames += 1            
        while ( name in MakeVar.__pdf_names ) or ( name in MakeVar.__var_names ) or ( not name ) or regname.known ( name ) : 
            
            part = ''.join ( ( random.choice ( ascii_letters ) for i in range ( 6 ) )  ) + suffix
            
            name = make_name ( prefix , part , suffix )
            
            MakeVar.__numnames += 1

        return name
            
    # =============================================================================
    ## create/modify  the variable
    #  Helper function for creation/modification/adjustment of variable
    #  @code
    #  v = self.make_var ( 10   , 'myvar' , 'mycomment' )
    #  v = self.make_var ( 10   , 'myvar' , 'mycomment' , '' ,     -1 , 1 )
    #  v = self.make_var ( 10   , 'myvar' , 'mycomment' , '' , 0 , -1 , 1 )
    #  v = self.make_var ( None , 'myvar' , 'mycomment' , '' , 0 , -1 , 1 )
    #  v = self.make_var ( None , 'myvar' , 'mycomment' , 10 , 0 , -1 , 1 )
    #  v = self.make_var ( v    , 'myvar' , 'mycomment' , 10 )
    #  @endcode
    #  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
    #  @date 2013-12-01
    def make_var ( self           ,
                   var            ,
                   name           , ## name 
                   comment = ''   , ## title 
                   fix     = None , *args ) :
        """Make/modify  the variable:
        
        v = self.make_var ( 10   , 'myvar' , 'mycomment' )
        v = self.make_var ( 10   , 'myvar' , 'mycomment' , '' ,     -1 , 1 )
        v = self.make_var ( 10   , 'myvar' , 'mycomment' , '' , 0 , -1 , 1 )
        v = self.make_var ( None , 'myvar' , 'mycomment' , '' , 0 , -1 , 1 )
        v = self.make_var ( None , 'myvar' , 'mycomment' , 10 , 0 , -1 , 1 )
        v = self.make_var ( v    , 'myvar' , 'mycomment' , 10 )        
        """
        # var = ( value )
        # var = ( min , max )
        # var = ( value , min , max )

        if   isinstance   ( var , tuple ) :
            assert name and isinstance ( name , string_types ) , "make_var: invalid name '%s'" % name
            var     = ROOT.RooRealVar ( self.var_name ( name ) , comment , *var )
            self.debug ( 'Create variable/1:  %s' % var ) 
            self.aux_keep.append ( var ) ##  ATTENTION: store newly created variable

        ## if only name is specified :
        if   isinstance  ( var , string_types ) and 1 <= len ( args ) <= 3 :
            assert name and isinstance ( name , string_types ) , "make_var: invalid name '%s'" % name
            var     = ROOT.RooRealVar( self.var_name ( var ) , name + comment , *args )
            self.debug ( 'Created variable/2: %s' % var ) 
            self.aux_keep.append ( var ) ##  ATTENTION: store newly created variable
            
        # var = value 
        if isinstance   ( var , num_types ) :
            assert name and isinstance ( name , string_types ) , "make_var: invalid name '%s'" % name
            if   not    args       :
                var = ROOT.RooRealVar ( self.var_name ( name ) , comment , var             )
                self.aux_keep.append ( var ) ##  ATTENTION: store newly created variable                           
                self.debug ( 'Create variable/3: %s' % var ) 
            elif 2 == len ( args ) :
                var = ROOT.RooRealVar ( self.var_name ( name ) , comment , var , *args     )
                self.aux_keep.append ( var ) ##  ATTENTION: store newly created variable
                self.debug ( 'Create variable/4: %s' % var ) 
            elif 3 == len ( args ) :
                var = ROOT.RooRealVar ( self.var_name ( name ) , comment , var , *args[1:] )
                self.aux_keep.append  ( var ) ##  ATTENTION: store newly created variable
                self.debug ( 'Create variable/5: %s' % var ) 

        ## create the variable from parameters 
        if not isinstance ( var , ROOT.RooAbsReal ) : 
            assert name and isinstance ( name , string_types ) , "make_var: invalid name '%s'" % name
            var = ROOT.RooRealVar ( self.var_name ( name ) , comment , *args )
            self.aux_keep.append ( var ) ##  ATTENTION: store newly created variable
            self.debug ( 'Create variable/6: %s' % var ) 
            
        ## fix it, if needed
        if   isinstance ( fix , bool       ) : pass 
        elif isinstance ( fix , num_types  ) :
                
            if hasattr ( var , 'getMin)') and fix < var.getMin() and hasattr ( var , 'setMin' ) :                                                                             
                self.warning ( "make_var: min-value for %s is redefined to be %s" % ( var.GetName () , fix ) )
                var.setMin ( fix )
            
            if hasattr ( var , 'getMax' ) and fix > var.getMax() and hasattr ( var , 'setMax' ) : 
                self.warning ( "make_var: max-value for %s is redefined to be %s" % ( var.GetName () , fix ) )
                var.setMax ( fix )
            
            if not var.isConstant () and hasattr ( var , 'fix'    ) : var.fix    ( fix )
            elif                         hasattr ( var , 'setVal' ) : var.setVal ( fix )

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
        
        for k , a in items_loop ( kwargs ) :
            
            klow = k.lower ().replace('_','')
            kup  = k.upper ().replace('_','')
            
            ## skip "drawing" options 
            if   klow in drawing_options                            : continue 
            if   klow in ( 'draw'            ,
                           'drawoption'      ,
                           'drawoptions'     ) : continue 
            
            if   isinstance ( a , ROOT.RooCmdArg ) : _args.append ( a )
            
            elif kup in ( 'VERBOSE' ,        ) and isinstance ( a , bool ) :
                
                if not verbose is None :
                    if a != verbose : 
                        logger.warning ( 'parse_args: Redefine VERBOSE to %s' %  a ) 
                        verbose = a                        
                if not silent is None :
                    if a == silent :
                        logger.warning ( 'parse_args: confusing VERBOSE/SILENT %s/%s' % ( a , silent ) )
                        silent = not a 
                _args.append ( ROOT.RooFit.Verbose (     a ) )
            elif kup in ( 'SILENT'           ,
                          'SILENCE'          ) and isinstance ( a , bool ) :
                if not silent is None :
                    if a != silent : 
                        logger.warning ( 'parse_args: Redefine SILENT to %s' %  a ) 
                        verbose = a                        
                if not verbose is None :
                    if a == verbose :
                        logger.warning ( 'parse_args: confusing SILENT/VERBOSE %s/%s' % ( a , verbose ) )
                        verbose = not a
                _args.append ( ROOT.RooFit.Verbose ( not a ) ) 
            elif kup in ( 'STRATEGY'         , 
                          'MINUITSTRATEGY'   ,
                          'STRATEGYMINUIT'   ) and isinstance ( a , integer_types ) and 0 <= a <= 2 : 
                _args.append ( ROOT.RooFit.Strategy (    a ) ) 
            elif kup in ( 'PRINTLEVEL'       ,
                          'MINUITPRINT'      ,
                          'MINUITLEVEL'      ) and isinstance ( a , integer_types ) and -1 <= a <= 3 :
                _args.append ( ROOT.RooFit.PrintLevel ( a ) ) 
            elif kup in ( 'PRINTEVALERRORS'  ,
                          'PRINTERRORS'      ,
                          'ERRORSPRINT'      ) and isinstance ( a , integer_types ) and -1 <= a :
                _args.append ( ROOT.RooFit.PrintEvalErrors ( a ) )                
            elif kup in ( 'TIMER'            ,
                          'TIMING'           ) and isinstance ( a , bool ) :
                _args.append ( ROOT.RooFit.Timer    ( a ) ) 
            elif kup in ( 'WARNING'          ,
                          'WARNINGS'         ) and isinstance ( a , bool ) :
                _args.append ( ROOT.RooFit.Warnings ( a ) ) 
            
            elif kup in ( 'SUMW2'            ,
                          'SUMW2ERR'         ,
                          'SUMW2ERROR'       ,
                          'SUMW2ERRORS'      ) and isinstance ( a , bool ) :
                
                if   a and dataset and     dataset.isWeighted()           : pass 
                elif a and dataset and not dataset.isWeighted()           :
                    self.warning ('parse_args: SumW2-flag is True  for non-weighted dataset')
                elif       dataset and not dataset.isWeighted() and not a : pass 
                elif       dataset and     dataset.isWeighted() and not a :
                    self.warning ('parse_args: SumW2-flag is False for     weighted dataset')                    

                _args.append (  ROOT.RooFit.SumW2Error( a ) )
                                    
            elif kup in ( 'ASYMPTOTIC'       ,
                          'ASYMPTOTICERR'    ,
                          'ASYMPTOTICERROR'  ,
                          'ASYMPTOTICERRORS' ) and isinstance ( a , bool ) and 61900 <= root_version_int :
                
                if   a and dataset and     dataset.isWeighted()           : pass 
                elif a and dataset and not dataset.isWeighted()           :
                    self.warning ('parse_args: AsymptoticError-flag is True  for non-weighted dataset')
                elif       dataset and not dataset.isWeighted() and not a : pass 
                elif       dataset and     dataset.isWeighted() and not a :
                    self.warning ('parse_args: AsymptoticError-flag is False for     weighted dataset')                    

                if a and root_version_int < 62006 :
                    self.warning ("``Asymptotic=True'' will crash if Title!=Name (ROOT-10668)")
                    
                _args.append (  ROOT.RooFit.AsymptoticError ( a ) )
                    
            elif kup in ( 'BATCH'            ,
                          'BATCHMODE'        ) and isinstance ( a , bool ) and 62000 <= root_version_int :
                _args.append (  ROOT.RooFit.BatchMode ( a ) )                                
            elif kup in ( 'EXTENDED' ,       ) and isinstance ( a , bool ) :
                _args.append   (  ROOT.RooFit.Extended ( a ) )                
            elif kup in ( 'CPU'              ,
                          'CPUS'             ,
                          'NCPU'             ,
                          'NCPUS'            ,
                          'NUMCPU'           ,
                          'NUMCPUS'          ) and isinstance ( a , int ) and 1<= a : 
                _args.append   (  ROOT.RooFit.NumCPU( a  ) ) 
            elif kup in ( 'CPU'              ,
                          'CPUS'             ,
                          'NCPU'             ,
                          'NCPUS'            ,
                          'NUMCPU'           ,
                          'NUMCPUS'          ) and \
                 isinstance ( a , list_types ) and 2 == len ( a )  and \
                 isinstance ( a[0] , integer_types ) and 1 <= a[1] and \
                 isinstance ( a[1] , integer_types ) and 0 <= a[1] <=3 :
                _args.append   (  ROOT.RooFit.NumCPU( a[0] ,  a[1] ) ) 
                
            elif kup in ( 'RANGE'            ,
                          'FITRANGE'         ,
                          'RANGES'           ,
                          'FITRANGES'        ) and isinstance ( a , string_types ) :
                _args.append   (  ROOT.RooFit.Range ( a ) )  
            elif kup in ( 'RANGE'            ,
                          'FITRANGE'         ) and isinstance ( a , list_types   ) \
                 and isinstance ( a[0] ,  num_types ) \
                 and isinstance ( a[1] ,  num_types ) \
                 and a[0] < a[1]  : 
                _args.append   (  ROOT.RooFit.Range ( a[0] , a[1] ) )
            elif kup in ( 'MINIMIZER'  ,     ) and isinstance ( a , list_types   ) \
                 and isinstance ( a[0] ,  string_types ) \
                 and isinstance ( a[1] ,  string_types ) :
                _args.append   (  ROOT.RooFit.Minimizer ( a[0] , a[1] ) )                 
            elif kup in  ( 'HESSE'    ,      ) and isinstance ( a , bool ) :
                _args.append   (  ROOT.RooFit.Hesse ( a )  )
            elif kup in  ( 'INITIALHESSE'    ,
                           'INITHESSE'       ,
                           'HESSEINIT'       ,
                           'HESSEINITIAL'    ) and isinstance ( a , bool ) :
                _args.append   (  ROOT.RooFit.InitialHesse ( a )  )
            elif kup in ( 'OPTIMIZE'         ,
                          'OPTIMISE'         ) and isinstance ( a , integer_types  ) :
                _args.append   (  ROOT.RooFit.Optimize     ( a )  )
            elif kup in ( 'MINOS'    ,       ) and isinstance ( a , bool           ) :
                _args.append   (  ROOT.RooFit.Minos        ( a )  )
            elif kup in ( 'MINOS'    ,       ) and isinstance ( a , ROOT.RooArgSet ) :
                _args.append   (  ROOT.RooFit.Minos        ( a )  )
            elif kup in ( 'MINOS'    ,       ) and isinstance ( a , string_types   ) \
                     and hasattr  ( self , 'params' ) and a in self.params ( dataset ) :                
                _v = self.params()[ a ]
                _s = ROOT.RooArgSet ( _v )
                self.aux_keep.append ( _s ) 
                _args.append   (  ROOT.RooFit.Minos        ( _s )  )                
            elif kup in ( 'MINOS'    ,       ) and not isinstance ( a , string_types ) :

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
                                
            elif kup in ( 'SAVE'     ,       ) and isinstance ( a , bool           ) :
                _args.append   (  ROOT.RooFit.Save         ( a )  )
            elif kup in ( 'CLONE'            ,
                          'CLONEDATA'        ) and isinstance ( a , bool           ) :
                _args.append   (  ROOT.RooFit.CloneData    ( a )  )
            elif kup in ( 'OFFSET'           ) and isinstance ( a , bool           ) :
                _args.append   (  ROOT.RooFit.Offset       ( a )  )
            elif kup in ( 'FITOPTIONS'       ,
                          'FITOPTION'        ) and isinstance ( a , string_types ) :
                _args.append   (  ROOT.RooFit.FitOptions   ( a )  )
                
            elif kup in ( 'CONSTRAINT'       ,
                          'CONSTRAINTS'      ,
                          'PARS'             ,
                          'PARAMS'           ,
                          'PARAMETER'        ,
                          'PARAMETERS'       ) :
                c = self.parse_constraints ( a )
                if c is None : self.error ('parse_args: Invalid constraint specification: %s/%s' % ( a , type ( a ) ) )
                else         : _args.append ( c ) 
                    
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
        if dataset :
            
            weighted = dataset.isWeighted ()            
            sw2      = check_arg ( 'SumW2Error'      , *_args )
            aer      = check_arg ( 'AsymptoticError' , *_args )

            if 62500 <= root_version_int                 and \
               isinstance ( dataset , ROOT.RooDataHist ) and ( not sw2 ) and ( not aer ) :
                _args.append ( ROOT.RooFit.SumW2Error ( True ) )                
            elif sw2 and aer :
                logger.warning ( "parse_args: Both ``SumW2Error'' and ``AsymptoticError'' are specified" )                
            elif weighted   and sw2 :
                value = bool ( sw2.getInt( 0 ) )
                if not value : logger.warning ("parse_args: 'SumW2=False' is specified for the weighted  dataset!")
            elif weighted and aer : 
                value = bool ( aer.getInt( 0 ) )
                if not value : logger.warning ("parse_args: 'AsymptoticError=False' is specified for the weighted  dataset!")
            elif weighted :                
                logger.warning ( "parse_args: Neither ``SumW2Error'' and ``AsymptoticError'' are specified for weighted dataset! ``SumW2=True'' is added" )
                _args.append ( ROOT.RooFit.SumW2Error ( True ) )                
            elif not weighted and sw2 :
                logger.warning ( "parse_args:``SumW2Error'' is specified for non-weighted dataset" )
            elif not weighted and aer :
                logger.warning ( "parse_args:``AsymptoticError'' is specified for non-weighted dataset" )

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
        assert isinstance ( var , ROOT.RooAbsRealLValue ) , 'Invalid type of ``var'' %s' % type ( var )
        
        if not hasattr ( var ,  'setVal' ) :
            raise ValueError ( "No value can be set for %s/%s" % ( var , type ( var ) ) )  

        ## convert to float 
        value = float ( value )

        ## check for the range, if defined 
        minmax = var.minmax ()
        if minmax :
            mn , mx = minmax
            if not ( mn <= value <= mx or isequal ( mn , value ) or isequal ( mx , value ) ) :
                raise ValueError ( "Value %s is outside of the [%s,%s] region" % ( value , mn , mx ) ) 
            
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
            self.warning ("component setter: unknown type for ``value'':%s/%s" % ( str( value) , type ( value ) ) )
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
    #  var3 = xxx.vars_multiply ( var2 )
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
            return ROOT.RooRealConstant.value ( res ) 
        elif f1 :
            # shortcut 
            if   1 == var1 : return var2                             ## SHORTCUT
            elif 0 == var1 : return ROOT.RooRealConstant.value ( 0 ) ## SHORTCUT
            # 
            var1 = ROOT.RooRealConstant.value ( var1 ) 
            return self.vars_multiply ( var1 , var2 , name , title )
        elif f2 : 
            # shortcut 
            if   1 == var2 : return var1                             ## SHORTCUT
            elif 0 == var2 : return ROOT.RooRealConstant.value ( 0 ) ## SHORTCUT
            # 
            var2 = ROOT.RooRealConstant.value ( var2 ) 
            return self.vars_multiply ( var1 , var2 , name , title )
        
        self.aux_keep.append ( var1 )
        self.aux_keep.append ( var2 )

        result = Ostap.MoreRooFit. Product  ( var1 , var2 )
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
            return ROOT.RooRealConstant.value ( 0 ) 
        elif 0 == c1 : 
            return self.vars_multiply ( var2 , c2 , name = name , title = title )
        elif 0 == c2 :
            return self.vars_multiply ( var1 , c1 , name = name , title = title )
        
        f1 = isinstance ( var1 , num_types )
        f2 = isinstance ( var2 , num_types )

        if f1 and f2 :
            res  = float ( var1 ) * float ( c1 ) + float ( var2 ) * float ( c2 )
            return ROOT.RooRealConstant.value ( res ) 
        elif f1 :
            ## shortcut 
            if 0 == var1 :
                return self.var_multiply ( var2 , c2 , name = name , title = title )  ## SHORTCUT 
            #
            var1 = ROOT.RooRealConstant.value ( float ( var1 ) * float ( c1 ) )                         
            return self.vars_add ( var1 , var2 , name = name , title = title )
        elif f2 :
            ## shortcut 
            if 0 == var2 :
                return self.var_multiply ( var1 , c1 , name = name , title = title )  ## SHORTCUT 
            #
            var2 = ROOT.RooRealConstant.value ( float ( var2 ) * float ( c2 ) ) 
            return self.vars_add ( var1 , var2 , name =name , title = title  )
        
        self.aux_keep.append ( var1 )
        self.aux_keep.append ( var2 )

        if c1 == 1 and c2 == 1 : 
            result = Ostap.MoreRooFit.Addition  (                var1 , var2           )
        else :
            result = Ostap.MoreRooFit.Addition  ( name , title , var1 , var2 , c1 , c2 )
                
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
            return ROOT.RooRealConstant.value ( res ) 
        elif f1 :
            ## 
            var1 = ROOT.RooRealConstant.value ( var1 )                         
            return self.vars_subtract ( var1 , var2 , name , title )
        elif f2 :
            ## shortcut 
            if 0 == var2 : return var1                      ## SHORTCUT
            #
            var2 = ROOT.RooRealConstant.value ( var2 ) 
            return self.vars_subtract ( var1 , var2 , name , title )

        self.aux_keep.append ( var1 )
        self.aux_keep.append ( var2 )

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
            return ROOT.RooRealConstant.value ( res ) 
        elif f1 :
            var1 = ROOT.RooRealConstant.value ( var1 ) 
            return self.vars_divide   ( var1 , var2 , name , title )
        elif f2 :
            return self.vars_multiply ( var1 , 1.0/var2 , name , title )
        
        self.aux_keep.append ( var1 )
        self.aux_keep.append ( var2 )

        result = Ostap.MoreRooFit.Division  ( var1 , var2 )
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
            return ROOT.RooRealConstant.value ( res ) 
        elif f1 :
            ## shortcut 
            if 0 == var1  : return ROOT.RooRealConstant.value ( 0 ) ## SHORTCUT
            #
            var1 = ROOT.RooRealConstant.value ( var1 ) 
            return self.vars_fraction ( var1 , var2 , name , title )
        elif f2 :
            ## shortcut
            if 0 == var2  : return ROOT.RooRealConstant.value ( 1 ) ## SHORTCUT
            #
            var2 = ROOT.RooRealConstant.value ( var2 ) 
            return self.vars_fraction ( var1 , var2 , name , title )
        
        self.aux_keep.append ( var1 )
        self.aux_keep.append ( var2 )

        result = Ostap.MoreRooFit.Fraction  ( var1 , var2 )
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
            return ROOT.RooRealConstant.value ( res ) 
        elif f1 :
            ## shortcut 
            if 0 == var1 : return ROOT.RooRealConstant.value ( -1.0 * scale ) ## shortcut
            #
            var1 = ROOT.RooRealConstant.value ( var1 ) 
            return self.vars_asymmetry ( var1 , var2 , scale = scale , name = name , title = title )
        elif f2 :
            ## shortcut 
            if 0 == var2 : return ROOT.RooRealConstant.value (  1.0 * scale ) ## shortcut
            #
            var2 = ROOT.RooRealConstant.value ( var2 ) 
            return self.vars_asymmetry ( var1 , var2 , scale = scale , name = name , title = title )
        
        self.aux_keep.append ( var1 )
        self.aux_keep.append ( var2 )

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
            return ROOT.RooRealConstant.value ( res ) 
        elif f1 :
            ## shortcut 
            if 1 == var1 : return ROOT.RooRealConstant.value ( 1 ) ## shortcut
            #
            var1 = ROOT.RooRealConstant.value ( var1 ) 
            return self.vars_power ( var1 , var2 , name , title )
        elif f2 :
            ## shortcut 
            if 0 == var2 : return ROOT.RooRealConstant.value (  1 ) ## shortcut
            #
            var2 = ROOT.RooRealConstant.value ( var2 ) 
            return self.vars_power ( var1 , var2 , name , title )
        
        self.aux_keep.append ( var1 )
        self.aux_keep.append ( var2 )

        result = Ostap.MoreRooFit.Power ( var1 , var2 )
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
            return ROOT.RooRealConstant.value ( res ) 
        elif f1 :
            ## shortcut 
            if 0 == var1 : return ROOT.RooRealConstant.value ( 1 ) ## shortcut
            #
            var1 = ROOT.RooRealConstant.value ( var1 ) 
            return self.vars_exp ( var1 , var2 , name , title )
        elif f2 :
            ## shortcut 
            if 0 == var2 : return ROOT.RooRealConstant.value ( 1 ) ## shortcut
            #
            var2 = ROOT.RooRealConstant.value ( var2 ) 
            return self.vars_exp ( var1 , var2 , name , title )
        
        self.aux_keep.append ( var1 )
        self.aux_keep.append ( var2 )

        result = Ostap.MoreRooFit.Exp ( var1 , var2 )
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
        
        ## rfv = ROOT.RooFormulaVar ( self.var_name ( name ) , title , formula_ , vlst )
        rfv = Ostap.FormulaVar ( self.var_name ( name ) , title , formula_ , vlst )
        
        self.aux_keep.append ( vlst )
        self.aux_keep.append ( rvf  )
        
        return rfv 

    # =========================================================================
    ## make a specific combination of variables:  alpha*var1*(bet+gamma*var2)    
    #  \f$ r = \alpha v_1 ( \beta + \gamma * v_2 ) \    
    def vars_combination ( self ,
                           var1 ,
                           var2 ,
                           alpha  = 1   ,
                           beta   = 1   ,
                           gamma  = 1   ,
                           name   = ''  , 
                           title  = ''  ) :
        
        f1 = isinstance ( var1 , num_types )
        f2 = isinstance ( var2 , num_types )

        assert isinstance ( alpha , num_types ) , "vars_combination: ``alpha'' must be numeric types!"
        assert isinstance ( beta  , num_types ) , "vars_combination: ``alpha'' must be numeric types!"
        assert isinstance ( gamma , num_types ) , "vars_combination: ``alpha'' must be numeric types!"

        if   0 == alpha               : return ROOT.RooRealConstant.value ( 0 ) 
        elif 0 == beta and 0 == gamma : return ROOT.RooRealConstant.value ( 0 ) 

        alpha = float ( alpha )
        beta  = float ( beta  )
        gamma = float ( gamma )
        
        
        if f1 and f2 :

            res  = beta   + gamma * float ( var2 )
            res *=          alpha * float ( var1 )
            
            return ROOT.RooRealConstant.value ( 0 )
        
        elif f1 :
            
            return self.sum_add      ( float ( var1 ) * alpha * beta        ,
                                       var2                                 ,
                                       c2    = alpha * gamma * float ( v1 ) , 
                                       name  = name                         ,
                                       title = title                        ) 
        elif f2 :
            
            return self.sum_multiply ( var1 ,
                                       alpha * ( beta + gamma * float ( v2 ) ) ,
                                       name  = name ,
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
    ## Convert pair of variables into ``sum'' & ``asymmetry'' pair:
    # - ``sum''       :  sum_scale  * ( var1 + var2  )
    # - ``asymmetry'' :  asym_scale * ( var1 - var2 ) / (var1+ + var2 ) 
    # @code
    # var1 = ...
    # var2 = ...
    # hsum , asum = self.vars_to_asymmetry ( var1 , var2 ) 
    # @endcode 
    def vars_to_asymmetry ( self             ,
                            var1             ,   ## the first variable 
                            var2             ,   ## the second variable
                            asym_scale = 1   ,   ## scale factor asymmetry
                            sum_scale  = 0.5 ,   ## scale factor for ``sum''
                            asym_name  = ''  ,   ## name for asymmetry variable
                            sum_name   = ''  ,   ## name for ``sum'' variable                              
                            asym_title = ''  ,   ## title for asymmetry variable 
                            sum_title  = ''  ) : ## title for ``sum'' variable 
        """Convert pair of variables into ``sum'' & ``asymmetry'' pair
        - ``sum''       :  sum_scale  * ( var1 + var2  )
        - ``asymmetry'' :  asym_scale * ( var1 - var2 ) / (var1+ + var2 ) 
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
    ## convert a pair of variables ``half-sum''&``asymmetry'' into ``var1'', ``var2''
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
        """Convert a pair of variables ``half-sum''&``asymmetry'' into ``var1'', ``var2''
        >>> halfsum   = ...
        >>> asymmetry = ...
        >>> var1 , var2 = self.vars_from_asymmetry ( halfsum , asymmetry )
        """

        if hsumvar is None :
            return hsumvar , hsumvar
        
        if isinstance ( hsumvar , ROOT.RooAbsArg ) and isinstance ( hsumvar , ROOT.RooAbsArg )  :
            s = hsumvare.name
            a = asymvare.name
            if not v1name  : v1name  = '%sL' % s
            if not v2name  : v2name = '%sR' % s            
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
    ## Helper function to prepare ``soft'' Gaussian constraint for the variable
    #  @attention the constraint is prepared, but not applied!
    #  @code
    #  sigma      = ...
    #  constraint = pdf.soft_constraint( sigma , VE ( 0.15 , 0.01**2 ) )
    #  @endcode 
    def soft_constraint ( self , var , value , name = '' , title = '' ) :
        """Prepare ``soft'' Gaussian constraint for the variable
        -  consraint is prepared but not applied!
        >>> sigma      = ...
        >>> constraint = pdf.make_constraint( sigma , VE ( 0.15 , 0.01**2 ) )
        """
        
        assert isinstance ( var   , ROOT.RooAbsReal ) ,\
               "Invalid ``v'': %s/%s"  % ( var , type ( var ) )               
        assert isinstance ( value , VE ),\
               "Invalid ``value'': %s/%s"  % ( value , type ( value ) )

        assert 0 < value.cov2() , 'Invalid error for %s' % value
        
        name  = name  if name  else 'Gauss_%s_%s'                      % ( var.GetName() , self.name ) 
        title = title if title else 'Gaussian Constraint(%s,%s) at %s' % ( var.GetName() , self.name , value )
        
        # value & error as RooFit objects: 
        val = ROOT.RooFit.RooConst ( value.value () )
        err = ROOT.RooFit.RooConst ( value.error () )
        
        # Gaussian constrains 
        gauss = ROOT.RooGaussian ( self.var_name ( name ) , title , var , val , err )
        
        # keep all the created technical stuff  
        self.aux_keep.append ( val   )
        self.aux_keep.append ( err   )
        self.aux_keep.append ( gauss )

        self.info ('Constraint is created %s=%s' % ( var.name , value ) )
        return  gauss 

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
               "Invalid ``value'': %s/%s"  % ( value , type ( value ) )
        
        if 1 < abs ( value.value() )  :
            return self.soft_ratio_constraint ( b , a , 1.0 / value )

        fa = isinstance ( a , ROOT.RooAbsReal )
        fb = isinstance ( b , ROOT.RooAbsReal )

        if   fa and fv :
            
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
               "Invalid ``value'': %s/%s"  % ( value , type ( value ) )

        fa = isinstance ( a , ROOT.RooAbsReal )
        fb = isinstance ( b , ROOT.RooAbsReal )

        if   fa and fv :
            
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
               "Invalid ``value'': %s/%s"  % ( value , type ( value ) )

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

        if   fa and fv :
            
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
               "Invalid ``value'': %s/%s"  % ( value , type ( value ) )
        
        fa = isinstance ( a , ROOT.RooAbsReal )
        fb = isinstance ( b , ROOT.RooAbsReal )

        if   fa and fv :
            
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
        """Create ready-to-use ``soft'' gaussian constraint for the variable
        
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
    ## sample ``random'' positive number of events
    #  @code
    #  n =  pdf.gen_sample ( 10            ) ## get poissonian 
    #  n =  pdf.gen_sample ( VE ( 10 , 3 ) ) ## get gaussian stuff
    #  @endcode
    def gen_sample ( self , nevents ) :
        """Sample ``random'' positive number of events
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

        assert is_good_number ( vmin ), 'Invalid type of ``min'' %s/%s' % ( vmin , type ( vmin ) )
        assert is_good_number ( vmax ), 'Invalid type of ``max'' %s/%s' % ( vmin , type ( vmin ) )
        assert vmin < vmax, 'Invalid min/max range: %s/%s' % ( vmin , vmax )
        
        return vmin , vmax

  
# =============================================================================


# =============================================================================
## simple convertor of 1D-histo to data set
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
class H1D_dset(MakeVar) :
    """Simple convertor of 1D-histogram into data set
    >>> h   = ...
    >>> dset = H1D_dset ( h )
    One can create `binned` (default) or `weighted` data set
    >>> wset = H1D_dset ( h , weighted = True  )
    >>> bset = H1D_dset ( h , weighted = False ) ## default 
    """
    
    w_min = -1.e+100 
    w_max =  1.e+100
    
    def __init__ ( self             , 
                   histo            ,
                   xaxis    = None  , ## predefined axis/variable ? 
                   density  = False ,
                   weighted = False , ## weighted or binned? 
                   silent   = False ) :
        
        import ostap.histos.histos
        
        #
        ## use mass-variable
        #
        assert isinstance ( histo , ROOT.TH1 ) and 1 == histo.dim () , "``histo'' is not ROOT.TH1"
        self.__histo      = histo 
        self.__histo_hash = hash ( histo )
        

        name           = histo.GetName()
        self.__xaxis   = self.make_var ( xaxis , 'x_%s' % name , 'x-axis(%s)' % name , None , *(histo.xminmax()) )
        
        self.__density = True if density else False 
        self.__silent  = silent 

        self.__wvar    = None

        with roo_silent ( self.silent ) :  

            ## create weighted dataset ?
            if weighted :
                
                wname = weighted if isinstance ( weighted , string_types ) else 'h1weight'
                
                self.__wvar = ROOT.RooRealVar  ( wname , "weight-variable" , 1 , self.w_min , self.w_max ) 
                self.__vset = ROOT.RooArgSet   ( self.__xaxis ,  self.__wvar )
                self.__wset = ROOT.RooArgSet   ( self.__wvar      )
                self.__warg = ROOT.RooFit.WeightVar ( self.__wvar ) , ROOT.RooFit.StoreError ( self.__wset )
                self.__dset = ROOT.RooDataSet  (
                    rootID ( 'whds_' )  , "Weighted data set for the histogram '%s'" % histo.GetTitle() ,
                    self.__vset , *self.__warg )

                xvar = self.__xaxis
                wvar = self.__wvar 
                with SETVAR ( xvar ) :
                    for i, x , v in histo.items () :
                        xvar.setVal     ( x.value () )
                        self.__dset.add ( self.__vset , v.value() , v.error() ) 
                        
            ## create binned dataset 
            else :
                
                self.__vlst = ROOT.RooArgList    ( self.xaxis )
                self.__vimp = ROOT.RooFit.Import ( self.histo , self.density )
                self.__dset = ROOT.RooDataHist   (
                    rootID ( 'bhds_' ) , "Binned data set for histogram '%s'" % histo.GetTitle() ,
                    self.__vlst  ,
                    self.__vimp  )
                
    @property     
    def xaxis ( self ) :
        """The histogram x-axis variable"""
        return self.__xaxis
    @property
    def histo ( self ) :
        """The  histogram itself"""
        return self.__histo
    @property
    def density( self ) :
        """Treat the histo as ``density'' histogram?"""
        return self.__density    
    @property
    def silent( self ) :
        """Use the silent mode?"""
        return self.__silent
    @property
    def dset ( self ) :
        """``dset'' : ROOT.RooDataHist object"""
        return self.__dset
    @property
    def histo_hash ( self ) :
        """Hash value for the histogram"""
        return self.__histo_hash
    @property
    def weight ( self ) :
        """``weight'' : get weight variable if defined, None otherwise"""
        return self.__wvar

# =============================================================================
## simple convertor of 2D-histo to data set
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H2D_dset(MakeVar) :
    """Simple convertor of 2D-histogram into data set
    """
    
    w_min = -1.e+100 
    w_max =  1.e+100
    
    def __init__ ( self             ,
                   histo            ,
                   xaxis    = None  ,
                   yaxis    = None  ,
                   density  = False ,
                   weighted = False ,
                   silent   = False ) :
        #
        import ostap.histos.histos
        
        assert isinstance ( histo , ROOT.TH2 ) and 2 == histo.dim() , "``histo'' is not ROOT.TH2"
        self.__histo      =        histo
        self.__histo_hash = hash ( histo )

        ## use mass-variable
        #
        name                 = histo.GetName()
        if not xaxis : xaxis = histo.xminmax() 
        if not yaxis : yaxis = histo.yminmax() 
        
        self.__xaxis = self.make_var ( xaxis , 'x_%s' % name , 'x-axis(%s)' % name ,
                                       xaxis , *(histo.xminmax()) )
        self.__yaxis = self.make_var ( yaxis , 'y_%s' % name , 'y-axis(%s)' % name ,
                                       yaxis , *(histo.yminmax()) )
        
        self.__density = True if density else False 
        self.__silent  = silent
        
        self.__wvar    = None

        with roo_silent ( silent ) : 

            ## create weighted dataset 
            if weighted :
                
                wname = weighted if isinstance ( weighted , string_types ) else 'h2weight'
                
                self.__wvar = ROOT.RooRealVar  ( wname , "weight-variable" , 1 , self.w_min , self.w_max )
                self.__vset = ROOT.RooArgSet   ( self.__xaxis ,  self.__yaxis , self.__wvar )
                self.__warg = ROOT.RooFit.WeightVar ( self.__wvar )
                self.__dset = ROOT.RooDataSet (
                    rootID ( 'whds_' )  , "Weighted data set for the histogram '%s'" % histo.GetTitle() ,
                    self.__vset ,
                    self.__warg )

                xvar = self.__xaxis
                yvar = self.__yaxis
                wvar = self.__wvar 
                with SETVAR ( xvar ) :
                    for i, x , y , v in histo.items () :
                        xvar.setVal     ( x.value () )
                        yvar.setVal     ( y.value () )
                        self.__dset.add ( self.__vset , v.value() , v.error() ) 

            ## create binned dataset 
            else : 
                self.__vlst  = ROOT.RooArgList    ( self.__xaxis , self.__yaxis )
                self.__vimp  = ROOT.RooFit.Import ( histo , density )
                self.__dset  = ROOT.RooDataHist   (
                    rootID ( 'bhds_' ) , "Binned sata set for histogram '%s'" % histo.GetTitle() ,
                    self.__vlst  ,
                    self.__vimp  )
            
    @property     
    def xaxis  ( self ) :
        """The histogram x-axis variable"""
        return self.__xaxis
    @property     
    def yaxis  ( self ) :
        """The histogram y-axis variable"""
        return self.__yaxis
    @property
    def histo ( self ) :
        """The  histogram itself"""
        return self.__histo
    @property
    def density( self ) :
        """Treat the histo as ``density'' histogram?"""
        return self.__density    
    @property
    def silent( self ) :
        """Use the silent mode?"""
        return self.__silent
    @property
    def dset ( self ) :
        """``dset'' : ROOT.RooDataHist object"""
        return self.__dset
    @property
    def histo_hash ( self ) :
        """Hash value for the histogram"""
        return self.__histo_hash
    @property
    def weight ( self ) :
        """``weight'' : get weight variable if defined, None otherwise"""
        return self.__wvar
    
# =============================================================================
## simple convertor of 3D-histo to data set
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H3D_dset(MakeVar) :
    """Simple convertor of 3D-histogram into data set
    """
    
    w_min = -1.e+100 
    w_max =  1.e+100
    
    def __init__ ( self             ,
                   histo            ,
                   xaxis    = None  ,
                   yaxis    = None  ,
                   zaxis    = None  ,
                   density  = False ,
                   weighted = False ,
                   silent   = False ) :
        
        import ostap.histos.histos

        assert isinstance ( histo , ROOT.TH3 ) and 3 == histo.dim () , "``histo'' is not ROOT.TH3"
        self.__histo      =        histo
        self.__histo_hash = hash ( histo )
        #
        ## use mass-variable
        #
        name                 = histo.GetName() 

        if not xaxis : xaxis = histo.xminmax()
        if not yaxis : yaxis = histo.yminmax()
        if not zaxis : zaxis = histo.zminmax()
        
        self.__xaxis = self.make_var ( xaxis , 'x_%s' % name , 'x-axis(%s)' % name ,
                                       xaxis , *(histo.xminmax()) )
        self.__yaxis = self.make_var ( yaxis , 'y_%s' % name , 'y-axis(%s)' % name ,
                                       yaxis , *(histo.yminmax()) )
        self.__zaxis = self.make_var ( zaxis , 'z_%s' % name , 'z-axis(%s)' % name ,
                                       zaxis , *(histo.zminmax()) )
        
        self.__density = True if density else False 
        self.__silent  = silent
        
        self.__wvar    = None
        
        with roo_silent ( silent ) : 

            ## create weighted dataset 
            if weighted :
                
                wname = weighted if isinstance ( weighted , string_types ) else 'h2weight'
                
                self.__wvar = ROOT.RooRealVar  ( wname , "weight-variable" , 1 , self.w_min , self.w_max )
                self.__vset = ROOT.RooArgSet   ( self.__xaxis ,  self.__yaxis , self.__zaxis, self.__wvar )
                self.__warg = ROOT.RooFit.WeightVar ( self.__wvar )
                self.__dset = ROOT.RooDataSet  (
                    rootID ( 'whds_' )  , "Weighted data set for the histogram '%s'" % histo.GetTitle() ,
                    self.__vset ,
                    self.__warg )

                xvar = self.__xaxis
                yvar = self.__yaxis
                zvar = self.__zaxis
                wvar = self.__wvar 
                with SETVAR ( xvar ) :
                    for i, x , y , z , v in histo.items () :
                        xvar.setVal     ( x.value () )
                        yvar.setVal     ( y.value () )
                        zvar.setVal     ( z.value () )
                        self.__dset.add ( self.__vset , v.value() , v.error() ) 

            ## create binned dataset 
            else : 
                                
                self.__vlst  = ROOT.RooArgList    ( self.__xaxis , self.__yaxis , self.__zaxis )
                self.__vimp  = ROOT.RooFit.Import ( histo , density )
                self.__dset  = ROOT.RooDataHist   (
                    rootID ( 'bhds_' ) , "Binned data set for histogram '%s'" % histo.GetTitle() ,
                    self.__vlst  ,
                    self.__vimp  )
            
    @property     
    def xaxis  ( self ) :
        """The histogram x-axis variable"""
        return self.__xaxis
    @property     
    def yaxis  ( self ) :
        """The histogram y-axis variable"""
        return self.__yaxis    
    @property     
    def zaxis  ( self ) :
        """The histogram z-axis variable"""
        return self.__zaxis
    @property
    def histo ( self ) :
        """The  histogram itself"""
        return self.__histo
    @property
    def density( self ) :
        """Treat the histo as ``density'' histogram?"""
        return self.__density    
    @property
    def silent( self ) :
        """Use the silent mode?"""
        return self.__silent
    @property
    def dset ( self ) :
        """``dset'' : ROOT.RooDataHist object"""
        return self.__dset
    @property
    def histo_hash ( self ) :
        """Hash value for the histogram"""
        return self.__histo_hash
    @property
    def weight ( self ) :
        """``weight'' : get weight variable if defined, None otherwise"""
        return self.__wvar
    
# =============================================================================
## @class XVar
#  Helper class to keep all properties of the x-variable
class XVar(MakeVar) :
    """Helper class to keep all properteis the x-variable
    """
    def __new__( cls, *args, **kwargs):
        if  python_version.major > 2 : obj = super(XVar, cls).__new__( cls )
        else                         : obj = super(XVar, cls).__new__( cls , *args , **kwargs )
        ##
        obj.__xvar_init = False  
        return obj

    def  __init__ ( self , xvar ):

        if self.__xvar_init : return
        else                : self.__xvar_init = False
        
        self.__xvar = None

        if   isinstance ( xvar , ROOT.TH1   ) : xvar = xvar.xminmax()
        elif isinstance ( xvar , ROOT.TAxis ) : xvar = xvar.GetXmin() , xvar.GetXmax()
        
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
            raise TypeError( "``x-variable''is not specified properly %s/%s" % ( xvar , type ( xvar ) ) )
                
    @property 
    def xvar ( self ) :
        """``x''-variable  (same as ``x'')"""
        return self.__xvar
    @property 
    def x    ( self ) :
        """``x''-variable (same as ``xvar'')"""
        return self.xvar
    def xminmax ( self ) :
        """Min/max values for x-variable"""
        return self.xvar.minmax()

    ## get the proper xmin/xmax range
    def xmnmx    ( self , xmin , xmax ) :
        """Get the proper xmin/xmax range
        """
        return self.vmnmx ( self.xvar , xmin , xmax )

# =============================================================================
## @class YVar
#  Helper class to keep all properties of the y-variable
class YVar(MakeVar) :
    """Helper class to keep all properteis the y-variable
    """
    def __new__( cls, *args, **kwargs):
        if  python_version.major > 2 : obj = super(YVar, cls).__new__( cls )
        else                         : obj = super(YVar, cls).__new__( cls , *args , **kwargs )
        ##
        obj.__yvar_init = False  
        return obj

    def  __init__ ( self , yvar ):
         
        if self.__yvar_init : return
        else                : self.__yvar_init = False

        self.__yvar = None

        if   isinstance ( yvar , ROOT.TH1   ) : yvar = yvar.xminmax()
        elif isinstance ( yvar , ROOT.TAxis ) : yvar = yvar.GetXmin() , yvar.GetXmax()
        
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
            raise TypeError( "``y-variable''is not specified properly %s/%s" % ( yvar , type ( yvar ) ) )
                
    @property 
    def yvar ( self ) :
        """``y''-variable (same as ``y'')"""
        return self.__yvar
    @property 
    def y    ( self ) :
        """``y''-variable (same as ``yvar'')"""
        return self.yvar
    def yminmax ( self ) :
        """Min/max values for y-variable"""
        return self.yvar.minmax()

    ## get the proper ymin/ymax range
    def ymnmx    ( self , ymin , ymax ) :
        """Get the proper ymin/ymax range
        """
        return self.vmnmx ( self.yvar , ymin , ymax )


# =============================================================================
## @class ZVar
#  Helper class to keep all properties of the z-variable
class ZVar(MakeVar) :
    """Helper class to keep all properteis the z-variable
    """
    def __new__( cls, *args, **kwargs):
        if  python_version.major > 2 : obj = super(ZVar, cls).__new__( cls )
        else                         : obj = super(ZVar, cls).__new__( cls , *args , **kwargs )
        ##
        obj.__zvar_init = False  
        return obj
        
    def  __init__ ( self , zvar ):
         
        if self.__zvar_init : return
        else                : self.__zvar_init = False

        self.__zvar = None

        if   isinstance ( zvar , ROOT.TH1   ) : zvar = zvar.xminmax()
        elif isinstance ( zvar , ROOT.TAxis ) : zvar = zvar.GetXmin() , yvar.GetXmax()
        
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
            raise TypeError( "``z-variable'' is not specified properly %s/%s" % ( zvar , type ( zvar ) ) )
                
    @property 
    def zvar ( self ) :
        """``y''-variable (same as ``z'')"""
        return self.__zvar
    @property 
    def z    ( self ) :
        """``z''-variable (same as ``zvar'')"""
        return self.zvar
    def zminmax ( self ) :
        """Min/max values for y-variable"""
        return self.zvar.minmax()

    ## get the proper xmin/xmax range
    def zmnmx    ( self , xmin , xmax ) :
        """Get the proper zmin/zmax range
        """
        return self.vmnmx ( self.zvar , zmin , zmax )

# =============================================================================
## @class Phases
#  helper class to build/keep the list of ``phi''-arguments
#   - needed e.g. for polynomial functions
class Phases(MakeVar) :
    """Helper class to build/keep the list of ``phi''-arguments,
    (needed e.g. for polynomial functions)
    """
    ## Create vector of phases (needed for various polynomial forms)
    def __init__( self  , power , the_phis = None ) :
        """Create vector of phases (needed for various polynomial forms)
        """

        ## check  the arguments 
        assert isinstance ( power , num_types ) and  int ( power ) == power and 0 <= power, \
               "Phases: invalid type/value for ``power''-parameter: %s/%s"  % (  power , type(power) )
        power = int ( power ) 

        if  isinstance ( the_phis , Phases ) : 
            self.__phis     = [ i for i in the_phis.phis ]  
            self.__phi_list = the_phis.phi_list            
            assert power == len( self.__phis ) , "Phases: Invalid length of ``phis''  %d/%s" %  ( power , len ( self.__phis ) )            
            return                                                   ## RETURN
        elif the_phis and isinstance ( the_phis , ROOT.RooArgList ) :
            self.__phis     = [ i for i in the_phis]  
            self.__phi_list = the_phis 
            assert power == len( self.__phis ) , "Phases: Invalid length of ``phis''  %d/%s" %  ( power , len ( self.__phis ) )            
            return                                                   ##  RETURN      
        elif the_phis and isinstance ( the_phis , (tuple,list) ) :
            self.__phis     = [ i for i in the_phis]  
            self.__phi_list = ROOT.RooArgList()
            for phi in the_phis : self.__phi_list.add ( phi )
            assert power == len( self.__phis ) , "Phases: Invalid length of ``phis''  %d/%s" %  ( power , len ( self.__phis ) )            
            return                                                   ## RETURN
        elif the_phis :
            self.warning("unknown type for ``the_phis'' %s/%s, skip it" % ( the_phis , type(the_phis) ) )

        self.__phis     = []
        self.__phi_list = ROOT.RooArgList()
        from math import pi
        for i in range( 0 , power ) :
            phi_i = self.make_var ( None ,
                                    'phi%d_%s'      % ( i , self.name )  ,
                                    '#phi_{%d}(%s)' % ( i , self.name )  ,
                                    None , 0 ,  -1.55 * pi  , 3.55 * pi  )
            self.__phis    .append ( phi_i ) 
            self.__phi_list.add    ( phi_i )
        
    ## set all phis to be 0
    def reset_phis ( self ) :
        """Set all phases to be zero
        >>> pdf = ...
        >>> pdf.reset_phis() 
        """
        for f in self.__phis : f.setVal(0)
        
    @property
    def phis ( self ) :
        """The list/tuple of ``phases'', used to parameterize various polynomial-like shapes
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
        """The list/ROOT.RooArgList of ``phases'', used to parameterize polynomial-like shapes
        """
        return self.__phi_list
    @property
    def phis_lst ( self ) :
        """The list/ROOT.RooArgList of ``phases'', used to parameterize polynomial-like shapes
        """
        return self.__phi_list

    ## ## number of parameters 
    ## def __len__     ( self ) :
    ##     """get the number of phi-parameters"""
    ##     return  len ( self.__phis )
    
    ## ## Get certain phi value(s)
    ## def __getitem__ ( self , index ) :
    ##     """Get certain phi-value by index
    ##     >>> obj = ...
    ##     >>> print(obj[3])
    ##     >>> print(obj[3:9])
    ##     """
    ##     return self.__phis [ index ]
    
    ## ## Change certain phi-value(s)
    ## #  @code
    ## #  >>> obj = ...
    ## #  >>> obj[3]   = 15
    ## #  >>> obj[1:3] = 2, 0   
    ## #  @endcode 
    ## def __setitem__ ( self , index , values ) :
    ##     """Change certain phi values
    ##     >>> obj = ...
    ##     >>> obj[3]   = 15
    ##     >>> obj[1:3] = 2, 0   
    ##     """
    ##     my_phis = self.__phis [ index ]
    ##     if isinstance ( my_phis , ROOT.RooAbsReal ) :
    ##         my_phis.setVal ( float ( values ) )
    ##     else  :            
    ##         for p , v in zip (  my_phis , values ) : p.setVal ( float ( v ) )
        

# =============================================================================
## @class ParamsPoly
#  Helper base class to implement polynomials 
class ParamsPoly(MakeVar) :
    """Helper base class to implement polynomials 
    """
    def __init__ ( self , power = 1 , pars = None ) :
        
        ## initialize the base class 
        MakeVar.__init__ ( self )
        
        assert pars or ( isinstance ( power , integer_types ) and 0 <= power ) ,\
               'Inconsistent power/npars setting'

        self.__pars     = [] 
        self.__pars_lst = ROOT.RooArgList()
        
        params = []
        if pars :
            for i , p in enumerate ( pars ) :                
                pp = self.make_var ( p ,
                                     'par%d_%s'            % ( i , self.name ) ,
                                     'parameter %d for %s' % ( i , self.name ) , None , p )
                if not isinstance ( p , ROOT.RooAbsReal ) : pp.setConstant ( False ) 
                params.append ( pp )
        else :
            for i in range ( power + 1 ) :                
                pp = self.make_var ( 0.0 ,
                                     'par%d_%s'            % ( i , self.name ) ,
                                     'parameter %d for %s' % ( i , self.name ) , None , 0.0 )
                pp.setConstant ( False ) 
                params.append ( pp )

        for p  in params :
            self.__pars.append  ( p )
            self.__pars_lst.add ( p )
            
        self.__pars = tuple ( self.__pars ) 

        assert self.pars , 'Invalid number of parameters!'

        self.config = {
            'name' : self.name ,
            'xvar' : self.xvar ,
            'pars' : self.pars }

    @property
    def pars  ( self ) :
        """``pars'' : the polynomial coefficients/parameters"""
        return self.__pars    
    @pars.setter
    def pars ( self , values ) :
        self.component_setter ( self.__pars , values )
            
    def reset_pars ( self , value = 0 ) :
        """Set all pars to be value 
        >>> pdf = ...
        >>> pdf.reset_pars() 
        """
        for f in self.__pars : f.setVal( value )

    @property
    def pars_lst ( self ) :
        """``pars_lst'' : the polynomial coefficients/parameters as RooArgList"""
        return self.__pars_lst
    
    @property
    def power ( self ) :
        """``power''  : polynomial degree """
        return len  ( self.pars ) - 1

# =============================================================================
## @class ShiftScalePoly
#  Helper base class to implement polynomials
#  \f$ f(x) = a + b P(x) \f$,
#  where \f$P(x)\f$ some special polynomial 
class ShiftScalePoly ( Phases ) :
    """Helper fbase class to implemnet polynomials 
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
                                    'bias/shift for %s' % self.name , a )
        ## parameter a 
        self.__b  = self.make_var ( b ,
                                    'b_%s'              % self.name ,
                                    'scale for %s'      % self.name , b )
        
        self.config = {
            'name' : self.name ,
            'xvar' : self.xvar ,
            'pars' : self.pars ,
            'a'    : self.a    ,
            'b'    : self.b    }

    @property
    def a ( self ) :
        """``a'' : bias parameter for polynomial:  f(x) = a + b*M(x)"""
        return self.__a
    @a.setter
    def a ( self , value ) :
        vv = float ( value )
        if self.__a.minmax () and not vv in self.__a  :
            self.error ("Value %s is outside the allowed region %s for %s"  % ( vv , self.__a.minmax() , self.__a.name ) )
        self.__a.setVal ( vv )

    @property
    def b ( self ) :
        """``scale'' : bias parameter for polynomial:  f(x) = a + b*M(x)"""
        return self.__b
    @b.setter
    def b ( self , value ) :
        vv = float ( value )
        if self.__b.minmax () and not vv in self.__b  :
            self.error ("Value %s is outside the allowed region %s for %s"  % ( vv , self.__b.minmax() . self.__b.name ) )
        self.__b.setVal ( vv )
    
    @property
    def shift  ( self ) :
        """``shift'' : bias/shift parameter  (same as ``a'')"""
        return self.__a
    @shift.setter
    def shift  ( self , value ) : self.a = value 
    @property
    def bias  ( self ) :
        """``bias'' : bias/shift parameter  (same as ``a'')"""
        return self.__a
    @bias.setter
    def bias   ( self , value ) : self.a = value
    
    @property
    def scale ( self ) :
        """``scale'' : bias/shift parameter  (same as ``b'')"""
        return self.__b
    @scale.setter
    def scale ( self , value ) : self.b = value

    @property
    def pars  ( self ) :
        """``pars'' :  polynomial parameters (same as ``phis'')"""
        return self.phis
    @pars.setter 
    def pars  ( self , values ) :
        self.phis = values 
    @property
    def pars_lst ( self ) :
        """``pars_lst'' :  polynomial parameters as RooArgList"""
        return self.phis_lst
    

# ==============================================================================
## Should one use ``similar'' component?
def component_similar ( same ) :
    """Should one use ``similar'' component?
    """
    if   same is Ellipsis           : return True
    elif same is NotImplemented     : return True
    elif isinstance ( same , str  ) \
         and same.strip().lower() in ( 'ditto' , 'similar' ) : return True
    return False

# =============================================================================
## Should one      ``clone''  component?
def component_clone  ( same ) :
    """Should one use ``cloned'' component?
    """    
    if isinstance ( same , str ) \
       and same.strip().lower() in ( 'clone' , 'cloned' , 'same' ) : return True        
    return False 
# =============================================================================

    
# =============================================================================
##  get <code>i</code>-th component from <code>what</code>
def get_i ( what , i , default = None ) :
    """
    """
    if   isinstance ( what , ROOT.RooArgList ) and i in what             : return what[i]
    elif isinstance ( what , ROOT.RooAbsReal ) and 0 == i                : return what 
    elif isinstance ( what , num_types       ) and 0 == i                : return what  
    elif isinstance ( what , list_types      ) and 0 <= i < len ( what ) : return what[i] 

    return default
        
        
# =============================================================================
## consruct MsgTopic
#  @see RooFit::MsgTopic
#  @code
#  topic = msgTopic ( ROOT.RooFit.Fitting ) 
#  topic = msgTopic ( ROOT.RooFit.Fitting , ROOT.RooFit.Caching )
#  topic = msgTopic ( 'Fitting' , 'Caching' )
#  @endcode
def msg_topic ( *topics ) :
    """onsruct MsgTopic
    >>> topic = msgTopic ( ROOT.RooFit.Fitting ) 
    >>> topic = msgTopic ( ROOT.RooFit.Fitting , ROOT.RooFit.Caching )
    >>> topic = msgTopic ( 'Fitting' , 'Caching' )
    """
    topic = 0
    for i in  topics : 
        if   isinstance ( i , integer_types )  : topic |= i
        elif isinstance ( i , string_types  )  :
            ii = i.lower() 
            if   ii == 'generation'            : topic |=  ROOT.RooFit.Generation 
            elif ii == 'minimization'          : topic |=  ROOT.RooFit.Minimization
            elif ii == 'minization'            : topic |=  ROOT.RooFit.Minimization
            elif ii == 'plotting'              : topic |=  ROOT.RooFit.Plotting
            elif ii == 'fitting'               : topic |=  ROOT.RooFit.Fitting 
            elif ii == 'integration'           : topic |=  ROOT.RooFit.Integration 
            elif ii == 'linkstatemgmt'         : topic |=  ROOT.RooFit.LinkStateMgmt
            elif ii == 'eval'                  : topic |=  ROOT.RooFit.Eval
            elif ii == 'caching'               : topic |=  ROOT.RooFit.Caching
            elif ii == 'optimization'          : topic |=  ROOT.RooFit.Optimization
            elif ii == 'optimisation'          : topic |=  ROOT.RooFit.Optimization
            elif ii == 'objecthandling'        : topic |=  ROOT.RooFit.ObjectHandling
            elif ii == 'inputarguments'        : topic |=  ROOT.RooFit.InputArguments
            elif ii == 'tracing'               : topic |=  ROOT.RooFit.Tracing
            elif ii == 'contents'              : topic |=  ROOT.RooFit.Contents
            elif ii == 'datahandling'          : topic |=  ROOT.RooFit.DataHandling
            elif ii == 'numintegration'        : topic |=  ROOT.RooFit.NumIntegration
            elif ii == 'numericintegration'    : topic |=  ROOT.RooFit.NumIntegration
            elif ii == 'numericalintegration'  : topic |=  ROOT.RooFit.NumIntegration
            elif ii == 'fastevaluations'       : topic |=  ROOT.RooFit.FastEvaluations
            else : logger.error ( 'MsgTopic/1: unknown topic %s, skip' % i )
        else : logger.error ( 'MsgTopic/2: unknown topic %s/%s, skip' % ( i , type ( i ) ) )
        
    return topic 
    
# =============================================================================
# upgraded constructor for class Ostap::Utils::RemoveTopicsd
# @code
# with RemoveTopic ( [ 'Fitting' , 'Plotting' ] ) :
#    ... do something ...
# with RemoveTopic ( ROOT.RooFit.Plotting | ROOT.RooFit.Fitting ) :
#    ... do something ...
# @endcode
# @see Ostap::Utils::AddTopic
# @see Ostap::Utils::RemoveTopic
def _rt_new_init_ ( self , topics , level = ROOT.RooFit.INFO , streams = -1  ) :
    """ Upgraded constructor for class Ostap::Utils::RemoveTopics
    >>> with RemoveTopic ( [ 'Fitting' , 'Plotting' ] ) :
    ...    ... do something ...
    >>> with RemoveTopic ( ROOT.RooFit.Plotting | ROOT.RooFit.Fitting ) :
    ...    ... do something ...
    - see Ostap::Utils::AddTopic
    - see Ostap::Utils::RemoveTopic
    """

    if isinstance ( topics , integer_types ) and 0 < topics and topics <= 2**16 :
        return self._old_init_ ( topics , level , streams )    
    if isinstance ( topics , string_types  ) : topics = topics.split()
    topic = msg_topic ( *topics )
    return self._old_init_ ( topic , level , streams )

if not hasattr ( Ostap.Utils.RemoveTopic , '_old_init_' ) :
    Ostap.Utils.RemoveTopic._old_init_ = Ostap.Utils.RemoveTopic.__init__
    Ostap.Utils.RemoveTopic.__init__   = _rt_new_init_
    Ostap.Utils.RemoveTopic.__enter__  = lambda s : s
    Ostap.Utils.RemoveTopic.__exit__   = lambda s,*_ : s.exit() 
    
# =============================================================================
# upgraded constructor for class Ostap::Utils::AddTopic
# @code
# with AddTopic ( [ 'Fitting' , 'Plotting' ] ) :
#    ... do something ...
# with AddTopic ( ROOT.RooFit.Plotting | ROOT.RooFit.Fitting ) :
#    ... do something ...
# @endcode
# @see Ostap::Utils::AddTopic
# @see Ostap::Utils::RemoveTopic
def _at_new_init_ ( self , topics , streams = -1  ) :
    """ Upgraded constructor for class Ostap::Utils::AddTopics
    >>> with RemoveTopic ( [ 'Fitting' , 'Plotting' ] ) :
    ...    ... do something ...
    >>> with RemoveTopic ( ROOT.RooFit.Plotting | ROOT.RooFit.Fitting ) :
    ...    ... do something ...
    - see Ostap::Utils::AddTopic
    - see Ostap::Utils::RemoveTopic
    """

    if isinstance ( topics , integer_types ) and 0 < topics and topics <= 2**16 :
        return self._old_init_ ( topics , streams )
    
    if isinstance ( topics , string_types  ) : topics = [ topics ]

    topic = msg_topic ( *topics )
    return self._old_init_ ( topic , streams )
                
if not hasattr (  Ostap.Utils.AddTopic , '_old_init_' ) :
    Ostap.Utils.AddTopic._old_init_ = Ostap.Utils.AddTopic.__init__
    Ostap.Utils.AddTopic.__init__   = _at_new_init_
    Ostap.Utils.AddTopic.__enter__  = lambda s : s
    Ostap.Utils.AddTopic.__exit__   = lambda s,*_ : s.exit() 
    


# ================================================================================
## remove topic from Roofit message streams
#  @see RooMsgService
#  @code
#  with remove_topic ( ROOT.RooFit.Fitting ) :
#    ...
#  with remove_topic ( ROOT.RooFit.Fitting | ROOT.RooFit.Plotting ) :
#    ...
#  with remove_topic ( [ 'Fitting' , 'Plotting' ] ) :
#    ...
#  @endcode
#  @see Ostap::Utils::RemoveTopic
#  @see Ostap::Utils::AddTopic
def remove_topic ( topics , level = ROOT.RooFit.INFO , stream  = -1 ) :
    """Remove topic from Roofit message streams
    - see RooMsgService
    >>> with remove_topic ( ROOT.RooFit.Fitting ) :
    ...  ...
    >>> with remove_topic ( ROOT.RooFit.Fitting | ROOT.RooFit.Plotting ) :
    ... ...
    >>> with remove_topic ( [ 'Fitting' , 'Plotting' ] ) :
    ... ...
    - see Ostap::Utils::RemoveTopic
    - see Ostap::Utils::AddTopic
    """
    return Ostap.Utils.RemoveTopic ( topics , level , stream ) 


# ================================================================================
## add topic from RooFit message streams
#  @see RooMsgService
#  @code
#  with add_topic ( ROOT.RooFit.Fitting ) :
#    ...
#  with add_topic ( ROOT.RooFit.Fitting | ROOT.RooFit.Plotting ) :
#    ...
#  with add_topic ( [ 'Fitting' , 'Plotting' ] ) :
#    ...
#  @endcode
#  @see Ostap::Utils::RemoveTopic
#  @see Ostap::Utils::AddTopic
def add_topic ( topics , stream  = -1 ) :
    """Add topic to RooFit message streams
    - see RooMsgService
    >>> with add_topic ( ROOT.RooFit.Fitting ) :
    ...  ...
    >>> with add_topic ( ROOT.RooFit.Fitting | ROOT.RooFit.Plotting ) :
    ... ...
    >>> with add_topic ( [ 'Fitting' , 'Plotting' ] ) :
    ... ...
    - see Ostap::Utils::RemoveTopic
    - see Ostap::Utils::AddTopic
    """
    return Ostap.Utils.AddTopic ( topics , level , stream ) 


# =============================================================================
## suppress certain message topics
#  @code
#  suppress_topics ( 'Fitting'  , 'Caching' ) 
#  @endcode 
def suppress_topics ( *topics ) :
    """suppress certain message topics
    >>> suppress_topics ( 'Fitting'  , 'Caching' ) 
    """
    if topics and 1 == len( topics ) :
        t = str ( topics [ 0 ] ).lower()
        if 'config' == t : return suppress_topics() 

    if not topics :
        newtopics = [] 
        import ostap.core.config as CONFIG
        if 'RooFit' in CONFIG.config :
            import string
            ws     = string.whitespace 
            node   = CONFIG.config [ 'RooFit' ]
            data   = node.get('RemoveTopics','(,)' )
            topics = tuple ( i.strip ( ws ) for i in data.split ( ',' ) if i.strip ( ws ) ) 
            
    if topics : 
        svc = ROOT.RooMsgService.instance()
        svc.saveState () 
        topic = msg_topic ( *topics ) 
        num   = svc.numStreams()
        for i in range ( num ) : ok = Ostap.Utils.remove_topic ( i , topic ) 

# =============================================================================
## and finally suppress exra RooFit topics! 
suppress_topics ()


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
#                                                                       The END 
# =============================================================================
