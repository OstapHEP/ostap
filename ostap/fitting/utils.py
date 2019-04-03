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
    #
    "H1D_dset"          , ## 1D-histogram to RooDataHist converter 
    "H2D_dset"          , ## 2D-histogram to RooDataHist converter 
    "H3D_dset"          , ## 3D-histogram to RooDataHist converter
    #
    'component_similar' , ## Should one use ``similar'' component?
    'component_clone'   , ## Should one use ``cloned'' component?
    # 
    'numcpu'            , ## number of CPUs
    'ncpu'              , ## fuction to builf ROOT.RooFit.NumCPU 
    )
# =============================================================================
import ROOT, math
import ostap.fitting.variables 
import ostap.fitting.roocollections
from   ostap.core.core     import rootID, VE, items_loop
from   ostap.core.types    import num_types, list_types, integer_types
from   ostap.logger.utils  import roo_silent
from   sys                 import version_info as python_version 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.utils' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
## MINUIT covariance matrix status:
# - status = -1 :  not available (inversion failed or Hesse failed)
# - status =  0 : available but not positive defined
# - status =  1 : covariance only approximate
# - status =  2 : full matrix but forced pos def
# - status =  3 : full accurate matrix
_cov_qual_ = {
    -1 :  '-1/not available (inversion failed or Hesse failed)' ,
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
    ## @attention ensure that important attributes are available even before __init__
    def __new__( cls, *args, **kwargs):
        if  python_version.major > 2 : obj = super(MakeVar, cls).__new__( cls )
        else                         : obj = super(MakeVar, cls).__new__( cls , *args , **kwargs )
            
        obj.__aux_keep = []                      ## ATTENTION!        
        obj.__name     = None                    ## ATTENTION!
        return obj
    
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
        """``aux_keep'' -  the list of objects to be kept by this PDF"""
        return self.__aux_keep
    @property
    def logger   ( self ) :
        """``logger'': get the local Logger object"""
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
        assert isinstance ( value , str ) , "``name'' must  be a string, %s/%s is given" % ( value , type(value) ) 
        self.__name = value

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
    def make_var ( self , var , name , comment , fix = None , *args ) :
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
            var = ROOT.RooRealVar ( name , comment , *var )

        # var = value 
        if isinstance   ( var , num_types ) :
            if   not    args       : var = ROOT.RooRealVar ( name , comment , var             )
            elif 2 == len ( args ) : var = ROOT.RooRealVar ( name , comment , var , *args     )
            elif 3 == len ( args ) : var = ROOT.RooRealVar ( name , comment , var , *args[1:] )
        
        ## create the variable from parameters 
        if not isinstance ( var , ROOT.RooAbsReal ) : 
            var = ROOT.RooRealVar ( name , comment , *args )
            self.aux_keep.append ( var ) ##  ATTENTION: store newly created variable
        
        ## fix it, if needed
        if   isinstance ( fix , bool       ) : pass 
        elif isinstance ( fix , num_types  ) :
                
            if hasattr ( var , 'getMin)') and fix < var.getMin() and hasattr ( var , 'setMin' ) :                                                                             
                self.warning ( "make_var: min-value for %s is redefined to be %s " % ( var.GetName () , fix ) )
                var.setMin ( fix )
            
            if hasattr ( var , 'getMax' ) and fix > var.getMax() and hasattr ( var , 'setMax' ) : 
                self.warning ( "make_var: max-value for %s is redefined to be %s " % ( var.GetName () , fix ) )
                var.setMax ( fix )
            
            if not var.isConstant () and hasattr ( var , 'fix'    ) : var.fix    ( fix )
            elif                         hasattr ( var , 'setVal' ) : var.setVal ( fix )

        return var

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
    ## technical function to parse arguments for <code>fitTo</code>  function
    def parse_args ( self ,  dataset = None , *args , **kwargs ) :
        """Technical function to parse arguments for fitTo function
        """
        _args = []
        for a in args :
            if not isinstance ( a , ROOT.RooCmdArg ) :
                self.warning ( 'parse_args: unknown argument type %s/%s, skip' % ( a , type ( a ) ) )
            else : _args.append ( a ) 

        from ostap.plotting.fit_draw import keys  as drawing_options
        
        for k , a in items_loop ( kwargs ) :
            
            klow = k.lower ()
            kup  = k.upper ()
            
            ## skip "drawing" options 
            if k.lower() in drawing_options                            : continue 
            if k.lower() in ( 'draw' , 'draw_option', 'draw_options' ) : continue 
            
            if   isinstance ( a , ROOT.RooCmdArg ) : _args.append ( a )
            
            elif kup in ( 'VERBOSE' ,        ) and isinstance ( a , bool ) :
                _args.append ( ROOT.RooFit.Verbose (     a ) ) 
            elif kup in ( 'SILENT'           ,
                           'SILENCE'         ) and isinstance ( a , bool ) :
                _args.append ( ROOT.RooFit.Verbose ( not a ) ) 
            elif kup in ( 'STRATEGY'         , 
                          'MINUITSTRATEGY'   ,
                          'MINUIT_STRATEGY'  ) and isinstance ( a , integer_types ) and 0 <= a <= 2 : 
                _args.append ( ROOT.RooFit.Strategy (    a ) ) 
            elif kup in ( 'PRINTLEVEL'       ,
                          'PRINT_LEVEL'      ,
                          'MINUITPRINT'      ,
                          'MINUIT_PRINT'     ,
                          'MINUITLEVEL'      ,
                          'MINUIT_LEVEL'     ) and isinstance ( a , integer_types ) and -1 <= a <= 3 :
                _args.append ( ROOT.RooFit.PrintLevel ( a ) ) 
            elif kup in ( 'PRINTEVALERRORS'  ,
                          'PRINT_EVAL_ERRORS',
                          'PRINTEVAL_ERRORS' ,
                          'PRINT_EVALERRORS' ,
                          'PRINTERRORS'      ,
                          'PRINT_ERRORS'     ,
                          'ERRORSPRINT'      ,
                          'ERRORS_PRINT'     ) and isinstance ( a , integer_types ) and -1 <= a :
                _args.append ( ROOT.RooFit.PrintEvalErrors ( a ) )                
            elif kup in ( 'TIMER'            ,
                          'TIMING'           ) and isinstance ( a , bool ) :
                _args.append ( ROOT.RooFit.Timer    ( a ) ) 
            elif kup in ( 'WARNING'          ,
                          'WARNINGS'         ) and isinstance ( a , bool ) :
                _args.append ( ROOT.RooFit.Warnings ( a ) ) 
            
            elif kup in ( 'WEIGHTED'         ,
                          'SUMW2'            ,
                          'SUMW2ERROR'       ) and isinstance ( a , bool ) :
                
                if   a and dataset and     dataset.isWeighted()           : pass 
                elif a and dataset and not dataset.isWeighted()           :
                    self.warning ('parse_args: SumW2-flag is True  for non-weighted dataset')
                elif       dataset and not dataset.isWeighted() and not a : pass 
                elif       dataset and     dataset.isWeighted() and not a :
                    self.warning ('parse_args: SumW2-flag is False for     weighted dataset')                    

                _args.append (  ROOT.RooFit.SumW2Error( a ) )
                    
            elif kup in ( 'EXTENDED' ,       ) and isinstance ( a , bool ) :
                _args.append   (  ROOT.RooFit.Extended ( a ) )                
            elif kup in ( 'NCPU'             ,
                          'NCPUS'            ,
                          'NUMCPU'           ,
                          'NUMCPUS'          ) and isinstance ( a , int ) and 1<= a : 
                _args.append   (  ROOT.RooFit.NumCPU( a  ) ) 
            elif kup in ( 'NCPU'             ,
                          'NCPUS'            ,
                          'NUMCPU'           ,
                          'NUMCPUS'          ) and \
                 isinstance ( a , list_types ) and 2 == len ( a )  and \
                 isinstance ( a[0] , integer_types ) and 1 <= a[1] and \
                 isinstance ( a[1] , integer_types ) and 0 <= a[1] <=3 :
                _args.append   (  ROOT.RooFit.NumCPU( a[0] ,  a[1] ) ) 
                
            elif kup in ( 'RANGE'            ,
                          'FITRANGE'         ,
                          'FIT_RANGE'        ,
                          'RANGES'           ,
                          'FITRANGES'        ,
                          'FIT_RANGES'       ) and isinstance ( a , string_types ) :
                _args.append   (  ROOT.RooFit.Range ( a ) )  
            elif kup in ( 'RANGE'            ,
                          'FITRANGE'         ,
                          'FIT_RANGE'        ) and isinstance ( a , list_types   ) \
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
                           'INITIAL_HESSE'   ,
                           'INIT_HESSE'      ,
                           'HESSE_INIT'      ,
                           'HESSE_INITIAL'   ) and isinstance ( a , bool ) :
                _args.append   (  ROOT.RooFit.InitialHesse ( a )  )
            elif kup in ( 'OPTIMIZE'         ,
                          'OPTIMISE'         ) and isinstance ( a , integer_types  ) :
                _args.append   (  ROOT.RooFit.Optimize     ( a )  )
            elif kup in ( 'MINOS'    ,       ) and isinstance ( a , bool           ) :
                _args.append   (  ROOT.RooFit.Minos        ( a )  )
            elif kup in ( 'MINOS'    ,       ) and isinstance ( a , ROOT.RooArgSet ) :
                _args.append   (  ROOT.RooFit.Minos        ( a )  )
            elif kup in ( 'SAVE'     ,       ) and isinstance ( a , bool           ) :
                _args.append   (  ROOT.RooFit.Save         ( a )  )
            elif kup in ( 'FITOPTIONS'       ,
                          'FITOPTION'        ,
                          'FIT_OPTIONS'      ,
                          'FIT_OPTION'       ) and isinstance ( a , string_types ) :
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

        keys       = [ a.GetName() for a in _args ]        
        if not 'NumCPU' in keys :
            if  dataset and not isinstance ( dataset , ROOT.RooDataHist ) :
                _args.append ( ncpu ( len ( dataset ) ) )
            else :
                nc = numcpu()
                if  1 < nc : _args.append ( ROOT.RooFit.NumCPU ( nc ) ) 

        keys = [ str ( a ) for a in _args ]
        keys.sort () 
        self.info ( 'parse_args: Parsed arguments %s' % keys )

        return tuple ( _args )
    
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
        
        name  = name  if name  else 'Gauss_%s_%s'                       % ( var.GetName() , self.name ) 
        title = title if title else 'Gauissian Constraint(%s,%s) at %s' % ( var.GetName() , self.name , value )
        
        # value & error as RooFit objects: 
        val = ROOT.RooFit.RooConst ( value.value () )
        err = ROOT.RooFit.RooConst ( value.error () )
        
        # Gaussian constrains 
        gauss = ROOT.RooGaussian ( name , title , var , val , err )
        
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
        
        if isinstance ( a , ROOT.RooAbsReal ) and isinstance ( a , ROOT.RooAbsReal ) :

            vlist   = ROOT.RooArgList    ( a , b )
            
            formula = '(%s)/(%s)'   % ( a.name , b.name )
            vname   = 'Ratio_%s_%s' % ( a.name , b.name )
            vtitle  = 'Ratio %s/%s' % ( a.name , b.name )
            
            var     = ROOT.RooFormulaVar ( vname , vtitle , formula , vlist )
            
            self.aux_keep.append ( vlist  )
            self.aux_keep.append ( var    )
            
            return self.soft_constraint ( var , value , name , title )

        elif isinstance ( a , ROOT.RooAbsReal ) and isinstance ( b , num_types + (VE,) ) :
            
            return self.soft_constraint ( a , value * b , name , title )
        
        elif isinstance ( b , ROOT.RooAbsReal ) and isinstance ( a , num_types + (VE,) ) :
            
            return self.soft_constraint ( b , a / v , name , title )
        
        raise TypeError('Unknown types a&b: %s/%s' % ( type ( a ) , type ( b ) ) )
    
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
        to the ratio of the variables 
        >>> N1 = ...
        >>> N2 = ...
        >>> rC = pdf.soft_sum_constraint ( N1 , N2 , VE (10,1**2) )
        """
        assert isinstance ( value , VE ) and 0 < value.cov2() ,\
               "Invalid ``value'': %s/%s"  % ( value , type ( value ) )
        
        if isinstance ( a , ROOT.RooAbsReal ) and isinstance ( a , ROOT.RooAbsReal ) :

            vlist   = ROOT.RooArgList    ( a , b )
            
            formula = '(%s)+(%s)' % ( a.name , b.name )
            vname   = 'Sum_%s_%s' % ( a.name , b.name )
            vtitle  = 'Sum %s+%s' % ( a.name , b.name )
            
            var     = ROOT.RooFormulaVar ( vname , vtitle , formula , vlist )
            
            self.aux_keep.append ( vlist  )
            self.aux_keep.append ( var    )
            
            return self.soft_constraint ( var , value , name , title )

        elif isinstance ( a , ROOT.RooAbsReal ) and isinstance ( b , num_types + (VE,) ) :
            
            return self.soft_constraint ( a , value - b , name , title )
        
        elif isinstance ( b , ROOT.RooAbsReal ) and isinstance ( a , num_types + (VE,) ) :
            
            return self.soft_constraint ( b , value - a  , name , title )

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
                
        if isinstance ( a , ROOT.RooAbsReal ) and isinstance ( a , ROOT.RooAbsReal ) :

            vlist   = ROOT.RooArgList    ( a , b )
            
            formula = '(%s)*(%s)'     % ( a.name , b.name )
            vname   = 'Product_%s_%s' % ( a.name , b.name )
            vtitle  = 'Product %s+%s' % ( a.name , b.name )
            
            var     = ROOT.RooFormulaVar ( vname , vtitle , formula , vlist )
            
            self.aux_keep.append ( vlist  )
            self.aux_keep.append ( var    )
            
            return self.soft_constraint ( var , value , name , title )

        elif isinstance ( a , ROOT.RooAbsReal ) and isinstance ( b , num_types + (VE,) ) :
            
            return self.soft_constraint ( a , value / b , name , title )
        
        elif isinstance ( b , ROOT.RooAbsReal ) and isinstance ( a , num_types + (VE,) ) :
            
            return self.soft_constraint ( b , 1 / value , name , title )
        
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

        return self.soft_ratio_constraint ( a , b , 1.0/value - 1.0 , name , value )
        
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
        
        >>> var     = ...                            ## the variable 
        >>> extcntr = var.constaint( VE(1,0.1**2 ) ) ## create constrains 
        >>> model.fitTo ( ... , extcntr )            ## use it in the fit 
        """
        
        ## create the gaussian constraint
        gauss  = self.soft_constraint ( var , value , name ,  title ) 
        
        cnts   = ROOT.RooArgSet ( gauss )
        
        result = ROOT.RooFit.ExternalConstraints ( cnts )
        
        self.aux_keep.append ( cnts   )
        
        return result 

# =============================================================================
## simple convertor of 1D-histo to data set
#  @code
#  h   = ...
#  dset = H1D_dset ( h )
#  @endcode 
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H1D_dset(MakeVar) :
    """Simple convertor of 1D-histogram into data set
    >>> h   = ...
    >>> dset = H1D_dset ( h )
    """
    def __init__ ( self            , 
                   histo           ,
                   xaxis   = None  ,
                   density = True  ,
                   silent  = False ) :
        #
        ## use mass-variable
        #
        assert isinstance ( histo , ROOT.TH1 ) , "``histo'' is not ROOT.TH1"
        self.__histo = histo 

        name           = histo.GetName()
        self.__xaxis   = self.make_var ( xaxis , 'x_%s' % name , 'x-axis(%s)' % name , None , *(histo.xminmax()) )
        
        self.__density = True if density else False 
        self.__silent  = silent 
        
        with roo_silent ( self.silent ) :  
            
            self.__vlst = ROOT.RooArgList    ( self.xaxis )
            self.__vimp = ROOT.RooFit.Import ( self.histo , self.density )
            self.__dset = ROOT.RooDataHist   (
                rootID ( 'hds_' ) ,
                "Data set for histogram '%s'" % histo.GetTitle() ,
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

# =============================================================================
## simple convertor of 2D-histo to data set
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H2D_dset(MakeVar) :
    """Simple convertor of 2D-histogram into data set
    """
    def __init__ ( self            ,
                   histo           ,
                   xaxis   = None  ,
                   yaxis   = None  ,
                   density = True  ,
                   silent  = False ) :
        #
        assert isinstance ( histo , ROOT.TH2 ) , "``histo'' is not ROOT.TH2"
        self.__histo = histo

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
        
        with roo_silent ( silent ) : 

            self.__vlst  = ROOT.RooArgList    ( self.xaxis , self.yaxis )
            self.__vimp  = ROOT.RooFit.Import ( histo , density )
            self.__dset  = ROOT.RooDataHist   (
                rootID ( 'hds_' ) ,
                "Data set for histogram '%s'" % histo.GetTitle() ,
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

# =============================================================================
## simple convertor of 3D-histo to data set
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H3D_dset(MakeVar) :
    """Simple convertor of 3D-histogram into data set
    """
    def __init__ ( self            ,
                   histo           ,
                   xaxis   = None  ,
                   yaxis   = None  ,
                   zaxis   = None  ,
                   density = True  ,
                   silent  = False ) :
        
        assert isinstance ( histo , ROOT.TH3 ) , "``histo'' is not ROOT.TH3"
        self.__histo = histo
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
                
        with roo_silent ( silent ) : 

            self.__vlst  = ROOT.RooArgList    ( self.xaxis , self.yaxis , self.zaxis )
            self.__vimp  = ROOT.RooFit.Import ( histo , density )
            self.__dset  = ROOT.RooDataHist   (
                rootID ( 'hds_' ) ,
                "Data set for histogram '%s'" % histo.GetTitle() ,
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

# =============================================================================
## @class Phases
#  helper class to build/keep the list of ``phi''-arguments
#   - needed e.g. for polynomial functions
class Phases(MakeVar) :
    """Helper class to build/keep the list of ``phi''-arguments (needed e.g. for polynomial functions)
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
                                    None , 0 ,  -0.85 * pi  , 1.55 * pi  )
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

        from ostap.core.types import num_types , list_types
        ##
        if   isinstance ( values , num_types          ) : values = [ values           ]
        elif isinstance ( values , VE                 ) : values = [ values.value()   ]
        elif isinstance ( values , ROOT.RooAbsReal    ) : values = [ float ( values ) ] 
        elif isinstance ( values , list_types         ) : pass
        elif isinstance ( values , ROOT.RooArgList    ) : pass
        else :
            raise TypeError("Unknown type for ``values'' %s/%s" % (  values , type ( values ) ) )

        for s , v in  zip ( self.__phis , values ) :
            vv = float ( v  )
            if s.minmax() and not vv in s :
                logger.error ("Value %s is outside the allowed region %s"  % ( vv , s.minmax() ) )                 
            s.setVal   ( vv )
        nphi = len ( self.__phis )

    @property
    def phi_list ( self ) :
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
    if isinstance   ( what , ROOT.RooArgList ) and i in what             : return what[i]
    elif isinstance ( what , ROOT.RooAbsReal ) and 0 == i                : return what 
    elif isinstance ( what , num_types       ) and 0 == i                : return what  
    elif isinstance ( what , list_types      ) and 0 <= i < len ( what ) : return what[i] 

    return default
        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
# The END 
# =============================================================================
