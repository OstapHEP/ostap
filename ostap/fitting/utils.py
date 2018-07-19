#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file utils.py
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
    'fitArgs'           , ## Function to parse 'fitTo'-arguments
    #
    'RangeVar'          , ## Helper class to temporary change a range for the variable 
    'MakeVar'           , ## Helper bade class that allow storage of newly created RooFit objects
    #
    'Adjust1D'          , ## adjust 1D-pdf to avois zeroes
    #
    "H1D_dset"          , ## 1D-histogram to RooDataHist converter 
    "H2D_dset"          , ## 2D-histogram to RooDataHist converter 
    "H3D_dset"          , ## 3D-histogram to RooDataHist converter
    #
    'component_similar' , ## Should one use ``similar'' component?
    'component_clone'   , ## Should one use ``cloned'' component?
    )
# =============================================================================
import ROOT, math
from   ostap.core.core     import rootID 
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
_nemax = 20000  ## number of events per CPU-core 
_ncmax =     6  ## maximal number of CPUs: there are some problems with >= 7
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
def ncpu (  events ) :
    """Prepare 'NumCPU' argument with reasonable choice of #cpu, depending on
    the number of events in dataset 
    """
    #
    n  = events // _nemax
    if n       <= 1 : return ROOT.RooFit.Save () ## fake!!! 
    # 
    n_cores = numcpu() 
    if n_cores <= 1 : return ROOT.RooFit.Save () ## fake!!! 
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
## "parse" common 'fitTo'-arguments
def fitArgs ( name , dataset = None , *args , **kwargs ) :
    """ Parse 'fitTo'-arguments 
    """
    _args = []
    for a in args :
        if not isinstance ( a , ROOT.RooCmdArg ) :
            logger.warning ( '%s unknown argument type %s, skip it ' % ( name , type ( a ) ) ) 
            continue
        _args.append ( a )

    from ostap.plotting.fit_draw import keys  as drawing_options 
    ncpu_added = False 
    for k , a in kwargs.iteritems() :
        
        ## skip "drawing" options 
        if k.lower() in drawing_options                            : continue 
        if k.lower() in ( 'draw' , 'draw_option', 'draw_options' ) : continue 
        
        if isinstance ( a , ROOT.RooCmdArg ) :
            logger.debug   ( '%s add keyword argument %s' % ( name , k ) )  
            _args.append ( a )
        elif k.upper() in ( 'WEIGHTED'   ,
                            'SUMW2'      ,
                            'SUMW2ERROR' ) and isinstance ( a , bool ) and dataset.isWeighted() :
            _args.append   (  ROOT.RooFit.SumW2Error( a ) )
            logger.debug   ( '%s add keyword argument %s/%s' % ( name , k , a ) )                 
        elif k.upper() in ( 'EXTENDED' , ) and isinstance ( a , bool ) :
            _args.append   (  ROOT.RooFit.Extended ( a ) )
            logger.debug   ( '%s add keyword argument %s/%s' % ( name , k , a ) )                 
        elif k.upper() in ( 'NCPU'       ,
                            'NCPUS'      ,
                            'NUMCPU'     ,
                            'NUMCPUS'    ) and isinstance ( a , int ) and 1<= a : 
            _args.append   (  ROOT.RooFit.NumCPU( a  ) ) 
            logger.debug   ( '%s add keyword argument %s/%s' % ( name , k , a ) )
            ncpu_added = True
        elif k.upper() in ( 'CONSTRAINT'  ,
                            'CONSTRAINTS' ,
                            'PARS'        ,
                            'PARAMS'      ,
                            'PARAMETER'   ,
                            'PARAMETERS'  ) :
            if   isinstance ( a , ROOT.RooCmdArg ) : _args.append ( a )
            elif isinstance ( a , (tuple,list)   ) :
                for ia in a :
                    if isinstance ( ia , ROOT.RooCmdArg ) : _args.append ( ia )
                    else : logger.warning( '%s skip keyword argument [%s] : %s' % ( name , k , a ) )
            else : logger.warning( '%s skip keyword argument %s: %s' % ( name , k , a ) )
                                        
        else : 
            logger.warning ( '%s unknown/illegal keyword argument type %s/%s, skip it ' % ( name , k , type ( a ) ) )
            continue            
        
    if not ncpu_added :
        logger.debug  ( '%s: NCPU is added ' % name ) 
        _args.append  (  ncpu ( len ( dataset ) ) )
            
    return tuple ( _args )

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
        obj = super(MakeVar, cls).__new__( cls , *args , **kwargs )
        obj.__aux_keep = []                      ## ATTENTION!        
        obj.__name     = None                    ## ATTENTION!
        return obj
    
    ##  produce ERROR    message using the local logger 
    def error   ( self , message , *args , **kwargs ) : return self.logger.error   ( message , *args , **kwargs )
    ##  produce WARNIING message using the local logger 
    def warning ( self , message , *args , **kwargs ) : return self.logger.warning ( message , *args , **kwargs )
    ##  produce INFO     message using the local logger 
    def info    ( self , message , *args , **kwargs ) : return self.logger.info    ( message , *args , **kwargs )
    ##  produce DEBUG    message using the local logger 
    def debug   ( self , message , *args , **kwargs ) : return self.logger.debug   ( message , *args , **kwargs )
    ##  produce VERBOSE  message using the local logger 
    def verbose ( self , message , *args , **kwargs ) : return self.logger.verbose ( message , *args , **kwargs )
    
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
        if isinstance   ( var , ( float , int , long ) ) :
            if   not    args       : var = ROOT.RooRealVar ( name , comment , var             )
            elif 2 == len ( args ) : var = ROOT.RooRealVar ( name , comment , var , *args     )
            elif 3 == len ( args ) : var = ROOT.RooRealVar ( name , comment , var , *args[1:] )
        
        ## create the variable from parameters 
        if not isinstance ( var , ROOT.RooAbsReal ) : 
            var = ROOT.RooRealVar ( name , comment , *args )
            self.aux_keep.append ( var )        ##  ATTENTION: store newly created variable
        
        ## fix it, if needed
        if   isinstance ( fix , bool                   ) : pass 
        elif isinstance ( fix , ( float , int , long ) ) :
                
            if hasattr ( var , 'getMin)') and fix < var.getMin() and hasattr ( var , 'setMin' ) :                                                                             
                self.warning ( "make_var: min-value for %s is redefined to be %s " % ( var.GetName () , fix ) )
                var.setMin ( fix )
            
            if hasattr ( var , 'getMax' ) and fix > var.getMax() and hasattr ( var , 'setMax' ) : 
                self.warning ( "make_var: max-value for %s is redefined to be %s " % ( var.GetName () , fix ) )
                var.setMax ( fix )
            
            if not var.isConstant () and hasattr ( var , 'fix'    ) : var.fix    ( fix )
            elif                         hasattr ( var , 'setVal' ) : var.setVal ( fix )

        return var


# =============================================================================
## simple convertor of 1D-histo to data set
#  @code
#  h1   = ...
#  dset = H1D_dset ( h1 )
#  @endcode 
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H1D_dset(MakeVar) :
    """Simple convertor of 1D-histogram into data set
    >>> h1   = ...
    >>> dset = H1D_dset ( h1 )
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
## simple class to adjust certaint PDF to avoid zeroes 
class Adjust1D(MakeVar) :
    """Simple class to ``adjust'' certain PDF to avoid zeroes
    - a small flat component is added and the compound PDF is constructed
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,
                   xvar             , 
                   pdf              ,
                   value    = 1.e-5 ) : 

        assert isinstance ( pdf  , ROOT.RooAbsPdf  ) , "``pdf''  must be ROOT.RooAbsPdf"
        assert isinstance ( xvar , ROOT.RooAbsReal ) , "``xvar'' must be ROOT.RooAbsReal"
        
        self.name      = name 
        self.__old_pdf = pdf

        from otap.fitting.basic import Flat1D 
        self.__flat    = Flat1D  ( xvar  , name = 'flat_' + name ) 
        self.__frac    = self.make_var ( value , 'fracA_%s'                     % name ,
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
## simple class to adjust certaint PDF to avoid zeroes 
class Adjust2D(MakeVar) :
    """Simple class to ``adjust'' certain PDF to avoid zeroes
    - a small flat component is added and the compound PDF is constructed
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,
                   xvar             , 
                   yvar             , 
                   pdf              ,
                   value    = 1.e-5 ) : 
        
        assert isinstance ( pdf  , ROOT.RooAbsPdf  ) , "``pdf''  must be ROOT.RooAbsPdf"
        assert isinstance ( xvar , ROOT.RooAbsReal ) , "``xvar'' must be ROOT.RooAbsReal"
        assert isinstance ( yvar , ROOT.RooAbsReal ) , "``yvar'' must be ROOT.RooAbsReal"
        
        self.name      = name 
        self.__old_pdf = pdf
        
        from otap.fitting.fit2d import Flat2D 
        self.__flat    = Flat2D  ( xvar , yvar , name = 'flat_' + name )
        self.__frac    = self.make_var ( value , 'fracA_%s'                     % name ,
                                   'small  fraction of flat component %s' % name ,
                                   value , 1.e-4 , 0 , 1 )
        
        self.__alist1  = ROOT.RooArgList ( self.__flat.pdf , self.__old_pdf )        
        self.__alist2  = ROOT.RooArgList ( self.__frac     )
        #
        ## final PDF
        # 
        self.__pdf     = ROOT.RooAddPdf  ( "adjust_"    + name ,
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
## simple class to adjust certaint PDF to avoid zeroes 
class Adjust3D(MakeVar) :
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
        
        from otap.fitting.fit3d import Flat3D        
        ##
        self.__flat    = Flat3D  ( xvar , yvar ,  xvar , name = 'flat_' + name )
        self.__frac    = self.make_var ( value , 'fracA_%s'                     % name ,
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
## @class Phases
#  helper class to build/keep the list of ``phi''-arguments (needed e.g. for polynomial functions)
class Phases(MakeVar) :
    """Helper class to build/keep the list of ``phi''-arguments (needed e.g. for polynomial functions)
    """
    ## Create vector of phases (needed for various polynomial forms)
    def __init__( self  , power , the_phis = None ) :
        """Create vector of phases (needed for various polynomial forms)
        """

        ## check  the arguments 
        assert isinstance ( power , (int ,long)) and 0 <= power, \
               "Phases: invalid type/value for ``power''-parameter: %s/%s"  % (  power , type(power) )

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

        nphi = len(self.__phis )
        
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

    @property
    def phi_list ( self ) :
        """The list/ROOT.RooArgList of ``phases'', used to parameterize polynomial-like shapes
        """
        return self.__phi_list
 
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
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
# The END 
# =============================================================================
