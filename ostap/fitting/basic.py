#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file basic.py
#  Set of useful basic utilities to build various fit models 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
# =============================================================================
"""Set of useful basic utilities to build various fit models"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'makeVar'       , ## helper function to create the proper RooRealVar
    'makeBkg'       , ## helper function to create smooth ``background''
    'H1D_dset'      , ## convertor of 1D-histo to RooDataHist 
    'H1D_pdf'       , ## convertor of 1D-histo to RooHistPdf 
    ##
    'PDF'           , ## useful base class for 1D-models
    'MASS'          , ## useful base class to create "signal" PDFs for mass-fits
    'RESOLUTION'    , ## useful base class to create "resolution" PDFs 
    'Fit1D'         , ## the model for 1D-fit: signal + background + optional components  
    ##
    'Adjust'        , ## addjust PDF to avoid zeroes (sometimes useful)
    'Convolution'   , ## helper utility to build convolution
    'Phases'        , ## helper utility to build/keep list of phases 
    ##
    'Generic1D_pdf' , ## wrapper over imported RooFit (1D)-pdf  
    )
# =============================================================================
import ROOT, math
from   ostap.core.core     import cpp , Ostap , VE , hID , rootID
from   ostap.histos.histos import h1_axis , h2_axes
from   ostap.logger.utils  import roo_silent 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.basic' )
else                       : logger = getLogger ( __name__         )
# =============================================================================
_nemax = 20000  ## number of events per CPU-core 
_ncmax =     6  ## maximal number of CPUs: there are some problems with >= 7
                ## @see https://sft.its.cern.ch/jira/browse/ROOT-4897
_ncpus = []
# =============================================================================
## MINUIT covarinace matrix status:
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
## MINUIT covarinace matrix status:
# - status = -1 :  not available (inversion failed or Hesse failed)
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
## prepare "NumCPU" argument with reasonable choice of number of cpus, depending on
#  number of events in dataset 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-03-31
def ncpu (  events ) :
    """Prepare 'NumCPU' argument with reasonable choice of N-cpu, depending on
    the number of events in dataset 
    """
    #
    n  = events // _nemax
    if n       <= 1 : return ROOT.RooFit.Save () # fake!!! 
    # 
    import multiprocessing
    n_cores = multiprocessing.cpu_count()
    if n_cores <= 1 : return ROOT.RooFit.Save () # fake!!! 
    #
    num = min ( n , n_cores , _ncmax )
    if not _ncpus : _ncpus.append ( num )   
    #
    return ROOT.RooFit.NumCPU ( num )

# =============================================================================
## create/modify  the variable
#  Helper function for creation/modification/adjustment of variable
#  @code
#    v = makeVar ( 10   , 'myvar' , 'mycomment' )
#    v = makeVar ( 10   , 'myvar' , 'mycomment' , '' ,     -1 , 1 )
#    v = makeVar ( 10   , 'myvar' , 'mycomment' , '' , 0 , -1 , 1 )
#    v = makeVar ( None , 'myvar' , 'mycomment' , '' , 0 , -1 , 1 )
#    v = makeVar ( None , 'myvar' , 'mycomment' , 10 , 0 , -1 , 1 )
#    v = makeVar ( v    , 'myvar' , 'mycomment' , 10 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-12-01
def makeVar ( var , name , comment , fix = None , *args ) :
    """Make/modify  the variable:
    
    v = makeVar ( 10   , 'myvar' , 'mycomment' )
    v = makeVar ( 10   , 'myvar' , 'mycomment' , '' ,     -1 , 1 )
    v = makeVar ( 10   , 'myvar' , 'mycomment' , '' , 0 , -1 , 1 )
    v = makeVar ( None , 'myvar' , 'mycomment' , '' , 0 , -1 , 1 )
    v = makeVar ( None , 'myvar' , 'mycomment' , 10 , 0 , -1 , 1 )
    v = makeVar ( v    , 'myvar' , 'mycomment' , 10 )
    """
    #
    # var = ( value )
    # var = ( min , max )
    # var = ( value , min , max ) 
    if   isinstance   ( var , tuple ) :
        var = ROOT.RooRealVar ( name , comment , *var )

    # var = value 
    if isinstance   ( var , ( float , int , long ) ) :
        if   not    args  : var = ROOT.RooRealVar ( name , comment , var             )
        elif 2==len(args) : var = ROOT.RooRealVar ( name , comment , var , *args     )
        elif 3==len(args) : var = ROOT.RooRealVar ( name , comment , var , *args[1:] )
        
    ## create the variable from parameters 
    if not isinstance ( var , ROOT.RooAbsReal ) : 
        var = ROOT.RooRealVar ( name , comment , *args )
        
    ## fix it, if needed:
        
    if   fix is False : pass
    elif fix is True  : var.fix ( var.getVal() )
    elif isinstance ( fix , ( float , int , long ) ) :
        
        if hasattr ( var , 'getMin' ) and fix < var.getMin() :
            logger.warning("Min-value for %s is redefined to be %s " % ( var.GetName() , fix ) )
            var.setMin ( fix )
            
        if hasattr ( var , 'getMax' ) and fix > var.getMax() :
            logger.warning("Max-value for %s is redefined to be %s " % ( var.GetName() , fix ) )
            var.setMax ( fix )
            
        if not var.isConstant () : var.fix    ( fix )
        else                     : var.setVal ( fix )

    return var


# =============================================================================
## make a RooArgList of variables/fractions 
def makeFracs ( N , pname , ptitle , model , fractions = True )  :
    """Make a RooArgList of variables/fractions 
    """
    if not isinstance ( N , (int,long) ) : raise TypeError('Invalid N type %s' % type(N) )
    elif   N < 2                         : raise TypeError('Invalid N=%d'      %      N  )
    ##
    fracs = ROOT.RooArgList()
    n     = (N-1) if fractions else N
    lst   = []
    for i in range(1,n+1) :
        if fractions : fi = makeVar ( None , pname  % i , ptitle % i , None , 1.0/N , 0 , 1     )
        else         : fi = makeVar ( None , pname  % i , ptitle % i , None , 1     , 0 , 1.e+6 )
        fracs.add  ( fi )
        setattr ( model , fi.GetName() , fi ) 
    ##    
    return fracs

# =============================================================================
## helper function to build composite (non-extended) PDF from components 
def addPdf ( pdflist , name , title , fname , ftitle , model , recursive = True ) :
    """Helper function to build composite (non-extended) PDF from components 
    """
    ##
    pdfs  = ROOT.RooArgList() 
    for pdf in pdflist : pdfs.add  ( pdf )
    fracs = makeFracs ( len ( pdfs ) , fname , ftitle , fractions = True , model = model ) 
    pdf   = ROOT.RooAddPdf ( name , title , pdfs, fracs , recursive )
    ##
    setattr ( model , '__addPdf_'       + name , pdf   )
    setattr ( model , '__addPdf_lst_'   + name , pdfs  )
    setattr ( model , '__addPdf_fracs_' + name , fracs )
    ##
    return pdf, fracs , pdfs
  
# =============================================================================
## helper class to temporary change a range for the variable 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class RangeVar(object) :
    """ Helper class to temporary change a range for the variable 
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
## "parse" common arguments for fit 
def fitArgs ( name , dataset = None , *args , **kwargs ) :
    """ ``parse'' common arguments for fit 
    """
    _args = []
    for a in args :
        if not isinstance ( a , ROOT.RooCmdArg ) :
            logger.warning ( '%s unknown argument type %s, skip it ' % ( name , type ( a ) ) ) 
            continue
        _args.append ( a )
        
    from ostap.plotting.fit_draw import keys     
    ncpu_added = False 
    for k,a in kwargs.iteritems() :

        ## skip "drawing" options 
        if k.lower() in keys                                       : continue 
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
## @class PDF
#  The helper base class for implementation of 1D-pdfs 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-08-21
class PDF (object) :
    """Useful helper base class for implementation of PDFs for 1D-fit
    """
    def __init__ ( self , name , xvar = None ) :
        self.name         = name
        self._signals     = ROOT.RooArgList ()
        self._backgrounds = ROOT.RooArgList ()
        self._components  = ROOT.RooArgList ()
        ## take care about sPlots 
        self._splots      = []
        self._properties  = {}
        
        if   isinstance ( xvar, ROOT.TH1   ) :
            xvar = xvar.xminmax()
        elif isinstance ( xvar , ROOT.TAxis ) :
            xvar = xvar.GetXmin() , xvar.GetXmax()

        self.__xvar = None 
        ## create the variable 
        if isinstance ( xvar , tuple ) and 2 == len(xvar) :  
            self.__xvar = makeVar ( xvar               , ## var 
                                    'x'                , ## name 
                                    'x-variable(mass)' , ## title/comment
                                    *xvar              , ## min/max 
                                    fix = None         ) ## fix ? 
        elif isinstance ( xvar , ROOT.RooAbsReal ) :
            self.__xvar = makeVar ( xvar               , ## var 
                                    'x'                , ## name 
                                    'x-variable/mass'  , ## title/comment
                                    fix = None         ) ## fix ? 
        else :
            ##logger.warning('x-variable is not specified (yet)')
            self.__xvar = makeVar( xvar , 'x' , 'x-variable' )
        

    @property 
    def xvar ( self ) :
        """``x''-variable for the fit (same as ``x'')"""
        return self.__xvar

    @property 
    def x    ( self ) :
        """``x''-variable for the fit (same as ``xvar'')"""
        return self.__xvar
    
    @property
    def xminmax ( self ) :
        """Min/max values for x-variable"""
        return self.__xvar.minmax()
        
    ## get all declared components 
    def components  ( self ) : return self._components
    ## get all declared signals 
    def signals     ( self ) : return self._signals
    ## get all declared backgrounds 
    def backgrounds ( self ) : return self._backgrounds 
    ## can be useful 
    def setPars     ( self ) :

        ## check own PDF 
        if hasattr ( self , 'pdf' ) :
            c = self.pdf
            if hasattr ( c , 'setPars' ) : c.setaPars()
            
        ## check other PDFs 
        cs = set() 
        for att in ( '_signals'     , '_backgrounds'     , '_components'     ,
                     'all_signals'  , 'all_backgrounds'  , 'all_components'  ,
                     'more_signals' , 'more_backgrounds' , 'more_components' ,
                     '_sigs'        , '_bkgs'            , '_cmps'           ) :
            if hasattr ( self , att ) :
                for a in getattr ( self , att ) : 
                    if hasattr ( a , 'setPars' ) and not a in cs :  
                        a.setPars()
                        cs.add ( a )
        del cs

    ## adjust PDF a little bit to avoid zeroes 
    def adjust ( self , value ) :
        """Adjust PDF a little bit to avoid zeroes 
        """
        if hasattr ( self , 'adjusted' ) :
            logger.warning ( "PDF is already adjusted!")
            return
        
        self.adjusted = Adjust ( self.name , self.mass , self.pdf , value ) 
        self.old_pdf  = self.adjusted.old_pdf
        self.pdf      = self.adjusted.pdf
        

    # =========================================================================
    ## make the actual fit (and optionally draw it!)
    #  @code
    #  r,f = model.fitTo ( dataset )
    #  r,f = model.fitTo ( dataset , weighted = True )    
    #  r,f = model.fitTo ( dataset , ncpu     = 10   )    
    #  r,f = model.fitTo ( dataset , draw = True , nbins = 300 )    
    #  @endcode 
    def fitTo ( self           ,
                dataset        ,
                draw   = False ,
                nbins  = 100   ,
                silent = False ,
                refit  = False , *args , **kwargs ) :
        """Perform the actual fit (and draw it)
        >>> r,f = model.fitTo ( dataset )
        >>> r,f = model.fitTo ( dataset , weighted = True )    
        >>> r,f = model.fitTo ( dataset , ncpu     = 10   )    
        >>> r,f = model.fitTo ( dataset , draw = True , nbins = 300 )    
        """
        if isinstance ( dataset , ROOT.TH1 ) :
            density = kwargs.pop ( 'density' , True   )
            chi2    = kwargs.pop ( 'chi2'    , False  )
            return self.fitHisto ( dataset , draw , silent , density , chi2 , *args , **kwargs ) 
        #
        ## treat the arguments properly
        #
        opts = fitArgs ( "PDF(%s).fitTo:" % self.name , dataset , *args , **kwargs )
        
        #
        ## define silent context
        with roo_silent ( silent ) :
            result =  self.pdf.fitTo ( dataset , ROOT.RooFit.Save () , *opts ) 
            if hasattr ( self.pdf , 'setPars' ) : self.pdf.setPars() 
                         
        st = result.status() if result else 9999    
        if 0 != st and silent and refit : 
            logger.warning ( 'PDF(%s).fitTo: status is %s. Refit in non-silent regime ' % ( self.name , fit_status ( st ) ) )
            return self.fitTo ( dataset , draw , nbins , False , refit , *args , **kwargs )
        
        for_refit = False
        if 0 != st   :
            for_refit = 'status' 
            logger.warning ( 'PDF(%s).fitTo: Fit status is %s ' % ( self.name , fit_status ( st ) ) )
        
        qual = result.covQual() if result else 9999 
        if   -1 == qual and dataset.isWeighted() : pass
        elif  3 != qual :
            for_refit = 'covariance' 
            logger.warning ( 'PDF(%s).fitTo: covQual    is %s ' % ( self.name , cov_qual ( qual ) ) ) 

        #
        ## check the integrals (when possible)
        #
        if hasattr ( self , 'nums' ) and self.nums :
            
            nsum = VE()            
            for i in self.nums :
                nsum += i.as_VE() 
                if hasattr ( i , 'getMax' ) and i.getVal() > i.getMax() - 0.05 * ( i.getMax() - i.getMin() ) :
                    logger.warning ( 'PDF(%s).fitTo: Variable %s == %s [too close to maximum %s]'
                                     % ( self.name , i.GetName() , i.getVal () , i.getMax () ) )
                    
            if not dataset.isWeighted() :
                if 0 < nsum.cov2() : 
                    nl = nsum.value() - 0.80 * nsum.error()
                    nr = nsum.value() + 0.80 * nsum.error()
                    if not nl <= len ( dataset ) <= nr :
                        logger.warning ( 'PDF(%s).fitTo: is problematic:  sum %s != %s ' % ( self.name , nsum , len( dataset ) ) )
                        for_refit = 'integral'

        #
        ## call for refit if needed
        #
        if refit and for_refit :
            logger.info ( 'PDF(%s).fitTo: call for refit:  %s/%s'  % ( self.name , for_refit , refit ) ) 
            if isinstance ( refit , ( int , long ) )  : refit -= 1
            else                                      : refit  = False
            return  self.fitTo ( dataset , draw , nbins , silent , refit , *args , **kwargs ) 


        ## draw it if requested
        from ostap.plotting.fit_draw import draw_options
        draw_opts = draw_options ( **kwargs )
        if draw_opts and not draw     : draw = draw_opts
        if isinstance ( draw , dict ) : draw_opts.update( draw )
        
        frame = self.draw ( dataset , nbins = nbins , silent = silent , **draw_opts ) if draw else None 

        if hasattr ( self.pdf , 'setPars' ) : self.pdf.setPars()
            
        for s in self.components  () : 
            if hasattr ( s , 'setPars' ) : s.setPars()
        for s in self.backgrounds () :  
            if hasattr ( s , 'setPars' ) : s.setPars() 
        for s in self.signals     () : 
            if hasattr ( s , 'setPars' ) : s.setPars() 

        ## 
        return result, frame 

    # =========================================================================
    ## ROOT-like function for fitting
    #  @code
    #  dataset = ...
    #  model   = ...
    #  model.Fit ( dataset , ... ) 
    #  @endcode 
    def Fit ( self , dataset , *args , **kwargs ) :
        """ROOT-like function for fitting
        >>> dataset = ...
        >>> model   = ...
        >>> model.Fit ( dataset , ... ) 
        """
        return self.fitTo ( dataset , *args , **kwargs ) 
        
    ## helper method to draw set of components 
    def _draw ( self , what , frame , options , base_color ) :
        """ Helper method to draw set of components """
        i = 0 
        for cmp in what : 
            cmps = ROOT.RooArgSet( cmp )
            if 0 <= base_color : 
                self.pdf .plotOn ( frame ,
                                   ROOT.RooFit.Components ( cmps                        ) ,
                                   ROOT.RooFit.LineColor  ( base_color + i ) , *options )
            else :
                self.pdf .plotOn ( frame ,
                                   ROOT.RooFit.Components ( cmps ) , *options )
                
            i += 1
        
    # ================================================================================
    ## draw fit results
    #  @code
    #  r,f = model.fitTo ( dataset )
    #  model.draw ( dataset , nbins = 100 ) 
    #  @endcode
    #  @param dataset  dataset to be drawn 
    #  @param nbins    binning scheme for frame/RooPlot 
    #  @param silent   silent mode ?
    #  @param data_options          drawing options for dataset
    #  @param signal_options        drawing options for `signal'     components    
    #  @param background_options    drawing options for `background' components 
    #  @param component_options     drawing options for 'other'      components 
    #  @param fit_options           drawing options for fit curve    
    #  @param base_signal_color     base color for signal components 
    #  @param base_background_color base color for background components
    #  @param base_component_color  base color for other components 
    #  @param data_overlay          draw points atop of fitting curves?  
    #  @see ostap.plotting.fit_draw 
    def draw ( self ,
               dataset               = None  ,
               nbins                 = 100   ,   ## Frame binning
               silent                = False ,   ## silent mode ?
               **kwargs                      ) :
        """Visualize the fits results
        >>> r,f = model.draw ( dataset )
        >>> model.draw ( dataset , nbins = 100 )
        >>> model.draw ( dataset , base_signal_color  = ROOT.kGreen+2 )
        >>> model.draw ( dataset , data_options = (ROOT.RooFit.DrawOptions('zp'),) )

        Produce also residual & pull:
        
        >>> f,r,h = model.draw ( dataset , nbins = 100 , residual = True , pull = True )
        
        Drawing options:
        - data_options            ## drawing options for dataset  
        - background_options      ## drawing options for background component(s)
        - crossterm1_options      ## drawing options for crossterm1 component(s)
        - crossterm2_options      ## drawing options for crossterm2 component(s)
        - signal_options          ## drawing options for signal     component(s)
        - component_options       ## drawing options for other      component(s)
        - total_fit_options       ## drawing options for the total fit curve
        
        Colors:                  
        - base_background_color   ## base color for background component(s)
        - base_crossterm1_color   ## base color for crossterm1 component(s)
        - base_crossterm2_color   ## base color for crossterm2 component(s)
        - base_signal_color       ## base color for signal     component(s)
        - base_component_color    ## base color for other      component(s)
        
        Other options:
        -  residual               ## make also residual frame
        -  pull                   ## make also residual frame
 
        For default values see ostap.plotting.fit_draw       
        
        """
        #
        import ostap.plotting.fit_draw as FD
        #
        ## again the context
        # 
        with roo_silent ( silent ) :
            
            if hasattr ( self , 'draw_var' ) and self.draw_var : drawvar = self.draw_var
            else                                               : drawvar = self.mass
            
            if nbins :  frame = drawvar.frame ( nbins )
            else     :  frame = drawvar.frame ()
            
            #
            ## draw invizible data (for normalzation of fitting curves)
            #
            data_options = kwargs.pop ( 'data_options' , FD.data_options )
            if dataset : dataset .plotOn ( frame , ROOT.RooFit.Invisible() , *data_options )
            
            ## draw various ``background'' terms
            boptions = kwargs.pop (     'background_options' , FD.   background_options )
            bbcolor  = kwargs.pop (  'base_background_color' , FD.base_background_color ) 
            if self.backgrounds () :
                self._draw( self.backgrounds() , frame , boptions , bbcolor )
                
            ## ugly :-(
            ct1options = kwargs.pop (     'crossterm1_options' , FD.   crossterm1_options )
            ct1bcolor  = kwargs.pop (  'base_crossterm1_color' , FD.base_crossterm1_color ) 
            if hasattr ( self , 'crossterms1' ) and self.crossterms1() : 
                self._draw( self.crossterms1() , frame , ct1options , ct1bcolor )

            ## ugly :-(
            ct2options = kwargs.pop (     'crossterm2_options' , FD.   crossterm1_options )
            ct2bcolor  = kwargs.pop (  'base_crossterm2_color' , FD.base_crossterm1_color )                 
            if hasattr ( self , 'crossterms2' ) and self.crossterms2() :
                self._draw( self.crossterms2() , frame , ct2options , ct2bcolor )

            ## draw ``other'' components
            coptions   = kwargs.pop (      'component_options' , FD.    component_options  )
            cbcolor    = kwargs.pop (   'base_component_color' , FD. base_component_color  ) 
            if self.components () :
                self._draw( self.components() , frame , coptions , cbcolor )

            ## draw ``signal'' components
            soptions   = kwargs.pop (         'signal_options' , FD.      signal_options  )
            sbcolor    = kwargs.pop (      'base_signal_color' , FD.   base_signal_color  ) 
            if self.signals    () :
                self._draw( self.signals() , frame , soptions , sbcolor )

            #
            ## the total fit curve
            #
            self.pdf .plotOn ( frame , *kwargs.pop ( 'total_fit_options' , FD. total_fit_options  ) )
            
            #
            ## draw data once more
            #
            if dataset : dataset  .plotOn ( frame , *data_options )            

            #
            ## suppress ugly axis labels
            #
            frame.SetXTitle ( '' )
            frame.SetYTitle ( '' )
            frame.SetZTitle ( '' )            
            
            #
            ## Draw the frame!
            #
            frame.Draw()
            
            residual =  kwargs.pop ( 'residual' , False )
            if residual and not  dataset :
                logger.warning("Can't produce residual without data")
                residual = False
                
            pull     =  kwargs.pop ( 'pull'     , False ) 
            if pull     and not  dataset :
                logger.warning("Can't produce residual without data")
                residual = False
                
            if kwargs :
                logger.warning("Ignored unknown options: %s" % kwargs.keys() )

            if not residual and not pull:
                return frame

            rframe =  None 
            if residual  :
                if   residual is True               : residual =      "P" ,
                elif isinstance  ( residual , str ) : residual = residual ,
                rframe  = frame.emptyClone ( rootID ( 'residual_' ) )
                rh      = frame.residHist()
                rframe.addPlotable ( rh , *residual ) 
                rframe.SetXTitle ( '' )
                rframe.SetYTitle ( '' )
                rframe.SetZTitle ( '' )
                
            pframe = None 
            if pull      : 
                if   pull is True               : pull =   "P",
                elif isinstance  ( pull , str ) : pull = pull ,
                pframe  = frame.emptyClone ( rootID ( 'pull_' ) )
                ph      = frame.pullHist()
                pframe.addPlotable ( ph , *pull ) 
                pframe.SetXTitle ( '' )
                pframe.SetYTitle ( '' )
                pframe.SetZTitle ( '' )
                
            return frame, rframe, pframe  

    # =========================================================================
    ## fit the histogram (and draw it)
    #  @code
    #  histo = ...
    #  r,f = model.fitHisto ( histo , draw = True ) 
    #  @endcode 
    def fitHisto ( self , histo , draw = False , silent = False , density = True , chi2 = False , *args , **kwargs ) :
        """Fit the histogram (and draw it)

        >>> histo = ...
        >>> r,f = model.fitHisto ( histo , draw = True )         
        """
        with RangeVar( self.mass , *(histo.xminmax()) ) : 
            
            ## convert it! 
            hdset     = H1D_dset ( histo , self.mass , density , silent )
            self.hset = hdset.dset
            
            if chi2 : return self.chi2fitTo ( self.hset , draw , None , silent , density , *args , **kwargs )
            else    : return self.fitTo     ( self.hset , draw , None , silent ,           *args , **kwargs )

    # =========================================================================
    ## make chi2-fit for binned dataset or histogram
    #  @code
    #  histo = ...
    #  r,f = model.chi2FitTo ( histo , draw = True ) 
    #  @endcode
    #  @todo add proper parsing of arguments for RooChi2Var 
    def chi2fitTo ( self,  dataset , draw = False , silent = False , density = True , *args , **kwargs ) :
        """ Chi2-fit for binned dataset or histogram
        >>> histo = ...
        >>> result , frame  = model.chi2FitTo ( histo , draw = True ) 
        """
        hdataset = dataset
        histo    = None 
        if isinstance  ( dataset , ROOT.TH1 ) :
            # if histogram, convert it to RooDataHist object:
            xminmax = dataset.xminmax() 
            with RangeVar( self.mass , *xminmax ) :                
                hdset     = H1D_dset ( dataset , self.mass , density , silent )
                self.hset = hdset.dset
                hdataset  = self.hset  
                histo     = dataset 
                
        with roo_silent ( silent ) : 

            lst1 = list ( fitArgs ( self.name , hdataset , *args , **kwargs ) )
            lst2 = []
            
            if       self.pdf.mustBeExtended () : lst2.append ( ROOT.RooFit.Extended ( True  ) )
            elif not self.pdf.canBeExtended  () : lst2.append ( ROOT.RooFit.Extended ( False ) )
            
            if not silent : lst2.append ( ROOT.RooFit.Verbose  () )
            if histo :
                if histo.natural() : lst2.append ( ROOT.RooFit.DataError ( ROOT.RooAbsData.Poisson ) )
                else               : lst2.append ( ROOT.RooFit.DataError ( ROOT.RooAbsData.SumW2   ) )  

            args_ = tuple ( lst2 + lst1  )
            #
            chi2 = ROOT.RooChi2Var ( rootID ( "chi2_" ) , "chi2(%s)" % self.name  , self.pdf , hdataset , *args_ )
            m    = ROOT.RooMinuit  ( chi2 ) 
            m.migrad   () 
            m.hesse    ()
            result = m.save ()
            
        from ostap.plotting.fit_draw import draw_options
        draw_opts = draw_options ( **kwargs )
        if draw_opts and not draw : draw = draw_opts 
        if isinstance ( draw , dict ) : draw_opts.update( draw )

        if not draw :
            return result, None 

        return result, self.draw ( hdataset , nbins = None , silent = silent , **draw_opts )

    # =========================================================================
    ## perform sPlot-analysis 
    #  @code
    #  r,f = model.fitTo ( dataset )
    #  model.sPlot ( dataset ) 
    #  @endcode 
    def sPlot ( self , dataset , silent = False ) : 
        """ Make sPlot analysis
        >>> r,f = model.fitTo ( dataset )
        >>> model.sPlot ( dataset ) 
        """
        if not hasattr ( self , 'alist2' ) :
            logger.error ('PDF(%s) has not attribute "alist2", no sPlot is possible' % self.name ) 
            raise AttributeError('PDF(%s) his not equipped for sPlot'                % self.name )

        with roo_silent ( silent ) :
            
            splot = ROOT.RooStats.SPlot ( rootID( "sPlot_" ) ,
                                          "sPlot"            ,
                                          dataset            ,
                                          self.pdf           ,
                                          self.alist2        )
        
            self._splots += [ splot ]
            
            return splot 

    ## simple 'function-like' interface
    #  @code
    #  pdf = ...
    #  x   = 0.45
    #  print 'Value of PDF at x=%f is %f' % ( x , pdf ( x ) ) 
    #  @endcode
    def __call__ ( self , x ) :
        """ Simple ``function-like// interface
        >>> pdf = ...
        >>> x = 0.45
        >>> print 'Value of PDF at x=%f is %f' % ( x , pdf ( x ) ) 
       """
        if isinstance ( self.mass , ROOT.RooRealVar ) :
            from ostap.fitting.roofit import SETVAR
            if x in self.mass :
                with SETVAR( self.mass ) :
                    self.mass.setVal ( x )
                    return self.pdf.getVal()
            else : return 0.0
            
        raise AttributeError, 'something wrong goes here'

    # ========================================================================
    # some generic stuff 
    # ========================================================================
    ## helper  function to implement some math stuff 
    def _get_stat_ ( self , funcall , *args , **kwargs ) :
        """Helper  function to implement some math stuff 
        """
        pdf         = self.pdf
        xmin , xmax = self.mass.minmax()
        
        if hasattr ( pdf , 'function' ) :    
            fun = pdf.function()
            if   hasattr ( pdf  , 'setPars'   ) : pdf.setPars()             
        else :            
            from ostap.fitting.roofit import PDF_fun
            fun = PDF_fun ( pdf , self.mass , xmin , xmax )
            
        return funcall (  fun , xmin , xmax , *args , **kwargs ) 
        
    
    ## get the effective RMS 
    def rms ( self ) :
        """Get the effective RMS
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'RMS: %s ' % pdf.rms()
        """
        pdf  = self.pdf 
        if   hasattr ( pdf , 'rms'            ) : return pdf.rms()        
        elif hasattr ( pdf , 'Rms'            ) : return pdf.Rms()        
        elif hasattr ( pdf , 'RMS'            ) : return pdf.RMS()        
        elif hasattr ( pdf , 'function'       ) :
            
            fun = pdf.function()
            if   hasattr ( pdf  , 'setPars'   ) : pdf.setPars() 
            if   hasattr ( fun , 'rms'        ) : return fun.rms()
            elif hasattr ( fun , 'variance'   ) : return fun.variance   ()**0.5  
            elif hasattr ( fun , 'dispersion' ) : return fun.dispersion () **5 
            
        from ostap.stats.moments import rms as _rms
        return  self._get_stat_ ( _rms )

    ## get the effective Skewness
    def skewness ( self ) :
        """Get the effective Skewness
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'SKEWNESS: %s ' % pdf.skewness()
        """
        ## use generic machinery 
        from ostap.stats.moments import skewness as _skewness
        return self._get_stat_ ( _skewness )


    ## get the effective Kurtosis
    def kurtosis ( self ) :
        """Get the effective Kurtosis
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'KURTOSIS: %s ' % pdf.kurtosis()
        """
        ## use generic machinery 
        from ostap.stats.moments import kurtosis as _kurtosis
        return self._get_stat_ ( _kurtosis ) 

    
    ## get the effective Full Width at Half Maximum
    def fwhm ( self ) :
        """Get the effective Full Width at  Half Maximum
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'FWHM: %s ' % pdf.fwhm()
        """
        ## use generic machinery 
        from ostap.stats.moments import width as _width
        w = self._get_stat_ ( _width )
        return w[1]-w[0]
    
    ## get the effective mode 
    def mode ( self ) :
        """Get the effective mode
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'MODE: %s ' % pdf.mode()
        """
        from ostap.stats.moments import mode as _mode
        return self._get_stat_ ( _mode )

    ## get the effective median
    def median ( self ) :
        """Get the effective median
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'MEDIAN: %s ' % pdf.median()
        """
        from ostap.stats.moments import median as _median
        return self._gets_stat_ ( _median )

    ## get the effective mean
    def get_mean ( self ) :
        """Get the effective Mean
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'MEAG: %s ' % pdf.get_mean()
        """
        from ostap.stats.moments import mean as _mean
        return self._get_stat_ ( _mean )

    ## get the effective moment for the distribution
    def moment ( self , N ) :
        """Get the effective moment
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'MOMENT: %s ' % pdf.moment( 10 )
        """
        ## use generic machinery 
        from ostap.stats.moments import moment as _moment
        return self._get_stat_ ( _moment , N ) 
    
    ## get the effective central moment for the distribution
    def central_moment ( self , N ) :
        """Get the effective central moment
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'MOMENT: %s ' % pdf.moment( 10 )
        """
        from ostap.stats.moments import central_moment as _moment
        return self._get_stat_ ( _moment , N ) 

    ## get the effective quantile 
    def quantile ( self , prob  ) :
        """Get the effective quantile
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'QUANTILE: %s ' % pdf.quantile ( 0.10 )
        """
        from ostap.stats.moments import quantile as _quantile
        return self._get_stat_ ( quantile , prob ) 
    
    ## get the symmetric confidence interval 
    def cl_symm ( self , prob , x0 =  None ) :
        """Get the symmetric confidence interval 
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'CL :  ',  pdf.cl_symm ( 0.10 )
        """
        from ostap.stats.moments import cl_symm as _cl
        return self._get_sstat_ ( _cl , prob , x0 ) 

    ## get the asymmetric confidence interval 
    def cl_asymm ( self , prob ) :
        """Get the asymmetric confidence interval 
        >>>  pdf = ...
        >>>  pdf.fitTo ( ... )
        >>>  print 'CL :  ',  pdf.cl_asymm ( 0.10 )
        """
        from ostap.stats.moments import cl_asymm as _cl
        return self._get_sstat_ ( _cl , prob ) 

    ## get the integral between xmin and xmax 
    def integral ( self , xmin , xmax ) :
        """Get integral between xmin and xmax
        >>> pdf = ...
        >>> print pdf.integral ( 0 , 10 )
        """
        ## check limits 
        if hasattr ( self.mass , 'getMin' ) :
            xmin = max ( xmin , self.mass.getMin() )
            xmax = max ( xmax , self.mass.getMin() )
        if hasattr ( self.mass , 'getMax' ) :
            xmin = min ( xmin , self.mass.getMax() )
            xmax = min ( xmax , self.mass.getMax() )
            
        if xmin == xmax : return 0 
            
        ## make a try to use the analytical integral 
        if hasattr ( self , 'pdf' ) :
            _pdf = self.pdf 
            if hasattr ( _pdf , 'setPars'  ) : _pdf.setPars() 
            try: 
                if hasattr ( _pdf , 'function' ) :
                    _func = _pdf.function() 
                    if hasattr ( _func , 'integral' ) :
                        return _func.integral ( xmin , xmax )
            except:
                pass
            
        ## use numerical ingeration
        from ostap.math.intergal import integral as _integral
        return _integral ( self , xmin , xmax )

    ## get the derivative at  point x 
    def derivative ( self , x ) :
        """Get derivative at point x 
        >>> pdf = ...
        >>> print pdf.derivative ( 0 ) 
        """
        ## check limits 
        if hasattr ( self.mass , 'getMin' ) and x < self.mass.getMin() : return 0.
        if hasattr ( self.mass , 'getMax' ) and x > self.mass.getMax() : return 0.

        ## make a try to use analytical derivatives 
        if hasattr ( self , 'pdf' ) :
            _pdf = self.pdf 
            if hasattr ( _pdf , 'setPars'  ) : _pdf.setPars() 
            try: 
                if hasattr ( _pdf , 'function' ) :
                    _func = _pdf.function() 
                    if hasattr ( _func , 'integral' ) :
                        return _func.derivative ( x )
            except:
                pass
            
        ## use numerical ingeration
        from ostap.math.derivative import derivative as _derivatve
        return _derivative ( self , x )


# =============================================================================
##  helper utilities to imlement resolution models.
# =============================================================================
class _CHECKMEAN(object) :
    check = True
def checkMean() :
    return True if  _CHECKMEAN.check else False
class Resolution(object) :    
    def __init__  ( self , resolution = True ) :
        self.check = False if resolution else True 
    def __enter__ ( self ) :
        self.old         = _CHECKMEAN.check 
        _CHECKMEAN.check =  self.check
    def __exit__  ( self , *_ ) :
        _CHECKMEAN.check =  self.old 
# =============================================================================
## helper base class for implementation  of various helper pdfs 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class MASS(PDF) :
    """Helper base class for implementation of various 1D-pdfs
    (presumably models for mass distributions)
    - it adds more functionality  to class PDF
    - 
    """
    def __init__ ( self            ,
                   name            ,
                   mass     = None ,
                   mean     = None ,
                   sigma    = None ) : 


        ##  if not mn is None : logger.warning('Ignore mn-argument')
        ## if not mx is None : logger.warning('Ignore mx-argument')
        
        m_name  = "m_%s"     % name
        m_title = "mass(%s)" % name

        if   isinstance ( mass , ROOT.TH1   ) :
            mass    = mass.xminmax()
            m_title = name.GetTitle()            
        elif isinstance ( mass , ROOT.TAxis ) :
            mass    = mass.GetXmin() , mass.GetXmax()

        ## create the variable 
        if isinstance ( mass , tuple ) and 2 == len(mass) :  
            mass = makeVar ( mass    , ## var 
                             m_name  , ## name 
                             m_title , ## title/comment
                             *mass   , ## min/max 
                             fix = None ) ## fix ? 
        elif isinstance ( mass , ROOT.RooAbsReal ) :
            mass = makeVar ( mass    , ## var 
                             m_name  , ## name 
                             m_title , ## title/comment
                             fix = None  ) ## fix ? 
        else :
            raise AttributeError("Unknown type of ``mass'' parameter %s/%s" % ( type ( mass ) , mass ) ) 

        ## intialize the base 
        PDF.__init__ ( self , name , xvar = mass )
        
        
        self._mn = self.mass.getMin ()
        self._mx = self.mass.getMax ()
        #
        _dm      = self._mx - self._mn 
        
        #
        ## mean-value
        # 
        self.__mean = makeVar ( mean              ,
                                "mean_%s"  % name ,
                                "mean(%s)" % name , mean ,  self._mn  , self._mx )
        ## 
        if checkMean () :
            
            if self.mean.isConstant() :
                if not self._mn <= self.mean.getVal() <= self._mx :
                    raise AttributeError ( 'MASS(%s): Fixed mass %s is not in mass-range (%s,%s)' % ( name , self.mean.getVal() , self._mn , self._mx ) )
            elif hasattr ( self.mean , 'setMin' ) and hasattr( self.mean , 'setMax' ) : 
                self.mean.setMin ( max ( self.mean.getMin () , self.mass.getMin() - 0.1 * _dm ) )
                self.mean.setMax ( min ( self.mean.getMax () , self.mass.getMax() + 0.1 * _dm ) )
                logger.debug ( 'MASS(%s) Mean range is redefined to be (%s,%s)' % ( name , self.mean.getMin() , self.mean.getMax() ) ) 
                    
        #
        ## sigma
        #
        sigma_max  = 2.0 * _dm / math.sqrt ( 12 )
        self.__sigma = makeVar ( sigma              ,
                                 "sigma_%s"  % name ,
                                 "sigma(%s)" % name , sigma , 0.01 * sigma_max , 0 , sigma_max )        
    
    @property 
    def mass ( self ) :
        """``mass''-variable for the fit (the same as ``x'' or ``xvar'')"""
        return self.xvar
    
    @property
    def mean ( self ):
        """This is my doc for property mean"""
        return self.__mean
    @mean.setter
    def mean ( self , value ) :
        self.mean.setVal ( value )
        return self.mean.getVal()
        
    @property
    def sigma ( self ):
        """This is doc for property sigma"""
        return self.__sigma
    @sigma.setter
    def sigma ( self , value ) :
        self.sigma.setVal ( value )
        return self.sigma.getVal()
    
    
# =============================================================================
## @class RESOLUTION
#  helper base class  to parameterize the resolution
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2017-07-13
class RESOLUTION(MASS) :
    """Helper base class  to parameterize the resolution
    """
    def __init__ ( self            ,
                   name            ,
                   mass     = None ,
                   sigma    = None , 
                   mean     = None ) : 
        from ostap.fitting.basic import Resolution
        with Resolution() :
            super(RESOLUTION,self).__init__ ( name  = name  ,
                                              mass  = mass  ,
                                              sigma = sigma ,
                                              mean  = mean  )

# =============================================================================
## simple convertor of 1D-histo to data set
#  @code
#  h1   = ...
#  dset = H1D_dset ( h1 )
#  @endcode 
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H1D_dset(object) :
    """Simple convertor of 1D-histogram into data set
    >>> h1   = ...
    >>> dset = H1D_dset ( h1 )
    """
    def __init__ ( self            , 
                   histo           ,
                   mass    = None  ,
                   density = True  ,
                   silent  = False ) :
        #
        ## use mass-variable
        #
        name = histo.GetName() 
        if mass : self.__haxis = mass 
        else    : self.__haxis = makeVar ( mass , 'm_%s' % name , 'mass(%s)' % name , None , *(histo.xminmax()) )

        self.impDens = density 
        
        with roo_silent ( silent ) :  
            
            self.var1    = self.mass
            self.x       = self.mass
            self.vlst    = ROOT.RooArgList    ( self.haxis )
            self.vimp    = ROOT.RooFit.Import ( histo     , density )
            self.dset    = ROOT.RooDataHist   (
                rootID ( 'hds_' ) ,
                "Data set for histogram '%s'" % histo.GetTitle() ,
                self.vlst  ,
                self.vimp  )
    @property
    def haxis  (self ):
        """The histogram axis"""
        return self.__haxis 
    
# =============================================================================
## simple convertor of 1D-histogram into PDF
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class H1D_pdf(H1D_dset,PDF) :
    """Simple convertor of 1D-histogram into PDF
    """
    def __init__ ( self            ,
                   name            ,
                   histo           ,
                   mass    = None  ,
                   density = True  ,
                   silent  = False ) :
        
        H1D_dset.__init__ ( self , histo , mass , density , silent )
        PDF     .__init__ ( self , name  , self.haxis )
        
        with roo_silent ( silent ) : 
            #
            ## finally create PDF :
            self.vset  = ROOT.RooArgSet  ( self.haxis )        
            self.pdf   = ROOT.RooHistPdf (
                'hpdf_%s'            % name ,
                'HistoPDF(%s/%s/%s)' % ( name , histo.GetName() , histo.GetTitle() ) , 
                self.vset  , 
                self.dset  )
            
        ## and declare it be be a "signal"
        self.signals().add ( self.pdf ) 
        

# =============================================================================
## @class Fit1D
#  The actual model for 1D ``mass-like'' fits
#  @code
#  ## signal components 
#  signal       = ... ## MUST be "Ostap"-object 
#  signal_2     = ... ## can be Ostap-object or bare RooFit pdf 
#  signal_3     = ... ## ditto
#  ...
#  signal_K     = ... ## ditto
#  ## background components 
#  background   = ... ## can be Ostap-object or bare RooFit pdf or just an integer number 
#  background_1 = ... ## can be Ostap-object or bare RooFit pdf
#  background_2 = ... ## can be Ostap-object or bare RooFit pdf
#  ...
#  background_M = ... ## can be Ostap-object or bare RooFit pdf
#  ## ``other'' components 
#  component_1  = ... ## can be Ostap-object or bare RooFit pdf 
#  component_2  = ... ## can be Ostap-object or bare RooFit pdf
#  ...
#  component_N  = ... ## can be Ostap-object or bare RooFit pdf
#  model = Fit1D (
#    signal           = signal     ,
#    othersignals     = [ signal_1 , signal_2 , ... , signal_K ] ,
#    background       = background ,
#    otherbackgrounds = [ background_1 , background_2 , ... , background_M ] ,
#    others           = [ component_1  , component_2  , ... , component_N  ] ,
#    ...
#  )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-08-02
class Fit1D (PDF) :
    """The actual fit-model for generic 1D ``mass-like'' fits
    ## 1) signal components 
    >>> signal       = ... ## main signal: it MUST be Ostap-object 
    >>> signal_2     = ... ## can be Ostap-object or bare RooFit pdf 
    >>> signal_3     = ... ## ditto
    ...
    >>> signal_K     = ... ## ditto
    ## 2) background components 
    >>> background   = ... ## main ``background'': it can be Ostap-object or bare RooFit pdf or just an integer number 
    >>> background_1 = ... ## can be Ostap-object or bare RooFit pdf
    >>> background_2 = ... ## can be Ostap-object or bare RooFit pdf
    ...
    >>> background_M = ... ## can be Ostap-object or bare RooFit pdf
    ## 3) ``other'' components (optional) 
    >>> component_1  = ... ## can be Ostap-object or bare RooFit pdf 
    >>> component_2  = ... ## can be Ostap-object or bare RooFit pdf
    ...
    >>> component_N  = ... ## can be Ostap-object or bare RooFit pdf
    ## build a model:
    model = Fit1D (
    ... signal           = signal     ,
    ... othersignals     = [ signal_1 , signal_2 , ... , signal_K ] ,
    ... background       = background ,
    ... otherbackgrounds = [ background_1 , background_2 , ... , background_M ] ,
    ... others           = [ component_1  , component_2  , ... , component_N  ] ,
    )
    """
    def __init__ ( self                          , 
                   signal                        ,  ## the main signal 
                   background          = None    ,  ## the main background 
                   othersignals        = []      ,  ## additional signal         components
                   otherbackgrounds    = []      ,  ## additional background     components
                   others              = []      ,  ## additional non-classified components 
                   suffix              = ''      ,  ## the suffix 
                   name                = ''      ,  ## the name 
                   extended            = True    ,  ## extended fits ?
                   combine_signals     = False   ,  ## combine signal PDFs into single "SIGNAL"     ? 
                   combine_backgrounds = False   ,  ## combine signal PDFs into single "BACKGROUND" ?            
                   combine_others      = False   ,  ## combine signal PDFs into single "COMPONENT"  ?             
                   recursive           = True    ): ## recursive fractions for NON-extended models?

        #
        if not name : 
            if   signal.name and suffix : name = signal.name + suffix
            elif signal.name            : name = signal.name + '1D' 
            #
            
        #
        self.extended   = extended 
        self.suffix     = suffix 
        self.signal     =      signal

        
        # Init base class 
        PDF.__init__ ( self , name + suffix , signal.xvar )

        self.mass       = self.xvar 

        #
        self.background = makeBkg ( background , 'Background' + suffix , self.mass )
        #
        ## update the lists of PDFs
        #
        self.signals     ().add ( self.signal     .pdf )
        self.backgrounds ().add ( self.background .pdf )

        self.more_signals       = othersignals
        self.more_backgrounds   = otherbackgrounds
        self.more_components    = others
        #
        ## treat additional signals
        #
        for c in self.more_signals : 
            if   isinstance ( c , ROOT.RooAbsPdf ) : self.signals() .add ( c     ) 
            elif hasattr    ( c , 'pdf' )          : self.signals() .add ( c.pdf ) 
            else : logger.error('Fit1D(%s): Unknown signal component type %s, skip it!' % ( self.name , type(c) ) ) 
        #
        ## treat additional backgounds 
        #
        for c in self.more_backgrounds :             
            if   isinstance ( c ,  ROOT.RooAbsPdf ) : self.backgrounds () .add ( c     ) 
            elif hasattr    ( c ,'pdf' )            : self.backgrounds () .add ( c.pdf ) 
            else : logger.error('Fit1D(%s): Unknown background component type %s, skip it!' % ( self.name , type(c) ) ) 
        #
        ## treat additional components
        #
        for c in self.more_components : 
            if   isinstance ( c ,  ROOT.RooAbsPdf ) : self.components () .add ( c     ) 
            elif hasattr    ( c ,'pdf' )            : self.components () .add ( c.pdf ) 
            else : logger.error('Fit1D(%s): Unknown additional component type %s, skip it!' % ( self.name , type(c) ) ) 

        #
        ## build PDF
        #
        
        ## all fit components 
        self.alist1 = ROOT.RooArgList ()

        ## all fit fractions (or yields for extended fits) 
        self.alist2 = ROOT.RooArgList ()
        
        
        self.all_signals     = self.signals     ()
        self.all_backgrounds = self.backgrounds ()
        self.all_components  = self.components  ()
        
        self.save_signal     = self.signal
        self.save_background = self.background
        
        ## combine signal components into single signal  (if needed) 
        if combine_signals and 1 < len( self.signals() ) :
            
            sig , fracs , sigs = addPdf ( self.signals()        ,
                                          'signal_'    + suffix ,
                                          'signal(%s)' % suffix ,
                                          'fS_%%d%s'   % suffix ,
                                          'fS(%%d)%s'  % suffix , recursive = True , model = self )
            ## new signal
            self.signal      = Generic1D_pdf   ( sig , self.mass , 'SIGNAL_' + suffix )
            self.all_signals = ROOT.RooArgList ( sig )
            self._sigs       = sigs 
            self.signals_fractions = fracs 
            for fi in fracs : setattr ( self , fi.GetName() , fi )
            logger.verbose('Fit1D(%s): %2d signals     are combined into single SIGNAL'     % ( self.name , len ( sigs ) ) ) 
            
        ## combine background components into single backhround (if needed ) 
        if combine_backgrounds and 1 < len( self.backgrounds() ) :
            
            bkg , fracs , bkgs = addPdf ( self.backgrounds()        ,
                                          'background_'    + suffix ,
                                          'background(%s)' % suffix ,
                                          'fB_%%d%s'       % suffix ,
                                          'fB(%%d)%s'      % suffix , recursive = True , model = self )
            ## new background
            self.background      = Generic1D_pdf   ( bkg , self.mass , 'BACKGROUND_' + suffix )
            self.all_backgrounds = ROOT.RooArgList ( bkg )
            self._bkgs           = bkgs 
            self.backgrounds_fractions = fracs 
            for fi in fracs : setattr ( self , fi.GetName() , fi ) 
            logger.verbose ('Fit1D(%s): %2d backgrounds are combined into single BACKGROUND' % ( self.name , len ( bkgs ) ) ) 

        ## combine other components into single component (if needed ) 
        if combine_others and 1 < len( self.components() ) :
            
            cmp , fracs , cmps = addPdf ( self.components()    ,
                                          'other_'    + suffix ,
                                          'other(%s)' % suffix ,
                                          'fC_%%d%s'  % suffix ,
                                          'fC(%%d)%s' % suffix , recursive = True , model = self )
            ## save old background
            self.other          = Generic1D_pdf   ( cmp , self.mass , 'COMPONENT_' + suffix )
            self.all_components = ROOT.RooArgList ( cmp )
            self.components_fractions = fracs 
            for fi in fracs : setattr ( self , fi.GetName() , fi ) 
            logger.verbose('Fit1D(%s): %2d components  are combined into single COMPONENT'    % ( self.name , len ( cmps ) ) )


        self.nums   = [] 

        ## build models 
        if extended :
            
            ns = len ( self.all_signals )
            if 1 == ns :
                sf = makeVar ( None , "S"+suffix , "Signal"     + suffix , None , 1 , 0 , 1.e+6 )
                self.alist1.add  ( self.all_signals[0]  )
                self.alist2.add  ( sf ) 
                self.s = sf
                self.S = sf
                self.S_name = self.s.GetName()
                self.nums.append ( self.s ) 
            elif 2 <= ns : 
                fis = makeFracs ( ns , 'S_%%d%s' % suffix ,  'S(%%d)%s'  % suffix , fractions  = False , model = self )
                for s in self.all_signals : self.alist1.add ( s )
                for f in fis :
                    self.alist2.add  ( f )
                    self.nums.append ( f ) 

            nb = len ( self.all_backgrounds )
            if 1 == nb :
                bf = makeVar ( None , "B"+suffix , "Background" + suffix , None , 1 , 0 , 1.e+6 )
                self.alist1.add  ( self.all_backgrounds[0]  )
                self.alist2.add  ( bf ) 
                self.b = bf
                self.B = bf
                self.B_name = self.b.GetName()
                self.nums.append ( self.b ) 
            elif 2 <= nb :
                fib = makeFracs ( nb , 'B_%%d%s' % suffix ,  'B(%%d)%s'  % suffix , fractions  = False , model = self )
                for b in self.all_backgrounds : self.alist1.add ( b )
                for f in fib :
                    self.alist2.add  ( f )
                    self.nums.append ( f ) 
                    
            nc = len ( self.all_components )
            if 1 == nc :
                cf = makeVar ( None , "C"+suffix , "Component" + suffix , None , 1 , 0 , 1.e+6 )
                self.alist1.add  ( self.all_components[0]  )
                self.alist2.add  ( cf ) 
                self.c = cf
                self.C = cf
                self.C_name = self.c.GetName()
                self.nums.append ( self.c ) 
            elif 2 <= nc : 
                fic = makeFracs ( nc , 'C_%%d%s' % suffix ,  'C(%%d)%s'  % suffix , fractions  = False , model = self )
                for c in self.all_components : self.alist1.add ( c )
                for f in fic :
                    self.alist2.add  ( f )
                    self.nums.append ( f ) 

        else :

            ns = len ( self.all_signals     )
            nb = len ( self.all_backgrounds )
            nc = len ( self.all_components  )
            
            for s in self.all_signals     : self.alist1.add ( s )
            for b in self.all_backgrounds : self.alist1.add ( b )
            for c in self.all_components  : self.alist1.add ( c )
            
            fic = makeFracs ( ns + nb + nc , 'f_%%d%s' % suffix ,  'f(%%d)%s'  % suffix , fractions  = True , model = self )
            for f in fic :
                self.alist2.add ( f )
                setattr ( self , f.GetName() , f ) 

        #
        ## The final PDF
        #       
        if   self.name       : title = "model(%s)"     % self.name
        else                 : title = "model(Fit1D)"
        
        name     = name if name else ( 'model_' + self.name )
        
        pdfargs  = name , title , self.alist1 , self.alist2
        if not extended : pdfargs = pdfargs + ( True if recursive else False , ) ## RECURSIVE ? 
        self.pdf = ROOT.RooAddPdf ( *pdfargs )
        
        if extended : 
            logger.debug ( "Extended     model ``%s'' with %s/%s components"  % ( self.pdf.GetName() , len( self.alist1) , len(self.alist2) ) )
        else : 
            logger.debug ( "Non-extended model ``%s'' with %s/%s components"  % ( self.pdf.GetName() , len( self.alist1) , len(self.alist2) ) )



        
# =============================================================================
## simple class to adjust certaint PDF to avoid zeroes 
class Adjust(object) :
    """
    Simple class to adjust certain PDF to avoid zeroes 
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,
                   mass             , 
                   pdf              ,
                   value    = 1.e-3 ) : 

        self.name    = name 
        self.mass    = mass
        self.old_pdf = pdf
        
        self.p0_pdf  = ROOT.RooPolynomial( 'p0_%s'     % name ,
                                           'poly0(%s)' % name , self.mass ) 

        
        self.num_f   = makeVar ( None  , 'valueT_%s'   % name ,
                                 'value/true(%s)'      % name ,
                                 None  ,
                                 value , 
                                 0     , 1 )
        
        self.alist1 = ROOT.RooArgList (
            self.p0_pdf  ,   
            self.old_pdf 
            )
        
        self.alist2 = ROOT.RooArgList (
            self.num_f    ,
            )
        #
        ## final PDF
        # 
        self.pdf  = ROOT.RooAddPdf  ( "adjust_"    + name ,
                                      "Adjust(%s)" % name ,
                                      self.alist1 ,
                                      self.alist2 )
        
# =============================================================================
## @class Convolution
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-07-13
class Convolution(object):
    """
    Helper class to make a convolution PDF 
    """
    def __init__ ( self           ,
                   name           ,
                   pdf            ,
                   mass           , 
                   convolution    ,
                   useFFT  = True ) :
        
        ## store "old" pdf 
        self.orig_pdf = pdf
        self.mass     = mass

        if   isinstance ( convolution ,   ROOT.RooAbsPdf       ) : self.convolution = convolution
        elif isinstance ( convolution , ( float , int , long ) ) and not isinstance ( convolution , bool ) :
            self.cnv_mean  = makeVar (
                0.0  ,
                'CnvMean'       + name ,
                'cnv_mean (%s)' % name , 
                0.0  , 0 ) 
            self.cnv_sigma = makeVar (
                None ,
                'CnvSigma'      + name ,
                'cnv_sigma(%s)' % name ,
                convolution ,
                convolution ,
                0.05 * convolution , 10 * convolution )
            self.cnv_gauss = ROOT.RooGaussian (
                'CnvGauss'     + name , 
                'CnvGauss(%s)' % name ,
                self.mass , self.cnv_mean , self.cnv_sigma )
            self.convolution = self.cnv_gauss
        else :
            raise AttributeError ( 'Unknown convolution type %s ' % convolution )

        #
        if useFFT :
            
            nb = 40000
            if hasattr ( self , 'cnv_sigma' ) :
                dm  = mass.getMax() - mass.getMin()
                dm /= self.cnv_sigma
                nb  = max ( nb , 100 * int (  dm ) ) 
                logger.debug('Convolution: choose #bins %d' % nb )
                
            self.mass.setBins ( nb , 'cache' ) 
            self.pdf = ROOT.RooFFTConvPdf ( 'FFT'     + name ,
                                            'FFT(%s)' % name ,
                                            self.mass , self.orig_pdf , self.convolution )
            self.pdf.setBufferFraction ( 0.25 )
            
        else      :
            self.pdf = ROOT.RooNumConvPdf ( 'CNV'     + name ,
                                            'CNV(%s)' % name ,
                                            self.mass , self.orig_pdf , self.convolution )
            if isinstance ( self.convolution , ROOT.RooGaussian ) :
                if hasattr ( self , 'cnv_mean' ) and hasattr ( self , 'cnv_sigma' ) :
                    self.pdf.setConvolutonWindow ( self.cnv_mean , self.cnv_sigma , 7 )


# =============================================================================
## helper class to build/keep the list of ``phi''-arguments (needed e.g. for polynomial functions)
class Phases(object) :
    """Helper class to build/keep the list of ``phi''-arguments (needed e.g. for polynomial functions)
    """
    ## Create vector of phases (needed for various polynomial forms)
    def __init__( self  , num , the_phis = None ) :
        """Create vector of phases (needed for various polynomial forms)
        """
        if the_phis :
            self.__phis     = [ i for i in the_phis.phis ]  
            self.__phi_list = the_phis.phi_list
            return
        
        self.__phis     = []
        self.__phi_list = ROOT.RooArgList()
        from math import pi 
        for i in range( 1 , num + 1 ) :
            phi_i = makeVar ( None ,
                              'phi%d_%s'      % ( i , self.name )  ,
                              '#phi_{%d}(%s)' % ( i , self.name )  ,
                              None , 0 ,  -0.75 * pi  , 1.25 * pi  )
            self.__phis.append  ( phi_i ) 
            self.phi_list.add ( phi_i )

    ## set all phis to be 0
    def reset_phis ( self ) :
        """Set all phases to be zero
        """
        for f in self.__phis : f.setVal(0)
        
    @property
    def phis ( self ) :
        """The list of ``phases'', used to parameterize various polynomial-like shapes
        """
        return tuple ( [ i for i in self.__phis ] )

    @property
    def phi_list ( self ) :
        """The list/ROOT.RooArgList of ``phases'', used to parameterize polynomial-like shapes
        """
        return self.__phi_list
    

# =============================================================================
## @class Generic1D_pdf
#  "Wrapper" over generic RooFit (1D)-pdf
#  @code
#     
#  raw_pdf = RooGaussian  ( ...     )
#  pdf     = Generic1D_pdf ( raw_pdf , varx = x )  
# 
#  @endcode 
#  If more functionality is required , more actions are possible:
#  @code
#  ## for sPlot 
#  pdf.alist2 = ROOT.RooArgList ( n1 , n2 , n3 ) ## for sPlotting 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-03-29
class Generic1D_pdf(PDF) :
    """Wrapper for generic RooFit pdf
    
    >>> raw_pdf = RooGaussian   ( ...     )
    >>> pdf     = Generic1D_pdf ( raw_pdf , varx = x )
    """
    ## constructor 
    def __init__ ( self , pdf , xvar , name = None ) :
        if not name : name = pdf.GetName()
        PDF  . __init__ ( self , name , xvar  )
        self.pdf  = pdf
        self.signals().add ( self.pdf )

# =============================================================================
## create simple background model
#  As basic the default - product of exponential function and  positive polynomial
#  @see  Bkg_pdf 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-04-03
def makeBkg ( bkg , name , xvar , **kwargs ) :
    """Helper function to create background models (around Bkg_pdf)
    
    >>> x =   .. ## the variable
    
    ## non-negative integer, construct PDF using Bkg_pdf 
    >>> bkg1  = makeBkg ( 3      , 'B1' , x ) ## use Bkg_pdf ( 'B1' , x , power = 3 )
    
    ## generic RooAbsPdf: use this PDF
    >>> pdf   = RooPolynomial( ... )
    >>> bkg2  = makeBkg ( pdf    , 'B2' , x ) ## use Generic1D_pdf ( pdf , x , 'B2' )
    
    ## some Ostap-based model, use it as it is  
    >>> model = Convex_pdf ( ... )
    >>> bkg3  = makeBkg ( models , 'B3' , x )  
    
    ## some RooAbsReal, use is as exponenial slope for Bkg_pdf
    >>> tau  = RooRealVar( ....  )
    >>> bkg4 = makeBkg ( tau     , 'B4' , x ) 
    
    """

    ## a  kind of  default ...
    if bkg is None : bkg = 0
    
    ## Regular case: degree of polynom in Bkg_pdf 
    if isinstance ( bkg , ( int , long ) ) and 0 <= bkg :

        import ostap.fitting.background as OFB 
        model = OFB.Bkg_pdf ( name , power = bkg , mass = xvar , **kwargs )
        return model
    
    ## native RooFit pdf ? 
    elif isinstance ( bkg , ROOT.RooAbsPdf ) :
        
        model = Generic1D_pdf ( bkg , varx = xvar , name = name ) 
        return model
    
    ## some Ostap-based model ?
    elif hasattr    ( bkg , 'pdf' ) and isinstance ( bkg.pdf , ROOT.RooAbsPdf ) :

        if hasattr ( bkg , 'mass' ) and not bkg.mass is xvar :
            logger.warning('Something weird is here') 
            try :
                logger.info ('make   try to create model') 
                model = bkg.__class__  ( name = name ,
                                         mass = mass , **kwargs ) 
                logger.info ('Model is  created!')
                return model 
            except:
                logger.warning('failure to create model') 
                
        model = bkg 
        return model

    ## interprete it as exponential slope for Bkg-pdf 
    elif isinstance ( bkg , ROOT.RooAbsReal ) :

        import ostap.fitting.background as OFB 
        model = OFB.Bkg_pdf ( name , mass = xvar , tau = bkg , **kwargs )
        return model
    
    raise  TypeError("Wrong type of bkg object: %s/%s " % ( bkg , type(bkg) ) )

        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
