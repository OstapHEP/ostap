#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/rooreduce.py
#  Reduction/serialization/deserialization for some Ostap/ROOT classes  
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
""" Reduction/serialization/deserialization for some Ostap/ROOT classes
- see RooAbsReal
- see RooRealVar
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    ) 
# =============================================================================
from   builtins                     import range
from   ostap.core.meta_info         import root_info 
from   ostap.math.base              import doubles
from   ostap.core.core              import Ostap
from   ostap.math.reduce            import root_factory
import ostap.fitting.variables
import ostap.fitting.roocollections
import ostap.histos.graph_reduce    as     GR 
import ROOT, random, array, ctypes, math, itertools    
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.rooreduce' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug( 'Some useful decorations for RooFit variables')
# =============================================================================
_new_methods_ = []
# =============================================================================
## Factory for deserialization of generic objects
#  @attention it stores the constructor kaprameters as local attributes
def root_store_factory ( klass , *params ) :
    """Factory for deserialization of generic object
    - attention: it stores the constructor kaprameters as local attributes
    """
    ## create the objects 
    obj = root_factory ( klass , *params )
    ## keep argumets with the newly created obnject  
    obj.__store = params    ## Attention - keep argumetns with newly crfeated object!
    return obj 


# =============================================================================
## Some native ROOT/RooFit objects 
# =============================================================================

    
# =============================================================================
## Dedicated unpickling factory for RooRealVar
def _rrv_factory ( args , errors , binnings , fixed , *attrs ) :
    """ Dedicated unpickling factory for `ROOT.RooRealVar`
    """

    ## create it
    rrv = ROOT.RooRealVar ( *args )

    ## set errors if needed 
    if   errors and  2 == len ( errors ) : rrv.setAsymError ( *errors ) 
    elif errors and  1 == len ( errors ) : rrv.setError     ( *errors )

    for b in binnings :
        if b : rrv.setBinning ( b , b.GetName ()) 

    rrv.setConstant ( fixed )
    
    if attrs :
        battrs = attrs[0] 
        for n , a in battrs : rrv.setAttribute          ( n , a )
        if 1 < len ( attrs ) :
            sattrs = attrs[1] 
            for n , a in sattrs : rrv.setStringAttribute    ( n , a )
            if 2 < len ( attrs ) :
                tattrs = attrs[2]             
                for n , a in tattrs : rrv.setTransientAttribute ( n , a ) 
    
    return rrv

# =============================================================================
## Reducing of <code>RooRealVar</code> for pickling/unpickling 
#  @see RooRooRealVar 
def _rrv_reduce ( rrv ) :
    """ Reducing of `ROOT.RooRealVar` for pickling 
    - see ROOT.RooRooRealVar 
    """
    name    = rrv.name 
    title   = rrv.title
    value   = rrv.getVal () 

    has_min = rrv.hasMin ()
    has_max = rrv.hasMax ()

    ## constructor arguments 
    if has_min and has_max :
        args = name , title , value ,  rrv.getMin () , rrv.getMax ()              , rrv.getUnit ()
    elif has_min : 
        args = name , title , value ,  rrv.getMin () , ROOT.RooNumber.infinity () , rrv.getUnit () 
    elif has_max : 
        args = name , title , value , -ROOT.RooNumber.infinity () , rrv.getMax () , rrv.getUnit () 
    else :
        args = name , title , value ,  rrv.getUnit () 

    ## errors 
    if   rrv.hasAsymError () :
        errors = rrv.getAsymErrorLo () , rrv.getAsymErrorHi ()
    elif rrv.hasError   () :
        errors = rrv.getError() ,
    else :
        errors = () 

    ## binings 
        
    binnings = tuple (   rrv.getBinning ( n , False )       for n in rrv.getBinningNames     () )

    ## attributes:
    
    battrs   = tuple ( ( n , rrv.getAttribute          ( n ) ) for   n       in rrv.attributes          () ) 
    sattrs   = tuple ( ( n , a                               ) for ( n , a ) in rrv.stringAttributes    () ) 
    tattrs   = tuple ( ( n , rrv.getTransientAttribute ( n ) ) for   n       in rrv.transientAttributes () ) 

    ## fixed ? 
    fixed = True if rrv.isConstant() else False 

    if tattrs : 
        content = args , errors , binnings , fixed , battrs , sattrs , tattrs
    elif sattrs :
        content = args , errors , binnings , fixed , battrs , sattrs
    elif battrs : 
        content = args , errors , binnings , fixed , battrs 
    else :        
        content = args , errors , binnings , fixed 
    
    return _rrv_factory , content 

ROOT.RooRealVar.__reduce__ = _rrv_reduce


# =============================================================================
## Reduce <code>RooConstVar</code>
#  @see RooConstVar 
def _rconst_reduce ( var ) :
    """ Reduce `ROOT.RooConstVar`
    - see ROOT.RooConstVar
    """
    return root_factory , ( type ( var )  ,
                            var.name      ,
                            var.title     ,
                            float ( var ) )  

ROOT.RooConstVar.__reduce__ = _rconst_reduce 

# ==============================================================================
## reduce <code>RooLinearVar</code>
#  @see RooLinearVar
def _rlinv_reduce ( var ) :
    """Reduce `RooLinearVar`
    - see `ROOT.RooLinearVar`
    """ 
    return root_factory , ( type ( var )  ,
                            var.name      ,
                            var.title     ,
                            var.variable  ,
                            var.slope     ,
                            var.offset    ,
                            var.unit      )
    
ROOT.RooLinearVar.__reduce__ = _rlinv_reduce 

# ==============================================================================
## factory for RooCategory objects
#  @see RooCategory
def _rcat_factory_  ( klass , name , title , index  , items ) :
    """factory for `ROOT.RooCategory` objects
    - see `ROOT.RooCategory`
    """
    cat = klass ( name , title )
    for label , index  in items : cat.defineType ( label , index )
    cat.setIndex ( index ) 
    return cat

# =============================================================================
## reduce RooCategory instance 
#  @see RooCategory
def _rcat_reduce_ ( cat ) :
    """reduce RooCategory instance 
    - see `ROOT.RooCategory`
    """
    items = tuple ( (l,i) for (l,i) in cat.items() )
    return _rcat_factory_ , ( type ( cat ) , cat.GetName() , cat.GetTitle() , cat.getIndex() , items ) 

ROOT.RooCategory   .__reduce__    = _rcat_reduce_ 

# =============================================================================
## factory for unpickling of <code>RooFormulaVar</code> and
#  <code>Ostap::FormulaVar</code>
#  @see RooFormualVar 
#  @see Ostap::FormualVar 
def _rfv_factory ( klass , args , vars ) :
    """Factory for unpickling of `RooFormulaVar` and `Ostap.FormulaVar`
    - see ROOT.RooFormulaVar 
    - see Ostap.FormualVar 
    """

    lst = ROOT.RooArgList ()
    for v in vars : lst.add ( v ) 

    margs = list  ( args  )
    margs.append  ( lst   )
    margs = tuple ( margs ) 
    
    rfv = klass   ( *margs )
    
    rfv.__args = args, vars, lst

    return rfv
    
# =============================================================================
## Reduce <code>RooFormulaVar</code> and <code>Ostap::FormulaVar</code> for pickling
#  @see RooFormulaVar 
#  @see Ostap::FormulaVar 
def _rfv_reduce ( rfv ) : 
    """Reduce `RooFormulaVar` and `Ostap::FormulaVar` for pickling
    - see RooFormulaVar 
    - see Ostap.FormulaVar 
    """

    name       = rfv.GetName     ()
    title      = rfv.GetTitle    ()
    
    rform      = rfv.formula     ()    
    expression = rform.GetTitle  () 

    vars       = tuple ( d for d in rform.actualDependents() ) 
    args       = name , title , expression
    
    return _rfv_factory , ( type ( rfv ) , args , vars ) 


if (6,22) <= root_info : ROOT.RooFormulaVar.__reduce__  = _rfv_reduce
else                   : Ostap.FormulaVar.__reduce__    = _rfv_reduce

# =============================================================================
## unpickle RooUniformBinning object
#  @see RooUniformBinnig  
def _rub_factory ( *args ) :
    """unpickle RooUniformBinning object
    -see ROOT.RooUniformBinnig
    """
    return ROOT.RooUniformBinning ( *args )

# =============================================================================
## reduce uniform binning scheme
#  @see RoUniformBinnig 
def _rub_reduce_ ( rub ) :
    """Reduce RooUniformBininkg Object
    - see ROOT.RooUniformBinning
    """
    nbins = rub.numBoundaries()
    if nbins : nbins -= 1
    content = rub.lowBound () , rub.highBound(), nbins , rub.GetName()
    
    return _rub_factory, content

# =============================================================================
## unpickle RooBinning object
#  @see RooBinnig  
def _rb_factory ( data , name  ) :
    """unpickle RooBinning object
    -see ROOT.RooUniformBinnig
    """
    return ROOT.RooBinning ( len ( data ) - 1 , data [ 0 ] , name )
# =============================================================================
## reduce RooBinning object
#  @see RooBinning 
def _rb_reduce_ ( rb  )  :
    """Reduce RooBinning object
    - see ROOT.RooBinning 
    """
    if rb.isUniform() : return _rub_reduce_ ( rb )

    nb    = rb.numBoundaries() 
    ab    = rb.array ()
    data  = array.array ( 'd' , [ 1.0 * ab[i] for i in range ( nb ) ] )

    content = data, rb.GetName()  
    return _rb_factory, content

# ==========================================================================
## unpickle RooRangeBinning object
#  @see RooRangeBinnig 
def _rrb_factory ( low , high , name ) :
    """unpickle RooRangeBinning object"""
    return ROOT.RooRangeBinning( low , high , name )
# ============================================================================
## reduce RooRangeBinnig object
#  @see RooRangeBinnig 
def _rrb_reduce_ ( rrb ) :
    """Reduce RooRangeBinnig object"""
    return _rrb_factory ,  ( rrb.lowBound() , rrb.highBound() , rrb.GetName() ) 

ROOT.RooBinning       .__reduce__ = _rb_reduce_
ROOT.RooUniformBinning.__reduce__ = _rub_reduce_
ROOT.RooRangeBinning  .__reduce__ = _rrb_reduce_
    
if not hasattr ( ROOT.RooGaussian , 'getX' ) :
    def _rgau_x_ ( pdf ) :
        """Get x-observable"""
        return Ostap.MoreRooFit.getX ( pdf )
    ROOT.RooGaussian.getX = _rgau_x_
    _new_methods_ += [ ROOT.RooGaussian.getX ]

if not hasattr ( ROOT.RooGaussian , 'getMean' ) :
    def _rgau_mean_ ( pdf ) :
        """Get x-observable"""
        return Ostap.MoreRooFit.getMean ( pdf )
    ROOT.RooGaussian.getMean = _rgau_mean_
    _new_methods_ += [ ROOT.RooGaussian.getMean ]

if not hasattr ( ROOT.RooGaussian , 'getSigma' ) :
    def _rgau_sigma_ ( pdf ) :
        """Get sigma"""
        return Ostap.MoreRooFit.getSigma ( pdf )
    ROOT.RooGaussian.getSigma = _rgau_sigma_
    _new_methods_ += [ ROOT.RooGaussian.getSigma ]

# ==========================================================================--
## reduce  RooGaussian object 
def _rgau_reduce_ ( pdf ) :
    """Reduce `RooGaussian` object"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.getX     () ,
                                  pdf.getMean  () , 
                                  pdf.getSigma () )

ROOT.RooGaussian.__reduce__ = _rgau_reduce_ 


if not hasattr ( ROOT.RooMultiVarGaussian , 'observables' ) :
    def _rmvgau_observables_ ( pdf ) :
        """Get observables"""
        return Ostap.MoreRooFit.observables( pdf )
    ROOT.RooMultiVarGaussian.observables = _rmvgau_observables_
    _new_methods_ += [ ROOT.RooMultiVarGaussian.observables ]

if not hasattr ( ROOT.RooMultiVarGaussian , 'mu_vec' ) :
    def _rmvgau_mu_vec_ ( pdf ) :
        """Get mu-vec"""
        return Ostap.MoreRooFit.mu_vec ( pdf )
    ROOT.RooMultiVarGaussian.mu_vec = _rmvgau_mu_vec_
    _new_methods_ += [ ROOT.RooMultiVarGaussian.mu_vec ]

# ================================================================================
## deserialize RooMultiVarGaussian object
#  @see RooMultiVarGaussian
def _rmvgau_factory_ ( klass , *args ) :
    """Deserialize `ROOT.RooMultiVarGaussian` 
    - see `ROOT.RooMultiVarGaussian`
    """

    if 5 != len ( args ) : raise pickle.UnpicklingError("RooMultiVarGaussian: ivalid record length!")

    muvec  = args [3]
    covm   = args [4]

    N = len ( muvec )
    if N * N != len ( covm ) : raise pickle.UnpicklingError("RooMultiVarGaussian: cannot reconstruct covmatrix")
    
    muvec = ROOT.TVectorD    ( N , muvec )
    covm  = ROOT.TMatrixDSym ( N , covm  )
    
    newargs = args [:3]  + ( muvec, covm ) 
    return root_store_factory ( klass , *newargs )  
    
# ==========================================================================--
## reduce  RooMultiVarGaussian object 
def _rmvgau_reduce_ ( pdf ) :
    """Reduce `RooMultiVarGaussian` object"""
    return _rmvgau_factory_ , ( type ( pdf )                                  ,
                                pdf.name                                      ,
                                pdf.title                                     ,
                                pdf.observables      ()                       ,
                                array.array ( 'd' , pdf.mu_vec ()           ) , 
                                array.array ( 'd' , pdf.covarianceMatrix () ) )

ROOT.RooMultiVarGaussian.__reduce__ = _rmvgau_reduce_ 


# ==========================================================================--
## Get the original fractions from the <code>RooAddPdf</code>
#  @code
#  addpdf = ...
#  fractions , recursive = addpdf.orig_fracs() 
#  @endcode 
#  @see RooAddPdf 
def _raddpdf_fractions ( radd ) :
    """ Get the original fractions from the <code>RooAddPdf</code>
    >>> addpdf = ...
    >>> fractions , recursive = addpdf.orig_fracs() 
    """
    rec   = ctypes.c_bool ()
    fracs = Ostap.MoreRooFit.fractions ( radd , rec )
    return fracs, rec.value
    
# ==========================================================================--
ROOT.RooAddPdf  .orig_fracs = _raddpdf_fractions

_new_methods_ += [ ROOT.RooAddPdf  .orig_fracs ] 

# ==========================================================================--
## reduce  RooAddPdf object 
def _raddpdf_reduce_ ( pdf ) :
    """Reduce `RooAddPdf` object"""
    content = type ( pdf ) , pdf.name , pdf.title , pdf.pdfList()
    pars    = pdf.coefList ()
    if 1 <= len ( pars ) :
        content   = content + pdf.orig_fracs () 
    return root_store_factory , content

ROOT.RooAddPdf.__reduce__ = _raddpdf_reduce_ 

# =============================================================================
## Is this <code>RooProdPdf</code> object conditional ?
#  @see RooProdPdf 
def _rprodpdf_cond_ ( pdf ) :
    """Is this `RooProdPdf` object conditional ?
    - see `ROOT.RooProdPdf` 
    """
    for p in pdf.pdfList() :
        vset = pdf.findPdfNSet ( p )
        if not set             : continue
        elif 1 <= len ( vset ) : return True 
    return False

ROOT.RooProdPdf.conditional = _rprodpdf_cond_ 
_new_methods_ += [ ROOT.RooProdPdf.conditional ]

# =============================================================================
## reduce RooProdPdf 
#  @see RooProdPdf
def _rprodpdf_reduce_ ( pdf ) :
    """Reduce `ROOT.RooProdPdf`
    - see `ROOT.RooProdPdf`
    """
    if pdf.conditional () :
        import pickle
        raise pickle.PicklingError("RooProdPdf is conditional (cannot be picked)")
    
    content = type ( pdf ) , pdf.name , pdf.title , pdf.pdfList()
    return root_store_factory , content

ROOT.RooProdPdf.__reduce__ = _rprodpdf_reduce_ 

# =============================================================================
## Factory for RooFFTConfPdf 
#  @see RooFFTConfPdf 
def _rfft_factory_ ( klass , args , params ) :
    """Factory for `ROOT.RooFFTConvPdf` 
    - see `ROOT.RooFFTConvPdf`
    """
    pdf = klass ( *args )
    pdf.__args = args 
    bs , bf , s1 , s2 = params
    pdf.setBufferStrategy ( bs )
    pdf.setBufferFraction ( bf )
    pdf.setShift ( s1 , s2 ) 
    return pdf

# =============================================================================
## reduce RooFFTConvPdf
#  @see RooFFTConvPdf 
def _rfft_reduce_ ( pdf ) :

    s1 = ctypes.c_double ( 0 )
    s2 = ctypes.c_double ( 0 )

    pars   = Ostap.MoreRooFit.fft_pars ( pdf , s1 , s2 )
    args   = ( pdf.name , pdf.title )  + \
             tuple ( p for p in pars ) + \
             ( pdf.getInterpolationOrder() , ) 
    
    params =  pdf.bufferStrategy() , pdf.bufferFraction() , s1.value , s2.value

    return _rfft_factory_ , ( type ( pdf ) , args , params )
    
ROOT.RooFFTConvPdf.__reduce__ = _rfft_reduce_ 

# =============================================================================
## Factory for RooSimultaneous
#  @see RooSimultaneous 
def _rsim_factory_ ( klass , args , catlst ) :
    """Factory for `ROOT.RooSimultaneous` 
    - see `ROOT.Simultaneous`
    """
    pdf = klass ( *args )
    for l , p in catlst : pdf.addPdf ( p , l )
    pdf.__args   = args 
    pdf.__catlst = catlst 
    return pdf

# ================================================================================
## Reduce RooSimultaneous object
#  @see RooSimultaneous 
def _rsim_reduce_ ( pdf ) :
    """Reduce RooSimultaneous object
    -see RooSimultaneous 
    """
    cat    = pdf.indexCat()
    labels = cat.labels()
    catlst = tuple ( ( l , pdf.getPdf ( l ) ) for l in labels )
    args   = pdf.name , pdf.title , cat 
    return _rsim_factory_ , ( type ( pdf ) , args , catlst ) 

ROOT.RooSimultaneous.__reduce__  = _rsim_reduce_ 

# =============================================================================
## access to underlying efficiency function from RooEfficiency
#  @see RooEfficiency 
#  @see Ostap::MorERooFit::get_eff 
def _reff_efficiency_  ( pdf )  :
    """Access to underlying efficiency function from RooEfficiency
    - see `ROOT.RooEfficiency`
    - see `Ostap.MoreRooFit.get_eff`
    """
    return Ostap.MoreRooFit.get_eff ( pdf )

# =============================================================================
## access to underlying accept/reject category from RooEfficiency
#  @see RooEfficiency 
#  @see Ostap::MorERooFit::get_cat 
def _reff_category_  ( pdf )  :
    """Access to underlying accept/reject category from RooEfficiency
    - see `ROOT.RooEfficiency`
    - see `Ostap.MoreRooFit.get_cat`
    """
    return Ostap.MoreRooFit.get_cat ( pdf )

# =============================================================================
## access to accept category from RooEfficiency
#  @see RooEfficiency 
#  @see Ostap::MorERooFit::get_acc
def _reff_accept_  ( pdf )  :
    """Access to accept category from RooEfficiency
    - see `ROOT.RooEfficiency`
    - see `Ostap.MoreRooFit.get_acc`
    """
    return Ostap.MoreRooFit.get_acc ( pdf )

ROOT.RooEfficiency. efficiency = _reff_efficiency_
ROOT.RooEfficiency. category   = _reff_category_
ROOT.RooEfficiency. accept     = _reff_accept_

_new_methods_ += [
    ROOT.RooEfficiency. efficiency , 
    ROOT.RooEfficiency. category   , 
    ROOT.RooEfficiency. accept     , 
    ]

# =============================================================================
## reduce RooEfficiency
#  @see RooEfficiency
def _reff_reduce_ ( pdf ) :
    """Reduce `ROOT.RooEfficiency`
    - see `ROOT.RooEfficiency`
    """
    return root_store_factory , ( type ( pdf )     ,
                                  pdf.name         ,
                                  pdf.title        ,
                                  pdf.efficiency() , 
                                  pdf.category  () , 
                                  pdf.accept    () )

ROOT.RooEfficiency.__reduce__  = _reff_reduce_ 

# ================================================================================
## get the list of coefficients from <code>RooPolyVar</code> & <code>RooPolynomial</code> & 
#  @see RooPolyVar
#  @see RooPolynomial
def _rpv_coefficients_ ( var ) :
    """Get the list of coefficients from `ROOT.RooPolyVar` & `ROOT.RooPolynomial` &
    - see `ROOT.RooPolyVar`
    - see `ROOT.RooPolynomial`
    """
    return Ostap.MoreRooFit.coefficients ( var ) 
# ================================================================================
## get the variable from <code>RooPolyVar</code> & <code>RooPolynomial</code> & 
#  @see RooPolyVar
#  @see RooPolynomial
def _rpv_variable_ ( var ) :
    """Get the variable from `ROOT.RooPolyVar` & `ROOT.RooPolynomial` &
    - see `ROOT.RooPolyVar`
    - see `ROOT.RooPolynomial`
    """
    return Ostap.MoreRooFit.get_variable ( var ) 
# ================================================================================
## get the lowest order from <code>RooPolyVar</code> & <code>RooPolynomial</code> & 
#  @see RooPolyVar
#  @see RooPolynomial
def _rpv_lowest_order_ ( var ) :
    """Get the lowest order from `ROOT.RooPolyVar` & `ROOT.RooPolynomial` &
    - see `ROOT.RooPolyVar`
    - see `ROOT.RooPolynomial`
    """
    return Ostap.MoreRooFit.lowest_order( var ) 



ROOT.RooPolyVar   . coefficients = _rpv_coefficients_
ROOT.RooPolynomial. coefficients = _rpv_coefficients_
ROOT.RooPolyVar   . get_variable = _rpv_variable_
ROOT.RooPolynomial. get_variable = _rpv_variable_
ROOT.RooPolyVar   . lowest_order = _rpv_lowest_order_
ROOT.RooPolynomial. lowest_order = _rpv_lowest_order_


_new_methods_ += [
    ROOT.RooPolyVar   . coefficients ,
    ROOT.RooPolynomial. coefficients ,
    ROOT.RooPolyVar   . get_variable ,
    ROOT.RooPolynomial. get_variable ,
    ROOT.RooPolyVar   . lowest_order ,
    ROOT.RooPolynomial. lowest_order ,
    ]

# =============================================================================
## reduce RooPolyVar & RooPolynomial
#  @see RooPolyVar
#  @see RooPolynomial
def _rpv_reduce_ ( var ) :
    """Reduce `ROOT.RooPolyVar` & `ROOT.RooPolynomial`
    - see `ROOT.RooPolyVar`
    - see `ROOT.RooPolynomial`
    """
    return root_store_factory , ( type ( var )        ,
                                  var.name            ,
                                  var.title           ,
                                  var.get_variable () , 
                                  var.coefficients () ,
                                  var.lowest_order () )
                                  
ROOT.RooPolyVar   . __reduce__  = _rpv_reduce_
ROOT.RooPolynomial. __reduce__  = _rpv_reduce_


# ================================================================================
## deserialize RooFitResult
#  @see RooFitResult
#  @see Ostap::Utils::FitResults
def _rrfr_factory_ ( klass , *args ) :
    """Deserialize RooFitResult
    - see `ROOT.RooFitResult`
    - see `Ostap.Utils.FitResults`
    """
    ## create

    newargs = args [ : 10 ] 
    covargs = args [   10 ] ## covariance related arguments 
    history = args [   11 ] ## the last argument, history 

    import pickle

    if 3 == len ( covargs ) :
        
        gcc , corm , covm = covargs 
        D    = len ( gcc )
        
        if D * D != len ( corm ) : raise pickle.UnpicklingError("FitResult(1): cannot reconstruct cormatrix")
        if D * D != len ( covm ) : raise pickle.UnpicklingError("FitResult(1): cannot reconstruct covmatrix")

        gcc      = doubles ( gcc ) 
        corm     = ROOT.TMatrixDSym ( D , corm )
        covm     = ROOT.TMatrixDSym ( D , covm )
        
        newargs +=  gcc , corm , covm
        
    else :
        
        covm = covargs 
        
        l2 = len ( covm )
        D  = int ( math.sqrt ( l2 ) ) 
        if D * D != l2 :
            raise pickle.UnpicklingError("FitResult(2): cannot reconstruct covmatrix")
        covm = ROOT.TMatrixDSym ( D , covm )

        newargs += covm ,
            
    tmpres = klass ( *newargs )
    
    for label , code in history : tmpres.add_to_history ( label , code )
    
    result = ROOT.RooFitResult ( tmpres )
    result.__args = args
    
    return result 
   
# ================================================================================
## Reduce RooFitResult
#  @see RooFitResult
def _rrfr_reduce_ ( res ) :
    """ Reduce `RooFitResult`
    - see `ROOT.RooFitResult`
    """
    nr      = Ostap.Utils.FitResults ( res )
    content = type ( nr )      , res.name , res.title , \
              res.constPars () , res.floatParsInit()  , res.floatParsFinal() , \
              nr.status     () , nr.covQual()         , nr.minNll()          , \
               nr.edm       () , nr.numInvalidNLL ()

    gcc  = nr.global_cc()
    N    = len ( gcc )
    covm = res.covarianceMatrix  ()
    corm = res.correlationMatrix ()
    if gcc and 1 <= N and covm.IsValid() and corm.IsValid() and N == covm.GetNrows() and N == corm.GetNrows() : 
        covargs = ( array.array ( 'd' , gcc  ) ,
                    array.array ( 'd' , corm ) ,
                    array.array ( 'd' , covm ) )
    else :
        covargs = array.array ( 'd' , covm )  
        
    history = tuple ( ( nr.statusLabelHistory(i) ,nr.statusCodeHistory(i) ) \
                      for i in range ( nr.numStatusHistory () ) )
    
    content += covargs , history 
  
    return _rrfr_factory_ , content 

ROOT.RooFitResult.__reduce__  = _rrfr_reduce_ 

# =============================================================================
## reconstruct/deserialize <code>RooPlot</code> object
def _rplot_factory_ ( klass , xmin , xmax , ymin , ymax , items )  :
    """Reconstruct/deserialize `ROOT.RooPlot` object
    """
    plot = klass  ( xmin , xmax , ymin , ymax )
    ##
    for ( obj , options , invisible ) in items :
        if   isinstance ( obj  , ROOT.RooPlotable ) :            
            plot.addPlotable ( obj , options , invisible )
        elif isinstance ( obj  , ROOT.TH1 ) and 1 == obj.GetDimension() :
            plot.addTH1      ( obj , options , invisible )
        else :
            plot.addObject   ( obj , options , invisible )            
    ## 
    plot.__store = items  ## ATTENTION!! keep the items! 
    return plot 
                

# =============================================================================
## reduce RooPlot object
def _rplot_reduce_ ( plot ) :
    """Reduce `ROOT.RooPlot` object"""
    return _rplot_factory_ , ( type ( plot )   ,
                               plot.GetXaxis().GetXmin() ,
                               plot.GetXaxis().GetXmax() ,
                               plot.GetMinimum () ,
                               plot.GetMaximum () ,
                               tuple ( i for i in plot.items() ) )

ROOT.RooPlot.__reduce__  = _rplot_reduce_ 

# ===============================================================================
## reconstruct/deserialize/unpickle <code>TGraph</code> object
#  @see RooCurve
#  @see RooPlotable
#  @see TGraph
def _rcurv_factory_ ( klass , *args ) : 
    """Reconstruct/deserialize/unpickle `ROOT.RooCurve` object
    -see `ROOT.RooCurve`
    -see `ROOT.RooPlotable`
    -see `ROOT.TGraph`
    """
    graph     = GR.graph_factory ( klass , *args [:-1] )
    ## 
    rplotatts = args[-1]
    graph.setYAxisLabel  (  rplotatts [0 ] )
    graph.setYAxisLimits ( *rplotatts [1:] )
    ## 
    return graph 

# =============================================================================
## Reduce/serialize/pickle  simple <code>TGraph</code> object
#  @see ROOT.TGraph
def _rcurv_reduce_ ( graph ) :
    """Reduce/serialize/pickle simple `ROOT.RooCurve` object\
    - see `ROOT.RooCurve`
    - see `ROOT.RooEllipse`
    - see `ROOT.RooPlotable`
    - see `ROOT.TGraph`
    """
    _ , content = GR.graph_reduce ( graph )
    rplotatts   = graph.getYAxisLabel () , graph.getYAxisMin (), graph.getYAxisMax ()
    return _rcurv_factory_ , content + ( rplotatts , ) 

ROOT.RooCurve  .__reduce__ = _rcurv_reduce_
ROOT.RooEllipse.__reduce__ = _rcurv_reduce_

# ===============================================================================
## reconstruct/deserialize/unpickle <code>RooHist</code> object
#  @see RooHist
#  @see RooPlotable
#  @see TGraphAsymmErrors
def _rhist_factory_ ( klass , *args ) :
    """Reconstruct/deserialize/unpickle `ROOT.RooHist` object
    -see `ROOT.RooHist`
    -see `ROOT.RooPlotable`
    -see `ROOT.TGraphAsymmErrors`
    """
    graph     = GR.graph_asymerrors_factory ( klass , *args [:-1] )
    rplotatts = args[-1]
    graph.setYAxisLabel  (  rplotatts [0 ] )
    graph.setYAxisLimits ( *rplotatts [1:] )
    ## 
    return graph 

# =============================================================================
## Reduce/serialize/pickle  simple <code>RooHist</code> object
#  @see ROOT.RooHist
#  @see ROOT.TGraphAsymmErorrs
def _rhist_reduce_ ( graph ) :
    """Reduce/serialize/pickle simple `ROOT.RooHist` object
    - see `ROOT.RooHist`
    - see `ROOT.TGraphAsymmErrors`
    """
    _ , content = GR.graph_asymerrors_reduce ( graph )
    rplotatts   = graph.getYAxisLabel (), graph.getYAxisMin (), graph.getYAxisMax () 
    return _rhist_factory_ , content + ( rplotatts , ) 

ROOT.RooHist.__reduce__ = _rhist_reduce_







# =============================================================================
## Reduce <code>Ostap::MoreRooFit::TwoVars</code> objects
#  @see Ostap::MoreRooFit.TwoVars
def _r2v_reduce ( var ) :
    """Reduce `Ostap::MoreRooFit::TwoVars` objects
    - see Ostap.MoreRooFit.TwoVars
    """
    return root_store_factory , ( type ( var ) , var.name , var.title , var.x() , var.y() )

Ostap.MoreRooFit.TwoVars.__reduce__  = _r2v_reduce

# ===================================================================
## Reduce <code>Ostap::MoreRooFit::Addition</code> objects
#  @see Ostap::MoreRooFit.Addition
def _radd1_reduce ( var ) :
    """Reduce `Ostap.MoreRooFit.Addition` objects
    - see Ostap.MoreRooFit.Addition
    """
    content = type ( var ) , var.name , var.title , var.x () , var.y ()
    return root_store_factory , content 

# ===================================================================
## Reduce <code>Ostap::MoreRooFit::Addition3</code> objects
#  @see Ostap::MoreRooFit.Addition3
def _radd2_reduce ( var ) :
    """Reduce `Ostap.MoreRooFit.Addition` objects
    - see Ostap.MoreRooFit.Addition2
    """
    content = type ( var ) , var.name , var.title , var.x () , var.y () , var.c1 () , var.c2 () 
    return root_store_factory , content

# =============================================================================
## reduce <code>Ostap::MoreRooFit::Id<code> object
#  @see Ostap.MoreRooFit::Id
def _rid_reduce ( var ) :
    """Reduce `Ostap.MoreRooFit.Id` object
    - see Ostap.MoreRooFit.Id
    """
    return root_store_factory , ( type ( var ) ,
                                  var.name     ,
                                  var.title    ,
                                  var.x ()     )

# =============================================================================
## reduce <code>Ostap::MoreRooFit::AddDeps<code> object
#  @see Ostap.MoreRooFit::AddDeps
def _radddep_reduce ( var ) :
    """Reduce `Ostap.MoreRooFit.AddDeps` object
    - see Ostap.MoreRooFit.AddDeps
    """
    return root_store_factory , ( type ( var ) ,
                                  var.name     ,
                                  var.title    ,
                                  var.x ()     ,
                                  var.vlst()   )



Ostap.MoreRooFit.Addition    .__reduce__  = _radd1_reduce
Ostap.MoreRooFit.Addition2   .__reduce__  = _radd2_reduce

Ostap.MoreRooFit.Subtraction .__reduce__  = _r2v_reduce
Ostap.MoreRooFit.Product     .__reduce__  = _r2v_reduce
Ostap.MoreRooFit.ProductPdf  .__reduce__  = _r2v_reduce

Ostap.MoreRooFit.Id          .__reduce__  = _rid_reduce
Ostap.MoreRooFit.AddDeps     .__reduce__  = _radddep_reduce

# ===================================================================
## Reduce <code>Ostap::MoreRooFit::Combination</code> objects
#  @see Ostap::MoreRooFit.Combination
def _rcomb_reduce ( var ) :
    """Reduce `Ostap.MoreRooFit.Combination` objects
    - see Ostap.MoreRooFit.Combination
    """
    return root_store_factory , ( type ( var ) ,
                                  var.name     ,
                                  var.title    ,
                                  var.x()      ,
                                  var.y()      ,
                                  var.alpha () ,
                                  var.beta  () ,
                                  var.gamma () )

Ostap.MoreRooFit.Combination.__reduce__  = _rcomb_reduce

# ===================================================================
## Reduce <code>Ostap::MoreRooFit::Asymmetry</code> objects
#  @see Ostap::MoreRooFit.Asymmetry
def _rasym_reduce ( var ) :
    """Reduce  `Ostap::MoreRooFit::Asymmetry` objects
    - see Ostap.MoreRooFit.Asymmetry
    """
    return root_store_factory , ( type ( var ) ,
                                  var.name , var.title ,
                                  var.x()  , var.y()   ,
                                  var.scale () )

Ostap.MoreRooFit.Asymmetry.  __reduce__  = _rasym_reduce

# =============================================================================
## Reduce <code>Ostap::MoreRooFit::Constant</code>
def _rconst2_reduce ( var ) :
    """ Reduce `ROOT.RooConstVar`
    """
    return root_store_factory , ( type ( var )  ,
                                  var.name      ,
                                  var.title     ,
                                  float ( var ) ,
                                  var.vlst()    )

Ostap.MoreRooFit.Constant.  __reduce__  = _rconst2_reduce

# =============================================================================
## Reduce <code>Ostap::MoreRooFit::Bernstein</code> object 
#  @see Ostap::MoreRooDit::Bernstein
def _rbern_reduce ( var ) :
    """Reduce `Ostap.MoreRooFit.Bernstein` object 
    - see Ostap.MoreRooDit.Bernstein
    """
    return root_store_factory , ( type ( var )  ,
                                  var.name      ,
                                  var.title     ,
                                  var.xvar   () ,
                                  var.pars   () ,                              
                                  var.xmin   () ,
                                  var.xmax   () )

Ostap.MoreRooFit.Bernstein.  __reduce__  = _rbern_reduce

# =============================================================================
## Reduce <code>Ostap::MoreRooFit::Monotonic</code> object 
#  @see Ostap::MoreRooDit::Monotonic
def _rmono_reduce ( var ) :
    """Reduce `Ostap.MoreRooFit.Monotonic` object 
    - see Ostap.MoreRooDit.Monotonic
    """
    return root_store_factory , ( type ( var )  ,
                                  var.name      ,
                                  var.title     ,
                                  var.xvar ()   ,
                                  var.pars ()   ,                               
                                  var.monotonic().increasing() , 
                                  var.xmin ()   ,
                                  var.xmax ()   ,
                                  var.a    ()   ,
                                  var.b    ()   )

Ostap.MoreRooFit.Monotonic.  __reduce__  = _rmono_reduce

# =============================================================================
## Reduce <code>Ostap::MoreRooFit::Convex</code> object 
#  @see Ostap::MoreRooDit::Monotonic
def _rconv1_reduce ( var ) :
    """Reduce `Ostap.MoreRooFit.Monotonic` object 
    - see Ostap.MoreRooDit.Monotonic
    """
    return root_store_factory , ( type ( var )  ,
                                  var.name      ,
                                  var.title     ,
                                  var.xvar ()   ,
                                  var.pars ()   ,                               
                                  var.convex ().increasing () , 
                                  var.convex ().convex     () , 
                                  var.xmin ()   ,
                                  var.xmax ()   ,
                                  var.a    ()   ,
                                  var.b    ()   )

Ostap.MoreRooFit.Convex.  __reduce__  = _rconv1_reduce

# =============================================================================
## Reduce <code>Ostap::MoreRooFit::ConvexOnly</code> object 
#  @see Ostap::MoreRooDit::ConvexOnly
def _rconv2_reduce ( var ) :
    """Reduce `Ostap.MoreRooFit.ConvexOnly` object 
    - see Ostap.MoreRooDit.ConvexOnly
    """
    return root_store_factory , ( type ( var )  ,
                                  var.name      ,
                                  var.title     ,
                                  var.xvar ()   ,
                                  var.pars ()   ,                               
                                  var.convex ().convex () , 
                                  var.xmin ()   ,
                                  var.xmax ()   ,
                                  var.a    ()   ,
                                  var.b    ()   )

Ostap.MoreRooFit.ConvexOnly.  __reduce__  = _rconv2_reduce

# =============================================================================
## unpickle <code>Ostap::MoreRooFit::BSpline</code> objects
#  @see Ostap::MoreRooFit::BSpline
def _rbspl_factory ( klass , name , title , xvar , knots , pars ) :
    """unpickle `Ostap.MoreRooFit.BSpline` objects
    - see Ostap.MoreRooFit.BSpline
    """
    plst  = ROOT.RooArgList()
    for v in pars : plst.add ( v )
    knots = doubles ( knots )
    #
    obj        = klass ( name , title , xvar , knots , pars )
    obj.__args = pars ## ATTENTION HERE! keep the args! 
    return obj 

# =============================================================================
## Reduce <code>Ostap::MoreRooFit::BSpline</code> object 
#  @see Ostap::MoreRooDit::BSpline
def _rbspl_reduce ( spl ) :
    """Reduce `Ostap.MoreRooFit.BSpline` object 
    - see Ostap.MoreRooDit.BSpline
    """
    return _rbspl_factory, ( type ( spl ) ,
                             spl.name     ,
                             spl.title    ,
                             spl.xvar  () ,
                             array.array( 'd' , spl.knots() ) ,
                             spl.pars()   )

Ostap.MoreRooFit.BSpline.  __reduce__  = _rbspl_reduce

# ============================
## reduce <code>Ostap::Models::Uniform<code> object
#  @see Ostap::Models.Uniform 
def _runi_reduce_ ( uni ) :
    """reduce <code>Ostap::Models::Uniform<code> object
    - see Ostap::Models.Uniform 
    """
    tail = ()
    ##
    d = uni.dim() 
    if   1 == d : tail = uni.x() ,
    elif 2 == d : tail = uni.x() , uni.y()
    elif 3 == d : tail = uni.x() , uni.y() , uni.z()
    ##
    return root_store_factory , ( type ( uni ) ,
                                  uni.name     ,
                                  uni.title    ) + tail

Ostap.Models.Uniform.__reduce__ = _runi_reduce_



# =============================================================================
## reduce BreitWigner
def _rbw_reduce_ ( pdf ):
    """Reduce BreitWigner"""
    return root_store_factory , ( type ( pdf )       ,
                                  pdf.name           ,
                                  pdf.title          ,
                                  pdf.x      ()      , 
                                  pdf.mass   ()      ,
                                  pdf.widths ()      ,
                                  pdf.breit_wigner() ) 

Ostap.Models.BreitWigner.__reduce__ = _rbw_reduce_ 

# =============================================================================
## reduce BreitWignerMC
def _rbwmc_reduce_ ( pdf ):
    """Reduce BreitWignerMC"""
    return root_store_factory , ( type ( pdf )          ,
                                  pdf.name              ,
                                  pdf.title             ,
                                  pdf.x      ()         , 
                                  pdf.mass   ()         ,
                                  pdf.widths ()         ,
                                  pdf.breit_wigner_MC() ) 

Ostap.Models.BreitWignerMC.__reduce__ = _rbwmc_reduce_ 

# =============================================================================
## reduce BWI
def _rbwi_reduce_ ( pdf ):
    """Reduce BWI"""
    return root_store_factory , ( type ( pdf )     ,
                                  pdf.name         ,
                                  pdf.title        ,
                                  Ostap.Models.BreitWigner ( pdf , 'pkl' + pdf.name ) ,
                                  pdf.magnitude () ,
                                  pdf.phase     () ,
                                  pdf.scale1    () ,
                                  pdf.scale2    () )

Ostap.Models.BWI.__reduce__ = _rbwi_reduce_ 

# =============================================================================
## reduce Flatte
def _rflatte_reduce_ ( pdf ):
    """Reduce Flatte"""
    return root_store_factory , ( type ( pdf )     ,
                                  pdf.name         ,
                                  pdf.title        ,
                                  pdf.x      ()    , 
                                  pdf.mass   ()    ,
                                  pdf.widths ()[0] ,
                                  pdf.widths ()[1] ,
                                  pdf.widths ()[2] ,
                                  pdf.flatte  ()   ) 

Ostap.Models.Flatte.__reduce__ = _rflatte_reduce_ 

# =============================================================================
## reduce FlatteBugg
def _rflattebugg_reduce_ ( pdf ):
    """Reduce FlatteBugg"""
    return root_store_factory , ( type ( pdf )     ,
                                  pdf.name         ,
                                  pdf.title        ,
                                  pdf.x      ()    , 
                                  pdf.mass   ()    ,
                                  pdf.widths ()[0] ,
                                  pdf.widths ()[1] ,
                                  pdf.widths ()[2] ,
                                  pdf.flatte_bugg () ) 

Ostap.Models.FlatteBugg.__reduce__ = _rflattebugg_reduce_ 

# =============================================================================
## reduce LASS
def _rlass_reduce_ ( pdf ):
    """Reduce LASS"""
    return root_store_factory , ( type ( pdf )     ,
                                  pdf.name         ,
                                  pdf.title        ,
                                  pdf.x      ()    , 
                                  pdf.mass   ()    ,
                                  pdf.widths ()[0] ,
                                  pdf.a      ()    , 
                                  pdf.b      ()    , 
                                  pdf.e      ()    , 
                                  pdf.lass   ()    )

Ostap.Models.LASS.__reduce__ = _rlass_reduce_ 

# =============================================================================
## reduce BWPS
def _rbwps_reduce_ ( pdf ):
    """Reduce BWPS"""
    return root_store_factory , ( type ( pdf )     ,
                                  pdf.name         ,
                                  pdf.title        ,
                                  pdf.x      ()    , 
                                  pdf.mass   ()    ,
                                  pdf.widths ()    ,
                                  pdf.phis   ()    ,
                                  pdf.bwps   ()    )

Ostap.Models.BWPS.__reduce__ = _rbwps_reduce_ 

# =============================================================================
## reduce BW3L
def _rbw3l_reduce_ ( pdf ):
    """Reduce BW3L"""
    return root_store_factory , ( type ( pdf )     ,
                                  pdf.name         ,
                                  pdf.title        ,
                                  pdf.x      ()    , 
                                  pdf.mass   ()    ,
                                  pdf.widths ()    ,
                                  pdf.bw3l   ()    )

Ostap.Models.BW3L.__reduce__ = _rbw3l_reduce_ 


# =============================================================================
## reduce Voigt
def _rvoigt_reduce_ ( pdf ):
    """Reduce Voigt"""
    return root_store_factory , ( type ( pdf )  ,
                                  pdf.name      ,
                                  pdf.title     ,
                                  pdf.x      () , 
                                  pdf.m0     () ,
                                  pdf.gamma  () ,
                                  pdf.sigma  () )

Ostap.Models.Voigt      .__reduce__ = _rvoigt_reduce_ 
Ostap.Models.PseudoVoigt.__reduce__ = _rvoigt_reduce_ 

# =============================================================================
## reduce CrystalBall
def _rcb_reduce_ ( pdf ):
    """Reduce CristalBall"""
    return root_store_factory , ( type ( pdf )  ,
                                  pdf.name      ,
                                  pdf.title     ,
                                  pdf.x      () , 
                                  pdf.m0     () ,
                                  pdf.sigma  () ,
                                  pdf.alpha  () ,
                                  pdf.n      () )

Ostap.Models.CrystalBall   .__reduce__ = _rcb_reduce_ 
Ostap.Models.CrystalBallRS .__reduce__ = _rcb_reduce_ 

# =============================================================================
## reduce CrystalBallDS
def _rcb2_reduce_ ( pdf ):
    """Reduce CristalBallDS"""
    return root_store_factory , ( type ( pdf )  ,
                                  pdf.name      ,
                                  pdf.title     ,
                                  pdf.x      () , 
                                  pdf.m0     () ,
                                  pdf.sigma  () ,
                                  pdf.alphaL () ,
                                  pdf.nL     () ,
                                  pdf.alphaR () ,
                                  pdf.nR     () )

Ostap.Models.CrystalBallDS .__reduce__ = _rcb2_reduce_ 

# =============================================================================
## reduce Needham
def _rneedham_reduce_ ( pdf ):
    """Reduce Needham"""
    return root_store_factory , ( type ( pdf )  ,
                                  pdf.name      ,
                                  pdf.title     ,
                                  pdf.x      () , 
                                  pdf.m0     () ,
                                  pdf.sigma  () ,                            
                                  pdf.a0     () ,
                                  pdf.a1     () ,
                                  pdf.a2     () )

Ostap.Models.Needham .__reduce__ = _rneedham_reduce_ 

# =============================================================================
## reduce Apollonious
def _rapo_reduce_ ( pdf ):
    """Reduce Apollonios"""
    return root_store_factory , ( type ( pdf )  ,
                                  pdf.name      ,
                                  pdf.title     ,
                                  pdf.x      () , 
                                  pdf.m0     () ,
                                  pdf.sigma  () ,                            
                                  pdf.alpha  () ,                            
                                  pdf.n      () ,
                                  pdf.b      () )

Ostap.Models.Apollonios.__reduce__ = _rapo_reduce_ 

# =============================================================================
## reduce Apollonious2
def _rapo2_reduce_ ( pdf ):
    """Reduce Apollonios2"""
    return root_store_factory , ( type ( pdf )  ,
                                  pdf.name      ,
                                  pdf.title     ,
                                  pdf.x      () , 
                                  pdf.m0     () ,
                                  pdf.sigmaL () ,                            
                                  pdf.sigmaR () ,                            
                                  pdf.beta   () )

Ostap.Models.Apollonios2.__reduce__ = _rapo2_reduce_ 

# =============================================================================
## reduce BifurcatedGauss
def _rgbf_reduce_ ( pdf ):
    """Reduce BifurcatedGauss"""
    return root_store_factory , ( type ( pdf )  ,
                                  pdf.name      ,
                                  pdf.title     ,
                                  pdf.x      () , 
                                  pdf.peak   () ,
                                  pdf.sigmaL () ,                            
                                  pdf.sigmaR () )

Ostap.Models.BifurcatedGauss.__reduce__ = _rgbf_reduce_ 


# =============================================================================
## reduce GenGaussV1
def _rggv1_reduce_ ( pdf ):
    """Reduce GenGaussV1"""
    return root_store_factory , ( type ( pdf )  ,
                                  pdf.name      ,
                                  pdf.title     ,
                                  pdf.x      () , 
                                  pdf.mu     () ,
                                  pdf.alpha  () ,                            
                                  pdf.beta   () )

Ostap.Models.GenGaussV1.__reduce__ = _rggv1_reduce_ 


# =============================================================================
## reduce GenGaussV2
def _rggv2_reduce_ ( pdf ):
    """Reduce GenGaussV2"""
    return root_store_factory , ( type ( pdf )  ,
                                  pdf.name      ,
                                  pdf.title     ,
                                  pdf.x      () , 
                                  pdf.xi     () ,
                                  pdf.alpha  () ,                            
                                  pdf.kappa  () )

Ostap.Models.GenGaussV2.__reduce__ = _rggv2_reduce_ 

# =============================================================================
## reduce SkewGauss
def _rskg_reduce_ ( pdf ):
    """Reduce SkewGauss"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () , 
                                  pdf.xi       () ,
                                  pdf.omega    () ,                            
                                  pdf.alpha    () )

Ostap.Models.SkewGauss.__reduce__ = _rskg_reduce_ 

# =============================================================================
## reduce ExGauss
def _rexg_reduce_ ( pdf ):
    """Reduce ExGauss"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () , 
                                  pdf.mu       () ,
                                  pdf.varsigma () ,                            
                                  pdf.k        () )

Ostap.Models.ExGauss.__reduce__ = _rexg_reduce_ 

# =============================================================================
## reduce NormalLaplace
def _rnl_reduce_ ( pdf ):
    """Reduce NormalLaplace"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () , 
                                  pdf.mu       () ,
                                  pdf.varsigma () ,                            
                                  pdf.kL       () ,
                                  pdf.kL       () )

Ostap.Models.NormalLaplace.__reduce__ = _rnl_reduce_ 

# =============================================================================
## reduce Novosibirsk
def _rnovo_reduce_ ( pdf ):
    """Reduce Novisibirsk"""
    return root_store_factory , ( type ( pdf )  ,
                                  pdf.name      ,
                                  pdf.title     ,
                                  pdf.x      () , 
                                  pdf.peak   () ,
                                  pdf.sigma  () ,                            
                                  pdf.tau    () )

Ostap.Models.Novosibirsk.__reduce__ = _rnovo_reduce_ 

# =============================================================================
## reduce Bukin
def _rbukin_reduce_ ( pdf ):
    """Reduce Bukin"""
    return root_store_factory , ( type ( pdf )  ,
                                  pdf.name      ,
                                  pdf.title     ,
                                  pdf.x      () , 
                                  pdf.peak   () ,
                                  pdf.sigma  () ,                            
                                  pdf.xi     () ,
                                  pdf.rhoL   () ,
                                  pdf.rhoR   () )

Ostap.Models.Bukin.__reduce__ = _rbukin_reduce_ 

# =============================================================================
## reduce StudentT
def _rstt_reduce_ ( pdf ):
    """Reduce StudentT"""
    return root_store_factory , ( type ( pdf )  ,
                                  pdf.name      ,
                                  pdf.title     ,
                                  pdf.x      () , 
                                  pdf.mu     () ,
                                  pdf.sigma  () ,                            
                                  pdf.n      () )

Ostap.Models.StudentT.__reduce__ = _rstt_reduce_ 


# =============================================================================
## reduce BifurcatedStudentT
def _rbstt_reduce_ ( pdf ):
    """Reduce BifurcatedStudentT"""
    return root_store_factory , ( type ( pdf )  ,
                                  pdf.name      ,
                                  pdf.title     ,
                                  pdf.x      () , 
                                  pdf.mu     () ,
                                  pdf.sigmaL () ,                            
                                  pdf.sigmaR () ,                            
                                  pdf.nL     () ,
                                  pdf.nR     () )

Ostap.Models.BifurcatedStudentT.__reduce__ = _rbstt_reduce_ 

# =============================================================================
## reduce PearsonIV
def _rp4_reduce_ ( pdf ):
    """Reduce PearsonIV"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () , 
                                  pdf.mu       () ,
                                  pdf.varsigma () ,                            
                                  pdf.n        () ,                            
                                  pdf.kappa    () )

Ostap.Models.PearsonIV.__reduce__ = _rp4_reduce_ 

# =============================================================================
## reduce SkewGenT
def _rsgt_reduce_ ( pdf ):
    """Reduce SkewGenT"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () , 
                                  pdf.mu       () ,
                                  pdf.sigma    () ,                            
                                  pdf.xi       () ,                            
                                  pdf.r        () ,                            
                                  pdf.zeta     () )

Ostap.Models.SkewGenT.__reduce__ = _rsgt_reduce_ 

# =============================================================================
## reduce SkewGenError
def _rsge_reduce_ ( pdf ):
    """Reduce SkewGenError"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () , 
                                  pdf.mu       () ,
                                  pdf.sigma    () ,                            
                                  pdf.xi       () ,                            
                                  pdf.p        () )

Ostap.Models.SkewGenError.__reduce__ = _rsge_reduce_ 

# =============================================================================
## reduce GramCharlierA
def _rgca_reduce_ ( pdf ):
    """Reduce GramCharlierA"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () , 
                                  pdf.mu       () ,
                                  pdf.sigma    () ,                            
                                  pdf.kappa3   () ,
                                  pdf.kappa4   () )

Ostap.Models.GramCharlierA.__reduce__ = _rgca_reduce_ 

# =============================================================================
## reduce PhaseSpace2
def _rps2_reduce_ ( pdf ):
    """Reduce PhaseSpace2"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () , 
                                  pdf.m1       () ,
                                  pdf.m2       () )

Ostap.Models.PhaseSpace2.__reduce__ = _rps2_reduce_ 

# =============================================================================
## reduce PhaseSpaceLeft
def _rpsl_reduce_ ( pdf ):
    """Reduce PhaseSpaceLeft"""
    return root_store_factory , ( type ( pdf )     ,
                                  pdf.name         ,
                                  pdf.title        ,
                                  pdf.x         () , 
                                  pdf.threshold () ,
                                  pdf.scale     () ,
                                  pdf.left      () )


Ostap.Models.PhaseSpaceLeft.__reduce__ = _rpsl_reduce_ 

# =============================================================================
## reduce PhaseSpaceRight
def _rpsr_reduce_ ( pdf ):
    """Reduce PhaseSpaceRight"""
    return root_store_factory , ( type ( pdf )     ,
                                  pdf.name         ,
                                  pdf.title        ,
                                  pdf.x         () , 
                                  pdf.threshold () ,
                                  pdf.L         () ,
                                  pdf.N         () )


Ostap.Models.PhaseSpaceRight.__reduce__ = _rpsr_reduce_ 

# =============================================================================
## reduce PhaseSpaceNL
def _rpsnl_reduce_ ( pdf ):
    """Reduce PhaseSpaceNL"""
    return root_store_factory , ( type ( pdf )     ,
                                  pdf.name         ,
                                  pdf.title        ,
                                  pdf.x         () , 
                                  pdf.low       () ,
                                  pdf.high      () ,
                                  pdf.N         () ,
                                  pdf.L         () )

Ostap.Models.PhaseSpaceNL.__reduce__ = _rpsnl_reduce_ 

# =============================================================================
## reduce PhaseSpacePol
def _rpspol_reduce_ ( pdf ):
    """Reduce PhaseSpacePol"""
    return root_store_factory , ( type ( pdf )     ,
                                  pdf.name         ,
                                  pdf.title        ,
                                  pdf.x         () , 
                                  pdf.psNL      () ,
                                  pdf.phis      () )

Ostap.Models.PhaseSpacePol.__reduce__ = _rpspol_reduce_ 


# =============================================================================
## reduce PhaseSpaceLeftExpoPol
def _rpslepol_reduce_ ( pdf ):
    """Reduce PhaseSpaceLeftExpoPol"""
    return root_store_factory , ( type ( pdf )     ,
                                  pdf.name         ,
                                  pdf.title        ,
                                  pdf.x         () , 
                                  pdf.psleft    () ,
                                  pdf.tau       () ,
                                  pdf.scale     () ,
                                  pdf.phis      () )

Ostap.Models.PhaseSpaceLeftExpoPol.__reduce__ = _rpslepol_reduce_ 


## PhaseSpace23L - not done .... why?


# =============================================================================
## reduce PolyPositive 
def _rpolpos_reduce_ ( pdf ):
    """Reduce PolyPositive"""
    return root_store_factory , ( type ( pdf )     ,
                                  pdf.name         ,
                                  pdf.title        ,
                                  pdf.x         () , 
                                  pdf.phis      () ,
                                  pdf.xmin      () ,
                                  pdf.xmax      () )

Ostap.Models.PolyPositive.__reduce__ = _rpolpos_reduce_

# =============================================================================
## reduce PolyPositiveEven 
def _rpolpose_reduce_ ( pdf ):
    """Reduce PolyPositiveEven"""
    return root_store_factory , ( type ( pdf )     ,
                                  pdf.name         ,
                                  pdf.title        ,
                                  pdf.x         () , 
                                  pdf.phis      () ,
                                  pdf.xmin      () ,
                                  pdf.xmax      () )

Ostap.Models.PolyPositiveEven.__reduce__ = _rpolpose_reduce_ 

# =============================================================================
## reduce PolyMonotonic
def _rpolmon_reduce_ ( pdf ):
    """Reduce PolyMonotonic"""
    return root_store_factory , ( type ( pdf )      ,
                                  pdf.name          ,
                                  pdf.title         ,
                                  pdf.x          () , 
                                  pdf.phis       () ,
                                  pdf.xmin       () ,
                                  pdf.xmax       () ,
                                  pdf.increasing () )

Ostap.Models.PolyMonotonic.__reduce__ = _rpolmon_reduce_ 

# =============================================================================
## reduce PolyConvex
def _rpolcon_reduce_ ( pdf ):
    """Reduce PolyConvex"""
    return root_store_factory , ( type ( pdf )      ,
                                  pdf.name          ,
                                  pdf.title         ,
                                  pdf.x          () , 
                                  pdf.phis       () ,
                                  pdf.xmin       () ,
                                  pdf.xmax       () ,
                                  pdf.increasing () ,
                                  pdf.convex     () )

Ostap.Models.PolyConvex.__reduce__ = _rpolcon_reduce_ 

# =============================================================================
## reduce PolyConvexOnly
def _rpolcono_reduce_ ( pdf ):
    """Reduce PolyConvexOnly"""
    return root_store_factory , ( type ( pdf )      ,
                                  pdf.name          ,
                                  pdf.title         ,
                                  pdf.x          () , 
                                  pdf.phis       () ,
                                  pdf.xmin       () ,
                                  pdf.xmax       () ,
                                  pdf.convex     () )

Ostap.Models.PolyConvexOnly.__reduce__ = _rpolcono_reduce_ 

# =============================================================================
## reduce ExpoPositive
def _rexppos_reduce_ ( pdf ):
    """Reduce ExpoPoisitive"""
    return root_store_factory , ( type ( pdf )      ,
                                  pdf.name          ,
                                  pdf.title         ,
                                  pdf.x          () , 
                                  pdf.tau        () ,
                                  pdf.phis       () ,
                                  pdf.xmin       () ,
                                  pdf.xmax       () )

Ostap.Models.ExpoPositive.__reduce__ = _rexppos_reduce_ 

# =============================================================================
## reduce PolySigmoid
def _rpolsigm_reduce_ ( pdf ):
    """Reduce PoLySigmoid"""
    return root_store_factory , ( type ( pdf )      ,
                                  pdf.name          ,
                                  pdf.title         ,
                                  pdf.x          () , 
                                  pdf.phis       () ,
                                  pdf.xmin       () ,
                                  pdf.xmax       () ,
                                  pdf.alpha      () ,
                                  pdf.x0         () )

Ostap.Models.PolySigmoid.__reduce__ = _rpolsigm_reduce_ 

# =============================================================================
## reduce TwoExpoPositive
def _r2exppos_reduce_ ( pdf ):
    """Reduce TwoExpoPositive"""
    return root_store_factory , ( type ( pdf )      ,
                                  pdf.name          ,
                                  pdf.title         ,
                                  pdf.x          () ,                            
                                  pdf.alpha      () ,
                                  pdf.delta      () ,
                                  pdf.x0         () ,
                                  pdf.phis       () ,
                                  pdf.xmin       () ,
                                  pdf.xmax       () )

Ostap.Models.TwoExpoPositive.__reduce__ = _r2exppos_reduce_ 

# =============================================================================
## reduce GammaDist
def _rgamdist_reduce_ ( pdf ):
    """Reduce GammaDist"""
    return root_store_factory , ( type ( pdf )      ,
                                  pdf.name          ,
                                  pdf.title         ,
                                  pdf.x          () ,                            
                                  pdf.k          () ,
                                  pdf.theta      () )

Ostap.Models.GammaDist     .__reduce__ = _rgamdist_reduce_ 
Ostap.Models.LogGammaDist  .__reduce__ = _rgamdist_reduce_ 
Ostap.Models.Log10GammaDist.__reduce__ = _rgamdist_reduce_ 

# =============================================================================
## reduce GenGammaDist
def _rggamdist_reduce_ ( pdf ):
    """Reduce GenGammaDist"""
    return root_store_factory , ( type ( pdf )      ,
                                  pdf.name          ,
                                  pdf.title         ,
                                  pdf.x          () ,                            
                                  pdf.k          () ,
                                  pdf.theta      () ,
                                  pdf.p          () ,
                                  pdf.low        () )

Ostap.Models.GenGammaDist.__reduce__ = _rggamdist_reduce_ 

# =============================================================================
## reduce Amoroso
def _ramoroso_reduce_ ( pdf ):
    """Reduce Amoroso"""
    return root_store_factory , ( type ( pdf )      ,
                                  pdf.name          ,
                                  pdf.title         ,
                                  pdf.x          () ,                            
                                  pdf.theta      () ,
                                  pdf.alpha      () ,
                                  pdf.beta       () ,
                                  pdf.a          () )

Ostap.Models.Amoroso.__reduce__ = _ramoroso_reduce_ 

# =============================================================================
## reduce LogGamma
def _rloggam_reduce_ ( pdf ):
    """Reduce LogGamma"""
    return root_store_factory , ( type ( pdf )      ,
                                  pdf.name          ,
                                  pdf.title         ,
                                  pdf.x          () ,                            
                                  pdf.nu         () ,
                                  pdf.lambd      () ,
                                  pdf.alpha      () )

Ostap.Models.LogGamma   .__reduce__ = _rloggam_reduce_ 

# =============================================================================
## reduce BetaPrime
def _rbetap_reduce_ ( pdf ):
    """Reduce BetaPrime"""
    return root_store_factory , ( type ( pdf )      ,
                                  pdf.name          ,
                                  pdf.title         ,
                                  pdf.x          () ,                            
                                  pdf.alpha      () ,
                                  pdf.beta       () ,
                                  pdf.scale      () ,
                                  pdf.shift      () )

Ostap.Models.BetaPrime.__reduce__ = _rbetap_reduce_ 

# =============================================================================
## reduce Landau
def _rlandau_reduce_ ( pdf ):
    """Reduce Landau"""
    return root_store_factory , ( type ( pdf )      ,
                                  pdf.name          ,
                                  pdf.title         ,
                                  pdf.x          () ,                            
                                  pdf.scale      () ,
                                  pdf.shift      () )

Ostap.Models.Landau.__reduce__ = _rlandau_reduce_ 

# =============================================================================
## reduce SinhAsinh
def _rshash_reduce_ ( pdf ):
    """Reduce SinhAsinh"""
    return root_store_factory , ( type ( pdf )   ,
                                  pdf.name       ,
                                  pdf.title      ,
                                  pdf.x       () ,                            
                                  pdf.mu      () ,
                                  pdf.sigma   () ,
                                  pdf.epsilon () ,
                                  pdf.delta   () )

Ostap.Models.SinhAsinh.__reduce__ = _rshash_reduce_ 

# =============================================================================
## reduce JohnsonSU
def _rjsu_reduce_ ( pdf ):
    """Reduce JohnsonSU"""
    return root_store_factory , ( type ( pdf )   ,
                                  pdf.name       ,
                                  pdf.title      ,
                                  pdf.x       () ,                            
                                  pdf.xi      () ,
                                  pdf.lambd   () ,
                                  pdf.delta   () ,
                                  pdf.gamma   () )

Ostap.Models.JohnsonSU.__reduce__ = _rjsu_reduce_ 

# =============================================================================
## reduce ATLAS
def _ratlas_reduce_ ( pdf ):
    """Reduce Atlas"""
    return root_store_factory , ( type ( pdf )   ,
                                  pdf.name       ,
                                  pdf.title      ,
                                  pdf.x       () ,                            
                                  pdf.mu      () ,
                                  pdf.sigma   () )

Ostap.Models.Atlas    .__reduce__ = _ratlas_reduce_ 
Ostap.Models.Sech     .__reduce__ = _ratlas_reduce_ 
Ostap.Models.Logistic .__reduce__ = _ratlas_reduce_ 
Ostap.Models.Hat      .__reduce__ = _ratlas_reduce_ 
Ostap.Models.Up       .__reduce__ = _ratlas_reduce_ 

# =============================================================================
## reduce BatesShae  
def _rbats_reduce_ ( pdf ):
    """Reduce BatesShape"""
    return root_store_factory , ( type ( pdf )   ,
                                  pdf.name       ,
                                  pdf.title      ,
                                  pdf.x       () ,                            
                                  pdf.mu      () ,
                                  pdf.sigma   () ,
                                  pdf.n       () )

Ostap.Models.BatesShape .__reduce__ = _rbats_reduce_ 

# =============================================================================
## reduce Slash
def _rslash_reduce_ ( pdf ):
    """Reduce Slasj"""
    return root_store_factory , ( type ( pdf )   ,
                                  pdf.name       ,
                                  pdf.title      ,
                                  pdf.x       () ,                            
                                  pdf.mu      () ,
                                  pdf.scale   () )

Ostap.Models.Slash    .__reduce__ = _rslash_reduce_ 

# =============================================================================
## reduce ARGUS
def _rargus_reduce_ ( pdf ):
    """Reduce ARGUS"""
    return root_store_factory , ( type ( pdf )   ,
                                  pdf.name       ,
                                  pdf.title      ,
                                  pdf.x       () ,                            
                                  pdf.mu      () ,
                                  pdf.c       () ,
                                  pdf.chi     () )

Ostap.Models.Argus  .__reduce__ = _rargus_reduce_ 

# =============================================================================
## reduce GenArgus
def _rgargus_reduce_ ( pdf ):
    """Reduce GenArgus"""
    return root_store_factory , ( type ( pdf )   ,
                                  pdf.name       ,
                                  pdf.title      ,
                                  pdf.x       () ,                            
                                  pdf.mu      () ,
                                  pdf.c       () ,
                                  pdf.chi     () ,
                                  pdf.dp     () )

Ostap.Models.GenArgus .__reduce__ = _rgargus_reduce_ 

# =============================================================================
## reduce Losev
def _rlosev_reduce_ ( pdf ):
    """Reduce Losev"""
    return root_store_factory , ( type ( pdf )   ,
                                  pdf.name       ,
                                  pdf.title      ,
                                  pdf.x       () ,                            
                                  pdf.mu      () ,
                                  pdf.alpha   () ,
                                  pdf.beta    () )

Ostap.Models.Losev.__reduce__ = _rlosev_reduce_ 

# =============================================================================
## reduce AsymmetricLaplace
def _ralap_reduce_ ( pdf ):
    """Reduce AsymmetricLaplace"""
    return root_store_factory , ( type ( pdf )   ,
                                  pdf.name       ,
                                  pdf.title      ,
                                  pdf.x       () ,                            
                                  pdf.mu      () ,
                                  pdf.lambdaL () ,
                                  pdf.lambdaR () )

Ostap.Models.AsymmetricLaplace.__reduce__ = _ralap_reduce_ 

# =============================================================================
## reduce FupN
def _rfupn_reduce_ ( pdf ):
    """Reduce FupN"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.N        () ,                            
                                  pdf.mu       () ,
                                  pdf.varsigma () )

Ostap.Models.FupN.__reduce__ = _rfupn_reduce_ 

# =============================================================================
## reduce Tsallis
def _rtsal_reduce_ ( pdf ):
    """Reduce Tsallis"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.pt       () ,                            
                                  pdf.n        () ,                            
                                  pdf.T        () ,
                                  pdf.mass     () )

Ostap.Models.Tsallis.__reduce__ = _rtsal_reduce_ 

# =============================================================================
## reduce QGSM
def _rqgsm_reduce_ ( pdf ):
    """Reduce QGSM"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.b        () ,                            
                                  pdf.mass     () )

Ostap.Models.QGSM.__reduce__ = _rqgsm_reduce_ 

# =============================================================================
## reduce Hagedorn
def _rhage_reduce_ ( pdf ):
    """Reduce Hagedorn"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.beta     () ,                            
                                  pdf.mass     () )

Ostap.Models.Hagedorn.__reduce__ = _rhage_reduce_ 

# =============================================================================
## reduce Tsallis2
def _rtsal2_reduce_ ( pdf ):
    """Reduce Tsallis2"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.pt       () ,                            
                                  pdf.y        () ,                            
                                  pdf.mass     () ,                            
                                  pdf.T        () ,
                                  pdf.q        () ,
                                  pdf.mu       () )

Ostap.Models.Tsallis2.__reduce__ = _rtsal2_reduce_ 

# =============================================================================
## reduce TwoExpos
def _r2expo_reduce_ ( pdf ):
    """Reduce TwoExpos"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.alpha    () ,                            
                                  pdf.delta    () ,
                                  pdf.x0       () )

Ostap.Models.TwoExpos.__reduce__ = _r2expo_reduce_ 

# =============================================================================
## reduce DoubleGauss
def _r2gau_reduce_ ( pdf ):
    """Reduce DoubleGauss"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.sigma    () ,                            
                                  pdf.fraction () ,
                                  pdf.scale    () ,
                                  pdf.mean     () )

Ostap.Models.DoubleGauss.__reduce__ = _r2gau_reduce_ 

# =============================================================================
## reduce Gumbel
def _rgumbel_reduce_ ( pdf ):
    """Reduce Gumbel"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.mu       () ,                            
                                  pdf.beta     () )

Ostap.Models.Gumbel.__reduce__ = _rgumbel_reduce_ 

# =============================================================================
## reduce Weibull
def _rweibull_reduce_ ( pdf ):
    """Reduce Weibull"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.scale    () ,                            
                                  pdf.shape    () ,
                                  pdf.shift    () )

Ostap.Models.Weibull.__reduce__ = _rweibull_reduce_ 

# =============================================================================
## reduce Rice
def _rrice_reduce_ ( pdf ):
    """Reduce Rice"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name         ,
                                  pdf.title        ,
                                  pdf.x         () ,                            
                                  pdf.nu        () ,                            
                                  pdf.varshigma () ,
                                  pdf.shift     () )

Ostap.Models.Rice.__reduce__ = _rrice_reduce_ 

# =============================================================================
## reduce RasingCosine 
def _rrcos_reduce_ ( pdf ):
    """Reduce RaisingCosine"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.mean     () ,                            
                                  pdf.scale    () )

Ostap.Models.RaisingCosine.__reduce__ = _rrcos_reduce_ 

# =============================================================================
## reduce q-Gaussian
def _rqgau_reduce_ ( pdf ):
    """Reduce QGaussian"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.mean     () ,                            
                                  pdf.scale    () ,
                                  pdf.q        () )

Ostap.Models.QGaussian.__reduce__ = _rqgau_reduce_ 

# =============================================================================
## reduce k-Gaussian
def _rkgau_reduce_ ( pdf ):
    """Reduce KGaussian"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.mean     () ,                            
                                  pdf.scale    () ,
                                  pdf.kappa    () )

Ostap.Models.KGaussian.__reduce__ = _rkgau_reduce_ 

# =============================================================================
## reduce Hyperbolic
def _rhyp_reduce_ ( pdf ):
    """Reduce Hyperbolic"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.mu       () ,                            
                                  pdf.sigma    () ,
                                  pdf.zeta     () ,
                                  pdf.kappa    () )

Ostap.Models.Hyperbolic.__reduce__ = _rhyp_reduce_ 

# =============================================================================
## reduce GenHyperbolic
def _rghyp_reduce_ ( pdf ):
    """Reduce GenHyperbolic"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.mu       () ,                            
                                  pdf.sigma    () ,
                                  pdf.zeta     () ,
                                  pdf.kappa    () ,
                                  pdf.lambd    () )

Ostap.Models.GenHyperbolic.__reduce__ = _rghyp_reduce_ 

# =============================================================================
## reduce Das
def _rdas_reduce_ ( pdf ):
    """Reduce Das"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.mu       () ,                            
                                  pdf.sigma    () ,
                                  pdf.kL       () ,
                                  pdf.kR       () )

Ostap.Models.Das.__reduce__ = _rdas_reduce_ 

# =============================================================================
## reduce GenInvGauss
def _rgig_reduce_ ( pdf ):
    """Reduce GenInvGauss"""
    return root_store_factory , ( type ( pdf )  ,
                                  pdf.name      ,
                                  pdf.title     ,
                                  pdf.x      () ,                            
                                  pdf.theta  () ,                            
                                  pdf.eta    () ,
                                  pdf.p      () ,
                                  pdf.shift  () )

Ostap.Models.GenInvGauss.__reduce__ = _rgig_reduce_ 

# =============================================================================
## reduce PositiveSpline 
def _rsplpos_reduce_ ( pdf ):
    """Reduce PositiveSpline"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.spline   () ,                            
                                  pdf.phis     () )

Ostap.Models.PositiveSpline   .__reduce__ = _rsplpos_reduce_ 
Ostap.Models.MonotonicSpline  .__reduce__ = _rsplpos_reduce_ 
Ostap.Models.ConvexSpline     .__reduce__ = _rsplpos_reduce_ 
Ostap.Models.ConvexOnlySpline .__reduce__ = _rsplpos_reduce_ 

# =============================================================================
## reduce HORNSdini
def _rdini_reduce_ ( pdf ):
    """Reduce HORNSdini"""
    return root_store_factory , ( type ( pdf )  ,
                                  pdf.name      ,
                                  pdf.title     ,
                                  pdf.x      () ,                            
                                  pdf.a      () ,                            
                                  pdf.delta  () ,
                                  pdf.phi    () )

Ostap.Models.HORNSdini.__reduce__ = _rdini_reduce_ 
Ostap.Models.HILLdini .__reduce__ = _rdini_reduce_ 

# =============================================================================
## reduce Histo1D
def _rh1d_reduce_ ( pdf ) :
    """reduce Histo1D"""
    return root_store_factory , ( type ( pdf )  ,
                                  pdf.name      ,
                                  pdf.title     ,
                                  pdf.x      () ,                            
                                  pdf.histo  () )

Ostap.Models.Histo1D.__reduce__ = _rh1d_reduce_ 

# =============================================================================
## reduce Histo2D
def _rh2d_reduce_ ( pdf ) :
    """reduce Histo2D"""
    return root_store_factory , ( type ( pdf )  ,
                                  pdf.name      ,
                                  pdf.title     ,
                                  pdf.x      () ,                            
                                  pdf.y      () ,                            
                                  pdf.histo  () )

Ostap.Models.Histo2D.__reduce__ = _rh2d_reduce_ 

# =============================================================================
## reduce Histo3D
def _rh3d_reduce_ ( pdf ) :
    """reduce Histo3D"""
    return root_store_factory , ( type ( pdf )  ,
                                  pdf.name      ,
                                  pdf.title     ,
                                  pdf.x      () ,                            
                                  pdf.y      () ,                            
                                  pdf.z      () ,                            
                                  pdf.histo  () )

Ostap.Models.Histo3D.__reduce__ = _rh3d_reduce_ 

# =============================================================================
## reduce CutoffGauss
def _rcutgau_reduce_ ( pdf ):
    """Reduce CutOffGauss"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.x0       () ,                            
                                  pdf.sigma    () ,
                                  pdf.cutoff   () )

Ostap.Models.CutOffGauss.__reduce__ = _rcutgau_reduce_ 

# =============================================================================
## reduce CutoffStudent
def _rcutstt_reduce_ ( pdf ):
    """Reduce CutOffStudent"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.x0       () ,                            
                                  pdf.nu       () ,
                                  pdf.sigma    () ,
                                  pdf.cutoff   () )

Ostap.Models.CutOffStudent.__reduce__ = _rcutstt_reduce_ 

# =============================================================================
##  2D PDFs
# =============================================================================

# =============================================================================
## reduce Poly2DPositive
def _rpol2d_reduce_ ( pdf ):
    """Reduce Poly2DPositive"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,                            
                                  pdf.nX       () ,                            
                                  pdf.nY       () ,
                                  pdf.phis     () )

Ostap.Models.Poly2DPositive.__reduce__ = _rpol2d_reduce_ 

# =============================================================================
## reduce Poly2DSymPositive
def _rpol2ds_reduce_ ( pdf ):
    """Reduce Poly2DSymPositive"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,                            
                                  pdf.n        () ,                            
                                  pdf.phis     () )

Ostap.Models.Poly2DSymPositive.__reduce__ = _rpol2ds_reduce_ 

# =============================================================================
## reduce PS2DPol
def _rps2dpol_reduce_ ( pdf ):
    """Reduce PS2DPol"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,
                                  pdf.psX      () ,                            
                                  pdf.psY      () ,                            
                                  pdf.nX       () ,                            
                                  pdf.nX       () ,                            
                                  pdf.phis     () )

Ostap.Models.PS2DPol.__reduce__ = _rps2dpol_reduce_ 

# =============================================================================
## reduce PS2DPolSym
def _rps2dpols_reduce_ ( pdf ):
    """Reduce PS2DPolSym"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,
                                  pdf.psX      () ,                            
                                  pdf.n        () ,                            
                                  pdf.phis     () )

Ostap.Models.PS2DPolSym.__reduce__ = _rps2dpols_reduce_ 

# =============================================================================
## reduce PS2DPol2
def _rps2dpol2_reduce_ ( pdf ):
    """Reduce PS2DPol2"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,
                                  pdf.psX      () ,                            
                                  pdf.psY      () ,                            
                                  pdf.mmax     () ,                            
                                  pdf.nX       () ,                            
                                  pdf.nY       () ,                            
                                  pdf.phis     () )

Ostap.Models.PS2DPol2.__reduce__ = _rps2dpol2_reduce_ 
Ostap.Models.PS2DPol3.__reduce__ = _rps2dpol2_reduce_ 

# =============================================================================
## reduce PS2DPol2Sym
def _rps2dpol2s_reduce_ ( pdf ):
    """Reduce PS2DPol2Sym"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,
                                  pdf.psX      () ,                            
                                  pdf.mmax     () ,                            
                                  pdf.nX       () ,                            
                                  pdf.phis     () )

Ostap.Models.PS2DPol2Sym.__reduce__ = _rps2dpol2s_reduce_ 
Ostap.Models.PS2DPol3Sym.__reduce__ = _rps2dpol2s_reduce_ 

# =============================================================================
## Reduce Expo2DPol
def _rexp2d_reduce_ ( pdf ):
    """Reduce Expo2DPol"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,
                                  pdf.taux     () ,                            
                                  pdf.tauy     () ,                            
                                  pdf.nX       () ,                            
                                  pdf.nY       () ,                            
                                  pdf.phis     () )

Ostap.Models.Expo2DPol.__reduce__ = _rexp2d_reduce_ 

# =============================================================================
## Reduce Expo2DPolSym
def _rexp2ds_reduce_ ( pdf ):
    """Reduce Expo2DPolSym"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,
                                  pdf.tau      () ,                            
                                  pdf.nX       () ,                            
                                  pdf.phis     () )

Ostap.Models.Expo2DPolSym.__reduce__ = _rexp2ds_reduce_ 

# =============================================================================
## Reduce ExpoPS2DPol
def _rexpps2d_reduce_ ( pdf ):
    """Reduce ExpoPS2DPol"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,
                                  pdf.tau      () ,                            
                                  pdf.psY      () ,                            
                                  pdf.nX       () ,                            
                                  pdf.nY       () ,                            
                                  pdf.phis     () )

Ostap.Models.ExpoPS2DPol.__reduce__ = _rexpps2d_reduce_ 

# =============================================================================
## reduce Spline2D
def _rspl2d_reduce_ ( pdf ):
    """Reduce Spline2D"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,                            
                                  pdf.spline   () ,                            
                                  pdf.phis     () )

Ostap.Models.Spline2D .__reduce__ = _rspl2d_reduce_ 

# =============================================================================
## reduce Spline2DSym
def _rspl2ds_reduce_ ( pdf ):
    """Reduce Spline2DSym"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,                            
                                  pdf.spline   () ,                            
                                  pdf.phis     () )

Ostap.Models.Spline2DSym.__reduce__ = _rspl2ds_reduce_ 

# =============================================================================
## reduce Gauss2D
def _rgauss2d_reduce_ ( pdf ):
    """Reduce Gauss2D"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,                            
                                  pdf.muX      () ,                            
                                  pdf.muY      () ,                            
                                  pdf.sigmaX   () ,                            
                                  pdf.sigmaY   () ,                            
                                  pdf.theta    () )

Ostap.Models.Gauss2D.__reduce__ = _rgauss2d_reduce_ 

# =============================================================================
## reduce Poly3DPositive
def _rpol3d_reduce_ ( pdf ):
    """Reduce Poly3DPositive"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,                            
                                  pdf.z        () ,                            
                                  pdf.nX       () ,                            
                                  pdf.nX       () ,                            
                                  pdf.nY       () ,                            
                                  pdf.nZ       () ,                            
                                  pdf.phis     () )

Ostap.Models.Poly3DPositive.__reduce__ = _rpol3d_reduce_ 

# =============================================================================
## reduce Poly3DSymPositive
def _rpol3ds_reduce_ ( pdf ):
    """Reduce Poly3DSymPosiitive"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,                            
                                  pdf.z        () ,                            
                                  pdf.nX       () ,                            
                                  pdf.phis     () )

Ostap.Models.Poly3DSymPositive.__reduce__ = _rpol3ds_reduce_ 

# =============================================================================
## reduce Poly3DMixPositive
def _rpol3dm_reduce_ ( pdf ):
    """Reduce Poly3DMixPosiitive"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,                            
                                  pdf.z        () ,                            
                                  pdf.nX       () ,                            
                                  pdf.nZ       () ,                            
                                  pdf.phis     () )

Ostap.Models.Poly3DMixPositive.__reduce__ = _rpol3dm_reduce_ 

# =============================================================================
## reduce Gauss3D
def _rgauss3d_reduce_ ( pdf ):
    """Reduce Gauss3D"""
    return root_store_factory , ( type ( pdf )    ,
                                  pdf.name        ,
                                  pdf.title       ,
                                  pdf.x        () ,                            
                                  pdf.y        () ,                            
                                  pdf.z        () ,                            
                                  pdf.muX      () ,                            
                                  pdf.muY      () ,                            
                                  pdf.muZ      () ,                            
                                  pdf.sigmaX   () ,                            
                                  pdf.sigmaY   () ,                            
                                  pdf.sigmaZ   () ,                            
                                  pdf.phi      () ,
                                  pdf.theta    () ,
                                  pdf.psi      () )

Ostap.Models.Gauss3D.__reduce__ = _rgauss3d_reduce_ 


# =============================================================================
## reduce Ostap::Functions::FuncRooTH1 
def _rfth1_reduce_ ( fun ):
    """Reduce Ostap.Functions.FuncTH1"""
    return root_factory , ( type ( fun ) ,
                            fun.histo () ,
                            fun.x     () )
# =============================================================================
## reduce Ostap::Functions::FuncRooTH3 
def _rfth2_reduce_ ( fun ):
    """Reduce Ostap.Functions.FuncRooTH2"""
    return root_factory , ( type ( fun ) ,
                            fun.histo () ,
                            fun.x     () , 
                            fun.y     () )
# =============================================================================
## reduce Ostap::Functions::FuncRooTH3
def _rfth3_reduce_ ( fun ):
    """Reduce Ostap.Functions.FuncRooTH3"""
    return root_factory , ( type ( fun ) ,
                            fun.histo () ,
                            fun.x     () ,
                            fun.y     () ,
                            fun.z     () )
# =============================================================================
## reduce Ostap::Functions::FuncRooFormula
def _rfff_reduce_ ( fun ):
    """Reduce Ostap.Functions.FuncRooFormula"""
    return root_factory , ( type ( fun ) , fun.expression() )

Ostap.Functions.FuncRooTH1     . __reduce__ = _rfth1_reduce_
Ostap.Functions.FuncRooTH2     . __reduce__ = _rfth2_reduce_
Ostap.Functions.FuncRooTH3     . __reduce__ = _rfth3_reduce_
Ostap.Functions.FuncRooFormula . __reduce__ =  _rfff_reduce_




# =============================================================================

_decorated_classes_ = (
    ## ROOT/RooFit classes 
    ROOT.RooBinning                    , 
    ROOT.RooUniformBinning             , 
    ROOT.RooRangeBinning               , 
    ROOT.RooArgSet                     , 
    ROOT.RooArgList                    , 
    ROOT.RooGaussian                   , 
    ROOT.RooMultiVarGaussian           , 
    ROOT.RooAddPdf                     , 
    ROOT.RooProdPdf                    , 
    ROOT.RooFFTConvPdf                 , 
    ROOT.RooSimultaneous               , 
    ROOT.RooEfficiency                 , 
    ROOT.RooPolyVar                    , 
    ROOT.RooPolynomial                 , 
    ROOT.RooFitResult                  ,
    ROOT.RooPlot                       ,
    ROOT.RooCurve                      ,
    ROOT.RooEllipse                    ,
    ROOT.RooHist                       ,
    ## Ostap classes 
    Ostap.MoreRooFit.TwoVars           , 
    Ostap.MoreRooFit.Addition          , 
    Ostap.MoreRooFit.Addition2         , 
    Ostap.MoreRooFit.Subtraction       , 
    Ostap.MoreRooFit.Product           , 
    Ostap.MoreRooFit.ProductPdf        , 
    Ostap.MoreRooFit.Id                , 
    Ostap.MoreRooFit.AddDeps           , 
    Ostap.MoreRooFit.Combination       , 
    Ostap.MoreRooFit.Asymmetry         , 
    Ostap.MoreRooFit.Constant          , 
    Ostap.MoreRooFit.Bernstein         , 
    Ostap.MoreRooFit.Monotonic         , 
    Ostap.MoreRooFit.Convex            , 
    Ostap.MoreRooFit.ConvexOnly        , 
    Ostap.MoreRooFit.BSpline           ,
    ##
    Ostap.Models.Uniform               ,
    ## BW & friends 
    Ostap.Models.BreitWigner           , 
    Ostap.Models.BreitWignerMC         , 
    Ostap.Models.BWI                   , 
    Ostap.Models.Flatte                ,
    Ostap.Models.LASS                  ,
    Ostap.Models.BWPS                  ,
    Ostap.Models.BW3L                  ,
    ## others 
    Ostap.Models.Uniform               ,
    Ostap.Models.Voigt                 , 
    Ostap.Models.PseudoVoigt           , 
    Ostap.Models.CrystalBall           , 
    Ostap.Models.CrystalBallRS         , 
    Ostap.Models.CrystalBallDS         , 
    Ostap.Models.Needham               ,
    Ostap.Models.Apollonios            , 
    Ostap.Models.Apollonios2           , 
    Ostap.Models.BifurcatedGauss       , 
    Ostap.Models.GenGaussV1            , 
    Ostap.Models.GenGaussV2            , 
    Ostap.Models.SkewGauss             , 
    Ostap.Models.ExGauss               , 
    Ostap.Models.NormalLaplace         , 
    Ostap.Models.Novosibirsk           , 
    Ostap.Models.Bukin                 , 
    Ostap.Models.StudentT              , 
    Ostap.Models.BifurcatedStudentT    , 
    Ostap.Models.PearsonIV             ,
    Ostap.Models.SkewGenT              ,
    Ostap.Models.GramCharlierA         , 
    Ostap.Models.PhaseSpace2           , 
    Ostap.Models.PhaseSpaceLeft        , 
    Ostap.Models.PhaseSpaceRight       , 
    Ostap.Models.PhaseSpaceNL          , 
    Ostap.Models.PhaseSpacePol         , 
    Ostap.Models.PhaseSpaceLeftExpoPol , 
    ##
    Ostap.Models.PolyPositive          , 
    Ostap.Models.PolyPositiveEven      , 
    Ostap.Models.PolyMonotonic         , 
    Ostap.Models.PolyConvex            , 
    Ostap.Models.PolyConvexOnly        , 
    Ostap.Models.ExpoPositive          , 
    Ostap.Models.PolySigmoid           , 
    Ostap.Models.TwoExpoPositive       , 
    Ostap.Models.GammaDist             , 
    Ostap.Models.LogGammaDist          , 
    Ostap.Models.Log10GammaDist        , 
    Ostap.Models.GenGammaDist          , 
    Ostap.Models.Amoroso               , 
    Ostap.Models.LogGamma              , 
    Ostap.Models.BetaPrime             , 
    Ostap.Models.Landau                , 
    Ostap.Models.SinhAsinh             , 
    Ostap.Models.JohnsonSU             , 
    Ostap.Models.Atlas                 , 
    Ostap.Models.Sech                  , 
    Ostap.Models.Logistic              ,
    Ostap.Models.Hat                   , 
    Ostap.Models.Up                    , 
    Ostap.Models.Slash                 , 
    Ostap.Models.Argus                 , 
    Ostap.Models.GenArgus              , 
    Ostap.Models.Losev                 , 
    Ostap.Models.AsymmetricLaplace     , 
    Ostap.Models.FupN                  , 
    Ostap.Models.Tsallis               , 
    Ostap.Models.QGSM                  , 
    Ostap.Models.TwoExpos              , 
    Ostap.Models.DoubleGauss           , 
    Ostap.Models.Gumbel                , 
    Ostap.Models.Weibull               , 
    Ostap.Models.QGaussian             , 
    Ostap.Models.KGaussian             , 
    Ostap.Models.Hyperbolic            , 
    Ostap.Models.GenHyperbolic         , 
    Ostap.Models.Das                   , 
    Ostap.Models.GenInvGauss           , 
    Ostap.Models.PositiveSpline        , 
    Ostap.Models.MonotonicSpline       , 
    Ostap.Models.ConvexSpline          , 
    Ostap.Models.ConvexOnlySpline      , 
    Ostap.Models.HORNSdini             , 
    Ostap.Models.HILLdini              , 
    Ostap.Models.Histo1D               , 
    Ostap.Models.Histo2D               , 
    Ostap.Models.Histo3D               , 
    Ostap.Models.CutOffGauss           , 
    Ostap.Models.CutOffStudent         , 
    ## 2D models 
    Ostap.Models.Poly2DPositive        , 
    Ostap.Models.Poly2DSymPositive     , 
    Ostap.Models.PS2DPol               , 
    Ostap.Models.PS2DPol2              , 
    Ostap.Models.PS2DPol3              , 
    Ostap.Models.PS2DPolSym            , 
    Ostap.Models.PS2DPol2Sym           , 
    Ostap.Models.PS2DPol3Sym           , 
    Ostap.Models.Expo2DPol             , 
    Ostap.Models.Expo2DPolSym          , 
    Ostap.Models.ExpoPS2DPol           , 
    Ostap.Models.Spline2DSym           , 
    Ostap.Models.Gauss2D               , 
    Ostap.Models.Tsallis2              , 
    ## 3D models 
    Ostap.Models.Poly3DPositive        , 
    Ostap.Models.Poly3DSymPositive     , 
    Ostap.Models.Poly3DMixPositive     , 
    Ostap.Models.Gauss3D               ,
    ##     
    Ostap.Functions.FuncRooTH1         , 
    Ostap.Functions.FuncRooTH2         , 
    Ostap.Functions.FuncRooTH3         , 
    Ostap.Functions.FuncRooFormula     ,  
    )

for t in _decorated_classes_ :
    _new_methods_.append ( t.__reduce__  )

_new_methods_ = tuple ( _new_methods_ ) 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
