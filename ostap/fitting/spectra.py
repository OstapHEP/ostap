#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file spectra.py
#
# A set of predefined ready-to-use
# functions and shapes for fitting of pt-spectra 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
# 
# =============================================================================
"""A set of predefined ready-to-use shapes and PDFs
for fitting pT-spectra"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2015-07-11"
__all__     = (
    #
    'Tsallis_pdf'   , ## Tsallis   PDF 
    'QGSM_pdf'      , ## QGSM      PDG 
    'GammaDist_pdf' , ## GammaDist PDF
    #
    'TsallisFun'    , ## Tsallis   function object 
    'QGSMFun'       , ## QGSM      function object 
    'GammaDistFun'  , ## GammaDist function object 
    #
    'tsallisTF1'    , ## Tsallis   create TF1
    'qgsmTF1'       , ## QGSM      create TF1 
    'gammaDistTF1'  , ## GammaDist create TF1 
    #
    )
# =============================================================================
import ROOT, math
# =============================================================================
from   ostap.logger.logger     import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.spectra' )
else                       : logger = getLogger ( __name__           )
# =============================================================================
from   ostap.core.core  import cpp, Ostap, VE , funID 
# =============================================================================
models = [] 
# =============================================================================
from ostap.fitting.distributions import   Tsallis_pdf
models.append (   Tsallis_pdf )
from ostap.fitting.distributions import   QGSM_pdf
models.append (      QGSM_pdf )
from ostap.fitting.distributions import   GammaDist_pdf
models.append ( GammaDist_pdf )

# ==============================================================================
## @class PtFitBase
#  helper object  
class PtFitBase(object) :
    """Helper class fro fitting pT-spectra
    """
    def __init__ ( self , ptmax , ptmin = 0 ) :
        self._integral = -1
        self._ptmin    = ptmin
        self._ptmax    = ptmax
        if ptmax <= ptmin :
            raise ArrtibuteError("PtFitBase: wrong ptmin/ptmax: %s/%s" % ( ptmin ,ptmax ) )

    ## get the mean-value 
    def mean     ( self ) : return self._fun.mean     ( self._ptmin , self._ptmax )
    ## get variance
    def variance ( self ) : return self._fun.variance ( self._ptmin , self._ptmax )
    ## get rms 
    def rms      ( self ) : return self._fun.rms      ( self._ptmin , self._ptmax )

    ## make final evaluation of the result 
    def result   ( self , pt , norm , changed  ) :
        ## recalculate integral (if needed) 
        if changed or 0 >= self._integral :
            self._integral = self._fun.integral ( self._ptmin , self._ptmax ) 
            
        return norm * self._fun ( pt ) / self._integral 
    
# ==============================================================================
## @class TsallisFun
#  helper object  
class TsallisFun(PtFitBase) :
    """Helper object  to wrap Tsallis function
    """ 
    def __init__ ( self , ptmax , ptmin = 0 ) :        
        ## initialize the base 
        PtFitBase.__init__ ( self , ptmax , ptmin ) 
        ## create Tsallis function
        self._fun = Ostap.Math.Tsallis(0,10,1.1)
        
    def __call__ ( self , x , pars ) :
        
        pt   = x    [ 0 ]
        
        norm = pars [ 0 ]  ## normalization (almost arbitrary) 
        m0   = pars [ 1 ]  ## 
        n    = pars [ 2 ]
        T    = pars [ 3 ]
        
        changed = False 
        if self._fun.setMass ( m0 ) : changed = True 
        if self._fun.setN    ( n  ) : changed = True 
        if self._fun.setT    ( T  ) : changed = True 

        return self.result ( pt , norm , changed ) 

# ==============================================================================
## @class QGSMFun
#  helper object  
class QGSMFun(PtFitBase) :
    """Helper object  to wrap QGSM-function
    """
    def __init__ ( self , ptmax , ptmin = 0 ) :
        ## initialize the base 
        PtFitBase.__init__ ( self , ptmax , ptmin ) 
        ## create QGSM function 
        self._fun = Ostap.Math.QGSM(0,1)
        
    ## get the mean-value 
    def mean     ( self ) :
        return self._fun.mean ( self._ptmin , self._ptmax )
    
    def __call__ ( self , x , pars ) :
        
        pt = x    [ 0 ]
        
        norm = pars [ 0 ]  ## normalization (almost arbitrary) 
        m0   = pars [ 1 ]  ## 
        b    = pars [ 2 ]
        
        changed = False 
        if self._fun.setMass ( m0 ) : changed = True 
        if self._fun.setB    ( b  ) : changed = True 

        return self.result ( pt , norm , changed ) 

# ==============================================================================
## @class GammaDistFun
#  helper object  
class GammaDistFun(PtFitBase) :
    """Helper object  to wrap Gamma-distribution
    """ 
    def __init__ ( self , ptmax , ptmin = 0 ) :
        ## initialize the base 
        PtFitBase.__init__ ( self , ptmax , ptmin ) 
        ## create QGSM function 
        self._fun = Ostap.Math.GammaDist(2,1)
        
    def __call__ ( self , x , pars ) :
        
        pt = x    [ 0 ]
        
        norm  = pars [ 0 ]  ## normalization (almost arbitrary) 
        k     = pars [ 1 ]  ## 
        theta = pars [ 2 ]
        
        changed = False 
        if self._fun.setK     ( k     ) : changed = True 
        if self._fun.setTheta ( theta ) : changed = True 
        
        return self.result ( pt , norm , changed ) 


# =============================================================================
## create Tsallis function (as TF1) for fitting pT-spectra
#  @code
#  >>> histo =
#  >>> fun   = tsallisTF1 ( ptmax = 50 )
#  >>> histo.Fit( fun , 'SI' )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def tsallisTF1( ptmax , ptmin = 0 , mass = None , name = '' ) :
    """Create Tsallis function (as TF1) for fitting pt-spectra:
    
    >>> histo =
    >>> fun   = tsallisFun ( ptmax = 50 )
    >>> histo.Fit( fun , 'SI' )    
    """
    if ptmax <= ptmin :
        raise ArrtibuteError("Tsallis/TF1: wrong ptmin/ptmax: %s/%s" % ( ptmin ,ptmax ) )
    
    obj       = TsallisFun ( ptmax , ptmin )
    if not name : name =  funID() 
    func      = ROOT.TF1( name , obj , ptmin , ptmax , 4 )
    func._obj = obj
    func.SetParNames (
        'norm' , ## 0
        'm'    , ## 1 
        'n'    , ## 2
        'T'    ) ## 3
    ##
    if isinstance ( mass , float ) and 0 <= mass : func.FixParameter ( 1 , mass )
    else                                         : func.SetParameter ( 1 , 0    )
    ## 
    func.SetParameter ( 0 , 1   ) ## norm 
    func.SetParameter ( 2 , 10  )
    func.SetParameter ( 3 , 1.1 )
    ##
    return func

# =============================================================================
## create QGSM function (as TF1) for fitting pT-spectra
#  @code
#  >>> histo =
#  >>> fun   = qgsmTF1 ( ptmax = 50 )
#  >>> histo.Fit( fun , 'SI' )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def qgsmTF1 ( ptmax , ptmin = 0 , mass = None , name = '') :
    """Create QGSM functin (as TF1) for fitting pt-spectra
    
    >>> histo =
    >>> fun   = qgsmTF1 ( ptmax = 50 )
    >>> histo.Fit( fun , 'SI' )
    """
    if ptmax <= ptmin :
        raise ArrtibuteError("QGSM/TF1: wrong ptmin/ptmax: %s/%s" % ( ptmin ,ptmax ) )
    
    obj       = QGSMFun  (                 ptmax , ptmin     )
    if not name : name =  funID() 
    func      = ROOT.TF1 ( name , obj , ptmin , ptmax , 3 )
    func._obj = obj
    func.SetParNames (
        'norm' , ## 0
        'm'    , ## 1 
        'b'    ) ## 2
    ##
    if isinstance ( mass , float ) and 0 <= mass : func.FixParameter ( 1 , mass )
    else                                         : func.SetParameter ( 1 , 0    )
    ## 
    func.SetParameter ( 0 , 1 )  ## norm 
    func.SetParameter ( 2 , 1 )
    ##
    return func

# =============================================================================
## create GammaDist function (as TF1) for fitting pT-spectra
#  @code
#  >>> histo =
#  >>> fun   = gammaDistTF1 ( ptmax = 50 )
#  >>> histo.Fit( fun , 'SI' )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-07-11
def gammaDistTF1 ( ptmax , ptmin = 0 , mass = None , name = '' ) :
    """Create GammaDist function (as TF1) for fitting pt-spectra
    
    >>> histo =
    >>> fun   = gammaDistTF1 ( ptmax = 50 )
    >>> histo.Fit( fun , 'SI' )
    """
    if ptmax <= ptmin :
        raise ArrtibuteError("GammaDist/TF1: wrong ptmin/ptmax: %s/%s" % ( ptmin ,ptmax ) )
    
    obj       = GammaDistFun (                 ptmax , ptmin     )
    if not name : name =  funID() 
    func      = ROOT.TF1     ( name , obj , ptmin , ptmax , 3 )
    func._obj = obj
    func.SetParNames (
        'norm'  , ## 0
        'k'     , ## 1 
        'theta' ) ## 2
    ##
    func.SetParameter ( 0 , 1 )  ## norm 
    func.SetParameter ( 1 , 2 )
    func.SetParameter ( 2 , 2 )
    ##
    return func

    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )
 
# =============================================================================
# The END 
# =============================================================================
