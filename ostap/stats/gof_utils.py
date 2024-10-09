#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/gof_utils.py
#  Set of utulities for goodness-of-fit studies 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2023-12-06
# =============================================================================
""" Simple utulities for goodness-of-fit studies 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2023-12-06"
__all__     = (
    'mean_var'  , ## mean and variance for (weighted) arrays
    'nEff'      , ## get number of effective entries
    'normalize' , ## "normalize" variables in dataset/structured array
    )
# =============================================================================
from   ostap.core.meta_info     import python_info 
from   ostap.core.core          import VE, Ostap
from   ostap.core.ostap_types   import string_types
from   ostap.stats.counters     import EffCounter
from   ostap.utils.basic        import numcpu
from   ostap.utils.utils        import splitter 
from   ostap.utils.progress_bar import progress_bar
import ROOT, sys, warnings  
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import numpy as np
    _np_floats = np.float16, np.float32, np.float64, np.float128
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    np = None
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.gof_utils' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Simple utilities for goodness-of-fit studies' )
# ============================================================================
## Get the mean and variance for (1D) data array with optional (1D) weight array
#  @code
#  ds = ... ## dataste as structured array
#  mean, cov2 = mean_var ( ds ['x'] )
#  @endcode
#  - with weight 
#  @code
#  ds = ... ## dataset as structured array with weight 
#  mean, cov2 = mean_var ( ds ['x'] , ds['weight'] )
#  @endcode 
def mean_var ( data , weight = None ) :
    """ Get the mean and variance for 1D-data array with optional 1D-weight array

    >>> ds = ... ## dataset as structured array
    >>> mean, cov2 = mean_var ( ds ['x'] )
    
    - with weight 
    
    >>> ds = ... ## dataset as structured array with weight 
    >>> mean, cov2 = mean_var ( ds ['x'] , ds['weight'] )
    """
    #
    if weight is None :
        mean = np.mean ( data , axis = 0 , dtype = np.float64 ) 
        var  = np.var  ( data , axis = 0 , dtype = np.float64 ) 
        return mean , var
    # 
    mean  = np.average (   data               , weights = weight , axis = 0 )
    var   = np.average ( ( data - mean ) ** 2 , weights = weight , axis = 0 )
    #
    return mean , var 

# =============================================================================
## Get the effectibe number of entries for 1D-array
#  \f{ n_{eff} = \frac{  \left\langle x \right\rangle^2}
#                     { \left\langle x^2 \right\rangle } \f}
def nEff ( weights ) :
    """ Get the effective number of entries for 1D-array
    n_eff = ( sum ( x )  ) ^2 / sum ( x^2 )
    """

    s1 = np.sum ( weights      , dtype = np.float64 )
    s2 = np.sum ( weights ** 2 , dtype = np.float64 )
    
    return s1 * s1 / s2

# =============================================================================
## Get the "normalized" input datasets
#  All floating felds  are calculated as
#  \f[ x = \frac{x - \left\langle x \right\rangle}{\sigma} \f]
#  where \f$ \left\langle x \right\rangle\f$ is mena value
#  and \f$ \sigma \f$ is a standard deviation.
# 
#  @code
#  ds      = ... # data set as structured array
#  dsn, J  = normalize ( ds ) 
#  @endcode
#
#  - If several datasets are specified, all floating names must be the same
#  and the mean and sigma are either taken either from the first dataset,
#  if <code>first=True</code> or as combined through all datasets otherwise 
#
#  @code
#  ds1 = ... # data set as structured array
#  ds2 = ... # data set as structured array
#  ds3 = ... # data set as structured array
#  ds1n, ds2n, ds3n = normalize ( ds1 , ds2 , ds3 , first = True )
#  @endcode
#
#  - If <code>weight</code> is specified, this floating column is considered
#  as the weight
#
#  @code
#  ds  = ... # data set as structured array with weight 
#  dsn = normalize ( ds , weight = 'weight' ) 
#  @endcode
#
#  @code
#  ds1 = ... # data set as structured array without weight 
#  ds2 = ... # data set as structured array with weight 
#  datasets    = normalize ( ds1 , ds2 , weight = ( None , 'weight'  ) )
#  ds1n , ds2n = datasets 
#  @endcode
#
#  @attention Only the floating point columns are transformed! 
#  @attention Input datasets are expected to be numpy structured arrays
#
#  @code
#  ds = ... # data set as structured array
#  dsn = normalize ( ds ) 
#  @endcode
def normalize2 ( datasets , weight = () , first = True ) :
    """ Get the `normalized' input datasets
    All floating fields  are calculated as
    
    x = (x - <x>)/sigma
    
    - <x> is a mean value 
    - is a standard deviation.
    
    - If several datasets are specified, all floating names must be the same
    and the mean and sigma are either taken either from the first dataset,
    if `first=True` or as combined through all datasets, otherwise  
    
    - If `weight` is specified, this floating column is concidered
    as the weight
    
    - attention Only the floating point columns are transformed! 
    - attention Input datasets are expected to be numpy structured arrays 
    """
    
    nd = len ( datasets ) 
    if not weight                             : weight = nd * [ ''     ]
    elif isinstance ( weight , string_types ) : weight = nd * [ weight ]
    
    assert ( len ( weight ) == nd ) and \
           all ( ( not w ) or isinstance ( w , string_types ) for w in weight ) , \
           'Invalid specification of weight!'
        
    weight = list ( weight )
    for i , w in enumerate ( weight ) :
        if not w : weight [ i ] = '' 
    weight = tuple ( weight )

    ds     = datasets [ 0  ]
    others = datasets [ 1: ]
    
    ## collect the floating columns 
    columns = []
    w0      = weight [ 0 ] 
    for n,t in ds.dtype.fields.items () :
        if t[0] in _np_floats  and n != w0 : columns.append ( n ) 
        
    vmeans  = [] 
    for i , c in enumerate ( columns ) :
        mean, var    = mean_var ( ds [ c ] , None if not w0 else ds [ w0 ] )
        vmeans.append ( VE ( mean , var ) )
        
    ## Number of events/effective entries 
    nevents = 1.0 * ds.shape [ 0 ] if not w0 else nEff ( ds [ w0 ] )
    
    if not first and others : 
        nevents = ds.shape[0] 
        for k , dd in enumerate ( others ) :
            
            wk = weight [ k + 1 ] 
            nn = 1.0 * dd.shape [ 0 ] if not wk else nEff ( dd [ wk ] )                
            
            for i , c in enumerate ( columns ) :
                
                mean, var = mean_var ( dd [ c ] , None if not wk else dd [ wk ] )                    
                vv = VE ( mean , var ) 
                vmeans [ i ] = Ostap.Math.two_samples ( vmeans [ i ] , nevents , vv , nn )
                
            nevents += nn
                
    result = []  
    for d in datasets :
        
        nds = d.copy ()
        for ic , c in enumerate ( columns ) :
            vv        = vmeans [ ic ]
            mean      = vv.value ()
            sigma     = vv.error ()                 
            a         = nds [ c ]
            nds [ c ] =  ( a - mean ) / sigma
            
        result.append ( nds )
            
    return tuple ( result ) 

# =============================================================================
if ( 3 , 0 ) <= sys.version_info : 
    code3 = """
def normalize ( ds , *others , weight = () , first = True ) :
    result = normalize2 ( ( ds , *others ) , weight = weight , first = first )      
    return result [ 0 ] if not others else result 
    """
    exec ( code3 )
else :
    code2 = """
def normalize ( ds , others = () , weight = () , first = True ) :
    datasets = [ ds ] + [ d for d in others ]
    result = normalize2 ( datasets , weight = weight , first = first )          
    return result [ 0 ] if not others else result 
    """
    exec ( code2 )
normalize.__doc__ = normalize2.__doc__ 


# =============================================================================
jl = None 
# =============================================================================
try : # =======================================================================
    # =========================================================================
    if ( 3 , 0 ) <= python_info :
        with warnings.catch_warnings(): 
            warnings.simplefilter ( "ignore" , category = DeprecationWarning  )
            import joblib as jl
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    jl = None
# =============================================================================
## @class PERMUTATOR
#  Helper class that allow to run permutattion test in parallel 
class PERMUTATOR(object) :
    """ Helper class that allow to run permutation test in parallel 
    """
    def __init__ ( self, gof, t_value , ds1 , ds2 ) :
        self.gof     = gof
        self.ds1     = ds1
        self.ds2     = ds2
        self.t_value = t_value
        
    # =========================================================================
    ## run N-permutations 
    def __call__ ( self , N , silent = True ) :
        np.random.seed()
        n1      = len ( self.ds1 )
        pooled  = np.concatenate ( [ self.ds1 , self.ds2 ] )
        counter = EffCounter()
        for i in progress_bar ( N , silent = silent , description = 'Permutations:') : 
            np.random.shuffle ( pooled )            
            tv       = self.gof.t_value ( pooled [ : n1 ] , pooled [ n1: ] )
            counter += bool ( self.t_value < tv  )
            
        del pooled
        return counter

# =============================================================================
jl = None 
# =============================================================================
if ( 3 , 0 ) <= python_info : # ===============================================
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        with warnings.catch_warnings(): 
            warnings.simplefilter ( "ignore" , category = DeprecationWarning  )
            import joblib as jl
        # =====================================================================
        ## Run NN-permutations in parallel using joblib 
        def joblib_run ( self , NN , silent = True ) :
            """ Run NN-permutations in parallel using joblib """
            nj    = 3 * numcpu () + 4
            lst   = splitter ( NN , nj )
            ## 
            conf  = { 'n_jobs' : -1 , 'verbose' : 0 }
            if    '1.3.0' <= jl.__version__ < '1.4.0' : conf [ 'return_as' ] = 'generator'           
            elif  '1.4.0' <= jl.__version__           : conf [ 'return_as' ] = 'unordered_generator'
            ##
            input   = ( jl.delayed (self)( N ) for N in lst )
            counter = EffCounter()
            with warnings.catch_warnings(): 
                warnings.simplefilter ( "ignore" ) ## , category = DeprecationWarning  )
                results = jl.Parallel ( **conf ) ( input )
                for r in progress_bar ( results , max_value = nj , silent = silent , description = 'Permutations:') :
                    counter += r
            ## 
            return counter
        # =====================================================================
        PERMUTATOR.run = joblib_run        
        # =====================================================================
        logger.debug ( 'Joblib will be  used foe parallel permuations')
        # =====================================================================        
    except ImportError : # ====================================================
        # =====================================================================
        jl = None
        
# =============================================================================
if not jl : # =================================================================
    # =========================================================================
    ## Run NN-permutations in parallel using WorkManager
    def pp_run ( self , NN , silent = True ) :
        """ Run NN-permutations in parallel using WorkManager"""
        nj    = 3 * numcpu () + 4
        lst   = splitter ( NN , nj )
        ##
        from ostap.parallel.parallel import WorkManager
        manager = WorkManager ( silent = silent )
        counter = EffCounter()
        ## 
        ## use the bare interface 
        for result in manager.iexecute ( self , lst , progress = not silent  , njobs = nj ) :
            counter += result 
        # 
        return counter
    # =========================================================================
    logger.debug ( 'Parallel will be  used for parallel permuations')
    # =====================================================================        
    PERMUTATOR.run = pp_run
    # =========================================================================

# =============================================================================
## @class TOYS
#  Helper class to tun toys for Goodness-of-Fit studies 
class TOYS(object) :
    """ Helper class that allow to run permutation test in parallel 
    """
    def __init__ ( self, gof , t_value , pdf , Ndata , sample = False ) :
        self.gof     = gof
        self.pdf     = pdf
        self.Ndata   = Ndata 
        self.t_value = t_value
        self.sample  = gof.sample 
        self.silent  = gof.silent
        
    # =========================================================================
    ## run N-toys 
    def __call__ ( self , N , silent = True ) :
        """ Run N-toys """
        ROOT.gRandom.SetSeed() 
        ROOT.RooRandom.randomGenerator().SetSeed()
        counter = EffCounter()
        for i in range ( N ) :
            
            ds = self.pdf.generate ( self.Ndata , sample = self.sample )
            tv = self.gof ( self.pdf , ds )
            counter += bool ( self.t_value < tv  )
           
            if isinstance ( ds , ROOT.RooDataSet ) :
                ds = Ostap.MoreRooFit.delete_data ( ds )
            del ds
            
        return counter 

    # =========================================================================
    ## Run N-toys in parallel using WorkManager
    def run ( self , NN , silent = False ) :
        """ Run NN-permutations in parallel using WorkManager"""
        nj    = 3 * numcpu () + 4
        lst   = splitter ( NN , nj )
        ##
        from ostap.parallel.parallel import WorkManager
        manager = WorkManager ( silent = silent )
        counter = EffCounter()
        ## 
        ## use the bare interface 
        for result in manager.iexecute ( self , lst , progress = not silent  , njobs = nj ) :
            counter += result 
        # 
        return counter
    
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    if not np : logger.warning ("Numpy  is not avalaille!") 
    if not jl : logger.warning ("Joblib is not avalaille!") 

# =============================================================================
##                                                                      The END 
# =============================================================================


    
