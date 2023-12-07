#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/gof.py
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
    'nEff'      , ## get number of effective entries
    'mean_var'  , ## mean and variance for (weighted) arrays
    'normalize' , ## "normalize" variables in dataset/structured array
    )
# =============================================================================
from   ostap.core.core        import VE, Ostap
from   ostap.core.ostap_types import string_types
import sys 
# =============================================================================
try :    
    import numpy as np
    _np_floats = np.float16, np.float32, np.float64, np.float128 
except ImportError :
    np = None
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.gof' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug( 'Simple utilities for goodness-of-fit studies')


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
    """Get the mean and variance for 1D-data array with optional 1D-weight array

    >>> ds = ... ## dataste as structured array
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
    """Get the effectibe number of entries for 1D-array

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
#  ds = ... # data set as structured array
#  dsn = normalize ( ds ) 
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
#  ds1n, ds2n, ds3n  = normalize ( ds1 , ds2 , ds3 , first = True ) 
#  @endcode
#
#  - If <code>weight</code> is specified, this floating column is considered
#  as the weight
#
#  @code
#  ds = ... # data set as structured array with weight 
#  dsn = normalize ( ds , weight = 'weight' ) 
#  @endcode
#
#  @code
#  ds1 = ... # data set as structured array without weight 
#  ds2 = ... # data set as structured array with weight 
#  ds1n , ds2n = normalize ( ds1 , ds2 , weight = ( None , 'weight'  ) ) 
#  @endcode
#
#  @attention Only the floating point columns are transformed! 
#  @attention Input datasets are expected to be numpy structured arrays
#
#  @code
#  ds = ... # data set as structured array
#  dsn = normalize ( ds ) 
#  @endcode

if (3,0) <= sys.version_info :
    
    def normalize ( ds , *others , weight = () , first = True ) :
        """ Get the `normalized' input datasets
        All floating felds  are calculated as
        
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

        nd = len ( others ) + 1 
        if not weight                             : weight = nd * [ ''     ]
        elif isinstance ( weight , string_types ) : weight = nd * [ weight ]
            
        assert ( len ( weight ) == nd ) and \
               all ( ( not w ) or isinstance ( w , string_types ) for w in weight ) , \
               'Invalid specification of weight!'
        
        weight = list ( weight )
        for i , w in enumerate ( weight ) :
            if not w : weight [ i ] = '' 
        weight = tuple ( weight )
        
        data    = ( ds , ) + others  
        result  = []  
    
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
                    vmean  [ i ] = Ostap.Math.two_samples ( vmean [ i ] , nevents , vv , nn )
                    
                nevents += nn
                
        result = []  
        for d in ( ds , *others ) :
            
            nds = d.copy ()
            for ic , c in enumerate ( columns ) :
                vv        = vmeans [ ic ]
                mean      = vv.value ()
                sigma     = vv.error ()                 
                a         = nds [  c ]
                nds [ c ] =  ( a - mean ) / sigma

            result.append ( nds )
            
        return result [ 0 ] if not others else tuple ( result ) 
        

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END 
# =============================================================================


    
