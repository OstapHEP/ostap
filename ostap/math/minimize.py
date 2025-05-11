#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/minimize.py
#  Module with some useful utilities for minimization of scalar functiom
#  - a kind of replacement for scipy.minimize.minimize_scalar when scipy is not accessible
#  - the actual code is copied from scipy.minimize  0.18.11
#
#  The main entry point is a function <code>minimizescalar</code>.
#  - a copy from scipy 0.18.11 
# =============================================================================
""" Module with some useful utilities for minimization of scalar functiom
- a kind of replacement for scipy.minimize.minimize_scalar when scipy is not accessible
- the actual code is copied from scipy.minimize  0.18.11

The main entry point is a function <code>minimize_scalar</code>.
- a copy from scipy 0.18.11
"""
# =============================================================================
__version__ = "$Revision:$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2018-10-05"
__all__     = (
    'minimize_scalar' , ## the main entry
    )
# =============================================================================
import warnings, math
from   ostap.math.base import scipy, numpy  
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.minimize' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
if scipy : # ==================================================================
    # =========================================================================
    minimize_scalar = scipy.optimize.minimize_scalar 
    # =========================================================================
else : # ======================================================================
    # =========================================================================
    from ostap.math.local_minimize import scalar_minimize as minimize_scalar 

# =============================================================================
if numpy and scipy : # ========================================================
    # =========================================================================
    ## get a minimum for 1D-function
    #  @code
    #  model = ...
    #  x = model.minimum() 
    #  @endcode 
    def sp_minimum_1D ( fun , xmin , xmax , x0 = None , *args ) :
        """ Get a minimum for 1D-function
        >>> model = ...
        >>> x = model.minimum () 
        >>>
        """
        if x0 == None : x0 = 0.5 * ( xmin + xmax )
        
        x0     = numpy.array ( [ x0 ] )
        
        bounds = [ ( xmin , xmax ) ]
        
        import scipy.optimize as spo
        res    = spo.minimize ( fun , x0 = x0 , bounds = bounds )
        if not res.success :
            logger.error ( "Can't minimize the function: %s" % res.message )
        return res.x[0]
        
    # =========================================================================
    ## get a maximum for 1D-function
    #  @code
    #  model = ...
    #  x = model.maximum() 
    #  @endcode 
    def sp_maximum_1D ( fun , xmin , xmax , x0 = None , *args ) :
        """ Get a maximum for 1D-function
        >>> model = ...
        >>> x = model.maximum () 
        >>>
        """
        funmin = lambda x , *a : -1.0 * ( float ( fun ( x , *a ) ) )
        return sp_minimum_1D ( funmin , xmin ,  xmax , x0 , *args )
    
    # =========================================================================
    ## get a minimum for 2D-function
    #  @code
    #  model2 = ...
    #  x , y = model2.minimum () 
    #  @endcode 
    def sp_minimum_2D ( fun  ,
                        xmin , xmax ,
                        ymin , ymax , x0 = () , *args ) :
        """ Get a maximum for 2D-function
        >>> model2 = ...
        >>> x , y = model2.maximum() 
        >>>
        """
        if not x0 :  x0 = 0.5 * ( xmin + xmax ) , 0.5 * ( ymin + ymax ) 

        x0     = numpy.array ( *x0  )
        
        bounds = [ ( xmin , xmax ) , ( ymin , ymax ) ]
        
        import scipy.optimize as spo
        res    = spo.minimize ( fun , x0 = x0 , bounds = bounds )
        if not res.success :
            logger.error ( "Can't minimize the function: %s" % res.message )
        return res.x[0] , res.x[1] 
            
    # =========================================================================
    ## get a maximum for 2D-function
    #  @code
    #  model2 = ...
    #  x , y = model2.maximum() 
    #  @endcode 
    def sp_maximum_2D ( fun ,
                        xmin , xmax ,
                        ymin , ymax , x0 = () , *args ) :
        """ Get a maximum for 2D-function
        >>> model2 = ...
        >>> x , y = model2.maximum () 
        >>>
        """
        funmin = lambda x , y , *a : -1.0 * ( float ( fun ( x , y , *a ) ) )
        return sp_minimum_2D  ( funmin ,
                                xmin , xmax ,
                                ymin , ymax , x0 , *args )
    
    # =========================================================================
    ## get a minimum for 3D-function
    #  @code
    #  model3 = ...
    #  x , y , z = model2.minimum () 
    #  @endcode 
    def sp_minimum_3D ( fun  ,
                        xmin , xmax ,
                        ymin , ymax ,
                        zmin , zmax , x0 = () , *args ) :
        """ Get a minimum for 3D-function
        >>> model3 = ...
        >>> x , y , z = model3.minimum() 
        >>>
        """
        if not x0 :  x0 = 0.5 * ( xmin + xmax ) , 0.5 * ( ymin + ymax ) , 0.5 * ( zmin + zmax ) 

        x0     = numpy.array ( *x0  )
        
        bounds = [ ( xmin , xmax ) , ( ymin , ymax ) , ( zmin , zmax ) ]
        
        import scipy.optimize as spo
        res    = spo.minimize ( fun , x0 = x0 , bounds = bounds )
        if not res.success :
            logger.error ( "Can't minimize the function: %s" % res.message )
        return res.x[0] , res.x[1] , res.x[2] 
        
    # =========================================================================
    ## get a maximum for 3D-function
    #  @code
    #  model3 = ...
    #  x , y , z = model3.maximum() 
    #  @endcode 
    def sp_maximum_3D ( fun ,
                        xmin , xmax ,
                        ymin , ymax ,
                        zmin , zmax , x0 = () , *args ) :
        """ Get a maximum for 3D-function
        >>> model3 = ...
        >>> x, y , z  = model3.maximum () 
        >>>
        """
        funmin = lambda x , y , z , *a : -1.0 * ( float ( fun ( x , y , z , *a ) ) )
        return sp_minimum_3D  ( funmin ,
                                xmin ,  xmax ,
                                ymin ,  ymax ,
                                zmin ,  zmax , x0 , *args )

    __all__ = __all__ +  (
        'sp_minimum_1D' ,
        'sp_minimum_2D' ,
        'sp_minimum_3D' ,
        'sp_maximum_1D' ,
        'sp_maximum_2D' ,
        'sp_maximum_3D' )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
