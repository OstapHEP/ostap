#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostath/integrator.py
#  Decoration modlel for C++ class Ostap::Math::Integrator
#  @see Ostap::Math::Integrator
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2022-12-18
# =============================================================================
""" Decoration modlel for C++ class Ostap::Math::Integrator
- see `Ostap.Math.Integrator`
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2022-12-18"
__all__     = (
    "Integrator" , ## C++ numerical integrator  
    )
# =============================================================================
from ostap.math.base      import Ostap, doubles 
from ostap.core.meta_info import root_info 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.integrator' )
else                       : logger = getLogger ( __name__                )

# =============================================================================
## actual C++ interator 
Integrator = Ostap.Math.Integrator


# =============================================================================
if root_info < (6,18) :
    
    _PyC1 = Ostap.Functions.PyCallable
    _PyC2 = Ostap.Functions.PyCallable2
    _PyC3 = Ostap.Functions.PyCallable3

    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_integrate 
    #  @see Ostap::Functions::.PyCallable        
    def _pyi_integrate ( integrator       ,
                         fun              ,
                         xmin             ,
                         xmax             ,
                         tag        = 0   ,
                         rescale    = 0   ,
                         aprecision = 0.0 ,
                         rprecision = 0.0 ,
                         rule       = 0   ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_integrate
        - see Ostap.Functions.PyCallable        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun = fun if isinstance ( fun , _PyC1 ) else _PyC1 ( fun , True ) 
        return Ostap.Math.py_integrate (
            integrator  ,
            pyfun       ,
            xmin , xmax , 
            tag  , rescale ,
            aprecision  ,
            rprecision  ,
            rule        )
    
    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_integrate_infinity 
    #  @see Ostap::Functions::.PyCallable        
    def _pyi_integrate_infinity  ( integrator       ,
                                   fun              ,
                                   tag        = 0   ,
                                   aprecision = 0.0 ,
                                   rprecision = 0.0 ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_integrate_infinity
        - see Ostap.Functions.PyCallable        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun = fun if isinstance ( fun , _PyC1 ) else _PyC1 ( fun , True ) 
        return Ostap.Math.py_integrate_infinity (
            integrator  ,
            pyfun       ,
            tag         , 
            aprecision  ,
            rprecision  )

    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_integrate_to_infinity 
    #  @see Ostap::Functions::.PyCallable        
    def _pyi_integrate_to_infinity  ( integrator       ,
                                      fun              ,
                                      xmin             , 
                                      tag        = 0   ,
                                      aprecision = 0.0 ,
                                      rprecision = 0.0 ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_integrate_to_infinity
        - see Ostap.Functions.PyCallable        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun = fun if isinstance ( fun , _PyC1 ) else _PyC1 ( fun , True ) 
        return Ostap.Math.py_integrate_to_infinity (
            integrator  ,
            pyfun       ,
            xmin        ,
            tag         ,
            aprecision  ,
            rprecision  )

    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_integrate_from_infinity 
    #  @see Ostap::Functions::.PyCallable        
    def _pyi_integrate_from_infinity  ( integrator       ,
                                        fun              ,
                                        xmax             , 
                                        tag        = 0   ,
                                        aprecision = 0.0 ,
                                        rprecision = 0.0 ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_integrate_from_infinity
        - see Ostap.Functions.PyCallable        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun = fun if isinstance ( fun , _PyC1 ) else _PyC1 ( fun , True ) 
        return Ostap.Math.py_integrate_to_infinity (
            integrator  ,
            pyfun       ,
            xmax        ,
            tag         ,
            aprecision  ,
            rprecision  )

    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_cauchy_pv 
    #  @see Ostap::Functions::.PyCallable        
    def _pyi_cauchy_pv ( integrator       ,
                         fun              ,
                         c                , 
                         xmin             , 
                         xmax             , 
                         tag        = 0   ,
                         rescale    = 0   , 
                         aprecision = 0.0 ,
                         rprecision = 0.0 ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_cauchy_pvy
        - see Ostap.Functions.PyCallable        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun = fun if isinstance ( fun , _PyC1 ) else _PyC1 ( fun , True ) 
        return Ostap.Math.py_cauchy_pv (
            integrator  ,
            pyfun       ,
            c           ,
            xmin        ,
            xmax        ,
            tag         ,
            rescale     , 
            aprecision  ,
            rprecision  )

    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_cauchy_pv_to_infinity 
    #  @see Ostap::Functions::.PyCallable        
    def _pyi_cauchy_pv_to_infinity ( integrator       ,
                                     fun              ,
                                     c                , 
                                     xmin             , 
                                     tag        = 0   ,
                                     rescale    = 0   , 
                                     aprecision = 0.0 ,
                                     rprecision = 0.0 ,
                                     width      = 0.0 ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_cauchy_pv_to_infinity
        - see Ostap.Functions.PyCallable        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun = fun if isinstance ( fun , _PyC1 ) else _PyC1 ( fun , True ) 
        return Ostap.Math.py_cauchy_pv_to_infinity (
            integrator  ,
            pyfun       ,
            c           ,
            xmin        ,
            tag         ,
            rescale     , 
            aprecision  ,
            rprecision  ,
            width       )

    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_cauchy_pv_from_infinity 
    #  @see Ostap::Functions::.PyCallable        
    def _pyi_cauchy_pv_from_infinity ( integrator       ,
                                       fun              ,
                                       c                , 
                                       xmax             , 
                                       tag        = 0   ,
                                       rescale    = 0   , 
                                       aprecision = 0.0 ,
                                       rprecision = 0.0 ,
                                       width      = 0.0 ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_cauchy_pv_from_infinity
        - see Ostap.Functions.PyCallable        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun = fun if isinstance ( fun , _PyC1 ) else _PyC1 ( fun , True ) 
        return Ostap.Math.py_cauchy_pv_from_infinity (
            integrator  ,
            pyfun       ,
            c           ,
            xmax        ,
            tag         ,
            rescale     , 
            aprecision  ,
            rprecision  ,
            width       )    

    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_cauchy_pv_infinity 
    #  @see Ostap::Functions::.PyCallable        
    def _pyi_cauchy_pv_infinity ( integrator       ,
                                  fun              ,
                                  c                , 
                                  tag        = 0   ,
                                  rescale    = 0   , 
                                  aprecision = 0.0 ,
                                  rprecision = 0.0 ,
                                  width      = 0.0 ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_cauchy_pv_infinity
        - see Ostap.Functions.PyCallable        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun = fun if isinstance ( fun , _PyC1 ) else _PyC1 ( fun , True ) 
        return Ostap.Math.py_cauchy_pv_infinity (
            integrator  ,
            pyfun       ,
            c           ,
            tag         ,
            rescale     , 
            aprecision  ,
            rprecision  ,
            width       )

    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_kramers_kronig  
    #  @see Ostap::Functions::.PyCallable        
    def _pyi_kramers_kronig ( integrator       ,
                              fun              ,
                              s                ,
                              xmin             ,
                              n          = 0   , 
                              tag        = 0   ,
                              rescale    = 0   , 
                              aprecision = 0.0 ,
                              rprecision = 0.0 ,
                              width      = 0.0 ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_kramers_kronig
        - see Ostap.Functions.PyCallable        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun = fun if isinstance ( fun , _PyC1 ) else _PyC1 ( fun , True ) 
        return Ostap.Math.py_kramers_kronig (
            integrator  ,
            pyfun       ,
            s           ,
            xmin        ,
            n           , 
            tag         ,
            rescale     , 
            aprecision  ,
            rprecision  ,
            width       )
        
    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_integrate_singular
    #  @see Ostap::Functions::.PyCallable        
    def _pyi_integrate_singular  ( integrator       ,
                                   fun              ,
                                   xmin             , 
                                   xmax             ,
                                   points           , 
                                   tag        = 0   ,
                                   aprecision = 0.0 ,
                                        rprecision = 0.0 ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_integrate_singular 
        - see Ostap.Functions.PyCallable        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun   = fun if isinstance ( fun , _PyC1 ) else _PyC1 ( fun , True )
        spoints = doubles ( points ) 
        return Ostap.Math.py_integrate_singular (
            integrator  ,
            pyfun       ,
            xmin        ,
            xmax        ,
            spoints     , 
            tag         ,
            aprecision  ,
            rprecision  )

    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_integrate_cquad
    #  @see Ostap::Functions::.PyCallable        
    def _pyi_integrate_cquad  ( integrator       ,
                                fun              ,
                                xmin             , 
                                xmax             ,
                                tag        = 0   ,
                                rescale    = 0   , 
                                aprecision = 0.0 ,
                                rprecision = 0.0 ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_integrate_cquad
        - see Ostap.Functions.PyCallable        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun   = fun if isinstance ( fun , _PyC1 ) else _PyC1 ( fun , True )
        return Ostap.Math.py_integrate_cquad (
            integrator  ,
            pyfun       ,
            xmin        ,
            xmax        ,
            tag         ,
            rescale     , 
            aprecision  ,
            rprecision  )

    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_integrate_romberg
    #  @see Ostap::Functions::.PyCallable        
    def _pyi_integrate_romberg ( integrator       ,
                                 fun              ,
                                 xmin             , 
                                 xmax             ,
                                 tag        = 0   ,
                                 rescale    = 0   , 
                                 aprecision = 0.0 ,
                                 rprecision = 0.0 ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_integrate_romberg
        - see Ostap.Functions.PyCallable        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun   = fun if isinstance ( fun , _PyC1 ) else _PyC1 ( fun , True )
        return Ostap.Math.py_integrate_romberg (
            integrator  ,
            pyfun       ,
            xmin        ,
            xmax        ,
            tag         ,
            rescale     , 
            aprecision  ,
            rprecision  )
    
    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_integrate2 
    #  @see Ostap::Functions::.PyCallable2        
    def _pyi_integrate2 ( integrator       ,
                          fun              ,
                          xmin             ,
                          xmax             ,
                          ymin             ,
                          ymax             ,
                          tag        = 0   ,
                          aprecision = 0.0 ,
                          rprecision = 0.0 ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_integrate2
        - see Ostap.Functions.PyCallable2        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun = fun if isinstance ( fun , _PyC2 ) else _PyC2 ( fun , True ) 
        return Ostap.Math.py_integrate2 (
            integrator  ,
            pyfun       ,
            xmin , xmax , 
            ymin , ymax , 
            tag  , 
            aprecision  ,
            rprecision  )

    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_integrate2X 
    #  @see Ostap::Functions::.PyCallable2        
    def _pyi_integrate2X ( integrator       ,
                           fun              ,
                           y                , 
                           xmin             ,
                           xmax             ,
                           tag        = 0   ,
                           rescale    = 0   , 
                           aprecision = 0.0 ,
                           rprecision = 0.0 ,
                           rule       = 0   ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_integrate2X
        - see Ostap.Functions.PyCallable2        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun = fun if isinstance ( fun , _PyC2 ) else _PyC2 ( fun , True ) 
        return Ostap.Math.py_integrate2X (
            integrator  ,
            pyfun       ,
            y           , 
            xmin , xmax , 
            tag         ,
            rescale     , 
            aprecision  ,
            rprecision  ,
            rule        )
    
    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_integrate2Y 
    #  @see Ostap::Functions::.PyCallable2        
    def _pyi_integrate2Y ( integrator       ,
                           fun              ,
                           x                , 
                           ymin             ,
                           ymax             ,
                           tag        = 0   ,
                           rescale    = 0   , 
                           aprecision = 0.0 ,
                           rprecision = 0.0 ,
                           rule       = 0   ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_integrate2Y
        - see Ostap.Functions.PyCallable2        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun = fun if isinstance ( fun , _PyC2 ) else _PyC2 ( fun , True ) 
        return Ostap.Math.py_integrate2Y (
            integrator  ,
            pyfun       ,
            x           , 
            ymin , ymax , 
            tag         ,
            rescale     , 
            aprecision  ,
            rprecision  ,
            rule        )
    
    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_integrate3 
    #  @see Ostap::Functions::.PyCallable3        
    def _pyi_integrate3 ( integrator       ,
                          fun              ,
                          xmin             ,
                          xmax             ,
                          ymin             ,
                          ymax             ,
                          zmin             ,
                          zmax             ,
                          tag        = 0   ,
                          aprecision = 0.0 ,
                          rprecision = 0.0 ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_integrate3
        - see Ostap.Functions.PyCallable3        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun = fun if isinstance ( fun , _PyC3 ) else _PyC3 ( fun , True ) 
        return Ostap.Math.py_integrate3 (
            integrator  ,
            pyfun       ,
            xmin , xmax , 
            ymin , ymax , 
            zmin , zmax , 
            tag         , 
            aprecision  ,
            rprecision  )

    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_integrate3XY
    #  @see Ostap::Functions::.PyCallable3        
    def _pyi_integrate3XY ( integrator       ,
                            fun              ,
                            z                ,
                            xmin             ,
                            xmax             ,
                            ymin             ,
                            ymax             ,
                            tag        = 0   ,
                            aprecision = 0.0 ,
                            rprecision = 0.0 ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_integrate3XY
        - see Ostap.Functions.PyCallable3        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun = fun if isinstance ( fun , _PyC3 ) else _PyC3 ( fun , True ) 
        return Ostap.Math.py_integrate3XY (
            integrator  ,
            pyfun       ,
            z           , 
            xmin , xmax , 
            ymin , ymax , 
            tag         , 
            aprecision  ,
            rprecision  )

    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_integrate3XZ
    #  @see Ostap::Functions::.PyCallable3        
    def _pyi_integrate3XZ ( integrator       ,
                            fun              ,
                            y                ,
                            xmin             ,
                            xmax             ,
                            zmin             ,
                            zmax             ,
                            tag        = 0   ,
                            aprecision = 0.0 ,
                            rprecision = 0.0 ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_integrate3XZ
        - see Ostap.Functions.PyCallable3        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun = fun if isinstance ( fun , _PyC3 ) else _PyC3 ( fun , True ) 
        return Ostap.Math.py_integrate3XZ (
            integrator  ,
            pyfun       ,
            y           , 
            xmin , xmax , 
            zmin , zmax , 
            tag         , 
            aprecision  ,
            rprecision  )

    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_integrate3YZ
    #  @see Ostap::Functions::.PyCallable3        
    def _pyi_integrate3YZ ( integrator       ,
                            fun              ,
                            z                ,
                            ymin             ,
                            ymax             ,
                            zmin             ,
                            zmax             ,
                            tag        = 0   ,
                            aprecision = 0.0 ,
                            rprecision = 0.0 ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_integrate3YZ
        - see Ostap.Functions.PyCallable3        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun = fun if isinstance ( fun , _PyC3 ) else _PyC3 ( fun , True ) 
        return Ostap.Math.py_integrate3YZ (
            integrator  ,
            pyfun       ,
            x           , 
            ymin , ymax , 
            zmin , zmax , 
            tag         , 
            aprecision  ,
            rprecision  )

    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_integrate3X
    #  @see Ostap::Functions::.PyCallable3        
    def _pyi_integrate3X  ( integrator       ,
                            fun              ,
                            y                ,
                            z                ,
                            xmin             ,
                            xmax             ,
                            tag        = 0   ,
                            rescale    = 0   , 
                            aprecision = 0.0 ,
                            rprecision = 0.0 ,
                            rule       = 0   ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_integrate3X
        - see Ostap.Functions.PyCallable3        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun = fun if isinstance ( fun , _PyC3 ) else _PyC3 ( fun , True ) 
        return Ostap.Math.py_integrate3X  (
            integrator  ,
            pyfun       ,
            y    , z    , 
            xmin , xmax , 
            tag         ,
            rescale     , 
            aprecision  ,
            rprecision  ,
            rule        )

    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_integrate3Y
    #  @see Ostap::Functions::.PyCallable3        
    def _pyi_integrate3Y  ( integrator       ,
                            fun              ,
                            x                ,
                            z                ,
                            ymin             ,
                            ymax             ,
                            tag        = 0   ,
                            rescale    = 0   , 
                            aprecision = 0.0 ,
                            rprecision = 0.0 ,
                            rule       = 0   ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_integrate3Y
        - see Ostap.Functions.PyCallable3        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun = fun if isinstance ( fun , _PyC3 ) else _PyC3 ( fun , True ) 
        return Ostap.Math.py_integrate3Y  (
            integrator  ,
            pyfun       ,
            x    , z    , 
            ymin , ymax , 
            tag         ,
            rescale     , 
            aprecision  ,
            rprecision  ,
            rule        )

    # =========================================================================
    ## Integrate python callable using Ostap.Math.Integrator
    #  @see Ostap::Math::py_integrate3Z
    #  @see Ostap::Functions::.PyCallable3        
    def _pyi_integrate3Z  ( integrator       ,
                            fun              ,
                            x                ,
                            y                ,
                            zmin             ,
                            zmax             ,
                            tag        = 0   ,
                            rescale    = 0   , 
                            aprecision = 0.0 ,
                            rprecision = 0.0 ,
                            rule       = 0   ) :
        """Integrate python callable using Ostap.Math.Integrator
        - see Ostap.Math.py_integrate3Z
        - see Ostap.Functions.PyCallable3        
        """
        assert callable ( fun ) , 'Function must be callable!'
        pyfun = fun if isinstance ( fun , _PyC3 ) else _PyC3 ( fun , True ) 
        return Ostap.Math.py_integrate3Z  (
            integrator  ,
            pyfun       ,
            x    , y    , 
            zmin , zmax , 
            tag         ,
            rescale     , 
            aprecision  ,
            rprecision  ,
            rule        )

    
    Integrator.integrate               = _pyi_integrate
    Integrator.integrate_infinity      = _pyi_integrate_infinity 
    Integrator.integrate_to_infinity   = _pyi_integrate_to_infinity 
    Integrator.integrate_from_infinity = _pyi_integrate_from_infinity 
    Integrator.cauchy_pv               = _pyi_cauchy_pv 
    Integrator.cauchy_pv_infinity      = _pyi_cauchy_pv_infinity 
    Integrator.cauchy_pv_to_infinity   = _pyi_cauchy_pv_to_infinity 
    Integrator.cauchy_pv_from_infinity = _pyi_cauchy_pv_from_infinity 
    Integrator.kramers_kronig          = _pyi_kramers_kronig   
    Integrator.integrate_singular      = _pyi_integrate_singular
    
    Integrator.integrate_cquad         = _pyi_integrate_cquad
    Integrator.integrate_romberg       = _pyi_integrate_romberg
    
    Integrator.integrate2              = _pyi_integrate2
    Integrator.integrate2X             = _pyi_integrate2X
    Integrator.integrate2Y             = _pyi_integrate2Y


    Integrator.integrate3              = _pyi_integrate3
    Integrator.integrate3XY            = _pyi_integrate3XY
    Integrator.integrate3XZ            = _pyi_integrate3XZ
    Integrator.integrate3YZ            = _pyi_integrate3YZ
    Integrator.integrate3X             = _pyi_integrate3X
    Integrator.integrate3Y             = _pyi_integrate3Y 
    Integrator.integrate3Z             = _pyi_integrate3Z

    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                     The EMD   
# =============================================================================

