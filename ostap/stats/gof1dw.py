#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/gof1dw.py
#  Set of utulities for two-sample tests and Godness-of-fit studies for (weighted) 1D-fits
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2024-09-16
# =============================================================================
""" Set of utulities for two-sample test and Godness-of-fit studies for (weighted) 1D-fits
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2024-09-16"
__all__     = (
    'KolmogorovSmirnov'  , ## Kolmogorov-Sminov GoF estimator 
    'Kuiper'             , ## Kuiper            GoF estimator 
    'AndersonDarling'    , ## Anderson-Darling  GoF estimator 
    'CramerVonMises'     , ## Cramer-von Mises  GoF estimator 
    'BerkJones'          , ## Berk-Jones        GoF estimator     
    'ZK'                 , ## ZK               GoF estimator
    'ZA'                 , ## ZA               GoF estimator
    'ZC'                 , ## ZC               GoF estimator
)
# =============================================================================
from   ostap.stats.gof_np      import GoFnp
from   ostap.stats.gofnd       import GoF
import ostap.stats.twosamples2 as     TS2 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.twosamples2' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Two-sample & GoF 1D-weighted tests' )
# =============================================================================
## @class BootstrapGoF 
#  Implementation of 1D GoFnp base class for getting p-value via boostrapping 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class BootstrapGoF ( GoFnp ) :
    """ 1D GoFnp base class to get p-value via BootStrap
    """
    def __init__ ( self            ,
                   method          ,
                   estimator       ,
                   nToys     = 400 , 
                   **kwargs        ) :
        
        self.__estimator = estimator        
        kwargs [ 'normalize' ] = False
        
        GoFnp.__init__ ( self ,
                         method = method ,
                         nToys  = nToys  , **kwargs )
    
    @property
    def estimator ( self ) :
        """`estimator` : actual t-value estimator
        """
        return self.__estimator
    
    # ==================================================================================
    @property
    def config ( self ) :
        """`config` : get all configuration parameters"""
        conf = {} 
        conf.update ( super().config  )
        conf [ 'estimator' ] = self.estimator
        return conf
            
    # =========================================================================
    ## Are weights supported by this GoF estimator?
    @property 
    def weights_supported ( self ) :
        """`weghts_supported`: Are weights supported by this estimator?
        """
        return True 

    # =========================================================================
    ## Are weights supported by this GoF estimator?
    @property 
    def negative_weights_supported ( self ) :
        """`negative_weghts_supported`: Are negative weights supported by this estimator?
        """
        return True
    
    # =========================================================================
    ## Good for two-samples comparison?
    #  Can this estimator be used for comparison of two samples?
    @property 
    def two_samples ( self ) :
        """`two_samples`: Can this estimator be used for comparison of two samples?
        """
        return True 
       
    # =========================================================================
    ## Calculate T-value for Goodness-of-Git
    #  @code
    #  data1   = ...
    #  data2   = ...
    #  weight1 = ...
    #  weight2 = ...
    #  tvalue  = gof.tvalue ( data1 , data2 , weight1 , weight2 )
    #  @endcode
    def tvalue ( self      ,
                 data1     ,
                 data2     ,
                 weight1   = None  ,
                 weight2   = None  ,
                 normalize = False ) :
        """ Calculate T-value for Goodness-of-Fit
        >>> data1   = ...
        >>> data2   = ...
        >>> weight1 = ...
        >>> weight2 = ...
        >>>> tvalue  = gof.tvalue ( data1 , data2 , weight1 , weight2 )
        """
        return self.estimator ( data1   = data1   ,
                                data2   = data2   ,
                                weight1 = weight1 ,
                                weight2 = weight2 )
    
    # =========================================================================
    ## Calculate the t & p-values
    #  @code
    #  gof = ...
    #  data1 , data2 = ...
    #  t , p = gof.pvalue ( data1 , data2 , normalize = False ) 
    #  @endcode 
    def pvalue ( self           , 
                 data1          ,
                 data2          , * ,
                 tvalue  = None , 
                 weight1 = None ,
                 weight2 = None ) : 
                
        """ Calculate the t & p-values
        >>> gof  = ...
        >>> data1 , data2 = ...
        >>> t   , p = gof.pvalue ( ds1 , ds2 , normalize = True ) 
        """

        ## transform ?
        uds1 , uds2 = self.unpack ( data1 , data2 ) 
        
        ## normalize ? 
        if self.normalize :
            uds1 , uds2 = self.normalize_pooled ( uds1 , uds2 ) 

        ## calculate t-value if not specified 
        t_value    = tvalue if not tvalue is None else self.tvalue ( uds1      ,
                                                                     uds2      ,
                                                                     weight1   = weight1 ,
                                                                     weight2   = weight2 ,
                                                                     normalize = False   )        
        ## use bootstrapping to get the p-value
        ## from ostap.stats.pvalue import BOOTSTRAPPER as RESAMPLER 
        from ostap.stats.pvalue import PERMUTATOR as RESAMPLER 
        resampler = RESAMPLER ( self                   ,
                                t_value                , 
                                uds1                   ,
                                uds2                   ,
                                weight1 = weight1      ,
                                weight2 = weight2      ) ;
        
        if self.parallel and resampler.run : counter , _ = resampler.run ( self.nToys , progress = self.progress )            
        else                               : counter , _ = resampler     ( self.nToys , progress = self.progress )
        
        ## get the efficiency/p-value from the counter
        p_value      = counter.eff

        self.ecdf    = resampler.ecdf
        self.counter = counter
        
        self.t_value = t_value 
        self.p_value = p_value 
        
        return self.t_value , self.p_value
    
# =============================================================================
## @class KolmogorovSmirnov
#  Two (weighted) sample test using Kolmogorov-Smirnov statistics
class KolmogorovSmirnov(BootstrapGoF) :
    """ Two (weighted) sample test using Kolmogorov-Smirnov' statistics """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using Kolmogorov-Smirnov statistics
        """
        BootstrapGoF.__init__ ( self ,
                                method    = 'Kolmogorov-Smirnov'   ,
                                estimator = TS2.kolmogorov_smirnov , **kwargs )

# =============================================================================
## @class Kuiper
#  Two (weighted) sample test using Kuiperstatistics
class Kuiper(BootstrapGoF) :
    """ Two (weighted) sample test using Kuiper' statistics """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using Kuiper' statistics
        """
        BootstrapGoF.__init__ ( self ,
                                method    = 'Kuiper'   ,
                                estimator = TS2.kuiper , **kwargs )

# =============================================================================
## @class CramerVonMises
#  Two (weighted) sample test using Kuiper' statistics
class CramerVonMises(BootstrapGoF) :
    """ Two (weighted) sample test using Cramer-von Mises statistics """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using Cramer-von Mises' statistics
        """
        BootstrapGoF.__init__ ( self ,
                                method    = 'Cramer-von Mises'   ,
                                estimator = TS2.cramer_von_mises , **kwargs )

# =============================================================================
## @class AndersonDarling
#  Two (weighted) sample test using Anderson-Darling' statistics
class AndersonDarling(BootstrapGoF) :
    """ Two (weighted) sample test using Anderson-Darling statistics """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using Anderson-Darling' statistics
        """
        BootstrapGoF.__init__ ( self ,
                                method    = 'Anderson-Darling'   ,
                                estimator = TS2.anderson_darling, **kwargs )
        
# =============================================================================
## @class BerkJones
#  Two (weighted) sample test using Berk-Jones' statistics
class BerkJones(BootstrapGoF) :
    """ Two (weighted) sample test using Berk-Jones statistics """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using Berk-Jones statistics
        """
        BootstrapGoF.__init__ ( self ,
                                method    = 'Berk-Jones'   ,
                                estimator = TS2.berk_jones, **kwargs )

# =============================================================================
## @class ZA
#  Two (weighted) sample test using Zhang' ZA statistics
class ZA(BootstrapGoF) :
    """ Two (weighted) sample test using Zhang' ZA statistics """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using Zhang's ZA statistics
        """
        BootstrapGoF.__init__ ( self ,
                                method    = "Zhang' ZA" ,
                                estimator = TS2.ZA      , **kwargs )

# =============================================================================
## @class ZC
#  Two (weighted) sample test using Zhang' ZC statistics
class ZC(BootstrapGoF) :
    """ Two (weighted) sample test using Zhang' ZC statistics """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using Zhang's ZC statistics
        """
        BootstrapGoF.__init__ ( self ,
                                method    = "Zhang' ZC" ,
                                estimator = TS2.ZC      , **kwargs )

# =============================================================================
## @class ZK
#  Two (weighted) sample test using Zhang' ZK statistics
class ZK(BootstrapGoF) :
    """ Two (weighted) sample test using Zhang' ZK statistics """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using Zhang's ZK statistics
        """
        BootstrapGoF.__init__ ( self ,
                                method    = "Zhang' ZK" ,
                                estimator = TS2.ZK      , **kwargs )
        
# =============================================================================

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
