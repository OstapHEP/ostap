#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/gof1dw.py
#  Set of utulities for classical Two-Samples & Godness-of-Fit tests 
#  for 1D-weighted case 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2024-09-16
# =============================================================================
""" Set of utulities for classical Two-Samples & Godness-of-Fit tests for 1D-weighted case 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2024-09-16"
__all__     = (
    # =========================================================================
    'KolmogorovSmirnov'  , ## Kolmogorov-Smitnov GoF/Two-Samples estimator 
    'Kuiper'             , ## Kuiper             GoF/Two-Samples estimator 
    'AndersonDarling'    , ## Anderson-Darling   GoF/Two-Samples estimator 
    'CramerVonMises'     , ## Cramer-von Mises   GoF/Two-Samples estimator 
    'BerkJones'          , ## Berk-Jones         GoF/Two-Samples estimator     
    'ZA'                 , ## ZA                 GoF/Two-Samples estimator
    'ZC'                 , ## ZC                 GoF/Two-Samples estimator
    'ZK'                 , ## ZK                 GoF/Two-Samples estimator
    'Wasserstein1'       , ## Wasserstein-1      GoF/Two-Samples estimator 
    'Wasserstein2'       , ## Wasserstein-2      GoF/Two-Samples estimator
    'Rosenblatt'         , ## Rosenblatt         GoF/Two-Samples estimator
    'Watson'             , ## Watson-U2          GoF/Two-Samples estimator
    ## Rank-based, *NON* weighted
    'MannWhitney'        , ## Mann-Whitney       GoF/Two-Samples estimator
    'Savage'             , ## Savage             GoF/Two-Samples estimator
    # =========================================================================
) # ===========================================================================
# =============================================================================
from   ostap.utils.core        import typename, classprop 
from   ostap.stats.gof_np      import GoFnp
import ostap.stats.twosamples2 as     TS2 
import ostap.stats.utils       as     U 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.gof1dw' )
else                       : logger = getLogger ( __name__ )
# =============================================================================
logger.debug ( 'Two-sample & GoF 1D-weighted tests' )
# =============================================================================\
## @class ClassicalGoF 
#  Implementation of helper base class that converts
#  the classical (weighted-1D) Two-Sample tests into
#  Ostap GoF/Two-Sample infrastructure  
#  - p-value is obtained via permutation/bootstrapping technique  
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class ClassicalGoF (GoFnp) :
    """ Helper base class that converts
    the classical (weighted-1D) Two-Sample tests into
    Ostap GoF/Two-Sample infrastructure
    
    - p-value is obtained via permutation/bootstrapping technique    
    """
    def __init__ ( self                  ,
                   estimator             ,
                   nToys     = 1000 , **kwargs ) :

        if not callable ( estimator ) :
            raise TypeError ( "Invalid `estimator' type %s" % typename ( estimator ) )
        
        self.__estimator = estimator        
        kwargs [ 'normalize' ] = False
        
        GoFnp.__init__ ( self , nToys  = nToys  , **kwargs )
    
    @property
    def estimator ( self ) :
        """`estimator` : actual classic t-value estimator
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
        """`weights_supported`: Are weights supported by this estimator?
        """
        return True 

    # =========================================================================
    ## Good for two-samples comparison?
    #  Can this estimator be used for Two Samples comparison ?
    @property 
    def two_samples ( self ) :
        """`two_samples`: Can this estimator be used for Two Samples comparison ?
        """
        return True 
       
    # =========================================================================
    ## Calculate T-value for Goodness-of-Fit
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
                       
# =============================================================================
## @class KolmogorovSmirnov
#  Two (weighted) sample test using Kolmogorov-Smirnov statistics
class KolmogorovSmirnov(ClassicalGoF) :
    """ T(weighted 1D) Two-Sample test using Kolmogorov-Smirnov' statistics """
    def __init__ ( self , **kwargs ) :
        """ (weighted 1D) Two-Sample test using Kolmogorov-Smirnov statistics
        """
        super().__init__ ( estimator = TS2.kolmogorov_smirnov , **kwargs )

    # =========================================================================
    ## Symbol/tag associated with this estimator 
    @classprop
    def symbol ( self ) :
        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_KS
                
    # =========================================================================
    ## Actual method implemented by the estimator 
    @property 
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_KS 

# =============================================================================
## @class Kuiper
#  Two (weighted) sample test using Kuiperstatistics
class Kuiper(ClassicalGoF) :
    """ Two (weighted) sample test using Kuiper' statistics """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using Kuiper' statistics
        """
        super().__init__ ( estimator = TS2.kuiper , **kwargs )

    # =========================================================================
    ## Symbol/tag associated with this estimator 
    @classprop
    def symbol ( self ) :
        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_K
                
    # =========================================================================
    ## Actual method implemented by the estimator 
    @property 
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_K 

# =============================================================================
## @class CramerVonMises
#  Two (weighted) sample test using Kuiper' statistics
class CramerVonMises(ClassicalGoF) :
    """ Two (weighted) sample test using Cramer-von Mises statistics """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using Cramer-von Mises' statistics
        """
        super().__init__ ( estimator = TS2.cramer_von_mises , **kwargs )

    # =========================================================================
    ## Symbol/tag associated with this estimator 
    @classprop
    def symbol ( self ) :
        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_CVM
                
    # =========================================================================
    ## Actual method implemented by the estimator 
    @property 
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_CVM 
        
# =============================================================================
## @class AndersonDarling
#  Two (weighted) sample test using Anderson-Darling' statistics
class AndersonDarling(ClassicalGoF) :
    """ Two (weighted) sample test using Anderson-Darling statistics """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using Anderson-Darling' statistics
        """
        super().__init__ ( estimator = TS2.anderson_darling, **kwargs )
        
    # =========================================================================
    ## Symbol/tag associated with this estimator 
    @classprop
    def symbol ( self ) :
        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_AD 
                
    # =========================================================================
    ## Actual method implemented by the estimator 
    @property 
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_AD  

# =============================================================================
## @class BerkJones
#  Two (weighted) sample test using Berk-Jones' statistics
class BerkJones(ClassicalGoF) :
    """ Two (weighted) sample test using Berk-Jones statistics """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using Berk-Jones statistics
        """
        super().__init__ ( estimator = TS2.berk_jones, **kwargs )

    # =========================================================================
    ## Symbol/tag associated with this estimator 
    @classprop
    def symbol ( self ) :
        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_BJ
                
    # =========================================================================
    ## Actual method implemented by the estimator 
    @property 
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_BJ  
        
# =============================================================================
## @class ZA
#  Two (weighted) sample test using Zhang' ZA statistics
class ZA(ClassicalGoF) :
    """ Two (weighted) sample test using Zhang' ZA statistics """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using Zhang's ZA statistics
        """
        super().__init__ ( estimator = TS2.ZA      , **kwargs )
        
    # =========================================================================
    ## Symbol/tag associated with this estimator 
    @classprop
    def symbol ( self ) :
        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_ZA
                
    # =========================================================================
    ## Actual method implemented by the estimator 
    @property 
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_ZA  

# =============================================================================
## @class ZC
#  Two (weighted) sample test using Zhang' ZC statistics
class ZC(ClassicalGoF) :
    """ Two (weighted) sample test using Zhang' ZC statistics """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using Zhang's ZC statistics
        """
        super ().__init__ ( estimator = TS2.ZC      , **kwargs )

    # =========================================================================
    ## Symbol/tag associated with this estimator 
    @classprop
    def symbol ( self ) :
        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_ZC
                
    # =========================================================================
    ## Actual method implemented by the estimator 
    @property 
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_ZC

# =============================================================================
## @class ZK
#  Two (weighted) sample test using Zhang' ZK statistics
class ZK(ClassicalGoF) :
    """ Two (weighted) sample test using Zhang' ZK statistics """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using Zhang's ZK statistics
        """
        super().__init__ ( estimator = TS2.ZK , **kwargs )

    # =========================================================================
    ## Symbol/tag associated with this estimator 
    @classprop
    def symbol ( self ) :
        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_ZK
                
    # =========================================================================
    ## Actual method implemented by the estimator 
    @property 
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_ZK

# =============================================================================
## @class Wasserstein1 
#  Two (weighted) sample test using the fitrst Wasserstein distance 
class Wasserstein1(ClassicalGoF) :
    """ Two (weighted) sample test using the First Wasserstein distance """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using the First Wasserstein distance
        """
        super () .__init__ ( estimator = TS2.wasserstein1 , **kwargs )
        
    # =========================================================================
    ## Symbol/tag associated with this estimator 
    @classprop
    def symbol ( self ) :
        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_W1
                
    # =========================================================================
    ## Actual method implemented by the estimator 
    @property 
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_W1

# =============================================================================
## @class Wasserstein2
#  Two (weighted) sample test using the second Wasserstein distance 
class Wasserstein2(ClassicalGoF) :
    """ Two (weighted) sample test using the second Wasserstein distance """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using the second Wasserstein distance
        """
        super () .__init__ ( estimator = TS2.wasserstein2 , **kwargs )

    # =========================================================================
    ## Symbol/tag associated with this estimator 
    @classprop
    def symbol ( self ) :
        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_W2
                
    # =========================================================================
    ## Actual method implemented by the estimator 
    @property 
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_W2

# =============================================================================
## @class Rosenblatt 
#  Two (weighted) sample test using the Rosenblatt 
class Rosenblatt(ClassicalGoF) :
    """ Two (weighted) sample test using Rosenblatt distance """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using the Rosenblatt distance
        """
        super () .__init__ ( estimator = TS2.rosenblatt , **kwargs )

    # =========================================================================
    ## Symbol/tag associated with this estimator 
    @classprop
    def symbol ( self ) :
        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_R
                
    # =========================================================================
    ## Actual method implemented by the estimator 
    @property 
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_R

# =============================================================================
## @class Watson
#  Two (weighted) sample test using the Watson U2
class Watson(ClassicalGoF) :
    """ Two (weighted) sample test using Watson U2 distance """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using the Watson U2 distance
        """
        super () .__init__ (  estimator = TS2.watson_u2 , **kwargs )
        
    # =========================================================================
    ## Symbol/tag associated with this estimator 
    @classprop
    def symbol ( self ) :
        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_U2
                
    # =========================================================================
    ## Actual method implemented by the estimator 
    @property 
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_U2

# =============================================================================
## @class MannWhitney 
#  Two (*NON*-weighted) sample test using the Mann-Whitney
class MannWhitney(ClassicalGoF) :
    """ Two (*NON*-weighted) sample test using Mann-Whitney distance """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using the Mann-Whitney distance
        """
        super () .__init__ (  estimator = TS2.mann_whitney , **kwargs )
        
    # =========================================================================
    ## Symbol/tag associated with this estimator 
    @classprop
    def symbol ( self ) :
        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_MW
                
    # =========================================================================
    ## Actual method implemented by the estimator 
    @property 
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_MW
    
    # =========================================================================
    ## Are weights supported by this GoF estimator?
    @property 
    def weights_supported ( self ) :
        """`weights_supported`: Are weights supported by this estimator?
        """
        return False 

# =============================================================================
## @class Savage
#  Two (*NON*-weighted) sample test using the Savage score 
class Savage(ClassicalGoF) :
    """ Two (*NON*-weighted) sample test using Savage distance """
    def __init__ ( self , **kwargs ) :
        """ Two (weighted) sample test using the Mann-Whitney distance
        """
        super () .__init__ (  estimator = TS2.savage , **kwargs )
        
    # =========================================================================
    ## Symbol/tag associated with this estimator 
    @classprop
    def symbol ( self ) :
        """`symbol`: Symbol/tag associated with this estimator
        """
        return U.symbol_S
                
    # =========================================================================
    ## Actual method implemented by the estimator 
    @property 
    def method ( self ) :
        """`method`: Actual method implemented by the estimator
        """
        return U.method_S
    
    # =========================================================================
    ## Are weights supported by this GoF estimator?
    @property 
    def weights_supported ( self ) :
        """`weights_supported`: Are weights supported by this estimator?
        """
        return False 

    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
