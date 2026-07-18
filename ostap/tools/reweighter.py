#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with utilities for reweigting using GBReweighetr 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2021-09-22
# =============================================================================
""" Module with utilities for reweighting using GBReweighter
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'Reweighter' , 
    ) 
# =============================================================================
from   ostap.math.math_base  import weight_trivial
import numpy, os  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.tools.reweighter' )
else                       : logger = getLogger( __name__ )
# =============================================================================
try : # =======================================================================
    # =========================================================================
    varname = 'NUMEXPR_MAX_THREADS'
    if not varname in os.environ :
        from ostap.utils.basic import numcpu 
        os.environ[ varname ]  = '%d'% numcpu ()
    # =========================================================================
    ## some patch 
    import sklearn.tree._classes as _cl
    if hasattr ( _cl , 'CRITERIA_REG' ) :
        if 'mse' in _cl.CRITERIA_REG and 'squared_error' not in _cl.CRITERIA_REG:
            _cl.CRITERIA_REG [ 'squared_error' ] = _cl.CRITERIA_REG [ 'mse']
        elif 'squared_error' in _cl.CRITERIA_REG and 'mse' not in _cl.CRITERIA_REG:
            _cl.CRITERIA_REG [ 'mse'] = _cl.CRITERIA_REG [ 'squared_error' ]

    ## 
    from hep_ml.reweight import GBReweighter as GBRW 
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    logger.error ( "Cannot import `GBReweighter` from `hep_ml.reweight`" , exc_info = True )
    GBRW = None 
    __all__ = () 
# =============================================================================
if GBRW : # ===================================================================
    # =========================================================================
    ## if not hasattr ( numpy , 'float' ) :
    ##    logger.warning ( 'No `numpy.float`... add it!')
    ##    numpy.float = numpy.float64
    # =========================================================================
    ## @class Reweighter
    #  Helper class for reweighting using <code>GBReweighter</code>
    #  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
    #  @date   2021-09-22 
    class Reweighter(object) :
        """ Helper class for reweighting using `GBReweighter`
        - see hep_ml.reweight.GBReweighter
        """
        def __init__ ( self                    , * , 
                       original                ,
                       target                  ,
                       original_weight = None  ,
                       target_weight   = None  , 
                       silent          = False , **params  ) :
            """ Helper class for reweighting using `GBReweighter`
            - see hep_ml.reweight.GBReweighter
            """
            cnf = { 'n_estimators'     : 60  , 
                    'learning_rate'    : 0.1 , 
                    'max_depth'        : 3   , 
                    'min_samples_leaf' : 100 }
            
            params .pop ( 'n_splits' , None )
            params .pop ( 'parallel' , None )
            params .pop ( 'progress' , None )
            cnf.update ( params ) 
            
            self.__params     = cnf 
            self.__silent     = True if silent else False 
            self.__reweighter = GBRW ( **self.params  ) 
            shape1 = original.shape
            shape2 = target.shape 
            
            assert len ( shape1 ) == len ( shape2 ) ,  "Inconsistent shapes: %s vs %s " % ( shape1 , shape2 )   
            assert len ( shape1 ) == 1 or shape1[-1] == shape2[-1] , \
                "Inconsistent shapes: %s vs %s " % ( shape1[-1] , shape2[-1] )   
                       
            assert ( original_weight is None ) or len ( original_weight ) == len ( original ) , \
                "Invalid length of ``original weights''"       
                 
            assert ( target_weight   is None ) or len ( target_weight   ) == len ( target   ) , \
                "Invalid length of ``target weights''"
        
            ## MC weights are used in training?
            self.__weight_used = not weight_trivial ( original_weight )
                
            if not silent and self.params : 
                title = '(GB)Reweighter configuration'
                table = self.table ( title = title  , prefix = '# ' )
                logger.info ( '%s:\n%s' %  ( title , table ) )

            if silent : 
                from ostap.logger.logger import logAttention as context
            else      : 
                from ostap.utils.basic   import NoContext    as context 
            
            with context () :
                self.__reweighter.fit ( original        ,
                                        target          ,
                                        original_weight = original_weight , 
                                        target_weight   = target_weight   )
                         

        @property
        def method ( self ) :
            """`method` : underlying method/engine"""
            return "GBReweighter"
    
        @property
        def params ( self ) :
            """`kwargs`: actual configuration usef for clreation of `GBReweighter`"""
            return self.__params
        
        @property
        def weight_used ( self ) :
            """`weight_used` : was orignal weight used for training?"""
            return self.__weight_used 
        
        @property
        def silent ( self ) :
            """`silent` : silent processing?"""
            return self.__silent
        
        @property 
        def config ( self ) :
            """`config` : Reweighter configuraton"""
            conf = {}
            conf.update ( self.__params )
            conf [ 'method'      ] = self.method
            conf [ 'silent'      ] = self.silent 
            conf [ 'weight_used' ] = self.weight_used
            return conf 
                    
        @property
        def reweighter ( self ) :
            """`reweighter' : get the underlying reweighter object"""
            return self.__reweighter
         
        # =========================================================================
        ## self-print get the configuration 
        def table (  self , title = '' , prefix = '# ') : 
            """ print configuration """
            from ostap.logger.utils import map2table_ex
            title = title if title else "%s configuration " % typename ( self ) 
            return map2table_ex ( self.config , 
                                header      = ( 'Parameter' , 'type' , 'value' ) ,
                                  ailgnment   = 'rcw'  , 
                                  prefix      = prefix ,
                                  title       = title  )
        
        def __str__  ( self ) : return self.table ( prefix = '' )
        def __repr__ ( self ) : return self.__str__ () 
        
        # =========================================================================
        ## Get/predict new weights for (new) original
        def weight ( self                   ,
                     original               ,
                     original_weight = None ) :
            """ Get/predict  new weights for (new) original
            """
            if not self.silent and self.weight_used and weight_trivial ( original_weight ) :
                logger.warning ( "Reweighter: `original-weight' was used for training but not provided for evaluation" ) 
                
            return self.reweighter.predict_weights (
                original        = original        ,
                original_weight = original_weight )

        weights     = weight 
        new_weight  = weight
        new_weights = weight
    
        # ==========================================================================
        ## Get/predict new weights for (new) original
        def __call__ ( self                   ,
                       original               ,
                       original_weight = None ) :
            """ Get/predict  new weights for (new) original
            """
            return self.weight ( original , original_weight = original_weight )

    ## extend a  bit the documentation string
    __doc__                     += '\n\nGBReweighter documentation:\n%s' % GBRW.__init__.__doc__ 
    Reweighter.__doc__          += '\n\nGBReweighter documentation:\n%s' % GBRW.__init__.__doc__ 
    Reweighter.__init__.__doc__ += '\n\nGBReweighter documentation:\n%s' % GBRW.__init__.__doc__ 
     
     
# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
        
# =============================================================================
##                                                                      The END 
# =============================================================================
