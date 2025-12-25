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
    if not hasattr ( numpy , 'float' ) :
        logger.warning ( 'No `numpy.float`... add it!')
        numpy.float = numpy.float64
    # =========================================================================
    ## @class Reweighter
    #  Helper class for reweighting using <code>GBReweighter</code>
    #  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
    #  @date   2021-09-22 
    class Reweighter(object) :
        """ Helper class for reweighting using `GBReweighter`
        - see hep_ml.reweight.GBReweighter
        """
        def __init__ ( self ,
                       original                ,
                       target                  ,
                       original_weight = None  ,
                       target_weight   = None  , 
                       silent          = False , **kwargs ) :
            """ Helper class for reweighting using `GBReweighter`
            - see hep_ml.reweight.GBReweighter
            """
            self.__kwargs     = kwargs 
            self.__reweighter = GBRW ( **kwargs ) 
            shape1 = original.shape
            shape2 = target.shape 
            
            assert len ( shape1 ) == len ( shape2 ) ,  "Inconsistent shapes: %s vs %s " % ( shape1 , shape2 )   
            assert len ( shape1 ) == 1 or shape1[-1] == shape2[-1] , \
                "Inconsistent shapes: %s vs %s " % ( shape1[-1] , shape2[-1] )   
                       
            assert ( original_weight is None ) or len ( original_weight ) == len ( original ) , \
                "Invalid length of ``original weights''"       
                 
            assert ( target_weight   is None ) or len ( target_weight   ) == len ( target   ) , \
                "Invalid length of ``target weights''"
        
            if not silent : 
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
        def kwargs ( self ) :
            """`kwargs`: actual configuration usef for clreationn of `GBReweighter`"""
            return self.__kwargs 
        
        @property
        def reweighter ( self ) :
            """`reweighter' : get the underlying reweighter object"""
            return self.__reweighter
         
        # =========================================================================
        def table    ( self , title = '' , prefix = '' ) :
            title = title if title else 'Reweighter configuration'
            from ostap.logger.utils import print_args 
            return print_args ( title = title , prefix = prefix , **self.__kwargs )
        
        def __str__  ( self ) : return self.table () 
        def __repr__ ( self ) : return self.table () 
        
        # =========================================================================
        ## Get/predict new weights for (new) original
        def weight ( self                   ,
                     original               ,
                     original_weight = None ) :
            """ Get/predict  new weights for (new) original
            """
            return self.reweighter.predict_weights (
                original        = original        ,
                original_weight = original_weight )
    
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
 
# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
        
# =============================================================================
##                                                                      The END 
# =============================================================================
