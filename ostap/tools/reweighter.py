#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with utilities for reweigting using GBReweighetr 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2021-09-22
# =============================================================================
"""Module with utilities for reweighting using GBReweighter"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'Reweighter' , 
    ) 
# =============================================================================
import warnings
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.tools.reweighter' )
else                       : logger = getLogger( __name__ )
# =============================================================================
## @class Reweighter
#  Helper class for reweighting using <code>GBReweighter</code>
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2021-09-22 
class Reweighter(object) :
    """ Helper class for reweighting using `GBReweighter`
    - see hep_ml.reweight.GBReweighter
    """
    
    def __init__ ( self , **kwargs )  :

        from hep_ml.reweight import GBReweighter as GBRW 
        self.__reweighter = GBRW ( **kwargs ) 
        self.__variables  = ()
        self.__nvars      = 0 
        
    @property
    def reweighter ( self ) :
        """``reweighter'' : get the underlying reweighter object"""
        return self.__reweighter

    # =========================================================================
    ## the main method to train the underlying BDT
    #  @param original  2D-array of data that needs to be reweighted ("MC")
    #  @param target    2D-array of data that is a rewighting target ("data")
    #  @param original_weight 1D array of input weights for "MC" dataset
    #  @param target_weight   1D array of input weights for "data" dataset
    #  @see hep_ml.reweight.GBReweighter
    def reweight  ( self                       ,    
                    original                   ,
                    target                     ,
                    original_weight     = None ,
                    target_weight       = None ) :
        """The main method to train the underlying BDT
        - see hep_ml.reweight.GBReweighter        
        """
        
        assert original.shape[-1] == target.shape[-1] , \
               "Inconsistent shapes: %s vs %s " % ( target.shape , original.shape )          
        assert ( not original_weight ) or len ( original_weight ) == len ( original ) , \
               "Invalid length of ``original weights''"        
        assert ( not target_weight   ) or len ( target_weight   ) == len ( target   ) , \
               "Invalid length of ``target weights''"
        
        self.__nvars = original.shape[-1]
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")    
            self.reweighter.fit ( original ,
                                  target   ,
                                  original_weight = original_weight if original_weight else None , 
                                  target_weight   = target_weight   if target_weight   else None )
            
    def weight ( self                   ,
                 original               ,
                 original_weight = None ) :

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")    
            return self.reweighter.predict_weights (
                original        = original         ,
                original_weight = original_weight  ) 
        
# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
