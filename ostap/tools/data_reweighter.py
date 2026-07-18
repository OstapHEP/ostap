#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with utilities for reweigting using GBReweighetr
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
    'DataReweighter' , 
    ) 
# =============================================================================
from   ostap.utils.basic        import typename
from   ostap.math.math_base     import weight_trivial
from   ostap.stats.statvars     import data_slice 
from   ostap.trees.cuts         import vars_and_cuts 
from   ostap.utils.progress_bar import progress_bar 
import ostap.trees.trees 
import ROOT 
# ============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.tools.reweighter' )
else                       : logger = getLogger( __name__ )
# =============================================================================
## @class DataReweighter
# Helper class to perform the reweighting using tree/source-interface
class DataReweighter(object) : 
    """ Helper class to perform (GB)reweighting using tree/source-interface
    """
    def __init__ ( self                      , 
                   reweighter_type           , * , 
                   original                  ,     
                   target                    ,  
                   target_variables          ,
                   original_variables = None , 
                   target_weight      = ''   , 
                   original_weight    = ''   , 
                   silent             = True , 
                   progress           = True , **params ) :
        
        ## 
        tvars, tweight, _ = vars_and_cuts ( target_variables , target_weight )
            
        ovars = original_variables if original_variables else tvars
        ovars, oweight, _ = vars_and_cuts ( ovars            , original_weight )
            
        self.__silent   = True if silent else False
        self.__progress = True if progress else False   
        self.__tvars    = tuple ( tvars )  
        self.__ovars    = tuple ( ovars )
        self.__tweight  = tweight
        self.__oweight  = oweight 

        ## 
        tdata, tw = data_slice ( target  , tvars, tweight , structured = False , progress = self.progress )
        odata, ow = data_slice ( original, ovars, oweight , structured = False , progress = self.progress ) 

        if not weight_trivial ( tw ) : tw /= numpy.sum ( tw )
        if not weight_trivial ( ow ) : ow /= numpy.sum ( ow )
        
        ## create the actual reweighter
        self.__rw = reweighter_type ( original        = odata , 
                                      target          = tdata , 
                                      original_weight = ow    , 
                                      target_weight   = tw    , 
                                      silent          = self.silent , **params  )

    @property
    def reweighter ( self ) :
        """`reweighter` : get the ctual Reweighter object
        """
        return self.__rw
    
    @property 
    def config ( self ) :
        """`config` : Reweighter configuraton"""
        conf = {} 
        conf.update ( self.__rw.config )
        conf [ 'progress'           ] = progress 
        conf [ 'Reweighter type'    ] = typename ( self.__rw )
        conf [ 'target-variables'   ] = self.__tvars
        conf [ 'original-variables' ] = self.__ovars
        if self.__tweight : conf [ 'target-weight'   ] = self.__tweight 
        if self.__oweight : conf [ 'originam-weight' ] = self.__oweight 
                   
    @property
    def silent ( self ) : 
        """`silent`: silent processing?
        """
        return self.__silent
    
    @property
    def progress ( self ) : 
        """`progress` : show progress bar?
        """
        return self.__progress 

    # =======================================================================================
    ## Get the weights for original, and add the weight back to original           
    def reweight ( self      ,
                   original  , * ,
                   name      = 'weight' ,
                   variables = None     ) :
        """ Get the weights for original, and add the weight back to original
        """
            
        assert isinstance ( original ,  ( ROOT.TTree , ROOT.RooDataSet ) ) , \
            "Invalid `original` type %s" % typename ( original ) 

        
        ## (1) processing  TChain with several files?    
        if isinstance ( original , ROOT.TChain ) and 1 < original.nFiles : 
            files         = original.files
            cname         = original.name  
            
            _progress = self.progress
            _silent   = self.silent 
            
            branches = ( set ( original.branches() ) | set ( original.leaves() ) ) if not self.silent else set()

            chain_progress = self.progress and ( 5 <= len ( files ) )
            if chain_progress : 
                self.__progrees = False 
                self.__silent   = True  
                
            for f in progress_bar ( files , silent = not chain_progress , description = 'Files:' ) : 
                ch = ROOT.TChain ( original.name , files = f )
                ch = self.reweight ( ch , name = name )
                
            chain = ROOT.TChain ( cname , files = files )
            if chain_progress : 
                self.__progress = _progress
                self.__silent   = _silent 
                
            if not self.silent :
                new_branches = sorted ( ( set ( chain.branches() ) | set ( chain.leaves() ) ) - branches )
                if new_branches :
                    n = len ( new_branches )
                    if 1 >= n : title = "Added %s branch to TChain(%s)"   % ( n , cname )
                    else      : title = "Added %s branches to TChain(%s)" % ( n , cname )
                    table = chain.table ( new_branches , title = title , prefix = '# ' )
                    logger.info ( '%s:\n%s' % ( title , table ) )
                else :
                    logger.warning ( "No branches are added to TChain(%s)" % cname )
                
            return chain 
            
        ## (2) processing single tree/source :
        
        ovars = self.__ovars 
        if not variables is None :
            ovars , _ , _ = vars_and_cuts ( variables , ''  )
            assert len ( ovars ) == len ( self.__ovars ) , "Invalid variables!"
            if not silent and not ovars == self.__ovars :
                logger.warning  (  "Variables to be used: [%s]" % ','.join ( ovars ) )
                
        ## ATTENTION: no weight here!
        odata, _ = data_slice ( original , ovars, '' , structured = False , progress = self.progress ) 
            
        the_weight = self.__rw ( odata ) 
        return original.add_new_buffer ( the_weight                 ,
                                         name     = name            ,
                                         report   = not self.silent , 
                                         progress = self.progress and not self.silent ) 

    # ===============================================================================================
    ## Calculate the weights for "original" to add them to the tree/source
    def __call__  ( self , original , name = 'weight' )  : 
        """ Calculate weights for `original' add add tem to the tree/source
        """
        return self.reweight ( original , name = name  )
    
