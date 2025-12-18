#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/parallel/parallel_add_rwweighting.py
#  (parallel) Add reweighting informaton  to looong TChain
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2023-06-15
# =============================================================================
"""(parallel) Add new branch to looong TChain
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2023-06-12"
__all__     = (
    'add_reweighting'  , ## add new branch to loooong TChain in parallel
    ) 
# =============================================================================
from   ostap.utils.basic       import numcpu, typename  
from   ostap.parallel.parallel import Task  , WorkManager
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.add_reweighting' )
else                       : logger = getLogger ( __name__     )
# =============================================================================
## @class AddReweighting
#  parallel adding of reweighting info for looong TChains
#  @see ostap/trees/trees.py
class AddReweighting(Task) :
    """ Add reweighting information to loooong TChain in parallel
    """
    def __init__          ( self , branch_name  , reweighter ) :
        self.branch_name = branch_name
        self.reweighter  = reweighter
        self.__output    = ()
    
    def initialize_local  ( self )                : self.__output = () 
    def initialize_remote ( self , jobid = -1   ) : self.__output = () 
    def process           ( self , jobid , tree ) :

        import ostap.trees.trees
        import ostap.tools.reweight
        
        files = []  
        tree.chain.add_reweighting ( weighter = self.reweighter  ,
                                     name     = self.branch_name ,
                                     verbose  = False            ,
                                     report   = False            ) 
        for f in tree.files :
            if not f in files : files.append ( f )

        ## list of processed  files
        self.__output = tuple ( files )

        return self.__output 
        
    ## merge results/datasets 
    def merge_results( self , result , jobid = -1 ) :
        
        if not self.__output : self.__output = result
        else                 :
            processed = sorted ( self.__output + result ) 
            self.__output = tuple ( processed ) 
            
    ## get the results 
    def results ( self ) : return self.__output

# =============================================================================================
## @class AddReweightingDS
#  parallel adding of reweighting info for dataset 
#  @see ostap/trees/trees.py
class AddReweightingDS(Task) :
    """ Add reweighting information to dataset
    """
    def __init__          ( self , branch_name  , reweighter ) :
        self.branch_name = branch_name
        self.reweighter  = reweighter
        self.__output    = ()
    
    def initialize_local  ( self )                : self.__output = () 
    def initialize_remote ( self , jobid = -1   ) : self.__output = () 
    def process           ( self , jobid , dset ) :

        import ostap.tools.reweight
        import ostap.fitting.dataset
        
        files = []  
        dset.add_reweighting ( weighter = self.reweighter  ,
                               name     = self.branch_name ,
                               progress = False            , 
                               report   = False            ,
                               parallel = False            )

        ## list of processed  files
        self.__output = dset 

        return self.__output 
        
    ## merge results/datasets 
    def merge_results( self , result , jobid = -1 ) :
        
        if not self.__output : self.__output = result
        else                 :
            self.__output += result 
            result.clear ()     ## ATTTENTION!!! 
            del result
            
    ## get the results 
    def results ( self ) : return self.__output
    
# =================================================================================
## Add reweighting to TChain in parallel
#  @see ROOT.TTree.add_new_branch
#  @code
#  chain = ....
#  chain.padd_reweigting ( weighter , name = 'pty_weight' ) 
#  @endcode
def add_reweighting ( chain            ,
                      weighter         , 
                      name             ,
                      verbose  = True  ,
                      report   = True  , **kwargs ) :
    """ Add reweihting info for loong chain in parallel
    - see ROOT.TTree.add_reweighting
    >>> chain = ....
    >>> chain.padd_reweighting ( 'new_branch' , 'px*py' )     
    """
    from ostap.trees.utils    import Chain
    from ostap.tools.reweight import tree_add_reweighting as _add_reweighting_ 
    
    if   isinstance ( chain , ROOT.TChain ) and 1 < chain.nFiles : pass 
    elif isinstance ( chain , ROOT.TTree  ) : 
        return _add_reweighting_ ( chain , weighter , name , verbose = verbose , report = report  ) 
    
    ch       = Chain ( chain )
    cname    = ch.name 
    branches = set   ( chain.branches() ) | set ( chain.leaves() ) 
    
    task     = AddReweighting ( name  ,  weighter )
    wmgr     = WorkManager    ( silent    = not verbose , **kwargs )
    trees    = ch.split       ( max_files = 1  )
    
    wmgr.process ( task , trees )
    
    nc = ROOT.TChain ( chain.name )
    for f in ch.files :  nc.Add ( f )

    if report : 
        new_branches = set ( nc.branches () ) | set ( nc.leaves() )
        new_branches = sorted ( new_branches - branches )
        if new_branches : 
            n = len ( new_branches )  
            if 1 == n  : title = 'Added %s branch to TChain(%s)'   % ( n , cname ) 
            else       : title = 'Added %s branches to TChain(%s)' % ( n , cname ) 
            table = nc.table ( new_branches , title = title , prefix = '# ' )
            logger.info ( '%s:\n%s' % ( title , table ) ) 
            
    return nc 

ROOT.TTree .padd_reweighting = add_reweighting 
ROOT.TChain.padd_reweighting = add_reweighting 

ROOT.RooDataSet .padd_reweighting = add_reweighting 

# =============================================================================
_decorated_classes_ = (
    ROOT.TTree      ,
    ROOT.TChain     ,
    ROOT.RooDataSet    
    )

_new_methods_       = (
    ROOT.TTree .padd_reweighting      , 
    ROOT.TChain.padd_reweighting      , 
    ROOT.RooDataSet .padd_reweighting )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
#                                                                       The END 
# =============================================================================
