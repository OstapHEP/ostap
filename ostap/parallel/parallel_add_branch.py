#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/parallel/parallel_add_branch.py
#  (parallel) Add new branch to looong TChain
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
# =============================================================================
"""(parallel) Add new branch to looong TChain
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'add_new_branch'  , ## add new branch to loooong TChain in parallel
    ) 
# =============================================================================
from   ostap.core.meta_info            import root_info, python_info 
from   ostap.parallel.parallel         import Task, WorkManager
import ostap.parallel.parallel_statvar
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.add_branch' )
else                       : logger = getLogger ( __name__     )
# =============================================================================
## @class AddBranch
#  parallel adding of new branch for looong TChains
#  @see ostap/trees/trees.py
class AddBranch(Task) :
    """Add new branch to loooong TChain in parallel
    """
    def __init__          ( self , branch_name  , function ) :
        self.branch_name = branch_name
        self.function    = function 
        self.__output    = ()
    
    def initialize_local  ( self )                : self.__output = () 
    def initialize_remote ( self , jobid = -1   ) : self.__output = () 
    def process           ( self , jobid , tree ) :

        import ostap.trees.trees

        files = []  
        tree.chain.add_new_branch ( self.branch_name , self.function , verbose = False , report = False ) 
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
    
# =================================================================================
## Add new branch  to TChain in parallel
#  @see ROOT.TTree.add_new_branch
#  @code
#  chain = ....
#  chain.padd_new_branch ( 'new_branch' , 'px*py' ) 
#  @endcode
def add_new_branch ( chain          ,
                     branch_name    ,
                     function       , 
                     verbose = True ,
                     report  = True , **kwargs ) :
    """Add new branch for loong chain in parallel
    - see ROOT.TTree.add_new_branch
    >>> chain = ....
    >>> chain.padd_new_branch ( 'new_branch' , 'px*py' )     
    """
    from ostap.trees.trees import Chain
    from ostap.trees.trees import add_new_branch as _add_branch_ 

    if ( 6 , 18 ) <= root_info < ( 6, 19) and python_info < ( 3 , 0 ) :
        if verbose : logger.info ( 'Switch to sequential processing...' ) 
        return _add_branch_ ( chain , branch_name , function , verbose = verbose , report = report  ) 
    
    if   isinstance ( chain , ROOT.TChain ) and 1 < len ( chain.files () ) : pass 
    elif isinstance ( chain , ROOT.TTree  ) : 
        return _add_branch_ ( chain , branch_name , function , verbose = verbose , report = report  ) 

    cname    = chain.name 
    ch       = Chain ( chain ) 
    branches = set   ( chain.branches() ) | set ( chain.leaves() ) 
    
    keep_it  = branch_name , function

    task     = AddBranch   ( branch_name ,  function  )
    wmgr     = WorkManager ( silent = not verbose , **kwargs )
    trees    = ch.split    ( max_files = 1  )

    import ostap.parallel.parallel as P
    print ( 'PARALLEL'  , P.worker    , P.WorkManager , P.dill , P.DILL_PY3_issue )
    print ( 'PICKLES/F' , function    , P.pickles ( function    ) )
    print ( 'PICKLES/B' , branch_name , P.pickles ( branch_name ) )
    print ( 'PICKLES/T' , task        , P.pickles ( task  ) )
    
    wmgr.process ( task , trees )

    nc = ROOT.TChain ( cname )
    for f in ch.files :  nc.Add ( f )

    if report : 
        new_branches = set ( nc.branches () ) | set ( nc.leaves() )
        new_branches = sorted ( new_branches - branches )
        if new_branches : 
            n = len ( new_branches )  
            if 1 == n  : title = 'Added %s branch to TChain'   % n
            else       : title = 'Added %s branches to TChain' % n
            table = nc.table ( new_branches , title = title , prefix = '# ' )
            logger.info ( '%s:\n%s' % ( title , table ) ) 
            
    return nc 

ROOT.TTree .padd_new_branch = add_new_branch 
ROOT.TChain.padd_new_branch = add_new_branch

# =============================================================================
_decorated_classes_ = (
    ROOT.TTree  ,
    ROOT.TChain ,    
    )

_new_methods_       = (
    ROOT.TTree .padd_new_branch ,
    ROOT.TChain.padd_new_branch ,
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
#                                                                       The END 
# =============================================================================
