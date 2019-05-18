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
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.add_branch' )
else                       : logger = getLogger ( __name__     )
# =============================================================================
import ROOT
from   ostap.parallel.parallel import Task, WorkManager
# =============================================================================
## @class AddBranch
#  parallel adding of new branch for looong TChains
#  @see ostap/trees/trees.py
class AddBranch(Task) :
    """Add new branch  to loooong TChain in parallel
    """

    def __init__          ( self , branch_name  , function ) :
        self.branch_name = branch_name
        self.function    = function 
        self.output = ()
    
    def initializeLocal   ( self ) : self.output = () 
    def process           ( self , tree ) :

        import ostap.trees.trees

        if isinstace ( self.function ,  str ) : function = self.function
        else :
            ftype    = type ( function )
            function = ftype ( function )

        files = set() 
        tree.chain.add_new_branch ( self.branch_name , function , verbose = False ) 
        for f in tree.files : files.add ( f )

        ## list of processed  files 
        self.output = list ( files )
        
    ## merge results/datasets 
    def _mergeResults( self , result) :
        if not  self.output : self.output = result
        else :
            s = set()
            for r in self.output : s.add ( r )
            for r in      result : s.add ( r )
            s = list ( s )
            s.sort()
            self.output = tuple( s ) 


# =================================================================================
## Add new branch  to TChain in parallel
#  @see ROOT.TTree.add_new_branch
#  @code
#  chain = ....
#  chain.padd_new_branch ( 'new_branch' , 'px*py' ) 
#  @endcode
def add_new_branch ( chain , branch_name , function , verbose = True ) :
    """Add new branch for loong chain in parallel
    - see ROOT.TTree.add_new_branch
    >>> chain = ....
    >>> chain.padd_new_branch ( 'new_branch' , 'px*py' )     
    """
    from ostap.trees.trees import Chain
    from ostap.trees.trees import add_new_branch as _add_branch_ 
    
    if   isinstance ( chain , ROOT.TChain ) and 1 < len ( chain.files () ) : pass 
    elif isinstance ( chain , ROOT.TTree  ) : 
        return _add_branch_ ( chain , branch_name , function , verbose = False ) 
    
    ch    = Chain ( chain ) 
    
    task  = AddBranch ( branch_name ,  function  )
    wmgr  = Parallel.WorkManager ( silent = not verbose  )
    trees = ch.split ( max_files = 1  )
    
    wmgr.process( task , trees )
    
    nc = ROOT.TChain ( chain.name )
    for f in ch.files :  nc.Add ( f )
    
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
# The END 
# =============================================================================
