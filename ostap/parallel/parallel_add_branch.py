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
from   ostap.core.core                 import valid_pointer 
from   ostap.utils.basic               import loop_items 
from   ostap.parallel.parallel         import Task, WorkManager, Checker 
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
#  parallel adding of new (prepared) branch for looong TChains
#  @see ostap/trees/trees.py
class AddNewBranch(Task) :
    """ Add new (prepared) branch to loooong TChain in parallel
    """
    def __init__ ( self , *recipe ) :
        self.recipe_args = recipe
        self.__output    = ()
    
    def initialize_local  ( self )                : self.__output = () 
    def initialize_remote ( self , jobid = -1   ) : self.__output = () 
    def process           ( self , jobid , tree ) :
        
        from ostap.trees.trees import push_2chain
        
        files = []
        chain = tree.chain
        push_2chain ( chain , *self.recipe_args , progress = False , report = False )
        for f in tree.files :
            if not f in files : files.append ( f )

        ## list of processed  files
        self.__output = tuple ( sorted ( files ) ) 

        return self.__output 
        
    ## merge results/datasets 
    def merge_results( self , result , jobid = -1 ) :
        
        if not self.__output : self.__output = result
        else                 :
            processed = sorted ( self.__output + result ) 
            self.__output = tuple ( processed ) 
            
    ## get the results 
    def results ( self ) : return self.__output
    
# =============================================================================
## @class AddUnpreparedBranch
#  parallel adding of new (prepared) branch for looong TChains
#  @see ostap/trees/trees.py
class AddUnpreparedBranch(AddNewBranch) :
    """ Add new (prepared) branch to loooong TChain in parallel
    """
    def __init__ ( self , branch , **kwargs ) :
        AddNewBranch.__init__ ( self ) 
        self.branch_arg    = branch 
        self.recipe_kwargs = kwargs.copy() 

    def process           ( self , jobid , tree ) :
        
        from ostap.trees.trees import prepare_branches
        chain = tree.chain
        args, expected , kw , keep = prepare_branches ( chain , self.branch_arg , **self.recipe_kwargs )
        self.recipe_args = args
        return AddNewBranch.process ( self , jobid , tree ) 
    
# =================================================================================
## Add new branch  to TChain in parallel
#  @see ROOT.TTree.add_new_branch
#  @code
#  chain = ....
#  chain.padd_new_branch ( 'new_branch' , 'px*py' ) 
#  @endcode
def add_new_branch ( chain           ,
                     branch          ,
                     progress = True ,
                     report   = True ,
                     backup   = True , **kwargs ) :
    """ Add new branch for loong chain in parallel
    - see ROOT.TTree.add_new_branch
    >>> chain = ....
    >>> chain.padd_new_branch ( 'new_branch' , 'px*py' )     
    """
    
    assert valid_pointer ( chain ) and isinstance ( chain , ROOT.TTree ) , \
        "TTree* is invalid!"

    from ostap.trees.trees import Chain, prepare_branches, push_2chain  

    if   isinstance ( chain , ROOT.TChain ) and 1 < len ( chain.files () ) : pass 
    elif isinstance ( chain , ROOT.TTree  ) :        
        from ostap.trees.trees import add_new_branch as _add_branch_ 
        return _add_branch_ ( chain , branch , progress = progress , report = report , **kwargs ) 

    ## check (un)picleability 
    check = Checker()
    
    verbose = kwargs.pop ( 'verbose' , False )

    if verbose : 
        title = 'All input arguments'
        logger.info ( '%s:\n%s' % ( title , check.pickling_table ( branch , prefix = '# ' , **kwargs ) ) )

    # ========================================================================
    ## process all rguments  
    args , expected , kw , keeps = prepare_branches ( chain , branch , **kwargs ) 


    if verbose : 
        title = 'All processed arguments'
        table = check.pickling_table ( *args  , prefix = '# ' , **kwargs ) 
        logger.info ( '%s:\n%s' % ( title , table ) )

    # =========================================================================
    ## perfect! 
    if   check.pickles_all ( *args )  :
        task = AddNewBranch ( *args )
    ## some 'args' are not pickable! 
    elif check.pickles_all ( branch , **kwargs ) :
        task = AddUnpreparedBranch ( branch , **kwargs )
    elif backup :         
        # =====================================================================
        table = check.pickling_table ( branch , *args , prefix = '# ' , **kwargs )
        logger.warning ( 'Not all arguments are pickable:\n%s' % table )        
        logger.warning ( 'Switch to (SLOW) sequential processing' )
        ## 
        chain      = push_2chain ( chain , *args , progress = progress , report = report )            
        missing    = sorted ( branch for branch in expected if not branch in chain  )
        if missing : logger.warning ( 'Missing expected brnaches: %s' % ( ', '.join ( m for m in missing ) ) )
        logger.attention ( "A (SLOQ) sequentional processing was used..." ) 
        return chain                           
    else  : # =================================================================
        # =====================================================================
        table = check.pickling_table ( branch , prefix = '# ' , **kwargs ) 
        logger.error ( 'Arguments are *NOT* pickleable:\n%s' % table )
        return
    
    files    = chain.files () 
    cname    = chain.name
    ch       = Chain ( chain ) 
    branches = set   ( chain.branches() ) | set ( chain.leaves() ) 

    wmgr     = WorkManager ( silent = not progress , **kw )
    trees    = ch.split    ( max_files = 1  )
    
    wmgr.process ( task , trees )

    nc = ROOT.TChain ( cname )
    for f in files : nc.Add ( f )

    if report : 
        new_branches = sorted ( ( set ( nc.branches () ) | set ( nc.leaves() ) ) - set ( branches ) ) 
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
