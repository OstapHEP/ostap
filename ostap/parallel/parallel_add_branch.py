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
        push_2chain ( chain , self.recipe_args , progress = False , report = False )
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
    
    logger.info ( 'PARALLEL-ADD/0' )

    assert valid_pointer ( chain ) and isinstance ( chain , ROOT.TTree ) , \
        "TTree* is invalid!"

    logger.info ( 'PARALLEL-ADD/1' )
    
    from ostap.trees.trees import Chain, prepare_branches, push_2chain  

    logger.info ( 'PARALLEL-ADD/2' )

    if   isinstance ( chain , ROOT.TChain ) and 1 < len ( chain.files () ) : pass 
    elif isinstance ( chain , ROOT.TTree  ) :        
        from ostap.trees.trees import add_new_branch as _add_branch_ 
        return _add_branch_ ( chain , branch , progress = progress , report = report  ) 

    logger.info ( 'PARALLEL-ADD/3' )

    check = Checker()

    logger.info ( 'PARALLEL-ADD/4' )
    
    verbose = kwargs.pop ( 'verbose' , False )

    if True : ## verbose : 
        title = 'All input arguments'
        try : 
            logger.info ( '%s:\n%s' % ( title , check.pickling_table ( branch , prefix = '# ' , **kwargs ) ) )
        except :
            logger.fatal ( 'UNCAUGHT exception/1' , exc_info = True )

    logger.info ( 'PARALLEL-ADD/5' )
        
    # ========================================================================
    ## process all rguments  
    args , expected , kw , keeps = prepare_branches ( chain , branch , **kwargs ) 

    logger.info ( 'PARALLEL-ADD/6-1' )
    from ostap.logger.utils import print_args
    logger.info ( 'ARGUMENTS-HERE\n%s' % print_args ( *args , **kw ) ) 
     
    logger.info ( 'PARALLEL-ADD/6-2' )
    
    if True : ## verbose : 
        title = 'All processed arguments'
        try : 
            table = check.pickling_table ( *args  , prefix = '# ' , **kwargs ) 
            logger.info ( '%s:\n%s' % ( title , table ) )
        except :
            logger.fatal ( 'UNCAUGHT exception/2' , exc_info = True )
            

    logger.info ( 'PARALLEL-ADD/7' )
        
    # =========================================================================
    ## perfect! 
    if   check.pickles_all ( *args )  :
        logger.info ( 'PARALLEL-ADD/8' )
        task = AddNewBranch ( *args )
        logger.info ( 'PARALLEL-ADD/9' )
    ## some 'args' are not pickable! 
    elif check.pickles_all ( branch , **kwargs ) :
        logger.info ( 'PARALLEL-ADD/10' )
        task = AddUnpreparedBranch ( branch , **kwargs )
        logger.info ( 'PARALLEL-ADD/11' )        
    elif backup :         
        # =====================================================================
        logger.info( 'PARALLEL-ADD/12' )
        table = check.pickling_table ( branch , *args , prefix = '# ' , **kwargs )
        logger.info ( 'PARALLEL-ADD/13' )

        logger.warning ( 'Not all arguments are pickable:\n%s' % table )        
        logger.warning ( 'Switch to sequential (SLOW) processing' )
        ## 
        chain      = push_2chain ( chain , args , progress = progress , report = report )            
        missing    = sorted ( branch for branch in expected if not branch in chain  )
        if missing : logger.warning ( 'Missing expected brnaches: %s' % ( ', '.join ( m for m in missing ) ) )
        logger.attention ( "A (SLOQ) sequentional processing was used..." ) 
        return chain                           
    else  : # =================================================================
        # =====================================================================
        print ( 'PARALLEL-ADD/15' )
        table = check.pickling_table ( branch , prefix = '# ' , **kwargs ) 
        logger.error ( 'Arguments are pickable:\n%s' % table )
        return

    print ( 'PARALLEL-ADD/15' )
    
    files    = chain.files () 
    cname    = chain.name
    ch       = Chain ( chain ) 
    branches = set   ( chain.branches() ) | set ( chain.leaves() ) 

    wmgr     = WorkManager ( silent = not progress , **kw )
    trees    = ch.split    ( max_files = 1  )
    
    print ( 'PARALLEL-ADD/16' )
    
    wmgr.process ( task , trees )

    print ( 'PARALLEL-ADD/17' )
    
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

    print ( 'PARALLEL-ADD/18' )
            
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
