#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/trees/base.py
#  Couple of simple utilities forTTree/TChain
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
""" Couple of simple utilities fo rTTree/TChain
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'tree_nonzero'  , ## TTree is valid and non-empty 
    'tree_valid'    , ## Ttree is valid and non-empty
    'tree_path'     , ## full path of TTree in the TFile 
    'tree_branches' , ## list of branches 
    'tree_leaves'   , ## list of leaves 
    'chain_files'   , ## list of files for the TChain 
    'chain_nfiles'  , ## number of files for the TChain 
    'chain_nFiles'  , ## number of files for the TChain
    'tree_branch'   , ## get the branch by name 
    'tree_leaf'     , ## get the leaf by name 
) 
# =============================================================================
from   ostap.core.core        import valid_pointer
from   ostap.core.ostap_types import string_types 
import ROOT 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.trees.base' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Couple of simple basic utilities for TTree/TChain')
# =============================================================================
## check validity/emptiness  of TTree/TChain
#  require non-zero poniter and non-empty Tree/Chain
def tree_nonzero ( tree ) :
    """ Check validity/emptiness  of TTree/TChain
    - require non-zero pointer and non-empty Tree/Chain
    """
    return valid_pointer ( tree ) and 0 < tree.GetEntries() 

ROOT.TTree .__nonzero__ = tree_nonzero 
ROOT.TChain.__nonzero__ = tree_nonzero 
ROOT.TTree .__bool__    = tree_nonzero 
ROOT.TChain.__bool__    = tree_nonzero 

## ditto 
tree_valid = tree_nonzero

# =============================================================================
## Number of entries in the tree/chain
#  @code
#  tree = ...
#  len ( tree ) 
#  @endcode 
def _tree_len_ ( tree ) :
    """ Number of entries in the tree/chain 
    >>> tree = ... 
    >>> len ( tree ) 
    """
    return tree.GetEntries() if tree_valid ( tree ) else 0 

ROOT.TTree .__len__ = _tree_len_
ROOT.TChain.__len__ = _tree_len_ 

# ==============================================================================
## get a full path of the TTree object 
#  @code
#  rdir = ...
#  path = rdir.full_path 
#  path = rdir.fullpath 
#  @endcode
def tree_path ( tree ) :
    """ Get a full path of the directory
    >>> tree = ...
    >>> path = tree.full_path 
    >>> path = tree.fullpath 
    """
    assert valid_pointer ( tree ) and isinstance ( tree , ROOT.TTree ) , \
        "tree_path: Invalid `tree' argument!"
    ## 
    rdir = tree.GetDirectory()
    if valid_pointer ( rdir ) :
        path      = rdir.GetPath ()
        h , s , p = path.rpartition(':/')
        if p : return  '/'.join ( ( p , tree.GetName () ) )
    return tree.GetName()

    ## if not rdir : return tree.GetName()
    ## return os.path.join ( rdir.full_path , tree.GetName() ) 

## ROOT.TTree.full_path = property ( tree_path , None , None , tree_path . __doc__ ) 
## ROOT.TTree.top_dir   = property ( top_dir   , None , None , top_dir   . __doc__ )
## ROOT.TTree.fullpath  = property ( tree_path , None , None , tree_path . __doc__ ) 
## ROOT.TTree.topdir    = property ( top_dir   , None , None , top_dir   . __doc__ )

# =============================================================================
## get list of files used for the given chain
#  @code
#  >>> chain = ... ## get the files 
#  >>> chain_files ( chain )
#  >>> chain.files () 
#  @endcode  
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-02-04
def chain_files ( chain ) :
    """ Get the list of files used for the chain    
    >>> chain = ..
    >>> chain_files( chain )
    >>> chain.files()
    """
    if tree_valid ( chain ) :
        if   isinstance ( chain , ROOT.TChain ) :
            return tuple ( i.GetTitle() for i in chain.GetListOfFiles() )
        elif isinstance ( chain , ROOT.TTree  ) :
            td = chain.topdir
            if isinstance ( td , ROOT.TDirectoryFile ) :
                return td.GetName () ,             
    return ()

ROOT.TTree.  files  = chain_files 
ROOT.TChain. files  = chain_files 

# =============================================================================
## get number of files used for the chain 
#  @code
#  >>> chain = ... ## get the files 
#  >>> n = chain_nfiles() 
#  >>> n = chain.nFiles() 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-02-04
def chain_nfiles ( chain ) :
    """ Get number of of files used for the chain
    
    >>> chain = ... ## get the files 
    >>> n = chain_nfiles ( chain ))
    >>> n = chain.nFiles()
    """
    if tree_valid ( chain ) and isinstance ( chain , ROOT.TTree  ) :
        if  isinstance ( chain , ROOT.TChain ) : return len ( chain.GetListOfFiles() )
        else : 
            td = chain.topdir
            if td and isinstance ( td , ROOT.TDirectoryFile ) : return 1       
    return 0

ROOT.TChain. nFiles = chain_nfiles
ROOT.TChain. nfiles = chain_nfiles
chain_nFiles        = chain_nfiles

# =============================================================================
## get the branches for the given tree/chain
#  @see TTree
#  @code
#  >>> tree = ...
#  >>> lst = tree.branches()
#  >>> for b in lst : print b
#  >>> lst = tree.branches( '.*(Muon).*' , re.I )
#  >>> for b in lst : print b
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-02-04
def tree_branches ( tree , pattern = '' , *args ) :
    """ Get the list of branch names
    
    >>> tree = ...
    >>> lst = tree.branches()
    >>> for b in lst : print b
    >>> lst = tree.branches( '.*(Muon).*' , re.I )
    >>> for b in lst : print b
    
    """
    assert valid_pointer ( tree ) and isinstance ( tree , ROOT.TTree ) , \
        "tree_branches: Invalid `tree' argument!"
    ##
    vlst = tuple ( sorted ( b.GetName() for b in tree.GetListOfBranches() ) )
    if not vlst or not pattern : return vlst 
    
    import re
    
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        c    = re.compile ( pattern , *args )
        lst  = sorted ( v for v in vlst if c.match ( v ) ) 
        return lst 
    except : # ================================================================
        # =====================================================================
        logger.error ('branches: exception is caught, skip it' , exc_info = True ) 

    return vlst 

ROOT.TTree.branches = tree_branches

# =============================================================================
## get the leaves for the given tree/chain
#  @see TTree
#  @code
#
#  >>> tree = ...
#  >>> lst = tree_leaves( tree )
#  >>> for l in lst : print l
#
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-02-04
def tree_leaves ( tree , pattern = '' , *args ) :
    """ Get the list of leaves names        
    
    >>> tree = ...
    >>> lst = tree_leaves( tree )
    >>> for l in lst : print l
    >>> lst = tree.leaves( '.*(muon).*', re.I )
    >>> for l in lst : print l
    """
    assert valid_pointer ( tree ) and isinstance ( tree , ROOT.TTree ) , \
        "tree_leaves: Invalid `tree' argument!"
    
    vlst =  tuple ( sorted ( v.GetName()  for v in tree.GetListOfLeaves() ) ) 
    if not vlst or not pattern : return vlst 

    if isinstance ( pattern , string_types ) : pattern  = [ pattern ]

    lst = set()
    for p in pattern :
        import re
        # =====================================================================
        try : # ===============================================================
            # =================================================================
            c    =  re.compile ( p , *args )
            vars = [ v for v in vlst if c.match ( v ) ]
            lst  = lst | set ( vars  ) 
        except :
            logger.error ('leaves("%s"): exception is caught, use all ' % p  , exc_info = True )
            
    return tuple ( sorted ( lst ) ) 

ROOT.TTree.leaves   = tree_leaves

# ==============================================================================
## Get the leaf with the certain name 
def tree_leaf ( tree , leaf ) :
    """ Get the leaf with certain name:
    >>> tree = ...
    >>> l = tree_leaf('pt') 
    """
    for v in tree.GetListOfLeaves() : 
        if leaf == v.GetName() : return v
    return None

ROOT.TTree.leaf   = tree_leaf

# ==============================================================================
## Get the branch with the certain name 
def tree_branch( tree , branch ) :
    """ Get the branchwith certain name:
    >>> tree = ...
    >>> l = tree_branch ('pt') 
    """
    for v in tree.GetListOfBranches() : 
        if branch== v.GetName() : return v
    return None

ROOT.TTree.branch = tree_branch

# ===============================================================================
_decorated_classes_ = (
    ROOT.TTree  ,
    ROOT.TChain ,   
)

_new_methods_ = (
    ##
    tree_nonzero            , ## TTree is valid and non-empty 
    tree_valid              , ## Ttree is valid and non-empty
    tree_path               , ## full path of TTree in TFile 
    chain_files             , ## list of files for the TChain 
    chain_nfiles            , ## number of files for the TChain 
    chain_nFiles            , ## number of files for the TChain
    tree_branches           , ## list of branches for the TTree
    tree_leaves             , ## list of leaves for the TTree
    tree_branch             , ## get the branch by name 
    tree_leaf               , ## get the leaf by name 
    ##
    ROOT.TTree .__nonzero__ , 
    ROOT.TChain.__nonzero__ , 
    ROOT.TTree .__bool__    , 
    ROOT.TChain.__bool__    ,
    ROOT.TTree .__len__     , 
    ROOT.TChain.__len__     ,
    ##
    ROOT.TChain. files      ,
    ROOT.TChain.nfiles      ,
    ROOT.TChain.nFiles      ,    
    ## 
    ROOT.TTree.branches     ,
    ROOT.TTree.leaves       ,
    ROOT.TTree.branch       ,
    ROOT.TTree.leaf         ,
)
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
