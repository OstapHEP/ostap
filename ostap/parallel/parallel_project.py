#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/parallel/parallel_project.py
#  Paralllel "project" from loooong chain/trees objects   
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
# =============================================================================
""" Paralllel `project' from loooong chain/trees objects   
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'cproject' , ## parallel project from looong TChain
    'tproject' , ## parallel project from looong TTree 
    ) 
# =============================================================================
from   ostap.parallel.parallel_statvar parallel_project
FIRST_ENTRY, LAST_ENTRY 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.parallel_project' )
else                       : logger = getLogger ( __name__     )
# =============================================================================  
## make a projection of the loooooooong chain into histogram using
#  multiprocessing functionality for per-file parallelisation
#  @code
#  >>> chain = ... ## large chain
#  >>> histo = ... ## histogram template 
#  >>> project        ( chain , histo , 'mass' , 'pt>10' )
#  >>> chain.pproject ( histo , 'mass' , 'pt>0' ) ## ditto 
#  >>> chain.cproject ( histo , 'mass' , 'pt>0' ) ## ditto 
#  @endcode
#  For 12-core machine, clear speedup factor of about 8 is achieved 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
def  cproject ( chain                    ,
                target                   ,
                expressions              ,
                cuts       = ''          ,
                first      = FIRST_ENTRY , 
                last       = LAST_ENTRY  ,
                progress   = True    , 
                use_frame  = True    ,                 
                chunk_size = -1      ,
                max_files  =  1      ,
                prohress   = True    , 
                silent     = False   , **kwargs ) :
    """ Make a projection of the loooong chain into histogram
    >>> chain = ... ## large chain
    >>> histo = ... ## histogram template 
    >>> cproject        ( chain , histo , 'mass' , 'pt>10' )
    >>> chain.ppropject ( histo , 'mass' , 'pt>0' ) ## ditto 
    >>> chain.cpropject ( histo , 'mass' , 'pt>0' ) ## ditto     
    For 12-core machine, clear speedup factor of about 8 is achieved     
    """
    return parallel_project ( chain                   ,
                              target                  ,
                              expressions             ,
                              cuts       = cuts       ,
                              first      = first      ,
                              last       = last       ,
                              progress   = progress   , 
                              use_frame  = use_frame  ,
                              chunk_size = chunk_size ,
                              max_files  = max_files  , 
                              silent     = silent     , **kwargs )

# =============================================================================  
## make a projection of the loooooooong tree into histogram using
#  multiprocessing functionality for per-file parallelisation
#  @code
#
#  >>> tree  = ... ## large tree 
#  >>> histo = ... ## histogram template 
#  >>> tproject ( tree , histo , 'mass' , 'pt>10' , maxentries = 1000000 )
#  >>> tree.pproject ( histo , 'mass' , 'pt>10' ) ## ditto 
#  @endcode
#  - significant gain can be achieved for very large ttrees with complicated expressions and cuts
#  - <code>maxentries</code> parameter should be rather large
#  @param tree       the tree
#  @param histo      the histogram
#  @param what       variable/expression/varlist to be projected
#  @param cuts       selection/weighting criteria 
#  @param nentries   number of entries to process  (>0: all entries in th tree)
#  @param first      the first entry to process
#  @param maxentries chunk size for parallel processing 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-09-23
def  tproject ( tree                 ,   ## the tree 
                target               ,   ## histogram 
                expressions          ,   ## variable/expression/list to be projected 
                cuts       = ''      ,   ## selection/weighting criteria
                first      = FIRST_ENTRY , 
                last       = LAST_ENTRY  ,
                progress   = True    , 
                use_frame  = True    ,
                chunk_size = 1000000 ,   ## chunk size 
                max_files  = 1       ,   ## not-used .... 
                silent     = False   , **kwargs ) : ## silent processing 
    """ Make a projection of the loooong tree into histogram
    >>> tree  = ... ## large chain
    >>> histo = ... ## histogram template 
    >>> tproject ( tree , histo , 'mass' , 'pt>10' )    
    >>> tree.pproject ( histo , 'mass' , 'pt>10' )    ## ditto 
    - significant gain can be achieved for very large TTrees with complicated expressions and cuts
    - maxentries parameter should be rather large
    Arguments:
    - tree       the tree
    - histo      the histogram
    - what       variable/expression/varlist to be projected
    - cuts       selection/weighting criteria 
    - nentries   number of entries to process  (>0: all entries in th tree)
    - first      the first entry to process
    - maxentries chunk size for parallel processing 
    """
    return parallel_poroject ( chain                   ,
                               target                  ,
                               expressions             ,
                               cuts       = cuts       ,
                               first      = first      ,
                               last       = last       ,
                               progress   = progress   , 
                               use_frame  = use_frame  ,
                               chunk_size = chunk_size ,
                               max_files  = max_files  , 
                               silent     = silent     , **kwargs )

ROOT.TTree .tproject = tproject
ROOT.TTree .pproject = tproject
ROOT.TChain.tproject = cproject
ROOT.TChain.pproject = cproject

# =============================================================================
_decorated_classes_ = (
    ROOT.TTree  ,
    ROOT.TChain ,    
    )

_new_methods_       = (
    ROOT.TTree .tproject ,
    ROOT.TTree .pproject ,
    ROOT.TChain.cproject ,
    ROOT.TChain.pproject ,     
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
#                                                                       The END 
# =============================================================================
