#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==========================================================================================
## @file ostap/utils/parallel_copy.py
#  Copy files in parallel using GNU paralllel or xargs (if/when available)
#  @see https://www.gnu.org/software/parallel/
#  @date   2024-04-18
#  @author Vanya  BELYAEV Ivan.Belyaev@itep.ru
# =============================================================================
"""  Copy files in parallel using GNU paralllel or xargs (if/when available) 
- see https://www.gnu.org/software/parallel/
"""
# =============================================================================
__author__  = 'Vanya BELYAEV  Ivan.Belyaev@itep.ru'
__date__    = "2020-01-18"
__version__ = '$Revision$'
__all__     = (
    'copy_files' , ## copy files in parallel usins GNU parallel
)
# =============================================================================
from   ostap.utils.utils import which
import os 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.utils.parallel.copy' )
else                       : logger = getLogger ( __name__                    )
# =============================================================================

# ============================================================================
## Copy files in parallel using GNU paralell or xargs 
#  @see https://www.gnu.org/software/parallel/
#  - switch to multiprocess-based parallelisation when GNU parallel is not avialable 
def copy_files ( file_pairs      ,
                 progress = True ,
                 copier   = None ,
                 copy_cmd = ''   , **kwargs ) :
    """Copy files in parallel using GNU parallel or xargs 
    - see https://www.gnu.org/software/parallel/
    - switch to multiprocess-based parallelisation when GNU parallel is not available 
    """
    if not copier :
        from ostap.utils.basic import copy_file
        copier = copy_file
        
    if not copy_cmd : copy_cmd = 'cp'

    ## keep arguments 
    pairs = tuple ( p for p in file_pairs )

    # ======================================================================
    ## sequential processing 
    # ======================================================================
    nfiles = len ( pairs ) 
    from ostap.utils.basic import numcpu 
    if nfiles <= 1 or numcpu () <=1  :
        from ostap.utils.progress_bar import progress_bar
        silent = nfiles <= 1 or not progress 
        copied = [] 
        for f, nf in progress_bar ( pairs , silent = silent ) :
            output = copier ( f , nf , progress = progress and nfiles <=1 )
            result = f, output 
            copied.append ( result ) 
        return tuple ( copied )

    parallel = which ( 'parallel' )
    xargs    = '' if parallel else which ( 'xargs' )

    # =====================================================================
    ## Neither GNU parallel nor xargs is available: 
    if not parallel and not xargs :
        from ostap.parallel.parallel_copy import copy_files as cpfiles
        return cpfiles ( file_pairs          ,
                         progress = progress ,
                         copier   = copier   , **kwargs ) 

    # =====================================================================
    ## (1) collect all directories
    ddirs = set() 
    for f, nf in pairs :
        ddir = os.path.dirname ( nf )
        if not os.path.exists ( ddir ) : 
            ddirs.add ( ddir )

    ## (2) recreate all required directories
    from ostap.utils.basic import make_dirs 
    while ddirs : make_dirs ( ddirs.pop() , exist_ok = True ) 

    if parallel :
        command = 'parallel --bar :::' if progress else 'parallel :::'
        thecmd  = '%s %%s %%s\n' % copy_cmd
    else        :
        command = 'xargs -P%d -L1 %s ' % ( numcpu() , copy_cmd ) 
        if progress and 'cp' == copy_cmd : command += ' -v ' 
        thecmd  = '%s %s\n' 

    ## (3) prepare the input file with the commands
    import ostap.utils.cleanup as CU
    tmpfile = CU.CleanUp.tempfile( suffix = '.lst' )  
    with open ( tmpfile , 'w' ) as cmd :
        for f, nf in pairs :
            cmd.write ( thecmd % ( f , nf ) )

    ## (4) finally, call GNU parallel/xargs command 
    import subprocess, shlex 
    with open ( tmpfile , 'r' ) as input :
        subprocess.check_call ( shlex.split ( command ) , stdin = input )

    ## (5) check the final results 
    results = set() 
    for f , nf in pairs :
        if os.path.exists ( nf ) and os.path.isfile ( nf ) :
            result = f , os.path.abspath ( nf ) 
            results.add ( result )
        else :

            logger.warning ( "copy_files: no expected output '%s'" % nf ) 

    return tuple ( results ) 


# ============================================================================
## Sync/copy files in parallel using GNU paralell or xargs 
#  @see https://www.gnu.org/software/parallel/
#  - switch to multiprocess-based parallelisation when GNU parallel is not avialable 
def sync_files ( file_pairs            ,
                 progress = True       ,
                 copier   = None       ,
                 copy_cmd = ''         , **kwargs ) :
    """Sync/copy files in parallel using GNU parallel or xargs 
    - see https://www.gnu.org/software/parallel/
    - switch to multiprocess-based parallelisation when GNU parallel is not available 
    """
    if not which ( 'rsync' ) :
        return copy_files ( file_pairs          ,
                            progress = progress ,
                            copier   = None     ,
                            copy_cmd = ''       , **kwargs )
    
    if not copy_cmd :
        copy_cmd = 'rsync -a'
    
    if not copier :
        from ostap.utils.basic import sync_file
        copier =  sync_file
        
    return copy_files ( file_pairs          ,
                        progress = progress ,
                        copier   = copier   ,
                        copy_cmd = copy_cmd , **kwargs )

# =============================================================================
if not which ( 'parallel' ) and not which ( 'xargs' ) :
    from ostap.parallel.parallel_copy import copy_files

# =============================================================================
if '__main__' == __name__ :

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END 
# =============================================================================
