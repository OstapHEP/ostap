# File: GaudiMP/Parallel.py
# Author: Pere Mato (pere.mato@cern.ch)

""" GaudiMP.Parallel module.
    This module provides 'parallel' processing support for GaudiPyhton.
    It is adding some sugar on top of public domain packages such as
    the 'multiprocessing' or the 'pp' packages. The interface can be made
    independent of the underlying implementation package.
    Two main class are defined: Task and WorkManager
"""
from   __future__        import print_function
__all__ = [ 'Task','WorkManager' ]
excluded_varnames = ['HOSTNAME', 'SSH_CLIENT', 'SSH_CONNECTION', 'DISPLAY']

import sys, os, time, copy


from ostap.utils.progress_bar  import progress_bar
from ostap.logger.logger       import getLogger
from ostap.parallel.task       import Task, Statistics ,  StatMerger 
logger  = getLogger('ostap.parallel.mp_gaudi')

# =============================================================================
vi = sys.version_info
if 3 <= vi.major and 6 <= vi.minor :
    import multiprocessing     as MP 
else : 
    try: 
        import multiprocess    as MP
        import dill
    except ImportError :
        import multiprocessing as MP 
# =============================================================================

        
def _prefunction( f, task , jobid , item) :
    return f( ( task , jobid , item ) )
def _ppfunction ( args ) :
    #--- Unpack arguments
    task, jobid , item = args
    with Statistics() as stat : 
        task.initialize_remote ( jobid )
        result = task.process  ( jobid , item )
        stat.stop()
        return result , stat

class WorkManager(object) :
    """ Class to in charge of managing the tasks and distributing them to
        the workers. They can be local (using other cores) or remote
        using other nodes in the local cluster """

    def __init__( self, ncpus='autodetect', ppservers=None , silent = False , **kwargs ) :
        
        if ncpus == 'autodetect' : self.ncpus = MP.cpu_count()
        else :                     self.ncpus = ncpus
        
        self.pool  = MP.Pool(self.ncpus)
        self.stats = StatMerger()
        
        self.silent = True if silent  else False 

    def __del__(self):
        if hasattr(self,'server') : self.server.destroy()

    def process(self, task, items, timeout=90000):
        if not isinstance(task,Task) :
            raise TypeError("task argument needs to be an 'Task' instance")
        # --- Call the Local initialialization
        task.initialize_local ()
        # --- Schedule all the jobs ....
            
        start = time.time()
        from itertools import repeat , count        
        
        jobs  = self.pool.imap_unordered ( _ppfunction, zip ( repeat ( task ) , count () , items ) )
        for result , stat in progress_bar ( jobs , max_value = len ( items ) , description = "# Job execution:" ,  silent = self.silent ) :
            task.merge_results    ( result )
            self.stats += stat   

        ## jobs  = self.pool.map_async ( _ppfunction, zip ( repeat ( task ) , count () , items ) )
        ## with ProgressBar ( max_value = len ( items ) , description = "# Job execution:" ,  silent = self.silent ) as bar :              
        ##    for result, stat in  jobs.get(timeout) :
        ##        task.merge_results    ( result )
        ##        self.stats += stat   
        ##        bar += 1
                    
        end = time.time()
        if not self.silent : 
            self.print_statistics()
            logger.info ( 'Time elapsed since server creation %f' % ( end - start ) ) 
        # --- Call the Local Finalize
        task.finalize()
        return task.results()
    
    def print_statistics(self):
        self.stats.print_stats ()


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )    
    logger.info ("Module ``%s'' is used for multiprocessing" % MP.__name__ )
    
    
# == EOF ====================================================================================
