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
import multiprocessing

from ostap.utils.progress_bar import ProgressBar
from ostap.logger.logger      import getLogger
from ostap.parallel.task      import Task, Statistics 
logger  = getLogger('ostap.parallel.mp_gaudi')

def _prefunction( f, task, item) :
    return f((task,item))
def _ppfunction( args ) :
    #--- Unpack arguments
    task, item = args
    stat = Statistics()
    #--- Initialize the remote side (at least once)
    if not task.__class__._initialized :
        for k,v in task.environ.items() :
            if k not in excluded_varnames : os.environ[k] = v
        task.initialize_remote()
        task.__class__._initialized = True
    #--- Reset the task output
    task.reset_output()
    #--- Call processing
    task.process(item)
    #--- Collect statistics
    stat.stop()
    return (copy.deepcopy(task.output), stat)


class WorkManager(object) :
    """ Class to in charge of managing the tasks and distributing them to
        the workers. They can be local (using other cores) or remote
        using other nodes in the local cluster """

    def __init__( self, ncpus='autodetect', ppservers=None , silent = False , **kwargs ) :
        
        if ncpus == 'autodetect' : self.ncpus = multiprocessing.cpu_count()
        else :                     self.ncpus = ncpus
        
        self.pool  = multiprocessing.Pool(self.ncpus)
        self.mode  = 'multicore'
        self.stats = {}
        
        self.silent = True if silent  else False 

    def __del__(self):
        if hasattr(self,'server') : self.server.destroy()

    def process(self, task, items, timeout=90000):
        if not isinstance(task,Task) :
            raise TypeError("task argument needs to be an 'Task' instance")
        # --- Call the Local initialialization
        task.initialize_local ()
        # --- Schedule all the jobs ....
        if self.mode == 'multicore' :
            start = time.time()
            jobs  = self.pool.map_async(_ppfunction, zip([task for i in items] , items ))
            
            with ProgressBar ( max_value = len ( items ) , description = "# Job execution:" ,  silent = self.silent ) as bar :              
                for result, stat in  jobs.get(timeout) :
                    task.merge_results    ( result )
                    self.merge_statistics ( stat   )
                    bar += 1
                    
            end = time.time()
            if not self.silent : 
                self.print_statistics()
                logger.info ( 'Time elapsed since server creation %f' % ( end - start ) ) 
        # --- Call the Local Finalize
        task.finalize()
    def print_statistics(self):
        njobs = 0
        for stat in self.stats.values():
            njobs += stat.njob
        logger.info ( 'Job execution statistics:' ) 
        logger.info ( 'job count | % of all jobs | job time sum | time per job | job server' ) 
        for name, stat  in self.stats.items():
            logger.info ( '   %6d |        %6.2f |   %10.3f |   %10.3f | %s' % (stat.njob, 100.*stat.njob/njobs, stat.time, stat.time/stat.njob, name) ) 

    def merge_statistics(self, stat):
        if stat.name not in self.stats : self.stats[stat.name] = Statistics()
        s = self.stats[stat.name]
        s.time += stat.time
        s.njob += 1

# == EOF ====================================================================================
