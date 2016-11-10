# File: GaudiMP/Parallel.py
# Author: Pere Mato (pere.mato@cern.ch)

""" GaudiMP.Parallel module.
    This module provides 'parallel' processing support for GaudiPyhton.
    It is adding some sugar on top of public domain packages such as
    the 'multiprocessing' or the 'pp' packages. The interface can be made
    independent of the underlying implementation package.
    Two main class are defined: Task and WorkManager
"""

__all__ = [ 'Task','WorkManager' ]
excluded_varnames = ['HOSTNAME', 'SSH_CLIENT', 'SSH_CONNECTION', 'DISPLAY']

import sys, os, time, copy
import multiprocessing

def _prefunction( f, task, item) :
    return f((task,item))
def _ppfunction( args ) :
    #--- Unpack arguments
    task, item = args
    stat = Statistics()
    #--- Initialize the remote side (at least once)
    if not task.__class__._initializeDone :
        for k,v in task.environ.items() :
            if k not in excluded_varnames : os.environ[k] = v
        task.initializeRemote()
        task.__class__._initializeDone = True
    #--- Reset the task output
    task._resetOutput()
    #--- Call processing
    task.process(item)
    #--- Collect statistics
    stat.stop()
    return (copy.deepcopy(task.output), stat)

class Statistics(object):
    def __init__(self):
        self.name  = os.getenv('HOSTNAME')
        self.start = time.time()
        self.time  = 0.0
        self.njob  = 0
    def stop(self):
        self.time = time.time() - self.start

class Task(object) :
    """ Basic base class to encapsulate any processing that is going to be porcessed in parallel.
        User class much inherit from it and implement the methods initializeLocal,
        initializeRemote, process and finalize.   """
    _initializeDone = False
    def __new__ ( cls, *args, **kwargs ):
        task = object.__new__( cls )
        task.output = ()
        task.environ = {}
        for k,v in os.environ.items(): task.environ[k] = v
        task.cwd = os.getcwd()
        return task
    def initializeLocal(self):
        pass
    def initializeRemote(self):
        pass
    def process(self, item):
        pass
    def finalize(self) :
        pass
    def _mergeResults(self, result) :
        if type(result) is not type(self.output) :
            raise TypeError("output type is not same as obtained result")
        #--No iteratable---
        if not hasattr( result , '__iter__' ):
            if hasattr(self.output,'Add') : self.output.Add(result)
            elif hasattr(self.output,'__iadd__') : self.output += result
            elif hasattr(self.output,'__add__') : self.output = self.output + result
            else : raise TypeError('result cannot be added')
        #--Dictionary---
        elif type(result) is dict :
            if self.output.keys() <= result.keys(): minkeys = self.output.keys()
            else: minkeys = result.keys()
            for key in result.keys() :
                if key in self.output :
                    if hasattr(self.output[key],'Add') : self.output[key].Add(result[key])
                    elif hasattr(self.output[key],'__iadd__') : self.output[key] += result[key]
                    elif hasattr(self.output[key],'__add__') : self.output[key] = self.output[key] + result[key]
                    else : raise TypeError('result cannot be added')
                else :
                    self.output[key] = result[key]
        #--Anything else (list)
        else :
            for i in range( min( len(self.output) , len(result)) ):
                if hasattr(self.output[i],'Add') : self.output[i].Add(result[i])
                elif hasattr(self.output[i],'__iadd__') : self.output[i] += result[i]
                elif hasattr(self.output[i],'__add__') : self.output[i] = self.output[i] + result[i]
                else : raise TypeError('result cannot be added')
    def _resetOutput(self):
        output =  (type(self.output) is dict) and self.output.values() or self.output
        for o in output :
            if hasattr(o, 'Reset'): o.Reset()


class WorkManager(object) :
    """ Class to in charge of managing the tasks and distributing them to
        the workers. They can be local (using other cores) or remote
        using other nodes in the local cluster """

    def __init__( self, ncpus='autodetect', ppservers=None) :
        if ncpus == 'autodetect' : self.ncpus = multiprocessing.cpu_count()
        else :                     self.ncpus = ncpus
        if ppservers :
            import pp
            self.ppservers = ppservers
            self.sessions = [ SshSession(srv) for srv in ppservers ]
            self.server = pp.Server(ncpus=self.ncpus, ppservers=self.ppservers)
            self.mode = 'cluster'
        else :
            self.pool = multiprocessing.Pool(self.ncpus)
            self.mode = 'multicore'
        self.stats = {}

    def __del__(self):
        if hasattr(self,'server') : self.server.destroy()

    def process(self, task, items, timeout=90000):
        if not isinstance(task,Task) :
            raise TypeError("task argument needs to be an 'Task' instance")
        # --- Call the Local initialialization
        task.initializeLocal()
        # --- Schedule all the jobs ....
        if self.mode == 'cluster' :
            jobs = [self.server.submit(_prefunction, (_ppfunction, task, item), (), ('GaudiMP.Parallel','time')) for item in items]
            for job in jobs :
                result, stat = job()
                task._mergeResults(result)
                self._mergeStatistics(stat)
            self._printStatistics()
            self.server.print_stats()
        elif self.mode == 'multicore' :
            start = time.time()
            jobs = self.pool.map_async(_ppfunction, zip([task for i in items] , items ))
            for result, stat in  jobs.get(timeout) :
                task._mergeResults(result)
                self._mergeStatistics(stat)
            end = time.time()
            self._printStatistics()
            print 'Time elapsed since server creation %f' %(end-start)
        # --- Call the Local Finalize
        task.finalize()
    def _printStatistics(self):
        njobs = 0
        for stat in self.stats.values():
            njobs += stat.njob
        print 'Job execution statistics:'
        print 'job count | % of all jobs | job time sum | time per job | job server'
        for name, stat  in self.stats.items():
            print '       %d |        %6.2f |     %8.3f |    %8.3f | %s' % (stat.njob, 100.*stat.njob/njobs, stat.time, stat.time/stat.njob, name)

    def _mergeStatistics(self, stat):
        if stat.name not in self.stats : self.stats[stat.name] = Statistics()
        s = self.stats[stat.name]
        s.time += stat.time
        s.njob += 1


class SshSession(object) :
    def __init__(self, hostname):
        import pyssh
        import pp
        self.host = hostname
        ppprefix =  os.path.dirname(os.path.dirname(pp.__file__))
        self.session = pyssh.Ssh(host=hostname)
        self.session.open()
        self.session.read_lazy()
        self.session.write('cd %s\n' % os.getcwd())
        self.session.read_lazy()
        self.session.write('setenv PYTHONPATH %s\n' % os.environ['PYTHONPATH'])
        self.session.read_lazy()
        self.session.write('setenv LD_LIBRARY_PATH %s\n' % os.environ['LD_LIBRARY_PATH'])
        self.session.read_lazy()
        self.session.write('setenv ROOTSYS %s\n' % os.environ['ROOTSYS'])
        self.session.read_lazy()
        self.session.write('%s %s/scripts-%s/ppserver.py \n'%(sys.executable, ppprefix, sys.version.split()[0] ))
        self.session.read_lazy()
        self.session.read_lazy()
        print 'started ppserver in ', hostname
    def __del__(self):
        self.session.close()
        print 'killed ppserver in ', self.host

# == EOF ====================================================================================
