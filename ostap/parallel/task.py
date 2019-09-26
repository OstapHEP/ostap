#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/parallel/task.py
#
#  Useful utilities for multiprocessing and parallel processing for Ostap
#  Actualy it is just a little bit upgraded version of original
#  GaudiMP.Parallel module developed by Pere MATO Pere.Mato@cern.ch
#
#  The original module relies on some obsolete python modules (e.g. pyssh)
#  and it has limitations related to the pickling.
#
#  The upgraded module relies on <code>pathos</code> suite that
#  has very attratcive functionality and solve pickling issues
#  @see https://github.com/uqfoundation/pathos
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-02-23
# =============================================================================
"""Useful utilities for multiprocessing and parallel processing for Ostap

Actualy it is just a little bit upgraded version of original
GaudiMP.Parallel module developed by Pere MATO Pere.Mato@cern.ch

The original module relies on some obsolete python modules (e.g. pyssh)
and it has limitations related to the pickling.

The upgraded module relies on <code>pathos</code> suite that
has very attractive functionality and solve the issues with pickling
@see https://github.com/uqfoundation/pathos

"""
# =============================================================================
__version__ = '$Revision$'
__author__  = 'Vanya BELYAEV Ivan.Belyaev@itep.ru'
__date__    = '2016-02-23'
__all__     = (
    'Task'          , ## the base class for task
    'GenericTask'   , ## the generic ``templated'' task
    'Statistics'    , ## helper class to collect statistics 
    'StatMerger'    , ## helper class to merge   statistics
    'TaskMerger'    , ## simple merger for task results
    'task_executor' , ## helper function to execute Task  
    'func_executor' , ## helper function to execute callable 
    )
# =============================================================================
import operator, abc  
from   ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.paralllel.task' )
else                      : logger = getLogger ( __name__               ) 
# ==============================================================================
## @class Task
#  Basic base class to encapsulate any processing
#  that is going to be porcessed in parallel.
#  User class much inherit from it and implement the methods:
#  - <code>initializeLocal</code>
#  - <code>initializeRemote</code>
#  - <code>process</code>
#  - <code>finalize</code>
#  @author Pere MATO Pere.Meto@cern.ch
class Task ( object ) :
    """ Basic base class to encapsulate any processing that is
    going to be processed in parallel.
    User class must inherit from it and implement the methods
    - initializeLocal
    - initializeRemote
    - process
    - finalize.
    """
    __metaclass__ = abc.ABCMeta
    
    ## Local initialization:  invoked once on localhost for master task
    def initialize_local  ( self )          :
        """Local initialization:  invoked once on localhost for master task"""
        pass
    
    ## Remote initialization: invoked for each slave task on remote host
    def initialize_remote ( self )          :
        """Remote initialization: invoked for each slave task on remote host
        -  default: use the same as ``local initialization''
        """
        return self.initialize_local () 
    
    ## Finalization action: invoked once on local host for the master task (after merge)
    def finalize          ( self )          :
        """Finalization action: invoked once on local host for the master task (after merge)"""
        pass
    
    ## Process the (slave) task on remote host
    #  @attention must return the results
    @abc.abstractmethod 
    def process           ( self , item   ) :
        """Process the (slave) task on remote host
        - It must return the results!
        """
        return None 

    ## Collect and merge the resutls (invoked at local host)
    @abc.abstractmethod 
    def merge_results     ( self , result ) :
        """Collect and merge the results (invoked at local host)"""
        pass

    ## get the final merged task results 
    @abc.abstractmethod 
    def results           ( self )          :
        """Get the final merged task results"""
        return None 

    ## shortcut of <code>process</code> method 
    def __call__          ( self , item   ) :
        """Shortcut of process method"""
        return self.process ( item ) 


# =============================================================================
## @class GenericTask
#  Generic ``temlated'' task for Parallel processing  
#    One needs to  define three functions/functors:
#    - processor   :<code>        output = processor   ( item )               </code>
#    - merger      :<code>updated_output = merger ( old_output , new_output ) </code>
#    - initializer :<code>        output = initializer (      )               </code> 
class GenericTask(Task) :
    """Generic ``temlated'' task for parallel processing.
    One needs to  define three functions/functors:
    - processor   :         output = processor   ( item ) 
    - merger      : updated_output = merger ( old_output , new_output )
    - initializer :         output = initializer (      )  
    """
    # =========================================================================
    def __init__ ( self                ,
                   processor           ,
                   merger      = None  ,
                   initializer = tuple ) :
        """Generic task for parallel processing. One needs to  define three functions/functors
        - processor   :         output = processor   ( item ) 
        - merger      : updated_output = merger      ( old_output , new_output )
        - initializer :         output = initializer (      )  
        """

        if not merger :
            import operator
            merger = operator.iadd
            
        self.__processor   = processor
        self.__merger      = merger
        self.__initializer = initializer
        self.__output      = None 
        
    # =========================================================================
    ## local initialization (executed once in parent process)
    def initialize_local   ( self ) :
        """Local initialization (executed once in parent process)"""
        self.__output = self.initializer () if self.initializer else None 
        
    # =========================================================================
    ## the actual processing of the single item 
    def process  ( self , item ) :
        """The actual processing of the single item"""
        return self.processor ( item )
        
    # =========================================================================
    ## merge results 
    def merge_results ( self , result ) :
        """Merge processing results"""
        self.__output = self.merger ( self.__output , result )

    # =========================================================================
    ## get the final  results
    def results ( self ) :
        """Get the final(merged) results"""
        return self.__output

    # =========================================================================
    @property
    def processor  ( self ) :
        """``processor'' : the actual function for each subprocess
        - Signature: output = processor ( item ) 
        """
        return self.__processor
    @property
    def merger     ( self ) :
        """``merger'' : the actual fuction to merge results
        - Signature: updated_output = merger ( old_output , new_output )         
        """
        return self.__merger
    @property
    def initializer ( self ) :
        """``initializer'' : the actual fuction to initialize local output  
        - Signature: output = initializer() 
        """
        return self.__initializer

# =============================================================================
## @class Statistics
#  helper class to collect statistics 
#  @author Pere MATO Pere.Meto@cern.ch
class Statistics(object):
    """Helper class to collect statistics
    """
    def __init__ ( self , host = None ) :
        import time
        if not host :
            import socket 
            self.__host = socket.getfqdn ()
        else :
            self.__host = host 
        self.__start = time.time ( )
        self.time  = 0.0
        self.njobs = 0
        
    def stop ( self ) :
        import time
        self.time   = time.time () - self.__start
        self.njobs += 1

    def __enter__ ( self      ) : return self
    def __exit__  ( self , *_ ) : self.stop()
    
    @property 
    def host ( self  ) :
        """``host'' :  the host where the jobs executed"""
        return self.__host

    def __repr__  ( self ) :
        return "Statistics(%s,time=%.5g,njobs=%d)" % ( self.host , self.time , self.njobs )
    
    __str__ = __repr__

# =============================================================================
## Simple class to merge the job execution statistics 
class StatMerger(object) :
    """Simple class to merge the job execution statistics
    """
    def __init__ ( self        ) : self.__merged = {}        
    def __iadd__ ( self , stat ) :

        if isinstance  ( stat , StatMerger ) :
            for host , se in stat :
                if host in self.__merged :
                    ss  = self.__merged [ host ]
                    ss += se 
                else :
                    self.__merged [ host ] = se 
            return self 
        
        if not stat.host in self.__merged :
            self.__merged  [ stat.host ] = Statistics()
            
        se  = self.__merged [ stat.host ]
        se.time  += stat.time
        se.njobs += stat.njobs 
        
        return self

    def __len__ ( self ) : return len ( self.__merged )
    # =========================================================================
    ## iterator over merged statictics :
    #  @code
    #  merged = ...
    #  for host , stat in merged :
    #     ...
    #  @endcode 
    def __iter__  ( self ) :
        """Iterator over merged statictics :
        >>> merged = ...
        >>> for host , stat in merged :
        ...
        """
        for h in self.__merged :
            yield h, self.__merged [ h ]
            
    @property
    def merged ( self ) :
        """``merged'' : get the full merged statistic"""
        return self.__merged 

    ## Print the job execution statistics
    #  @code 
    #  merged = ...
    #  merged.print_stat ()
    #  @endcode 
    def print_stats ( self , prefix = '' , cputime = None ) :
        """Print job execution sstatistics
        >>> merged = ...
        >>> merged.print_stats () 
        """
        suffix = ''
        if cputime and 0 < cputime :

            sumtime = 0
            for host in self.__merged :
                se       = self.__merged[host]
                sumtime += se.time
                
            if 0 < sumtime :
                
                h1 , r1 = divmod ( cputime , 3600 )
                m1 , s1 = divmod ( r1      ,   60 )
                
                h2 , r2 = divmod ( sumtime , 3600 )
                m2 , s2 = divmod ( r2      ,   60 )
                
                h1 = int ( h1 ) 
                h2 = int ( h2 )
                
                m1 = int ( m1 ) 
                m2 = int ( m2 )
                
                s1 = int ( s1 )
                s2 = int ( s2 )
                
                if   h1     : suffix  = ' %02d:%02d:%02ds'    % ( h1 , m1 , s1 )
                elif m1     : suffix  = ' %02d:%02ds'         % (      m1 , s1 )
                else        : suffix  = ' %02ds'              %             s1 
                
                if   h2     : suffix += ' vs %02d:%02d:%02ds' % ( h2 , m2 , s2 )
                elif m2     : suffix += ' vs %02d:%02ds'      % (      m2 , s2 )
                else        : suffix += ' vs %02ds'           %             s2 
                
                gain = float ( sumtime ) / cputime
                suffix += ' (Gain: %.1f)'   % gain

                
        for line in self.__str__  ( prefix , suffix ).replace('\n#','\n').split('\n') : logger.info ( line ) 
                
    ## standard printout 
    def __str__  ( self , prefix = '' , suffix = '' ) :

        text = [] 
        text.append ( '%sJob execution statistics:%s'% ( prefix , suffix ) )
        text.append ( ' #jobs |   %   | total time |  time/job  | job server' )        
        fmt1 = '%6d | %5.1f | %10.4g | %10.4g | %-s'
        fmt0 = '%6d |       |            |            | %-s'
        njobs = self.njobs
        
        keys  = self.__merged.keys() 
        for host in sorted ( keys ) :
            se   = self.__merged [ host ]
            nj   = se.njobs
            time = se.time
            
            mean = time / nj if 1 <= nj else 0.0
            
            if 1 <= nj :
                text.append ( fmt1 % ( nj , 100. * nj / njobs , time , mean , host ) )
            else :
                text.append ( fmt0 % ( nj , host ) )
                                            
        return '\n# '.join  ( text ) 

    @property 
    def njobs ( self ) :
        """``njobs'' : total number of jobs"""
        return sum ( s.njobs for s in self.__merged.values() ) 
    
    __repr__ = __str__


# =============================================================================
## Merge task results 
#  @code
#  merger = TaskMerger()
#  jobs   = pool.uimap  ( .... )
#  for result , stat in jobs :
#      merger += result
#  merged = merger.result 
#  @encode 
class TaskMerger(object) :
    """Merge task resuls
    >>> merger = TaskMerger()
    >>> jobs   = pool.uimap  ( .... )
    >>> for result , stat in jobs :
    ...    merger += result 
    ... merged = merger.result
    """
    def __init__ ( self , merger = operator.iadd , init = None  ) :
        
        self.__merger = merger
        self.__result = init    
        
    # ========================================================================= 
    ## Merge task results
    #  @code
    #  merger = TaskMerger()
    #  jobs   = pool.uimap  ( .... )
    #  for result , stat in jobs :
    #      merger += result 
    #  merged = merger.result
    #  @encode 
    def __iadd__ ( self , result ) :
        """Merge task resuls
        >>> merger = TaskMerger()
        >>> jobs   = pool.uimap  ( .... )
        >>> for result , stat in jobs :
        ...    merger += result 
        ... merged = merger.result
        """
        self.merge  ( result ) 
        return self
    
    # ========================================================================= 
    ## Merge task resuls
    #  @code
    #  merger = TaskMerger()
    #  jobs   = pool.uimap  ( .... )
    #  for result , stat in jobs :
    #      merger.merge ( result )
    #  merged = merger.result
    #  @encode 
    def merge    ( self , result ) :
        """Merge task results
        >>> merger = TaskMerger()
        >>> jobs   = pool.uimap  ( .... )
        >>> for result , stat in jobs :
        ...    merger.merge ( result )
        ... merged = merger.result
        """
        
        if   self.__result is None : self.__result = result 
        elif self.__merger  :
            self.__result = self.__merger ( self.__result , result )
        elif hasattr ( self.__result , '__iadd__' ) :
            self.__result += result 
        elif hasattr ( self.__result , '__add__'  ) or hasattr ( result , '__radd__' ) :
            self.__result  = self.__result +  result
        else :
            raise TypeError ( 'TaskMerger: no merge is defined for %s and %s' % ( type ( self.__result ) , type ( result ) ) )
        
        return self
        
    @property
    def result  ( self ) :
        """``result'' : get the merged results"""
        return self.__result
    
# =============================================================================
## helper function to execute the task and collect stattistic
#  (unfornately due to limitation of <code>parallel python</code> one cannot
#  use decorators here :-(
def task_executor ( item ) :
    """Helper function to execute the task and collect job execution statistic
    - unfornately due to limitation of ``parallel python'' one cannot
    use python decorators here :-(
    """

    ## unpack
    task = item [ 0  ] 
    args = item [ 1: ] 

    ## perform remote  inialization (if needed) 
    task.initialize_remote() 
        
    with Statistics ()  as stat :    
        result = task.process ( *args ) 
        return result , stat

# =============================================================================
## helper function to execute the function and collect stattistic
#  (unfornately due to limitation of <code>parallel python</code> one cannot
#  use decorators here :-(
def func_executor ( item ) :
    """Helper function to execute the task and collect job execution statistic
    - unfornately due to limitation of ``parallel python'' one cannot
    use python decorators here :-(
    """

    ## unpack
    fun  = item [ 0  ] 
    args = item [ 1: ] 
    
    with Statistics ()  as stat :
        return fun ( *args ) , stat 
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )    
# =============================================================================
#                                                                       The END 
# =============================================================================

