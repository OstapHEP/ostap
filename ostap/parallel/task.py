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
#  - <code>initialize_local</code>
#  - <code>initialize_remote</code>
#  - <code>process</code>
#  - <code>finalize</code>
# 
#  One can specify following attributes:
#  - <code>directory</code>: the working directory for the job
#  - <code>environment</code>: additional environmental variables 
#  - <code>append_to</code>: append some path-like enviroment varibales 
#  - <code>prepend_to</code>: prepend some path-like enviroment varibales 
#  - <code>dot_in_path</code>: shoud the '.' be added to sysy.path?
#  @author Pere MATO Pere.Meto@cern.ch
class Task ( object ) :
    """ Basic base class to encapsulate any processing that is
    going to be processed in parallel.
    User class must inherit from it and implement the methods
    - initialize_local
    - initialize_remote
    - process
    - finalize.
    One can specify following attributes:
    - <directory : the working directory for the job
    - environment : additional environmental variables 
    - append_to : append some path-like enviroment varibales 
    - prepend_to : prepend some path-like enviroment varibales 
    - dot_in_path : shoud the '.' be added to sys.path?
    """
    __metaclass__ = abc.ABCMeta
    
    ## @attention ensure that the important attributes are available even before __init__
    def __new__( cls , *args , **kwargs):
        obj = super( Task , cls).__new__( cls )
        ## define the local trash 
        
        obj.__directory   = None
        obj.__environment = {}
        obj.__prepend_to  = {}
        obj.__append_to   = {}
        obj.__dot_in_path = None
        
        return obj

    ## Local initialization:  invoked once on localhost for master task
    def initialize_local  ( self )          :
        """Local initialization:  invoked once on localhost for master task"""
        pass
    
    ## Remote initialization: invoked for each slave task on remote host
    def initialize_remote ( self , jobid = -1 )  :
        """Remote initialization: invoked for each slave task on remote host
        - default: run ``local initialization''
        """
        return self.initialize_local () 
    
    ## Finalization action: invoked once on local host for the master task (after merge)
    def finalize          ( self )          :
        """Finalization action: invoked once on local host for the master task (after merge)"""
        pass
    
    ## Process the (slave) task on remote host
    #  @attention must return the result
    #  @param jobid jobid 
    @abc.abstractmethod 
    def process           ( self , jobid , *params ) :
        """Process the (slave) task on remote host
        - It must return the results!
        """
        return None 

    ## Collect and merge the results (invoked at local host)
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
    #  @param jobid jobid 
    #  @param item  item 
    def __call__          ( self , jobid , *params ) :
        """Shortcut of process method"""
        return self.process ( jobid , *params ) 

    @property
    def output  ( self ) :
        """``output'' : get a task output"""
        return self.results()

    @property
    def directory ( self ) :
        """``Directory'' : directory where job starts"""
        return self.__directory
    
    @directory.setter 
    def directory ( self , value ) :
        self.__directory = value 

    @property
    def environment ( self ) :
        """``environment'' : additional environment for the job"""
        return self.__environment
    
    @property
    def append_to ( self ) :
        """``append_to'' : a dictionary of enbvironemtn varibales to be appended"""
        return self.__append_to

    @property
    def prepend_to ( self ) :
        """``prepend_to'' : a dictionary of enbvironemtn varibales to be appended"""
        return self.__prepend_to
    
    @property
    def dot_in_path ( self ) :
        """``dot in path'' : has a dot in sys.path?"""
        return  self.__dot_in_path
    
    @dot_in_path.setter 
    def dot_in_path ( self , value ) :
        self.__dot_in_path = value
        
# =============================================================================
## @class GenericTask
#  Generic ``templated'' task for Parallel processing  
#  One needs to  define three functions/functors:
#    - processor   :<code>        output = processor ( jobid , item )         </code>
#    - merger      :<code>updated_output = merger ( old_output , new_output ) </code>
#    - initializer :<code>        output = initializer (      )               </code> 
#    - environment : additional environment for the job 
#    - append_to   : additional variables to be ''appended''
#    - prepend_to  : additional variables to be ''prepended''
class GenericTask(Task) :
    """Generic ``temlated'' task for parallel processing.
    One needs to  define three functions/functors:
    - processor   :         output = processor   ( jobid , item ) 
    - merger      : updated_output = merger ( old_output , new_output )
    - initializer :         output = initializer (      )  
    - initializer :         output = initializer (      )
    - directory   : change to this directory  (if it exists)
    - environment : additional environment for the job 
    - append_to   : additional variables to be ''appended''
    - prepend_to  : additional variables to be ''prepended''
    """
    # =========================================================================
    def __init__ ( self                ,
                   processor           ,
                   merger      = None  ,
                   initializer = tuple ,
                   directory   = None  ,
                   environment = {}    ,
                   append_to   = {}    ,
                   prepend_to  = {}    ) :
        """Generic task for parallel processing. One needs to define three functions/functors
        - processor   :         output = processor   ( jobid , item ) 
        - merger      : updated_output = merger      ( old_output , new_output )
        - initializer :         output = initializer (      )  
        - directory   : change to this directory  (if it exists) 
        - environment : additional environment for the job 
        - append_to   : additional variables to be ''appended''
        - prepend_to  : additional variables to be ''prepended''
        """

        if not merger :
            import operator
            merger = operator.iadd

        self.__processor   = processor
        self.__merger      = merger
        self.__initializer = initializer
        self.__output      = None
        
        self.directory     = directory
        self.environment  . update ( environment ) 
        self.append_to    . update ( append_to   ) 
        self.prepend_to   . update ( prepend_to  ) 
        
    # =========================================================================
    ## local initialization (executed once in parent process)
    def initialize_local   ( self ) :
        """Local initialization (executed once in parent process)"""
        self.__output = self.initializer () if self.initializer else None 
        
    # =========================================================================
    ## the actual processing of the single item 
    def process  ( self , jobid , *params ) :
        """The actual processing of the single item"""
        return self.processor ( jobid , *params )
        
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
        """``host'' :  the host where the job executed"""
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
                else        : suffix  = ' %2ds'               %             s1 
                
                if   h2     : suffix += ' vs %02d:%02d:%02ds' % ( h2 , m2 , s2 )
                elif m2     : suffix += ' vs %02d:%02ds'      % (      m2 , s2 )
                else        : suffix += ' vs %2ds'            %             s2 
                
                gain = float ( sumtime ) / cputime
                suffix += ' (Gain: %.1f)'   % gain

                
        title = prefix + 'Job execution statistics' + suffix
        logger.info ( title + "\n%s" % self.table ( title = title , prefix = "# " ) )
                
    ## standard printout as table 
    def table  ( self , title = 'Jobs execution statistics' , prefix = '' ) :

        text =  [ (' #jobs ' , '%' , ' total  time' , 'time/job' , 'job server') ]
        
        njobs = self.njobs        
        keys  = self.__merged.keys()
        
        for host in sorted ( keys ) :
            se   = self.__merged [ host ]
            nj   = se.njobs
            time = se.time
            
            mean = time / nj if 1 <= nj else 0.0
            
            if 1 <= nj :
                line = ( "%6d "     % nj                    ,
                         " %5.1f "  % ( 100. * nj / njobs ) ,
                         " %10.4g " % time ,
                         " %10.4g " % mean ,
                         " %-s"     % host )
            else :
                line = "%6d "% nj , '', '' , '' , " %-s" % host

            text.append ( line )
            
        import ostap.logger.table as T
        return T.table ( text , title = title , prefix = prefix )

    ## standard printout 
    def __str__  ( self ) :
        return self.table  ( title = "Job execution statistics" )

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
## helper function to execute the task and collect statistic
#  (unfornately due to limitation of <code>parallel python</code> one cannot
#  use decorators here :-(
#  @see Task 
def task_executor ( item ) :
    """Helper function to execute the task and collect job execution statistic
    - unfortunately due to limitation of ``parallel python'' one cannot
    use python decorators here :-(
    - see Task 
    """
    
    ## unpack
    task  = item [ 0  ]
    jobid = item [ 1  ] 
    args  = item [ 2: ] 

    import os, re, sys  

    what      =  r'(?<!\\)\$[A-Za-z_][A-Za-z0-9_]*' 
    expandvars = lambda item : re.sub ( what , '' , os.path.expandvars ( item ) )

    ## change the current working directory 
    if task.directory :
        directory = expandvars ( task.directory )
        if os.path.exists ( directory ) and os.path.isdir( directory ) :
            logger.debug ( 'Task %s: change the current directory to %s' %  ( jobid , directory ) )
            os.chdir ( directory ) 
            
    ## modify/update environment variables, if needed 
    for key in task.environment :
        item  = expandvars  ( task.environment [ key ] )
        logger.debug ( 'Task %s: modify the environment variable %s : %s ' % ( jobid , key , item ) )
        os.environ [ key ] =  item
        
    ## 2. prepend paths 
    for key in task.prepend_to   :
        item  = expandvars ( task.prepend_to [ key ] )
        ncmps = item.split ( os.pathsep )
        hask  = os.environ.get ( key , None )
        if hask is None : cmps =                             ncmps
        else            : cmps = hask.split ( os.pathsep ) + ncpms
        #
        path = os.pathsep.join ( cmps )
        os.environ [ key ] = path 
        logger.debug ( 'Task %s: prepend path %s : %s ' % ( jobid , key , path ) )

    ## 3. append paths 
    for key in task.append_to   :
        item  = expandvars ( task.append_to [ key ] )
        ncmps = item.split ( os.pathsep )
        hask  = os.environ.get ( key , None )
        if hask is None : cmps = ncmps
        else            : cmps = ncmps + hask.split ( os.pathsep )
        #
        path = os.pathsep.join ( cmps )
        os.environ [ key ] = path 
        logger.debug ( 'Task %s: append  path %s : %s ' % ( jobid , key , path ) )
        
    ## 4. Is current directory in the path? 
    if task.dot_in_path and not '.' in sys.path :
        sys.path  = ['.'] + sys.path
        logger.debug ( "Task %s: '.' is added to sys.path" % jobid )
        
    ## perform remote  inialization (if needed) 
    task.initialize_remote ( jobid ) 
        
    with Statistics ()  as stat :    
        result = task.process ( jobid , *args ) 
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
    fun   = item [ 0  ]
    jobid = item [ 1  ] 
    args  = item [ 2: ] 
    
    with Statistics ()  as stat :
        return fun ( jobid , *args ) , stat 
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )    
# =============================================================================
#                                                                       The END 
# =============================================================================

