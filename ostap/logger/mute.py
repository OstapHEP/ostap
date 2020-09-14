#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with some simple but useful utilities
#   - suppression of stdout/stderr 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
#  
# =============================================================================
"""Module with some simple but useful utilities"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    #
    'tee_py'             , ## tee for Python's printouts
    'tee_cpp'            , ## tee for C++'s    printouts
    'output'             , ## redirect stdout/stderr into the file 
    'mute_py'            , ## suppress stdout/strerr Python printout 
    'silence_py'         , ## ditto 
    'mute'               , ## context manager to suppress stdout/strerr printout 
    'silence'            , ## ditto
    ##
    'TeeCpp'             , ## context manager (t ee   for C/C++  code) 
    'TeePy'              , ## context manager (tee    for python code) 
    'MuteC'              , ## context manager (mute   for C/C++  code) 
    'MutePy'             , ## context manager (mute   for python code) 
    'OutputC'            , ## context manager (output for C/C++  code) 
    )
# =============================================================================
import sys, os ## attention here!!
# =============================================================================
## @class MutePy
#  Very simple context manager to suppress python printout 
class MutePy(object):
    """A context manager for doing a ``deep suppression'' of stdout and stderr in 
    Python, i.e. will suppress all print, even if the print originates in a 
    compiled C/Fortran sub-function.
    This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).      
    
    stallen from  
    http://stackoverflow.com/questions/11130156/suppress-stdout-stderr-print-from-python-functions
    """
    def __init__( self , out = True , err = False ):
        self._out = out
        self._err = err
        
    def __enter__(self):
        #
        ## helper class to define empty stream 
        class Silent(object):
            def write(self,*args,**kwards) : pass

        self.stdout = sys.stdout
        self.stderr = sys.stderr
        
        if self._out : sys.stdout = Silent() 
        if self._err : sys.stderr = Silent() 

        return self
    
    def __exit__(self, *_):
        
        sys.stdout = self.stdout
        sys.stderr = self.stderr
        
# ============================================================================
## @class MuteC
#  context manager to suppress pythion prinout
#  the actual code is stallen from
#  http://stackoverflow.com/questions/11130156/suppress-stdout-stderr-print-from-python-functions
#  A fix is added for "IOError: [Errno 24] Too many open files" :
#  original code leaks the file descriptors
class MuteC(object):
    """A context manager for doing a ``deep suppression'' of stdout and stderr in 
    Python, i.e. will suppress all print, even if the print originates in a 
    compiled C/Fortran sub-function.
    This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).      
    
    stallen from  
    http://stackoverflow.com/questions/11130156/suppress-stdout-stderr-print-from-python-functions
    """
    #
    ## class variables: dev-null device & instance counter 
    _devnull = 0
    _cnt     = 0
    
    def __init__( self , out = True , err = False ):
        
        self._out = out
        self._err = err

        # increment instance counter 
        self.__class__._cnt += 1

        # create dev-null if not done yet 
        if not self.__class__._devnull :
            self.__class__._devnull = os.open ( os.devnull , os.O_WRONLY )            

    def __del__  ( self ) :
        
        # decrement instance counter 
        self.__class__._cnt -= 1
        
        # close dev-null if not done yet 
        if self.__class__._cnt <= 0 and self.__class__._devnull : 
            os.close ( self.__class__._devnull  )
            self.__class__._devnull = 0
            
    ## context-manager 
    def __enter__(self):
        
        ## Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds =  os.dup(1), os.dup(2)  # leak was here !!!
        
        ## mute it!
        if self._out : os.dup2 ( self.__class__._devnull , 1 )  ## C/C++
        if self._err : os.dup2 ( self.__class__._devnull , 2 )  ## C/C++

        return self
    
    ## context-manager 
    def __exit__(self, *_):
        
        # Re-assign the real stdout/stderr back to (1) and (2)  (C/C++)
        if self._err : os.dup2 ( self.save_fds[1] , 2 )
        if self._out : os.dup2 ( self.save_fds[0] , 1 )
        
        # fix the  file descriptor leak
        # (there were no such line in example, and it causes
        #      the sad:  "IOError: [Errno 24] Too many open files"
        
        os.close ( self.save_fds[1] ) 
        os.close ( self.save_fds[0] )

# =============================================================================
## dump all stdout/stderr information (including C/C++) into separate file
#  @code
#  with output ('output.txt') :
#           print 'ququ!'
#  @endcode 
#  @see MuteC 
class OutputC(object) :
    """Dump all stdout/stderr information into separate file:    
    >>>  with output ('output.txt') :
    ...             print 'ququ!'    
    """
    ## constructor: file name 
    def __init__ ( self , filename , out = True , err = False ) : 
        """Constructor
        """
        self._out  = out 
        self._err  = err
        self._file = open ( filename , 'w' ) 
            
    ## context-manager 
    def __enter__(self):
        
        if self._out : sys.stdout.flush()
        if self._err : sys.stderr.flush()
        
        self._file.flush()
        
        self._file.__enter__ () 
        ## Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds =  os.dup(1), os.dup(2)  # leak was here !!!
        
        ## mute it!
        if self._out : os.dup2 ( self._file.fileno() , 1 )  ## C/C++
        if self._err : os.dup2 ( self._file.fileno() , 2 )  ## C/C++

        return self
    
    ## context-manager 
    def __exit__( self , *_ ):

        if self._out : sys.stdout.flush()
        if self._err : sys.stderr.flush()
        
        self._file.flush()
        
        # Re-assign the real stdout/stderr back to (1) and (2)  (C/C++)
        if self._err : os.dup2 ( self.save_fds[1] , 2 )
        if self._out : os.dup2 ( self.save_fds[0] , 1 )
        
        # fix the  file descriptor leak
        # (there were no such line in example, and it causes
        #      the sad:  "IOError: [Errno 24] Too many open files"        
        os.close ( self.save_fds[1] ) 
        os.close ( self.save_fds[0] )
        
        self._file.__exit__ ( *_ )
        
        sys.stdout.flush()
        sys.stderr.flush()
        
# =============================================================================
## very simple context manager to duplicate Python-printout into file ("tee")
#  into separate file
#  @code
#  with tee('tee.txt') :
#           print 'ququ!'
#  @endcode
#  @attention: only Python printouts are grabbed 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  date 2012-07-06
class TeePy(object) :
    """Very simple context manager to duplicate Python-printout into file (``tee'')
    into separate file
    
    >>>  with tee('tee.txt') :
    ...        print 'ququ!'
    
    Unfortunately only Python printouts are grabbed 
    """
    ## constructor 
    def __init__( self , filename ):
        
        self._file = open ( filename , 'w' ) 

    ## context manager 
    def __enter__(self):
        
        self._file . __enter__ ()
        
        ## helper class to define empty stream 
        class _Tee(object):
            def __init__ ( self , the_file , the_stream ) :
                
                self._stream = the_stream 
                self._log    = the_file
                
            def write(self,*args) :
                
                self._stream .write ( *args ) 
                self._log    .write ( *args )
                
        self.stdout =  sys.stdout        
        sys.stdout  = _Tee ( self._file , self.stdout ) 

        return self
    
    ## context manager 
    def __exit__(self, *_):

        self._file.flush  ()
        self.stdout.flush ()
        
        sys.stdout = self.stdout
        
        self._file.__exit__ ( *_ )

# =============================================================================
## very simple context manager to duplicate C++-printout into file ("tee")
#  into separate file
#  @code
#  >>> with tee_cpp('tee.txt') :
#  ...         some_cpp_function() 
#  @endcode
#  @see Ostap::Utils::Tee
#  @attention: Python&C-printouts probably  are not affected 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  date 2012-07-07
class TeeCpp(object) :
    """Very simple context manager to duplicate C++-printout into file 
    into separate file
    
    >>> with tee_cpp('tee.txt') :
    ...         some_cpp_function()
    
    """
    def __init__ ( self , fname ) :
        sys.stdout.flush ()
        sys.stderr.flush ()
        from ostap.core.core import cpp 
        self.__tee = cpp.Ostap.Utils.Tee ( fname ) 
        
    ## context manager
    def __enter__ ( self      ) :
        sys.stdout.flush ()
        sys.stderr.flush ()
        self.__tee.enter ()
        return self
    
    ## context manager
    def __exit__  ( self , *_ ) :
        self.__tee.exit  ()
        del self.__tee
        sys.stdout.flush ()
        sys.stderr.flush ()


# =============================================================================
## very simple context manager to duplicate Python-printout into file ("tee")
#  into separate file
#  @code
#  >>> with tee_py ('tee.txt') :
#  ...         print 'ququ!'
#  @endcode
#  @attention: only Python prinouts are grabbed 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  date 2012-07-06
def tee_py ( filename ) :
    """Very simple context manager to duplicate Python-printout into file ("tee")
    into separate file
    >>> with tee('tee.txt') :
    ...        print 'ququ!'
    Unfortunately only Python printouts are grabbed 
    """
    return TeePy ( filename ) 
    
# =============================================================================
## very simple context manager to duplicate C++-printout into file ('tee')
#  into separate file
#  @code
#  >>> with tee_cpp ('tee.txt') : some_cpp_code() 
#  @endcode
#  @attention: only C/C++ printouts are grabbed 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  date 2012-07-06
def tee_cpp ( filename ) :
    """Very simple context manager to duplicate C++-printout into file ('tee')
    into separate file
    >>> with tee_cpp('tee.txt') : some_cpp_code()
    Unfortunately only C/C++ printouts are grabbed 
    """
    return TeeCpp ( filename ) 


# =============================================================================
## simple context manager to redirect all (C/C++/Python) printout 
#  into separate file
#  @code
#  >>> with output ('output.txt') :
#  ...         print 'ququ!'
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  date 2012-07-06
def output ( fname , cout = True , cerr = False ) :
    """ Simple context manager to redirect all (C/C++/Python) printotu
    
    >>> with output ('output.txt') :
    ...               print 'ququ!'
    
    """
    return OutputC ( fname  , cout , cerr )

# =============================================================================
## simple context manager to suppress C/C++-printout
#
#  @code
#  >>> with mute () :
#  ...        <some code here>
#  @endcode 
def mute ( cout = True , cerr = False )   :
    """Simple context manager to suppress C/C++ printout
    
    >>> with mute () :
    ...     <some code here>
    """
    return MuteC ( cout , cerr )

# =============================================================================
## simple context manager to suppress Python-printout
#
#  @code
#  >>> with mute_py () :
#  ...        <some code here>
#  @endcode 
def mute_py ( cout = True , cerr = False )   :
    """Simple context manager to suppress python printouts
    
    >>> with mute_py () :
    ...    <some code here>    
    """
    return MutePy ( cout , cerr )

# ==============================================================================
## ditto 
silence_py  = mute_py  # ditto
silence     = mute     # ditto


# =============================================================================
if '__main__' == __name__ :
    
    from   ostap.logger.logger    import getLogger
    logger = getLogger ('ostap.logger.mute')
    
    from ostap import banner
    logger.info ( __file__  + '\n' + banner )
    logger.info ( 80*'*'   )
    logger.info ( __doc__  )
    logger.info ( 80*'*' )
    logger.info ( ' Author  : %s' %         __author__    ) 
    logger.info ( ' Version : %s' %         __version__   ) 
    logger.info ( ' Date    : %s' %         __date__      )
    logger.info ( ' Symbols : %s' %  list ( __all__     ) )
    logger.info ( 80*'*' ) 
    
# =============================================================================
##                                                                     The END 
# =============================================================================
