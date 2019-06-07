#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/utils/cleanup.py
#  Module with some simple but useful utilities
#  to deal  with tempoary files and directories 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10  
# =============================================================================
"""Module with some simple but useful utilities
to deal  with tempoary files and directories 
Module with some simple but useful utilities for
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
__all__     = (
    'CleanUp'  , ## Simple (base) class to deal with temporary files and directories
    'TempFile' , ## Simple placeholder for temproary file  
    )
# =============================================================================
import os
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils.cleanup' )
else                       : logger = getLogger( __name__              )
del getLogger
# =============================================================================
## @class CleanUp
#  Simple (base) class to get temporary files
#  and directories and to remove them at exit
class  CleanUp(object) :
    """Simple (base) class to get temporary files
    and directories and to remove them at exit
    """
    _tmpfiles = set()
    _tmpdirs  = set()

    @property
    def tmpdir ( self ) :
        """``tmpdir'' : return the name of temporary managed directory
        - the managed directory will be cleaned-up and deleted at-exit
        >>> o    = CleanUp()
        >>> tdir = o.tmpdir 
        """
        tdir = CleanUp.tempdir () 
        return tdir
    
    @property
    def tmpdirs  ( self ) :
        """``tmpdirs'' - list of currently registered managed temporary directories"""
        return tuple ( self._tmpdirs )
    
    @tmpdirs.setter
    def tmpdirs ( self, values ) :
        if instance ( values , str ) : values = [ values ]
        for o in values :
            if o and isinstance ( o ,  str ) :
                self._tmpdirs.add ( o )
                logger.debug ( 'temporary directory     added %s' % o )
                
    @property 
    def tmpfiles ( self ) :
        """``tempfiles'' : list of registered managed temporary files"""
        return list ( self._tmpfiles )
    
    @tmpfiles.setter
    def tmpfiles ( self , other ) :
        if isinstance ( other , str ) : other = [ other ]
        for o in other :
            if o and isinstance ( o , str ) : 
                self._tmpfiles.add ( o )
                logger.debug ( 'temporary file          added %s' % o )

    @staticmethod
    def tempdir ( suffix = '' , prefix = 'tmp_' ) :
        """Get the name of the temporary directory.
        The directory will be cleaned-up and deleted at-exit.
        >>> dirname = CleanUp.tempdir() 
        """
        import tempfile
        tmp = tempfile.mkdtemp ( suffix = suffix , prefix = prefix ) 
        CleanUp._tmpdirs.add ( tmp )
        logger.debug ( 'temporary directory requested %s' % tmp   )
        return tmp        
    
    @staticmethod
    def tempfile ( suffix = '' , prefix = 'tmp_' , dir = None ) :
        """Get the name of the temporary file. The file will be deleted at-exit
        >>> fname = CleanUp.tempfile() 
        """
        import tempfile
        _file = tempfile.NamedTemporaryFile( suffix = suffix ,
                                             prefix = prefix ,
                                             dir    = dir    , 
                                             delete = False  )
        fname = _file.name
        _file.close()
        os.unlink(fname)
        assert not os.path.exists ( fname )
        CleanUp._tmpfiles.add ( fname )
        logger.debug ( 'temporary file      requested %s' % fname )
        return fname

    @staticmethod
    def remove_file ( fname ) :
        """Remove the (temporary) file
        """
        if os.path.exists ( fname ) and os.path.isfile ( fname ) :
            logger.verbose ( 'remove temporary file : %s' % fname )
            try    : os.remove ( fname )
            except : pass
        if os.path.exists  ( fname ) and os.path.isfile ( fname ) :
            logger.error   ( 'failed to remove file : %s' %  fname  )
            return False 
        return True 

    @staticmethod
    def remove_dir ( fdir ) :
        """Remove the (temporary) directory
        """
        if os.path.exists ( fdir ) and os.path.isdir ( fdir ) :
            logger.verbose ( 'remove temporary dir : %s' % fdir  )
            ## 1: collect all files & subdirectories 
            for root, subdirs, files in os.walk ( fdir  , topdown = False ):
                ## 2: remove all files 
                for ff in files : CleanUp.remove_file ( os.path.join ( root , ff  ) ) 
                ## 3: remove all directories 
                for dd in subdirs :
                    dd = os.path.join ( root , dd )
                    logger.verbose ( 'remove subdirectory %s in temporary directory %s ' % ( dd , fdir ) )
                    try    : os.rmdir   ( dd  )
                    except : pass
                    if os.path.exists ( dd ) and os.path.isdir ( dd )   :
                        logger.error ( 'failed to remove %s in temporary directory %s ' % ( dd , fdir ) )                        
            ## 4: finally remove the root
            try    : os.rmdir ( fdir )
            except : pass 
        if os.path.exists ( fdir ) and os.path.isdir ( fdir ) :
            logger.error ( 'failed to  remove : %s' % fdir  )
            
    @staticmethod
    def remove ( fname ) :
        """Remove temporary object (if any) 
        """
        if   os.path.exists ( fname ) and os.path.isdir  ( fname ) :
            return CleanUp.remove_dir  ( fname )
        elif os.path.exists ( fname ) and os.path.isfile ( fname ) :
            return CleanUp.remove_file ( fname )
        
# =============================================================================
import atexit
@atexit.register
def _cleanup_ () :

    ## 1. clean up the files 
    tmp_files  = CleanUp._tmpfiles
    logger.debug ( 'remove temporary files: %s' % list ( tmp_files ) )
    while tmp_files :
        f = tmp_files.pop()
        CleanUp.remove_file ( f )

    ## 2. clean up  the directories 
    tmp_dirs = CleanUp._tmpdirs
    logger.debug ( 'remove temporary directories: %s' % list ( tmp_dirs ) )
    while tmp_dirs :
        f = tmp_dirs.pop()
        CleanUp.remove_dir ( f )

# =============================================================================
## @class TempFile
#  Base class-placeholder for the temporary  file 
class TempFile(object) :
    """Base class-placeholder for the temporary  file 
    """
    
    def __init__ ( self , suffix = '' , prefix = '' , dir = None ) :
        self.__filename = CleanUp.tempfile (  suffix = suffix ,
                                              prefix = prefix ,
                                              dir    = dir    )
        
    @property
    def filename ( self ) :
        """``filename'': the actual name of temporary file"""
        return self.__filename

    ## context  manager enter: no action 
    def __enter__ ( self      ) :
        """Context  manager enter: no action"""
        return self

    ## Context manager  exit: delete the file 
    def __exit__  ( self , *_ ) :
        """Context manager  exit: delete the file
        """
        if self.__filename and os.path.exists ( self.__filename ) :
            CleanUp.remove_file ( self.__filename )

    ## delete the file 
    def __del__  ( self ) :
        """Delete the temporary file
        """        
        if self.__filename and os.path.exists ( self.__filename ) :
            CleanUp.remove_file ( self.__filename )
            
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    logger.info ( 80*'*' ) 
            
# =============================================================================
# The END 
# =============================================================================
