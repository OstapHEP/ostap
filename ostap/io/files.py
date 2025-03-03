#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file files.py 
#  Module with useful tool to handle set of files
#  Usually it is used as a base class for other tools 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
# =============================================================================
""" Module with useful tool to handle set of files
Usually it is used as a base class for other tools 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    'Files'       , ## base class tool to handle set of files
    'copy_files'  , ## a bit specific copy of set of files
    'sync_files'  , ## synchonize the files
    'sync_dirs'   , ## synchonize the directories     
)
# =============================================================================
from   ostap.core.ostap_types import integer_types, path_types, sized_types  
from   ostap.parallel.task    import Task
import os, glob, random, math   
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.io.files' )
else                       : logger = getLogger( __name__         )
del getLogger
# =============================================================================
## @class TreatFileTask
#  Helper class for parallel processing of long list of files 
#  It defiend the basic teatemtn of the, calling `Files.treatFile` method
class TreatFileTask(Task) :
    """ Helper class for parallel processing of long list of files 
    It defiens the basic teatemtn of the, calling `Files.treatFile` method
    """
    def __init__ ( self , data ) :
        self.__data  = data.clone ( files = [] ) 
    def initialize_local  ( self              ) : pass 
    def initialize_remote ( self , jobid = -1 ) : pass
    ## actgual processing  
    def process ( self , jobid , files ) :
        """ Actual processing"""
        for f in files :  self.__data.treatFile ( f ) 
        return self.__data 
    ## get the results 
    def results (  self ) : return self.__data 
    ## merge the results 
    def merge_results  ( self , results , jobid = -1 ) :    
        self.__data += results

# ===========================================================================
## get the size with units for the file-size units
def fsize_unit ( size ) :
    """ Get the size with units for the file-size units"""
    assert isinstance ( size , integer_types ) , "unit: invalid 'size' type"
    if   ( 1024 ** 5 ) <= size : return  size // ( 1024 ** 5 ) , 'PB'
    elif ( 1024 ** 4 ) <= size : return  size // ( 1024 ** 4 ) , 'TB'
    elif ( 1024 ** 3 ) <= size : return  size // ( 1024 ** 3 ) , 'GB'
    elif ( 1024 ** 2 ) <= size : return  size // ( 1024 ** 2 ) , 'MB'
    elif ( 1024      ) <= size : return  size // ( 1024      ) , 'kB'
    return size , 'B'

# =============================================================================
## @class Files
#  Simple utility to pickup the list of files
#  - collect file the file
#  - keeps and store them
#  - union, difference, intersection, subtraction, etc...
#  - copy and move collection of files
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @author Alexander BARANOV a.baranov@cern.ch
#  @date   2014-06-08  
class Files(object):
    """ Simple utility to pickup the list of files     

    - collect file the file
    - keeps and store them
    - union, difference, intersection, subtraction, etc...
    - copy and move collection of files

    >>> data  = Files( '*.root' )
    >>> files = data.files
    """
    def __init__( self                ,
                  files               ,
                  description = ""    ,
                  maxfiles    = -1    ,
                  silent      = False ,
                  sorted      = False ,    ## sort the list of files  ? 
                  parallel    = False ) :  ## collect them in parallel?
        #
        if   isinstance ( files , str   ) : files = [ files ]
        elif isinstance ( files , Files ) : files = files.files   
        #
        
        self.__description  = description
        self.__silent       = silent 

        pats = list ( set ( f for f in files ) )
        pats.sort() 
        self.__patterns     = tuple ( pats ) 
                
        self.__maxfiles     = maxfiles
        self.__sorted       = True if sorted else False 
        self.__files        = []        
        self.__bad_files    = set()
 
        assert isinstance ( maxfiles , int ) , "Invalid type for 'maxfiles'!"
        
        #  =====================================================================
        # convert list of patterns into the list of files 
        # =====================================================================
        
        _files   = self.the_files ()

        if not self.silent :
            logger.info ('Loading: %s  #patterns/files: %s/%d' % ( self.description   ,
                                                                   len(self.patterns) , 
                                                                   len( _files )    ) )            
        ## can use parallel processing here
        nfiles = len ( _files )

        chunk_size = min ( 20 , max ( nfiles // 3 , 2 ) )         
        if parallel and chunk_size < nfiles and ( maxfiles < 0 or nfiles <= maxfiles ) :
            
            jobs = []
            from   ostap.utils.utils import chunked 
            for chunk in chunked ( _files , chunk_size ) : jobs.append ( chunk )
                
            psilent = self.silent
            self.silent = True 
            
            task = TreatFileTask ( self )
            from ostap.parallel.parallel import WorkManager 
            wmgr = WorkManager ( silent = True , progress = not self.silent )

            wmgr.process ( task , jobs )

            results      = task.results()            
            self        += results             
            self.silent  = psilent
            
        else : 

            self.__add_files ( _files , maxfiles )
        
        if not self.silent : logger.info ('Loaded: %s' % self )
            
    @property 
    def files     ( self ) :
        """'files' : the list(tuple) of files"""
        return tuple ( self.__files )
    
    @property
    def description ( self ) :
        """'description': description of this collection"""
        return self.__description
    @description.setter
    def description ( self , value ) :
        self.__description = str ( value )
        
    @property
    def patterns  ( self ) :
        """'patterns' : the list of patterns"""
        return tuple ( self.__patterns )
    
    @property
    def silent    ( self ) :
        """'silent' : silent processing?"""
        return self.__silent
    @silent.setter
    def silent    ( self , value ) :
        self.__silent = True if value else False

    @property
    def verbose ( self ) :
        """'verbose' : verbose processing?"""
        return not self.silent
    @verbose.setter 
    def verbose ( self , value ) :
        self.__silent = False if value else True 

    @property
    def maxfiles ( self ) :
        """'maxfiles' : maximal number of files to collect"""
        return self.__maxfiles
    
    @property
    def sorted ( self ) :
        """'sorted' : keep the list of files sorted?"""
        return self.__sorted

    @property
    def bad_files ( self ) :
        """'bad_files' : collection of bad/problematic files """
        return self.__bad_files 

    # =========================================================================
    ## get the disk size of files, only existing files are counted 
    def getsize ( self ) :
        """ Get the disk size of files, only existing/valid files are counted"""
        s = 0
        for f in self.__files :
            s += max ( 0 , self.get_file_size ( f ) ) 
        return s 
    
    # ==========================================================================
    ## Get the size of the file
    #  @param fname input file name
    #  @return ROOT file sizeor -1 for non-existingt/invalid files
    def get_file_size ( self , fname ) :
        """ Get the size of the file
        - fname : input file name
        - return file size or -1 for non-existing/invalid files
        """
        if os.path.exists ( fname ) and os.path.isfile ( fname ) :
            return os.path.getsize ( fname )

        return -1
    # =========================================================================
    ## check if the file is a part of collection
    #  @code
    #  files = ...
    #  if the_file in files : ....
    #  @endcode 
    def __contains__ ( self , the_file ) :
        """ Check if the file is a part of collection
        >>> files = ...
        >>> if the_file in files : ....
        """
        return the_file in self.__files 
        
    # =========================================================================
    ## Get the list of files from the patterns 
    def the_files ( self ) :
        """ Get the list of files from the patterns"""
        
        _files = set ()
        
        for pattern in  self.patterns :
            
            for f in glob.iglob ( pattern ) :
                f = os.path.abspath  ( f )
                f = os.path.normpath ( f )                    
                _files.add ( f )
                
        return tuple ( sorted ( _files ) ) 

    # =========================================================================
    ## add files 
    def __add_files ( self , files , max_files = -1 ) :
        """ Add files/patterns to data collector
        """
        
        if isinstance ( files  , str ) : files  = [ files  ]

        ## eliminate duplicates and sort 
        files = tuple ( sorted ( set ( files ) ) )
        
        nfiles    = len ( files )
        max_files = max_files if 0 <= max_files <= nfiles else nfiles 
        
        from ostap.utils.progress_bar import progress_bar
        for f in progress_bar ( files , silent = self.silent ) :
            
            if max_files <= len ( self.files ) :
                logger.debug ('Max-files limit is reached %s ' % max_files )
                break

            ## treat the file
            self.treatFile ( f )

        if self.sorted :
            self.__files.sort () 

    ## the specific action for each file 
    def treatFile ( self, the_file ) :
        if not the_file in self.__files :
            self.__files.append ( the_file )

    # ===============================================================================
    ## clone it! 
    def  clone ( self               ,
                 files       = None ,
                 description = None ) :
        """ Clone the object
        """
        import copy
        result = copy.copy ( self )
        
        if not files       is None :
            if isinstance ( files , str ) : files = files , 
            result.__files        =  [ f for f in files ] 
            result.__patterns     = ()
            if result.sorted : result.__files.sort()            
        if not description is None :
            result.descrpiption = str ( description )
            
        return result
    
    # =========================================================================
    ## check operations
    def check_ops ( self , other ):
        """Check operations
        """
        return isinstance ( other , Files ) 

    # =========================================================================
    ## check that the list of files is the same
    #  @code
    #  ds1 = ...
    #  ds2 = ...\
    #  ds1 == ds2 
    #  @endcode 
    def __eq__ ( self , other ) :
        """Check that the list of fiels is the same
        >>> ds1 = ...
        >>> ds2 = ...
        >>> ds1 == ds2 
        """
        if not self.check_ops ( other ) : return False 
        ##
        return set ( self.files ) == set ( other.files ) 

    ## merge two sets together
    def __or__ ( self , other ) :
        """ Merge two sets together
        >>> ds1 = ...
        >>> ds2 = ...
        >>> ds  = ds1 + ds2 
        >>> ds  = ds1 | ds2 ## ditto
        """
        if not self.check_ops ( other ) : return NotImplemented

        files       = self.files + tuple ( f for f in other.files if not f in self.files ) 
        description = "|".join ( "(%s)" for s in ( self.description , other.description ) )
        
        result = self.clone ( files = files , description = description )        
        for f in other.bad_files : result.bad_files.add ( f )
        return result 
    
    ## get an intersection of two datasets 
    def __and__ (  self , other ) :
        """ Get intersection of two sets
        >>> ds1 = ...
        >>> ds2 = ...
        >>> ds  = ds1 & ds2 ## get intersection 
        >>> ds  = ds1 * ds2 ## ditto
        """
        if not self.check_ops ( other ) : return NotImplemented

        files       = tuple ( f for f in self.files if f in other.files )
        description = "&".join ( "(%s)" for s in ( self.description , other.description ) )
        
        result = self.clone ( files = files , description = description )        
        for f in other.bad_files : result.bad_files.add ( f )
        return result 

    ## get an exclusive OR for two datasets 
    def __xor__ (  self , other ) :
        """ Get an exclusive OR for  two sets
        >>> ds1 = ...
        >>> ds2 = ...
        >>> ds  = ds1 ^ ds2 ## get an exclusive OR 
        """
        if not self.check_ops ( other ) : return NotImplemented

        files1 = [ f for f in self .files if not f in other.files ]
        files2 = [ f for f in other.files if not f in self.files  ]
        files  = tuple ( files1 + files2 ) 
        
        description = "^".join ( "(%s)" for s in ( self.description , other.description ) )
        
        result = self.clone ( files = files , description = description )        
        for f in other.bad_files : result.bad_files.add ( f )
        return result 

    __add__  = __or__    
    __mul__  = __and__

    ## append with another dataset
    def __ior__ ( self , other ) :
        """ Append with another dataset 
        >>> ds1 = ...
        >>> ds2 = ...
        >>> ds  |= ds2 
        >>> ds  += ds2 ## ditto
        """
        if not self.check_ops ( other ) : return NotImplemented

        files       = self.files + tuple ( f for f in other.files if not f in self.files )
        
        self.__files       = list ( files ) 
        self.__description = "|".join ( "(%s)" for s in ( self.description , other.description ) )

        for f in other.bad_files : self.bad_files.add ( f )        
        return self

    __iadd__  = __ior__    
    
    ## get union of two datasets 
    def union ( self , other ) :
        """ Union of two datasets"""
        return self | other
    
    ## get intersection of two datasets 
    def intersection ( self , other ) :
        """ Intersection of two datasets"""
        return self & other

    ## get difference of two datasets 
    def difference  ( self , other ) :
        """ Differfence of two datasets"""
        return self - other

    ## get symmetric difference of two datasets 
    def symmetric_difference  ( self , other ) :
        """ Symmetric Differfence of two datasets"""
        return self ^ other

    ## subtraction for datasets 
    def __sub__ (  self , other ) :
        """ Get subtraction of two sets
        >>> ds1 = ...
        >>> ds2 = ...
        >>> ds  = ds1 - ds2 ## get subtraction 
        """
        if not isinstance ( other , Files ) : return NotImplemented

        files       = tuple ( f for f in self.files if not f in other.files )
        description = "-".join ( "(%s)" for s in ( self.description , other.description ) ) 

        return self.clone ( files = files , description = description )

    ## remove the files from  another dataset 
    def __isub__ (  self , other ) :
        """ Remove file from another dataset 
        >>> ds1 = ...
        >>> ds2 = ...
        >>> ds1 -= ds2 ## get subtraction 
        """
        if not isinstance ( other , Files ) : return NotImplemented

        files       = tuple ( f for f in self.files if not f in other.files )

        self.__files       = list ( files ) 
        self.__description = "-".join ( "(%s)" for s in ( self.description , other.description ) )
        
        for f in other.bad_files : self.bad_files.add ( f )        
        return self

    # =========================================================================
    ## get a sample of at most  n-elements (if n is integer and >=1 )  or n-fraction 
    def sample_files ( self ,  n , sort ) :
        """ Get a sample of at most  n-elements (if n is integer and >=1 )  or n-fraction 
        """
        if   isinstance ( n , int   ) and 1 <= n <= len ( self.files ) :
            files = random.sample ( self.files , n )
            if sort : files.sort()
            return tuple ( files ) 
        elif isinstance ( n , float ) and 0 < n < 1 :
            n = int ( math.ceil ( len ( self.files ) * n ) )
            return self.sample_files ( n , sort ) 
        else  :
            raise IndexError ( "Invalid fraction: %s/%s" % ( n , len( self.files ) ) )

    # =========================================================================
    ##  Get a sub-sample
    #   @code
    #   files = ...
    #   f1 = files.sample ( 5   ) ##  5     files
    #   f2 = files.samlpe ( 0.1 ) ## 10% of files 
    #   @endcode
    def sample ( self , n , sort = True , **kwargs ) :
        """ Get a sub-sample
        >>> files = ...
        >>> f1 = files.sample ( 5   ) ##  5     files
        >>> f2 = files.sample ( 0.1 ) ## 10% of files 
        """
        files       = self.sample_files ( n , sort )
        description = "Sample(%d): %s" % ( n , self.description )
        
        return self.clone ( files = files , description = description )

    # =========================================================================
    ##  Get an element or slice 
    #   @code
    #   files = ...
    #   f1 = files[5] 
    #   f2 = files[4:10]
    #   @endcode
    def __getitem__ ( self , item ) :
        """ Get a sub-sample
        >>> files = ...
        >>> f1 = files[5] 
        >>> f2 = files[4:10]
        """
        
        files       = self.files [ item ]
        if isinstance ( item , int ) : 
            description = "Item(%d): %s" % ( item , self.description )
        else :
            description = "%s: %s" % ( item , self.description )
        
        return self.clone ( files = files , description = description  )
        
    ## printout 
    def __str__(self):
        """ The specific printout
        """
        return "<#files: %4d>" % len ( self.files ) 
    
    def __repr__    ( self ) : return self.__str__()
    def __nonzero__ ( self ) : return bool ( self.files )
    def __bool__    ( self ) : return bool ( self.files )
    def __len__     ( self ) : return len  ( self.files )

    # =========================================================================
    ## Print collection of files as table
    #  @code
    #  files = ...
    #  print ( files.table() )    
    #  @endcode
    def table ( self , title = 'Files' , prefix = '' ) :
        """ Print collection of files as table
        """
        rows  = [ ( '#' , 'size' , 'name' ) ]
        files = self.files
        bad   = sorted ( self.bad_files ) 
        nn    = max ( len ( files ) , len ( bad ) ) 
        nfmt  = '%%%dd' % ( math.floor ( math.log10 ( nn ) ) + 1 )

        from itertools  import chain
        total_size =  0 
        for i , f in enumerate ( chain ( files , bad ) , start = 1 ) : 

            row   = [ nfmt % i ]
            fsize = self.get_file_size ( f )
            if 0 <= fsize :
                total_size += fsize 
                vv , unit   = fsize_unit ( fsize )
                row.append ( '%3d %s' % ( vv , unit ) ) 
            else : 
                row.append ( '???' )
                
            row .append ( f   ) 
            rows.append ( row )

        ## the last row: summary
        from ostap.logger.colorized import infostr 
        vv , unit  = fsize_unit ( total_size  )
        row   = '' , \
            infostr ( '%3d %s' % ( vv , unit ) ) , \
            infostr ( self.commonpath )  
        rows.append ( row )
        
        import ostap.logger.table as T
        return T.table ( rows , title = title , prefix = prefix , alignment = 'rrw' ) 

    # =========================================================================
    ## Split object into several chunks of smaller size
    #  @code
    #  data   = ...
    #  chunks = data.split_chunks ( 10 ) 
    #  @endcode 
    def split_chunks ( self , chunk_size = 10 ) :
        """ Split object into several chunks of smaller size
        >>> data   = ...
        >>> chunks = data.split_chunks ( 10 ) 
        """
        from ostap.utils.utils import chunked 
        return tuple ( self.clone ( files = chunk ) for chunk in chunked ( self.files , chunk_size ) )

    # =========================================================================
    ## Split object into several chunks 
    #  @code
    #  data   = ...
    #  chunks = data.split_groups ( 10 ) 
    #  @endcode 
    def split_groups ( self , groups = 10 ) :
        """ Split the object into several groups
        >>> data   = ...
        >>> groups = data.split_groups ( 10 ) 
        """
        from ostap.utils.utils import divide
        result = []
        for group in divide ( groups , self.files ) :
            files = list ( group ) 
            if not files : continue
            result.append ( self.clone ( files = files ) )
            
        return tuple  ( result )
    
    # =========================================================================
    ## get a common path (prefix) for all files in collection
    @property 
    def commonpath ( self ) :
        """'commonpath': common path (prefix) for all files in collection
        """
        ##
        
        files = self.files        
        if not files : return '' 

        from ostap.utils.basic import commonpath 
        cp = commonpath ( files )
        return cp if os.path.isdir ( cp ) else os.path.dirname ( cp ) 

    # =========================================================================
    ## copy all the files to new directory
    #  - new directory will be created (if needed)
    #  - common path (prefix) for all files will be replaced by new directory
    def copy_files ( self , new_dir = None , parallel = False , also_bad = False ) :
        """ Copy all the files to new directory
        - new directory will be created (if needed)
        - common path (prefix) for all files will be replaced by new directory
        """

        files_to_copy = set ( self.__files )
        if also_bad and self.bad_files :
            files_to_copy |= set ( self.bad_files )

        ## use generic function 
        copied = copy_files ( files_to_copy              ,
                              new_dir  = new_dir         ,   
                              parallel = parallel        ,
                              progress = not self.silent )
        
        if not self.silent :
            from ostap.utils.basic import commonpath 
            cp = commonpath ( copied  )
            logger.info ( "copy_files: #%d files are copied to '%s'" %  ( len ( copied ) , cp ) )

        return self.clone ( files = sorted ( copied ) ) 


    # =========================================================================
    ## Sync/copy all the files to new directory
    #  - new directory will be created (if needed)
    #  - common path (prefix) for all files will be replaced by new directory
    def sync_files ( self , new_dir = None , parallel = False , also_bad = False ) :
        """ Sync/copy all the files to new directory
        - new directory will be created (if needed)
        - common path (prefix) for all files will be replaced by new directory
        """

        files_to_copy = set ( self.__files )
        if also_bad and self.bad_files :
            files_to_copy |= set ( self.bad_files )

        ## use generic function 
        copied = sync_files ( files_to_copy              ,
                              new_dir  = new_dir         ,   
                              parallel = parallel        ,
                              progress = not self.silent )
        
        if not self.silent :
            from ostap.utils.basic import commonpath 
            cp = commonpath ( copied  )
            logger.info ( "sync_files: #%d files are sync/copied to '%s'" %  ( len ( copied ) , cp ) )

        return self.clone ( files = sorted ( copied ) ) 
        
# =========================================================================
## Copy set of files into new directory.
#  Files are copied in a way to preserve they name uniqueness:
#  - a.txt                -> NEWDIR/a.txt
#  - subdira.txt          -> NEWDIR//subdira.txt
#  - subdir/subdir/a.txt  -> NEWDIR/subdir/subdir.txt
#  Essentially a common prefix for all input files will be replaced
#  by the destination directory 
#  @param file_to_copy sequence of files to be copied or sequence of (source,destination) pairs 
#  @param new_dir destination directory, for None temproary directory will be used. 
#  @param copier the low-level copy routine ot be used 
#  @param progress show progrees bar if possible
#  @return list of copied files 
def copy_files ( files_to_copy           ,
                 new_dir   = None        ,
                 parallel  = False       ,
                 copier    = None        ,
                 copy_cmd  = ''          , 
                 progress  = False       ) :
    """ Copy set of files into new directory.

    Files are copied in a way to preserve they name uniqueness:

    - a.txt                -> NEWDIR/a.txt
    - subdira.txt          -> NEWDIR//subdira.txt
    - subdir/subdir/a.txt  -> NEWDIR/subdir/subdir.txt

    Essentially a common prefix for all input files is replaced  
    by the destination directory 

    - file_to_copy sequence of files to be copied or seuqence of (source,destination) pairs 
    - new_dir      destination directory, for None temproary directory wil lbe used
    - copier       low-level copy routine ot be used 
    - progress     show progrees bar if possible
    
    A list of copied files is returned 
    """
    
    the_files = [ f for f in files_to_copy ] 

    ## 1) no input 
    if   not the_files : pairs = [] 
    ## 1) input: list of files :
    elif all ( isinstance ( f , path_types ) and os.path.exists ( f ) and os.path.isfile ( f ) for f in the_files ) :

        ## use  temporary directory
        if new_dir is None or new_dir == '' :
            import ostap.utils.cleanup as CU
            new_dir = CU.CleanUp.tempdir()
            
        ## create directory if needed 
        if not os.path.exists ( new_dir ) : os.makedirs ( new_dir )
        
        from ostap.utils.basic  import writeable
        assert writeable ( new_dir ), \
            "copy_files: the destination directory `%s' is not writable!" % new_dir 
        
        nd = os.path.abspath  ( new_dir )
        nd = os.path.normpath ( nd      ) 
        
        ## get the common path.prefix for all input files 
        from ostap.utils.basic import commonpath 
        cp = commonpath ( os.path.abspath ( f ) for f in the_files )
        cp = cp if os.path.isdir ( cp ) else os.path.dirname ( cp ) 
        
        pairs = []
        for f in files_to_copy :
            fs   = os.path.abspath  ( f )   
            nf   = fs.replace       ( cp , nd ) 
            nf   = os.path.abspath  ( nf )
            pair = f , nf
            pairs.append ( pair )
    
    ## 2) input: list of pairs 
    elif  all ( isinstance ( f , sized_types )          \
                and 2 == len ( f )                      \
                and isinstance ( f [ 0 ] , path_types ) \
                and isinstance ( f [ 1 ] , path_types ) \
                and os.path.exists ( f [ 0 ] )          \
                and os.path.isfile ( f [ 0 ] ) for f in the_files ) :

        pairs = the_files
        
    else :
        raise TypeError ( "Invalid input file list" ) 
                                  

    if not copier :
        from ostap.utils.basic  import copy_file
        copier = copy_file

    from   ostap.utils.basic      import numcpu
    from   ostap.utils.utils      import which
    
    nfiles = len ( pairs  )
    if parallel and 1 < nfiles and 1 < numcpu () :
        
        from ostap.utils.parallel_copy import copy_files as parallel_copy        
        copied = parallel_copy ( pairs , maxfiles = 1 , copier = copier , copy_cmd = copy_cmd , progress = progress  )
        copied = tuple ( f [ 1 ] for f in copied )

    else :
        
        ## sequential copy file-by-file
        from ostap.utils.progress_bar import progress_bar
        silent = nfiles <= 1 or not progress 
        copied = [] 
        for f, nf in progress_bar ( pairs , silent = silent ) :
            result = copier ( f , nf , progress = progress and nfiles <=1 ) 
            copied.append ( result ) 
            
    return tuple ( copied )

# =========================================================================
## Sync/copy set of files into new directory.
#  Files are copied in a way to preserve they name uniqueness:
#  - a.txt                -> NEWDIR/a.txt
#  - subdira.txt          -> NEWDIR//subdira.txt
#  - subdir/subdir/a.txt  -> NEWDIR/subdir/subdir.txt
#  Essentially a common prefix for all input files will be replaced
#  by the destination directory 
#  @param file_to_copy sequence of files to be copied or sequence of (source,destination) pairs 
#  @param new_dir destination directory, for None temproary directory wil lbe used
#  @param copier the low-level copy routine ot be used 
#  @param progress show progrees bar if possible
#  @return list of copied files 
def sync_files ( files_to_copy          ,
                 new_dir   = None       ,
                 parallel  = False      ,
                 copier    = None       ,
                 copy_cmd  = ''         , 
                 progress  = False      ) :
    
    """ Sync/copy set of files into new directory.

    Files are copied in a way to preserve they name uniqueness:

    - a.txt                -> NEWDIR/a.txt
    - subdira.txt          -> NEWDIR//subdira.txt
    - subdir/subdir/a.txt  -> NEWDIR/subdir/subdir.txt

    Essentially a common prefix for all input files is replaced  
    by the destination directory 

    - file_to_copy sequence of files to be copied or sequence of (source,destination) pairs 
    - new_dir      destination directory, for None temproary directory wil lbe used
    - copier       low-level copy routine ot be used 
    - progress      show progrees bar if possible
    
    A list of copied files is returned 
    """

    from   ostap.utils.utils  import which
    if not which ( 'rsync' ) :
        from ostap.utils.basic import copy_file
        return copy_files ( files_to_copy        ,
                            new_dir  = new_dir   ,
                            parallel = parallel  ,
                            copier   = copy_file ,
                            copy_cmd = ''        ,
                            progress = progress  )

    if not copy_cmd : copy_cmd = 'rsync -a'
    if not copier : 
        from ostap.utils.basic  import sync_file
        copier = sync_file
        
    return copy_files ( files_to_copy        ,
                        new_dir  = new_dir   ,
                        parallel = parallel  ,
                        copier   = copier    ,
                        copy_cmd = copy_cmd  ,
                        progress = progress  )

# =============================================================================
## Sync/copy two directories 
#  Files are copied in a way to preserve they name uniqueness:
#  - a.txt                -> NEWDIR/a.txt
#  - subdira.txt          -> NEWDIR//subdira.txt
#  - subdir/subdir/a.txt  -> NEWDIR/subdir/subdir.txt
#  Essentially a common prefix for all input files will be replaced
#  by the destination directory 
#  @param file_to_copy sequence of files to be copied
#  @param new_dir destination directory, for None temproary directory wil lbe used
#  @param copier the low-level copy routine ot be used 
#  @param progress show progrees bar if possible
#  @return list of copied files 
def sync_dirs  ( source_dir             ,
                 new_dir                ,
                 parallel  = False      ,
                 copier    = None       ,
                 copy_cmd  = ''         , 
                 progress  = False      ) :
    """ Sync/copy two directories 

    Files are copied in a way to preserve they name uniqueness:

    - a.txt                -> NEWDIR/a.txt
    - subdira.txt          -> NEWDIR//subdira.txt
    - subdir/subdir/a.txt  -> NEWDIR/subdir/subdir.txt

    Essentially a common prefix for all input files is replaced  
    by the destination directory 

    - file_to_copy sequence of files to be copied
    - new_dir      destination directory, for None temproary directory wil lbe used
    - copier       low-level copy routine ot be used 
    - progress      show progrees bar if possible
    
    A list of copied files is returned 
    """

    assert source_dir and os.path.exists ( source_dir ) and os.path.isdir ( source_dir ) , \
        'Invalid source directory %s' % source_dir 

    ## create directory if needed 
    if not os.path.exists ( new_dir ) : os.makedirs ( new_dir )
    
    from   ostap.utils.basic      import writeable
    assert os.path.exists ( new_dir ) and os.path.isdir ( new_dir ) and writeable ( new_dir ) , \
        'Invalid new directory %s' % new_dir 

    source = os.path.abspath  ( source_dir )
    source = os.path.normpath ( source     )

    the_files  = set() 
    ## collect all files in the source directory
    for root, dirs, files in os.walk ( source ) :
        root = os.path.abspath  ( root )
        root = os.path.normpath ( root )
        for f in files : the_files.add ( os.path.join ( root , f ) ) 

    newdir = os.path.abspath  ( new_dir    )
    newdir = os.path.normpath ( newdir     )

    pairs  = [ ( f, f.replace ( source , newdir ) ) for f in the_files ]
    pairs  = tuple ( pairs ) 

    return sync_files ( pairs               ,
                        new_dir  = new_dir  , 
                        parallel = parallel ,
                        copier   = copier   ,
                        copy_cmd = copy_cmd ,
                        progress = progress )

# =============================================================================



# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    logger.info ( 80*'*' ) 
            
# =============================================================================
##                                                                      The END 
# =============================================================================
