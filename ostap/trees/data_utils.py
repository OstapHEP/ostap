#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/trees/data.py 
#
#  Useful utilities to get access to datafiles & chains
#  Actualy it is just a little bit modified version (with globbing) of
#  original ``Ostap'' code by Alexander BARANOV
#
#  @code
#
#  >>> data  = Data('Bc/MyTree', '*.root' )
#  >>> chain = data.chain
#  >>> flist = data.files 
#
#  @endcode
# 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @author Alexander BARANOV a.baranov@cern.ch
#  @date   2014-06-08
#
# =============================================================================
""" Useful utilities to get access to datafiles & chains

Actualy it is a little bit modified version (with globbing) of
the original ``Ostap'' code by Alexander BARANOV

>>> data  = Data('Bc/MyTree', '*.root' )
>>> chain = data.chain
>>> flist = data.files 

>>> data  = Data('Bc/MyTree', 'a.root' )
>>> chain = data.chain
>>> flist = data.files 

>>> data  = Data('Bc/MyTree', [ 'a.root' , 'b.root' ] )
>>> chain = data.chain
>>> flist = data.files 

"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-06-08"
__all__     = (
    'Files'       , ## collect files  
    'Data'        , ## collect files and create     TChain
    'Data2'       , ## collect files and create two TChain objects 
    )
# =============================================================================
from   ostap.core.core import rootError, rootWarning
import ROOT, os, glob, random, math
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.trees.data' )
else                       : logger = getLogger ( __name__     )
# =============================================================================
if not hasattr ( ROOT.TTree , '__len__' ) :  
    ROOT.TTree. __len__ = lambda s : s.GetEntries()

# =============================================================================
## protocols for remote (ROOT) files 
protocols = (
    'root:'  ,
    'xroot:' ,
    'http:'  ,
    'https:' ,    
    )
# =============================================================================
## Does the filenamse starts with protocol?
def has_protocol ( fname ) :
    """Does the filenamse starts with protocol? """
    for p in protocols :
        if fname.startswith ( p ) : return True
    return False
# =============================================================================
## Strip protocl from the file name 
def strip_protocol ( fname ) :
    """Strip protocl from the file name"""
    for p in protocols :
        if fname.startswith ( p ) : return fname.replace ( p , '' ) 
    return fname

# =============================================================================
class DataProcessor (object) :
    def __init__   ( self , data  )  :
        self.data = data
    def initialize ( self ) :
        return self.data.clone()
    def __call__ ( self , items )  :
        data = self.data.clone() 
        data.add_files ( items )
        return data

from ostap.parallel.task import   Task

# =============================================================================
## @class DataTask
class DataTask(Task) : 
    def __init__ ( self , data ) :
        self.__data  = data.clone ( files = [] ) 
    def initialize_local  ( self              ) : pass 
    def initialize_remote ( self , jobid = -1 ) : pass
    ## actgual processing  
    def process ( self , jobid , files ) :
        """Actual processing"""
        for f in files :  self.__data.treatFile ( f ) 
        return self.__data 
    ## get the results 
    def results (  self ) : return self.__data 
    ## merge the results 
    def merge_results  ( self , results , jobid = -1 ) :    
        self.__data += results
        
# =============================================================================
## @class Files
#  Simple utility to pickup the list of files 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @author Alexander BARANOV a.baranov@cern.ch
#  @date   2014-06-08  
class Files(object):
    """Simple utility to pickup the list of files     
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
            
            task = DataTask ( self )
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
    ## check if the file is a part of collection
    #  @code
    #  files = ...
    #  if the_file in files : ....
    #  @endcode 
    def __contains__ ( self , the_file ) :
        """Check if the file is a part of collection
        >>> files = ...
        >>> if the_file in files : ....
        """
        return the_file in self.__files 
        
    # =========================================================================
    ## Get the list of files from the patterns 
    def the_files ( self ) :
        """Get the list of files from the patterns"""
        
        _files = set ()
        
        for pattern in  self.patterns :
            
            _added = False  
            for p in protocols :
                if pattern.startswith ( p ) :
                    _files.add ( pattern )
                    _added = True 
                    break
                
            if not _added : 
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
    ## check that the list of fiels is the same
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
        """ get intersection of two sets
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
        """ get an exclusive OR for  two sets
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
        """Union of two datasets"""
        return self | other
    
    ## get intersection of two datasets 
    def intersection ( self , other ) :
        """Intersection of two datasets"""
        return self & other

    ## get difference of two datasets 
    def difference  ( self , other ) :
        """Differfence of two datasets"""
        return self - other

    ## get symmetric difference of two datasets 
    def symmetric_difference  ( self , other ) :
        """Symmetric Differfence of two datasets"""
        return self ^ other

    ## subtraction for datasets 
    def __sub__ (  self , other ) :
        """ get subtraction of two sets
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
        """ Remove file forom anothere dataset 
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
        """get a sample of at most  n-elements (if n is integer and >=1 )  or n-fraction 
        """
        if   isinstance ( n , int   ) and 1 <= n <= len ( self.files ) :
            files = random.sample ( self.files , n )
            if sort : files.sort()
            return tuple ( files ) 
        elif isinstance ( n , float ) and 0 < n < 1 :
            n = int ( math.floor ( len ( self.files ) * n ) )
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
        """Get a sub-sample
        >>> files = ...
        >>> f1 = files.sample ( 5   ) ##  5     files
        >>> f2 = files.sample ( 0.1 ) ## 10% of files 
        """
        files       = self.sample_files ( n , sort )
        description = "Sample(%d): %s" % ( n , self.description )
        
        return self.clone ( files = files , description = description , patterns = () )
    
    # =========================================================================
    ##  Get an element or slice 
    #   @code
    #   files = ...
    #   f1 = files[5] 
    #   f2 = files[4:10]
    #   @endcode
    def __getitem__ ( self , item ) :
        """Get a sub-sample
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
        """The specific printout
        """
        return "<#files: %4d>" % len ( self.files ) 
    
    def __repr__    ( self ) : return self.__str__()
    def __nonzero__ ( self ) : return bool ( self.files )
    def __bool__    ( self ) : return bool ( self.files )
    def __len__     ( self ) : return len  ( self.files )

    # =========================================================================
    ## merge all files using <code>hadd</code> script from ROOT
    #  @param output  name of the output merged file, if None,
    #                 the temporary name will be generated,
    #                 that will be deleted at the end of the session
    #  @param opts   options for command <code>hadd</code>
    #  @return the name of the merged file
    # OPTIONS:
    # -a                                   Append to the output
    # -k                                   Skip corrupt or non-existent files, do not exit
    # -T                                   Do not merge Trees
    # -O                                   Re-optimize basket size when merging TTree
    # -v                                   Explicitly set the verbosity level: 0 request no output, 99 is the default
    # -j                                   Parallelize the execution in multiple processes
    # -dbg                                 Parallelize the execution in multiple processes in debug mode (Does not delete partial files stored inside working directory)
    # -d                                   Carry out the partial multiprocess execution in the specified directory
    # -n                                   Open at most 'maxopenedfiles' at once (use 0 to request to use the system maximum)
    # -cachesize                           Resize the prefetching cache use to speed up I/O operations(use 0 to disable)
    # -experimental-io-features            Used with an argument provided, enables the corresponding experimental feature for output trees
    # -f                                   Gives the ability to specify the compression level of the target file(by default 4) 
    # -fk                                  Sets the target file to contain the baskets with the same compression
    #                                      as the input files (unless -O is specified). Compresses the meta data
    #                                      using the compression level specified in the first input or the
    #                                      compression setting after fk (for example 206 when using -fk206)
    # -ff                                  The compression level use is the one specified in the first input
    # -f0                                  Do not compress the target file
    # -f6                                  Use compression level 6. (See TFile::SetCompressionSettings for the support range of value.)    
    def hadd ( self , output = None , opts = "-ff" ) :
        """<erge all files using <code>hadd</code> script from ROOT
        - `output`  name of the output merged file
        - `opts`   options for command <code>hadd</code>
        It returns the name of the merged file
        
        If no output file name is specified, the temporary name
        will be generate and the temporary file will be deleted
        at the end of the session

        OPTIONS:
        # -a                                   Append to the output
        # -k                                   Skip corrupt or non-existent files, do not exit
        # -T                                   Do not merge Trees
        # -O                                   Re-optimize basket size when merging TTree
        # -v                                   Explicitly set the verbosity level: 0 request no output, 99 is the default
        # -j                                   Parallelize the execution in multiple processes
        # -dbg                                 Parallelize the execution in multiple processes in debug mode (Does not delete partial files stored inside working directory)
        # -d                                   Carry out the partial multiprocess execution in the specified directory
        # -n                                   Open at most 'maxopenedfiles' at once (use 0 to request to use the system maximum)
        # -cachesize                           Resize the prefetching cache use to speed up I/O operations(use 0 to disable)
        # -experimental-io-features            Used with an argument provided, enables the corresponding experimental feature for output trees
        # -f                                   Gives the ability to specify the compression level of the target file(by default 4) 
        # -fk                                  Sets the target file to contain the baskets with the same compression
        #                                      as the input files (unless -O is specified). Compresses the meta data
        #                                      using the compression level specified in the first input or the
        #                                      compression setting after fk (for example 206 when using -fk206)
        # -ff                                  The compression level use is the one specified in the first input
        # -f0                                  Do not compress the target file
        # -f6                                  Use compression level 6. (See TFile::SetCompressionSettings for the support range of value.)                            
        """
        from ostap.utils.utils import hadd as hadd_
        return hadd_ ( self.files , ouput = output , opts = opts  )
    
    # =========================================================================
    ## get a common path (prefix) for all files in collection
    #  - protocols are ignored 
    @property 
    def commonpath ( self ) :
        """'commonpath': common path (prefix) for all files in collection
        - protocols are ignored 
        """
        from ostap.utils.basic import commonpath 
        ##
        if any ( has_protocol ( f ) for f in self.__files ) :
            files = []
            for f in self.__files :
                files .append ( strip_protocol ( f ) ) 
        else : files = self.__files 
        cp = commonpath ( files )
        return cp if os.path.isdir ( cp ) else os.path.dirname ( cp ) 
    # =========================================================================

    # =========================================================================
    ## copy all the files to new directory
    #  - new directory will be created (if needed)
    #  - common path (prefix) for all files will be replaced by new directory
    def copy_files ( self , new_dir , parallel = False , also_bad = False ) :
        """copy all the files to new directory
        - new directory will be created (if needed)
        - common path (prefix) for all files will be replaced by new directory
        """
        
        from ostap.utils.basic  import writeable,    copy_file
        from ostap.io.root_file import copy_file as copy_root_file 

        ## create directory if needed 
        if not os.path.exists ( new_dir ) : os.makedirs ( new_dir )
        
        assert writeable ( new_dir ), \
               "New directory ``%s'' is not writable!" % new_dir 

        nd = os.path.abspath  ( new_dir )
        nd = os.path.normpath ( nd      ) 
        nd = os.path.realpath ( nd      )
        
        cp = self.commonpath

        files_to_copy = set ( self.__files )
        if also_bad and self.bad_files :
            files_to_copy |= set ( self.bad_files )
            
        if parallel :

            regular_files  = [] 
            root_files     = [] 
            for f in files_to_copy :
                fs = os.path.normpath ( strip_protocol ( f ) ) 
                nf = fs.replace ( cp , nd ) 
                nf = os.path.normpath ( nf )
                pair = f , nf 
                if has_protocol ( f ) : root_files   .append ( pair )
                else                  : regular_files.append ( pair ) 
                
            from ostap.parallel.parallel_copy import copy_files as parallel_copy
            copied1 = []
            copied2 = []
            
            if regular_files :
                copied1 = parallel_copy ( regular_files , maxfiles = 1 , copier = copy_file      , progress = not self.silent )
            if root_files :
                copied2 = parallel_copy (    root_files , maxfiles = 1 , copier = copy_root_file , progress = not self.silent )
                
            copied  = [ d[1] for d in copied1 ] + [ d[1] for d in copied2 ] 

        else :
            
            copied = []
            from ostap.utils.progress_bar import progress_bar
            nf = len ( files_to_copy ) 
            for f in progress_bar ( files_to_copy , silent = self.silent or nf <=1 ) :
                
                fs = os.path.normpath ( strip_protocol ( f ) ) 
                nf = fs.replace ( cp , nd ) 
                nf = os.path.normpath ( nf )
                
                if not has_protocol ( f ) :
                    result = copy_file      ( f , nf , progress = ( 1 == nf ) and self.verbose ) 
                else                      :
                    result = copy_root_file ( f , nf , progress = ( 1 == nf ) and self.verbose ) 
                    
                copied.append ( result )
                
        copied = tuple ( copied )
        return self.clone ( files = copied ) 
    
# =============================================================================
## @class Data
#  Simple utility to access to certain chain in the set of ROOT-files
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @author Alexander BARANOV a.baranov@cern.ch
#  @date   2014-06-08  
class Data(Files):
    """Simple utility to access to certain chain in the set of ROOT-files
    >>> data  = Data('Bc/MyTree', '*.root' )
    >>> chain = data.chain
    >>> flist = data.files 
    """
    def __init__( self                 ,
                  chain                ,
                  files        = []    ,
                  description  = ''    , 
                  maxfiles     = -1    ,
                  check        = True  , 
                  silent       = False ,
                  sorted       = True  , 
                  parallel     = False ) : 

        ## we will need Ostap machinery for trees&chains here
        import ostap.trees.trees 

        ## decorate files 
        if isinstance ( files , str ) : files = [ files ]

        self.__check = True if check else False
        
        ## chain name 
        self.__chain_name = chain if isinstance ( chain , str ) else chain.name 

        ## update description
        if not description : description = "ROOT.TChain(%s)" % self.chain_name 
        
        ## initialize the  base class 
        Files.__init__( self , files , description  , maxfiles , silent = silent , sorted = sorted , parallel = parallel )

    @property
    def validate ( self ) :
        """'check' : make check for `TTree`/`TChain` structures
        """
        return self.__check
    
    @property
    def chain_name ( self ) :
        """'chain_name' : the name of TTree/TChain object
        """
        return self.__chain_name
    
    @property
    def chain ( self ) :
        """'chain' : (re)built and return `TChain` object"""
        ch = ROOT.TChain( self.chain_name )
        for f in self.files : ch.Add ( f )
        return ch
    
    @property
    def check ( self ) :
        """'check' : make check for `TTree`/`TChain` structures
        """
        return self.__check

    ## check the content of the two trees
    @staticmethod 
    def check_trees ( tree1 , tree2 , the_file = '' ) :

        if tree1 and tree2 :
            
            branches1 = set ( tree1.branches () )                
            leaves1   = set ( tree1.leaves   () )
            
            branches2 = set ( tree2.branches () )                
            leaves2   = set ( tree2.leaves   () )

            if branches1 != branches2 :
                missing = list ( sorted ( branches1 - branches2 ) ) 
                extra   = list ( sorted ( branches2 - branches1 ) ) 
                if missing : logger.warning ( "Tree('%s'): %3d missing branches: %s in %s" %  ( tree1.GetName() , len ( missing ) , missing , the_file ) )
                if extra   : logger.warning ( "Tree('%s'): %3d extra   branches: %s in %s" %  ( tree1.GetName() , len ( extra   ) , extra   , the_file ) )
                
            if ( ( branches1 != leaves1 ) or ( branches2 != leaves2 ) ) and leaves1 != leaves2 :
                missing = list ( sorted ( leaves1 - leaves2 ) )
                extra   = list ( sorted ( leaves2 - leaves1 ) ) 
                if missing : logger.warning ( "Tree('%s'): %3d missing leaves:   %s in %s" %  ( tree1.GetName() , len ( missing ) , missing , the_file ) )
                if extra   : logger.warning ( "Tree('%s'): %3d extra   leaves:   %s in %s" %  ( tree1.GetName() , len ( extra   ) , extra   , the_file ) )

                
    ## the specific action for each file 
    def treatFile ( self, the_file ) :
        """Add the file to TChain
        """

        ## suppress Warning/Error messages from ROOT 
        with rootError() :

            ## new temporary chain/tree for this file 
            tree  = ROOT.TChain ( self.chain_name )
            tree.Add ( the_file )

            if 1 <= len ( tree ) and 1 <= len ( tree.branches () ) :
                
                if self.check:
                    files = self.files 
                    if files :
                        chain = ROOT.TChain ( self.chain_name )
                        chain.Add ( files [ 0 ] ) 
                        self.check_trees ( tree , chain , the_file )
                        del chain
                        
                Files.treatFile ( self , the_file )
                
            else :
                
                self.bad_files.add ( the_file )
                if not self.silent : 
                    logger.warning ( "No/empty chain  '%s' in file '%s'" % ( self.chain_name , the_file ) )
                    
            del tree

            
    ## printout 
    def __str__(self):
        with rootWarning() :
            nf = len ( self.files     )
            nc = '??'
            if not self.silent :
                chain = self.chain
                nc    = len ( chain )
                del chain 
            ne = len ( self.bad_files )
        return "<#files: {}; Entries: {}>"              .format ( nf , nc ) if not self.bad_files else \
               "<#files: {}; Entries: {}; No/empty: {}>".format ( nf , nc  ,  ne )
    
    # =========================================================================
    ## check operations
    def check_ops ( self , other ):
        """Check operations
        """
        return isinstance ( other , Data ) and self.chain_name == other.chain_name

    # =========================================================================
    ## get DataFrame for the chain
    #  @see ROOT::RDataFrame
    @property
    def frame ( self ) :
        """'frame': get the DataFrame for the chain"""
        from   ostap.frames.frames import DataFrame, frame_progress 
        f = DataFrame ( self.chain )
        if not self.filent:
            pb = frame_progress ( f , len ( self.chain ) ) 
        return f 
    
# =============================================================================
## @class Data2
#  Simple utility to access two chains in the set of ROOT-files
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @author Alexander BARANOV a.baranov@cern.ch
#  @date   2014-06-08  
class Data2(Data):
    """Simple utility to access to certain chain in the set of ROOT-files    
    >>> data  = Data('Bc/MyTree', '*.root' )
    >>> chain = data.chain
    >>> flist = data.files 
    """
    
    def __init__( self                ,
                  chain1              ,
                  chain2              , 
                  files       = []    ,
                  description = ''    ,
                  maxfiles    = -1    ,
                  check       = True  , 
                  silent      = False ,
                  sorted      = True  , 
                  parallel    = False ) :

        ## decorate files 
        if isinstance ( files , str ) : files = [ files ]


        ## chain name 
        self.__chain2_name = chain2 if isinstance ( chain2 , str ) else chain2.name 

        if not description :
            description = chain1.GetName() if hasattr ( chain1 , 'GetName' ) else str(chain1)
            description = "%s&%s" % ( description , self.chain2.GetName() )

        Data.__init__( self                      ,
                       chain       = chain1      ,
                       files       = files       ,
                       description = description ,
                       maxfiles    = maxfiles    ,
                       check       = check       ,
                       silent      = silent      ,
                       sorted      = sorted      , 
                       parallel    = parallel    ) 
        
    @property 
    def files1    ( self ) :
        """'files1' : the list of files fore the first chain (same as 'files') """
        return self.files 

    @property 
    def files2    ( self ) :
        """'files2' : the list of files for the second chain (the same)"""
        return self.files

    @property
    def chain1_name ( self ) :
        """'chain1_name' : the name of the first TTree/TChain object (same as 'chain_name')
        """
        return self.chain_name

    @property
    def chain2_name ( self ) :
        """'chain2_name' : the name of the second TTree/TChain object
        """
        return self.__chain2_name

    @property
    def chain1 ( self ) :
        """'chain1' : (re)built and return the first `TChain` object (same as `chain')
        """
        return self.chain
    
    @property
    def chain2 ( self ) :
        """'chain2' : (re)built and return the second `TChain` object
        """
        ch = ROOT.TChain( self.chain2_name )
        for f in self.files2 : ch.Add ( f )
        return ch
    
    ## the specific action for each file 
    def treatFile ( self, the_file ) :
        """Add the file to TChain
        """
        ## suppress Warning/Error messages from ROOT 
        with rootError() :
            
            tree1 = ROOT.TChain ( self.chain1_name  )
            tree1.Add ( the_file )
            
            tree2 = ROOT.TChain ( self.chain2_name  )
            tree2.Add ( the_file )

            ok1 = 1 <= len ( tree1 ) and 1 <= len ( tree1.branches () ) 
            if self.check and ok1 :
                files = self.files
                if files :                    
                    ch1 = ROOT.TChain( self.chain1_name )
                    ch1.Add ( files [ 0 ] ) 
                    self.check_trees ( tree1 , ch1 , the_file )
                    del ch1
                
            ok2 = 1 <= len ( tree2 ) and 1 <= len ( tree2.branches () )
            if self.check and ok2 :
                files = self.files
                if files :                    
                    ch2 = ROOT.TChain( self.chain2_name )
                    ch2.Add ( files [ 0 ] ) 
                    self.check_trees ( tree2 , ch2 , the_file )
                    del ch2

            if  ok1 and ok2      :
                
                Files.treatFile ( self , the_file ) 
                
            elif ok2 :
                
                self.bad_files.add ( the_file  )
                if not self.silent : 
                    logger.warning ( "No/empty chain1 '%s'      in file '%s'" % ( self.chain1_name , the_file ) )

            elif ok1 :
                
                self.bad_files.add ( the_file )
                if not self.silent : 
                    logger.warning ( "No/empty chain2 '%s'      in file '%s'" % ( self.chain2_name , the_file ) ) 

            else :
                
                self.bad_files.add ( the_file )
                if not self.silent :                 
                    logger.warning ( "No/empty chains '%s'/'%s' in file '%s'" % ( self.chain1_name ,
                                                                                  self.chain2_name , the_file ) )
            del tree1
            del tree2
            
    ## printout 
    def __str__(self):

        with rootWarning() :
            nf  = len ( self.files      )
            nf2 = len ( self.files2     )
            
            nc   = '??'
            nc2  = '??'
            
            if not self.silent :
                chain1 = self.chain1
                nc     = len ( chain1 )
                del  chain1
                chain2 = self.chain2
                nc     = len ( chain2 )
                del  chain2
                
            ne  = len ( self.bad_files1 )
            
        sf  =  set ( self.files ) == set ( self.files2 )
        
        if not self.bad_files :
            return "<#files: {}; Entries: {}/{}>"   .format ( nf, nc , nc2 ) if sf else \
                   "<#files: {}/{}; Entries: {}/{}>".format ( nf , nf2 , nc , nc2 )
        else :
            return "<#files: {}; Entries: {}/{}; No/empty :{}>"   .format ( nf , nc , nc2 , ne ) if sf else \
                   "<#files: {}/{}; Entries: {}/{}; No/empty :{}>".format ( nf , nf2 , nc , nc2 , ne )
        

    def __nonzero__ ( self )  :
        return bool ( self.files ) and bool ( self.files2 )


    # =========================================================================
    ## check operations
    def check_ops ( self , other ):
        """Check operations
        """
        return isinstance ( other , Data2 )          and \
               self.chain1_name == other.chain1_name and \
               self.chain2_name == other.chain2_name

    # =========================================================================
    ## get DataFrame for the first chain
    #  @see ROOT::RDataFrame
    @property
    def frame1 ( self ) :
        """'frame1': Get the DataFrame for the chain (same as 'frame')
        """
        return self.frame
    
    # =========================================================================
    ## get DataFrame for the chain
    #  @see ROOT::RDataFrame
    @property
    def frame2 ( self ) :
        """'frame2': Get the DataFrame for the second chain"""
        from   ostap.frames.frames import DataFrame, frame_progress 
        f = DataFrame ( self.chain2  )
        if not self.silent :
            pb = frame_progres ( f , len ( self.chain2 ) ) 
        return f
    
# =============================================================================

# =============================================================================
if '__main__' == __name__ :

    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
