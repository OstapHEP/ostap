#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/io/rootfiles.py
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'RootFiles'    , ## utility to collect and deal wth ROOT files 
) 
# =============================================================================
from   ostap.io.files     import Files
from   ostap.io.root_file import copy_file, copy_root_file 
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.io.rootfiles' )
else                       : logger = getLogger( __name__ )
# =============================================================================
## @class RootFiles
#  Some specialzation of gheneric class Files for ROOT files
class RootFiles(Files) :
    """ Specialzation of gheneric class Files for ROOT files
    - it is aware of the ROOT protocols as file names. 
    """
    ## protocols for remote (ROOT) files 
    root_protocols = (
        'root:'  ,
        'xroot:' ,
        'http:'  ,
        'https:' ,    
    )
    # =========================================================================
    ## Has ROOT protocol in th efile name? 
    @staticmethod 
    def has_protocol ( fname ) :
        """ Has ROOT protocol in th efile name?
        """
        for p in RootFiles.root_protocols :
            if p and fname.startswith ( p ) : return p
            return False
        
    # =========================================================================
    ## Strip protocol from the file name    
    @staticmethod 
    def strip_protocol ( fname ) :
        """ Strip ROOT protocol from the file name"""
        for p in RootFiles.root_protocols :
            if p and fname.startswith ( p ) : return fname.replace ( p , '' ) 
        return fname

    # =========================================================================
    ## Get the list of files from the patterns 
    def the_files ( self ) :
        """ Get the list of files from the patterns"""

        # 1. get the regular files
        files = set ( Files.the_files ( self ) )

        # 2. add explicit ROOT patterns 
        for pattern in  self.patterns :
            if RootFiles.has_protocol ( pattern ) :
                files.add ( pattern )
        
        return tuple ( sorted ( files ) ) 

    # =========================================================================
    ## get a common path (prefix) for all files in collection
    #  - protocols are ignored 
    @property 
    def commonpath ( self ) :
        """'commonpath': common path (prefix) for all files in collection
        - protocols are ignored 
        """
        ##
        if any ( RootFiles.has_protocol ( f ) for f in self.files ) :
            files = []
            for f in self.files :
                files .append ( RootFiles.strip_protocol ( f ) ) 
        else : files = self.files
        
        if not files : return '' 

        from ostap.utils.basic import commonpath 
        cp = commonpath ( os.path.abspath ( f ) for f in files )
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
        
        ## use the temporary directory
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

        cp = self.commonpath
        
        files_to_copy = set ( self.files )
        if also_bad and self.bad_files :
            files_to_copy |= set ( self.bad_files )

        for_copy = [] 
        for f in files_to_copy :
            fs   = os.path.abspath ( RootFiles.strip_protocol ( f ) ) 
            nf   = fs.replace ( cp , nd ) 
            nf   = os.path.abspath ( nf )
            pair = f , nf
            for_copy.append ( pair )

        nfiles = len ( for_copy ) 

        ## parallel copy 
        from   ostap.utils.basic      import numcpu
        if parallel and 1 < nfiles and 1 < numcpu () :
            
            from   ostap.utils.utils  import which
            if which ( 'parallel' ) : from ostap.utils.parallel_copy    import copy_files as parallel_copy
            else                    : from ostap.parallel.parallel_copy import copy_files as parallel_copy
            
            copied = parallel_copy ( for_copy , maxfiles = 1 , copier = copy_root_file , progress = not self.silent )
            copied = [ f [ 1 ] for f in copied ]

        ## sequential copy
        else :
            
            copied = []
            from ostap.utils.progress_bar import progress_bar
            for f, nf in progress_bar ( for_copy , silent = self.silent or nfiles <=1 ) :
                copied.append ( copy_root_file ( f , nf , progress = ( 1 == nfiles ) and self.verbose ) ) 
                
        if not self.silent :
            logger.info ( "copy_files: #%d files are copied to '%s'" %  ( len ( copied ) , nd ) )
            
        return self.clone ( files = copied )

    # ==========================================================================
    ## Get the size of the file
    #  @param fname input file name
    #  @return ROOT file sizeor -1 for non-existingt/invalid files
    def get_file_size ( self , fname ) :
        """ Get the size of the file
        - fname : input file name
        - return file size or -1 for non-existing/invalid files
        """
        
        ## (1) try as a regular file: 
        fsize = Files.get_file_size ( self , fname )
        if 0 <= fsize : return fsize                      ## RETURN

        ## (2) try to open it as ROOT file
        # ====================================================================+
        try : # ==============================================================+
            # ================================================================+
            with ROOT.TFile.Open ( fname , 'r' , exception = True ) as r :
                return r.GetSize()                        ## RETURN
            # ================================================================+
        except ( OSError , IOError ) : # =====================================+
            # ================================================================+
            pass

        return -1

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
    def hadd ( self , output = None , opts = "-ff -O" , **kwargs ) :
        """ Merge all files using <code>hadd</code> script from ROOT
        - `output`  name of the output merged file
        - `opts`   options for command <code>hadd</code>
        It returns the name of the merged file
        
        If no output file name is specified, the temporary name
        will be generate and the temporary file will be deleted
        at the end of the session

        OPTIONS:
        # -a                         Append to the output
        # -k                         Skip corrupt or non-existent files, do not exit
        # -T                         Do not merge Trees
        # -O                         Re-optimize basket size when merging TTree
        # -v                         Explicitly set the verbosity level: 
        #                            0 request no output, 
        #                            99 is the default
        # -j                         Parallelize the execution in multiple processes
        # -dbg                       Parallelize the execution in multiple processes in debug mode 
        #                            (Does not delete partial files stored inside working directory)
        # -d                         Carry out the partial multiprocess execution 
        #                            in the specified directory
        # -n                         Open at most 'maxopenedfiles' at once 
        #                            (use 0 to request to use the system maximum)
        # -cachesize                 Resize the prefetching cache use to speed up 
        #                            I/O operations(use 0 to disable)
        # -experimental-io-features  Used with an argument provided, enables the corresponding 
        #                            experimental feature for output trees
        # -f                         Gives the ability to specify the compression level of 
        #                            the target file(by default 4) 
        # -fk                        Sets the target file to contain the baskets with the same compression
        #                            as the input files (unless -O is specified). Compresses the meta data
        #                            using the compression level specified in the first input or the
        #                            compression setting after fk (for example 206 when using -fk206)
        # -ff                        The compression level use is the one specified in the first input
        # -f0                        Do not compress the target file
        # -f6                        Use compression level 6. 
        #                            (See TFile::SetCompressionSettings for the support range of value.)                            
        """
        from ostap.utils.utils import hadd as hadd_
        return hadd_ ( self.files , output = output , opts = opts  , **kwargs )

    # =========================================================================
    ## Split data into severals chunks of smaller size and merge (using <code>hadd</code>) each chunk
    #  @code
    #  data   = ..
    #  merged = data.merge_chunks ( 10 ) 
    #  @endcode 
    def __merge_chunks (  self , chunks , opts = '-ff -O' , parallel = True ) :
        """ Split data into severals chunks of smaller size and merge (using <code>hadd</code>) each chunk        
        >>> data   = ..
        >>> merged = data.merge_chunks ( 10 ) 
        """

        import ostap.utils.cleanup as  CU
        cu     = CU.CleanUp()
        tmpdir = cu.tmpdir

        if parallel and chunks and 2 <= len ( chunks[0].files ) :
            
            from   ostap.parallel.parallel import WorkManager
            from   ostap.utils.utils       import hadd2       as hadd_
            
            if     '-j' in opts : opts = opts.replace('-j','')
            if not '-v' in opts : opts = '%s -v 0' % opts

            pargs  = [ ( c.files , None , tmpdir , opts ) for c in chunks ] 
            
            wm     = WorkManager ( silent = True , progress = not self.silent )
            
            merged = [ o for o in wm.iexecute ( hadd_ , pargs , progress = not self.silent ) ]

            return self.clone ( files = merged ) 
            
        if not '-j' in opts : opts = '%s -j' % opts

        output = []

        if self.silent and not '-v' in opts : opts = '%s -v 0' % opts
        for c in progress_bar ( chunks , silent = self.silent or not chunks ) :
            output.append (  c.hadd ( opts = opts , dir = tmpdir ) ) 
                    
        return self.clone ( output )

    # =========================================================================
    ## Split/merge data into severals chunks of smaller size and merge (using <code>hadd</code>) each chunk
    #  @code
    #  data   = ..
    #  merged = data.merge_chunks ( 10 ) 
    #  @endcode 
    def merge_chunks (  self , chunk_size , opts = '-ff -O' , parallel = True ) :
        """ Split data into severals chunks of smaller size and merge (using <code>hadd</code>) each chunk        
        >>> data   = ..
        >>> merged = data.merge_chunks ( 10 ) 
        """
        chunks = self.split_chunks ( chunk_size )
        return self.__merge_chunks ( chunks , opts = opts , parallel = parallel )
        
    # =========================================================================
    ## Split/merge data into several group and merge (using  <code>hadd</code>) each group
    #  @code
    #  data   = ..
    #  merged = data.merge_groups ( 10 ) 
    #  @endcode 
    def merge_groups (  self , groups , opts = '-ff -O' , parallel = True ) :
        """ Split data into severals groups and merge (using <code>hadd</code>) each group
        >>> data   = ..
        >>> merged = data.merge_groups ( 10 ) 
        """
        chunks = self.split_groups ( groups )  
        return self.__merge_chunks ( chunks , opts = opts , parallel = parallel )

# =============================================================================



# =============================================================================
##                                                                      The END 
# =============================================================================
