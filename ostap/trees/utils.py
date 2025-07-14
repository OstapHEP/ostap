#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/trees/utils.py
#  Utilities to make TTree/TChain suitable for multiprocessing:
#  - pickling
#  - unpickling 
#  - spltitting 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
""" Utilities to make TTree/TChain suitable for multiprocessing:
- pickling
- unpickling 
- spltitting 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'Chain' , ## wrapper for Chain 
    'Tree'  , ## wrapper for TTree
) 
# =============================================================================
from   itertools              import accumulate 
from   ostap.core.ostap_types import string_types, integer_types 
from   ostap.utils.cleanup    import CleanUp
from   ostap.core.core        import valid_pointer, rootException
from   ostap.math.base        import FIRST_ENTRY, LAST_ENTRY, evt_range, all_entries
from   ostap.utils.basic      import file_info
from   ostap.utils.utils      import split_range 
import ostap.trees.base   
import ROOT, os 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.trees.utils' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Utilities to make TTree/TChain suitable for multiptocessing' )
# =============================================================================

## @class Chain
#  Simple class to make TChain suitable for multiprocessing:
# - pickling
# - unpickling 
# - spltitting
class Chain(CleanUp) :
    """ Simple class to make TChain suitable for multiprocesisng:
    - pickling
    - unpickling 
    - spltitting
    """
    # =========================================================================
    ##  Chain -> pickled state 
    def __getstate__  ( self ) :
        """ Chain -> pickled state 
        """
        ## full file name 
        def fullfn ( f ) :
            if os.path.exists ( f ) and os.path.isfile ( f ) :
                ff = os.path.abspath ( f  )
                if os.path.samefile  ( ff , f ) : return ff
            return f

        ## 
        file_infos = tuple ( ( fullfn ( f ) , file_info ( f ) ) for f in self.files )         
        ## the the host name 
        import socket 
        host = socket.getfqdn ().lower()        
        ## 
        return { 'name'     : self.__name  ,
                 'first'    : self.__first ,            
                 'last'     : self.__last  ,            
                 'files'    : file_infos   ,
                 'host'     : host         }
    # =========================================================================
    ## pickled state -> Chain 
    def __setstate__  ( self , state ) :
        """ Pickled state -> Chain 
        """
        
        self.__name    = state [ 'name'    ]
        self.__first   = state [ 'first'   ]
        self.__last    = state [ 'last'    ]
        #
        origin         = state [ 'host'    ]
        file_infos     = state [ 'files'   ]
        ##
        import socket 
        host    = socket.getfqdn ().lower()
        ##
        files = [] 
        
        same_host = origin == host
        if same_host : files = [ f [ 0 ] for f in file_infos ]
        else         :
            for fname , finfo in file_infos :
                if file_info ( fname )  == finfo :
                    ## hosts are different but the files are the same (shared file system?)
                    files.append ( fname )
                else :
                    # =========================================================
                    ## the file needs to be copied locally
                    ## @todo implement here  the parallel copy? 
                    from ostap.utils.scp_copy import scp_copy
                    full_name  = '%s:%s' % ( origin , fname ) 
                    copied , t = scp_copy  ( full_name )
                    if copied :
                        cinfo = file_info ( copied )
                        c = '%s --> %s' % ( full_name , "%s:%s" % ( host , copied ) )
                        if cinfo [ : 4 ] == finfo [ : 4 ] :
                            size = cinfo [ 1 ]
                            s = size / float ( 1024 ) / 1025 ##  MB 
                            v = s    / t                     ## MB/s 
                            logger.debug ( 'Chain.__setstate__ : file copied %s :  %.3f[MB] %.2f[s] %.3f[MB/s]' % ( c , s , t , v ) )
                            files.append ( copied ) ## APPEND 
                        else :
                            logger.error ( 'Chain.__setstate__ : Something wrong with the copy %s : %s vs %s ' % ( c , cinfo[:4] , finfo[:4] ) )
                    else : logger.error  ( "Chain.__setstate__: Failure to scp_copy the file: %s"  % full_name )

        if len ( files ) != len ( file_infos ) :
            logger.error ( "Chain.__setstate__: Some files are missing!" )
            
        self.__files = tuple ( files )

    # ======================================================================================
    ## create Chain object:
    #  - either from the real TTree/TChain:
    #  @code
    #  chain = ...  ## ROOT.TChain
    #  ch = Chain ( chain ) 
    #  @endcode
    #  - or from description:
    #  @code
    #  ch = Chain ( name = 'n', files = [ ... ] )
    #  @endcode 
    def __init__ ( self                  ,  
                   tree    = None        ,
                   name    = None        ,
                   files   =  []         ,
                   first   = FIRST_ENTRY ,
                   last    = LAST_ENTRY  ) :         
        """ Create Chain object 
        
        - either from the real TTree/TChain:
        
        >>> chain = ...  ## ROOT.TChain
        >>> ch = Chain ( chain )
        
        - or from description:
        
        >>> ch = Chain ( name = 'n', files = [ ... ] )
        """
        
        if   isinstance ( files , string_types ) : files = files ,

        ## 3 valid cases: 
        assert ( isinstance ( tree , Chain ) and tree.chain ) or \
            ( name and files                                ) or \
            ( isinstance ( tree , ROOT.TTree                ) and valid_pointer ( tree ) ) , \
            "Invalid tree/name/files combination: %s/%s%s" % ( tree , name , files    )
        
        ## copy-like, ignore other arguments  
        if isinstance  ( tree , Chain ) and tree.chain :

            if name  and name != tree.name : logger.warning ( 'Chain: explicitley specified name is inored!' )
            if files and set ( files ) != set ( tree.files ) :
                logger.warning ( 'Chain: explicitely specified files are inored!' )
                
            self.__name  = tree.name 
            self.__files = tree.files
            self.__last  = tree.last 
            
        elif name and files :
            
            self.__name    = name  
            self.__files   = files
            self.__first   = first 
            self.__last    = last 
            
            ## check that object exists for each file:
            for i , f in enumerate ( self.__files ) :
                with ROOT.TFile.Open ( f , 'read' , exception = False ) as rf : 
                    assert rf                , "Tree: TFile (%s) is invalid or non-existing " % f 
                    assert self.__name in rf , "Tree: TFile (%s) ihas no `%s` object"         % ( f , self.__name  ) 

            ## reconstruct the chain 
            chain = self.chain
            assert chain , "Chain:  invalid reconstructed TChain!"
            
        elif isinstance ( tree , ROOT.TTree  ) and valid_pointer ( tree ) :
            
            ## ATTENTION: full path is specified here! 
            self.__name = tree.full_path
            
            ## get the files from Tree/TChain 
            files = tree.files () ## get the files from Tree/TChain 
            
            if files : self.__files = tuple ( files ) 
            else     : 
                ## Memory resident tree? put it into the temporary file. keep ?
                tmpfile = Cleanup.tmpfile ( suffix = '.root' , prefix = 'ostap-tree-' , keep = True ) 
                logger.info ( "TTree is memory resident! Write it into the temporary TFile:%s" % tmpfile )
                with ROOT.TFile ( tmpfile , 'NEW' ) as rf : rf [ self.name ] = tree    
                self.__files = tmpfile ,

            self.__first   = first
            self.__last    = last 
                
        ## reconstruct the chain
        ch = self.chain 
        assert ch , "Chain: invalid reconstructed TChain!"

        ## check & adjust first/last setting 
        first , last = evt_range ( ch , first , last ) 
        assert 0 <= first <= last , "Chain: invalid first/last setting!"
        
        self.__first = first
        self.__last  = last 

        nfiles  = self.nFiles
        nevents = len ( ch )
        
        ## remove the files that are outside of the range 
        if 2 <= nfiles and ( 0 < self.first or self.last < nevents ) :

            first, last = self.first, self.last 
            files       = list ( self.files )
            
            asizes = [ s for s in self.accumulated_sizes() ]
            
            while 2 <= len ( asizes ) and last < asizes [ -2 ]  :
                asizes.pop ()
                files .pop ()

            while 2 <= len ( asizes ) and asizes [ 0 ] <= first :
                first -= asizes [ 0 ]
                last  -= asizes [ 0 ]
                asizes.pop ( 0 )
                files .pop ( 0 )
                
            ## redefine files, first & last
                
            self.__nfiles = tuple ( files )
            self.__first  = first
            self.__last   = last 
            
    # =========================================================================
    ## Generator to Get the sizes for all files
    #  @code
    #  ch = Chain ( ... )
    #  for s in ch.sizes () : 
    #  @endcode 
    def sizes ( self ) :
        """ Generator to get the sizes for all files
        >>> ch = Chain ( ... )
        >>> for s in ch.sizes () : 
        """
        for f in self.__files :
            ch = ROOT.TChain ( self.name , f )
            yield len ( ch )
    # ===========================================================================
    ## Generator to get accumulated sized for all files
    #  @code
    #  ch = Chain ( ... )
    #  for s in ch.accumulated_sizes () : 
    #  @endcode 
    def accumulated_sizes ( self ) :
        """ Generator to get accumulated sized for all files
        >>> ch = Chain ( ... )
        >>> for s in ch.accumulated_sizes () : 
        """
        for ss in accumulate ( self.sizes () ) : yield ss
        
    # =========================================================================
    ## *GENERATOR* to split the chain for several chains with at most
    #  chunk_size entries and at most max_files
    #  @code
    #  chain = Chain ( ... )
    #  for i in chain.split ( chunk_size = 1000 ) :
    #  ... 
    #  @endcode 
    def split ( self , chunk_size = -1 , max_files = -1  ) :
        """ GENERATOR Split the tree for several trees 
        with at most  `chunk_size` entries
        >>> chain  = Chain ( .... ) 
        >>> for t in tree.split ( chunk_size = 1000000 )  : 
        >>> ...
        """
        assert isinstance ( chunk_size , integer_types ) , \
            'Chain.split : invalid chunk_size %s' % chunk_size

        if not self.nFiles  : return                                   ## RETURN 
        
        if 1 == self.nFiles :
            tree = Tree ( name  = self.name         ,
                          file  = self.files [ 0 ]  ,
                          first = self.first        ,  
                          last  = self.last         )
            for t in tree.split ( chunk_size = chunk_size ) : yield t  ## YIELD
            return                                                     ## RETURN 
        
        asizes = tuple ( s for s in self.accumulated_sizes()  )

        ## the first tree: can be incomplete 
        first  = self.first
        tree   = Tree ( name = self.name , file = self.files [ 0 ] , first = first )
        for t in tree.split ( chunk_size = chunk_size ) : yield t      ## YIELD  
        
        ## full trees (all of them are complete 
        for f in self.files [1:-1] :
            tree = Tree ( name = self.name , file = f  )
            for t in tree.split ( chunk_size = chunk_size ) : yield t  ## YIELD 
            
        ## the last tree: can be incomplete
        last   = self.last - asizes [ -2 ] ## ATTENTION HERE!
        tree   = Tree ( name = self.name , file = self.files[-1] , last = last  )
        for t in tree.split ( chunk_size = chunk_size ) : yield t      ## YIELD 

    @property
    def chain ( self ) :
        """`chain' : get the underlying tree/chain"""
        chain = ROOT.TChain ( self.name )
        for f in self.__files : chain.Add ( f ) 
        return chain 

    @property
    def name    ( self ) :
        """`name'   : TTree/TChain name"""
        return self.__name
    
    @property
    def files   ( self ) :
        """`files'   : the files"""
        return self.__files
    
    @property    
    def nFiles   ( self ) :
        """`nFiles'   : the numer of files"""
        return len(self.__files)
    
    @property    
    def first   ( self ) :
        """`first' : the first event to process"""
        return self.__first
    
    @property
    def last    ( self ) :
        """`last' : the last event (not-inclusive)"""
        return self.__last
    
    @property
    def nevents ( self ) :
        """`nevents' : number of events to process ( == last - first)"""
        return self.__last - self.first 

    # ============================================================================
    ## get DataFrame
    #  @code
    #  tree = ....
    #  f  =  tree.frame () ## get
    #  f1 =  tree.frame ('px', 'py' , 'pz') ## get frame with default branches
    #  f2 =  tree.frame ( tx = 'px/pz' , ty = 'py/pz') ## define new variables
    #  @endcode 
    def  frame ( self , *vars , **newvars ) :
        """`frame' : get ROOT.RDataFrame for the given chain/tree
        >>> tree = ....
        >>> f  = tree.frame () ## get
        >>> f1 = tree.frame ('px', 'py' , 'pz') ## get frame with default branches
        >>> f2 = tree.frame ( tx = 'px/pz' , ty = 'py/pz') ## define new variables 
        """
        chain  = self.chain
        assert all_entries ( chain , self.first , self.last ), \
            "Chain.frame: can not use RDataFrame due to first/last setting!"
        
        from ostap.frames.frames import DataFrame
        from ostap.core.core     import strings    
        fnames = strings    ( *self.files )
        vnames = strings    ( *vars       )
        df     = DataFrame  ( self.name , fnames , vnames )
        for k in new_vars : df =  df.Define ( k , new_vars  [ k ] )
        return  df                              

    # =========================================================================
    ## delegate all other attributes to the underlying chain object 
    def __getattr__  ( self , attr ) :
        """ Delegate all other attributes to the underlying chain object"""
        return getattr  ( self.chain , attr )

# =============================================================================
## @class Tree
#  Simple class to make TTree suitable for multiprcessing: 
# - pickling
# - unpickling 
# - spltitting
class Tree(Chain) :
    """ Simple class to make TTree suitable for multiprocessing:
    - pickling
    - unpickling 
    - spltitting
    """
    ## pickle 
    def __getstate__  ( self )         : return Chain.__getstate__  ( self )
    ## unpickle 
    def __setstate__  ( self , state ) :        Chain.__setstate__  ( self , state ) 
    ## constructor 
    def __init__      ( self                   ,  
                        tree    =  None        ,
                        name    =  None        ,
                        file    =  ''          ,
                        first   =  FIRST_ENTRY ,
                        last    =  LAST_ENTRY  ) :
        
        if name and file :            
            assert isinstance ( file , string_types ), 'Tree: `file` should be single file name!'
            
        elif isinstance ( tree , Chain  )  : 
            assert 1 ==  tree.nFiles , "Tree: only single files are OK!"
             
        else :
            assert isinstance    ( tree , ROOT.TTree ) , 'Tree: tree  is not TTree!'
            assert valid_pointer ( tree )              , "Tree: TTree is invalid!" 
            if isinstance ( tree , ROOT.TChain ) :
                assert 1 == tree.nFiles() , 'Tree: only single file is allowed!'
              
        Chain.__init__ ( self , tree = tree , name = name , files = [ file ] , first = first , last = last  )
        
        assert 1 == self.nFiles , 'Invalid number of files!'

    # ==========================================================================
    ## *GENERATOR* split the tree into smaller trees
    #  @code
    #  tree = Tree  ( ... )
    #  for i in tree.split ( chunk_size = 1000 ) :
    #  ... 
    #  @endcode 
    def split ( self , chunk_size = -1 , max_files = -1  ) :
        """ *GENERATOR* split the tree into smaller Tree objects 
        with at most  `chunk_size` entries
        >>> tree  = Tree ( .... ) 
        >>> for t in tree.split ( chunk_size = 1000000 )  : 
        >>> ...
        """        
        assert isinstance ( chunk_size , integer_types ) , \
            'Tree.split: invalid chunk_size %s' % chunk_size

        ## no split 
        if chunk_size <= 0 :
            yield self                     ## YIELD            
            return                         ## RETUN
        
        the_file = self.file
        for first, last in split_range ( self.first , self.last , chunk_size ) :
            yield Tree ( name  = self.name ,
                         file  = the_file  ,
                         first = first     ,
                         last  = last      )
    @property
    def file ( self ) :
        """`file'   : the file name """
        fs = self.files 
        assert 1 == len  ( fs ) , 'Tree.file: Invalid number of files %s' % len ( fs ) 
        return fs [ 0 ] 
    
    # =========================================================================
    ## delegate all other attributes to the underlying chain object 
    def __getattr__  ( self , attr ) :
        """ Delegate all other attributes to the underlying tree/chain object"""
        return getattr  ( self.chain , attr )
    
    @property
    def tree ( self ) :
        """`tree' : get the underlying TTree/TChain"""
        return self.chain

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END
# =============================================================================
