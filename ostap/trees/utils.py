#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/trees/utils.py
#  Simple utilities for TTree/TChain treatment in multiprocesisng 
#  - pickling
#  - unpickling 
#  - spltitting 
#  machinery for splitting
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
""" Simple utilities for TTree/TChain treatment in multiprocesisng 
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
    'Ttree' , ## wrapper for TTree
) 
# =============================================================================
from   itertools              import accumulate 
from   ostap.core.ostap_types import string_types 
from   ostap.utils.cleanup    import CleanUp
from   ostap.core.core        import valid_pointer, rootException
from   ostap.math.base        import FIRST_ENTRY, LAST_ENTRY, evt_range, all_entries
from   ostap.utils.basic      import file_info
import ostap.trees.base   
import ROOT, os 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.trees.utils' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Utilities for TTree/TChain pickling/unpickling and splitting' )
# =============================================================================
## @class Chain
#  Simple class to make TChain "pickleable"
#  It is needed for multiprcessing:
# - pickling
# - unpickling 
# - spltitting
class Chain(CleanUp) :
    """ Simple class to make TTree/TChain "pickleable 
    It is needed mainly for multiprocessing: 
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
        #
        ##
        import socket 
        host    = socket.getfqdn ().lower()
        ##
        self.__chain   = None
        self.__files   = [] 
        
        for fname , finfo in file_infos :
            if   same_host : self.__files.append ( fname )
            else :
                fnew = file_info ( fname )
                if fnew == finfo :
                    ## hosts are different but the files are the same (shared file system?)
                    self.__files.append ( fname )
                else :
                    # =========================================================
                    ## the file needs to be copied locally
                    ## @todo implement the parallel copy! 
                    from ostap.utils.scp_copy import scp_copy
                    full_name  = '%s:%s' % ( origin , fname ) 
                    copied , t = scp_copy  ( full_name )
                    if copied :
                        cinfo = file_info ( copied )
                        c = '%s -> %s' % ( full_name , "%s:%s" % ( self.__host , copied ) )
                        if cinfo[:4] == finfo [:4] :
                            size = cinfo[1]
                            s =  cinfo [ 1 ] / float ( 1024 ) / 1025 ##  MB 
                            v = s / t                  ## MB/s 
                            logger.debug ( 'File copied %s :  %.3f[MB] %.2f[s] %.3f[MB/s]' % ( c , s , t , v ) )
                            self.__files..append ( copied ) ## APPEND 
                        else :
                            logger.error ( 'Something wrong with the copy %s : %s vs %s ' % ( c , cinfo[:4] , finfo[:4] ) )
                    else : logger.error ("Cannot copy the file %s"  % full_name )

        self.__files = tuple ( self.__files )

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

        assert ( isinstance ( tree , Chain ) and tree.chain ) or \
            ( name and files                                ) or \
            ( isinstance ( tree , ROOT.TTree                ) and valid_pointer ( tree ) ) , \
            "Invalid tree/name/files combination: %s/%s%s" % ( tree , name , files    )
        
        ## copy-like, ignore other arguments  
        if isinstance  ( tree , Chain ) and tree.chain :

            if name  and name != tree.name : logger.warning ( 'Chain: explicitley specified name is inored!' )
            if files and set ( files ) != set ( ttee.files ) :
                logger.warning ( 'Chain: explicitely specified files are inored!' )
                
            self.__name  = tree.name 
            self.__files = tree.files
            self.__last  = tree.last 
            
        elif name and files :
            
            self.__name    = tree.name 
            self.__files   = tree.files
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
            
            input_chain = tree
            
            ## ATTENTION: full path is specified here! 
            self.__name = tree.full_path
            
            ## get the files from Tree/TChain 
            files = tree.files () ## get the files from Tree/TChain 
            
            if files : self.__files = tuple ( files ) 
            else     : 
                ## Memory resident tree? put it into the temporary file. keep ?
                tmpfile = Cleanup.tmpfile ( suffix = '.root' , prefix = 'ostap-tree-' , keep = True ) 
                logger.info ( "TTree is memory resident! Write in into the temporary TFile:%s" % tmpfile )
                with ROOT.TFile ( tmpfile , 'NEW' as rf ) : rf [ self.name ] = tree    
                self.__files = tmpfile ,

            self.__first   = first
            self.__last    = last 
                
        ## reconstruct the chain
        ch = self.chain 
        assert ch , "Chain: invalid reconstructed TChain!"

        ## check & adjust first/last setting 
        first , last = evt_range ( ch , first , last ) 
        assert 0 <= first <= last , "Chain: invalid first/last setting!" )
        
        self.__first = first
        self.__last  = last 

        nfiles  = self.nFiles
        nevents = len ( ch )
        
        ## remove the files that are outsize of the range 
        if 2 <= self.nFiles and ( 0 < self.first or self.last < nevents ) :

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
    ## GENERATOR to split the chain for several chains with at most chunk_size entries and at most max_files 
    def split ( self , chunk_size = -1 , max_files = -1  ) :
        """ GENERATOR Split the tree for several trees with chunk_size entries
        >>> tree = ....
        >>> trees = tree.split ( chunk_size = 1000000 ) 
        """        
        if chunk_size <= 0 : chunk_size = 100000
        if max_files  <= 0 : max_files  = 1 

        assert isinstance ( chunk_size , integer_types ) and 0 < chunk_size , \
            'Chain.split : invalid chunk_size %s' % chunk_size
        
        if 1 == self.nFiles :
            tree = Tree ( self.name                 ,
                          file  = self.files [ 0 ]  ,
                          first = self.first        ,  
                          last  = self.last         )
            for t in tree.split ( chunk_size = chunk_size ) : yield t
            
        asizes = tuple ( s for s in self.accumulated_sizes()  )
        first  = self.first
        last   = self.last = asizes [ -2 ] ## ATTENTION HERE!
        
        ## the first tree, can be incomplete 
        tree = Tree ( self.name , file = self.files[0 ] , first = first )
        for t in tree.split ( chunk_size = chubnk_size ) : yield t
        
        ## full trees (all of them are complete 
        for f in self.files [1:-1] :
            tree = Tree ( self.name , file = f  )
            for t in tree.split ( chunk_size = chunk_size ) : yield t
            
        ## the last tree (can be incomplete) 
        tree = Tree ( self.name , files = self.files[-1] , last = last  )
        for t in tree.split ( chunk_size = chunk_size ) : yield t 

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
#  Simple class to make TTree "pickleable"
#  It is needed for multiprcessing: 
# - pickling
# - unpickling 
# - spltitting
class Tree(Chain) :
    """ Simple class to make TTree "pickleable"
    It is needed for multiprcessing: 
    - pickling
    - unpickling 
    - spltitting
    """
    ## pickle 
    def __getstate__  ( self )         : return Chain.__getstate__  ( self )
    ## unpickle 
    def __setstate__  ( self , state ) :        Chain.__setstate__  ( self , state ) 
    ## constructir 
    def __init__      ( self                   ,
                        tree    =  None        ,
                        name    =  None        ,
                        file    =  ''          ,
                        first   =  FIRST_ENTRY ,
                        last    =  LAST_ENTRY  ) :
        
        if name and file :            
            assert isinstance ( file , string_type  ), 'Tree: `file` should be single file name!'
            
        elif valid_pointer    ( tree ) :
            assert isinstance ( tree , ROOT.TTree ) , 'Tree: tree is not TTree!'
            
            if isinstance ( tree , ROOT.TChain ) :
                assert 1 == tree.nFiles() , 'Tree: only single file is allowed!'
                
        Chain.__init__ ( self , tree , name  , files = [ file ] , first = first , last = last  )
        
        assert 1 == self.nFiles , 'Invalid number of files!'

    # ==========================================================================
    ## Generator to split the tree into smaller trees
    def split ( self , chunk_size = -1 , max_files = -1  ) :
        """ GENERATOR to split the tree into smaller Tree objects 
        """
        if chunk_size <= 0 : chunk_size = 1000000

        assert isinstance ( chunk_size , integer_types ) and 0 < chunk_size , \
            'Tree.split : invalid chunk_size %s' % chunk_size

        ## small enough ?
        if self.last - self.first < chunk_size ) :  yield self,
        
        ## split it!
        the_file = self.file
        for first, last in split_range ( self.first , self.last , chubnk_size ) :
            yield Tree ( name  = self.name ,
                         file  = the_file  ,
                         first = first     ,
                         last  = last      )
            
    @property
    def file ( self ) :
        """`file'   : the file name """
        fs = self.files 
        assert 1 == len  ( fs ) , 'Tree: Invalid number of files %s' % len ( fs ) 
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
