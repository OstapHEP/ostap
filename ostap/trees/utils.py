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
from   collections            import namedtuple 
from   itertools              import accumulate 
from   ostap.core.ostap_types import string_types, integer_types, path_types  
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
### helper type to represenrt the file/tree  with correct/valid  number of entries 
FileItem        = namedtuple ( 'FileItem' , ( 'file_name' , 'size' , 'hash' ) )
file_item_types = path_types + ( FileItem , )
# =============================================================================
## full file name 
def fullfn ( f ) :
    if os.path.exists ( f ) and os.path.isfile ( f ) :
        ff = os.path.abspath ( f  )
        if os.path.samefile  ( ff , f ) : return ff
    return f
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
        ## 
        file_infos = tuple ( ( item , file_info ( item.file_name ) ) for item in self.__items )         
        ## the the host name 
        import socket 
        host = socket.getfqdn ().lower()        
        ## 
        return { 'name'     : self.__name  ,
                 'first'    : self.__first ,            
                 'last'     : self.__last  ,            
                 'items'    : file_infos   ,
                 'host'     : host         }
    # =========================================================================
    ## pickled state -> Chain 
    def __setstate__  ( self , state ) :
        """ Pickled state -> Chain 
        """
        
        self.__name    = state [ 'name'  ]
        self.__first   = state [ 'first' ]
        self.__last    = state [ 'last'  ]
        #
        origin         = state [ 'host'  ]
        file_infos     = state [ 'items' ]
        ##
        import socket 
        host    = socket.getfqdn ().lower()
        ##
        items   = [] 
        
        same_host = origin == host
        if same_host : items = [ f [ 0 ] for f in file_infos ]
        else         :
            for item , finfo in file_infos :
                if file_info ( item.file_name ) == finfo :
                    ## hosts are different but the files are the same (shared file system?)
                    item.append ( entry )
                else :
                    # =========================================================
                    ## the file needs to be copied locally
                    ## @todo implement here  the parallel copy? 
                    from ostap.utils.scp_copy import scp_copy
                    full_name  = '%s:%s' % ( origin , entry.file_name ) 
                    copied , t = scp_copy  ( full_name )
                    if copied :
                        cinfo = file_info ( copied )
                        c = '%s --> %s' % ( full_name , "%s:%s" % ( host , copied ) )
                        if cinfo [ : 4 ] == finfo [ : 4 ] :
                            size = cinfo [ 1 ]
                            s = size / float ( 1024 ) / 1025 ##  MB 
                            v = s    / t                     ## MB/s 
                            logger.debug ( 'Chain.__setstate__ : file copied %s :  %.3f[MB] %.2f[s] %.3f[MB/s]' % ( c , s , t , v ) )
                            ## 
                            items.append ( self.__mke_item (  copied ) ) ## APPEND 
                        else :
                            logger.error ( 'Chain.__setstate__ : Something wrong with the copy %s : %s vs %s ' % ( c , cinfo[:4] , finfo[:4] ) )
                    else : logger.error  ( "Chain.__setstate__: Failure to scp_copy the file: %s"  % full_name )

        if len ( items ) != len ( file_infos ) :
            logger.error ( "Chain.__setstate__: Some files are missing!" )
            
        self.__items = tuple ( items )

    # ======================================================================================
    def __make_item ( self , file_name ) :
        
        if isinstance ( file_name , string_types ) :

            file_name  = fullfn ( file_name ) 
            ch         = ROOT.TChain ( self.name )
            ch.Add ( file_name )
            size       = len ( ch )
            size_hash  = self.__make_hash ( self.name , file_name , size )
            # 
            return FileItem ( file_name , size , size_hash )
        
        elif isinstance ( file_name , path_types ) : return self.__make_item ( str ( file_name ) )
        elif isinstance ( file_name , FileItem   ) :
            return file_name if self.__good_item ( file_name ) else self.__make_item ( file_name.file_name )

        raise TypeError ( "Unknown/invalid `file_name` type %s" % typename ( file_name ) ) 

    # ======================================================================================
    def __make_hash  ( self , tree_name , file_name , size ) :
        return hash ( ( self.name , file_name , size ) ) 

    
    ## of this is a good/valid entry? 
    def __good_item  ( self , item  ) :
        return isinstance ( item , FileItem ) and \
            item.hash == self.__make_hash ( self.name , item.file_name , item.size  )
    
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
                   files   = []          ,
                   first   = FIRST_ENTRY ,
                   last    = LAST_ENTRY  ) :         
        """ Create Chain object 
        
        - either from the real TTree/TChain:
        
        >>> chain = ...  ## ROOT.TChain
        >>> ch = Chain ( chain )
        
        - or from description:
        
        >>> ch = Chain ( name = 'n', files = [ ... ] )
        """

        if   isinstance ( files , file_item_types ) : files = files ,
        
        ## "copy" constructor 
        if isinstance ( tree , Chain ) :
            
            if name  and name  != tree.name  : logger.warning ( "explicitly specified name `%s` is ignored!" % name  )
            if files :
                if   files == tree.items : pass
                elif files == tree.files : pass
                else : logger.warning ( "explicitly  specified list of `files` is ignored!" )
                
            tfirst , tlast = evt_range  ( sum ( s for s in tree.sizes () ) , first , last )
            if tfirst != tree.first or tlast != tree.last : logger.warning ( "explicitly specified `first/last` %s/%s are ignored!"  % ( first , last ) )
            
            self.__name  = tree.name
            self.__items = tree.items            
            self.__first = tree.first
            self.__last  = tree.last  

        elif isinstance ( tree , ROOT.TTree  ) and valid_pointer ( tree ) and 1 <= tree.nFiles : 
            
            self.__name  = tree.fullpath
            self.__items = tuple ( self.__make_item ( fname ) for fname in tree.files ) 
            self.__first = first
            self.__last  = last 
                    
        elif isinstance ( tree , ROOT.TTree ) and valid_pointer ( tree ) : 

            self.__name  = tree.fullpath
            
            ## Memory resident tree? put it into the temporary file. keep ?
            tmpfile = Cleanup.tmpfile ( suffix = '.root' , prefix = 'ostap-tree-' , keep = True ) 
            logger.attention ( "TTree is a memory resident! Write it into the temporary TFile:%s" % tmpfile )
            with ROOT.TFile ( tmpfile , 'NEW' ) as rf : rf [ self.name ] = tree
            
            self.__items = self.__make_item ( tmpfile ) , 
            self.__first = first
            self.__last  = last 
            
        elif not tree is None :
            raise TypeError ( 'Unknown/invalid type of `tree` : %s' % typename ( tree ) )
        
        elif not name or not isinstance ( name , string_types ) :            
            raise TypeError ( 'Unknown/invalid type of `name` : %s' % typename ( name ) )
        
        elif files and all ( isinstance ( e , file_item_types ) for e in files ) :
            
            self.__name  = name
            self.__items = tuple ( self.__make_item  ( e ) for e in files ) 
            self.__first = first
            self.__last  = last 
            
        else :
            
            raise TypeError ( 'Inconsistent tree/name/path: %s/%s/%s' % ( typename ( tree  ) ,
                                                                          typename ( name  ) ,
                                                                          typename ( files ) ) )        

        
        asizes  = [ s for s in self.accumulated_sizes() ]
        nevents = asizes [ -1 ]
        
        ## check & adjust first/last setting 
        first , last = evt_range ( nevents  , first , last ) 
        assert 0 <= first <= last , "Chain: invalid first/last setting!"
        
        self.__first = first
        self.__last  = last
        
        ## remove the files that are outside of the range 
        if 2 <= self.nFiles and ( 0 < self.first or self.last < nevents ) :

            first, last = self.first, self.last 
            items       = list ( self.items )
                        
            while 2 <= len ( asizes ) and last < asizes [ -2 ]  :
                asizes.pop ()
                items.pop  ()

            while 2 <= len ( asizes ) and asizes [ 0 ] <= first :
                first -= asizes [ 0 ]
                last  -= asizes [ 0 ]
                asizes  .pop    ( 0 )
                items   .pop    ( 0 )
                
            ## redefine files, first & last
                
            self.__items  = tuple ( items  )
            self.__first  = first
            self.__last   = last 
            
    # =========================================================================
    ## Generator to get the sizes for all files
    #  @code
    #  ch = Chain ( ... )
    #  for s in ch.sizes () : 
    #  @endcode 
    def sizes ( self ) :
        """ Generator to get the sizes for all files
        >>> ch = Chain ( ... )
        >>> for s in ch.sizes () : 
        """
        for item in self.__items : yield item.size

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
            tree = Tree ( name  = self.name        ,
                          file  = self.items [ 0 ] ,
                          first = self.first       ,  
                          last  = self.last        )
            for t in tree.split ( chunk_size = chunk_size ) : yield t  ## YIELD
            return                                                     ## RETURN 
        
        asizes = tuple ( s for s in self.accumulated_sizes()  )

        ## the first tree: can be incomplete 
        first  = self.first
        tree   = Tree ( name = self.name , file = self.items [ 0 ] , first = first )
        for t in tree.split ( chunk_size = chunk_size ) : yield t      ## YIELD  
        
        ## full trees (all of them are complete 
        for entry in self.items [ 1 : -1 ] :
            tree = Tree ( name = self.name , file = entry  )
            for t in tree.split ( chunk_size = chunk_size ) : yield t  ## YIELD 
            
        ## the last tree: can be incomplete
        last   = self.last - asizes [ -2 ] ## ATTENTION HERE!
        tree   = Tree ( name = self.name , file = self.items [ -1 ] , last = last  )
        for t in tree.split ( chunk_size = chunk_size ) : yield t      ## YIELD 

    @property
    def chain ( self ) :
        """`chain' : get the underlying tree/chain"""
        chain = ROOT.TChain ( self.name )
        for item in self.items  : chain.Add ( item.file_name ) 
        return chain 

    @property
    def name    ( self ) :
        """`name'   : TTree/TChain name"""
        return self.__name
    
    @property
    def items  ( self ) :
        """`items'   : a tuple of items """
        return self.__items 
    
    @property
    def files   ( self ) :
        """`files'   : the files"""
        return tuple ( f.file_name  for f in self.items ) 
    
    @property    
    def nFiles   ( self ) :
        """`nFiles'   : the numer of files"""
        return len ( self.items )
    
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
                        tree    =  None        , * , 
                        name    =  None        ,
                        file    =  ''          ,
                        first   =  FIRST_ENTRY ,
                        last    =  LAST_ENTRY  ) :

        assert tree is None or isinstance ( tree , Tree ) or \
            ( isinstance ( tree , ROOT.TTree ) and valid_pointer ( tree ) ) , \
            "Tree: Unknown/invalid type for `tree`:%s" % typename ( tree )
        
        assert not file or isinstance ( file , file_item_types ) , \
            "Tree: `file` must be either None or single file/entry! %s" % typename ( file ) 
        
        Chain.__init__ ( self             ,
                         tree  = tree     ,
                         name  = name     ,
                         files = [ file ] ,
                         first = first    ,
                         last  = last     )
        
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
        
        for first, last in split_range ( self.first , self.last , chunk_size ) :
            yield Tree ( name  = self.name        ,
                         file  = self.items [ 0 ] ,
                         first = first            ,
                         last  = last             )
    @property
    def file ( self ) :
        """`file'   : the file name """
        assert 1 == self.nFiles , 'Tree.file: Invalid number of files %s' % self.nFiles 
        return self.items [ 0 ].file_name 
    
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
