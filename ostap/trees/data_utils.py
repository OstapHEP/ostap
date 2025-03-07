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
    'RootFiles'   , ## collect ROOT files  
    'Data'        , ## collect ROOT files and create TChain(s)
    )
# =============================================================================
from   collections            import defaultdict 
from   ostap.core.ostap_types import string_types 
from   ostap.core.core        import rootError, rootWarning
from   ostap.io.files         import Files, fsize_unit 
from   ostap.io.root_file     import RootFiles
import ROOT, os, math  
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.trees.data' )
else                       : logger = getLogger ( __name__     )
# =============================================================================
if not hasattr ( ROOT.TTree , '__len__' ) :  
    ROOT.TTree. __len__ = lambda s : s.GetEntries()

chain_types  = string_types + ( ROOT.TChain, ROOT.TTree ) 
# =============================================================================
## @class Data
#  Simple utility to access to certain chain in the set of ROOT-files
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @author Alexander BARANOV a.baranov@cern.ch
#  @date   2014-06-08  
class Data(RootFiles):
    """ Simple utility to access to certain chain in the set of ROOT-files
    >>> data  = Data('Bc/MyTree', '*.root' )
    >>> chain = data.chain
    >>> flist = data.files 
    """
    def __init__( self                 ,
                  files                , / , 
                  *chains              , 
                  description  = ''    , 
                  maxfiles     = -1    ,
                  check        = True  , 
                  silent       = False ,
                  parallel     = False ) : 

        ## we will need Ostap machinery for trees&chains here
        import ostap.trees.trees
        
        if   isinstance ( chains , ROOT.TTree   ) : chains = chains.path , 
        elif isinstance ( chains , string_types ) : chains = chains      ,

        assert chains and all ( isinstance ( c , chain_types ) for c in chains ) , \
            "Invald type of `chains` %s" % str ( chains )

        self.__check = True if check else False
        
        ## chain names
        self.__chain_names = tuple ( ( c.path if isinstance ( c , ROOT.TTree ) else c ) for c in chains ) 

        ## update description
        if not description :
            description = "Chains:{%s}" %  ( ', '.join ( self.__chain_names ) ) 

        self.__bad_files = defaultdict(set) 
        
        ## initialize the  base class 
        RootFiles.__init__( self                      ,
                            files                     ,
                            description = description ,
                            maxfiles    = maxfiles    ,
                            silent      = silent      ,
                            parallel    = parallel    )
        
    # =========================================================================    
    @property
    def validate ( self ) :
        """'check' : make check for `TTree`/`TChain` structures
        """
        return self.__check
    
    @property
    def chain_names ( self ) :
        """'chain_names' : the names of TTree/TChain objects
        """
        return self.__chain_names

    # ========================================================================
    ## get the chain by index  
    def get_chain ( self , index ) :
        """ Get the chain by index
        """
        if not 0 <= index < len ( self.__chain_names ) :
            raise IndexError ( "Index out the range!" )        
        cname     = self.chain_names [ index ] 
        ch        = ROOT.TChain ( cname )
        files     = self.files
        bad_files = self.bad_files[cname] 
        for f in files :
            if not f in bad_files : ch.Add ( f ) 
        return ch

    @property
    def nchains ( self ) :
        """`nchains` : number of chains 
        """
        return len ( self.chain_names )
    
    @property
    def chains ( self ) :
        """'chains' : (re)built and return all `TChain` objects"""
        result = []
        for index in range ( self.nchains ) : 
            result.append ( self.get_chain ( index ) ) 
        return tuple ( result )
        
    @property
    def chain ( self ) :
        """'chain' : (re)built and return the first `TChain` object , same as `chain1`"""
        return self.get_chain ( 0 )
    
    @property
    def chain1 ( self ) :
        """'chain1' : (re)built and return the first `TChain` object , same as `chain`"""
        return self.get_chain ( 0 )

    @property
    def chain2 ( self ) :
        """'chain2' : (re)built and return the second `TChain` object"""
        return self.get_chain ( 1 )
    
    @property
    def check ( self ) :
        """'check' : make check for `TTree`/`TChain` structures
        """
        return self.__check

    @property
    def bad_files ( self ) :
        """`bad_files` : dictionaryof bad files """
        return self.__bad_files
    
    ## check the content of the two trees
    @staticmethod 
    def check_trees ( tree1 , tree2 , the_file = '' ) :

        if tree1 and tree2 :
            
            branches1 = set ( tree1.branches () )                
            leaves1   = set ( tree1.leaves   () )
            
            branches2 = set ( tree2.branches () )                
            leaves2   = set ( tree2.leaves   () )

            if branches1 != branches2 :
                missing = sorted ( branches1 - branches2 ) 
                extra   = sorted ( branches2 - branches1 )
                if missing : logger.warning ( "TTree('%s'): %3d missing branches: {%s} in %s" % ( tree1.GetName() , len ( missing ) , ', '.join ( missing ) , the_file ) )
                if extra   : logger.warning ( "TTree('%s'): %3d extra   branches: {%s} in %s" % ( tree1.GetName() , len ( extra   ) , ', '.join ( extra   ) , the_file ) )
                
            if ( ( branches1 != leaves1 ) or ( branches2 != leaves2 ) ) and leaves1 != leaves2 :
                missing = sorted ( leaves1 - leaves2 ) 
                extra   = sorted ( leaves2 - leaves1 ) 
                if missing : logger.warning ( "TTree('%s'): %3d missing leaves:   {%s} in %s" % ( tree1.GetName() , len ( missing ) , ', '.join ( missing ) , the_file ) )
                if extra   : logger.warning ( "TTree('%s'): %3d extra   leaves:   {%s} in %s" % ( tree1.GetName() , len ( extra   ) , ', '.join ( extra   ) , the_file ) )
                
    # ===================================================================================================
    ## the specific action for each file 
    def treatFile ( self, the_file ) :
        """ Add the file to TChain
        """
        
        files    = self.files
        has_tree = False ## has at leas opne valid tree ?
        
        for i, cname in enumerate ( self.chain_names ) :

            bad_files = self.bad_files[cname]
            
            ## (1) suppress Warning/Error messages from ROOT
            with rootError() :
            
                ## new temporary chain/tree for this file 
                tree  = ROOT.TChain ( cname )
                tree.Add ( the_file )

                if not tree or 0 == len ( tree.branches () ) :
                    self.bad_files.add ( the_file )
                    if not self.silent : 
                        logger.warning ( "No/empty chain '%s' in file '%s'" % ( cname , the_file ) )
                    bad_files.add ( the_file ) 
                    continue
                
                has_tree = True 
                if self.check and files :
                    last_file = '' 
                    if cname in self.bad_files :
                        bf = self.bad_files[cname] 
                        for f in reversed ( files ) :
                            if f in self.bad_files : continue 
                            
                            
                    chain = ROOT.TChain ( cname )
                    chain.Add ( files [ 0 ] ) 
                    self.check_trees ( tree , chain , the_file )
                    del chain

                del tree
                
        if has_tree : Files.treatFile ( self , the_file )
    
    # ===================================================================================================
    ## printout 
    def __str__(self):
        chains = self.chains 
        result = 'Data(%s): #files=%d #entries={%s}' % (
            ', '.join ( self.chain_names )  ,
            len( self.files ) ,
            ', '.join ( str ( len ( c ) ) for c in chains ) )
        if any ( self.bad_files[c] for c in self.chain_names ) :
            result += ' #bad={%s}' % ( ', '.join ( str ( len ( self.bad_files[c] ) )  for c in self.chain_names ) )
        return result 
    
    # =========================================================================
    ## Print collection of files as table
    #  @code
    #  files = ...
    #  print ( files.table() )    
    #  @endcode
    def table ( self , title = '' , prefix = '' , style = '' ) :
        """ Print collection of files as table
        """
        rows  = [ ( '#' , '#entries' , 'size' , 'name' ) ]
        files = self.files
        nn    = len ( files )
        nfmt  = '%%%dd' % ( math.floor ( math.log10 ( nn ) ) + 1 )

        total_size    = 0
        total_entries = 0
        
        from itertools          import chain 
        for i , f in enumerate ( chain ( files , bad ) , start = 1 ) : 

            row   = [ nfmt % i ]
            fsize = self.get_file_size ( f )
            
            if 0 <= fsize :
                
                ch = ROOT.TChain ( self.chain_name ) 
                ch.Add ( f )
                entries = len ( ch )
                
                row.append ( '%d' % entries )

                vv , unit = fsize_unit ( fsize )
                row.append ( '%3d %s' % ( vv , unit) ) ## file size  

                total_size    += fsize
                total_entries += entries
                
            else :
                
                row.append ( '???' ) ## entries 
                row.append ( '???' ) ## filezies
                
            row .append ( f   ) 
            rows.append ( row )

        ## summary row
        from ostap.logger.colorized import infostr 
        vv , unit  = fsize_unit ( total_size )
        row   = ''                               , \
            infostr ( '%d' % total_entries     ) , \
            infostr ( '%3d %s' % ( vv , unit ) ) , \
            infostr ( self.commonpath )
        
        rows.append ( row )

        title = title if title else "Data(chain='%s')" % self.chain_name 
        import ostap.logger.table as T
        return T.table ( rows , title = title , prefix = prefix , alignment = 'rrrw' , style = style ) 

    # =========================================================================
    ## check operations
    def check_ops ( self , other ):
        """ Check operations
        """
        return isinstance ( other , Data ) and self.chain_names == other.chain_names

    ## Get RDatataFrame for the given chain
    def get_frame ( self , index ) :
        """ Get RDatataFrame for the given chain
        """
        from ostap.frames.frames import DataFrame, frame_progress
        ch = self.get_chain ( index  ) 
        fr = DataFrame ( ch )
        if not self.filent: pb = frame_progress ( fr , len ( ch ) ) 
        return fr 
    
    # =========================================================================
    ## get DataFrame for the chain
    #  @see ROOT::RDataFrame
    @property
    def frames ( self ) :
        """'frames': get all DataFrame for the chain """
        result = []
        for index in range ( self.nchains ) : 
            result.append ( self.get_frame ( index ) )
        return tuple ( result ) 

    # =========================================================================
    ## get DataFrame for the chain
    #  @see ROOT::RDataFrame
    @property
    def frame ( self ) :
        """'frame': get the DataFrame for the first chain (same as `frame`)"""
        return self.get_frame ( 0 ) ;
    
    # =========================================================================
    ## get DataFrame for the chain
    #  @see ROOT::RDataFrame
    @property
    def frame1 ( self ) :
        """'frame1': get DataFrame for the first chain (same ad `frame`)"""
        return self.get_frame ( 0 ) ;
    
    # =========================================================================
    ## get DataFrame for the chain
    #  @see ROOT::RDataFrame
    @property
    def frame2 ( self ) :
        """'frame2': get DataFrame for the second chain"""
        return self.get_frame ( 1 ) ;

# =============================================================================
if '__main__' == __name__ :

    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
