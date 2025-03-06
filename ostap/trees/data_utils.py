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
    'Data'        , ## collect ROOT files and create     TChain
    'Data2'       , ## collect ROOT files and create two TChain objects 
    )
# =============================================================================
from   ostap.core.core     import rootError, rootWarning
from   ostap.io.files      import Files, fsize_unit 
from   ostap.io.root_file  import RootFiles
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
                  chains               , 
                  files        = []    ,
                  description  = ''    , 
                  maxfiles     = -1    ,
                  check        = True  , 
                  silent       = False ,
                  sorted       = True  , 
                  parallel     = False ,
                  transform    = None  ) : 

        ## we will need Ostap machinery for trees&chains here
        import ostap.trees.trees
        
        if   isinstance ( chains , ROOT.TTree   ) : chains = chains.path , 
        elif isinstance ( chains , string_types ) : chains = chains      ,

        assert chains and all ( isinsatnce ( c , chain_types ) for c in chains ) , \
            "Invald type of `chains` %s" % str ( chains )

        self.__check = True if check else False
        
        ## chain names
        self.__chain_names = tuple ( ( c.path if isinstance ( c , ROOT.TTree ) else c ) for c in chains ) 

        ## update description
        if not description :
            description = "Chains:{%s}" %  ( ', '.join ( self.__chain_names ) ) 
            
        ## initialize the  base class 
        RootFiles.__init__( self                       ,
                            files        = files       ,
                            descxription = description ,
                            maxfiles     = maxfiles    ,
                            silent       = silent      ,
                            sorted       = sorted      ,
                            parallel     = parallel    )
        
    # =========================================================================    
    @property
    def validate ( self ) :
        """'check' : make check for `TTree`/`TChain` structures
        """
        return self.__check
    
    @property
    def chain_names ( self ) :
        """'chain_name' : the name of TTree/TChain object
        """
        return self.__chain_names

    # ========================================================================
    ## get the chain by index  
    def get_chain ( self , index ) :
        """ Get the chain by index
        """
        if not 0 <= index < len ( self.__chain_names ) :
            raise IndexError ( "Index out the range!" )
        ch    = ROOT.TChain ( self.__chain_names [ index ] )
        files = self.files 
        for f in files : ch.Add ( f ) 
        return ch
    
    @property
    def chains ( self ) :
        """'chains' : (re)built and return all `TChain` objects"""
        result = []
        for i , _ in enunerate ( self.chain_names ) :
            result.append ( self.get_chain ( i ) ) 
        return tuple ( results )
        
    @property
    def chain ( self ) :
        """'chain' : (re)built and return the fist `TChain` object , same as `chain1`"""
        return self.get_chain ( 0 )
    
    @property
    def chain1 ( self ) :
        """'chain1' : (re)built and return the fist `TChain` object , same as `chain`"""
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
                if missing : logger.warning ( "TTree('%s'): %3d missing branches: %s in %s" %  ( tree1.GetName() , len ( missing ) , missing , the_file ) )
                if extra   : logger.warning ( "TTree('%s'): %3d extra   branches: %s in %s" %  ( tree1.GetName() , len ( extra   ) , extra   , the_file ) )
                
            if ( ( branches1 != leaves1 ) or ( branches2 != leaves2 ) ) and leaves1 != leaves2 :
                missing = list ( sorted ( leaves1 - leaves2 ) )
                extra   = list ( sorted ( leaves2 - leaves1 ) ) 
                if missing : logger.warning ( "TTree('%s'): %3d missing leaves:   %s in %s" %  ( tree1.GetName() , len ( missing ) , missing , the_file ) )
                if extra   : logger.warning ( "TTree('%s'): %3d extra   leaves:   %s in %s" %  ( tree1.GetName() , len ( extra   ) , extra   , the_file ) )
                
    # ===================================================================================================
    ## the specific action for each file 
    def treatFile ( self, the_file ) :
        """ Add the file to TChain
        """
        
        files = self.files 
        for i, cname in enumerate ( self.chain_names ) :
            
            ## (1) suppress Warning/Error messages from ROOT
            with rootError() :
            
                ## new temporary chain/tree for this file 
                tree  = ROOT.TChain ( cname )
                tree.Add ( the_file )

                if not tree or 0 == len ( tree.branches () ) :
                    self.bad_files.add ( the_file )
                    if not self.silent : 
                        logger.warning ( "No/empty chain '%s' in file '%s'" % ( cname , the_file ) )
                    continue
            
                if self.check and files : 
                    
                    chain = ROOT.TChain ( cname )
                    chain.Add ( files [ 0 ] ) 
                    self.check_trees ( tree , chain , the_file )
                    del chain
                            
                Files.treatFile ( self , the_file )
                
                del tree
    
    # ===================================================================================================
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
        bad   = sorted ( self.bad_files )
        nn    = max ( len ( files ) , len ( bad ) ) 
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
        return isinstance ( other , Data ) and self.chain_name == other.chain_name

    ## Get RDatataFrame for the given chain
    def get_frame ( self , index ) :\
        """ Get RDatataFrame for the given chain
        """
        from ostap.frames.frames import DataFrame, frame_progress
        ch = self.get_chain ( index  ) 
        fr = DataFrame ( ch )
        if not self.filent:
            pb = frame_progress ( fr , len ( ch ) ) 
        return fr 

    # =========================================================================
    ## get DataFrame for the chain
    #  @see ROOT::RDataFrame
    @property
    def frames ( self ) :
        """'frames': get all DataFrame for the chain """
        result = []
        for index, _ inenumerate ( self.chain_names ) : 
            result.append ( self.get_frame ( index ) )
        return tuple ( result ) 

    # =========================================================================
    ## get DataFrame for the chain
    #  @see ROOT::RDataFrame
    @property
    def frame ( self ) :
        """'frame': get the DataFrame for the chain"""
        return self.get_frame ( 0 ) ;
    
    # =========================================================================
    ## get DataFrame for the chain
    #  @see ROOT::RDataFrame
    @property
    def frame1 ( self ) :
        """'frame1': get DataFrame for the first chain"""
        return self.get_frame ( 1 ) ;
    
    # =========================================================================
    ## get DataFrame for the chain
    #  @see ROOT::RDataFrame
    @property
    def frame2 ( self ) :
        """'frame2': get DataFrame for the second chain"""
        return self.get_frame ( 2 ) ;


# =============================================================================
## @class Data2
#  Simple utility to access two chains in the set of ROOT-files
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @author Alexander BARANOV a.baranov@cern.ch
#  @date   2014-06-08  
class Data2(Data):
    """ Simple utility to access to certain chain in the set of ROOT-files    
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
                  parallel    = False ,
                  transform   = None  ) : 

        Data.__init__ ( self ,
                        chains      = ( chain1 , chain2 ) ,
                        files       = files               ,
                        description = description         ,
                        maxfiles    = maxfiles            ,
                        check       = check               ,
                        silent      = silent              ,
                        sorted      = sorted              ,
                        parallel    = parallel            ,
                        transform   = transform           )
        
    # =========================================================================
    ## check operations
    def check_ops ( self , other ):
        """ Check operations
        """
        return isinstance ( other , Data2 )          and \
               self.chain1_name == other.chain1_name and \
               self.chain2_name == other.chain2_name

# =============================================================================
if '__main__' == __name__ :

    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
