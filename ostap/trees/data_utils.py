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
                  chain                ,
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

        ## decorate files 
        if isinstance ( files , str ) : files = [ files ]

        self.__check = True if check else False
        
        ## chain name 
        self.__chain_name = chain if isinstance ( chain , str ) else chain.name 

        ## update description
        if not description : description = "ROOT.TChain(%s)" % self.chain_name 
        
        ## initialize the  base class 
        RootFiles.__init__( self , files , description  , maxfiles , silent = silent , sorted = sorted , parallel = parallel )

    # =========================================================================    
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
        ## (1) suppress Warning/Error messages from ROOT 
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
    def table ( self , title = '' , prefix = '' ) :
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
        return T.table ( rows , title = title , prefix = prefix , alignment = 'rrrw' ) 

    # =========================================================================
    ## check operations
    def check_ops ( self , other ):
        """ Check operations
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

        ## decorate files 
        if isinstance ( files , str ) : files = [ files ]

        ## chain name 
        self.__chain2_name = chain2 if isinstance ( chain2 , str ) else chain2.name 

        if not description :
            description = chain1.GetName() if hasattr ( chain1 , 'GetName' ) else str ( chain1 )
            description = "%s&%s" % ( description , self.__chain2_name  )

        Data.__init__( self                      ,
                       chain       = chain1      ,
                       files       = files       ,
                       description = description ,
                       maxfiles    = maxfiles    ,
                       check       = check       ,
                       silent      = silent      ,
                       sorted      = sorted      , 
                       parallel    = parallel    ,
                       transform   = transform   ) 
        
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
        """ Add the file to TChain
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
                
                RootFiles.treatFile ( self , the_file ) 
                
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
                nc2    = len ( chain2 )
                del  chain2
                
            ne  = len ( self.bad_files )
            
        sf  =  set ( self.files ) == set ( self.files2 )
        
        if not self.bad_files :
            return "<#files: {}; Entries: {}/{}>"   .format ( nf ,      nc  , nc2 ) if sf else \
                   "<#files: {}/{}; Entries: {}/{}>".format ( nf , nf2 , nc , nc2 )
        else :
            return "<#files: {}; Entries: {}/{}; No/empty :{}>"   .format ( nf ,       nc , nc2 , ne ) if sf else \
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

    # =========================================================================
    ## Print collection of files as table
    #  @code
    #  files = ...
    #  print ( files.table() )    
    #  @endcode
    def table ( self , title = '' , prefix = '' ) :
        """ Print collection of files as table
        """
        
        rows  = [ ( '#' , '#entries1' , '#entries2' , 'size' , 'name' ) ]
        
        files = self.files
        bad   = sorted ( self.bad_files )
        nn    = max ( len ( files ) , len ( bad ) ) 
        nfmt  = '%%%dd' % ( math.floor ( math.log10 ( nn ) ) + 1 )

        total_entries1 = 0
        total_entries2 = 0
        total_size     = 0
        
        from itertools          import chain 
        for i , f in enumerate ( chain ( files , bad ) , start = 1 ) : 

            row   = [ nfmt % i ]
            fsize = self.get_file_size ( f )
            
            if 0 <= fsize :
                
                ch1 = ROOT.TChain ( self.chain1_name ) 
                ch1.Add ( f )
                entries1 = len ( ch1 )                
                row.append ( '%d' % entries1 )

                ch2 = ROOT.TChain ( self.chain2_name ) 
                ch2.Add ( f )
                entries2 = len ( ch2 )                
                row.append ( '%d' % entries2 )
                
                vv , unit = fsize_unit ( fsize ) 
                row.append ( '%3d %s' % ( vv , unit ) ) ## value 
                
                total_size     += fsize
                total_entries1 += entries1
                total_entries2 += entries2

            else :
                
                row.append ( '???' ) ## entries1
                row.append ( '???' ) ## entries2 
                row.append ( '???' ) ## file size
                
            row .append ( f   ) 
            rows.append ( row )

        ## summary row
        from ostap.logger.colorized import infostr 
        vv , unit  = fsize_unit ( total_size )
        row   = ''                               , \
            infostr ( '%d' % total_entries1  )   , \
            infostr ( '%d' % total_entries2  )   , \
            infostr ( '%3d %s' % ( vv , unit ) ) , \
            infostr ( self.commonpath ) 
        rows.append ( row )

        title = title if title else "Data2(chani1='%s',chani2='%s')" % ( self.chain1_name , self.chain2_name ) 
        import ostap.logger.table as T
        return T.table ( rows , title = title , prefix = prefix , alignment = 'rrrrw' ) 


# =============================================================================
if '__main__' == __name__ :

    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
