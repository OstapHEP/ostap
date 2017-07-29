#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file data.py 
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

>>> data  = DataAndLumi('Bc/MyTree', '*.root' )
>>> chain = data.chain
>>> flist = data.files 
>>> lumi  = data.lumi
>>> print data.getLumi() 

>>> data  = DataAndLumi('Bc/MyTree', [ 'A/*.root' , 'B/B/D/*.root' ] )
>>> chain = data.chain
>>> flist = data.files 
>>> lumi  = data.lumi
>>> print data.getLumi() 

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
import ROOT, glob 
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
protocols = (
    'root:/' ,
    'http:/' ,
    'https:/',    
    )
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
    def __init__( self                  ,
                  files                 ,
                  description = ""      ,
                  maxfiles    = 1000000 ,
                  silent      = False   ) :  
        #
        from copy import deepcopy
        #
        self.files        = []
        self.patterns     = deepcopy ( files ) 
        self.description  = description
        self.maxfiles     = maxfiles
        self.silent       = silent 
        #
        self.add_files ( deepcopy ( files ) ) 

    def __getstate__ ( self ) :

        return {
            'files'       : self.files       ,
            'patterns'    : self.patterns    ,
            'description' : self.description ,
            'maxfiles'    : self.maxfiles    ,
            'silent'      : self.silent      ,            
            }
    
    def __setstate__ ( self , state ) :
        self.files       = state.get('files'      , []    ) 
        self.patterns    = state.get('patterns'   , []    )
        self.description = state.get('description', ''    )
        self.maxfiles    = state.get('maxfiles'   , 1e+6  )
        self.silent      = state.get('silent'     , False )
        
    ## add files 
    def add_files ( self , patterns ) :
        """ Add files/patterns to data collector
        """
        
        if isinstance ( patterns , str ) : patterns = [ patterns ]
        
        _files  = set ()
        for pattern in patterns :

            _added = False  
            for p in protocols :
                if p in pattern : 
                    if not pattern in self.files :
                        _files.add  ( pattern )
                    _added = True
                    break
                
            if not _added : 
                for f in glob.iglob ( pattern ) :
                    if not f in self.files :
                        _files.add ( f )
                        
        if not self.silent :
            logger.info ('Loading: %s  #patterns/files: %s/%d' % ( self.description ,
                                                                   len(patterns)    , 
                                                                   len( _files )    ) )
        ## update list of patterns  
        self.patterns += patterns
        
        from ostap.utils.progress_bar import ProgressBar 
        with ProgressBar ( max_value = len(_files) , silent = self.silent ) as bar :
            self.progress = bar 
            for f in _files : 
                if len ( self.files ) < self.maxfiles : self.treatFile  ( f )
                else :
                    logger.warning ('Maxfiles limit is reached %s ' % self.maxfiles )
                    break

        if not self.silent :
            logger.info ('Loaded: %s' % self )
            
    ## the specific action for each file 
    def treatFile ( self, the_file ) :
        self.files.append ( the_file )
        self.progress += 1
        
    ## merge two sets togather
    def add_data ( self , another ) :
        """ Merge two datasets
        """
        return self.add_files ( another.patterns )    

    ## clone it! 
    def  clone ( self ) :
        """ Clone the object
        """
        from copy import deepcopy
        return deepcopy ( self )

    ## merge two sets togather
    def __iadd__ ( self , another ) :
        """ Merge two datasets
        """
        self.add_data ( another )    
        return self

    
    ## merge two sets togather
    def __add__ ( self , another ) :
        """ Merge two sets together
        """
        result  = self.clone() 
        result += another 
        return result
        
    ## clone !
    def clone  ( self ) :
        """ Clone the  object
        """
        result = Files( files        = []               ,
                        description  = self.description ,
                        maxfiles     = self.maxfiles    ,
                        silent       = True             )
        
        from copy import deepcopy
        result.files    = deepcopy ( self.files    )
        result.patterns = deepcopy ( self.patterns )
        result.silent   = self.silent 
                                     
        return result 


    ## printout 
    def __str__(self):
        """The specific printout
        """
        return "<#files: %4d>" % len ( self.files ) 
    
    def __repr__    ( self ) : return self.__str__()
    def __nonzero__ ( self ) : return bool(self.files)
    def __len__     ( self ) : return len (self.files)
    
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
    
    def __init__( self                  ,
                  chain                 ,
                  files       = []      ,
                  description = ''      , 
                  maxfiles    = 1000000 ,
                  silent      = False   ) :  
        
        self.e_list1    = set()  
        self.chain      = ROOT.TChain ( chain ) if isinstance ( chain , str ) else chain 
        if not description :
            description = self.chain.GetName()
        Files.__init__( self , files , description  , maxfiles , silent )


    def __getstate__ ( self ) :

        state = Files.__getstate__( self )
        state [ 'e_list1' ] = self.e_list1
        state [ 'chain'   ] = self.chain.GetName()
        
        return state

    def __setstate__ ( self , state ) :

        Files.__setstate__ ( self , state )
        
        self.e_list1     = state.get('e_list1'    , set() )
        self.chain       = ROOT.TChain( state['chain'] )
        for f in  self.files :
            self.chain.Add ( f ) 
        
    ## the specific action for each file 
    def treatFile ( self, the_file ) :
        """Add the file to TChain
        """
        ## suppress Warning/Error messages from ROOT 
        from ostap.logger.utils import rootError
        with rootError() :
            
            tmp1 = ROOT.TChain ( self.chain .GetName() )
            tmp1.Add ( the_file )

            if  tmp1.ok() :
                Files.treatFile ( self     , the_file ) 
                self.chain .Add ( the_file )
            else : 
                self.e_list1.add ( the_file )
                logger.warning ( "No/empty chain  '%s'      in file '%s'" % ( self.chain .GetName() ,
                                                                              the_file )            ) 
                
    ## printout 
    def __str__(self):
        return  "<#files: {}; Entries: {}>".format ( len ( self.files ) ,
                                                     len ( self.chain ) )
    
    def __nonzero__ ( self )  :
        return bool(self.files) and bool(len(self.chain))

    ## merge two sets togather
    def add_data ( self , another ) :
        """ Merge two datasets
        """
        self.patterns += another.patterns
        for f in another.e_list1 :self.e_list1.add ( f ) 
        for f in another.files :
            if not f in self.files : 
                self.chain. Add ( f )
                
        return self 

    ## clone !
    def clone  ( self ) :
        """ Clone the  object
        """
        result         = Data ( chain       = self.chain.GetName() ,
                                files       = []                   ,
                                description = self.description     , 
                                maxfiles    = self.maxfiles        ,
                                silent      = True                 )
        
        result.chain   = ROOT.TChain ( self.chain.GetName() )
        for f in self.files : result.chain.Add ( f ) 
        
        from copy import deepcopy
        result.files    = deepcopy ( self.files    ) 
        result.patterns = deepcopy ( self.patterns ) 
        result.e_list1  = deepcopy ( self.e_list1  ) 
        
        return result 
    
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
    
    def __init__( self                  ,
                  chain1                ,
                  chain2                , 
                  files       = []      ,
                  description = ''      ,
                  maxfiles    = 1000000 ,
                  silent      = False   ) :  
        
        self.e_list2 = set()
        self.chain2  = ROOT.TChain ( chain2 ) if isinstance ( chain2 , str ) else chain2 
        if not description :
            description = chain1.GetName() if hasattr ( chain1 , 'GetName' ) else str(chain1)
            description = "%s&%s" % ( description , self.chain2.GetName() )
        Data.__init__( self , chain1 , files , description , maxfiles , silent )
        self.chain1  = self.chain 

    def __getstate__ ( self ) :

        state = Data.__getstate__( self )
        state [ 'e_list2' ] = self.e_list2
        state [ 'chain2'  ] = self.chain2.GetName()
        
        return state

    def __setstate__ ( self , state ) :

        Files.__setstate__ ( self , state )
        
        self.e_list1     = state.get   ('e_list1'    , set() )
        self.e_list2     = state.get   ('e_list2'    , set() )
        self.chain       = ROOT.TChain ( state['chain' ] )
        self.chain2      = ROOT.TChain ( state['chain2'] )
        for f in  self.files :
            self.chain .Add ( f ) 
            self.chain2.Add ( f )
            
        self.chain1  = self.chain 


    ## the specific action for each file 
    def treatFile ( self, the_file ) :
        """Add the file to TChain
        """
        ## suppress Warning/Error messages from ROOT 
        from ostap.logger.utils import rootError
        with rootError() :
            tmp1 = ROOT.TChain ( self.chain .GetName() )
            tmp1.Add ( the_file )            
            tmp2 = ROOT.TChain ( self.chain2.GetName() )
            tmp2.Add ( the_file )

            if  tmp1.ok() and tmp2.ok() :
                Files.treatFile ( self     , the_file ) 
                self.chain .Add ( the_file )
                self.chain2.Add ( the_file )
            elif tmp2.ok() and not tmp1.ok() : 
                self.e_list1.add ( the_file  )
                logger.warning ( "No/empty chain1 '%s'      in file '%s'" % ( self.chain .GetName() ,
                                                                              the_file )            ) 
            elif tmp1.ok() and not tmp2.ok() : 
                self.e_list2.add ( the_file )
                logger.warning ( "No/empty chain2 '%s'      in file '%s'" % ( self.chain2.GetName() ,
                                                                          the_file )            ) 
            else :
                self.e_list1.add ( the_file )
                self.e_list2.add ( the_file )
                logger.warning ( "No/empty chains '%s'/'%s' in file '%s'" % ( self.chain .GetName() ,
                                                                          self.chain2.GetName() ,
                                                                          the_file )            )
    ## printout 
    def __str__(self):    
        return  "<#files: {}; Entries1: {}; Entries2: {}>".format ( len ( self.files  ) ,
                                                                    len ( self.chain  ) ,
                                                                    len ( self.chain2 ) )
    
    def __nonzero__ ( self )  :
        return bool(self.files) and bool(len(self.chain)) and bool(len(self.chain2))

    ## merge two sets togather
    def add_data ( self , another ) :
        """ Merge two datasets
        """
        self.patterns += another.patterns
        for f in another.e_list1 :self.e_list1.add ( f ) 
        for f in another.e_list2 :self.e_list2.add ( f ) 
        for f in another.files :
            if not f in self.files : 
                self.chain1. Add  ( f )
                self.chain2. Add  ( f )
                self.files.append ( f ) 
        return self 
    
    ## clone !
    def clone  ( self ) :
        """ Clone the  object
        """
        result         = Data2 ( chain1      = self.chain .GetName() ,
                                 chain2      = self.chain2.GetName() ,
                                 files       = []                    ,
                                 description = self.description      , 
                                 maxfiles    = self.maxfiles         ,
                                 silent      = True                  )
        
        result.chain   = ROOT.TChain ( self.chain .GetName() )
        result.chain2  = ROOT.TChain ( self.chain2.GetName() )
        for f in self.files :
            result.chain .Add ( f ) 
            result.chain2.Add ( f )
            
        result.chain1  = result.chain

        from copy import deepcopy
        result.files    = deepcopy ( self.files    ) 
        result.patterns = deepcopy ( self.patterns ) 
        result.e_list1  = deepcopy ( self.e_list1  ) 
        result.e_list2  = deepcopy ( self.e_list2  ) 
        
        return result 
    
# =============================================================================
## @class DataAndLumi
#  Simple utility to access to certain chain in the set of ROOT-files
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @author Alexander BARANOV a.baranov@cern.ch
#  @date   2014-06-08  
class DataAndLumi(Data2):
    """Simple utility to access to certain chain in the set of ROOT-files    
    >>> data  = DataAndLumi('Bc/MyTree', '*.root' )
    >>> chain = data.chain
    >>> flist = data.files 
    >>> lumi  = data.lumi
    >>> print data.getLumi() 
    """

    def __init__( self               ,
                  chain              ,
                  files       = []   ,
                  description = ''   , 
                  lumi_chain  = 'GetIntegratedLuminosity/LumiTuple' , 
                  maxfiles    = 1000000                             ,
                  silent      = False ) :  

        if not description :
            description = chain.GetName() if hasattr ( chain , 'GetName' ) else str(chain)
        Data2.__init__ ( self , chain , lumi_chain , files , description , maxfiles  , silent ) 
        self.lumi = self.chain2 

    def __getstate__ ( self ) :
        return Data2.__getstate__( self )

    def __setstate__ ( self , state ) :
        Data2.__setstate__ ( self , state )                   
        self.lumi  = self.chain2 
        
    ## get the luminosity 
    def getLumi ( self ):
        """Get the luminosity
        """
        ## suppress Warning/Error messages from ROOT 
        from ostap.logger.utils import rootError
        with rootError() :
            from   ostap.contrib.lhcb.lumi import getLumi
            return getLumi ( self.chain2  )
        
    ## printout 
    def __str__(self):        
        return "<Luminosity: {}pb-1; #files: {}; Entries: {}>".format(
            self.getLumi()     ,
            len ( self.files ) ,
            len ( self.chain ) )
    
    def __nonzero__ ( self )  :
        return bool(self.files) and bool(len(self.chain)) and bool(len(self.chain2))

    ## clone !
    def clone  ( self ) :
        """ Clone the  object
        """
        result          = DataAndLumi (
            chain       = self.chain .GetName() ,
            lumi_chain  = self.chain2.GetName() ,
            files       = []                    ,
            description = self.description      , 
            maxfiles    = self.maxfiles         ,
            silent      = True                  )
        
        result.chain   = ROOT.TChain ( self.chain .GetName() )
        result.chain2  = ROOT.TChain ( self.chain2.GetName() )
        for f in self.files :
            result.chain .Add ( f ) 
            result.chain2.Add ( f )
            
        result.chain1  = result.chain
        result.lumi    = result.chain2

        from copy import deepcopy
        result.files    = deepcopy ( self.files    ) 
        result.patterns = deepcopy ( self.patterns ) 
        result.e_list1  = deepcopy ( self.e_list1  ) 
        result.e_list2  = deepcopy ( self.e_list2  ) 
        
        return result 
    
# =============================================================================

# =============================================================================
if '__main__' == __name__ :

    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
