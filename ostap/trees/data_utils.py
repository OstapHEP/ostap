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
protocols = (
    'root:/' ,
    'http:/' ,
    'https:/',    
    )
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
                  maxfiles    = -1      ,
                  silent      = False   ) :
        #
        if isinstance ( files , str ) :  files =  [ files ]
        #
        #
        
        self.description  = description
        self.maxfiles     = maxfiles
        self.__silent     = silent 

        from copy import deepcopy
        self.__patterns   = list ( set ( deepcopy ( files ) ) )
        self.__patterns.sort()
        
        self.__files      = []        

        # =====================================================================
        # convert list of patterns into the list of files 
        # =====================================================================
        
        _files   = self.the_files ()

        if not self.silent :
            logger.info ('Loading: %s  #patterns/files: %s/%d' % ( self.description   ,
                                                                   len(self.patterns) , 
                                                                   len( _files )    ) )        
        self.add_files ( _files )

        if not self.silent :
            logger.info ('Loaded: %s' % self )

    ## append the file to the list of the files 
    def _append_to_files (  self , f ) :
        """Append the file to the list of the files"""
        self.__files.append ( f )
        
    @property 
    def files     ( self ) :
        """``files'' : the list of files"""
        return tuple ( self.__files    )

    @property
    def patterns  ( self ) :
        """``patterns'' : the list of patterns"""
        return tuple ( self.__patterns )
    
    @property
    def silent    ( self ) :
        """``silent'' : silent processing?"""
        return self.__silent
    @silent.setter
    def silent    ( self , value ) :
        self.__silent = True if value else  False
        
    # =========================================================================
    ## Get the list of files form the patterns 
    def the_files ( self ) :
        """Get the list of files form the patterns"""
        
        _files = set ()
        
        _prots = 0 
        _norms = 0 
        for pattern in  self.patterns :
            
            _added = False  
            for p in protocols :
                if pattern.startswith ( p ) or p in pattern : 
                    if not pattern in _files :
                        _files.add  ( pattern )
                        _prots +=1 
                    _added = True
                    break
                
            if not _added : 
                for f in glob.iglob ( pattern ) :
                    f = os.path.abspath  ( f )
                    f = os.path.normpath ( f )                    
                    if not f in _files :
                        _norms += 1 
                        _files.add ( f )

        ## final sort 
        _files = list ( _files )
        _files.sort   ()
        
        return tuple ( _files )
        
    # =========================================================================
    ## set files 
    def set_files ( self , files ) :
        """Set files"""
        assert isinstance ( files , ( list , tuple , set ) ),\
               'Invalid list of files %s' % files
        self.__files = list ( set ( files ) ) 
        self.__files.sort()

    # =========================================================================
    ## set patterns
    def set_patterns ( self , patterns ) :
        """Set files"""
        assert isinstance ( patterns , ( list , tuple , set ) ),\
               'Invalid list of patterns %s' % patterns
        self.__patterns = list ( set ( patterns ) ) 
        self.__patterns.sort()

    def __getstate__ ( self ) :

        return {
            'files'       : self.files       ,
            'patterns'    : self.patterns    ,
            'description' : self.description ,
            'maxfiles'    : self.maxfiles    ,
            'silent'      : self.silent      ,            
            }
    
    def __setstate__ ( self , state ) :
        self.__files     = state.get ( 'files'      , []    ) 
        self.__patterns  = state.get ( 'patterns'   , []    )
        self.description = state.get ( 'description', ''    )
        self.maxfiles    = state.get ( 'maxfiles'   , -1    )
        self.silent      = state.get ( 'silent'     , False )

    ## add files 
    def add_files ( self , files ) :
        """ Add files/patterns to data collector
        """
        
        if isinstance ( files  , str ) : files  = [ files  ]

        ## eliminate duplicates and sort 
        files = list  ( set ( files ) )
        files.sort ()
        nf    = len ( files )
        max_value = nf if 0 >= self.maxfiles else min ( nf , self.maxfiles )
        from ostap.utils.progress_bar import ProgressBar 
        with ProgressBar ( max_value = max_value , silent = self.silent ) as bar :
            self.progress = bar 
            for f in files :
                if   0 >= self.maxfiles                 : self.treatFile ( f ) 
                elif len ( self.files ) < self.maxfiles : self.treatFile ( f )
                else :
                    logger.debug ('Maxfiles limit is reached %s ' % self.maxfiles )
                    break
                
    ## the specific action for each file 
    def treatFile ( self, the_file ) :
        self.__files.append ( the_file )
        self.progress += 1
 
     ## merge two sets together
    def add_data ( self , other ) :
        """ Merge two datasets
        """
        #
        f1 = set ( self .files    )
        f2 = set ( other.files    )
        p1 = set ( self .patterns )
        p2 = set ( other.patterns )
        
        self.set_files    ( f1 | f2 )
        self.set_patterns ( p1 | p2 )
        
        return self
    
    ## clone it! 
    def  clone ( self ) :
        """ Clone the object
        """
        from copy import deepcopy
        return deepcopy ( self )
    
    ## merge two sets together
    def __iadd__ ( self , other ) :
        """ Merge two datasets
        >>> ds   = ...
        >>> ds1  = ...
        >>> ds  += ds1
        >>> ds  |= ds1 ## ditto 
        """
        if not isinstance ( other , Files ) :
            return NotImplemented
        return self.add_data ( other )    
    
    ## merge two sets together
    def __add__ ( self , other ) :
        """ Merge two sets together
        >>> ds1 = ...
        >>> ds2 = ...
        >>> ds  = ds1 + ds2 
        >>> ds  = ds1 | ds2 ## ditto
        """
        if not isinstance ( other , Files ) :
            return NotImplemented
        result  = self.clone() 
        result += other 
        return result
    
    __or__  = __add__
    __ror__ = __add__
    __ior__ = __iadd__

    ## get an intersection of two datasets 
    def __and__ (  self , other ) :
        """ get intersection of two sets
        >>> ds1 = ...
        >>> ds2 = ...
        >>> ds  = ds1 & ds2 ## get intersection 
        >>> ds  = ds1 * ds2 ## ditto
        """
        if  not isinstance  ( other , Files ) : return NotImplemented

        result = Files( files        = []               ,
                        description  = self.description ,
                        maxfiles     = self.maxfiles    ,
                        silent       = True             )
        
        f1 = set ( self .files    )
        f2 = set ( other.files    )
        p1 = set ( self .patterns )
        p2 = set ( other.patterns )
        
        result.files    = list ( f1 & f2 )
        result.patterns = list ( p1 | p2 )
        result.silent   = self.silent
        
        return result
    
    __rand__ = __and__
    __mul__  = __and__
    __rmul__ = __and__
    
    ## clone !
    def clone  ( self ) :
        """ Clone the  object
        """
        result = Files( files        = []               ,
                        description  = self.description ,
                        maxfiles     = self.maxfiles    ,
                        silent       = True             )
        
        from copy import deepcopy
        result.set_files    ( deepcopy ( self.files    ) ) 
        result.set_patterns ( deepcopy ( self.patterns ) )
        result.silent     = self.silent 
                                     
        return result 

    ## get a sample of at most  n-elements (if n is integer and >=1 )  or n-fraction 
    def sample_files ( self ,  n , sort ) :
        
        if   isinstance ( n , int   ) and 1 <= n <= len  ( self.files ) :
            files = random.sample ( self.files , n )
            if sort : files.sort()
            return files 
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
        files = self.sample_files ( n , sort )
        return Files ( files       = files ,
                       description = kwargs.get( 'description' , self.description ) ,
                       maxfiles    = kwargs.get( 'maxfiles'    , self.maxfiles    ) ,
                       silent      = kwargs.get( 'silent'      , self.silent      ) )
    
    ##  reload!
    def reload ( self , silent = True ) :
        
        prev_silent   = self.silent
        self.__silent = silent
        
        self.__files  = []
        self.add_files ( self.the_files ()  )
        
        self.__silent = prev_silent
        
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
        files = self.files[ item ]
        return Files ( files       = files            ,
                       description = self.description ,
                       maxfiles    = self.maxfiles    ,
                       silent      = self.silent      )
    ## printout 
    def __str__(self):
        """The specific printout
        """
        return "<#files: %4d>" % len ( self.files ) 
    
    def __repr__    ( self ) : return self.__str__()
    def __nonzero__ ( self ) : return bool ( self.files )
    def __len__     ( self ) : return len  ( self.files )


    ## merge all files using <code>hadd</code> script from ROOT
    #  @param output  name of the output merged file, if None,
    #                 the temporary name will be generated,
    #                 that will be deleted at the end of the session
    #  @param opts   options for command <code>hadd</code>
    #  @return the name of the merged file 
    def hadd ( self , output = None , opts = "-ff" ) :
        """<erge all files using <code>hadd</code> script from ROOT
        - `output`  name of the output merged file
        - `opts`   options for command <code>hadd</code>
        It returns the name of the merged file
        
        If no output file name is specified, the temporary name
        will be generate and the temporary file will be deleted
        at the end of the session
        """
        if not output :
            import ostap.utils.cleanup as CU
            output = CU.tempfile ( suffix = '.root' )
            

        import subprocess
        
        args    = [ 'hadd' ] + opts.split() + [ output ] + [ f for f in self.files ]
        subprocess.check_call ( args )
        
        if os.path.exists ( output ) and os.path.isfile ( output ) :
            return output 
        
        raise IOError ( "The output file %s does not exist!" % output )


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
                  maxfiles    = -1      ,
                  silent      = False   ,
                  quick       = False   ) :  

        ## we will need Ostap machinery for trees&chains here
        import ostap.trees.trees 

        self.__quick = True if quick else False
        
        ## decorate files 
        if isinstance ( files , str ) : files = [ files ]

        if quick and 0 < maxfiles :
            logger.warning("``Quick'' and ``maxfiles>0'' are not compatible, switch to ``quick=False''")
            quick = False 
            
        if quick :
            quick = self._use_quick_(  files )
            if not quick :
                logger.warning("Patterns are not compatible with quick search, switch to ``quick=False''")
        
        self.e_list1    = set()        
        self.chain      = ROOT.TChain ( chain ) if isinstance ( chain , str ) else chain 
        if not description : description = self.chain.GetName()
        ##
        if not quick : 
            Files.__init__( self , files , description  , maxfiles , silent = silent )
        else :
            Files.__init__( self , []    , description  , maxfiles , silent = True   )
            self.silent = silent
            self.files  = self._quick_add_ ( self.chain ,  files )
            if not self.silent : logger.info ('Loaded: %s' % self )

    @property
    def  quick ( self ) :
        """``quick'' :  quick processing?"""
        return self.__quick
    
    ## can ``quick-add'' be applied for certain patterns?
    def _use_quick_ ( self , files ) :
        """Can ``quick-add'' be applied for certain patterns?
        """
        for f in files :
            _d,_f = os.path.split ( f )
            if _d :
                if 0 <= _d.find ( '*' )                         : return False
                if 0 <= _d.find ( '[' ) or 0 <= _d.find ( ']' ) : return False
        return True


    ## quick add  files to the chain 
    def _quick_add_ ( self , chain , files ) :
        """Quick add  files to the chain 
        """
        if isinstance ( files , str ) :  files = [  files ]
        _files = set ( files )
        _used  = set ( chain.files () ) 
        for f in _files :
            if f in _used : continue 
            n = chain.Add ( f ) 
            if n <= 0 : logger.info ( 'Quick-add: no files for ``%s'' pattern are found!' )
        return chain.files()[:] 
            
    ## picking 
    def __getstate__ ( self ) :

        state = Files.__getstate__( self )
        state [ 'e_list1' ] = self.e_list1
        state [ 'chain'   ] = self.chain.GetName()
        stat  [ 'quick'   ] = self.quick 
        return state

    ## unpickling 
    def __setstate__ ( self , state ) :

        Files.__setstate__ ( self , state )
        self.__quick     = state.get('quick'      , False )
        self.e_list1     = state.get('e_list1'    , set() )
        self.chain       = ROOT.TChain( state['chain'] )
        for f in  self.files :   self.chain.Add ( f ) 

    ## check the content of the two trees 
    def check_trees ( self , tree1 , tree2 , the_file = '' ) :

        if tree1 and tree2 :
            
            branches1 = set ( tree1.branches () )                
            leaves1   = set ( tree1.leaves   () )
            
            branches2 = set ( tree2.branches () )                
            leaves2   = set ( tree2.leaves   () )

            if branches1 != branches2 :
                missing = list ( branches1 - branches2 )
                missing . sort ()
                extra   = list ( branches2 - branches1 )
                extra   . sort ()                    
                logger.warning ( "Tree('%s'): missing/extra branches %s/%s in %s" %  ( tree1.GetName() , missing , extra , the_file ) )
                
            if ( ( branches1 != leaves1 ) or ( branches2 != leaves2 ) ) and leaves1 != leaves2 :
                missing = list ( leaves1 - leaves2 )
                missing . sort ()
                extra   = list ( leaves2 - leaves1 )
                extra   . sort ()                    
                logger.warning ( "Tree('%s'): missing/extra leaves   %s/%s in %s" %  ( tree1.GetName() , missing , extra , the_file ) )

                
    ## the specific action for each file 
    def treatFile ( self, the_file ) :
        """Add the file to TChain
        """
        ## suppress Warning/Error messages from ROOT 
        from ostap.logger.utils import rootError
        with rootError() :
            
            tmp1 = ROOT.TChain ( self.chain .GetName() )
            tmp1.Add ( the_file )
            
            if  tmp1 :                
                self.check_trees ( tmp1 , self.chain , the_file )                
                Files.treatFile ( self     , the_file ) 
                self.chain .Add ( the_file )
            else : 
                self.e_list1.add ( the_file )
                if not self.silent : 
                    logger.warning ( "No/empty chain  '%s' in file '%s'" % ( self.chain .GetName() ,
                                                                             the_file )            ) 
                    
    ## printout 
    def __str__(self):
        from ostap.logger.utils import rootWarning
        with rootWarning() :
            nf = len ( self.files   )
            nc = len ( self.chain   )
            ne = len ( self.e_list1 )
        return "<#files: {}; Entries: {}>"              .format ( nf , nc ) if not self.e_list1 else \
               "<#files: {}; Entries: {}; No/empty: {}>".format ( nf , nc  ,  ne ) 


    ## conversion to boolean
    def __nonzero__ ( self )  :
        return bool(self.files) and bool(self.chain)

    ## merge two sets together
    def __iadd__ ( self , other ) :
        """ Merge two datasets
        >>> ds   = ...
        >>> ds1  = ...
        >>> ds  += ds1
        >>> ds  |= ds1 ## ditto 
        """
        if not isinstance ( other , Data ) : return NotImplemented
        return self.add_data ( other )    

    ## merge two sets together
    def add_data ( self , other ) :
        """ Merge two datasets
        """
        if self.chain.GetName() != other.chain.GetName() :
            raise NotImplementedError("Can't merge different chains %s and %s" % (
                self .chain.GetName() , 
                other.chain.GetName() , 
                ) )

        self.set_patterns ( self.patterns + other.patterns )
        for f in other.e_list1 : self.e_list1.add ( f ) 
        for f in other.files :
            if not f in self.files : 
                self.chain. Add ( f )
                
        return self
    
    ## get an intersection of two datasets 
    def __and__ (  self , other ) :
        """ get intersection of two sets
        >>> ds1 = ...
        >>> ds2 = ...
        >>> ds  = ds1 & ds2 ## get intsersecion
        >>> ds  = ds1 * ds2 
        """
        if not isinstance ( other , Data ) : return NotImplemented
        
        if self.chain.GetName() != other.chain.GetName() :
            raise NotImplementedError("Can't intersect different chains %s and %s" % (
                self .chain.GetName() , 
                other.chain.GetName() , 
                ) )
        
        result         = Data ( chain       = self.chain.GetName() ,
                                files       = []                   ,
                                description = self.description     , 
                                maxfiles    = self.maxfiles        ,
                                silent      = True                 )


        f1 = set ( self .files    )
        f2 = set ( other.files    )
        p1 = set ( self .patterns )
        p2 = set ( other.patterns )
        
        result.files    = list ( f1 & f2 )
        result.patterns = list ( p1 | p2 )
        result.e_list1  = self.e_list1 | other.e_list1
        result.silent   = self.silent
        
        for f in result.files : result.chain.Add ( f )
            
        return result
    
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
        result.set_files    ( deepcopy ( self.files    ) )
        result.set_patterns ( deepcopy ( self.patterns ) )
        result.e_list1      = deepcopy ( self.e_list1  ) 
        
        return result 
    
    ##  reload!
    def reload ( self , silent = True ) :
        self.chain   = ROOT.TChain ( self.chain.GetName() )
        self.e_list1 = set() 
        ##
        Files.reload ( self , silent ) 

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
        files = self.sample_files ( n , sort )
        return Data ( files       = files ,
                      chain       = kwargs.get ( 'chain'       , self.chain.GetName() ) , 
                      description = kwargs.get ( 'description' , self.description     ) ,
                      maxfiles    = kwargs.get ( 'maxfiles'    , self.maxfiles        ) ,
                      silent      = kwargs.get ( 'silent'      , self.silent          ) ,
                      quick       = kwargs.get ( 'quick'       , self.quick           ) )
    
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
        files = self.files[ item ]
        return Data ( files       = files                 ,
                      chain       = self.chain.GetName () , 
                      description = self.description      ,
                      maxfiles    = self.maxfiles         ,
                      silent      = self.silent           ,
                      quick       = self.quick            )

    # =========================================================================
    ## get DataFrame for the chain
    #  @see ROOT::RDataFrame
    @property
    def frame ( self ) :
        """Get the DataFrame for the chain"""
        from   ostap.frames.frames import DataFrame
        f = DataFrame ( self.chain )
        if not self.filent:
            pb = f.ProgressBar ( len ( self.chain ) )
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
    
    def __init__( self                  ,
                  chain1                ,
                  chain2                , 
                  files       = []      ,
                  description = ''      ,
                  maxfiles    = -1      ,
                  silent      = False   ,
                  quick       = False   ,
                  missing1st  = True    ,   ##  allow 1st missing chain 
                  missing2nd  = True    ) : ##  allow 2nd missing chain

        ## decorate files 
        if isinstance ( files , str ) : files = [ files ]

        if quick and 0 < maxfiles :
            logger.warning("``Quick'' and ``maxfiles>0'' are not compatible, switch to ``quick=False''")
            quick = False 
            
        if quick :
            quick = self._use_quick_(  files )
            if not quick :
                logger.warning("Patterns are not compatible with quick search, switch to ``quick=False''")

        self.e_list2 = set()
        self.chain2  = ROOT.TChain ( chain2 ) if isinstance ( chain2 , str ) else chain2 
        if not description :
            description = chain1.GetName() if hasattr ( chain1 , 'GetName' ) else str(chain1)
            description = "%s&%s" % ( description , self.chain2.GetName() )

        if quick and ( not missing1st or not missing2nd ) : 
            logger.warning ("Data2: Can't combine ``quick'' and ``missing''  options, set ``quick=False''") 
            quick = False

        self.missing1st = missing1st
        self.missing2nd = missing2nd
        
        self.__files2 = []
        if not quick : 
            Data.__init__( self , chain1 , files , description , maxfiles , silent , quick = False )
            self.chain1  = self.chain 
            self.set_files  ( self.chain .files()[:] ) 
            self.set_files2 ( self.chain2.files()[:] )
        else :
            Data.__init__( self , chain1 , []    , description , maxfiles , silent =  True , quick = True )
            self.chain1  = self.chain 
            self.silent  = silent
            self.set_files  ( self._quick_add_ ( self.chain  ,  files ) ) 
            self.set_files2 ( self._quick_add_ ( self.chain2 ,  files ) ) 
            self.set_files  ( self.chain .files()[:] ) 
            self.set_files2 ( self.chain2.files()[:] ) 
            if not self.silent :
                logger.info ('Loaded: %s' % self )

    @property 
    def files2    ( self ) :
        """``files2'' : the list of files"""
        return tuple ( self.__files2   ) 

    # =========================================================================
    ## set files 
    def set_files2 ( self , files ) :
        """Set files2"""
        assert isinstance ( files , ( list , tuple , set ) ),\
               'Invalid list of files %s' % files
        self.__files2 = list ( set ( files ) ) 
        self.__files2.sort()
        
    def __getstate__ ( self ) :

        state = Data.__getstate__( self )
        state [ 'e_list2'    ] = self.e_list2
        state [ 'chain2'     ] = self.chain2.GetName()
        state [ 'files2'     ] = self.__files2 
        state [ 'missing1st' ] = self.missing1st
        state [ 'missing2nd' ] = self.missing2nd
        
        return state

    def __setstate__ ( self , state ) :

        Data.__setstate__ ( self , state )
        
        self.e_list2     = state.get   ( 'e_list2'    , set() )
        self.chain2      = ROOT.TChain (  state['chain2'] )
        self.__files2    = state.get   ( 'files2'     , []    )
        self.missing1st  = state.get   ( 'missing1st' , True  )
        self.missing2nd  = state.get   ( 'missing2nd' , True  )
        
        for f in  self.files2 : self.chain2.Add ( f ) 
        
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

            if tmp1 : self.check_trees ( tmp1 , self.chain  , the_file )
            if tmp2 : self.check_trees ( tmp2 , self.chain2 , the_file )
  
            if  tmp1 and tmp2      : 
                Files.treatFile ( self     , the_file ) 
                self.chain .Add      ( the_file )
                self.chain2.Add      ( the_file )
                self.__files2.append ( the_file ) 
            elif tmp2 and not tmp1 and self.missing1st :
                self.e_list1.add ( the_file  )
                if not self.silent : 
                    logger.warning ( "No/empty chain1 '%s'      in file '%s'" % ( self.chain .GetName() ,
                                                                                  the_file )            )
                self.chain2.Add ( the_file )
                self.__files2.append ( the_file ) 
            elif tmp1 and not tmp2 and self.missing2nd :  
                self.e_list2.add ( the_file )
                if not self.silent : 
                    logger.warning ( "No/empty chain2 '%s'      in file '%s'" % ( self.chain2.GetName() ,
                                                                                  the_file )            ) 
                self.chain .Add ( the_file )            
                self.set_files  ( self.files + ( the_file , ) )  
            else :
                self.e_list1.add ( the_file )
                self.e_list2.add ( the_file )
                if not self.silent :                 
                    logger.warning ( "No/empty chains '%s'/'%s' in file '%s'" % ( self.chain .GetName() ,
                                                                                  self.chain2.GetName() ,
                                                                                  the_file )            )
    ## printout 
    def __str__(self):

        from ostap.logger.utils import rootWarning
        with rootWarning() :
            nf  = len( self.files   )
            nf2 = len( self.files2  )
            nc  = len( self.chain   )
            nc2 = len( self.chain2  )
            ne  = len( self.e_list1 )
            ne2 = len( self.e_list2 )
            
        sf  =  set(self.files) == set(self.files2)
        
        if not self.e_list1 and not self.e_list2 :
            return "<#files: {}; Entries: {}/{}>"   .format ( nf, nc , nc2 ) if sf else \
                   "<#files: {}/{}; Entries: {}/{}>".format ( nf , nf2 , nc , nc2 )
        else :
            return "<#files: {}; Entries: {}/{}; No/empty :{}/{}>"   .format ( nf , nc , nc2 , ne , ne2 ) if sf else \
                   "<#files: {}/{}; Entries: {}/{}; No/empty :{}/{}>".format ( nf , nf2 , nc , nc2 , ne , ne2 )
        
    def __nonzero__ ( self )  :
        return bool(self.files) and bool(self.files2) and bool(self.chain) and bool(self.chain2)

    ## merge two sets together
    def __iadd__ ( self , other ) :
        """ Merge two datasets
        >>> ds   = ...
        >>> ds1  = ...
        >>> ds  += ds1
        >>> ds  |= ds1 ## ditto 
        """
        if not isinstance ( other , Data2 ) : return NotImplemented        
        return self.add_data ( other )    

    ## merge two sets together
    def add_data ( self , other ) :
        """ Merge two datasets
        """
        if self.chain1.GetName() != other.chain1.GetName() :
            raise NotImplementedError("Can't merge different chains %s and %s" % (
                self .chain1.GetName() , 
                other.chain1.GetName() , 
                ) )
        
        if self.chain2.GetName() != other.chain2.GetName() :
            raise NotImplementedError("Can't merge different chains %s and %s" % (
                self .chain2.GetName() , 
                other.chain2.GetName() , 
                ) )

        self.set_patterns ( self.patterns + other.patterns )
        for f in other.e_list1 :self.e_list1.add ( f ) 
        for f in other.e_list2 :self.e_list2.add ( f )
        #
        for f in other.files :
            if not f in self.files : 
                self.chain1. Add   ( f )
                self._append_to_files ( f ) 
        for f in other.files2 :
            if not f in self.files2 : 
                self.chain2. Add   ( f )
                self.__files2.append ( f ) 
        return self 

    ## get an intersection of two datasets 
    def __and__ (  self , other ) :
        """ get intersection of two sets
        >>> ds1 = ...
        >>> ds2 = ...
        >>> ds  = ds1 & ds2 ## get intersection
        >>> ds  = ds1 * ds2 ## ditto
        """
        if not isinstance ( other , Data2 ) : return NotImplemented
        
        if self.chain1.GetName() != other.chain1.GetName() :
            raise NotImplementedError("Can't intersect different chains %s and %s" % (
                self .chain1.GetName() , 
                other.chain1.GetName() , 
                ) )
        
        if self.chain2.GetName() != other.chain2.GetName() :
            raise NotImplementedError("Can't intersect different chains %s and %s" % (
                self .chain2.GetName() , 
                other.chain2.GetName() , 
                ) )

        result         = Data2 ( chain1      = self.chain .GetName() ,
                                 chain2      = self.chain2.GetName() ,
                                 files       = []                    ,
                                 description = self.description      , 
                                 maxfiles    = self.maxfiles         ,
                                 silent      = True                  )
        
        f1  = set ( self .files    )
        f1o = set ( other.files    )
        f2  = set ( self .files2   )
        f2o = set ( other.files2   )
        p1  = set ( self .patterns )
        p2  = set ( other.patterns )
        
        result.set_files    ( f1 & f1o )
        result.set_files2   ( f2 & f2o )        
        result.set_patterns ( p1 | p2  )
        result.e_list1  = self.e_list1 | other.e_list1
        result.e_list2  = self.e_list2 | other.e_list2
        result.silent   = self.silent
        
        for f in result.files  : result.chain .Add ( f ) 
        for f in result.files2 : result.chain2 .Add ( f ) 
            
        return result

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
        for f in self.files  : result.chain .Add ( f ) 
        for f in self.files2 : result.chain2.Add ( f ) 
            
        result.chain1  = result.chain

        from copy import deepcopy
        result.set_files    ( self.files    ) 
        result.set_files2   ( self.files2   ) 
        result.set_patterns ( self.patterns ) 
        result.e_list1      = deepcopy ( self.e_list1    ) 
        result.e_list2      = deepcopy ( self.e_list2    ) 
        result.missing1st   = deepcopy ( self.missing1st ) 
        result.missing2nd   = deepcopy ( self.missing2nd ) 
        
        return result 

    ##  reload!
    def reload ( self , silent = True ) :
        self.__files2  = []
        ##
        self.chain2    = ROOT.TChain ( self.chain2.GetName() )
        self.e_list2   = set ()
        ##
        Data.reload ( self , silent ) 
        ##
        self.chain1  = self.chain 

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
        files = self.sample_files ( n , sort )
        return Data2 ( files       = files ,
                       chain1      = kwargs.get ( 'chain1'      , self.chain1.GetName() ) , 
                       chain2      = kwargs.get ( 'chain2'      , self.chain2.GetName() ) , 
                       description = kwargs.get ( 'description' , self.description      ) ,
                       maxfiles    = kwargs.get ( 'maxfiles'    , self.maxfiles         ) ,
                       silent      = kwargs.get ( 'silent'      , self.silent           ) ,
                       quick       = kwargs.get ( 'quick'       , self.quick            ) ,
                       missing1st  = kwargs.get ( 'missing1st'  , self.missing1st       ) ,
                       missing2nd  = kwargs.get ( 'missing2nd'  , self.missing2nd       ) )
    
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
        files = self.files[ item ]
        return Data2 ( files       = files                  ,
                       chain1      = self.chain1.GetName () , 
                       chain2      = self.chain2.GetName () , 
                       description = self.description       ,
                       maxfiles    = self.maxfiles          ,
                       silent      = self.silent            ,
                       quick       = self.quick             ,
                       missing1st  = self.missing1st        ,
                       missing2nd  = self.missing2nd        )
 
    # =========================================================================
    ## get DataFrame for the chain
    #  @see ROOT::RDataFrame
    @property
    def frame2 ( self ) :
        """Get the DataFrame for the chain"""
        from   ostap.frames.frames import DataFrame
        f = DataFrame ( self.chain2 )
        if not self.silent :
            pb = f.ProgressBar ( len ( self.chain2 ) ) 
        return f
    
# =============================================================================

# =============================================================================
if '__main__' == __name__ :

    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
