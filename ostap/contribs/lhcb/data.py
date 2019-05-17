#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/contribs/lhcb/data.py 
#
#  Useful utilities to get access to datafiles & chains
#  Actualy it is just a little bit modified version (with globbing) of
#  original ``Ostap'' code by Alexander BARANOV
#
#  @code
#  data  = DataAndLumi('Bc/MyTree', '*.root' )
#  chain = data.chain
#  flist = data.files 
#  lumi  = data.lumi
#  print data.getLumi() 
#
#  data  = DataAndLumi('Bc/MyTree', [ 'A/*.root' , 'B/B/D/*.root' ] )
#  chain = data.chain
#  flist = data.files 
#  lumi  = data.lumi
#  print data.getLumi()
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
    'DataAndLumi' , ## collect files and create two TChain objects (one for Lumi)
    )
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger    import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.contribs.lhcb.data' )
else                       : logger = getLogger ( __name__     )
# =============================================================================
from copy                     import deepcopy
from ostap.logger.utils       import rootError
from ostap.trees.data_utils   import Data2 
from ostap.contribs.lhcb.lumi import getLumi
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

    def __init__( self                ,
                  chain               ,
                  files       = []    ,
                  description = ''    , 
                  lumi_chain  = 'GetIntegratedLuminosity/LumiTuple' , 
                  maxfiles    = 1000000                             ,
                  silent      = False ,
                  missing     = True  ) :  

        if not description :
            description = chain.GetName() if hasattr ( chain , 'GetName' ) else str(chain)
        Data2.__init__ ( self , chain , lumi_chain , files , description , maxfiles  , silent , quick = False , missing1st = missing , missing2nd = False ) 
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
        with rootError() :
            try :
                return getLumi ( self.chain2  )
            except ImportError :
                logger.error('DataAndLumi:getLumi is not available!')
                return -1.e+6
                
    ## printout 
    def __str__(self):
        
        l   = self.getLumi()
        nf  = len ( self.files   )
        nf2 = len ( self.files2  )
        nc  = len ( self.chain   )
        ne  = len ( self.e_list1 )
        ne2 = len ( self.e_list2 )
        
        sf  = set(self.files) == set(self.files2)
        
        if not self.e_list1 and not self.e_list2 :            
            return "<Luminosity: {}pb-1; #files: {}; Entries: {};>"   .format ( l , nf ,       nc ) if sf else \
                   "<Luminosity: {}pb-1; #files: {}/{}; Entries: {};>".format ( l , nf , nf2 , nc )
        else :            
            return "<Luminosity: {}pb-1; #files: {}; Entries: {}; No/empty: {}/{}>"   .format( l , nf ,       nc , ne , ne2 ) if sf else \
                   "<Luminosity: {}pb-1; #files: {}/{}; Entries: {}; No/empty: {}/{}>".format( l , nf , nf2 , nc , ne , ne2 )


    ## get an intersection of two datasets 
    def __and__ (  self , other ) :
        """ get intersection of two sets
        >>> ds1 = ...
        >>> ds2 = ...
        >>> ds  = ds1 & ds2 ## get intersection
        >>> ds  = ds1 & ds2 ## ditto
        """
        if not isinstance ( other , DataAndLumi ) : return NotImplemented
        
        result          = DataAndLumi (
            chain       = self.chain .GetName() ,
            lumi_chain  = self.chain2.GetName() ,
            files       = []                    ,
            description = self.description      , 
            maxfiles    = self.maxfiles         ,
            silent      = True                  )
        
        f1  = set ( self .files   )
        f1o = set ( other.files   )
        f2  = set ( self .files2  )
        f2o = set ( other.files2  )
        p   = set ( self .patterns )
        po  = set ( other.patterns )
        
        result.files    = list ( f1 & f1o )
        result.files2   = list ( f2 & f2o )        
        result.patterns = list ( p  | po  )
        result.e_list1  = self.e_list1 | other.e_list1
        result.e_list2  = self.e_list2 | other.e_list2
        result.silent   = self.silent
        
        for f in result.files  : result.chain .Add ( f ) 
        for f in result.files2 : result.chain2.Add ( f ) 
            
        return result



    ## clone it !
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
        
        for f in self.files  : result.chain .Add ( f ) 
        for f in self.files2 : result.chain2.Add ( f ) 
            
        result.chain1  = result.chain
        result.lumi    = result.chain2

        result.files    = deepcopy ( self.files    ) 
        result.files2   = deepcopy ( self.files2   ) 
        result.patterns = deepcopy ( self.patterns ) 
        result.e_list1  = deepcopy ( self.e_list1  ) 
        result.e_list2  = deepcopy ( self.e_list2  ) 
        
        return result 

    ##  reload!
    def reload ( self ) :
        self.files   = [] 
        self.files2  = [] 
        self.chain   = ROOT.TChain ( self.chain .GetName() )
        self.chain2  = ROOT.TChain ( self.chain2.GetName() )
        self.e_list1 = set () 
        self.e_list2 = set () 
        ## 
        self.add_files ( deepcopy ( self.patterns ) )
        ##
        self.chain1  = self.chain 
        self.lumi    = self.chain2 
 
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
        return DataAndLumi ( files       = files ,
                             chain       = kwargs.get ( 'chain'       , self.chain1.GetName() ) ,
                             lumi_chain  = kwargs.get ( 'lumi_chain'  , self.chain2.GetName() ) ,                             
                             description = kwargs.get ( 'description' , self.description      ) ,
                             maxfiles    = kwargs.get ( 'maxfiles'    , self.maxfiles         ) ,
                             silent      = kwargs.get ( 'silent'      , self.silent           ) ,
                             missing     = kwargs.get ( 'missing'     , self.missing1st       ) )

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
        return DataAndLumi ( files       = files                  ,
                             chain       = self.chain1.GetName () ,
                             lumi_chain  = self.chain2.GetName () ,                             
                             description = self.description       ,
                             maxfiles    = self.maxfiles          ,
                             silent      = self.silent            ,
                             missing     = self.missing1st        )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap import banner
    logger.info ( __file__ + '\n' + banner  )
    logger.info ( 80*'*' )
    logger.info ( __doc__  )
    logger.info ( 80*'*' )
    logger.info ( ' Author  : %s' %         __author__    ) 
    logger.info ( ' Version : %s' %         __version__   ) 
    logger.info ( ' Date    : %s' %         __date__      )
    logger.info ( ' Symbols : %s' %  list ( __all__     ) )
    logger.info ( 80*'*' ) 
    
# =============================================================================
# The END 
# =============================================================================
