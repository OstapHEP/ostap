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
    'Lumi'        , ## collect files and create TChain object for Lumi
    'DataAndLumi' , ## collect files and create two TChain objects (one for Lumi)
    )
# =============================================================================
from copy                     import deepcopy
from ostap.core.core          import rootError
from ostap.trees.data_utils   import Data
from ostap.contribs.lhcb.lumi import getLumi
from ostap.utils.basic        import typename
from ostap.logger.symbols     import light_bulb as lumi_symbol 
from ostap.logger.symbols     import folder     as folder_symbol 
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger    import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.contribs.lhcb.data' )
else                       : logger = getLogger ( __name__     )
# =============================================================================
## @class  Lumi
#  Simple utility to collect luministy information form the standard trees 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-08  
class Lumi(Data):
    """ Simple utility to access to certain chain in the set of ROOT-files    
    >>> data  = DataAndLumi('Bc/MyTree', '*.root' )
    >>> chain = data.chain
    >>> flist = data.files 
    >>> lumi  = data.lumi
    >>> print data.getLumi() 
    """

    def __init__( self                 ,
                  files                , / , 
                  lumi        = 'GetIntegratedLuminosity/LumiTuple' , 
                  description = ''     , 
                  maxfiles    = 100000 ,
                  check       = True   ,
                  silent      = False  , 
                  parallel    = False  ) :  
        
        Data.__init__ ( self , files , lumi ,  
                        description = description ,
                        maxfiles    = maxfiles    ,
                        check       = check       , 
                        silent      = silent      , 
                        parallel    = parallel    )
    # =========================================================================
    @property
    def lumi ( self ) :
        """'lumi': recreate and return luminosity chain
        """
        return self.chain
    
    ## get the luminosity 
    def getLumi ( self ):
        """ Get the luminosity
        """
        ## suppress Warning/Error messages from ROOT 
        with rootError() : return getLumi ( self.chain  )
        
    @property
    def luminosity ( self ) :
        """'luminosity': get the Luminosity"""
        return self.getLumi() 
        
    ## printout 
    def __str__(self):
        """ The specific printout
        """
        l = self.getLumi() 
        result = '%s%s: %s/pb %d%s' % (
            lumi_symbol             ,
            typename ( self )       ,
            str ( l )               , 
            len      ( self.files ) , 
            file_symbol if file_symbol else 'files' )
        
        return result

# =============================================================================
## @class DataAndLumi
#  Simple utility to access to certain chain in the set of ROOT-files
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @author Alexander BARANOV a.baranov@cern.ch
#  @date   2014-06-08  
class DataAndLumi(Data):
    """ Simple utility to access to certain chain in the set of ROOT-files    
    >>> data  = DataAndLumi( '*.root' , 'Bc/MyTree' )
    >>> chain = data.chain
    >>> flist = data.files 
    >>> lumi  = data.lumi
    >>> print data.getLumi() 
    """

    def __init__( self                 ,
                  files                , / , 
                  *chains              ,
                  description = ''     ,
                  lumi        = 'GetIntegratedLuminosity/LumiTuple' , 
                  maxfiles    = 100000 ,
                  check       = True   ,
                  silent      = False  , 
                  parallel    = False  ) :  

        allchains = chains + ( lumi , ) 
        Data.__init__ ( self  , files , *allchains ,
                        description = description  ,
                        maxfiles    = maxfiles     ,
                        check       = check        , 
                        silent      = silent       ,
                        parallel    = parallel     )
        
    # =========================================================================
    @property
    def lumi ( self ) :
        """'lumi': recreate and return luminosity (the last) chain
        """
        return self.get_chain ( self.nchains - 1 ) 
    
    ## get the luminosity 
    def getLumi ( self ):
        """ Get the luminosity
        """
        ## suppress Warning/Error messages from ROOT 
        with rootError() : return getLumi ( self.lumi  )

    @property
    def luminosity ( self ) :
        """'luminosity': get the Luminosity"""
        return self.getLumi() 
        
    ## printout 
    def __str__(self):
        """ The specific printout
        """
        l = self.getLumi() 
        result = '%s%s(%s): %s/pb %d%s  %d%s' % (
            folder_symbol + lumi_symbol        ,
            typename ( self )                  , 
            ','.join ( self.chain_names[:-1] ) ,
            str ( l )                    , 
            len ( self.chain_names ) - 1 ,
            tree_symbol + chain_symbol if tree_symbol and chain_symbol else 'chains' ,            
            len ( self.files       ) , 
            file_symbol  if file_symbol  else 'files'  )
        
        return result 
    
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
##                                                                      The END 
# =============================================================================
