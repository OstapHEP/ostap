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
from copy                     import deepcopy
from ostap.core.core          import rootError
from ostap.trees.data_utils   import Data2 
from ostap.contribs.lhcb.lumi import getLumi
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger    import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.contribs.lhcb.data' )
else                       : logger = getLogger ( __name__     )
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

    def __init__( self                 ,
                  chain                ,
                  files       = []     ,
                  description = ''     , 
                  lumi_chain  = 'GetIntegratedLuminosity/LumiTuple' , 
                  maxfiles    = 100000 ,
                  check       = True   ,
                  silent      = False  ,
                  sort        = True   , 
                  parallel    = False  ) :  

        chain      =      chain if isinstance (      chain , str ) else      chain.name
        lumi_chain = lumi_chain if isinstance ( lumi_chain , str ) else lumi_chain.name
        
        if not description :
            description = "ROOT.TChain(%s)+lumi %s" % ( chain , lumi_chain )
            
        Data2.__init__ ( self                      ,
                         chain1      = chain       ,
                         chain2      = lumi_chain  ,
                         files       = files       ,
                         description = description ,
                         maxfiles    = maxfiles    ,
                         check       = check       , 
                         silent      = silent      ,
                         sort        = sort        , 
                         parallel    = parallel    ) 

    # =========================================================================
    @property
    def lumi ( self ) :
        """``lumi'': recreate and return luminosity chain (same as ``chain2'')
        """
        return self.chain2
    
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
    @property
    def luminosity ( self ) :
        """``luminosity'': get the Luminosity"""
        return self.getLumi() 
        
    ## printout 
    def __str__(self):
        
        l   = self.getLumi()
        nf  = len ( self.files      )
        nf2 = len ( self.files2     )
        nc  = len ( self.chain      )
        ne  = len ( self.bad_files  )
        
        sf  = set ( self.files ) == set ( self.files2 )
        
      #  if not self.e_list1 and not self.e_list2 :            
      #      return "<Luminosity: {}pb-1; #files: {}; Entries: {};>"   .format ( l , nf ,       nc ) if sf else \
      #             "<Luminosity: {}pb-1; #files: {}/{}; Entries: {};>".format ( l , nf , nf2 , nc )
      #  else :            
        return "<Luminosity: {}pb-1; #files: {}; Entries: {}; No/empty: {}>"   .format( l , nf ,       nc , ne ) if sf else \
                   "<Luminosity: {}pb-1; #files: {}/{}; Entries: {}; No/empty: {}>".format( l , nf , nf2 , nc , ne )

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
