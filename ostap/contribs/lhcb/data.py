#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
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
#  >>> data  = DataAndLumi('Bc/MyTree', '*.root' )
#  >>> chain = data.chain
#  >>> flist = data.files 
#  >>> lumi  = data.lumi
#  >>> print data.getLumi()
#
#  @endcode
# 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @author Alexander BARANOV a.baranov@cern.ch
#  @date   2014-06-08
#
#                    $Revision$
#  Last modification $Date$
#  by                $Author$
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
    'DataAndLumi' )
# =============================================================================
import ROOT, glob 
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.contrib.lhcb.data' )
else                       : logger = getLogger ( __name__     )
# =============================================================================
import ostap.trees.data as DATA
Files = DATA.Files
Data  = DATA.Data 
Data2 = DATA.Data2 
            
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
        
    ## get the luminosity 
    def getLumi ( self ):
        """Get the luminosity
        """
        from   Ostap.GetLumi import getLumi
        return getLumi ( self.chain2  )

    ## printout 
    def __str__(self):
        
        return "<Luminosity: {}pb-1; #files: {}; Entries: {}>".format(
            self.getLumi()     ,
            len ( self.files ) ,
            len ( self.chain ) ) 


# =============================================================================
## add it to original module 
DATA.DataAndLumi = DataAndLumi 
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.logger.line import line 
    logger.info ( __file__ + '\n' + line  )
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
