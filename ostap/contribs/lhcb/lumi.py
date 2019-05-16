#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/contrib/lhcb/lumi.py
#  Helper function to extract luminosity (LHCb specific)
#
#  @code
#
#  >>> l1 = getLumi ( 'myfile.root' )
#  >>> l2 = getLumi ( tree  )
#  >>> l3 = getLumi ( chain )
#  >>> l4 = getLumi ( file  )
#  >>> l5 = getLumi ( [ any sequence of above ]  )
#
#  @endcode
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-10-16
# =============================================================================
"""Helper function to extract luminosity (LHCb specific)

Get lumi :

>>> l1 = getLumi ( 'myfile.root' )
>>> l2 = getLumi ( tree  )
>>> l3 = getLumi ( chain )
>>> l4 = getLumi ( file  )
>>> l5 = getLumi ( [ any sequence of above ]  )

"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2012-10-16"
# =============================================================================
__all__     = (
    'getLumi'     ,  ## get the lumi
    )
# =============================================================================
import ROOT, os  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger   import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.contribs.lhcb.lumi' )
else                       : logger = getLogger ( __name__        )
# ==============================================================================
import ostap.trees.trees 
import ostap.io.root_file
from   ostap.logger.utils  import rootError 
from   ostap.core.core     import VE, hID
# ==============================================================================
lumi_tree = 'GetIntegratedLuminosity/LumiTuple'
lumi      = 'IntegratedLuminosity'
lumi_err  = 'IntegratedLuminosityErr'
lumi_cuts = '0<=%s && 0<=%s' % ( lumi , lumi_err )
# ==============================================================================
## get luminosity from Lumi tuple
#
#  @param data  (INPUT) tree, chain, file, filename or sequence
#  @return the luminosity
#  @attention Linear addition of uncertainties is used here 
#
#  @code
#
#  >>> l1 = getLumi ( 'myfile.root' )
#  >>> l2 = getLumi ( tree  )
#  >>> l3 = getLumi ( chain )
#  >>> l4 = getLumi ( file  )
#  >>> l5 = getLumi ( [ any sequence of above ]  )
#
#  @endcode
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-09-20
def getLumi ( data , *args ) :
    """Get lumi :
    
    >>> l1 = getLumi ( 'myfile.root' )
    >>> l2 = getLumi ( tree  )
    >>> l3 = getLumi ( chain )
    >>> l4 = getLumi ( file  )
    >>> l5 = getLumi ( [ any sequence of above ]  )
    """
    #
    if args :
        data = [ data ]
        for a in args : data.append ( a )
        return getLumi ( data ) 
        
    if isinstance ( data , str ) :
        ## expand the actual file name
        data = os.path.expandvars ( data )
        data = os.path.expanduser ( data )
        data = os.path.expandvars ( data )
        data = os.path.expandvars ( data )
        
        try :    
            tree = ROOT.TChain ( lumi_tree ) 
            tree.Add ( data )   
            return getLumi ( tree )
        except :
            logger.error('Unable to get lumi/1 for %s' % data )
            return VE()
        
        #
    if isinstance ( data , ROOT.TFile ) :
        try :
            tree = data.Get( lumi_tree ) 
            return getLumi ( tree ) 
        except:
            logger.error('Unable to get lumi/2 for %s' % data.GetName() )
            return VE()
        
    if isinstance ( data , ROOT.TTree ) :

        ## try :
        with rootError() : ## suppress errors from ROOT
            
            ## if hasattr ( data , 'pstatVar' ) :
            ##    stat = data.pstatVar ( [ lumi , lumi_err ] , lumi_cuts , chunk_size = -1 , max_files = 10 )
            ## else  :
            stat = data. statVar ( [ lumi , lumi_err ] , lumi_cuts )

            ##
            s1 = stat[ lumi     ]
            s2 = stat[ lumi_err ]
            ##
            return VE ( s1.sum() , s2.sum() **2 )
        
        ##except :
        ##    logger.error('Unable to get lumi/3 for %s' % data.GetName() )
        ##    return VE()
        
    l = VE() 
    for i in data :
        k = getLumi ( i )
        ## @attention: linear addition of uncertainties: 
        l = VE ( l.value() + k.value() , ( l.error() + k.error () ) ** 2 ) 

    return l 


# =============================================================================
if '__main__' == __name__ :
    
    from ostap import banner
    logger.info ( __file__  + '\n' + banner ) 
    logger.info ( 80*'*'   )
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
