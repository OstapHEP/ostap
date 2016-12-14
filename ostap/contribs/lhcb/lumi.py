#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
# =============================================================================
## @file
#  Helper function to extract luminosity 
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
#  
#                    $Revision$
#  Last modification $Date$
#  by                $Author$
# =============================================================================
"""Helper function to extract luminosity 

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
    'getLumi'  ,  ## get the lumi
    )
# =============================================================================
import ROOT 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'Ostap.contrib.lhcb.lumi' )
else                       : logger = getLogger ( __name__        )
# =============================================================================
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
    tree_name = 'GetIntegratedLuminosity/LumiTuple'
    #
    from   ostap.core.core      import VE, hID
    import ostap.io.root_file
    import ostap.trees.trees 
    #
    if args :
        data = [ data ]
        for a in args : data.append ( a )
        return getLumi ( data ) 
        
    if isinstance ( data , str ) :
        ## expand the actual file name
        import os 
        data = os.path.expandvars ( data )
        data = os.path.expanduser ( data )
        data = os.path.expandvars ( data )
        data = os.path.expandvars ( data )
        
        try :    
            tree = ROOT.TChain ( tree_name ) 
            tree.Add ( data )   
            lumi = getLumi ( tree )
            return lumi
        except :
            logger.error('Unable to get lumi(1) for %s' % data )
            return VE()
        
        #
    if isinstance ( data , ROOT.TFile ) :
        try :
            tree = data.Get( tree_name ) 
            return getLumi ( tree ) 
        except:
            logger.error('Unable to get lumi(2) for %s' % data.GetName() )
            return VE()
        
    if isinstance ( data , ROOT.TTree ) :

        ## print data
        from ostap.logger.utils import rootError 
        
        try:
            #
            
            with rootError() : ## suppress errors from ROOT
                
                ## @attention here we are using sumVar! 
                l1 = data.sumVar ( '1.0*IntegratedLuminosity+0.0*IntegratedLuminosityErr' , '0<=IntegratedLuminosity' )
                l2 = data.sumVar ( '1.0*IntegratedLuminosity+1.0*IntegratedLuminosityErr' , '0<=IntegratedLuminosity' )
                l3 = data.sumVar ( '1.0*IntegratedLuminosity-1.0*IntegratedLuminosityErr' , '0<=IntegratedLuminosity' )            
                #
                l1.setError ( 0.5 * abs ( l2.value () - l3.value () ) )
                #
                l0 = data.sumVar ( 'IntegratedLuminosity' , '0 >IntegratedLuminosity'      )
                if 0 != l0.value() : logger.error( 'Something weird happens with Lumi/1: %s' % l0 )  
                l0 = data.sumVar ( 'IntegratedLuminosity' , 'IntegratedLuminosity>100000'  )
                if 0 != l0.value() : logger.error( 'Something weird happens with Lumi/2: %s' % l0 )  
                l0 = data.sumVar ( 'IntegratedLuminosity' , '0>IntegratedLuminosityErr'    )
                if 0 != l0.value() : logger.error( 'Something weird happens with Lumi/3: %s' % l0 )  
                # 
                return l1
        except :
            logger.error('Unable to get lumi(3) for %s' % data.GetName() )
            return VE()
        
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
