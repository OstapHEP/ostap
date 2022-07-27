#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/core/known_issues.py
#  Known issues for ostap
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2020-09-30
# =============================================================================
"""Known issues"""
# =============================================================================
__version__ = "$Revision:$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2020-09-30"
__all__     = (
    'ROOT_issue_6470' , ## virtual C++ methods overriden in python are ignored for subclasses
    'ROOT_issue_6370' , ## dill cannot serialize ROOT objects for python3 
    'DILL_issue_356'  , ## dill cannot serialize ROOT objects for python3 
    'DILL_ROOT_issue' , ## dill cannot serialize ROOT objects for python3 
)
# =============================================================================
import ROOT, sys 
groot = ROOT.ROOT.GetROOT()
# =============================================================================
## dill has problems with serialization of ROOT objects for python3
#  @see https://github.com/root-project/root/issues/6370
#  @see https://github.com/uqfoundation/dill/issues/356
ROOT_issue_6370 = (3,6) <= sys.version_info and groot.GetVersionInt() < 62406 
DILL_issue_356  = (3,6) <= sys.version_info and groot.GetVersionInt() < 62406 
DILL_ROOT_issue = ROOT_issue_6370 or DILL_issue_356



# =============================================================================
## Virtual C++ methods overriden in python are ignored for subclasses 
#  @see https://github.com/root-project/root/issues/6470
ROOT_issue_6470   = groot.GetVersionInt() > 62300


# =============================================================================
if '__main__' == __name__ :
    
    # =========================================================================
    from ostap.logger.logger import getLogger 
    if '__main__' ==  __name__ : logger = getLogger( 'ostap.core.known_issues' )
    else                       : logger = getLogger( __name__     )
    # =========================================================================
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    issue = 0 
    if DILL_ROOT_issue :
        issue +=1 
        logger.info (" %d: DILL-ROOT-issue" %   issue ) 
        logger.info ('dill fails serialization of ROOT objects for 3.6<=python')
        logger.info ('It affects parallelization (pathos,multiprocess,ppft)'   )
        logger.info (' - see https://github.com/root-project/root/issues/6370' )
        logger.info (' - see https://github.com/uqfoundation/dill/issues/356'  )
        logger.info (80*'*')

    if ROOT_issue_6470 :
        issue +=1 
        logger.info (" %2d: Inheritance problem in PyROOT" % issue )
        logger.info ('C++ virtial methods overriden in python, ignored in subclasses')
        logger.info ('It affect `pyselectors` for *NEW* PyROOT' )
        logger.info (' - see https://github.com/root-project/root/issues/6470' )
        logger.info (80*'*')
        
# =============================================================================
##                                                                      The END 
# =============================================================================
    
