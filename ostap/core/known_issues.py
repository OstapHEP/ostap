#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/core/known_issues.py
#  Known issues for ostap
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2020-09-30
# =============================================================================
"""Known issues 
"""
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



# =============================================================================
## dill has problems with serialization of ROOT objects for python3
#  @see https://github.com/root-project/root/issues/6370
#  @see https://github.com/uqfoundation/dill/issues/356
ROOT_issue_6370 = 3 <= sys.version_info.major and 6 <= sys.version_info.minor
DILL_issue_356  = 3 <= sys.version_info.major and 6 <= sys.version_info.minor
DILL_ROOT_issue = ROOT_issue_6370 or DILL_issue_356


groot = ROOT.ROOT.GetROOT()

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
 
# =============================================================================
##                                                                      The END 
# =============================================================================
    
