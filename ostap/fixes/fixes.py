#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file fixes.py 
#  Couple of minor fixes for Ostap
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-02-23
# =============================================================================
"""Couple of minor fixes for Ostap
"""
# =============================================================================
__version__ = '$Revision$'
__author__  = 'Vanya BELYAEV Ivan.Belyaev@itep.ru'
__date__    = '2016-02-23'
__all__     = () ## noting to import 
# =============================================================================
import ROOT, cppyy, os  
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.fixes.fixes')
else                      : logger = getLogger ( __name__           ) 
# =============================================================================
## gbl = cppyy.gbl 
## suppress welcome message from RooFit
## from ostap.logger.utils import mute
## with mute() :
##    v = ROOT.RooRealVar()
##    del v
    
## try :
##     enabled = ROOT.ROOT.IsImplicitMTEnabled ()
##     if enabled : logger.debug ("ImplicitMT is  enabled")
##     else       : logger.debug ("ImplicitMT is disabled")
## except AttributeError :
##     ROOT.ROOT.IsImplicitMTEnabled  = lambda *_ : False
##     ROOT.ROOT.EnableImplicitMT     = lambda *_ : False
##     ROOT.ROOT.DisableImplicitMT    = lambda *_ : False 
##     logger.info ("``Enable/Disable''Implicit MT is not available") 


print('ROOT:', ROOT, dir(ROOT) ) 

## # =============================================================================
## # Include path for ACLiC:
## # =============================================================================
## opath = ROOT.gSystem.GetIncludePath()
## logger.debug ( 'Old include ath: %s' % opath  )
## opath = opath.replace ( '-I' , ' ' ) . split ()
## ## add gsl ? 
## opath.append ( '$OSTAPDIR/include' ) 

## npath  = []
## npath_ = []
## for item in opath :

##     if not item : continue
    
##     if   item[0] == '"' and item[-1]== '"' : item = item[1:-1]
##     elif item[0] == "'" and item[-1]== "'" : item = item[1:-1]
    
##     nitem = os.path.expandvars (  item ) 
##     nitem = os.path.expandvars ( nitem ) 
##     nitem = os.path.expandvars ( nitem ) 
##     nitem = os.path.expanduser ( nitem ) 
##     nitem = os.path.expandvars ( nitem ) 

##     if os.path.exists ( nitem ) and os.path.isdir ( nitem ) :
##         if not item in npath and not nitem in npath_ :
##             if ' ' in item : npath.append ( '"' + item + '"')
##             else           : npath.append (       item      )
##             npath_.append ( nitem ) 
##             ## for CINT
##             groot = ROOT.ROOT.GetROOT()
##             groot.ProcessLine ('.include  %s' % nitem )

## print('ROOT:', ROOT, dir(ROOT) ) 
## if npath : npath[0] = '-I'+npath[0]
## npath = ' -I'.join(npath)
## ROOT.gSystem.SetIncludePath( npath ) 
## logger.debug ( 'New include path: %s' % npath  )
## # =============================================================================

    
# =============================================================================
if '__main__' == __name__ : 
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger ) 

# =============================================================================
##                                                                      The END
# =============================================================================
