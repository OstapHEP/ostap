#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file fixes.py 
#  Couple of minor fixes for Ostap
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-02-23
# =============================================================================
""" Couple of minor fixes for Ostap """
# =============================================================================
__version__ = '$Revision$'
__author__  = 'Vanya BELYAEV Ivan.Belyaev@cern.ch'
__date__    = '2016-02-23'
__all__     = () ## noting to import 
# =============================================================================
import os, sys 
# ============================================================================
## if ( 3 , 13 ) <= sys.version_info :
##    import builtins
##    sys.modules['__builtin__'] = builtins   
##    if not hasattr ( builtins , "__main__" ) :
##        import __main__
##        builtins.__main__ = __main__
# ============================================================================
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
groot = ROOT.ROOT.GetROOT()
if groot and groot.GetVersionInt() < 62800 :
    ROOT.gEnv.SetValue ( "RooFit.Banner" , 0 )    
    _ = ROOT.RooRealVar()

# =============================================================================
# Include path for ACLiC:
# =============================================================================
old_path     = ROOT.gSystem.GetIncludePath()
tokens       = old_path.replace ( '-I', ' ' ) .split()
unique_items = [ '$OSTAP_DIR/include' ]
seen_items   = set()

for token in tokens :
    
    item = token.strip ()
    if not item : continue
    if item [ 0 ] == item [ -1 ] :        
        if   item [ 0 ] == "'" : item = item [ 1 : -1 ]
        elif item [ 0 ] == '"' : item = item [ 1 : -1 ]

    item = os.path.expandvars ( item )
    item = os.path.expandvars ( item )
    item = os.path.expanduser ( item )
    item = os.path.expandvars ( item )
    item = os.path.normpath   ( item )
    
    if   not os.path.exists ( item ) : continue
    elif not os.path.isdir  ( item ) : continue
    
    if not item in seen_items  :
        seen_items .add ( item )
        
        if ' ' in item : unique_items.append ( '"%s"' % item )
        else           : unique_items.append (          item )

new_path = ' -I'.join ( unique_items )        
ROOT.gSystem.SetIncludePath( new_path )
# =============================================================================
    
# =============================================================================
if '__main__' == __name__ : 

    from   ostap.logger.logger    import getLogger
    logger = getLogger ('ostap.logger.fixes')

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger ) 

    new_path = ROOT.gSystem.GetIncludePath() 
    
    logger.info ( 'OLD include path:\n%s' % '\n'.join ( old_path.split() ) ) 
    logger.info ( 'NEW include path:\n%s' % '\n'.join ( new_path.split() ) ) 
    
# =============================================================================
##                                                                      The END
# =============================================================================
