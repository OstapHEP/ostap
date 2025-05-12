#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file fixes.py 
#  Couple of minor fixes for Ostap
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-02-23
# =============================================================================
""" Couple of minor fixes for Ostap
"""
# =============================================================================
__version__ = '$Revision$'
__author__  = 'Vanya BELYAEV Ivan.Belyaev@itep.ru'
__date__    = '2016-02-23'
__all__     = () ## noting to import 
# =============================================================================
import os 
# ============================================================================
## @class MuteC
#  context manager to suppress pythion prinout
#  the actual code is stallen from
#  http://stackoverflow.com/questions/11130156/suppress-stdout-stderr-print-from-python-functions
#  A fix is added for "IOError: [Errno 24] Too many open files" :
#  original code leaks the file descriptors
class MuteC(object):
    """ A context manager for doing a ``deep suppression'' of stdout and stderr in 
    Python, i.e. will suppress all print, even if the print originates in a 
    compiled C/Fortran sub-function.
    This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).      
    
    stallen from  
    http://stackoverflow.com/questions/11130156/suppress-stdout-stderr-print-from-python-functions
    """
    #
    ## class variables: dev-null device & instance counter 
    _devnull = 0
    _cnt     = 0
    
    def __init__( self , out = True , err = False ):
        
        self._out = out
        self._err = err

        # increment instance counter 
        self.__class__._cnt += 1

        # create dev-null if not done yet 
        if not self.__class__._devnull :
            self.__class__._devnull = os.open ( os.devnull , os.O_WRONLY )            

    def __del__  ( self ) :
        
        # decrement instance counter 
        self.__class__._cnt -= 1
        
        # close dev-null if not done yet 
        if self.__class__._cnt <= 0 and self.__class__._devnull : 
            os.close ( self.__class__._devnull  )
            self.__class__._devnull = 0
            
    ## context-manager 
    def __enter__(self):
        
        ## Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds =  os.dup(1), os.dup(2)  # leak was here !!!
        
        ## mute it!
        if self._out : os.dup2 ( self.__class__._devnull , 1 )  ## C/C++
        if self._err : os.dup2 ( self.__class__._devnull , 2 )  ## C/C++

        return self
    
    ## context-manager 
    def __exit__(self, *_):
        
        # Re-assign the real stdout/stderr back to (1) and (2)  (C/C++)
        if self._err : os.dup2 ( self.save_fds[1] , 2 )
        if self._out : os.dup2 ( self.save_fds[0] , 1 )
        
        # fix the  file descriptor leak
        # (there were no such line in example, and it causes
        #      the sad:  "IOError: [Errno 24] Too many open files"
        
        os.close ( self.save_fds[1] ) 
        os.close ( self.save_fds[0] )

# =============================================================================
with MuteC ( True , True ) : 
    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions = True

# =============================================================================
# Include path for ACLiC:
# =============================================================================
opath = ROOT.gSystem.GetIncludePath()
## logger.debug ( 'Old include ath: %s' % opath  )
opath = opath.replace ( '-I' , ' ' ) . split ()
## add gsl ? 
opath.append ( '$OSTAPDIR/include' ) 

npath  = []
npath_ = []
for item in opath :

    if not item : continue
    
    if   item[0] == '"' and item[-1]== '"' : item = item[1:-1]
    elif item[0] == "'" and item[-1]== "'" : item = item[1:-1]
    
    nitem = os.path.expandvars (  item ) 
    nitem = os.path.expandvars ( nitem ) 
    nitem = os.path.expandvars ( nitem ) 
    nitem = os.path.expanduser ( nitem ) 
    nitem = os.path.expandvars ( nitem ) 

    if os.path.exists ( nitem ) and os.path.isdir ( nitem ) :
        if not item in npath and not nitem in npath_ :
            if ' ' in item : npath.append ( '"' + item + '"')
            else           : npath.append (       item      )
            npath_.append ( nitem ) 
            ## for CINT
            groot = ROOT.ROOT.GetROOT()
            groot.ProcessLine ('.include  %s' % nitem )

if npath : npath[0] = '-I'+npath[0]
npath = ' -I'.join(npath)
ROOT.gSystem.SetIncludePath( npath ) 
npath = ROOT.gSystem.GetIncludePath()
# =============================================================================
    
# =============================================================================
if '__main__' == __name__ : 

    from   ostap.logger.logger    import getLogger
    logger = getLogger ('ostap.logger.fixes')

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger ) 

# =============================================================================
##                                                                      The END
# =============================================================================
