#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/utils/scp_copy.py
#  Copy files using scp 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Decoration of Tree/Chain objects for efficient use in python"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (  
    'scp_copy' , ## copy files using scp 
  ) 
# =============================================================================
import ROOT, sys , os , time
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils.scp_copy' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( "Copy files using SCP")
# =============================================================================
from   ostap.core.core     import valid_pointer
import ostap.trees.cuts
# =============================================================================
_scp_copy = None
# =============================================================================
if sys.version_info.major < 3  :
    # try to import pathos 
    try :                
        # =================================================================
        import pathos.secure 
        # =================================================================
        ## copy <code>source</code> to <code>destination</code> using
        #  <code>pathos.secure.Copier</code
        def _scp_copy ( source , destination ) :
            """Copy ``source'' to ``destination'' using pathos.secure.Copier
            - see pathos.secure.Copier
            """
            copier = pathos.secure.Copier ( 'SSH-copier' )
            copier ( source = source , destination = destination ) 
            copier.launch ()
            r = copier.response()
            if r : logger.warning ('Response from SSH-copier: %s' % r )
            
            logger.debug ( 'pathos.secure.Copier will be used for scp-copy')
            
    except ImportError : ## pathos is not acessible, use subprocess 
        _scp_copy = None
# =============================================================================                
if not _scp_copy :
    # =========================================================================
    ## copy <code>source</code> to <code>destination</code> using
    #  <code>subprocess.check_call</code> + <code>scp</code>
    def _scp_copy ( source , destination ) :
        """Copy ``source'' to ``destination'' using
        subprocess.check_call & scp
        - see subprocess.check_call 
        """
        import subprocess                 
        try :
            r = subprocess.check_call ( "scp %s %s" % ( source , destination ) ,
                                        shell              = True ,
                                        close_fds          = True ,
                                        universal_newlines = True )            
        except subprocess.CalledProcessError as e :
            logger.error  ( str ( e ) ) 

# =====================================================================================
# Utility to copy files via ssh 
# =====================================================================================
## Copy remote file to local host
#  @code 
#  result , time = scp_to_local ( 'lxplus701.cern.ch:a.cpp' , 'b.cpp' )
#  @endcode
#  - Use either <code>pathos.secure.Copier</code> or <code>subprocess.check_call</code>
def scp_copy ( source  , destination = None ) :
    """Copy remote file to local host:
    
    >>> result , time = scp_copy ( source = 'lxplus701.cern.ch:a.cpp', destination='b.cpp' )
    
    - Use either pathos.secure.Copier or subprocess.check_call
    """
    
    if not destination :
        _ , fext = os.path.splitext ( source )  
        fext  = fext if fext else '.root'
        import ostap.utils.cleanup as CU
        destination = CU.CleanUp.tempfile ( suffix = fext , prefix = 'copy-' )
        
    start = time.time () 
    r = _scp_copy ( source , destination )
    stop  = time.time () 
    delta = stop - start
    
    return ( destination, delta )  if \
           os.path.exists ( destination ) and \
           os.path.isfile ( destination ) and \
           os.access      ( destination , os.R_OK ) else   ( None , delta ) 


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
