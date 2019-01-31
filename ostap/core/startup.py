#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## The trivial startup sctript for Ostap session
#
#  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
#  @date   2006-10-08
#
# =============================================================================
"""This is a trivial startup script for Ostap session
 
"""
# =============================================================================
__author__  = 'Vanya BELYAEV Ivan.Belyaev@itep.ru'
__date__    = "2006-10-08"
__version__ = "$Revision$"
__all__     = () 
# =============================================================================
## logging
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' == __name__ : logger = getLogger ( 'ostap.core.startup' )
else                      : logger = getLogger ( __name__ )
# =============================================================================
from ostap.utils.basic import with_ipython
    
# =============================================================================
## check if the file is actually "empty"
def _empty_ ( fname ) :
    
    try :

        import os 
        if not os.path.exists ( fname ) or 0 == os.path.getsize( fname ) : return True
                
        ## find at least one non-empty line not starting from '#'
        with open ( fname , 'r' ) as f :
            for l in f :
                if not l  or  '#' == l [0] : continue
                l1 = l.strip()
                if     l1 and '#' != l1[0] : return False 
            
    except IOError : pass

    return True
     
# =============================================================================
import datetime, time 
start_time  = datetime.datetime.now()
start_time_ = time.time() 
logger.info ( 'Ostap  session started %s' % start_time.strftime('%c')  )

import os,sys
## define the name of the history file
__history__ = os.path.curdir + os.sep + '.ostap_history'

def _rename_ ( base , version = 0  ) :
    
    if version <= 0 :
        version = 0 
        fname   = '%s'    %    base 
        fnamep1 = '%s.%d' %  ( base , 1 )
    else :
        fname   = '%s.%d' %  ( base , version     )
        fnamep1 = '%s.%d' %  ( base , version + 1 )
        
    if os.path.exists ( fname ) and _empty_ ( fname ) :
        os.remove ( fname )
        return
    
    if os.path.exists ( fnamep1 ) :
        if  _empty_   ( fnamep1 ) : os.remove ( fnamep1 )
        else : _rename_ ( base , version + 1 )
        
    if os.path.exists ( fname ) :
        if _empty_    ( fname ) : os.remove ( fname )
        else : os.rename ( fname , fnamep1 )
        
                
# =============================================================================
## write history at the end 
def _prnt_() :
    end_time = datetime.datetime.now()
    delta    = time.time() - start_time_
    logger.info ( 'Ostap  session   ended %s [%.1fs]' %  ( end_time.strftime('%c') , delta ) )

# =============================================================================
## line completer 
import rlcompleter
import readline
readline.clear_history() 
readline.parse_and_bind("tab: complete")

# =============================================================================
## write history
def write_history ( fname ) :
    """Write history file 
    """
    ## first, delete old history file

    try :        
        _rename_ ( __history__  )
    except :
        logger.warning ( "Can't erase old history file(s)", exc_info = True ) 

    end_time = datetime.datetime.now()   
    command  = [ a for a in sys.argv ]
    if command : command[0] = os.path.basename( command[0] )
    command  = ' '.join(command)
    curdir   = os.getcwd() 
    
    delta    = time.time() - start_time_

    written  = False 
    if with_ipython() :

        try :
            
            import IPython, getpass
            ip  = IPython.get_ipython()
            me  = getpass.getuser()
            with open ( fname , 'w' ) as f :
                f.write( '# Ostap  session by %s started at %s\n' % ( me , start_time.strftime('%c' ) ) ) 
                f.write( '# Command from CWD=%s \n# %s\n' % ( curdir , command  ) ) 
                for record in ip.history_manager.get_range() :
                    f.write( record[2] + '\n' )
                f.write( '# Ostap  session by %s   ended at %s [%.1fs]\n' % ( me ,   end_time.strftime('%c' ) , delta ) )
                
            if os.path.exists( fname ) and os.path.isfile ( fname ) :
                if not _empty_ ( fname ) :
                    logger.info ( 'Ostap  history file: %s' % __history__ )
                    return                
            written = False            
        except:
            written = False 
            
    ## use 'old-style' history
    try :
        readline.write_history_file ( fname )
        written = True 
    except :
        written = False
    
    if written and os.path.exists( fname ) and os.path.isfile ( fname ) and not _empty_ ( fname ) : 
        logger.info ( 'Ostap  history file: %s' % __history__ )
    
# =============================================================================
import atexit
atexit.register ( _prnt_ )
atexit.register ( write_history , __history__ )    


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
