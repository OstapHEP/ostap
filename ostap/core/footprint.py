#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file footprint.py
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-03-10
# =============================================================================
""" Ostap footproint 
""" 
# =============================================================================
__author__  = "Vanya BELYAEV  Ivan.Belyaev@itep.ru"
__date__    = "2014-03-10"
__version__ = "$Revision$"
__all__     = (
    'footprint_files' ,  
)
# =============================================================================
from ostap.core.meta_info import python_info, root_info, ostap_info, user 
from ostap.core.cache_dir import cache_dir
import os, sys, datetime
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.core.footprint'  )
else                       : logger = getLogger ( __name__                )
# =============================================================================
start_time = datetime.datetime.now()
# =============================================================================
#
footprint_files = []
## central file, if exists and writeables 
fpfile = '$OSTAPDIR/.footprints' 
fpfile_ = os.path.expanduser ( os.path.expandvars ( fpfile ) )
if os.path.exists ( fpfile_ ) and \
   os.path.isfile ( fpfile_ ) and \
   os.access      ( fpfile_ , os.W_OK ) : footprint_files.append ( fpfile_ )
## luser file in (sriteavle) cache directory
footprint_files .append ( os.path.join ( cache_dir , '.footprints' ) )
footprint_files = tuple ( footprint_files ) 

print ( footprint_files )

# =============================================================================
import atexit 
@atexit.register 
def add_footprint () :
    # =========================================================================
    for fp_file in footprint_files : # ========================================
        # =====================================================================
        try : # ================================================================
            # =================================================================
            end_time = datetime.datetime.now ()
            ts       = start_time.strftime ( '%Y-%m-%d %H:%M:%S' )
            te       = end_time  .strftime ( '%Y-%m-%d %H:%M:%S' )
            with open ( fp_file , 'at' ) as fp : # =========================
                fp.write ( " + ostap session %s\n" % os.getpid () )
                fp.write ( "   - started        : %s\n"            % ts          )
                fp.write ( "   - ended          : %s\n"            % te          )             
                fp.write ( "   - USER           : %s\n"            % user        )
                fp.write ( "   - CWD            : %s\n"            % os.getcwd() )
                fp.write ( "   - argv           : %s\n"            % sys.argv    )
                fp.write ( "   - ostap  version : %s.%s.%s.%s\n"   % ostap_info  )
                fp.write ( "   - ROOT   version : %s.%s/%s\n"      % root_info   )
                fp.write ( "   - python version : %s.%s.%s.%s%s\n" % python_info )
                logger.debug ( "Footprint is added to %s"          % fp_file ) 
            # =====================================================================
        except : # ================================================================
            pass 

# =============================================================================
if '__main__' == __name__ : 
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
