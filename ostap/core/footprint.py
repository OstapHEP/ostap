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
__all__     = () 
# =============================================================================
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
## central file, if exists and writeables 

# =============================================================================
import atexit 
@atexit.register 
def add_footprint () :
    # =========================================================================
    from ostap.core.meta_info import python_info, root_info, ostap_info, user 
    ## list of footprint-files 
    files = []
    ## (1) central file. common for everybody (is exists and writeable) 
    footprints_file = os.path.join ( '$OSTAPDIR' , '.footprints' )
    footprints_file = os.path.expanduser ( os.path.expandvars ( footprints_file ) )
    if os.path.exists ( footprints_file ) and \
       os.path.isfile ( footprints_file ) and \
       os.access      ( footprints_file , os.W_OK ) :
        files.append ( footprints_file )
    ## (2) user file in (writeable) cache directory
    from   ostap.core.cache_dir import cache_dir
    if cache_dir : files . append ( os.path.join ( cache_dir , '.footprints' ) )
    ## 
    nfp = len ( files ) 
    for i , fp_file in enumerate ( files ) : # ===============================
        # =====================================================================
        try : # ===============================================================
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
                if ( i + 1 == nfp ) : 
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
