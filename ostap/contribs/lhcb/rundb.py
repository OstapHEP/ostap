#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Collection of simple utilites to deal with LHCb RunDB 
#  via http://lbrundb.cern.ch/api
#
#  @thanks Alex PEARCE 
# 
#  Run this script to get example of run/fill information that could be
#  extracted using these utilities
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#
#                    $Revision$
#  Last modification $Date$
#  by                $Author$
# =============================================================================
"""Collection of simple utilites to deal with LHCb RunDB 
via http://lbrundb.cern.ch/api

Run this script to get example of run/fill information that could be
extracted using these utilities

Thanks to Alex PEARCE
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2015-01-12"
__all__     = (
    'run_info'   , ## get run information from RunDB  
    'fill_info'  , ## get fill information from RunDB 
    'fill_number', ## get fill number from given run-number 
    'run_url'    , ## pattern for run-information in LHCb RunDB 
    'fill_url'   , ## pattern for fill-information in LHCb RunDB 
    ) 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.contrib.lhcb.rundb' )
else                       : logger = getLogger ( __name__           )
# =============================================================================
logger.debug ( 'Collection of utilities to deal with LHCb RunDB')
# =============================================================================
import json,urllib
# =============================================================================
## local caches to minimize accesses to RunDB 
# =============================================================================
_fills_  = {} ## mapping: #run  -> #fill
_finfos_ = {} ## mapping: #fill -> fill-info
_rinfos_ = {} ## mapping: #run  -> run-info
# =============================================================================
## LHCb RunDB url format 
run_url  = 'http://lbrundb.cern.ch/api/run/{0}/'  ## pattern for runs  
fill_url = 'http://lbrundb.cern.ch/api/fill/{0}/' ## pattern for fills 
# =============================================================================
## get run info for the given run from LHCb runDB
#  @code
#  run  =  169064
#  rinfo = run_info ( run ) 
#  print 'Run Info %s ' %  rinfo
#  print 'Magnet: %s' % rinfo['magnetState']
#  print 'Velo  : %s' % rinfo['veloPosition']
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-01-12
#  @param run_number  run number
#  @return run information from LHCb RunDB 
#  @code
def run_info ( run_num ) :
    """Get run info for the given run from LHCb runDB
    >>> run  =  169064
    >>> rinfo = run_info ( run ) 
    >>> print 'Run Info %s ' %  rinfo 
    >>> print 'Magnet: %s' % rinfo['magnetState']
    >>> print 'Velo  : %s' % rinfo['veloPosition']
    """
    
    global _rinfos_
    rinfo  = _rinfos_.get ( run_num , None )
    if rinfo : return rinfo 
    
    try :
        
        #
        url   = run_url.format  ( run_num  )
        _obj  = urllib.urlopen  ( url      )
        rinfo = json.load       ( _obj     )

        rinfo = rinfo if rinfo else None
        _rinfos_ [ run_num ] = rinfo 
        return rinfo
        
    except:
        return None 

    return None

# =============================================================================
## get fill info for the given run from LHCb runDB
#  @code
#  fill  =  4691
#  finfo = fill_info ( fill ) 
#  print 'Fill Info %s ' %  finfo
#  print 'Colliding bunches : %s' % finfo ['nCollidingBunches']
#  print 'Beam-1 bunches    : %s' % finfo ['nBunchesB1']
#  print 'Beam-2 bunches    : %s' % finfo ['nBunchesB2']
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-01-12
#  @param fill_num  fill number
#  @return fill information from LHCb RunDB 
#  @code
def fill_info ( fill_num ) :
    """Get fill info for the given run from LHCb runDB
    >>> fill  =  4691
    >>> finfo = fill_info ( fill ) 
    >>> print 'Fill Info %s ' %  finfo
    >>> print 'Colliding bunches : %s' % finfo ['nCollidingBunches']
    >>> print 'Beam-1 bunches    : %s' % finfo ['nBunchesB1']
    >>> print 'Beam-2 bunches    : %s' % finfo ['nBunchesB2']    
    """
    
    global _finfos_
    finfo  = _finfos_.get ( fill_num , None )
    if finfo : return finfo 
    
    try :
        #
        url   = fill_url.format ( fill_num  )
        _obj  = urllib.urlopen  ( url      )
        finfo = json.load       ( _obj     )
        
        finfo = finfo if finfo else None 
        
        _finfos_ [ fill_num ] = finfo 
        return finfo
    
    except:
        return None 

    return None
    
# ===============================================================================
## get #fill from #run using LHCb RunDB
#  The function has been posted by Alex Pearce
#  at lhcb-davinci mailing list at 2015-01-12
#  @code
#  run  =  169064
#  fill = fill_number ( run )
#  print 'Run/Fill# %s/%s ' % ( run , fill
#  @endcode 
#  @author Alex PEARCE   alex.pearce@cern.ch
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-01-12
#  @param run_number  run number
#  @return fill number  or -1 
def fill_number ( run_number ) :
    """ Get #fill#  from #run via LHCb RunDB
    >>> run  =  169064
    >>> fill = fill_number ( run )
    >>> print 'Run/Fill# %s/%s ' % ( run , fill)
    """
    global _fills_
    fill_num = _fills_.get( run_number , None )
    if fill_num : return fill_num

    rinfo = run_info ( run_number )

    fill_num = rinfo.get ( 'fillid' , -1 ) if rinfo else -1
    _fills_ [ run_number ] = fill_num
    return  fill_num
    
# =============================================================================
if '__main__' == __name__  :
    
    from ostap import banner
    logger.info ( 80*'*' )        
    logger.info ( __file__ + '\n' + banner )
    logger.info ( 80*'*' )
    logger.info ( __doc__  )
    logger.info ( 80*'*' )
    logger.info ( ' Author  : %s' %         __author__    ) 
    logger.info ( ' Version : %s' %         __version__   ) 
    logger.info ( ' Date    : %s' %         __date__      )
    logger.info ( ' Symbols : %s' %  list ( __all__     ) )
    logger.info ( 80*'*' )    


    logger.info ( 80*'*' )        
    logger.info ( 'Run  info for #run=%d ' %  169064 )
    logger.info ( 80*'*' )        
    rinfo = run_info  ( 169064 )
    for k in rinfo : logger.info ( "  %30s : %-35s " % ( k , rinfo [ k ] ) )

    logger.info ( 80*'*' )        
    logger.info ( 'Fill info for #fill=%d ' % 4691 )
    logger.info ( 80*'*' )
    
    finfo = fill_info (   4691 ) 
    for k in finfo : logger.info ( "  %30s : %-35s " % ( k , finfo [ k ] ) )
    
    logger.info ( 80*'*' )    

    from ostap.utils.utils import timing    
    runs =  [ 0,  1 , 169064 , 5 , 6 , 98241980 , 169064  , 2334 , 2334 , 524387 ]
    for run in runs :
        with timing() :
            fill = fill_number ( run )
            logger.info ( 'Run/Fill#:%12d/%-10d ' % ( run , fill ) )
            
    logger.info ( 80*'*' )    

# =============================================================================
# The END 
# =============================================================================
