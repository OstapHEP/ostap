#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file Utilities for getting the environmen varibakles 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
# =============================================================================
""" Utilities for getting the environmen varibakles 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    ## 
    'has_env'           , ## case-insensitive check for environment variable   
    'get_env'           , ## case-insensitive access to environment variable
    ##
    'OSTAP_BATCH'       , ## environment varianle for batch processing
    'OSTAP_BUILD_DIR'   , ## environment varianle for build-dir
    'OSTAP_DISPLAY'     , ## Ostap display 
    'OSTAP_CONFIG'      , ## Ostap config 
    'OSTAP_CACHE_DIR'   , ## Ostap cache directory
    'OSTAP_TMP_DIR'     , ## ostap TMP dir 
    'OSTAP_PROTOCOL'    , ## Ostap pickling protocol     
    'OSTAP_TABLE'       , ## Ostap table style
    'OSTAP_PARALLEL'    , ## Ostap parallel worker 
    'OSTAP_NCPUS'       , ## max number of paralell workers 
    ## 
)
# =============================================================================
import os
# =============================================================================
OSTAP_BATCH         = 'OSTAP_BATCH'
OSTAP_BUILD_DIR     = 'OSTAP_BUILD_DIR'
OSTAP_DISPLAY       = 'OSTAP_DISPLAY'
OSTAP_CONFIG        = 'OSTAP_CONFIG'
OSTAP_DIR           = 'OSTAP_DIR'
OSTAP_CACHE_DIR     = 'OSTAP_CACHE_DIR'   ## Ostap cache dir 
OSTAP_TMP_DIR       = 'OSTAP_TMP_DIR'     ## Ostap TMP   dir 
OSTAP_PROTOCOL      = 'OSTAP_PROTOCOL'    ## pickling protocol 
OSTAP_TABLE         = 'OSTAP_TABLE'       ## table style 
OSTAP_PARALLEL      = 'OSTAP_PARALLEL'    ## Ostap parallel worker 
OSTAP_NCPUS         = 'OSTAP_NCPUS'       ##Max number of parallel workers 
# =============================================================================
## transformation:  no blanks, no understores, no dashes 
#  - case-insensitive
#  - space ignored
#  - underline ignored
#  - dashes ignored 
transform = lambda v : v.replace(' ','').replace('_','').replace('-','').lower()
# =============================================================================
## case-insensitive check for existence of the environment variable
#  - case-insensitive
#  - space ignored
#  - underline ignored
#  @code
#  has_env ( 'Ostap_Table_Style' ) 
#  @endcode
def has_env ( variable ) :
    """ Case-insensitive check for existence of the environment variable    
    - case-insensitive
    - space ignored
    - underline ignored    
    - dashes ignored 
    >>> has_env ( 'Ostap_Table_Style' )    
    """
    new_var = transform ( variable )
    for key in os.environ  :
        if transform ( key ) == new_var  : return True
    return False
# =============================================================================
## case-insensitive access for the environment variable
#  - case-insensitive
#  - space ignored
#  - underline ignored
#  - dashes ignored 
#  @code
#  var = get_env ( 'Ostap_Table_Style' , '' ) 
#  @endcode
#  In ambiguous case warning message is printed and the last value is returned 
def get_env ( variable , default , silent = False ) :
    """ Case-insensitive access for the environment variable    
    - case-insensitive
    - space ignored
    - underline ignored    
    >>> var = has_env ( 'Ostap_Table_Style' , 'empty' )
    In ambiguous case  warning message is printed and the last value is returned 
    """
    new_var   = transform ( variable )
    found     = [] 
    for key in os.environ :
        if transform ( key ) == new_var  :
            value = os.environ [ key ]
            item  = key, value
            found.append  ( item )
            
    if not found : return default

    if 1 < len ( found )  and not silent :
        rows = [ ( 'Variable' , 'value' ) ]
        for k, v in found  :
            row = '%s' % k , '%s' % v
            rows.append ( row )
        title  = "'%s' matches" % variable
        import ostap.logger.table as T 
        table = T.table ( rows,  title = title  , prefix = '# ' , alignment = 'll' )
        from ostap.logger.logger import getLogger
        logger = getLogger ( 'ostap.get_env' ) 
        title2 = "Found %s matches for '%s', the last is taken" % ( len ( found ) , variable )  
        logger.warning ( '%s\n%s' %  ( title2 , table ) ) 
    
    return found [ -1] [ 1 ] 

# =============================================================================
## print enviroments as a table 
def show_env ( title = 'Environment' ) :
    """ Print enviroments as a table 
    """
    rows = [  ( '#' , 'Variable' , 'Value' ) ]

    i = 0
    keys = sorted ( os.environ.keys() )
    for key in keys :
        value = os.environ [ key ] 
        i += 1
        if 'PATH' in key :
            vv  = value.split(':')
            vv  =  '\n'.join ( vv ) 
            row = '%d' % i , key , vv 
        else : 
            row = '%d' % i , key , value 
        rows.append ( row )
    import ostap.logger.table as T
    return T.table ( rows , title = title , prefix = '# ' , alignment = 'rlw' )
    
# =============================================================================
if __name__ == '__main__' :

    from ostap.logger.logger import getLogger
    logger = getLogger ( 'ostap.utils.basic' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    title = 'Environment'
    logger.info ( '%s:\n%s' % ( title , show_env ( title = title ) ) )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
