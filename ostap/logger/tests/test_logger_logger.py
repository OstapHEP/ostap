#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file ostap/logger/tests/test_logger_logger.py
# Copyright (c) Ostap developpers.
# =============================================================================
""" Test module or soem core utils 
"""
# =============================================================================
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'test_logger_logger'  )
else                       : logger = getLogger ( __name__          )
# =============================================================================
from   ostap.logger.logger import *


def test_logger_logger () :
    
    logger.verbose  ( 'This is VERBOSE  message'  ) 
    logger.debug    ( 'This is DEBUG    message'  ) 
    logger.info     ( 'This is INFO     message'  ) 
    logger.warning  ( 'This is WARNING  message'  ) 
    logger.error    ( 'This is ERROR    message'  ) 
    logger.fatal    ( 'This is FATAL    message'  ) 
    logger.critical ( 'This is CRITICAL message'  ) 
    
    with logColor() : 
        
        logger.verbose  ( 'This is VERBOSE  message'  ) 
        logger.debug    ( 'This is DEBUG    message'  ) 
        logger.info     ( 'This is INFO     message'  )
        logger.warning  ( 'This is WARNING  message'  ) 
        logger.error    ( 'This is ERROR    message'  ) 
        logger.fatal    ( 'This is FATAL    message'  ) 
        logger.critical ( 'This is CRITICAL message'  )
        
        with noColor () : 
            logger.verbose  ( 'This is VERBOSE  message'  ) 
            logger.debug    ( 'This is DEBUG    message'  ) 
            logger.info     ( 'This is INFO     message'  )
            logger.warning  ( 'This is WARNING  message'  ) 
            logger.error    ( 'This is ERROR    message'  ) 
            logger.fatal    ( 'This is FATAL    message'  ) 
            logger.critical ( 'This is CRITICAL message'  )
            
        logger.verbose  ( 'This is VERBOSE  message'  ) 
        logger.debug    ( 'This is DEBUG    message'  ) 
        logger.info     ( 'This is INFO     message'  )
        logger.warning  ( 'This is WARNING  message'  ) 
        logger.error    ( 'This is ERROR    message'  ) 
        logger.fatal    ( 'This is FATAL    message'  ) 
        logger.critical ( 'This is CRITICAL message'  )
        
        with keepColor() :
            
            logger.verbose  ( 'This is VERBOSE  message'  ) 
            logger.debug    ( 'This is DEBUG    message'  ) 
            logger.info     ( 'This is INFO     message'  ) 
            logger.warning  ( 'This is WARNING  message'  )
            
            make_colors()
            
            logger.error    ( 'This is ERROR    message'  ) 
            logger.fatal    ( 'This is FATAL    message'  ) 
            logger.critical ( 'This is CRITICAL message'  ) 

    logger.verbose  ( 'This is VERBOSE  message'  ) 
    logger.debug    ( 'This is DEBUG    message'  ) 
    logger.info     ( 'This is INFO     message'  ) 
    logger.warning  ( 'This is WARNING  message'  ) 
    logger.error    ( 'This is ERROR    message'  ) 
    logger.fatal    ( 'This is FATAL    message'  ) 
    logger.critical ( 'This is CRITICAL message'  ) 
    logger.info ( 80*'*'  )

# =============================================================================
if __name__ == '__main__' :

    test_logger_logger ()
    
# =============================================================================
##                                                                     The END 
# =============================================================================

