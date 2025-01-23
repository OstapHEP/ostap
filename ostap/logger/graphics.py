#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/logger/graphics.py 
#  Pseudographics symbols, picked from terminaltables 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2025-01-23
# =============================================================================
""" Pseudographics symbols, picked from terminaltables 
Module with decoration of TH* objects for efficient use in python
"""
# =============================================================================
__version__ = "$Revision: 207180 $"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
# =============================================================================
__all__     = ( 'CHAR_F_INNER_HORIZONTAL'      , 
                'CHAR_F_INNER_INTERSECT'       , 
                'CHAR_F_INNER_VERTICAL'        , 
                'CHAR_F_OUTER_LEFT_INTERSECT'  ,
                'CHAR_F_OUTER_LEFT_VERTICAL'   ,
                'CHAR_F_OUTER_RIGHT_INTERSECT' ,
                'CHAR_F_OUTER_RIGHT_VERTICAL'  ,
                'CHAR_H_INNER_HORIZONTAL'      ,
                'CHAR_H_INNER_INTERSECT'       ,
                'CHAR_H_INNER_VERTICAL'        ,
                'CHAR_H_OUTER_LEFT_INTERSECT'  ,
                'CHAR_H_OUTER_LEFT_VERTICAL'   ,
                'CHAR_H_OUTER_RIGHT_INTERSECT' ,
                'CHAR_H_OUTER_RIGHT_VERTICAL'  ,
                'CHAR_INNER_HORIZONTAL'        ,
                'CHAR_INNER_INTERSECT'         ,
                'CHAR_INNER_VERTICAL'          ,
                'CHAR_OUTER_BOTTOM_HORIZONTAL' , 
                'CHAR_OUTER_BOTTOM_INTERSECT'  ,
                'CHAR_OUTER_BOTTOM_LEFT'       ,  
                'CHAR_OUTER_BOTTOM_RIGHT'      , 
                'CHAR_OUTER_LEFT_INTERSECT'    , 
                'CHAR_OUTER_LEFT_VERTICAL'     , 
                'CHAR_OUTER_RIGHT_INTERSECT'   , 
                'CHAR_OUTER_RIGHT_VERTICAL'    , 
                'CHAR_OUTER_TOP_HORIZONTAL'    , 
                'CHAR_OUTER_TOP_INTERSECT'     , 
                'CHAR_OUTER_TOP_LEFT'          , 
                'CHAR_OUTER_TOP_RIGHT'         )
# =============================================================================    
import sys 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.logger.graphics' )
else                       : logger = getLogger( __name__                )
# =============================================================================
if 'win32' == sys.platform : # ================================================
    # =========================================================================
    CHAR_F_INNER_HORIZONTAL      = b'\xc4'.decode('ibm437')
    CHAR_F_INNER_INTERSECT       = b'\xc5'.decode('ibm437')
    CHAR_F_INNER_VERTICAL        = b'\xb3'.decode('ibm437')
    CHAR_F_OUTER_LEFT_INTERSECT  = b'\xc3'.decode('ibm437')
    CHAR_F_OUTER_LEFT_VERTICAL   = b'\xb3'.decode('ibm437')
    CHAR_F_OUTER_RIGHT_INTERSECT = b'\xb4'.decode('ibm437')
    CHAR_F_OUTER_RIGHT_VERTICAL  = b'\xb3'.decode('ibm437')
    CHAR_H_INNER_HORIZONTAL      = b'\xc4'.decode('ibm437')
    CHAR_H_INNER_INTERSECT       = b'\xc5'.decode('ibm437')
    CHAR_H_INNER_VERTICAL        = b'\xb3'.decode('ibm437')
    CHAR_H_OUTER_LEFT_INTERSECT  = b'\xc3'.decode('ibm437')
    CHAR_H_OUTER_LEFT_VERTICAL   = b'\xb3'.decode('ibm437')
    CHAR_H_OUTER_RIGHT_INTERSECT = b'\xb4'.decode('ibm437')
    CHAR_H_OUTER_RIGHT_VERTICAL  = b'\xb3'.decode('ibm437')
    CHAR_INNER_HORIZONTAL        = b'\xc4'.decode('ibm437')
    CHAR_INNER_INTERSECT         = b'\xc5'.decode('ibm437')
    CHAR_INNER_VERTICAL          = b'\xb3'.decode('ibm437')
    CHAR_OUTER_BOTTOM_HORIZONTAL = b'\xc4'.decode('ibm437')
    CHAR_OUTER_BOTTOM_INTERSECT  = b'\xc1'.decode('ibm437')
    CHAR_OUTER_BOTTOM_LEFT       = b'\xc0'.decode('ibm437')
    CHAR_OUTER_BOTTOM_RIGHT      = b'\xd9'.decode('ibm437')
    CHAR_OUTER_LEFT_INTERSECT    = b'\xc3'.decode('ibm437')
    CHAR_OUTER_LEFT_VERTICAL     = b'\xb3'.decode('ibm437')
    CHAR_OUTER_RIGHT_INTERSECT   = b'\xb4'.decode('ibm437')
    CHAR_OUTER_RIGHT_VERTICAL    = b'\xb3'.decode('ibm437')
    CHAR_OUTER_TOP_HORIZONTAL    = b'\xc4'.decode('ibm437')
    CHAR_OUTER_TOP_INTERSECT     = b'\xc2'.decode('ibm437')
    CHAR_OUTER_TOP_LEFT          = b'\xda'.decode('ibm437')
    CHAR_OUTER_TOP_RIGHT         = b'\xbf'.decode('ibm437')
    # =========================================================================
else : # ======================================================================
    # =========================================================================    
    CHAR_F_INNER_HORIZONTAL      = '\033(0\x71\033(B'
    CHAR_F_INNER_INTERSECT       = '\033(0\x6e\033(B'
    CHAR_F_INNER_VERTICAL        = '\033(0\x78\033(B'
    CHAR_F_OUTER_LEFT_INTERSECT  = '\033(0\x74\033(B'
    CHAR_F_OUTER_LEFT_VERTICAL   = '\033(0\x78\033(B'
    CHAR_F_OUTER_RIGHT_INTERSECT = '\033(0\x75\033(B'
    CHAR_F_OUTER_RIGHT_VERTICAL  = '\033(0\x78\033(B'
    CHAR_H_INNER_HORIZONTAL      = '\033(0\x71\033(B'
    CHAR_H_INNER_INTERSECT       = '\033(0\x6e\033(B'
    CHAR_H_INNER_VERTICAL        = '\033(0\x78\033(B'
    CHAR_H_OUTER_LEFT_INTERSECT  = '\033(0\x74\033(B'
    CHAR_H_OUTER_LEFT_VERTICAL   = '\033(0\x78\033(B'
    CHAR_H_OUTER_RIGHT_INTERSECT = '\033(0\x75\033(B'
    CHAR_H_OUTER_RIGHT_VERTICAL  = '\033(0\x78\033(B'
    CHAR_INNER_HORIZONTAL        = '\033(0\x71\033(B'
    CHAR_INNER_INTERSECT         = '\033(0\x6e\033(B'
    CHAR_INNER_VERTICAL          = '\033(0\x78\033(B'
    CHAR_OUTER_BOTTOM_HORIZONTAL = '\033(0\x71\033(B'
    CHAR_OUTER_BOTTOM_INTERSECT  = '\033(0\x76\033(B'
    CHAR_OUTER_BOTTOM_LEFT       = '\033(0\x6d\033(B'
    CHAR_OUTER_BOTTOM_RIGHT      = '\033(0\x6a\033(B'
    CHAR_OUTER_LEFT_INTERSECT    = '\033(0\x74\033(B'
    CHAR_OUTER_LEFT_VERTICAL     = '\033(0\x78\033(B'
    CHAR_OUTER_RIGHT_INTERSECT   = '\033(0\x75\033(B'
    CHAR_OUTER_RIGHT_VERTICAL    = '\033(0\x78\033(B'
    CHAR_OUTER_TOP_HORIZONTAL    = '\033(0\x71\033(B'
    CHAR_OUTER_TOP_INTERSECT     = '\033(0\x77\033(B'
    CHAR_OUTER_TOP_LEFT          = '\033(0\x6c\033(B'
    CHAR_OUTER_TOP_RIGHT         = '\033(0\x6b\033(B'
    # =========================================================================    


# =============================================================================
if '__main__' == __name__ :
    
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger ) 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
