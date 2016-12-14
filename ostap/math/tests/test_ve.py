# Copyright (c) Ostap developpers.
"""
Test module for ostap/math/ve.py.
"""

from ostap.utils.docme  import docme
from ostap.math.ve      import VE

from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.ve' )
else                       : logger = getLogger ( __name__        )

docme ( __name__ , logger = logger )

def test_ve():
    a = VE(100,100)
    b = VE(400,400)
    
    logger.info ( 'a=%s, b=%s' % ( a , b ) )
    
    logger.info ( 'a+b         %s' % ( a + b ) )
    logger.info ( 'a-b         %s' % ( a - b ) )
    logger.info ( 'a*b         %s' % ( a * b ) )
    logger.info ( 'a/b         %s' % ( a / b ) )
    logger.info ( 'a/(a+b)     %s' % ( a.frac ( b ) ) )
    logger.info ( '(a-b)/(a+b) %s' % ( a.asym ( b ) ) )
