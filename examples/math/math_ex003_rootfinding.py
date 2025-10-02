#!/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# Copyright Ostap developers
# =============================================================================
## @file examples/math/math_ex003_rootfinding.py
#  Simple example of (local) root-finding in ostap
# =============================================================================
"""Simple example of (local) root finding in Ostap 
"""
# =============================================================================
from   __future__          import print_function
from ostap.math.rootfinder import find_root  ##  local utility

K    = 1
N    = 2 * K + 1

fun  = lambda x :   1.0     *           ( x - 1.0 )**N
der1 = lambda x :   1.0 * N *           ( x - 1.0 )**(N-1)
der2 = lambda x :   1.0 * N * ( N-1 ) * ( x - 1.0 )**(N-2)

kwargs = { 'full_output' : True , 'disp' : False }

r1 = find_root ( fun , -0.5 , 10 , deriv1 = der1 , deriv2 = der2 , **kwargs )
r2 = find_root ( fun , -0.5 , 10 , deriv1 = der1                 , **kwargs )
r3 = find_root ( fun , -0.5 , 10                                 , **kwargs )

print ( 'Halley:%s' % str(r1) ) 
print ( 'Newton:%s' % str(r2) ) 
print ( 'Plain :%s' % str(r3) ) 

try :
    import scipy.optimize as SP
    r4 = SP.brentq ( fun , -0.5 , 10 , **kwargs) 
    print ( 'Brent :%s' % str(r4) )
except :
    pass 



# =============================================================================
# The END 
# =============================================================================
