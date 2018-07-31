#!/bin/env python
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
from ostap.math.rootfinder import find_root  ##  local utility

K    = 1
N    = 2 * K + 1

fun  = lambda x :   1.0     *           ( x - 1.0 )**N
der1 = lambda x :   1.0 * N *           ( x - 1.0 )**(N-1)
der2 = lambda x :   1.0 * N * ( N-1 ) * ( x - 1.0 )**(N-2)

kwargs = { 'full_output' : True , 'disp' : False }

print     'Halley:' , find_root ( fun , -0.5 , 10 , deriv1 = der1 , deriv2 = der2 , **kwargs )
print     'Newton:' , find_root ( fun , -0.5 , 10 , deriv1 = der1                 , **kwargs )
print     'Plain :' , find_root ( fun , -0.5 , 10                                 , **kwargs ) 

try :
    import scipy.optimize as SP
    print 'Brent :' , SP.brentq ( fun , -0.5 , 10 , **kwargs)
except :
    pass 



# =============================================================================
# The END 
# =============================================================================
