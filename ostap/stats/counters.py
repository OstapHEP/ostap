#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
# =============================================================================
## @file
#  @see Ostap::StatEntity
#  @see Ostap::WStatEntity
#  @see Ostap::NStatEntity
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-09-12
# =============================================================================
"""Simple counters 
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@nikhef.nl"
__date__    = "2009-09-12"
__version__ = ""
# =============================================================================
__all__     = (
    'SE'  , ## simple smart counter (1-bi histo)  C++ Ostap::StatEntity 
    'SWE' , ## simple smart counter with weight :     Ostap::WStatEntity 
    'NWE' , ## simple smart running counter     :     Ostap::NStatEntity 
    ) 
# =============================================================================
import ROOT, cppyy
cpp   = cppyy.gbl
Ostap = cpp.Ostap

# =============================================================================
SE            = Ostap.StatEntity 
WSE           = Ostap.WStatEntity 
NSE           = Ostap.NStatEntity 

# =============================================================================
# temporary trick, to be removed 
# =============================================================================

SE .__repr__ = lambda s : 'Stat: '+ s.toString()
SE .__str__  = lambda s : 'Stat: '+ s.toString()

# =============================================================================
# minor decoration for StatEntity 
# ============================================================================= 
if not hasattr ( SE , '_orig_sum'  ) : 
    _orig_sum    = SE.sum
    SE._orig_sum = _orig_sum
    
if not hasattr ( SE , '_orig_mean' ) : 
    _orig_mean    = SE.mean
    SE._orig_mean = _orig_mean

from ostap.math.ve import VE 

SE. sum     = lambda s : VE ( s._orig_sum  () , s.sum2()       )
SE. minmax  = lambda s :    ( s.flagMin()     , s.flagMax()    ) 
SE. mean    = lambda s : VE ( s._orig_mean () , s.meanErr()**2 )

# =============================================================================
# minor decoration for WStatEntity 
# ============================================================================= 
if not hasattr ( WSE , '_orig_sum'  ) : 
    _orig_sum     = WSE.sum
    WSE._orig_sum = _orig_sum

if not hasattr ( WSE , '_orig_mean' ) : 
    _orig_mean_wse = WSE.mean
    WSE._orig_mean = _orig_mean_wse
    
WSE. sum     = lambda s : VE ( s._orig_sum  () , s.sum2()       )
WSE. mean    = lambda s : VE ( s._orig_mean () , s.meanErr()**2 )
WSE. minmax  = lambda s :            s.values  ().minmax() 
WSE.__repr__ = lambda s : 'WStat: '+ s.toString()
WSE.__str__  = lambda s : 'WStat: '+ s.toString()

# =============================================================================
if '__main__' == __name__  :
    
    print '*'*120
    print                      __doc__
    print ' Author  : %s ' %   __author__
    print ' Version : %s ' %   __version__
    print ' Date    : %s ' %   __date__
    print '*'*120
    
    import random
    
    cnt = SE() 
    for i in xrange(10000) :
        cnt += random.gauss(1,1)
        
    print 'Counter:', cnt
    print 'Mean : ', cnt.mean()
    print 'RMS  : ', cnt.rms ()
    
    print '*'*120
    
# =============================================================================
# The END 
# =============================================================================
