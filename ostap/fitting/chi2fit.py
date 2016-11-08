#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
# =============================================================================
## @file
#  @see Ostap::Math::Chi2Fit
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-09-12
# =============================================================================
"""Simple chi2-fit
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@nikhef.nl"
__date__    = "2009-09-12"
__version__ = ""
# =============================================================================
__all__     = (
    'C2FIT' , ## simple chi2-fit 
    ) 
# =============================================================================
import ROOT, cppyy
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.chi2fit' )
else                       : logger = getLogger ( __name__                )
# =============================================================================
from ostap.math.ve import VE,cpp 
cpp   = cppyy.gbl
Ostap = cpp.Ostap
C2FIT = Ostap.Math.Chi2Fit

C2FIT . __str__  = lambda s : s.toString ()
C2FIT . __repr__ = lambda s : s.toString ()
# =============================================================================
## chi2-probabilty
def _c2_prob_  ( s ) :
    """Chi2 probabiilty
    >>> r = h.hfit ( ... )
    >>> r.Prob()
    """
    dofs = s.points() - s.size()
    return ROOT.TMath.Prob ( s.chi2() , dofs )

C2FIT . Prob        = _c2_prob_
C2FIT . Probability = _c2_prob_
C2FIT . prob        = _c2_prob_
C2FIT . probability = _c2_prob_
C2FIT . __len__     = lambda s     : s.size  (   )
C2FIT . __getitem__ = lambda s , i : s.param ( i )

# =============================================================================
if '__main__' == __name__  :
    
    from ostap.logger.line import line 
    logger.info ( __file__  + '\n' + line  ) 
    logger.info ( 80*'*'   )
    logger.info ( __doc__  )
    logger.info ( 80*'*' )
    logger.info ( ' Author  : %s' %         __author__    ) 
    logger.info ( ' Version : %s' %         __version__   ) 
    logger.info ( ' Date    : %s' %         __date__      )
    logger.info ( ' Symbols : %s' %  list ( __all__     ) )
    logger.info ( 80*'*' ) 

    
# =============================================================================
# The END 
# =============================================================================
