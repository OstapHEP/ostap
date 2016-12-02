#!/usr/bin/env ipython 
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$ 
# =============================================================================
## @file load_ostap.py
#  
#     .oooooo.                .                        
#    d8P'  `Y8b             .o8                        
#   888      888  .oooo.o .o888oo  .oooo.   oo.ooooo.  
#   888      888 d88(  "8   888   `P  )88b   888' `88b 
#   888      888 `"Y88b.    888    .oP"888   888   888 
#   `88b    d88' o.  )88b   888 . d8(  888   888   888 
#    `Y8bood8P'  8""888P'   "888" `Y888""8o  888bod8P' 
#                                            888       
#                                           o888o      
#                                                    
#  Simple interactive analysis environment to provide access to zillions
#  useful decorators for ROOT (and not only ROOT) objects&classes  
# 
#  This file is a part of 
#  <a href="http://cern.ch/lhcb-comp/Analysis/Bender/index.html">Bender project</a>
#  <b>``Python-based Interactive Environment for Smart and Friendly Physics Analysis''</b>
#
#  The package has been designed with the kind help from
#  Pere MATO and Andrey TSAREGORODTSEV. 
#  And it is based on the 
#  <a href="http://cern.ch/lhcb-comp/Analysis/LoKi/index.html">LoKi project:</a>
#  <b>``C++ ToolKit for Smart and Friendly Physics Analysis''</b>
#
#  By usage of this code one clearly states the disagreement 
#  with the smear campaign of Dr.O.Callot et al.: 
#  ``No Vanya's lines are allowed in LHCb/Gaudi software''
#
#  @date   2012-02-15
#  @author Vanya BELYAEV Ivan.Belyaevitep.ru
#
#                    $Revision$
#  Last modification $Date$
#                 by $Author$
# =============================================================================
""" Simple interactive PyRoot-based analysis environment
to provide access to zillions useful decorators for ROOT (and not only ROOT) objects&classes
    
    This file is a part of BENDER project:

  ``Python-based Interactive Environment for Smart and Friendly Physics Analysis''

The project has been designed with the kind help from Pere MATO and Andrey TSAREGORODTSEV. 

And it is based on the LoKi project:
 
   ``C++ ToolKit for Smart and Friendly Physics Analysis''

By usage of this code one clearly states the disagreement with the smear campaign 
of Dr.O.Callot et al.:

   ``No Vanya's lines are allowed in LHCb/Gaudi software''
"""
# =============================================================================
__author__  = 'Vanya BELYAEV Ivan.Belyaev@itep.ru'
__date__    = "2012-09-10"
__version__ = '$Revision$'
# =============================================================================
import ROOT, sys
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
__logger = getLogger ( 'ostap.core.load_ostap' )
# =============================================================================
## 1) load Ostap-style file
# =============================================================================
import ostap.plotting.style 
# =============================================================================
# The Heart 
# =============================================================================
## load zillions of decorations for ROOT-objects 
import ostap.core.pyrouts        ## NB: the most important line!
import ostap.io.zipshelve as DBASE

# =============================================================================
## minor decoration for default shelve module 
import ostap.io.shelve_ext

# =============================================================================
## import useful context managers
from ostap.logger.utils import *
from ostap.utils.utils  import *
# ============================================================================= 
## prepend the path 
if '.' not in sys.path :
    __logger.debug('Prepend sys.path with $PWD')
    sys.path = ['.'] + sys.path 
    
# =============================================================================
from ostap.core.core import cpp, Ostap, VE, SE, WSE, hID 
# =============================================================================
if '__main__' == __name__ :

    logger = __logger
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = __logger )

# =============================================================================
# The END 
# =============================================================================

