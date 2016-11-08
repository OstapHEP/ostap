#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
# =============================================================================
## @file
#
#  Simple file to provide "easy" access in python for
#  the basic ROOT::Math classes
#  @see Ostap/Point3DTypes.h
#  @see Ostap/Vector3DTypes.h
#  @see Ostap/Vector4DTypes.h
#  @see Ostap/GenericVectorTypes.h
#
#  The usage is fairly trivial:
#
#  @code
#
#  import ostap.math.base
#
#  @endcode
#
#  Important: All types are defined in corresponding C++ namespaces
#
#  @code
#
#  import ostap.math.base 
#
#  import cppyy
#  cpp   = cppyy.gbl                           ## global C++ namespace 
#  Ostap = cpp.Ostap                           ## get C++ namespace Ostap
#
#  p3 = Ostap.XYZPoint(0,1,2)               ## use C++ type Ostap::XYZPoint
#
#  @endcode
#
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-09-12
# =============================================================================
"""Simple file to provide 'easy' access in python for the basic ROOT::Math classes

  see Ostap/Point3DTypes.h
  see Ostap/Vector3DTypes.h
  see Ostap/Vector4DTypes.h
  see Ostap/GenericVectorTypes.h
  see Ostap/Line.h

  The lines and planes are decorated:
     see Ostap/GeomFun.h

  The usage is fairly trivial:

  >>> import ostap.math.base 

  Important: All types are defined in corresponding
               C++ namespaces: Gaudi & Gaudi::Math

  >>> import LHCbMath.Types
  >>> from GaudiPython.Bindings import gbl as cpp ## get global C++ namespace
  >>> Gaudi = cpp.Gaudi                           ## get C++ namespace Gaudi
  >>> p3 = Gaudi.XYZPoint(0,1,2)                  ## use C++ type Gaudi::XYZPoint

  >>> dir( Gaudi.Math )
  >>> dir( Gaudi      )

  Last modification $Date$
                 by $Author$

"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@nikhef.nl"
__date__    = "2009-09-12"
__version__ = "Version$Revision$"
# =============================================================================
__all__     = () ## nothing to be imported !
# =============================================================================
import ROOT, cppyy



# =============================================================================
# The END 
# =============================================================================
