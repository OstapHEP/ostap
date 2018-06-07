#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file pypdf.py
#  Very specific helper class PyPDF to implement "pythonic" PDF for RooFit
#  - Typical usage:
#  @code
#  import math
#  from ostap.fitting.basic import MASS
#  from ostap.fitting.pypdf import PyPDF 
#  class PyGauss(MASS,PyPDF) :
#       norm = 1.0 / math.sqrt ( 2 * math.pi )
#       def __init__ ( self         ,
#                      name         ,
#                      xvar         ,    // the obeservable 
#                      mean         ,    // parameter: mean-value 
#                      sigma        ,    // parameter: sigma
#                      pdf   = None ) :  // MANDATORY ARGUMENT
#           MASS .__init__ ( self , name      , xvar , mean , sigma ) 
#           PyPDF.__init__ ( self , self.name , ( self.xvar  ,
#                                                 self.mean  ,
#                                                 self.sigma ) , pdf = pdf )
#           ## MANDATORY: 
#           self.config = {  
#               'name'  : self.name ,
#               'xvar'  : self.xvar ,
#               'mean'  : self.mean ,
#               'sigma' : self.mean ,
#               'pdf'   : None        ## ATTENTION!
#              }
#        ## MANDATORY:
#        def evaluate ( self ) : ## The main method 
#           varlist = self.varlist
#           x  = float ( varlist[0] ) 
#           m  = float ( varlist[1] ) 
#           s  = float ( varlist[2] ) 
#           dx = ( x - m ) / x
#           return math.exp ( -0.5 * dx * dx ) * self.norm / s
#   @endcode
#
# Note:
#  - 1. The double inheritance pattern: It is not mandatory, but it allows to
#       get all benefits from the left bas-class, <code>MASS</code> in this  case
#  - 2. The <code>__init__</code> method *must* have keyword argument
#       <code>pdf</code> (default is <code>None</code>)
#  - 3. One *must* specify the <code>config</code> dictionary
#
#  The latter two items are mandatory for the proper implemtation of RooAbsPdf::clone
#  - @see Ostap::Models::PyPdf 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2018-06-07
# =============================================================================
"""Very specific helper class to implement ``pythonic'' PDF for RooFit

    Typical usage:
    
    ... import math
    ... from ostap.fitting.basic import MASS
    ... from ostap.fitting.pypdf import PyPDF 
    ... class PyGauss(MASS,PyPDF) :
    ...    norm = 1.0 / math.sqrt ( 2 * math.pi )
    ...    def __init__ ( self         ,
    ...                   name         ,
    ...                   xvar         ,    // the obeservable 
    ...                   mean         ,    // parameter: mean-value 
    ...                   sigma        ,    // parameter: sigma
    ...                   pdf   = None ) :  // MANDATORY ARGUMENT
    ...        MASS .__init__ ( self , name      , xvar , mean , sigma ) 
    ...        PyPDF.__init__ ( self , self.name , ( self.xvar  ,
    ...                                              self.mean  ,
    ...                                              self.sigma ) , pdf = pdf )
    ...        ## MANDATORY: 
    ...        self.config = {  
    ...          'name'  : self.name ,
    ...          'xvar'  : self.xvar ,
    ...          'mean'  : self.mean ,
    ...          'sigma' : self.mean ,
    ...          'pdf'   : None        ## ATTENTION!
    ...         }
    ...
    ...    ##MANDATORY:
    ...    def evaluate ( self ) :  ## the main method 
    ...        varlist = self.varlist     ## <--- ATTENTION! get variables 
    ...        x  = float ( varlist[0] ) 
    ...        m  = float ( varlist[1] ) 
    ...        s  = float ( varlist[2] ) 
    ...        dx = ( x - m ) / s
    ...        return math.exp ( -0.5 * dx * dx ) * self.norm / s 
    
    Note:
    
    1. the double inhetiance pattern: It is not mandatory, but it allows to
    get all benefits from the left bas-class, MASS in this  case
    2. The __init__ method *must* have keyword argument `pdf` (default is `None`)
    3. one *must* specify the `config` dictionary
    
    The latter two items are mandatory for the proper implemtation of RooAbsPdf::clone
     - see Ostap::Models::PyPdf 
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'PyPDF' , ## ``pythonic'' PDF for RooFit 
    )
# =============================================================================
import ROOT, math
from   ostap.core.core      import Ostap
import ostap.fitting.roofit 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.pypdf' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
## @class PyPDF
#  helper base class to implement ``pure-python'' PDF
#  - Typical usage:
#  @code
#  import math 
#  from ostap.fitting.basic import MASS
#  from ostap.fitting.pypdf import PyPDF 
#  class PyGauss(MASS,PyPDF) :
#       norm = 1.0 / math.sqrt ( 2 * math.pi )
#       def __init__ ( self         ,
#                      name         ,
#                      xvar         ,    // the obeservable 
#                      mean         ,    // parameter: mean-value 
#                      sigma        ,    // parameter: sigma
#                      pdf   = None ) :  // MANDATORY ARGUMENT
#           MASS .__init__ ( self , name      , xvar , mean , sigma ) 
#           PyPDF.__init__ ( self , self.name , ( self.xvar  ,
#                                                 self.mean  ,
#                                                 self.sigma ) , pdf = pdf )
#           ## MANDATORY: 
#           self.config = {  
#               'name'  : self.name ,
#               'xvar'  : self.xvar ,
#               'mean'  : self.mean ,
#               'sigma' : self.mean ,
#               'pdf'   : None        ## ATTENTION!
#              }
#        ## MANDATORY:
#        def evaluate ( self ) : ## The main method 
#           varlist = self.varlist
#           x  = float ( varlist[0] ) 
#           m  = float ( varlist[1] ) 
#           s  = float ( varlist[2] ) 
#           dx = ( x - m ) / s
#           return math.exp ( -0.5 * dx * dx ) * self.norm / s
#   @endcode
#
# Note:
#  - 1. The double inheritance pattern: It is not mandatory, but it allows to
#       get all benefits from the left bas-class, <code>MASS</code> in this  case
#  - 2. The <code>__init__</code> method *must* have keyword argument
#       <code>pdf</code> (default is <code>None</code>)
#  - 3. One *must* specify the <code>config</code> dictionary
#
#  The latter two items are mandatory for the proper implemtation of RooAbsPdf::clone
#  - @see Ostap::Models::PyPdf 
class PyPDF (object) :
    """Helper base class to implement ``pure-python'' PDF
    
    Typical usage:
    
    ... import math
    ... from ostap.fitting.basic import MASS
    ... from ostap.fitting.pypdf import PyPDF 
    ... class PyGauss(MASS,PyPDF) :
    ...    norm = 1.0 / math.sqrt ( 2 * math.pi )
    ...    def __init__ ( self         ,
    ...                   name         ,
    ...                   xvar         ,    // the obeservable 
    ...                   mean         ,    // parameter: mean-value 
    ...                   sigma        ,    // parameter: sigma
    ...                   pdf   = None ) :  // MANDATORY ARGUMENT
    ...        MASS .__init__ ( self , name      , xvar , mean , sigma ) 
    ...        PyPDF.__init__ ( self , self.name , ( self.xvar  ,
    ...                                              self.mean  ,
    ...                                              self.sigma ) , pdf = pdf )
    ...        ## MANDATORY: 
    ...        self.config = {  
    ...          'name'  : self.name ,
    ...          'xvar'  : self.xvar ,
    ...          'mean'  : self.mean ,
    ...          'sigma' : self.mean ,
    ...          'pdf'   : None        ## ATTENTION!
    ...         }
    ...
    ...    ##MANDATORY:
    ...    def evaluate ( self ) :  ## the main method 
    ...        varlist = self.varlist     ## <--- ATTENTION! get variables 
    ...        x  = float ( varlist[0] ) 
    ...        m  = float ( varlist[1] ) 
    ...        s  = float ( varlist[2] ) 
    ...        dx = ( x - m ) / s
    ...        return math.exp ( -0.5 * dx * dx ) * self.norm / s 
    
    Note:
    
    1. the double inhetiance pattern: It is not mandatory, but it allows to
    get all benefits from the left bas-class, MASS in this  case
    2. The __init__ method *must* have keyword argument `pdf` (default is `None`)
    3. one *must* specify the `config` dictionary
    
    The latter two items are mandatory for the proper implemtation of RooAbsPdf::clone
     - see Ostap::Models::PyPdf 
    """
    def __init__ ( self          ,
                   name          ,
                   vars          , ## all variables (observables&parameters) 
                   title  = ''   ,
                   pdf    = None ) :

        ## convert to RooArgList if needed 
        if not isinstance ( vars , ROOT.RooAbsCollection ) :
            vv = ROOT.RooArgList ()
            for v in vars :
                if not v in vv :  vv.add ( v )
            vars = vv 
            
        self.__pyvars = vars
        
        self.__pyname   =  name 
        if not title : title = "PyPDF(%s)"  % name 
        
        assert ( not pdf ) or isinstance ( pdf , Ostap.Models.PyPdf  ) , \
               "Invalid type of ``pdf'': %s/%s" % ( pdf , type ( pdf ) )
        
        
        if not pdf : pdf = Ostap.Models.PyPdf ( self , name , title , self.__pyvars  )

        logger.debug ( 'PyPDF: use %s/%s as PDF' % (  pdf , type(pdf) ) )
        
        self.pdf = pdf
    
    @property
    def vars ( self ) :
        """``vars'' - all variables for PyPDF
        - it should *NOT* be used innside `evaluate`-function!
        """
        return self.__pyvars
    
    @property
    def varlist ( self ) :
        """``varlist'' : get all variables (as RooArgList)  from Ostap::Models::PyPdf
        - all variables from evalauet function *MUST* be accessed via this collection
        """
        return self.pdf.varlist ()
    
    @property
    def varset ( self ) :
        """``varset''  : get all variables (as RooArgSet)   from Ostap::Models::PyPdf
        - all variables from evalauet function *MUST* be accessed via this collection
        """
        return self.pdf.varset () 

    ## clone the object (needed by C++ partner class)
    def clone ( self , **kwargs ) :
        """Clone the object
        -  attention: existing PDF is ``copied'', unless specified via kwargs (by C++)
        """
        conf =  {
            'name'    : self.__pyname ,
            'vars'    : self.vars     ,
            'pdf'     : self.pdf      ## ATTENTION: clone also PDF!
            }
        
        conf.update ( kwargs )
        
        KLASS = self.__class__
        return KLASS ( **conf )
 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
# The END 
# =============================================================================
