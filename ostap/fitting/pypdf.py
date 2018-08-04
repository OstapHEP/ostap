#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/pypdf.py
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
#  @endcode
#
# Note:
#  - 1. The double inheritance pattern: It is not mandatory, but it allows to
#       get all benefits from the left base-class, <code>MASS</code> in this  case
#  - 2. The <code>__init__</code> method *must* have keyword argument
#       <code>pdf</code> (default is <code>None</code>)
#  - 3. One *must* specify the <code>config</code> dictionary
#
#  The latter two items are mandatory for the proper implemtation of RooAbsPdf::clone
#  - @see Ostap::Models::PyPdf
#
#
# Also analytical integrals can be specified using two methods:
#  - <code>get_analytical_integral</code>, python' partner of <code>RooAbsPdf::getAnalyticalIntegral</code>
#  - <code>analytical_integral</code>, python' partner of <code>RooAbsPdf::analyticalIntegral</code>
#  @code
#    ## declare analytical integral 
#    def get_analytical_integral ( self ) :
#        """Declare the analytical integral"""
#        
#        x  = self.varlist[0]
#        
#        if self.matchArgs ( x ) : return 1 ## get the integration code
#        
#        return 0
#    
#    ## calculate the analytical integral 
#    def analytical_integral ( self ) :
#        """Calculate the analytical integral"""
#
#        assert 1 == self.intCode , 'Invalid integration code!'
#        
#        vlist = self.varlist
#
#        rn    = self.rangeName        
#        xv    = vlist [ 0 ]        
#        xmax  = xv.getMax ( rn )
#        xmin  = xv.getMin ( rn )
#        
#        m     = float ( vlist [ 1 ] ) 
#        s     = float ( vlist [ 2 ] )        
#        
#        return CDF ( xmax , m , s  ) - CDF ( xmin , m , s  )
# @endcode
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
    
    1. the double inheritance pattern: It is not mandatory, but it allows to
    get all benefits from the left base-class, MASS in this  case
    2. The __init__ method *must* have keyword argument `pdf` (default is `None`)
    3. one *must* specify the `config` dictionary
    
    The latter two items are mandatory for the proper implemtation of RooAbsPdf::clone
     - see Ostap::Models::PyPdf

     Also analytical integrals can be specified using two methods:
     - get_analytical_integral : python' partner of RooAbsPdf::getAnalyticalIntegral
     - analytical_integral     : python' partner of RooAbsPdf::analyticalIntegral
     ... ## declare analytical integral 
     ... def get_analytical_integral ( self ) :
     ...     x  = self.varlist[0]
     ...     if self.matchArgs ( x ) : return 1 ## get the integration code
     ...     return 0
     ... ## calculate the analytical integral 
     ... def analytical_integral ( self ) :
     ...     assert 1 == self.intCode , 'Invalid integration code!'
     ...     vlist = self.varlist
     ...     rn    = self.rangeName        
     ...     xv    = vlist [ 0 ]        
     ...     xmax  = xv.getMax ( rn )
     ...     xmin  = xv.getMin ( rn )
     ...     m     = float ( vlist [ 1 ] ) 
     ...     s     = float ( vlist [ 2 ] )        
     ... return CDF ( xmax , m , s  ) - CDF ( xmin , m , s  )
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
        """``vars'' - all variables used fro PyPDF creation
        - it should *NOT* be used inside `evaluate`-function!
        """
        return self.__pyvars
    
    @property
    def varlist ( self ) :
        """``varlist'' : get all variables (as RooArgList)  from Ostap::Models::PyPdf
        - all variables *MUST* be accessed via this collection!
        """
        return self.pdf.varlist ()
    
    @property
    def varset ( self ) :
        """``varset''  : get all variables (as RooArgSet)   from Ostap::Models::PyPdf
        - all variables *MUST* be accessed via this collection!
        """
        return self.pdf.varset () 

    @property
    def allDeps ( self ) :
        """``allDeps'' - the first parameter  of ``getAnalyticalIntegra/get_analytical_integral''-method 
        - attention:  it coudl be a null pointer! check it before usage! 
        """
        return self.pdf.allDeps()
    
    @property
    def analDeps ( self ) :
        """``analDeps'' - the second parameter  of ``getAnalyticalIntegral/get_analystical_integral''-method 
        - attention:  it could be a null pointer! check it before usage!
        """
        return self.pdf.allDeps()
    @property
    
    def rangeName ( self )  :
        """``rangeName'' -  the ``rangeName''  parameter for ``(get_)analytical_integral''-method
        - attention - it could be a null pointer! check it before usage! 
        """
        return self.pdf.rangeName ()
    
    @property
    def intCode   ( self )  :
        """``intCode'' -  the ``integration code'' parameter for ``analytic_integral''-method
        """
        return self.pdf.intCode ()

    # =================================================================================
    ## Safe shortcut for <code>RooAbsPdf.matchArgs(allDeps,analDeps,*vars)</code> 
    def matchArgs ( self , *vars ) :
        """Safe shortcut for RooAbsPdf.matchArgs ( allDeps , analDeps , *vars )
        """
        vv = ROOT.RooArgSet() 
        for v in vars : vv.add ( v )
        return self.pdf.matchArgs ( vv ) 
    
                                     
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
