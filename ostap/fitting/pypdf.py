#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/pypdf.py
#  Very specific helper class PyPDF to implement "pythonic" PDF for RooFit
#  - Typical usage:
#  @code
#  import math
#  from ostap.fitting.basic import PEAK
#  from ostap.fitting.pypdf import PyPDF 
#  class PyGauss(PEAK,PyPDF) :
#       norm = 1.0 / math.sqrt ( 2 * math.pi )
#       def __init__ ( self         ,
#                      name         ,
#                      xvar         ,    // the obeservable 
#                      mean         ,    // parameter: mean-value 
#                      sigma        ,    // parameter: sigma
#                      pdf   = None ) :  // MANDATORY ARGUMENT
#           PEAK .__init__ ( self , name      , xvar , mean , sigma ) 
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
#       get all benefits from the left bas-class, <code>PEAK</code> in this  case
#  - 2. The <code>__init__</code> method *must* have keyword argument
#       <code>pdf</code> (default is <code>None</code>)
#  - 3. One *must* specify the <code>config</code> dictionary
#
#  The latter two items are mandatory for the proper implemtation of RooAbsPdf::clone
#  - @see Ostap::Models::PyPdf 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2018-06-07
# =============================================================================
"""Very specific helper class to implement 'pythonic' PDF for RooFit

    Typical usage:
    
    ... import math
    ... from ostap.fitting.basic import PEAK
    ... from ostap.fitting.pypdf import PyPDF 
    ... class PyGauss(PEAK,PyPDF) :
    ...    norm = 1.0 / math.sqrt ( 2 * math.pi )
    ...    def __init__ ( self         ,
    ...                   name         ,
    ...                   xvar         ,    // the obeservable 
    ...                   mean         ,    // parameter: mean-value 
    ...                   sigma        ,    // parameter: sigma
    ...                   pdf   = None ) :  // MANDATORY ARGUMENT
    ...        PEAK .__init__ ( self , name      , xvar , mean , sigma ) 
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
    get all benefits from the left base-class, PEAK in this  case
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
    'PyPDF'      , ## 'pythonic' PDF for RooFit 
    'PyPDFLight' , ## 'pythonic' PDF for RooFit 
    )
# =============================================================================
from   ostap.core.ostap_types import sequence_types 
from   ostap.core.core        import Ostap
from   ostap.utils.basic      import typename, prntrf  
import ostap.fitting.roofit
import ROOT, math, abc 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.pypdf' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
# @class PyPDF
# Very simple "Pythonic" PDF
# One must define:
#  - <code>clone</code> method 
#  - proper constructor that takes keyword <code>clone</code>
#  - <code>evaluate</code> method 
#  - method <code>__reduce__</code> for serialization (if needed)
# 
# Twomore methods arte needed if one needs analystical integrals :
#  - `get_analytical_integral`
#  - `analytical_integral`
class PyPDF(Ostap.Models.PyPdf) :
    """ Very simple `pythonic' PDF 
    One must define: 
    - `clone` method 
    - proper constructor that takes keyword <code>clone</code>
    - `evaluate` method 
    - method `reduce` for serialization (if needed) 
    """
    def __init__ ( self , name , title = '' , variables = () , clone = None ) :
        """ Constructor that accepts `clone` argument 
        """

        assert not clone or isinstance ( clone , PyPDF ), \
            "PyPDF: invalid `clone` type:%s " % typename ( clone ) 
        
        assert clone or ( variables                                and \
                          isinstance ( variables, sequence_types ) and \
                          all ( isinstance ( v , ROOT.RooAbsReal ) for v in variables ) ) , \
                          "PyPDF: invalid `variables`: %s/%s" %  ( typename ( variables ) , str ( variables ) ) 

                
        if clone :
            super ( PyPDF , self ) .__init__ ( clone , name if name else clone.name )            
        else     :
            vv = ROOT.RooArgList ( v for v in variables ) 
            super ( PyPDF, self ) .__init__ ( name , title if title else 'PyPDf(%s)' % name , vv )

        self._keep = variables, 
        if clone : self._keep += clone._keep
        
    ## redefine 
    def matchArgs ( self , *vars ) :
        return self.match_args ( ROOT.RooArgSet ( v for v in vars ) ) 
    ## redefine 
    def matchArg  ( self ,  var  ) : return self.match_arg  ( var  ) 
    
    @abc.abstractmethod 
    def clone ( self  , newname = '' ) :
        raise NotImplementedError("PyPDF.clone must be implemented")

    @abc.abstractmethod
    def evaluate ( self ) :
        raise NotImplementedError("PyPDF.evaluate must be implemented")

    @property
    def variables ( self ) :
        """`variables` : get list of variables (same as `varlist()` 
        - see `Ostap.Models.PyPdf.varlist` 
        """
        return self.varlist()

# =============================================================================
## Very simple `ready-to-use' pythonic PDF
#  @see Ostap.Modeld.PyPdfLite 
def PyPDFLite ( name            ,
                function        ,
                variables       ,
                title     = ''  ) : 
    """ Very simple `ready-to-use' pythonic PDF
    - see `Ostap.Modeld.PyPdfLite` 
    """

    assert callable ( function ) , \
        "PyPDFLite: invalid `function`: %s/%s" %  ( typename ( function ) , prntrf ( function ) )
    
    assert variables and isinstance ( variables, sequence_types )   and \
        all ( isinstance ( v , ROOT.RooAbsReal ) for v in variables ) , \
        "PyPDFLite: invalid `variables`: %s/%s" %  ( typename ( variables ) , str ( variables ) ) 
    
    title = title if title else 'Python PDF %s with %s/%s' % ( name                  ,
                                                               typename ( function ) ,
                                                            prntrf   ( function ) )
    vv = ROOT.RooArgList ( v for v in variables )
    return Ostap.Models.PyPdfLite ( name , title , function , vv )


def __pdfl_str__    ( self ) :
    return '%s(%s,%s,%s/%s,#%d)' % ( typename ( self ) ,
                                     self.name  ,
                                     self.title ,
                                     typename ( self.function () ) ,
                                     prntrf   ( self.function () ) ,
                                     self.numrefs() )

Ostap.Models.PyPdfLite.__str__  = __pdfl_str__
Ostap.Models.PyPdfLite.__repr__ = __pdfl_str__

def ppdfl_factory ( config ) : return PyPDFLite ( **config )
def _ppdfl_reduce_ ( self ) :
    config = { 'name'      : self.name       ,
               'function'  : self.function() ,
               'variables' : self.varlist()  , 
               'title'     : self.title      } 
    return  ppdfl_factory , ( config , )


Ostap.Models.PyPdfLite.__reduce__ = _ppdfl_reduce_

old_PyROOT = False 
# =============================================================================
if old_PyROOT :

    __all__ = (
        'PyPDF'  , ## 'pythonic' PDF for RooFit 
        'PyPDF2' , ## 'pythonic' PDF for RooFit 
        )
    # =========================================================================
    ## @class PyPDF
    #  helper base class to implement 'pure-python' PDF
    #  - Typical usage:
    #  @code
    #  import math 
    #  from ostap.fitting.basic import PEAK
    #  from ostap.fitting.pypdf import PyPDF 
    #  class PyGauss(PEAK,PyPDF) :
    #       norm = 1.0 / math.sqrt ( 2 * math.pi )
    #       def __init__ ( self         ,
    #                      name         ,
    #                      xvar         ,    // the obeservable 
    #                      mean         ,    // parameter: mean-value 
    #                      sigma        ,    // parameter: sigma
    #                      pdf   = None ) :  // MANDATORY ARGUMENT
    #           PEAK .__init__ ( self , name      , xvar , mean , sigma ) 
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
    #       get all benefits from the left base-class, <code>PEAK</code> in this  case
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
        """Helper base class to implement 'pure-python' PDF
        
        Typical usage (for very OLD PyROOT):
        
        ... import math
        ... from ostap.fitting.basic import PEAK
        ... from ostap.fitting.pypdf import PyPDF 
        ... class PyGauss(PEAK,PyPDF) :
        ...    norm = 1.0 / math.sqrt ( 2 * math.pi )
        ...    def __init__ ( self         ,
        ...                   name         ,
        ...                   xvar         ,    // the obeservable 
        ...                   mean         ,    // parameter: mean-value 
        ...                   sigma        ,    // parameter: sigma
        ...                   pdf   = None ) :  // MANDATORY ARGUMENT
        ...        PEAK .__init__ ( self , name      , xvar , mean , sigma ) 
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
        get all benefits from the left base-class, PEAK in this  case
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
        _storage = []
        
        def __init__ ( self          ,
                       name          ,
                       vars          , ## all variables (observables&parameters) 
                       title  = ''   ,
                       pypdf  = None ) :
            
            assert ( not pypdf ) or isinstance ( pypdf , Ostap.Models.PyPdf  ) , \
                   "Invalid type of 'pypdf': %s/%s" % ( pypdf , type ( pypdf ) )
            
            if not pypdf :
                
                ## convert to RooArgList if needed 
                if not isinstance ( vars , ROOT.RooArgList ) :
                    vv = ROOT.RooArgList ()
                    for v in vars :
                        if not v in vv :  vv.add ( v )
                    vars = vv    
                        
                self.__pyvars = vars
                
                if not title : title = "PyPDF(%s)"  % name 
                
                pypdf = Ostap.Models.PyPdf ( self , name , title , self.__pyvars  )
                
            ## take care on ROOT ownership pf the object
            ROOT.SetOwnership ( pypdf , False )
            
            ## define the pdf 
            self.__pypdf = pypdf
            
            self.config =  {
                'name'    : self.pypdf.GetName   () ,
                'title'   : self.pypdf.GetTitle  () ,
                'vars'    : self.pypdf.variables () ,  
                'pypdf'   : None 
                }

        @property
        def pypdf ( self ) :
            """'pypdf' : get the actual Ostap::Models::PyPdf object"""
            return self.__pypdf
        @pypdf.setter
        def pypdf ( self , value ) :
            assert value is None  or isinstance (  value , ROOT.RooAbsPdf ),\
                   'Invalid type of pypdf:%s' % type ( value )        
            self.__pypdf == value 
            
        @property
        def varlist ( self ) :
            """'varlist' : get all variables (as RooArgList)  from Ostap::Models::PyPdf
            - all variables *MUST* be accessed via this collection!
            """
            return self.pypdf.varlist ()
        
        @property
        def allDeps ( self ) :
            """'allDeps' - the first parameter  of 'getAnalyticalIntegra/get_analytical_integral'-method 
            - attention:  it could be a null pointer! check it before usage! 
            """
            return self.pypdf.allDeps()
        
        @property
        def analDeps ( self ) :
            """'analDeps' - the second parameter  of 'getAnalyticalIntegral/get_analystical_integral'-method 
            - attention:  it could be a null pointer! check it before usage!
            """
            return self.pypdf.allDeps()
        
        @property
        def rangeName ( self )  :
            """'rangeName' -  the 'rangeName'  parameter for `(get_)analytical_integral'-method
            - attention - it could be a null pointer! check it before usage! 
            """
            return self.pypdf.rangeName ()
        
        @property
        def intCode   ( self )  :
            """'intCode' -  the 'integration code' parameter for 'analytic_integral'-method
            """
            return self.pypdf.intCode ()
        
        # =====================================================================
        ## The method  that MUST be implemented
        #  @code
        #  pdf = ...
        #  pdf.evaluate() 
        #  @endcode
        def evaluate ( self ) :
            """The method  that MUST be implemented
            >>> pdf = ...
            >>> pdf.evaluate() 
            """
            raise NotImplementedError("PyPDF: 'evaluate' is not implemented!")
        
        # =====================================================================
        ## get a value of certain variale
        #  @code
        #  pdf = ...
        #  a = pdf.variable ( 'a' )
        #  b = pdf.variable (  2  )
        #  @endcode 
        def variable ( self , tag ) :
            """Get a value of certain variale
            >>> pdf = ...
            >>> a = pdf.variable ( 'a' )
            >>> b = pdf.variable (  2  )
            """
            return self.pypdf.variable ( tag ) 
        # =========================================================================
        ## Safe shortcut for <code>RooAbsPdf.matchArgs(allDeps,analDeps,*vars)</code> 
        def matchArgs ( self , *vars ) :
            """Safe shortcut for RooAbsPdf.matchArgs ( allDeps , analDeps , *vars )
            """
            vv = ROOT.RooArgSet() 
            for v in vars : vv.add ( v )
            return self.pypdf.matchArgs ( vv ) 
        
        ## clone the object (needed by C++ partner class)
        def clone ( self , **kwargs ) :
            """Clone the object
            """
            
            conf = {}
            conf.update ( self.config )
            
            ## modify the name if the name is in config  
            if 'name' in conf : conf['name'] += '_copy'
            
            conf.update ( kwargs )
            
            KLASS = self.__class__
            return KLASS ( **conf )

        
# =============================================================================
## @class PyPDF2
#  'Light' version of 'pythonic-Pdf'
#  
#  @see Ostap::Models::PyPdf
#  @see Ostap::Models::PyPdf2
class PyPDF2(object) :
    """  'Light' version of 'pythonic-Pdf'
    """
    def __init__ (  self       ,
                    name       ,
                    function   , ##  python callable (var1,var2,...,varn)-> double  
                    vars       , 
                    title = '' ) :

        ## function must be valid function! 
        assert function and callable ( function ) , "'function' is not callable!"

        if not title : title = 'PyPDF2(%s)' % name

        self.__pyfunction = function

        self.__argvars    = vars 
        if not isinstance ( vars , ROOT.RooArgList )  :
            self.__pyvars = vars 
            vv = ROOT.RooArgList()
            for v in vars : vv.add ( v ) 
            vars = vv
            
        self.__vars  = vars
        
        ## create the actual RooAbsPdf 
        pypdf        = Ostap.Models.PyPdf2 ( name          ,
                                             title         ,
                                             self.function ,
                                             self.__vars   )
        ROOT.SetOwnership ( pypdf , False )

        self.__pypdf = pypdf


        ## finally define PDF 
        self.pdf     = self.pypdf

        ## store configuration
        self.config = {
            'name'     : self.pypdf.GetName  () ,
            'title'    : self.pypdf.GetTitle () ,
            'function' : self.function          ,
            'vars'     : self.variables                                           
            }

    @property
    def pypdf ( self ) :
        """'pypdf'   : get the actual Ostap::Models::PyPdf2 object"""
        return  self.__pypdf 

    @property
    def function ( self ) :
        """'function' : get the actual python function/callable"""
        return self.__pyfunction
    
    @property
    def variables  ( self ) :
        """'variables' : list(ROOT.RooArgList) of all variables"""
        return self.__vars
    
    @property
    def vars       ( self ) :
        """'vars' : tuple of all variables"""
        return  tuple ( [ i for i in self.__vars ] ) 
    
    ## clone the object (needed by C++ partner class)
    def clone ( self , **kwargs ) :
        """Clone the object
        """
        
        conf = {}
        conf.update ( self.config )
        
        ## modify the name if the name is in config  
        if 'name' in conf : conf['name'] += '_copy'

        conf.update ( kwargs )
        
        KLASS = self.__class__
        return KLASS ( **conf )
        
        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END 
# =============================================================================
