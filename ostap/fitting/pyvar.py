#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/pyvar.py
#  Very specific helper class PyVAR to implement "pythonic" RooAbsReal for RooFit
#  - Typical usage:
#  @code
#  import math
#  from ostap.fitting.pyvar import PyVAR
#  class MyVar(PyVAR) :
#     
#     def evaluate ( self ) :
# 
#       vlist = self.varlist
#       x     = float ( vlist[0] )
#       x0    = float ( vlist[1] )
#       sigma = float ( vlist[2] )
#
#       dx =  (x-x0)/sigma
#       return  math.exp ( 0.5* dx ** 2 ) / math.sqrt ( 2*math.pi*sigma )
#
#   ## varibales 
#   x     = RooRealVar( ... )
#   m0    = RooRealVar( ... )
#   sigma = RooRealVar( ... )
#
#   ## create the formula 
#   myVar = MyVar ( name = 'my_var' , vars = ( x,m0,sigma) )
#
#   for v in ( 0.1 , 0.2  , 0.3 , 0.4 , ... ) :
#      x.setVal ( v )
#      print 'function value is %s ', myVar()
#   @endcode
#  - @see Ostap::Functions::PyVar
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2018-06-07
# =============================================================================
"""Very specific helper class to implement ``pythonic'' RooAbsReal for RooFit

- Typical usage:

    ... import math
    ... from ostap.fitting.pyvar import PyVAR
    ... class MyVar(PyVAR) :
    ...
    ...     def evaluate ( self ) :
    ...
    ...        vlist = self.varlist
    ...        x     = float ( vlist[0] )
    ...        x0    = float ( vlist[1] )
    ...        sigma = float ( vlist[2] )
    ...
    ...        dx =  (x-x0)/sigma
    ...
    ...        return  math.exp ( 0.5* dx ** 2 ) / math.sqrt ( 2*math.pi*sigma )
    
    ... ## variables 
    ... x     = RooRealVar( ... )
    ... m0    = RooRealVar( ... )
    ... sigma = RooRealVar( ... )
    
    ... ## create the formula 
    ... myVar = MyVar ( name = 'my_var' , vars = ( x,m0,sigma) )
    
    ... for v in ( 0.1 , 0.2  , 0.3 , 0.4 , ... ) :
    ...     x.setVal ( v )
    ...     print 'function value is %s ', myVar()
    
    - see Ostap::Functions::PyVar 
    """
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'PyVAR' , ## ``pythonic'' RooAbsReal for RooFit 
    )
# =============================================================================
import ROOT, math
from   ostap.core.core      import Ostap
import ostap.fitting.roofit 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.pyvar' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
## @class PyVAR
#  Very specific helper class PyVAR to implement "pythonic" RooAbsReal for RooFit
#  - Typical usage:
#  @code
#  import math
#  from ostap.fitting.pyvar import PyVAR
#  class MyVar(PyVAR) :
#     
#     def evaluate ( self ) :
# 
#       vlist = self.varlist
#       x     = float ( vlist[0] )
#       x0    = float ( vlist[1] )
#       sigma = float ( vlist[2] )
#
#       dx =  (x-x0)/sigma
#       return  math.exp ( 0.5* dx ** 2 ) / math.sqrt ( 2*math.pi*sigma )
#
#   ## varibales 
#   x     = RooRealVar( ... )
#   m0    = RooRealVar( ... )
#   sigma = RooRealVar( ... )
#
#   ## create the formula 
#   myVar = MyVar ( name = 'my_var' , vars = ( x,m0,sigma) )
#
#   for v in ( 0.1 , 0.2  , 0.3 , 0.4 , ... ) :
#      x.setVal ( v )
#      print 'function value is %s ', myVar()
# @endcode
class PyVAR (object) :
    """Helper base class to implement ``pure-python'' RooAbsVar
    
    - Typical usage:

    ... import math
    ... from ostap.fitting.pyvar import PyVAR
    ... class MyVar(PyVAR) :
    ...
    ...     def evaluate ( self ) :
    ...
    ...        vlist = self.varlist
    ...        x     = float ( vlist[0] )
    ...        x0    = float ( vlist[1] )
    ...        sigma = float ( vlist[2] )
    ...
    ...        dx =  (x-x0)/sigma
    ...
    ...        return  math.exp ( 0.5* dx ** 2 ) / math.sqrt ( 2*math.pi*sigma )
    
    ... ## variables 
    ... x     = RooRealVar( ... )
    ... m0    = RooRealVar( ... )
    ... sigma = RooRealVar( ... )
    
    ... ## create the formula 
    ... myVar = MyVar ( name = 'my_var' , vars = ( x,m0,sigma) )
    
    ... for v in ( 0.1 , 0.2  , 0.3 , 0.4 , ... ) :
    ...     x.setVal ( v )
    ...     print 'function value is %s ', myVar()
    
    - see Ostap::Functions::PyVar 

    """
    def __init__ ( self          ,
                   name          ,   ## the name 
                   vars   = []   ,   ## all variables 
                   title  = ''   ,   ## the title 
                   pyvar  = None ) : ## helper 

        assert vars or ( pyvar and isinstance ( pyvar , Ostap.Functions.PyVar ) ),\
               "Invaild configuration of PyVAR: %s/%s" % ( vars , pypdf ) 
        
        if not pyvar :
            
            ## convert to RooArgList if needed 
            self.__pyvars  = vars        
            if not title : title = "PyVAR(%s)"  % name
        
            if not isinstance ( vars , ROOT.RooArgList ) :
                vv = ROOT.RooArgList ()
                for v in vars :
                    if not v in vv :  vv.add ( v )
                vars = vv 
            self.__pyvars = vars 
            pyvar = Ostap.Functions.PyVar ( self , name , title , self.__pyvars  )
            ROOT.SetOwnership ( pyvar , False )
            
        ## finally set the variable!
        self.__var   = pyvar

    @property
    def var     ( self ) :
        """``var'': the actual C++ RooAbsReal object"""
        return self.__var
    
    @property
    def name    ( self ) :
        """``name'' : the name of the varibale"""
        return self.var.GetName()
    
    @property
    def title   ( self ) :
        """``title'' : the title of the variable"""
        return self.var.GetTitle()
    
    @property
    def varlist ( self ) :
        """``variables'' : get all variables (as RooArgList)  from Ostap::Functions::PyVar
        """
        return self.var.varlist ()
    
    # =========================================================================
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
        return self.pdf.variable ( tag ) 

    ## clone the object (needed by C++ partner class)
    def clone ( self , **kwargs ) :
        """Clone the object
        -  attention: existing PDF is ``copied'', unless specified via kwargs (by C++)
        """
        conf =  {
            'name'    : self.name    ,
            'title'   : self.title   , 
            'vars'    : self.varlist ,
            'pyvar'   : None
            }
        
        conf.update ( kwargs )

        KLASS = self.__class__
        return KLASS ( **conf )

    # ==========================================================================
    ## the method  that MUST be implemented
    #  @code
    #  var = ...
    #  var.evaluate() 
    #  @endcode
    def evaluate ( self ) :
        """The method  that MUST be implemented
        >>> var = ...
        >>> var.evaluate() 
        """
        raise NotImplementedError("PyVAR: ``evaluate'' is not implemented!")
    

    # =========================================================================
    ## get the value of the function
    #  @code
    #  var = ...
    #  print var() 
    #  @endcode
    def __call__  ( self , local = True ) :
        """Get the value of the function
        >>> var = ...
        >>> var() 
        """
        return self.evaluate() if local else self.getVal () 

    # =========================================================================
    ## get the value of the RooAbsReal
    #  @code
    #  var = ...
    #  var.getVal()
    #  @endcode
    def getVal    ( self ) :
        """Get the value of RootAbsReal
        >>> var = ...
        >>> var.getVal()
        """
        return self.var.getVal() 
       

# =============================================================================
## @class PyVAR2
#  Light version of PyVAR
class PyVAR2 (object) :
    """Helper base class to implement ``pure-python'' RooAbsVar
    
    - Typical usage:

    ... import math
    ... from ostap.fitting.pyvar import PyVAR
    ... class MyVar(PyVAR) :
    ...
    ...     def evaluate ( self ) :
    ...
    ...        vlist = self.varlist
    ...        x     = float ( vlist[0] )
    ...        x0    = float ( vlist[1] )
    ...        sigma = float ( vlist[2] )
    ...
    ...        dx =  (x-x0)/sigma
    ...
    ...        return  math.exp ( 0.5* dx ** 2 ) / math.sqrt ( 2*math.pi*sigma )
    
    ... ## variables 
    ... x     = RooRealVar( ... )
    ... m0    = RooRealVar( ... )
    ... sigma = RooRealVar( ... )
    
    ... ## create the formula 
    ... myVar = MyVar ( name = 'my_var' , vars = ( x,m0,sigma) )
    
    ... for v in ( 0.1 , 0.2  , 0.3 , 0.4 , ... ) :
    ...     x.setVal ( v )
    ...     print 'function value is %s ', myVar()
    
    - see Ostap::Functions::PyVar 

    """
    def __init__ ( self          ,
                   name          ,   ## the name
                   function      ,   ## the function 
                   vars          ,   ## all variables 
                   title  = ''   ):  ## the title 
        
        ## function must be valid function! 
        assert function and callable ( function ) , "``function'' is not callable!"

        if not title : title = "PyVAR2(%s)"  % name
        
        self.__pyfunction = function
        self.__argvars    = vars

        self.__argvars    = vars 
        if not isinstance ( vars , ROOT.RooArgList )  :
            self.__pyvars = vars 
            vv = ROOT.RooArgList()
            for v in vars : vv.add ( v ) 
            vars = vv
            
        self.__vars  = vars
        
        ## create the actual RooAbsPdf 
        pyvar        = Ostap.Models.PyVar2 ( name          ,
                                             title         ,
                                             self.function ,
                                             self.__vars   )
        ROOT.SetOwnership ( pyvar , False )
    
        ## finally set the variable!
        self.__var   = pyvar

    @property
    def var       ( self ) :
        """``var'': the actual C++ RooAbsReal object"""
        return self.__var
    
    @property
    def name      ( self ) :
        """``name'' : the name of the varibale"""
        return self.var.GetName()
    
    @property
    def title     ( self ) :
        """``title'' : the title of the variable"""
        return self.var.GetTitle()
    
    @property
    def variables ( self ) :
        """``variables'' : list(ROOT.RooArgList) of all variables"""
        return self.__vars
    
    @property
    def vars      ( self ) :
        """``vars'' : tuple of all variables"""
        return  tuple ( [ i for i in self.__vars ] ) 

    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END 
# =============================================================================
