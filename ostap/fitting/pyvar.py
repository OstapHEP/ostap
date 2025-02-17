#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/pyvar.py
#
#  For *OLD* PyROOT:
#  Very specific helper class PyVAR to implement "pythonic" RooAbsReal for RooFit
#
#  @see Ostap::Functions::PyVar
#  @see Ostap::Functions::PyVar2
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2018-06-07
# =============================================================================
"""Very specific helper class to implement 'pythonic' RooAbsReal for RooFit

- For *OLD* PyROOT only

- see Ostap.Functions.PyVar
- see Ostap.Functions.PyVarLite
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    'PyVAR'     , ## simple utility to build "pythonic" RooAbsReal 
    'PyVARLite' , ## simple utility to build "pythonic" RooAbsReal 
    ) 
# =============================================================================
from   ostap.core.ostap_types       import sequence_types 
from   ostap.core.core              import Ostap
from   ostap.utils.basic            import typename, prntrf  
import ostap.fitting.roocollections
import ostap.fitting.variables     
import ROOT, abc, math
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.pyvar' )
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
class PyVAR(Ostap.Functions.PyVar) :
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
        assert not clone or isinstance ( clone , PyVAR), \
            "PyVAR: invalid `clone` type:%s " % typename ( clone ) 
        
        assert clone or ( variables                                and \
                          isinstance ( variables, sequence_types ) and \
                          all ( isinstance ( v , ROOT.RooAbsReal ) for v in variables ) ) , \
                          "PyPDF: invalid `variables`: %s/%s" %  ( typename ( variables ) , str ( variables ) )      
        if clone :
            super ( PyVAR , self ) .__init__ ( clone , name if name else clone.name )            
        else     : 
            vv = ROOT.RooArgList () ;
            for v for v in variables : vv.add ( v ) 
            super ( PyVAR, self ) .__init__ ( name , title if title else 'PyPDf(%s)' % name , vv )

        self._keep = variables, 
        if clone : self._keep += clone._keep
        
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
#  @see Ostap.Functions.PyVarLite 
def PyVARLite ( name            ,
                function        ,
                variables       ,
                title     = ''  ) : 
    """ Very simple `ready-to-use' pythonic variable
    - see `Ostap.Functions.PyVarLite` 
    """

    assert callable ( function ) , \
        "PyVARLite: invalid `function`: %s/%s" %  ( typename ( function ) , prntrf ( function ) )
    
    assert variables and isinstance ( variables, sequence_types )   and \
        all ( isinstance ( v , ROOT.RooAbsReal ) for v in variables ) , \
        "PyVARLite: invalid `variables`: %s/%s" %  ( typename ( variables ) , str ( variables ) ) 
    
    title = title if title else 'Python PDF %s with %s/%s' % ( name                  ,
                                                               typename ( function ) ,
                                                               prntrf   ( function ) )
    vv = ROOT.RooArgList ( v for v in variables )
    return Ostap.Functions.PyVarLite ( name , title , function , vv )

# =============================================================================
## printout of PyVarLite object 
def _pvarl_str_    ( self ) :
    """ printout of PyVarLite object """
    return '%s(%s,%s,%s/%s,#%d)' % ( typename ( self ) ,
                                     self.name  ,
                                     self.title ,
                                     typename ( self.function () ) ,
                                     prntrf   ( self.function () ) ,
                                     self.numrefs() )

Ostap.Functions.PyVarLite.__str__   = _pvarl_str_
Ostap.Functions.PyVarLite.__repr__  = _pvarl_str_

## The factory to de-serialize the PyVarLine onjects 
def pvarl_factory ( config ) : return PyVARLite ( **config )
## Reduce/serialize  PyPdfLite objects 
def _pvarl_reduce_ ( self ) :
    """ Reduce/serialize  PyPdfLite objects """
    config = { 'name'      : self.name       ,
               'function'  : self.function() ,
               'variables' : self.varlist()  , 
               'title'     : self.title      } 
    return  pvarl_factory , ( config , )

Ostap.Functions.PyVarLite.__reduce__  = _pvarl_reduce_
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
