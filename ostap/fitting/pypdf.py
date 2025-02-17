#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/pypdf.py
#  Very specific helper classes  to implement "pythonic" PDF for RooFit
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2018-06-07
# =============================================================================
""" Very specific helper class to implement 'pythonic' PDF for RooFit
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'PyPDF'     , ## 'pythonic' PDF for RooFit 
    'PyPDFLite' , ## 'pythonic' PDF for RooFit 
    )
# =============================================================================
from   ostap.core.ostap_types       import sequence_types 
from   ostap.core.core              import Ostap
from   ostap.utils.basic            import typename, prntrf  
import ostap.fitting.roocollections
import ostap.fitting.variables 
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

# =============================================================================
## printout of PyPdfLite object 
def _ppdfl_str_    ( self ) :
    """ printout of PyPdfLite object """
    return '%s(%s,%s,%s/%s,#%d)' % ( typename ( self ) ,
                                     self.name  ,
                                     self.title ,
                                     typename ( self.function () ) ,
                                     prntrf   ( self.function () ) ,
                                     self.numrefs() )

Ostap.Models.PyPdfLite.__str__  = _ppdfl_str_
Ostap.Models.PyPdfLite.__repr__ = _ppdfl_str_

## The factory to de-serialize the PyPdfLine onjects 
def ppdfl_factory ( config ) : return PyPDFLite ( **config )
## Reduce/serialize  PyPdfLite objects 
def _ppdfl_reduce_ ( self ) :
    """ Reduce/serialize  PyPdfLite objects """
    config = { 'name'      : self.name       ,
               'function'  : self.function() ,
               'variables' : self.varlist()  , 
               'title'     : self.title      } 
    return  ppdfl_factory , ( config , )

Ostap.Models.PyPdfLite.__reduce__ = _ppdfl_reduce_

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END 
# =============================================================================
