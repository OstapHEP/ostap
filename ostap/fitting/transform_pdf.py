#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostap/fitting/transform_pdf.py
#  PDF for transformed variables 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2020-05-11
# =============================================================================
"""PDF for transformed variables
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2020-05-11"
__all__     = (
    ##
    'TrPDF'         , ## 1D-transformed PDF
    ##
    )
# =============================================================================
import ostap.fitting.roofit
from   ostap.fitting.variables import var_mul  
from   ostap.fitting.funbasic  import AFUN1 
from   ostap.fitting.pdfbasic  import PDF1, make_pdf 
import ROOT
# =============================================================================
from   ostap.logger.logger     import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.transform_pdf' )
else                       : logger = getLogger ( __name__                      )
# =============================================================================
## @class TrPDF
#  PDF of the transformed variable.

#  - e.g. gaussian in log10 scale:
#  @code
#  
#  x  = ROOT.RooRealVar ( 'x'  ,''        ,1,100000 ) ## original/old variable
#  lx = ROOT.RooRealVar ( 'lx' ,'log10(x)',0,5      ) ## new_variable
#  LX = Fun1D( lx , lx )
#  NX = 10 ** LX   ## old variable as a function of new variable 
#
#  ## PDF as function of old variable:
#  g1 = Gauss_pdf ( 'G1' , xvar = x , mean = 10000 , sigma = 10000 )
#
#  ## PDF as function of new variable 
#  g2 = TrPDF ( pdf = g1, new_var = NX )
#  @endcode
# 
#  Optionally the absolute value of the jacobian can be specified:
#  @code
#  J = math.log(10) * ( 10**LX)
#  ## PDF as function of new variable 
#  g2 = TrPDF ( pdf = g1, new_var = NX , jacob = J )
#  @endcode
class TrPDF(PDF1) :
    """ PDF of the transformed variable.
    
    - e.g. gaussian in log10 scale:
    >>> x  = ROOT.RooRealVar ( 'x'  ,''        ,1,100000 ) ## original/old variable
    >>> lx = ROOT.RooRealVar ( 'lx' ,'log10(x)',0,5      ) ## new_variable
    >>> LX = Fun1D( lx , lx )
    >>> NX = 10 ** LX                      ## old variable as function of new variable
    - PDF as function of old variable:
    >>> g1 = Gauss_pdf ( 'G1' , xvar = x , mean = 10000 , sigma = 10000 )
    - PDF as function of new variable 
    >>> g2 = TrPDF ( pdf = g1, new_var = NX )
    
    Optionally the absolute value of the jacobian can be specified:
    >>> J = math.log(10) * ( 10**LX) 
    >>> g2 = TrPDF ( pdf = g1, new_var = NX , jacob = J )
    """
    
    def __init__ ( self         ,
                   pdf          ,   ## template PDF of "old" variable  
                   new_var      ,   ## old variable as function of a new variable
                   jacob = None ,   ## absolute value of the Jacobian: |d(old)/d(new)|
                   name  = ''   ) : ## proposed name   
        
        assert pdf     and isinstance ( pdf     ,  PDF1 ) , 'Invalid PDF type %s'     % type ( pdf     )
        assert new_var and isinstance ( new_var , AFUN1 ) , 'Invalid new_var type %s' % type ( new_var )
        
        xvar = new_var.xvar

        name = name if name else "Transform_%s" % pdf.name 
        PDF1.__init__ ( self , name , xvar = xvar )
        
        if not jacob : jacob = abs ( new_var.dFdX () )

        assert isinstance ( jacob , AFUN1 ) , 'Invalid Jacobian %s' % type ( jacob )
        if not xvar in jacob : self.warning ( 'Jacobian does not depend on xvar!')  

        self.__jacob   = jacob
        self.__new_var = new_var 
        self.__ori_pdf = pdf
        
        self.__clo_pdf = pdf.clone ( xvar = new_var.fun )

        ## new PDF as a function
        self.__new_fun = var_mul ( jacob.fun  , self.__clo_pdf.pdf )

        self.__new_pdf = make_pdf ( self.__new_fun , [ xvar ] , self.name + '_' )
        
        ## finally the new PDF:
        self.pdf  = self.__new_pdf.pdf
        
        self.config = {
            'name'    : self.name     ,
            'pdf'     : self.orig_pdf ,
            'new_var' : self.new_var  ,
            'jacob'   : self.jacob    ,  
            }
        
        self.checked_keys.add  ( 'pdf'      )
        self.checked_keys.add  ( 'new_var'  )
        self.checked_keys.add  ( 'jacob'    )
        
    @property
    def orig_pdf ( self ) :
        """'orig_pdf': original, not-transformed PDF"""
        return self.__ori_pdf

    @property
    def new_pdf ( self ) :
        """'new_pdf':  transformed PDF"""
        return self.__new_pdf
    
    @property
    def new_var ( self ) :
        """'new_var'  : new/old variable """
        return self.__new_var

    @property
    def jacob( self ) :
        """'jacob'  : absolute value of the Jacobian |d(old)/d(new)| """
        return self.__jacob 

        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
