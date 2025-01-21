#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostap/fitting/constraint.py
#  Function to use/add/apply constrains as a product of RooAbsPdfs 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2024-01-20
# =============================================================================
""" Function to use/add/apply constrains as a product of RooAbsPdfs 
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    'make_constrained' , ## helper function to create constraied PDFs
)
# =============================================================================
from   ostap.core.meta_info         import root_info 
from   ostap.core.core              import typename 
import ostap.fitting.variables
import ostap.fitting.roocollections
import ROOT 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.constraint' )
else                       : logger = getLogger ( __name__                   )
# =============================================================================
## make constrained PDF 
def make_constrained ( pdf , *constraints ) :
    """ Make constrained PDF 
    """
    assert isinstance ( pdf , ROOT.RooAbsPdf ) , "Invalid pdf type: %s" % typename ( pdf )
    assert all ( isinstance ( c , ROOT.RooAbsPdf ) for c in constraints ) , \
        "Invalid type of constraint(s):[%s]" % ( ','.join ( typename ( c ) for c in constraints ) )
    
    name  = 'Constrained_%s[%s]' % ( pdf.name  , '&'.join( c.name for c in constraints ) )
    title = 'Constrained_%s[%s]' % ( pdf.title , '&'.join( c.name for c in constraints ) )

    if isinstance ( pdf , ROOT.RooSimultaneous ) :
        category = pdf.indexCat()
        result   = ROOT.RooSimultaneous ( name , title , category )
        keep_it  = [] 
        for label, index in category.items () :
            cmp           = pdf.getPdf ( label )
            newcmp , tail = make_constrained ( cmp , *constraints )
            keep_it.append ( newcmp ) 
            if tail : keep_it.append ( tail ) 
            result.addPdf ( newcmp , str ( label ) )
        return result, tuple ( keep_it )   

    if root_info < ( 6, 24 ) :
        plist = ROOT.RooArgSet( pdf )
        for p in constraints  : plist.add ( p ) 
        return ROOT.RooProdPdf ( name , title , plist ) , () 
        
    return ROOT.RooProdPdf ( name , title , ( pdf , ) + constraints ) , () 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END 
# =============================================================================
