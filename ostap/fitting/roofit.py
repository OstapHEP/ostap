#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/roofit.py
#  Module with decoration of some RooFit objects for efficient use in python
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Decoration of some RooFit objects for efficient use in python
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'setStorage'    , ## define the default storage for  RooDataStore 
    'useStorage'    , ## define (as context) the default storage for  RooDataStore 
    'PDF_fun'       , ## wrapper of PDF to "simple" function 
    'SETVAR'        , ## context manager to preserve the current value for RooRealVar
    'FIXVAR'        , ## context manager to fix/unfix the variable 
    'var_from_name' , ## "convert" name/expression into variable/formula
    ) 
# =============================================================================
from   ostap.core.core              import Ostap, VE, valid_pointer
from   ostap.core.ostap_types       import integer_types, string_types 
from   ostap.fitting.variables      import SETVAR, FIXVAR  
from   ostap.fitting.dataset        import setStorage, useStorage
from   ostap.utils.basic            import typename 
import ostap.fitting.roocollections
import ostap.fitting.rooreduce 
import ostap.fitting.roofitresult
import ostap.fitting.printable
import ostap.fitting.roocmdarg   
import ROOT, random, math, copy  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.rootfit' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug( 'Some useful decorations for RooFit objects')
# =============================================================================
_new_methods_ = []

# =============================================================================
## product of two PDFs 
def _pdf_mul_ ( pdf1 , pdf2 ) :
    """ Easy contruct for the product of two PDFs:
    
    >>> pdf1 = ...
    >>> pdf2 = ...
    
    >>> product = pdf1 * pdf2 
    """
    return Ostap.Models.Product (
        '%s*%s'             % ( pdf1.GetName  () , pdf2.GetName  () ) ,
        'Product: %s & %s ' % ( pdf1.GetTitle () , pdf2.GetTitle () ) ,
        pdf1 , pdf2 )
# ============================================================================
ROOT.RooAbsPdf . __mul__  = _pdf_mul_ 

_new_methods_ += [
    ROOT.RooAbsPdf.__mul__  , 
    ]

# ============================================================================
## "convert" name/expression into variable/formula
def var_from_name ( w , varset ) :
    """ Convert name/expression into variable/formula
    """
    w = w.strip() 
    if    0 <= w.find('(') < what.find(')') : pass
    elif  0 <  w.find('*')                  : pass
    elif  0 <  w.find('/')                  : pass
    elif  0 <  w.find('%')                  : pass 
    elif  0 <  w.find('+')                  : pass
    elif  0 <  w.find('-')                  : pass
    else :
        v = varset[w]
        return v
    ##
    
    vlst = ROOT.RooArgList()
    for s in varset : vlst.add ( s )
    #
    f = Ostap.FormulaVar( w , w , vlst )
    return f 

# =============================================================================
## @class PDF_fun
#  Helper class to wrap PDF as 'function'
#  can be helpful for some pure math-operations
#  @code
#  pdf,var = ....
#  fun     = PDF( fun , var , xmin=0 , xmax=1 )
#  from ostap.stats.moments import mean, mode, median, CL   
#  print 'MEAN    : %s' % mean    ( fun , 0 , 1 )
#  print 'MODE    : %s' % mode    ( fun , 0 , 1 )
#  print 'MEDIAN  : %s' % median  ( fun , 0 , 1 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-03-29
class PDF_fun(object):
    """ Helper class to wrap PDF as 'function'
    >>> pdf,var = ....
    >>> fun     = PDF( pdf , var , xmin=0 , xmax=1 )
    >>> print fun(0.1),fun(0.5) 
    >>> from ostap.stats.moments import mean, mode, median  
    >>> print 'MEAN    : %s' % mean    ( fun , 0 , 1 )
    >>> print 'MODE    : %s' % mode    ( fun , 0 , 1 )
    >>> print 'MEDIAN  : %s' % median  ( fun , 0 , 1 )
    """
    ##
    def __init__ ( self , pdf , xvar , xmin = None , xmax = None , norm = 1 , normset = None ) :
        
        self.__pdf        = pdf
        self.__norm       = float ( norm )

        self.__normset    = None 
        if normset : 
            self.__normset    = ROOT.RooArgSet ()
            for v in normset : self.__normset.add ( v )
            
        ## ostap stuff: 
        if not isinstance ( pdf , ROOT.RooAbsPdf ) :
            if hasattr ( self.pdf , 'pdf' ) :
                self.__pdf_ =   pdf 
                self.__pdf  = pdf.pdf

        self.__xvar   = xvar
        self.__xmin   = None 
        self.__xmax   = None
        
        if not xmin is None : self.__xmin = xmin 
        if not xmax is None : self.__xmax = xmax

        if hasattr ( xvar , 'getMin' ) :
            if self.__xmin is None : self.__xmin = xvar.getMin()
            else                   : self.__xmin = max ( self.__xmin , xvar.getMin() )
            
        if hasattr ( xvar , 'getMax' ) :
            if self.__xmax is None : self.__xmax = xvar.getMax()
            else                   : self.__xmax = min ( self.__xmax , xvar.getMax() )
            
        if self.__xmin is None :
            raise AttributeError ( "xmin can't be deduced from  input arguments" )

        if self.__xmax is None :
            raise AttributeError ( "xmax can't be deduced from  input arguments" )
        
        if self.__xmin > self.__xmax :
            self.__xmin , self.__xmax = self.__xmax , self.__xmin

    @property
    def pdf ( self ) :
        """'pdf' : get the actual RooAbsPdf
        """
        return self.__pdf
    @property
    def norm ( self ) :
        """'norm' : additional normalization factor
        """
        return self.__norm

    @property
    def xvar ( self ) :
        """'xvar': x-variable
        """
        return self.__xvar
    
    def xmin     ( self ) : return self.__xmin
    def xmax     ( self ) : return self.__xmax
    
    ## the main method 
    def __call__ ( self , x , pars = [] ) :

        ## for ROOT.TF1
        if   hasattr ( x , '__len__' ) and hasattr( x , 'getitem' ) and not len( x )   : x = x [ 0 ]
        elif hasattr ( x , '__len__' ) and hasattr( x , 'getitem' ) and 0 != x.size()  : x = x [ 0 ]
        
        ## try to be efficient 
        if not self.__xmin <= x <= self.__xmax : return 0 
        
        with SETVAR ( self.__xvar ) :
            self.__xvar.setVal ( x )
            if self.__normset : 
                return self.__pdf.getVal ( self.__normset ) * self.__norm
            else : 
                return self.__pdf.getVal (                ) * self.__norm
            
# ============================================================================
ROOT.RooPlot.__len__ =  lambda s : math.floor ( s.numItems() ) 

# =============================================================================
def _rp_contains_ ( plot , what ) :
    return isinstance ( what , integer_types ) and 0 <= what < len ( plot )

ROOT.RooPlot.__contains__ =  _rp_contains_

# =============================================================================
## Get <code>RooPlot</code> componet 
#  @code
#  frame = ...
#  o = frame[2] 
#  @endcode
def _rp_getitem_  ( plot , index ) :
    """ Get `ROOT.RooPlot` componet 
    >>> frame = ...
    >>> o = frame[2] 
    """
    
    if not isinstance ( index, integer_types ) or not index in plot : 
        raise IndexError('Index %s in not in ROOT.RooPlot' % index )

    return plot.getObject ( index )

ROOT.RooPlot.__getitem__  =  _rp_getitem_

# =============================================================================
## Iterator over <code>RooPlot</code> componnet 
#  @code
#  frame = ...
#  for o in frame : ...
#  @endcode
def _rp_iter_  ( plot ) :
    """ Iterator over `RooPlot` componnets 
    >>> frame = ...
    >>> for o in frame : ...
    """

    n = len ( plot ) 
    for i in range ( n ) :
        yield plot.getObject ( i ) 

ROOT.RooPlot.__iter__  =  _rp_iter_

# =============================================================================
## Iterator over <code>RooPlot</code> items 
#  @code
#  frame = ...
#  for (obj,options,inviisble) in frame.items()  : ...
#  for (obj,options,invisible) in frame.iteritems()  : ... ## ditto 
#  @endcode
def _rp_items_  ( plot ) :
    """ Iterator over <code>RooPlot</code> items  
    >>> frame = ...
    >>> for (obj,options,inviisble) in frame.items()  : ...
    >>> for (obj,options,invisible) in frame.iteritems()  : ... ## ditto 
    """
    n = len ( plot )
    names = set() 
    for i in range ( n ) :
        obj       = plot.getObject      ( i )
        name      = plot.nameOf         ( i )  
        options   = plot.getDrawOptions ( name )
        # Trick! 
        invisible = plot.getInvisible   ( name )
        if invisible and name in names : invisible = False
        elif invisible                 : names.add ( name ) 
        #
        yield obj , str ( options ) , invisible

ROOT.RooPlot.items      =  _rp_items_
ROOT.RooPlot.iteritems  =  _rp_items_
  
# =============================================================================
## format <code>RooPlot</code> as a table
#  @code
#  frame = ...
#  print ( frame.table( title = 'Title', prefix = '# ' ) 
#  @endcode 
#  @see RooPlot 
def _rp_table_ ( plot , prefix = '' , title = '' ) :
    """Format <code>RooPlot</code> as a table
    >>> frame = ...
    >>> print ( frame.table( title = 'Title', prefix = '# ' ) 
    - see `RooPlot`
    """
    
    def _name ( obj  ) :
        n = typename ( obj )
        p = n.find ( 'cppyy.gbl.' ) 
        return n [ p + 10 : ] if 0 < p else n
    
    if not title  :
        title = 'RooPlot %s' % plot.name
        from ostap.logger.symbols import palette as plot_symbol
        if plot_symbol : title = '%s %s' % ( plot_symbol , title )
        
    table = [ ( 'Index' , 'Type' , 'Option' , 'Binw' , 'Line' , 'Marker' , 'Fill' , 'Name' ) ]

    names = set() 
    for index , obj in enumerate ( plot )  :
        
        name      = plot.nameOf ( index )
        options   = str ( plot.getDrawOptions ( name ) )
        invisible = plot.getInvisible ( name )
        ## Trick! 
        if invisible :
            if not name in names : options = options + ":I"
            names.add  ( name )
            
        row  = [ '%2d' % index , _name ( obj ) , options ]

        if not isinstance ( obj , ROOT.RooHist ) : row.append ( '' )
        else : row.append ( '%g' % obj.getNominalBinWidth() ) 
            

        ## line attributes 
        if isinstance ( obj , ROOT.TAttLine ) : 
            row.append ( '%2d/%d/%d' % ( obj.GetLineStyle () ,
                                         obj.GetLineColor () ,
                                         obj.GetLineWidth () ) ) 
        else : row.append ('')

        ## marker attributes 
        if isinstance ( obj , ROOT.TAttMarker ) : 
            row.append ( '%2d/%d/%s' % ( obj.GetMarkerStyle () ,
                                         obj.GetMarkerColor () ,
                                         obj.GetMarkerSize  () ) )
        else : row.append ('')

        ## Fill attributes 
        if isinstance ( obj , ROOT.TAttFill ) : 
            row.append ( '%4d/%-3d' % ( obj.GetFillStyle () ,
                                        obj.GetFillColor () ) )
        else : row.append ('')

        
        ## 
        row.append  ( name )
        
        row = tuple ( row  )
        
        table.append ( row )

    import ostap.logger.table as T
    return T.table ( table , title = title, prefix = prefix , alignment = 'clcclllw' )

ROOT.RooPlot.table    =  _rp_table_
ROOT.RooPlot.__str__  =  _rp_table_
ROOT.RooPlot.__repr__ =  _rp_table_


# ================================================================================
## Friendly print of RooSimultaneous
def _rsim_table_ ( pdf , title = '' , prefix = '' , style = '' ) :
    """ Friendly print of RooSimultaneous
    """
    assert isinstance ( pdf , ROOT.RooSimultaneous ) , \
        'PDF is not RooSimultaneopus!'
    
    rows = [ ( 'Category' , 'Index' , 'Type' , 'Name' , 'Title' ) ] 
    cat    = pdf.indexCat() 
    for label, index in cat.items() :
        cpdf = pdf.getPdf ( label )
        if not valid_pointer ( cpdf ) :
            print ('INVALID POINTER' , label )
        else :
            row  = '%s' % label , '%+d' % index , typename ( cpdf ) , cpdf.name , cpdf.title
            rows.append ( row )

    import ostap.logger.table as T
    title = title if title else '%s(%s,%s)' % ( typename ( pdf ) , pdf.name , pdf.title )
    return T.table ( rows, title = title , prefix = prefix , alignment = 'lllll' )
     
ROOT.RooSimultaneous.table    = _rsim_table_
ROOT.RooSimultaneous.__str__  = _rsim_table_
ROOT.RooSimultaneous.__repr__ = _rsim_table_

_new_methods_ += [
    ROOT.RooSimultaneous.table    , 
    ROOT.RooSimultaneous.__str__  ,
    ROOT.RooSimultaneous.__repr__ , 
    ]

# ================================================================================
## copy RooPlot object
def _rp_copy_ ( plot ) :
    """ Copy RooPlot object"""
    return copy.deepcopy ( plot )

ROOT.RooPlot.copy      =  _rp_copy_ 

_new_methods_ += [
    ROOT.RooPlot.__len__      ,
    ROOT.RooPlot.__contains__ ,
    ROOT.RooPlot.__getitem__  ,
    ROOT.RooPlot.__iter__     ,
    ROOT.RooPlot.__iter__     ,
    ROOT.RooPlot.__str__      ,
    ROOT.RooPlot.__repr__     ,
    ROOT.RooPlot.table        ,    
    ROOT.RooPlot.items        , 
    ROOT.RooPlot.iteritems    , 
    ]

# =============================================================================
_decorated_classes_ = (
    ROOT.RooPlot         , 
    ROOT.RooSimultaneous , 
    )

_new_methods_ = tuple ( _new_methods_ ) 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                     The END 
# =============================================================================
