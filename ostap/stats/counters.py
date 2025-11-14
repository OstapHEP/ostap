#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  @see Ostap::StatEntity
#  @see Ostap::WStatEntity
#  @see Ostap::NStatEntity
#  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
#  @date 2009-09-12
# =============================================================================
"""Simple counters 
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@nikhef.nl"
__date__    = "2009-09-12"
__version__ = ""
# =============================================================================
__all__     = (
    'SE'             , ## simple smart counter (1-bi histo)  C++ Ostap::StatEntity 
    'WSE'            , ## simple smart counter with weight :     Ostap::WStatEntity 
    'NSE'            , ## simple smart running counter     :     Ostap::NStatEntity
    'table_counters' , ## make table of counters
    'EffCounter'     , ## simple counter for efficiency
    ) 
# =============================================================================
from   ostap.math.ve          import Ostap, VE
from   ostap.math.base        import std, isequal, isequalf, iszero
from   ostap.utils.basic      import typename 
from   ostap.core.ostap_types import ( dictlike_types , sequence_types ,
                                       integer_types  , sized_types    )
import ROOT, cppyy, math, sys  
# =============================================================================
_new_methods_ = []
# =============================================================================
SE    = Ostap.StatEntity 
WSE   = Ostap.WStatEntity 
NSE   = Ostap.NStatEntity 

# =============================================================================
# minor decoration for StatEntity 
# ============================================================================= 
if not hasattr ( SE , '_orig_sum'  ) : 
    _orig_sum    = SE.sum
    SE._orig_sum = _orig_sum
    
if not hasattr ( SE , '_orig_mean' ) : 
    _orig_mean    = SE.mean
    SE._orig_mean = _orig_mean

SE. sum     = lambda s : VE ( s._orig_sum  () , s.sum2 ()      )
SE. minmax  = lambda s :    ( s.min        () , s.max  ()      ) 
SE. mean    = lambda s : VE ( s._orig_mean () , s.meanErr()**2 )

# ==============================================================================
## Update the counter
#  @code
#  >>> cnt  = ...
#  >>> cnt += value
#  @endcode
def _se_iadd_ ( self , other ) :
    """ Update the counter
    >>> cnt  = ...
    >>> cnt += value 
    """
    return self.add ( other )

# ==============================================================================
## Update the counter
#  @code
#  >>> cnt  = ...
#  >>> cnt -= value
#  @endcode
def _se_isub_ ( self , other ) :
    """ Update the counter
    >>> cnt  = ...
    >>> cnt -= value 
    """
    return self.add ( -other )

# ==============================================================================
## add two counters
#  @code
#  >>> cnt1 = ...
#  >>> cnt2 = ...
#  >>> cnt  = cnt1 +  cnt2 
#  @endcode
def _se_add_ ( self , other ) :
    """ Update the counter
    >>> cnt  = ...
    >>> cnt += value 
    """
    r  = SE ( self )
    r.add ( other )
    return r

# ==============================================================================
## add two counters
#  @code
#  >>> cnt1 = ...
#  >>> cnt2 = ...
#  >>> cnt  = cnt1 +  cnt2 
#  @endcode
def _wse_add_ ( self , other ) :
    """ Update the counter
    >>> cnt  = ...
    >>> cnt += value 
    """
    r = WSE ( self )
    r.add ( other )
    return r
# =============================================================================

SE  .__iadd__ =  _se_iadd_
SE  .__add__  =  _se_add_
WSE .__add__  = _wse_add_

# ==============================================================================
## numbers are close enough?
def closef ( a , b , N = 1000 ) :
    return ( a == b ) or isequalf ( a , b ) or abs ( a - b ) * N <= abs ( a ) + abs ( b )
    
# =============================================================================
## equal counters?
def _se_eq_ ( s1 , s2 ) :
    """ Numerically equal counters?"""
    ##
    if not isinstance ( s2 , SE ) : return NotImplemented
    ## 
    ## closef   ( s1.mu2 () , s2.mu2 () ) and \
    return s1.n     () == s2.n   ()           and \
           isequalf ( s1.mu  () , s2.mu  () ) and \
           isequalf ( s1.mu2 () , s2.mu2 () ) and \
           s1.min () == s2.min ()             and \
           s1.max () == s2.max ()
# =============================================================================
## non-equal counters?
def _se_ne_ ( s1 , s2 ) :
    """ Numerically non-equal counters?
    """
    if not isinstance ( s2 , SE ) : return NotImplemented 
    return not ( s1 == s2 )

SE.__eq__ = _se_eq_
SE.__ne__ = _se_ne_

_new_methods_ += [
    SE.sum       ,
    SE.mean      ,
    SE.minmax    ,
    SE.__add__   ,
    SE.__iadd__  ,    
    SE.__isub__  ,    
    SE.__eq__    , 
    SE.__ne__    ,
    ]

# =============================================================================
## equal counters?
def _wse_eq_ ( s1 , s2 ) :
    """ Numerically equal counters?"""
    ##
    if not isinstance ( s2 , WSE ) : return NotImplemented
    ##
    ## closef     ( s1.mu2 () , s2.mu2 () ) and \
    return s1.n       () == s2.n       ()       and \
           isequalf   ( s1.mu  () , s2.mu  () ) and \
           isequalf   ( s1.mu2 () , s2.mu2 () ) and \
           s1.weights () == s2.weights ()       and \
           s1.values  () == s2.values  ()

# =============================================================================
## non-equal counters?
def _wse_ne_ ( s1 , s2 ) :
    """ Numerically non-equal counters?"""
    ##
    if not isinstance ( s2 , WSE ) : return NotImplemented
    ##    
    return not ( s1 == s2 )

WSE.__eq__ = _wse_eq_
WSE.__ne__ = _wse_ne_

# =============================================================================
# minor decoration for WStatEntity 
# ============================================================================= 
if not hasattr ( WSE , '_orig_sum'  ) : 
    _orig_sum     = WSE.sum
    WSE._orig_sum = _orig_sum

if not hasattr ( WSE , '_orig_mean' ) : 
    _orig_mean_wse = WSE.mean
    WSE._orig_mean = _orig_mean_wse
    
WSE. sum     = lambda s : VE ( s._orig_sum  () , s.sum2()       )
WSE. mean    = lambda s : VE ( s._orig_mean () , s.meanErr()**2 )
WSE. minmax  = lambda s : ( s.min () , s.max () ) 

_new_methods_ += [
    WSE.sum         ,
    WSE.mean        ,
    WSE.minmax      ,
    WSE.__add__     ,
    WSE.__eq__      ,
    WSE.__ne__      ,
]

# =============================================================================
## Count iterable
#  @code
#  c = SE.count ( [1,2,3,4] )  
#  @endcode
def se_count ( values ) :
    """ Count iterable
    >>> c = SE.count ( [1,2,3,4] )  
    """
    cnt = SE ()
    for v in values : cnt.add ( v ) 
    return cnt
SE.count = staticmethod ( se_count )

_new_methods_ += [
    SE.count ,
    ]

# ============================================================================== 
## Build the single row from the counter 
def cnt_row ( counter , key = '' ) :
    """ Build the single row from the SE-counter 
    """
    mean        = counter.mean   ()
    rms         = counter.rms    ()
    minv, maxv  = counter.minmax () 
    vsum        = counter.sum    ()    
    nEff        = counter.nEff   ()

    from ostap.logger.pretty  import pretty_float, fmt_pretty_values

    if isinstance ( nEff , integer_types ) : nEff , expo0 = '%d' % nEff , 0 
    else                                   : nEff , expo0 = pretty_float ( nEff )

    vsum , expo1  = vsum.pretty_print ( parentheses = False ) 
    mean , expo2  = mean.pretty_print ( parentheses = False ) 
    rms  , expo3  = pretty_float      ( rms )
    fmtx , expo4  = fmt_pretty_values ( minv , maxv , with_sign = True )
    
    mnmx          = '%s/%s' % ( fmtx , fmtx )
    if expo4 : scale = 1.0/(10**expo4)
    else     : scale = 1 
    mnmx = mnmx % ( minv * scale , maxv * scale  )

    from   ostap.logger.symbols   import times
    
    return '%s' % key , \
        nEff  , '%s10^{%+d}' % ( times , expo0 ) if expo0 else ''  , \
        vsum  , '%s10^{%+d}' % ( times , expo1 ) if expo1 else ''  , \
        mean  , '%s10^{%+d}' % ( times , expo2 ) if expo2 else ''  , \
        rms   , '%s10^{%+d}' % ( times , expo3 ) if expo3 else ''  , \
        mnmx  , '%s10^{%+d}' % ( times , expo4 ) if expo4 else '' 

# =============================================================================
## Build 1 or 3 rows in the table for the counter :
def cnt_rows ( cnt , key = '' ) :
    """ Build 1 or 3 rows in the table for the counter """
    
    if   isinstance ( cnt , WSE ) :
        return \
            cnt_row ( cnt             , ( '%s'         % key ) if key else ''       ) , \
            cnt_row ( cnt.values   () , ( '%s/values'  % key ) if key else 'values' ) , \
            cnt_row ( cnt.weights  () , ( '%s/weights' % key ) if key else 'weight' ) ,
    
    elif isinstance ( cnt , NSE ) :
        return \
            cnt_row ( cnt             , ( '%s'         % key ) if key else ''       ) , \
            cnt_row ( cnt.counter1 () , ( '%s/first'   % key ) if key else 'first'  ) , \
            cnt_row ( cnt.counter2 () , ( '%s/second'  % key ) if key else 'second' ) ,
    
    return cnt_row ( cnt  , '%s' % key  ) , 
    
VSE  = std.vector(SE)
VWSE = std.vector(WSE)
counter_types = SE, WSE, NSE 
# =============================================================================
## Make table of counters
#  @code
#  counters = .... ## sequence or mapping for counters
#  table    = table_counters ( counters , prrefix = '# ' )
#  logger.info ( 'Table is \n%s' % table )
#  @endcode
def table_counters ( counters , prefix = '' , title = '' , style = None ) :
    """ Make table of counters
    >>> counters = .... ## sequence or mapping for counters
    >>> table     = counters_table ( counters , prrefix = '# ' )
    >>> logger.info ( 'Table is \n%s' % table )
    """
                              
    if   isinstance ( counters , dictlike_types ) and \
         all ( isinstance ( v  , counter_types  ) for v in counters.values() ) : pass
    
    elif isinstance ( counters , sequence_types ) and \
         isinstance ( counters , sized_types )    and \
         all ( isinstance ( v  , counter_types  ) for v in counters  ) :
        
        cnts = {}
        n = len ( counters )
        from ostap.logger.symbols import labels 
        for l , c in zip ( labels ( n ) , counters ) : cnts [ l ] = c
        counters = cnts 
    
    elif isinstance ( counters , sequence_types ) and \
         all ( isinstance ( v  , counter_types  ) for v in counters  ) :
        
        cnts = {}
        for c , i  in enumerate ( counters , start = 1 ) : cnts [ i ] = c
        counters = cnts 
        
    elif isinstance ( counters , counter_types ) :
        
        counters = { 1 : counters }
        
    ## vectors of counters 
    elif  isinstance ( counters , ( VSE , VWSE ) ) :        
        cnts = tuple ( v for v in counters )
        return table_counters ( cnts , prefix = prefix , title = title , style = style )
    
    else :
        
        raise TypeError ( "counter_table: Invalid type for 'counters' %s" % typename ( counters ) )

    from   ostap.logger.symbols  import sum_symbol, rms_symbol  
    rows = [ ( ''          ,    \
               '#'        , '' , 
               sum_symbol , '' , 
               'mean'     , '' ,
               rms_symbol , '' , 
               'min/max'  , '' ) ]
    
    for key in counters :
        cnt = counters [ key ]        
        for row in cnt_rows ( cnt , key ) : rows.append ( row )
        
    import ostap.logger.table as T
    rows = T.remove_empty_columns ( rows ) 
    if not title : title = 'Table of %d counters' 
    table = T.table ( rows , prefix = prefix , title = title , alignment = "lcccccccccccc" , style = style )
    #
    return table 

# =============================================================================
SE   .__repr__ = lambda s : table_counters ( { '' : s } , title = 'Counter' )
SE   .__str__  = lambda s : table_counters ( { '' : s } , title = 'Counter' )  
NSE  .__repr__ = lambda s : table_counters ( { '' : s } , title = 'Counter' )
NSE  .__str__  = lambda s : table_counters ( { '' : s } , title = 'Counter' )  

_new_methods_ += [
    SE .__repr__   ,
    SE .__str__    ,
    NSE.__repr__   ,
    NSE.__str__    ,
    table_counters ]

        
# =============================================================================
WSE .__repr__ = lambda s : table_counters ( { '' : s } , title = 'Counter' )
WSE .__str__  = lambda s : table_counters ( { '' : s } , title = 'Counter' )  

VSE .__str__  = lambda s : table_counters ( s , title = 'Counters' )
VSE .__repr__ = lambda s : table_counters ( s , title = 'Counters' )
VWSE.__str__  = lambda s : table_counters ( s , title = 'Counters' )
VWSE.__repr__ = lambda s : table_counters ( s , title = 'Counters' )

_new_methods_ += [
    WSE.__repr__    ,
    WSE.__str__     ,
    ##
    VWSE.__repr__   ,
    VWSE.__str__    ,
    VSE .__repr__   ,
    VSE .__str__    ]

# ==============================================================================
## @class EffCounter
#  A primitive "effciciency" counter
#  - It keeps number of `accepted` and `rejected` entries
#  and allows to calcuale the (binomial) efficiency
#  @code
#  import random
#  cnt = EffCounter
#  for i in range ( 1000 ) :
#      x = random.uniform ( 0 , 1 ) 
#      cnt += ( x < 0.5 ) 
#  eff = cnt.eff 
#  @endcode 
class EffCounter(object):
    """ A primitive `effciency' counter
    - It keeps number of `accepted` and `rejected` entries
    and allows to calcuale the (binomial) efficiency

    >>> import random
    >>> cnt = EffCounter
    >>> for i in range ( 1000 ) :
    ...     x = random.uniform ( 0 , 1 ) 
    ...     cnt += ( x < 0.5 )
    >>> eff = cnt.eff 
    """
    __slots__ = '__A' , '__R'
    # =========================================================================
    def __init__ ( self ) :
        ## accepted 
        self.__A = 0 ## acccepted
        ## rejected          
        self.__R = 0 ## rejected 
    # =========================================================================
    @property
    def accepted ( self ) :
        """`accepted` : number of `accepted' cases (same as `success`)"""
        return self.__A
    # =========================================================================
    @property
    def rejected ( self ) :
        """`rejected` : number of `rejected' cases"""
        return self.__R
    @property
    def success  ( self ) :
        """`success` : number of `success' cases (same as `accepted')"""
        return self.accepted 
    # =========================================================================
    @property
    def failure  ( self ) :
        """`failure` : number of `failed' cases (same as `rejected')"""
        return self.rejected 
    # =========================================================================
    @property
    def total    ( self ) :
        """`total`    : total number of entries """
        return self.__A + self.__R
    # ==========================================================================
    ## Increment the counter 
    def add ( self , other ) :
        """ Increment the counter """
        self += other    
        return self 
    # =========================================================================
    ## Increment the counter 
    def __iadd__  ( self , other ) :
        """ Increment the counter """
        if isinstance ( other , EffCounter ) :
            self.__A += other.accepted
            self.__R += other.rejected
            return self
        if other : self.__A += 1
        else     : self.__R += 1 
        return self
    # ========================================================================
    ## Add two conuters 
    def __add__ ( self , other ) :
        """ Add two counters """
        if isinstance ( other , EffCounter ) :
            newc  = EffCounter()
            newc += self
            newc += other 
            return newc
        return NotImplemented
    # ==========================================================================
    ## Get (binomial) efficiency
    #  @see Ostap::Math::binomEff 
    @property 
    def efficiency ( self ) :
        """`efficiency` : get (binomial) efficiency
        - see `Ostap.Math.binomEff`         
        """
        return Ostap.Math.binomEff ( self.accepted , self.total )
    # ==========================================================================
    @property 
    def eff        ( self ) :
        """`eff` : get (binomial) efficiency"""
        return self.efficiency
    # ==========================================================================
    def __getstate__ ( self ) : return self.__A , self.__R 
    def __setstate__ ( self , state ) :
        self.__A , self.__R = state
    # ==========================================================================
    
# ==============================================================================
## REDUCE 
# ==============================================================================
from ostap.math.reduce import root_factory

# =============================================================================
## Reduce Ostap::StatEntity 
def  _se_reduce_ ( cnt ) :
    """ Reduce `Ostap.StatEntity` object
    """
    return  root_factory , ( type ( cnt ) ,
                             cnt.n   () ,
                             cnt.mu  () ,
                             cnt.mu2 () ,
                             cnt.min () ,
                             cnt.max () )

# =============================================================================
## Reduce Ostap::WStatEntity 
def  _wse_reduce_ ( cnt ) :
    """ "Reduce `Ostap.WStatEntity` object
    """
    return  root_factory , ( type ( cnt )   ,
                             cnt.mu      () ,
                             cnt.mu2     () ,
                             cnt.values  () ,
                             cnt.weights () )
# =============================================================================
## Reduce Ostap::NStatEntity 
def  _nse_reduce_ ( cnt ) :
    """ Reduce `Ostap.NStatEntity` object
    """
    return  root_factory , ( type ( cnt ) ,
                             cnt.N    ()  ,
                             cnt.cnt1 ()  , 
                             cnt.cnt2 () )

SE .__reduce__ =  _se_reduce_
WSE.__reduce__ = _wse_reduce_
NSE.__reduce__ = _nse_reduce_

_new_methods_ += [
    SE .__reduce__ , 
    WSE.__reduce__ , 
    NSE.__reduce__ , 
]

# ==============================================================================
_new_methods_ = tuple ( _new_methods_ )

# =============================================================================
_decorated_classes_ = (
    SE , WSE , NSE 
    )

# =============================================================================
if '__main__' == __name__  :

    from ostap.logger.logger import getLogger
    logger = getLogger( 'ostap.stats.counters' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    import random
    
    cnt = SE() 
    for i in range(10000) :
        cnt += random.gauss(1,1)
        
    logger.info ( 'Counter: %s' % cnt        ) 
    logger.info ( 'Mean   : %s' % cnt.mean() ) 
    logger.info ( 'RMS    : %s' % cnt.rms () ) 
    
    logger.info (80*'*')
    
# =============================================================================
##                                                                      The END 
# =============================================================================
