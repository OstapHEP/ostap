#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file valerrors.py
#  Helper structures (point-with-errors) for (T)Graphs
#
#  - These objects are *NOT* for math! 
#  - The objects are used for graphs 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07 
# =============================================================================
""" Helper structures (point-with-erorrs) for (T)Graphs 
- These objects are *NOT* for math!
- The objects are used for graphs 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'AsymErrors'         ,
    'ValWithErrors'      ,
    'ValWithMultiErrors' ,
    'VAE'                , ## shortcut for ValWithErrors 
    'VME'                , ## shortcut for ValWithMultiErrors     
    ) 
# =============================================================================
from   ostap.core.ostap_types import num_types, sized_types, sequence_types   
from   ostap.core.core        import VE
from   ostap.math.base        import iszero, isequal  
import math 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils.valerrors' )
else                       : logger = getLogger( __name__                )
# =============================================================================
logger.debug ( "Helper structures (point-with-errors) for (T)Graphs")
# =============================================================================
## helper constants  for effective varibacle of the split normal distribution 
_C1 = ( 1.0 - 2.0 / math.pi ) 
## helper constants  for effective varibacle of the split normal distribution 
_C2 = math.sqrt ( 2.0 / math.pi )
# =============================================================================
#  @class AsymErrors
#  Represenattion of `asymmetric errors`
#  @code
#  ae = AsymErorrs ( -0.5 ,  1.0 ) 
#  ae = AsymErorrs ( +1.0 , -0.5 ) ## ditto
#  ae = AsymErorrs (  0.5 ,  1.0 ) ## ditto
#  @endcode 
class AsymErrors (object) :
    """`Asymmetric Errors`
    >>> ae = AsymErorrs ( -0.5 ,  1.0 ) 
    >>> ae = AsymErorrs ( +1.0 , -0.5 ) ## ditto
    >>> ae = AsymErorrs (  0.5 ,  1.0 ) ## ditto
    """
    __slots__ = ( '__positive' , '__negative' )
    ##
    def __init__ ( self , negative = 0 , positive = 0 ) :

        assert isinstance ( negative , num_types ) , "Invalid type of 'negative' error"
        assert isinstance ( positive , num_types ) , "Invalid type of 'positive' error"
        
        pos = float ( positive )
        neg = float ( negative )
        
        if iszero ( pos ) : pos = 0.0
        if iszero ( neg ) : neg = 0.0        
        
        if   neg <= 0 <= pos       : self.__negative , self.__positive =  neg , pos 
        elif pos <= 0 <= neg       : self.__negative , self.__positive =  pos , neg   
        elif 0 <= neg and 0 <= pos : self.__negative , self.__positive = -neg , pos 
        else :
            raise TypeError ( "Invalid setting of errors: %s/%s" % ( negative , positive ) )

        if iszero ( self.__negative ) : self.__negative = 0.0
        if iszero ( self.__positive ) : self.__positive = 0.0

        assert self.__negative <= 0 , 'Negative error %+g is not non-positive!' % ( self.__negative )  
        assert self.__positive >= 0 , 'Positive error %+g is not non-negative!' % ( self.__positive ) 
        
    @property
    def negative ( self ) :
        """'negative' : get the negative (low) error"""
        return self.__negative 
    @property
    def positive ( self ) :
        """'positive' : get the positive (high) error"""
        return self.__positive

    # =========================================================================
    ## An `effective error' (as split normal distribution)
    #  @see https://en.wikipedia.org/wiki/Split_normal_distribution
    @property
    def eff_error ( self ) :
        """ An `effective error' (as split normal distribution)
        - see https://en.wikipedia.org/wiki/Split_normal_distribution
        """
        s1 = abs ( self.negative )
        s2 = abs ( self.positive )
        variance = _C1 * ( ( s2 - s1 )**2 ) + s1 * s2
        return math.sqrt ( variance )
    
    # =========================================================================
    ## An `effective bias' (as split normal distribution)
    #  @see https://en.wikipedia.org/wiki/Split_normal_distribution
    @property
    def eff_bias ( self ) :
        """ An `effective bias' (as split normal distribution)
        - see https://en.wikipedia.org/wiki/Split_normal_distribution
        """
        s1 = abs ( self.negative )
        s2 = abs ( self.positive )
        return _C2 * ( s2 - s1 ) 
    
    # =========================================================================
    ## picklings/unpickling 
    def __getstate__ ( self ) :
        return self.__negative, self.__positive    
    def __setstate__ ( self , state ) :
        self.__negative = state [ 0 ]
        self.__positive = state [ 1 ]

    # =========================================================================
    ## update (via quadratic sum) with another AsymErrors object 
    def __iadd__ ( self , other ) :
        """update (via quadratic sum) with another AsymErrors object"""
        if not isinstance ( other , AsymErrors ) : return NotImplemented
        neg = other.negative
        pos = other.positive
        if   neg and self.__negative : self.__negative = -math.hypot ( self.__negative , neg )
        elif neg                     : self.__negative =  neg 
        if   pos and self.__positive : self.__positive =  math.hypot ( self.__positive , pos )
        elif pos                     : self.__positive =  pos 
        ##
        return self 

    # =========================================================================
    ## Quadratic sum of two AsymErrors object 
    def __add__ ( self , other ) :
        """Quadratic sum of two AsymErrors object"""
        if not isinstance ( other , AsymErrors ) : return NotImplemented
        result = AsymErrors ( self.__negative , slef.__positive )
        result += other
        return result

    # =========================================================================
    ## equality 
    def __eq__  ( self , other ) :
        if not isinstance ( other , AsymErrors ) : return NotImplememnted
        return  isequal ( self.negative , other.negative ) \
            and isequal ( self.positive , other.positive ) 

    # =========================================================================
    ## non-equality 
    def __ne__  ( self , other ) :
        if not isinstance ( other , AsymErrors ) : return NotImplememnted
        return not ( self == other )

    # =========================================================================
    ## tuple-like 
    def __len__ ( self ) : return 2

    # =========================================================================
    ## get the error by index
    #  0 : negative error
    #  1 : positive error 
    def __getitem__ ( self , index ) :
        """ get the error by index
        - 0 : negative error
        - 1 : positive error 
        """
        if not 0 <= index < 2 : raise IndexError( "Invalid index %s!" % index ) 
        return self.__negative if 0 == index else self.__positive

    # ==========================================================================
    ## conversion to string 
    def toString ( self , format = '( -/%.5g +/%-.5g ) ' ) :
        """Conversion to string
        """
        return format % ( abs ( self.__negative ) , self.__positive )
    
    def __str__  ( self ) : return self.toString () 
    def __repr__ ( self ) : return self.__str__ () 

# =============================================================================
## @class ValWithErrors
#  Value with asymmetric error
#  @code
#  value = ValWihErrors ( 1.0 ,  0.5 ,  0.7 )
#  value = ValWihErrors ( 1.0 , -0.5 ,  0.7 )
#  value = ValWihErrors ( 1.0 ,  0.7 , -0.5 )
#  value = ValWihErrors ( 1.0 ,  AsymErrors ( 0.5 ,0.7 ) 
#  value = ValWihErrors ( VE ( 1 , 0.5**2 ) , AsymErrors ( 0.5 ,0.7  ) 
#  @endcode 
class ValWithErrors(object) : 
    """Value with asymmetric errors
    >>> value = ValWihErrors ( 1.0 ,  0.5 ,  0.7 )
    >>> value = ValWihErrors ( 1.0 , -0.5 ,  0.7 )
    >>> value = ValWihErrors ( 1.0 ,  0.7 , -0.5 )
    >>> value = ValWihErrors ( 1.0 ,  AsymErrors ( 0.5 ,0.7 ) 
    >>> value = ValWihErrors ( VE ( 1 , 0.5**2 ) , AsymErrors ( 0.5 ,0.7  ) 
    """    
    __slots__ = (  '__value' , '__errors' )
    ##

    # =========================================================================
    ## Value with asymmetric error
    #  @code
    #  value = ValWihErrors ( 1.0 ,  0.5 ,  0.7 )
    #  value = ValWihErrors ( 1.0 , -0.5 ,  0.7 )
    #  value = ValWihErrors ( 1.0 ,  0.7 , -0.5 )
    #  value = ValWihErrors ( 1.0 ,  AsymErrors ( 0.5 ,0.7 ) 
    #  value = ValWihErrors ( VE ( 1 , 0.5**2 ) , AsymErrors ( 0.5 ,0.7  ) 
    #  @endcode 
    def __init__  ( self , value = 0 , *errors ) :
        """Value with asymmetric errors
        >>> value = ValWihErrors ( 1.0 ,  0.5 ,  0.7 )
        >>> value = ValWihErrors ( 1.0 , -0.5 ,  0.7 )
        >>> value = ValWihErrors ( 1.0 ,  0.7 , -0.5 )
        >>> value = ValWihErrors ( 1.0 ,  AsymErrors ( 0.5 ,0.7 ) 
        >>> value = ValWihErrors ( VE ( 1 , 0.5**2 ) , AsymErrors ( 0.5 ,0.7  ) 
        """            
        if   isinstance ( value , ValWithMultiErrors ) :
            self.__value    = value.value
            self.__errors   = AsymErrors ( value.neg_error , value.pos_error )
        elif isinstance ( value , ValWithErrors ) :
            self.__value    = value.value
            self.__errors   = value.errors            
        elif isinstance ( value , VE ) :
            if 0 <= value.cov2() : error = value.error ()
            else                 : error = 0            
            self.__value  = float      ( value   )            
            self.__errors = AsymErrors ( -error , error )
            
        elif isinstance ( value , sized_types )        and \
                 3 == len ( value )                    and \
                 all ( isinstance ( v  , num_types ) for v in value ) :

            self.__value  = float      (  value [ 0  ] )
            self.__errors = AsymErrors ( *value [ 1: ] )
            
        elif isinstance ( value , sized_types )        and \
                 2 == len ( value )                    and \
                 isinstance ( value[0] , num_types   ) and \
                 isinstance ( value[1] , sized_types ) and \
                 2 == len ( value[1] )                 and \
                 all ( isinstance ( v  , num_types ) for v in value[1] ) :

            self.__value  = float      (  value [ 0 ] )
            self.__errors = AsymErrors ( *value [ 1 ] ) 

        elif isinstance ( value , num_types ) :
            
            self.__value  = float      ( value   )
            self.__errors = AsymErrors ()

        else  :
            raise TypeError ( "Invalid 'value' %s/%s" % ( type ( value ) , str ( value ) ) ) 
                
        nerrs = len ( errors )
        if   not errors : pass 
        elif 1 == nerrs and isinstance ( errors [ 0 ] , AsymErrors  ) :
            self.__errors += errors [ 0 ] 
        elif 1 == nerrs and isinstance ( errors [ 0 ] , sized_types ) \
                 and 2 == len ( errors [ 0] ) \
                 and all ( isinstance ( v  , num_types ) for v in errors [ 0 ]  ) :
            self.__errors += AsymErrors ( *errors[0] )
        elif 2 == nerrs  and all ( isinstance ( v  , num_types ) for v in errors ) :
            self.__errors += AsymErrors ( *errors )
        else            : raise TypeError ( "Invalid 'errors' %s" % str ( errors ) ) 
            
    @property
    def value ( self ) :
        """'value' : the value """
        return self.__value
    @property
    def errors ( self ) :
        """'errors' : the asymmetric errors"""
        return self.__errors
    
    @property
    def neg_error ( self ) :
        """'neg_error' : the negative (low) error"""
        return self.__errors.negative
    @property
    def pos_error ( self ) :
        """'pos_error' : the positive (high) error"""
        return self.__errors.positive

    # ========================================================================
    ## conversion to float 
    def __float__ ( self ) :
        """conversion to float"""
        return self.__value

    # ========================================================================
    ## (numerical) equality
    def __eq__ ( self , other ) :
        """(Numerical) equality of object 
        """
        if   isinstance ( other , ValWithErrors ) :
            return isequal ( self.value , other.value   ) and self.errors == other.errors        
        elif isinstance ( other , VE ) :            
            if not isequal ( self.value , other.value() ) : return False            
            oerr = other.error() 
            return  isequal (       self.pos_error   , oerr ) \
                and isequal ( abs ( self.neg_error ) , oerr )        
        return NotImplemented
    
    # ========================================================================
    ## (numerical) non equality
    def __ne__ ( self , other ) :
        """(Numerical) equality of object 
        """
        if   isinstance ( other , ValWithErrors ) : return not ( self == other )         
        elif isinstance ( other , VE )            : return not ( self == other ) 
        return NotImplemented
    
    # ========================================================================
    ## Conversion to VE, as split normal distribution with optional bias 
    #  @see https://en.wikipedia.org/wiki/Split_normal_distribution
    def asVE ( self , bias = False ) :
        """Conversion to VE with optional bias 
        - see https://en.wikipedia.org/wiki/Split_normal_distribution
        """
        vv = self.value if not bias else self.value + self.__errors.eff_bias 
        return VE ( vv , self.__errors.eff_error ** 2 )

    # ========================================================================
    ## picklings/unpickling 
    def __getstate__ ( self ) :
        return self.__value, self.__errors 
    def __setstate__ ( self , state ) :
        self.__value  = state [ 0 ]
        self.__errors = state [ 1 ]


    # ==========================================================================
    ## conversion to string 
    def toString ( self , format = '( %+.5g -/%.5g +/%-.5g ) ' ) :
        """Conversion to string
        """
        return format % ( self.value , abs ( self.neg_error ) , self.pos_error )

    def __str__  ( self ) : return self.toString () 
    def __repr__ ( self ) : return self.__str__  () 

# ===========================================================================
## decode error data into the flat list of errors 
def flat_errors ( *items ) :
    """decode error data into the flat list of errors"""
    flat = [] 
    for item in items :
        even = ( 0 == len ( flat ) % 2 )
        if   isinstance ( item , num_types ) :
            ## trivial case 
            flat.append ( float ( item ) )
        elif even and isinstance ( item , AsymErrors  ) :
            flat.append ( item.negative )
            flat.append ( item.positive )
        elif even and isinstance ( item , sized_types ) and all ( isinstance ( v , num_types ) for v in item  ) :
            ee = AsymErrors ( *item ) 
            flat.append ( ee.negative ) 
            flat.append ( ee.positive ) 
        elif even and isinstance ( item , sequence_types ) :
            ## recursion
            for e in flat_errors ( *item ) : flat.append ( e ) 
        else : raise TypeError ( "Invalid 'items' %s" % str ( items ) )
    ##
    assert 0 == len ( flat ) % 2 , "Invalid 'items' %s" % str ( items ) 
    return tuple ( flat ) 
            
# =============================================================================
## @class ValWithMultiErrors
#  Value with multiple asymmetric error 
class ValWithMultiErrors(object) : 
    """Value with multiple asymmetric error"""    
    __slots__ = (  '__value' , '__errors' )
    ## 
    def __init__  ( self , value = 0 , *errors ) :

        self.__errors = [] 
        if   isinstance ( value , ValWithMultiErrors ) :
            
            self.__value   = value.value
            self.__errors += list ( value.errors )
            
        elif   isinstance ( value , ValWithErrors ) :
            
            self.__value   = value.value
            self.__errors  .append ( value.errors )
            
        elif isinstance ( value , VE ) :
            
            if 0 <= value.cov2() : error = value.error ()
            else                 : error = 0                        
            self.__value   = float ( value )
            self.__errors  .append ( AsymErrors ( -error , error ) )
            
        elif isinstance ( value , sized_types )        and \
                 3 == len ( value )                    and \
                 all ( isinstance ( v  , num_types ) for v in value ) :

            self.__value  = float      (  value [ 0  ] )
            self.__errors .append ( AsymErrors ( *value [ 1: ] ) )
            
        elif isinstance ( value , sized_types )        and \
                 2 == len ( value )                    and \
                 isinstance ( value[0] , num_types   ) and \
                 isinstance ( value[1] , sized_types ) and \
                 2 == len ( value[1] )                 and \
                 all ( isinstance ( v  , num_types ) for v in value[1] ) :

            self.__value  = float      (  value [ 0 ] )
            self.__errors .append ( AsymErrors ( *value [ 1 ] ) )
            
        elif isinstance ( value , num_types )  : self.__value  = float ( value )            
            
        else  : raise TypeError ( "Invalid 'value' %s" % str ( value ) ) 

        ## decode errors parameters 
        errs = flat_errors ( *errors )
        ## convert flat list into the list of AsymErrors 
        es   = ( AsymErrors (  errs [ i ] , errs [ i + 1 ] ) for i in range ( 0 , len ( errs ) , 2 ) )
        ## loop over additional errors        
        for ee in es  : self.__errors.append ( ee )
        
        ## ensure that the list of errors is not empty 
        if not self.__errors : self.__errors = [ AsymErrors () ]
        
        ## make it immutable 
        self.__errors   = tuple ( self.__errors ) 
        
    @property
    def value ( self ) :
        """'value' : the value """
        return self.__value
    @property
    def errors ( self ) :
        """'errors' : the list of asymmetric errors"""
        return self.__errors
    @property
    def nerrors ( self ) :
        """'nerrors' : numer of associated asymmetrical errors (pairs)"""
        return len ( self.__errors ) 
    
    ## conversion to float 
    def __float__ ( self ) :
        """conversion to float"""
        return self.__value
        
    ## picklings/unpickling 
    def __getstate__ ( self ) :
        return self.__value, self.__errors 
    def __setstate__ ( self , state ) :
        self.__value  = state [ 0 ]
        self.__errors = state [ 1 ]

    @property
    def neg_error ( self ) :
        """'neg_error' : the effective (sum squared) negative (low)error """
        return math.sqrt ( sum ( e.negative**2 for e in self.__errors ) ) 
    @property
    def pos_error ( self ) :
        """'pos_error' : the effective (sum squared) positive (high) error"""
        return math.sqrt ( sum ( e.positive**2 for e in self.__errors ) ) 


    # ========================================================================
    ## (numerical) equality
    def __eq__ ( self , other ) :
        """(Numerical) equality of object 
        """
        if   isinstance ( other , ValWithMultiErrors  ) :
            return isequal ( self.value , other.value ) and self.errors       == other.errors
        elif isinstance ( other , ValWithErrors       ) and 1 == len ( self.errors ) :
            return isequal ( self.value , other.value ) and self.errors [ 0 ] == other.errors
        elif isinstance ( other , VE                  ) and 1 == len ( self.errors ) :
            if not isequal ( self.value , other.value() ) : return False
            verr = self.errors [ 0 ]
            oerr = other.error () 
            return  isequal ( abs ( verr.negative )  , oerr ) \
                and isequal (       verr.positive    , oerr )

        return NotImplemented

    # ========================================================================
    ## (numerical) non equality
    def __ne__ ( self , other ) :
        """(Numerical) equality of object 
        """
        if   isinstance ( other , ValWithMultiErrors )                              : return not ( self == other ) 
        elif isinstance ( other , ValWithErrors      ) and 1 == len ( self.errors ) : return not ( self == other ) 
        elif isinstance ( other , VE                 ) and 1 == len ( self.errors ) : return not ( self == other ) 
        return NotImplemented

    # ========================================================================
    ## Combine multiple errors and create ValWithErrors object
    #  @code
    #  mve = ...
    #  vae = mv.asVAE () 
    #  @endcode 
    def asVAE ( self ) :
        """ Combine multiple errors and create ValWithErrors object
        >>> mve = ...
        >>> vae = mv.asVAE () 
        """
        sume = self.errors [ 0 ]
        sume = AsymErrors ( negative = sume.negative , positive = sume.positive )
        for ae in self.errors [ 1 : ] : sume += ae 
        return ValWithErrors ( self.value , sume ) 
    
    # ========================================================================
    ## Conversion to VE, as split normal distribution with optional bias 
    #  @see https://en.wikipedia.org/wiki/Split_normal_distribution
    def asVE ( self , bias = False ) :
        """Conversion to VE with optional bias 
        - see https://en.wikipedia.org/wiki/Split_normal_distribution
        """
        vae = self.asVAE ()
        return vae.asVE ( bias = bias )

    # ==========================================================================
    ## conversion to string 
    def toString ( self ,
                   format  = '( %+.5g %s) '    ,   ## global format 
                   format2 = ' -/%-.5g +/%.5g' ) : ## AsymErrors formar 
        """Conversion to string
        """

        errs = [ e.toString ( format2 ) for e in self.errors ] 
        errs = ','.join ( errs ) 
        return format % ( self.value , errs )

    def __str__  ( self ) : return self.toString () 
    def __repr__ ( self ) : return self.__str__  () 

# =============================================================================
## shortcut 
AE  = AsymErrors

# =============================================================================
## shortcut 
VAE = ValWithErrors

# =============================================================================
## shortcut 
VME = ValWithMultiErrors

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
