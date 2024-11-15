#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file valerrors.py
#  Helper structures (point-with-errors) for (T)Graphs
#
#  - These objects are *NOT* for math! 
#  - The objects are used for graphs&plots only!  
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07 
# =============================================================================
""" Helper structures (point-with-erorrs) for (T)Graphs 
- These objects are *NOT* for math!
- The objects are used for graphs&plots only! 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'AsymErrors'         ,
    'ValWithErrors'      ,
    'ValWithMultiErrors' ,    
    'AE'                 , ## shortcut for AsymErrors 
    'VAE'                , ## shortcut for ValWithErrors 
    'VME'                , ## shortcut for ValWithMultiErrors     
) 
# =============================================================================
from   ostap.core.ostap_types import num_types, sized_types, sequence_types   
from   ostap.core.core        import VE
from   ostap.math.base        import iszero, isequal
from   ostap.logger.utils     import fmt_pretty_errs, fmt_pretty_err2  
import math, copy  
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
    ## make a copy 
    def __copy__ ( self ) :
        """Make a copy"""
        return AsymErrors ( negative = self.__negative , positive = self.__positive )
    
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
        variance = _C1 * ( ( s2 - s1 ) ** 2 ) + s1 * s2
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
        """ Update (via quadratic sum) with another AsymErrors object"""
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
    ## Quadratic sum of two AsymErrors objects 
    def __add__ ( self , other ) :
        """ Quadratic sum of two AsymErrors objects
        """
        if not isinstance ( other , AsymErrors ) : return NotImplemented
        r  = copy.copy ( self ) 
        r += other
        return r

    # =========================================================================
    ## self-scaling
    def __imul__ ( self , other ) :
        """ Self-scaling 
        """
        if not isinstance ( other , num_types ) : return NotImplemented
        ## 
        a = self.__negative * other
        b = self.__positive * other
        self.__negative = min ( a , b )
        self.__positive = max ( a , b )        
        return self 
    # =========================================================================
    ## self-scaling
    def __idiv__ ( self , other ) :
        """ Self-scaling 
        """
        if not isinstance ( other , num_types ) : return NotImplemented
        ##
        a = self.__negative / other
        b = self.__positive / other
        self.__negative = min ( a , b )
        self.__positive = max ( a , b )        
        return self 
    # ==========================================================================
    ## Self-scaling
    __itruediv__ = __idiv__    
    # ==========================================================================
    ## Multiplication/scaling
    def __mul__ ( self , other ) :
        """ Multiplication/scaling
        """
        if not isinstance ( other , num_types ) : return NotImplemented
        r  = copy.copy ( self ) 
        r *= other
        return r  
    # ==========================================================================
    ## division/scaling 
    def __div__ ( self , other ) :
        """ Division/scaling 
        """
        if not isinstance ( other , num_types ) : return NotImplemented
        r  = copy.copy ( self ) 
        r /= other
        return r 
    # ==========================================================================
    ## right multuiplication 
    __rmul__    = __mul__
    __truediv__ = __div__ 
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
        """ Get the error by index
        - 0 : negative error
        - 1 : positive error 
        """
        if not 0 <= index < 2 : raise IndexError( "Invalid index %s!" % index ) 
        return self.__negative if 0 == index else self.__positive

    # ==========================================================================
    ## conversion to string 
    def toString ( self , format = '( -/%.5g +/%-.5g ) ' ) :
        """C onversion to string
        """
        return format % ( abs ( self.__negative ) , self.__positive )
    
    # =========================================================================
    ## pretty print
    #  @code
    #  errors = ...
    #  result, expo = errors.pretty_print () 
    #  @endcode
    def pretty_print  ( self               ,
                        width       = 6    ,
                        precision   = 4    ,
                        parentheses = True ) : 
        """ Pretty print
        >>> errors = ...
        >>> result, expo = errors.pretty_print () 
        """
        return pretty_ae ( self                      ,
                           width       = width       ,
                           precision   = precision   ,
                           parentheses = parentheses )
    
    def __str__  ( self ) : return self.toString () 
    def __repr__ ( self ) : return self.__str__ () 

# =============================================================================
the_types = num_types + ( VE , ) 
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
    """ Value with asymmetric errors
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
        """ Value with asymmetric errors
        >>> value = ValWihErrors ( 1.0 ,  0.5 ,  0.7 )
        >>> value = ValWihErrors ( 1.0 , -0.5 ,  0.7 )
        >>> value = ValWihErrors ( 1.0 ,  0.7 , -0.5 )
        >>> value = ValWihErrors ( 1.0 ,  AsymErrors ( 0.5 ,0.7 ) 
        >>> value = ValWihErrors ( VE ( 1 , 0.5**2 ) , AsymErrors ( 0.5 ,0.7  ) 
        """            
        if   isinstance ( value , ValWithMultiErrors ) :
            self.__value    = float ( value.value ) 
            self.__errors   = AsymErrors ( value.neg_error , value.pos_error )
        elif isinstance ( value , ValWithErrors ) :
            self.__value    = float ( value.value ) 
            self.__errors   = AsymErrors ( value.neg_error , value.pos_error )
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

    # =========================================================================
    ## make a copy 
    def __copy__ ( self ) :
        """ Make a copy """ 
        return ValWithErrors ( self ) 
        
    # =========================================================================
    ## self-scaling
    #  @code
    #  vae = ...
    #  vae *= 5 
    #  @endcode 
    def __imul__ ( self , other ) :
        """ Self-scaling
        >>> vae = ...
        >>> vae *= 5 
        """
        if isinstance ( other , VE ) and other.cov2() <= 0 : other = float ( other ) 
        ## 
        if isinstance ( other , num_types ) :
            self.__value  *= other
            self.__errors *= other
            return self
        elif isinstance ( other , VE ) : 
            a  = self.__value
            ae = self.__errors
            b  = other.value()
            be = other.error() 
            be = AE ( -be , be )
            self.__value  = a  * b 
            self.__errors = ae * b + be * a 
            return self
        
        return NotImplemented
    
    # =========================================================================
    ## self-scaling
    #  @code
    #  vae = ...
    #  vae /= 5 
    #  @endcode 
    def __idiv__ ( self , other ) :
        """ Self-scaling
        >>> vae = ...
        >>> vae /= 5 
        """
        if isinstance ( other , VE ) and other.cov2() <= 0 : other = float ( other )
        ## 
        if isinstance ( other , num_types ) :
            if iszero ( other ) : return NotImplemented            
            self.__value  /= other
            self.__errors /= other
            return self
        elif isinstance ( other , VE ) :
            if iszero ( other.value() ) : return NotImplemented
            self *= ( 1.0 / other ) 
            return self 
        
        return NotImplemented         
    # =========================================================================
    ## Self-scaling
    __itruediv__ = __idiv__    
    # ==========================================================================
    ## Multiplication/scaling
    #  @code
    #  vae = ...
    #  vae * 5
    #  5 * vae 
    #  @endcode
    def __mul__ ( self , other ) :
        """ Multiplication/scaling
        >>> vae = ...
        >>> vae * 5
        >>> 5 * vae 
        """
        if isinstance ( other , the_types ) :
            r  = copy.copy ( self ) 
            r *= other
            return r
        
        return NotImplemented 
            
    # ==========================================================================
    ## Division/scaling
    #  @code
    #  vae = ...
    #  vae / 5
    #  @endcode
    def __div__ ( self , other ) :
        """ Division/scaling
        >>> vae = ...
        >>> vae * 5
        >>> 5 * vae 
        """
        if not isinstance ( other , the_types ) : return NotImplemented
        r  = copy.copy ( self ) 
        r /= other
        return r 
    # ==========================================================================
    __rmul__    = __mul__
    __truediv__ = __div__ 

    # ==========================================================================
    # add scalar
    def __iadd__ ( self , other ) :
        """ Additivbe update
        """
        if   isinstance ( other , VAE       ) :
            self.__value  += other.value
            self.__errors += other.errors
            return self
        elif isinstance ( other , VE        ) : 
            self.__value  += other.value() 
            if 0 < other.cov2() : 
                oe = other.error() 
                self.__errors += AE ( -oe , oe )
            return self
        elif isinstance ( other , num_types ) : 
            self.__value += other
            return self
            
        return NotImplemented 

    # ===========================================================================
    # subtract scalar
    def __isub__ ( self , other ) :
        """ Subtraction update 
        """
        if   isinstance ( other , VAE       ) :
            self.__value  -=        other.value
            self.__errors += ( -1 * other.errors )
            return self
        elif isinstance ( other , VE        ) : 
            self.__value  -= other.value() 
            if 0 < other.cov2() : 
                oe = other.error() 
                self.__errors += AE ( -oe , oe )
            return self
        elif isinstance ( other , num_types ) : 
            self.__value -= other
            return self

        return NotImplemented 

    # ==========================================================================
    # add two objecte 
    def __add__ ( self , other ) :
        """ Sum of two objects 
        """
        if not isinstance ( other , the_types ) : return NotImplemented        
        r  = copy.copy ( self ) 
        r += other
        return r 
    # ==========================================================================
    # subract scalar
    def __sub__ ( self , other ) :
        """ Difference of two objects 
        """
        if not isinstance ( other , the_types ) : return NotImplemented        
        r  = copy.copy ( self ) 
        r -= other
        return r     
    # ==========================================================================
    # add scalar
    def __radd__ ( self , other ) :
        """ Sum of two objects 
        """
        if not isinstance ( other , the_types ) : return NotImplemented        
        r  = copy.copy ( self ) 
        r += other
        return r 
    # ==========================================================================
    # subtract scalar
    def __rsub__ ( self , other ) :
        """ Difference of two objects 
        """
        if not isinstance ( other , the_types ) : return NotImplemented        
        r  = copy.copy ( self )
        r *= -1
        r += other 
        return r 
    
    # ========================================================================
    ## conversion to float 
    def __float__ ( self ) :
        """ Conversion to float
        """
        return self.__value

    # ========================================================================
    ## (numerical) equality
    def __eq__ ( self , other ) :
        """ (Numerical) equality of object 
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
        """ (Numerical) equality of object 
        """
        if   isinstance ( other , ValWithErrors ) : return not ( self == other )         
        elif isinstance ( other , VE )            : return not ( self == other ) 
        return NotImplemented
    
    # =========================================================================
    ## An `effective error' (as split normal distribution)
    #  @see https://en.wikipedia.org/wiki/Split_normal_distribution
    @property
    def eff_error ( self ) :
        """ An `effective error' (as split normal distribution)
        - see https://en.wikipedia.org/wiki/Split_normal_distribution
        """
        return self.__errors.eff_error
    
    # =========================================================================
    ## An `effective bias' (as split normal distribution)
    #  @see https://en.wikipedia.org/wiki/Split_normal_distribution
    @property
    def eff_bias ( self ) :
        """ An `effective bias' (as split normal distribution)
        - see https://en.wikipedia.org/wiki/Split_normal_distribution
        """
        return self.__errors.eff_bias 

    # =========================================================================
    ## An `effective value' == value + bias (as split normal distribution)
    #  @see https://en.wikipedia.org/wiki/Split_normal_distribution
    @property
    def eff_value  ( self ) :
        """ An `effective value' == value + bias (as split normal distribution)
        - see https://en.wikipedia.org/wiki/Split_normal_distribution
        """
        return self.__value + self.__errors.eff_bias 
    
    # =========================================================================
    ## An `effective error' (as split normal distribution)
    #  @see https://en.wikipedia.org/wiki/Split_normal_distribution
    @property
    def eff_error ( self ) :
        """ An `effective error' (as split normal distribution)
        - see https://en.wikipedia.org/wiki/Split_normal_distribution
        """
        return self.__errors.eff_error 

    # ========================================================================
    ## Conversion to VE, as split normal distribution with optional bias 
    #  @see https://en.wikipedia.org/wiki/Split_normal_distribution
    def asVE ( self , bias = False ) :
        """ Conversion to VE with optional bias
        - see https://en.wikipedia.org/wiki/Split_normal_distribution
        """
        vv = self.value if not bias else self.eff_value 
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
        """ Conversion to string
        """
        return format % ( self.value , abs ( self.neg_error ) , self.pos_error )

    # =========================================================================
    ## pretty print
    #  @code
    #  value = ...
    #  result, expo = value.pretty_printt() 
    #  @endcode
    def pretty_print ( self               ,
                       width       = 6    ,
                       precision   = 4    ,
                       parentheses = True ) :
        """ Pretty print
        >>> value = ...
        >>> result, expo = value.pretty_print () 
        """
        return pretty_vae ( self                      ,
                            width       = width       ,
                            precision   = precision   ,
                            parentheses = parentheses )

    def __str__  ( self ) : return self.toString () 
    def __repr__ ( self ) : return self.__str__  () 

# ===========================================================================
## decode error data into the flat list of errors 
def flat_errors ( *items ) :
    """ Decode error data into the flat list of errors
    """
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
    """ Value with multiple asymmetric error """    
    __slots__ = (  '__value' , '__errors' )
    ## 
    def __init__  ( self , value = 0 , *errors ) :

        self.__errors = [] 
        if   isinstance ( value , ValWithMultiErrors ) :
            
            self.__value   = float ( value.value  ) 
            self.__errors += list  ( value.errors )
            
        elif   isinstance ( value , ValWithErrors ) :
            
            self.__value   = flaot ( value.value  ) 
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
    
    # =========================================================================
    ## make a copy 
    def __copy__ ( self ) :
        """ Make a copy
        """ 
        return ValWithMultiErrors ( self ) 
        
    # =========================================================================
    ## self-scaling
    #  @code
    #  vae = ...
    #  vae *= 5 
    #  @endcode 
    def __imul__ ( self , other ) :
        """ Self-scaling
        >>> vae = ...
        >>> vae *= 5 
        """
        if not isinstance ( other , num_types ) : return NotImplemented
        self.__value  *= other
        self.__errors  = tuple ( e*other for e in self.errors ) 
        return self    
    # =========================================================================
    ## self-scaling
    #  @code
    #  vae = ...
    #  vae /= 5 
    #  @endcode 
    def __idiv__ ( self , other ) :
        """ Self-scaling
        >>> vae = ...
        >>> vae /= 5 
        """
        if not isinstance ( other , num_types ) : return NotImplemented
        self.__value  /= other
        self.__errors  = tuple ( e/other for e in self.errors ) 
        return self    
    # =========================================================================
    ## Self-scaling
    __itruediv__ = __idiv__    
    # ==========================================================================
    ## Multiplication/scaling
    #  @code
    #  vae = ...
    #  vae * 5
    #  5 * vae 
    #  @endcode
    def __mul__ ( self , other ) :
        """ Multiplication/scaling
        >>> vae = ...
        >>> vae * 5
        >>> 5 * vae 
        """
        if not isinstance ( other , num_types ) : return NotImplemented
        r  = copy.copy ( self )
        r *= other
        return r 
    # ==========================================================================
    ## Division/scaling
    #  @code
    #  vae = ...
    #  vae / 5
    #  @endcode
    def __div__ ( self , other ) :
        """ Division/scaling
        >>> vae = ...
        >>> vae * 5
        >>> 5 * vae 
        """
        if not isinstance ( other , num_types ) : return NotImplemented
        r  = copy.copy ( self )
        r /= other
        return r
    
    # ==========================================================================
    # add scalar
    def __iadd__ ( self , other ) :
        """ Additive update 
        """
        if not isinstance ( other , num_types ) : return NotImplemented
        self.__value += other
        return self
    # ===========================================================================
    # subtract scalar
    def __isub__ ( self , other ) :
        """ Subtraction update 
        """
        if not isinstance ( other , num_types ) : return NotImplemented
        self.__value -= other
        return self
    # ==========================================================================
    # add scalar
    def __add__ ( self , other ) :
        """ Sum of two objects 
        """
        if not isinstance ( other , num_types ) : return NotImplemented        
        r  = copy.copy ( self ) 
        r += other
        return r 
    # ==========================================================================
    # subract scalar
    def __sub__ ( self , other ) :
        """ Difference of two objects 
        """
        if not isinstance ( other , num_types ) : return NotImplemented        
        r  = copy.copy ( self ) 
        r -= other
        return r     
    # ==========================================================================
    # add scalar
    def __radd__ ( self , other ) :
        """ Sum of two objects 
        """        
        if not isinstance ( other , num_types ) : return NotImplemented        
        r  = copy.copy ( self ) 
        r += other
        return r 
    # ==========================================================================
    # subtract scalar
    def __rsub__ ( self , other ) :
        """ Difference of two objects 
        """
        if not isinstance ( other , num_types ) : return NotImplemented        
        r  = copy.copy ( self )
        r *= -1
        r += other 
        return r 

    # ==========================================================================
    __rmul__    = __mul__
    __truediv__ = __div__ 
    # ==========================================================================
    ## conversion to float 
    def __float__ ( self ) :
        """ Conversion to float"""
        return self.__value
        
    ## picklings/unpickling 
    def __getstate__ ( self ) :
        return self.__value, self.__errors 
    def __setstate__ ( self , state ) :
        self.__value  = state [ 0 ]
        self.__errors = state [ 1 ]

    @property
    def neg_error ( self ) :
        """'neg_error' : the effective (sum squared) negative/low error """
        return -math.sqrt ( sum ( e.negative**2 for e in self.__errors ) ) 
    @property
    def pos_error ( self ) :
        """'pos_error' : the effective (sum squared) positive/high error"""
        return  math.sqrt ( sum ( e.positive**2 for e in self.__errors ) ) 

    # ========================================================================
    ## (numerical) equality
    def __eq__ ( self , other ) :
        """ (Numerical) equality of object 
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
        """ (Numerical) equality of object 
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
        """ Conversion to VE with optional bias 
        - see https://en.wikipedia.org/wiki/Split_normal_distribution
        """
        vae = self.asVAE ()
        return vae.asVE ( bias = bias )

    # ==========================================================================
    ## conversion to string 
    def toString ( self ,
                   format  = '( %+.5g %s) '    ,   ## global format 
                   format2 = ' -/%-.5g +/%.5g' ) : ## AsymErrors formar 
        """ Conversion to string
        """
        errs = [ e.toString ( format2 ) for e in self.errors ] 
        errs = ','.join ( errs ) 
        return format % ( self.value , errs )

    # =========================================================================
    ## pretty print
    #  @code
    #  value  = ...
    #  result, expo = value.pretty_print () 
    #  @endcode
    def pretty_print   ( self               ,
                         width       = 6    ,
                         precision   = 4    ,
                         parentheses = True ) :
        """ Pretty print
        >>> value = ...
        >>> result, expo = value.pretty_print () 
        """
        return pretty_vme ( self                      ,
                            width       = width       ,
                            precision   = precision   ,
                            parentheses = parentheses )
    
    def __str__  ( self ) : return self.toString () 
    def __repr__ ( self ) : return self.__str__  () 

# =============================================================================
## Pretty prints
# =============================================================================
## Format for  pretty print of asymmetric errors
#  @code
#  ae = AsymErrors ( -5 , 10000 )
#  fmte , expo = fmt_pretty_ae ( ae , width = 8, precision = 6 ) 
#  @endcode
def fmt_pretty_ae ( errors             ,
                    width       = 6    ,
                    precision   = 4    ,
                    parentheses = True )  :
    """ Format for  pretty print of asymmetric errors
    >>> ae = AsymErrors ( -5 , 10000 )
    >>> fmte , expo = fmt_pretty_ae ( ae , width = 8, precision = 6 ) 
    """
    assert isinstance ( errors , AsymErrors ), \
        "Invalid type of `errors`: %s" % type ( errors )

    errs = errors.negative , errors.positive 
    _ , fmte , expo = fmt_pretty_errs ( value     = 0.0       ,
                                        errors    = errs      ,
                                        width     = width     ,
                                       precision = precision )

    fmt = '-/%s +/%s' % ( fmte , fmte )
    ## 
    if parentheses : fmt = '( ' + fmt + ' )'
    return fmt , expo 

# =============================================================================
## Formats for  pretty print of value with asymmetric errors
#  @code
#  vae = ValWithErrors ( -5 , -1, 20 ) 
#  fmt , fmtv , fmte , expo = fmt_pretty_vae ( vae , width = 8, precision = 6 ) 
#  @endcode
def fmt_pretty_vae (  value              ,
                      width       = 6    ,
                      precision   = 4    ,
                      parentheses = True )  :
    """ Formats for  pretty print of value with asymmetric errors
    >>> vae = ValWithErrors ( -5 , -1, 20 ) 
    >>> fmt , fmtv , fmte , expo = fmt_pretty_vae ( vae , width = 8, precision = 6 ) 
    """
    assert isinstance ( value , ValWithErrors ), \
        "Invalid type of `value`: %s" % type ( value)

    v , errlow , errhigh = value.value , value.neg_error , value.pos_error
    return fmt_pretty_err2 ( value       = v           ,
                             errlow      = errlow      ,
                             errhigh     = errhigh     ,
                             width       = width       ,
                             precision   = precision   ,
                             parentheses = parentheses )

# =============================================================================
## Formats for  pretty print of value with multiple asymmetric errors
#  @code
#  vme = ValWithMultiErrors ( ... ) 
#  fmt , fmtv , fmte , expo = fmt_pretty_vme ( vme , width = 8, precision = 6 ) 
#  @endcode
def fmt_pretty_vme ( value              ,
                     width       = 6    ,
                     precision   = 4    ,
                     parentheses = True )  :
    """ Formats for  pretty print of value with multipl asymmetric errors
    >>> vme = ValWithMultipliErrors ( ... ) 
    >>> fmt , fmtv , fmte , expo = fmt_pretty_vme ( vme , width = 8, precision = 6 ) 
    """
    assert isinstance ( value , ValWithMultiErrors ), \
        "Invalid type of `value`: %s" % type ( value )
    
    errors  = tuple ( e.negative for e in value.errors ) 
    errors += tuple ( e.positive for e in value.errors ) 

    fmtv , fmte , expo = fmt_pretty_errs ( value     = value.value ,
                                           errors    = errors      , 
                                           width     = width       ,
                                           precision = precision   )
    fmt = " %s" % fmtv
    for e in value.errors : fmt += " -/%s +/%s"
    if parentheses : fmt = '( ' + fmt + ' )'
    return fmt , fmtv , fmte , expo

# =============================================================================
## pretty print for asymemtric errors
#  @code
#  ae = AsymErrors ( -5 , 10000 )
#  s , expo = pretty_ae ( ae , width = 8, precision = 6 ) 
#  @endcode
def pretty_ae ( errors             ,
                width       = 6    ,
                precision   = 4    ,
                parentheses = True )  :
    """ Pretty print for asymemtric errors
    >>> ae = AsymErrors ( -5 , 10000 )
    >>> s , expo = pretty_ae ( ae , width = 8, precision = 6 ) 
    """
    assert isinstance ( errors , AsymErrors ), \
        "Invalid type of `errors`: %s" % type ( errors )
    
    fmt , expo = fmt_pretty_ae ( errors                    ,
                                 width       = width       ,
                                 precision   = precision   ,
                                 parentheses = parentheses ) 

    values = abs ( errors.negative )  , errors.positive
    
    if expo:
        scale  = 10 ** expo 
        values = tuple ( v / scale for v in values )

    return fmt % values , expo 

# =============================================================================
## pretty print for value with asymmetric errors
#  @code
#  vae = ValWithErrors ( 1 , -5 , 10000 )
#  s , expo = pretty_vae ( vae , width = 8, precision = 6 ) 
#  @endcode
def pretty_vae ( value              ,
                 width       = 6    ,
                 precision   = 4    ,
                 parentheses = True )  :
    """ Pretty print for value  with asymemtric errors
    >>> vae = ValWithErrors ( ... )
    >>> s , expo = pretty_vae ( vae , width = 8, precision = 6 ) 
    """    
    assert isinstance ( value , ValWithErrors ), \
        "Invalid type of `value`: %s" % type ( value )
    
    fmt, _ , _ , expo = fmt_pretty_vae ( value                     , 
                                         width       = width       ,
                                         precision   = precision   ,
                                         parentheses = parentheses ) 
    
    values  = value.value , abs ( value.neg_error ) , value.pos_error    
    if expo:
        scale  = 10 ** expo 
        values = tuple ( v / scale for v in values )

    return fmt % values , expo 

# =============================================================================
## pretty print for value with multiple asymmetric errors
#  @code
#  vme = ValWithMltiErrors ( ... )
#  s , expo = pretty_vme ( vme , width = 8, precision = 6 ) 
#  @endcode
def pretty_vme ( value              ,
                 width       = 6    ,
                 precision   = 4    ,
                 parentheses = True )  :
    """ Pretty print for value  with multiple asymmetric errors
    >>> vme = ValWithMultiErrors ( ... )
    >>> s , expo = pretty_vme ( vme , width = 8, precision = 6 ) 
    """    
    assert isinstance ( value , ValWithErrors ), \
        "Invalid type of `value`: %s" % type ( value)

    fmt, _ , _ , expo = fmt_pretty_vme ( value                     , 
                                         width       = width       ,
                                         precision   = precision   ,
                                         parentheses = parentheses ) 
    
    values = [ vme.value ]
    for e in mve.errors :
        values += [ abs ( e.negative ) , e.possitive ]
        
    values  = tuple ( values )
    if expo:
        scale  = 10 ** expo 
        values = tuple ( v / scale for v in values )

    return fmt % values , expo 

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
