#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file valerrors.py
#  Helper structures (point-with-errors) for (T)Graphs
#
#  These objects are *NOT* for math!
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07 
# =============================================================================
""" Helper structures (point-with-erorrs) for (T)Graphs 
- These objects are *NOT*  for math!
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
#  @class AsymErrors
#  Asymmetric errors 
class AsymErrors (object) :
    """Asymmetric Errors
    """
    __slots__ = ( '__positive' , '__negative' )
    ##
    def __init__ ( self , negative = 0 , positive = 0 ) :

        assert isinstance ( negative , num_types ) , "Invalid type of 'negative'"
        assert isinstance ( positive , num_types ) , "Invalid type of 'positive'"
        
        pos = float ( positive )
        neg = float ( negative )
        
        if   neg <= 0 <= pos       : self.__negative , self.__positive =  neg , pos 
        elif pos <= 0 <= neg       : self.__negative , self.__positive =  pos , neg   
        elif 0 <= neg and 0 <= pos : self.__negative , self.__positive = -neg , pos 
        else :
            raise TypeError ( "Invalid setting of errors: %s/%s" % ( negative , positive ) )
        
    @property
    def negative ( self ) :
        """'negative' : get the negative (low) error"""
        return self.__negative 
    @property
    def positive ( self ) :
        """'positive' : get the positive (high) error"""
        return self.__positive
    
    ## picklings/unpickling 
    def __getstate__ ( self ) :
        return self.__negative, self.__positive    
    def __setstate__ ( self , state ) :
        self.__negative = state [ 0 ]
        self.__positive = state [ 1 ]

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

    ## Quadratic sum of two AsymErrors object 
    def __add__ ( self , other ) :
        """Quandratic sum of two AsymErrors object"""
        if not isinstance ( other , AsymErrors ) : return NotImplemented
        result = AsymErrors ( self.__negative , slef.__positive )
        result += other
        return result
    
    ## tuple-like 
    def __len__ ( self ) : return 2
    ## get the error by index
    #  0 : negative error
    #  1 : positive error 
    def __getitem__ ( self , index ) :
        """ get the error by index
        - 0 : negative error
        - 1 : positive error 
        """
        if not 0 <= index <2 : raise IndexError("Invalid index %s!" % index ) 
        return self.__negative if 0 == index else self.__positive
    
    def __str__ ( self ) :
        return  "( -/%.5g +/%.5g )" % ( self.__negative , self.__positive )
    ##def __repr__ ( self ) :
    ##    return  "AsymErrors( negative=%+.6g, positive=%+.6g )" % ( self.__negative , self.__positive )
    def __repr__ ( self ) : return self.__str__ () 

# =============================================================================
## @class ValWithErrors
#  Value with asymmetric error 
class ValWithErrors(object) : 
    """Value with asymmetric error"""    
    __slots__ = (  '__value' , '__errors' )
    ## 
    def __init__  ( self , value = 0 , *errors ) :
        
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

        else  : raise TypeError ( "Invalid 'value' %s" % str ( value ) ) 
                
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

    def __str__ ( self ) :
        return  "( %.5g -/%.5g +/%.5g )" % \
               ( self.__value , abs ( self.neg_error ) , self.pos_error )
    ## def __repr__ ( self ) :
    ##    return  "ValWithErrors(%+.6g, ( %+.6g , %+.6g ))" % \
    ##           ( self.__value , self.neg_error, self.pos_error )
    def __repr__ ( self ) : return self.__str__ () 


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


        ## decore errors parameters 
        errs = flat_errors ( *errors )
        ## convert flat list into the list of AsymErrors 
        es   = ( AsymErrors (  errs[i] , errs[i+1] ) for i in range ( 0 , len ( errs ) , 2 ) )
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

    def __str__ ( self ) :

        fragments = [ '( %.5g' % self.__value ]
        for e in self.__errors :
            fragments += [  '-/%.5g +/%.5g' % ( e.negative , e.positive ) ]
        fragments.append ( " )" )                    
        return  ''.join ( fragments ) 

    ## def __repr__ ( self ) :

    ##     fragments = [] 
    ##     for e in self.__errors :
    ##         fragments += [  '(%+.6g,%+.6g)' % ( e.negative , e.positive ) ]
    ##     errors = ','.join ( fragments ) 

    ##     fragments = [ 'ValWithMultiErrors( %+.6g, [ ' % self.__value ]
    ##     fragments.append ( errors )                    
    ##     fragments.append ( " ] )" )
        
    ##     return  ''.join ( fragments ) 


## shortcut 
VAE = ValWithErrors


## shortcut 
VME = ValWithMultiErrors

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
