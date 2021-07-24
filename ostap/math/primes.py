#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/primes.py
#  Get prime numbers
#  @code
#  np = primes  ( 1000 )##  get all prime numbers that are smaller than 1000
#  @endcode
#  The function <code>primes</code> use sieve algorithm to get
#  the prime numbers and it could be relatively slow.e.g. for N>2**32 
#  A bit more effecient version with reusage of the prime numbers
#  is by using of the class Primes:
#  @code 
#  p  = Primes  ( 1000000 )
#  for n in p :  ...
#  @endcode
#  The class <code>Primes</code> stores and reuses the prime numbers when possible,
#  thus allowing a bit more CPU efficient treatment.
# 
#  @see https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n/3035188#3035188
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2016-07-15
# =============================================================================
""" Get prime numbers 
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2009-09-12"
__version__ = "Version $Revision:$"
# =============================================================================
__all__     = (
    'primes'            , ## function to get an array of prime numbers ,
    'Primes'            , ## class to store/reuse/handle the prime numbers
    'idivisors'         , ## generator to get all divisors for the given number    
    'divisors'          , ##                  all divisors for the given number
    'all_prime_factors' , ## get all prime factors for the given number 
    'prime_factors'     , ## get all non-repetitive prime factors for the given number 
    )
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger    import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.primes' )
else                       : logger = getLogger ( __name__            )
# =============================================================================
import bisect , random  
from   builtins import range
# =============================================================================
try :
    
    import numpy
    npvers   = tuple ( numpy.version.version.split('.') )
    npb_type = bool if ( '1', '20') <= npvers else numpy.bool 
                     
    atype = numpy.uint64
    # =========================================================================
    ## get the array of prime numbers that do not exceed <code>n</code>
    #  - <code>numpy</code> (fast?) version 
    #  @code
    #  nums = primes  ( 1000) 
    #  @endcode
    #  @param n the upper edge
    #  @return the array o fprime numbers that do not exceed  <code>n</code>
    #  @see https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n/3035188#3035188
    def primes ( n ):
        """ Get the array of prime numbers that do not exceed n
        - numpy (fast?) version
        >>> nums = primes  ( 1000) 
        - see https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n/3035188#3035188
        """
        if   n == 2 : return numpy.array ([] , dtype = atype ) 
        elif n <= 3 : numpy.array ([2]       , dtype = atype ) 
        elif n <= 5 : numpy.array ([2,3]     , dtype = atype ) 

        sieve = numpy.ones(n//3 + (n%6==2), dtype = npb_type )
        for i in range(1,int(n**0.5)//3+1):
            if sieve[i]:
                k=3*i+1|1
                sieve[       k*k//3     ::2*k] = False
                sieve[k*(k-2*(i&1)+4)//3::2*k] = False
        return numpy.r_[2,3,((3*numpy.nonzero(sieve)[0][1:]+1)|1)]
    
except ImportError :

    # =========================================================================
    ## get the array of prime numbers that do not exceed <code>n</code>
    #  - pure python (slow?) version 
    #  @code
    #  nums = primes  ( 1000) 
    #  @endcode
    #  @param n the upper edge
    #  @return the array o fprime numbers that do not exceed  <code>n</code>
    #  @see https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n/3035188#3035188
    def primes ( n ):
        """ Get the array of prime numbers that do not exceed n
        - pure python  (slow?) version
        >>> nums = primes  ( 1000 ) 
        - see https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n/3035188#3035188
        """
        if   n == 2 : return []
        elif n <= 3 : return [2]
        elif n <= 5 : return [2,3]
        
        n, correction = n-n%6+6, 2-(n%6>1)
        sieve = [True] * (n//3)
        for i in range(1,int(n**0.5)//3+1):
            if sieve[i]:
                k=3*i+1|1
                sieve[      k*k//3      ::2*k] = [False] * ((n//6-k*k//6-1)//k+1)
                sieve[k*(k-2*(i&1)+4)//3::2*k] = [False] * ((n//6-k*(k-2*(i&1)+4)//6-1)//k+1)
        return [2,3] + [3*i+1|1 for i in range(1,n//3-correction) if sieve[i]]

    logger.warning ("primes: Numpy can't be imported, pure python/(slow?) version will be used")

# =============================================================================
## @class Primes
#  produce and store the primary numbers
#  @code
#  p = Primes(1000)
#  for q in p                     : print q  ## loop over all primes 
#  for q in p.range ( 500 , 700 ) : print q  ## loop over primes betee 500 and 700
#  qs  = p[:20]    ##  get teh first twenty primes
#  qs  = p[10:20]  ##  get the range of primes
#  q   = p.choice (  300 , 800 )  ## get random prime betwewn 300 and 800
#  @endcode
class Primes(object) :
    """Produce and store the primary numbers
    >>> p = Primes(1000)
    >>> for q in p                     : print q  ## loop over all primes 
    >>> for q in p.range ( 500 , 700 ) : print q  ## loop over primes betee 500 and 700
    >>> qs  = p[:20]    ##  get teh first twenty primes
    >>> qs  = p[10:20]  ##  get the range of primes
    >>> q   = p.choice (  300 , 800 )  ## get random prime betwewn 300 and 800
    """
    
    __primes = [] ##  the storage of prime numbers 
    __last   = 0  ##  the (current) upper edge 

    # =========================================================================
    ## create the instance
    #  Note, that the actual storage of prime numbers is shared between
    #  the instances,. allowing to reuse already calculated numbers 
    def __init__ ( self , n = 1000000 ) :

        assert isinstance ( n , int ) and 1 < n , 'Invalid upper edge for prime numbers '
        
        if not self.__primes or self.__last < n : 
            self.__primes = primes ( n )
            self.__last   = n
            
        self.__N = n  

    # =========================================================================
    @property
    def primes ( self ) :
        """``primes'' :  currently known array of primes"""
        return self.__primes

    # =========================================================================
    @property
    def upper  ( self ) :
        """``last''  : get the upper edge for produced prime numbers"""
        return self.__last

    # =========================================================================
    @property
    def N      ( self ) :
        """``N'' :  get the current boundary of the prime numbers"""
        return self.__N

    # =========================================================================
    ## Is this number already in the list of primes ?
    #  @code
    #  primes = Primes  ( ... )
    #  if 234133 in primes : ... 
    #  @endcode
    def __contains__ ( self ,  n ) :
        """Is this number already in the list of primes ?
        >>> primes = Primes  ( ... )
        >>> if 234133 in primes : ... 
        """
        if not isinstance ( n , int ) or n <= 1 : return False

        ## extend the list of primes if needed 
        if n > self.__last :
            self.__primes = primes ( 2 * n )
            self.__last   = n

        i = bisect.bisect_left ( self.__primes , n )
        #
        return  i != len ( self.__primes ) != i and n == self.__primes[i]
        
    # =========================================================================
    ## iterator over primary numbers
    #  @code
    #  primes = Primes ( 1000 )
    #  for n in primes : ...
    #  @endcode
    def __iter__ ( self ) :
        """Iterator  over primary numbers between [2,N]        
        >>> primes = Primes ( 1000 )
        >>> for n in primes : ...
        """
        for n in self.__primes :
            if n < self.__N : yield n
            else            : break 

    # =========================================================================
    ## iterator over the range of primary numbers in the range (min,max)
    #  @code
    #  primes = Primes ( 10000 )
    #  for n in primes.range ( 300 , 500 ) : ...
    #  @endcode
    def range  (  self , min , max ) :
        """Iterator  over primary numbers between [min,max]
        >>> primes = Primes ( 10000 )
        >>> for n in primes.range ( 300 , 500 ) : ...
        """
        bmin = bisect.bisect_left  ( self.__primes , min )
        bmax = bisect.bisect_left  ( self.__primes , max )
        for i  in range ( bmin , bmax ) :
            yield self.__primes[i]

    # =========================================================================
    ## Get random primary number in the specified range
    #  @code
    #  primes = Primes ( 10000 )
    #  p      = primes.choice()
    #  p      = primes.choice( min = 300 )
    #  p      = primes.choice( min = 300 , max = 500 )    
    #  @endcode
    def choice ( self , min = None , max = None ) :
        """Get random primary number in the specified range
        >>> primes = Primes ( 10000 )
        >>> p      = primes.choice()
        >>> p      = primes.choice( min = 300 )
        >>> p      = primes.choice( min = 300 , max = 500 )    
        """
        if ( min is None or min <= self.__primes[0] ) and \
               ( max is None or max > self.__primes [-1] ) :
            return random.choice ( self.__primes ) 
        
        if min is None : bmin = 0
        else           : bmin = bisect.bisect_left ( self.__primes , min )
        
        bmax = bisect.bisect_left ( self.__primes , max )
        
        return random.choice ( self [ bmin : bmax] ) 
        
    # =========================================================================
    ## get primary number with certain index or slice
    #  @code
    #  p = Primes ( 4000 )
    #  n  = p[100]
    #  s  = p[-5:] 
    #  @endcode
    def __getitem__ ( self , index ) :
        """Get primary number with certain index or slice
        >>> primes = Primes ( 4000 )
        >>> n  = primes [100]
        >>> s  = primes [-5:]
        """
        return self.__primes [ index ]

# =============================================================================
## Generator to get ALL divisors for the given number
#  @code
#  for d in idivisors ( 100 ) : print d 
#  @endcode
#  The code is taken from:
#  @see https://stackoverflow.com/questions/171765/what-is-the-best-way-to-get-all-the-divisors-of-a-number
def idivisors ( n ) :
    """ Generator to get ALL divisors for the given number
    >>> for d in idivisors ( 100 ) :
    ...   print d 
    The code is taken from:
    - see https://stackoverflow.com/questions/171765/what-is-the-best-way-to-get-all-the-divisors-of-a-number
    """
    # get factors and their counts
    factors = {}
    nn = n
    i = 2
    imax    = math.sqrt ( n ) 
    while i <= imax :
        while nn % i == 0:
            factors[i] = factors.get(i, 0) + 1
            nn //= i
        i += 1
    if nn > 1:
        factors[nn] = factors.get(nn, 0) + 1

    primes_ = list ( factors.keys() )

    # generates factors from primes[k:] subset
    def generate(k):
        if k == len(primes_):
            yield 1
        else:
            rest = generate(k+1)
            prime = primes_[k]
            for factor in rest:
                prime_to_i = 1
                # prime_to_i iterates prime**i values, i being all possible exponents
                for _ in range(factors[prime] + 1):
                    yield factor * prime_to_i
                    prime_to_i *= prime

    # in python3, `yield from generate(0)` would also work
    for factor in generate(0):
        yield factor

# =============================================================================
## Get ALL divisors for the given number
#  @code
#  d100 =  divisors ( 100 )
#  @endcode
#  The code is taken from:
#  @see https://stackoverflow.com/questions/171765/what-is-the-best-way-to-get-all-the-divisors-of-a-number
def divisors ( n ) :
    """ Get ALL divisors for the given number
    >>> d100 = idivisors ( 100 ) 
    The code is taken from:
    - see https://stackoverflow.com/questions/171765/what-is-the-best-way-to-get-all-the-divisors-of-a-number
    """
    return tuple ( idivisors ( n ) )


# =============================================================================
## get all prime factors for the given number
#  @code
#  factors = all_prime_factors ( 12345 ) 
#  @endcode
#  It uses precalculated list of prime numbers, thefore for the
#  large values and the first call it coudl be very slow 
def all_prime_factors ( n ) :
    """ Get all prime factors for the given number:
    >>> factors = all_prime_factors ( 12345 ) 
    - It uses the precalculated list of prime numbers, thefore for the
    large values and the first call it could be rather slow 
    """
    
    _primes = Primes ( n )

    factors = []

    _n  = n 
    for p in _primes :
        
        if   p == _n :
            factors.append ( p )
            break 
        elif p >  _n : break
        
        a , b = divmod ( _n , p ) 
        while 0 == b :
            factors.append ( p )
            _n = a
            a , b = divmod ( _n , p )

    return tuple ( factors ) 
            
            
        
# =============================================================================
## get all non-repetitive prime factors for the given number
#  @code
#  factors = prime_factors ( 12345 ) 
#  @endcode
#  It uses precalculated list of prime numbers, thefore for the
#  large values and the first call it coudl be very slow 
def prime_factors ( n ) :
    """ Get all non-repetitive prime factors for the given number:
    >>> factors = all_prime_factors ( 12345 ) 
    - It uses the precalculated list of prime numbers, thefore for the
    large values and the first call it could be rather slow 
    """
    
    _primes = Primes ( n )

    factors = []

    _n = n 
    for p in _primes :        

        if   p == _n :
            factors.append ( p )
            break 
        elif p >  _n : break
        
        a , b = divmod ( _n , p ) 
        if 0 == b : factors.append ( p )
        while 0 == b :
            _n = a
            a , b = divmod ( _n , p )
            

    return tuple ( factors ) 
            
            
        
    
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================

    
    
