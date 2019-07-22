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
    'primes'  , ## function to get an array of prime numbers ,
    'Primes'  , ## class to store/reuse/handle the prime numbers 
    )
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger    import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.primes' )
else                       : logger = getLogger ( __name__            )
# =============================================================================
import bisect , random  
from   builtins import  range
# =============================================================================
try :
    
    import numpy
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
    def primes(n):
        """ Get the array of prime numbers that do not exceed n
        - numpy (fast?) version
        >>> nums = primes  ( 1000) 
        - see https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n/3035188#3035188
        """
        if   n == 2 : return numpy.array ([] , dtype = atype ) 
        elif n <= 3 : numpy.array ([2]       , dtype = atype ) 
        elif n <= 5 : numpy.array ([2,3]     , dtype = atype ) 

        sieve = numpy.ones(n//3 + (n%6==2), dtype=numpy.bool)
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
        >>> nums = primes  ( 1000) 
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
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================

    
    
