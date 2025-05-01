#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
# Utilies  to play woth random seeds 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-02-10
# =============================================================================
""" Utilies  to play woth random seeds 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2013-02-10"
# =============================================================================
__all__     = (
    ##
    'RandomSeed'         , ## context manager to set/keep seed for the PYTOHN random generator
    'randomSeed'         , ## context manager to set/keep seed for the PYTHON random generator
    'random_seed'        , ## context manager to set/keep seed for the PYTHON random generator
    ## 
    'RootRandomSeed'     , ## context manager to set/keep seed for the ROOT random generator
    'rootRandomSeed'     , ## context manager to set/keep seed for the ROOT random generator
    'root_random_seed'   , ## context manager to set/keep seed for the ROOT random
    ## 
    'RooRandomSeed'      , ## context manager to set/keep seed for the RooFit random generator
    'rooRandomSeed'      , ## context manager to set/keep seed for the RooFit random generator
    'roo_random_seed'    , ## context manager to set/keep seed for the RooFit random generator
    ## 
    )
# =============================================================================
import random 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.core.random_seed' )
else                       : logger = getLogger( __name__                 )
# =============================================================================
## @class RandomSeed 
#  Context manager to set/keep seed/state of PYTHON random generator
#  @code
#  import random 
#  with RandomSeed ( 'this is seed' ) :
#  ... for i in range ( 10 ) : print ( i, random.uniform (0,1) )
#  @endcode
#  Seed can be None, int,float,string, bytes, .
class RandomSeed(object) :
    """ Context manager to set/keep seed/state of PYTHON random generator
    >>> import random 
    >>> with RandomSeed ( 'this is seed' ) :
    >>> ... for i in range ( 10 ) : print ( i, random.uniform (0,1) )
    - Seed can be None, int,float,string, bytes, .
    """
    def __init__ ( self , seed = None )  :
        
        self.__seed  = seed
        self.__state = None  

    def __enter__  ( self ) :
        
        self.__state = random.getstate ()
        random.seed ( self.__seed )
        return self
    
    def __exit__   ( self , *_ ) :
        random.setstate ( self.__state )
        
    @property
    def seed ( self ) :
        """'seed' : the seed for random
        """
        return self.__seed 
    
# ==============================================================================
## Context manager to set/keep seed/state of random generator
#  @code
#  import random 
#  with randomSeed ( 'this is seed' ) :
#  ... for i in range ( 10 ) : print ( i, random.uniform (0,1) )
#  @endcode
#  Seed can be None, int,float,string, bytes, ... 
def randomSeed ( seed = None ) :
    """Context manager to set/keep seed/state of random generator
    >>> import random 
    >>> with randomSeed ( 'this is seed' ) :
    >>> ... for i in rnage ( 10 ) : print ( i, random.uniform (0,1) )
    - Seed can be None, int,float,string, bytes, ...
    """
    return RandomSeed ( seed )

# ==============================================================================
## Context manager to set/keep seed/state of random generator
#  @code
#  import random 
#  with random_seed ( 'this is seed' ) :
#  ... for i in range ( 10 ) : print ( i, random.uniform (0,1) )
#  @endcode
#  Seed can be None, int,float,string, bytes, ...
random_seed = randomSeed     

# ==============================================================================
## @class RootRandomSeed 
#  Context manager to set/keep seed/state of ROOT random generator
#  @code
#  import random 
#  with RootRandomSeed ( 'this is seed' ) :
#  ... for i in range ( 10 ) : print ( i, ROOT.gRandom.Rndm () )
#  @endcode
#  Seed can be None, int,float,string, bytes, .
class RootRandomSeed(object) :
    """ Context manager to set/keep seed/state of ROOT random generator
    >>> import random 
    >>> with RootRandomSeed ( 'this is seed' ) :
    >>> ... for i in range ( 10 ) : print ( i, ROOT.gRandom.Rndm() )
    - Seed can be None, int,float,string, bytes or any hashable object 
    """
    def __init__ ( self , seed = None )  :
        
        if   seed is None              : self.__seed = 0
        elif isinstance ( seed , int ) : self.__seed = seed
        else                           : self.__seed = hash ( seed ) 
        
        self.__old  = 0   
        
    def __enter__  ( self ) :

        if ROOT.gRandom :
            self.__old = ROOT.gRandom.GetSeed()
            ROOT.gRandom.SetSeed ( self.__seed ) 
            
        return self
    
    def __exit__   ( self , *_ ) :
        if ROOT.gRandom :
            ROOT.gRandom.SetSeed ( self.__old ) 
        
    @property
    def seed ( self ) :
        """'seed' : the actual seed for random
        """
        return self.__seed 

# ==============================================================================
## Context manager to set/keep seed/state of ROOT random generator
#  @code
#  import random 
#  with rootRandomSeed ( 'this is seed' ) :
#  ... for i in range ( 10 ) : print ( i, ROOT.gRandom.Rndm () )
#  @endcode
#  Seed can be None, int,float,string, bytes, .
def rootRandomSeed ( seed = None ) :
    """Context manager to set/keep seed/state of ROOT random generator
    >>> import random 
    >>> with rootRandomSeed ( 'this is seed' ) :
    >>> ... for i in range ( 10 ) : print ( i, ROOT.gRandom.Rndm() )
    - Seed can be None, int,float,string, bytes or any hashable object 
    """
    return RootRandomSeed ( seed ) 

# ==============================================================================
## Context manager to set/keep seed/state of ROOT random generator
#  @code
#  import random 
#  with rootRandomSeed ( 'this is seed' ) :
#  ... for i in range ( 10 ) : print ( i, ROOT.gRandom.Rndm () )
#  @endcode
#  Seed can be None, int,float,string, bytes or any hashable object 
root_random_seed = rootRandomSeed



# ==============================================================================
## @class RooRandomSeed 
#  Context manager to set/keep seed/state of RooFit random generator
#  @code
#  import random 
#  with RooRandomSeed ( 'this is seed' ) :
#  ... for i in range ( 10 ) : print ( i, ROOT.RooRandom.uniform () )
#  @endcode
#  Seed can be None, int,float,string, bytes, .
class RooRandomSeed(object) :
    """ Context manager to set/keep seed/state of ROOT random generator
    >>> import random 
    >>> with RooRandomSeed ( 'this is seed' ) :
    >>> ... for i in range ( 10 ) : print ( i, ROOT.RooRandom.uniform() )
    - Seed can be None, int,float,string, bytes or any hashable object
    """

    def __init__ ( self , seed = None )  :
        
        if   seed is None              : self.__seed = 0
        elif isinstance ( seed , int ) : self.__seed = seed
        else                           : self.__seed = hash ( seed ) 
        
        self.__old  = 0   
        

    def __enter__  ( self ) :
        
        rg = ROOT.RooRandom.randomGenerator() 
        if rg :
            self.__old = rg.GetSeed()
            rg.SetSeed ( self.__seed ) 
            
        return self
    
    def __exit__   ( self , *_ ) :
        
        rg = ROOT.RooRandom.randomGenerator() 
        if rg : rg.SetSeed ( self.__old ) 
        
    @property
    def seed ( self ) :
        """'seed' : the actual seed for random
        """
        return self.__seed 

# ==============================================================================
## Context manager to set/keep seed/state of RooFit random generator
#  @code
#  import random 
#  with rooRandomSeed ( 'this is seed' ) :
#  ... for i in range ( 10 ) : print ( i, ROOT.RooRandom.uniform () )
#  @endcode
#  Seed can be None, int,float,string, bytes, .
def rooRandomSeed ( seed = None ) :
    """ Context manager to set/keep seed/state of RooFit random generator
    >>> import random 
    >>> with rooRandomSeed ( 'this is seed' ) :
    >>> ... for i in range ( 10 ) : print ( i, ROOT.RooRandom.uniform() )
    - Seed can be None, int,float,string, bytes or any hashable object 
    """
    return RooRandomSeed ( seed ) 

# ==============================================================================
## Context manager to set/keep seed/state of RooFit random generator
#  @code
#  import random 
#  with roo_random_seed ( 'this is seed' ) :
#  ... for i in range ( 10 ) : print ( i, ROOT.RooRandom.uniform () )
#  @endcode
#  Seed can be None, int,float,string, bytes or any hashable object 
roo_random_seed = rooRandomSeed


# =============================================================================
if '__main__' == __name__ :

    # logging 
    from ostap.logger.logger import getLogger 
    logger = getLogger( 'ostap.fitting.random_seed' )
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    logger.info ( 80*'*' ) 
            
# =============================================================================
##                                                                      The END 
# =============================================================================
