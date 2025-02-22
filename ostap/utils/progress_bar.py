#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Progress bars for Ostap
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date   2011-12-01
#
#  - ProgressBar
# 
#  Use as class:
# 
#  @code
#  count = 0
#  total = 100000
#  bar   = ProgressBar(count, total, 77, mode='fixed', char='#')
#  while count <= total:
#    count += 1
#    bar   += 1 
#  @endcode
#  Better to use as context manager:
#  @code
#  count = 0
#  total = 100000
#  with ProgressBar(count, total, 77, mode='fixed', char='#') as bar  :
#    while count <= total:
#      count += 1
#      bar   += 1 
#  @endcode
#  With helper function:
#  @code 
#  for i in progress_bar  ( range(10000 ) ) :
#      .. do something here ...
#  @endcode 
#
#  - RunningBar
#  @code
#  with RunningBar () as bar :
#    for i in range(2000) :
#           ...
#           bar += 1 
#  @endcode
#  With helper function:
#  @code
#  for i in running_bar  ( range(10000 ) ) :
#      .. do something here ...
#  @endcode 
#
#  ProgressBar is an improvement from the original found at:
#    http://code.activestate.com/recipes/168639/
#
# =============================================================================
""" Simple utilities for `progrees bar'

- ProgressBar

Use as class:

>>> count = 0
>>> total = 100000
>>> bar   = ProgressBar(count, total, 77, mode='fixed', char='#')
>>> while count <= total:
...   count += 1
...   bar   += 1 

Better to use as context manager:

>>> with ProgressBar(count, total, 77, mode='fixed', char='#') as bar  :
>>>    while count <= total:
...       count += 1
...       bar   += 1 

Use with helper function:

>>> for i in progress_bar  ( range(10000 ) ) :
...       <do something here>

- RunningBar

>>> with RunningBar () as bar :
>>>    for i in range(2000) :
...       bar += 1

With helper function:

>>> for i in running_bar  ( range(10000 ) ) :
...    <do something here>


This class is an improvement from the original found at:
@see http://code.activestate.com/recipes/168639/
"""
# =============================================================================
__author__   = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__version__  = "$Revision$"
__date__     = "2011-12-01"
__all__      = (
    "ProgressBar"   , ## Progress bas as it is 
    "RunningBar"    , ## Running bar 
    "progress_bar"  , ## helper function fro ProgressBar 
    "running_bar"     ## helper function for RunningBar 
    )
# =============================================================================
from   ostap.core.ostap_types import sized_types 
from   ostap.utils.basic      import isatty, terminal_size
from   ostap.logger.symbols   import clock_ticks, runner, finish, clock, arrow_right     
from   ostap.logger.colorized import allright, infostr 
import sys , os, time, datetime 
# =============================================================================
## get number of columns for xterm
#  @code
#  ncols = columns()
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
def columns () :
    """ Get number of columns for xterm
    """
    width , height = terminal_size()
    return width

# =============================================================================
## @class ProgressBar
#
#  This class is an improvement from the original found at:
#    http://code.activestate.com/recipes/168639/
#
# =============================================================================
# A Python Library to create a Progress Bar.
# Copyright (C) 2008  BJ Dierkes <wdierkes@5dollarwhitebox.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# =============================================================================
#
#  - ProgressBar
# 
#  Use as class:
#  @code
#  counnt = 0 
#  total  = 100000
#  bar    = ProgressBar(total, 77, mode='fixed', char='#')
#  while count <= total:
#    ...
#    count += 1 
#    bar   += 1 
#  @endcode
#  Or
#  @code
#  total = 100000
#  bat = ProgressBar(total, 77, mode='fixed', char='#') as bar  :
#  while bar : 
#    ...
#    bar += 1 
#  @endcode
#  Better to use as context manager:
#  @code
#  count = 0
#  total = 100000
#  with ProgressBar(count, total, 77, mode='fixed', char='#') as bar  :
#    while count <= total:
#      count += 1
#      bar   += 1 
#  @endcode
#  With helper function:
#  @code 
#  for i in progress_bar  ( range ( 10000 ) ) :
#      .. do something here ...
#  @endcode 
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
class ProgressBar(object):
    """ ProgressBar
    
    Use as class:
    
    >>> count = 0
    >>> total = 100000
    >>> bar   = ProgressBar(total, 77, mode='fixed', char='#')
    >>> while count <= total:
    ...   count += 1
    ...   bar   += 1 


    >>> total = 100000
    >>> bar   = ProgressBar(total, 77, mode='fixed', char='#')
    >>> while bar :
    ...   bar   += 1 

    Better to use as context manager:
    
    >>> with ProgressBar(total, 77, mode='fixed', char='#') as bar  :
    >>>    while count <= total:
    ...       count += 1
    ...       bar   += 1 
    
    Use with helper function:
    
    >>> for i in progress_bar  ( range(10000 ) ) :
    ...       ... do something here ... 
    """
    def __init__( self                   , 
                  max_value = 100        ,
                  width     = 110        ,
                  min_value = 0          ,
                  output    = sys.stdout , **kwargs ) :

        ## tty ? 
        tty = isatty ( output )
        
        self.silent   = kwargs.get ( 'silent' , not tty )
        self.r        = '\r' ## '\r' if tty else '\n'
        
        self.char = kwargs.get ( 'char' , '#'       ) ##
        self.mode = kwargs.get ( 'mode' , 'fixed'   ) ## fixed or dynamic
        if not self.mode in ( 'fixed' , 'dynamic' ) :
            self.mode = 'fixed'
            
        self.bar      = ''
        self.min      = min_value
        self.max      = max_value
        self.span     = max ( max_value - min_value , 1 )
        self.last     = '' 
        ##
        self.output   = output
        ## 
        ncols         = columns () - 12
        self.width    = ncols if ncols > 10 else width
        
        self.prefix   = kwargs.get ( 'description' , '' )  ## description
        self.width    = self.width - len ( self.prefix )
        
        ##  self.amount   = 0    
        self.amount   = self.min ## ??    

        self._hashes  = -1 
        self.__end    = None 
        
        self.update_amount ( self.min )
        self.build_bar ()
        self.show      ()
        self.__start  = time.time ()

    # =========================================================================
    ## Is this bar active?
    def __nonzero__ ( self ) :
        """ Is this bar active? """
        return self.min <= self.amount < self.max
    
    # =========================================================================
    ## Is this bar active?
    def __bool__    ( self ) :
        """ Is this bar active? """
        return self.min <= self.amount < self.max

    # =========================================================================
    ## Incremeemt self.amount
    def increment_amount ( self , add_amount = 1 ) :
        """ Incremeemt self.amount """ 
        return self.update_amount ( self.amount + add_amount )

    # =========================================================================
    ## Update self.amount with 'new_amount', and then rebuild & show the bar
    def update_amount ( self , new_amount = None ):
        """ Update self.amount with 'new_amount', and then rebuild& show the bar
        """
        if new_amount is None   : new_amount = self.amount
        ##
        self.amount = min ( max ( new_amount , self.min ) , self.max )
        ##))
        if self.silent : return self
        ##
        if self.max == self.amount and self.__end is None :self.__end = time.time() 
        ##  
        if self.build_bar() : self.show()
        else :
            dmin  = self.amount - self.min
            dmax  = self.max    - self.amount
            dspan = max ( 20 , 0.001 * self.span ) 
            if dmin < dspan or dmax < dspan : self.show ()
        ##
        return self

    # =========================================================================
    ## Figure new percent complete, and rebuild the bar string base on self.amount.
    def build_bar ( self ) :
        """ Figure new percent complete, and rebuild the bar string base on self.amount.
        """
        diff         = float ( self.amount - self.min )
        done         =  ( diff / float ( self.span ) ) * 100.0
        percent_done = int ( round ( done ) )

        # figure the proper number of 'character' make up the bar 
        all_full     = self.width - 2
        num_hashes   = int ( round ( ( percent_done * all_full ) / 100 ) )

        if 100 <= done and self.__end is None :            
            self.__end = time.time ()
            
        if self.__end is None and num_hashes == self._hashes : return False 
        
        eta  = ''
        leta = len ( eta )
        
        if  self.__end is None and 7 < num_hashes and 1 < done :
            
            now   = time.time ()
            feta  = int ( ( 100 - done ) *  ( now -  self.__start ) / done )
            h , _ = divmod ( feta            , 3600 )
            m , s = divmod ( feta - h * 3600 ,   60 )
            if   h     : eta = 'ETA %02dh%02d\'%02d" ' % ( h , m , s )
            elif m     : eta = 'ETA %02d\'%02d" '      % (     m , s )
            elif s >=1 : eta = 'ETA %02d" '            %           s 
            leta = len ( eta )
            
        elif ( not self.__end is None ) and 6 < num_hashes :
            
            now   = self.__end 
            feta  = int ( now -  self.__start ) 
            h , _ = divmod ( feta            , 3600 )
            m , s = divmod ( feta - h * 3600 ,   60 )
            if   h     : eta = '%02dh%02d\'%02d" ' % ( h , m , s )
            elif m     : eta = '%02d\'%02d" '      % (     m , s )
            elif s >=1 : eta = '%d" '              %           s 
            leta = len ( eta )

        if finish   and 100 <= done : 
            eta   = finish + ' ' + eta
            leta += 2
        elif runner and done < 100 :
            eta   = runner + ' ' + eta
            leta += 2
        
        if self.mode == 'dynamic':
            # build a progress bar with self.char (to create a dynamic bar
            # where the percent string moves along with the bar progress.            
            if eta and leta < num_hashes : 
                self.bar = eta + allright ( self.char * ( num_hashes - leta ) )
            else :
                self.bar = allright ( self.char * num_hashes ) 
        else:
            # build a progress bar with self.char and spaces (to create a 
            # fixed bar (the percent string doesn't move)
            if eta and leta + 1 < num_hashes :
                self.bar = eta + allright ( self.char * ( num_hashes - leta ) ) + ' ' * ( all_full - num_hashes )
            else : 
                self.bar = allright ( self.char * num_hashes ) + ' ' * ( all_full - num_hashes )
 
        percent_str  = str ( percent_done ) + "%"
        
        self.bar     = '[ ' + self.bar + ' ] ' + infostr ( percent_str ) 
        
        self._hashes  = num_hashes
        self._done    = done 

        return True

    # =========================================================================
    ## Increment amount
    def __iadd__ ( self , i ) :
        """ Incremen amount """ 
        return self.increment_amount ( i )

    # ==========================================================================
    ## Show the bar 
    def __str__(self):
        """ Show hthe bar """
        return str  ( self.bar )

    # ==========================================================================
    ## Show the  bar 
    def show ( self ) :
        """ Show the bar """
        if not self.silent and self.bar != self.last :   
            if self.prefix : self.output.write( self.prefix ) 
            self.output.write ( self.bar + self.r ) 
            self.last = self.bar  
            self.output.flush ()
            
    # =========================================================================
    ## The end 
    def end  ( self  ) :
        if not self.silent :
            if self.__end is None : self.__end = time.time () 
            self.build_bar()
            if self.prefix : self.output.write( self.prefix ) 
            self.output.write ( self.bar + '\n' ) 
        self.output.flush()
        self.silent = True

    # =========================================================================
    ## Context manager: ENTER 
    def __enter__ ( self      ) : 
        """ Context manager: ENTER 
        """
        self.__start  = time.time ()
        self.show() 
        return self
        
    # =========================================================================
    ## Context manager EXIT 
    def __exit__  ( self , *_ ) :
        """ Context manager EXIT 
        """
        self.end ()

    def __del__   ( self      ) : self.end ()

# =============================================================================
_bar_  = tuple ( ( runner + ' ' + tick + ' '  + '\r' ) for tick in clock_ticks ) 
_lbar_ = len ( _bar_ )
_done_ = finish + ' ' + infostr  ( 'Done #%-d' ) 
# =============================================================================
## @class RunningBar 
#  - RunningBar
#  @code
#  with RunningBar () as bar :
#    for i in xrange(2000) :
#           ...
#           bar += 1 
#  @endcode
#  With helper function:
#  @code
#  for i in running_bar  ( range(10000 ) ) :
#      .. do something here ...
#  @endcode 
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
class RunningBar(object):
    """ RunningBar
    
    >>> with RunningBar () as bar :
    >>>    for i in xrange(2000) :
    ...       bar += 1
    
    With helper function:
    
    >>> for i in running_bar  ( xrange(10000 ) ) :
    ...    do something here
    
    """
    def __init__ ( self, *fargs , **kwargs ) :
        
        self.silent   = kwargs.get( 'silent' , False ) or not isatty() 

        self.amount   =  0 
        self.freq     =  int ( kwargs.get ( 'frequence' , 100 ) )
        self.prefix   =  kwargs.get ( 'description' , ''     )
        self.__start  =  time.time ()
        self.__last   =  None 
        self.__shown  =  0
        
        self.update_amount() 

        
    def increment_amount(self, add_amount = 1):
        return self if self.silent else self.update_amount ( self.amount + add_amount )

    def update_amount(self, new_amount = None ):
        """Update self.amount with 'new_amount', and then rebuild the bar string.
        """
        if self.silent : return self ## really silent 
        #
        if not new_amount: new_amount = self.amount
        self.amount = new_amount 
        ##
        if not self.silent : self.show() 
        ##
        return self
    
    def __iadd__ ( self , i ) :
        return self if self.silent else self.increment_amount ( i )
    
    def __str__(self):
        return str(self.bar)

    def show ( self ) :

        if self.silent : return

        if   0 < self.freq          and 0 == self.amount % self.freq : return self.show_ ()
        elif self.amount <=      10                                  : return self.show_ () 
        elif self.amount <=     100 and 0 == self.amount % 10        : return self.show_ () 
        elif self.amount <=    1000 and 0 == self.amount % 50        : return self.show_ () 
        elif self.amount <=   10000 and 0 == self.amount % 100       : return self.show_ () 
        elif self.amount <= 1000000 and 0 == self.amount % 500       : return self.show_ ()
        elif 0 == self.amount % 1000                                 : return self.show_ () 

    def show_ ( self , now = None ) :

        if   self.silent : return        
        elif now is None : now = time.time() 

        ## avoid very frequent printout        
        if self.__last and now - self.__last < 0.1 : return
        
        ## index 
        bar = _bar_ [ self.__shown % _lbar_ ]
        
        sys.stdout.write ( self.prefix + _bar_ [ self.__shown % _lbar_ ] )
        sys.stdout.flush ()

        self.__shown += 1
        self.__last   = now 

        return
    
    def end  ( self ) :
        if self.silent : return
        start   = datetime.datetime.fromtimestamp ( self.__start )
        end     = datetime.datetime.now ()
        delta   = end - start

        hours    , reminder = divmod ( delta.seconds , 3600 )
        minutes  , seconds  = divmod ( reminder      ,   60 )
        days     = delta.days 
        secs     = seconds + 1.e-6 * delta.microseconds
        
        if   days    : timebar = ' %d days %02dh%02d\'%04.1f"'  % ( delta.days , hours , minutes , secs )
        elif hours   : timebar = ' %02dh%02d\'%04.1f"'          % ( hours , minutes , secs )
        elif minutes : timebar = ' %02d\'%04.1f"'               % (         minutes , secs )
        else         : timebar = ' %.2f"'                       %                     secs 

        timebar += ' ' + start.strftime ('%a %d %b %Y, %I:%M%p')
        timebar += ' ' + arrow_right 
        timebar += ' ' + end.strftime   ('%a %d %b %Y, %I:%M%p')
        
        sys.stdout.write ( self.prefix + ( _done_ % self.amount ) + ' ' + clock + timebar + '\n' ) 
        sys.stdout.flush ()
        self.silent = True
        
    ## CONTEXT manager ENTER 
    def __enter__ ( self      ) :
        self.__start  = time.time ()
        self.__last   = None 
        self.__shown  = 0 
        return self
    
    ## CONTEXT manager EXIT 
    def __exit__  ( self , *_ ) : self.end()
    def __del__   ( self , *_ ) : self.end()

# ==============================================================================
## helper function to display running bar 
#  @code
#  for i in running_bar  ( range(10000 ) ) :
#      .. do something here ...
#  @endcode 
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
def running_bar ( iterable , frequency = 100 , description = '' , **kwargs ) :
    """ Helper function to display runnning bar 
    >>> for i in running_bar  ( xrange(10000 ) ) :
    ...    do something here
    """
    with RunningBar ( frequency = frequency , description = description , **kwargs ) as bar :
        for i in iterable :
            bar += 1
            yield i

# =============================================================================
## helper function to display progress bar
#  @code 
#  for i in progress_bar  ( xrange(10000 ) ) :
#      .. do something here ...
#  @endcode 
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
def progress_bar ( iterable , max_value = None , **kwargs ) :
    """ Helper function to display progress bar 
    
    >>> for i in progress_bar  ( range ( 10000 ) ) :
    ...      do something here
    """
    if isinstance ( iterable , int ) and  0 <= iterable  :
        iterable = range ( iterable )
        
    if   max_value is None and isinstance ( iterable , sized_types ) :
        max_value = len ( iterable ) 
    elif max_value is None and hasattr    ( iterable , '__len__'   ) :
        max_value = len ( iterable )
    elif max_value is None and hasattr    ( iterable , 'size'      ) :
        max_value = iterable.size()

    if   max_value is None : bar = RunningBar  ( **kwargs )
    elif max_value <  1    : bar = RunningBar  ( **kwargs )
    else                   : bar = ProgressBar ( max_value = max_value , **kwargs ) 

    with bar :
        bar.show () 
        for i in iterable :
            yield i
            bar += 1
                        
# =============================================================================
## simple test 
def test_bars ():

    limit = 1000 
    
    import time
    
    print('Example 1: Fixed Bar')
    with ProgressBar ( limit ,  mode = 'fixed' ) as bar : 
        for i in range ( limit + 1 ):
            bar += 1 
            time.sleep ( 0.005 )
 
    print('Example 2: Dynamic Bar')
    with ProgressBar( limit , mode = 'dynamic' , char = '-' ) as bar : 
        for i in range ( limit + 1 ):
            bar += 1 
            time.sleep ( 0.005 )

    for i in progress_bar ( range ( 5000 ) , description = "Doing something ") : 
        time.sleep(0.001)
        
    for i in running_bar  ( range ( 5000 ) , description = "Empty looping ") :         
        time.sleep(0.001)

# ==============================================================================
if __name__ == '__main__':

    from   ostap.logger.logger import getLogger
    logger = getLogger('ostap.utils.progress_bar')
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
  
    test_bars ()
    logger.info ( 80*'*' ) 
    
# =============================================================================
##                                                                      The END 
# =============================================================================

