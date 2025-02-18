#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/minuit.py
#  Module with decoration of some (T)Minuit functions for efficient use in python
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Decoration for some (T)Minuit functions for efficient use in python
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = () 
# =============================================================================
from   ostap.core.ostap_types import integer_types, string_types
from   ostap.core.core        import cpp, VE
from   ostap.utils.valerrors  import VAE 
from   ostap.logger.colorized import allright, attention
import ROOT, math, ctypes
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger      import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.minuit' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Some useful decorations for (T)Minuit functions')
# =============================================================================
partypes = integer_types
# =============================================================================
## return codes from MINUIT commands
return_codes = {
    0  : allright  ( 'command executed normally'                         ) ,
    1  : attention ( 'command is blank, ignored'                         ) ,
    2  : attention ( 'command line unreadable, ignored'                  ) ,
    3  : attention ( 'unknown command, ignored'                          ) ,
    4  : attention ( 'abnormal termination (e.g., MIGRAD not converged)' ),
    5  : 'command is a request to read PARAMETER definitions' , 
    6  : "'SET INPUT' command"   ,
    7  : "'SET TITLE' command"   ,
    8  : "'SET COVAR' command"   ,
    9  : 'reserved'              ,
    10 : 'END command'           ,
    11 : 'EXIT or STOP command'  ,
    12 : 'RETURN command'        ,
    }
# =============================================================================
## get the parameter from (T)Minuit
#  @code
#  mn = ...             # TMinuit object
#
#  p1 = mn[0]                # get the parameter by index 
#  p1 = mn.par(0)            # ditto 
#  p1 = mn.parameter(0)      # ditto
#
#  p2 = mn['par2']           # get the parameter by name 
#  p2 = mn.par('par2')       # ditto 
#  p2 = mn.parameter('par2') # ditto
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-09-28
def _mn_par_ ( self , i ) :
    """Get the parameter from (T)Minuit

    >>> mn = ...             # TMinuit object
    >>> p1 = mn[0]                # get the parameter by index 
    >>> p1 = mn.par(0)            # ditto 
    >>> p1 = mn.parameter(0)      # ditto
    
    >>> p2 = mn['par2']           # get the parameter by name 
    >>> p2 = mn.par('par2')       # ditto 
    >>> p2 = mn.parameter('par2') # ditto
    
    """
    if not i in self : raise IndexError
    if isinstance ( i , string_types ) : i = _mn_index_ ( self , i )
    
    #
    val = ctypes.c_double ( 0 )  
    err = ctypes.c_double ( 0 ) 
    #
    res = self.GetParameter ( i , val , err )
    #
    val = float ( val.value )
    err = float ( err.value )
    #
    return VE ( val , err * err )

# =============================================================================
## Is given parameter knows to MINUIT ?
#  @code
#  minuit = ...
#  if 5    in minuit : ...
#  if 'p6' in minuit : ... 
#  @endcode 
def _mn_contains_ (  self , p ) :
    """Is given parameter known to MINUIT ?
    >>> minuit = ...
    >>> if 5    in minuit : ...
    >>> if 'p6' in minuit : ... 
    """

    if isinstance ( p, partypes ) :
        return 0 <= p < self.GetNumPars()

    if isinstance ( p ,   string_types ) :
        
        val = ctypes.c_double ( 0 )
        err = ctypes.c_double ( 0 )
        low = ctypes.c_double ( 0 )
        up  = ctypes.c_double ( 0 )
        idx = ctypes.c_int    ( 0 )

        for i in self :
            name = ROOT.TString() 
            self.mnpout ( i , name , val , err , low , up , idx )
            if 0 <= idx.value and str ( name ) == p : return True 
            
    return False 
    
ROOT.TMinuit . __contains__ = _mn_contains_ 
ROOT.TMinuit . __len__      = lambda s : s.GetNumPars() 

ROOT.TMinuit . par         = _mn_par_
ROOT.TMinuit . parameter   = _mn_par_
ROOT.TMinuit . __getitem__ = _mn_par_
ROOT.TMinuit . __call__    = _mn_par_

# =============================================================================
## get the parameter name
#  @code
#  minuit = ...
#  name = minuit.par_name ( 4 ) 
#  @endcode
def _mn_parname_ ( self , index ) :
    """Get the parameter name
    >>> minuit = ...
    >>> name = minuit.par_name ( 4 ) 
    """
    if not index in self :
        raise IndexError ( "No parameter with index %s" % index )
    
    val  = ctypes.c_double ( 0 )
    err  = ctypes.c_double ( 0 )
    low  = ctypes.c_double ( 0 )
    up   = ctypes.c_double ( 0 )
    idx  = ctypes.c_int    ( 0 )
    name = ROOT.TString() 
    self.mnpout ( index  , name , val , err , low , up , idx )
    if 0 <= idx.value : return str ( name ) 

    raise IndexError ( "No parameter with index %s" % index )
     

ROOT.TMinuit . par_name = _mn_parname_

# =============================================================================
## Get the index for parameter with the name
#  code
#  minuit = ...
#  index  = minuit.par_index ( 'p4' )
#  @endcode
def _mn_index_ ( self , name ) :
    """Get the index for parameter with the name 
    >>> minuit = ...
    >>> index  = minuit.par_index ( 'p4' )
    """
    val = ctypes.c_double ( 0 )
    err = ctypes.c_double ( 0 )
    low = ctypes.c_double ( 0 )
    up  = ctypes.c_double ( 0 )
    idx = ctypes.c_int    ( 0 )
    
    for i in self :
        
        pname = ROOT.TString() 
        self.mnpout ( i , pname , val , err , low , up , idx )
        if 0 <= idx.value and str ( pname ) == name : return i 
        
    raise IndexError ( "No parameter with name %s" % name )

ROOT.TMinuit . par_index = _mn_index_
ROOT.TMinuit . index     = _mn_index_

# =============================================================================
## iterator over TMinuit parameter indices
#  @code
#  m = ... #TMinuit object
#  for i in m : print m[i]
#  @endcode
def _mn_iter_ ( self ) :
    """Iterator for TMinuit indices:

    >>> m = ... #TMinuit object
    >>> for i in m : print m[i]
    
    """
    i = 0
    while i < len ( self )  :
        yield i
        i += 1

ROOT.TMinuit . __iter__ = _mn_iter_

# =============================================================================
## Simple wrapper for <code>ROOT.TMinuit.mnexcm</code> function
#  @see TMinuit::mnexcm
#  Execute MINUIT command
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-04-01
def _mn_exec_ ( self , command , *args ) :
    """Execute MINUIT  command
    Simple wrapper for `ROOT.TMinuit.mnexcm` function
    - see ROOT.TMinuit.mnexcm    
    """

    if args :
        
        from array import array
        arglist = array ( 'd' , [ i for i in args ]  )
        #
        ierr = ctypes.c_int ( 0 )
        #
        self.mnexcm ( command , arglist , len ( arglist ) , ierr )
        result = int ( ierr.value )
        
    else :
        
        result = self.Command ( command )
        
    if result and result in return_codes and result < 10 :
        lst = [ str ( i ) for i in args ]
        lst = ' '.join  ( lst )
        logger.warning ( "Command %s -> %s:%s" % ( command + ' ' + lst     ,
                                                   result                  ,
                                                   return_codes [ result ] ) )
        
    return result 

_mn_exec_ . __doc__  += '\n' + ROOT.TMinuit.mnexcm . __doc__

ROOT.TMinuit.execute = _mn_exec_

# =============================================================================
## excute MINUIT "SHOW" command
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-04-01
def _mn_show_ ( self , what = 'PAR' , *args ) :
    """Execute MINUIT  'SHOW'-command
    """
    return _mn_exec_ ( self , "SHOW %s" % what.upper() ,  *args )


ROOT.TMinuit.show = _mn_show_

# =============================================================================
## set the parameter  (and optionally fix it!)
#  @code
#  minuit = ...
#
#  minuit.setPar        ( 0 , 10.0       ) ## par(0) == 10.0 
#  minuit.setPar        ( 0 , 10.0, True ) ## ... and fix it!
#
#  minuit.set_par       ( 0 , 10.0       ) ## par(0) == 10.0 
#  minuit.set_par       ( 0 , 10.0, True ) ## ... and fix it!
# 
#  minuit.setParameter  ( 0 , 10.0       ) ## par(0) == 10.0 
#  minuit.setParameter  ( 0 , 10.0, True ) ## ... and fix it!
#
#  minuit.set_parameter ( 0 , 10.0       ) ## par(0) == 10.0 
#  minuit.set_parameter ( 0 , 10.0, True ) ## ... and fix it!
#
#  minuit [ 0 ] = 10                       ## ditto! 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-09-28
def _mn_set_par_ ( self , i , val , fix = False ) :
    """Set MINUIT parameter for some value and optionally fix it!

    >>> minuit = ...
    
    >>> minuit.setPar        ( 0 , 10.0       ) ## par(0) == 10.0 
    >>> minuit.setPar        ( 0 , 10.0, True ) ## ... and fix it!
    
    >>> minuit.set_par       ( 0 , 10.0       ) ## par(0) == 10.0 
    >>> minuit.set_par       ( 0 , 10.0, True ) ## ... and fix it!
    
    >>> minuit.setParameter  ( 0 , 10.0       ) ## par(0) == 10.0 
    >>> minuit.setParameter  ( 0 , 10.0, True ) ## ... and fix it!
    
    >>> minuit.set_parameter ( 0 , 10.0       ) ## par(0) == 10.0 
    >>> minuit.set_parameter ( 0 , 10.0, True ) ## ... and fix it!
    
    >>> minuit [ 0 ] = 10                       ## ditto! 
    
    """
    if not i in self :
        raise IndexError ( "Invalid parameter index %s!" % i )
    
    if isinstance ( i , string_types ) : i = _mn_index_ ( self , i )
        
    #
    if hasattr ( val , 'value' ) : val = val.value()
    #
    ierr =  _mn_exec_ ( self , "SET PAR" , i + 1 , val )
    #
    if fix : self.FixParameter ( i ) 
    #
    return ierr 

ROOT.TMinuit . setPar        = _mn_set_par_
ROOT.TMinuit . setParameter  = _mn_set_par_
ROOT.TMinuit . set_par       = _mn_set_par_
ROOT.TMinuit . set_parameter = _mn_set_par_

ROOT.TMinuit . fixPar        = lambda s,i,v: _mn_set_par_ ( s , i , v , True )
ROOT.TMinuit . fix_par       = lambda s,i,v: _mn_set_par_ ( s , i , v , True )
ROOT.TMinuit . fixParameter  = lambda s,i,v: _mn_set_par_ ( s , i , v , True )
ROOT.TMinuit . fix_parameter = lambda s,i,v: _mn_set_par_ ( s , i , v , True )
ROOT.TMinuit . fix           = lambda s,i,v: _mn_set_par_ ( s , i , v , True )


ROOT.TMinuit . __setitem__ = _mn_set_par_ 


# =============================================================================
## release the parameter
#  @code
#  mn = ... # TMinuit  obejct
#  mn.release ( 1 ) 
#  mn.rel     ( 1 ) ## ditto! 
#  mn.rel_par ( 1 ) ## ditto! 
#  mn.release ( 1 , 2 , 3  ) 
#  mn.rel     ( 1 , 2 , 3  ) ## ditto! 
#  mn.rel_par ( 1 , 2 , 3  ) ## ditto! 
#  mn.release ( 'p1' ) 
#  mn.rel     ( 'p1' ) ## ditto! 
#  mn.rel_par ( 'p1' ) ## ditto! 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-09-28
def _mn_rel_par_ ( self , *pars ) :
    """Release MINUIT parameter for some value

    >>> mn = ... # TMinuit  obejct
    >>> mn.release ( 1 ) 
    >>> mn.rel     ( 1 ) ## ditto! 
    >>> mn.rel_par ( 1 ) ## ditto!
    
    >>> mn.release ( 1 , 2 , 3 ) 
    >>> mn.rel     ( 1 , 2 , 3 ) ## ditto! 
    >>> mn.rel_par ( 1 , 2 , 3 ) ## ditto!
    
    >>> mn.release ( 'p1' ) 
    >>> mn.rel     ( 'p1' ) ## ditto! 
    >>> mn.rel_par ( 'p1' ) ## ditto!
    """
    for i in pars :
        
        if not i in self :
            raise IndexError ( "Invalid parameter index %s!" % i )
        #
        if isinstance ( i , string_types ) : i = _mn_index_ ( self , i )
        #
        _mn_exec_ ( self , "REL" , i + 1 )
        
    
ROOT.TMinuit . release = _mn_rel_par_ 
ROOT.TMinuit . rel     = _mn_rel_par_ 
ROOT.TMinuit . rel_par = _mn_rel_par_ 

# ===========================================================
## Perform the actual MINUIT minimization:
#  @code 
#  m = ... #
#  m.fit()       ## run migrade! 
#  m.migrade ()  ## ditto
#  m.fit ( method = 'MIN' )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-09-28
def _mn_min_ ( self                  ,
               maxcalls  = 5000      ,
               tolerance = 0.01       ,
               method    = 'MIGRADE' ) :
    """Perform the actual MINUIT minimization:

    >>> m = ... #
    >>> m.fit()       ## run migrade! 
    >>> m.migrade ()  ## ditto
    >>> m.fit ( method = 'MIN' ) 
    
    """
    #
    return _mn_exec_ ( self , method , maxcalls , tolerance ) 

ROOT.TMinuit . migrade  = _mn_min_
ROOT.TMinuit . migrad   = _mn_min_
ROOT.TMinuit . fit      = _mn_min_

ROOT.TMinuit . hesse    = lambda s,*a : _mn_exec_ ( s , 'HESSE'   , *a )
ROOT.TMinuit . minimize = lambda s,*a : _mn_exec_ ( s , 'MIN'     , *a )
ROOT.TMinuit . seek     = lambda s,*a : _mn_exec_ ( s , 'SEEK'    , *a )
ROOT.TMinuit . simplex  = lambda s,*a : _mn_exec_ ( s , 'SIMPLEX' , *a )

#=============================================================================
## Print MINUIT information
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-09-28
def _mn_str_ ( self , l = 3 , v = 0.0 ) :
    """Print MINUIT information:

    >>> m = ...
    >>> print m
    
    """
    #
    self.mnprin ( l , v )
    return '\n'

ROOT.TMinuit . __str__   =  _mn_str_ 
ROOT.TMinuit . __repr__  =  _mn_str_ 

# =============================================================================
## define/add parameter to TMinuit 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-09-28
def _mn_add_par_ ( self       ,
                   name       ,
                   start      , 
                   step  = -1 ,
                   low   = 0  ,
                   high  = 0  ) :
    """Define/add parameter to MUNUIT
    >>> m.addPar ( 'ququ' , 10 , 0.1 )
    
    """
    if hasattr ( start , 'value' ) : start = start . value()
    if hasattr ( step  , 'value' ) : step  = step  . value()
    ## 
    if step < 0 : 
        if low < high : step = 1.e-3 *     ( high - low ) 
        elif start    : step = 1.e-3 * abs ( start      ) 
    ##
    ipar = len          ( self )
    ##
    ierr = ctypes.c_int ( 0    ) 
    ##
    self.mnparm ( ipar , name ,  start , step , low , high , ierr )
    #
    result = int ( ierr.value )
    if result : logger.error ("Error from TMinuit.mnparm %s" % result)
    ## 
    return result 

ROOT.TMinuit . addpar           = _mn_add_par_
ROOT.TMinuit . addPar           = _mn_add_par_
ROOT.TMinuit . add_par          = _mn_add_par_
ROOT.TMinuit . add_parameter    = _mn_add_par_
ROOT.TMinuit . defpar           = _mn_add_par_
ROOT.TMinuit . def_par          = _mn_add_par_
ROOT.TMinuit . define_parameter = _mn_add_par_
ROOT.TMinuit . newpar           = _mn_add_par_
ROOT.TMinuit . newPar           = _mn_add_par_
ROOT.TMinuit . new_par          = _mn_add_par_
ROOT.TMinuit . new_parameter    = _mn_add_par_

# =============================================================================
## get MINOS errors
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-09-28 
def _mn_minerr_ ( self , i ) :
    """Get MINOS errors for parameter:

    >>> m = ...       # TMinuit object
    >>> pos , neg = m.minosErr( 0 )
    
    """
    #
    if not i in self : raise IndexError
    #
    if isinstance ( i , string_types ) : i = _mn_index_ ( self , i )
    #
    eplus  = ctypes.c_double ( 0 ) 
    eminus = ctypes.c_double ( 1 ) 
    epara  = ctypes.c_double ( 2 ) 
    gcc    = ctypes.c_double ( 3 ) 
    #
    self.mnerrs ( i , eplus , eminus , epara , gcc )
    #
    return float ( eplus.value ) , float ( eminus.value )  

ROOT.TMinuit .   minErr     = _mn_minerr_ 
ROOT.TMinuit . minosErr     = _mn_minerr_ 
ROOT.TMinuit . minos_err    = _mn_minerr_ 
ROOT.TMinuit . minos_errors = _mn_minerr_ 
ROOT.TMinuit .   minErrs    = _mn_minerr_ 
ROOT.TMinuit . minosErrs    = _mn_minerr_ 

# =============================================================================
## run MINOS
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-09-28 
def _mn_minos_ ( self , *args ) :
    """Get MINOS errors for parameter:
    
    >>> m = ...       # TMinuit object
    >>> result = m.minos ( 1 , 2  )
    
    """
    ipars  = []
    for i in args :
        if not i in self : raise IndexError
        if isinstance ( i , string_types ) : i = _mn_index_ ( self , i )
        ipars.append  ( i + 1 ) ## note + 1 here! 

    return _mn_exec_ ( self , 'MINOS' , 500 , *tuple(ipars) ) 

ROOT.TMinuit . minos = _mn_minos_
# =============================================================================

# =============================================================================
## get current Minuit statistics 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2013-04-01
def _mn_stat_  ( self ) :
    """Get current Minuit status

    >>> mn   = ... # TMinuit object
    >>> stat = mn.stat()
    
    
    Returns concerning the current status of the minimization
    *-*      =========================================================
    *-*       User-called
    *-*          Namely, it returns:
    *-*        FMIN: the best function value found so far
    *-*        FEDM: the estimated vertical distance remaining to minimum
    *-*        ERRDEF: the value of UP defining parameter uncertainties
    *-*        NPARI: the number of currently variable parameters
    *-*        NPARX: the highest (external) parameter number defined by user
    *-*        ISTAT: a status integer indicating how good is the covariance
    *-*           matrix:  0= not calculated at all
    *-*                    1= approximation only, not accurate
    *-*                    2= full matrix, but forced positive-definite
    *-*                    3= full accurate covariance matrix
    *
    
    """
    
    fmin    = ctypes.c_double ( 1 ) 
    fedm    = ctypes.c_double ( 2 )
    errdef  = ctypes.c_double ( 3 ) 
    
    npari   = ctypes.c_int    ( 1 )
    nparx   = ctypes.c_int    ( 2 )
    istat   = ctypes.c_int    ( 0 )
    #
    self . mnstat( fmin, fedm, errdef, npari , nparx , istat )
    #
    return { 'FMIN'   : float ( fmin   . value ) ,
             'FEDM'   : float ( fedm   . value ) ,
             'ERRDEF' : float ( errdef . value ) ,
             'NPARI'  : int   ( npari  . value ) ,
             'NPARX'  : int   ( nparx  . value ) ,
             'ISTAT'  : int   ( istat  . value ) } 

_mn_stat_ . __doc__  += '\n' + ROOT.TMinuit.mnstat . __doc__

ROOT.TMinuit.stat    = _mn_stat_ 
# =============================================================================
## get UP-parameter for err-def 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2013-04-01
def _mn_get_errdef_ ( self ) :
    """Get UP-parameter used to define the uncertainties

    >>> mn = ... # TMoniut object
    >>> up  = mn.GetErrorDef()
    
    """
    return _mn_stat_ ( self ) ['ERRDEF']

ROOT.TMinuit.errDef           = _mn_get_errdef_ 
ROOT.TMinuit.err_def          = _mn_get_errdef_ 
ROOT.TMinuit.error_def        = _mn_get_errdef_ 
ROOT.TMinuit.error_definition = _mn_get_errdef_ 
ROOT.TMinuit.GetErrorDef      = _mn_get_errdef_ 

# =============================================================================
## create N-sigma contour 
#  @code
#  mn    = ... # TMinuit object
#  graph = mn.contour ( 100 , 'par1' , 'par2' , 1  )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2011-06-07
def _mn_contour_ ( self , npoint , par1  , par2 , nsigma = 1 ) :
    """Create n-sigma contour for par1 vs par2

    >>> mn = ... # TMinuit object
    >>> graph = mn.contour ( 100 , 1 , 2 )
    
    """
    if     npoint < 4   : raise ValueError ( 'contour: npoint (%s) must be >= 4'  % npoint )
    if not par1 in self : raise ValueError ( 'contour: par1(%s) is not in Minuit' % par1   )
    if not par2 in self : raise ValueError ( 'contour: par2(%s) is not in Minuit' % par2   )
    
    if isinstance ( par1 , string_types ) : par1 = _mn_index_ ( self , par1 )
    if isinstance ( par2 , string_types ) : par2 = _mn_index_ ( self , par2 )
    
    if     par1 == par2 : raise ValueError ( 'contour: par1 == par2(%s) '         % par2   )
    #

    name1 = _mn_par_name_  ( self , par1 )
    name2 = _mn_par_name_  ( self , par2 )
    
    ## save old error defintion
    old_err_def = self.GetErrorDef()

    ## set new error definition
    self.SetErrorDef ( nsigma * nsigma )
    
    graph  = self.Contour ( npoint , par1 , par2 )
    
    logger.debug ( "Build 2D-contour %s vs %s at %s sigma" % ( name2 ,   name1 , nsigma ) )
    
    ## restore old error defininion
    self.SetErrorDef ( old_err_def ) 

    status = self.GetStatus()
    #
    if graph and 0 == status : return graph
    
    logger.error ( 'TMinuit::Contour: status %i' % status ) 

    return graph


_mn_contour_ . __doc__ += '\n' + ROOT.TMinuit.Contour . __doc__ 

ROOT.TMinuit . contour = _mn_contour_


# =============================================================================
## create N-sigma contours
#  @code
#  mn = ... # TMinuit object
#  graphs = mn.contours ( 100 , 'par1' , 'par2' , 1 ,  2 , 3 )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2011-06-07
def _mn_contours_ ( self , npoints , par1  , par2 , *nsigmas ) :
    """ Create n-sigma contours for par1 vs par2

    >>> mn = ... # TMinuit object
    >>> graphs = mn.contours ( 100 , 'par1' , 'par2' , 1 ,  2 , 3 )
    """
    
    if     npoints < 4  : raise ValueError ( 'contour: npoint (%s) must be >= 4'  % npoints )
    if not par1 in self : raise ValueError ( 'contour: par1(%s) is not in Minuit' % par1    )
    if not par2 in self : raise ValueError ( 'contour: par2(%s) is not in Minuit' % par2    )
    
    if isinstance ( par1 , string_types ) : par1 = _mn_index_ ( self , par1 )
    if isinstance ( par2 , string_types ) : par2 = _mn_index_ ( self , par2 )
    
    if     par1 == par2 : raise ValueError ( 'contour: par1 == par2(%s) '         % par2   )
    #
    
    graphs = []
    
    for nsigma in sigmas :

        g = _mn_contour_ ( self , npoints ,  par1  , par2 , nsigma )
        graphs.append ( g )

    return tuple ( graphs ) 
        
    
ROOT.TMinuit . contours = _mn_contours_

# =============================================================================
## get the covariance matrix from TMinuit
def _mn_cov_ ( self , size = -1 , root = False ) :
    """Get the covariance matrix from TMinuit

    >>> mn  = ... # TMinuit object
    >>> cov = mn.cov() 
    
    """
    #
    if size <= 0 : size = len ( self )
    size = min ( size , len ( self ) ) 
    #
    from array import array
    matrix = array ( 'd' , [ 0 for i in range ( 0 , size * size ) ]  )
    self.mnemat ( matrix , size )
    #
    import ostap.math.linalg
    from   ostap.core.core import Ostap 
    mtrx = Ostap.Math.SymMatrix ( size )() 
    for i in range ( 0 , size ) :
        for j in range ( i , size ) :            
            mtrx [ i , j ] = matrix [ i * size + j ]
            
    return mtrx

# =============================================================================
## get the correlation matrix from TMinuit
def _mn_cor_ ( self , size = -1 , root  = False ) :
    """ Get the correlation matrix from TMinuit

    >>> mn  = ... # TMinuit object
    >>> cor = mn.cor() 
        
    """
    #
    cov = self.cov ( size , root )
    #
    if   isinstance ( cov , ROOT.TMatrix ) :

        size  = cov.GetNrows()
        root  = True
        
    else : size = cov.kRows

    ## use ROOT matrices 
    if root : cor = ROOT.TMatrix  ( size , size )
    else    : cor = cov.__class__ () 

    for i in range ( 0 , size ) :
        
        d_i = cov ( i , i )
        cor [ i , i ] = 1 if 0 < d_i else 0
        
        for j in range ( i + 1 , size  ) :
            
            d_j = cov ( j , j )
            
            if 0 != cov ( i , j ) and 0 < d_i and 0 < d_j  :                
                cor [ i ,   j ] = cov ( i , j ) / math.sqrt ( d_i * d_j )                
            else :
                cor [ i ,   j ] = 0

    return cor
            
_mn_cor_ . __doc__ += '\n' + ROOT.TMinuit.mnemat . __doc__ 


ROOT.TMinuit . cov  = _mn_cov_
ROOT.TMinuit . cor  = _mn_cor_
ROOT.TMinuit . corr = _mn_cor_


# =============================================================================
## Build parameter table from (T)Minuit
def _mn_table_ ( self , title =  '' , prefix = '' ) :
    """ Build parameter table from (T)Minuit
    """
    
    header = ( 'Parameter' , '' , 'Value' )
    rows   = [] 

    from ostap.fitting.utils import fit_status  , cov_qual
    from ostap.logger.pretty import pretty_float, pretty_ve, pretty_vae 

    status = self.GetStatus()
    if status :
        status = fit_status ( status ) 
        row = '  Status' , '' , status 
        rows.append ( row )

    stat   = _mn_stat_ (  self )

    istat  = stat.pop ( 'ISTAT'  , None )
    if not istat is None : 
        cq = ''
        if   -1 == istat             : cq =              cov_qual ( istat )
        elif  3 == istat             : cq = allright   ( cov_qual ( istat ) )
        elif  istat in ( 0 , 1 , 2 ) : cq = attentiont ( cov_qual ( istat ) )
        else                         : cq =              cov_qual ( istat )
        row = 'Covariance matrix quality' , '' , cq 
        rows.append  ( row ) 
    
    fmin   = stat.pop ( 'FMIN'   , None )
    if not fmin is None :
        s , n = pretty_float ( fmin ) 
        if n : n = '[10^%+d]' % n
        else : n = '' 
        row = 'Minimized FCN value' ,  n , s  
        rows.append ( row ) 
    
    fedm   = stat.pop ( 'FEDM'   , None )
    if not fedm is None :
        s , n = pretty_float ( fedm ) 
        if n : n = '[10^%+d]' % n
        else : n = '' 
        row = 'Estimated distance to minimum' ,  n , s 
        rows.append ( row ) 
    
    errdef = stat.pop ( 'ERRDEF' , None ) 
    ## needed ? 

    has_limits = False
    has_minos  = False

    val  = ctypes.c_double ( 0 )
    err  = ctypes.c_double ( 0 )
    low  = ctypes.c_double ( 0 ) 
    high = ctypes.c_double ( 0 ) 
    idx  = ctypes.c_int    ( 0 )
    
    dct_pars  = {} 
    ## loop over all parameters 
    for i in  self :
        
        name = ROOT.TString() 
        self.mnpout ( i , name , val , err , low , high , idx )
        if not 0 <= idx.value : continue
        
        dct = {}
        
        dct [ 'name' ] = '#%-2d: %s' % ( i , str ( name ).strip() ) 
        dct [ 'value'] = val.value
        
        if low.value < high.value :
            dct [ 'low' ] = low .value
            dct [ 'high'] = high.value
            has_limits = True
            
        if 0 <= err.value :
            dct [ 'error' ] = err.value
            mn_plus , mn_minus = _mn_minerr_ ( self , str ( name ) )
            if ( 0 <= mn_plus and mn_minus <= 0 ) :
                error = err.value 
                if abs ( mn_minus ) != mn_plus or abs ( mn_minus ) != error or mn_plus != error : 
                    dct [ 'minos+' ] = mn_plus 
                    dct [ 'minos-' ] = mn_minus 
                    has_minos = True
                
        dct_pars[ i ] = dct
    
    if has_minos :
        ## some parameters have MINOS errors, add columns
        header = ( header ) + ( 'neg-minos' , 'pos-minos' ) 
        rows   = [ r + ( '','' ) for r in rows ] 
        
    if has_limits :
        ## some parameters have LIMITS, add columns 
        header = ( header ) + ( 'low limit' , 'high limit' ) 
        rows   = [ r + ( '','' ) for r in rows ]         

    for p in dct_pars :
        
        pdict  = dct_pars [ p ]

        row = []
        
        row.append ( pdict.pop ( 'name' )      )
        
        value    = pdict.pop ( 'value'         )
        error    = pdict.pop ( 'error'  , None ) 
        mn_plus  = pdict.pop ( 'minos+' , None )
        mn_minus = pdict.pop ( 'minos-' , None )

        if   ( not mn_plus is None and not mn_minus is None ) :
            s , n = pretty_vae ( VAE ( value , mn_minus , mn_plus ) , parentheses = False ) 
        elif not error is None :
            s , n = pretty_ve  ( VE  ( value , error    * error   )  , parentheses = False ) 
        else : 
            s , n = pretty_float ( value )
            s = "%s(fixed)" % s
            
        if n : row.append ( '[10^%+d]' % n )
        else : row.append ( ''             )
            
        row.append ( s )
                                 
        if False and has_minos :
            mn_plus  = pdict.pop ( 'minos+' , None ) 
            mn_minus = pdict.pop ( 'minos-' , None )
            
            if mn_plus  is None : mn_plus  = ''
            else                : mn_plus  = '%8f' % ( mn_plus   * 10** n )
            if mn_minus is None : mn_minus = ''
            else                : mn_minus = '%8f' % ( mn_minus  * 10** n )

            row.append ( mn_minus )  
            row.append ( mn_plus  )  

        if has_limits :
            low  = pdict.pop ( 'low'  , None ) 
            high = pdict.pop ( 'high' , None )
            if low  is None : low  = ''
            else            : low  = '%8f' %  ( low  * 10 ** n ) 
            if high is None : high = ''
            else            : high = '%8f' %  ( high * 10 ** n ) 
            row.append ( low ) 
            row.append ( high )
            
        row = tuple ( row )
        rows.append ( row ) 
        
    rows = [ header ] + rows
    
    import ostap.logger.table as T 
    if rows : rows = T.remove_empty_columns ( rows )
    return T.table ( rows ,  title = title  , prefix = prefix ) 


ROOT.TMinuit .table = _mn_table_
              
_decorates_classes_ = (
    ROOT.TMinuit ,
)

_new_methdos_ = (
    ## 
    ROOT.TMinuit . __contains__ , 
    ROOT.TMinuit . __len__      , 
    ## 
    ROOT.TMinuit . par          , 
    ROOT.TMinuit . parameter    , 
    ROOT.TMinuit . __getitem__  , 
    ##     
    ROOT.TMinuit . par_name    , 
    ROOT.TMinuit . par_index   , 
    ROOT.TMinuit . index       , 
    ## 
    ROOT.TMinuit . __iter__  , 
    ## 
    ROOT.TMinuit.execute  , 
    ROOT.TMinuit.show     , 
    ## 
    ROOT.TMinuit . setPar        , 
    ROOT.TMinuit . setParameter  , 
    ROOT.TMinuit . set_par       , 
    ROOT.TMinuit . set_parameter , 
    ## 
    ROOT.TMinuit . fixPar        , 
    ROOT.TMinuit . fix_par       , 
    ROOT.TMinuit . fixParameter  , 
    ROOT.TMinuit . fix_parameter , 
    ROOT.TMinuit . fix           , 
    ## 
    ROOT.TMinuit . __setitem__   , 
    ##         
    ROOT.TMinuit . release  , 
    ROOT.TMinuit . rel      , 
    ROOT.TMinuit . rel_par  , 
    ##     
    ROOT.TMinuit . migrade  , 
    ROOT.TMinuit . migrad   , 
    ROOT.TMinuit . fit      , 
    ## 
    ROOT.TMinuit . hesse    , 
    ROOT.TMinuit . minimize , 
    ROOT.TMinuit . seek     , 
    ROOT.TMinuit . simplex  , 
    ## 
    ROOT.TMinuit . __str__          , 
    ROOT.TMinuit . __repr__         , 
    ## 
    ROOT.TMinuit . addpar           , 
    ROOT.TMinuit . addPar           , 
    ROOT.TMinuit . add_par          , 
    ROOT.TMinuit . add_parameter    , 
    ROOT.TMinuit . defpar           , 
    ROOT.TMinuit . def_par          , 
    ROOT.TMinuit . define_parameter , 
    ROOT.TMinuit . newpar           , 
    ROOT.TMinuit . newPar           , 
    ROOT.TMinuit . new_par          , 
    ROOT.TMinuit . new_parameter    , 
    ##     
    ROOT.TMinuit .   minErr         , 
    ROOT.TMinuit . minosErr         , 
    ROOT.TMinuit . minos_err        , 
    ROOT.TMinuit . minos_errors     , 
    ROOT.TMinuit .   minErrs        , 
    ROOT.TMinuit . minosErrs        , 
    ##     
    ROOT.TMinuit . minos            , 
    ## 
    ROOT.TMinuit.errDef             ,
    ROOT.TMinuit.err_def            ,
    ROOT.TMinuit.error_def          ,
    ROOT.TMinuit.error_definition   ,
    ROOT.TMinuit.GetErrorDef        ,
    ## 
    ROOT.TMinuit . contours         , 
    ROOT.TMinuit . contour          , 
    ## 
    ROOT.TMinuit . cov              , 
    ROOT.TMinuit . cor              , 
    ROOT.TMinuit . corr             ,  
    ROOT.TMinuit .table             ,
    )

    
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
