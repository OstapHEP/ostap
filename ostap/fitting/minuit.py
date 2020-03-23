#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/minuit.py
#  Module with decoration of some (T)Minuit functions for efficient use in python
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
#
# =============================================================================
"""Decoration for some (T)Minuit functions for efficient use in python"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = () 
# =============================================================================
import ROOT, ctypes 
from   sys                    import version_info as python_version
from   ostap.core.core        import cpp, VE
from   ostap.core.ostap_types import integer_types, string_types 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.minuit' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Some useful decorations for (T)Minuit functions')
# =============================================================================
partypes = integer_types
# =============================================================================
## get the parameter from Minuit 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-09-28
def _mn_par_ ( self , i ) :
    """Get the parameter from minuit

    >>> mn = ...             # TMinuit object
    >>> p1 = mn[0]           # get the parameter 
    >>> p1 = mn.par(0)       # ditto 
    >>> p1 = mn.parameter(0) # ditto 
    """
    if not i in self : raise IndexError
    #
    ## val = ROOT.Double ( 0 )
    ## err = ROOT.Double ( 0 )
    ##
    val = ctypes.c_double ( 0 )  
    err = ctypes.c_double ( 0 ) 
    #
    res = self.GetParameter ( i , val , err )
    #
    val = float ( val.value )
    err = float ( err.value )
    #
    return VE ( val , err*err )

# =============================================================================
def _mn_contains_ (  self , p ) :
    return  isinstance ( p, partypes )  and 0 <= p < self.GetNumPars()
    
ROOT.TMinuit . __contains__ = _mn_contains_ 
ROOT.TMinuit . __len__      = lambda s : s.GetNumPars() 

ROOT.TMinuit . par         = _mn_par_
ROOT.TMinuit . parameter   = _mn_par_
ROOT.TMinuit . __getitem__ = _mn_par_
ROOT.TMinuit . __call__    = _mn_par_

# =============================================================================
## iterator over TMinuit indices 
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
## @see TMinuit::mnexcm
## excute MINUIT command
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-04-01
def _mn_exec_ ( self , command , *args ) :
    """Execute MINUIT  command
    Simple wrapper for ROOT.TMinuit.mnexcm function
    - see ROOT.TMinuit.mnexcm    
    """
    if not  args : args = 0 ,
    
    from array import array
    arglist = array ( 'd' , [ i for i in args ]  )
    #
    ierr = ctypes.c_int ( 0 ) 
    ##
    self.mnexcm ( command , arglist , len ( arglist ) , ierr )
    #
    return int ( ierr.value )

_mn_exec_ . __doc__  += '\n' + ROOT.TMinuit.mnexcm . __doc__

ROOT.TMinuit.execute = _mn_exec_

# =============================================================================
## excute MINUIT "SHOW" command
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-04-01
def _mn_show_ ( self , what = 'PAR' , *args ) :
    """Execute MINUIT  command
    """
    if not args : args = [ 0 ]
    ##
    what = what.upper()
    whar = what.replace ( 'SHOW',' ' )
    return _mn_exec_ ( self , 'SHOW ' + what , *args )

ROOT.TMinuit.show = _mn_show_

# =============================================================================
## set the parameter 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-09-28
def _mn_set_par_ ( self , i , val , fix = False ) :
    """Set MINUIT parameter for some value
    """
    if not i in self : raise IndexError
    #
    if hasattr ( val , 'value' ) : val = val.value()
    #
    ierr =  _mn_exec_ ( self , "SET PAR" , i + 1 , val )
    #
    if fix : self.FixParameter ( i ) 
    #
    return ierr 

ROOT.TMinuit . setPar       = _mn_set_par_
ROOT.TMinuit . setParameter = _mn_set_par_

ROOT.TMinuit . fixPar       = lambda s,i,v: _mn_set_par_ ( s , i , v , True )
ROOT.TMinuit . fixParameter = lambda s,i,v: _mn_set_par_ ( s , i , v , True )
ROOT.TMinuit . fix          = lambda s,i,v: _mn_set_par_ ( s , i , v , True )


ROOT.TMinuit . __setitem__ = _mn_set_par_ 


# =============================================================================
## release the parameter 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-09-28
def _mn_rel_par_ ( self , i ) :
    """Release MINUIT parameter for some value

    >>> mn = ... # TMinuit  obejct
    >>> mn.release ( 1 ) 
    """
    if not i in self : raise IndexError
    #
    return _mn_exec_ ( self , "REL" , i + 1 )
    #

ROOT.TMinuit . release = _mn_rel_par_ 

# ===========================================================
## set the parameter 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-09-28
def _mn_min_ ( self                  ,
               maxcalls  = 5000      ,
               tolerance = 0.1       ,
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

# =============================================================================
## set the parameter 
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

ROOT.TMinuit . Print     =  _mn_str_ 
ROOT.TMinuit . __str__   =  _mn_str_ 
ROOT.TMinuit . __repr__  =  _mn_str_ 

# =============================================================================
## define/add parameter to TMinuit 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-09-28
def _mn_add_par_ ( self    , name      ,
                   start   , step = -1 ,
                   low = 0 , high = 0  ) :
    """Define/add parameter to MUNUIT
    >>> m.addPar ( 'ququ' , 10 , 0.1 )
    
    """
    if hasattr ( start , 'value' ) : start = start . value()
    if hasattr ( step  , 'value' ) : step  = step  . value()
    ## 
    if step < 0 : step = abs ( 0.01 * start ) 
    ##
    from array import array
    starts  = array ( 'd' , 1 * [ start ] )
    steps   = array ( 'd' , 1 * [ step  ] )
    #
    ipar    = len         ( self )
    ##
    ierr = ctypes.c_int ( 0 ) 
    ##
    self.mnparm ( ipar , name ,  start , step , low , high , ierr )
    #
    return int ( ierr.value )

ROOT.TMinuit . addpar = _mn_add_par_
ROOT.TMinuit . addPar = _mn_add_par_
ROOT.TMinuit . defpar = _mn_add_par_
ROOT.TMinuit . defPar = _mn_add_par_
ROOT.TMinuit . newpar = _mn_add_par_
ROOT.TMinuit . newPar = _mn_add_par_

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
    eplus  = ROOT.Double ( 0 ) 
    eminus = ROOT.Double ( 0 ) 
    epara  = ROOT.Double ( 0 ) 
    gcc    = ROOT.Double ( 0 ) 
    #
    self.mnerrs ( i , eplus , eminus , epara , gcc )
    #
    return eplus,eminus 

ROOT.TMinuit .   minErr  = _mn_minerr_ 
ROOT.TMinuit . minosErr  = _mn_minerr_ 
ROOT.TMinuit .   minErrs = _mn_minerr_ 
ROOT.TMinuit . minosErrs = _mn_minerr_ 

# =============================================================================
## run MINOS
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2012-09-28 
def _mn_minos_ ( self , *args ) :
    """Get MINOS errors for parameter:
    
    >>> m = ...       # TMinuit object
    >>> result = m.minos( 1 , 2  )
    
    """
    ipars  = []
    for i in args :
        if not i in self : raise IndexError
        ipars.append ( i )

    return _mn_exec_ ( self , 'MINOS' , 200 , *tuple(ipars) ) 

ROOT.TMinuit . minos = _mn_minos_
# =============================================================================

# =============================================================================
## get current Minuit statistics 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2013-04-01
def _mn_stat_  ( self ) :
    """Get current Minuit status

    >>> mn   = ... # TMoniut object
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
    # fmin    = ROOT.Double ( )
    # fedm    = ROOT.Double ( )
    # errdef  = ROOT.Double ( )
    
    fmin    = ctypes.c_double () 
    fedm    = ctypes.c_double ()
    errdef  = ctypes.c_double () 
    
    npari   = ctypes.c_int ( 1 )
    nparx   = ctypes.c_int ( 2 )
    istat   = ctypes.c_int ( 0 )
    #
    self . mnstat( fmin, fedm, errdef, npari , nparx , istat )
    #
    return { 'FMIN'   : float( fmin   . value ) ,
             'FEDM'   : float( fmin   . value ) ,
             'ERRDEF' : float( errdef . value ) ,
             'NPARI'  : int  ( npari  . value ) ,
             'NPARX'  : int  ( nparx  . value ) ,
             'ISTAT'  : int  ( istat  . value ) } 

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

ROOT.TMinuit.errDef      = _mn_get_errdef_ 
ROOT.TMinuit.GetErrorDef = _mn_get_errdef_ 

# =============================================================================
## create N-sigma contour 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2011-06-07
def _mn_contour_ ( self , npoint , par1  , par2 , nsigma = 1 ) :
    """Create n-sigma contour for par1 vs par2

    >>> mn = ... # TMinuit object
    >>> graph = mn.contour( 100 , 1 , 2 )
    
    """
    if     npoint < 4   : raise ValueError ( 'contour: npoint (%s) must be >= 4'  % npoint )
    if not par1 in self : raise ValueError ( 'contour: par1(%s) is not in Minuit' % par1   )
    if not par2 in self : raise ValueError ( 'contour: par2(%s) is not in Minuit' % par2   )
    if     par1 == par2 : raise ValueError ( 'contour: par1 == par2(%s) '         % par2   )
    #
    ## save old error defintion
    #
    old_err_def = self.GetErrorDef()
    #
    ## set new error definition
    #
    self.SetErrorDef ( nsigma * nsigma )
    
    graph  = self.Contour ( npoint , par1 , par2 )

    #
    ## restore old error defininion
    #
    status = self.GetStatus()
    self.SetErrorDef ( old_err_def ) 
    #
    if graph and 0 == status : return graph
    logger.error ( 'TMinuit::Contour: status %i' % status ) 
    return graph 

_mn_contour_ . __doc__ += '\n' + ROOT.TMinuit.Contour . __doc__ 

ROOT.TMinuit . contour = _mn_contour_


_rv = ROOT.gROOT.GetVersionInt() // 10000
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
    matrix = array ( 'd' , [ 0 for i in range(0, size * size) ]  )
    self.mnemat ( matrix , size )
    #
    from ostap.math.linalg import Ostap 
    mtrx = Ostap.Math.SymMatrix ( size )() 
    for i in range ( 0 , size ) :
        for j in range ( i , size ) :            
            mtrx [ i , j ] = matrix [ i * size + j ]
            
    return mtrx

# =============================================================================
## get the correlation matrix from TMinuit
def _mn_cor_ ( self , size = -1 , root  = False ) :
    """Get the correlation matrix from TMinuit

    >>> mn  = ... # TMinuit object
    >>> cor = mn.cor() 
        
    """
    #
    cov = self.cov ( size , root )
    #
    from math import sqrt
    #
    if   isinstance ( cov , ROOT.TMatrix ) :

        size  = cov.GetNrows()
        root  = True
        
    else : size = cov.kRows

    ## use ROOT matrices 
    if root : cor = ROOT.TMatrix  ( size , size )
    else    : cor = cov.__class__ () 

    for i in range(0, size ) :
        
        d_i = cov ( i , i )
        cor [ i , i ] = 1 if 0 < d_i else 0
        
        for j in range ( i + 1 , size  ) :
            
            d_j = cov ( j , j )
            
            if 0 != cov ( i , j ) and 0 < d_i and 0 < d_j  :
                
                if root and _rv < 6  : cor [ i ] [ j ] = cov ( i , j ) / sqrt ( d_i * d_j )
                else                 : cor [ i ,   j ] = cov ( i , j ) / sqrt ( d_i * d_j )
                
            else :
                
                if root and _rv < 6  : cor [ i ] [ j ] = 0 
                else                 : cor [ i ,   j ] = 0

    return cor
            
_mn_cor_ . __doc__ += '\n' + ROOT.TMinuit.mnemat . __doc__ 


ROOT.TMinuit . cov  = _mn_cov_
ROOT.TMinuit . cor  = _mn_cor_
ROOT.TMinuit . corr = _mn_cor_

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
