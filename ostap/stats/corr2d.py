#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/corr2d.py 
#  Simple 2D-decorrelation transformation 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-08
# =============================================================================
""" Simple 2D-decorrelation transformation """
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-06-08"
__all__     = ( 'Corr2D', )
# =============================================================================
from   ostap.core.core    import cpp , WSE, Ostap 
from   ostap.trees.cuts   import expression_types 
import ostap.math.linalg
import ROOT,math
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.corr2D' )
else                       : logger = getLogger ( __name__ )
# =============================================================================
## get error-function 
if  not hasattr ( math , 'erf' ) :
    from ostap.math.math_ve import erf as _erf 
    math.erf = _erf
    logger.debug ( 'ostap.math.erf       is added to math' )
# =============================================================================
## error function
#  @see http://en.wikipedia.org/wiki/Error_function
erf = math.erf
# ============================================================================
## get Gaussian cdf 
if not hasattr ( math , 'gauss_cdf' ) :
    from ostap.math.math_ve import gauss_cdf as _gauss_cdf
    math.gauss_cdf = _gauss_cdf
    logger.debug ( 'ostap.math.gauss_cdf is added to math' )
# =============================================================================
gauss_cdf = math.gauss_cdf
# =============================================================================
## invoke decorations of  matrices & vectors 
C2 = Ostap.Math.SymMatrix ( 2 )
V2 = Ostap.Math.Vector    ( 2 )
# =============================================================================
## Simple 2D-decorrelation transformation 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-08
class Corr2D(object) :
    """ Simple 2D-decorrelation transformation 
    """
    # =========================================================================
    ## constructor
    #  @param chain     (INPUT) tree/chain or dataset 
    #  @param var1      (INPUT) expression for the first variable 
    #  @param var2      (INPUT) expression for the second variable
    #  @param selection (INPUT) cuts, if needed
    #  @param first     (INPUT) the first event
    #  @param last      (INPUT) the last event     
    def __init__  ( self           ,
                    dataset        ,
                    exprs1         ,
                    exprs2         ,
                    selection = '' , *args )  :

        ## extra arguments 
        self.__args = args
        
        assert exprs1 and isinstance ( exprs1 , expression_types ) , \
            "Invalid type of `exprs1`:%s" % type ( exprs1 )
        assert exprs2 and isinstance ( exprs2 , expression_types ) , \
            "Invalid type of `exprs2`:%s" % type ( exprs2 )
        assert not selection or isinstance ( selection , expression_types ) , \
            "Invalid type of `selection`: %s" % type ( selection )
        
        self.__var1      = str ( exprs1    ).strip () 
        self.__var2      = str ( exprs2    ).strip () 
        self.__selection = str ( selection ).strip () if selection else ''

        ## get the statistics&covariances 
        self.__wcov = Ostap.StatVar.statCov ( dataset        ,
                                              self.var1      ,
                                              self.var2      ,
                                              self.selection ,
                                              *self.args     ) 
        
        self.__counter1 = self.__wcov.counter1 ()  
        self.__counter2 = self.__wcov.counter2 ()

        ## the covariance matrix 
        self.__cov2     = Ostap.Math.covariance  (  self.__wcov )
        
        ## the correlation matrix 
        self.__corr     = Ostap.Math.correlation (  self.__wcov )

        ## eigensystem: eigenvalues and transformation 
        self.__evalues, self.__evectors = self.__cov2.eigenVectors( sorted = True , ascending = False )

        ## the first and second eigenvectors
        self.__evector1 = self.__evectors.column ( 0 ) 
        self.__evector2 = self.__evectors.column ( 1 ) 
        
        ## get the first an second rows of transformation matrix
        self.__row1 = self.__evectors.row ( 0 ) 
        self.__row2 = self.__evectors.row ( 1 )
        
        ## global bias 
        self.__delta = V2()
        self.__delta [ 0 ] = float ( self.__counter1.mean() ) 
        self.__delta [ 1 ] = float ( self.__counter2.mean() ) 

        nmean1 = -1 * self.__delta [ 0 ] 
        nmean2 = -1 * self.__delta [ 1 ] 
        
        dv1 = '((%s)%+.16g)' % ( self.var1 , nmean1 )
        dv2 = '((%s)%+.16g)' % ( self.var2 , nmean2 )

        self.__decorrelated1 = "(%+.16g*%s%+.16g*%s)" % ( self.__row1 [ 0 ] , dv1 , self.__row1 [ 1 ] , dv2 )
        self.__decorrelated2 = "(%+.16g*%s%+.16g*%s)" % ( self.__row2 [ 0 ] , dv1 , self.__row2 [ 1 ] , dv2 )
        
        scale1 = math.sqrt ( self.__evalues [ 0 ] )
        scale2 = math.sqrt ( self.__evalues [ 1 ] )

        ## get the first and second rows of transformation matrix and scale them 
        self.__srow1 = self.__row1 / scale1 
        self.__srow2 = self.__row2 / scale2 

        self.__normalized1 = "(%+.16g*%s%+.16g*%s)" % ( self.__srow1 [ 0 ] , dv1 , self.__srow1 [ 1 ] , dv2 )
        self.__normalized2 = "(%+.16g*%s%+.16g*%s)" % ( self.__srow2 [ 0 ] , dv1 , self.__srow2 [ 1 ] , dv2 )
        
        ## get the first and second rows of transformation matrix and scale them
        sqrt2 = math.sqrt ( 2.0 )
        self.__ssrow1 = self.__srow1  / sqrt2 
        self.__ssrow2 = self.__srow2  / sqrt2 
        
        self.__u1 = "(%+.16g*%s%+.16g*%s)" % ( self.__ssrow1 [ 0 ] , dv1 , self.__ssrow1 [ 1 ] , dv2 )
        self.__u2 = "(%+.16g*%s%+.16g*%s)" % ( self.__ssrow2 [ 0 ] , dv1 , self.__ssrow2 [ 1 ] , dv2 )
        
        self.__uniform1 = "(0.5+0.5*TMath::Erf(%s))" % self.__u1 
        self.__uniform2 = "(0.5+0.5*TMath::Erf(%s))" % self.__u2


    @property
    def var1 ( self ) :
        """`var1' : the first expression, same as `expr1`"""
        return self.__var1

    @property
    def var2 ( self ) :
        """`var2' : the second expression, same as `expr2`"""
        return self.__var2

    @property
    def expr1 ( self ) :
        """`expr1' : the first expression, same as `var1`"""
        return self.var1
    
    @property
    def expr2 ( self ) :
        """`expr2' : the second expression, same as `var2`"""
        return self.var2

    @property
    def selection ( self ) :
        """`selection` : selection/weight  expression (same as `cuts` or `weight`)"""
        return self.__selection

    @property
    def cuts ( self ) :
        """`cuts` : selection/weight  expression (same as `selection` or `weight`)"""
        return self.selection

    @property
    def weight ( self ) :
        """`weight` : selection/weight  expression (same as `selection` or `cuts`)"""
        return self.selection

    @property
    def args  ( self ) :
        """`args` : extra arguments for `Ostap.StatVar.statCov` call"""
        return self.__args
    
    @property
    def covariance ( self ) :
        """`covariance' :  `Ostap::Math::(W)Covariance` object"""
        return self.__wcov

    @property
    def counter1 ( self ) : 
        """`conter1' : counetr for the first expression"""
        return self.covariance.counter1()

    @property
    def counter2 ( self ) : 
        """`conter2' : counetr for the second expression"""
        return self.covariance.counter2()

    @property
    def correlation ( self ) :
        """`correlation` : get the correlation coefficient """
        return self.covariance.correlation()
        
    @property
    def cov2 ( self ) :
        """`cov2` : covariance matrix"""
        return self.__cov2 

    @property
    def corr ( self ) :
        """`corr` : correlation matrix"""
        return self.__corr 

    @property
    def eigenvalues  ( self ) :
        """`eigenvalues` : eigenvalues of covariance matrix"""
        return self.__evalues

    @property
    def eigenvectors ( self ) :
        """`eigenvectors` : matrix of eigenvectors: each coluimn is eigenvector, each eigenvector is normalized to 1"""
        return self.__evectors 

    @property
    def decorrelated1 ( self ) :
        """`decorrelated1` : the first decorrelated variable  """
        return self.__decorrelated1
        
    @property
    def decorrelated2 ( self ) :
        """`decorrelated2` : the second decorrelated variable  """
        return self.__decorrelated2
    
    @property
    def normalized1 ( self ) :
        """`normalized1` : the first decorrelated&normalized variable  """
        return self.__normalized1

    @property
    def normalized2 ( self ) :
        """`normalized2` : the second decorrelated&normalized variable  """
        return self.__normalized2

    @property
    def uniform1 ( self ) :
        """`uniform1` : the first decorrelated/normalized&uniform variable  """
        return self.__uniform1
    
    @property
    def uniform2 ( self ) :
        """`uniform2` : the second decorrelated/normalized&uniform variable  """
        return self.__uniform2

    # =========================================================================
    ## Make transformation from original to decorrelated variables
    #  @code
    #  x  , y  = ...
    #  xd , yd = corr2d.decorrelated ( x , y ) 
    #  @endcode 
    def decorrelated  ( self , x , y ) :
        """ Make a transformation from original to decorrelated variables
        >>> x  , y  = ...
        >>> xd , yd = corr2d.decorrelated ( x , y ) 
        """
        dv = V2 ( x , y ) - self.__delta
        ## 
        xd = self.__row1 * dv
        yd = self.__row2 * dv
        ## 
        return xd , yd 

    # =========================================================================
    ## Make transformation from original to decorrelated&normalized variables
    #  @code
    #  x  , y  = ...
    #  xn , yn = corr2d.normalized ( x , y )
    #  @endcode 
    def normalized  ( self , x , y ) :
        """ Make a transformation from original to decorrelated variables
        >>> x  , y  = ...
        >>> xn , yn = corr2d.normalized ( x , y )
        """
        dv = V2 ( x , y ) - self.__delta
        ## 
        xn = self.__srow1 * dv
        yn = self.__srow2 * dv
        ## 
        return xn , yn
    
    # =========================================================================
    ## Make transformation from original to decorrelated,normalized&uniform variables
    #  @code
    #  x  , y  = ...
    #  xn , yn = corr2d.uniform ( x , y )
    #  @endcode 
    def uniform ( self , x , y ) :
        """ Make a transformation from original to decorrelated variables
        >>> x  , y  = ...
        >>> xn , yn = corr2d.uniform ( x , y (   
        """
        dv = V2 ( x , y ) - self.__delta
        ## 
        xu = gauss_cdf ( self.__srow1 * dv ) 
        yu = gauss_cdf ( self.__srow2 * dv ) 
        ## 
        return xu , yu
    
# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                     The END 
# =============================================================================
