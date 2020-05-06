#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/corr2d.py 
#  Simple 2D-decorrelation transformation 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-08
# =============================================================================
"""Simple 2D-decorrelation transformation """
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2014-06-08"
__all__     = ( 'Corr2D', )
# =============================================================================
import ROOT,math
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'Ostap.Corr2D' )
else                       : logger = getLogger( __name__ )
from   ostap.core.core import cpp , WSE
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
## get Gaussian cdf 
if not hasattr ( math , 'gauss_cdf' ) :
    from ostap.math.math_ve import gauss_cdf as _gauss_cdf
    math.gauss_cdf = _gauss_cdf
    logger.debug ( 'ostap.math.gauss_cdf is added to math' )
# =============================================================================
gauss_cdf = math.gauss_cdf
# =============================================================================
## simple 2D-decorrelation transformation 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-06-08
class Corr2D(object) :
    """Simple 2D-decorrelation transformation 
    """
    ## constructor
    #  @param chain     (INPUT) tree/chain or dataset 
    #  @param var1      (INPUT) expression for the first variable 
    #  @param var2      (INPUT) expression for the second variable
    #  @param selection (INPUT) cuts, if needed
    #  @param first     (INPUT) the first event
    #  @param last      (INPUT) the last event     
    def __init__  ( self              ,
                    chain             ,
                    var1              ,
                    var2              , 
                    selection = None  ,
                    first     = 0     ,
                    last      = 2**63 ) :

        import ostap.math.linalg
        
        SV    =  Ostap.StatVar
        ##
        self.stat1 = WSE () 
        self.stat2 = WSE ()
        self.cov2  = Ostap.SymMatrix(2) ()
        ##
        self.var1 = var1
        self.var2 = var2
        
        ## 
        if isinstance ( chain , ROOT.RooDataSet ) :
            chain = chain.store().tree() 
        ## 
        if selection :
            self.num   = SV.statCov ( chain      ,
                                      var1       ,
                                      var2       ,
                                      selection  ,
                                      self.stat1 ,
                                      self.stat2 ,                                 
                                      self.cov2  ,
                                      first      , last ) 
        else :
            self.num   = SV.statCov ( chain      ,
                                      var1       ,
                                      var2       ,
                                      self.stat1 ,
                                      self.stat2 ,                                 
                                      self.cov2  ,
                                      first      , last ) 


        self.corr = Ostap.SymMatrix(2)()

        rms1 = self.stat1.rms()
        rms2 = self.stat2.rms()
        
        self.corr[0,0] = self.cov2(0,0) / ( rms1 * rms1 )
        self.corr[0,1] = self.cov2(0,1) / ( rms1 * rms2 )
        self.corr[1,1] = self.cov2(1,1) / ( rms2 * rms2 )
        
        self.correlation = self.corr(0,1)
        
        ## get eigen vectors 
        self.e = self.cov2.eigenVectors()
        
        self.vct0  = self.e[1][0]
        self.vct1  = self.e[1][1]

        import math

        ## normalized eigen vectors 
        self.vct0 /= math.sqrt ( self.e[0][0] ) 
        self.vct1 /= math.sqrt ( self.e[0][1] )
        
        self.m1 = -1*self.stat1.mean().value()
        self.m2 = -1*self.stat2.mean().value()
        
        self.nvar1 = "(%+g*((%s)%+g)%+g*((%s)%+g))" % ( self.vct0[0] , self.var1 , self.m1 , self.vct0[1] , self.var2 , self.m2 ) 
        self.nvar2 = "(%+g*((%s)%+g)%+g*((%s)%+g))" % ( self.vct1[0] , self.var1 , self.m1 , self.vct1[1] , self.var2 , self.m2 ) 

        logger.info ( 'Correlation %.3f%%' % ( 100 * self.correlation )  )
        
        logger.info ( 'The 1st decorrelated variable:\n %s ' % self.nvar1 )
        logger.info ( 'The 2nd decorrelated variable:\n %s ' % self.nvar2 )

        ## normalize eigenvectors for 1/sqrt(2), just for convinency of erf. 
        sqr2i         = math.sqrt(0.5) 
        self.nvct0    = Ostap.Vector2() 
        self.nvct1    = Ostap.Vector2()
        
        self.nvct0[0] = self.vct0[0]*sqr2i
        self.nvct0[1] = self.vct0[1]*sqr2i
        self.nvct1[0] = self.vct1[0]*sqr2i
        self.nvct1[1] = self.vct1[1]*sqr2i

        self.qvar1 = "0.5+0.5*TMath::Erf( %+g*((%s)%+g) %+g*((%s)%+g) )" % ( self.nvct0[0] , self.var1 , self.m1 , self.nvct0[1] , self.var2 , self.m2 ) 
        self.qvar2 = "0.5+0.5*TMath::Erf( %+g*((%s)%+g) %+g*((%s)%+g) )" % ( self.nvct1[0] , self.var1 , self.m1 , self.nvct1[1] , self.var2 , self.m2 ) 

        logger.info ( 'The 1st decorrelated normalized variable:\n %s ' % self.qvar1 )
        logger.info ( 'The 2nd decorrelated normalized variable:\n %s ' % self.qvar2 )

    
    def   fvar1   ( self ) :
        #
        v1 = str ( self.qvar1 )
        v1 = v1.replace ( '0.5+0.5*TMath::Erf' , 'gauss_cdf' )
        v1 = v1.replace ( self.var1    ,  's.' + self.var1   )
        v1 = v1.replace ( self.var2    ,  's.' + self.var2   )
        #
        return eval ( 'lambda s: ' + v1 )
    
    def   fvar2   ( self ) :
        #
        v2 = str ( self.qvar2 )
        v2 = v2.replace ( '0.5+0.5*TMath::Erf' , 'gauss_cdf' )
        v2 = v2.replace ( self.var1    ,  's.' + self.var1   )
        v2 = v2.replace ( self.var2    ,  's.' + self.var2   )
        #
        return eval ( 'lambda s: ' + v2 )
    
    def   fun2D  ( self ) :
        #
        v1 = str ( self.qvar1 )
        v1 = v1.replace ( '0.5+0.5*TMath::Erf' , 'gauss_cdf' )
        v1 = v1.replace ( self.var1    ,  's.' + self.var1   )
        v1 = v1.replace ( self.var2    ,  's.' + self.var2   )
        #
        v2 = str ( self.qvar2 )
        v2 = v2.replace ( '0.5+0.5*TMath::Erf' , 'gauss_cdf' )
        v2 = v2.replace ( self.var1    ,  's.' + self.var1   )
        v2 = v2.replace ( self.var2    ,  's.' + self.var2   )
        #        #
        return eval ( 'lambda s: ( ' + v1 + ',' + v2 + ' ) ' )
    
    def __repr__ ( self ) :

        result  =   'Events                       %s'       % self.num
        result += '\nStat(var1) [mean/RMS/minmax] %s/%s/%s' % ( self.stat1.mean   () ,
                                                                self.stat1.rms    () ,
                                                                self.stat1.minmax () ) 
        result += '\nStat(var2) [mean/RMS/minmax] %s/%s/%s' % ( self.stat2.mean   () ,
                                                                self.stat2.rms    () ,
                                                                self.stat2.minmax () ) 
        result += '\nCovariance  \n%s' % self.cov2 
        
        result += '\nStat(var1)           %s' % self.stat1
        result += '\nStat(var1) [values ] %s' % self.stat1.values  () 
        result += '\nStat(var1) [weights] %s' % self.stat1.weights ()
        result += '\nStat(var2)           %s' % self.stat1
        result += '\nStat(var2) [values ] %s' % self.stat1.values  () 
        result += '\nStat(var2) [weights] %s' % self.stat1.weights ()
        
        result += '\nCorrelation: %.3f%%\n%s' % ( self.corr(0,1)*100 , self.corr )
        
        result += '\nThe first  eigenvector %s '                      % self.vct0
        result += '\nThe second eigenvector %s '                      % self.vct1

        result += '\nThe 1st decorrelated variable:\n %s '            % self.nvar1 
        result += '\nThe 2nd decorrelated variable:\n %s '            % self.nvar2 
        
        result += '\nThe 1st decorrelated normalized variable:\n %s ' % self.qvar1 
        result += '\nThe 2nd decorrelated normalized variable:\n %s ' % self.qvar2 

        return result

    __str__ = __repr__ 


# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
