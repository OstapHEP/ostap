#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/average.py
#  Simple averages for several (independent) inconsistend measurmenets
#  @see  M.Trassinelli and M.Maxton,
#        ``A minimalistic and general weighted averaging method for inconsistent data'',
#        https://arxiv.org/abs/2406.08293
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2024-06-24
# =============================================================================
""" Simple averages for several (independent) inconsistent measurmenet
 - see  M.Trassinelli and M.Maxton,
        ``A minimalistic and general weighted averaging method for inconsistent data'',
        https://arxiv.org/abs/2406.08293
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2024-06-24"
__version__ = ""
# =============================================================================
__all__     = (
    ## as objects 
    'Standard'               , ## standard inverse variance weighted average 
    'Birge'                  , ## Birge's/scaled inverse-varinace weighted average
    'Conservative'           , ## Conservative average 
    'Jeffreys'               , ## JEffreys's prior average
    ## as functions 
    'standard_average'       , ## standard     average
    'birge_average'          , ## scaled/Birge average
    'conservative_average'   , ## ``conservative'' average
    'jeffreys_average'       , ## Jeffrey limit  average
)
# =============================================================================
from   ostap.math.ve          import VE
from   ostap.math.base        import iszero
import ROOT, math, abc 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.average' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
## trivial helper function to expand the arguments 
def ve_expand ( *values ) :
    result = []
    for item in values :
        if isinstance ( item , VE ) : result .  append    ( item )
        else :                        result += ve_expand ( item )
    return result

# =============================================================================
## @class Average
#  Helper abstract base class to implement varios averages
class Average(object):
    """Helper abstarct base class to implement average for (inconsistent) data
    """
    def __init__ ( self , *data ) :
        
        self.__data = tuple ( ve_expand ( *data ) )
        
        assert 2 <= len ( self.__data ) and all ( ( isinstance ( v , VE ) and 0 < v.cov2() ) for v in self.__data ) , \
            "Invalid input data: either not ValueWithjError or invalid uncertainty!"

        self.__xmin   = min ( ( v.value() - 3 * v.error() ) for v in self.__data )
        self.__xmax   = max ( ( v.value() + 3 * v.error() ) for v in self.__data )
        self.__smin   = min ( v.cov2()                      for v in self.__data ) ** 0.5 
        
        self.__result = None
        self.__pvalue = None
        self.__graph  = None
        
    # =========================================================================
    ## Get p-value from the chi2-distribution
    @property
    def pvalue ( self ) :
        """'pvalue' : p-value (from chi2-distribution)"""
        if self.__pvalue is None :
            average       = self.weighted_average()
            value         = average.value() 
            self.__pvalue = self.p_value ( value ) 
        return self.__pvalue

    # =======================================================================
    ## Calcuaet the p-value form chi2 distributon 
    def p_value ( self , mu ) :
        """Calculate the p-value"""        
        ## get chi2 value 
        c2      = self.chi2 ( mu )
        k       = len ( self ) 
        return Ostap.Math.gamma_inc_P ( k / 2.0 , c2 / 2.0 )  
        
        
    @property
    def values  ( self ) :
        """'values' : input data, array of `ValueWithError` objects"""
        return self.__data
    @property
    def data    ( self ) :
        """'data' : input data, array of `ValueWithError` objects"""
        return self.__data

    # ============================================================================
    ## Get the associated graph    
    @property
    def graph ( self ) :
        """'graph' : get the associated graph"""
        if self.__graph is None :
            import ostap.histos.graphs
            self.__graph = ROOT.TGraphErrors()
            for i,v in enumerate ( self.values ) :
                x = v
                y = VE ( i + 0.5 , 0 )
                self.__graph.append ( x , y )
        return self.__graph 

    # ============================================================================
    ## make the NLL graph 
    def draw_nll ( self , opts = '' , xmin = None , xmax = None , **kwargs ) :
        from ostap.math.models import f1_draw
        if xmin is None : xmin = self.__xmin
        if xmax is None : xmax = self.__xmax
        return f1_draw ( self , opts , xmin = xmin , xmax = xmax , **kwargs ) 
        
    # ============================================================================
    ## len of data 
    def __len__  ( self ) : return len ( self.__data )

    # ============================================================================
    ## iterator over the data 
    def __iter__ ( self ) :
        """"iterator over the data"""        
        for v in self.__data : yield v
        
    # ==========================================================================
    ## Get the value of negative log-likelihood 
    @abc.abstractmethod
    def NLL ( self , mu ) :
        """Get the value of negative log-likelihood"""
        return None 

    @property 
    def result ( self ) :
        """'result' : the actual average"""
        if self.__result is None :
            self.__result = self.get_average ()
        return self.__result

    # ============================================================================
    ## Get the negatile log-likelihood 
    def __call__ ( self , mu ) :
        """Get the negatile log-likelihood"""
        return self.NLL ( mu ) 
        
    # ==========================================================================
    ## calculate the aaverage 
    def get_average ( self ) :
        """Calculate the average"""
        
        from ostap.math.minimize import minimize_scalar
        result = minimize_scalar ( self , bracket = ( self.__xmin , self.__xmax ) )
        x = result.x

        from ostap.math.derivative import Derivative2
        N  = len ( self ) 
        d2 = Derivative2 ( self , step = self.__smin / ( 100 * N ) ) 
        e2 = d2 ( x ) 
        
        c2 = 1.0 / e2 
        return VE ( x , c2 )

    ## get a simple (inverse-variance) weighted average 
    def weighted_average ( self ) :
        """Get a simple (inverse-variance) weighted average"""        
        wtotal = 0.0
        wvalue = 0.0
        for v in self.values :
            weight = 1.0 / v.cov2()
            wvalue += v.value() * weight
            wtotal +=             weight
        wvalue /= wtotal        
        return VE ( wvalue , 1.0 / wtotal ) 

    # ==============================================================================
    ## get the chi2 value 
    def chi2 ( self , x  = None ) :
        """Get the chi2 value 
        """
        if x is None : x = self.weighted_average().value()  
        c2 = 0.0 
        for v in self.values :
            c2 += ( v.value() - x ) ** 2 / v.cov2 ()
        return c2
            
# ===================================================================================
## @class Standard
#  The standard average (inverse variance)
#  This one just for cmpleteness and comparisons 
class Standard(Average) :
    """The standard average (inverse variance)
      This one just for cmpleteness and comparisons 
    """

    ## get the negative log-likelihood 
    def NLL ( self , mu ) :
        """Get the negative log-likelihood"""
        return 0.5 * self.chi2 ( mu )

    ## get an average 
    def get_average ( self ) :
        return self.weighted_average ()

# =============================================================================
## @class Birge
#  Birge (scaled inverse variance) average fo rinconcistent measurements
#  @see Birge, Raymond T
#       ``The Calculation of Errors by the Method of Least Squares'',
#        Phys. Rev., 40 (1932) 207--227, doi = 10.1103/PhysRev.40.207},
#  @see https://link.aps.org/doi/10.1103/PhysRev.40.207}
#  @see M.Trassinelli and M.Maxton,
#       "A minimalistic and general weighted averaging method for inconsistent data",
#  @see https://arxiv.org/abs/2406.08293
class Birge(Standard) :
    """ Birge (scaled inverse-variance) average for inconcistent  measurements
    - see Birge, Raymond T
       ``The Calculation of Errors by the Method of Least Squares'',
        Phys. Rev., 40 (1932) 207--227, doi = 10.1103/PhysRev.40.207},
    - see https://link.aps.org/doi/10.1103/PhysRev.40.207}
    - see M.Trassinelli and M.Maxton,
           "A minimalistic and general weighted averaging method for inconsistent data",
    - see https://arxiv.org/abs/2406.08293
    """    
    def __init__ ( self , *data ) :        
        super(Birge,self).__init__( *data )        
        self.__scale = None
        
    @property 
    def scale ( self ) :
        """'scale' : Birge's scale factor"""
        if self.__scale is None :
            self.get_average() 
        return self.__scale

    ## get an average 
    def get_average ( self ) :
        """Gen a scaled average
        """
        
        result       = self.weighted_average ()
        value , cov2 = result.value() , result.cov2()
        
        N            = len ( self ) 
        chi2         = self.chi2 ( value ) 
        ratio        = chi2 / ( N - 1 ) 
        self.__scale = ratio 
        
        return VE ( value , ratio * cov2 )

# ==============================================================================
## @class Conservative
## Get the negative log-likelihood for ``conservative`` estimate for the
#  average of several inconsistent measurments 
#  @see Sivia D.S. and Skilling J,
#      ``Data analysis: a Bayesian tutorial'',
#       Oxford University Press. 2006
#  @see https://books.google.it/books/about/Data_Analysis.html?id=6O8ZAQAAIAAJ&redir_esc=y
#  @see M.Trassinelli and M.Maxton,
#       "A minimalistic and general weighted averaging method for inconsistent data",
#  @see https://arxiv.org/abs/2406.08293
class Conservative(Average) :
    """Get the negative log-likelihood for ``conservative`` estimate for the
       average of several inconsistent measurments 
    - see Sivia D.S. and Skilling J,
          ``Data analysis: a Bayesian tutorial'',
        Oxford University Press. 2006
    - see https://books.google.it/books/about/Data_Analysis.html?id=6O8ZAQAAIAAJ&redir_esc=y
    - see M.Trassinelli and M.Maxton,
           "A minimalistic and general weighted averaging method for inconsistent data",
    - see https://arxiv.org/abs/2406.08293
    """
    _c1 = math.log ( 2.0 * math.sqrt ( 2.0 * math.pi ) ) 

    ## get the negative log-likelihood 
    def NLL ( self , mu ) :
        """Get the negative log-likelihood"""
        
        ##  \f$ \frac{ 1- e^{-x}}{x} \f$
        def term ( x2 ) :
            """ Expression ( 1 - exp(-x) ) / x """            
            if   abs ( x2 ) < 1.e-15 : return  1.0 - 0.5 * x2 
            elif abs ( x2 ) < 1      : return -math.expm1(-x2)/x2 
            return ( 1.0 - math.exp ( -x2 ) ) / x2
        
        nll = 0.0
        for v in self.values :
            vv , cov2 = v.value() , v.cov2() 
            x2   = 0.5 * ( vv - mu ) ** 2 / cov2
            dd   = math.log ( term ( x2 ) ) - 0.5 * math.log ( cov2 ) - self._c1 
            nll -= dd
            
        return nll 

# =============================================================================
## @class Jeffreys
#  Jeffreys limit for the avergae of inconsistent data
#  @see M.Trassinelli and M.Maxton,
#       "A minimalistic and general weighted averaging method for inconsistent data",
#  @see https://arxiv.org/abs/2406.08293
class Jeffreys(Average) :
    """Jeffreys limit for the avergae of inconsistent data
    - see M.Trassinelli and M.Maxton,
           "A minimalistic and general weighted averaging method for inconsistent data",
    - see https://arxiv.org/abs/2406.08293
    """
    _c2 = 2.0 / math.sqrt ( math.pi )
    _c3 = 1.0 / math.sqrt ( 2.0     )
    _c4 = 2.0 * math.sqrt ( 2.0     ) 

    ## get the negative log-likelihood 
    def NLL ( self , mu ) :
        """Get the negative log-likelihood"""
        
        ## expression erf(x)/x 
        def term ( z ) :
            """Expression erfc(x)/x"""        
            if   abs ( z ) < 1.e-10 : return self._c2 * ( 1 - z*z/3 )
            return                           self._c2 * ( math.erf  ( z ) / z ) 
        
        nll = 0.0
        for v in self.values :        
            vv , e = v.value() , v.error() 
            x    = self._c3 * ( vv - mu ) / e         
            dd   = math.log ( term ( x ) ) - math.log ( self._c4 * e  )  
            nll -= dd
            
        return nll 
            
## ============================================================================
## Make a simple average of several measurements
#  @code
#  val1  = ...
#  val2  = ...
#  ave   = standard_average ( va1 , val2 ) 
#  @endcode
def standard_average ( *values ) :
    """Make a simple average of severla measurements
    >>> val1  = ...
    >>> val2  = ...
    >>> val3  = ...
    >>> ave   = standard_average ( va1 , val2 , val3 ) 
    """
    return Standard ( *values ).result

## ============================================================================
## Make a scaled/Birge average of several inconsistent measurements
#  @code
#  val1            = ...
#  val2            = ...
#  average , scale = birge_average ( va1 , val2 ) 
#  @endcode
#  @see Birge, Raymond T
#       ``The Calculation of Errors by the Method of Least Squares'',
#        Phys. Rev., 40 (1932) 207--227, doi = 10.1103/PhysRev.40.207},
#  @see https://link.aps.org/doi/10.1103/PhysRev.40.207}
#  @see M.Trassinelli and M.Maxton,
#       "A minimalistic and general weighted averaging method for inconsistent data",
#  @see https://arxiv.org/abs/2406.08293
def birge_average ( *values ) :
    """Make a scaled/Birge average of several inconsisten measurements
    >>> val1            = ...
    >>> val2            = ...
    >>> average , scale = birge_average ( va1 , val2 ) 
    -  @see Birge, Raymond T
            ``The Calculation of Errors by the Method of Least Squares'',
            Phys. Rev., 40 (1932) 207--227, doi = 10.1103/PhysRev.40.207},
    - see https://link.aps.org/doi/10.1103/PhysRev.40.207}
    - see  M.Trassinelli and M.Maxton,
           ``A minimalistic and general weighted averaging method for inconsistent data'',
    - see https://arxiv.org/abs/2406.08293
    """

    b      = Birge ( *values )
    result = b.result
    scale  = b.scale
    return result , math.sqrt ( scale ) 

## ============================================================================
## Make a ``conservative'' average of several inconsistent measurements
#  @code
#  val1  = ...
#  val2  = ...
#  ave   = conservative_average ( va1 , val2 ) 
#  @endcode
#  @see Sivia D.S. and Skilling J,
#       ``Data analysis: a Bayesian tutorial'',
#        Oxford University Press. 2006
#  @see https://books.google.it/books/about/Data_Analysis.html?id=6O8ZAQAAIAAJ&redir_esc=y
#  @see M.Trassinelli and M.Maxton,
#       "A minimalistic and general weighted averaging method for inconsistent data",
#        https://arxiv.org/abs/2406.08293
def conservative_average ( *values ) :
    """Make a ``conservative'' average of several inconsistent measurements
    >>> val1 = ...
    >>> val2 = ...
    >>> ave  = conservative_average ( va1 , val2 ) 
    - see Sivia D.S. and Skilling J,
        ``Data analysis: a Bayesian tutorial'',
         Oxford University Press. 2006
    - see https://books.google.it/books/about/Data_Analysis.html?id=6O8ZAQAAIAAJ&redir_esc=y
    - see  M.Trassinelli and M.Maxton,
    ``A minimalistic and general weighted averaging method for inconsistent data'',
    - see https://arxiv.org/abs/2406.08293
    """
    
    return Conservative ( *values ).result 

## ============================================================================
## Make ``Jeffrey's limit'' average of several inconsistent measurements
#  @code
#  val1  = ...
#  val2  = ...
#  ave   = jeffreys_average ( va1 , val2 ) 
#  @endcode
#  @see https://books.google.it/books/about/Data_Analysis.html?id=6O8ZAQAAIAAJ&redir_esc=y
#  @see M.Trassinelli and M.Maxton,
#       "A minimalistic and general weighted averaging method for inconsistent data",
#        https://arxiv.org/abs/2406.08293
def jeffreys_average ( *values ) :
    """Make ``Jeffrey's limit'' average of several inconsistent measurements
    >>> val1 = ...
    >>> val2 = ...
    >>> ave  = jeffreys_average ( va1 , val2 ) 
    - see https://books.google.it/books/about/Data_Analysis.html?id=6O8ZAQAAIAAJ&redir_esc=y
    - see  M.Trassinelli and M.Maxton,
    ``A minimalistic and general weighted averaging method for inconsistent data'',
    - see https://arxiv.org/abs/2406.08293
    """
    return Jeffreys ( *values ).result
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme  import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END
# =============================================================================
