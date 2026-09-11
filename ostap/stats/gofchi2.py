#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/gof1dw.py
#  Set of utulities for two-sample tests and Godness-of-fit studies for (weighted) 1D-fits
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2024-09-16
# =============================================================================
""" Set of utulities for two-sample test and Godness-of-fit studies for (weighted) 1D-fits
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2024-09-16"
__all__     = (
    'Chi2'               , ## binned Chi2       GoF estimator
)
# =============================================================================
from   ostap.core.ostap_types import sequence_types 
from   ostap.core.core        import Ostap
from   ostap.utils.core       import typename 
from   ostap.stats.gofnp      import GoFnp
from   ostap.stats.utils      import ( weight_trivial ,
                                       check_all      , 
                                       num_features   ,
                                       num_samples    , nEff ) 
from   ostap.histos.axes      import axis_from_edges, h1_axis, h2_axes, h3_axes
from   ostap.logger.symbols   import chi2 as chi2_symbol
import ROOT, numpy, math   
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.gofchi2' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Two-sample & GoF (binned) Chi2 test' )
# =============================================================================
## @class Chi2
#  Binned Chi2 for Two-Sample/Goodness-of-Fit Test
#  - it works for 1,2&3 dimensions
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
class Chi2 ( GoFnp ) :
    """ Binned Chi2 for Two-Sample/Goodness-of-Fit Test 
    - it works for 1,2&3 dimensions
    """
    def __init__ ( self             , * , 
                   binning   = ()   , 
                   nToys     = 1000 , 
                   **kwargs         ) :

        if   isinstance ( binning , ROOT.TH3   ) and 3 == binning.GetDimension () :
            self.__axes = binning .GetXaxis()  , binning.GetYaxis() , binning.GetZaxis()
        elif isinstance ( binning , ROOT.TH2   ) and 2 == binning.GetDimension () :
            self.__axes = binning .GetXaxis()  , binning.GetYaxis() 
        elif isinstance ( binning , ROOT.TH1   ) and 1 == binning.GetDimension () :
            self.__axes = binning .GetXaxis()  , 
        elif isinstance ( binning , sequence_types ) and \
             1 <= len ( binning ) <= 3               and \
             all ( isinstance ( b , ROOT.TAxis )    for b in binning ) : self.__axes = tuple ( binning ) 
        elif isinstance ( binning , sequence_types ) and \
             1 <= len ( binning ) <= 3               and \
             all ( isinstance ( b , int ) and 1 < b for b in binning ) : self.__axes = tuple ( binning ) 
        elif isinstance ( binning , ROOT.TAxis )                       : self.__axes = binning , 
        elif isinstance ( binning , int        ) and 1 < binning       : self.__axes = binning ,
        elif isinstance ( binning , sequence_types ) and not binning   : self.__axes = () 
        else :
            raise TypeError ( "Invalid `binning` structure: %s/%s" % ( typename ( bininng ) , binning ) )
                        
        kwargs [ 'normalize' ] = False
        
        GoFnp.__init__ ( self ,
                         method = chi2_symbol ,
                         nToys  = nToys       , **kwargs )
        
    @property
    def dimension ( self ) :
        """`dimension` : actual chi2-dimension 
        """
        return len ( self.__axes )

    @property
    def axes ( self ) :
        """`axes` : actual binning schemes
        """
        return self.__axes

    # ==================================================================================
    @property
    def config ( self ) :
        """`config` : get all configuration parameters"""
        conf = {} 
        conf.update ( super().config  )
        conf [ 'dimension' ] = self.dimension
        conf [ 'axes'      ] = self.axes
        return conf
            
    # =========================================================================
    ## Are weights supported by this GoF estimator?
    @property 
    def weights_supported ( self ) :
        """`weights_supported`: Are weights supported by this estimator?
        """
        return True 

    # =========================================================================
    ## Good for two-samples comparison?
    #  Can this estimator be used for Two Samples comparison ?
    @property 
    def two_samples ( self ) :
        """`two_samples`: Can this estimator be used for Two Samples comparison ?
        """
        return True 
       
    # =========================================================================
    ## Calculate T-value for Goodness-of-Fit
    #  @code
    #  data1   = ...
    #  data2   = ...
    #  weight1 = ...
    #  weight2 = ...
    #  tvalue  = gof.tvalue ( data1 , data2 , weight1 , weight2 )
    #  @endcode
    def tvalue ( self      ,
                 data1     ,
                 data2     ,
                 weight1   = None  ,
                 weight2   = None  ,
                 normalize = False ) :
        """ Calculate T-value for Goodness-of-Fit
        >>> data1   = ...
        >>> data2   = ...
        >>> weight1 = ...
        >>> weight2 = ...
        >>> tvalue  = gof.tvalue ( data1 , data2 , weight1 , weight2 )
        """
        
        ## transform ?
        uds1 , uds2 = self.unpack ( data1 , data2 )

        ## check consistency of input data 
        check_all ( uds1 , uds2 , weight1 , weight2 , typename ( self ) )

        nf  = num_features ( uds1 )
        if not 1 <= nf <= 3 : raise ValueError ("Chi2: #num_features must be between 1 and 3")

        dim = self.dimension 
        
        if not dim :
            
            n1  = nEff ( uds1 , weight1 )
            n2  = nEff ( uds2 , weight2 )
            nt  = min  ( n1   , n2      )
            
            ## indicative number of events per bin 
            ne  = 16
            nn  = math.ceil ( nt ** ( 1 / nf ) ) 

            if   1 == nf :
                nb  = max ( 20 , nn ) 
                self.__axes = nb ,                
            elif 2 == nf :
                nb  = max ( 10 , nn ) 
                self.__axes = nb , nb 
            elif 3 == nf :
                nb  = max (  5 , nn ) 
                self.__axes = nb , nb , nb 
            
        dim = self.dimension                                  
        if nf != dim : raise ValueError ( "Mismatch in #num_features=%s & #dimension=%d" % ( nf , dim ) ) 
        
        w1_trivial = weight_trivial ( weight1 )
        w2_trivial = weight_trivial ( weight2 )

        w1 = 1.0 if w1_trivial else numpy.ascontiguousarray ( weight1.ravel () , dtype = numpy.float64 )
        w2 = 1.0 if w2_trivial else numpy.ascontiguousarray ( weight2.ravel () , dtype = numpy.float64 )

        ## loop over known binnings/axes 
        for i, axis in enumerate ( self.axes ) :

            if isinstance ( axis , int ) :

                d1 = uds1 if 1 == nf else uds1 [ : , i ]
                d2 = uds2 if 1 == nf else uds2 [ : , i ]
                
                from   ostap.stats.counters import WECDF
                from   ostap.math.math_base import data2vct
                
                if w1_trivial : wecdf = WECDF ( data2vct ( d1 ) , 1.0             , False ) 
                else          : wecdf = WECDF ( data2vct ( d1 ) , data2vct ( w1 ) , False )
                
                wsum1  = num_samples ( uds1 ) if w1_trivial else numpy.sum ( weight1 )
                wsum2  = num_samples ( uds2 ) if w2_trivial else numpy.sum ( weight2 )
                wscale = wsum1 / wsum2 
                
                if w2_trivial : wecdf.add ( data2vct ( d2 ) ,                 wscale   )
                else          : wecfd.add ( data2vct ( d2 ) , data2vct ( w2 * wscale ) )
                
                N  = axis
                quantiles   = wecdf.quantiles_[N-1] ()                
                axis        = axis_from_edges ( quantiles )
                if not self.silent : logger.info ( '%s: choose the binning scheme: #%d/%d %s' % ( typename ( self ) , i , nf , axis ) ) 

                del wecdf
                
                axes = list ( self.axes )
                axes [ i ]  = axis 
                self.__axes = axes 
              
        if   3 == dim : histo1 = h3_axes  ( *self.axes , double = ROOT.TH3D )
        elif 2 == dim : histo1 = h2_axes  ( *self.axes , double = ROOT.TH2D )
        elif 1 == dim : histo1 = h1_axis  ( *self.axes , double = ROOT.TH1D )

        histo2 = histo1.clone() 

        if 1 == dim :

            d1 = numpy.ascontiguousarray ( uds1 .ravel () , dtype = numpy.float64 )
            d2 = numpy.ascontiguousarray ( uds2 .ravel () , dtype = numpy.float64 )

            self.fill_1D ( histo1 , d1 , weight = w1 )
            self.fill_1D ( histo2 , d2 , weight = w2 )

        elif 2 == dim :  

            d1 = [ numpy.ascontiguousarray ( uds1 [:, i] , dtype = numpy.float64 ) for i in range ( dim ) ]
            d2 = [ numpy.ascontiguousarray ( uds2 [:, i] , dtype = numpy.float64 ) for i in range ( dim ) ]

            self.fill_2D ( histo1 , *d1 , weight = w1 )
            self.fill_2D ( histo2 , *d2 , weight = w2 )
            
        elif 3 == dim : 
                
            d1 = [ numpy.ascontiguousarray ( uds1 [:, i] , dtype = numpy.float64 ) for i in range ( dim ) ]
            d2 = [ numpy.ascontiguousarray ( uds2 [:, i] , dtype = numpy.float64 ) for i in range ( dim ) ]
            
            self.fill_3D ( histo1 , *d1 , weight = w1 )
            self.fill_3D ( histo2 , *d2 , weight = w2 )

        ## now we can calculate chi2
        h1 = histo1.density ()
        h2 = histo2.density ()
        
        ## loop over histogram bins and calculate chi2
        chi2 = 0 
        for i , j in zip ( h1 , h2 ) :
            v1   = h1 [ i ]
            v2   = h2 [ i ]
            cov2 = v1.cov2 () + v2.cov2 ()
            if 0 < cov2 :
                dv    = v1.value () - v2.value ()
                chi2 += dv * dv / cov2 

        del h1, h2, histo1, histo2
        
        return chi2 
    
    ## fill 1D histogram 
    def fill_1D ( self  , h1 , x , weight = 1.0 ) :
        n  = len ( x  )  
        sc = Ostap.fill_TH1 ( h1 , n , x , weight )
        if sc.isFailure() : raise ValueError ( "Error from Ostap.fill_TH1 %s" % sc ) 
    
    ## fill 2D histogram 
    def fill_2D ( self  , h2 , x , y , weight = 1.0 ) :
        nx = len ( x )  
        ny = len ( y )  
        if nx != ny       : raise ValueError ( "Mismatch in array sizes!" ) 
        sc = Ostap.fill_TH2 ( h2 , nx , x , y , weight )
        if sc.isFailure() : raise ValueError ( "Error from Ostap.fill_TH2 %s" % sc ) 
        
    ## fill 3D histogram 
    def fill_3D ( self  , h3 , x , y , z , weight = 1.0 ) :
        nx = len ( x )  
        ny = len ( y )  
        nz = len ( z )  
        if nx != ny or ny != nz  : raise ValueError ( "Mismatch in array sizes!" ) 
        sc = Ostap.fill_TH3 ( h3 , nx , x , y , z , weight )
        if sc.isFailure() : raise ValueError ( "Error from Ostap.fill_TH3 %s" % sc ) 
                    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
