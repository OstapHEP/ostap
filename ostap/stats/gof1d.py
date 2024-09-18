#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/gof_1d.py
#  Set of utulities for goodness-of-fit studies for 1D-fits
#  @see https://doi.org/10.1111/1467-9868.00337 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2024-09-16
# =============================================================================
""" Simple utulities for goodness-of-fit studies for 1D-fits  
- see https://doi.org/10.1111/1467-9868.00337
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2024-09-16"
__all__     = (
    )
# =============================================================================
from   collections            import defaultdict, namedtuple
from   ostap.core.meta_info   import root_info 
from   ostap.fitting.funbasic import AFUN1
from   ostap.fitting.pdfbasic import PDF1
from   ostap.core.core        import SE, VE, Ostap
from   ostap.math.base        import doubles, axis_range  
from   ostap.math.models      import f1_draw 
from   ostap.math.math_ve     import significance
import ostap.fitting.ds2numpy 
import ostap.fitting.roofit
import ROOT, math  
# =============================================================================
try :    
    import numpy as np
except ImportError :
    np = None
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.gof1d' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Simple utilities for goodness-of-1D-fit studies' )
# =============================================================================
if  ( 6 , 32 ) <= root_info : data2vct = lambda s : s
else                        : data2vct = lambda s : doubles ( s ) 
# =============================================================================
## Get Kolmogorov-Smirnov statistis KS
#  @code
#  cdf_data =...
#  ks2  = kolmogorov_smirnov ( cdf_data )
#  @endcode
#  @param cdf_data sorted array of F0(X_i) - values of CDF at X data points
#  @return Kolmogorov-Smirnov statistisc KS
def kolmogorov_smirnov ( cdf_data ) :
    """ Get Kolmogorov-Smirnov statistis  KS
    - `cdf_data` : sorted array of F0(X_i) - values of CDF at X data points
    >>> cdf_data =...
    >>> ks2  = kolmogorov_smirnov ( cdf_data )
    """
    n      = len ( cdf_data ) 
    result = max ( max ( ( i + 1.0 ) / n - Fi , Fi - float ( i ) / n ) for ( i, Fi )  in enumerate ( cdf_data )  ) 
    return math.sqrt ( result ) 
# =============================================================================
## Get Anderson-Darling  statistiscs AD^2
#  @code
#  cdf_data =...
#  ad2      = anderson_darling ( cdf_data )
#  @endcode
#  @param cdf_data sorted array of F0(X_i) - values of CDF at X data points
#  @return Anderson-Darling statistisc AD^2
def anderson_darling ( cdf_data ) :
    """ Get Anderson-Darling statistiscs AD^2
    - `cdf_data` : sorted array of F0(X_i) - values of CDF at X data points
    >>> cdf_data =...
    >>> ad2      = anderson_darling ( cdf_data )    
    """
    n       = len ( cdf_data ) 
    flog    = math.log
    result  = sum ( ( i + 0.5 ) * flog ( Fi ) + ( n - i -  0.5 ) * flog ( 1 - Fi ) for ( i , Fi )  in enumerate ( cdf_data ) ) 
    result *= -2.0 / n
    result -= n
    return result 
# =============================================================================
## Get Cramer-von Mises statistics CM^2
#  @code
#  cdf_data = ...
#  cm2      = cramer_von_mises ( cdf_data )
#  @endcode
#  @param cdf_data sorted array of F0(X_i) - values of CDF at X data points
#  @return Cramer-von Mises statistisc CM^2
def cramer_von_mises ( cdf_data  ) :
    """ Get Cramer-von Mises statistics CM^2
    - `cdf_data` : sorted array of F0(X_i) - values of CDF at X data points
    >>> cdf_data =...
    >>> cm2     = cramer_von_mises ( cdf_data )    
    """
    n       = len ( cdf_data ) 
    result  = sum ( ( Fi - ( i + 0.5 ) / n ) ** 2 for ( i, Fi ) in enumerate ( cdf_data ) ) 
    result += 12.0 / n
    return result
# =============================================================================
## Get ZK statististics
#  @code
#  cdf_data = ...
#  zk       = ZK ( cdf_data )
#  @endcode
#  @param cdf_data sorted array of F0(X_i) - values of CDF at X data points
#  @return ZK statistisc ZK 
def ZK  ( cdf_data ) :
    """ Get ZK statististics 
    - `cdf_data` : sorted array of F0(X_i) - values of CDF at X data points
    >>> cdf_data =...
    >>> zk       = ZK ( cdf_data )    
    """
    n      = len ( cdf_data ) 
    flog   = math.log 
    result = max ( ( i     + 0.5 ) * flog ( ( i + 0.5     ) / ( n *       Fi   ) ) +
                   ( n - i - 0.5 ) * flog ( ( n - i - 0.5 ) / ( n * ( 1 - Fi ) ) ) for ( i , Fi ) in enumerate ( cdf_data ) )
    return result 
# =============================================================================
## Get ZA statististics
#  @code
#  cdf_data = ...
#  za       = ZA ( cdf_data )
#  @endcode
#  @param cdf_data sorted array of F0(X_i) - values of CDF at X data points
#  @return ZA statistisc ZA 
def ZA  ( cdf_data ) :
    """ Get ZA statististics 
    - `cdf_data` : sorted array of F0(X_i) - values of CDF at X data points
    >>> cdf_data =...
    >>> za       = ZA ( cdf_data )    
    """
    n       = len ( cdf_data ) 
    flog    = math.log
    result  = sum ( flog ( Fi) / ( n - i - 0.5 ) + flog ( 1 - Fi ) / ( i + 0.5 ) for ( i , Fi )  in enumerate ( cdf_data ) )
    result *= -1 
    return result    
# =============================================================================
## Get ZC statististics
#  @code
#  cdf_data = ...
#  zc       = ZC ( cdf_data )
#  @endcode
#  @param cdf_data sorted array of F0(X_i) - values of CDF at X data points
#  @return ZC statistisc ZC
def ZC  ( cdf_data ) :
    """ Get ZC statististics 
    - `cdf_data` : sorted array of F0(X_i) - values of CDF at X data points
    >>> cdf_data =...
    >>> zc       = ZC ( cdf_data )    
    """
    n      = len ( cdf_data ) 
    flog   = math.log
    result = sum ( ( flog ( ( 1.0 / Fi - 1 ) / ( ( n - 0.5 ) / ( i + 0.25 ) - 1 ) ) ) ** 2 for ( i , Fi ) in enumerate ( cdf_data ) )
    return result 
    
# =============================================================================
## @class GoF1D
#  Goodness of 1D-fits 
#  @see https://doi.org/10.1111/1467-9868.00337
class GoF1D(object) :
    """ Goodness of 1D-fits 
    - see https://doi.org/10.1111/1467-9868.00337
    """
    def __init__ ( self    ,
                   pdf     ,
                   dataset ) :
        
        assert isinstance ( pdf     , PDF1            ) , 'Invalid type of `PDF`:%s'     % type ( pdf     )
        assert isinstance ( dataset , ROOT.RooDataSet ) , 'Invalid type of `daatset`:%s' % type ( dataset )
        
        vars = pdf.vars
        assert 1 == len ( vars )   , 'GoF1D: Only 1D-pdfs are allowed!'
        assert pdf.xvar in dataset , 'GoF1D: `xvar`:%s is not in dataset!' % ( self.xvar.name ) 
        
        cdf = pdf.cdf ()
        if hasattr ( pdf.pdf , 'function' ) :
            fun = pdf.pdf.function()
            if hasattr ( fun , 'cdf' ) :
                cdf = fun.cdf() 

        self.__cdf   = cdf
        self.__xmnmx = pdf.xminmax()
        
        ## vectorized version of CDF 
        vct_cdf = np.vectorize ( cdf )
        
        ## data in form of numpy sructured array
        varname = pdf.xvar.name 
        data    = dataset.tonumpy ( varname ) [ varname ] 

        ## sorted data 
        self.__data     = np.sort  ( data )

        ## empirical CDF function 
        self.__ecdf     = Ostap.Math.ECDF ( data2vct ( self.__data ) )

        ## evalute CDF for sorted data 
        self.__cdf_data = vct_cdf ( self.__data )
 
        self.__KS = kolmogorov_smirnov ( self.__cdf_data )
        self.__AD = anderson_darling   ( self.__cdf_data )
        self.__CM = cramer_von_mises   ( self.__cdf_data )
        self.__ZK = ZK                 ( self.__cdf_data )
        self.__ZA = ZA                 ( self.__cdf_data )
        self.__ZC = ZC                 ( self.__cdf_data )
       
        
        del vct_cdf 

    # =========================================================================
    ## size of dataset
    @property 
    def N  ( self ) :
        """`N` : size of dataset"""
        return len ( self.__data ) 

    # =========================================================================
    ## empirical CDF for data
    def ecdf ( self ) :
        """`ecdf` : empirical CDF for data"""
        return self.__ecdf
    
    # =========================================================================
    ## Get Kolmogorov-Smirnov statistiscs 
    @property 
    def kolmogorov_smirnov ( self ) :
        """ Get Kolmogorov-Smirnov statistiscs KS
        """
        return self.__KS 

    # ===============================================
    ## Get Anderson-Darling  statistiscs 
    @property 
    def anderson_darling ( self ) :
        """ Get Anderson-Darling statistiscs 
        """
        return self.__AD 

    # =========================================================================
    ## Get Cramer-von Mises statistics 
    @property 
    def cramer_von_mises ( self ) :
        """ Get Cramer-von Mises statistics 
        """
        return self.__CM 
    
    # =========================================================================
    ## Get ZK statististics 
    @property 
    def ZK  ( self ) :
        """ Get ZK statistics
        """
        return self.__ZK 
        
    # =========================================================================
    ## Get ZA statististics 
    @property 
    def ZA  ( self ) :
        """ Get ZA statistics
        """        
        return self.__ZA
    
    # =========================================================================
    ## Get ZC statististics 
    @property 
    def ZC  ( self ) :
        """ Get ZC statistics
        """        
        return self.__ZC 
            
    @property
    def data ( self ) :
        """`data` : sorted  input data (as numpy array)"""
        return self.__data

    @property
    def cdf_data ( self ) :
        """`cdata` : vector of CDF(x) data (as numpy array)"""
        return self.__cdf_data 
    
    # ==========================================================================
    ## Print the summary as Table
    def table ( self , title = '' , prefix = '' , width = 5 , precision = 3 ) :
        """ Print the summary as Table """
        rows = [ ( 'Statistics' , 'Value' ) ] 
        from   ostap.logger.pretty import pretty_float
        import ostap.logger.table  as     T 

        ks , expo = pretty_float ( self.kolmogorov_smirnov , width = width , precision = precision )
        if expo : row = 'Kolmogorov-Smirnov' , ks , '10^%+d' % expo
        else    : row = 'Kolmogorov-Smirnov' , ks              
        rows.append ( row )
        
        ad , expo = pretty_float ( self.anderson_darling , width = width , precision = precision )
        if expo : row = 'Anderson-Darling' , ad , '10^%+d' % expo
        else    : row = 'Anderson-Darling' , ad              
        rows.append ( row )

        cm , expo = pretty_float ( self.cramer_von_mises  , width = width , precision = precision )
        if expo : row = 'Cramer-von Mises' , cm , '10^%+d' % expo
        else    : row = 'Cramer-von Mises' , cm              
        rows.append ( row )

        zk , expo = pretty_float ( self.ZK , width = width , precision = precision )
        if expo : row = 'ZK' , zk, '10^%+d' % expo
        else    : row = 'ZK' , zk             
        rows.append ( row )

        za , expo = pretty_float ( self.ZA , width = width , precision = precision )
        if expo : row = 'ZA' , za, '10^%+d' % expo
        else    : row = 'ZA' , za             
        rows.append ( row )

        zc , expo = pretty_float ( self.ZC , width = width , precision = precision )
        if expo : row = 'ZC' , zc, '10^%+d' % expo
        else    : row = 'ZC' , zc             
        rows.append ( row )

        title = title if title else 'Goodness of 1D-fit' 
        return T.table ( rows , title = title , prefix = prefix , alignment = 'lcl' )

    __repr__ = table
    __str__  = table
    
    # =========================================================================
    ## Draw fit CDF & empirical ECDF 
    def draw  ( self , *args , **kwargs ) :
        """ Draw fit CDF & empirical CDF
        """
        cdf  = self.__cdf
        ecdf = self.__ecdf
        
        xmin, xmax =  ecdf.xmin () , ecdf.xmax ()
        
        if hasattr ( cdf , 'xminmax' ) :
            xmn , xmx =  cdf.xminmax()
            xmin = min ( xmin , xmn )
            xmax = max ( xmax , xmx )
        else :
            if hasattr ( cdf , 'xmin' ) : xmin = min ( xmin , cdf.xmin () )
            if hasattr ( cdf , 'xmax' ) : xmax = max ( xmax , cdf.xmax () )

        xmin , xmax = axis_range ( xmin , xmax , delta = 0.05 )
        
        xmin = kwargs.pop ( 'xmin' , xmin )
        xmax = kwargs.pop ( 'xmax' , xmax )

        if isinstance ( cdf , AFUN1 ) :
            self.__frame = cdf.draw ()
            self.__frame.draw ()
        else : 
            self.__draw_fun = lambda x : cdf ( x ) 
            f1_draw   ( self.__draw_fun , color = ROOT.kOrange + 1 , linewidth = 3 , xmin = xmin , xmax = xmax , **kwargs ) 

        ecdf.draw ( 'same' , color = 2 , linewidth = 3 , xmin = xmin , xmax = xmax , **kwargs )

# =============================================================================
## @class GoF1DToys
#  Check Goodness of 1D-fits using toys 
class GoF1DToys(object) :
    """ Check Goodness of 1D-fits ysing toys 
    """
    ## result of GoGptoys 
    Result = namedtuple ( 'Result' , 'statistics counter pvalue nsigma' )
    # =========================================================================
    def __init__ ( self           ,
                   pdf            ,
                   dataset        ,
                   nToys  = 1000  ,
                   silent = False ) :
        
        assert isinstance ( pdf     , PDF1            ) , 'Invalid type of `PDF`:%s'     % type ( pdf     )
        assert isinstance ( dataset , ROOT.RooDataSet ) , 'Invalid type of `daatset`:%s' % type ( dataset )
        
        vars = pdf.vars
        assert 1 == len ( vars )   , 'GoF1D: Only 1D-pdfs are allowed!'
        assert pdf.xvar in dataset , 'GoF1D: `xvar`:%s is not in dataset!' % ( self.xvar.name ) 
        
        cdf = pdf.cdf ()
        if hasattr ( pdf.pdf , 'function' ) :
            fun = pdf.pdf.function()
            if hasattr ( fun , 'cdf' ) :
                cdf = fun.cdf() 
                
        ## vectorized version of CDF 
        vct_cdf = np.vectorize ( cdf )
        
        ## data in form of numpy sructured array
        varname = pdf.xvar.name 
        data    = dataset.tonumpy ( varname ) [ varname ] 

        ## sorted data 
        self.__data     = np.sort  ( data )

        ## evalute CDF for sorted data 
        self.__cdf_data = vct_cdf ( self.__data )

        self.__vct_cdf = vct_cdf
        self.__cdf     = cdf
        self.__pdf     = pdf
                
        self.__ecdfs   = {} 

        self.__KS = kolmogorov_smirnov ( self.__cdf_data )
        self.__AD = anderson_darling   ( self.__cdf_data )
        self.__CM = cramer_von_mises   ( self.__cdf_data )
        self.__ZK = ZK                 ( self.__cdf_data )
        self.__ZA = ZA                 ( self.__cdf_data )
        self.__ZC = ZC                 ( self.__cdf_data )

        self.__KS_cnt = SE ()
        self.__AD_cnt = SE ()
        self.__CM_cnt = SE ()
        self.__ZK_cnt = SE ()
        self.__ZA_cnt = SE ()
        self.__ZC_cnt = SE ()

        self.__nToys  = 0 
        if 0 < nToys : self.run ( nToys , silent = silent ) 

    # ===============================================================================
    ## run toys 
    def run ( self , nToys = 1000 , silent = False ) :
        """ Run toys 
        """ 
        assert isinstance ( nToys , int ) and 0 < nToys , "Invalid `nToys` argument!"

        N       = len ( self.__data )
        varname = self.__pdf.xvar.name
        
        results = defaultdict(list)

        from ostap.utils.progress_bar import progress_bar 
        for i in progress_bar ( nToys , silent = silent ) :

            dset     = self.__pdf.generate ( N , sample = True )
            data     = dset.tonumpy ( varname ) [ varname ] 
            data     = np.sort ( data )
            cdf_data = self.__vct_cdf ( data )
            
            ks       = kolmogorov_smirnov ( cdf_data )
            ad       = anderson_darling   ( cdf_data )
            cm       = cramer_von_mises   ( cdf_data )
            zk       = ZK                 ( cdf_data )
            za       = ZA                 ( cdf_data )
            zc       = ZC                 ( cdf_data )
            
            self.__KS_cnt += ks 
            self.__AD_cnt += ad 
            self.__CM_cnt += cm 
            self.__ZK_cnt += zk 
            self.__ZA_cnt += za 
            self.__ZC_cnt += zc 

            results [ 'KS'  ].append ( ks )    
            results [ 'AD'  ].append ( ad ) 
            results [ 'CM'  ].append ( cm ) 
            results [ 'ZK'  ].append ( zk ) 
            results [ 'ZA'  ].append ( za ) 
            results [ 'ZC'  ].append ( zc ) 

            ## delete data
            if isinstance ( dset , ROOT.RooDataSet ) : 
                dset = Ostap.MoreRooFit.delete_data ( dset )
                
            del dset
            del data
            del cdf_data            

        ## accumulate number of toys 
        self.__nToys += nToys 

        ECDF = Ostap.Math.ECDF 
        for key in results :
            data = results [ key ]
            if not data : continue
            if not key in self.__ecdfs : self.__ecdfs [ key ] = ECDF ( data , True ) ## complementary ECDF!
            else                       : self.__ecdfs [ key ]  .add  ( data2vct ( data ) ) 
            
        del results 
    
    # =========================================================================
    ## number of toys 
    @property
    def nToys ( self ) :
        """`nToys` : number of toys"""
        return self.__nToys
    
    # =========================================================================
    ## ECDFs
    @property 
    def ecdfs ( self ) :
        """`ecdfs` : toy results as empirical cumulative distribution functions"""
        return self.__ecdfs

    # =========================================================================
    ## Helper methof to get result
    def result ( self , value , counter , label ) :
        """ Helper method to get result 
        """
        if self.__ecdfs and label in self.__ecdfs :
            pvalue = self.__ecdfs [ label ] . estimate ( value )
            nsigma = significance ( pvalue )            
        else : 
            pvalue = VE ( -1 , 0 )
            nsigma = VE ( -1 , 0 )
        ## 
        return self.Result ( value   ,
                             counter ,
                             pvalue  ,
                             nsigma  )
            
    # =========================================================================
    ## Get Kolmogorov-Smirnov statistiscs 
    @property 
    def kolmogorov_smirnov ( self ) :
        """ Get Kolmogorov-Smirnov statistiscs KS
        """
        return self.result ( self.__KS , self.__KS_cnt , 'KS' ) 
    
    # ===============================================
    ## Get Anderson-Darling  statistiscs 
    @property 
    def anderson_darling ( self ) :
        """ Get Anderson-Darling statistiscs 
        """
        return self.result ( self.__AD , self.__AD_cnt , 'AD' ) 
        
    # =========================================================================
    ## Get Cramer-von Mises statistics 
    @property 
    def cramer_von_mises ( self ) :
        """ Get Cramer-von Mises statistics 
        """
        return self.result ( self.__CM , self.__CM_cnt , 'CM' ) 
        
    # =========================================================================
    ## Get ZK statististics 
    @property 
    def ZK  ( self ) :
        """ Get ZK statistics 
        """
        return self.result ( self.__ZK , self.__ZK_cnt , 'ZK' ) 
        
    # =========================================================================
    ## Get ZA statististics 
    @property 
    def ZA  ( self ) :
        """ Get ZA statistics
        """        
        return self.result ( self.__ZA , self.__ZA_cnt , 'ZA' ) 
    
    # =========================================================================
    ## Get ZC statististics 
    @property 
    def ZC  ( self ) :
        """ Get ZC statistics 
        """        
        return self.result ( self.__ZC , self.__ZC_cnt , 'ZC' ) 
            
    # =========================================================================
    ## format a row in the table
    def _row  ( self , what , result , width = 5 , precision = 3 ) :
        """ Format a row in the table
        """

        value      = result.statistics
        counter    = result.counter
        pvalue     = result.pvalue
        nsigma     = result.nsigma
        
        mean       = counter.mean   ()
        rms        = counter.rms    () 
        vmin, vmax = counter.minmax () 
        
        mxv        = max ( value , mean.value() , mean.error() , rms , vmin , vmax )
        
        from   ostap.logger.pretty import fmt_pretty_ve 
        
        fmt, fmtv , fmte , expo = fmt_pretty_ve ( VE ( mxv ,  mean.cov2() ) ,
                                                  width       = width       ,
                                                  precision   = precision   , 
                                                  parentheses = False       )
        
        if expo : scale = 10**expo
        else    : scale = 1
        fmt2 = '%s/%s' % ( fmtv , fmtv ) 
        
        return ( what  ,
                 fmtv  %  ( value / scale )                    ,
                 ( mean / scale ).toString ( fmt )             ,
                 fmtv  %  ( rms  / scale )                     ,
                 fmt2  %  ( vmin / scale , vmax / scale )      ,
                 ( '10^%+d' % expo  if expo else '' )          ,                  
                 ( 100 * pvalue ) .toString ( '%.2f +/- %-.2f' ) , 
                 ( nsigma       ) .toString ( '%.1f +/- %-.1f' ) ) 
    
    # =========================================================================
    ## Make a summary table
    def table ( self , title = '' , prefix = '' , width = 5 , precision = 3 ) :
        """ Make a summary table
        """
        
        import ostap.logger.table  as     T 
                
        rows = [ ( 'Statistics' , 'value' , 'mean' , 'rms' , 'min/max' , 'factor' , 'p-value [%]' , '#sigma' ) ] 

        rows.append ( self._row ( 'Kolmogorov-Smirnov' , self.kolmogorov_smirnov , width = width , precision = precision ) )
        rows.append ( self._row ( 'Anderson-Darling'   , self.anderson_darling   , width = width , precision = precision ) )
        rows.append ( self._row ( 'Cramer-von Mises'   , self.cramer_von_mises   , width = width , precision = precision ) )
        rows.append ( self._row ( 'ZK'                 , self.ZK                 , width = width , precision = precision ) )
        rows.append ( self._row ( 'ZA'                 , self.ZA                 , width = width , precision = precision ) )
        rows.append ( self._row ( 'ZC'                 , self.ZC                 , width = width , precision = precision ) )

        ## skip empty columns        
        has_expo = False 
        for row in rows[1:] :
            r = list ( row )
            if r[-3] :
                has_expo = True
                break

        if not has_expo :
            new_rows = []
            for row in rows :
                r = list ( row )
                del r [ -3 ]
                new_rows.append ( r ) 
            rows = new_rows 

            
        if   not title and self.nToys :
            title = 'Goodness of 1D-fit with #%d toys' % self.nToys  
        elif not title :
            title = 'Goodness of 1D-fit'
        
        return T.table ( rows , title = title , prefix = prefix , alignment = 'lccccccccccc' )

    __repr__ = table
    __str__  = table
            
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================


