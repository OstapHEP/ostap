#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/twosamples.py
#  Two-Sample Tests for 1D non-weighted datasets 
#  @see https://www.jstor.org/stable/25471118
#  For dataset with weights, see 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2024-10-16
# =============================================================================
""" Two-Sample Tests for 1D non-weighted datasets 
 - see https://www.jstor.org/stable/25471118
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2024-09-16"
__all__     = (
    ##
    "TwoSamplesTest"     , ## Two-Samples Test
    "TwoSamplesToys"     , ## toys for Two-samples Test 
    ## 
)
# =============================================================================
from   ostap.core.core         import VE, Ostap
from   ostap.math.math_base    import axis_range, data2vct, pos_infinity
from   ostap.utils.cidict      import cidict_fun
from   ostap.logger.pretty     import pretty_float, format_pow10 
from   ostap.stats.gof_utils   import Labels, Keys, clip_pvalue 
from   ostap.stats.utils       import full_header, t_header  
from   ostap.math.math_ve      import significance
from   ostap.math.ve           import fmt_pretty_ve
from   ostap.stats.counters    import SE, EffCounter
from   ostap.stats.utils       import data2ecdf, weight_trivial 
from   ostap.stats.pvalue      import make_permutations
from   ostap.logger.symbols    import plus_minus, infinity_pos 
from   collections             import defaultdict, namedtuple
import ostap.stats.twosamples2 as TS2  
import ostap.logger.table      as     T
import ostap.math.models 
import ostap.fitting.ds2numpy 
import ostap.fitting.roofit
import ROOT, math, numpy    
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
logger = getLogger ( 'ostap.stats.twosamples' )
# =============================================================================
logger.debug ( 'Two-Samples Test for 1D non-weighted datasets' )
# =============================================================================
## @class TwoSamplesTest 
#  Two-Samples Test
class TwoSamplesTest ( object ) :
    """ Two-sample test 
    """
    def __init__ ( self        ,
                   data1       ,
                   data2       , 
                   weight1     = None ,
                   weight2     = None , 
                   frequencies = TS2.epps_singleton_frequencies () ) :
        
        ## transform input data 
        ecdf1 = data2ecdf ( data1 , weight1 )
        ecdf2 = data2ecdf ( data2 , weight2 )
        
        self.__frequencies = tuple  ( frequencies )
            
        self.__ecdf1   = ecdf1
        self.__ecdf2   = ecdf2
        self.__weight1 = weight1 
        self.__weight2 = weight2 
        self.__ecdf    = ecdf1.merge ( ecdf2 ) 
        
        self.__estimators = {            
            'KS'  : TS2.kolmogorov_smirnov ( self.ecdf1 , self.ecdf2 ) , 
            'K'   : TS2.kuiper             ( self.ecdf1 , self.ecdf2 ) , 
            'AD'  : TS2.anderson_darling   ( self.ecdf1 , self.ecdf2 ) , 
            'CM'  : TS2.cramer_von_mises   ( self.ecdf1 , self.ecdf2 ) , 
            'BJ'  : TS2.berk_jones         ( self.ecdf1 , self.ecdf2 ) , 
            'ZK'  : TS2.ZK                 ( self.ecdf1 , self.ecdf2 ) , 
            'ZA'  : TS2.ZA                 ( self.ecdf1 , self.ecdf2 ) , 
            'ZC'  : TS2.ZC                 ( self.ecdf1 , self.ecdf2 ) , 
            'W1'  : TS2.wasserstein1       ( self.ecdf1 , self.ecdf2 ) , 
            'W2'  : TS2.wasserstein2       ( self.ecdf1 , self.ecdf2 ) , 
            'R'   : TS2.rosenblatt         ( self.ecdf1 , self.ecdf2 ) , 
            'U2'  : TS2.watson_u2          ( self.ecdf1 , self.ecdf2 ) ,
            'AB'  : TS2.ansari_bradley     ( self.ecdf1 , self.ecdf2 ) , 
            'ST'  : TS2.siegel_tukey       ( self.ecdf1 , self.ecdf2 ) , 
            'MD'  : TS2.mood               ( self.ecdf1 , self.ecdf2 ) , 
            'VDW' : TS2.van_der_waerden    ( self.ecdf1 , self.ecdf2 ) , 
            'BM'  : TS2.brunner_munzel     ( self.ecdf1 , self.ecdf2 ) ,
            'ES'  : TS2.epps_singleton     ( self.ecdf1 , self.ecdf2 , frequencies = self.frequencies ) ,   
        }
        if  weight_trivial ( weight1 ) and weight_trivial ( weight2 ) :
            self.__estimators [ 'MW' ] =  TS2.mann_whitney ( self.ecdf1 , self.ecdf2 )
            self.__estimators [ 'S'  ] = TS2.savage        ( self.ecdf1 , self.ecdf2 )
        
    @property
    def ecdf1 ( self ) :
        """`ecdf1` : empirical CDF for the 1st dataset
        """
        return self.__ecdf1 
    @property
    def ecdf2 ( self ) :
        """`ecdf2` : empirical CDF for the 2nd dataset
        """
        return self.__ecdf2
    
    @property
    def ecdf ( self ) :
        """`ecdf` : total/combined CDF for both datasets
        """
        return self.__ecdf
    
    @property
    def weight1  ( self  ) :
        """`weight1`: weight for the 1st dataset"""
        return self.__weight1 
    
    @property
    def weight2  ( self  ) :
        """`weight2`: weight for the 2nd dataset"""
        return self.__weight2 
        
    @property
    def frequencies ( self ) :
        """`frequencies`: frequences for Epps-Singleton method"""
        return self.__frequencies 
    
    # =========================================================================
    ## Print the summary as Table
    def table ( self             , * , 
                title     = ''   ,
                prefix    = ''   ,
                precision = 4    , 
                width     = 6    ,
                style     = None ) :
        """ Print the summary as Table
        """
        ##
        header = ( 'Statistics' , )  + t_header 
        rows   = [] 
        
        for label , value  in self.estimators.items () :
            
            the_label = Labels.get ( label , label )
            
            result , expo = pretty_float ( value , width = width , precision = precision )            
            row = the_label , result , ( '[%s]' % format_pow10 ( expo ) )  if expo else ''
            
            rows.append ( row )
            
        rows  = [ header ] + sorted ( rows ) 
        rows  = T.remove_empty_columns ( rows )
        title = title if title else 'Two-Samples Test'    
        return T.table ( rows , title = title , prefix = prefix , alignment = 'lcl' , style = style  )

    ## print estimator as table 
    __repr__ = table
    __str__  = table
    
    # =========================================================================
    @property
    def estimators ( self ) :
        """`estimators` : get a dictionary of all estimators"""
        return self.__estimators

    # =========================================================================
    ## draw all involved empirical CDFS
    def draw  ( self , option = '' , *options , **kwargs ) :
        """ Draw fit CDF & empirical CDF
        """

        ecdf1 = self.ecdf1
        ecdf2 = self.ecdf2
        ecdf  = self.ecdf

        xmin  , xmax  = ecdf .xmin (), ecdf .xmax()
        xmin1 , xmax1 = ecdf1.xmin (), ecdf1.xmax()
        xmin2 , xmax2 = ecdf2.xmin (), ecdf2.xmax()

        xmin = min ( xmin , xmin1 , xmin2 )
        xmax = max ( xmax , xmax1 , xmax2 )

        xmin , xmax = axis_range ( xmin , xmax , delta = 0.10 )

        xmin = kwargs.pop ( 'xmin' , xmin )
        xmax = kwargs.pop ( 'xmax' , xmax )
        
        option  = option.strip() 
        optsame = 'same' if not option else '%s %s' % ( 'same' , option ) 
        r  = ecdf .draw ( option  , *options , color = 8 , linewidth = 3 , xmin = xmin , xmax = xmax , **kwargs )
        r1 = ecdf1.draw ( optsame , *options , color = 2 , linewidth = 3 , xmin = xmin , xmax = xmax , **kwargs )
        r2 = ecdf2.draw ( optsame , *options , color = 4 , linewidth = 3 , xmin = xmin , xmax = xmax , **kwargs )
        
        return r, r1, r2

# =============================================================================
## @class TwoSamplesToys
#  Toys for Two-Samples Test
class TwoSamplesToys(TwoSamplesTest):
    """ Toys for Two-sample test 
    """
    ## result of 2-samples test
    Result = namedtuple ( 'Result' , 'statistics counter pvalue nsigma' )
    # =========================================================================
    ## Initialize TSToys object 
    def __init__ ( self            ,
                   data1           ,
                   data2           ,
                   weight1 = None  ,
                   weight2 = None  , 
                   nToys   = 100   , 
                   silent  = False ) :
        
        ## initialize the base
        super().__init__ (  data1   = data1   , 
                            data2   = data2   , 
                            weight1 = weight1 , 
                            weight2 = weight2 ) 

        self.__counters = defaultdict(SE) 
        self.__ecdfs    = {} 
        self.__silent   = True if silent else False
        
        self.__counter  = EffCounter() 
        self.__counter += True 
        
        self.__nToys  = 0 
        if 0 < nToys : self.run ( nToys , silent = silent ) 
    
    # =========================================================================
    ## number of toys/permutations 
    @property
    def nToys ( self ) :
        """`nToys` : number of toys/permutations"""
        return self.__nToys
    
    # =========================================================================
    ## ECDFs
    @property 
    def ecdfs ( self ) :
        """`ecdfs` : toys/permutations results as empirical cumulative distribution functions"""
        return self.__ecdfs
    
    # =========================================================================
    ## Counters  
    @property 
    def counters ( self ) :
        """`counters` : toy results as counters"""
        return self.__counters
    
    # ===============================================================================
    ## run toys 
    def run ( self , nToys = 1000 , silent = False ) :
        """ Run toys 
        """ 
        assert isinstance ( nToys , int ) and 0 < nToys , "Invalid `nToys` argument!"

        n1   = len ( self.ecdf1 )
        n2   = len ( self.ecdf1 )
        
        data1 , weight1 = self.ecdf1.raw_data ()
        data2 , weight2 = self.ecdf2.raw_data ()
         
        results    = defaultdict(list)
        counters   = self.counters
        estimators = self.estimators 
        
        no_weight = weight_trivial ( weight1 ) and weight_trivial ( weight2 )
        
        ref_ks  = estimators.get ( 'KS'  , pos_infinity )
        ref_k   = estimators.get ( 'K'   , pos_infinity )
        ref_ad  = estimators.get ( 'AD'  , pos_infinity )
        ref_cm  = estimators.get ( 'CM'  , pos_infinity )
        ref_bj  = estimators.get ( 'BJ'  , pos_infinity )
        ref_zk  = estimators.get ( 'ZK'  , pos_infinity )
        ref_za  = estimators.get ( 'ZA'  , pos_infinity )
        ref_zc  = estimators.get ( 'ZC'  , pos_infinity )
        ref_w1  = estimators.get ( 'W1'  , pos_infinity )
        ref_w2  = estimators.get ( 'W2'  , pos_infinity )
        ref_r   = estimators.get ( 'R'   , pos_infinity )
        ref_u2  = estimators.get ( 'U2'  , pos_infinity )
        ref_ab  = estimators.get ( 'AB'  , pos_infinity )
        ref_st  = estimators.get ( 'ST'  , pos_infinity )
        ref_md  = estimators.get ( 'MD'  , pos_infinity ) 
        ref_vdw = estimators.get ( 'VDW' , pos_infinity )
        ref_bm  = estimators.get ( 'BM'  , pos_infinity )
        ref_es  = estimators.get ( 'ES'  , pos_infinity )
        ref_mw  = estimators.get ( 'MW'  , pos_infinity )
        ref_s   = estimators.get ( 'S'   , pos_infinity ) 
                             
                                
        for ds1, ds2, w1, w2 in make_permutations ( nToys , 
                                                      data1   = data1   , 
                                                      data2   = data1   , 
                                                      weight1 = weight1 , 
                                                      weight2 = weight2 ) :
        
            ## transform input data 
            ecdf1 = data2ecdf ( ds1 , w1  )
            ecdf2 = data2ecdf ( ds2 , w2  )
                
            ks  = TS2.kolmogorov_smirnov ( ecdf1 , ecdf2 )
            k   = TS2.kuiper             ( ecdf1 , ecdf2 )
            ad  = TS2.anderson_darling   ( ecdf1 , ecdf2 )            
            cm  = TS2.cramer_von_mises   ( ecdf1 , ecdf2 )             
            bj  = TS2.berk_jones         ( ecdf1 , ecdf2 )             
            zk  = TS2.ZK                 ( ecdf1 , ecdf2 )
            za  = TS2.ZA                 ( ecdf1 , ecdf2 )
            zc  = TS2.ZC                 ( ecdf1 , ecdf2 )
            w1  = TS2.wasserstein1       ( ecdf1 , ecdf2 )             
            w2  = TS2.wasserstein2       ( ecdf1 , ecdf2 )             
            r   = TS2.rosenblatt         ( ecdf1 , ecdf2 )             
            u2  = TS2.watson_u2          ( ecdf1 , ecdf2 )             
            ab  = TS2.ansari_bradley     ( ecdf1 , ecdf2 )
            st  = TS2.siegel_tukey       ( ecdf1 , ecdf2 )
            md  = TS2.mood               ( ecdf1 , ecdf2 )
            vdw = TS2.van_der_waerden    ( ecdf1 , ecdf2 )
            bm  = TS2.brunner_munzel     ( ecdf1 , ecdf2 )
            es  = TS2.epps_singleton     ( ecdf1 , ecdf2 , frequencies = self.frequencies )
            
            
            counters [ 'KS'  ] += ks 
            counters [ 'K'   ] += k 
            counters [ 'AD'  ] += ad
            counters [ 'CM'  ] += cm
            counters [ 'BJ'  ] += bj
            counters [ 'ZK'  ] += zk
            counters [ 'ZA'  ] += za
            counters [ 'ZC'  ] += zc
            counters [ 'W1'  ] += w1 
            counters [ 'W2'  ] += w2 
            counters [ 'R'   ] += r
            counters [ 'U2'  ] += u2
            counters [ 'AB'  ] += ab 
            counters [ 'ST'  ] += st  
            counters [ 'MD'  ] += md 
            counters [ 'VDW' ] += vdw 
            counters [ 'BM'  ] += bm 
            counters [ 'ES'  ] += es 
            
            results  [ 'KS'  ] .append ( ks  )    
            results  [ 'K'   ] .append ( k   ) 
            results  [ 'AD'  ] .append ( ad  ) 
            results  [ 'CM'  ] .append ( cm  ) 
            results  [ 'BJ'  ] .append ( bj  ) 
            results  [ 'ZK'  ] .append ( zk  ) 
            results  [ 'ZA'  ] .append ( za  ) 
            results  [ 'ZC'  ] .append ( zc  ) 
            results  [ 'W1'  ] .append ( w1  ) 
            results  [ 'W2'  ] .append ( w2  ) 
            results  [ 'R'   ] .append ( r   ) 
            results  [ 'U2'  ] .append ( u2  ) 
            results  [ 'AB'  ] .append ( ab  ) 
            results  [ 'ST'  ] .append ( st  ) 
            results  [ 'MD'  ] .append ( md  ) 
            results  [ 'VDW' ] .append ( vdw ) 
            results  [ 'BM'  ] .append ( bm  ) 
            results  [ 'ES'  ] .append ( es  ) 
            
            in_excess =   ( ref_ks  <= ks  or 
                            ref_k   <= k   or 
                            ref_ad  <= ad  or 
                            ref_cm  <= cm  or 
                            ref_bj  <= bj  or 
                            ref_zk  <= zk  or 
                            ref_za  <= za  or 
                            ref_zc  <= zc  or 
                            ref_w1  <= w1  or 
                            ref_w2  <= w2  or  
                            ref_r   <= r   or
                            ref_u2  <= u2  or  
                            ref_ab  <= ab  or 
                            ref_st  <= st  or 
                            ref_md  <= md  or 
                            ref_vdw <= vdw or 
                            ref_bm  <= bm  or 
                            ref_es  <= es  ) 
                                                                                                            
            if no_weight : 
                mw = TS2.mann_whitney  ( ecdf1 , ecdf2 )             
                s  = TS2.savage        ( ecdf1 , ecdf2 )     
                counters [ 'MW'  ] += mw
                counters [ 'S'   ] += s 
                results  [ 'MW'  ] .append ( mw ) 
                results  [ 'S'   ] .append ( s  ) 
                in_excess = in_excess or ref_mw <= mw or ref_s <= s
                
            self.__counter += in_excess 
            self.__nToys   += 1 
            
        ## accumulate number of to

        ECDF = Ostap.Math.ECDF 
        for key in results :
            data = results [ key ]
            if not data : continue
            if not key in self.__ecdfs : self.__ecdfs [ key ] = ECDF ( data , True ) ## complementary ECDF!
            else                       : self.__ecdfs [ key ]  .add  ( data2vct ( data ) ) 
            
        del results

    # ========================================================================
    ## Helper method to get the result
    def result ( self , label ) :
        """ Helper method to get the result 
        """
        if not label in self.estimators : return None
        if not label in self.ecdfs      : return None
        if not label in self.counters   : return None
        ##
        value   = self.estimators   [ label ]
        ecdf    = self.ecdfs        [ label ] 
        counter = self.counters     [ label ] 
        ##
        pvalue = ecdf. estimate ( value  ) ## estimate the p-value
        #
        pv = clip_pvalue ( pvalue , 0.5 ) 
        nsigma = significance ( pv ) ## convert  it to significace
        
        return self.Result ( value   ,
                             counter ,
                             pvalue  ,
                             nsigma  )
    
    # =========================================================================
    ## format a row in the summary table
    def row  ( self , what , result , width = 6 , precision = 4 ) :
        """ Format a row in the sumamry table
        """
        value      = result.statistics
        counter    = result.counter
        pvalue     = result.pvalue
        nsigma     = result.nsigma
        
        mean       = counter.mean   ()
        rms        = counter.rms    () 
        vmin, vmax = counter.minmax () 
        
        mxv = max ( abs ( value        ) ,
                    abs ( mean.value() ) ,
                    mean.error()         , rms ,
                    abs ( vmin )  , abs ( vmax ) ) 
        
        fmt, fmtv , fmte , expo = fmt_pretty_ve ( VE ( mxv ,  mean.cov2() ) ,
                                                  width       = width       ,
                                                  precision   = precision   , 
                                                  parentheses = False       )
        
        if expo : scale = 10**expo
        else    : scale = 1
        
        fmt2 = '%s/%s' % ( fmtv , fmtv ) 

        vs  = value / scale
        vm  = mean  / scale
        vr  = rms   / scale
        vmn = vmin  / scale
        vmx = vmax  / scale
        
        pvalue = str ( ( 100 * pvalue ) .toString ( '%% 5.2f %s %%-.2f' % plus_minus ) )
        nsigma = str ( nsigma.toString ( '%%.2f %s %%-.2f' % plus_minus ) if float ( nsigma ) < 1000 else infinity_pos  ) 

        return ( what  ,
                 fmtv  % vs ,
                 fmt   % ( vm.value() , vm.error() ) ,
                 fmtv  % vr                          ,
                 fmtv  % vmn                         ,
                 fmtv  % vmx                         ,
                 ( '[%s]' % format_pow10 ( expo ) )  if expo else '' , pvalue , nsigma )
    
    # =========================================================================
    ## Make a summary table
    def table ( self , title = '' , prefix = '' , width = 6 , precision = 4 , style = None ) :
        """ Make a summary table
        """
        import ostap.logger.table  as     T                 
        header = ( 'Statistics'  , ) + full_header
        rows   = [] 
        
        for label in self.ecdfs :
            
            result  = self.result ( label )
            if not result : continue

            the_label = Labels.get ( label , label )
            row = self.row ( the_label , result , width = width , precision = precision )
            rows.append ( row ) 
                    
        if   not title and self.nToys : title = 'Goodness of 1D-fit with #%d toys' % self.nToys  
        elif not title                : title = 'Goodness of 1D-fit'

        rows = [ header ] + sorted ( rows ) 
        rows = T.remove_empty_columns ( rows ) 
        return T.table ( rows , title = title , prefix = prefix , alignment = 'lccccccc' , style = style )

    __repr__ = table
    __str__  = table

    # =========================================================================
    ## Draw ECDF for toys & statistical estimator 
    def draw  ( self , what , option = '' , *options , **kwargs ) :
        """ Draw ECDF for toys & statistical estgimator 
        """
        key = cidict_fun ( what ) 
        if   key in Keys [ 'KS' ] and 'KS' in self.ecdfs :             
            result = self.result ( 'KS' )
            ecdf   = self.ecdfs  [ 'KS' ]
        elif key in Keys [ 'K'  ] and 'K'  in self.ecdfs : 
            result = self.result ( 'K' )
            ecdf   = self.ecdfs  [ 'K' ]
        elif key in Keys [ 'AD' ] and 'AD' in self.ecdfs :             
            result = self.result ( 'AD' )
            ecdf   = self.ecdfs  [ 'AD' ]
        elif key in Keys [ 'CM' ] and 'CM' in self.ecdfs : 
            result = self.result  ( 'CM' )
            ecdf   = self.ecdfs   [ 'CM' ]
        elif key in Keys [ 'BJ' ] and 'BJ' in self.ecdfs : 
            result = self.result  ( 'BJ' )
            ecdf   = self.ecdfs   [ 'BJ' ]            
        elif key in Keys [ 'ZK' ]  and 'ZK' in self.ecdfs : 
            result = self.result  ( 'ZK' )
            ecdf   = self.ecdfs   [ 'ZK' ]
        elif key in Keys [ 'ZA' ]   and 'ZA' in self.ecdfs :  
            result = self.result  ( 'ZA' )
            ecdf   = self.ecdfs   [ 'ZA' ]
        elif key in Keys [ 'ZC' ] and 'ZC' in self.ecdfs : 
            result = self.result  ( 'ZC' )
            ecdf   = self.ecdfs   [ 'ZC' ]
        elif key in Keys [ 'W1' ] and 'W1' in self.ecdfs : 
            result = self.result  ( 'W1' )
            ecdf   = self.ecdfs   [ 'W1' ]
        elif key in Keys [ 'W2' ] and 'W2' in self.ecdfs : 
            result = self.result  ( 'W2' )
            ecdf   = self.ecdfs   [ 'W2' ]
        elif key in Keys [ 'R' ] and 'R' in self.ecdfs : 
            result = self.result  ( 'R' )
            ecdf   = self.ecdfs   [ 'R' ]
        elif key in Keys [ 'U2' ] and 'U2' in self.ecdfs : 
            result = self.result  ( 'U2' )
            ecdf   = self.ecdfs   [ 'U2' ]
        elif key in Keys [ 'AB' ] and 'AB' in self.ecdfs : 
            result = self.result  ( 'AB' )
            ecdf   = self.ecdfs   [ 'AB' ]
        elif key in Keys [ 'ST' ] and 'ST' in self.ecdfs : 
            result = self.result  ( 'ST' )
            ecdf   = self.ecdfs   [ 'ST' ]
        elif key in Keys [ 'MD' ] and 'MD' in self.ecdfs : 
            result = self.result  ( 'MD' )
            ecdf   = self.ecdfs   [ 'MD' ] 
        elif key in Keys [ 'VDW' ] and 'VDW' in self.ecdfs : 
            result = self.result  ( 'VDW' )
            ecdf   = self.ecdfs   [ 'VDW' ]    
        elif key in Keys [ 'BM' ] and 'BM' in self.ecdfs : 
             result = self.result  ( 'BM' )
             ecdf   = self.ecdfs   [ 'BM' ] 
        elif key in Keys [ 'ES' ] and 'ES' in self.ecdfs : 
             result = self.result  ( 'ES' )
             ecdf   = self.ecdfs   [ 'ES' ]    
        elif key in Keys [ 'MW' ] and 'MW' in self.ecdfs : 
            result = self.result  ( 'MW' )
            ecdf   = self.ecdfs   [ 'MW' ]
        elif key in Keys [ 'S' ] and 'S' in self.ecdfs : 
            result = self.result  ( 'S' )
            ecdf   = self.ecdfs   [ 'S' ]
        else :
            raise KeyError (  "draw: Invalid `what`:%s" % what )

        xmin , xmax = ecdf.low_edge () , ecdf.high_edge ()
        value       = result.statistics
        xmin        = min ( xmin , value )
        xmax        = max ( xmax , value )
        xmin , xmax = axis_range ( xmin , xmax , delta = 0.20 )

        kwargs [ 'xmin' ] = kwargs.get ( 'xmin' , xmin ) 
        kwargs [ 'xmax' ] = kwargs.get ( 'xmax' , xmax )

        result    = ecdf.draw  ( option , *options , **kwargs ) 
        line      = ROOT.TLine ( value  , 1e-3 , value , 1 - 1e-3 )
        ## 
        line.SetLineWidth ( 4 ) 
        line.SetLineColor ( 8 ) 
        line.draw ( 'same' )
        ##
        self._line = line
        return result, line
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================




