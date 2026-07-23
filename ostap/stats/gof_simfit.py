#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/gof_simfit.py
#  Set of utulities for goodness-of-fit studies for SimFit 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2024-09-16
# =============================================================================
""" Simple utulities for goodness-of-fit studies for Simfit 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2024-09-16"
__all__     = (
    'GoFSimFit1D'      , ## helper utility for GoF estimate 
    'GoFSimFit1DToys'  , ## helper utility for GoF estimate with toys
    'GoFSimFit'        , ## helper utility for generic GoF estimate 
    'GoFSimFitType'    , ## helper utility for generic GoF estimate 
)
# =============================================================================
from   ostap.core.ostap_types   import string_types, dictlike_types 
from   ostap.fitting.pdfbasic   import PDF1, APDF1
from   ostap.core.core          import VE, Ostap
from   ostap.math.math_base     import axis_range, product 
from   ostap.utils.cidict       import cidict_fun
from   ostap.utils.basic        import loop_items, typename   
from   ostap.stats.counters     import SE, EffCounter 
from   ostap.logger.pretty      import pretty_float
from   ostap.math.ve            import fmt_pretty_ve
from   ostap.math.math_ve       import significance
from   ostap.logger.symbols     import plus_minus, times, greek_lower_sigma, clock, toys 
from   ostap.logger.colorized   import infostr
from   ostap.stats.gof_utils    import Labels, Keys, clip_pvalue, combine_pvalues 
from   ostap.fitting.simfit     import SimFit 
from   ostap.stats.gof          import AGoF
from   ostap.stats.gof1d        import ( GoF1D , vct_clip      ,
                                         kolmogorov_smirnov    ,
                                         anderson_darling      ,
                                         cramer_von_mises      ,
                                         kuiper , berk_jones   ,
                                         ZK , ZA , ZC          )
from   ostap.plotting.color     import Navy, DarkGreen
from   ostap.utils.timing       import timing 
from   collections              import defaultdict, namedtuple
import ostap.stats.gofnd        as     GoFnD 
import ostap.logger.table       as     T
import ostap.fitting.roofit
import ROOT, math, numpy, copy  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.gof_simfit' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Simple utilities for goodness-of-SimFit studies' )
# =============================================================================
class GoFSimFitBase(object) :
    """ Base class for Goodness-of-fit for Simultaneous Fits
    """
    def __init__ ( self               ,
                   pdf                ,
                   dataset            ,
                   parameters  = None ) :
        
        assert isinstance ( pdf     , SimFit          ) , 'Invalid type of `pdf`:%s'     % typename ( pdf     )
        assert isinstance ( dataset , ROOT.RooDataSet ) , 'Invalid type of `dataset`:%s' % typename ( dataset )

        assert pdf.sample in dataset , "Sample category `%s` not in dataset!" % pdf.sample.name 

        self.__pdf        = pdf
        self.__parameters = parameters if parameters else pdf.parameters ( dataset )
        
        ## load parameters here
        self.__pdf.load_params ( self.__parameters , silent = True ) 
        
        self.__gofs = {}
        self.__N    = {}

    @property
    def pdf ( self ) :
        """`pdf`: SimFit/PDF for simultaneous fit
        """
        return self.__pdf
    
    @property
    def sample ( self ) :
        """`sample`: sample/category  variable for simultaneous fit
        """
        return self.pdf.sample 

    @property
    def parameters ( self ) :
        """`parameters' : fit parameters, e.g. fit-result or dictionary or ...
        """
        return self.__parameters 
        
    @property
    def gofs ( self ) :
        """`gofs` : individual GoF estimators for Simfit components
        """
        return self.__gofs 

    @property
    def N ( self ) :
        """`N` : dictionary { sample : #events}
        """
        return self.__N 

    # =============================================================================
    ## serialize the object 
    def __getstate__ ( self ) :
        """ Serialize the object
        """
        if self.parameters : self.pdf.load_params ( self.parameters , silent = True )
        return { 'pdf'         : self.pdf        ,
                 'parameters'  : self.parameters , 
                 'gofs'        : self.gofs       ,
                 'N'           : self.N          }
    # ============================================================================
    ## De-serialize the object 
    def __setstate__ ( self , state ) :
        """ De-serialize the object\
        """         
        self.__pdf         = state.pop ( 'pdf'  )
        self.__parameters  = state.pop ( 'parameters' , {} ) 
        self.__gofs        = state.pop ( 'gofs'       , {} )
        self.__N           = state.pop ( 'N'          , {} )
        ## 
        if self.parameters : self.pdf.load_params ( self.parameters , silent = True )
        
# =============================================================================
## @class GoFSimFit1D 
#  Goodness-of-fit for simultaneous 1D-fits
#  - All components of the simultaneous fit must be 1D-components 
#  - GoF is estimated for each 1D-component
class GoFSimFit1D(GoFSimFitBase) :
    """ Goodness-of-fit for Simultaneous 1D-fits
    - All components of the Simultaneous fit must be 1D-components!
    - GoF is estimated for each 1D-component
    """
    def __init__ ( self               ,
                   pdf                ,
                   dataset            ,
                   parameters  = None ) : 
        
        ## initialize the base class 
        GoFSimFitBase.__init__ ( self                    ,
                                 pdf        = pdf        ,
                                 dataset    = dataset    ,
                                 parameters = parameters )

        name = self.sample.name
        for key , cmp  in self.pdf.categories.items ()  :
            
            assert isinstance ( cmp , PDF1 ) , "Component `%s` is not PDF1` %s" % ( key , typename ( cmp ) )
            obs      = cmp.pdf.getObservables ( dataset )
            category = '%s==%s::%s' % ( name , name , key )
             
            ds       = dataset.subset ( variables = obs , cuts = category )
            
            gof      = GoF1D ( cmp , ds )
            
            self.gofs [ key ] = gof 
            self.N    [ key ] = len ( ds ) 
            
    ## serialize the object 
    def __getstate__ ( self ) :
        """ Serialize the object 
        """
        ## (1) serialize the base 
        return GoFSimFitBase.__getstate__ ( self )
    
    ## De-serialize the object 
    def __setstate__ ( self , state ) :
        """ De-serialize the object """        
        ## (1) de-serialize the base 
        return GoFSimFitBase.__setstate__ ( self , state )

    # =========================================================================
    ## all estimators together
    @property 
    def estimators ( self ) :
        """`estimators` : get all statistical estimators
        """
        return { key : v.estimators for key , v in self.gofs.items () }
        
    @property
    def kolmogorov_smirnov_estimator ( self ) :
        """ Get Kolmogorov-Smirnov statistics' KS
        """
        return { key : v.kolmogorov_smirnov_estimator for key , v in self.gofs.items () }
    
    @property
    def kuiper_estimator ( self ) :
        """ Get Kuiper' statistics K
        """
        return { key : v.kuiper_estimator for key , v in self.gofs.items () }

    @property
    def anderson_darling_estimator ( self ) :
        """ Get Anderson-Darling' statistics AD
        """
        return { key : v.anderson_darling_estimator for key , v in self.gofs.items () }
    
    @property
    def cramer_von_mises_estimator ( self ) :
        """ Get Cramer-on Mises' statistics CM
        """
        return { key : v.cramer_von_mises_estimator for key , v in self.gofs.items () }
    
    @property
    def ZK_estimator ( self ) :
        """ Get Zhang's statistics ZK
        """
        return { key : v.ZK_estimator for key , v in self.gofs.items () }
    
    @property
    def ZA_estimator ( self ) :
        """ Get Zhang's statistics ZA
        """
        return { key : v.ZA_estimator for key , v in self.gofs.items () }
    
    @property
    def ZC_estimator ( self ) :
        """ Get Zhang's statistics ZC
        """
        return { key : v.ZC_estimator for key , v in self.gofs.items () }

    @property
    def berk_jones_estimator ( self ) :
        """ Get Berk-Jones 's statistics BJ
        """
        return { key : v.berk_jones_estimator for key , v in self.gofs.items () }
    
    # ====================================================================================
    ## Print the summary as Table  (for simfit)
    def table ( self             , * ,
                title     = ''   ,
                prefix    = ''   ,
                precision = 4    , 
                width     = 6    ,
                style     = None ) :
        """ Print the summary as Table  (for SimFit)
        """
        ##
        
        keys = tuple ( self.gofs.keys() )
        
        header = "Statistics",  
        for k in keys : header += ( str ( k ) , '' )
        
        rows = []

        ## get the estimators from the fist GOF 
        estimators = self.gofs [ keys [ 0 ] ].estimators.keys()
        for label in estimators : 
            
            the_label = Labels.get ( label , label )

            if   'KS' == label : values = self.kolmogorov_smirnov_estimator
            elif 'K'  == label : values = self.kuiper_estimator
            elif 'AD' == label : values = self.anderson_darling_estimator
            elif 'CM' == label : values = self.cramer_von_mises_estimator
            elif 'ZK' == label : values = self.ZK_estimator
            elif 'ZA' == label : values = self.ZA_estimator
            elif 'ZC' == label : values = self.ZC_estimator
            elif 'BJ' == label : values = self.berk_jones_estimator
            else               : continue
            
            row = the_label ,
            
            for k, v in values.items () :
                result , expo = pretty_float ( v  , width = width , precision = precision )                
                if expo : row += ( result , '10^%+d' % expo )
                else    : row += ( result , ''              )
                
            rows.append ( row )

        rows  = [ header ] + sorted ( rows ) 
        title = title if title else 'Goodness of 1D (Sim)Fit'
        rows  = T.remove_empty_columns ( rows )
        return T.table ( rows , title = title , prefix = prefix , alignment = 'lccccccccccccccccccc' , style = style  )

    ## print estimators as table 
    __repr__ = table
    __str__  = table
    report   = table
    
    ## Draw fit CDF & empirical ECDF 
    def draw  ( self , sample , option = '' , **kwargs ) :
        """ Draw fit CDF & empirical CDF
        """
        gof = self.gofs.get ( sample , None )
        if gof is None : raise KeyError ( "Invalid sample `%s`" % sample )
        return gof.draw ( option , **kwargs )

    
# =============================================================================
## @class GoFSimFitToys
#  Check Goodness of 1D (Sim)Fits using toys 
class GoFSimFit1DToys(GoFSimFit1D) :
    """ Check Goodness-of-Fit with toys (SimFit case)
    """
    ## result of GoF-toys 
    Result = namedtuple ( 'Result' , 'statistics counter pvalue nsigma' )
    # =========================================================================
    ## Initialize GoF1D toys object :
    #  @code
    #  gof  = GoFSimfit     ( ... ) 
    #  toys = GoFSimFitToys ( gof )
    #  @endcode 
    def __init__ ( self , gof ) :
        """ Initialize GoF1D toys object :
        >>> gof  = GoFSimFit     ( ... ) 
        >>> toys = GoFSimFitToys ( gof ) 
        """
        assert isinstance ( gof , GoFSimFit1D ) , \
            "Invalid `gof`-parameter: %s" % typename ( gof ) 

        ## mimic the copy-constructor for the base class 
        state = GoFSimFit1D.__getstate__ ( gof ) 
        GoFSimFit1D.__setstate__ ( self , state )

        self.__counters = { k : defaultdict(SE) for k in self.gofs }
        self.__ecdfs    = { k : {}              for k in self.gofs }
        self.__total    = defaultdict(EffCounter) 
        self.__nToys    = 0

    ## serialize the object 
    def __getstate__ ( self ) :
        """ Serialize the object 
        """
        ## (1) serialize the base 
        state = GoFSimFit1D.__getstate__ ( self )
        # 
        state [ 'counters' ] = self.__counters
        state [ 'ecdfs'    ] = self.__ecdfs 
        state [ 'total'    ] = self.__total 
        state [ 'nToys'    ] = self.__nToys
        # 
        return state 
    
    ## De-serialize the object 
    def __setstate__ ( self , state ) :
        """ De-serialize the object 
        """        
        ## (1) de-serialize the base 
        GoFSimFit1D.__setstate__ ( self , state )
        # 
        self.__counters   = state.pop ( 'counters'  )
        self.__ecdfs      = state.pop ( 'ecdfs'     )
        self.__total      = state.pop ( 'total'     )
        self.__nToys      = state.pop ( 'nToys' , 0 )

    # ===============================================================================
    ## run toys 
    def run ( self ,
              nToys    = 1000  , * ,
              parallel = False ,
              fitconf  = {}    , 
              silent   = False ,
              nSplit   = 0     ) :
        """ Run toys 
        """
        assert isinstance ( nToys , int ) and 1 <= nToys , "Invalid `nToys` argument!"

        if parallel :
            
            from ostap.parallel.parallel_gof import parallel_goftoys as parallel_toys 
            result = parallel_toys ( gof      = self       ,
                                     nToys    = nToys      ,
                                     nSplit   = nSplit     ,
                                     fitconf  = fitconf    , 
                                     silent   = True       ,
                                     progress = not silent )
            if result : self+= result  
            return self 
        
        results  = { k : defaultdict(list) for k in self.gofs } 
        counters = self.counters 

        from ostap.utils.progress_bar import progress_bar
        for i in progress_bar ( nToys , silent = silent , description = 'Toys:') :

            self.pdf.load_params ( self.parameters , silent = True )            
            dset     = self.pdf.generate ( self.N  , sample = True )
            sample   = self.pdf.sample

            total_KS = True  
            total_K  = True 
            total_AD = True 
            total_CM = True 
            total_ZK = True 
            total_ZA = True 
            total_ZC = True 
            total_BJ = True 
                      
            for key , gof  in self.gofs.items ()  :

                obs      = gof.pdf.pdf.getObservables ( dset )
                category = '%s==%s::%s' % ( sample.name , sample.name , key )
                ds       = dset.subset ( variables = obs , cuts = category )
                varname  = gof.pdf.xvar.name

                ## convert to numpy & sort it! 
                data     = ds.tonumpy ( varname ) [ varname ]
                data     = numpy.sort ( data )
                
                vct_cdf  = gof.vcdf

                ## CLIP... does one need it? 
                cdf_data = vct_clip ( vct_cdf ( data ) )
            
                ks       = kolmogorov_smirnov ( cdf_data )
                k        = kuiper             ( cdf_data )
                ad       = anderson_darling   ( cdf_data )
                cm       = cramer_von_mises   ( cdf_data )
                zk       = ZK                 ( cdf_data )
                za       = ZA                 ( cdf_data )
                zc       = ZC                 ( cdf_data )
                bj       = berk_jones         ( cdf_data )
                
                total_KS = total_KS and gof .kolmogorov_smirnov_estimator <= ks
                total_K  = total_K  and gof .            kuiper_estimator <= k
                total_AD = total_AD and gof .  anderson_darling_estimator <= ad
                total_CM = total_CM and gof .  cramer_von_mises_estimator <= cm
                total_ZK = total_ZK and gof .                ZK_estimator <= zk 
                total_ZA = total_ZA and gof .                ZA_estimator <= za 
                total_ZC = total_ZC and gof .                ZC_estimator <= zc 
                total_BJ = total_BJ and gof .        berk_jones_estimator <= bj
                
                cnts = counters [ key ]                
                cnts [ 'KS' ] += ks
                cnts [ 'K'  ] += k
                cnts [ 'AD' ] += ad
                cnts [ 'CM' ] += cm
                cnts [ 'ZK' ] += zk
                cnts [ 'ZA' ] += za
                cnts [ 'ZC' ] += zc
                cnts [ 'BJ' ] += bj

                res  = results [ key ]
                res [ 'KS'  ].append ( ks )    
                res [ 'K'   ].append ( k  )
                res [ 'AD'  ].append ( ad ) 
                res [ 'CM'  ].append ( cm ) 
                res [ 'ZK'  ].append ( zk ) 
                res [ 'ZA'  ].append ( za ) 
                res [ 'ZC'  ].append ( zc ) 
                res [ 'BJ'  ].append ( bj ) 

                ## delete data
                if isinstance ( ds , ROOT.RooDataSet ) : ds.clear () 
                del ds
                del data
                del cdf_data
                
            self.__total [ 'KS' ] += total_KS
            self.__total [ 'K'  ] += total_K
            self.__total [ 'AD' ] += total_AD
            self.__total [ 'CM' ] += total_CM
            self.__total [ 'ZK' ] += total_ZK
            self.__total [ 'ZA' ] += total_ZA
            self.__total [ 'ZC' ] += total_ZC
            self.__total [ 'BJ' ] += total_BJ
                
            ## delete data
            if isinstance ( dset , ROOT.RooDataSet ) : dset.clear () 
            del dset
            
        ## accumulate number of toys 
        self.__nToys += nToys 

        ECDF = Ostap.Math.ECDF 
        for key, vv  in results.items()  :
            for e , data in vv.items() : 
                if not data : continue
                ecdfs = self.__ecdfs [ key ]
                if not e in ecdfs : ecdfs [ e ] = ECDF ( data , True ) ## complementary ECDF!
                else              : ecdfs [ e ] .add   ( data2vct ( data ) ) 
                
        del results 
        return self
    
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
    ## Counters  
    @property 
    def counters ( self ) :
        """`counters` : toy results as counters"""
        return self.__counters

    # =========================================================================
    ## Total/global counters
    @property 
    def total ( self ) :
        """`total` : total/global counters
        """
        return self.__total
    
    # =========================================================================
    ## Get Kolmogorov-Smirnov statistics 
    @property 
    def kolmogorov_smirnov ( self ) :
        """ Get Kolmogorov-Smirnov statistiscs KS
        """
        return self.result ( 'KS' ) 
    
    # ===============================================
    ## Get Anderson-Darling  statistics
    @property 
    def anderson_darling ( self ) :
        """ Get Anderson-Darling statistiscs 
        """
        return self.result ( 'AD' ) 
        
    # =========================================================================
    ## Get Cramer-von Mises statistics 
    @property 
    def cramer_von_mises ( self ) :
        """ Get Cramer-von Mises statistics 
        """
        return self.result ( 'CM' ) 
        
    # =========================================================================
    ## Get ZK statistics
    @property 
    def ZK  ( self ) :
        """ Get ZK statistics 
        """
        return self.result ( 'ZK' ) 
        
    # =========================================================================
    ## Get ZA statistics
    @property 
    def ZA  ( self ) :
        """ Get ZA statistics
        """        
        return self.result ( 'ZA' ) 
    
    # =========================================================================
    ## Get ZC statistics
    @property 
    def ZC  ( self ) :
        """ Get ZC statistics 
        """        
        return self.result ( 'ZC' ) 

    # =========================================================================
    ## Get Kuiper statistics
    @property 
    def kuiper ( self ) :
        """ Get Kuiper statistics 
        """        
        return self.result ( 'K' ) 

    # =========================================================================
    ## Get Berk-Jonesstatistics
    @property 
    def berk_jones ( self ) :
        """ Get Berk-Jones statistics 
        """        
        return self.result ( 'BJ' ) 
    
    # ========================================================================
    ## Helper method to get the result
    def result ( self , sample , label ) :
        """ Helper method to get the result 
        """
        ##
        value   = self.estimators   [ sample ] [ label ]
        ecdf    = self.ecdfs        [ sample ] [ label ] 
        counter = self.counters     [ sample ] [ label ] 
        ##
        pvalue = ecdf. estimate ( value  ) ## estimate the p-value
        #
        pv     = clip_pvalue  ( pvalue , 0.5 )
        nsigma = significance ( pv ) ## convert  it to significace
        
        return self.Result ( value   ,
                             counter ,
                             pvalue  ,
                             nsigma  )

    # =========================================================================
    ## format a row in the summary table
    def row  ( self , sample , what , result , width = 6 , precision = 4 ) :
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
        sigma  = str ( nsigma.toString ( '%%.2f %s %%-.2f' % plus_minus ) if float ( nsigma ) < 1000 else '+inf' ) 

        return ( what ,
                 sample ,
                 fmtv   % vs ,
                 fmt    % ( vm.value() , vm.error() ) ,
                 fmte   % vr                          ,
                 fmt2   %  ( vmn , vmx )              ,
                 ( '%s10^%+d' %  ( times , expo )  if expo else '' )   ,                  
                 pvalue , sigma )
    
    # =========================================================================
    ## Make a summary table
    def table ( self , title = '' , prefix = '' , width = 6 , precision = 4 , style = None ) :
        """ Make a summary table
        """
        import ostap.logger.table  as     T                 
        header = ( 'Statistics'        ,
                   'Sample'            , 
                   'value'             , 
                   'mean'              ,
                   'rms'               ,
                   'min/max'           ,
                   'factor'            ,
                   'p-value [%]'       ,
                   '#%s' % greek_lower_sigma ) 
        
        rows = []

        ## collect p-values to calculate the combined p-value 
        to_combine = defaultdict(list)
        
        for sample , ecdfs in self.ecdfs.items()  :
            for label in ecdfs :
                result  = self.result ( sample , label )
                if not result : continue
                the_label = Labels.get ( label , label )
                row = self.row ( sample , the_label , result , width = width , precision = precision )
                rows.append ( row ) 
                to_combine[ the_label ].append ( result.pvalue ) 
                
        if   not title and self.nToys : title = 'Goodness of 1D-fit with #%d toys' % self.nToys  
        elif not title                : title = 'Goodness of 1D-fit'

        # ======================================================================
        ## get the total from the toys 
        for e , cnt in self.total.items() :

            ## get the binomial efficiency 
            pvalue = cnt.efficiency

            pv     = clip_pvalue ( pvalue , 0.5 )
            nsigma = significance ( pv ) ## convert  it to significace
            
            pvalue = str ( ( 100 * pvalue ) .toString ( '%% 5.2f %s %%-.2f' % plus_minus ) )
            sigma  = str ( nsigma.toString ( '%%.2f %s %%-.2f' % plus_minus ) if float ( nsigma ) < 1000 else '+inf' ) 

            label = Labels.get( e , e ) 
            row   = label , '*COMBINED*' , '' , '' , '' , '' , ''  , pvalue , sigma 
            rows.append ( row ) 

        for label , pvs in to_combine.items () :
            
            clipped  = [ clip_pvalue ( p , 0.5 ) for p in pvs ]                 
            for method in  ( 'fisher' , 'pearson' , 'tippett' , 'stouffer' ) :
                
                combined = combine_pvalues ( clipped , method = method )
                
                pvalue   = combined 
                pv       = clip_pvalue ( pvalue , 0.5 )
                
                nsigma  = significance ( pv ) ## convert  it to significace            
                pvalue  = str ( ( 100 * pvalue ) .toString ( '%% 5.2f %s %%-.2f' % plus_minus ) )
                sigma   = str ( nsigma.toString ( '%%.2f %s %%-.2f' % plus_minus ) if float ( nsigma ) < 1000 else '+inf' ) 

                mm      = method.upper() 
                row     = label , '*%s*' % mm  , '' , '' , '' , '' , ''  , pvalue , sigma 
                rows.append ( row )
        
            minimal = min ( clipped , key = lambda p : float ( p ) )
            pvalue  = minimal  
            pv       = clip_pvalue ( pvalue , 0.5 )
                
            nsigma  = significance ( pv ) ## convert  it to significace            
            pvalue  = str ( ( 100 * pvalue ) .toString ( '%% 5.2f %s %%-.2f' % plus_minus ) )
            sigma   = str ( nsigma.toString ( '%%.2f %s %%-.2f' % plus_minus ) if float ( nsigma ) < 1000 else '+inf' ) 

            mm      = 'MINIMAL' 
            row     = label , '*%s*' % mm  , '' , '' , '' , '', '' , pvalue , sigma 
            rows.append ( row )
            
        rows  = [ header ] + sorted ( rows )
        rows  = T.remove_empty_columns ( rows ) 
        return T.table ( rows , title = title , prefix = prefix , alignment = 'lcccccccc' , style = style )
        
    __repr__ = table
    __str__  = table
    report   = table
    
    # =========================================================================
    ## Draw ECDF for toys & statistical estimator 
    def draw  ( self , sample , what , option = '' , **kwargs ) :
        """ Draw ECDF for toys & statistical estgimator 
        """
        if not sample in self.ecdfs :
            raise KeyError (  "draw: Invalid `sample`:%s" % sample  )

        ecdfs = self.ecdfs[ sample ]

        key = cidict_fun ( what )
        if   key in Keys [ 'KS' ] and 'KS' in ecdfs :             
            result = self.result ( sample , 'KS' )
            ecdf   = ecdfs [ 'KS' ]
            ## logger.info ( 'Toy results for %s/Kolmogorov-Smirnov estimate' % sample ) 
        elif key in Keys [ 'K'  ] and 'K'  in ecdfs : 
            result = self.result ( sample , 'K' )
            ecdf   = ecdfs [ 'K' ]
            ## logger.info ( 'Toy results for %s/Kuiper estimate'             % sample ) 
        elif key in Keys [ 'AD' ] and 'AD' in ecdfs :             
            result = self.result ( sample , 'AD' )
            ecdf   = ecdfs  [ 'AD' ]
            ## logger.info ( 'Toy results for %s/Anderson-Darling estimate'   % sample ) 
        elif key in Keys [ 'CM' ] and 'CM' in ecdfs : 
            result = self.result  ( sample , 'CM' )
            ecdf   = ecdfs   [ 'CM' ]
            ## logger.info ( 'Toy results for %s/Cramer-von Mises  estimate'  % sample ) 
        elif key in Keys [ 'ZK' ]  and 'ZK' in ecdfs : 
            result = self.result  ( sample , 'ZK' )
            ecdf   = ecdfs   [ 'ZK' ]
            ## logger.info ( 'Toy results for %s/Zhang/ZK estimate'           % sample ) 
        elif key in Keys [ 'ZA' ]   and 'ZA' in ecdfs :  
            result = self.result  ( sample , 'ZA' )
            ecdf   = ecdfs   [ 'ZA' ]
            ## logger.info ( 'Toy results for %s/Zhang/ZA estimate'           % sample ) 
        elif key in Keys [ 'ZC' ] and 'ZC' in ecdfs : 
            result = self.result  ( sample , 'ZC' )
            ecdf   = ecdfs   [ 'ZC' ]
            ## logger.info ( 'Toy results for %s/Zhang/ZC estimate'           % sample ) 
        elif key in Keys [ 'BJ' ] and 'BJ' in ecdfs : 
            result = self.result  ( sample , 'BJ' )
            ecdf   = ecdfs   [ 'BJ' ]
            ## logger.info ( 'Toy results for %s/Zhang/ZC estimate'           % sample ) 
        else :
            raise KeyError (  "draw: Invalid `sample/what`: %s/%s" % ( sample , what ) )

        xmin , xmax = ecdf.xmin () , ecdf.xmax ()
        value       = result.statistics
        xmin        = min ( xmin , value )
        xmax        = max ( xmax , value )
        xmin , xmax = axis_range ( xmin , xmax , delta = 0.10 )
        ## 
        kwargs [ 'xmin' ] = kwargs.get ( 'xmin' , xmin ) 
        kwargs [ 'xmax' ] = kwargs.get ( 'xmax' , xmax )
        ##
        result    = ecdf.draw  ( option , **kwargs ) 
        line1     = ROOT.TLine ( value  , 1e-3 , value , 1 - 1e-3 )
        
        ## horisontal line 
        xmin      = kwargs [ 'xmin' ]
        xmax      = kwargs [ 'xmax' ]
        dx        = ( xmax - xmin ) / 100 
        e         = ecdf ( value )
        line2     = ROOT.TLine ( xmin + dx , e , xmax - dx , e )
        ## 
        line2.SetLineWidth ( 2 ) 
        line2.SetLineColor ( Navy       ) 
        line2.SetLineStyle ( 9 ) 
        ## 
        line1.SetLineWidth ( 4 ) 
        line1.SetLineColor ( DarkGreen  )
        ##
        line2.draw ( 'same' )
        line1.draw ( 'same' )
        ##
        self._line1 = line1
        self._line2 = line2
        ##
        return result, line1, line2   
    
    # =========================================================================
    ## merge two objects:
    def merge ( self , other ) :
        self += other
        return self
    
    ## merge two objects (needed for parallel execution):
    def __iadd__ ( self , other ) :
        """ Merge two GoF-toys objects
        - needed for paralell execution 
        """        
        if not isinstance ( other , GoFSimFit1DToys ) : return NotImplemented 

        ## (1) merge ECDFs 
        for key, content in loop_items ( other.ecdfs   ) :
            if not key in self.__ecdfs : self.__ecdfs [ key ] = content
            else :
                ecdfs = self.__ecdfs [ key ]
                for e , ecdf in content.items () :
                    if   e in ecdfs : ecdfs [ e ] += ecdf
                    else            : ecdfs [ e ]  = ecdf
                    
        ## (2) merge counters 
        for key, content in loop_items ( other.counters ) :
            if not key in self.__counters : self.__counters [ key ] = content
            else :
                counters = self.__counters [ key ]
                for e , cnt in content.items () : counters [ e ] += cnt

        ## (3) merge total
        for key, content in loop_items ( other.total ) :
            if not key in self.__total : self.__total [ key ]  = content
            else                       : self.__total [ key ] += content

        ## (4) number of toys
        self.__nToys += other.nToys
                
        return self 

# =============================================================================
## @class GoFSimFit
#  Goodness-of-fit for simultaneous fits
#  - It is a bit more  general and (a bit less efficient)  GoF estimator for 
#    simultaneous fits
class GoFSimFit(GoFSimFitBase) :
    """ Generic Goodness-of-fit for simultaneous fits
    - It is a bit more generic and (a bit less CPU efficient) GoF estimator for 
    simultaneous fits
    """
    def __init__ ( self               ,
                   pdf                ,
                   dataset            , * , 
                   estimators         , 
                   parameters  = None , **runconfig ) :

        # ========================================================================
        ## initialize the base class 
        GoFSimFitBase.__init__ ( self                    ,
                                 pdf        = pdf        ,
                                 dataset    = dataset    ,
                                 parameters = parameters )

        assert isinstance ( estimators , dictlike_types ) , \
            "Invalid type for `estimators`: %s" % ( typename ( estimators ) ) 
        assert estimators and all ( k in self.sample for k in estimators ) , \
            "Invalid keys in `estimators`: %s"  % ( ',%s' % ( k for k in estimators ) ) 
        assert all ( l in estimators for l in self.sample.labels() ) , \
            "Missing keys in `estimators`: %s"  % ( ',%s' % ( k for k in estimators ) ) 
        
        ## self.__dataset  = dataset
        
        self.__tvalues       = {} 
        self.__counters      = defaultdict(SE) 
        self.__ecdfs         = {} 
        self.__total         = EffCounter ()
        self.__nToys         = 0 
        self.__gofs_parallel = False 
        self.__time          = 0
        
        name = self.sample.name
        for key , cmp in self.pdf.categories.items ()  :
            
            assert isinstance ( cmp , APDF1 ) , "Component `%s` is not APDF1` %s" % ( key , typename ( cmp ) )
            obs      = cmp.pdf.getObservables ( dataset )
            category = '%s==%s::%s' % ( name , name , key )
            ds       = dataset.subset ( variables = obs , cuts = category )
            
            gof = estimators [ key ] 
            if isinstance ( gof , string_types ) : pass 
            
            ## get t-value for this component 
            self.__tvalues [ key ] = gof.tvalue ( cmp , ds )
                                                                
            self.gofs  [ key  ] = gof 
            self.N     [ key  ] = len ( ds ) 
        
        ## all GoFs can run in parallel ? 
        self.__gofs_parallel = all ( hasattr ( gof , 'parallel' ) and gof.parallel for gof in self.gofs.values () ) 
            
        ## run toys if requested
        if runconfig : self.run ( **runconfig )

    ## serialize the object 
    def __getstate__ ( self ) :
        """ Serialize the object 
        """
        ## (1) serialize the base 
        state = GoFSimFitBase.__getstate__ ( self )
        # 
        state [ 'tvalues'       ] = self.__tvalues 
        state [ 'counters'      ] = self.__counters
        state [ 'ecdfs'         ] = self.__ecdfs 
        state [ 'total'         ] = self.__total 
        state [ 'nToys'         ] = self.__nToys
        state [ 'gofs_parallel' ] = self.__gofs_parallel 
        state [ 'time'          ] = self.__time 
        # 
        return state 
    
    ## De-serialize the object 
    def __setstate__ ( self , state ) :
        """ De-serialize the object """        
        ## (1) de-serialize the base 
        GoFSimFitBase.__setstate__ ( self , state )
        # 
        self.__tvalues       = state.pop ( 'tvalues'       , {}    )
        self.__counters      = state.pop ( 'counters'              )
        self.__ecdfs         = state.pop ( 'ecdfs'                 )
        self.__total         = state.pop ( 'total'                 )
        self.__nToys         = state.pop ( 'nToys'         , 0     )
        self.__gofs_parallel = state.pop ( 'gofs_parallel' , False )
        self.__time          = state.pop ( 'time'          , 0     )
        
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
    ## Counters  
    @property 
    def counters ( self ) :
        """`counters` : toy results as counters"""
        return self.__counters

    # =========================================================================
    ## Total/global counters
    @property 
    def total ( self ) :
        """`total` : total/global counters
        """
        return self.__total

    # =========================================================================
    ## get all t-values for fit-components
    def tvalues ( self ) : 
        """`tvalues` : dictionary { component : t-value }
        """
        return self.__tvalues

    # ==========================================================================
    @property
    def gofs_parallel ( self ) :
        """`gofs_parallel`: all gofs can run parallel?
        """
        return self.__gofs_parallel

    # ==========================================================================
    @property
    def time  ( self ) :
        """`time` : total CPU time used [s] """
        return self.__time
        
    # =========================================================================
    ## merge two objects (needed for parallel execution):
    def merge ( self , other ) :
        """ Merge two GoF-toys objects
        - needed for parallel execution 
        """        
        self += other
        return self

    # =========================================================================
    ## merge two objects (needed for parallel execution):
    def __iadd__ ( self , other ) :
        """ Merge two GoF-toys objects
        - needed for parallel execution 
        """        
        if not isinstance ( other , GoFSimFit  ) : return NotImplemented 

        ## (1) merge ECDFs 
        for key, content in loop_items ( other.ecdfs   ) :
            if not key in self.__ecdfs    : self.__ecdfs [ key ]  = content
            else                          : self.__ecdfs [ key ] += content 
                    
        ## (2) merge counters 
        for key, content in loop_items ( other.counters ) :
            if not key in self.__counters : self.__counters [ key ]  = content
            else                          : self.__counters [ key ] += content 

        ## (3) merge total
        self.__total += other.total 
        
        ## (4) number of toys
        self.__nToys += other.nToys
                
        return self 

    # ===========================================================================
    ## run toys to get p-value
    def run ( self     ,
              nToys    = 500    ,
              parallel = False  ,
              fitconf  = {}     , 
              silent   = False  ,
              nSplit   = 0      ) :
        """ Run toys to get the p-value
        """        
        assert  isinstance ( nToys , int ) and 0 < nToys , "Invalid number of toys: %s" % nToys

        if parallel and not self.gofs_parallel : 
            logger.warning  ( "Parallel processing is switched OFF!" ) 
            parallel  = False 
            
        if parallel :
            
            from ostap.parallel.parallel_gof import parallel_goftoys as parallel_toys 
            result = parallel_toys ( gof      = self       ,
                                     nToys    = nToys      ,
                                     nSplit   = nSplit     ,
                                     fitconf  = {}         , 
                                     silent   = True       ,
                                     progress = not silent )
            
            if result : self += result
            
            return self 

        ## run with timer 
        with timing ( logger = None ) as timer : result = self.__run ( nToys , silent )
            
        if 0 < timer.delta : self.__time += timer.delta
        return result
    
    ## run toys to get p-value
    def __run ( self              ,
                nToys    = 500    ,
                silent   = False  ) : 
        """ Run toys to get the p-value
        """        
        assert  isinstance ( nToys , int ) and 0 < nToys , "Invalid number of toys: %s" % nToys

        print ( 'RUN:' , nToys )
        
        from ostap.utils.progress_bar import progress_bar
        for t in progress_bar ( nToys , silent = silent , description = 'Toys:' ) : 
              
            name = self.sample.name
            new_dataset = self.pdf.generate ( self.N , sample = True )

            ttv   = self.tvalues()
            above = True 
            ## above    = False 
         
            for sample , component in self.pdf.categories.items ()  :
                
                obs      = component.pdf.getObservables ( new_dataset )
                category = '%s==%s::%s' % ( name , name , sample )
                ds       = new_dataset.subset ( variables = obs , cuts = category )
                
                gof = self.gofs [ sample ]                
                tv  = gof.tvalue ( component , ds )

                if not sample in self.__ecdfs : self.__ecdfs [ sample ] = Ostap.Math.ECDF ( tv , True )
                else                          : self.__ecdfs [ sample ].add ( tv ) 
              
                self.__counters [ sample ] += tv 
                
                above = above and ttv [ sample ] <= tv 
                ## above = above or  ttv [ sample ] <= tv
                
                ds.clear ()
                del ds 
                
            ## the the global counter       
            self.__total += above 
                    
            self.__nToys += 1 
     
            new_dataset.clear() 
            del new_dataset 
            
    # ============================================================================
    ## get dictionary of t-values & dictionary of p-values
    def pvalues  ( self ) :
        """ Get dictionaries of t and p-values 
        """
        
        tvs = self.tvalues ()
        pvs = {}  
        for key, ecdf in self.__ecdf.items() : 
            tv = tvs.get ( key , None )
            if not tv is None : pvs [ key ] = ecdf.estimate ( tv )
            
        ## global  p-value 
        pvs [ '*SIMFIT*' ] = self.counter.eff
            
        return tvs, pvs 

    # ==================================================================================
    ## Helper method to create the informative title for the summary table
    @property 
    def title ( self ) :
        """ Helper method to create the informative title for the summary table
        """
        the_types = defaultdict(list)
        for k , g in self.gofs.items() : the_types [ typename ( g ) ].append ( k )            
        name  = ';'.join ( '[%s]:%s' % ( ','.join ( k ) , t ) for t, k in the_types.items () )
        name  = 'GoFSimFit(%s)' % name            
        title = '%s summary' % name 
        if self.nToys : title  = '%s for #%d toys' % ( title , self.nToys )
        if self.time  :
            if clock : title  = '%s %s %.1fs' % ( title , clock , self.time )
            else     : title  = '%s [%.1fs]'  % ( title         , self.time )
        return title 
        
    # ==================================================================================
    ## output as table      
    def table  ( self             , 
                 title     = ''   , 
                 prefix    = ''   , 
                 precision = 4    , 
                 width     = 6    , 
                 style     = None ) :
        """ The table 
        """ 
        
        keys = tuple ( k for k in self.sample.labels () ) 

        if not title : title = self.title

        ## no toys, no p-values, ...
        if not self.__ecdfs or not self.__counters : 
            rows = [ ('Category' , 't-value' , 'Factor' ) ] 
            for  key, tv in self.tvalues().items() :
                tt , expo = pretty_float ( tv , precision = precision , width = width ) 
                row = key, tt , '10^{%d}' % expo if expo else ''
                rows.append ( row )
            rows  = T.remove_empty_columns ( rows )
            return T.table  ( rows , title  = title , prefix = prefix , alignment = 'lcc' )

        rows    = []
        pvalues = [] 
        for key, tv in self.tvalues().items()  :
            
            cnt  = self.counters.get ( key , None )
            ecdf = self.__ecdfs .get ( key , None  ) 
      
            if ecdf is None or cnt is None : continue 
          
            pvalue = ecdf.estimate ( tv )
            pv     = clip_pvalue   ( pvalue , 0.5 )
            nsigma = significance  ( pv ) ## convert  it to significace

            mean       = cnt.mean ()
            rms        = cnt.rms    () 
            vmin, vmax = cnt.minmax () 

            mxv = max ( abs ( tv   ) , abs ( mean.value() ) ,
                        abs ( vmin ) , abs ( vmax         ) )
            mxe = max ( rms , mean.error() )

            fmt, fmtv , fmte , expo = fmt_pretty_ve ( VE ( mxv , mxe * mxe )  ,
                                                      width       = width     ,
                                                      precision   = precision ,
                                                      parentheses = False     )
            
            if expo : scale = 10 ** expo
            else    : scale = 1 
             
            smean = mean / scale 
            spv   = pv * 100 
            
            spv   = str ( ( 100 * pv  ) .toString ( '%% 5.2f %s %%-.2f' % plus_minus ) )
            nsig  = str ( nsigma.toString ( '%%.2f %s %%-.2f' % plus_minus ) if float ( nsigma ) < 1000 else '+inf' ) 

            row = ( key ,
                    fmtv % ( tv    / scale )  , 
                    fmt  % ( smean . value () , smean . error () ) ,
                    fmte % ( rms   / scale )  , 
                    fmtv % ( vmin  / scale )  , 
                    fmtv % ( vmax  / scale )  ,
                    '10^{%d}' % expo if expo else '' , 
                    spv  , 
                    nsig )  

            rows.append ( row )
            pvalues.append ( pvalue )

        rows    = sorted ( rows )

        ## global statistics
        pvalue  = self.total.eff        
        pv      = clip_pvalue ( pvalue , 0.5 )        
        nsigma  = significance ( pv ) ## convert  it to significance            
        pvalue  = str ( ( 100 * pvalue ) .toString ( '%% 5.2f %s %%-.2f' % plus_minus ) )
        sigma   = str ( nsigma.toString ( '%%.2f %s %%-.2f' % plus_minus ) if float ( nsigma ) < 1000 else '+inf' ) 
        row     = '*COMBINED*' , '' , '' , '' , '', '' , ''  , pvalue , sigma  
        rows.append ( row )
        
        if pvalues and len ( pvalues ) + 1 == len ( rows ) :
            
            clipped  = [ clip_pvalue ( p , 0.5 ) for p in pvalues ]            
            for method in  ( 'fisher' , 'pearson' , 'tippett' , 'stouffer' ) :
                
                combined = combine_pvalues ( clipped , method = method )
                
                pvalue   = combined 
                pv       = clip_pvalue ( pvalue , 0.5 )
                
                nsigma  = significance ( pv ) ## convert  it to significace            
                pvalue  = str ( ( 100 * pvalue ) .toString ( '%% 5.2f %s %%-.2f' % plus_minus ) )
                sigma   = str ( nsigma.toString ( '%%.2f %s %%-.2f' % plus_minus ) if float ( nsigma ) < 1000 else '+inf' ) 

                mm      = method.upper()                
                row    = '*%s*' % mm  , '' , '' , '' , '', '' , ''  , pvalue , sigma 
                rows.append ( row )

            minimal = min ( clipped , key = lambda p : float ( p ) )
            pvalue  = minimal  
            pv       = clip_pvalue ( pvalue , 0.5 )
                
            nsigma  = significance ( pv ) ## convert  it to significace            
            pvalue  = str ( ( 100 * pvalue ) .toString ( '%% 5.2f %s %%-.2f' % plus_minus ) )
            sigma   = str ( nsigma.toString ( '%%.2f %s %%-.2f' % plus_minus ) if float ( nsigma ) < 1000 else '+inf' ) 

            mm      = 'MINIMAL' 
            row    = '*%s*' % mm  , '' , '' , '' , '', '' , ''  , pvalue , sigma 
            rows.append ( row )
            
            rows [ -6 : ] = sorted ( rows [ -6 : ] )
        
        header = ( 'Category' , 't-value' , 't-mean' , 't-rms' , 't-min' , 't-max' , 'Factor' , 'p-value [%]' , '#%s' % greek_lower_sigma )
        rows   = [ header ] + rows 
        rows   = T.remove_empty_columns ( rows ) 
        title = title if title else 'SimFit-GoF statistics #%d' % self.nToys
        return T.table ( rows , title = title , prefix = prefix )     

    __str__  = table 
    __repr__ = table
    report   = table 
    # =========================================================================

# ==============================================================================
## Another helper class for Goodness-of-fit estimator for SimFit
class GoFSimFitType(GoFSimFit) :
    """ Another helper class for Goodness-of-Fit estimator for SimFit
    """
    def __init__ ( self              ,
                   GOF               , 
                   pdf               ,
                   dataset           , * , 
                   parameters = None ,
                   gofconfig  = {}   , **runconfig ) :

        ## proper instance ? 
        assert isinstance ( GOF ,  AGoF ) , 'Invalid type of `GOF`:%s' % typename ( GOF  )

        ## proper PDF type ?
        assert isinstance ( pdf , SimFit ) , 'Invalid type of `pdf`:%s' % typename ( pdf )
        
        estimators = {}        
        gof        = None 
        for k in pdf.categories :
            
            if   gof is None               : gof = GOF
            elif hasattr ( GOF , 'clone' ) : gof = GOF.clone () 
            elif hasattr ( GOF , 'copy'  ) : gof = GOF.copy  () 
            else                           : gof = copy.deepcopy ( GOF )
            ##
            estimators [ k ] = gof
            
        ## initialize the base class 
        GoFSimFit.__init__ ( self,
                             pdf        = pdf        ,
                             dataset    = dataset    ,
                             parameters = parameters ,
                             estimators = estimators , **runconfig )
        
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
