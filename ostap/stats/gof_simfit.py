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
    'GoFSimFit'      , ## helper utility for GoF estimate 
    'GoFSimFitToys'  , ## helper utility for GoF estimate with toys 
)
# =============================================================================
from   ostap.fitting.pdfbasic   import PDF1
from   ostap.core.core          import VE, Ostap
from   ostap.math.base          import axis_range, np2raw    
from   ostap.utils.cidict       import cidict_fun
from   ostap.utils.basic        import loop_items, typename   
from   ostap.stats.counters     import SE, EffCounter 
from   ostap.logger.pretty      import pretty_float
from   ostap.math.ve            import fmt_pretty_ve
from   ostap.math.math_ve       import significance
from   ostap.logger.symbols     import plus_minus, times, greek_lower_sigma
from   ostap.logger.colorized   import infostr
from   ostap.stats.gof_utils    import Labels, Keys, clip_pvalue
from   ostap.fitting.simfit     import SimFit 
from   ostap.stats.gof1d        import ( GoF1D , vct_clip     ,
                                         kolmogorov_smirnov   ,
                                         anderson_darling     ,
                                         cramer_von_mises     ,
                                         kuiper , ZK , ZA, ZC )
from   collections              import defaultdict, namedtuple
import ostap.stats.gofnd        as     GoFnD 
import ostap.logger.table       as     T
import ostap.fitting.roofit
import ROOT, math, numpy 
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
        """sample`: sample/category  variable for simultaneous fit
        """
        return self.pdf.sample 

    @property
    def parameters ( self ) :
        """`parameters' : fit parameters, e.g. fit-resutl or dictinoary or ...
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

    ## serialize the object 
    def __getstate__ ( self ) :
        """ Serialize the object"""
        self.__pdf.load_params ( self.parameters , silent = True )
        return { 'pdf'        : self.pdf        ,
                 'parameters' : self.parameters , 
                 'gofs'       : self.gofs       ,
                 'N'          : self.N          }
    
    ## De-serialize the object 
    def __setstate__ ( self , state ) :
        """ De-serialize the object """         
        self.__pdf        = state.pop ( 'pdf'  )
        self.__parameters = state.pop ( 'parameters' , {} ) 
        self.__gofs       = state.pop ( 'gofs' )
        self.__N          = state.pop ( 'N'    )
        self.__pdf.load_params ( self.__parameters , silent = True )
    
# =============================================================================
## @class GoFSimFit
#  Goodness-of-fit for simultaneous 1D-fits
#  - All components of the simultaneous fit must be 1D-components 
#  - GoF is estimated for each 1D-component
class GoFSimFit(GoFSimFitBase) :
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
            ds       = dataset.subset ( variables = obs ,cuts = category )
            gof      = GoF1D ( cmp , ds )
            
            self.gofs [ key ] = gof 
            self.N    [ key ] = len ( ds ) 
            
    # =========================================================================
    ## all estimators togather
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

    ## Draw fit CDF & empirical ECDF 
    def draw  ( self , sample , opts = '' , *args , **kwargs ) :
        """ Draw fit CDF & empirical CDF
        """
        gof = self.gofs.get ( sample , None )
        if gof is None : raise KeyError ( "Invalid sample `%s`" % sample )
        return gof.draw ( opts , *args , **kwargs )
    
# =============================================================================
## @class GoFSimFitToys
#  Check Goodness of 1D (Sim)Fits using toys 
class GoFSimFitToys(GoFSimFit) :
    """ Check Goodness-of-Fit with toys (Simfit caase)
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
        assert isinstance ( gof , GoFSimFit ) , "Invalid `gof`-parameter"

        ## mimic the copy-constructor for the base class 
        state = GoFSimFit.__getstate__ ( gof ) 
        GoFSimFit.__setstate__ ( self , state )

        self.__counters = { k : defaultdict(SE) for k in self.gofs }
        self.__ecdfs    = { k : {}              for k in self.gofs }
        self.__total    = defaultdict(EffCounter) 
        self.__nToys    = 0

    ## serialize the object 
    def __getstate__ ( self ) :
        """ Serialize the object 
        """
        #
        ## (1) serialize the base 
        state = GoFSimFit.__getstate__ ( self )
        # 
        state [ 'counters' ] = self.__counters
        state [ 'ecdfs'    ] = self.__ecdfs 
        state [ 'total'    ] = self.__total 
        state [ 'nToys'    ] = self.__nToys
        # 
        return state 
    
    ## De-serialize the object 
    def __setstate__ ( self , state ) :
        """ De-serialize the object """
        
        ## (1) de-serialize the base 
        GoFSimFit.__setstate__ ( self , state )
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
              silent   = False ,
              nSplit   = 0     ) :
        """ Run toys 
        """
        assert isinstance ( nToys , int ) and 0 < nToys , "Invalid `nToys` argument!"

        if parallel :
            
            from ostap.parallel.parallel_gof1d import parallel_gof1dtoys as parallel_toys 
            self += parallel_toys ( gof      = self       ,
                                    nToys    = nToys      ,
                                    nSplit   = nSplit     ,
                                    silent   = True       ,
                                    progress = not silent )
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
                
                total_KS = total_KS and gof .kolmogorov_smirnov_estimator <= ks
                total_K  = total_K  and gof .            kuiper_estimator <= k
                total_AD = total_AD and gof .  anderson_darling_estimator <= ad
                total_CM = total_CM and gof .  cramer_von_mises_estimator <= cm
                total_ZK = total_ZK and gof .                ZK_estimator <= zk 
                total_ZA = total_ZA and gof .                ZA_estimator <= za 
                total_ZC = total_ZC and gof .                ZC_estimator <= zc 
                
                cnts = counters [ key ]                
                cnts [ 'KS' ] += ks
                cnts [ 'K'  ] += k
                cnts [ 'AD' ] += ad
                cnts [ 'CM' ] += cm
                cnts [ 'ZK' ] += zk
                cnts [ 'ZA' ] += za
                cnts [ 'ZC' ] += zc

                res  = results [ key ]
                res [ 'KS'  ].append ( ks )    
                res [ 'K'   ].append ( k  )
                res [ 'AD'  ].append ( ad ) 
                res [ 'CM'  ].append ( cm ) 
                res [ 'ZK'  ].append ( zk ) 
                res [ 'ZA'  ].append ( za ) 
                res [ 'ZC'  ].append ( zc ) 

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
                 fmtv   % vr                          ,
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
        for sample , ecdfs in self.ecdfs.items()  :
            for label in ecdfs :
                result  = self.result ( sample , label )
                if not result : continue
                the_label = Labels.get ( label , label )
                row = self.row ( sample , the_label , result , width = width , precision = precision )
                rows.append ( row ) 
                    
        if   not title and self.nToys : title = 'Goodness of 1D-fit with #%d toys' % self.nToys  
        elif not title                : title = 'Goodness of 1D-fit'

        for e , cnt in self.total.items() :

            ## get the binomial efficiency 
            pvalue = cnt.efficiency

            pv     = clip_pvalue ( pvalue , 0.5 )
            nsigma = significance ( pv ) ## convert  it to significace
            
            pvalue = str ( ( 100 * pvalue ) .toString ( '%% 5.2f %s %%-.2f' % plus_minus ) )
            sigma  = str ( nsigma.toString ( '%%.2f %s %%-.2f' % plus_minus ) if float ( nsigma ) < 1000 else '+inf' ) 

            label = Labels.get( e , e ) 
            row   = label , infostr ( '***' ) , '---' , '---' , '---' , '---' , ''  , pvalue , sigma 
            rows.append ( row ) 

        rows  = [ header ] + sorted ( rows )
        rows  = T.remove_empty_columns ( rows ) 
        return T.table ( rows , title = title , prefix = prefix , alignment = 'lcccccccc' , style = style )
        
    __repr__ = table
    __str__  = table
    
    # =========================================================================
    ## Draw ECDF for toys & statistical estimator 
    def draw  ( self , sample , what , opts = '' , *args , **kwargs ) :
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
        else :
            raise KeyError (  "draw: Invalid `sample/what`: %s/%s" % ( sample , what ) )
            
        xmin , xmax = ecdf.xmin () , ecdf.xmax ()
        value       = result.statistics
        xmin        = min ( xmin , value )
        xmax        = max ( xmax , value )
        xmin , xmax = axis_range ( xmin , xmax , delta = 0.20 )

        kwargs [ 'xmin' ] = kwargs.get ( 'xmin' , xmin ) 
        kwargs [ 'xmax' ] = kwargs.get ( 'xmax' , xmax )

        result    = ecdf.draw  ( opts , *args , **kwargs ) 
        line1     = ROOT.TLine ( value , 1e-3 , value , 1 - 1e-3 )
        
        ## horisontal line 
        xmin      = kwargs['xmin']
        xmax      = kwargs['xmax']
        dx        = ( xmax - xmin ) / 100 
        e         = ecdf ( value )
        line2     = ROOT.TLine ( xmin + dx , e , xmax - dx , e )
        ## 
        line2.SetLineWidth ( 2 ) 
        line2.SetLineColor ( 4 ) 
        line2.SetLineStyle ( 9 ) 
        ## 
        line1.SetLineWidth ( 4 ) 
        line1.SetLineColor ( 8 )
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
    
    ## merge two objects:
    def __iadd__ ( self , other ) :
        """ Merge two GoF-toys objects
        """        
        if not isinstance ( other , GoFSimFitToys ) : return NotImplemented 

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

# ==============================================================================
## Another helper base class for Goodness-of-fit estimator for SimFit
class GoFSimFitType(GoFSimFitBase) :
    """ Another helper base class for Goodness-of-fit estimator for SimFit
    """
    def __init__ ( self               ,
                   GOF_type           , 
                   pdf                ,
                   dataset            ,
                   parameters  = None , **config ) :
        
        ## initialize the base class 
        GoFSimFitBase.__init__ ( self,
                                 pdf        = pdf        ,
                                 dataset    = dataset    ,
                                 parameters = parameters )
        
        
        self.__cmp = {} 
        ## create the GoF-objects for each component 
        name = self.sample.name
        for key , cmp  in self.pdf.categories.items ()  :
            
            obs      = cmp.pdf.getObservables ( dataset )
            category = '%s==%s::%s' % ( name , name , key ) 
            dset     = dataset.subset ( variables = obs ,cuts = category )

            self.__cmp [ key ] = cmp , dset

            ## the method 
            self.gofs  [ key ] = GOF_type ( **config )
            self.N     [ key ] = len ( dset )
            
        self.__tvalues = {}
        self.__pvalues = {}
        
    ## serialize the object 
    def __getstate__ ( self ) :
        """ Serialize the object"""
        state = GoFSimFitBase.__getstate__ ()
        state [ 'components' ] = self.__cmp
        return state 
    
    ## De-serialize the object 
    def __setstate__ ( self , state ) :
        """ De-serialize the object """
        self.__cmp        = state.pop  ( 'components' )
        GoFSimFitBase.__setstate__ ( state )
        
    # =========================================================================
    ## Get dictionary of t-values (for each SimFit component)
    def __call__ ( self ) :
        """ Get dictionary of t-values (for each SimFit component)
        """
        if not self.__tvalues :
            ## evaluate t-values for each component 
            for key , value in self.__cmp.items()  :
                cmp, data = value
                gof = self.gofs [ key ]
                self.__tvalues [ key ] = gof ( cmp , data )
        return self.__tvalues
    
    # =========================================================================
    ## Get dictionary of t-values (for each SimFit component)
    tvalues = __call__
    
    # =========================================================================
    ## Get dictionary of t&p-values (for each SimFit component)
    def pvalues ( self ) :
        """ Get dictionary of t-values (for each SimFit component)
        """
        if not self.__pvalues : 
            ## evaluate t&p-values for each component 
            for key , value in self.__cmp.items()  :
                cmp, data = value
                gof = self.gofs [ key ]
                if not gof.silent :
                    logger.info ( 'GoF.pvalues[%s]: processing sample "%s"' % ( typename ( gof ) , key ) )
                self.__pvalues [ key ] = gof.pvalue ( cmp , data )                
        return self.__pvalues
    
    # =========================================================================
    ## Draw ECDF for toys & statistical estimator 
    def draw  ( self , sample , opts = '' , *args , **kwargs ) :
        """ Draw ECDF for toys & statistical estgimator 
        """
        gof = self.gofs.get ( sample , None ) 
        if not gof : raise KeyError ( 'Unknown sample "%s"!' % sample )
        ## draw it
        t = self.tvalues()[sample]
        return gof.draw ( tvalue = t , opts = opts , *args , **kwargs )
    
    # =========================================================================
    def table ( self             , * , 
                title     = ''   ,
                prefix    = ''   ,
                precision = 4    , 
                width     = 6    ,
                style     = None ) :
        """ Print the summary as Table  (for SimFit)
        """
        gof_type = typename ( self )             
        tvalues  = self.tvalues () 
        pvalues  = self.pvalues ()
        title    = title if title else 'GoF:%s' % gof_type
        
        tvalues  = self.tvalues () 
        pvalues  = self.pvalues () 

        rows     = [ ( 'Sample'            , 
                       't-value'           , 
                       't-mean'            ,
                       't-rms'             ,
                       't-min/max'         ,
                       'factor'            ,
                       'p-value [%]'       ,
                       '#%s' % greek_lower_sigma ) ] 
                     
        for sample in sorted ( pvalues ) :
            
            tvalue , pvalue = pvalues [ sample ]
            tv     , te     = pretty_float ( tvalue , width = width , precision = precision  , with_sign = True ) 
            
            if te : row = sample , tv , '%s10^%+d' % ( times , te  ) 
            else  : row = sample , tv , ''
            
            ## ECDF ?
            gof  = self.gofs [ sample ]
            ecdf = gof.ecdf

            if ecdf :
                counter    = ecdf   .counter ()
                tmean      = counter.mean    ()
                trms       = counter.rms     ()
                tmin, tmax = counter.minmax  ()
                mxv = max ( abs ( tvalue )        ,
                            abs ( tmean.value() ) ,
                            tmean.error()         , trms ,
                            abs ( tmin )  , abs ( tmax ) )                        
                fmt, fmtv , fmte , expo = fmt_pretty_ve ( VE ( mxv ,  tmean.cov2() ) ,
                                                          width       = width       ,
                                                          precision   = precision   , 
                                                          parentheses = False       )                    
                if expo : scale = 10**expo
                else    : scale = 1
                
                tv    = tvalue / scale
                tmean = tmean  / scale 
                trms  = trms   / scale 
                tmin  = tmin   / scale
                tmax  = tmax   / scale
                
                row = ( sample        ,
                        fmtv % tvalue ,
                        fmt  % ( tmean.value() , tmean.error() ) , 
                        fmtv % trms   , 
                        '%s/%s' % ( fmtv % tmin , fmtv % tmax ) )
                
                if expo : row +=  '%s10^%+d' % ( times , expo ) ,
                else    : row +=  '' , 
                
            
            pv     = clip_pvalue ( pvalue , 0.5 ) 
            nsigma = significance ( pv ) ## convert  it to significace
            
            p = pvalue * 100 
            pvalue = '% 5.2f %s %-.2f' % ( p.value() , plus_minus , p.error () )
            n = nsigma 
            nsigma = '%.2f %s %-.2f'   % ( n.value() , plus_minus , n.error () ) if float ( nsigma ) < 1000 else '+inf'
            
            row +=  pvalue , nsigma
            
            rows.append ( row ) 

        rows  = T.remove_empty_columns ( rows )
        title = title if title else 'Goodness-Of-Fit %s for simfit' % gof_type  
        return T.table ( rows                 ,
                         title     = title    ,
                         prefix    = '# '     ,
                         alignment = 'llcccc' ) 
        
# =============================================================================
## Goodness-of-fit estimator for SimFit using PPD method
class PPDSimFit(GoFSimFitType) :
    """ Goodness-of-fit estimator for SimFit using PPD method
    """
    def __init__ ( self               ,
                   pdf                ,
                   dataset            ,
                   parameters  = None , **config ) :
        
        GoFSimFitType.__init__ ( self ,
                                 GOF_type   = GoFnD.PPD  , 
                                 pdf        = pdf        ,
                                 dataset    = dataset    ,
                                 parameters = parameters , **config ) 

# =============================================================================
## Goodness-of-fit estimator for SimFit using DNN method
class DNNSimFit(GoFSimFitType) :
    """ Goodness-of-fit estimator for SimFit using DNN method
    """    
    def __init__ ( self               ,
                   pdf                ,
                   dataset            ,
                   parameters  = None , **config ) :
        
        GoFSimFitType.__init__ ( self ,
                                 GOF_type   = GoFnD.DNN  , 
                                 pdf        = pdf        ,
                                 dataset    = dataset    ,
                                 parameters = parameters , **config ) 

    # ==========================================================================
    ## Get dictionary of histograms with u-value distribution
    def uvalues ( self ) :
        """ Get dictionary of histograms with u-value distribution
        """
        uvalues = {}
        ## evaluate t&p-values for each component 
        for key , value in self.__cmp.items()  :
            cmp, data = value
            gof = self.gofs [ key ]
            uvalues [ key ] = gof.histo 
        return uvalues

# =============================================================================
## Goodness-of-fit estimator for SimFit using USTAT method
class USTATSimFit(GoFSimFitType) :
    """ Goodness-of-fit estimator for SimFit using USTAT method
    """    
    def __init__ ( self               ,
                   pdf                ,
                   dataset            ,
                   parameters  = None , **config ) :
        
        GoFSimFitType.__init__ ( self ,
                                 GOF_type   = GoFnD.USTAT , 
                                 pdf        = pdf        ,
                                 dataset    = dataset    ,
                                 parameters = parameters , **config ) 
        
    # ==========================================================================
    ## Get dictionary of histograms with u-value distribution
    def uvalues ( self ) :
        """ Get dictionary of histograms with u-value distribution
        """
        uvalues = {}
        ## evaluate t&p-values for each component 
        for key , value in self.__cmp.items()  :
            cmp, data = value
            gof = self.gofs [ key ]
            uvalues [ key ] = gof.histo 
        return uvalues
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================


