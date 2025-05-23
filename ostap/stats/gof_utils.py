#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/stats/gof_utils.py
#  Set of utulities for goodness-of-fit studies 
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2023-12-06
# =============================================================================
""" Simple utulities for goodness-of-fit studies 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2023-12-06"
__all__     = (
    'mean_var'   , ## mean and variance for (weighted) arrays
    'nEff'       , ## get number of effective entries
    'normalize'  , ## "normalize" variables in dataset/structured array
    'Estimators' , ## helper mixin class to print statistical estimators 
    'Summary'    , ## helper mixin class to print statistical estimators 
    'GoFSummary' , ## helper class to print summary of GoF methods 
)
# =============================================================================
from   collections              import namedtuple
from   ostap.core.core          import VE, Ostap
from   ostap.utils.cidict       import cidict_fun
from   ostap.core.ostap_types   import string_types
from   ostap.math.base          import axis_range, numpy   
from   ostap.math.math_ve       import significance
from   ostap.stats.counters     import EffCounter
from   ostap.utils.basic        import numcpu, loop_items 
from   ostap.utils.utils        import splitter
from   ostap.utils.memory       import memory_enough 
from   ostap.utils.progress_bar import progress_bar
import ROOT, sys, warnings, math  
# =============================================================================
## float types in pumy 
np_floats = ( numpy.float16  ,
              numpy.float32  ,
              numpy.float64  ,
              numpy.float128 ) if numpy else () 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.gof_utils' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Simple utilities for goodness-of-fit studies' )
# ============================================================================
## Get the mean and variance for (1D) data array with optional (1D) weight array
#  @code
#  ds = ... ## dataste as structured array
#  mean, cov2 = mean_var ( ds ['x'] )
#  @endcode
#  - with weight 
#  @code
#  ds = ... ## dataset as structured array with weight 
#  mean, cov2 = mean_var ( ds ['x'] , ds['weight'] )
#  @endcode 
def mean_var ( data , weight = None ) :
    """ Get the mean and variance for 1D-data array with optional 1D-weight array

    >>> ds = ... ## dataset as structured array
    >>> mean, cov2 = mean_var ( ds ['x'] )
    
    - with weight 
    
    >>> ds = ... ## dataset as structured array with weight 
    >>> mean, cov2 = mean_var ( ds ['x'] , ds['weight'] )
    """
    #
    if weight is None :
        mean = numpy.mean ( data , axis = 0 , dtype = numpy.float64 ) 
        var  = numpy.var  ( data , axis = 0 , dtype = numpy.float64 ) 
        return mean , var
    # 
    mean  = numpy.average (   data               , weights = weight , axis = 0 )
    var   = numpy.average ( ( data - mean ) ** 2 , weights = weight , axis = 0 )
    #
    return mean , var 

# =============================================================================
## Get the effectibe number of entries for 1D-array
#  \f{ n_{eff} = \frac{  \left\langle x \right\rangle^2}
#                     { \left\langle x^2 \right\rangle } \f}
def nEff ( weights ) :
    """ Get the effective number of entries for 1D-array
    n_eff = ( sum ( x )  ) ^2 / sum ( x^2 )
    """

    s1 = numpy.sum ( weights      , dtype = numpy.float64 )
    s2 = numpy.sum ( weights ** 2 , dtype = numpy.float64 )
    
    return s1 * s1 / s2

# =============================================================================
## Get the "normalized" input datasets
#  All floating felds  are calculated as
#  \f[ x = \frac{x - \left\langle x \right\rangle}{\sigma} \f]
#  where \f$ \left\langle x \right\rangle\f$ is mena value
#  and \f$ \sigma \f$ is a standard deviation.
# 
#  @code
#  ds      = ... # data set as structured array
#  dsn, J  = normalize ( ds ) 
#  @endcode
#
#  - If several datasets are specified, all floating names must be the same
#  and the mean and sigma are either taken either from the first dataset,
#  if <code>first=True</code> or as combined through all datasets otherwise 
#
#  @code
#  ds1 = ... # data set as structured array
#  ds2 = ... # data set as structured array
#  ds3 = ... # data set as structured array
#  ds1n, ds2n, ds3n = normalize ( ds1 , ds2 , ds3 , first = True )
#  @endcode
#
#  - If <code>weight</code> is specified, this floating column is considered
#  as the weight
#
#  @code
#  ds  = ... # data set as structured array with weight 
#  dsn = normalize ( ds , weight = 'weight' ) 
#  @endcode
#
#  @code
#  ds1 = ... # data set as structured array without weight 
#  ds2 = ... # data set as structured array with weight 
#  datasets    = normalize ( ds1 , ds2 , weight = ( None , 'weight'  ) )
#  ds1n , ds2n = datasets 
#  @endcode
#
#  @attention Only the floating point columns are transformed! 
#  @attention Input datasets are expected to be numpy structured arrays
#
#  @code
#  ds = ... # data set as structured array
#  dsn = normalize ( ds ) 
#  @endcode
def normalize ( ds , *others , weight = () , first = True ) :
    """ Get the `normalized' input datasets
    All floating fields  are calculated as
    
    x = (x - <x>)/sigma
    
    - <x> is a mean value 
    - is a standard deviation.
    
    - If several datasets are specified, all floating names must be the same
    and the mean and sigma are either taken either from the first dataset,
    if `first=True` or as combined through all datasets, otherwise  
    
    - If `weight` is specified, this floating column is concidered
    as the weight
    
    - attention Only the floating point columns are transformed! 
    - attention Input datasets are expected to be numpy structured arrays 
    """

    datasets = ( ds , ) + others
    
    nd = len ( datasets ) 
    if not weight                             : weight = nd * [ ''     ]
    elif isinstance ( weight , string_types ) : weight = nd * [ weight ]
    
    assert ( len ( weight ) == nd ) and \
           all ( ( not w ) or isinstance ( w , string_types ) for w in weight ) , \
           'Invalid specification of weight!'
        
    weight = list ( weight )
    for i , w in enumerate ( weight ) :
        if not w : weight [ i ] = '' 
    weight = tuple ( weight )

    ds     = datasets [ 0  ]
    others = datasets [ 1: ]
    
    ## collect the floating columns 
    columns = []
    w0      = weight [ 0 ] 
    for n,t in ds.dtype.fields.items () :
        if t [ 0 ] in np_floats  and n != w0 : columns.append ( n ) 
        
    vmeans  = [] 
    for i , c in enumerate ( columns ) :
        mean, var    = mean_var ( ds [ c ] , None if not w0 else ds [ w0 ] )
        vmeans.append ( VE ( mean , var ) )
        
    ## Number of events/effective entries 
    nevents = 1.0 * ds.shape [ 0 ] if not w0 else nEff ( ds [ w0 ] )
    
    if not first and others : 
        nevents = ds.shape[0] 
        for k , dd in enumerate ( others ) :
            
            wk = weight [ k + 1 ] 
            nn = 1.0 * dd.shape [ 0 ] if not wk else nEff ( dd [ wk ] )                
            
            for i , c in enumerate ( columns ) :
                
                mean, var = mean_var ( dd [ c ] , None if not wk else dd [ wk ] )                    
                vv = VE ( mean , var ) 
                vmeans [ i ] = Ostap.Math.two_samples ( vmeans [ i ] , nevents , vv , nn )
                
            nevents += nn
                
    result = []  
    for d in datasets :
        
        nds = d.copy ()
        for ic , c in enumerate ( columns ) :
            vv        = vmeans [ ic ]
            mean      = vv.value ()
            sigma     = vv.error ()                 
            a         = nds [ c ]
            nds [ c ] =  ( a - mean ) / sigma
            
        result.append ( nds )
            
    return tuple ( result ) 


# =============================================================================
## @class PERMUTATOR
#  Helper class that allow to run permutattion test in parallel 
class PERMUTATOR(object) :
    """ Helper class that allow to run permutation test in parallel 
    """
    def __init__ ( self, gof, t_value , ds1 , ds2 ) :
        self.gof     = gof
        self.ds1     = ds1
        self.ds2     = ds2
        self.t_value = t_value
        
    # =========================================================================
    ## run N-permutations 
    def __call__ ( self , N , silent = True ) :
        numpy.random.seed()
        n1      = len ( self.ds1 )
        pooled  = numpy.concatenate ( [ self.ds1 , self.ds2 ] )
        counter = EffCounter()
        for i in progress_bar ( N , silent = silent , description = 'Permutations:') : 
            numpy.random.shuffle ( pooled )            
            tv       = self.gof.t_value ( pooled [ : n1 ] , pooled [ n1: ] )
            counter += bool ( self.t_value < tv  )            
        del pooled
        return counter

# =============================================================================
jl = None 
# =============================================================================
if False : # ==================================================================
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        with warnings.catch_warnings(): 
            warnings.simplefilter ( "ignore" , category = DeprecationWarning  )
            import joblib as jl
            jl_version = tuple ( int ( i ) for i in jl.__version__.split('.') )
        # =====================================================================
        ## Run NN-permutations in parallel using joblib 
        def joblib_run ( self , NN , silent = True ) :
            """ Run NN-permutations in parallel using joblib """
            me    = math.ceil ( memory_enough() ) + 1 
            nj    = min ( 2 * numcpu () + 3 , me ) 
            lst   = splitter ( NN , nj )
            if not silent : logger.info ( 'permutations: #%d parallel subjobs to be used' % nj ) 
            ## 
            conf  = { 'n_jobs' : -1 , 'verbose' : 0 }
            if    (1,3,0) <= jl_version < (1,4,0) : conf [ 'return_as' ] = 'generator'           
            elif  (1,4,0) <= jl_version           : conf [ 'return_as' ] = 'unordered_generator'
            ##
            input   = ( jl.delayed (self)( N ) for N in lst )
            counter = EffCounter()
            with warnings.catch_warnings(): 
                warnings.simplefilter ( "ignore" ) ## , category = DeprecationWarning  )
                results = jl.Parallel ( **conf ) ( input )
                for r in progress_bar ( results , max_value = nj , silent = silent , description = 'Permutations:') :
                    counter += r
            ## 
            return counter
        # =====================================================================
        PERMUTATOR.run = joblib_run        
        # =====================================================================
        logger.debug ( 'Joblib will be  used for parallel permutations')
        # =====================================================================        
    except ImportError : # ====================================================
        # =====================================================================
        jl = None

# =============================================================================
if not jl : # =================================================================
    # =========================================================================
    ## Run NN-permutations in parallel using WorkManager
    def pp_run ( self , NN , silent = True ) :
        """ Run NN-permutations in parallel using WorkManager"""
        me    = math.ceil ( memory_enough() ) + 1 
        nj    = min ( 2 * numcpu () + 3 , me ) 
        lst   = splitter ( NN , nj )
        ##
        if not silent : logger.info ( 'permutations: #%d parallel subjobs to be used' % nj ) 
        counter = EffCounter()
        ## 
        ## use the bare interface 
        from ostap.parallel.parallel import WorkManager
        with WorkManager ( silent = silent ) as manager : 
            for result in manager.iexecute ( self , lst , progress = not silent  , njobs = nj , description = 'Permutations:') :
                counter += result 
        # 
        return counter
    # =========================================================================
    logger.debug ( 'Parallel will be  used for parallel permutations')
    # =====================================================================        
    PERMUTATOR.run = pp_run
    # =========================================================================

# =============================================================================
## @class TOYS
#  Helper class to run toys for Goodness-of-Fit studies 
class TOYS(object) :
    """ Helper class that allow to run permutation test in parallel 
    """
    def __init__ ( self, gof , t_value , pdf , Ndata , sample = False ) :
        self.gof     = gof
        self.pdf     = pdf
        self.Ndata   = Ndata 
        self.t_value = t_value
        self.sample  = gof.sample 
        self.silent  = gof.silent
        
    # =========================================================================
    ## run N-toys 
    def __call__ ( self , N , silent = True ) :
        """ Run N-toys """
        ROOT.gRandom.SetSeed() 
        ROOT.RooRandom.randomGenerator().SetSeed()
        counter = EffCounter()
        for i in range ( N ) :
            
            ds = self.pdf.generate ( self.Ndata , sample = self.sample )
            tv = self.gof ( self.pdf , ds )
            counter += bool ( self.t_value < tv  )
           
            if isinstance ( ds , ROOT.RooDataSet ) :
                ds = Ostap.MoreRooFit.delete_data ( ds )
            del ds
            
        return counter 

    # =========================================================================
    ## Run N-toys in parallel using WorkManager
    def run ( self , NN , silent = False ) :
        """ Run NN-permutations in parallel using WorkManager"""
        me    = math.ceil ( memory_enough() ) + 1 
        nj    = min ( 2 * numcpu () + 3 , me ) 
        lst   = splitter ( NN , nj )
        if not silent : logger.info ( 'toys: #%d parallel subjobs to be used' % nj ) 
        ##
        counter = EffCounter()
        ## 
        ## use the bare interface 
        from ostap.parallel.parallel import WorkManager
        with WorkManager ( silent = silent ) as manager : 
            for result in manager.iexecute ( self , lst , progress = not silent , njobs = nj , description = 'Toys:' ) :
                counter += result 
        # 
        return counter

# =============================================================================
## Short labels for various statitical estimators 
Labels = {
    'KS' : 'Kolmogorov-Smirnov' ,
    'K'  : 'Kuiper'             ,
    'AD' : 'Anderson-Darling'   ,
    'CM' : 'Cramer-von Mises'   ,
    'ZK' : 'Zhang/ZK'           ,
    'ZA' : 'Zhang/ZA'           ,
    'ZC' : 'Zhang/ZC'           ,        
}
# ============================================================================
## @class Estimators
#  Helper mixin class to format the table of estimators 
class Estimators(object) :
    """ Helper mixin class to format the table of estimators 
    """
    # ==========================================================================
    ## Print the summary as Table
    def table ( self , title = '' , prefix = '' , width = 6 , precision = 4 , style = None ) :
        """ Print the summary as Table
        """
        import ostap.logger.table  as     T 
        from   ostap.logger.pretty import pretty_float
        ##
        rows = [ ( 'Statistics' , 'Value' ) ]
        for label , value  in loop_items ( self.estimators ) :
            
            the_label = Labels.get ( label , label )

            result , expo = pretty_float ( value , width = width , precision = precision )

            if expo : row = the_label , result , '10^%+d' % expo
            else    : row = the_label , result 
            rows.append ( row )
            
        title = title if title else 'Goodness of 1D-fit' 
        return T.table ( rows , title = title , prefix = prefix , alignment = 'lcl' , style = style  )

# ============================================================================
## @class Summary
#  Helper mixin class to keep (and format) the statistics
class Summary(object) :
    """ Helper mixin class to keep (and format) the statistics
    """
    # =============================================================================
    KS_keys = 'ks' , 'kolmogorov' , 'kolmogorovsmirnov' 
    K_keys  = 'k'  , 'kuiper'  
    AD_keys = 'ad' , 'anderson'   , 'andersondarling' 
    CM_keys = 'cm' , 'cramer'     , 'cramervonmises' 
    ZK_keys = 'zk' , 'zhangk'     , 'zhangzk'
    ZA_keys = 'za' , 'zhanga'     , 'zhangza'
    ZC_keys = 'zc' , 'zhangc'     , 'zhangzc'
    # =========================================================================
    ## result of toys/permutations 
    Result = namedtuple ( 'Result' , 'statistics counter pvalue nsigma' )
    # =========================================================================
    ## Helper method to get result
    def result ( self , label ) :
        """ Helper method to get result 
        """
        if not label in self.estimators : return None
        if not label in self.ecdfs      : return None
        if not label in self.counters   : return None
        ##
        value   = self.estimators   [ label ]
        ecdfs   = self.ecdfs        [ label ] 
        counter = self.counters     [ label ] 
        ##
        pvalue = ecdfs. estimate ( value ) ## get the p-value 
        nsigma = significance ( pvalue )   ## convert oit to significace 
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
        
        from   ostap.math.ve import fmt_pretty_ve         
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
        
        from ostap.logger.symbols import plus_minus, times         
        return ( what  ,
                 fmtv  % vs ,
                 fmt   % ( vm.value() , vm.error() ) ,
                 fmtv  % vr                          ,
                 fmt2  %  ( vmn , vmx )              ,
                 ( '%s10^%+d' %  ( times , expo )  if expo else '' )   ,                  
                 ( 100 * pvalue ) .toString ( '%% 5.2f %s %%-.2f' % plus_minus ) , 
                 ( nsigma       ) .toString ( '%%.2f %s %%-.2f' % plus_minus ) )
    
    # =========================================================================
    ## Make a summary table
    def table ( self , title = '' , prefix = '' , width = 6 , precision = 4 , style = None ) :
        """ Make a summary table
        """
        import ostap.logger.table  as     T                 
        rows = [ ( 'Statistics' , 'value' , 'mean' , 'rms' , 'min/max' , 'factor' , 'p-value [%]' , '#sigma' ) ]
        
        for label in self.ecdfs :
            
            result  = self.result ( label )
            if not result : continue

            the_label = Labels.get ( label , label )
            row = self.row ( the_label , result , width = width , precision = precision )
            rows.append ( row ) 
                    
        if rows : rows = T.remove_empty_columns ( rows )
        if   not title and self.nToys :
            title = 'Goodness of 1D-fit with #%d toys' % self.nToys  
        elif not title :
            title = 'Goodness of 1D-fit'

        return T.table ( rows , title = title , prefix = prefix , alignment = 'lccccccc' , style = style )

    # =========================================================================
    ## Draw ECDF for toys & statistical estgimator 
    def draw  ( self , what , opts = '' , *args , **kwargs ) :
        """ Draw ECDF for toys & statistical estgimator 
        """
        key = cidict_fun ( what ) 
        if   key in self.KS_keys and 'KS' in self.ecdfs :             
            result = self.result ( 'KS' )
            ecdf   = self.ecdfs  [ 'KS' ]
            logger.info ( 'Toy resuls for Kolmogorov-Smirnov estimate' ) 
        elif key in self.K_keys  and 'K'  in self.ecdfs : 
            result = self.result ( 'K' )
            ecdf   = self.ecdfs  [ 'K' ]
            logger.info ( 'Toy resuls for Kuiper estimate' ) 
        elif key in self.AD_keys and 'AD' in self.ecdfs :             
            result = self.result ( 'AD' )
            ecdf   = self.ecdfs  [ 'AD' ]
            logger.info ( 'Toy resuls for Anderson-Darling estimate' ) 
        elif key in self.CM_keys  and 'CM' in self.ecdfs : 
            result = self.result  ( 'CM' )
            ecdf   = self.ecdfs   [ 'CM' ]
            logger.info ( 'Toy resuls for Cramer-von Mises  estimate' ) 
        elif key in self.ZK_keys  and 'ZK' in self.ecdfs : 
            result = self.result  ( 'ZK' )
            ecdf   = self.ecdfs   [ 'ZK' ]
            logger.info ( 'Toy resuls for Zhang/ZK estimate' ) 
        elif key in self.ZA_keys  and 'ZA' in self.ecdfs :  
            result = self.result  ( 'ZA' )
            ecdf   = self.ecdfs   [ 'ZA' ]
            logger.info ( 'Toy resuls for Zhang/ZA estimate' ) 
        elif key in self.ZC_keys and 'ZC' in self.ecdfs : 
            result = self.result  ( 'ZC' )
            ecdf   = self.ecdfs   [ 'ZC' ]
            logger.info ( 'Toy resuls for Zhang/ZC estimate' ) 
        else :
            raise KeyError (  "draw: Invalid `what`:%s" % what )
            
        xmin , xmax = ecdf.xmin () , ecdf.xmax ()
        value     = result.statistics
        xmin      = min ( xmin , value )
        xmax      = max ( xmax , value )
        xmin , xmax  = axis_range ( xmin , xmax , delta = 0.10 )

        kwargs [ 'xmin' ] = kwargs.get ( 'xmin' , xmin ) 
        kwargs [ 'xmax' ] = kwargs.get ( 'xmax' , xmax )
        
        result    = ecdf.draw  ( opts , *args , **kwargs ) 
        line      = ROOT.TLine ( value , 0 , value , 1 )
        ## 
        line.SetLineWidth ( 4 ) 
        line.SetLineColor ( 8 ) 
        line.draw('same')
        ##
        self._line = line
        return result, line  

# =============================================================================
## @class GoFSummary
#  Helper class for format summary table 
class GoFSummary ( object) :
    
    def __init__ ( self ) :
        self.__header = ( 'Method' , 't-value' , '' , '#Toys' , 'p-value [%]' , '#sigma' ) 
        self.__items  = []

    @property 
    def items  ( self ) :
        """`items`: get all items"""
        return tuple ( self.__items )
        
    def add_row ( self , method , tvalue , pvalue , nToys ) :
        item = method, tvalue, pvalue, nToys 
        self.__items.append ( item ) 

    # =======================================================================
    ## Make a summary table 
    def table ( self , title = '' , prefix = '' , width = 5 , precision = 3 , style = '' ) :
        """ Make a summary table 
        """
        import ostap.logger.table  as     T 
        from   ostap.logger.pretty import pretty_float
        ## 
        rows = [ self.__header ]
        for item in self.__items :
            m, t, p, n = item
            sv = significance ( p ) 
            tv , texpo = pretty_float ( t , width = width , precision = precision )
            pv   = p * 100
            pval = '%5.3f +/- %.3f' %  ( pv.value() , pv.error() )
            sval = '%.2f +/- %.2f'  %  ( sv.value() , sv.error() )
            row  = m , tv , '[10^%+d]' % texpo if texpo else '' , '%d' % n , pval , sval 
            rows.append ( row )

        rows = T.remove_empty_columns ( rows )
        title = title if title else 'Goodness of 1D-fit' 
        return T.table ( rows , title = title , prefix = prefix , alignment = 'lclcc', style = style  )
                
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    if not numpy : logger.warning ("Numpy  is not avalaille!") 
    if not jl    : logger.warning ("Joblib is not avalaille!") 

# =============================================================================
##                                                                      The END 
# =============================================================================


    
