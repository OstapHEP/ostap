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
    'mean_var'    , ## mean and variance for (weighted) arrays
    'nEff'        , ## get number of effective entries
    'normalize'   , ## "normalize" variables in dataset/structured array
    'clip_pvalue' , ## clip-value 
)
# =============================================================================
from   ostap.core.meta_info     import root_info 
from   ostap.core.core          import VE, Ostap
from   ostap.utils.cidict       import cidict, cidict_fun
from   ostap.core.ostap_types   import string_types, num_types 
from   ostap.math.base          import doubles, axis_range
from   ostap.math.math_ve       import significance
from   ostap.math.ve            import fmt_pretty_ve 
from   ostap.stats.counters     import SE, WSE, EffCounter
from   ostap.utils.basic        import numcpu, loop_items 
from   ostap.utils.utils        import splitter
from   ostap.utils.memory       import memory_enough 
from   ostap.utils.progress_bar import progress_bar
from   ostap.logger.symbols     import times, plus_minus, greek_lower_sigma
from   ostap.logger.pretty      import pretty_float
from   ostap.plotting.color     import Orange, Green, Blue
import ostap.logger.table       as     T 
import ROOT, sys, math, numpy  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.stats.gof_utils' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Simple utilities for goodness-of-fit studies' )#
# =============================================================================
# float types in numpy
np_floats = ( numpy.float16  ,
              numpy.float32  ,
              numpy.float64  ,
              numpy.float128 )
# ============================================================================
if  ( 6 , 32 ) <= root_info : data2vct = lambda s : s
else                        : data2vct = lambda s : doubles ( s ) 
# =============================================================================
## Get the mean and variance for (1D) data array with optional (1D) weight array
#  @code
#  ds = ... ## dataset as structured array
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
## Get the effective number of entries for 1D-array
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
#  All floating fields  are calculated as
#  \f[ x = \frac{x - \left\langle x \right\rangle}{\sigma} \f]
#  where \f$ \left\langle x \right\rangle\f$ is a mean  value
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
## lower-case shortcuts:
Keys = {    
    'KS' : ( 'ks' , 'kolmogorov' , 'kolmogorovsmirnov' ) , 
    'K'  : ( 'k'  , 'kuiper'                           ) , 
    'AD' : ( 'ad' , 'anderson'   , 'andersondarling'   ) , 
    'CM' : ( 'cm' , 'cramer'     , 'cramervonmises'    ) , 
    'ZK' : ( 'zk' , 'zhangk'     , 'zhangzk'           ) , 
    'ZA' : ( 'za' , 'zhanga'     , 'zhangza'           ) , 
    'ZC' : ( 'zc' , 'zhangc'     , 'zhangzc'           ) , 
    }
# =============================================================================
assert Labels.keys() == Keys.keys() , "Mismatch between Labels & Keys structures!"
# =============================================================================
## clip p-value 
def clip_pvalue ( pvalue , clip = 0.5 ) :
    """ Clip p-value
    """
    pv = VE ( pvalue )
    ## everything is fine, no need to clip 
    if 0 < pv.value() < 1 : return pv  
    ##
    clip = min ( 0.5 , abs ( clip ) * pv.error () )
    ## 
    if   1 <= pv.value() : pv = VE ( 1 - clip , pv.cov2() )
    elif 0 >= pv.value() : pv = VE (     clip , pv.cov2() )
    ## 
    return pv 
    
# =============================================================================
## @class PERMUTATOR
#  Helper class that allow to run permutation test in parallel 
class PERMUTATOR(object) :
    """ Helper class that allow to run permutation test in parallel 
    """
    def __init__ ( self, gof, t_value , ds1 , ds2 ) :
        
        self.gof     = gof
        self.ds1     = ds1
        self.ds2     = ds2
        self.t_value = t_value
        self.__ecdf  = None
        
    # =========================================================================
    ## run N-permutations 
    def __call__ ( self , N , silent = True ) :
        
        counter, tvalues = self.run_toys ( N = N , silent = silent )
        
        if not self.__ecdf : self.__ecdf = Ostap.Math.ECDF ( tvalues , True )
        else               : self.__ecdf.add ( data2vct ( tvalues )  )
        
        return counter 

    # =========================================================================
    ## run N-toys 
    def run_toys ( self, N , silent = True ) :
        """ Run N-toys
        """
        
        numpy.random.seed()
        n1      = len ( self.ds1 )
        pooled  = numpy.concatenate ( [ self.ds1 , self.ds2 ] )
        counter = EffCounter()
        tvalues = [] 
        for i in progress_bar ( N , silent = silent , description = 'Permutations:') : 
            numpy.random.shuffle ( pooled )            
            tv       = self.gof.t_value ( pooled [ : n1 ] , pooled [ n1: ] )
            tvalues.append ( float ( tv ) )
            counter += bool ( self.t_value < tv  )            
        del pooled
        return counter, tuple ( tvalues )

    @property 
    def ecdf ( self ) :
        """`ecdf` : empirical CDF for t-values from permutations 
        """
        return self.__ecdf
    
    @ecdf.setter
    def ecdf ( self , value ) :
        self.__ecdf = value
        
# =============================================================================
jl = None 
# =============================================================================
if False : # ==================================================================
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        import joblib as jl
        jl_version = tuple ( int ( i ) for i in jl.__version__.split('.') )
        # =====================================================================
        ## Run NN-permutations in parallel using joblib 
        def joblib_run ( self , NN , silent = True ) :
            """ Run NN-permutations in parallel using joblib
            """
            me    = math.ceil ( memory_enough() ) + 1 
            nj    = min ( 2 * numcpu () + 3 , me ) 
            lst   = splitter ( NN , nj )
            njobs = len ( [ k for k in splitter ( NN , nj ) ] )
            if not silent : logger.info ( 'permutations: #%d parallel subjobs to be used with joblib' % njobs ) 
            ## 
            conf  = { 'n_jobs' : -1 , 'verbose' : 0 }
            if    (1,3,0) <= jl_version < (1,4,0) : conf [ 'return_as' ] = 'generator'           
            elif  (1,4,0) <= jl_version           : conf [ 'return_as' ] = 'unordered_generator'
            ##
            input   = ( jl.delayed (self.run_toys)( N ) for N in lst )
            counter = EffCounter()
            tvalues = () 
            results = jl.Parallel ( **conf ) ( input )
            for r in progress_bar ( results , max_value = njobs , silent = silent , description = 'Permutations:') :
                cnt , tvals = r 
                counter += cnt
                tvalues += tvals 
                ##
                
            if not self.ecdf : self.ecdf = Ostap.Math.ECDF ( tvalues , True )
            else             : self.ecdf.add ( data2vct ( tvalues )  )

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
        njobs = len ( [ k for k in splitter ( NN , nj ) ] )
        if not silent : logger.info ( 'permutations: #%d parallel subjobs to be used with WorkManager' % njobs ) 
        counter = EffCounter()
        tvalues = () 
        ## 
        ## use the bare interface 
        from ostap.parallel.parallel import WorkManager
        with WorkManager ( silent = silent ) as manager : 
            for result in manager.iexecute ( self.run_toys ,
                                             lst           ,
                                             progress    = not silent      ,
                                             njobs       = njobs           ,
                                             description = 'Permutations:' ) :
                cnt , tvals = result 
                counter += cnt
                tvalues += tvals 
        ##
        if not self.ecdf : self.ecdf = Ostap.Math.ECDF ( tvalues , True )
        else             : self.ecdf.add ( data2vct ( tvalues )  )
        ## 
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
    """ Helper class that allow to run toys in parallel 
    """
    def __init__ ( self    , 
                   gof     , * , 
                   t_value ,
                   pdf     ,
                   Ndata   , sample = False ) :
        
        self.gof     = gof
        self.pdf     = pdf
        self.Ndata   = Ndata 
        self.t_value = t_value
        self.sample  = gof.sample 
        self.silent  = gof.silent
        self.__ecdf  = None 
        
    # =========================================================================
    ## run N-toys 
    def __call__ ( self , nToys , silent = True ) :
        """ Run N-toys
        """
        counter , ecdf = self.run_toys ( nToys = nToys , silent = silent )
        return counter 
    
    # =========================================================================
    ## run N-toys 
    def run_toys ( self , nToys , silent = True ) :
        """ Run N-toys
        """
        ROOT.gRandom                     .SetSeed () 
        ROOT.RooRandom.randomGenerator() .SetSeed ()
        
        counter = EffCounter ()
        tvalues = [] 
        for i in range ( nToys ) :
            
            dset     = self.pdf.generate ( self.Ndata , sample = self.sample )
            tv       = self.gof ( self.pdf , dset )
            counter += bool ( self.t_value > tv   ) ## NOTE SIGN HERE!

            tvalues.append ( tv ) 
            if isinstance  ( dset , ROOT.RooDataSet ) : dset.clear () 
            del dset

        tvalues = tuple ( tvalues ) 
        if not self.ecdf : self.__ecdf = Ostap.Math.ECDF ( tvalues  , True ) 
        else             : self.ecdf.add ( data2vct ( tvalues )  )

        return counter, self.ecdf 

    # =========================================================================
    ## Run N-toys in parallel using WorkManager
    def run ( self , nToys  , silent = False ) :
        """ Run N-toys in parallel using WorkManager
        """
        me       = math.ceil ( memory_enough() ) + 1 
        nj       = min ( 2 * numcpu () + 3 , me ) 
        the_list = [ j for j in splitter ( nToys , nj ) ] 
        njobs    = len ( the_list )
        if not silent : logger.info ( 'toys: #%d parallel subjobs to be used' % njobs ) 
        ##
        counter = EffCounter()
        tvalues = ()
        ##        
        ## use the bare interface 
        from ostap.parallel.parallel import WorkManager
        with WorkManager ( silent = silent ) as manager :
            
            for result in manager.iexecute ( self.run_toys ,
                                             the_list      ,
                                             progress     = not silent ,
                                             njobs        = njobs      ,
                                             description  = 'Toys:'    ) :
                
                cnt , ecdf = result
                
                counter += cnt
                if not self.ecdf : self.__ecdf =   ecdf 
                else             : self.ecdf.add ( ecdf )                
         
        return counter 

    @property 
    def ecdf ( self ) :
        """`ecdf` : empirical CDF for t-values from toys/pseudoexperiments 
        """
        return self.__ecdf

# =============================================================================
## Format the row for GoF tables
#  @code
#  tvalue = ...
#  pvalue = ...
#  ecdf   = ... 
#  header , row = format_row ( tvalue = tvalue, pvalue = pvalue , ecdf = ecdf )
#  @endcode
def format_row ( tvalue    = None ,
                 pvalue    = None ,
                 ecdf      = None ,
                 counter   = None , 
                 precision = 4    ,
                 width     = 6    ) :
    """ Format the row for the GoF tables
    >>> tvalue = ...
    >>> pvalue = ...
    >>> ecdf   = ...
    >>> header , row = format_row ( tvalue = tvalue , pvalue = pvalue  , ecdf = ecdf )
    """
    
    has_tvalue  = not tvalue  is None and isinstance ( tvalue  , num_types ) 
    has_pvalue  = not pvalue  is None and isinstance ( pvalue  , VE        ) 
    has_ecdf    = not ecdf    is None and isinstance ( ecdf    , Ostap.Math.ECDF ) 
    has_counter = not counter is None and isinstance ( counter , ( SE , WSE )   )

    if has_ecdf  and not has_counter  :
        counter     = ecdf.counter ()
        has_counter = True 
        
    if has_tvalue and has_pvalue and has_counter :
        
        header = ( 't-value'    ,
                   't-mean'     ,
                   't-rms'      ,
                   't-min/max'  ,                
                   '%s[..]' % times , 'p-value [%]' , '#%s' % greek_lower_sigma ) 
        
        pv         = clip_pvalue  ( pvalue ) 
        nsigma     = significance ( pv ) ## convert  it to significance
        
        mean       = counter.mean   ()
        rms        = counter.rms    () 
        vmin, vmax = counter.minmax () 
        
        mxv        = max ( abs ( tvalue       ) ,
                           abs ( mean.value() ) ,
                           mean.error()         , rms ,
                           abs ( vmin )  , abs ( vmax ) ) 
        
        fmt, fmtv , fmte , expo = fmt_pretty_ve ( VE ( mxv ,  mean.cov2() ) ,
                                                  precision   = precision   ,
                                                  width       = width       , 
                                                  parentheses = False       )
        
        if expo : scale = 10**expo
        else    : scale = 1
    
        vs  = tvalue / scale
        vm  = mean   / scale
        vr  = rms    / scale
        vmn = vmin   / scale
        vmx = vmax   / scale
        
        fmt2 = '%s/%s' % ( fmtv , fmtv ) 
        
        pvalue  = pvalue * 100
        pvalue  = '%5.2f %s %.2f' % ( pvalue.value() , plus_minus , pvalue.error () )
        nsigma  = '%.2f %s %.2f'  % ( nsigma.value() , plus_minus , nsigma.error () )
        
        row = ( fmtv  % vs ,
                fmt   % ( vm.value() , vm.error() ) ,
                fmtv  % vr                          ,
                fmt2  % ( vmn , vmx )               ,
                ( '%s10^%+d' %  ( times , expo )  if expo else '' ) , pvalue , nsigma )
        
        return header , row

    elif has_tvalue and has_pvalue :

        
        header = ( 't-value'  , '%s[..]' % times , 'p-value [%]' , '#%s' % greek_lower_sigma ) 

        pv        = clip_pvalue  ( pvalue )
        nsigma    = significance ( pv     )
        tv , expo = pretty_float ( tvalue , precision = precision , width = width )
            
        pvalue  = pvalue * 100
        pvalue  = '%5.2f %s %.2f' % ( pvalue.value() , plus_minus , pvalue.error() )        
        nsigma  =  '%.2f %s %.2f' % ( nsigma.value() , plus_minus , nsigma.error () )        
        row     = tv , '%s10^%+d' % ( times , expo ) if expo else '' , pvalue , nsigma

        return header , row 

    elif has_tvalue :
        
        header    = ( 't-value'  , '%s[..]' % times )          
        tv , expo = pretty_float ( tvalue , precision = precision , width = width )
        row       = tv , '%s10^%+d' % ( times , expo ) if expo else ''        
        return header, row

    ## no data for table 
    return () , ()

# ==========================================================================
## Get results in form of the table 
def format_table ( tvalue    = None ,
                   pvalue    = None ,
                   ecdf      = None ,
                   counter   = None ,
                   precision = 4    ,
                   width     = 6    , 
                   title     = ''   ,
                   prefix    = ''   ,
                   style     = ''   ) :
    """ Get results in form of the table 
    """
    
    header , row = format_row ( tvalue    = tvalue    ,
                                pvalue    = pvalue    ,
                                ecdf      = ecdf      ,
                                counter   = counter   , 
                                precision = precision ,
                                width     = width     ) 
    rows = [ header ]
    rows.append ( row )
    rows = T.remove_empty_columns ( rows ) 
    return T.table ( rows               ,
                     title     = title  ,
                     prefix    = prefix ,
                     alignment = 10*'c' ,
                     style     = style  )

# =============================================================================
## Draw ECDF + 2 lines when/if t-value specified
#  @code
#  ecdf   = ...
#  tvalue = ...
#  result = draw_ecdf ( ecdf , tvalue = tvalue ) 
#  @endcode 
def draw_ecdf ( ecdf , tvalue = None , option = '' , **kwargs ) :
    """ Draw ECDF + 2 lines when/if t-value specified
    >>> ecdf   = ...
    >>> tvalue = ...
    >>> result = draw_ecdf ( ecdf , tvalue = tvalue ) 
    """

    has_tvalue = not tvalue is None and isinstance ( tvalue , num_types )
    
    xmin , xmax = ecdf.xmin () , ecdf.xmax ()
    
    if has_tvalue :
        tvalue  = float ( tvalue ) 
        xmin    = min ( xmin , tvalue )
        xmax    = max ( xmax , tvalue )

    xmin , xmax = axis_range ( xmin , xmax , delta = 0.20 )

    ## some transformation  
    kw = cidict ( transform = cidict_fun , **kwargs )

    kw [ 'xmin'      ] = kw.pop ( 'xmin'       , xmin   ) 
    kw [ 'xmax'      ] = kw.pop ( 'xmax'       , xmax   )
    kw [ 'color'     ] = kw.pop ( 'linecolor'  , Orange )
    kw [ 'linewidth' ] = kw.pop ( 'linewidth'  , 2      )
    kw [ 'maxvalue'  ] = kw.pop ( 'maxvalue'   , 1.1    )
    
    result = ecdf.draw  ( option = option , **kw )
    
    ## draw ECDF 
    if not  has_tvalue : return result
    
    ## vertical line 
    vline     = ROOT.TLine ( tvalue , 1e-3 , tvalue , 1 - 1e-3 )
    
    ## horisontal line 
    xmin      = kw [ 'xmin' ]
    xmax      = kw [ 'xmax' ]
    dx        = ( xmax - xmin ) / 100 
    e         = ecdf ( tvalue )
    hline     = ROOT.TLine ( xmin + dx , e , xmax - dx , e )

    ## 
    vline.SetLineWidth  ( 4     ) 
    vline.SetLineColor  ( Green )
    ## 
    hline.SetLineWidth  ( 2     ) 
    hline.SetLineColor  ( Blue  )  
    hline.SetLineStyle  ( 9     ) 
    ##
    vline.draw ( 'same' )
    hline.draw ( 'same' )
    ## 
    return result, vline, hline 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    if not jl : logger.warning ("Joblib is not avalaille!") 

# =============================================================================
##                                                                      The END 
# =============================================================================


    
