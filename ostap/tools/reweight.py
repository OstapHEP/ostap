#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with utilities for reweigting 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-05-10
# =============================================================================
"""Module with utilities for reweighting"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'Weight'          , ## the weighting object  
    'makeWeights'     , ## helper function for   weighting iterations 
    'WeightingPlot'   , ## helper object to define the weighting rule & target
    'ComparisonPlot'  , ## helper object to keep comparision plot
    'W2Tree'          , ## helper to add the calculated weight to ROOT.TTree
    'W2Data'          , ## helper to add the clacualted weight to ROOT.RooAbsData
    ) 
# =============================================================================
from   ostap.core.pyrouts     import VE, SE
from   ostap.math.base        import iszero
from   ostap.core.ostap_types import string_types, list_types, num_types  
from   ostap.math.operations  import Mul as MULT  ## needed for proper abstract multiplication
import ostap.io.zipshelve     as     DBASE ## needed to store the weights&histos
from   ostap.trees.funcs      import FuncTree, FuncData ## add weigth to TTree/RooDataSet
from   ostap.logger.utils     import pretty_ve 
import ostap.histos.histos 
import ostap.histos.compare 
import ostap.trees.trees
import ostap.fitting.dataset
import ROOT, operator
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger    import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.tools.reweight' )
else                       : logger = getLogger ( __name__               )
# =============================================================================
## @class AttrGetter
#  simple class to bypass <code>operator.attrgetter</code> that
#  has some problem with multiprocessing
class AttrGetter(object):
    """Simple class to bypass operator.attrgetter that
    has some problem with multiprocessing
    """
    def __init__ ( self , *attributes ) :
        self.__attributes = attributes 
    def __call__ ( self , obj ) :
        getter = operator.attrgetter( *self.__attributes )
        return getter ( obj )

    @property 
    def attributes ( self ) :
        """``attributes'': the actual attributes
        """
        return self.__attributes
    
# =============================================================================
## @class Weight
#  helper class for semiautomatic reweighting of data 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-05-10
class Weight(object) :
    """Helper class for semi-automatic reweighting of data

    It reads various ``weigth components'' from DBASE and calculates
    the ``global'' event weight via simple accumulation of weights
    
    Simplest case: calculate weights for one variable.
    
    factors = [ ( fun , name ) ] 
    where ``fun'' is a function to get variable from TTree/RooDataSet 
    and   ``name'' is address in data base with reweighting information 
    e.g.
    
    >>> factors = [ ( lambda s : s.pt  , 'pt-data' ) ]
    
    It assumes, that DBASE has an entry with name 'pt-data', such
    that
    
    >>> funcs = database['pt-data']
    
    Here ``func'' will be function/callable or list of functions/callables 
    with weights:
    
    >>> pt = fun ( entry ) 
    >>> w  = 1.0    
    >>> w *= func ( pt )              ## for the single function/callable
    >>> for f in funcs : w*= f ( pt ) ## for list of functions/callables  
    
    quantity ``w'' will be an event  weight
    
    if list of ``factors'' constains more than one entry,
    to each entry weigth is calculated independently and total weight
    is a product of individual weights.

    For simplistic case one can put some storable function or function object
    >>> db['pt-data'] = math.exp
    >>> db['pt-data'] = [ math.exp , math.cosh ] ## List of object also works
    
    For realistic case the weights usually stored as (a list of) 1D or 2D histograms
    >>>  db['pt-data'] = [h1,h2,h3,...,hn]
    where h1,...,hn are iteratively obtained histograms with corrections, such as
    the total correction is calculated by the product of them
    """
    ## constructor of  
    def __init__ ( self                   ,
                   dbase   = "weights.db" , ## the name of data base with the weights 
                   factors = []           ) :
        
        #
        ## make some statistic
        #
        self.__counter = SE ()
        self.__nzeroes = 0 

        self.__vars    = [] 
        if not factors : return
        self.__dbase = dbase
        
        ## open database


        self.__table = [ ( 'Reweighting' , 'accessor' , '#' , 'merged?' , 'skip') ] 
        rows = []
        
        with DBASE.open ( dbase , 'r' ) as db : ## READONLY

            logger.debug ( 'Reweigting database: \n%s' % db.table ( prefix = '# ' ) ) 
                
            ## loop over the weighting factors and build the function
            for wvar in factors :

                funval  = wvar.accessor  ## accessor to the variable 
                funname = wvar.address   ## address  in database 
                merge   = wvar.merge     ## merge sequence of callables?
                skip    = wvar.skip      ## skip   some of them?


                row = []
                
                row.append ( funname ) 

                if isinstance ( funval , str ) :
                    row.append ( funval ) 
                    ## funval = operator.attrgetter( funval ) 
                    funval = AttrGetter( funval )
                elif isinstance ( funval , AttrGetter ) :
                    atts = funval.attributes
                    if 1 == len ( atts ) : atts = atts[0]
                    row.append ( str ( atts ) )
                else :
                    row.append ( '' ) 

                ## 
                functions  = db.get ( funname , [] ) ## db[ funname ]
                if not functions :
                    logger.warning ( "No reweighting is available for ``%s'', skip it" % funname )
                    continue
                                
                if not isinstance (  functions , ( list , tuple ) ) :
                    functions = [ functions ]                    
                
                flen = len(functions) 
                if   0 < skip and skip      < flen :
                    logger.info  ("Use only %d first iterations for ``%s''" % ( skip , funname ) )
                    functions = functions[:skip] 
                elif 0 > skip and abs(skip) < flen :
                    logger.info  ("Skip last %d iterations for ``%s''" % ( skip , funname ) )
                    functions = functions[:skip] 
                elif 0 == skip :
                    pass
                else :
                    logger.error("Invalid ``skip'' parameter %s/%d for ``%s''" % ( skip , flen , funname ) )
                row.append ( '%d' % flen ) 
                
                ## nullify the uncertainties except for the last histogram
                _functions = []
                _first     = True 
                for f in reversed ( functions ) :
                    if isinstance ( f , ROOT.TH1 ) and _first : 
                        ff = f.clone()
                        for i in ff :
                            v     = float ( ff[i] )
                            ff[i] = VE(v,0)
                        _functions.append ( ff  )                        
                        _first = False 
                    else :
                        _functions.append ( f  )
                        
                _functions.reverse() 
                functions = _functions

                row.append  ( '+' if merge else '-' )
                row.append  ( '%s' % skip           )
                
                ## merge list of functions into single function 
                if merge and 1 < len ( functions)  : 
                            
                    ## single_func = functions[0] * functions [1] 
                    single_func = MULT ( functions[0] , functions [1] )
                    
                    for fun in functions [2:] :

                        ## multiply it                               
                        ## single_func *= fun
                        single_func = MULT ( single_func , fun )
                            
                    functions  = [ single_func ]
                    
                self.__vars += [ ( funname , funval , functions , SE() ) ]

                self.__table.append ( row )
                
        self.__vars = tuple ( self.__vars ) 

    @property
    def dbase ( self ) :
        """``dbase'' : th ename of the database woith reweigting information
        """
        return self.__dbase
    
    @property
    def variables ( self ) :
        """```variables'' - the structure of reweigjhting variables/machinery"""
        return self.__vars 
        
    @property 
    def stat    ( self ) :
        """``stat'' : get the statistic of used weights"""
        return self.__counter
    @property 
    def nZeroes ( self ) :
        """``nZeroes'': Number of null weights"""
        return self.__nzeroes
    
    ## calculate the weight for the given "event"
    def __call__ ( self , s ) :
        """   Calculate the weigth for the given ``event'' (==record in TTree/TChain or RooDataSet):
        >>> weight = Weight ( ... )
        >>> tree   = ...
        >>> w = weight ( tree )
        """

        ## initialize the weight 
        weight  = VE(1,0) 

        ## loop over functions 
        for i in self.__vars :
            
            funval    = i[1] ## accessor 
            functions = i[2] ## the functions 

            ##  get the weight arguments for the given event 
            v       = funval ( s )

            ww = VE(1.0)
            
            for f in functions :

                if isinstance ( v , tuple ) : w = f ( *v )
                else                        : w = f (  v )

                ww *= w # update the weight factor 

            ## keep the statistics
            cnt  = i[3]
            cnt += ww.value()

            ## update the global weight 
            weight *= ww
            
        vw = weight.value()
        
        self.__counter += vw 
        if iszero ( vw ) : self.__nzeroes += 1
            
        return vw

    ## format information as a table 
    def table ( self , title = None , prefix = ' ' ) :
        """Format weighting information as table
        """
        import ostap.logger.table as Table
        if title is None : title = 'Weighter(%s)' % self.__dbase 
        return Table.table ( self.__table , title = title , prefix = prefix ,
                             alignment = 'llcc' )

    # ===============================================================================
    ## get the dictionary of graphs
    #  for the regular case, it will be graphs that illustrates the convergency
    #  for reweighting iterations 
    def graphs ( self ) :
        """Get the dictionary of graphs
        """
        grphs = {}
        
        import ostap.histos.graphs 
        with DBASE.open ( self.dbase , 'r' ) as db : ## READONLY
            
            for i in self.__vars :
                
                address   = i[0]
                functions = db.get ( address , () )
                if not functions : continue
                
                graph = ROOT.TGraphAsymmErrors ( len ( functions )  )
                for n , w  in enumerate ( functions ) :
                    if not hasattr ( w  ,'stat' ) : continue

                    cnt = w.stat() 
                    wmn , wmx =   cnt.minmax()
                    graph [ n ] = n , 0 , 0 , 1 , 1 - wmn , wmx - 1  
                    
                grphs [ address ] = graph
                
        return grphs     
        
    # =========================================================================
    ## @class WeightingVar
    #  Helper class to keep information about single reweighting
    #  - accessor : an accessor function that extracts the variable(s)
    #               from  TTree/TChain/RooDataSet
    #  - address  : the  address in DBASE, where reweighting callable(s) is/are stored
    #  - merge    : merge list of callables from DB into the single callable ?
    #  - skip     : use only certain elements from the list of callables from DBASE
    #
    # Schematic data flow to get the weigth for the given event 
    #  - tree/chain/dataset -> accessor -> database(address) -> weight
    class  Var(object) :
        """Helper class to keep information about singe reweighting
        - ``accessor'' : an accessor function that extracts the variable(s) from  TTree/TChain/RooDataSet
        - ``address''  : the  address in DBASE, where reweigftjnig callable(s) is/are stored
        - ``merge''    : merge list of callables from DB into the single callable ?
        - ``skip''     : use only certain elements from the list of callables from DBASE
        
        Schematic data flow to get the weigth for the given event 
        - tree/chain/dataset -> accessor -> database(address) -> weight
        """
        ##  constructor
        #   @param sccessor  the accessor function:  tree -> variable(s) 
        #   @param address   the address of   reweighintg object in DBASE 
        #   @param merge     merge sequence of reweigthing objects ?
        #   @param skip      skip some reweigting objects ?
        def __init__ ( self ,
                       accessor         ,   ## accessor function:  tree -> variable(s) 
                       address          ,   ## the address of   reweighting object in DBASE 
                       merge     = True ,   ## merge sequence of reweighting objects ?
                       skip      = None ) : ## skip some reweigting objects ? 
            """Keep information about single reweighting
            - ``accessor'' : an accessor function that extracts the variable(s) from  TTree/TChain/RooDataSet
            - ``address''  : the  address in DBASE, where reweigftjnig callable(s) is/are stored
            - ``merge''    : merge list of callables from DB into the single callable ?
            - ``skip''     : use only certain elements from the list of callables from DBASE
            
            Schematic data flow to get the weigth for the given event 
            - tree/chain/dataset -> accessor -> database(address) -> weight
            """

            if   isinstance ( accessor , string_types ) :
                ## accessor = operator.attrgetter (  accessor )
                accessor = AttrGetter (  accessor )
            elif isinstance ( accessor , list_types   ) and accessor and \
                     all ( isinstance ( i , string_types ) for i in accessor ) : 
                ## accessor = operator.attrgetter ( *accessor )
                accessor = AttrGetter ( *accessor )
            
            assert callable ( accessor ) , \
                   "Invalid type of ``accessor'' %s/%s" % ( accessor , type( accessor ) )
            
            self.__accessor = accessor ,
            self.__address  = str ( address )
            self.__merge    = True if merge else False
            self.__skip     = skip if skip  else 0
            
        @property
        def accessor ( self ) :
            """``accessor'' - the accessor function to get the variable from TTree/TChain/RooDataSet
            e.g.  to get the variable ``x'' : 
            >>> xvar   = lambda s : s.x
            get to     variables for 2D-reweighting:
            >>> xyvars = lambda s : s.x,s.y         
            """
            return self.__accessor[0]
        
        @property
        def address  ( self ) :
            """``address''  - the address in DB with the reweighting information
            for the  given varibale(s).
            A callable (or list of  callables) is expected in database
            These callables accept their result of ``accessor'' as their argument
            """
            return self.__address
        
        @property
        def merge ( self ) :
            """``merge'' - merge list  of callables into the single  callable?"""
            return self.__merge
        
        @property
        def skip ( self ) :
            """``skip''  - use only certain elements from the list of callables from DBASE
            0 < skip ?   - use only the first  ``skip'' elements 
            0 > skip ?   - skip last  ``abs(skip)'' elements        
            """
            return self.__skip 

# =============================================================================
## default function to make MC projection
#  @param daatset dataset MC  dataset (typically TTree)
#  @param histo   histogram template
#  @param what    historgam varibales
#  @param how     historgam template 
def mc_data_projector  ( dataset , histo , what , how ) :
    """Default function to make MC projection
    """
    
    dataset.project ( histo , what , how )
    return histo 


# =============================================================================
## @class WeightingPlot
#  helper class to manage/keep ``weighting-plot'' 
class WeightingPlot(object) :    
    """Helper class to manage/keep ``weighting-plot''
    
    - ``what'' : the variable/expression to be used for ``weighting-plot''
    Used  as the second argument of ``dataset.project'' method to produce
    ``weighted-plot'':
    
    >>> dataset.project ( mchisto , WHAT , how , ... ) 
    
    - ``how'' : the ``weight'' expression to be used for ``weighed-plots''
    Used  as the third argument of ``dataset.project'' method to produce
    ``weighted-plot'':
    
    >>> dataset.project ( mchisto , what , HOW , ... )
    
    Typically it refers to ``weight'' variable in datasete
    
    - ``address'' : the addres in ``weighting-database''
    where to store the obtained weights
    
    - ``data'' : the ``data'' object, or  the ``target'' for the reweighting procedure
    Typically it is a histogram. But it could be any kind of callable 
    
    - ``mchisto'' : template/shape for the mc-historgam, to be used for reweighting.
    It is used as the  first argument of ``dataset.project'' method
    
    >>> dataset.project ( MCHISTO , what , how , ... )             
    """
    def __init__ ( self              ,
                   what              ,  
                   how               ,
                   address           ,
                   data              ,
                   mc_histo  = None  ,
                   w         = 1.0   ,
                   projector = None  ,
                   ignore    = False ) :
        """Helper class to manage/keep ``weighting-plot''
        
        - ``what'' : the variable/expression to be used for ``weighting-plot''
        Used  as the second argument of ``dataset.project'' method to produce
        ``weighted-plot'':
        
        >>> dataset.project ( mchisto , WHAT , how , ... ) 
        
        - ``how'' : the ``weight'' expression to be used for ``weighed-plots''
        Used  as the third argument of ``dataset.project'' method to produce
        ``weighted-plot'':
        
        >>> dataset.project ( mchisto , what , HOW , ... )
        
        Typically it refers to ``weight'' variable in datasete
        
        - ``address'' : the addres in ``weighting-database''
        where to store the obtained weights
        
        - ``data'' : the ``data'' object, or  the ``target'' for the reweighting procedure
        Typically it is a histotgram. But it could be any kind of callable 
        
        - ``mc_histo'' : template/shape for the mc-historgam, to be used for reweighting.
        It is used as the  first argument of ``dataset.project'' method
        
        >>> dataset.project ( MCHISTO , what , how , ... )         
        """
        
        self.__what      = str(what)     if isinstance ( what , str ) else what 
        self.__how       = str(how )     if isinstance ( how  , str ) else how 
        self.__address   = str(address) 
        self.__data      = data
        self.__mc        = mc_histo      if mc_histo else data.clone()
        assert isinstance ( self.__mc , ROOT.TH1 ), \
               "WPlot: invalid type of ``mchisto'' %s/%s"  % ( self.__mc , type ( self.__mc ) )
        self.__w         = w

        if not projector : projector = mc_data_projector
            
        self.__projector = projector
        
        assert self.projector and callable ( self.projector ) ,"``Projector'' must be callable!"
        self.__ignore    = True if ignore else False 
        
    @property
    def  what ( self ) :
        """``what'' : the variable/expression to be used for ``weighting-plot''
        Used  as the second argument of ``dataset.project'' method to produce
        ``weighted-plot'':
        >>> dataset.project ( mchisto , WHAT , how , ... ) 
        """
        return self.__what
    @property
    def  how  ( self ) :
        """``how'' : the ``weight'' expression to be used for ``weighed-plots''
        Used  as the third argument of ``dataset.project'' method to produce
        ``weighted-plot'':
        >>> dataset.project ( mchisto , what , HOW , ... )
        Typically it refers to ``weight'' variable in datasete
        """
        return self.__how 
    @property
    def  address ( self ) :
        """``address'' : the address in ``weighting-database''
        where to store the obtained weights
        """
        return self.__address
    @property
    def  data ( self ) :
        """``data'' : the ``data'' object, or  the ``target'' for the reweighting procedure
        Typically it is a histotgram. But it could be any kind of callable 
        """
        return self.__data 
    @property
    def  mc_histo ( self ) :
        """``mc_histo'' : template/shape for the mc-histogram, to be used for reweighting.
        It is used as the  first argument of ``dataset.project'' method 
        >>> dataset.project ( MCHISTO , what , how , ... )         
        """
        return self.__mc
    @property
    def projector ( self ) :
        """``projector'' :  callable function to build MC distribution:
        hmc = projector ( dataset , hmc ) 
        """
        return self.__projector
    
    @property
    def w  ( self )   :
        """``w''  - relative weight (relative importance is this variable)"""
        return self.__w 
    @property
    def ignore ( self ) :
        """``ignore'' : do nto use varibal ein reweights, but use for comparsion"""
        return self.__ignore
    
# =============================================================================
from collections import namedtuple 
ComparisonPlot = namedtuple ( 'ComparisonPlot' , ( 'what' , 'data' , 'mc' , 'weight' ) )

# =============================================================================
## draw comparison plots
def _cmp_draw_ ( self ) :
    """Draw comparison plots
    """
    hd   = self.data
    hmc  = self.mc
    hw   = self.weight 

    if isinstance ( hw , ROOT.TH1 ) :
        if hasattr ( hw     , 'red'   ) : hw .red   ()         
        hw.SetMinimum  ( 0 )
        hw.SetMaximum  ( max ( 1.3 , 1.2 * hw.ymax () ) ) 
        hw.draw ()
        hw.level ( 1.0 , linestyle = 1 , linecolor = ROOT.kRed )

    if isinstance ( hd , ROOT.TH1 ) and isinstance ( hmc , ROOT.TH1 ) :
        if hasattr ( hd     , 'green'  ) : hd .green   () 
        if hasattr ( hmc    , 'blue'   ) : hmc.blue   () 
        rmax   = 1.2 * max ( hd .GetMaximum () , hmc.GetMaximum() )
        hd .SetMinimum ( 0 )
        hmc.SetMinimum ( 0 )        
        scale  = ROOT.gPad.GetUymax() / rmax
        hd  .Scale ( scale  )
        hmc .Scale ( scale  )
        hd  .draw  ( 'same' )
        hmc .draw  ( 'same' )
    elif isinstance ( hd  , ROOT.TH1 ) :
        if hasattr ( hd     , 'green'  ) : hd .green   () 
        rmax   = 1.2 *       hd .GetMaximum ()         
        hd .SetMinimum ( 0 )
        scale  = ROOT.gPad.GetUymax() / rmax
        hd  .Scale ( scale  )
        hd  .draw  ( 'same' )
    elif isinstance ( hmc , ROOT.TH1 ) :
        if hasattr ( hmc    , 'blue'   ) : hmc.blue   () 
        rmax   = 1.2 *       hmc.GetMaximum ()         
        hmc.SetMinimum ( 0 )        
        scale  = ROOT.gPad.GetUymax() / rmax
        hmc .Scale ( scale  )
        hmc .draw  ( 'same' )
        

ComparisonPlot. draw = _cmp_draw_

# =============================================================================
## The main  function: perform one re-weighting iteration 
#  and reweight "MC"-data set to looks as "data"(reference) dataset
#  @code
#  results = makeWeights (
#   dataset           , ## data source to be  reweighted (DataSet, TTree, abstract source)
#   plots             , ## reweighting plots
#   database          , ## datadabse to store/update reweigting results
#   delta             , ## stopping criteria for "mean"    weight variation
#   minmax            , ## stopping criteria for "min/max" weight variation
#   power             , ## effective power to apply to the weigths
#   debug      = True , ## store debuig information in database
#   make_plots = True , ## produce useful comparison plots
#   tag        = 'RW' ) ## tag for better printout 
#  @endcode
#  If <code>make_plots = False</code> function returns the tuple of active reweitings:
#  @code
#  active        = makeWeights ( ... , make_plots = False , ... ) \
#  @endcode
#  otherwise it also returns the tuple of comparison plots 
#  @code 
#  active, cmp_plots = makeWeights ( ... , make_plots = True  , ... )
#  for item in  cmp_plots :
#     what    = item.what
#     hdata   = item.data
#     hmc     = item.mc
#     hweight = item.weight
#  @endcode
#  If no more rewighting iteratios required, <code>active</code> is an empty tuple 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-05-10
def makeWeights  ( dataset                    ,
                   plots      = []            , 
                   database   = "weights.db"  ,
                   compare    = None          , ## comparison function 
                   delta      = 0.01          , ## delta for ``mean''  weight variation
                   minmax     = 0.03          , ## delta for ``minmax'' weight variation
                   power      = None          , ## auto-determination
                   debug      = True          , ## save intermediate information in DB
                   make_plots = False         , ## make plots 
                   tag        = "Reweighting" ) :
    """The main  function: perform one re-weighting iteration 
    and reweight ``MC''-data set to looks as ``data''(reference) dataset
    >>> results = makeWeights (
    ... dataset           , ## data source to be  reweighted (DataSet, TTree, abstract source)
    ... plots             , ## reweighting plots
    ... database          , ## datadabse to store/update reweigting results
    ... delta             , ## stopping criteria for `mean`    weight variation
    ... minmax            , ## stopping criteria for `min/max` weight variation
    ... power             , ## effective power to apply to the weigths
    ... debug      = True , ## store debuig information in database
    ... make_plots = True , ## produce useful comparison plots
    ... tag        = 'RW' ) ## tag for better printout
    
    If `make_plots = False`,  it returns the tuple of active reweitings:
    >>> active        = makeWeights ( ... , make_plots = False , ... )
    
    Otherwise it also returns list of comparison plots 
    >>> active, cmp_plots = makeWeights ( ... , make_plots = True  , ... )
    >>> for item in  cmp_plots :
    ...    what    = item.what
    ...    hdata   = item.data
    ...    hmc     = item.mc
    ...    hweight = item.weight
    
    If no more rewighting iteratios required, <code>active</code> is an empty tuple 
    """

    assert 0 < delta  , "makeWeights(%s): Invalid value for ``delta''  %s" % ( tag , delta  )
    assert 0 < minmax , "makeWeights(%s): Invalid value for ``minmax'' %s" % ( tag , minmax ) 

    from ostap.logger.colorized   import allright , attention , infostr 
    from ostap.utils.basic        import isatty

    nplots  = len ( plots )
    ## if 1 < nplots :
    ##     import  math
    ##     fudge_factor = math.sqrt ( 1.0 / max ( 2.0 , nplots -  1.0 ) )
    ##     delta   = delta  * fudge_factor
    ##     minmax  = minmax * fudge_factor

    ## list of plots to compare 
    cmp_plots  = []
    ## reweighting summary table
    header       = ( 'Reweighting' , 'wmin/wmax' , 'OK?' , 'wrms[%]' , 'OK?' , 'chi2/ndf' , 'ww' , 'exp' )
    
    rows         = {}
    save_to_db = [] 
    ## number of active plots for reweighting
    for wplot in plots  :
        
        what      = wplot.what       ## variable/function to plot/compare 
        how       = wplot.how        ## weight and/or additional cuts 
        address   = wplot.address    ## address in database 
        hdata0    = wplot.data       ## original "DATA" object 
        hmc0      = wplot.mc_histo   ## original "MC"   histogram 
        ww        = wplot.w          ## relative weight
        projector = wplot.projector  ## projector for MC data
        ignore    = wplot.ignore     ## ignore for weigtht building?
        #
        # normalize the data
        #
        hdata = hdata0
        if isinstance ( hdata , ROOT.TH1 ) :  hdata = hdata.density ()

        # =====================================================================
        ## make a plot on (MC) data with the weight
        # =====================================================================
        hmc0 = projector ( dataset , hmc0 , what , how )
        
        st   = hmc0.stat()
        mnmx = st.minmax()
        if  iszero ( mnmx [0] ) :
            logger.warning ( "%s: statistic goes to zero %s/``%s''" % ( tag , st , address ) ) 
        elif mnmx [0] <= 0      :
            logger.warning ( "%s: statistic is negative  %s/``%s''" % ( tag , st , address ) ) 
            
        # =====================================================================
        ## normalize MC
        # =====================================================================
        hmc = hmc0.density() 

        # =====================================================================
        ## calculate  the reweighting factor : a bit conservative (?)
        #  this is the only important line
        # =====================================================================
        
        #  try to exploit finer binning if/when possible
        hboth = isinstance ( hmc   , ROOT.TH1 ) and isinstance ( hdata , ROOT.TH1 ) 
            
        if   hboth and 1 == hmc.dim () and 1 == hdata.dim () and \
               len ( hmc ) >= len( hdata ) :   
            w = ( 1.0 / hmc ) * hdata                                 ## NB! 
        elif hboth and 2 == hmc.dim () and 2 == hdata.dim () and \
                 ( hmc.binsx() >= hdata.binsx() ) and \
                 ( hmc.binsy() >= hdata.binsy() ) :
            w = ( 1.0 / hmc ) * hdata                                 ## NB! 
        elif hboth and 3 == hmc.dim () and 3 == hdata.dim () and \
                 ( hmc.binsx() >= hdata.binsx() ) and \
                 ( hmc.binsy() >= hdata.binsy() ) and \
                 ( hmc.binsz() >= hdata.binsz() ) :
            w = ( 1.0 / hmc ) * hdata                                 ## NB! 
        else                            : 
            w = hdata / hmc                                           ## NB!


        # =====================================================================
        ## scale & get the statistics of weights 
        w   /= w.stat().mean().value()
        cnt  = w.stat()
        #
        mnw , mxw = cnt.minmax()
        wvar  = cnt.rms()/cnt.mean()
        good1 = wvar.value()      <= delta
        good2 = abs ( mxw - mnw ) <= minmax 
        good  = good1 and good2  ## small variance?         
        #

        c2ndf  = 0
        for i in w : c2ndf += w[i].chi2 ( 1.0 )
        c2ndf /= ( len ( w ) - 1 ) 

        ## build  the row in the summary table 
        row = address  ,  \
              '%-5.3f/%5.3f' % ( cnt.minmax()[0]    , cnt.minmax()[1] ) , \
              allright ( '+' ) if good2 else attention ( '-' ) , \
              (wvar * 100).toString('%6.2f+-%-6.2f') , \
              allright ( '+' ) if good1 else attention ( '-' ) , '%6.2f' % c2ndf 

        ## make plots at the start of  each iteration? 
        if make_plots : 
            item = ComparisonPlot ( what , hdata , hmc , w ) 
            cmp_plots.append ( item )
            
        row = tuple ( list ( row ) + [ '%4.3f' % ww if 1 != ww  else '' ] ) 
        
        rows [ address ] = row
        
        #
        ## make decision based on the variance of weights 
        #
        mnw , mxw = cnt.minmax()
        if ( not good )  and ( not ignore ) : ## small variance?
            save_to_db.append ( ( address , ww , hdata0 , hmc0 , hdata , hmc , w ) )

        # =====================================================================
        ## make a comparison (if needed)
        # =====================================================================
        if compare : compare ( hdata0 , hmc0 , address )

    
    active  = tuple ( [ p[0] for p in save_to_db ] )  
    nactive = len ( active )  

    if power and callable ( power ) : 
        eff_exp = power ( nactive ) 
    elif isinstance ( power , num_types ) and 0 < power <= 1.5 :
        eff_exp = 1.0 * power 
    elif 1 == nactive and 1 < len ( plots ) :  
        eff_exp = 0.95 
    elif 1 == nactive  :  
        eff_exp = 1.00 
    else               : 
        eff_exp = 1.10 / max ( nactive , 1 ) 
            
    while database and save_to_db :

        entry = save_to_db.pop() 
        
        address , ww , hd0 , hm0 , hd , hm , weight = entry  

        cnt = weight.stat()
        mnw, mxw = cnt.minmax()

        ## avoid too large or too small  weights
        for i in weight :
            w = weight [ i ]
            if   w.value() < 0.5 :
                weight [ i ] = VE ( 0.5 , w.cov2 () ) 
            elif w.value() > 2.0 :
                weight [ i ] = VE ( 2.0 , w.cov2 () ) 
    
        if 1 < nactive and 1 != ww :
            eff_exp *= ww
            logger.info  ("%s: apply ``effective exponent'' of %.3f for ``%s''" % ( tag , eff_exp  , address ) )

        if 1 != eff_exp and 0 < eff_exp : 
            weight = weight ** eff_exp
            row = list ( rows [ address ] )
            row.append ( '%4.3f' % eff_exp )            
            rows [ address ] = tuple ( row ) 
                        
        with DBASE.open ( database ) as db :
            
            db [ address ] = db.get( address , [] ) + [ weight ]
            
            if debug :
                addr        = address + ':REWEIGHTING'
                db [ addr ] = db.get ( addr , [] ) + list ( entry[2:] )
                
        del hd0, hm0 , hd , hm , weight , entry 

    table = [  header ]
    for   row in   rows  : table.append ( rows[row] )
        
    import ostap.logger.table as Table
    logger.info ( '%s, active:#%d \n%s ' % ( tag , nactive , Table.table ( table , title = tag , prefix = '# ' , alignment = 'lccccccc' ) ) )

    cmp_plots = tuple ( cmp_plots ) 
    return ( active , cmp_plots ) if make_plots else active

# =============================================================================
## @class W2Tree
#  Helper class to add the weight into <code>ROOT.TTree</code>
#  @code
#  w    = Weight ( ... )      ## the weighting object
#  tree = ...                 ## The tree
#  wf   = W2Tree ( w , tree ) ## create the weighting function
#  tree.add_new_branch ( 'weight' , wf ) 
#  @endcode
class W2Tree(FuncTree) :
    """Helper class to add the weight into ROOT.TTree
    >>> w    = Weight ( ... )      ## the weighting object
    >>> tree = ...                 ## The tree
    >>> wf   = W2Tree ( w , tree ) ## create the weighting function
    >>> tree.add_new_branch ( 'weight' , wf )
    """
    def __init__ ( self , weight , tree = None ) :
        
        assert isinstance ( weight , Weight    )                 , 'Wrong type of weight!'
        assert tree is None or isinstance ( tree  , ROOT.TTree ) , 'Wrong type of tree!'
        
        FuncTree.__init__ ( self , tree )
        self.__weight = weight
        
    ## evaluate the weighter for the given TTree entry 
    def evaluate ( self ) : 
        """Evaluate the weighter for the given TTree entry"""
        t = self.tree ()
        w = self.__weight ( t ) 
        return w

    @property
    def weight ( self ) :
        """``weight'' : get the weighter object"""
        return  self.__weight

# =============================================================================
## @class W2Data
#  Helper class to add the weight into <code>RooDataSet</code>
#  @code
#  w    = Weight ( ... )      ## the weighting object
#  ds   = ...                 ## dataset 
#  wf   = W2Data ( w , tree ) ## create the weighting function
#  ds.add_new_var ( 'weight' , wf ) 
#  @endcode
class W2Data(FuncData) :
    """Helper class to add the weight into <code>RooDataSet</code>
    >>> w    = Weight ( ... )      ## the weighting object
    >>> ds   = ...                 ## dataset 
    >>> wf   = W2Data ( w , tree ) ## create the weighting function
    >>> ds.add_new_var ( 'weight' , wf ) 
    """
    
    def __init__ ( self , weight , data = None ) :

        assert isinstance ( weight , Weight    )                      , 'Wrong type of weight!'
        assert data is None or isinstance ( data  , ROOT.RooAbsData ) , 'Wrong type of data!'
        
        FuncData.__init__ ( self , data )
        self.__weight = weight
        
    ## evaluate the weighter for the given RooAbsData entry 
    def evaluate ( self ) : 
        """Evaluate the weighter for the given RooAbsData entry"""
        
        d = self.data ()
        w = self.__weight ( d ) 
        return w
    
    @property
    def weight ( self ) :
        """``weight'' : get the weighter object"""
        return  self.__weight
    
# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
