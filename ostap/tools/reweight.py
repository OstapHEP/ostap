#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with utilities for iterative reweigting 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-05-10
# =============================================================================
""" Module with utilities for iterative reweighting
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'Weight'            , ## the weighting object  
    'makeWeights'       , ## helper function for   weighting iterations 
    'WeightingPlot'     , ## helper object to define the weighting rule & target
    'ComparisonPlot'    , ## helper object to keep comparision plot
    'W2Tree'            , ## helper to add the calculated weight to ROOT.TTree
    'W2Data'            , ## helper to add the clacualted weight to ROOT.RooAbsData
    'backup_to_ROOT'    , ## backup the reweighting DBASE to ROOT file 
    'restore_from_ROOT' , ## restre the reweighting DBASE from ROOT file 
    ) 
# =============================================================================
from   ostap.core.meta_info   import root_info 
from   ostap.core.pyrouts     import VE, SE, Ostap, hID  
from   ostap.math.base        import iszero, axis_range 
from   ostap.utils.strings    import split_string 
from   ostap.core.ostap_types import string_types, list_types, num_types, sized_types, sequence_types    
from   ostap.math.operations  import Mul  as MULT       ## needed for proper abstract multiplication
from   ostap.trees.funcs      import FuncTree, FuncData ## add weight to TTree/RooDataSet
from   ostap.utils.utils      import CallThem, AttrGetter 
from   ostap.utils.strings    import is_formula
from   ostap.utils.basic      import typename, numcpu  
from   ostap.math.reduce      import root_factory
from   ostap.logger.symbols   import plus_minus, thumb_up , checked_no, checked_yes, number, chi2ndf   
from   ostap.logger.colorized import allright , attention
import ostap.io.zipshelve     as     DBASE              ## needed to store the weights&histos
import ostap.logger.table     as     T 
import ostap.core.core 
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
_new_methods_ = []
# =============================================================================
## converted iterations ?
OK  = thumb_up  
## not yet converging 
NOT = checked_no
## tag in DBASE for list of known reweightings
tag_reweightings = 'REWEIGHTINGS'
## tag in DBASE for collection of comparsoon plots
tag_comparisons  = 'COMPARISON_PLOTS'
## tag in DBASE for convergency graphs
tag_graphs       = 'CONVERGENCY_GRAPHS'
## suffix for reweighting address 
tag_rw_suffix    = 'REWEIGHTING'

# =============================================================================
## Try to normalize  
def normalize_data ( data , another = None  ) :
    """ Make a try to normalize the distribution 
    """
    if   isinstance ( data, ROOT.TH1   ) : return data.density()
    elif hasattr    ( data, 'density'  ) : return data.density()

    ## (1) make a try to use internal integration
    my_integral = None if not hasattr ( data , 'integral' ) else data.integral 

    if my_integral :
        try  :
            ii = my_integral()
            if 0 < ii : return data/ ii
        except TypeError : pass

    ## can we use numerical integration ? 
    num_integral  = callable ( data )
    
    ## (1) make a try to use internal integration 
    if my_integral or num_integral :
        
        if   hasattr ( data    , 'xmin' ) and hasattr ( data    , 'xmax' ) :            
            xmin, xmax = data.xmin()    , data.xmax() 
        elif hasattr ( another , 'xmin' ) and hasattr ( another , 'xmax' ) :
            xmin, xmax = another.xmin() , another.xmax()            
        else :
            
            return data                   ## RETURN!! 

        ## use interal integration 
        if my_integral :            
            try :
                ii = my_integral ( xmin , xmax )
                if 0 < ii : return data / ii
            except TypeError : pass

        ## use numerical integration as the last resort
        if num_integration :
            import ostap.math.integral as II
            ii = II.integral ( data , xmin , xmax )
            if 0 < ii : return data / ii

    ## nothing to do...
    return data 
        
# =============================================================================
## Set the proper axis range for the comparison plots
def adjust_histo_range ( histo , min_zero = True ) :
    """ Set the proper axis range for the comparison plots
    """
    
    vmin, vmax  = histo.minmax ( errors = True )
    vmax        = max ( vmax , 1.05 )
    vmin        = min ( vmin , 0.95 )
    
    if   0.95 <= vmin and vmax <= 1.05 :
        
        histo.SetMinimum  ( 0.95 )
        histo.SetMaximum  ( 1.05 )
        histo.SetContour  ( 50   )
        
    if 0.90 <= vmin and vmax <= 1.10 :
        
        histo.SetMinimum  ( 0.90 )
        histo.SetMaximum  ( 1.10 )
        histo.SetContour  ( 40   )
        
    elif 0.80 <= vmin and vmax <= 1.20 :
        
        histo.SetMinimum  ( 0.80 )
        histo.SetMaximum  ( 1.20 )
        histo.SetContour  ( 40   )
        
    elif 0.75 <= vmin and vmax <= 1.25 :
        
        histo.SetMinimum  ( 0.75 )
        histo.SetMaximum  ( 1.25 )
        histo.SetContour  ( 50   )
        
    elif 0.50 <= vmin and vmax <= 1.50 :
        
        histo.SetMinimum  ( 0.50 )
        histo.SetMaximum  ( 1.50 )
        histo.SetContour  ( 50   )
        
    else :
        
        vmn, vmx = axis_range ( vmin , vmax , delta = 0.20 )
        histo.SetMinimum  ( vmn  )
        histo.SetMaximum  ( vmx  )
        histo.SetContour  ( 50   )

    ## reset minimum again for 1D histograms 
    if min_zero and isinstance ( histo , ROOT.TH1 ) and 1 == histo.dim() : 
        histo.SetMinimum  ( min ( vmin , 0.0 ) )

    return histo 
    
# =============================================================================
## @class Weight
#  helper class for semiautomatic reweighting of data 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-05-10
class Weight(object) :
    """ Helper class for semi-automatic reweighting of data

    It reads various `weight components' from DBASE and calculates
    the `global' event weight via simple accumulation of weights
    
    Simplest case: calculate weights for one variable:
    
    factors = [ ( fun , name ) ] 

    where `fun' is a function to get variable from TTree/RooDataSet 
    and   `name' is address in data base with reweighting information 
    e.g.
    
    >>> factors = [ ( lambda s : s.pt  , 'pt-data' ) ]
    
    It assumes, that DBASE has an entry with name 'pt-data', such that
    
    >>> funcs = database['pt-data']
    
    Here `funcs' will be callable or list of callables 
    with weights:
    
    >>> pt = fun  ( entry ) 
    >>> w  = 1.0    
    >>> w *= func ( pt )              ## for the single function/callable
    >>> for f in funcs : w*= f ( pt ) ## for list of functions/callables  
    
    The quantity `w' will be an event  weight
    
    If list of `factors' constains more than one entry,
    to each entry weigth is calculated independently and total weight
    is a product of individual weights.

    For simplistic case one can put some storable function or function object

    >>> db['pt-data'] = math.exp
    >>> db['pt-data'] = [ math.exp , math.cosh ] ## List of object also works
    
    For realistic/practical cases the weights usually stored as (a list of) 1D or 2D histograms

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
        self.__dbase   = dbase
        self.__factors = tuple ( factors ) 
        
        ## open database

        self.__table = [ ( 'Reweighting' , 'accessor' , number , 'merged?' , 'skip') ] 
        rows = []

        # =============================================================================
        try : # =======================================================================
            with DBASE.open ( dbase , 'r' ) as db : pass
            # =========================================================================
        except Exception as e : # =====================================================
            # =========================================================================
            message = "DBASE `%s' can not be opened (or does not exist)" % dbase 
            logger.error  ( message , exc_info = True )
            raise OSError ( message ) 

        with DBASE.open ( dbase , 'r' ) as db : ## READONLY

            logger.debug ( 'Weight: Reweighting database: \n%s' % db.table ( prefix = '# ' ) ) 
                
            ## loop over the weighting factors and build the function
            for wvar in self.__factors :

                funval  = wvar.accessor  ## accessor to the variable(s) 
                funname = wvar.address   ## address  in database 
                merge   = wvar.merge     ## merge sequence of callables?
                skip    = wvar.skip      ## skip   some of them?

                row = []
                
                row.append ( funname ) 

                if isinstance ( funval , str ) :
                    row.append ( funval ) 
                    funval = AttrGetter( funval )
                ## elif isinstance ( funval , AttrGetter ) :
                ##    atts = funval.attributes
                ##    if 1 == len ( atts ) : atts = atts[0]
                ##    row.append ( str ( atts ) )
                else :
                    row.append ( str ( funval ) ) 

                ## get the reweighting infromation from database 
                functions  = db.get ( funname , [] ) ## db[ funname ]
                if not functions :
                    logger.warning ( "Weight: No reweighting is available for `%s', skip it" % funname )
                    continue
                                
                if not isinstance (  functions , sequence_types ) :
                    functions = [ functions ]                    
                
                flen = len ( functions ) 
                if   0 < skip and skip      < flen :
                    logger.info  ("Weight: Use only %d first iterations for `%s'" % ( skip , funname ) )
                    functions = functions[:skip] 
                elif 0 > skip and abs(skip) < flen :
                    logger.info  ("Weight: Skip last %d iterations for `%s'" % ( skip , funname ) )
                    functions = functions[:skip] 
                elif 0 == skip :
                    pass
                else :
                    logger.error("Weight: Invalid `skip' parameter %s/%d for `%s'" % ( skip , flen , funname ) )
                row.append ( '%d' % flen ) 
                
                ## nullify the uncertainties except for the last histogram
                _functions = []
                _first     = True 
                for f in reversed ( functions ) :
                    if isinstance ( f , ROOT.TH1 ) and _first : 
                        ff = f.clone()
                        for i in ff :
                            v        = float ( ff [ i ] )
                            ff [ i ] = VE ( v , 0 )
                        _functions.append ( ff  )                        
                        _first = False 
                    else :
                        _functions.append ( f   )
                        
                _functions.reverse() 
                functions = _functions

                row.append  ( checked_yes if merge else checked_no )
                if skip : row.append  ( '%s' % skip )
                
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
                
        self.__table = T.remove_empty_columns ( self.__table )
        self.__vars = tuple ( self.__vars ) 

    @property
    def dbase ( self ) :
        """'dbase' : the name of the database for storage of reweigting information
        """
        return self.__dbase
    
    @property
    def variables ( self ) :
        """'variables' - the structure of reweighting variables/machinery"""
        return self.__vars 
        
    @property 
    def stat    ( self ) :
        """'stat' : get the statistic of used weights"""
        return self.__counter
    @property 
    def nZeroes ( self ) :
        """'nZeroes': Number of null weights"""
        return self.__nzeroes

    # ==============================================================================
    ## non-empty weighter?
    #  return True if there is at leats one rweweigjtig recorf in database
    def __bool__    ( self ) :
        """ non-empty weighter?
        - return True if there is at least one weweigjtig recorf in database        
        """
        with DBASE.open ( self.dbase , 'r' ) as db : 
            rws = db.get( tag_reweightings , () )
            if not rws : return False
            for rw in rws :
                r = db.get ( rw , () )
                if r : return True
            return False
    __nonzero__ = __bool__ 

    ## calculate the weight for the given "event"
    def __call__ ( self , s ) :
        """ Calculate the weigth for the given `event' (==record in TTree/TChain or RooDataSet):
        >>> weight = Weight ( ... )
        >>> tree   = ...
        >>> w = weight ( tree )
        """

        ## initialize the weight 
        weight  = VE ( 1 , 0 ) 

        ## loop over functions 
        for i in self.__vars :
            
            funval    = i [ 1 ] ## accessor 
            functions = i [ 2 ] ## the functions 

            ##  get the weight arguments for the given event 
            v         = funval ( s )

            ww = VE ( 1.0 )
            
            for f in functions :

                if isinstance ( v , tuple ) : w = f ( *v )
                else                        : w = f (  v )

                ww *= w # update the weight factor 

            ## keep the statistics
            cnt  = i [ 3 ]
            cnt += ww.value()

            ## update the global weight 
            weight *= ww
            
        vw = weight.value()
        
        self.__counter += vw 
        if iszero ( vw ) : self.__nzeroes += 1
            
        return vw

    ## format information as a table 
    def table ( self , title = None , prefix = ' ' ) :
        """ Format weighting information as table
        """
        import ostap.logger.table as T
        if title is None : title = "Weight('%s')" % self.__dbase 
        return T.table ( self.__table , title = title , prefix = prefix , alignment = 'llcc' )
    
    # ===============================================================================
    ## get the dictionary of graphs
    #  for the regular case, it will be graphs that illustrates the convergency
    #  for reweighting iterations 
    def graphs ( self ) :
        """ Get the dictionary of graphs
        """
        grphs = {}

        import ostap.histos.graphs 
        with DBASE.open ( self.dbase , 'r' ) as db : ## READONLY
            
            for i in self.__vars :
                
                address   = i [ 0 ]
                functions = db.get ( address , () )
                if not functions : continue

                ## graph with the spread(min/max) weights 
                gr1 = ROOT.TGraphAsymmErrors ( len ( functions )  )
                ## graph with the rms weights 
                gr2 = ROOT.TGraphErrors      ( len ( functions )  )

                has_points = False 
                for n , w  in enumerate ( functions , start = 1 ) :
                    
                    if not hasattr ( w  ,'stat' ) : continue
                    
                    cnt         = w.stat()
                    cnt_mean    = cnt.mean  () 
                    wmean       = cnt_mean  .value()
                    wrms        = cnt.rms   ()
                    wmn , wmx   = cnt.minmax()
                    gr1 [ n - 1 ] =      n , 0 , 0 , wmean , wmn - wmean , wmx - wmean  
                    gr2 [ n - 1 ] = VE ( n , 0.0 ) , VE ( wmean , wrms ** 2) 
                    has_points = True

                if has_points :
                    
                    gr1.SetMarkerStyle ( 1  )
                    gr1.SetFillColor   ( ROOT.kOrange )
                    gr1.SetFillStyle   ( 1001 )
                    gr1.SetLineColor   ( ROOT.kOrange )
                    gr1.SetMarkerColor ( ROOT.kOrange )
                    
                    gr2.SetMarkerStyle ( 20 )
                    gr2.SetLineColor   ( 2  )
                    gr2.SetMarkerColor ( 2  )
                    
                    ## combine two graphs together 
                    graph = ROOT.TMultiGraph ()
                    graph.Add      ( gr1 , '3'  )
                    graph.Add      ( gr2 , 'pl' )
                    graph.SetTitle ( '%s;#iteration;weights' % address )
                    
                    grphs [ address ] = graph

        # =====================================================================
        ## store graphs in DATABASE:
        if grphs : # ==========================================================
            # =================================================================
            try : # ===========================================================
                # =============================================================
                with DBASE.open ( self.dbase ) as db :
                    db [ tag_graphs ] = grphs  
                    for key, gr in grphs.items() : db [ '%s:CONVERGENCY' % key ] = gr
                # ============================================================
            except : # =======================================================
                pass
            
        return grphs     

    def __str__    ( self ) : return self.table()
    def __repr__   ( self ) : return self.table()

    # =========================================================================
    ## reduce the object for serialization (for parallel processing)
    def __reduce__ ( self ) :
        """ Reduce the object for serialization (for parallel processing)
        """
        return root_factory , ( type ( self )  ,
                                self.__dbase   ,
                                self.__factors )
                                
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
        """ Helper class to keep information about single reweighting
        - `accessor'  : an accessor function that extracts the variable(s) from  TTree/TChain/RooDataSet
        - `address'   : the address in DBASE, where reweighting  callable(s) is/are stored
        - `merge'     : merge list of callables from DB into the single callable ?
        - `skip'      : use only certain elements from the list of callables from DBASE
        
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
            """ Keep information about single reweighting
            - `accessor' : an accessor function that extracts the variable(s) from  TTree/TChain/RooDataSet
            - `address'  : the  address in DBASE, where reweigftjnig callable(s) is/are stored
            - `merge'    : merge list of callables from DB into the single callable ?
            - `skip'     : use only certain elements from the list of callables from DBASE
            
            Schematic data flow to get the weigth for the given event 
            - tree/chain/dataset -> accessor -> database(address) -> weight
            """

            if isinstance ( accessor , string_types ) :
                accessor = split_string ( accessor , strip = True , respect_groups = True )
                if 1 == len ( accessor ) : accessor = accessor [ 0 ]

            if   accessor and isinstance ( accessor , string_types ) and is_formula ( accessor ) :
                ## use exporession here (suitable both for for TTree/RooAbsData)
                accessor = Ostap.Functions.Expression ( accessor )
                
            elif accessor and isinstance ( accessor , string_types ) :
                ## accessor = operator.attrgetter (  accessor )
                accessor = AttrGetter (  accessor )
                
            elif accessor and isinstance ( accessor , sequence_types   ) and \
                     all ( callable ( i ) for i in accessor ) :
                ## merge into single callable 
                accessor = CallThem ( accessor )
                
            elif accessor and isinstance ( accessor , sequence_types   )       and \
                     all ( isinstance ( i , string_types ) for i in accessor ) and \
                     any ( is_formula ( i ) for i in accessor  ) :

                ## convert to universal function for TTree/RooAbsData 
                accessor = tuple ( Ostap.Functions.Expression ( i ) for i in accessor )
                ## merge into single callable                 
                accessor = CallThem ( accessor )
                
            elif accessor and isinstance ( accessor , sequence_types   )       and \
                     all ( isinstance ( i , string_types ) for i in accessor ) : 
                
                ## accessor = operator.attrgetter ( *accessor )
                accessor = AttrGetter ( *accessor )

            assert accessor and callable ( accessor ) , \
                   "Invalid type of `accessor' %s/%s" % ( accessor , typename ( accessor ) )
            
            self.__accessor = accessor ,
            self.__address  = str ( address )
            self.__merge    = True if merge else False
            self.__skip     = skip if skip  else 0
            
        @property
        def accessor ( self ) :
            """`accessor' - the accessor function to get the variable from TTree/TChain/RooDataSet
            e.g.  to get the variable `x' : 
            >>> xvar   = lambda s : s.x
            get to     variables for 2D-reweighting:
            >>> xyvars = lambda s : s.x,s.y         
            """
            return self.__accessor[0]
        
        @property
        def address  ( self ) :
            """`address'  - the address in DB with the reweighting information
            for the  given variables(s).

            A callable (or list of  callables) is expected in database

            These callables accept their result of `accessor' as their argument
            """
            return self.__address
        
        @property
        def merge ( self ) :
            """`merge' - merge list  of callables into the single  callable?"""
            return self.__merge
        
        @property
        def skip ( self ) :
            """`skip'   - use only certain elements from the list of callables from DBASE
            0 < skip ?  - use only the first  `skip' elements 
            0 > skip ?  - skip last  `abs(skip)' elements        
            """
            return self.__skip 

# =============================================================================
## default function to make MC projection
#  @param daatset dataset MC  dataset (typically TTree)
#  @param histo   histogram template
#  @param what    histogram varibales
#  @param how     histogram template 
def mc_data_projector  ( dataset , histo , what , how ) :
    """ Default function to make MC projection
    """
    dataset.project ( histo , what , cuts = how )
    return histo 

# =============================================================================
## @class WeightingPlot
#  helper class to manage/keep `weighting-plot' 
class WeightingPlot(object) :    
    """ Helper class to manage/keep `weighting-plot'
    
    - `what' : the variable/expression to be used for `weighting-plot'
    Used  as the second argument of `dataset.project' method to produce
    `weighted-plot':
    
    >>> dataset.project ( mchisto , WHAT , how , ... ) 
    
    - `how' : the `weight' expression to be used for `weighed-plots'
    Used  as the third argument of `dataset.project' method to produce
    `weighted-plot':
    
    >>> dataset.project ( mchisto , what , HOW , ... )
    
    Typically it refers to `weight' variable in the dataset
    
    - `address' : the addres in `weighting-database'
    where to store the obtained weights
    
    - `data' : the `data' object, or  the `target' for the reweighting procedure
    Typically it is a histogram. But it could be any kind of callable 
    
    - `mchisto' : template/shape for the mc-histogram, to be used for reweighting.
    It is used as the  first argument of `dataset.project' method
    
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
        """ Helper class to manage/keep `weighting-plot'
        
        - `what' : the variable/expression to be used for `weighting-plot'
        Used  as the second argument of `dataset.project' method to produce
        `weighted-plot':
        
        >>> dataset.project ( mchisto , WHAT , how , ... ) 
        
        - `how' : the `weight' expression to be used for `weighed-plots'
        Used  as the third argument of `dataset.project' method to produce
        `weighted-plot':
        
        >>> dataset.project ( mchisto , what , HOW , ... )
        
        Typically it refers to `weight' variable in dataset
        
        - `address' : the addres in `weighting-database'
        where to store the obtained weights
        
        - `data' : the `data' object, or the `target' for the reweighting procedure
        Typically it is a histotgram. But it could be any kind of callable 
        
        - `mc_histo' : template/shape for the mc-histogram, to be used for reweighting.
        It is used as the  first argument of `dataset.project' method
        
        >>> dataset.project ( MCHISTO , what , how , ... )         
        """

        assert data and callable ( data ) , "WeightingPlot: `data' object must be callable!"
        
        self.__what      = str(what)     if isinstance ( what , str ) else what 
        self.__how       = str(how )     if isinstance ( how  , str ) else how 
        self.__address   = str(address)

        ## (1) make a try to normalize input data
        ndata = normalize_data  ( data ) 
        if not ndata is data :
            logger.info ( "WeightingPlot: DATA converted to normalized `density-like' form for %s/%s" % ( self.what ,
                                                                                                          self.address ) )
            data = ndata
            
        self.__data = data 
        ## (2) check for non-positive data when/if possible  
        if isinstance ( data , ROOT.TH1 ) :
            st     = data.stat()
            mn, mx = st.minmax()
            if iszero ( mn ) or mn < 0 :
                logger.warning ("Weighting(%s,%s): DATA Statistics is non-positive: min/max=%.3g/%.3g : " % ( self.what    ,
                                                                                                              self.address ,
                                                                                                              float ( mn ) ,
                                                                                                              float ( mx ) ) )
                
        self.__mc       = mc_histo if mc_histo else data.clone()
        self.__w        = w

        assert projector is None or callable ( projector ) , \
            "`Projector' must be None or callable!"        
        self.__projector = projector
        
        self.__ignore    = True if ignore else False 
        
    @property
    def  what ( self ) :
        """`what` : the variable/expression to be used for `weighting-plot'
        Used  as the second argument of `dataset.project' method to produce
        `weighted-plot':
        >>> dataset.project ( mchisto , WHAT , how , ... ) 
        """
        return self.__what
    @property
    def  how  ( self ) :
        """`how` : the `weight' expression to be used for `weighed-plots'
        Used  as the third argument of `dataset.project' method to produce
        `weighted-plot':
        >>> dataset.project ( mchisto , what , HOW , ... )
        Typically it refers to `weight' variable in datasete
        """
        return self.__how 
    @property
    def  address ( self ) :
        """`address` : the address in `weighting-database'
        where to store the obtained weights
        """
        return self.__address
    @property
    def  data ( self ) :
        """`data` : the `data' object, or  the `target' for the reweighting procedure
        Typically it is a histogram. But it could be any kind of callable/function 
        """
        return self.__data 
    @property
    def  mc_histo ( self ) :
        """`mc_histo` : template/shape for the mc-histogram, to be used for reweighting.
        It is used as the  first argument of `dataset.project' method 
        >>> dataset.project ( MCHISTO , what , how , ... )         
        """
        return self.__mc
    @mc_histo.setter
    def mc_histo ( self , value ) :
        self.__mc = value
        
    @property
    def projector ( self ) :
        """`projector` :  callable function to build MC distribution:
        >>> hmc = projector ( dataset , hmc ) 
        """
        return mc_data_projector if self.__projector is None else self.__projector
    
    @property
    def w  ( self )   :
        """`w`  - relative weight (relative importance of this variable)
        - the importance can be (de/in)creased, if needed. 
        """
        return self.__w
    
    @property
    def ignore ( self ) :
        """`ignore` : do not use variable in reweights; use for comparsion only
        - the reweighting results for this variable ar nto saved to database 
        """
        return self.__ignore
    
# =============================================================================
from collections import namedtuple 
ComparisonPlot = namedtuple ( 'ComparisonPlot' ,
                              ( 'what'   ,
                                'data'   ,
                                'mc'     ,
                                'weight' ) )

_store = set() 
# =============================================================================
## draw comparison plots
def _cmp_draw_ ( self ) :
    """ Draw comparison plots
    """
    hd   = self.data
    hmc  = self.mc
    hw   = self.weight 
        
    adjust_histo_range ( hw )

    
    if   isinstance ( hw , ROOT.TH3 ) and 3 == hw.dim () :

        hw.draw  ( copy = True )

        minx, miny, minz = hw.minimum_bin()
        maxx, maxy, maxz = hw.maximum_bin()
        
        ax   = hw.GetXaxis()
        ay   = hw.GetYaxis()
        az   = hw.GetZaxis()

        minv = hw [ minx, miny, minz ]
        maxv = hw [ maxx, maxy, maxz ]

        rows = [ ( '' , 'value' , 'x' , 'y' , 'z' ) ]
        
        x , y , z = ax.GetBinCenter ( minx ) , ay.GetBinCenter ( miny ) , az.GetBinCenter ( minz )        
        row  = 'Min-weight' , '%.2f' % minv , '%.4g' % x , '%.4g' % y , '%.4g' % z
        rows.append ( row )

        x , y , z = ax.GetBinCenter ( maxx ) , ay.GetBinCenter ( maxy ) , az.GetBinCenter ( maxz )        
        row  = 'Max-weight' , '%.2f' % maxv , '%.4g' % x , '%.4g' % y , '%.4g' % z
        rows.append ( row )

        title = 'Comparison plot for %s' % self.what
        table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lcccc' )
        logger.info ( '%s:\n%s' % ( title , table ) )
        
    elif isinstance ( hw , ROOT.TH2 ) and 2 == hw.dim () :

        zmin,zmax = hw.zminmax()

        nx = hw.binsx()
        ny = hw.binsy()

        from ostap.plotting.style import useStyle 
        with useStyle ( 'Z' , textformat = '.2f' ) :
            
            if   10 <  nx and 10 < ny : hw.draw  ( 'cont4z' , copy = True )
            elif 10 <  nx or  10 < ny : hw.draw  ( 'cont4z' , copy = True )
            else                      : hw.texte ( 'colz'   , copy = True )
    
        minx, miny = hw.minimum_bin()
        maxx, maxy = hw.maximum_bin()

        ax   = hw.GetXaxis()
        ay   = hw.GetYaxis()
        
        minv = hw [ minx, miny ]
        maxv = hw [ maxx, maxy ]

        pmin1 = ROOT.TMarker ( ax.GetBinCenter ( minx ) , ay.GetBinCenter ( miny ) , 43 )
        pmax1 = ROOT.TMarker ( ax.GetBinCenter ( maxx ) , ay.GetBinCenter ( maxy ) , 47 )        
        pmin1.SetMarkerColor ( 2 )
        pmax1.SetMarkerColor ( 2 )        
        pmin1.SetMarkerSize  ( 2 )
        pmax1.SetMarkerSize  ( 2 )
        
        pmin2 = pmin1.Clone  () 
        pmax2 = pmax1.Clone  () 
        pmin2.SetMarkerSize  ( pmin1.GetMarkerSize() + 1  )
        pmax2.SetMarkerSize  ( pmax1.GetMarkerSize() + 1  )
        pmin2.SetMarkerColor ( 0 )
        pmax2.SetMarkerColor ( 0 )

        for p in ( pmin2 , pmax2 , pmin1 , pmax1 ) :
            p.DrawClone ()
            _store.add ( p ) 
        
        rows = [ ( '' , 'value' , 'x' , 'y'  ) ]
        
        x , y = ax.GetBinCenter ( minx ) , ay.GetBinCenter ( miny ) 
        row  = 'Min-weight' , '%.2f' % minv , '%.4g' % x , '%.4g' % y 
        rows.append ( row )

        x , y = ax.GetBinCenter ( maxx ) , ay.GetBinCenter ( maxy ) 
        row  = 'Max-weight' , '%.2f' % maxv , '%.4g' % x , '%.4g' % y 
        rows.append ( row )

        title = 'Comparison plot for %s' % self.what
        table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lccc' )
        logger.info ( '%s:\n%s' % ( title , table ) )

    elif isinstance ( hw , ROOT.TH1 ) and 1 == hw.dim () :
        
        if hasattr ( hw , 'yellow' ) : hw .yellow ()         

        hw.draw  ( copy = True )
        hw.level ( 1.0 , linestyle = 1 , linecolor = ROOT.kRed  )

        minx = hw.minimum_bin()
        maxx = hw.maximum_bin()

        ax   = hw.GetXaxis()

        minv = hw [ minx ]
        maxv = hw [ maxx ]

        rows = [ ( '' , 'value' , 'x' ) ]
        
        x = ax.GetBinCenter ( minx ) 
        row  = 'Min-weight' , '%.2f' % minv , '%.4g' % x 
        rows.append ( row )

        x    = ax.GetBinCenter ( maxx ) 
        row  = 'Max-weight' , '%.2f' % maxv , '%.4g' % x 
        rows.append ( row )

        pmin = ROOT.TMarker ( ax.GetBinCenter ( minx ) , 1.0 , 43 )
        pmax = ROOT.TMarker ( ax.GetBinCenter ( maxx ) , 1.0 , 47 )
        
        pmin.SetMarkerColor ( 92 )
        pmax.SetMarkerColor ( 92 )
        pmin.SetMarkerSize  (  3 )
        pmax.SetMarkerSize  (  2 )
        pmin.DrawClone() 
        pmax.DrawClone()

        _store.add ( pmin )
        _store.add ( pmax )
                
        title = 'Comparison plot for %s' % self.what
        table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lcc' )
        logger.info ( '%s:\n%s' % ( title , table ) )
        
    if isinstance ( hd , ROOT.TH1 ) and isinstance ( hmc , ROOT.TH1 ) and 1 == hd.dim() :

        hd  = hd .clone()
        hmc = hmc.clone()

        if hasattr ( hd     , 'red'   ) : hd .red   () 
        if hasattr ( hmc    , 'blue'  ) : hmc.blue  () 
        rmax   = 1.2 * max ( hd .GetMaximum () , hmc.GetMaximum() )
        hd .SetMinimum ( 0 )
        hmc.SetMinimum ( 0 )
        pad = ROOT.Ostap.Utils.get_pad () 
        if pad :
            scale  = pad.GetUymax() / rmax        
            hd  .Scale ( scale  )
            hmc .Scale ( scale  )        
        hd  .draw  ( 'same' , copy = True )
        hmc .draw  ( 'same' , copy = True )
        
    elif isinstance ( hd  , ROOT.TH1 ) and 1 == hd.dim() :

        hd = hd.clone()
        if hasattr ( hd     , 'red'  ) : hd .red   () 
        rmax   = 1.2 *       hd .GetMaximum ()         
        hd .SetMinimum ( 0 )
        pad = ROOT.Ostap.Utils.get_pad () 
        if pad :        
            scale  = pad.GetUymax() / rmax
            hd  .Scale ( scale  )
        hd  .draw  ( 'same' , copy = True )

    elif isinstance ( hmc , ROOT.TH1 ) and 1 == hmc.dim() :

        hmc = hmc.clone()
        
        if hasattr ( hmc    , 'blue'   ) : hmc.blue   () 
        rmax   = 1.2 *       hmc.GetMaximum ()         
        hmc.SetMinimum ( 0 )        
        pad = ROOT.Ostap.Utils.get_pad () 
        if pad :        
            scale  = pad.GetUymax() / rmax
            hmc .Scale ( scale  )
        hmc .draw  ( 'same' , copy = True )

    return self

ComparisonPlot. draw = _cmp_draw_

# =============================================================================
## The main  function: perform one re-weighting iteration 
#  and reweight "MC"-data set to looks as "data"(reference) dataset
#  @code
#  results = makeWeights (
#   dataset              , ## data source to be  reweighted (DataSet, TTree, abstract source)
#   plots                , ## reweighting plots
#   database             , ## datadabse to store/update reweigting results
#   delta                , ## stopping criteria for "mean"    weight variation
#   minmax               , ## stopping criteria for "min/max" weight variation
#   power                , ## effective power to apply to the weigths
#   debug        = True  , ## store debuig information in database
#   make_plots   = True  , ## produce useful comparison plots
#   force_update = False , ## force DB update even for "good" results 
#   tag          = 'RW'  ) ## tag for better printout 
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
def makeWeights  ( dataset                      ,
                   plots        = []            , 
                   database     = "weights.db"  ,
                   compare      = None          , ## comparison function 
                   delta        = 0.01          , ## delta for `mean'  weight variation
                   minmax       = 0.03          , ## delta for `minmax' weight variation
                   maxchi2      = 1.00          , ## maximum chi2/ndf 
                   power        = None          , ## auto-determination
                   debug        = True          , ## save intermediate information in DB
                   make_plots   = True          , ## make comparison plots (and draw them)
                   wtruncate    = ( 0.5 , 1.5 ) , ## truncate too small or too large weights
                   force_update = False         , ## force DB update even for "good" results 
                   tag          = "Reweighting" ) :
    """ The main  function: perform one re-weighting iteration 
    and reweight `MC'-data set to looks as `data'(reference) dataset
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

    assert 0 < delta  , "makeWeights(%s): Invalid value for `delta'  %s" % ( tag , delta  )
    assert 0 < minmax , "makeWeights(%s): Invalid value for `minmax' %s" % ( tag , minmax ) 

    if wtruncate and isinstance ( wtruncate , sized_types ) and 2 == len ( wtruncate ) \
       and isinstance ( wtruncate [ 0 ] , num_types ) \
       and isinstance ( wtruncate [ 1 ] , num_types ) \
       and wtruncate [0] < 1.0 and 1.0 < wtruncate [1] : pass
    elif not wtruncate : wtruncate = ()  ## No truncation  
    else :        
        logger.warning  ( "Invalid 'wtruncate' parameter '%s', use ()" % ( str ( wtruncate ) ) )
        wtruncate = () 
    
    nplots  = len ( plots )

    ## #itemd that have been updated at this iterations and *not-good*
    for_update   = 0 
    ## list of plots to compare 
    cmp_plots    = []
    ## reweighting summary table
    header       = ( 'Reweighting' , 'wmin/wmax' , 'OK?' , 'wrms[%]' , 'OK?' , chi2ndf , 'OK?' , 'ww' , 'exp' )
    rows         = {}
    save_to_db   = [] 
    ## number of active plots for reweighting
    for wplot in plots  :
        
        what      = wplot.what       ## variable/function to plot/compare 
        how       = wplot.how        ## weight and/or additional cuts 
        address   = wplot.address    ## address in database 
        hdata0    = wplot.data       ## original "DATA" object 
        hmc0      = wplot.mc_histo   ## original "MC"   histogram 
        ww        = wplot.w          ## relative weight
        projector = wplot.projector  ## projector for MC data
        ignore    = wplot.ignore     ## ignore for weight building?
        #
        # normalize the data
        #
        hdata = hdata0
        if isinstance ( hdata , ROOT.TH1 ) and not hdata.is_density() : 
            hdata = hdata.density ()
            logger.atention ( "'Data' histogram converted to 'density' for '%s'" %  what )
            
        # =====================================================================
        ## make a plot on (MC) data with the weight
        # =====================================================================
        hmc0 = projector ( dataset , hmc0 , what , how )
        
        st   = hmc0.stat()
        mnmx = st.minmax()
        if  iszero ( mnmx [ 0 ] ) :
            logger.warning ( "%s: MC statistic goes to zero %s/`%s'" % ( tag , st , address ) )            
        elif mnmx [ 0 ] <= 0    :
            logger.warning ( "%s: MC statistic is negative  %s/`%s'" % ( tag , st , address ) ) 
            
        # =====================================================================

        ## make a try to normalize MC projection:
        hmc = normalize_data ( hmc0 , hdata )

        # =====================================================================
        ## calculate  the reweighting factor : a bit conservative (?)
        #  this is the only important line
        # =====================================================================
        
        #  try to exploit finer binning if/when possible
        hboth = isinstance ( hmc   , ROOT.TH1 ) and isinstance ( hdata , ROOT.TH1 ) and hmc.dim() == hdata.dim () 

        ## 
        if hboth and hmc.finer_bins ( hdata ) : w = ( 1.0 / hmc ) * hdata    ## NB!
        else                                  : w = hdata / hmc              ## NB!

        # =====================================================================
        ## Here we need to get the statistics of weights
        # =====================================================================
        ## scale & get the statistics of weights
        wstat  = w.stat      () 
        wmean  = wstat.mean  () 
        w     /= wmean.value ()
        
        wstat  = w.stat      ()
        wmean  = wstat.mean  ()

        ## RMS: 
        wrms   = wstat.rms   () / wmean
        
        ## MINMAX: 
        wmin , wmax = wstat.minmax()
        ## CHI2:
        meanv  = float ( wmean )
        ndf    = len ( w ) 
        c2ndf  = sum ( w [ i ].chi2 ( meanv ) for i in w ) 
        c2ndf /= ndf 

        # =====================================================================
        ## make decision
        # =====================================================================
        
        good1 = abs   ( wmax - wmin ) <= minmax
        good2 = float ( wrms        ) <= delta
        good3 = c2ndf                 <= maxchi2 
        good  = good1 and good2 and good3 ## small variance?         
        #
        
        # =====================================================================
        ## document the decision 
        # =====================================================================
                
        ## build  the row in the summary table 
        row = [ address  , '%-6.3f/%6.3f' % ( wmin , wmax ) ] 

        if   ignore : row.append ( ''  )
        elif good1  : row.append ( OK  )
        else        : row.append ( NOT )
        
        row.append ( ( wrms * 100 ).toString('%%+7.3f %s %%-7.3f' % plus_minus ) )

        if   ignore : row.append ( ''  )
        elif good2  : row.append ( OK  )
        else        : row.append ( NOT )
        
        row.append (  '%7.3f' % c2ndf  )
        
        if   ignore : row.append ( ''  )
        elif good3  : row.append ( OK  )
        else        : row.append ( NOT )

        ## apply weight truncation:
        if not ignore and wtruncate : 
            wmin , wmax = wtruncate
            truncated   = False 
            for i in w  :
                we = w [ i ]
                if   we.value() < wmin :
                    w [ i ] = VE ( wmin , we.cov2 () )
                    truncated    = True 
                elif we.value() > wmax :
                    w [ i ] = VE ( wmax , we.cov2 () ) 
                    truncated    = True
                    
            if truncated :
                logger.info ( "%s: weights for '%s' are truncated to be %.3f<w<%.3f" % ( tag , what , wmin , wmax ) )
            
        ## make plots at the start of  each iteration? 
        if make_plots : 
            item = ComparisonPlot ( what , hdata , hmc , w ) 
            cmp_plots.append ( item )

        if ignore : row.append ( '' ) 
        else      : row.append ( '%4.3f' % ww if 1 != ww  else '' ) 

        rows [ address ] = tuple ( row ) 
        
        # =====================================================================
        ## prepare the reweighting entry to datatbase  
        # =====================================================================
        if not ignore  :
            # =================================================================
            ## update DB only for "not-good" or "forced" entries
            # =================================================================
            if force_update or not good :
                if not good : for_update += 1
                save_to_db.append ( ( address , ww , hdata0 , hmc0 , hdata , hmc , w ) )

        # =====================================================================
        ## make a comparison (if needed)
        # =====================================================================
        if compare : compare ( hdata0 , hmc0 , address )

    active  = tuple ( [ p [ 0 ] for p in save_to_db ] )  
    nactive = len ( active )  

    if   power is None :
        ## Variables are 100% uncorelated? power can be 1.0 
        p1 = 1.0
        ## Variables are 100% correlated? power should be 1/nactive  
        p2 = 1.0 / ( nactive - 1 ) if 2 <= nactive  else 1.0
        ## 
        ## for realistic case take something intermediate 
        ## 
        alpha   = 0.75 
        ## average between 100% correlated and 100% uncorrelated 
        eff_exp = alpha * p1 + ( 1 - alpha ) * p2 
        
    elif power and callable ( power )                          : eff_exp = power ( nactive ) 
    elif isinstance ( power , num_types ) and 0 < power <= 2.0 : eff_exp = 1.0 * power 
    elif 1 == nactive and 1 < len ( plots )                    : eff_exp = 0.95 ## slight underreweighting
    elif 1 == nactive                                          : eff_exp = 1.01 ## slight overreweighting
    else                                                       : eff_exp = 1.05 / max ( nactive - 1 , 1 ) 

    # =============================================================================
    ## save the results into data base +  abit of black magid with effective exponent 
    # =============================================================================
    while database and save_to_db :

        entry = save_to_db.pop() 
        
        address , ww , hd0 , hm0 , hd , hm , weight = entry  

        cnt = weight.stat()
        mnw , mxw = cnt.minmax()

        if 1 < nactive and 1 != ww :            
            eff_exp *= ww
            logger.info  ("%s: apply `effective exponent' of %.3f for `%s'" % ( tag , eff_exp  , address ) )

        if 1 != eff_exp and 0 < eff_exp :
            
            weight = weight ** eff_exp

            row = list ( rows [ address ] )
            row.append ( '%4.3f' % eff_exp )            
            rows [ address ] = tuple ( row ) 
                        
        with DBASE.open ( database ) as db :

            rw = list ( db.get( address , [] ) )
            rw.append ( weight )
            db [ address ] = tuple ( rw )

            ## update the list of known reweightings
            
            tag       = tag_reweightings
            addresses = db.get( tag , () )
            addresses = set ( addresses ) 
            addresses.add ( address )
            addressed = sorted ( addresses ) 
            db [ tag ] = tuple ( addresses )
            
            ## save more information to database 
            if debug and True :
                addr        = '%s:%s' % ( address , tag_rw_suffix ) 
                rw = list  ( db.get ( addr , [] ) )
                e2 = tuple ( entry [ 2 : ] ) 
                rw.append ( e2 )                  
                db [ addr ] = tuple ( rw ) 
                        
        del hd0, hm0 , hd , hm , weight , entry 

    table = [  header ]
    for   row in   rows  : table.append ( rows [ row ] )
        
    import ostap.logger.table as Table
    table = Table.remove_empty_columns ( table ) 
    table = Table.table ( table , title = tag , prefix = '# ' , alignment = 'lccccccc' )
    logger.info ( '%s: %d active reweightings\n%s' % ( tag , nactive , table ) ) 

    cmp_plots = tuple ( cmp_plots )

    if make_plots and cmp_plots :
        groot = ROOT.ROOT.GetROOT() 
        if groot and not groot.IsBatch() : 
            from ostap.plotting.canvas import use_canvas 
            for item in cmp_plots :
                colz = \
                    isinstance ( item.data   , ROOT.TH1 ) and 1 == item.data  .dim () and \
                    isinstance ( item.mc     , ROOT.TH1 ) and 1 == item.mc    .dim () and \
                    isinstance ( item.weight , ROOT.TH1 ) and 1 == item.weight.dim ()                
                with use_canvas ( '%s/%s' % ( tag , item.what ) , colz = colz ) : item.draw()

        ## save comparison plots to database 
        if debug and True :
            tag = tag_comparisons
            with DBASE.open ( database ) as db :
                comparison_plots  = list ( db.get( tag , () ) )
                comparison_plots.append ( cmp_plots ) 
                db [ tag ] = tuple ( comparison_plots ) 
                
    ## return ( active , cmp_plots ) if make_plots else active
    return ( for_update , cmp_plots ) if make_plots else for_update 

# =============================================================================
## helper factory for proper (de)serialization 
def w2t_factory ( *args ) : return W2Tree ( *args ) 
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
    """ Helper class to add the weight into ROOT.TTree
    >>> w    = Weight ( ... )      ## the weighting object
    >>> tree = ...                 ## The tree
    >>> wf   = W2Tree ( w , tree ) ## create the weighting function
    >>> tree.add_new_branch ( 'weight' , wf )
    """
    store = set()
    
    def __init__ ( self , weight = None , tree = None , clone = None ) :

        assert not clone or isinstance ( clone , W2Tree ) , \
            "W2Tree: Invalid 'clone' type!"  % typename ( clone )        
        assert clone or isinstance ( weight , Weight ) , \
            "W2Tree: Invalid 'weight' type!" % typename ( weight )

        ## initialize the base 
        super(W2Tree,self).__init__ ( tree = tree , clone = clone )

        if clone : self.__weight = clone.weight 
        else     : self.__weight = weight
        
    ## clone function 
    def clone ( self , name = "" ) :
        """ Clone it! """
        cloned = W2Tree ( weight = self.weight , tree = self.the_tree , clone = self )
        ROOT.SetOwnership ( cloned , False )
        if ( 6 , 32 ) <= root_info : self.store.add ( cloned )         
        return cloned 

    ## evaluate the weighter for the given TTree entry 
    def evaluate ( self ) : 
        """ Evaluate the weighter for the given TTree entry """
        t = self.tree ()
        w = self.__weight ( t ) 
        return w
    
    ## REDUCE 
    def __reduce__  ( self ) : return w2t_factory , ( self.weight , ) 

    @property
    def weight ( self ) :
        """`weight` : get the weighter object"""
        return  self.__weight

# =============================================================================
## helper factory for proper (de)serialization 
def w2d_factory ( *args ) : return W2Data ( *args ) 
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
    """ Helper class to add the weight into <code>RooDataSet</code>
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
    
    ## REDUCE 
    def __reduce__  ( self ) : return w2d_factory , ( self.weight , ) 

    @property
    def weight ( self ) :
        """`weight' : get the weighter object"""
        return  self.__weight

# =============================================================================
## Add specific re-weighting information into <code>ROOT.TTree</code>
#  @see ostap.tools.reweight
#  @see ostap.tools.reweight.Weight 
#  @see ostap.tools.reweight.W2Tree 
#  @code
#  w    = Weight ( ... ) ## weighting object ostap.tools.reweight.Weight 
#  tree = ...
#  tree.add_reweighting ( w ) 
#  @endcode 
def tree_add_reweighting ( tree                 ,
                           weighter             ,
                           name      = 'weight' ,
                           progress  = True     ,
                           report    = True     ,
                           parallel  = False    ) :
    """ Add specific re-weighting information into ROOT.TTree
    
    >>> w    = Weight ( ... ) ## weighting object ostap.tools.reweight.Weight 
    >>> data = ...
    >>> data.add_reweighting ( w )
    - see ostap.tools.reweight
    - see ostap.tools.reweight.Weight 
    - see ostap.tools.reweight.W2Tree 
    """
    assert isinstance ( weighter , Weight     ), "Invalid type of `weighter'!"  
    assert isinstance ( tree     , ROOT.TTree ), "Invalid type of `tree`!"

    if not weighter :
        logger.attention ( 'Weighter object is empty! Add 1.0 to the data' )
        if parallel and isinstance ( tree , ROOT.TChain ) and 1 < tree.nFiles : 
            import ostap.parallel.parallel_add_branch
            return tree.padd_new_branch ( 1.0 , name = name , progress = progress , report = report )            
        return tree.add_new_branch ( 1.0 , name = name , progress = progress , report = report )
    
    ## create the weighting function 
    wfun = W2Tree ( weighter )
    
    if parallel and isinstance ( tree , ROOT.TChain ) and 1 < tree.nFiles and 1 < numcpu () :
        from ostap.parallel.parallel import Checker
        checker = Checker()
        if checker.pickles ( weighter , wfun ) :         
            import ostap.parallel.parallel_add_branch
            return tree.padd_new_branch ( wfun , name = name  , progress = progress , report = report ) 
    
    ## regular sequential processing 
    return tree.add_new_branch ( wfun , name = name , progress = progress , report = report ) 

ROOT.TTree.add_reweighting = tree_add_reweighting    

# =============================================================================
## Add specific re-weighting information into dataset
#  @see ostap.tools.reweight
#  @see ostap.tools.reweight.Weight 
#  @see ostap.tools.reweight.W2Data 
#  @code
#  w    = Weight ( ... ) ## weighting object ostap.tools.reweight.Weight 
#  data = ...
#  data.add_reweighting ( w ) 
#  @endcode 
def data_add_reweighting ( data                ,
                           weighter            ,
                           name     = 'weight' ,
                           progress = False    ,
                           report   = True     ,
                           parallel = False    ) :
    """ Add specific re-weighting information into dataset
    
    >>> w    = Weight ( ... ) ## weighting object ostap.tools.reweight.Weight 
    >>> data = ...
    >>> data.add_reweighting ( w )
    - see ostap.tools.reweight
    - see ostap.tools.reweight.Weight 
    - see ostap.tools.reweight.W2Data 
    """

    assert isinstance ( weighter , Weight          ), "Invalid type of `weighter'!"
    assert isinstance ( data     , ROOT.RooDataSet ), "Invalid type of `data`!"

    if not weighter :
        logger.attention ( 'Weighter object is empty! Add 1.0 to the data' )
        return data.add_new_var ( name , "1.0" , progress = progress )
        
    if parallel and 1 < numcpu () :
        from ostap.parallel.parallel import Checker
        checker = Checker()
        if checker.pickles ( weighter ) :
            return parallel_data_add_reweighting ( data                ,
                                                   weighter            ,
                                                   name     = name     ,
                                                   progress = progress ,
                                                   report   = report   )   

    ## create the weighting function 
    wfun = W2Data ( weighter )

    ## regular sequential processing 
    return data.add_new_var ( name , wfun , progress = progress ) 

# ===========================================================================
## Add reweighting to very large RooDataSet by splitting it into chunks and
#  processing the chunks in parallel 
def parallel_data_add_reweighting ( dataset          ,
                                    weighter         , 
                                    name             ,                                    
                                    progress = True  ,                                    
                                    report   = True  , **kwargs ) :
    """ Add reweighting to very large RooDataSet by splititng it into chunks and
     processing the chunks in parallel 
    """
    assert isinstance( dataset , ROOT.RooDataSet ) , \
        "Invald type for `dataset` : %s" % typename ( dataset ) 
    
    if numcpu() < 2 : return data_add_reweighting ( dataset  ,
                                                    weighter ,
                                                    name     = name     ,
                                                    progress = progress ,
                                                    report   = report   , 
                                                    parallel = False    ) 

    ## list of branches 
    branches = set ( dataset.branches() ) if report else set() 

    ## create the weighting function 
    wfun = W2Data ( weighter  )

    from ostap.parallel.parallel import Checker, WorkManager
    
    checker = Checker()
    if not checker.pickles ( weighter , wfun ) :
        return data_add_reweighting ( dataset , 
                                     weighter , 
                                     name     = name     , 
                                     progress = progress ,
                                     report   = report   , 
                                     parallel = False    ) 


    nchunks  = 3 * numcpu ()
        
    ## split dataset into `nchunks` chunks:
    from ostap.utils.utils       import split_n_range
    chunks = ( dataset [ i : j ] for i , j in split_n_range ( 0 , len ( dataset ) , nchunks ) )

    ## Task to add new variable
    from ostap.parallel.parallel_add_branch import AddNewVar 
    task     = AddNewVar   ( name            ,  wfun    )

    silent   = kwargs.pop ( 'silent'  , not progress )    
    wmgr     = WorkManager ( silent = silent , progress = progress , **kwargs )
    
    ## start processing 
    wmgr.process ( task , chunks )

    ## new dataset 
    newds = task.results()

    ## check consistency: 
    assert len ( dataset ) == len ( newds ) , "Mismatch in dataset sizes!"
    
    if report :
        new_branches = sorted ( set ( newds.branches() ) - branches )
        if new_branches :
            n = len ( new_branches )
            if 1 >= n : title = "Added %s variable to RooDataSet(%s)"  % ( n , dataset.GetName () ) 
            else      : title = "Added %s variables to RooDataSet(%s)" % ( n , dataset.GetName () )
            table = newds.table ( new_branches , title = title , prefix = '# ' )
            logger.info ( '%s:\n%s' % ( title , table ) )
        else :                                
            logger.error ( "No variables are added to RooDataSet(%s)" %  dataset.GetName () )


    dataset.add_new_var ( name , '1' )             
    dataset.clear ()

    dataset.append ( newds )

    newds.clear()
    del newds
    del task

    
    return dataset

##    del dataset
##    del task    
##    return newds 

ROOT.RooDataSet.add_reweighting = data_add_reweighting

_new_methods_ += [
    ROOT.RooDataSet .add_reweighting , 
    ROOT.TTree      .add_reweighting , 
    ]

# =========================================================================
## backup the database content to vanilla ROOT file *IF/WHEN* possible
#  @attention Only simplistic histogram-bases structures can be converted.
#  @code
#  dbname     = ...
#  root_file = backup_to_ROOT ( dbname , 'root_file.root' )
#  @endcode
#  If root file name is not specified, a temporary file will be created 
def backup_to_ROOT ( dbase , root_file = '' , * , reweightings = () ) :
    """ Backup the database content to vanilla ROOT file *IF/WHEN* possible 

    >>> dbname     = ...
    >>> backup_to_ROOT ( dbname , 'root_file.root' )

    If root file name is not specified, a temporary file will be created

    """
    if not root_file :
        from   ostap.utils.cleanup    import CleanUp        
        root_file = CleanUp.tempfile ( suffix = '.root' , prefix ='ostap-reweight-' )
        
    ## (re)create ROOT TFile 
    with ROOT.TFile ( root_file , 'c' ) : pass
    
    with DBASE.open ( dbase , 'r' ) as db : ## READONLY
        
        rws  = db.get ( tag_reweightings , () )
        
        rws  = set ( rws ) | set ( reweightings )
        rws  = sorted ( rws )
        
        good = []
        for r in rws :
            
            if not r in db :
                logger.warning ( 'Record %s not in database, skip' % r )  
                continue
            
            rw = db.get ( r , None )
            if not rw      :
                logger.warning ( 'Record %s not in database, skip' % r )  
                continue

            if isinstance ( rw , ROOT.TObject ) : rw = [ rw ]
            
            if not all ( isinstance ( o , ROOT.TObject ) and callable ( o )  for o in rw ) : 
                logger.warning ( 'Record %s not callable TObject or sequence of callable TObjects, skip!' % r )  
                continue
            
            with ROOT.TFile ( root_file , 'u' ) as rf : 
                for i , o in enumerate ( rw ) : rf [ '%s/%d' %  ( r , i ) ] = o 
                
            good.append ( r )

            
        tlst = ROOT.TList()
        for r in good :
            tstr = ROOT.TObjString ( r )
            tlst.Add  ( tstr ) 

    title = 'Reweighting DB -> ROOT'
    with ROOT.TFile ( root_file , 'u' ) as rf :
        rf [ tag_reweightings ] = tlst
        table = rf.ls_table ( title = title , prefix = '# ' )        
        logger.info ('%s:\n%s' % ( title , table ) )        
        rf.ls_tree   ( prefix = '# ' )
    
    return root_file 

# ==================================================================================
## Convert the ROOT-file into reweighting DBASE
#  @code
#  root_file_name = ...
#  respore_form_ROOT ( root_file_name , 'dbase.db' ) 
#  @endcode 
def restore_from_ROOT ( root_file , dbname = '' , * , reweightings = () ) :
    """ Convert the ROOT-file into proper dbase for reweighting 

    >>> root_file_name = ...
    >>> restore_from_ROOT ( root_file_name , 'dbase.db' ) 

    """
    if not dbname :
        from   ostap.utils.cleanup    import CleanUp        
        dbname = CleanUp.tempfile ( suffix = '.root' , prefix ='ostap-reweight-' )

    ## (Re)create the data base 
    with DBASE.open ( dbname , 'c' ) : pass

    ## read input ROOT file 
    with ROOT.TFile ( root_file , 'r' ) as rf :

        ## Read the entry 'tag_reweightings:  It *must* be TList with TObjStrings
        rws = rf.get ( tag_reweightings , () )
        if not rws                             or \
           not isinstance ( rws , ROOT.TList ) or \
           not all ( isinstance ( r , ROOT.TObjString ) for r in rws ) : 
            logger.warning ( 'No valid reweightings are found!' ) 
            return

        ## add 'reweightings' from arguments 
        rws  = set ( [ str ( r ) for r in rws ] ) | set ( reweightings )
        rws  = tuple ( sorted ( rws ) )

        ## now `rws` is alist of potential reweighting structures, placed in corresponsig directories 
        
        good = []
        
        ## loop over the specified directories/reweightings  
        for r in rws :

            ## root directory            
            rdir = rf.get ( r , None ) 
            if not rdir or not isinstance ( rdir , ROOT.TDirectory ) :
                logger.warning ( 'Not valid TDirectory %s, skip' % r )
                continue
            
            objects = []
            ## loop over the content of the directory
            for key, value in rdir.iteritems ( recursive = False , no_dir = True ) :                
                if value and callable ( value ) :
                    objects.append ( value )
                else : 
                    logger.warning ( 'Invalid callable %s/%S, skip' % ( r , key ) ) 
                    continue 

            good.append ( r )            
            with DBASE.open ( dbname , 'u' ) as db : db [ r ] = objects

    ## add the list of reweightings to 
    with DBASE.open ( dbname , 'u' ) as db :
        rr = db.get ( tag_reweightings , () )
        rr = set ( rr ) | set ( good )        
        db [ tag_reweightings ] = tuple ( sorted ( rr ) )  

    ## finally show the content of converted DBASE
    logger.info ( 'Reweighting DB from ROOT: %s -> %s' % ( root_file , dbname ) )     
    with DBASE.open ( dbname , 'r' ) as db :
        db.ls()
        rw = db.get ( tag_reweightings , () )
        if not rw : logger.error ( "DB has no '%s' entry" % tag_reweightings )
        for r in rw :
            if not r in db : logger.error ( "DB has no '%s' reweighting" % r )
                
    return dbname 

# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================

