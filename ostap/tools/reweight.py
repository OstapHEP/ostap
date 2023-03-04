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
from   ostap.core.pyrouts     import VE, SE, Ostap 
from   ostap.math.base        import iszero
from   ostap.core.core        import split_string 
from   ostap.core.ostap_types import string_types, list_types, num_types, sized_types, sequence_types    
from   ostap.math.operations  import Mul as MULT  ## needed for proper abstract multiplication
import ostap.io.zipshelve     as     DBASE ## needed to store the weights&histos
from   ostap.trees.funcs      import FuncTree, FuncData ## add weigth to TTree/RooDataSet
from   ostap.utils.utils      import CallThem, is_formula  
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
        """`attributes': the actual attributes
        """
        return self.__attributes
    # print attributes 
    def __str__ ( self ) :
        return ','.join ( self.__attributes )
    __repr__ = __str__ 
# =============================================================================
## @class Weight
#  helper class for semiautomatic reweighting of data 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-05-10
class Weight(object) :
    """Helper class for semi-automatic reweighting of data

    It reads various `weigth components' from DBASE and calculates
    the `global' event weight via simple accumulation of weights
    
    Simplest case: calculate weights for one variable.
    
    factors = [ ( fun , name ) ] 
    where `fun' is a function to get variable from TTree/RooDataSet 
    and   `name' is address in data base with reweighting information 
    e.g.
    
    >>> factors = [ ( lambda s : s.pt  , 'pt-data' ) ]
    
    It assumes, that DBASE has an entry with name 'pt-data', such
    that
    
    >>> funcs = database['pt-data']
    
    Here `func' will be function/callable or list of functions/callables 
    with weights:
    
    >>> pt = fun ( entry ) 
    >>> w  = 1.0    
    >>> w *= func ( pt )              ## for the single function/callable
    >>> for f in funcs : w*= f ( pt ) ## for list of functions/callables  
    
    quantity `w' will be an event  weight
    
    if list of `factors' constains more than one entry,
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

            logger.debug ( 'Weight: Reweigting database: \n%s' % db.table ( prefix = '# ' ) ) 
                
            ## loop over the weighting factors and build the function
            for wvar in factors :

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
                    
                ## 
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
    
    ## calculate the weight for the given "event"
    def __call__ ( self , s ) :
        """Calculate the weigth for the given `event' (==record in TTree/TChain or RooDataSet):
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
        import ostap.logger.table as T
        if title is None : title = "Weight('%s')" % self.__dbase 
        return T.table ( self.__table , title = title , prefix = prefix ,
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

                ## grpah with the spread(min/max) weights 
                gr1 = ROOT.TGraphAsymmErrors ( len ( functions )  )
                ## graph with the rms weights 
                gr2 = ROOT.TGraphErrors      ( len ( functions )  )

                has_points = False 
                for n , w  in enumerate ( functions , start = 1 ) :
                    
                    if not hasattr ( w  ,'stat' ) : continue
                    
                    cnt = w.stat()
                    wmean       = cnt.mean  ().value()
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
                    
                    ## combine two graphs togather 
                    graph = ROOT.TMultiGraph ()
                    graph.Add      ( gr1 , '3'  )
                    graph.Add      ( gr2 , 'pl' )
                    graph.SetTitle ( '%s;iteration;weights' % address )
                    
                    grphs [ address ] = graph

        return grphs     

    def __str__  ( self ) : return self.table()
    def __repr__ ( self ) : return self.table()

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
        - `accessor'  : an accessor function that extracts the variable(s) from  TTree/TChain/RooDataSet
        - `address'   : the  address in DBASE, where reweigftjnig callable(s) is/are stored
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
            """Keep information about single reweighting
            - `accessor' : an accessor function that extracts the variable(s) from  TTree/TChain/RooDataSet
            - `address'  : the  address in DBASE, where reweigftjnig callable(s) is/are stored
            - `merge'    : merge list of callables from DB into the single callable ?
            - `skip'     : use only certain elements from the list of callables from DBASE
            
            Schematic data flow to get the weigth for the given event 
            - tree/chain/dataset -> accessor -> database(address) -> weight
            """

            if isinstance ( accessor , string_types ) :
                accessor = split_string ( accessor , strip = True )
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
                   "Invalid type of `accessor' %s/%s" % ( accessor , type( accessor ) )
            
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
            for the  given varibale(s).
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
    """Default function to make MC projection
    """
    
    dataset.project ( histo , what , how )
    return histo 


# =============================================================================
## @class WeightingPlot
#  helper class to manage/keep `weighting-plot' 
class WeightingPlot(object) :    
    """Helper class to manage/keep `weighting-plot'
    
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
        """Helper class to manage/keep `weighting-plot'
        
        - `what' : the variable/expression to be used for `weighting-plot'
        Used  as the second argument of `dataset.project' method to produce
        `weighted-plot':
        
        >>> dataset.project ( mchisto , WHAT , how , ... ) 
        
        - `how' : the `weight' expression to be used for `weighed-plots'
        Used  as the third argument of `dataset.project' method to produce
        `weighted-plot':
        
        >>> dataset.project ( mchisto , what , HOW , ... )
        
        Typically it refers to `weight' variable in datasete
        
        - `address' : the addres in `weighting-database'
        where to store the obtained weights
        
        - `data' : the `data' object, or  the `target' for the reweighting procedure
        Typically it is a histotgram. But it could be any kind of callable 
        
        - `mc_histo' : template/shape for the mc-histogram, to be used for reweighting.
        It is used as the  first argument of `dataset.project' method
        
        >>> dataset.project ( MCHISTO , what , how , ... )         
        """

        assert data and callable ( data ) , "ReweightingPlot: `data' object must be callable!"
        
        self.__what      = str(what)     if isinstance ( what , str ) else what 
        self.__how       = str(how )     if isinstance ( how  , str ) else how 
        self.__address   = str(address)

        ## ATTENTION! make a ty to convert data to density!!
        if isinstance ( data , ROOT.TH1 ) :
            if not data.is_density() :
                logger.info ( "WeightingPlot: Convert 'data' histogram into 'density' for '%s'" % self.what )
            self.__data = data.density()
        else :
            self.__data = data

        self.__mc       = mc_histo      if mc_histo else data.clone()
        
        assert isinstance ( self.__mc , ROOT.TH1 ), \
               "WPlot: invalid type of `mchisto' %s/%s"  % ( self.__mc , type ( self.__mc ) )
        self.__w         = w

        if not projector : projector = mc_data_projector
            
        self.__projector = projector
        
        assert self.projector and callable ( self.projector ) ,"`Projector' must be callable!"
        self.__ignore    = True if ignore else False 
        
    @property
    def  what ( self ) :
        """'what' : the variable/expression to be used for `weighting-plot'
        Used  as the second argument of `dataset.project' method to produce
        `weighted-plot':
        >>> dataset.project ( mchisto , WHAT , how , ... ) 
        """
        return self.__what
    @property
    def  how  ( self ) :
        """'how' : the `weight' expression to be used for `weighed-plots'
        Used  as the third argument of `dataset.project' method to produce
        `weighted-plot':
        >>> dataset.project ( mchisto , what , HOW , ... )
        Typically it refers to `weight' variable in datasete
        """
        return self.__how 
    @property
    def  address ( self ) :
        """'address' : the address in `weighting-database'
        where to store the obtained weights
        """
        return self.__address
    @property
    def  data ( self ) :
        """'data' : the `data' object, or  the `target' for the reweighting procedure
        Typically it is a histogram. But it could be any kind of callable/function 
        """
        return self.__data 
    @property
    def  mc_histo ( self ) :
        """'mc_histo' : template/shape for the mc-histogram, to be used for reweighting.
        It is used as the  first argument of `dataset.project' method 
        >>> dataset.project ( MCHISTO , what , how , ... )         
        """
        return self.__mc
    @property
    def projector ( self ) :
        """'projector' :  callable function to build MC distribution:
        hmc = projector ( dataset , hmc ) 
        """
        return self.__projector
    
    @property
    def w  ( self )   :
        """'w'  - relative weight (relative importance is this variable)"""
        return self.__w 
    @property
    def ignore ( self ) :
        """'ignore' : do not use variable in reweights, but use for comparsion"""
        return self.__ignore
    
# =============================================================================
from collections import namedtuple 
ComparisonPlot = namedtuple ( 'ComparisonPlot' , ( 'what' , 'data' , 'mc' , 'weight' ) )

_store = set() 
# =============================================================================
## draw comparison plots
def _cmp_draw_ ( self ) :
    """Draw comparison plots
    """
    hd   = self.data
    hmc  = self.mc
    hw   = self.weight 

    if   isinstance ( hw , ROOT.TH3 ) and 3 == hw.dim () :
        hw.draw  ( copy = True )

        ## minx, miny, minz = hw.minimum_bin()
        ## maxx, maxy, maxz = hw.maximum_bin()
        
        ## ax = hw.GetXaxis()
        ## ay = hw.GetYaxis()
        ## az = hw.GetZaxis()

        ## minv = hw [ minx, miny, minz ]
        ## maxv = hw [ maxx, maxy, maxz ]
        
        ## logger.info ( "ComparisonPlot('%s'): min-weight %+.3f at (%.5g,%.5g,%.5g)" % ( self.what ,
        ##                                                                                minv      ,
        ##                                                                                ax.GetBinCenter ( minx ) ,
        ##                                                                                ay.GetBinCenter ( miny ) ,
        ##                                                                                az.GetBinCenter ( minz ) ) )
        ## logger.info ( "ComparisonPlot('%s'): max-weight %+.3f at (%.5g,%.5g,%.5g)" % ( self.what ,
        ##                                                                                maxv      ,
        ##                                                                                ax.GetBinCenter ( maxx ) ,
        ##                                                                                ay.GetBinCenter ( maxy ) ,
        ##                                                                                az.GetBinCenter ( maxz ) ) )
                
    elif isinstance ( hw , ROOT.TH2 ) and 2 == hw.dim () :

        zmin,zmax = hw.zminmax()

        if   0.95 <= zmin and zmax <= 1.05 :
            hw.SetMinimum  ( 0.95 )
            hw.SetMaximum  ( 1.05 )
            hw.SetContour  ( 50   ) 
        elif 0.90 <= zmin and zmax <= 1.10 :
            hw.SetMinimum  ( 0.90 )
            hw.SetMaximum  ( 1.10 )
            hw.SetContour  ( 40   ) 
        elif 0.80 <= zmin and zmax <= 1.20 :
            hw.SetMinimum  ( 0.80 )
            hw.SetMaximum  ( 1.20 )
            hw.SetContour  ( 40   ) 
        elif 0.75 <= zmin and zmax <= 1.25 :
            hw.SetMinimum  ( 0.75 )
            hw.SetMaximum  ( 1.25 )
            hw.SetContour  ( 50   ) 
        elif 0.50 <= zmin and zmax <= 1.50 :
            hw.SetMinimum  ( 0.50 )
            hw.SetMaximum  ( 1.50 )
            hw.SetContour  ( 50   )
        else :
            hw.SetMinimum  ( 0.50 )
            hw.SetMaximum  ( 2.00 )
            hw.SetContour  ( 50   ) 

        nx = hw.binsx()
        ny = hw.binsy()
        if   10 <  nx and 10 < ny : hw.draw  ( 'cont4z' ,                copy = True )
        elif 10 <  nx or  10 < ny : hw.draw  ( 'cont4z' ,                copy = True )
        else                      : hw.texte ( 'colz'   , fmt = '5.3f' , copy = True )

        minx, miny = hw.minimum_bin()
        maxx, maxy = hw.maximum_bin()

        ax = hw.GetXaxis()
        ay = hw.GetYaxis()

        pmin = ROOT.TMarker ( ax.GetBinCenter ( minx ) , ay.GetBinCenter ( miny ) , 43 )
        pmax = ROOT.TMarker ( ax.GetBinCenter ( maxx ) , ay.GetBinCenter ( maxy ) , 20 )
        pmin.SetMarkerColor ( 2 )
        pmax.SetMarkerColor ( 2 )
        pmin.SetMarkerSize  ( 2 )
        pmax.SetMarkerSize  ( 2 )
        pmin.DrawClone() 
        pmax.DrawClone()
        
        ## minv = hw [ minx, miny ]
        ## maxv = hw [ maxx, maxy ]
        
        ## logger.info ( "ComparisonPlot('%s'): min-weight %+.3f at (%.5g,%.5g)" % ( self.what ,
        ##                                                                           minv      ,
        ##                                                                           ax.GetBinCenter ( minx ) ,
        ##                                                                           ay.GetBinCenter ( miny ) ) ) 
        ## logger.info ( "ComparisonPlot('%s'): max-weight %+.3f at (%.5g,%.5g)" % ( self.what ,
        ##                                                                           maxv      ,
        ##                                                                           ax.GetBinCenter ( maxx ) ,
        ##                                                                           ay.GetBinCenter ( maxy ) ) ) 
        
    elif isinstance ( hw , ROOT.TH1 ) and 1 == hw.dim () :
        if hasattr ( hw     , 'red'   ) : hw .red   ()         
        hw.SetMinimum  ( 0 )
        hw.SetMaximum  ( max ( 1.3 , 1.2 * hw.ymax () ) ) 
        hw.draw ( copy = True )
        hw.level ( 1.0 , linestyle = 1 , linecolor = ROOT.kRed )
            
    if isinstance ( hd , ROOT.TH1 ) and isinstance ( hmc , ROOT.TH1 ) and 1 == hd.dim() :

        hd  = hd .clone()
        hmc = hmc.clone()
        
        if hasattr ( hd     , 'green'  ) : hd .green   () 
        if hasattr ( hmc    , 'blue'   ) : hmc.blue   () 
        rmax   = 1.2 * max ( hd .GetMaximum () , hmc.GetMaximum() )
        hd .SetMinimum ( 0 )
        hmc.SetMinimum ( 0 )        
        scale  = ROOT.gPad.GetUymax() / rmax        
        hd  .Scale ( scale  )
        hmc .Scale ( scale  )        
        hd  .draw  ( 'same' , copy = True )
        hmc .draw  ( 'same' , copy = True )
        
    elif isinstance ( hd  , ROOT.TH1 ) and 1 == hd.dim() :

        hd = hd.clone()
        if hasattr ( hd     , 'green'  ) : hd .green   () 
        rmax   = 1.2 *       hd .GetMaximum ()         
        hd .SetMinimum ( 0 )
        scale  = ROOT.gPad.GetUymax() / rmax
        hd  .Scale ( scale  )
        hd  .draw  ( 'same' , copy = True )
        
    elif isinstance ( hmc , ROOT.TH1 ) and 1 == hmc.dim() :

        hmc = hmc.clone()
        
        if hasattr ( hmc    , 'blue'   ) : hmc.blue   () 
        rmax   = 1.2 *       hmc.GetMaximum ()         
        hmc.SetMinimum ( 0 )        
        scale  = ROOT.gPad.GetUymax() / rmax
        hmc .Scale ( scale  )
        hmc .draw  ( 'same' , copy = True )

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
def makeWeights  ( dataset                    ,
                   plots        = []            , 
                   database     = "weights.db"  ,
                   compare      = None          , ## comparison function 
                   delta        = 0.01          , ## delta for `mean'  weight variation
                   minmax       = 0.03          , ## delta for `minmax' weight variation
                   power        = None          , ## auto-determination
                   debug        = True          , ## save intermediate information in DB
                   make_plots   = True          , ## make comparison plots (and draw them)
                   wtruncate    = ( 0.5 , 1.5 ) , ## truncate too small/large weights
                   force_update = False         , ## force DB update  ven for "good" results 
                   tag          = "Reweighting" ) :
    """The main  function: perform one re-weighting iteration 
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
    else :
        logger.warning  ( "Invalid 'wtruncate' parameter '%s', use (0.5,2.0)" % ( str ( wtruncate ) ) )
        wtruncate = 0.5 , 2.0 
    
    from ostap.logger.colorized   import allright , attention , infostr 
    from ostap.utils.basic        import isatty

    nplots  = len ( plots )

    for_update   = 0 
    ## list of plots to compare 
    cmp_plots    = []
    ## reweighting summary table
    header       = ( 'Reweighting' , 'wmin/wmax' , 'OK?' , 'wrms[%]' , 'OK?' , 'chi2/ndf' , 'ww' , 'exp' )
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
            logger.warning ( "'Data' histogram converted to 'density' for '%s'" %  what )
                                                                        
        # =====================================================================
        ## make a plot on (MC) data with the weight
        # =====================================================================
        hmc0 = projector ( dataset , hmc0 , what , how )
        
        st   = hmc0.stat()
        mnmx = st.minmax()
        if  iszero ( mnmx [0] ) :
            logger.warning ( "%s: statistic goes to zero %s/`%s'" % ( tag , st , address ) ) 
        elif mnmx [0] <= 0      :
            logger.warning ( "%s: statistic is negative  %s/`%s'" % ( tag , st , address ) ) 
            
        # =====================================================================
        ## normalize MC
        # =====================================================================
        if isinstance ( hmc0 , ROOT.TH1 ) or hasattr ( hmc , 'density' ) :
            hmc = hmc0.density()
        else : hmc = hmc0 

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
                 ( hmc.binsy() >= hdata.binsy() ) and \
                 ( hmc.binsx() >  hdata.binsx()   or  \
                   hmc.binsy() >  hdata.binsy() ) :            
            w = ( 1.0 / hmc ) * hdata                                 ## NB! 
        elif hboth and 3 == hmc.dim () and 3 == hdata.dim () and \
                 ( hmc.binsx() >= hdata.binsx() ) and \
                 ( hmc.binsy() >= hdata.binsy() ) and \
                 ( hmc.binsz() >= hdata.binsz() ) and \
                 ( hmc.binsx() >  hdata.binsx()   or  \
                   hmc.binsy() >  hdata.binsy()   or  \
                   hmc.binsz() >  hdata.binsz() ) :
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
        row = [ address  , '%-5.3f/%5.3f' % ( cnt.minmax()[0]    , cnt.minmax()[1] ) ] 

        if   ignore : row.append ( ''                )
        elif good2  : row.append ( allright  ( '+' ) )
        else        : row.append ( attention ( '-' ) )
        
        row.append ( (wvar * 100).toString('%6.2f+-%-6.2f') )

        if   ignore : row.append ( ''                )
        elif good1  : row.append ( allright  ( '+' ) )
        else        : row.append ( attention ( '-' ) )
        
        row.append (  '%6.2f' % c2ndf  ) 

        ## apply weight truncation:
        if not ignore : 
            wmin , wmax = wtruncate
            truncated = False 
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
        
        #
        ## make decision based on the variance of weights 
        #
        mnw , mxw = cnt.minmax()
        if not ignore  :
            ## update DB for "not-good" or "forced" entries 
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
        power   = lambda n : 0.5 * ( 1.0 / max ( n , 1 ) + 1 )
        ## average between 100% correlated and 100% uncorrelated 
        eff_exp = 0.5 * ( 1.0 / max ( nactive , 1 ) + 1.0 )
    elif power and callable ( power ) : 
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

        if 1 < nactive and 1 != ww :
            eff_exp *= ww
            logger.info  ("%s: apply `effective exponent' of %.3f for `%s'" % ( tag , eff_exp  , address ) )

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
    logger.info ( '%s: %d active reweightings\n%s' % ( tag , nactive , Table.table ( table , title = tag , prefix = '# ' , alignment = 'lccccccc' ) ) )

    cmp_plots = tuple ( cmp_plots )

    if make_plots and cmp_plots and not ROOT.gROOT.IsBatch() :
        from ostap.plotting.canvas import use_canvas 
        for item in cmp_plots :
            with use_canvas ( '%s/%s' % ( tag , item.what ) ) : item.draw()
                    
    ## return ( active , cmp_plots ) if make_plots else active
    return ( for_update , cmp_plots ) if make_plots else for_update 

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
        """`weight' : get the weighter object"""
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
        """`weight' : get the weighter object"""
        return  self.__weight
    
# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================

