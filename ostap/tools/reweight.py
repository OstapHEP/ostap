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
    'Weight'        , ## the weighting object  
    'makeWeights'   , ## helper function for   weighting iterations 
    'WeightingPlot' , ## helper object to define the weighting rule & target
    'W2Tree'        , ## helper to add the calculated weight to ROOT.TTree
    'W2Data'        , ## helper to add the clacualted weight to ROOT.RooAbsData
    ) 
# =============================================================================
import ROOT, operator
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger    import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.tools.reweight' )
else                       : logger = getLogger ( __name__               )
# =============================================================================
logger.info ( 'Set of utitilities for re-weigthing')
# =============================================================================
from   ostap.core.pyrouts     import VE, SE
from   ostap.math.base        import iszero
from   ostap.core.ostap_types import string_types, list_types 
from   ostap.math.operations  import Mul as MULT  ## needed for proper abstract multiplication
import ostap.io.zipshelve     as     DBASE ## needed to store the weights&histos
from   ostap.trees.funcs      import FuncTree, FuncData ## add weigth to TTree/RooDataSet
import ostap.trees.trees
import ostap.fitting.dataset
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

        ## open database 
        with DBASE.open ( dbase , 'r' ) as db : ## READONLY
            
            for k in db :
                e = db[k]
                if hasattr ( e , '__len__' ) :  
                    logger.debug( "DBASE ``%.15s'' key ``%.15s'' #%d" % ( dbase ,  k, len( e ) ) ) 
                
            ## loop over the weighting factors and build the function
            for wvar in factors :

                funval  = wvar.accessor  ## accessor to the variable 
                funname = wvar.address   ## address  in database 
                merge   = wvar.merge     ## merge sequence of callables?
                skip    = wvar.skip      ## skip   some of them?
                
                if isinstance ( funval , str ) :
                    ## funval = operator.attrgetter( funval ) 
                    funval = AttrGetter( funval ) 
                    
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
                    functions = functions[:-1*skip] 
                elif 0 == skip :
                    pass
                else :
                    logger.error("Invalid ``skip'' parameter %s/%d for ``%s''" % ( skip , flen , funname ) )
                
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
                
        self.__vars = tuple ( self.__vars ) 

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

            ##  get the weight arguments for given event 
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

    # =========================================================================
    ## @class WeightingVar
    #  Helper class to keep information about singe reweighting
    #  - accessor : an accessor function that extracts the variable(s)
    #               from  TTree/TChain/RooDataSet
    #  - address  : the  address in DBASE, where reweigftjnig callable(s) is/are stored
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
                       address          ,   ## the address of   reweighintg object in DBASE 
                       merge     = True ,   ## merge sequence of reweigthing objects ?
                       skip      = None ) : ## skip some reweigting objects ? 
            """Keep information about singe reweighting
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
            self.__address  = str( address )
            self.__merge    = True if merge else False
            self.__skip     = skip if skip  else 0
            
        @property
        def accessor ( self ) :
            """``accessor'' - the accessro function to get the variable from TTree/TChain/RooDataSet
            e.g.  to get the variable ``x'' : 
            >>> xvar   = lambda s : s.x
            get to     variables fo  2D-reweighting:
            >>> xyvars = lambda s : s.x,s.y         
            """
            return self.__accessor[0]
        @property
        def address  ( self ) :
            """``address''  - the address in DB with the reweighting information
            for the  given varibale(s).
            A callable (or list of  callables) is expected in database
            These callables accept ther result of ``accessor'' as their argument
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
    Typically it is a histotgram. But it could be any kind of callable 
    
    - ``mchisto'' : template/shape for the mc-historgam, to be used for reweighting.
    It is used as the  first argument of ``dataset.project'' method
    
    >>> dataset.project ( MCHISTO , what , how , ... )             
    """
    def __init__ ( self            ,
                   what            ,  
                   how             ,
                   address         ,
                   data            ,
                   mc_histo = None ,
                   w        = 1.0  ) :
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
        
        self.__what     = str(what)     if isinstance ( what , str ) else what 
        self.__how      = str(how )     if isinstance ( how  , str ) else how 
        self.__address  = str(address) 
        self.__data     = data
        self.__mc       = mc_histo      if mc_histo else data.clone()
        assert isinstance ( self.__mc , ROOT.TH1 ), \
               "WPlot: invalid type of ``mchisto'' %s/%s"  % ( self.__mc , type ( self.__mc ) )
        self.__w        = w
        
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
        """``address'' : the addres in ``weighting-database''
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
        """``mc_histo'' : template/shape for the mc-historgam, to be used for reweighting.
        It is used as the  first argument of ``dataset.project'' method 
        >>> dataset.project ( MCHISTO , what , how , ... )         
        """
        return self.__mc
    @property
    def w  ( self )   :
        """``w''  - relative weigtht (relative importance is this variable)"""
        return self.__w 

# =============================================================================
## The main  function: perform one re-weighting iteration 
#  and reweight "MC"-data set to looks as "data"(reference) dataset
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-05-10
def makeWeights  ( dataset                  ,
                   plots    = []            , 
                   database = "weights.db"  ,
                   compare  = None          , ## comparison function 
                   delta    = 0.001         , ## delta for ``mean''  weight variation
                   minmax   = 0.05          , ## delta for ``minmax'' weight variation
                   power    = 0             , ## auto-determination
                   debug    = True          , ## save intermediate information in DB
                   tag      = "Reweighting" ) : 

    assert 0 < delta  , "makeWeights(%s): Invalid value for ``delta''  %s" % ( tag , delta  )
    assert 0 < minmax , "makeWeights(%s): Invalid value for ``minmax'' %s" % ( tag , minmax ) 

    from ostap.logger.colorized   import allright , attention , infostr 
    from ostap.utils.basic        import isatty

    power   = power if power >= 1 else len ( plots ) 

    nplots  = len ( plots )
    if 1 < nplots :
        import  math
        fudge_factor = math.sqrt ( 1.0 / max ( 2.0 , nplots -  1.0 ) )
        delta   = delta  * fudge_factor
        minmax  = minmax * fudge_factor
        
    save_to_db = [] 
    ## number of active plots for reweighting
    for wplot in plots  :
        
        what    = wplot.what       ## variable/function to plot/compare 
        how     = wplot.how        ## weight and/or additional cuts 
        address = wplot.address    ## address in database 
        hdata0  = wplot.data       ## original "DATA" object 
        hmc0    = wplot.mc_histo   ## original "MC"   histogram 
        ww      = wplot.w          ## relative weight 
        #
        # normailze the data
        #
        hdata = hdata0
        if isinstance ( hdata , ROOT.TH1 ) :  hdata = hdata.density ()
        
        # =====================================================================
        ## make a plot on (MC) data with the weight
        # =====================================================================
        dataset.project ( hmc0 , what , how )
        
        st   = hmc0.stat()
        mnmx = st.minmax()
        if iszero ( mnmx[0] ) :
            logger.warning ( "Reweighting: statistic goes to zero %s/``%s''" % ( st , address ) ) 
            
        # =====================================================================
        ## normalize MC
        # =====================================================================
        hmc = hmc0.density() 
        
        # =====================================================================
        ## calculate  the reweighting factor : a bit conservative (?)
        #  this is the only important line
        # =====================================================================
        
        #  try to exploit finer binning if/when possible
        if isinstance ( hmc   , ROOT.TH1 ) and \
           isinstance ( hdata , ROOT.TH1 ) and \
           len ( hmc ) >= len( hdata )  : w =  ( 1.0 / hmc ) * hdata ## NB!      
        else                            : w =  hdata / hmc           ## NB!

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
        afunc1 = allright if good1 else attention 
        afunc2 = allright if good2 else attention
        #
        message  = "%s: %24s:" % ( tag , address ) 
        message += ' ' + 'mean=%12s' % cnt.mean().toString('(%4.2f+-%4.2f)') 
        message += ' ' + afunc2 ( 'min/max=%-5.3f/%5.3f' % ( cnt.minmax()[0] , cnt.minmax()[1] ) ) 
        message += ' ' + afunc1 ( 'rms=%s[%%]' % (wvar * 100).toString('(%4.2f+-%4.2f)') ) 
        logger.info  ( message ) 
        #
        ## make decision based on the variance of weights 
        #
        mnw , mxw = cnt.minmax()
        if good : ## small variance?
            message  = "%s: No more reweights for %s" % ( tag , address )
            message += ' ' + allright (  "min/max/rms=%+3.1f/%+3.1f/%3.1f[%%]" % ( ( mnw - 1 ) * 100 ,  ( mxw - 1 ) * 100 , 100 * wvar ) )
            logger.info ( message ) 
            del w , hdata , hmc 
        else :
            save_to_db.append ( ( address , ww , hdata0 , hmc0 , hdata , hmc , w ) )
        # =====================================================================
        ## make a comparison (if needed)
        # =====================================================================
        if compare : compare ( hdata0 , hmc0 , address )
    
    ## for single reweighting 
    if 1 == nplots : power = 1
    
    if power != nplots :
        logger.info ( "%s: ``power'' is %g/#%d"  % ( tag , power , nplots  ) )

    active = [ p[0] for p in save_to_db ]    
    all    = [ p.address for p in plots ]
    for i , a in enumerate ( all ) :
        if a in active :
            if isatty () : all[i] = attention ( a )
            else         : all[i] = '*' + a + '*'
        else :
            if isatty () : all[i] = allright  ( a )
            
    logger.info ( "%s: reweights are: %s" % ( tag ,  ( ', '.join ( all ) ) ) ) 
    
    if len ( active ) != nplots :
        if database and save_to_db : 
            power += ( nplots - len ( active ) ) 
            logger.info  ("%s: ``power'' is changed to %g" %  ( tag , power ) )
            
    nactive = len ( active )  
    while database and save_to_db :

        entry = save_to_db.pop() 
        
        address, ww , hd0, hm0, hd , hm , weight = entry  

        eff_exp = 1.0  / power
        
        eff_exp = 0.95 / ( 1.0 * nactive ) ** 0.5
        
        eff_exp = 0.5  if 1 < nactive else 1 
        
        if 1 != ww :
            eff_exp *= ww
            logger.info  ("%s: apply ``effective exponent'' of %.3f for ``%s''" % ( tag , eff_exp  , address ) )
            
        if 1 != eff_exp and 0 < eff_exp : 
            weight = weight ** eff_exp

        ## print 'WEIGHT stat', eff_exp, weight.stat()
        
        ## hmmmm... needed ? yes! 
        #if 1 < power : weight = weight ** ( 1.0 / power )
        
        ## relative importance
        #if 1 != ww :
        #    logger.info  ("%s: apply ``relative importance factor'' of %.3g for ``'%s'" % ( tag , ww , address ) )
        #    weight = weight ** ww 

        with DBASE.open ( database ) as db :
            
            db[address] = db.get( address , [] ) + [ weight ]
            
            if debug :
                addr        = address + ':REWEIGHTING'
                db [ addr ] = db.get ( addr , [] ) + list ( entry[2:] )
                
        del hd0, hm0 , hd , hm , weight , entry 
        
    return active 

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
# The END 
# =============================================================================
