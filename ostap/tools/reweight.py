#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with utilities for reweigting 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-05-10
#  
#                    $Revision$
#  Last modification $Date$
#  by                $Author$
# =============================================================================
"""Module with utilities for reweighting"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'Weight'      ,
    'makeWeights' ,
    ) 
# =============================================================================
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.tools.reweigh' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
logger.info ( 'Set of utitilities for re-weigthing')
from   ostap.core.pyrouts import VE, SE, iszero 
import ostap.io.zipshelve as     DBASE ## needed to store the weights&histos 
# =============================================================================
## @class Weight
#  helper class for semiautomatic reweighting of data 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-05-10
class Weight(object) :
    """Helper class for semi-automatic reweighting of data

    It reads various ``weigth components'' from DBASE and calculates
    the ``global'' event weight via smiple accumulation of weights
    
    Simplest case: calculate weights for one variable.
    
    factors = [ ( fun , name ) ] 
    where ``fun'' is a function to get variabel from TTree/RooDataSet 
    and   ``name'' is address in data base with reweighting information 
    e.g.
    
    # factors = [ ( lambda s : s.pt  , 'pt-data' ) ]
    
    It assumes, that DATA base has an entry with name 'pt-data', such
    that
    
    # funcs = database['pt-data']
    
    Here ``func'' will be function/callable or list of functions/calables
    with weights:
    
    # pt = fun ( entry ) 
    # w  = 1.0    
    # w *= func ( pt )              ## for the single function/callable
    # for f in funcs : w*= f ( pt ) ## for list of functions/callables  
    
    quantity ``w'' will be am event  weight
    
    if list ``factors'' constains more than one entry.,
    to each entry weigth is calculated independently and total weight
    is a product of individual weights..

    For smiplistic case one can put some storable function or function object
    # db['pt-data'] = math.exp
    # db['pt-data'] = [ math.exp , math.cosh ] ## List of object also works:
    
    For realistic case the weights usually stored as (a list of) 1D or 2D historgams
    # db['pt-data'] = [h1,h2,h3,...,hn]
    where h1,...,hn are iteratively obtained histograms with corrections, such as
    the total correction is calculated by the product of them
    """
    def __init__ ( self                   ,
                   dbase   = "weights.db" , ## the name of data base with the weights 
                   factors = []           ) :
        
        #
        ## make some statistic
        #
        self._counter = SE ()
        self._nzeroes = 0 

        self.vars = [] 
        if not factors : return

        ## open database 
        with DBASE.open ( dbase , 'r' ) as db : ## READONLY
            
            for k in db :
                e = db[k]
                if hasattr ( e , '__len__' ) :  
                    logger.debug( 'DBASE "%.15s" key "%.15s" #%d' % ( dbase ,  k, len( e ) ) ) 
                
            ## loop over the weighting factors and build the function
            for f in factors :

                funval  = f[0]  ## accessor to the variable 
                funname = f[1]  ## address  in database 

                if isinstance ( funval , str ) :
                    varnam = funval 
                    funval = lambda s : getattr ( s , varnam )
                    
                ## 
                functions  = db.get ( funname , [] ) ## db[ funname ]
                if not functions :
                    logger.warning('No reweighting is available for %s, skip it' % funname )

                merge = True
                if 2 < len ( f ) : merge = f[2] 
                
                if not isinstance (  functions , ( list , tuple ) ) :
                    functions = [ functions ]
                    
                ## merge list of functions into single function 
                if merge and 1 < len ( functions)  : 
                
                    single_func = functions[0] * functions [1] 
                    
                    for fun in functions [2:] :
                        single_func *= fun
                            
                    functions  = [ single_func ]
                    
                self.vars += [ ( funname , funval , functions , SE() ) ]  


    ## get the statistic of weights 
    def stat    ( self ) :
        "Get the statistic of used weights"
        return self._counter
    ## number of zero weights 
    def nZeroes ( self ) :
        "Number of null weights"
        return self._nzeroes
    
    ## calculate the weight for the given "event"
    def __call__ ( self , s ) :
        """
        Calculate the weigth for the given ``event''
        (record in TTree/TChain or RooDataSet)
        
        """

        ## initialize the weight 
        weight  = VE(1,0) 

        ## loop over functions 
        for i in self.vars :
            
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
        
        self._counter += vw 
        if iszero ( vw ) : self._nzeroes += 1
            
        return vw 

# =============================================================================
## make one re-weighting iteration 
#  and reweight "MC"-data set to looks as "data"(reference) dataset
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-05-10
def makeWeights  ( dataset                 ,
                   plots    = []           ,
                   database = "weights.db" ,
                   compare  = None         ,    ## comparison function 
                   delta    = 0.001        ,    ## delta for weigth variance 
                   debug    = True         ) :  ## save intermediate information in DB 

    more = False 
    ## loop over plots 
    for r in plots  :

        what    = r [0]         ## variable/function to plot/compare 
        how     = r [1]         ## weight or additional cuts 
        address = r [2]         ## address in database 
        hdata0  = r [3]                          .clone () ## original "DATA" histogram
        hmc0    = r [4] if 4 < len(r) else hdata0.clone () ## original "MC"   histogram 

        #
        ## black magic to take into account the difference in bins and normalizations
        #
        hdata = hdata0 
        if hasattr ( hdata , 'rescale_bins' ) : 
            hdata = hdata.rescale_bins ( 1.0   )
            
        ## normalize the data:
        hmean = None 
        if hasattr ( hdata , 'mean' ) and hasattr ( hdata , '__idiv__' ) :

            ## normalization point
            hmean  = hdata.mean()
            #
            if isinstance ( hdata , ROOT.TH2 ) : hdata /= hdata ( *hmean )
            else                               : hdata /= hdata (  hmean )

        #
        ## make a plot on (MC) data with the weight
        # 
        dataset.project ( hmc0 , what , how )

        st   = hmc0.stat()
        mnmx = st.minmax()
        if iszero ( mnmx[0] ) :
            logger.warning ( 'Statistic goes to zero %s/"%s"' % ( st , address ) ) 
            
        #
        ## black magic to take into account the difference in bins and normalizations
        # 
        hmc = hmc0.rescale_bins ( 1.0 )
        
        if hmean is None : pass
        else             :
            if isinstance ( hmc , ROOT.TH2 ) : hmc /= hmc ( *hmean )
            else                             : hmc /= hmc (  hmean )

        #
        ## calculate  the reweighting factor : a bit conservative (?)
        power = min ( 2.0 , len ( plots ) )                   ## NB!
        #  this is the only important line
        #  try to exploit finer binning if possible 
        if len ( hmc ) >= len( hdata )  : 
            w     = ( ( 1.0   / hmc ) * hdata ) ** ( 1.0 / power )  ## NB!
        else :
            w     = ( ( hdata / hmc )         ) ** ( 1.0 / power )  ## NB!
            
        #
        ## get the statistics of weights 
        #
        cnt  = w.stat()
        mnmx = cnt.minmax()
        if not mnmx [0] <= 1 <= mnmx[1] : w /= cnt.mean().value()
        cnt  = w.stat()
        #
        wvar = cnt.rms()/cnt.mean()
        logger.info ( 'Reweighting "%-.15s: Mean/minmax:%s/(%.4f,%.4f) Vars:%s[%%]' %
                      ( address         ,
                        cnt.mean()      ,
                        cnt.minmax()[0] ,
                        cnt.minmax()[1] , wvar * 100 ) ) 
        #
        ## make decision based on variance of weights 
        #
        if wvar.value() <= delta / len ( plots ) : ## small variance? 
            save = False
            logger.info("No more reweighting for %s [%.3f%%]" %  ( address , wvar * 100 ) ) 
        else            :
            save = True 

        #
        ## make a comparison (if needed)
        # 
        if compare :
            compare ( hdata0 , hmc0 , address )
        
        ## update data base 
        if save and database and address :
            with DBASE.open ( database ) as db :

                db[address] = db.get( address , [] ) + [ w ]
                
                if debug :
                    addr        = address + ':REWEIGHTING'
                    entry       = ( hdata0 , hmc0 , hdata , hmc , w ) 
                    db [ addr ] = db.get ( addr , [] ) + [ entry ]
                    
        ## 
        more = more or save

        del hdata0, hmc0, hdata, hmc, w  
        
    return more

# =============================================================================
## some simple comparsion 
def hCompare ( data , mc , title = '' , spline = True ) :

    
    if not isinstance ( data , ( ROOT.TH1D , ROOT.TH1F ) ) : return
    if not isinstance ( mc   , ( ROOT.TH1D , ROOT.TH1F ) ) : return

    data.cmp_prnt( mc )

    hd  = data.rescale_bins ( 1 ) 
    hm  = mc  .rescale_bins ( 1 ) 

        
# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
