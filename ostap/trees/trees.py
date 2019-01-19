#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with decoration of Tree/Chain objects for efficient use in python
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Decoration of Tree/Chain objects for efficient use in python"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (  
    'Chain' , ## helper class , needed for multiprocessing 
    'Tree'  , ## helper class , needed for multiprocessing 
  ) 
# =============================================================================
import ROOT
from ostap.core.core    import std , Ostap, VE, hID
from ostap.logger.utils import multicolumn
from ostap.utils.basic  import terminal_size, isatty 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger, allright,  attention 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.trees.trees' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Some useful decorations for Tree/Chain objects')
# =============================================================================
import ostap.trees.cuts
# =============================================================================
_large = 2**64
# =============================================================================
from ostap.core.core import valid_pointer
# =============================================================================
## check validity/emptiness  of TTree/TChain
#  require non-zero poniter and non-empty Tree/Chain
def _tt_nonzero_ ( tree ) :
    """Check validity/emptiness  of TTree/TChain
    - require non-zero poniter and non-empty Tree/Chain
    """
    return valid_pointer ( tree ) and 0 < len ( tree )
ROOT.TTree .__nonzero__ = _tt_nonzero_
ROOT.TChain.__nonzero__ = _tt_nonzero_
ROOT.TTree .__bool__    = _tt_nonzero_
ROOT.TChain.__bool__    = _tt_nonzero_

# =============================================================================
## Iterator over ``good events'' in TTree/TChain:
#  @code 
#    >>> tree = ... # get the tree
#    >>> for i in tree.withCuts ( 'pt>5' ) : print i.y
#  @endcode
#  @attention: TTree::GetEntry is already invoked for accepted events,
#              no need in second call
#  @see Analysis::PyIterator
#  @see Ostap::Formula
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-06
def _iter_cuts_ ( self , cuts , first = 0 , last = _large , progress = False ) :
    """Iterator over ``good events'' in TTree/TChain:
    
    >>> tree = ... # get the tree
    >>> for i in tree.withCuts ( 'pt>5' ) : print i.y
    
    Attention: TTree::GetEntry is already invoked for accepted events,
    no need in second call 
    """
    #
    last = min ( last , len ( self )  )
    
    pit = Ostap.PyIterator ( self , cuts , first , last )
    if not pit.ok() : raise TypeError ( "Invalid Formula: %s" % cuts )
    #
    from ostap.utils.progress_bar import ProgressBar 
    with ProgressBar ( min_value = first        ,
                       max_value = last         ,
                       silent    = not progress ) as bar :
        
        step = 13.0 * max ( bar.width , 101 ) / ( last - first ) 
        
        _t = pit.tree()
        _o = _t 
        while valid_pointer ( _t ) :

            yield _t
            _t      = pit.next()             ## advance to the next entry  

            if progress : 
                current = pit.current() - 1  ## get the current entry index 
                if not _t                          \
                       or  _t != _o                \
                       or current - first   < 120  \
                       or last    - current < 120  \
                       or 0 == current % 100000    \
                       or 0 == int ( step * ( current - first ) ) % 5  :
                    
                    ## show progress bar 
                    bar.update_amount( current )
                    _o = _t
                    
        if progress : bar.update_amount( last ) 

    del pit
    self.GetEntry(0)
    

ROOT.TTree .withCuts  = _iter_cuts_ 
ROOT.TChain.withCuts  = _iter_cuts_ 

ROOT.TTree. __len__   = lambda s : s.GetEntries()

# =============================================================================
## Iterator over ``good events'' in TTree/TChain:
#  @code 
#    >>> tree = ... # get the tree
#    >>> for i in tree( 0, 100, 'pt>5' ) : print i.y
#  @endcode
#  @see Ostap::PyIterator
#  @see Ostap::Formula
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-06
def _tc_call_ ( self , first = 0 , last = -1  , cuts = None , progress = False ) :
    """Iterator over ``good events'' in TTree/TChain:
    
    >>> tree = ... # get the tree
    >>> for i in tree(0, 100 , 'pt>5' ) : print i.y
    
    """
    #
    if last < 0 : last = ROOT.Tree.kMaxEntries
    
    last = min ( last , len ( self )  )

    from ostap.utils.progress_bar import ProgressBar 
    with ProgressBar ( min_value = first        ,
                       max_value = last         ,
                       silent    = not progress ) as bar :
        
        step = 13.0 * max ( bar.width , 101 ) / ( last - first ) 

        pit = 1 
        if cuts :
            
            pit = Ostap.PyIterator ( self , cuts , first , last )
            if not pit.ok() : raise TypeError ( "Invalid Formula: %s" % cuts )
            #
            
            _t = pit.tree()
            _o = _t 
            while valid_pointer ( _t ) :
                
                yield _t                         ## YIELD 
                _t      = pit.next()             ## advance to the next entry  
                
                if progress : 
                    current = pit.current() - 1  ## get the current entry index 
                    if not _t                          \
                           or  _t != _o                \
                           or current - first   < 120  \
                           or last    - current < 120  \
                           or 0 == current % 100000    \
                           or 0 == int ( step * ( current - first ) ) % 5  :
                        
                    ## show progress bar 
                        bar.update_amount( current )
                        _o = _t
        else :
            
            ## just explicit loop 
            for current in range ( first , last + 1 ) :
                
                if progress :
                    if     current - first   < 120  \
                           or last - current < 120  \
                           or 0 == current % 100000 \
                           or 0 == int ( step * ( current - first ) ) % 5  :
                        
                        bar.update_amount( current )
                        
                if 0 >= self.GetEntry ( current ) : break
                yield self                         ## YIELD! 
                
                    
        if progress : bar.update_amount( last ) 

    del pit
    self.GetEntry(0)
    

ROOT.TTree .__call__  = _tc_call_ 
ROOT.TChain.__call__  = _tc_call_

# =============================================================================
## help project method for ROOT-trees and chains 
#
#  @code 
#    >>> h1   = ROOT.TH1D(... )
#    >>> tree.Project ( h1.GetName() , 'm', 'chi2<10' ) ## standart ROOT 
#    
#    >>> h1   = ROOT.TH1D(... )
#    >>> tree.project ( h1.GetName() , 'm', 'chi2<10' ) ## ditto 
#    
#    >>> h1   = ROOT.TH1D(... )
#    >>> tree.project ( h1           , 'm', 'chi2<10' ) ## use histo
# 
#    ## make invididual projections of 'm1' and 'm2' and make a sum of distributions
#    >>> h1   = ROOT.TH1D(... )
#    >>> tree.project ( h1           , ['m1','m2'] , 'chi2<10' ) ## use histo
#
#    ## make invididual projections of 'm1' and 'm2' and make a sum of distributions
#    >>> h1   = ROOT.TH1D(... )
#    >>> tree.project ( h1           , "m1,m2"     , 'chi2<10' )
#    >>> tree.project ( h1           , "m1;m2"     , 'chi2<10' )
#  @endcode
#
#  @param tree   the tree
#  @param histo  the histogram or histogram name 
#  @param what variable/expression to be projected.
#              It could be a list/tuple of variables/expressions or just a comma-separated expression
#  @param cuts expression for cuts/weights
#  @see TTree::Project
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-07-06
def _tt_project_ ( tree               ,
                   histo              ,
                   what               ,
                   cuts       = ''    ,
                   options    = ''    ,
                   nentries   = -1    ,
                   firstentry =  0    ,
                   silent     = False ) :
    """Helper project method
    
    >>> tree = ...
    
    >>> h1   = ROOT.TH1D(... )
    >>> tree.Project ( h1.GetName() , 'm', 'chi2<10' ) ## standart ROOT 
    
    >>> h1   = ROOT.TH1D(... )
    >>> tree.project ( h1.GetName() , 'm', 'chi2<10' ) ## ditto 
    
    >>> h1   = ROOT.TH1D(... )
    >>> tree.project ( h1           ,  'm', 'chi2<10' ) ## use histo

    ## make invididual projections of m1 and m2 and make a sum of distributions
    >>> h1   = ROOT.TH1D(... )
    >>> tree.project ( h1           , ('m1','m2') , 'chi2<10' ) ## two variables 
    >>> tree.project ( h1           , 'm1,m2'     , 'chi2<10' ) ## ditto
    >>> tree.project ( h1           , 'm1;m2'     , 'chi2<10' ) ## ditto
    
    - tree  : the tree
    - histo : the histogram (or histogram name)
    - what  : variable/expression to project. It can be expression or list/tuple of expression or comma (or semicolumn) separated expression
    - cuts  : selection criteria/weights 
    """
    #
    
    ## if nentries < 0 :
    nentries = ROOT.TTree.kMaxEntries
        
    args = options , nentries , firstentry, silent
    ## 
    hname = histo 
    if   hasattr    ( histo , 'GetName' ) : hname = histo.GetName()
    ## elif isinstance ( histo , str       ) : 
    ##    h = ROOT.gROOT.FindObject ( hname )
    ##    if h : histo = h

    ## reset it!
    if histo and isinstance ( histo , ROOT.TH1  ) : histo.Reset()
    #
    if isinstance ( cuts  , ROOT.TCut ) : cuts = str(cuts) 
    if not what : return 0, histo
    #
    ## trivial 1-item list
    if hasattr ( what , '__len__' ) and 1 == len ( what ) and not isinstance ( what , (str, ROOT.TCut) ): 
        what = what[0]

    ## check for semicolumn-separated list of expressions:
    if isinstance ( what , str ) and ';' in what : 
        what = what.split(';')
        if 1 == len(what) : what = what[0]

    ## check for comma-separated list of expressions:
    if isinstance ( what , str ) and ',' in what :
        if '(' in what and ')' in what : pass 
        else :
            what = what.split(',')
            if 1 == len( what ) : what = what[0]

    #
    if   isinstance ( what  , str       ) : what =     what 
    elif isinstance ( what  , ROOT.TCut ) : what = str(what)  
    elif isinstance ( histo , ROOT.TH1  ) : 
        rr = 0 
        hh = histo.clone()
        for v in what :
            r , h  = _tt_project_ ( tree , hh , v , cuts , options , *args )
            rr    += r
            histo += h
        hh.Delete()
        del hh 
        return rr , histo
    elif isinstance ( histo , str ) :
        ## process the head of the list: the first call creates the histo... 
        rr, hh =  _tt_project_ ( tree , histo , what[0] , cuts , *args )
        histo  = hh
        if 1 == len ( what )   : return rr , histo
        # normal processing of the tail of the list using created historgam 
        hh      = histo.clone()
        r1 , h1 = _tt_project_ ( tree , hh , what[1:] , cuts , *args )
        rr     += r1
        histo  += h1
        hh.Delete()
        del hh, h1 
        return rr , histo

    ## the basic case 
    from ostap.core.core import ROOTCWD
    with ROOTCWD() :
        ROOT.gROOT.cd ()
        ## make projection 
        result = tree.Project ( hname , what , cuts , *args[:-1] )
        if   isinstance ( histo , ROOT.TH1 ) : return result, histo
        elif isinstance ( histo , str      ) :
            h = ROOT.gROOT.FindObject ( hname )
            if h : return result, h

    return result, histo

ROOT.TTree .project = _tt_project_
ROOT.TChain.project = _tt_project_

# =============================================================================
## get the statistic for certain expression in Tree/Dataset
#  @code
#  tree  = ... 
#  stat1 = tree.statVar( 'S_sw/effic' )
#  stat2 = tree.statVar( 'S_sw/effic' ,'pt>1000')
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-09-15
def _stat_var_ ( tree , expression , *cuts ) :
    """Get a statistic for the  expression in Tree/Dataset
    
    >>> tree  = ... 
    >>> stat1 = tree.statVar ( 'S_sw/effic' )
    >>> stat2 = tree.statVar ( 'S_sw/effic' ,'pt>1000')
    
    """

    if isinstance ( expression , str ) :
        from ostap.core.core import split_string
        expression = split_string ( expression , ',;:' ) 
        
    if 1 != len ( expression ) :
        return _stat_vars_ ( tree ,  expression , *cuts )
    
    expression = expression[0] 
    
    return Ostap.StatVar.statVar ( tree , expression , *cuts )
    
ROOT.TTree     . statVar = _stat_var_
ROOT.TChain    . statVar = _stat_var_

# =============================================================================
## get the statistic for certain expressions in Tree/Dataset
#  @code
#  tree  = ... 
#  stat1 = tree.statVars( [ 'S_sw/effic', 'pt1' , 'pt2' ] ) 
#  stat2 = tree.statVars( [ 'S_sw/effic', 'pt1' , 'pt2' ] , 'mass>10') 
#  @endcode
#  It is more efficient than getting statistics individually for each expression
#  @see Ostap::Math::StatVar 
#  @see Ostap::Math::StatVar::statVars 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2018-11-03
def _stat_vars_ ( tree , expressions , *cuts ) :
    """Get the statistic for certain expressions in Tree/Dataset
    >>> tree  = ... 
    >>> stat1 = tree.statVars( [ 'S_sw/effic', 'pt1' , 'pt2' ] ) 
    >>> stat2 = tree.statVars( [ 'S_sw/effic', 'pt1' , 'pt2' ] , 'mass>10') 
    - It is more efficient than getting statistics individually for each expression
    - see Ostap::Math::StatVar
    - see Ostap::Math::StatVar::statVars 
    """
    from ostap.core.core import std, strings, split_string, WSE 
    
    if isinstance ( expressions , str ) :
        expressions = split_string ( expressions , ',;:' ) 

    if not expressions : return {}    
    if 1 == len ( expressions ) :
        return _stat_var_ ( tree , expressions[0] , *cuts )
    
    vct = strings ( *expressions )
    res = std.vector(WSE)() 

    ll  = Ostap.StatVar.statVars ( tree , res , vct , *cuts )
    assert res.size() == vct.size(), 'Invalid size of structures!'

    N = res.size()
    results = {} 

    for i in range(N) :
        results[ vct[i] ] = WSE ( res[i] ) 

    return results 

ROOT.TTree     . statVars = _stat_vars_
ROOT.TChain    . statVars = _stat_vars_

# =============================================================================
## get the statistic for pair of expressions in Tree/Dataset
#  @code
#  tree  = ...
#  stat1 , stat2 , cov2 , len = tree.statCov ( 'x' , 'y' )
#  # apply some cuts 
#  stat1 , stat2 , cov2 , len = tree.statCov ( 'x' , 'y' , 'z>0' )
#  # use only subset of events
#  stat1 , stat2 , cov2 , len = tree.statCov ( 'x' , 'y' , 'z>0' , 100 , 10000 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-09-15
def _stat_cov_ ( tree        ,
                 expression1 ,
                 expression2 ,
                 cuts = ''   , *args  ) :
    """Get the statistic for pair of expressions in Tree/Dataset
    
    >>>  tree  = ...
    >>>  stat1 , stat2 , cov2 , len = tree.statCov( 'x' , 'y' )
    
    Apply some cuts:
    >>> stat1 , stat2 , cov2 , len = tree.statCov( 'x' , 'y' , 'z>0' )
    
    Use only subset of events
    >>> stat1 , stat2 , cov2 , len = tree.statCov( 'x' , 'y' , 'z>0' , 100 , 10000 )
    """
    import ostap.math.linalg 
    stat1  = Ostap.WStatEntity       ()
    stat2  = Ostap.WStatEntity       ()
    cov2   = Ostap.Math.SymMatrix2x2 ()

    if cuts : 
        length = Ostap.StatVar.statCov ( tree        ,
                                         expression1 ,
                                         expression2 ,
                                         cuts        ,
                                         stat1       ,
                                         stat2       ,
                                         cov2        , 
                                         *args       )
    else :
        length = Ostap.StatVar.statCov ( tree        ,
                                         expression1 ,
                                         expression2 ,
                                         stat1       ,
                                         stat2       ,
                                         cov2        ,
                                         *args       )
        
    return stat1 , stat2 , cov2, length

ROOT.TTree     . statCov = _stat_cov_
ROOT.TChain    . statCov = _stat_cov_

# =============================================================================
## get the statistic for the list of expressions 
#  @code
#  tree  = ...
#  stats , cov2 , len = tree.statCovs( ['x' , 'y'] )
#  # apply some cuts 
#  stats , cov2 , len = tree.statCovs( [ 'x' , 'y' , 'z'] , 'z>0' )
#  # use only subset of events
#  stats , cov2 , len = tree.statCovs( [ 'x' , 'y' , 'z' ], 'z>0' , 100 , 10000 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2017-02-19
def _stat_covs_ ( tree        ,
                  expressions ,
                  cuts = ''   , *args  ) :
    """Get the statistic for the list of expressions 
    >>> tree  = ...
    >>> stats , cov2 , len = tree.statCovs( ['x' , 'y'] )
    Apply some cuts 
    >>> stats , cov2 , len = tree.statCovs( [ 'x' , 'y' , 'z'] , 'z>0' )
    Use only subset of events
    >>> stats , cov2 , len = tree.statCovs( [ 'x' , 'y' , 'z' ], 'z>0' , 100 , 10000 )
    """
    ##
    if isinstance ( expressions , str ) : expressions = [ expressions ]
    ##
    
    import ostap.math.linalg 
    
    _SV    = std.vector('std::string')
    _vars  = _SV()
    vars   = expressions
    for e in vars : _vars.push_back( e )
    
    WSE    = Ostap.WStatEntity
    _WV    = std.vector( WSE )
    _stats = _WV()
    _DV    =  std.vector('double')
    _cov2  = _DV()

    if cuts : 
        length = Ostap.StatVar._statCov ( tree   ,
                                          _vars  ,
                                          cuts   ,
                                          _stats ,
                                          _cov2  ,
                                          *args  ) 
    else :
        length = Ostap.StatVar._statCov ( tree   ,
                                          _vars  ,
                                          _stats ,
                                          _cov2  ,
                                          *args  )
        
    l = len(_vars)
    if 0 == length : 
        return None , None , 0 
    elif l != len ( _stats ) or l*(l+1)/2 != len( _cov2 ):
        logger.error("statCovs: unexpected output %d/%s/%s" % ( l           ,
                                                                len(_stats) ,
                                                                len(_cov2 ) ) )
        return None, None, length


    ## get the statistics of variables
    stats = tuple ( [ WSE(s) for s in _stats ] ) 
    
    import ostap.math.linalg
    COV2 = Ostap.Math.SymMatrix ( l )
    cov2 = COV2 () 

    for i in range( l ) :
        for j in range ( i + 1 ) :
            ij = i * ( i + 1 ) / 2 + j
            cov2[i,j] = _cov2[ ij ]
            
    return stats, cov2 , length

ROOT.TTree     . statCovs = _stat_covs_
ROOT.TChain    . statCovs = _stat_covs_

# =============================================================================
## Get min/max for the certain variable in chain/tree
#  @code  
#  chain = ...
#  mn,mx = chain.vminmax('pt')
#  mn,mx = chain.vminmax('pt','y>3')
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-09-19
def _tc_minmax_ ( tree , var , cuts = '' , delta = 0.0 )  :
    """Get min/max for the certain variable in chain/tree
    >>> chain = ...
    >>> mn,mx = chain.vminmax('pt')
    >>> mn,mx = chain.vminmax('pt','y>3')
    """
    if hasattr ( tree , 'pstatVar' ) : 
        if cuts : s = tree.pstatVar ( var , cuts )
        else    : s = tree.pstatVar ( var )
    else :
        if cuts : s = tree.statVar  ( var , cuts )
        else    : s = tree.statVar  ( var )

    mn,mx = s.minmax()
    if mn < mn and 0.0 < delta :
        dx   = delta * 1.0 * ( mx - mn )  
        mx  += dx   
        mn  -= dx   
    return mn , mx

ROOT.TTree     . vminmax = _tc_minmax_
ROOT.TChain    . vminmax = _tc_minmax_

# =============================================================================
## @var _h_one_
#  special helper histogram for summation
_h_one_ = ROOT.TH1D( hID() , '' , 3 , -1 , 2 ) ; _h_one_.Sumw2()
# =============================================================================
## make a sum over expression in Tree/Dataset
#
#  @code
#
#  >>> dataset = ...
#  ## get corrected number of events 
#  >>> n_corr  = dataset.sumVar ( "S_sw/effic" )
#
#  @endcode
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-09-15
def _sum_var_old_ ( tree , expression ) :
    """Make a sum over expression in Tree/Dataset
    
    >>> dataset = ...
    ## get corrected number of signale events  
    >>> n_corr  = dataset.sumVar_( 'S_sw/effic' )
    
    """
    _h_one_.Reset() 
    tree.project ( _h_one_ , '1' , expression )
    return _h_one_.accumulate()

    
ROOT.TTree      . sumVar_ = _sum_var_old_
ROOT.TChain     . sumVar_ = _sum_var_old_

# =============================================================================
## make a sum over expression in Tree/Dataset
#
#  @code
#
#  >>> dataset = ...
#
#  ## get corrected number of events 
#  >>> n_corr     = dataset.sumVar ( "S_sw/effic" )
#
#  ## get corrected number of events 
#  >>> n_corr_pt  = dataset.sumVar ( "S_sw/effic" , 'pt>1')
#
#  @endcode
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-09-15
def _sum_var_ ( tree , expression , *cuts ) :
    """Make a sum over expression in Tree/Dataset
    
    >>> dataset = ...
    ## get corrected number of signal events  
    >>> n_corr     = dataset.sumVar ( 'S_sw/effic' )
    
    ## get corrected number of signal events  
    >>> n_corr_pt  = dataset.sumVar ( 'S_sw/effic' , 'pt>1')
    
    """
    ## if hasattr ( tree , 'pStatVar' ) : w = tree.pStatVar ( expression , *cuts )
    ## else                             : w = tree. statVar ( expression , *cuts )
    w = tree. statVar ( expression , *cuts )
    ##
    return VE ( w.sum() , w.sum2() )

ROOT.TTree      . sumVar = _sum_var_
ROOT.TChain     . sumVar = _sum_var_

# =============================================================================
## get the leaves for the given tree/chain
#  @see TTree
#  @code
#
#  >>> tree = ...
#  >>> lst = tree.leaves()
#  >>> for l in lst : print l
#
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-02-04
def _rt_leaves_ ( t , pattern = '' , *args ) :
    """ Get the list of leaves names        
    
    >>> tree = ...
    >>> lst = tree.leaves()
    >>> for l in lst : print l
    >>> lst = tree.leaves( '.*(muon).*', re.I )
    >>> for l in lst : print l
    """
    vlst =  t.GetListOfLeaves()

    if not vlst : return tuple()
    if pattern :        
        try : 
            import re
            c  =  re.compile ( pattern , *args )
            lst  = [ v.GetName() for v in vlst if c.match ( v.GetName () ) ]
            lst.sort()
            return tuple ( lst ) 
        except :
            logger.error ('leaves: exception is caught, skip it' , exc_info = True ) 
            
    lst  = [ v.GetName() for v in vlst  ]
    lst.sort()
    return tuple ( lst ) 


    
    _lst = [ l.GetName() for l in _lst ] 
    _lst.sort()
    return tuple( _lst ) 

ROOT.TTree.leaves   = _rt_leaves_

# ==============================================================================
## Get the leaf with the certain name 
def _rt_leaf_ ( tree , leaf ) :
    """Get the leaf with certain name:
    >>> tree = ...
    >>> l = tree.leaf('pt') 
    """
    lst = tree.GetListOfLeaves()
    for i in lst :
        if leaf == i.GetName() : return i
    return None

ROOT.TTree.leaf   = _rt_leaf_


# =============================================================================
## get the branches for the given tree/chain
#  @see TTree
#  @code
#  >>> tree = ...
#  >>> lst = tree.branches()
#  >>> for b in lst : print b
#  >>> lst = tree.branches( '.*(Muon).*' , re.I )
#  >>> for b in lst : print b
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-02-04
def _rt_branches_ ( t , pattern = '' , *args ) :
    """Get the list of branch names
    
    >>> tree = ...
    >>> lst = tree.branches()
    >>> for b in lst : print b
    >>> lst = tree.branches( '.*(Muon).*' , re.I )
    >>> for b in lst : print b
    
    """
    vlst =  t.GetListOfBranches()
    if not vlst : return tuple()

    if pattern :        
        try : 
            import re
            c  =  re.compile ( pattern , *args )
            lst  = [ v.GetName() for v in vlst if c.match ( v.GetName () ) ]
            lst.sort()
            return tuple ( lst ) 
        except :
            logger.error ('branches: exception is caught, skip it' , exc_info = True ) 
            
    lst  = [ v.GetName() for v in vlst  ]
    lst.sort()
    return tuple ( lst ) 


ROOT.TTree.branches = _rt_branches_


# =============================================================================
## simplified printout for TTree/TChain
#  @see TTree
#  @code
#
#  >>> tree = ...
#  >>> print tree
#
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-02-04
def _rt_print_ ( t ) :
    """Simplified print out for tree/chain

    >>> tree = ...
    >>> print tree
    """
    ##
    res = "Name: %s Enries/#%d" %  ( t.GetName() , t.GetEntries() ) 
    if hasattr ( t , 'GetNtrees' ) : res += " Chain/#%d " %       t.GetNtrees()
    ##
    _l          = t.leaves   ()
    res        += "\nLeaves:\n%s"    % multicolumn ( list( _l ) , indent = 2 , pad = 1 )

    ## collect non-trivial branches 
    _b          = t.branches ()
    
    _bs = set  ( _b )
    _ls = set  ( _l )
    _b  = list ( _bs - _ls ) 
    _b . sort() 
    if _b : res += "\nNon-trivial branches:\n%s" % multicolumn ( _b , indent = 2 ,  pad = 1 )

    return res.replace ('\n','\n# ') 

ROOT.TTree.__repr__ = _rt_print_
ROOT.TTree.__str__  = _rt_print_
ROOT.TTree.pprint   = _rt_print_


# =============================================================================

__std_ints  = ( 'char' , 'short'  , 'int' , 'long' , 'long long' )
__std_uints = tuple ( [ 'unsigned ' + i for i in __std_ints ] )
__std_types = ( 'bool' , 'float' , 'double' , 'long double' ) + __std_ints + __std_uints 
__scalars   = __std_types + ( 'Bool_t'    ,
                              'Char_t'    ,
                              'UChar_t'   ,
                              'Short_t'   ,
                              'UShort_t'  , 
                              'Int_t'     ,
                              'UInt_t'    ,
                              'Float_t'   ,
                              'Double_t'  ,
                              'Long64_t'  ,
                              'ULong64_t' )
__vectors  =  tuple ( [ 'vector<' + i + '>' for i in __scalars ] )
__types    =  list ( __scalars)  + list ( __vectors )
tmp = set()
for t in __types :
    while 0 <= t.find ( 2 * ' ' ) : t = t.replace ( 2 * ' ' , ' ' )
    tmp.add ( t )
__types    = tuple ( tmp )
del tmp 


def __in_types ( t ) :
    while 0 <= t.find ( 2 * ' ' ) : t = t.replace ( 2 * ' ' , ' ' )
    return t in __types 
    
# ==============================================================================
## print tree as table 
def _rt_table_0_ ( tree , pattern = None , cuts = '' , *args ) :
    """
    """
    ## get list of branches 
    brs = tree.leaves ( pattern )
    if 'TObject' in brs :
        brs = list  ( brs )
        brs.remove  ( 'TObject' ) 
        brs = tuple ( brs )

    ## collect information
    _vars = []
    if hasattr ( tree , 'pstatVar' ) : s0 = tree.pstatVar ( '1' , cuts , *args )
    else                             : s0 = tree. statVar ( '1' , cuts , *args )
    n0    = s0.nEntries  ()

    ## no entries passed the cuts 
    brs   = () if 0 == n0 else brs

    ## from ostap.utils.progress_bar import ProgressBar 
    ## progress = ProgressBar ( max_value = len ( brs ) )

    ## get the interesting branches:
    
    bbs     = []
    selvars = 0 
    for b in brs :
        
        l = tree.leaf ( b )

        if not l :
            logger.warning ("table: can't get the leaf  \"%s\"" % b )
            continue
        
        tn       = l.GetTypeName ()
        typename = tn

        selvars += 1 
        if not __in_types ( tn ) : continue
        
        bbs.append ( b ) 

    report  = '# %s("%s","%s"' % ( tree.__class__.__name__ , tree.GetName  () , tree.GetTitle () )
    if tree.GetDirectory() :  report += ',%s' % tree.GetDirectory().GetName()
    report += ')'
    if tree.topdir :  report += '\n# top-dir:%s' % tree.topdir.GetName()
    report += '\n# ' + allright ( '%d entries; %d/%d variables (selected/total)' % ( len ( tree  ) ,
                                                                                     selvars       ,
                                                                                     len ( tree.branches() ) ) )

    if hasattr ( tree , 'pstatVar' ) : bbstats = tree.pstatVar ( bbs , cuts , *args )
    else                             : bbstats = tree. statVar ( bbs , cuts , *args )

    from ostap.stats.counters import WSE 
    if isinstance ( bbstats , WSE )  : bbstats = { bbs[0] : bbstats } 
    
    for b in brs :

        
        l = tree.leaf ( b )

        if not l :
            logger.warning ("table: can't get the leaf  \"%s\"" % b )
            continue
        
        tn       = l.GetTypeName ()
        typename = tn
        
        br = l.GetBranch()
        if br :  
            n = br.GetTitle()
            typename = '%s %s '
            p = n.find('[')
            if  0 <= p :
                p2 = n.find( '/' , p + 1 )
                if p < p2 : typename = '%s %s' % ( tn , n[p:p2] )
                else      : typename = '%s %s' % ( tn , n[p:  ] )            
            else          : typename = '%s'    %   tn  

        typename = typename.replace ( 'Float_t'  , 'float'  ) 
        typename = typename.replace ( 'Double_t' , 'double' ) 
        typename = typename.replace ( 'Bool_t'   , 'bool'   )
        
        rr = [ b , typename ]
        
        stat = bbstats.get ( b , None  )
        if stat :  
            n    = stat.nEntries() 
            mnmx = stat.minmax ()
            mean = stat.mean   () 
            rms  = stat.rms    ()
            rr += [ ( '%+.5g' % mean.value() ).strip()                  , ## 2
                    ( '%.5g'  % rms          ).strip()                  , ## 3 
                    ( '%+.5g' % mnmx[0]      ).strip()                  , ## 4
                    ( '%+.5g' % mnmx[1]      ).strip()                  , ## 5
                    '' if  n == n0 else '%.3g' % ( float ( n ) / n0 ) ]   ## 6            
        else :
            ## logger.warning ("table: can't get info for the leaf \"%s\"" % b )
            rr +=  [ '-' , '-' , '-' , '-' , '' ]
            
        _vars.append ( tuple  ( rr ) )
        
    _vars.sort() 
    name_l  = len ( 'Variable' )  
    mean_l  = len ( 'mean' ) 
    rms_l   = len ( 'rms'  ) 
    min_l   = len ( 'min'  )  
    max_l   = len ( 'max'  )  
    num_l   = len ( '#'    )    
    type_l  = len ( 'type' )    
    for v in _vars :
        name_l = max ( name_l , len ( v[0] ) )
        type_l = max ( type_l , len ( v[1] ) )
        mean_l = max ( mean_l , len ( v[2] ) )
        rms_l  = max ( rms_l  , len ( v[3] ) )
        min_l  = max ( min_l  , len ( v[4] ) )
        max_l  = max ( max_l  , len ( v[5] ) )
        num_l  = max ( num_l  , len ( v[6] ) )

    
    __vars = []
    for v in _vars :
        if not ' ' in v[1]  :
            __vars.append ( v )
            continue 
        tn    = v [1]
        cl    = len ( tn )
        ml    = type_l
        vv    = list  ( v )
        vv[1] = tn.replace ( ' ' , ( ml + 1 - cl ) * ' ' , 1 ) 
        vv    = tuple ( vv ) 
        __vars.append ( vv )
        
    _vars = __vars 

    sep      = '# +%s+%s+%s+%s+%s+' % ( ( name_l       + 2 ) * '-' ,
                                        ( type_l       + 2 ) * '-' ,
                                        ( mean_l+rms_l + 5 ) * '-' ,
                                        ( min_l +max_l + 5 ) * '-' ,
                                        ( num_l        + 2 ) * '-' )
    fmt = '# | %%-%ds | %%-%ds | %%%ds / %%-%ds | %%%ds / %%-%ds | %%%ds |'  % (
        name_l ,
        type_l ,
        mean_l ,
        rms_l  ,
        min_l  ,
        max_l  ,
        num_l  )
    
    header  = fmt % ( 'Variable' ,
                      'type'     , 
                      'mean'     ,
                      'rms'      ,
                      'min'      ,
                      'max'      ,
                      '#'        )
    
    report += '\n' + sep
    report += '\n' + header
    report += '\n' + sep            
    for v in _vars :
        line    =  fmt % ( v[0] , v[1] , v[2] , v[3] , v[4] , v[5] , v[6] )
        report += '\n' + line  
    report += '\n' + sep
    
    return report, len ( sep )  


# ==============================================================================
## print rot-tree in a form of the table
#  @code
#  data = ...
#  print dat.table() 
#  @endcode
def _rt_table_ (  dataset ,  variables = [] ,   cuts = '' , *args ) :
    """print dataset in a form of the table
    >>> dataset = ...
    >>> print dataset.table()
    """
    return _rt_table_0_ ( dataset ,  variables , cuts , *args )[0]


# =============================================================================
## simplified printout for TTree/TChain
#  @see TTree
#  @code
#
#  >>> tree = ...
#  >>> print tree.pprint()
#
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-02-04
def _rt_print_ ( t ) :
    """Simplified print out for tree/chain

    >>> tree = ...
    >>> print tree.pprint() 
    """
    ##
    res = "Name: %s Enries/#%d" %  ( t.GetName() , t.GetEntries() ) 
    if hasattr ( t , 'GetNtrees' ) : res += " Chain/#%d " %       t.GetNtrees()
    ##
    _l          = t.leaves   ()
    res        += "\nLeaves:\n%s"    % multicolumn ( list( _l ) , indent = 2 , pad = 1 )

    ## collect non-trivial branches 
    _b          = t.branches ()
    
    _bs = set  ( _b )
    _ls = set  ( _l )
    _b  = list ( _bs - _ls ) 
    _b . sort() 
    if _b : res += "\nNon-trivial branches:\n%s" % multicolumn ( _b , indent = 2 ,  pad = 1 )

    return res.replace ('\n','\n# ') 


# =============================================================================
##  print DataSet
def _rt_print2_ ( data  ) :
    """Print TTree/TChain"""

    br = len ( data.branches() ) 
    l  = len ( data            )
    if 10000000 < br * l : return _rt_print_ ( data )
    
    if not isatty() : return _rt_table_ ( data )
    th  , tw   = terminal_size()
    rep , wid  = _rt_table_0_ ( data ) 
    if wid < tw  : return rep
    ##
    return _rt_print_ ( data )


ROOT.TTree.__repr__ = _rt_print2_
ROOT.TTree.__str__  = _rt_print2_
ROOT.TTree.table    = _rt_table_ 


# =============================================================================
## get list of files used for the given chain
#  @code
#  >>> chain = ... ## get the files 
#  >>> files = chain.files() 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-02-04
def _rc_files_ ( chain ) :
    """Get the list of files used for the chain
    
    >>> chain = ... ## get the files 
    >>> files = chain.files()
    """
    lst = chain.GetListOfFiles()
    return [ i.GetTitle() for i in lst ]


ROOT.TChain. files = _rc_files_

# =============================================================================
## get number of used for the given chain
#  @code
#  >>> chain = ... ## get the files 
#  >>> n = chain.nFiles() 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-02-04
def _rc_nfiles_ ( chain ) :
    """Get number of of files used for the chain
    
    >>> chain = ... ## get the files 
    >>> n = chain.nFiles()
    """
    lst = chain.GetListOfFiles()
    return lst.GetEntries()

ROOT.TChain. nFiles = _rc_nfiles_

# =============================================================================
## get the chain of reduced size (in terms of number of input files)
#  @code
#  chain = ...
#  new_chain = chain[1:3] ## keep only files 1-3
#  print len(chain), len(new_chain)
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-03-17
def _rc_getslice_ ( self , start , stop , *step ) :
    """ Get the chain of reduced size (in terms of number of input files) 
    >>> chain = ...
    >>> new_chain = chain[1:3] ## keep pnly files 1-3
    >>> print len(chain), len(new_chain)
    """
    _files = self.files()
    ## get slice 
    _files = _files[ slice(start,stop,*step) ] 
    _chain = ROOT.TChain( self.GetName() , self.GetTitle() )
    for _f in _files : _chain.Add ( _f )
    return _chain

ROOT.TChain.__getslice__ = _rc_getslice_

# =============================================================================
## Get the chain corresponding to the subset of files
#  @code
#  chain  = ...
#  chain1 = chain[2 ] ## the third file only  
#  chain2 = chain[:2] ## the first three files  
#  chain3 = chain[-1] ## the last file 
#  chain4 = chain[0:-1:2] ## every 2nd file 
#  @endcode
def _rc_getitem_ ( self , index ) :
    """ Get the chain corresponding to the subset of files
    >>> chain = ...
    >>> chain1 = chain[2 ] ## the third file only  
    >>> chain2 = chain[:2] ## the first three files  
    >>> chain3 = chain[-1] ## the last file 
    >>> chain4 = chain[0:-1:2] ## every 2nd file 
    """
    
    _files = self.files()
    
    if isinstance ( index , ( int , long ) ) :
        
        assert 0 <= index < len ( _files ), "Invalid index %s" % index
        _c = ROOT.TChain( self.GetName() , self.GetTitle() )
        _c.Add ( _files[ index] )
        return _c 
        
    if isinstance ( index , slice ) :
        _fs = _files [ index ] 
        _c  = ROOT.TChain( self.GetName() , self.GetTitle() )
        for _f in _fs :_c.Add ( _f  )
        return _c 

    raise TypeError("Invalid index type %s/%s" % ( index  , type ( index ) ) )

ROOT.TChain.__getitem__ = _rc_getitem_

# =============================================================================
## get "slice" from TTree in a form of numpy.array
#  @code
#  tree = ...
#  varr = tree.slice('Pt','eta>3')
#  print varr 
#  @endcode 
#  @see numpy.array 
#  @author Albert BURSCHE
#  @date 2015-07-08
def _rt_slice_ ( tree , varname , cut = '' ) :
    """ Get ``slice'' from TTree in a form of numpy.array
    ##
    >>> tree = ...
    >>> varr = tree.slice('Pt','eta>3')
    >>> print varr 
    """
    #
    ## decode the name (if needed)
    if isinstance ( varname , str ) :
        varname = varname.strip()
        varname = varname.replace ( ':' , ',' )
        varname = varname.replace ( ';' , ',' )
        varname = varname.replace ( ' ' , ',' )
        varname = varname.split   (       ',' )
        if 1 == len ( varname ) : varname = varname[0].strip()
        else :
            for i in range( 0 , len(varname) ) : 
                varname[i] = varname[i].strip()  
    #
    if       isinstance ( varname ,  ( list , tuple ) ) :
        ## forward to appropriate method 
        return tree.slices ( varname , cut )
    elif not isinstance ( varname , str ) :
        raise AttibuteError ( 'Invalid type %s' % varname )
    
    ##
    p1 = varname.find( '[')
    if 0 < p1 :
        p2 = varname.find( ']' , p1 + 1 )
        if p1 < p2 :
            raise AttributeError("TTree:slice: can't slice array-like variable '%s'" % varname )
            
    ge   = long( tree.GetEstimate() ) 
    tree.SetEstimate ( max ( len ( tree ) , ge ) )
    ##
    n    = tree.Draw ( varname , cut , "goff" )
    ##
    import numpy
    sl =   numpy.array ( numpy.frombuffer ( tree.GetV1() , count = n ) , copy = True )
    ##
    tree.SetEstimate ( ge ) 
    return sl 


# =============================================================================
## get "slices" from TTree in a form of numpy.array
#  @code
#  tree = ...
#  varrs1 = tree.slices ( ['Pt','eta'] , 'eta>3' )
#  print varrs1 
#  varrs2 = tree.slices (  'Pt , eta'  , 'eta>3' )
#  print varrs2
#  varrs3 = tree.slices (  'Pt : eta'  , 'eta>3' )
#  print varrs3
#  @endcode 
#  @see numpy.array 
#  @author Albert BURSCHE
#  @date 2015-07-08  
def _rt_slices_ ( tree , varnames , cut = '' ) :
    """ Get ``slices'' from TTree in a form of numpy.array
    
    >>> tree = ...
    
    >>> varrs1 = tree.slices( ['Pt' , 'eta'] ,'eta>3')
    >>> print varrs1
    
    >>> varrs2 = tree.slices( 'Pt,eta'  ,'eta>3')
    >>> print varrs2
    
    >>> varrs3 = tree.slices( 'Pt : eta' ,'eta>3')
    >>> print varrs3
    """
    #
    varname = varnames 
    ## decode the name (if needed)
    for sep in ( ',' , ':' , ';' ) :
        if isinstance ( varname , str ) :
            varname = varname.strip() 
            varname = varname.split( sep )
            if 1 == len ( varname ) : varname = varname[0].strip()
            else :
                for i in range( 0 , len(varname) ) : 
                    varname[i] = varname[i].strip()  
    #
    if       isinstance ( varname , str ) :
        ## forward to appropriate method 
        return tree.slice ( varname , cut )
    elif not isinstance ( varname ,  ( list , tuple ) ) :
        raise AttibuteError ( 'Invalid type %s' % varname )
    ##
    import numpy
    a = numpy.array ( [tree.slice(name, cut) for name in varname ] )
    a.sort()
    return a


ROOT.TTree .slice  = _rt_slice_
ROOT.TTree .slices = _rt_slices_

def _not_implemented_ ( self , method , *args , **kwargs ) :
    raise NotImplementedError('%s: the method "%s" is not implemented' % ( self.__class__ , method ) ) 

ROOT.TChain.slice  = lambda s,*x : _not_implemented_( s , 'slice'  , *x ) 
ROOT.TChain.slices = lambda s,*x : _not_implemented_( s , 'slices' , *x ) 

# =============================================================================
## extending the existing chain 
def _tc_iadd_ ( self ,  other ) :
    """ Add elements (files,  chains) to existing chain
    >>>  chain  = ...
    >>>  chain += 'myfile.root'
    >>>  chain += ( 'myfile1.root' , 'myfile2.root' )    
    """
    if   self == other                             : return self    
    elif isinstance ( o , ( list , tuple , set ) ) :        
        for f in other : _tc_iadd_ (  self , f )
        return  self
    
    elif isinstance ( other , ROOT.TChain ) :
        return _tc_iadd_ ( self , other.files() ) 
    
    elif isinstance ( other , str ) :        
        if not other in self.files () : self.Add ( other )
        return self
        
    return NotImplemented 

# =============================================================================
## summing two existing chains
def _tc_add_ ( self ,  other ) :
    """ Add two  chains together 
    >>>  chain1 = ...
    >>>  chain2 = ...
    >>>  chain3 =  chain1         + chain2
    >>>  chain4 =  chain1         + 'my_file.root'
    >>>  chain5 =  chain1         + ( 'my_file1.root' , 'my_file2.root' )
    >>>  chain5 =  'my_file.root' + chain2 
    """
    left  = ROOT.TChain ( self.GetName() )
    left += self
    left += other 
    return  left

ROOT.TChain.__iadd__ = _tc_iadd_
ROOT.TChain.__add__  = _tc_add_
ROOT.TChain.__radd__ = _tc_add_


from ostap.io.root_file import top_dir
ROOT.TTree.topdir = property ( top_dir , None , None ) 

# =============================================================================
from ostap.utils.utils import CleanUp
# =============================================================================
## @class Chain
#  simple class to keep pickable definitinon of tree/chain
#  it is needed for multiprcessing 
class Chain(CleanUp) :
    """Simple class to keep definition of tree/chain
    """
    def __getstate__  ( self ) :
        return { 'name'     : self.__name     ,
                 'files'    : self.__files    ,
                 'first'    : self.__first    ,            
                 'nevents'  : self.__nevents  }
    def __setstate__  ( self , state ) :        
        self.__name    = state [ 'name'    ]
        self.__files   = state [ 'files'   ]
        self.__first   = state [ 'first'   ]
        self.__nevents = state [ 'nevents' ]
        ## reconstruct the chain 
        self.__chain   = ROOT.TChain ( self.__name  )
        for f in self.__files  : self.__chain.Add ( f )
                    
    def __init__ ( self , tree = None ,  name = None , files = [] , first = 0 , nevents = -1  ) :

        assert ( name and files                  ) or \
               ( name and files and tree is True ) or \
               ( isinstance ( tree , ROOT.TTree  ) and valid_pointer ( tree  ) ) ,  \
               "Invalid tree/name/files combination: %s/%s%s" % ( tree , name , files    )
        assert isinstance ( first , int ) and  0 <= first     , \
               "Invalid ``first'' %s/%s"                      % ( first , type ( first ) ) 

        self.__first   = int ( first )  
        self.__nevents = nevents if 0 <= nevents < ROOT.TChain.kMaxEntries else -1 
        self.__chain   = None
        
        if files and isinstance ( files , str ) : files = files,

        if name and files :

            self.__name  = name
            self.__files = files

            if not tree is True : 
                chain = self.__create_chain() 
                assert valid_pointer ( chain ), 'Invalid TChain!'
                assert len ( files ) == len( chain.files() ) , 'Invalid length of files'
                self.__chain = chain 
            
        elif valid_pointer ( tree ) :
            
            self.__name = tree.GetName()
            
            if   isinstance ( tree ,  ROOT.TChain ) :
                
                self.__files = tuple ( tree.files() ) 
                self.__chain = tree
                
            elif isinstance ( tree ,  ROOT.TTree ) :
                
                topdir = tree.topdir
                
                if isinstance ( topdir , ROOT.TFile ) : self.__files = topdir.GetName() ,
                else :
                    
                    fname  = CleanUp.tempfile ( suffix = '.root' , prefix = 'tree_' )
                    from ostap.core.core import ROOTCWD
                    with ROOTCWD() : 
                        import ostap.io.root_file
                        with ROOT.TFile ( fname , 'NEW') as rfile :
                            rfile.cd()
                            tname = tree.GetName() 
                            rfile[ tname ] = tree 
                            rfile.ls()
                    self.__files   =   fname ,
                    self.tmpfiles += [ fname ]
                    
                chain = ROOT.TChain( tree.GetName() )
                chain.Add ( self.__files[0] )
                tmp = chain
                assert len ( chain ) == len ( tree ) , 'Something wrong happens here :-( '
                self.__chain = chain

        ## ##  last adjustment
        ## ll = len ( self )  
        ## if  ll < self.__first :
        ##     logger.warning ( 'Tree/Chain will be empty %s/%s' %   ( self.__first , ll ) )
        ##     self.__first   = ll
        ##     self.__nevents = 0
                        
        ## if number of events is specified: 
        ## if 0 < self.__nevents :
        ##    ll = len ( self )  
        ##    self.__nevents = min ( self.__nevents , ll - self.__first )

    ## split the chain for several chains  with at most chunk_size entries
    def slow_split ( self , chunk_size = 200000 ) :
        """Split the tree for several trees with chunk_size entries
        >>> tree = ....
        >>> trees = tree.split ( chunk_size = 1000000 ) 
        """

        if chunk_size <= 0 : chunk_size = ROOT.TChain.kMaxEntrie
        
        trees = []

        ievt = 0
        nevt = 0
        
        for f in self.__files :
            
            if 0 <= self.__nevents and self.__nevents <= nevt : break  ## BREAK

            ## get the length of the current tree 
            tt = Tree ( name  = self.name , file = f )
            ll = len  ( tt )
            if ievt + ll < self.__first : continue                     ## CONTINUE 
            
            ##  
            first   = self.__first - ievt if ievt <= self.__first else 0
            nevents = -1 
            if 0 <= self.__nevents and self.__nevents  < nevt + ll : 
                nevents  = self.__nevents - nevt
            t = Tree ( tt.chain , name = self.name  , file = f , first = first , nevents = nevents )
            trees += list ( t.split  ( chunk_size ) ) 
            
            ievt += ll
            nevt += nevents if 0 <= nevents else ll 

        return tuple ( trees ) 

    ## simple generator split lst into chunks 
    def get_slices ( self , first , last , chunk_size ) :
        """Split ``lst'' into  chunks
        """
        for i in range ( 0 , last - first , chunk_size ) :
            ## yield lst [ i : min ( i + chunk_size , list_size ) ] 
            yield slice ( first + i , first + min ( i + chunk_size , last -first  ) ) 
                       
    ## split the chain for several chains with at most chunk_size entries
    def split ( self , chunk_size = -1 , max_files = 10 ) :
        """Split the tree for several trees with chunk_size entries
        >>> tree = ....
        >>> trees = tree.split ( chunk_size = 1000000 ) 
        """
        if chunk_size <= 0 : chunk_size = ROOT.TChain.kMaxEntries
        if max_files  <= 0 : max_files  = 1 
        
        if 0 != self.first or 0 < self.__nevents :
            return self.slow_split ( chunk_size )
        
        ## first split on per-file basis
        fs     = self.files
        chains = [ Chain ( tree = True , name = self.name , files = fs[c] ) for c in self.get_slices ( 0 , len(fs) , max_files ) ]
        ## chains = [ Chain ( self.chain[c] ) for c in self.get_slices ( 0 , len(fs) , max_files ) ]

        if chunk_size < 0 : return tuple ( chains )
        
        result = []  
        for ch in chains :

            size = len ( ch ) 
            if   size > chunk_size and 1 == ch.nFiles : 
                tree    = Tree ( ch.chain , name = ch.name , file = ch.files[0] , first = ch.first , nevents = ch.nevents ) 
                result += tree.split ( chunk_size )
            elif size > chunk_size : result += ch.slow_split ( chunk_size )                     
            else                   : result.append ( ch ) 

        return  tuple ( result ) 


    ##  number of entries in the Tree/Chain
    def __len__ ( self ) :
        if self.__chain is None : self.__chain = self.__create_chain () 
        return len ( self.__chain )
    
    def __create_chain ( self ) :
        """``chain'' : get the underlying tree/chain"""
        c = ROOT.TChain ( self.__name )
        for f in self.__files  : c.Add ( f )
        return c

    @property
    def chain ( self ) :
        """``chain'' : get the underlying tree/chain"""
        if self.__chain is None : self.__chain = self.__create_chain () 
        return self.__chain

    @property
    def name    ( self ) :
        """``name''   : TTree/TChain name"""
        return self.__name
    @property
    def files   ( self ) :
        """``files''   : the files"""
        return self.__files
    @property
    def nFiles   ( self ) :
        """``nFiles''   : the numer of files"""
        return len(self.__files)
    @property
    def first   ( self ) :
        """``first'' : the first event to process"""
        return self.__first
    @property
    def last    ( self ) :
        """``last'' : the last event (not-inclusive)"""
        ll = len ( self ) 
        return self.__first + min ( ll , self.nevents )  
    @property
    def nevents ( self ) :
        """``nevents'' : number of events to process"""
        return self.__nevents

    ## get DataFrame 
    def  frame ( self , *vars ) :
        """``frame'' : get ROOT.RDataFrame for the given chain/tree
        >>> tree = ....
        >>> f  =  tree.frame () ## get
        >>> f1 =  tree.frame ('px', 'py' , 'pz') ## get frame with default branches
        """
        from ostap.frames.frames import DataFrame
        from ostap.core.core     import strings 
        fnames = strings  ( *self.files )
        vnames = strings  ( *vars  )
        return DataFrame ( self.name , fnames , vnames ) 
        
    def __str__ ( self ) :
        r = "Chain('%s',%s" % ( self.name , self.__files )
        if 0 != self.first or 0 <= self.__nevents : r += ",%s,%s" % ( self.first , self.__nevents )            
        return r + ")"
    __repr__ = __str__
    
# =============================================================================
## @class Tree
#  simple class to keep 'persistent' definition of the tree
#  it is needed for multiprocessing  
class Tree(Chain) :
    """Simple class to keep definition of tree/chain
    """
    def __getstate__  ( self )         : return Chain.__getstate__  ( self )   
    def __setstate__  ( self , state ) :        Chain.__setstate__  ( self , state ) 

    def __init__ ( self , tree = None ,  name = None , file = '' , first = 0 , nevents = -1 ) :

        if name and file :
            
            assert isinstance ( file , str ), '"File should be single file name!"'
            
        elif valid_pointer  ( tree )  :
            
            if isinstance ( tree , ROOT.TChain ) :
                assert 1 == len (  tree.files() ) , 'Tree is for ROOT.TTree only!'
                
        Chain.__init__ ( self , tree , name  , files = [ file ] , first = first , nevents = nevents )
        
        assert 1 == self.nFiles , 'Invalid number of files!'

    ## split the tree for several trees with max=chunk_size entries
    def split ( self , chunk_size = 200000  ) :
        """Split the tree for several trees with max=chunk_size entries
        >>> tree = ....
        >>> trees = tree.split ( chunk_size = 1000000 ) 
        """
        
        if 0 == self.nevents : return ()
        
        assert isinstance  ( chunk_size , int ) , "Illegal type of ``chunk_size'' %s" % chunk_size
        
        ## no splitting ?
        if 0 >= chunk_size : return  self,
        
        ll   = len ( self )
        last = min ( ll , self.first + self.nevents if 0 <= self.nevents else ROOT.TChain.kMaxEntries ) 

        result = [] 
        for s in self.get_slices ( self.first , last , chunk_size ) :
            start , stop , stride = s.indices ( ll )
            if start < stop : 
                t = Tree ( tree = self.chain , name = self.name , file = self.file , first = start , nevents = stop - start ) 
                result.append ( t ) 
                
        return tuple ( result ) 

    ## get a slice for the tree 
    def __getslice__ ( self , start , stop ) :
        """ Get a slice for the given tree 
        >>> tree  = ...
        >>> tree1 = tree[:1000]  ## the first 1000 events 
        """
        s  = slice ( start , stop )
        ll = len   ( self  )
        
        if self.first < ll : start , stop = 0 , 0

        last = self.first + self.nevents if 0 < self.nevents else ll
        last = min ( last , ll )
        
        start, stop , stride = s.indices ( last - self.first )
        
        start += self.first
        stop  += self.first
        
        return Tree ( name = self.name , file = self.file , first = start , nevents = stop )
        
    @property
    def file ( self ) :
        """``file''   : the file name """
        fs = self.files 
        assert 1 == len  ( fs ) , 'Invalid number of files %s' % len ( fs ) 
        return fs [0] 

    def __str__ ( self ) :
        r = "Tree('%s','%s'" % ( self.name , self.file )
        if 0 != self.first or 0 <= self.nevents : r += ",%s,%s" % ( self.first , self.nevents )     
        return r + ")"
    __repr__ = __str__

# =============================================================================
## Get the effective entries in data frame
#  @code
#  data = ...
#  neff = data.nEff('b1*b1')
#  @endcode
def _rt_nEff_  ( self , cuts = '' , *args ) :
    """Get the effective entries in data frame 
    >>> data = ...
    >>> neff = data.nEff('b1*b1')
    """
    return Ostap.StatVar.nEff ( self , cuts , *args )

ROOT.TTree.nEff = _rt_nEff_ 
# =============================================================================

from  ostap.stats.statvars import data_decorate as _dd
_dd ( ROOT.TTree )

# =============================================================================
_decorated_classes_ = (
    ROOT.TTree   ,
    ROOT.TChain     
    )
_new_methods_       = (
    #
    ROOT.TTree .withCuts  ,
    ROOT.TChain.withCuts  ,
    ROOT.TTree. __len__   ,
    #
    ROOT.TTree .__call__  ,
    ROOT.TChain.__call__  ,
    #
    ROOT.TTree .project   ,
    ROOT.TChain.project   ,
    #
    ROOT.TTree .statVar   ,
    ROOT.TChain.statVar   ,
    ROOT.TTree .statCov   ,
    ROOT.TChain.statCov   ,
    ROOT.TTree .statCovs  ,
    ROOT.TChain.statCovs  ,
    #
    ROOT.TTree .vminmax   ,
    ROOT.TChain.vminmax   ,
    #
    ROOT.TTree .sumVar_   ,
    ROOT.TChain.sumVar_   ,
    #
    ROOT.TTree .sumVar    ,
    ROOT.TChain.sumVar    ,
    #
    ROOT.TTree .branches  , 
    ROOT.TTree .__repr__  , 
    ROOT.TTree .__str__   ,
    #
    ROOT.TChain.files        ,
    ROOT.TChain.__getslice__ ,
    #
    ROOT.TTree.slice        ,
    ROOT.TTree.slices       ,
    #
    ROOT.TTree.nEff             , 
    ROOT.TTree.get_moment       , 
    ROOT.TTree.central_moment   , 
    ROOT.TTree.mean             ,
    ROOT.TTree.rms              ,
    ROOT.TTree.skewness         ,
    ROOT.TTree.kurtosis         ,
    ROOT.TTree.quantile         ,
    ROOT.TTree.median           ,
    ROOT.TTree.quantiles        ,
    ROOT.TTree.interval         ,
    ROOT.TTree.terciles         ,
    ROOT.TTree.quartiles        ,
    ROOT.TTree.quintiles        ,
    ROOT.TTree.deciles          ,

    )
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
