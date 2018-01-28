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
__all__     = () 
# =============================================================================
import ROOT
from ostap.core.core import cpp, VE
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.trees.trees' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Some useful decorations for Tree/Chain objects')
# =============================================================================
import ostap.trees.cuts
# =============================================================================
_large = 2**64
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
    
    pit = cpp.Ostap.PyIterator ( self , cuts , first , last )
    if not pit.ok() : raise TypeError ( "Invalid Formula: %s" % cuts )
    #
    from ostap.utils.progress_bar import ProgressBar 
    with ProgressBar ( min_value = first        ,
                       max_value = last         ,
                       silent    = not progress ) as bar :
        
        step = 13.0 * max ( bar.width , 101 ) / ( last - first ) 
        
        _t = pit.tree()
        _o = _t 
        while _t :

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
def _tc_call_ ( self , first = 0 , last = _large , cuts = None , progress = False ) :
    """Iterator over ``good events'' in TTree/TChain:
    
    >>> tree = ... # get the tree
    >>> for i in tree(0, 100 , 'pt>5' ) : print i.y
    
    """
    #
    last = min ( last , len ( self )  )

    from ostap.utils.progress_bar import ProgressBar 
    with ProgressBar ( min_value = first        ,
                       max_value = last         ,
                       silent    = not progress ) as bar :
        
        step = 13.0 * max ( bar.width , 101 ) / ( last - first ) 

        pit = 1 
        if cuts :
            
            pit = cpp.Ostap.PyIterator ( self , cuts , first , last )
            if not pit.ok() : raise TypeError ( "Invalid Formula: %s" % cuts )
            #
            
            _t = pit.tree()
            _o = _t 
            while _t :
                
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
def _tt_project_ ( tree , histo , what , cuts = '' , *args ) :
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

    ## check for comma-separated list of expressions:
    if isinstance ( what , str ) :
        what = what.split(',')
        if 1 == len(what) : what = what[0]

    ## check for semicolumn-separated list of expressions:
    if isinstance ( what , str ) :
        what = what.split(';')
        if 1 == len(what) : what = what[0]

    #
    if   isinstance ( what  , str       ) : what =     what 
    elif isinstance ( what  , ROOT.TCut ) : what = str(what)  
    elif isinstance ( histo , ROOT.TH1  ) : 
        rr = 0 
        hh = histo.clone()
        for v in what :
            r , h  = _tt_project_ ( tree , hh , v , cuts , *args )
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
        result = tree.Project ( hname , what , cuts , *args )
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
    return cpp.Ostap.StatVar.statVar ( tree , expression , *cuts )

ROOT.TTree     . statVar = _stat_var_
ROOT.TChain    . statVar = _stat_var_

# =============================================================================
## get the statistic for pair of expressions in Tree/Dataset
#  @code
#  tree  = ...
#  stat1 , stat2 , cov2 , len = tree.statCov( 'x' , 'y' )
#  # apply some cuts 
#  stat1 , stat2 , cov2 , len = tree.statCov( 'x' , 'y' , 'z>0' )
#  # use only subset of events
#  stat1 , stat2 , cov2 , len = tree.statCov( 'x' , 'y' , 'z>0' , 100 , 10000 )
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
    stat1  = cpp.Ostap.WStatEntity       ()
    stat2  = cpp.Ostap.WStatEntity       ()
    cov2   = cpp.Ostap.Math.SymMatrix2x2 ()

    if cuts : 
        length = cpp.Ostap.StatVar.statCov ( tree        ,
                                             expression1 ,
                                             expression2 ,
                                             cuts        ,
                                             stat1       ,
                                             stat2       ,
                                             cov2        , 
                                             *args       )
    else :
        length = cpp.Ostap.StatVar.statCov ( tree        ,
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
    
    _SV    = cpp.std.vector('std::string')
    _vars  = _SV()
    vars   = expressions
    for e in vars : _vars.push_back( e )
    
    WSE    = cpp.Ostap.WStatEntity
    _WV    = cpp.std.vector( WSE )
    _stats = _WV()
    _DV    = cpp.std.vector('double')
    _cov2  = _DV()

    if cuts : 
        length = cpp.Ostap.StatVar._statCov ( tree   ,
                                              _vars  ,
                                              cuts   ,
                                              _stats ,
                                              _cov2  ,
                                              *args  ) 
    else :
        length = cpp.Ostap.StatVar._statCov ( tree   ,
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
    COV2 = cpp.Ostap.Math.SymMatrix ( l )
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
    if cuts : s = tree.statVar ( var , cuts )
    else    : s = tree.statVar ( var )
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
_h_one_ = ROOT.TH1D( "_h_one_", '' , 3 , -1 , 2 ) ; _h_one_.Sumw2()
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
    w = tree.statVar ( expression , *cuts )
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
def _rt_leaves_ ( t ) :
    """ Get the list of leaves names        
    
    >>> tree = ...
    >>> lst = tree.leaves()
    >>> for l in lst : print l
    """
    _lst =  t.GetListOfLeaves()
    if not _lst : return tuple() 
    _lst = [ l.GetName() for l in _lst ] 
    _lst.sort()
    return tuple( _lst ) 

ROOT.TTree.leaves   = _rt_leaves_

# =============================================================================
## get the branches for the given tree/chain
#  @see TTree
#  @code
#
#  >>> tree = ...
#  >>> lst = tree.branches()
#  >>> for b in lst : print b
#
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-02-04
def _rt_branches_ ( t ) :
    """Get the list of branch names
    
    >>> tree = ...
    >>> lst = tree.branches()
    >>> for b in lst : print b
    
    """
    _lst =  t.GetListOfBranches()
    if not _lst : return tuple() 
    _lst = [ l.GetName() for l in _lst ] 
    _lst.sort()
    return tuple( _lst ) 

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
    _b          = t.branches ()
    res        +=        "\nBranches: %s" % list( _b )
    _l          = t.leaves   ()
    if _l != _b : res += "\nLeaves: %s"   % list( _l )
    return res

ROOT.TTree.__repr__ = _rt_print_
ROOT.TTree.__str__  = _rt_print_

# =============================================================================
## get lst of files used for the given chain
#  @code
#
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
    )
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
