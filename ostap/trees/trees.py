#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/trees/trees.py
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
    'Chain'           , ## helper class , needed for multiprocessing 
    'Tree'            , ## helper class , needed for multiprocessing
    'ActiveBranches'  , ## context manager to activate certain branches 
    'active_branches' , ## context manager to activate certain branches
  ) 
# =============================================================================
from   ostap.core.meta_info      import root_info
from   ostap.core.core           import ( std , Ostap , VE  , WSE , hID ,
                                          ROOTCWD , strings  , split_string ,
                                          valid_pointer           ) 
from   ostap.core.ostap_types    import ( integer_types  , long_type      ,
                                          string_types   , sequence_types ,
                                          sized_types    , num_types      ,
                                          dictlike_types , list_types     )
from   ostap.utils.utils         import chunked
from   ostap.utils.basic         import isatty, terminal_size, NoContext 
from   ostap.utils.scp_copy      import scp_copy
from   ostap.math.reduce         import root_factory
from   ostap.utils.progress_bar  import progress_bar
import ostap.trees.treereduce 
import ostap.histos.histos
import ostap.trees.param
import ROOT, os, math, array 
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
_large = ROOT.TVirtualTreePlayer.kMaxEntries
# =============================================================================


# =============================================================================
## check validity/emptiness  of TTree/TChain
#  require non-zero poniter and non-empty Tree/Chain
def _tt_nonzero_ ( tree ) :
    """Check validity/emptiness  of TTree/TChain
    - require non-zero poniter and non-empty Tree/Chain
    """
    return valid_pointer ( tree ) and 0 < tree.GetEntries() 

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
#  If only (small) fraction of branches is used in <code>cuts</code> and/or
#  only small fraction of branches wil lbe used in the loop,
#  the processing can be speed up siginificantly
#  by specification of "active" branches:
#  @code
#  tree = ...
#  sum_y = 0 
#  for i in tree.withCuts ( 'pt>10' , active = ( 'pt' , 'y') ) : sum_y += i.y 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-06
def _iter_cuts_ ( tree , cuts = '' , first = 0 , last = _large , progress = False , active = () ) :
    """Iterator over ``good events'' in TTree/TChain:
    
    >>> tree = ... # get the tree
    >>> for i in tree.withCuts ( 'pt>5' ) : print i.y
    
    Attention: TTree::GetEntry is already invoked for accepted events,
    no need in second call

    - If only (small) fraction of branches is used in 'cuts' and/or
    only small fraction of branches wil lbe used in the loop,
    the processing can be speed up siginificantly
    by specification of ``active'' branches:

    >>> tree = ...
    >>> sum_y = 0 
    >>> for i in tree.withCuts ( 'pt>10' , active = ( 'pt' , 'y') ) : sum_y += i.y 
    """
    #
    last  = min ( last , len ( tree ) )
    first = max ( 0    , first        ) 
    #
    if not cuts : cuts = '1'
    #
    if active :
        
        cvar    = tree.the_variables ( cuts , *active )
        abrs    = tuple ( set ( [ c for c in cvar ] ) ) 
        context = ActiveBranches  ( tree , *abrs )
        
    else : context = NoContext () 

    from   ostap.utils.progress_conf import progress_conf
    with context :
        
        if progress : pit = Ostap.PyIterator ( tree , progress_conf , cuts , first , last )
        else        : pit = Ostap.PyIterator ( tree ,                 cuts , first , last )
        ##
        if not pit.ok() : raise TypeError ( "Invalid Formula: %s" % cuts )
        
        ## the tree is advanced 
        mytree = pit.tree()
        while valid_pointer ( mytree ) :
            yield mytree
            mytree = pit.next()
            
        del pit

    ## set the tree at the initial position 
    ievt = tree.GetEntryNumber ( 0 )
    if 0 <= ievt : tree.GetEntry ( ievt )

    
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
def _tc_call_ ( tree , first = 0 , last = -1  , cuts = None , progress = False , active = () ) :
    """Iterator over ``good events'' in TTree/TChain:
    
    >>> tree = ... # get the tree
    >>> for i in tree(0, 100 , 'pt>5' ) : print i.y
    
    """
    #
    if last < 0 : last = _large
    
    last  = min ( last , len ( tree ) )
    first = max ( 0    , first        ) 
    
    if active :
        
        cvar    = tree.the_variables ( cuts , *active )
        abrs    = tuple ( set ( [ c for c in cvar ] ) ) 
        context = ActiveBranches  ( tree , *abrs )
        
    else :
        
        from ostap.utils.basic import NoContext 
        context = NoContext () 

    from   ostap.utils.progress_conf import progress_conf
    with context : 
        
        if cuts : ## use Ostap.PyIterator 
            
            if progress : pit = Ostap.PyIterator ( tree , progress_conf , cuts , first , last )
            else        : pit = Ostap.PyIterator ( tree ,                 cuts , first , last )
            
            if not pit.ok() : raise TypeError ( "Invalid Formula: %s" % cuts )
            
            mytree = pit.tree()
            while valid_pointer ( mytree ) :                
                yield mytree                    ## YIELD HERE                  
                mytree = pit.next()             ## advance to the next entry  
                
            del pit
            
        else : ## trivial loop 

            ## explicit loop over entries 
            for entry in progress_bar ( range ( firts , last ) , silent = not progress )  :
                
                ievt = tree.GetEntryNumber ( entry  )
                
                if ievt < 0 :
                    logger.error ( "Invald entry/1 %s< skip it!" % entry )
                    continue
                
                ievt = tree.GetEntry ( ievt )            
                if ievt < 0 :
                    logger.error ( "Invald entry/2 %s< skip it!" % entry )
                    continue
                
                yield tree                     ## YIELD HERE 

            
    ## set the tree at the initial position 
    ievt = tree.GetEntryNumber ( 0 )
    if 0 <= ievt : tree.GetEntry ( ievt )
    

ROOT.TTree .__call__  = _tc_call_ 
ROOT.TChain.__call__  = _tc_call_

try : 
    from numpy import frombuffer as _frombuffer 
    def get_result ( vct ) :
        return _frombuffer ( vct.data() , count = len ( vct ) , dtype = float )
        ## return _array ( data , dtype = float )
except ImportError :
    from array import array as _array 
    def get_result ( vct ) :
        return _array ( 'd' , vct )

# =============================================================================
##  Iterate over tree entries and get a row/array of values for each good entry
#   @code
#   tree = ...
#   for row in tree.rows ( 'a a+b/c sin(d)' , 'd>0' ) :
#      print ( row ) 
#   @code 
def _tt_rows_ ( tree , variables , cuts = '' , first = 0 , last = -1 , progress = False , active = () ) :
    """Iterate over tree entries and get a row/array of values for each good entry
    >>> tree = ...
    >>> for row in tree.rows ( 'a a+b/c sin(d)' , 'd>0' ) :
    >>>    print ( row ) 
    """
    
    if last < 0 : last = _large
    last  = min ( last , len ( tree ) )
    first = max ( 0    , first        ) 
    
    if isinstance ( variables , string_types ) : variables = split_string ( variables , ' ,;:' )
    vars = []
    for v in variables :
        vars += split_string ( v , ' ,;:' )

    if active :
        
        cvar    = tree.the_variables ( cuts , vars ,  *active )
        abrs    = tuple ( set ( [ c for c in cvar ] ) ) 
        context = ActiveBranches  ( tree , *abrs )
        
    else :
        
        from ostap.utils.basic import NoContext 
        context = NoContext () 
        
    vars = strings ( vars ) 

    from   ostap.utils.progress_conf import progress_conf
    with context :
        
        getter = Ostap.Trees.Getter ( tree , vars ) 
        result = std.vector('double')()
        
        if cuts : ## more efficient loop using Ostap.PyIterator 
            
            if progress : pit = Ostap.PyIterator ( self , progress_conf , cuts , first , last )
            else        : pit = Ostap.PyIterator ( self ,                 cuts , first , last )
            
            assert pit and pit.ok() , 'ROWS: Invalid formula %s' % cuts  
            
            mytree = pit.tree() 
            while valid_pointer ( mytree ) :
                
                sc = getter.eval ( result )
                if sc.isFailure () :
                    logger.error('ROWS: Error status from getter %s' % sc  ) 
                    break
                
                yield get_result ( result )                 
                mytree = pit.next()
                
            del pit
        
        else : ## trivial loop
        
            for event in progress_bar ( range ( first , last ) , silent = not progress ) :
                
                tt      = getter.tree()
                
                ievent  = tt.GetEntryNumber ( event  )
                if ievent < 0 :
                    logger.error('ROWS: Cannot read entry %s' % event ) 
                    break
                
                ientry  = tt.GetEntry ( ievent )
                if ientry < 0 :
                    logger.error('ROWS: Cannot get entry  %s' % event ) 
                    break
                
                sc = getter.eval ( result )
                if sc.isFailure () :
                    logger.error('ROWS: Error status from getter %s' % sc  ) 
                    break 
                
                yield get_result ( result )                 

    del getter
    
    ## set the tree at the initial position 
    ievt = tree.GetEntryNumber ( 0 )
    if 0 <= ievt : tree.GetEntry ( ievt )

            
ROOT.TTree .rows  = _tt_rows_ 

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
#  @endcode
#  @attention For 2D&3D cases if variables specifed as singel string, the order is Z,Y,X,
#             Otherwide the natural order is used.
#  @param tree       (INPUT) the tree
#  @param histo      (INPUT/UPDATE) the histogram or histogram name 
#  @param what       (INPUT) variable/expression to be projected.
#                            It could be a list/tuple of variables/expressions
#                            or just a comma or semicolumn-separated expression
#  @param cuts       (INPUT) expression for cuts/weights
#  @param options    (INPUT) options to be propagated to <code>TTree.Project</code>
#  @param first      (INPUT) first entry to process
#  @param last       (INPUT) last entry to process
#  @param use_frame  (INPUT) use DataFrame for processing?
#  @param silent     (INPUT) silent processing?
#  @see TTree::Project
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-07-06
def tree_project ( tree               ,
                   histo              ,
                   what               ,
                   cuts       = ''    ,
                   first      =  0    , 
                   last       = -1    , 
                   use_frame  = True  , ## use DataFrame ? 
                   progress   = False ) :
    """Helper project method
    
    >>> tree = ...
    
    >>> h1   = ROOT.TH1D(... )
    >>> tree.Project ( h1.GetName() , 'm', 'chi2<10' ) ## standart ROOT 
    
    >>> h1   = ROOT.TH1D(... )
    >>> tree.project ( h1.GetName() , 'm', 'chi2<10' ) ## ditto 
    
    >>> h1   = ROOT.TH1D(... )
    >>> tree.project ( h1           ,  'm', 'chi2<10' ) ## use histo

    - histo : the histogram (or histogram name)
    - what  : variable/expression to project. It can be expression or list/tuple of expression or comma (or semicolumn) separated expression
    - cuts  : selection criteria/weights
    
    - For 2D&3D cases if variables specifed as singel string, the order is Z,Y,X,
    - Otherwide the natural order is used.

    """
    #
    assert isinstance ( histo , ROOT.TH1 ) , "Invalid type of 'histo' %s" % type ( histo )
    histo.Reset()
    
    assert not cuts or \
           isinstance ( cuts , string_types    ) or \
           isinstance ( cuts , ROOT.TCut       ) ,  \
           "Invalid 'cuts' %s" % type ( cuts ) 
    
    first   = max ( first , 0 )

    ## chain helper object?
    if   isinstance ( tree , Chain           ) : tree = tree.chain 
    elif isinstance ( tree , ROOT.RooAbsData ) :
        from ostap.fitting.dataset import ds_project as _ds_project_
        result = _ds_project_ ( tree , histo , what , cuts , first , last , progress )
        return result.nEntries () , result 
    
    nevents = len ( tree )
    last    = nevents if last < 0 else min ( nevents , last )

    if isinstance ( cuts , ROOT.TCut ) : cuts = str ( cuts )
    elif not cuts                      : cuts = ''
    
    tail = first , last , use_frame , progress
    
    if isinstance ( what , string_types ) :
        
        ## ATTENTION reverse here! 
        what = tuple ( reversed ( [ v.strip() for v in split_string ( what , ':;,' ) ] ) ) 
        return tree_project ( tree , histo , what , cuts , *tail )

    assert isinstance ( what , list_types ) , "tree_project: invalid 'what' %s" % what 

    hdim = histo.dim()
    assert len ( what ) == hdim  or ( 1 == hdim and hdim < len ( what ) ) , \
           "tree_project: invalid 'what' %s" % what 
    
    ## special treatment for 1D histograms: several variables are summed/projected togather 
    if 1 == hdim and hdim < len ( what ) :
        htmp  = histo.clone()
        total = 0 
        for w in what :
            n   , htmp = tree_project ( dataset , htmp , w , cuts , *tail )
            total += n
            histo += htmp 
        del htmp
        return total , histo 

    ## use frame if requested and if/when possible 
    if use_frame and isinstance ( tree , ROOT.TTree ) and 0 == first and len ( tree ) <= last :
        import ostap.frames.frames as F 
        frame   = F.DataFrame ( tree )
        if progress : counter = F.frame_progress ( frame , len ( tree ) )
        else        : counter = frame.Count()
        fargs   = what + ( cuts , ) 
        hh      = F.frame_project ( frame , histo , *fargs )
        return counter.GetValue () , hh 

    ## use 'old' method since it is slightly more efficient 
    if not progress and isinstance ( tree , ROOT.TTree ) : 
        hname = histo.GetName() 
        what  = ' : '.join ( reversed ( what ) )
        num   = tree.Project ( hname , what , '' , '' , last - first , first )
        return num , histo 
        
    active = tree.the_variables ( cuts , *what )
    from   ostap.utils.progress_conf import progress_conf
    with ActiveBranches  ( tree , *active ) :    
        
        args = ( histo ,) + what + ( cuts , first , last ) 
        if   3 == hdim :
            if progress : sc = Ostap.HistoProject.project3 ( tree , progress_conf , *args )
            else        : sc = Ostap.HistoProject.project3 ( tree ,                 *args )
            if not sc.isSuccess () :
                logger.error ( "Error from Ostap.HistoProject.project3 %s" % sc )
                return 0 , None
            return histo.nEntrie () , histo
        elif 2 == hdim :
            if progress : sc = Ostap.HistoProject.project2 ( tree , progress_conf , *args )
            else        : sc = Ostap.HistoProject.project2 ( tree ,                 *args )
            if not sc.isSuccess() :
                logger.error ( "Error from Ostap.HistoProject.project2 %s" % sc )
                return 0 , None
            return histo.nEntries () , histo 
        elif 1 == hdim :
            if progress : sc = Ostap.HistoProject.project  ( tree , progress_conf , *args )
            else        : sc = Ostap.HistoProject.project  ( tree ,                 *args )
            if not sc.isSuccess() :
                logger.error ( "Error from Ostap.HistoProject.project1 %s" % sc )
                return 0 , None
            return histo.nEntries() , histo
        
    raise TypeError ( 'tree_project, invalid case' )

ROOT.TTree .project = tree_project
ROOT.TChain.project = tree_project


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
#  @endcode
#
#  @param tree       (INPUT) the tree
#  @param histo      (INPUT/UPDATE) the histogram or histogram name 
#  @param what       (INPUT) variable/expression to be projected.
#                            It could be a list/tuple of variables/expressions
#                            or just a comma or semicolumn-separated expression
#  @param cuts       (INPUT) expression for cuts/weights
#  @param options    (INPUT) options to be propagated to <code>TTree.Project</code>
#  @param nentries   (INPUT) number of entries to process
#  @param firstentry (INPUT) first entry to process
#  @param use_frame  (INPUT) use DataFrame for processing?
#  @param silent     (INPUT) silent processing?
#  @see TTree::Project
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-07-06
def tree_project_old ( tree               ,
                       histo              ,
                       what               ,
                       cuts       = ''    ,
                       options    = ''    ,
                       nentries   = -1    ,
                       firstentry =  0    ,
                       use_frame  = False , ## use DataFrame ? 
                       silent     = False ) :
    """Helper project method
    
    >>> tree = ...
    
    >>> h1   = ROOT.TH1D(... )
    >>> tree.Project ( h1.GetName() , 'm', 'chi2<10' ) ## standart ROOT 
    
    >>> h1   = ROOT.TH1D(... )
    >>> tree.project ( h1.GetName() , 'm', 'chi2<10' ) ## ditto 
    
    >>> h1   = ROOT.TH1D(... )
    >>> tree.project ( h1           ,  'm', 'chi2<10' ) ## use histo

    - histo : the histogram (or histogram name)
    - what  : variable/expression to project. It can be expression or list/tuple of expression or comma (or semicolumn) separated expression
    - cuts  : selection criteria/weights 
    """
    #

    if nentries < 0 : nentries = _large

    args = () 
    if options or 0 < firstentry or 0 < nentries < _large :
        args = options , nentries , firstentry


    if   isinstance ( histo , ROOT.TH1     ) : hname = histo.GetName()
    elif isinstance ( histo , string_types ) :
        hname = histo
        groot = ROOT.ROOT.GetROOT()
        h     = groot.FindObject ( hname )
        if h and isinstance ( h , ROOT.TH1 ) : histo = h
        else                                 : histo = None
        
    assert histo and isinstance ( histo , ROOT.TH1 ) ,\
           "project: invalid type of ``histo'': %s " % type ( histo )

    ## reset the histogram 
    histo.Reset() 

    if isinstance ( cuts , ROOT.TCut ) : cuts = str ( cuts )

    ## comma, column or semicolumn separated list
    if isinstance ( what , string_types ) :
        ## attention! note reversed here! 
        what = [ w.strip() for w in reversed ( split_string ( what , ',;:' ) ) ] ## attention! note reversed here! 
        
    assert isinstance ( what , sequence_types  ) and \
           isinstance ( what , sized_types     ) and 1 <= len ( what ) ,\
           "project: invalid ``what''/1: %s" % str ( what )
    
    what = [ w for w in what ]    
    assert all ( isinstance ( w , string_types ) for w in what ) , \
           "project: invalid ``what''/2: %s" % str ( what )
    what = [ w.strip()  for w in what ]    


    ## special treatment for 1D histograms 
    if 1 == histo.dim() and 1 < len ( what ) :
        nr = 0
        hh = histo.clone() 
        for i , w in enumerate ( what ) : 
            n , h = tree_project_old ( tree , hh , w , cuts , *args , use_frame = use_frame , silent = silent )
            histo += h
            nr    += n 
        del hh
        return nr , histo 

    ## number of variables must match the dimensionality of the histogram 
    assert len ( what ) == histo.dim(), \
           "project: dimension mismatch : ``what''/%d vs ``dim''/%d " % ( len ( what ) , histo.dim() ) 

    ## use frame if requested and if possible 
    if use_frame and isinstance ( tree , ROOT.TTree ) and not args :
        
        from ostap.frames.frames import DataFrame, frame_project
        frame   = DataFrame ( tree , enable = False )
        counter = frame.Filter( cuts ).Count()  if cuts else frame.Count() 
        hh      = frame_project ( frame , histo , what , cuts )
        return counter.GetValue () , hh 

    ## the basic case 
    with ROOTCWD() :
        
        groot = ROOT.ROOT.GetROOT() 
        groot.cd ()

        ## attention! note reversed here! 
        vars = ' : '.join ( ( '(%s)' % w for w in reversed ( what ) ) )  ## attention! note reversed here! 

        ## make temporary histogram
        uid   = vars , cuts , histo.GetName() , tree.GetName() , args , use_frame , silent 
        
        htemp = histo.clone ( prefix = 'htmp_%X_' % ( hash ( uid ) % 2**32 ) )
        hname = htemp.GetName  () 

        ## make projection to temporary histogram 
        result = tree.Project ( hname , vars  , cuts , *args )

        if   0 == result :
            
            del htemp
            return result , histo
        
        elif 0 < result and result == htemp.GetEntries() :
            
            histo += htemp
            
            del htemp 
            return result , histo 
        
        ## make a try extract the temporary histogram from ROOT memory 
        hroot = groot.FindObject ( hname )
        if hroot and ( hroot is htemp ) :
            
            histo += hroot 
            
            del hroot
            del htemp
            
            return result, histo
        
        else :
            
            logger.error( "project: cannot  access temporary histo")
            
            histo += htemp
            
            del hroot
            del htemp 
            
    return result, histo

ROOT.TTree .project_old = tree_project_old
ROOT.TChain.project_old = tree_project_old 


# =============================================================================
## check if object is in tree/chain  :
#  @code
#  tree = ...
#  if obj in tree : 
#  ...
#  @endcode
#  Operation is defiend:
#  - integer value is "in tree" if it corresponds to the valid entry number
#  - string  value is "in tree" if it corresponds to the name of branch or leaf 
def _rt_contains_ ( tree , obj ) :
    """Check if object is in tree/chain  :
    >>> tree = ...
    >>> if obj in tree : 
    ...
    Operation is defiend:
    - integer value is ``in tree'' if it corresponds to the valid entry number
    - string  value is ``in tree'' if it corresponds to the name of branch or leaf 
    """
    
    if   isinstance ( obj , integer_types ) :        
        return 0 <= obj < len ( tree )
    
    elif isinstance ( obj , string_types  ) :
        return ( obj in tree.branches () ) or ( obj in tree.leaves () )
    
    return False 

ROOT.TTree .__contains__ = _rt_contains_
ROOT.TChain.__contains__ = _rt_contains_

# =============================================================================
## get the statistic for certain expression(s) in Tree/Dataset
#  @code
#  tree  = ... 
#  stat1 = tree.statVar ( 'S_sw/effic' )
#  stat2 = tree.statVar ( 'S_sw/effic' , 'pt>1000' )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-09-15
def _stat_var_ ( tree , expression , *cuts ) :
    """Get a statistic for the  expression in Tree/Dataset
    
    >>> tree  = ... 
    >>> stat1 = tree.statVar ( 'S_sw/effic' )
    >>> stat2 = tree.statVar ( 'S_sw/effic' ,'pt>1000')
    
    """
    
    if isinstance ( expression , string_types ) :
        
        explist = split_string ( expression , ' ,;:' )         
        if 1 != len ( explist ) :
            return _stat_vars_ ( tree , explist , *cuts )  ## RETURRN
        
    else :
        
        return _stat_vars_ ( tree ,  expression , *cuts ) ## RETURN 
    
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

    if isinstance ( expressions , string_types ) :
        return _stat_var_ ( tree , expressions , *cuts ) 
    
    if not expressions : return {}
    
    vct = strings ( *expressions )
    res = std.vector(WSE)() 

    ll  = Ostap.StatVar.statVars ( tree , res , vct , *cuts )

    assert res.size() == vct.size(), 'stat_vars: Invalid size of structures!'

    N = res.size()
    results = {} 

    for i in range ( N ) :
        results[ vct [ i ] ] = WSE ( res[i] ) 

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
    cov2   = Ostap.Math.SymMatrix(2) ()

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
        length = Ostap.StatVar.statCov ( tree   ,
                                         _vars  ,
                                         cuts   ,
                                         _stats ,
                                         _cov2  ,
                                         *args  ) 
    else :
        length = Ostap.StatVar.statCov ( tree   ,
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

    if not pattern :

        lst  = [ v.GetName() for v in vlst  ]
        lst.sort ()
        return tuple ( lst ) 

    if isinstance ( pattern , string_types ) : pattern  = [ pattern ]

    lst = set()
    for p in pattern : 
        try : 
            import re
            c    =  re.compile ( p , *args )
            vars = [ v.GetName() for v in vlst if c.match ( v.GetName () ) ]
            lst  = lst | set ( vars  ) 
        except :
            logger.error ('leaves("%s"): exception is caught, use all ' % p  , exc_info = True ) 
            lst  = lst | set ( [ v.GetName() for v in vlst  ] )
            
    lst = list ( lst )
    lst.sort () 
    return tuple ( lst ) 

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

    from   ostap.logger.utils     import multicolumn

    res = "Name: %s Entries/#%d" %  ( t.GetName() , t.GetEntries() ) 
    if hasattr ( t , 'GetNtrees' ) : res += " Chain/#%d " %       t.GetNtrees()

    ##
    _l          = list ( set ( t.leaves () ) ) 
    _l . sort ()
    _lt = [ "%s:%s" % ( l , t.leaf(l).get_type_short().replace (' [','[' ) ) for l in _l ]
    res        += "\nLeaves:\n%s"    % multicolumn ( _lt , indent = 2 , pad = 1 )

    ## collect non-trivial branches 
    _b          = t.branches ()

    _bs = set  ( _b )
    _ls = set  ( _l )
    _b  = list ( _bs - _ls ) 
    _b . sort () 
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


def _in_types ( t ) :
    while 0 <= t.find ( 2 * ' ' ) : t = t.replace ( 2 * ' ' , ' ' )
    return t in __types 


# ==============================================================================
## print tree as table 
def _rt_table_0_ ( tree , pattern = None , cuts = '' , prefix = '' , title = '' , *args ) :
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
        if not _in_types ( tn ) : continue
        
        bbs.append ( b ) 

    if hasattr ( tree , 'pstatVar' ) : bbstats = tree.pstatVar ( bbs , cuts , *args )
    else                             : bbstats = tree. statVar ( bbs , cuts , *args )

    from ostap.stats.counters import WSE 
    if isinstance ( bbstats , WSE )  : bbstats = { bbs[0] : bbstats } 
    
    for b in brs :
        
        l = tree.leaf ( b )

        if not l :
            logger.warning ("table: can't get the leaf  \"%s\"" % b )
            continue
        
        typename = l.get_type()
        
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

    index_l =   int ( math.ceil ( math.log10( len ( _vars ) + 1 ) ) )
        
    fmt_name  = '%%%ds. %%-%ds' % ( index_l , name_l )
    fmt_type  = '%%-%ds'        % type_l
    fmt_mean  = '%%%ds'         % mean_l
    fmt_rms   = '%%-%ds'        % rms_l
    fmt_min   = '%%%ds'         % mean_l
    fmt_max   = '%%-%ds'        % rms_l
    fmt_num   = '%%%ds'         % num_l

    title_l = index_l + 2 + name_l 
    header = (
        ( '{:^%d}' % title_l ).format ( 'Variable' ) ,
        ( '{:^%d}' % type_l  ).format ( 'type'     ) ,
        ( '{:^%d}' % mean_l  ).format ( 'mean'     ) ,
        ( '{:^%d}' % rms_l   ).format ( 'rms'      ) ,
        ( '{:^%d}' % min_l   ).format ( 'min'      ) ,
        ( '{:^%d}' % max_l   ).format ( 'max'      ) ,
        ( '{:^%d}' % num_l   ).format ( '#'        ) )    
               
    table_data = [ header ] 
    for i , v in enumerate ( _vars ) :
        table_data.append ( ( fmt_name  % ( i + 1 , v [ 0 ] ) ,
                              fmt_type  %           v [ 1 ] ,
                              fmt_mean  %           v [ 2 ] ,
                              fmt_rms   %           v [ 3 ] ,
                              fmt_min   %           v [ 4 ] ,
                              fmt_max   %           v [ 5 ] ,
                              fmt_num   %           v [ 6 ] ) )

    if not title :
        
        tt = tree.GetTitle()
        if tt and tt != tree.GetName() : 
            title  = '%s("%s","%s") %d entries,' % ( tree.__class__.__name__ , tree.path , tt , len ( tree ) )
        else :
            title  = '%s("%s") %d entries,'      % ( tree.__class__.__name__ , tree.path ,      len ( tree ) )

        nb = len ( tree.branches () )
        title += '%d branches' % nb 
        nl = len ( tree.leaves   () )
        if nl != nb : title += '%d leaves' % nl
        
        if isinstance ( tree , ROOT.TChain ) :
            nfiles = len ( tree.files() )
            if 1 < nfiles : title += '/%d files ' % nfiles 
        
    import ostap.logger.table as T
    t  = T.table ( table_data , title , prefix = prefix )
    w  = T.table_width ( t )
    return t , w 
    

# ==============================================================================
## get a type of TLeaf object
#  @code
#  tree = ...
#  leaf = t.leaf ( 'QQQ' )
#  print leaf.get_type ( )
#  @endcode 
def _tl_type_ ( leaf ) :
    """Get a type for TLeaf object
    >>> tree = ...
    >>> leaf = t.leaf ( 'QQQ' )
    >>> print leaf.get_type ( )
    """
    
    if not leaf : return 'NULL'
    
    branch   = leaf.GetBranch   ()
    typename = leaf.GetTypeName () 
    
    name     = branch.GetTitle() 
    p1       = name. find ( '[' ) 
    p2       = name.rfind ( ']' )
    if   0 < p1 < p2 :
        typename = '%s [%s]' % ( typename , name [ p1 + 1 : p2 ] )
    elif 0 < p1 :
        typename = '%s [%s]' % ( typename , name [ p1 : ] )

    typename = typename.replace ( 'Float_t'  , 'float'  ) 
    typename = typename.replace ( 'Double_t' , 'double' ) 
    typename = typename.replace ( 'Bool_t'   , 'bool'   )
    
    return typename 


# =============================================================================
_short_types_ = {
    'Char_t'     : 'B' ,
    'UChar_t'    : 'b' ,
    'Short_t'    : 'S' ,
    'UShort_t'   : 's' ,
    'Int_t'      : 'I' , 
    'UInt_t'     : 'i' ,
    'Float_t'    : 'F' ,
    'Float16_t'  : 'f' ,
    'Double_t'   : 'D' ,
    'Double32_t' : 'd' ,
    'Long64_t'   : 'L' ,
    'ULong64_t'  : 'l' ,
    'Bool_t'     : 'O' ,
    ##
    'double'     : 'D' ,
    'float'      : 'F' , 
    'bool'       : 'O' ,
    }
# ==============================================================================
## get a type of TLeaf object
#  @code
#  tree = ...
#  leaf = t.leaf ( 'QQQ' )
#  print leaf.get_short_type ( )
#  @endcode 
def _tl_type_short_ ( leaf ) :
    """Get a type for TLeaf object
    >>> tree = ...
    >>> leaf = t.leaf ( 'QQQ' )
    >>> print leaf.get_type ( )
    """

    ts = leaf.get_type ()
    for k in reversed ( sorted ( _short_types_ ) ) :
        ts = ts.replace ( k , _short_types_[k] )
    return ts 

# ==============================================================================
ROOT.TLeaf . get_type       = _tl_type_
ROOT.TLeaf . get_type_short = _tl_type_short_
ROOT.TLeaf . get_short_type = _tl_type_short_

# ==============================================================================
## print rot-tree in a form of the table
#  @code
#  data = ...
#  print dat.table() 
#  @endcode
def _rt_table_ (  dataset ,  variables = [] ,   cuts = '' , prefix = '' , title = '' , *args ) :
    """print dataset in a form of the table
    >>> dataset = ...
    >>> print dataset.table()
    """
    return _rt_table_0_ ( dataset , variables , cuts , prefix , title , *args )[0]


# =============================================================================
##  print DataSet
def _rt_print2_ ( data  , prefix = '' ) :
    """Print TTree/TChain"""
    
    br = len ( data.branches () ) + len ( data.leaves() )  
    l  = len ( data             )
    if 10000000 < br * l : return _rt_print_ ( data )
    
    from   ostap.utils.basic      import terminal_size, isatty
    if not isatty() : return _rt_table_ ( data )
    th  , tw   = terminal_size()
    rep , wid  = _rt_table_0_ ( data , prefix = prefix ) 
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
## iterator over individual trees in the echain
#  @code
#  chain = ...
#  for tree in  chain.trees ()  :
#       print len(tree)
#  @endcode 
def _rc_itrees_   ( self ) :
    """Iterator over individual trees in the echain
    >>> chain = ...
    >>> for tree in  chain.trees ()  :
    ...     print len(tree)
    """

    _files = self.files()
    for _f in _files :

        c = ROOT.TChain ( self.GetName() )
        c.Add ( _f )
        yield c
        
    
    
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
    
    if isinstance ( index , integer_types ) :
        
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
#  print ( varr )  
#  @endcode 
#  @see numpy.array 
#  @author Albert BURSCHE
#  @date 2015-07-08
def _rt_slice_ ( tree , varname , cut = '' , weight = '' , transpose = False , first = 0 , last = _large ) :
    """ Get ``slice'' from TTree in a form of numpy.array
    ##
    >>> tree = ...
    >>> varr , _  = tree.slice('Pt','eta>3')
    >>> print ( varr )  
    """

    if isinstance ( varname , string_types ) : varname = split_string ( varname , ' ,;:' )
    names = []
    for v in varname :
        names += split_string ( v , ' ,;:' )

    if weight : names.append ( weight )
        
    if not names : return () 
    
    result = []
    
    import numpy    
    for chunk in chunked ( names , 10 ) : ## blocks up to 10 variables 
        
        l    = len ( chunk )
        vars = ':'.join ( chunk )

        ge   = tree.GetEstimate() 
        n    = tree.Draw ( vars , cut , "goff" , last , first )
        n    = tree.GetSelectedRows () 
        if 0 <= ge <= n + 1 :
            tree.SetEstimate ( max ( n + 1 , ge ) )
            n = tree.Draw ( vars , cut , "goff" )
            n = tree.GetSelectedRows () 
        
        for k in range ( l ) :
            result.append ( numpy.array ( numpy.frombuffer ( tree.GetVal ( k ) , count = n ) , copy = True ) )
            
        tree.SetEstimate ( ge ) 

    if not result :
        return None , None 
        
    if weight :
        weights = result[ -1]
        result  = result[:-1]
    else :
        weights = None 

    if not result :
        return result, weights 

    result = numpy.stack ( result )

    if transpose :
        result = result.transpose()
        
    return result, weights

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
def _rt_slices_ ( tree , varnames , cut = '' , weight = '' , transpose = False  , first = 0 , last = _large ) :
    """ Get ``slices'' from TTree in a form of numpy.array
    
    >>> tree = ...
    
    >>> vars1 = tree.slices( ['Pt' , 'eta'] ,'eta>3')
    >>> print vars1
    
    >>> vars2 = tree.slices( 'Pt,eta'  ,'eta>3')
    >>> print vars2
    
    >>> vars3 = tree.slices( 'Pt : eta' ,'eta>3')
    >>> print vars3
    """
    #
    return tree.slice ( varnames , cut , weight , transpose , first , last )


ROOT.TTree .slice  = _rt_slice_
ROOT.TTree .slices = _rt_slices_


# =============================================================================
## get "slices" from TChain in a form of numpy.array
#  @code
#  chain = ...
#  varrs1 = chain.slices ( ['Pt','eta'] , 'eta>3' )
#  print varrs1 
#  varrs2 = chain.slices (  'Pt , eta'  , 'eta>3' )
#  print varrs2
#  varrs3 = chain.slices (  'Pt : eta'  , 'eta>3' )
#  print varrs3
#  @endcode 
#  @see numpy.array 
def _rc_slice_ ( chain , varname , cut = '' , weight = '', transpose = False ) :
    """Get ``slices'' from TChain in a form of numpy.array
    >>> chain = ...
    >>> varrs1 = chain.slices ( ['Pt','eta'] , 'eta>3' )
    >>> print varrs1 
    >>> varrs2 = chain.slices (  'Pt , eta'  , 'eta>3' )
    >>> print varrs2
    >>> varrs3 = chain.slices (  'Pt : eta'  , 'eta>3' )
    >>> print varrs3
    - see numpy.array 
    """

    files = chain.files()
    
    import numpy    

    result , weights = () , () 
    
    for  i , f in enumerate ( files ) :
        
        t = ROOT.TChain( chain.GetName() )
        t.Add ( f )
        
        r , w = _rt_slice_ ( t , varname , cut , weight , transpose )
        if 0 == i :
            result  = r
            weights = w
        else :
            result  = numpy.concatenate ( ( result  , r ) , axis = 0 if transpose else 1 )
            
            if    weight is None or 0 == len ( weights ) : weights = w
            elif  w      is None or 0 == len ( w       ) : pass
            else :
                weights = numpy.concatenate ( ( weights , w ) )
            
        del t
        
    return result , weights


ROOT.TChain .slice  = _rc_slice_
ROOT.TChain .slices = _rc_slice_


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


# ==============================================================================
addbranch_types = string_types + num_types + ( ROOT.TH1 , Ostap.IFuncTree )
# ==============================================================================
## basic types of objects that can be used for <code>add_new_branch</code> methods
#  - string formula
#  - constant number
#  - <code>ROOT.TH1</code> 
#  - <code>ROOT.TH2</code> 
#  - <code>ROOT.TH3</code> 
#  - <code>Ostap.IFuncTree</code>
#  - array-like objects
#  - python callable
#  @see Ostap.IFuncTree
#  @see TH1 
def btypes ( obj ) :
    """Basic types of objects that can be used for `add_new_branch`methods
    - string formula
    - constant number
    - `ROOT.TH1`
    - `ROOT.TH2`
    - `ROOT.TH3`
    - `Ostap.IFuncTree`
    - array-like objects
    - python callable
    """

    if   isinstance ( obj , addbranch_types ) : return True
    elif btypes_array ( obj )                 : return True
    
    return callable ( obj ) 

# =============================================================================
## basic types of array-line objects that can be used for <code>add_new_branch</code> methods 
def btypes_array ( obj ) :
    """Basic types of array-line objects that can be used for <code>add_new_branch</code> methods 
    """
    
    if   isinstance ( obj, string_types )          : return False

    ## efficient treatment for ROOT versions from  6/24 
    elif isinstance ( obj , array.array ) and  (6,24)<=root_info and \
         obj.typecode in ( 'f' , 'd' ,'i' , 'l' )  : return True 

    ## efficient treatment for ROOT verisons from 6/24 
    elif numpy and (6,24)<= root_info       and \
         isinstance ( obj , numpy.ndarray ) and \
         obj.dtype  in ( numpy.float32 ,
                         numpy.float64 ,
                         numpy.int32   ,
                         numpy.int64   )           : return True
    
    ## generic case with array-like structure
    elif isinstance ( obj , sized_types    ) and \
         isinstance ( obj , sequence_types ) and \
         hasattr    ( obj , '__getitem__'  ) : return True

    return False

    
# ==============================================================================
## add new branch to the chain
#  @see Ostap::Trees::add_branch
#  @see Ostap::IFuncTree   
def _chain_add_new_branch ( chain , name , function , verbose = True , value = 0 ) :
    """ Add new branch to the tree
    - see Ostap::Trees::add_branch
    - see Ostap::IFuncTree 
    """
    assert isinstance ( chain , ROOT.TChain ), 'Invalid chain!'

    if len ( chain.files() ) <= 1 :
        return add_new_branch ( chain               ,
                                name     = name     ,
                                function = fnuction , 
                                verbose  = verbose  ,
                                value    = value    ) 
    
    if isinstance ( function , dictlike_types ) :
        assert name     is None , 'add_branch: when function is dict, name must be None!'
        name , function = function , None 
        
    names = name
    if isinstance ( names , string_types )  : names =  [ names ]    
    for n in names : 
        assert not n in chain.branches() ,'Branch %s already exists!' % n 
        
    assert ( isinstance ( name , dictlike_types ) and function is None ) or btypes ( function ) ,\
           "add_branch: invalid type of ``function'': %s/%s" % ( function , type ( function ) )  

    if   isinstance ( name  ,  dictlike_types ) and function is None : pass    
    elif isinstance ( function , addbranch_types )                   : pass 
    elif btypes_array ( function ) :    
        return _chain_add_new_branch_array ( chain                ,
                                             name      = name     ,
                                             the_array = function ,
                                             verbose   = verbose  ,
                                             value     = value    ) 

    files = chain.files   ()
    cname = chain.GetName () 
    

    tree_verbose  = verbose and      len ( files ) < 5
    chain_verbose = verbose and 5 <= len ( files )

    import ostap.io.root_file
    for fname in progress_bar ( files , len ( files ) , silent = not chain_verbose ) :
        
        logger.debug ('Add_new_branch: processing file %s' % fname )
        with ROOT.TFile.Open  ( fname , 'UPDATE' , exception = True ) as rfile :
            ## get the tree 
            ttree = rfile.Get ( cname )
            ## treat the tree 
            add_new_branch    ( ttree                   ,
                                name     = name         ,
                                function = function     ,
                                verbose  = tree_verbose ,
                                value    = value        )
            
    ## recollect the chain 
    newc = ROOT.TChain ( cname )
    for f in files : newc.Add ( f  )
    
    return newc 

# ==============================================================================
## add new branch to the chain
#  @see Ostap::Trees::add_branch
#  @see Ostap::IFuncTree   
def _chain_add_new_branch_array ( chain           ,
                                  name            ,
                                  the_array       ,
                                  verbose = True  ,
                                  value   = 0     ) : 
    """ Add new branch to the tree
    - see Ostap::Trees::add_branch
    - see Ostap::IFuncTree 
    """
    assert isinstance ( chain , ROOT.TChain ), 'Invalid chain!'

    names = name
    if isinstance ( names , string_types )  : names =  [ names ]    
    for n in names : 
        assert not n in chain.branches() ,'Branch %s already exists!' % n 

    assert isinstance ( the_array , sized_types    ) and \
           isinstance ( the_array , sequence_types ) and \
           hasattr    ( the_array , '__getitem__'  ) ,   \
           "Invalid type of ``the_array'' %s/%s" % ( the_array , type ( the_array ) ) 
    
    files = chain.files   ()
    cname = chain.GetName () 
    
    from ostap.utils.progress_bar import progress_bar

    verbose = verbose and 1 < len ( files )
    
    import ostap.io.root_file

    start = 0
    
    tree_verbose  = verbose and      len ( files ) < 5
    chain_verbose = verbose and 5 <= len ( files )

    for i , fname in enumerate ( progress_bar ( files , len ( files ) , silent = not chain_verbose  ) ) :
        
        logger.debug ('Add_new_branch: processing file %s' % fname )
        with ROOT.TFile.Open  ( fname , 'UPDATE' , exception = True ) as rfile :
            ## get the tree 
            ttree = rfile.Get ( cname )
            ## treat the tree
            size  = len ( ttree )
            end   = min ( start + size , len ( the_array ) ) 
                          
            if   0 == i       : what = the_array
            elif end <= start : what = ()
            else              : what = the_array [ start : end ]
            
            add_new_branch    ( ttree                   ,
                                name     = name         ,
                                function = what         ,
                                verbose  = tree_verbose ,
                                value    = value        ) 
            start += size

    ## recollect the chain 
    newc = ROOT.TChain ( cname )
    for f in files : newc.Add ( f  )
    
    return newc 


# =============================================================================
try : 
    import numpy, ctypes  
except ImportError :
    numpy = None
    
# ==============================================================================
## Add new branch to the tree
# 
#   - Using formula:
#   @code
#   >>> tree = ....
#   >>> tree.add_new_branch ( 'pt2' , 'pt*pt' ) ## use formula
#   @endcode
# 
#   - Sampling from 1D-histogram
#   @code 
#   >>> tree = ...
#   >>> h1   = ...  ## 1D histogram to be sampled 
#   >>> tree.add_new_branch ( 'ntracks' , h1 )
#   @endcode
# 
#   - Sampling from 2D-histogram
#   @code
#   >>> tree = ...
#   >>> h2   = ...  ## 2D histogram to be sampled 
#   >>> tree.add_new_branch ( [ 'pt', 'eta' ] ,  h2 )
#   @endcode
#
#   - Sampling from 3D-histogram
#   @code
#   >>> tree = ...
#   >>> h3   = ...  ## 3D histogram to be sampled 
#   >>> tree.add_new_branch ( [ 'pt', 'eta' , 'ntracks' ] ,  h3 )
#   @endcode
# 
#   - adding 1D-histogram as function:
#     for  each entry it gets the value of `mLb` variable
#     and stores the value from the histogram in new `S_sw` variable:
#   @code
#   >>> h  = ...  ## histogram, e.g. sPlot 
#   >>> fn = Ostap.Functions.FuncTH1 ( sph , 'mLb' )
#   >>> tree.add_new_branch ( 'S_sw' , fn )
#   @endcode
# 
#   - arbitrary function derived from Ostap.ITreeFunc:
#   @code
#   >>> fun = ...
#   >>> tree.add_new_branch ( 'Var' , fun )
#   @endcode
# 
#   -  Several variables can be filled at once:
#   @code
#   >>> tree = ....
#   >>> tree.add_new_branch ( { 'pt2' : 'pt*pt'  ,
#   ...                         'et2' : 'pt*pt+mass*mass' } , None ) ## use formulas
#   @endcode
#   @attention it makes a try to reopen the file with tree in UPDATE mode,
#              and it fails when it is not possible!
#
#  @see Ostap::Trees::add_branch
#  @see Ostap::IFuncTree 
def add_new_branch ( tree , name , function , verbose = True , value = 0 ) :
    """ Add new branch to the tree

    - Using formula:
    
    >>> tree = ....
    >>> tree.add_new_branch ( 'pt2' , 'pt*pt' ) ## use formula
    
    - Sampling from 1D-histogram
    
    >>> tree = ...
    >>> h1   = ...  ## 1D histogram to be sampled 
    >>> tree.add_new_branch ( 'ntracks' , h1 )
    
    - Sampling from 2D-histogram
    
    >>> tree = ...
    >>> h2   = ...  ## 2D histogram to be sampled 
    >>> tree.add_new_branch ( [ 'pt', 'eta' ] ,  h2 )
    
    - Sampling from 3D-histogram
    
    >>> tree = ...
    >>> h3   = ...  ## 3D histogram to be sampled 
    >>> tree.add_new_branch ( [ 'pt', 'eta' , 'ntracks' ] ,  h3 )

    - adding histogram as function: for  each entry it gets the value of `mLb` variable
    and stores the value from the histogram in new `S_sw` variable:
    
    >>> h  = ...  ## histogram, e.g. sPlot 
    >>> fn = Ostap.Functions.FuncTH1 ( sph , 'mLb' )
    >>> tree.add_new_branch ( 'S_sw' , fn )
    
    - arbitrary function derived from Ostap.ITreeFunc:
    
    >>> fun = ...
    >>> tree.add_new_branch ( 'Var' , fun )
    
    -  Several variables can be filled at once:
    >>> tree = ....
    >>> tree.add_new_branch ( { 'pt2' : 'pt*pt'  ,
    ...                         'et2' : 'pt*pt+mass*mass' } , None ) ## use formulas
    
    
    - ATTENTION: it makes a try to reopen the file with tree in UPDATE mode,
    and it fails when it is not possible!
    
    - see Ostap::Trees::add_branch
    - see Ostap::IFuncTree
    
    """
    if not tree :
        logger.error (  "add_branch: Invalid Tree!" )
        return
    elif isinstance ( tree  , ROOT.TChain ) and 1 < len ( tree.files() ) :
        return _chain_add_new_branch ( tree                ,
                                       name                ,
                                       function = function ,
                                       verbose  = verbose  ,
                                       value    = value    )
    
    if isinstance ( function , dictlike_types ) :
        assert name     is None , 'add_branch: when function is dict, name must be None!'
        name , function = function , None 

    names = name 
    if isinstance ( names , string_types ) : names = [ names ]
    names = [ n.strip() for n in names ] 
    for n in names : 
        assert not n in tree.branches() ,"``Branch'' %s already exists!" % n


    assert ( isinstance ( name , dictlike_types ) and function is None ) or btypes ( function ) ,\
           "add_branch: invalid type of ``function'': %s/%s" % ( function , type ( function ) )  


    if isinstance ( name  ,  dictlike_types ) and function is None :
        
        typeformula = False 
        for k in  name.keys() :
            
            assert not k in tree.branches() ,'Branch %s already exists!' % k
            v = name [ k ]
            if   isinstance ( v , string_types    ) : pass
            elif isinstance ( v , Ostap.IFuncTree ) : typeformula = True
            else : raise TypeError ('add_branch: Unknown brnach %s/%s for %s'  % ( v , type( v ) , k ) )
                    
        if typeformula : MMAP = Ostap.Trees.FUNCTREEMAP
        else           : MMAP = std.map  ( 'std::string'       , 'std::string' )
        
        mmap = MMAP () 
        for k in name.keys() :
            v = name [ k ]
            if typeformula and isinstance ( v , string_types ) :
                v = Ostap.Functions.FuncFormula ( v , tree )
                funcs.append ( v ) 
            ## mmap.insert ( PAIR ( k , v ) )
            mmap[ k ] = v 
            
        args = mmap ,
        
    elif isinstance ( function , addbranch_types ) :

        args = tuple ( [n  for n in names ] + [ function ] )


    ## efficient case with array 
    elif ( 6 , 24 ) <= root_info                   and \
             isinstance ( function , array.array ) and \
             function.typecode in ( 'f' , 'd' , 'h' , 'i' , 'l' , 'H' , 'I' , 'L' ) :
        
        data = function 
        args = tuple ( [ n  for n in names ] + [ data , len ( data ) , value ] )
        
    ## efficient case with array 
    elif numpy and (6,24) <= root_info            and \
         isinstance ( function , numpy.ndarray )  and \
         function.dtype  in ( numpy.float32 ,
                              numpy.float64 ,
                              numpy.int16   ,
                              numpy.int32   ,
                              numpy.int64   ,
                              numpy.uint16  ,
                              numpy.uint32  , 
                              numpy.uint64  ) :        
        data = function 
        dt   = data.dtype 
        ct   = numpy.ctypeslib._ctype_from_dtype( dt )
        buff = data.ctypes.data_as ( ctypes.POINTER( ct ) )  
        args = tuple ( [ n  for n in names ] + [ buff , len ( data ) , value ] )
        
    ## generic case with array-like structure 
    elif isinstance ( function , sized_types    ) and \
         isinstance ( function , sequence_types ) and \
         hasattr    ( function , '__getitem__'  ) :

        data = function 
        from ostap.trees.funcs import PyTreeArray as PTA
        args = tuple ( [n  for n in names ] + [ PTA ( data , value = value ) ] )

    elif callable ( function )  :
        
        from ostap.trees.funcs import PyTreeFunction as PTF
        args = tuple ( [n  for n in names ] + [ PTF ( function  ) ] )

    else :

        logger.warning ('addbranch: suspicion case name/function:  %s/%s %s/%s' % (
            name , type(name) , function, type(function) ) ) 
                        
        args = tuple ( [n  for n in names ] + [ function ] )
        
    tname = tree.GetName      ()
    tdir  = tree.GetDirectory ()
    tpath = tree.path

    from ostap.io.root_file        import REOPEN 
    from ostap.utils.progress_conf import progress_conf
    with ROOTCWD() , REOPEN ( tdir ) as tfile :
        
        tfile.cd() 
        ttree = tfile.Get ( tpath )

        ## add progress bar 
        if verbose : sc    = Ostap.Trees.add_branch ( ttree , progress_conf , *args )
        else       : sc    = Ostap.Trees.add_branch ( ttree ,                 *args )
        
        if   sc.isFailure     () : logger.error ( "Error from Ostap::Trees::add_branch %s" % sc )
        elif tfile.IsWritable () :
            tfile.Write( "" , ROOT.TObject.kOverwrite )
            logger.debug ('Write back TTree %s to %s' % ( tpath , tfile ) )            
        else :
            logger.error ("Can't write TTree %s back to the file %s" % ( tpath , tfile ) )

ROOT.TTree.add_new_branch = add_new_branch 

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
def add_reweighting ( tree , weighter , name = 'weight' , verbose = False ) :
    """Add specific re-weighting information into ROOT.TTree
    
    >>> w    = Weight ( ... ) ## weighting object ostap.tools.reweight.Weight 
    >>> data = ...
    >>> data.add_reweighting ( w )
    - see ostap.tools.reweight
    - see ostap.tools.reweight.Weight 
    - see ostap.tools.reweight.W2Tree 
    """
    
    import ostap.tools.reweight as W
    
    assert isinstance ( weighter , W.Weight ), "Invalid type of ``weighting''!"
    
    ## create the weighting function 
    wfun = W.W2Tree ( weighter )
    
    return tree.add_new_branch (  name , wfun , verbose = verbose ) 

ROOT.TTree.add_reweighting = add_reweighting
    

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
## get all variables needed to evaluate the expressions for the given tree
#  @code
#  tree = 
#  vars = tree.the_variables ( [ 'x>0&& y<13' , 'zzz*15' ] )   
#  vars = the_variables ( tree , [ 'x>0&& y<13' , 'zzz*15' ] ) ## ditto
#  @endcode 
def the_variables ( tree , expression , *args ) :
    """Get all variables needed to evaluate the expressions for the given tree
    >>> tree = 
    >>> vars = tree.the_variables ( tree , [ 'x>0&& y<13' , 'zzz' ]  )
    >>> vars =      the_variables (        [ 'x>0&& y<13' , 'zzz' ]  ) ##  ditto
    """
    from ostap.core.core import fID
    
    if isinstance  ( expression, ( list , tuple ) ) :
        exprs = list ( expression ) 
    else :
        exprs = [ expression ]
        
    for e in args :
        if isinstance  ( e , ( list , tuple ) ) :
            exprs = exprs + list ( e ) 
        else :
            exprs.append ( e )

    vars = set() 
    for e in exprs :

        if not e : continue
        
        tf = Ostap.Formula ( fID() , str ( e ) , tree )
        if not tf.ok()  :
            logger.error ('the_variables: Invalid formula "%s"' % e )
            del tf 
            return None
        
        i    =  0
        leaf = tf.GetLeaf ( i )
        while leaf :
            lname = leaf.GetName()
            vars.add ( lname )
            i += 1
            leaf = tf.GetLeaf ( i )                

        del tf

    vvars    = list ( vars )
    
    leaves   = tree.leaves   ()
    branches = tree.branches ()
    all      = set ( leaves + branches )

    ## If variable-sized vectors, add the lengths...

    sizes = set()  
    for v in vvars :

        l  = tree.GetLeaf ( v )
        t  = l.GetTitle ()
        p1 = t. find ( '[' )
        p2 = t.rfind ( ']' )
        if 0 < p1 < p2 :
            n = t [ p1 + 1 : p2 ]
            n = n.replace ( '][' , ',' )
            n = n.split   ( ',' ) 
            for i in n :
                if i in all : sizes.add ( i )
            
        b  = l   .GetBranch ( )
        t  = b   .GetTitle  ( )
        p1 = t. find ( '[' )
        p2 = t.rfind ( ']' )
        if 0 < p1 < p2 :
            n = t [ p1 + 1 : p2 ]
            n = n.replace ( '][' , ',' )
            n = n.split   ( ',' ) 
            for i in n :
                if i in all : sizes.add ( i )
                
    vars  = [ v for v in vars if not v in sizes ]
    vars.sort ()
    
    sizes = list ( sizes )
    sizes.sort ()

    ## SIZES MUST BE FIRST! 
    return  tuple ( sizes ) + tuple ( vars )  ## SIZES MUST BE FIRST! 


ROOT.TTree.the_variables = the_variables

# ===============================================================================
## Get all ``size''-variables
#  @code
#  tree  = ... 
#  sizes = tree.size_vars() 
#  @endcode 
def _rt_size_vars_ ( tree ) :
    """Get all ``size''-variables
    >>> tree  = ... 
    >>> sizes = tree.size_vars() 
    """
    leaves   = tree.leaves   ()
    branches = tree.branches ()
    all      = set ( leaves + branches )
    sizes    = set () 

    for leaf in leaves :
        
        l = tree.GetLeaf( leaf )
        t = l.GetTitle ()

        p1 = t. find ( '[' )
        p2 = t.rfind ( ']' )
        if 0 < p1 < p2 :
            n = t [ p1 + 1 : p2 ]
            n = n.replace ( '][' , ',' )
            n = n.split   ( ',' ) 
            for i in n :
                if i in all : sizes.add ( i )
            
        b  = l   .GetBranch ( )
        t  = b   .GetTitle  ( )
        p1 = t. find ( '[' )
        p2 = t.rfind ( ']' )
        if 0 < p1 < p2 :
            n = t [ p1 + 1 : p2 ]
            n = n.replace ( '][' , ',' )
            n = n.split   ( ',' ) 
            for i in n :
                if i in all : sizes.add ( i )

    return tuple ( sorted ( sizes ) )

ROOT.TTree.size_vars = _rt_size_vars_

# ===============================================================================
## Get all ``array-like''-variables
#  @code
#  tree  = ... 
#  arrays = tree.array_vars() 
#  @endcode 
def _rt_array_vars_ ( tree ) :
    """Get all ``array-like''-variables
    >>> tree  = ... 
    >>> arrays = tree.array_vars() 
    """
    leaves   = tree.leaves   ()
    branches = tree.branches ()
    all      = set ( leaves + branches )
    
    arrays   = set () 

    for leaf in leaves :
        
        l = tree.GetLeaf( leaf )
        t = l.GetTitle ()

        p1 = t. find ( '[' )
        p2 = t.rfind ( ']' )
        if 0 < p1 < p2 :
            arrays.add ( leaf )
            
        b  = l   .GetBranch ( )
        t  = b   .GetTitle  ( )
        p1 = t. find ( '[' )
        p2 = t.rfind ( ']' )
        if 0 < p1 < p2 :
            arrays.add ( b.GetName()  )

    return tuple ( sorted ( arrays ) )

ROOT.TTree.array_vars = _rt_array_vars_

# ===============================================================================
## @class ActiveBranches
#  Context manager to activate only certain branches in the tree.
#  It drastically speeds up the iteration over the tree.
#  @code
#  tree = ...
#  with ActiveBraches( tree , '*_Lb' , 'eta_Lc') :
#    for i in range(1000000) :
#        tree.GetEntry ( i )
#        print tree.pt_Lb, tree.eta_Lc 
#  @endcode
class ActiveBranches(object) :
    """Context manager to activate only certain branches in the tree.
    - It drastically speeds up the iteration over the tree.
    >>> tree = ...
    >>> with ActiveBraches( tree , '*_Lb', 'eta_Lc') :
    ...    for i in range(1000000) :
    ...    tree.GetEntry ( i )
    ...    print tree.pt_Lb, tree.eta_Lc
    """
    def __init__ ( self , tree , *vars ) :
        
        assert tree and vars , 'ActiveBrnaches: both tree and vars must be valid!'
        
        self.__tree = tree
        self.__vars = vars 
        
    ## context manager: ENTER 
    def __enter__ ( self ) :
        ## deactivate the all branches 
        self.__tree.SetBranchStatus ( '*' , 0 )     ## deactivate *ALL* branches
        for var in self.__vars :
            ##  activate certain branches 
            self.__tree.SetBranchStatus ( var , 1 )  ## activate only certain branches
            
        return self.__tree 
    
    ## context manager: EXIT 
    def __exit__ ( self , *_ ) :
        ## reactivate all branches again 
        self.__tree.SetBranchStatus ( '*' , 1 )     ## reactivate *ALL* branches
        
# ===============================================================================
## Context manager to activate only certain branches in the tree.
#  It drastically speeds up the iteration over the tree.
#  @code
#  tree = ...
#  with active_braches ( tree , '*_Lb' , 'eta_Lc') :
#    for i in range(1000000) :
#        tree.GetEntry ( i )
#        print tree.pt_Lb, tree.eta_Lc 
#  @endcode
def active_branches ( tree , *vars ) :
    """Context manager to activate only certain branches in the tree.
    - It drastically speeds up the iteration over the tree.
    >>> tree = ...
    >>> with active__braches ( tree , '*_Lb', 'eta_Lc') :
    ...    for i in range(1000000) :
    ...    tree.GetEntry ( i )
    ...    print tree.pt_Lb, tree.eta_Lc
    """
    return ActiveBranches ( tree , *vars ) 
    

# =============================================================================
## files and utilisties for TTree/TChain "serialization"
# =============================================================================
## get some file info for the given path 
def file_info ( fname ) :
    """Get some file info for the given path
    """
    p , s , f = fname.partition ( '://' )
    if p and s : return 'Protocol'
    if os.path.exists ( fname ) and os.path.isfile ( fname ) and os.access ( fname , os.R_OK ) :
        s = os.stat ( fname )
        return s.st_mode , s.st_size , s.st_uid, s.st_gid, s.st_atime , s.st_mtime , s.st_ctime
    return 'Invalid'
# =============================================================================
from ostap.utils.cleanup  import CleanUp
# =============================================================================
## @class Chain
#  simple class to keep pickable definition of tree/chain
#  it is needed for multiprcessing 
class Chain(CleanUp) :
    """Simple class to keep definition of tree/chain ``pickable''
    """
    def __getstate__  ( self ) :

        def fullfn ( f ) :
            if os.path.exists ( f ) and os.path.isfile ( f ) :
                ff = os.path.abspath ( f )
                if os.path.samefile ( ff , f ) : return ff
            return f
        
        file_infos = tuple ( ( fullfn ( f ) , file_info ( f ) ) for f in self.__files ) 
        
        return { 'name'     : self.__name    ,
                 'first'    : self.__first   ,            
                 'nevents'  : self.__nevents ,
                 'files'    : file_infos     ,
                 'host'     : self.__host    }

    def __setstate__  ( self , state ) :
        
        self.__name    = state [ 'name'    ]
        self.__first   = state [ 'first'   ]
        self.__nevents = state [ 'nevents' ]
        #
        origin         = state [ 'host'    ]
        file_infos     = state [ 'files'   ]
        #
        ##
        import socket 
        self.__host    = socket.getfqdn ().lower()
        ##
        same_host      = origin == self.__host
        files_         = []
        self.__lens    = ()
        self.__chain   = None
        
        for fname , finfo in file_infos :
            if   same_host : files_.append ( fname )
            else :
                fnew = file_info ( fname )
                if fnew == finfo :
                    ## hosts are different but the files are the same (shared file system?)
                    files_.append ( fname )
                else :
                    # =========================================================
                    ## the file need to be copied locally
                    full_name  = '%s:%s' % ( origin , fname ) 
                    copied , t = scp_copy   ( full_name )
                    if copied :
                        cinfo = file_info ( copied )
                        c = '%s -> %s' % ( full_name , "%s:%s" % ( self.__host , copied ) )
                        if cinfo[:4] == finfo [:4] :
                            files_.append ( copied )
                            size = cinfo[1]
                            s =  cinfo [ 1 ] / float ( 1024 ) / 1025 ##  MB 
                            v = s / t                  ## MB/s 
                            logger.debug ( 'File copied %s :  %.3f[MB] %.2f[s] %.3f[MB/s]' % ( c , s , t , v ) )
                        else :
                            logger.error ( 'Something wrong with the copy %s : %s vs %s ' % ( c , cinfo[:4] , finfo[:4] ) )
                    else : logger.error ("Cannot copy the file %s"  % full_name )
                        
        self.__files = tuple ( files_ )
        ## self.__files = tuple ( [ i   for i,j in file_infos ] )
        
        ## reconstruct the chain
        import ROOT
        import ostap.io.root_file

        ## for f in  self.__files :
        ##     rf = ROOT.TFile.Open ( f , 'read' , exception = False )
        ##     if rf and not rf.IsZombie()  :
        ##         logger.verbose ('Chain/setstate:file %s is OK ' % f ) 
        ##     else :
        ##         logger.error   ('Chain/settate: file %s is not OK ' % f ) 
        ##     del rf
            
        ## self.__chain   = ROOT.TChain ( self.__name  )
        ## for f in self.__files  : self.__chain.Add ( f )
                
    # ======================================================================================
    ## create Chain object:
    #  - either from the real TTree/TChain:
    #  @code
    #  chain = ...  ## ROOT.TChain
    #  ch = Chain ( chain ) 
    #  @endcode
    #  - or from description:
    #  @code
    #  ch = Chain ( name = 'n', files = [ ... ] )
    #  @endcode 
    def __init__ ( self , tree = None , name = None , files = [] , first = 0 , nevents = -1  ) :
        """ Create Chain object 
        
        - either from the real TTree/TChain:
        
        >>> chain = ...  ## ROOT.TChain
        >>> ch = Chain ( chain )
        
        - or from description:
        
        >>> ch = Chain ( name = 'n', files = [ ... ] )
        """
        assert   isinstance ( tree , Chain       ) or \
               ( name and files                  ) or \
               ( isinstance ( tree , ROOT.TTree  ) and valid_pointer ( tree  ) )   ,\
               "Invalid tree/name/files combination: %s/%s%s" % ( tree , name , files    )
        
        assert isinstance ( first , int ) and  0 <= first     , \
               "Invalid ``first'' %s/%s"                      % ( first , type ( first ) ) 
        
        self.__first   = int ( first )  
        self.__nevents = nevents if 0 <= nevents < _large else -1 
        self.__chain   = None
        self.__name    = 'Unknown!'
        
        import socket 
        self.__host    = socket.getfqdn ().lower()

        ## copy-like, ignore other arguments  
        if isinstance  ( tree , Chain ) : 
            name    = tree.name 
            files   = tree.files
            first   = tree.first
            nevents = tree.nevents 
            tree    = tree.chain

        if files and isinstance ( files , str ) : files = files,
                         
        if name and files :

            self.__name  = name
            self.__files = files

            if isinstance ( tree , ROOT.TTree ) and valid_pointer ( tree ) : 
                chain = self.__create_chain() 
                assert valid_pointer ( chain ), 'Invalid TChain!'
                assert len ( files ) == len ( chain.files() ) , 'Invalid length of files'
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
                    fname  = CleanUp.tempfile ( suffix = '.root' , prefix = 'ostap-tree-' )
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
                assert chain.GetEntries() == tree.GetEntries () , 'Something wrong happens here :-( '
                self.__chain = chain

        # =====================================================================
        ## The final adjustment
        self.__lens = () 
        if 0 < self.__first or 0 < nevents < _large :
            
            _first = self.__first
            _files = []
            total  = 0 
            for f in self.__files :
                
                t = ROOT.TChain ( self.name )
                t.Add ( f )
                clen = t.GetEntries() 

                if _first < clen :
                    
                    _files.append (  ( f , clen ) )
                    total += clen
                    
                else             :                    
                    _first -= clen
                    continue 

                if 0 < nevents and _first + nevents < total :
                    break 
                
            ## redefine quantities:
            self.__files = tuple ( ( f [0] for f in _files ) )
            self.__lens  = tuple ( ( f [1] for f in _files ) )
            self.__first = _first
            
            if 0 < nevents and _first + nevents < total : pass
            else                                        : nevents = -1 
            
            self.__nevents = nevents 

    # =========================================================================
    ## get/calculate the lengths of the individual trees 
    def calc_lens ( self ) :
        """Get/calculate the lengths of the individual trees
        """
        if self.__lens : return self.__lens

        lens = []
        for f in self.__files :
            t = ROOT.TChain ( self.name )
            t.Add ( f )
            lens.append ( t.GetEntries () )
            
        self.__lens = tuple ( lens )
        return self.__lens
        
    ## split the chain for several chains  with at most chunk_size entries
    def slow_split ( self , chunk_size = 200000 ) :
        """ Split the chain/tree for several chains/trees with at most chunk_size entries
        >>> tree = ....
        >>> trees = tree.split ( chunk_size = 1000000 ) 
        """

        if chunk_size <= 0 : chunk_size = _large 
        
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
        if chunk_size <= 0 : chunk_size = _large
        if max_files  <= 0 : max_files  = 1 
        
        if 0 != self.first or 0 < self.__nevents :
            return self.slow_split ( chunk_size )
        
        ## first split on per-file basis
        fs     = self.files
        chains = [ Chain ( tree = True , name = self.name , files = fs[c] ) for c in self.get_slices ( 0 , len ( fs ) , max_files ) ]
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

        if not self.__lens : self.calc_lens()
            
        if self.__lens  :
            total  = sum ( self.__lens )
            len1   = total - self.__first 
            return len1 if self.__nevents < 0 else min ( len1 , self.__nevents )
        
        if self.__chain is None : self.__chain = self.__create_chain () 
        return len ( self.__chain )
    
    def __create_chain ( self ) :
        """'chain' : get the underlying tree/chain"""
        import ROOT 
        c = ROOT.TChain ( self.name )
        for f in self.__files  : c.Add ( f )
        return c

    @property
    def chain ( self ) :
        """``chain'' : get the underlying tree/chain"""
        if self.__chain is None :
            self.__chain = self.__create_chain ()
            
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

    # ============================================================================
    ## get DataFrame
    #  @code
    #  tree = ....
    #  f  =  tree.frame () ## get
    #  f1 =  tree.frame ('px', 'py' , 'pz') ## get frame with default branches
    #  f2 =  tree.frame ( tx = 'px/pz' , ty = 'py/pz') ## define new variables
    #  @endcode 
    def  frame ( self , *vars , **newvars ) :
        """``frame'' : get ROOT.RDataFrame for the given chain/tree
        >>> tree = ....
        >>> f  = tree.frame () ## get
        >>> f1 = tree.frame ('px', 'py' , 'pz') ## get frame with default branches
        >>> f2 = tree.frame ( tx = 'px/pz' , ty = 'py/pz') ## define new variables 
        """
        from ostap.frames.frames import DataFrame
        from ostap.core.core     import strings 
        fnames = strings  ( *self.files )
        vnames = strings  ( *vars       )
        df = DataFrame  ( self.name , fnames , vnames )
        for k in new_vars : df =  df.Define ( k , new_vars [k] )
        return  df                              
    
    def __str__ ( self ) :
        r = "Chain('%s',%s" % ( self.name , self.__files )
        if 0 != self.first or 0 <= self.__nevents :
            r += ",%s,%s" % ( self.first , self.__nevents )            
        return r + ")"
    __repr__ = __str__

    # =========================================================================
    ## delegate all other attributes to the underlying chain object 
    def __getattr__  ( self , attr ) :
        """Delegate all other attributes to the underlying chain object"""
        return getattr  ( self.chain , attr )

    # =========================================================================
    ## add/merge two chains 
    def __add__ ( self , other ) :
        """Add/merge two chains
        """
        if  0 != self.first   : return NotImplemented 
        if -1 != self.nevents : return NotImplemented

        if   isinstance ( other , Chain ) and self.name == other.name and\
               0 == other.first and -1 == other.nevents :
            
            files1 = set  ( self.files  )
            files2 = set  ( other.files ) 
            files  = list ( files1 | files2 )
            files.sort()
            
            return Chain ( name = self.name , files = files )
        
        elif isinstance ( other , ROOT.TChain ) and self.name != other.name :
            
            files1 = set  ( self.files      )
            files2 = set  ( other.files()   )  
            files  = list ( files1 | files2 )
            files.sort ()
            
            return Chain ( name = self.name , files = files )
        
        return NotImplemented
        
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
        last = min ( ll , self.first + self.nevents if 0 <= self.nevents else _large ) 

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


    # =========================================================================
    ## delegate all other attributes to the underlying chain object 
    def __getattr__  ( self , attr ) :
        """Delegate all other attributes to the underlying chain object"""
        return getattr  ( self.chain , attr )
    
    @property
    def tree ( self ) :
        """``tree'' : get the underlying tree/chain"""
        return self.chain

# =============================================================================
_decorated_classes_ = (
    ROOT.TTree  ,
    ROOT.TChain ,   
    ROOT.TLeaf      
    )
_new_methods_       = (
    #
    ROOT.TTree .withCuts  ,
    ROOT.TChain.withCuts  ,
    ROOT.TTree. __len__   ,
    #
    ROOT.TTree. rows      ,
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
    ROOT.TTree.size_vars    , 
    ROOT.TTree.array_vars   , 
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
    #
    ROOT.TTree.the_variables    ,
    ROOT.TTree.add_new_branch   ,
    ROOT.TTree.add_reweighting  ,
    ##
    ROOT.TLeaf.get_type         ,
    ROOT.TLeaf.get_type_short   ,
    ROOT.TLeaf.get_short_type   ,
    )
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
